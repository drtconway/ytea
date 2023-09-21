use noodles::{
    bam::{self, IndexedReader},
    bgzf::Reader,
    core::{region::ParseError, Region},
    sam::Header,
};
use rand::Rng;
use rayon::prelude::{IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::sync::mpsc::channel;
use std::{
    cmp::min,
    io::{Error, ErrorKind, Result},
};
use threadpool::ThreadPool;

use crate::summarise::Summariser;

fn do_all_calculate_coverage_2(bam: &str, loci: &[(String, usize)]) -> Result<BasicInfo> {
    let mut tasks: Vec<(String, String, usize)> = Vec::new();
    for x in loci {
        tasks.push((bam.to_string(), x.0.to_string(), x.1));
    }

    let pool = ThreadPool::new(4);
    let (tx, rx) = channel();
    for task in tasks {
        let bam: String = task.0.to_string();
        let locus: String = task.1.to_string();
        let pos: usize = task.2;
        let tx = tx.clone();
        pool.execute(move || {
            let res = do_calculate_coverage(&bam, &locus, pos);
            tx.send(res).expect("send failed");
        });
    }

    let mut cov_stats = Summariser::new();
    let mut insert_stats = Summariser::new();
    let mut read_length = 0;

    for res in rx.iter().take(loci.len()) {
        if let Some((x_cov, x_stats, x_read_length)) = res {
            if x_stats.n > 0 {
                cov_stats.add(x_cov as f64);
                insert_stats.add_other(&x_stats);
            }
            if x_read_length > read_length {
                read_length = x_read_length;
            }
        }
    }
    println!("sampled {} positions successfully", cov_stats.n);
    let res = BasicInfo {
        bam: bam.to_string(),
        cov: cov_stats.mean(),
        read_length: read_length,
        insert_size_mean: insert_stats.mean(),
        insert_size_sd: insert_stats.sd(),
    };
    println!("{}", serde_json::to_string(&res).unwrap());
    Ok(res)
}

fn do_calculate_coverage(bam: &str, locus: &str, pos: usize) -> Option<(usize, Summariser, usize)> {
    let mut reader = bam::indexed_reader::Builder::default()
        .build_from_path(bam)
        .expect("reader failed");
    let header = reader.read_header().expect("read header failed");
    let res = calculate_coverage(&mut reader, &header, locus, pos).expect("coverage failed");
    res
}

fn do_all_calculate_coverage(bam: &str, loci: &[(String, usize)]) -> Result<BasicInfo> {
    let mut tasks: Vec<(String, String, usize)> = Vec::new();
    for x in loci {
        tasks.push((bam.to_string(), x.0.to_string(), x.1));
    }
    let mut cov_stats = Summariser::new();
    let mut insert_stats = Summariser::new();
    let mut read_length: usize = 0;
    let results: Vec<Option<(usize, Summariser, usize)>> = tasks
        .par_iter()
        .map(|x| do_calculate_coverage(&x.0, &x.1, x.2))
        .collect();
    for res in results.iter() {
        match res {
            Some((x_cov, x_stats, x_read_length)) => {
                if x_stats.n > 0 {
                    cov_stats.add(*x_cov as f64);
                    insert_stats.add_other(&x_stats);
                }
                if *x_read_length > read_length {
                    read_length = *x_read_length;
                }
            }
            None => {}
        }
    }
    let res = BasicInfo {
        bam: bam.to_string(),
        cov: cov_stats.mean(),
        read_length: read_length,
        insert_size_mean: insert_stats.mean(),
        insert_size_sd: insert_stats.sd(),
    };
    Ok(res)
}

fn calculate_coverage(
    reader: &mut IndexedReader<Reader<File>>,
    header: &Header,
    locus: &str,
    pos: usize,
) -> Result<Option<(usize, Summariser, usize)>> {
    let region: Region = locus
        .parse()
        .map_err(|err: ParseError| Error::new(ErrorKind::Other, err.to_string()))?;
    let query: bam::reader::Query<std::fs::File> = reader.query(&header, &region)?;
    let mut cov = 0;
    let mut insert_stats = Summariser::new();
    let mut read_length = 0;
    let mut n = 0;
    for res in query {
        let rec = res?;
        n += 1;
        if rec.flags().is_unmapped() {
            continue;
        }
        if rec.flags().is_duplicate() {
            continue;
        }
        if rec.flags().is_secondary() {
            continue;
        }
        if rec.sequence().len() > read_length {
            read_length = rec.sequence().len();
        }
        match (rec.alignment_start(), rec.alignment_end()) {
            (Some(start_pos), Some(end_pos)) => {
                let st = start_pos.get();
                let en = end_pos.get();
                if st <= pos && pos <= en {
                    cov += 1;
                }
            }
            _ => {}
        }
        match (
            rec.reference_sequence_id(),
            rec.mate_reference_sequence_id(),
        ) {
            (Some(this_id), Some(mate_id)) => {
                let l = rec.template_length().abs();
                if this_id == mate_id && l < 2000 {
                    insert_stats.add(l as f64);
                }
            }
            _ => {}
        }
    }
    if n < 10 {
        // no reads here
        return Ok(None);
    }
    Ok(Some((cov, insert_stats, read_length)))
}

#[derive(Serialize, Deserialize, Debug)]
pub struct BasicInfo {
    bam: String,
    cov: f64,
    read_length: usize,
    insert_size_mean: f64,
    insert_size_sd: f64,
}

pub fn gather_chromosome_info(bam: &str) -> Result<(Vec<String>, Vec<usize>)> {
    let mut reader = bam::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = reader.read_header()?;
    let mut chrom_names: Vec<String> = Vec::new();
    let mut chrom_lengths: Vec<usize> = Vec::new();
    for item in header.reference_sequences().iter() {
        let chrom_name = item.0.to_string();
        let chrom_length = item.1.length().get();
        if chrom_length < 4 * 1024 * 1024 {
            // ignore anything < 4Mb
            continue;
        }
        if chrom_name.contains("_") {
            // skip all the alts.
            continue;
        }
        chrom_names.push(chrom_name);
        chrom_lengths.push(chrom_length);
    }
    Ok((chrom_names, chrom_lengths))
}

pub fn chromosome_ranges(bam: &str) -> Result<Vec<String>> {
    let (chrom_names, chrom_lengths) = gather_chromosome_info(bam)?;
    let mut res = Vec::new();
    for i in 0..chrom_names.len() {
        res.push(format!("{}:{}-{}", chrom_names[i], 1, chrom_lengths[i]));
    }
    println!("chunks={:?}", res);
    Ok(res)
}

pub fn split_genome(bam: &str, chunk_size: usize) -> Result<Vec<String>> {
    let (chrom_names, chrom_lengths) = gather_chromosome_info(bam)?;
    let mut res = Vec::new();
    for i in 0..chrom_names.len() {
        let mut pos = 1;
        while pos < chrom_lengths[i] {
            let st = pos;
            let en = min(pos + chunk_size - 1, chrom_lengths[i]);
            res.push(format!("{}:{}-{}", chrom_names[i], st, en));
            pos = en + 1;
        }
    }
    println!("chunks={:?}", res);
    Ok(res)
}

impl BasicInfo {
    pub fn from_bam(bam: &str) -> Result<BasicInfo> {
        let (chrom_names, chrom_lengths) = gather_chromosome_info(bam)?;
        let total_length: usize = chrom_lengths.iter().sum();

        chromosome_ranges(bam)?;

        let n: usize = 3000;
        let mut rng = rand::thread_rng();
        let mut xs: Vec<(usize, usize, usize, usize, usize)> = Vec::new();
        while xs.len() < n {
            let u: f64 = rng.gen();
            let v = (total_length as f64) * u;
            let j0 = v as usize;
            let j1 = j0 + 500;
            let j = (j0 + j1) / 2;

            let mut t = 0;
            for k in 0..chrom_lengths.len() {
                let l = chrom_lengths[k];
                let next_t = t + l;
                if t <= j0 && j1 < next_t {
                    let x = (j, k, j0 - t, j1 - t, j - t);
                    xs.push(x);
                }
                t = next_t;
            }
        }
        xs.sort();

        if false {
            let ys: Vec<(String, usize)> = Vec::from_iter(xs.iter().map(|x| {
                let (_, k, st, en, m) = x;
                let locus = format!("{}:{}-{}", chrom_names[*k], st, en);
                (locus, *m)
            }));
            let res = do_all_calculate_coverage(bam, &ys)?;
            println!("{}", serde_json::to_string(&res).unwrap());
            return Ok(res);
        }

        let ys: Vec<(String, usize)> = Vec::from_iter(xs.iter().map(|x| {
            let (_, k, st, en, m) = x;
            let locus = format!("{}:{}-{}", chrom_names[*k], st, en);
            (locus, *m)
        }));
        let res = do_all_calculate_coverage_2(bam, &ys)?;
        println!("{}", serde_json::to_string(&res).unwrap());
        Ok(res)
    }
}
