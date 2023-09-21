use crate::{
    files::{open_reader, open_writer},
    kmers::Kmer,
};
use noodles::{
    core::Position,
    fasta::{self, record::Definition},
};
use regex::Regex;
use sha2::{Digest, Sha256};
use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader, Error, ErrorKind, Result},
};

pub enum Strand {
    Fwd,
    Rev,
}

pub struct RepeatDb {
    pub chroms: Vec<String>,
    pub begins: Vec<u32>,
    pub ends: Vec<u32>,
    pub strands: Vec<String>,
    pub labels: Vec<String>,
    pub classes: Vec<String>,
}

pub fn prepare_repeats(
    rpts: &str,
    genome: &str,
    dest_name: &str,
    opt_window: &Option<usize>,
    opt_class: Option<String>,
    opt_filter: Option<String>,
) -> Result<()> {
    let filter: Option<Regex> = opt_filter.map(|re| Regex::new(&re).unwrap());

    let db: RepeatDb = load_repeat_db(rpts)?;

    if opt_class == Some("_LIST_".to_string()) {
        let mut xs: HashSet<String> = HashSet::new();
        for class in db.classes.iter() {
            xs.insert(class.to_string());
        }
        let mut ys: Vec<String> = Vec::from_iter(xs.iter().map(|s| s.to_string()));
        ys.sort();
        for x in ys.iter() {
            println!("{}", x);
        }
        return Ok(());
    }

    let mut db_idx: HashMap<String, (usize, usize)> = HashMap::new();
    for i in 0..db.chroms.len() {
        let chrom_name: &str = &db.chroms[i];
        let mut entry = db_idx.entry(chrom_name.to_string()).or_insert((i, i));
        entry.1 = i;
    }
    let mut boxed = open_reader(genome)?;
    let buf = BufReader::new(boxed.as_mut());
    let mut genome_reader = fasta::Reader::new(buf);
    let mut fa_file = open_writer(dest_name)?;
    let mut fa_out = fasta::Writer::new(fa_file.as_mut());
    let window = match opt_window {
        None => 0,
        Some(w) => *w,
    };

    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut num_processed: usize = 0;

    let mut elements: HashMap<String, Vec<Vec<u8>>> = HashMap::new();

    for res in genome_reader.records() {
        let rec = res?;
        let chrom_name = rec.name();
        let chrom_sequence = rec.sequence();
        if let Some((j, k)) = db_idx.get(chrom_name) {
            for i in *j..=*k {
                let chrom_name: &str = &db.chroms[i];
                let begin = db.begins[i];
                let end = db.ends[i];
                let label: &str = &db.labels[i];
                let class: &str = &db.classes[i];
                if let Some(cls) = &opt_class {
                    if cls != class {
                        continue;
                    }
                }
                if let Some(re) = &filter {
                    if !re.is_match(label) {
                        continue;
                    }
                }
                num_processed += 1;
                if true {
                    let st = Position::try_from(begin as usize).expect("bad start");
                    let en = Position::try_from(end as usize).expect("bad end");
                    let seq: fasta::record::Sequence =
                        chrom_sequence.slice(st..=en).expect("bad interval");
                    let v: Vec<u8> = Kmer::frequency_vector(4, &seq).unwrap();
                    let item = elements.entry(label.to_string()).or_insert(Vec::new());
                    item.push(v);
                    let digest = Sha256::digest(seq);
                    let mut x: u64 = 0;
                    for l in 0..8 {
                        x = (x << 8) | (digest[l] as u64);
                    }
                    if !seen.insert(x) {
                        continue;
                    }
                }
                let seq_id = format!(
                    "{}:{}-{}~{}~{}~{}",
                    chrom_name, begin, end, window, class, label
                );

                let st0 = if (begin as usize) <= window {
                    1
                } else {
                    begin as usize - window
                };
                let en0 = if (end as usize) + window > chrom_sequence.len() {
                    chrom_sequence.len()
                } else {
                    (end as usize) + window
                };
                let st = Position::try_from(st0).expect("bad start");
                let en = Position::try_from(en0).expect("bad end");
                let seq = chrom_sequence.slice(st..=en).expect("bad interval");
                let fa = fasta::Record::new(Definition::new(seq_id, None), seq);
                fa_out.write_record(&fa)?;
            }
        }
    }

    println!(
        "{:?}",
        elements
            .iter()
            .map(|itm| (itm.0.to_string(), itm.1.len()))
            .collect::<Vec<(String, usize)>>()
    );
    println!(
        "{} unique sequences found out of {} scanned.",
        seen.len(),
        num_processed
    );
    Ok(())
}

pub fn load_repeat_db(filename: &str) -> Result<RepeatDb> {
    let mut boxed = open_reader(filename)?;
    let reader = BufReader::new(boxed.as_mut());
    let mut line_number = 0;
    let mut chroms: Vec<String> = Vec::new();
    let mut begins: Vec<u32> = Vec::new();
    let mut ends: Vec<u32> = Vec::new();
    let mut strands: Vec<String> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    let mut classes: Vec<String> = Vec::new();
    for res in reader.lines() {
        let line = res?;
        line_number += 1;
        if line_number <= 3 {
            continue;
        }
        let parts = Vec::from_iter(line.split_whitespace());
        if parts.len() != 15 {
            println!(
                "{}:{}: expected 15 fields, got {}",
                filename,
                line_number,
                parts.len()
            );
            continue;
        }
        let chrom = parts[4].to_string();
        chroms.push(chrom);
        let begin: u32 = parts[5]
            .parse::<u32>()
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
        begins.push(begin);
        let end: u32 = parts[6]
            .parse::<u32>()
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
        ends.push(end);
        let strand = match parts[8] {
            "+" => Ok("+".to_string()),
            "C" => Ok("-".to_string()),
            _ => Err(Error::new(
                ErrorKind::Other,
                format!("expected '+' or 'C'; got '{}'", parts[8]),
            )),
        }?;
        strands.push(strand);
        let label = parts[9].to_string();
        labels.push(label);
        let class = parts[10].to_string();
        classes.push(class);
    }
    Ok(RepeatDb {
        chroms,
        begins,
        ends,
        strands,
        labels,
        classes,
    })
}
