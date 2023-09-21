use std::{
    collections::{HashMap, HashSet},
    io::{Error, ErrorKind, Result},
    sync::mpsc::{channel, Sender},
    usize,
};

use noodles::{bam, core::Position, sam::Header};

use indicatif::{ProgressBar, ProgressStyle};
use threadpool::ThreadPool;

use crate::{bam_helpers::chromosome_levels, locus::Locus, regions::Regions, tsv::Table};

pub fn seqname(hdr: &Header, chrom_id: usize) -> String {
    hdr.reference_sequences()
        .get_index(chrom_id)
        .map(|x| x.0.to_string())
        .unwrap_or_default()
}

pub fn seqlen(hdr: &Header, chrom_id: usize) -> usize {
    hdr.reference_sequences()
        .get_index(chrom_id)
        .map(|x| x.1.length().get())
        .unwrap_or_default()
}

enum DiscReads {
    Done(String),
    Disc(String, HashMap<String, HashMap<u32, (u32, u32, u32)>>),
}

#[derive(Debug, Clone, Copy)]
enum Side {
    Left,
    Right,
}

pub fn get_disc(
    bam: &str,
    peaks_filename: &str,
    _repeats: &Option<Regions>,
    _verbose: bool,
) -> Result<()> {
    let peaks = Table::read_from_tsv(peaks_filename)?;

    let prog = ProgressBar::new(peaks.len() as u64);
    let sty = ProgressStyle::with_template(
        "{prefix} [{elapsed_precise}] [{wide_bar}] {percent}% ({pos}/{len})",
    )
    .unwrap();
    prog.set_style(sty);

    let chroms = chromosome_levels(bam)?;

    let mut loci: Vec<(u32, u32, u32)> = Vec::new();
    for row_num in 0..peaks.len() {
        let row = peaks.get(row_num);
        let chrom_name = row.get(0).to_string();
        let chrom_id = chroms.get(&chrom_name).unwrap();
        let peak_pos: usize = row
            .get(1)
            .parse::<usize>()
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
        let peak_sd: f64 = row
            .get(5)
            .parse::<f64>()
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
        let width = 150 + (3.0 * peak_sd) as usize;
        loci.push((chrom_id as u32, peak_pos as u32, width as u32));
    }

    let pool = ThreadPool::new(8);
    let (tx, rx) = channel();
    for some_loci in loci.chunks(5000) {
        let the_loci: Vec<(u32, u32, u32)> = some_loci.to_vec();
        let bam_name = bam.to_string();
        let tx: Sender<DiscReads> = tx.clone();
        let repeats = None;
        pool.execute(move || {
            get_disc_inner(&bam_name, the_loci, &tx, &repeats, &None).expect("get_clips failed");
        });
    }
    let mut todo: usize = peaks.len();
    loop {
        let msg = rx.recv().expect("recv failed");
        match msg {
            DiscReads::Done(_chrom_name) => {
                todo -= 1;
                if todo == 0 {
                    break;
                }
            }
            DiscReads::Disc(locus, stuff) => {
                for item in stuff.iter() {
                    let chrom_name = item.0;
                    for jtem in item.1.iter() {
                        println!(
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            locus, chrom_name, jtem.0, jtem.1 .0, jtem.1 .1, jtem.1 .2
                        );
                    }
                }
            }
        }
    }
    Ok(())
}

fn get_disc_inner(
    bam: &str,
    loci: Vec<(u32, u32, u32)>,
    tx: &Sender<DiscReads>,
    repeats: &Option<Regions>,
    _opt_prog: &Option<ProgressBar>,
) -> Result<()> {
    let chroms = chromosome_levels(bam)?;

    let mut reader = bam::indexed_reader::Builder::default().build_from_path(bam)?;
    let hdr = reader.read_header()?;

    let mut good_chroms: HashSet<usize> = HashSet::new();
    let mut bad_chroms: HashSet<usize> = HashSet::new();

    for (chrom_id, pos, width) in loci.iter() {
        let locus_left = Locus::new(*chrom_id, *pos - width, *pos);
        let locus_right = Locus::new(*chrom_id, *pos + 1, *pos + width);
        let mut counts: HashMap<String, HashMap<u32, (u32, u32, u32)>> = HashMap::new();
        let mut disc_count: usize = 0;

        for (locus, side) in [(locus_left, Side::Left), (locus_right, Side::Right)] {
            let first = Position::try_from(locus.first as usize).unwrap();
            let last = Position::try_from(locus.last as usize).unwrap();
            let region: noodles::core::Region =
                noodles::core::Region::new(locus.chrom_name(&chroms), first..=last);
            let chrom_name = locus.chrom_name(&chroms);

            let query: bam::reader::Query<std::fs::File> = reader.query(&hdr, &region)?;
            for res in query {
                let rec = res?;
                if rec.flags().is_unmapped() {
                    continue;
                }
                if rec.flags().is_duplicate() {
                    continue;
                }
                if rec.flags().is_secondary() {
                    continue;
                }
                if let Some(q) = rec.mapping_quality() {
                    if q.get() <= 12 {
                        continue;
                    }
                }

                match (
                    rec.reference_sequence_id(),
                    rec.alignment_start(),
                    rec.mate_reference_sequence_id(),
                    rec.mate_alignment_start(),
                ) {
                    (Some(chrom_id), Some(pos), Some(mate_chrom_id), Some(mate_pos)) => {
                        let px: usize = pos.get();
                        let mpx: usize = mate_pos.get();

                        let distance = if px > mpx { px - mpx } else { mpx - px };
                        if chrom_id == mate_chrom_id && distance < 2500 {
                            continue;
                        }
                        let mate_in_rep = if let Some(rpts) = repeats {
                            rpts.contains(&chrom_name, mate_pos.get())
                        } else {
                            false
                        };

                        if bad_chroms.contains(&mate_chrom_id) {
                            continue;
                        }
                        let mate_chrom_name = seqname(&hdr, mate_chrom_id);
                        if !good_chroms.contains(&mate_chrom_id) {
                            if mate_chrom_name.contains("_") {
                                bad_chroms.insert(mate_chrom_id);
                                continue;
                            } else {
                                good_chroms.insert(mate_chrom_id);
                            }
                        }

                        let chrom_entry: &mut HashMap<u32, (u32, u32, u32)> =
                            counts.entry(mate_chrom_name).or_insert(HashMap::new());
                        let mut pos_entry: &mut (u32, u32, u32) =
                            chrom_entry.entry(mpx as u32).or_insert((0, 0, 0));
                        match side {
                            Side::Left => {
                                pos_entry.0 += 1;
                            }
                            Side::Right => {
                                pos_entry.1 += 1;
                            }
                        }
                        if mate_in_rep {
                            pos_entry.2 += 1;
                        }
                        disc_count += 1;
                    }
                    _ => {}
                }
            }
        }

        if disc_count >= 2 {
            tx.send(DiscReads::Disc("".to_string(), counts))
                .expect("send failed");
        }
        tx.send(DiscReads::Done("".to_string()))
            .expect("send failed");
    }

    Ok(())
}
