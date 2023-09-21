use std::{
    collections::{BTreeMap},
    io::{Error, ErrorKind, Result, Write},
    sync::mpsc::{channel, Sender},
    usize,
};

use noodles::{
    bam,
    core::region::ParseError,
    fastq::{self, record::Definition},
    sam::{record::cigar::op::Kind, Header},
};

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use threadpool::ThreadPool;

use crate::{info::chromosome_ranges, files::open_writer, regions::Regions, bam_helpers::{with_cigar}};


pub enum Side {
    Left,
    Right,
}

pub struct ClipRec {
    left_count: usize,
    right_count: usize,
    mate_in_rep_count: usize,
}

impl ClipRec {
    pub fn new() -> ClipRec {
        ClipRec {
            left_count: 0,
            right_count: 0,
            mate_in_rep_count: 0,
        }
    }

    pub fn add_left(&mut self, mate_in_rep: bool) {
        self.left_count += 1;
        if mate_in_rep {
            self.mate_in_rep_count += 1;
        }
    }

    pub fn add_right(&mut self, mate_in_rep: bool) {
        self.right_count += 1;
        if mate_in_rep {
            self.mate_in_rep_count += 1;
        }
    }

    pub fn write(&self, out: &mut dyn Write, chrom_name: &str, clip_pos: usize) -> Result<()> {
        if self.left_count + self.right_count >= 2 {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                chrom_name, clip_pos, self.left_count, self.right_count, self.mate_in_rep_count
            )
        } else {
            Ok(())
        }
    }
}

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

#[derive(Debug)]
pub enum ClipData {
    Done(String),
    Fastq(Vec<(String, String, String)>),
    Clips(String, Vec<(usize, usize, usize, usize)>),
}

pub fn get_clips_2(
    bam: &str,
    fastq_filename: &Option<&str>,
    counts_filename: &Option<&str>,
    repeats: &Option<Regions>,
    verbose: bool,
) -> Result<()> {
    let mut opt_out_holder: Option<Box<dyn Write>> = match fastq_filename {
        None => None,
        Some(filename) => {
            let f = open_writer(filename)?;
            Some(f)
        }
    };
    let mut opt_out = match &mut opt_out_holder {
        None => None,
        Some(boxed) => {
            let writer = boxed.as_mut();
            Some(fastq::Writer::new(writer))
        }
    };

    let mut opt_csv: Option<Box<dyn Write>> = match counts_filename {
        None => None,
        Some(filename) => {
            let mut f = open_writer(filename)?;
            writeln!(
                f.as_mut(),
                "chrom\tpos\tleft_count\tright_count\tmate_in_rep"
            )?;
            Some(f)
        }
    };

    let multi = MultiProgress::new();
    let sty = ProgressStyle::with_template(
        "{prefix} [{elapsed_precise}] [{wide_bar}] {percent}% ({pos}/{len})",
    )
    .unwrap();

    let pool = ThreadPool::new(8);
    let (tx, rx) = channel();
    let chrom_ranges = chromosome_ranges(bam)?;
    for chrom_range in chrom_ranges.iter() {
        let opt_prog = if verbose {
            let prog = multi.add(ProgressBar::new(1));
            prog.set_style(sty.clone());
            Some(prog)
        } else {
            None
        };
        let locus = chrom_range.to_string();
        let locus_parts: Vec<&str> = locus.split(':').collect();
        let chrom_name = locus_parts[0].to_string();
        let bam_name = bam.to_string();
        let tx = tx.clone();
        let repeats = repeats.as_ref().map(|rpts| rpts.subset(&chrom_name));
        pool.execute(move || {
            get_clips_2_inner(&bam_name, &locus, &tx, &repeats, &opt_prog).expect("get_clips failed");
        });
    }
    let mut todo = chrom_ranges.len();
    loop {
        let msg = rx.recv().expect("recv failed");
        match msg {
            ClipData::Done(_chrom_name) => {
                todo -= 1;
                if todo == 0 {
                    break;
                }
            }
            ClipData::Fastq(items) => match &mut opt_out {
                None => {}
                Some(out) => {
                    for (rd_id, clip_seq, clip_qual) in items {
                        let fq =
                            fastq::Record::new(Definition::new(rd_id, ""), clip_seq, clip_qual);
                        out.write_record(&fq)?;
                    }
                }
            },
            ClipData::Clips(chrom_name, clips) => match &mut opt_csv {
                None => {}
                Some(csv) => {
                    for clip in clips {
                        writeln!(
                            csv,
                            "{}\t{}\t{}\t{}\t{}",
                            chrom_name, clip.0, clip.1, clip.2, clip.3
                        )?;
                    }
                }
            },
        }
    }
    Ok(())
}

pub fn get_clips_2_inner(
    bam: &str,
    locus: &str,
    tx: &Sender<ClipData>,
    repeats: &Option<Regions>,
    opt_prog: &Option<ProgressBar>,
) -> Result<()> {
    let mut fastq_data: Vec<(String, String, String)> = Vec::new();
    let mut clip_data: Vec<(usize, usize, usize, usize)> = Vec::new();

    let mut counts: BTreeMap<usize, ClipRec> = BTreeMap::new();

    let mut reader = bam::indexed_reader::Builder::default().build_from_path(bam)?;
    let hdr = reader.read_header()?;
    let region: noodles::core::Region = locus
        .parse()
        .map_err(|err: ParseError| Error::new(ErrorKind::Other, err.to_string()))?;
    let chrom_name = region.name().to_string();
    let interval = region.interval();
    let chrom_length =
        interval.end().expect("no end").get() - interval.start().expect("no start").get();
    if let Some(prog) = opt_prog {
        prog.set_prefix(chrom_name.to_string());
        prog.set_length(chrom_length as u64);
        prog.set_position(0);
    }
    let query: bam::reader::Query<std::fs::File> = reader.query(&hdr, &region)?;
    let frac = 2048.0 / (chrom_length as f64);
    let mut prev_pos = 0;
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

        match (rec.reference_sequence_id(), rec.alignment_start()) {
            (Some(_chrom_id), Some(pos)) => {

                if let Some(prog) = opt_prog {
                    let px: usize = pos.get();
                    if (prev_pos as f64 * frac) as i64 != (px as f64 * frac) as i64 {
                        prog.set_position(px as u64);
                    }
                    prev_pos = px;
                }

                let px: usize = pos.get();
                while let Some(item) = counts.first_key_value() {
                    let clip_pos: usize = *item.0;
                    if clip_pos >= px {
                        break;
                    }
                    let clip_rec = item.1;
                    if clip_rec.left_count + clip_rec.right_count >= 2 {
                        clip_data.push((
                            clip_pos,
                            clip_rec.left_count,
                            clip_rec.right_count,
                            clip_rec.mate_in_rep_count,
                        ));
                    }
                    counts.pop_first();
                }

                let mut mate_chrom_name: String = "*".to_string();
                let mut mate_pos: usize = 0;
                let mut mate_in_rep = false;
                match (rec.mate_reference_sequence_id(), rec.mate_alignment_start()) {
                    (Some(mate_chrom_id_0), Some(mate_pos_0)) => {
                        mate_chrom_name = seqname(&hdr, mate_chrom_id_0);
                        mate_pos = mate_pos_0.get();
                        if let Some(rpts) = repeats {
                            mate_in_rep = rpts.contains(&chrom_name, mate_pos);
                        }
                    }
                    _ => {}
                }

                let seg = if rec.flags().is_first_segment() {
                    "1"
                } else {
                    "2"
                };

                with_cigar(&rec, |_i, p, q, k, l| {
                    match k {
                        Kind::SoftClip => {
                            if q == 0 {
                                let clip_rec = counts.entry(p + l).or_insert(ClipRec::new());
                                clip_rec.add_left(mate_in_rep);
                            } else {
                                let clip_rec = counts.entry(p).or_insert(ClipRec::new());
                                clip_rec.add_right(mate_in_rep);
                            }
                            let seq: String = rec.sequence().to_string();
                            if seq.len() < 9 {
                                return Ok(());
                            }
                            let qual: String = rec.quality_scores().to_string();
                            let clip_seq = &seq[q..(q + l)];
                            let clip_qual = &qual[q..(q + l)];
                            let side = if q == 0 { "L" } else { "R" };
                            let rd_id_base =
                                rec.read_name().map(|r| r.to_string()).unwrap_or_default();
                            let rd_id = format!(
                                "{}~{}~{}~{}~{}~{}~{}",
                                rd_id_base, chrom_name, p, mate_chrom_name, mate_pos, side, seg
                            );
                            fastq_data.push((rd_id, clip_seq.to_string(), clip_qual.to_string()));
                        }
                        _ => {}
                    }
                    Ok(())
                })?;
            }
            _ => {}
        }

        if fastq_data.len() > 1000 {
            tx.send(ClipData::Fastq(fastq_data)).expect("send failed");
            fastq_data = Vec::new();
        }
        if clip_data.len() > 1000 {
            tx.send(ClipData::Clips(chrom_name.to_string(), clip_data))
                .expect("send failed");
            clip_data = Vec::new();
        }
    }

    while let Some(item) = counts.first_key_value() {
        let clip_pos: usize = *item.0;
        let clip_rec = item.1;
        if clip_rec.left_count + clip_rec.right_count >= 2 {
            clip_data.push((
                clip_pos,
                clip_rec.left_count,
                clip_rec.right_count,
                clip_rec.mate_in_rep_count,
            ));
        }
        counts.pop_first();
    }

if fastq_data.len() > 0 {
        tx.send(ClipData::Fastq(fastq_data)).expect("send failed");
    }
    if clip_data.len() > 0 {
        tx.send(ClipData::Clips(chrom_name.to_string(), clip_data))
            .expect("send failed");
    }
    if let Some(prog) = opt_prog {
        prog.finish_and_clear();
    }
    tx.send(ClipData::Done(chrom_name)).expect("send failed");

    Ok(())
}
