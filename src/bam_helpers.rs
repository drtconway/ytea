use std::{
    
    io::{Error, ErrorKind, Result},
};

use noodles::{
    bam,
    core::region::ParseError,
    sam::{self, record::cigar::op::Kind},
};

use crate::factor::Factor;

pub fn with_alignments<F>(filename: &str, mut f: F) -> Result<()>
where
    F: FnMut(&sam::Header, &sam::alignment::Record) -> Result<()>,
{
    let mut reader = bam::reader::Builder::default().build_from_path(filename)?;
    let header = reader.read_header()?;
    for res in reader.records(&header) {
        let rec = res?;
        f(&header, &rec)?;
    }
    Ok(())
}

pub fn with_mapped_alignments<F>(filename: &str, mut f: F) -> Result<()>
where
    F: FnMut(&sam::Header, &sam::alignment::Record) -> Result<()>,
{
    let mut reader = bam::reader::Builder::default().build_from_path(filename)?;
    let header = reader.read_header()?;
    for res in reader.records(&header) {
        let rec = res?;
        if rec.flags().is_unmapped() {
            continue;
        }
        if rec.flags().is_duplicate() {
            continue;
        }
        f(&header, &rec)?;
    }
    Ok(())
}

pub fn with_mapped_alignment_in_regions<F>(
    filename: &str,
    region_name: &str,
    mut f: F,
) -> Result<()>
where
    F: FnMut(&sam::Header, &sam::alignment::Record) -> Result<()>,
{
    let mut reader = bam::indexed_reader::Builder::default().build_from_path(filename)?;
    let header = reader.read_header()?;
    let region = region_name
        .parse()
        .map_err(|err: ParseError| Error::new(ErrorKind::Other, err.to_string()))?;
    let query: bam::reader::Query<std::fs::File> = reader.query(&header, &region)?;
    for res in query {
        let rec = res?;
        if rec.flags().is_unmapped() {
            continue;
        }
        if rec.flags().is_duplicate() {
            continue;
        }
        if let Some(q) = rec.mapping_quality() {
            if q.get() < 12 {
                continue;
            }
        }
        f(&header, &rec)?;
    }
    Ok(())
}

pub fn with_cigar<F>(rd: &sam::alignment::Record, mut f: F) -> Result<()>
where
    F: FnMut(usize, usize, usize, Kind, usize) -> Result<()>,
{
    let cig = rd.cigar();
    match rd.alignment_start() {
        None => {}
        Some(p0) => {
            let mut i: usize = 0;
            let mut p: usize = p0.get();
            let mut q: usize = 0;
            for op in cig.iter() {
                let n: usize = op.len();
                f(i, p, q, op.kind(), n)?;
                i += 1;
                match op.kind() {
                    Kind::Match => {
                        p += n;
                        q += n;
                    }
                    Kind::Insertion => {
                        q += n;
                    }
                    Kind::Deletion => {
                        p += n;
                    }
                    Kind::HardClip => {
                        // nothing!
                    }
                    Kind::SoftClip => {
                        q += n;
                    }
                    Kind::Skip => {
                        p += n;
                    }
                    Kind::Pad => {
                        // nothing!
                    }
                    Kind::SequenceMatch => {
                        p += n;
                        q += n;
                    }
                    Kind::SequenceMismatch => {
                        p += n;
                        q += n;
                    }
                }
            }
        }
    }

    return Ok(());
}

pub fn chromosome_levels(bam: &str) -> Result<Factor> {
    let mut reader = bam::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = reader.read_header()?;
    let mut items = Vec::new();
    let n = header.reference_sequences().len();
    for chrom_id in 0..n {
        let opt_chrom_name = header.reference_sequences().get_index(chrom_id);
        if let Some(chrom_name_ref) = opt_chrom_name {
            let chrom_name = chrom_name_ref.0.to_string();
            items.push((chrom_id, chrom_name));
        }
    }
    Ok(Factor::from_pairs(&items))
}