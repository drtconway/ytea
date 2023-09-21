use std::{cmp::max, collections::HashMap, path::Path, ffi::OsStr, fs::File, io::{BufReader, Result}};

use flate2::read::MultiGzDecoder;
use noodles::bed;

use crate::factor::Factor;

pub struct Regions {
    chroms: Factor,
    index: HashMap<usize, Vec<(usize, usize)>>,
}

impl Regions {
    pub fn new() -> Regions {
        Regions {
            chroms: Factor::new(),
            index: HashMap::new(),
        }
    }

    pub fn add(&mut self, chrom_name: &str, first: usize, last: usize) {
        let chrom_id = self.chroms.get_or_add(chrom_name);
        self.index
            .entry(chrom_id)
            .or_insert(Vec::new())
            .push((first, last));
    }

    pub fn simplify(&mut self) {
        for item in self.index.iter_mut() {
            item.1.sort();
            let mut i = 1;
            let mut j = 0;
            while i < item.1.len() {
                if item.1[j].1 >= item.1[i].0 {
                    item.1[j].1 = max(item.1[j].1, item.1[i].1);
                    i += 1;
                } else {
                    j += 1;
                    if j < i {
                        item.1[j] = item.1[i];
                    }
                    i += 1;
                }
            }
            while j + 1 < item.1.len() {
                item.1.pop();
            }
        }
    }

    pub fn contains(&self, chrom_name: &str, pos: usize) -> bool {
        match self.chroms.get(chrom_name).and_then(|chrom_id| self.index.get(&chrom_id)) {
            None => { false }
            Some(v) => {
                let mut lo = 0;
                let mut hi = v.len();
                while lo < hi {
                    let mid = (lo + hi) / 2;
                    let (first, last) = v[mid];
                    if pos < first {
                        hi = mid;
                    } else if pos > last {
                        lo = mid + 1;
                    } else {
                        return true;
                    }
                }
                lo < v.len() && v[lo].0 <= pos && pos <= v[lo].1
            }
        }
    }

    pub fn subset(&self, chrom_name: &str) -> Regions {
        let chroms = self.chroms.clone();
        let mut index = HashMap::new();
        if let Some(chrom_id) = chroms.get(chrom_name) {
            if let Some(v) = self.index.get(&chrom_id) {
                index.insert(chrom_id, v.to_vec());
            }
        }
        Regions { chroms, index }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_regions_1() {
        let mut r = Regions::new();
        r.add(&"chr1", 100, 110);
        assert_eq!(r.contains(&"chr1", 99), false);
        assert_eq!(r.contains(&"chr1", 100), true);
        assert_eq!(r.contains(&"chr1", 101), true);
        assert_eq!(r.contains(&"chr1", 110), true);
        assert_eq!(r.contains(&"chr1", 111), false);
    }

    #[test]
    fn test_regions_2() {
        let mut r = Regions::new();
        r.add(&"chr1", 100, 110);
        r.add(&"chr1", 120, 130);
        r.add(&"chr1", 125, 128);
        r.add(&"chr1", 140, 150);
        r.simplify();
        assert_eq!(r.index.get(&0), Some(&vec![(100, 110), (120, 130), (140, 150)]));
    }
}

fn with_bed_file<F>(filename: &str, mut f: F) -> Result<()>
where
    F: FnMut(&bed::Record<4>) -> Result<()>,
{
    let path = Path::new(filename);
    if path.extension() == Some(OsStr::new("gz")) {
        let mut reader = File::open(filename)
            .map(MultiGzDecoder::new)
            .map(BufReader::new)
            .map(bed::Reader::new)?;
        for rec in reader.records() {
            let r = rec?;
            f(&r)?;
        }
    } else {
        let mut reader = File::open(filename)
            .map(BufReader::new)
            .map(bed::Reader::new)?;
        for rec in reader.records() {
            let r = rec?;
            f(&r)?;
        }
    }
    Ok(())
}

pub fn regions_from_bed(filename: &str) -> Result<Regions> {
    let mut r = Regions::new();
    with_bed_file(filename, |rec| {
        let chrom_name = rec.reference_sequence_name();
        r.add(chrom_name, rec.start_position().get(), rec.end_position().get());
        Ok(())
    })?;
    r.simplify();
    Ok(r)
}