use std::{collections::HashMap};

use crate::stabby::{Interval, Stabby};

pub struct FeatureIndexBuilder {
    accumulator: HashMap<String, HashMap<Interval, Vec<String>>>,
}

impl FeatureIndexBuilder {
    pub fn new() -> FeatureIndexBuilder {
        FeatureIndexBuilder {
            accumulator: HashMap::new(),
        }
    }

    pub fn add(&mut self, chrom: &str, first: u32, last: u32, value: &str) {
        let ivl = Interval::new(first, last);
        if !self.accumulator.contains_key(chrom) {
            self.accumulator.insert(chrom.to_string(), HashMap::new());
        }
        let chrom_accumulator: &mut HashMap<Interval, Vec<String>> = self
            .accumulator
            .get_mut(chrom)
            .expect("vanishing chromosome");
        if !chrom_accumulator.contains_key(&ivl) {
            chrom_accumulator.insert(ivl, Vec::new());
        }
        let ivl_features: &mut Vec<String> =
            chrom_accumulator.get_mut(&ivl).expect("vanishing interval");
        let old_len = ivl_features.len();
        ivl_features.push(value.to_string());
        let new_len = ivl_features.len();
        if (old_len as f64).sqrt().floor() != (new_len as f64).sqrt().floor() {
            ivl_features.sort();
            ivl_features.dedup();
        }
    }
}

pub struct FeatureIndex {
    features: HashMap<String, HashMap<Interval, Vec<String>>>,
    index: HashMap<String,Stabby>
}

impl FeatureIndex {
    pub fn new(builder: &FeatureIndexBuilder) -> FeatureIndex {
        let mut features: HashMap<String, HashMap<Interval, Vec<String>>> = HashMap::new();
        let mut index: HashMap<String,Stabby> = HashMap::new();
        for chrom_entry in builder.accumulator.iter() {
            let chrom = chrom_entry.0;
            let mut chrom_features: HashMap<Interval, Vec<String>> = HashMap::new();
            let mut ivls: Vec<Interval> = Vec::new();
            for ivl_entry in chrom_entry.1.iter() {
                let ivl = ivl_entry.0;
                let mut items: Vec<String> = ivl_entry.1.to_vec();
                items.sort();
                items.dedup();
                chrom_features.insert(*ivl, items);
                ivls.push(*ivl);
            }

            features.insert(chrom.to_string(), chrom_features);
            index.insert(chrom.to_string(), Stabby::new(&ivls));
        }
        FeatureIndex { features, index }
    }

    pub fn find(&self, chrom: &str, pos: u32) -> Vec<&str> {
        let mut res: Vec<&str> = Vec::new();
        match (self.features.get(chrom), self.index.get(chrom)) {
            (None, None) => {}
            (Some(feats), Some(idx)) => {
                let hits: Vec<Interval> = idx.stab(pos);
                for ivl in hits.iter() {
                    match feats.get(ivl) {
                        None => {}
                        Some(values) => {
                            for v in values.iter() {
                                res.push(v);
                            }
                        }
                    }
                }
            }
            _ => unreachable!()
        }
        res
    }
}