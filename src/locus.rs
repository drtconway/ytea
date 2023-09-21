use crate::factor::Factor;


#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct Locus {
    pub chrom_id: u32,
    pub first: u32,
    pub last: u32,
}

impl Locus {
    pub fn new(chrom_id: u32, first: u32, last: u32) -> Locus {
        Locus { chrom_id, first, last }
    }

    pub fn to_string(&self, chrom_names: &Factor) -> String {
        format!("{}:{}-{}", chrom_names.levels[self.chrom_id as usize], self.first, self.last)
    }

    pub fn chrom_name(&self, chrom_names: &Factor) -> String {
        chrom_names.levels[self.chrom_id as usize].to_string()
    }
}

