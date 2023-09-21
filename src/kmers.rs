use num::{one, zero, Unsigned};
use rand::{seq::SliceRandom, thread_rng};

use crate::summarise::Summariser;

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Kmer(u64);

impl Kmer {
    pub fn make(seq: &str) -> Option<Kmer> {
        if seq.len() > 32 {
            return None;
        }
        let mut x = Kmer(0);
        for c in seq.chars() {
            let b = Kmer::base(c)?;
            x.0 = (x.0 << 2) | b.0;
        }
        Some(x)
    }

    pub fn make_many(k: usize, seq: &str) -> Option<Vec<Kmer>> {
        let mut xs = Vec::new();
        let mut i = 0;
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        for c in seq.chars() {
            match Kmer::base(c) {
                None => {
                    i = 0;
                    x.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    i += 1;
                    if i == k {
                        xs.push(Kmer(x.0 & msk));
                        i -= 1;
                    }
                }
            }
        }
        Some(xs)
    }

    pub fn with_many<F>(k: usize, seq: &str, mut f: F) -> ()
    where
        F: FnMut(&Kmer),
    {
        let mut i = 0;
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        for c in seq.chars() {
            match Kmer::base(c) {
                None => {
                    i = 0;
                    x.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    i += 1;
                    if i == k {
                        x.0 &= msk;
                        f(&x);
                        i -= 1;
                    }
                }
            }
        }
    }

    pub fn with_many_both<S, F>(k: usize, seq: &S, mut f: F) -> ()
    where
        S: AsRef<[u8]>,
        F: FnMut(&Kmer, &Kmer),
    {
        let shift = 2 * (k - 1);
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        let mut y = Kmer(0);
        let mut i = 0;
        for c in seq.as_ref() {
            match Kmer::byte(c) {
                None => {
                    i = 0;
                    x.0 = 0;
                    y.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    y.0 = (y.0 >> 2) | ((3 - b.0) << shift);
                    i += 1;
                    if i == k {
                        x.0 &= msk;
                        f(&x, &y);
                        i -= 1;
                    }
                }
            }
        }
    }

    pub fn base(c: char) -> Option<Kmer> {
        match c {
            'A' | 'a' => Some(Kmer(0)),
            'C' | 'c' => Some(Kmer(1)),
            'G' | 'g' => Some(Kmer(2)),
            'T' | 't' | 'U' | 'u' => Some(Kmer(3)),
            _ => None,
        }
    }

    pub fn byte(c: &u8) -> Option<Kmer> {
        match c {
            b'A' | b'a' => Some(Kmer(0)),
            b'C' | b'c' => Some(Kmer(1)),
            b'G' | b'g' => Some(Kmer(2)),
            b'T' | b't' | b'U' | b'u' => Some(Kmer(3)),
            _ => None,
        }
    }

    pub fn rev(&self, k: usize) -> Kmer {
        let m2 = 0x3333333333333333;
        let m3 = 0x0F0F0F0F0F0F0F0F;
        let m4 = 0x00FF00FF00FF00FF;
        let m5 = 0x0000FFFF0000FFFF;
        let m6 = 0x00000000FFFFFFFF;

        let mut x = Kmer(self.0);
        x.0 = ((x.0 >> 2) & m2) | ((x.0 & m2) << 2);
        x.0 = ((x.0 >> 4) & m3) | ((x.0 & m3) << 4);
        x.0 = ((x.0 >> 8) & m4) | ((x.0 & m4) << 8);
        x.0 = ((x.0 >> 16) & m5) | ((x.0 & m5) << 16);
        x.0 = ((x.0 >> 32) & m6) | ((x.0 & m6) << 32);
        x.0 >>= 64 - 2 * k;
        x
    }

    pub fn rev_comp(&self, k: usize) -> Kmer {
        Kmer(!self.0).rev(k)
    }

    pub fn render(&self, k: usize) -> String {
        let mut s = String::new();
        let mut y = self.rev(k);
        for _i in 0..k {
            match y.0 & 3 {
                0 => {
                    s.push('A');
                }
                1 => {
                    s.push('C');
                }
                2 => {
                    s.push('G');
                }
                3 => {
                    s.push('T');
                }
                _ => {
                    unreachable!();
                }
            }
            y.0 >>= 2;
        }
        s
    }

    pub fn ham(x: &Kmer, y: &Kmer) -> usize {
        let m1 = 0x5555555555555555;
        let z = x.0 ^ y.0;
        let v = (z | (z >> 1)) & m1;
        v.count_ones() as usize
    }

    pub fn frequency_vector<S, T>(k: usize, seq: &S) -> Option<Vec<T>>
    where
        S: AsRef<[u8]>,
        T: Unsigned + Clone + std::ops::AddAssign,
    {
        let n: usize = 1 << (2 * k);
        let mut v: Vec<T> = Vec::new();
        v.resize(n, zero::<T>());
        Kmer::with_many_both(k, seq, |x, y| {
            v[x.0 as usize] += one::<T>();
            v[y.0 as usize] += one::<T>();
        });
        Some(v)
    }
}

pub struct Cluster {
    centroid: Vec<f64>,
    members: Vec<usize>,
}

impl Cluster {
    pub fn new(v: &Vec<u8>) -> Cluster {
        let centroid: Vec<f64> = Vec::from_iter(v.iter().map(|x| (*x as f64)));
        let members: Vec<usize> = Vec::new();
        Cluster { centroid, members }
    }

    pub fn add(&mut self, i: usize, v: &Vec<u8>) {
        self.members.push(i);
        let n = self.members.len() as f64;
        for j in 0..self.centroid.len() {
            self.centroid[j] = ((n - 1.0) * self.centroid[j] + (v[j] as f64)) / n;
        }
    }

    pub fn euclidean(&self, v: &Vec<u8>) -> f64 {
        let mut d2: f64 = 0.0;
        for j in 0..self.centroid.len() {
            d2 += (self.centroid[j] - (v[j] as f64)).powi(2);
        }
        d2.sqrt()
    }

    pub fn recalculate_mean(&mut self, vectors: &Vec<Vec<u8>>) {
        let n: f64 = self.members.len() as f64;
        for j in 0..self.centroid.len() {
            let mut sx = 0.0;
            for i in self.members.iter() {
                sx += vectors[*i][j] as f64;
            }
            self.centroid[j] = sx / n;
        }
    }
}

pub struct KMeans {}

impl KMeans {
    pub fn cluster_n(vectors: &Vec<Vec<u8>>, num_clusters: usize) -> Vec<Cluster> {
        let mut perm: Vec<usize> = Vec::new();
        for i in 0..vectors.len() {
            perm.push(i);
        }
        let mut rng: rand::rngs::ThreadRng = thread_rng();
        perm.shuffle(&mut rng);
        let mut clusters = Vec::new();
        for j in 0..num_clusters {
            let i = perm[j];
            clusters.push(Cluster::new(&vectors[i]));
        }
        let mut assignment: Vec<usize> = Vec::new();
        assignment.resize(vectors.len(), 0);
        for w in 0..50 {
            for j in 0..num_clusters {
                clusters[j].members.clear();
            }
            perm.shuffle(&mut rng);
            for i in perm.iter() {
                KMeans::assign(&mut clusters, *i, &vectors[*i]);
            }
            for j in 0..num_clusters {
                clusters[j].recalculate_mean(vectors);
            }
            let mut reassigned = 0;
            let mut s = Summariser::new();
            for j in 0..num_clusters {
                for i in clusters[j].members.iter() {
                    if assignment[*i] != j {
                        reassigned += 1;
                    }
                    assignment[*i] = j;
                    let d = clusters[j].euclidean(&vectors[*i]);
                    s.add(d);
                }
            }
            println!("iteration {}, mean distance = {}, sd = {}", w, s.mean(), s.sd());
            println!("iteration {}, has {} reassignments", w, reassigned);
            if reassigned == 0 {
                break;
            }
        }
        clusters
    }

    fn assign(clusters: &mut [Cluster], i: usize, v: &Vec<u8>) {
        let num_clusters = clusters.len();
        let mut j = 0;
        let mut d = clusters[0].euclidean(v);
        for k in 1..num_clusters {
            let d0 = clusters[k].euclidean(v);
            if d0 < d || (d0 == d && clusters[k].members.len() < clusters[j].members.len()) {
                j = k;
                d = d0;
            }
        }
        clusters[j].members.push(i);
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_1() {
        {
            let x = Kmer::make("TCTG");
            assert_eq!(x, Some(Kmer(0xDE)));
        }
        {
            let x = Kmer::make("GATCGGT");
            assert_eq!(x, Some(Kmer(0x236B)));
        }
        {
            let x = Kmer::make("WIBBLE");
            assert_eq!(x, None);
        }
    }

    #[test]
    fn test_make_many_1() {
        let k: usize = 11;
        let seq = "ACTGTAGTCCTAGCTACTCGGGAGGCTGAGGCACAAGAATTGCTTGAACCCGGGAAGCAGAGGTTGGAGTGAACCA";
        let opt_xs = Kmer::make_many(k, seq);
        assert_ne!(opt_xs, None);
        let xs = opt_xs.unwrap();
        assert_eq!(seq.len(), 76);
        assert_eq!(xs.len(), 66);
        assert_eq!(xs[0].render(k), "ACTGTAGTCCT");
        assert_eq!(xs[7].render(k), "TCCTAGCTACT");
    }

    #[test]
    fn test_reverse_complement_1() {
        let x = Kmer::make("ACTGTAGTCCT").unwrap();
        assert_eq!(x.rev_comp(11).render(11), "AGGACTACAGT");
    }
}
