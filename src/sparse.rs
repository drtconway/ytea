use crate::ransel::Ransel;

pub struct Sparse {
    xs: Vec<u32>,
}

impl Sparse {
    pub fn new(xs: &[u32]) -> Sparse {
        Sparse { xs: Vec::from(xs) }
    }
}

impl Ransel for Sparse {
    type Position = u32;
    type Rank = usize;

    fn count(&self) -> usize {
        self.xs.len()
    }

    fn size(&self) -> u32 {
        self.xs.last().unwrap_or(&0) + 1
    }

    fn rank(&self, x: u32) -> usize {
        let mut lo: usize = 0;
        let mut hi: usize = self.xs.len();
        while lo < hi {
            let mid: usize = (lo + hi) / 2;
            if self.xs[mid] < x {
                lo = mid + 1;
            } else if self.xs[mid] > x {
                hi = mid;
            } else {
                return mid;
            }
        }
        return lo;
    }

    fn select(&self, i: usize) -> u32 {
        return self.xs[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sparse_1() {
        let xs: [u32; 5] = [2, 3, 7, 9, 11];
        let s = Sparse::new(&xs);
        assert_eq!(s.count(), 5);
        assert_eq!(s.size(), 12);
        assert_eq!(s.rank(0), 0);
        assert_eq!(s.rank(1), 0);
        assert_eq!(s.rank(2), 0);
        assert_eq!(s.rank(3), 1);
        assert_eq!(s.rank(4), 2);
        assert_eq!(s.rank(5), 2);
        assert_eq!(s.rank(6), 2);
        assert_eq!(s.rank(7), 2);
        assert_eq!(s.rank(8), 3);
        assert_eq!(s.rank(9), 3);
        assert_eq!(s.rank(10), 4);
        assert_eq!(s.rank(11), 4);
        assert_eq!(s.rank(12), 5);
        assert_eq!(s.rank(13), 5);
        assert_eq!(s.select(0), 2);
        assert_eq!(s.select(1), 3);
        assert_eq!(s.select(2), 7);
        assert_eq!(s.select(3), 9);
        assert_eq!(s.select(4), 11);
    }
}
