use crate::ransel::Ransel;
use bitintr::Popcnt;

pub struct Dense {
    x_max: u32,
    x_count: usize,
    words: Vec<u64>,
    toc: Vec<usize>,
}

impl Dense {
    pub fn new(xs: &[u32]) -> Dense {
        let mut x_max: usize = 0;
        for x in xs.iter() {
            let xu: usize = *x as usize;
            if xu > x_max {
                x_max = xu;
            }
        }
        let w_max = 1 + (x_max + 63) / 64;
        let mut words: Vec<u64> = Vec::new();
        words.resize(w_max, 0);
        for x in xs.iter() {
            let w: usize = (x >> 6) as usize;
            let b: u64 = (x & 63) as u64;
            words[w] |= 1u64 << b;
        }
        let mut toc: Vec<usize> = Vec::new();
        let mut cum: usize = 0;
        for word in words.iter() {
            toc.push(cum);
            cum += (*word).popcnt() as usize;
        }
        toc.push(cum);
        Dense {
            x_max: x_max as u32,
            x_count: xs.len(),
            words,
            toc,
        }
    }
}

impl Ransel for Dense {
    type Position = u32;
    type Rank = usize;

    fn size(&self) -> Self::Position {
        self.x_max + 1
    }

    fn count(&self) -> Self::Rank {
        self.x_count
    }

    fn rank(&self, x: Self::Position) -> Self::Rank {
        let w: usize = (x >> 6) as usize;
        let b: u64 = (x & 63) as u64;
        let m: u64 = (1 << b) - 1;
        self.toc[w] + (self.words[w] & m).popcnt() as usize
    }

    fn select(&self, i: Self::Rank) -> Self::Position {
        let mut k = 0;
        let mut c = self.toc.len();
        while c > 0 {
            let s = c / 2;
            let l = k + s;
            if self.toc[l] < i {
                k = l + 1;
                c -= s + 1;
            } else {
                c = s;
            }
        }
        if k > 0 && self.toc[k] > i {
            k -= 1;
        }
        let mut x0: Self::Position = (k << 6) as Self::Position;
        let mut j = i + 1 - self.toc[k];
        let mut w = self.words[k];
        assert!(w.popcnt() >= j as u64);
        loop {
            let v = w.trailing_zeros();
            x0 += v + 1;
            w >>= v + 1;
            j -= 1;
            if j == 0 {
                break;
            }
        }
        x0 - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dense_1() {
        let xs: [u32; 5] = [2, 3, 7, 9, 11];
        let s = Dense::new(&xs);
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

    #[test]
    fn test_dense_2() {
        let xs: [u32; 15] = [0, 7, 15, 31, 49, 63, 64, 81, 127, 128, 173, 254, 255, 256, 257];
        let s = Dense::new(&xs);
        assert_eq!(s.count(), 15);
        assert_eq!(s.size(), 258);
        for i in 0..xs.len() {
            assert_eq!(s.select(i), xs[i]);
            assert_eq!(s.rank(xs[i]), i);
        }
    }

}
