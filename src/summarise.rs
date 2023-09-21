#[derive(Clone, Copy, Debug)]
pub struct Summariser {
    pub n: usize,
    sx: f64,
    sx2: f64,
}

impl Summariser {
    pub fn new() -> Summariser {
        Summariser {
            n: 0,
            sx: 0.0,
            sx2: 0.0,
        }
    }

    pub fn add(&mut self, x: f64) {
        self.n += 1;
        self.sx += x;
        self.sx2 += x * x;
    }

    pub fn add_multiple(&mut self, x: f64, n: usize) {
        self.n += n;
        self.sx += (n as f64) * x;
        self.sx2 += (n as f64) * x * x;
    }

    pub fn add_other(&mut self, other: &Summariser) {
        self.n += other.n;
        self.sx += other.sx;
        self.sx2 += other.sx2;
    }

    pub fn mean(&self) -> f64 {
        self.sx / (self.n as f64)
    }

    pub fn var(&self) -> f64 {
        let m = self.mean();
        self.sx2 / (self.n as f64) - m * m
    }

    pub fn sd(&self) -> f64 {
        if self.n == 0 {
            -1.0
        } else {
            self.var().sqrt()
        }
    }
}

impl Default for Summariser {
    fn default() -> Self {
        Self::new()
    }
}

pub struct MedianOf5 {
    vecs: Vec<Vec<f64>>,
}

impl MedianOf5 {
    pub fn new() -> MedianOf5 {
        MedianOf5 { vecs: Vec::new() }
    }

    pub fn add(&mut self, x0: f64) {
        let mut x = x0;
        let mut i = 0;
        loop {
            while self.vecs.len() <= i {
                self.vecs.push(Vec::new());
            }
            self.vecs[i].push(x);
            if self.vecs[i].len() < 5 {
                break;
            }
            self.vecs[i].sort_by(|a, b| a.partial_cmp(b).unwrap());
            x = self.vecs[i][2];
            self.vecs[i].clear();
            i += 1;
        }
    }

    pub fn median(&self) -> Option<f64> {
        match self.vecs.last() {
            None => None,
            Some(v0) => {
                if v0.len() == 0 {
                    None
                } else {
                    let mut v: Vec<f64> = v0.clone();
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    Some(v[v.len() / 2])
                }
            }
        }
    }
}
