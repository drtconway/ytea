use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hasher},
    mem::swap,
};

pub struct Feistel {
    split: usize,
    mask: u64,
    keys: Vec<u64>,
}

impl Feistel {
    pub fn new(bits: usize, keys: Vec<u64>) -> Feistel {
        let split = bits / 2;
        let mask: u64 = (1 << split) - 1;
        Feistel { split, mask, keys }
    }

    pub fn encode(&self, x: u64) -> u64 {
        let n = self.keys.len();
        let mut x_l = x >> self.split;
        let mut x_r = x & self.mask;
        for i in 0..n {
            x_l ^= Feistel::hash(x_r, self.keys[i]) & self.mask;
            swap(&mut x_l, &mut x_r);
        }
        (x_r << self.split) | x_l
    }

    pub fn decode(&self, x: u64) -> u64 {
        let n = self.keys.len();
        let mut x_l = x >> self.split;
        let mut x_r = x & self.mask;
        for j in 0..n {
            let i = n - 1 - j;
            x_l ^= Feistel::hash(x_r, self.keys[i]) & self.mask;
            swap(&mut x_l, &mut x_r);
        }
        (x_r << self.split) | x_l
    }

    pub fn hash(x: u64, k: u64) -> u64 {
        let mut h = DefaultHasher::new();
        h.write_u64(k);
        h.write_u64(x);
        h.finish()
    }
}

pub struct Perm {
    n: u64,
    i: u64,
    f: Feistel,
}

impl Perm {
    pub fn new(n: u64, seed: u64) -> Perm {
        let m = n.next_power_of_two();
        let mut bits: usize = m.ilog2() as usize;
        if (bits & 1) == 1 {
            bits += 1;
        }
        let mut keys: Vec<u64> = Vec::new();
        for i in 0..4 {
            keys.push(Feistel::hash(i, seed))
        }
        Perm {
            n,
            i: 0,
            f: Feistel::new(bits, keys),
        }
    }
}

impl Iterator for Perm {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        if self.i == self.n {
            return None
        }
        let mut r = self.f.encode(self.i);
        while r >= self.n {
            r = self.f.encode(r);
        }
        self.i += 1;
        Some(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feistel_1() {
        let keys: Vec<u64> = Vec::from([0xfaf6d0f4374b8b1d, 0xc3c309110ce06759]);

        let f = Feistel::new(64, keys);
        {
            let x = 23;
            let y = f.encode(x);
            let z = f.decode(y);
            assert_eq!(x, z);
        }
    }

    #[test]
    fn test_feistel_1b() {
        let keys: Vec<u64> = Vec::from([0xfaf6d0f4374b8b1d, 0xc3c309110ce06759]);

        let f = Feistel::new(8, keys);
        {
            let x = 23;
            let y = f.encode(x);
            let z = f.decode(y);
            assert_eq!(x, z);
        }
    }

    #[test]
    fn test_feistel_2() {
        let keys: Vec<u64> = Vec::from([
            0xfaf6d0f4374b8b1d,
            0xc3c309110ce06759,
            0xca7c868a85660aff,
            0x4836e41c8df89cba,
            0x70e29455c9a1d91f,
            0xa6b9dd427a3f2a5b,
            0x0936dc7d62364b4d,
            0xcf6ceab6a74075a7,
            0xe0ba03a971e91f8a,
            0x728bbbc088e9bde7,
        ]);

        let f = Feistel::new(64, keys);
        for x in [
            0xe3977ad97bfc59b4,
            0xeb03b9543cd24c03,
            0x0b68c4cc54308296,
            0x5c2d70c160fd99ff,
            0x882216cf10f46983,
            0x6ebed7b4610dac51,
            0x1773b63d0386be9a,
            0xc29ecb14441c7964,
            0x4ee37f7b4178b78a,
            0xb1591652b76312a7,
            0xd656e65291ad0e3c,
            0xc58a7051dc90c3b3,
            0x98befbfdd6807101,
            0xf5a2f59f2d7d3b24,
            0x4473a4a6d62af064,
            0x9283b92a0f3bc00b,
            0xd7edb18740732bb2,
            0x393c10fe65928e9a,
            0xee2fe63dd35c3caf,
            0xae8417547a24ebdf,
            0x05aca77c57534dfd,
            0x049c06833075d5f4,
            0x261e9efb79ed9dc3,
            0xd8a2889ad606d3a0,
            0x967ace458658738e,
            0x858f56a220f86fac,
            0x6399819b954a6052,
            0x0b0acca02b3a0961,
            0x64470eae55038c89,
            0x747a92c9d824d965,
            0x0c9015d960297149,
            0x67fbb352cae3580e,
            0xdbce139c58a0817b,
            0x5baa688af81b1b9c,
            0x3c5ef0a9ae0ddaa9,
            0x0c4997a4baee48c0,
            0xd83321ca293a594a,
            0x5e283f9b025446f7,
            0xda9cdfc81253bf44,
            0xb8a8595b20a3dd00,
            0x473d52fe6fae93b2,
            0x7dea6d81c33f5237,
            0x7b2910d5f1c1605e,
            0x5bb6062d44b09e48,
            0xe7c8bc4aeab169b6,
            0x2725a3892f234c89,
            0xa847b0dcf0503b8d,
            0x8c6f07202eae1278,
            0x558cad385e4e7200,
            0x505f936cb98a17b7,
            0x38ecf09e3dbc4071,
            0xa0783fc7fae6c7e6,
            0x2f331be607e9aac9,
            0x47b481edd85410fa,
            0x6e80fe65240cf477,
            0x54ffce7298ec0ff9,
            0x5a2c1571bc1263cd,
            0x73a7459e0f26cf29,
            0x0cba57c4f4e47ffe,
            0xf2f35d37f3711a6c,
            0x259ec9b0e1b0fb52,
            0xad76da625193a631,
            0xa14ce2331b61af1c,
            0x0b23d888ac42ab98,
            0xe84036d8fd16457f,
            0x173650eb7cc58d65,
            0xac77916df86d6333,
            0x8529f39627c8845a,
            0x0eb0cb232ef681fa,
            0x79f24c5ebf1caf41,
            0xb34d303338ffa660,
            0x181633d376e6afa0,
            0xed017757d0e16c81,
            0x49cb962cf0208363,
            0x64cb4c3da1b67a7e,
            0xe8be0c632540a93e,
            0x443662ec90a09bd3,
            0x36ea4042bb167a4d,
            0x90fd21192ca343a3,
            0xfef1db89b8518577,
            0xec404573141dd6c9,
            0xc34919f29e1de58f,
            0xf71a293157d931ff,
            0x814cde3a4124b983,
            0x5e0146b575817810,
            0x8282dccb87f67ba2,
            0x62defd8874b1f449,
            0x0d63d5c0b5a2416e,
            0x8c069117cc60b481,
            0xfc5782298e092208,
            0x1f3abe690bf72112,
            0x27fbced64dd35d86,
            0x6e812de10d16c521,
            0x2af94e1bc28de6c0,
            0xec87035a20211731,
            0x3649331ad954c8d7,
            0x5a392ae7f4257070,
            0x2eab129ed3809a62,
            0x271a5eb46d918081,
            0x696997d1c21f0d40,
        ] {
            let y = f.encode(x);
            let z = f.decode(y);
            assert_eq!(x, z);
        }
    }

    use std::collections::HashSet;

    #[test]
    fn test_perm_1() {
        let n = 25;
        let mut seen = HashSet::new();
        for x in Perm::new(n, 19) {
            assert!(x < n);
            assert!(!seen.contains(&x));
            seen.insert(x);
        }
        assert_eq!(seen.len() as u64, n);
    }

    #[test]
    fn test_perm_2() {
        let n = 2500;
        let mut seen = HashSet::new();
        for x in Perm::new(n, 19) {
            assert!(x < n);
            assert!(!seen.contains(&x));
            seen.insert(x);
        }
        assert_eq!(seen.len() as u64, n);
    }
}