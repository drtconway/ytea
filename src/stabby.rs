use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, VecDeque},
    fmt::Display,
};

use crate::{
    listy::{Listy, ListyElement},
    ransel::Ransel,
    sparse::Sparse,
};

#[derive(
    Clone, Copy, Eq, PartialOrd, Ord, Default, Hash, PartialEq, Debug, Serialize, Deserialize,
)]
pub struct Interval {
    first: u32,
    last: u32,
}

impl Interval {
    pub fn new(first: u32, last: u32) -> Interval {
        Interval {
            first: first,
            last: last,
        }
    }

    pub fn zero() -> Interval {
        Interval { first: 0, last: 0 }
    }
}

impl Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}, {}]", self.first, self.last)
    }
}

pub fn make_smaller(items: &[Interval]) -> (Vec<Interval>, HashMap<Interval, Vec<Interval>>) {
    let n = items.len();
    let mut i: usize = 0;
    let mut basic: Vec<Interval> = Vec::new();
    let mut smaller: HashMap<Interval, Vec<Interval>> = HashMap::new();
    while i < n {
        let mut j: usize = i;
        while j + 1 < n && items[j].first == items[j + 1].first {
            j += 1;
        }
        if j > i {
            let v = &items[i..j];
            if v.len() > 0 {
                smaller.insert(items[j], Vec::from(v));
            }
        }
        basic.push(items[j]);
        i = j + 1;
    }
    (basic, smaller)
}

#[derive(Serialize, Deserialize)]
pub struct DenseStabby {
    smaller: HashMap<Interval, Vec<Interval>>,
    start: HashMap<u32, Interval>,
    start2: HashMap<u32, Interval>,
    parent: HashMap<Interval, Interval>,
    last: HashMap<Interval, Interval>,
    left: HashMap<Interval, Interval>,
}

impl DenseStabby {
    pub fn new(q_max: u32, items: &[Interval]) -> DenseStabby {
        let (basic, smaller) = make_smaller(items);
        let mut event: HashMap<u32, Vec<Interval>> = HashMap::new();
        for item in basic.iter() {
            event.entry(item.last).or_insert(Vec::new()).push(*item);
            event.entry(item.first).or_insert(Vec::new()).push(*item);
        }

        let mut start: HashMap<u32, Interval> = HashMap::new();
        let mut start2: HashMap<u32, Interval> = HashMap::new();
        let mut parent: HashMap<Interval, Interval> = HashMap::new();
        let mut last: HashMap<Interval, Interval> = HashMap::new();
        let mut left: HashMap<Interval, Interval> = HashMap::new();

        let mut l: Listy<Interval> = Listy::new();
        let mut saved: HashMap<Interval, ListyElement<Interval>> = HashMap::new();
        let mut rml: usize = 0;

        for q in 0..=q_max {
            match l.back() {
                None => {}
                Some(a) => {
                    start.insert(q, *a);
                }
            }
            match event.get(&q) {
                None => {}
                Some(v) => {
                    for a in v.iter().rev() {
                        match saved.get(a) {
                            None => {
                                if a.first == q {
                                    start.insert(q, *a);
                                    let ptr = l.push_back(*a);
                                    saved.insert(*a, ptr);
                                } else {
                                    // TODO
                                    println!("it happened");
                                }
                            }
                            Some(s) => {
                                assert_eq!(a.last, q);
                                let op = l.prev(s).map(|v| {   
                                l.get(&v)});
                                let mut p: Interval = Interval::zero();
                                match op {
                                    None => {}
                                    Some(pp) => {
                                        p = *pp;
                                    }
                                };
                                parent.insert(*a, p);
                                match last.get(&p) {
                                    None => {}
                                    Some(w) => {
                                        left.insert(*a, *w);
                                    }
                                }
                                last.insert(p, *a);
                                l.remove(s);
                                let mut s_count = 0;
                                for e in saved.iter() {
                                    if e.1 == s {
                                        s_count += 1;
                                        assert_eq!(e.0, a);
                                    }
                                }
                                assert_eq!(s_count, 1);
                                saved.remove(a);
                            }
                        }
                    }

                    while rml + 1 < basic.len() && basic[rml + 1].first <= q {
                        rml += 1;
                    }
                    if basic[rml].first <= q {
                        start2.insert(q, basic[rml]);
                    }
                }
            }
        }

        DenseStabby {
            smaller: smaller,
            start: start,
            start2: start2,
            parent: parent,
            last: last,
            left: left,
        }
    }

    pub fn stabs(&self, q: u32) -> bool {
        self.start.contains_key(&q)
    }

    pub fn stab(&self, q: u32) -> Vec<Interval> {
        let mut res: Vec<Interval> = Vec::new();
        let mut kew: VecDeque<Interval> = VecDeque::new();
        let mut ov: Option<&Interval> = self.start.get(&q);
        loop {
            match ov {
                None => {
                    break;
                }
                Some(v) => {
                    if *v == Interval::zero() {
                        break;
                    }
                    kew.push_front(*v);
                    ov = self.parent.get(v);
                }
            }
        }
        loop {
            match kew.pop_back() {
                None => {
                    break;
                }
                Some(a) => {
                    res.push(a);
                    match self.smaller.get(&a) {
                        None => {}
                        Some(s) => {
                            for r in s.iter().rev() {
                                if r.last < q {
                                    break;
                                }
                                res.push(*r);
                            }
                        }
                    }
                    let mut ot = self.left.get(&a);
                    loop {
                        match ot {
                            None => {
                                break;
                            }
                            Some(t) => {
                                if t.last < q {
                                    break;
                                }
                                kew.push_back(*t);
                                ot = self.last.get(t);
                            }
                        }
                    }
                }
            }
        }
        res.reverse();
        res
    }

    pub fn stab_interval(&self, qi: &Interval) -> Vec<Interval> {
        let lq = qi.first;
        let rq = qi.last;

        let mut ot: Option<&Interval> = None;
        match self.start.get(&lq) {
            None => {}
            Some(u) => {
                ot = Some(u);
            }
        }
        match self.start2.get(&rq) {
            None => {}
            Some(u) => match ot {
                None => {
                    ot = Some(u);
                }
                Some(t) => {
                    if t.first < u.first {
                        ot = Some(u);
                    }
                }
            },
        }

        let mut res: Vec<Interval> = Vec::new();
        match ot {
            None => {
                return res;
            }
            Some(t) => {
                if t.last < lq {
                    return res;
                }
            }
        }

        let mut kew: VecDeque<Interval> = VecDeque::new();
        loop {
            match ot {
                None => {
                    break;
                }
                Some(t) => {
                    if *t == Interval::zero() {
                        break;
                    }
                    kew.push_front(*t);
                    ot = self.parent.get(&t);
                }
            }
        }

        loop {
            match kew.pop_back() {
                None => {
                    break;
                }
                Some(a) => {
                    res.push(a);

                    match self.smaller.get(&a) {
                        None => {}
                        Some(s) => {
                            for r in s.iter().rev() {
                                if r.last < lq {
                                    break;
                                }
                                res.push(*r);
                            }
                        }
                    }

                    ot = self.left.get(&a);

                    loop {
                        match ot {
                            None => {
                                break;
                            }
                            Some(t) => {
                                if t.last < lq {
                                    break;
                                }
                                kew.push_back(*t);
                                ot = self.last.get(t);
                            }
                        }
                    }
                }
            }
        }
        res.reverse();
        res
    }
}

pub struct Stabby {
    domain: Sparse,
    dense: DenseStabby,
}

impl Stabby {
    pub fn new(xs: &[Interval]) -> Stabby {
        let domain = Self::make_domain(xs);
        let mut ys: Vec<Interval> = Vec::new();
        let mut y_max = 0;
        for x in xs.iter() {
            let (y_f,y_l) = domain.rank2(x.first, x.last);
            let y = Interval::new((y_f * 2) as u32, (y_l * 2) as u32);
            if y.last > y_max {
                y_max = y.last;
            }
            ys.push(y);
        }
        ys.sort();
        let dense = DenseStabby::new(y_max + 1, &ys);

        Stabby { domain, dense }
    }

    fn make_domain(xs: &[Interval]) -> Sparse {
        let mut ys: Vec<u32> = Vec::new();
        ys.push(0);
        for x in xs.iter() {
            ys.push(x.first);
            ys.push(x.last);
        }
        ys.sort();
        ys.dedup();
        Sparse::new(&ys)
    }

    pub fn stab(&self, q: u32) -> Vec<Interval> {
        let qd = self.sparse_to_dense(q);
        let ys = self.dense.stab(qd);
        let mut xs: Vec<Interval> = Vec::new();
        for y in ys.iter() {
            let y1 = self.dense_to_sparse(y.first);
            let y2 = self.dense_to_sparse(y.last);
            xs.push(Interval::new(y1, y2));
        }
        xs
    }

    pub fn stab_interval(&self, q: &Interval) -> Vec<Interval> {
        let qd = Interval::new(self.sparse_to_dense(q.first), self.sparse_to_dense(q.last));
        let ys = self.dense.stab_interval(&qd);
        let mut xs: Vec<Interval> = Vec::new();
        for y in ys.iter() {
            let y1 = self.dense_to_sparse(y.first);
            let y2 = self.dense_to_sparse(y.last);
            xs.push(Interval::new(y1, y2));
        }
        xs
    }

    fn sparse_to_dense(&self, x: u32) -> u32 {
        let (r1, r2) = self.domain.rank2(x, x + 1);
        let q = if r1 != r2 { r1 * 2 } else { r1 * 2 - 1 };
        q as u32
    }

    fn dense_to_sparse(&self, x: u32) -> u32 {
        self.domain.select((x / 2) as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_smaller() {
        let items: [Interval; 3] = [
            Interval { first: 1, last: 3 },
            Interval { first: 2, last: 5 },
            Interval { first: 2, last: 7 },
        ];
        let (basic, _smaller) = make_smaller(&items);
        assert_eq!(basic.len(), 2);
    }

    #[test]
    fn test_stabby_1() {
        let src: Vec<Interval> = Vec::from([
            Interval::new(1, 2),
            Interval::new(1, 4),
            Interval::new(3, 5),
        ]);
        let s = DenseStabby::new(5, &src);
        println!("s={}", serde_yaml::to_string(&s).unwrap());
        let w = Interval::new(1, 4);
        let v = Vec::from([Interval::new(1, 2)]);
        assert_eq!(s.smaller.len(), 1);
        assert_eq!(s.smaller.get(&w), Some(&v));
        assert_eq!(s.stabs(3), true);
        assert_eq!(s.stabs(0), false);
        assert_eq!(
            s.stab(3),
            Vec::from([Interval::new(1, 4), Interval::new(3, 5)])
        );
    }

    #[test]
    fn test_stabby_2() {
        let mut src: Vec<Interval> = Vec::from([
            Interval::new(7, 62),
            Interval::new(38, 57),
            Interval::new(22, 65),
            Interval::new(5, 88),
            Interval::new(23, 43),
            Interval::new(73, 75),
            Interval::new(36, 89),
            Interval::new(31, 55),
            Interval::new(17, 91),
            Interval::new(24, 80),
            Interval::new(40, 47),
            Interval::new(52, 67),
            Interval::new(16, 45),
            Interval::new(86, 90),
            Interval::new(46, 56),
            Interval::new(33, 84),
            Interval::new(21, 56),
            Interval::new(15, 93),
            Interval::new(48, 51),
            Interval::new(68, 96),
            Interval::new(35, 64),
            Interval::new(63, 78),
            Interval::new(71, 99),
            Interval::new(3, 4),
            Interval::new(18, 85),
            Interval::new(37, 83),
            Interval::new(10, 25),
            Interval::new(32, 66),
            Interval::new(9, 49),
            Interval::new(3, 13),
            Interval::new(6, 94),
            Interval::new(75, 98),
            Interval::new(8, 69),
            Interval::new(39, 74),
            Interval::new(11, 41),
            Interval::new(14, 59),
            Interval::new(24, 75),
            Interval::new(1, 2),
            Interval::new(50, 53),
            Interval::new(0, 42),
            Interval::new(49, 76),
            Interval::new(25, 55),
            Interval::new(60, 69),
            Interval::new(15, 18),
            Interval::new(52, 79),
            Interval::new(15, 74),
            Interval::new(30, 82),
            Interval::new(34, 81),
            Interval::new(39, 97),
            Interval::new(51, 81),
            Interval::new(29, 87),
            Interval::new(61, 92),
            Interval::new(12, 70),
            Interval::new(26, 54),
            Interval::new(12, 52),
            Interval::new(52, 54),
            Interval::new(12, 72),
            Interval::new(12, 58),
            Interval::new(44, 88),
        ]);
        src.sort();
        let s = DenseStabby::new(100, &src);
        for q in 0..=100 {
            let mut expected: Vec<Interval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stabs(q), expected.len() > 0);
            assert_eq!(s.stab(q), expected);
        }
    }

    #[test]
    fn test_stabby_3() {
        let mut src: Vec<Interval> = Vec::from([
            Interval::new(7, 62),
            Interval::new(38, 57),
            Interval::new(22, 65),
            Interval::new(5, 88),
            Interval::new(23, 43),
            Interval::new(73, 75),
            Interval::new(36, 89),
            Interval::new(31, 55),
            Interval::new(17, 91),
            Interval::new(24, 80),
            Interval::new(40, 47),
            Interval::new(52, 67),
            Interval::new(16, 45),
            Interval::new(86, 90),
            Interval::new(46, 56),
            Interval::new(33, 84),
            Interval::new(21, 56),
            Interval::new(15, 93),
            Interval::new(48, 51),
            Interval::new(68, 96),
            Interval::new(35, 64),
            Interval::new(63, 78),
            Interval::new(71, 99),
            Interval::new(3, 4),
            Interval::new(18, 85),
            Interval::new(37, 83),
            Interval::new(10, 25),
            Interval::new(32, 66),
            Interval::new(9, 49),
            Interval::new(3, 13),
            Interval::new(6, 94),
            Interval::new(75, 98),
            Interval::new(8, 69),
            Interval::new(39, 74),
            Interval::new(11, 41),
            Interval::new(14, 59),
            Interval::new(24, 75),
            Interval::new(1, 2),
            Interval::new(50, 53),
            Interval::new(0, 42),
            Interval::new(49, 76),
            Interval::new(25, 55),
            Interval::new(60, 69),
            Interval::new(15, 18),
            Interval::new(52, 79),
            Interval::new(15, 74),
            Interval::new(30, 82),
            Interval::new(34, 81),
            Interval::new(39, 97),
            Interval::new(51, 81),
            Interval::new(29, 87),
            Interval::new(61, 92),
            Interval::new(12, 70),
            Interval::new(26, 54),
            Interval::new(12, 52),
            Interval::new(52, 54),
            Interval::new(12, 72),
            Interval::new(12, 58),
            Interval::new(44, 88),
        ]);
        src.sort();
        let s = DenseStabby::new(100, &src);

        for qi in [
            Interval::new(39, 84),
            Interval::new(13, 44),
            Interval::new(2, 57),
            Interval::new(50, 75),
            Interval::new(2, 11),
        ] {
            let mut expected: Vec<Interval> = Vec::new();
            for ivl in src.iter() {
                if (ivl.first <= qi.last) && (ivl.last >= qi.first) {
                    expected.push(*ivl);
                }
            }
            println!("[{}, {}] -> {}", qi.first, qi.last, expected.len());
            assert_eq!(s.stab_interval(&qi), expected);
        }
    }

    #[test]
    fn test_stabby_4a() {
        let src: Vec<Interval> = Vec::from([
            Interval::new(0, 10),
            Interval::new(0, 58),
            Interval::new(2, 10),
            Interval::new(2, 44),
            Interval::new(4, 10),
            Interval::new(4, 58),
            Interval::new(6, 10),
            Interval::new(6, 58),
            Interval::new(8, 10),
            Interval::new(8, 36),
            Interval::new(12, 14),
            Interval::new(16, 18),
            Interval::new(20, 22),
            Interval::new(24, 26),
            Interval::new(28, 30),
            Interval::new(32, 34),
            Interval::new(32, 36),
            Interval::new(38, 40),
            Interval::new(42, 44),
            Interval::new(46, 50),
            Interval::new(48, 50),
            Interval::new(52, 54),
            Interval::new(56, 58),
            Interval::new(60, 64),
            Interval::new(60, 134),
            Interval::new(60, 136),
            Interval::new(60, 142),
            Interval::new(62, 64),
            Interval::new(62, 132),
            Interval::new(66, 68),
            Interval::new(70, 72),
            Interval::new(74, 76),
            Interval::new(78, 80),
            Interval::new(82, 84),
            Interval::new(86, 88),
            Interval::new(90, 92),
            Interval::new(94, 96),
            Interval::new(98, 100),
            Interval::new(102, 104),
            Interval::new(106, 108),
            Interval::new(110, 112),
            Interval::new(114, 116),
            Interval::new(118, 120),
            Interval::new(122, 124),
            Interval::new(126, 128),
            Interval::new(130, 132),
            Interval::new(130, 134),
            Interval::new(130, 136),
            Interval::new(138, 140),
            Interval::new(144, 146),
            Interval::new(144, 258),
            Interval::new(144, 260),
            Interval::new(148, 150),
            Interval::new(152, 154),
            Interval::new(156, 158),
            Interval::new(160, 162),
            Interval::new(164, 166),
            Interval::new(168, 170),
            Interval::new(172, 174),
            Interval::new(176, 178),
            Interval::new(180, 182),
            Interval::new(184, 186),
            Interval::new(188, 190),
            Interval::new(192, 194),
            Interval::new(196, 198),
            Interval::new(200, 202),
            Interval::new(204, 206),
            Interval::new(208, 210),
            Interval::new(212, 214),
            Interval::new(216, 218),
            Interval::new(220, 222),
            Interval::new(224, 226),
            Interval::new(228, 230),
            Interval::new(232, 234),
            Interval::new(236, 238),
            Interval::new(240, 244),
            Interval::new(242, 244),
            Interval::new(242, 256),
            Interval::new(246, 248),
            Interval::new(250, 252),
            Interval::new(254, 256),
            Interval::new(254, 258),
        ]);

        let s = DenseStabby::new(259, &src);

        for q in 0..259 {
            let mut expected: Vec<Interval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stab(q), expected);
        }
    }

    #[test]
    fn test_stabby_4() {
        let src: Vec<Interval> = Vec::from([
            Interval::new(55519604, 55519916),
            Interval::new(55519604, 55545079),
            Interval::new(55519646, 55519916),
            Interval::new(55519646, 55539543),
            Interval::new(55519690, 55519916),
            Interval::new(55519690, 55545079),
            Interval::new(55519718, 55519916),
            Interval::new(55519718, 55545079),
            Interval::new(55519731, 55519916),
            Interval::new(55519731, 55535914),
            Interval::new(55520405, 55520479),
            Interval::new(55522102, 55522166),
            Interval::new(55523721, 55523822),
            Interval::new(55528878, 55528992),
            Interval::new(55533873, 55533960),
            Interval::new(55535712, 55535763),
            Interval::new(55535712, 55535914),
            Interval::new(55537483, 55537585),
            Interval::new(55538765, 55539543),
            Interval::new(55543938, 55544074),
            Interval::new(55544025, 55544074),
            Interval::new(55544220, 55544369),
            Interval::new(55544907, 55545079),
            Interval::new(55547292, 55550006),
            Interval::new(55547292, 55617614),
            Interval::new(55547292, 55617622),
            Interval::new(55547292, 55618880),
            Interval::new(55548378, 55550006),
            Interval::new(55548378, 55617608),
            Interval::new(55558775, 55558968),
            Interval::new(55564313, 55564497),
            Interval::new(55564902, 55565041),
            Interval::new(55565703, 55565850),
            Interval::new(55568194, 55568363),
            Interval::new(55573619, 55573777),
            Interval::new(55577315, 55577356),
            Interval::new(55578247, 55578342),
            Interval::new(55579679, 55579781),
            Interval::new(55581567, 55581698),
            Interval::new(55585051, 55585167),
            Interval::new(55586618, 55586734),
            Interval::new(55588879, 55588956),
            Interval::new(55598416, 55599039),
            Interval::new(55603978, 55604076),
            Interval::new(55615451, 55615506),
            Interval::new(55617144, 55617608),
            Interval::new(55617144, 55617614),
            Interval::new(55617144, 55617622),
            Interval::new(55617869, 55618444),
            Interval::new(55634061, 55636392),
            Interval::new(55634061, 55693844),
            Interval::new(55634061, 55693863),
            Interval::new(55637552, 55637599),
            Interval::new(55640627, 55640705),
            Interval::new(55643158, 55643213),
            Interval::new(55643319, 55643425),
            Interval::new(55644637, 55644720),
            Interval::new(55645349, 55645432),
            Interval::new(55646259, 55646322),
            Interval::new(55646415, 55646486),
            Interval::new(55647347, 55647453),
            Interval::new(55654900, 55654953),
            Interval::new(55656131, 55656220),
            Interval::new(55656305, 55656371),
            Interval::new(55660157, 55660193),
            Interval::new(55661956, 55662026),
            Interval::new(55666991, 55667093),
            Interval::new(55667862, 55667958),
            Interval::new(55671319, 55671376),
            Interval::new(55671995, 55672046),
            Interval::new(55672893, 55673079),
            Interval::new(55679682, 55679795),
            Interval::new(55680712, 55680759),
            Interval::new(55680855, 55680918),
            Interval::new(55683785, 55683834),
            Interval::new(55684943, 55685048),
            Interval::new(55684994, 55685048),
            Interval::new(55684994, 55693823),
            Interval::new(55686370, 55686444),
            Interval::new(55687645, 55687705),
            Interval::new(55693663, 55693823),
            Interval::new(55693663, 55693844),
        ]);
        let s: Stabby = Stabby::new(&src);
        let mut q = 55519000;
        while q < 55700000 {
            let mut expected: Vec<Interval> = Vec::new();
            for ivl in src.iter() {
                if ivl.first <= q && q <= ivl.last {
                    expected.push(*ivl);
                }
            }
            assert_eq!(s.stab(q), expected);
            q += 150;
        }

    }
}
