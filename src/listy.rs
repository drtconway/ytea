use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fmt::{self, Debug, Display},
    marker::PhantomData,
};

#[derive(Clone, Copy, Eq, Default, Hash, PartialEq, Serialize, Deserialize, Debug)]
pub struct ListyElement<T>(u64, PhantomData<T>);

impl<T> Display for ListyElement<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "<{}>", self.0)
    }
}

#[derive(Serialize, Deserialize)]
pub struct Listy<T>
where
    T: std::fmt::Debug,
{
    counter: u64,
    items: HashMap<u64, T>,
    pred: HashMap<u64, u64>,
    succ: HashMap<u64, u64>,
    first: Option<u64>,
    last: Option<u64>,
}

impl<T> Listy<T>
where
    T: std::fmt::Debug,
{
    pub fn new() -> Listy<T> {
        Listy {
            counter: 0,
            items: HashMap::new(),
            pred: HashMap::new(),
            succ: HashMap::new(),
            first: None,
            last: None,
        }
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn front(&self) -> Option<&T> {
        match self.first {
            None => None,
            Some(x) => self.items.get(&x),
        }
    }

    pub fn front_ptr(&self) -> Option<ListyElement<T>> {
        match self.first {
            None => None,
            Some(x) => Some(ListyElement(x, PhantomData)),
        }
    }

    pub fn back(&self) -> Option<&T> {
        match self.last {
            None => None,
            Some(x) => self.items.get(&x),
        }
    }

    pub fn back_ptr(&self) -> Option<ListyElement<T>> {
        match self.last {
            None => None,
            Some(x) => Some(ListyElement(x, PhantomData)),
        }
    }

    pub fn push_front(&mut self, value: T) -> ListyElement<T> {
        let x = self.counter;
        self.counter += 1;
        self.items.insert(x, value);
        match (&self.first, &self.last) {
            (None, None) => {
                self.first = Some(x);
                self.last = Some(x);
            }
            (Some(fst), Some(_lst)) => {
                self.succ.insert(x, *fst);
                self.pred.insert(*fst, x);
                self.first = Some(x);
            }
            _ => unreachable!(),
        }
        ListyElement(x, PhantomData)
    }

    pub fn push_back(&mut self, value: T) -> ListyElement<T> {
        let x = self.counter;
        self.counter += 1;
        self.items.insert(x, value);
        match (&self.first, &self.last) {
            (None, None) => {
                self.first = Some(x);
                self.last = Some(x);
            }
            (Some(_fst), Some(lst)) => {
                self.succ.insert(*lst, x);
                self.pred.insert(x, *lst);
                self.last = Some(x);
            }
            _ => unreachable!(),
        }
        ListyElement(x, PhantomData)
    }

    pub fn pop_front(&mut self) -> Option<T> {
        match (&self.first, &self.last) {
            (None, None) => None,
            (Some(fst), Some(lst)) => {
                let res: u64 = *fst;
                if fst == lst {
                    self.first = None;
                    self.last = None;
                } else {
                    let w = self.succ.get(fst).cloned();
                    match w {
                        Some(nxt) => {
                            self.succ.remove(fst);
                            self.pred.remove(&nxt);
                            self.first = Some(nxt);
                        }
                        None => unreachable!(),
                    }
                }
                self.items.remove(&res)
            }
            _ => unreachable!(),
        }
    }

    pub fn pop_back(&mut self) -> Option<T> {
        match (&self.first, &self.last) {
            (None, None) => None,
            (Some(fst), Some(lst)) => {
                let res: u64 = *lst;
                if fst == lst {
                    self.first = None;
                    self.last = None;
                } else {
                    let w = self.pred.get(lst).cloned();
                    match w {
                        Some(prv) => {
                            self.succ.remove(&prv);
                            self.pred.remove(lst);
                            self.last = Some(prv);
                        }
                        None => unreachable!(),
                    }
                }
                self.items.remove(&res)
            }
            _ => unreachable!(),
        }
    }

    pub fn get(&self, ptr: &ListyElement<T>) -> &T {
        self.items
            .get(&ptr.0)
            .expect("attempt to dereference dead ListyElement")
    }

    pub fn prev(&self, ptr: &ListyElement<T>) -> Option<ListyElement<T>> {
        self.pred.get(&ptr.0).map(|a| ListyElement(*a, PhantomData))
    }

    pub fn next(&self, ptr: &ListyElement<T>) -> Option<ListyElement<T>> {
        self.succ.get(&ptr.0).map(|a| ListyElement(*a, PhantomData))
    }

    pub fn remove(&mut self, ptr: &ListyElement<T>) -> Option<T> {
        let x = ptr.0;
        let res = match (&self.first, &self.last) {
            (None, None) => None,
            (Some(fst), Some(lst)) => {
                match (self.pred.get(&x).cloned(), self.succ.get(&x).cloned()) {
                    (None, None) => {
                        self.first = None;
                        self.last = None;
                    }
                    (Some(w), Some(y)) => {
                        self.pred.remove(&x);
                        self.succ.remove(&x);
                        self.succ.insert(w, y);
                        self.pred.insert(y, w);
                        if x == *fst {
                            self.first = Some(y);
                        }
                        if x == *lst {
                            self.last = Some(w)
                        }
                    }
                    (Some(w), None) => {
                        self.pred.remove(&x);
                        self.succ.remove(&w);
                        self.last = Some(w);
                    }
                    (None, Some(y)) => {
                        self.succ.remove(&x);
                        self.pred.remove(&y);
                        self.first = Some(y);
                    }
                }
                self.items.remove(&x)
            }
            _ => unreachable!(),
        };
        res
    }

    #[allow(dead_code)]
    fn dump(&self) {
        println!("ends: {}, {}", self.first.unwrap_or(999), self.last.unwrap_or(999));
        print!("items: [");
        for i in self.items.iter() {
            print!("<{}>,", i.0);
        }
        println!("]");
        print!("pred: [");
        for e in self.pred.iter() {
            print!("<{}, {}>,", e.0, e.1);
        }
        println!("]");
        print!("succ: [");
        for e in self.succ.iter() {
            print!("<{}, {}>,", e.0, e.1);
        }
        println!("]");
        }

    #[allow(dead_code)]
    fn sanity_check(&self) {
        for e in self.pred.iter() {
            assert_eq!(self.succ.get(e.1), Some(e.0));
            assert!(self.items.contains_key(e.0));
            assert!(self.items.contains_key(e.1));
        }
        for e in self.succ.iter() {
            assert_eq!(self.pred.get(e.1), Some(e.0));
            assert!(self.items.contains_key(e.0));
            assert!(self.items.contains_key(e.1));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cons_1() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        l.push_back(42);
        assert_eq!(l.len(), 1);
        assert_eq!(l.pop_front(), Some(42));
        assert_eq!(l.len(), 0);
    }

    #[test]
    fn test_cons_2() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        l.push_back(42);
        assert_eq!(l.len(), 1);
        assert_eq!(l.pop_back(), Some(42));
        assert_eq!(l.len(), 0);
    }

    #[test]
    fn test_cons_3() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        l.push_front(42);
        assert_eq!(l.len(), 1);
        assert_eq!(l.pop_back(), Some(42));
        assert_eq!(l.len(), 0);
    }

    #[test]
    fn test_cons_4() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        for x in [42, 19, 23] {
            l.push_back(x);
        }
        assert_eq!(l.len(), 3);
        assert_eq!(l.pop_front(), Some(42));
        assert_eq!(l.pop_back(), Some(23));
        assert_eq!(l.pop_front(), Some(19));
        assert_eq!(l.pop_back(), None);
    }

    #[test]
    fn test_remove_1() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        l.push_back(42);
        l.push_back(19);
        let p = l.push_back(23);
        l.push_back(11);
        l.push_back(99);
        assert_eq!(l.len(), 5);
        assert_eq!(l.remove(&p), Some(23));
        assert_eq!(l.len(), 4);
        assert_eq!(l.pop_front(), Some(42));
        assert_eq!(l.pop_front(), Some(19));
        assert_eq!(l.pop_front(), Some(11));
        assert_eq!(l.pop_front(), Some(99));
    }

    #[test]
    fn test_remove_2() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        let p = l.push_back(42);
        l.push_back(19);
        l.push_back(23);
        l.push_back(11);
        l.push_back(99);
        assert_eq!(l.len(), 5);
        assert_eq!(l.remove(&p), Some(42));
        assert_eq!(l.len(), 4);
        assert_eq!(l.pop_front(), Some(19));
        assert_eq!(l.pop_front(), Some(23));
        assert_eq!(l.pop_front(), Some(11));
        assert_eq!(l.pop_front(), Some(99));
    }

    #[test]
    fn test_remove_3() {
        let mut l: Listy<u32> = Listy::new();
        assert_eq!(l.len(), 0);
        l.push_back(42);
        l.push_back(19);
        l.push_back(23);
        l.push_back(11);
        let p = l.push_back(99);
        assert_eq!(l.len(), 5);
        assert_eq!(l.remove(&p), Some(99));
        assert_eq!(l.len(), 4);
        assert_eq!(l.pop_front(), Some(42));
        assert_eq!(l.pop_front(), Some(19));
        assert_eq!(l.pop_front(), Some(23));
        assert_eq!(l.pop_front(), Some(11));
    }
}
