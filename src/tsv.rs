use std::io::{BufRead, BufReader, Result};

use crate::files::open_reader;

trait AsA {
    fn asa(src: &str) -> Self;
}

impl AsA for String {
    fn asa(src: &str) -> Self {
        src.to_string()
    }
}

impl AsA for u32 {
    fn asa(src: &str) -> Self {
        src.parse().expect("AsA failed")
    }
}

impl AsA for u64 {
    fn asa(src: &str) -> Self {
        src.parse().expect("AsA failed")
    }
}

impl AsA for usize {
    fn asa(src: &str) -> Self {
        src.parse().expect("AsA failed")
    }
}

impl AsA for f64 {
    fn asa(src: &str) -> Self {
        src.parse().expect("AsA failed")
    }
}

pub enum RowImpl {
    Unparsed(String),
    Parsed(Vec<String>)
}

pub struct Row {
    parts: Vec<String>
}

impl Row {
    pub fn new(line: String) -> Result<Row> {
        let parts: Vec<String> = Vec::from_iter(line.split_ascii_whitespace().map(|x| x.to_string()));
        Ok(Row {parts})
    }

    pub fn get(&self, idx: usize) -> &str {
        &self.parts[idx]
    }
}

pub struct Table {
    rows: Vec<Row>
}

impl Table {
    pub fn read_from_tsv(filename: &str) -> Result<Table> {
        let mut boxed = open_reader(filename)?;
        let reader = BufReader::new(boxed.as_mut());
        let mut rows = Vec::new();
        let mut first = true;
        for line_res in reader.lines() {
            let line = line_res?;
            if first {
                first = false;
                continue;
            }
            let row = Row::new(line)?;
            rows.push(row);
        }
        Ok(Table {rows})
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn get(&self, row_idx: usize) -> &Row {
        &self.rows[row_idx]
    }
}