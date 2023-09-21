use polars::prelude::*;
use std::{collections::HashMap, io::Result};

use crate::{files::open_writer, summarise::Summariser};

pub fn peaks(filename: &str, outname: &str) -> Result<()> {
    let mut chroms: HashMap<String, usize> = HashMap::new();
    let mut tbl: Vec<(String, Vec<(usize, usize, usize)>)> = Vec::new();
    {
        let df = CsvReader::from_path(filename)
            .unwrap()
            .with_delimiter(b'\t')
            .finish()
            .unwrap()
            .sort(["chrom", "pos"], false, false)
            .unwrap();
        println!("{:?}", df.shape());
        for i in 0..df.shape().0 {
            let row = df.get_row(i).unwrap().0;
            match (&row[0], &row[1], &row[2], &row[3]) {
                (
                    AnyValue::Utf8(chrom),
                    AnyValue::Int64(pos0),
                    AnyValue::Int64(left_count),
                    AnyValue::Int64(right_count),
                ) => {
                    let chrom_id = match chroms.get(*chrom) {
                        None => {
                            let chrom_id = chroms.len();
                            chroms.insert(chrom.to_string(), chrom_id);
                            tbl.push((chrom.to_string(), Vec::new()));
                            chrom_id
                        }
                        Some(chrom_id) => *chrom_id,
                    };
                    let pos = *pos0 as usize;
                    let lc = *left_count as usize;
                    let rc = *right_count as usize;
                    tbl[chrom_id].1.push((pos, lc, rc));
                }
                _ => {
                    panic!("badly formatted data")
                }
            }
        }
    }
    let mut out_box = open_writer(outname)?;
    let out = out_box.as_mut();
    writeln!(out, "chrom\tpos\tmax_both\tleft_sd\tright_sd\tboth_sd")?;

    let window_max: usize = 100;
    for item in tbl.iter() {
        let chrom_name: &str = &item.0;
        let counts: &[(usize, usize, usize)] = &item.1;

        let n = counts.len();
        let mut i = 0;
        while i < n {
            let mut j = i + 1;
            while j < n && counts[j].0 - counts[j - 1].0 <= window_max {
                j += 1;
            }
            if j == i + 1 {
                // singleton case.

                let p = counts[i].0;
                let l = counts[i].1;
                let r = counts[i].2;
                let both = l + r;
                writeln!(out, "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}", chrom_name, p, both, 0, 0, 0)?;
                i += 1;
            } else {
                let mut max_clip = 0;
                let mut candidate_pos = 0;
                let mut lsum: Summariser = Summariser::new();
                let mut rsum: Summariser = Summariser::new();
                let mut bsum: Summariser = Summariser::new();
                for k in i..j {
                    let p = counts[k].0;
                    let l = counts[k].1;
                    let r = counts[k].2;
                    let both = l + r;
                    if both > max_clip {
                        max_clip = both;
                        candidate_pos = p;
                    }
                    lsum.add_multiple(p as f64, l);
                    rsum.add_multiple(p as f64, r);
                    bsum.add_multiple(p as f64, both);
                }
                writeln!(
                    out,
                    "{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}",
                    chrom_name,
                    candidate_pos,
                    max_clip,
                    lsum.sd(),
                    rsum.sd(),
                    bsum.sd()
                )?;
                i = j;
            }
        }
    }
    Ok(())
}
