use std::io::{Error, ErrorKind, Result};
use std::sync::mpsc::{channel, Receiver, Sender};

use noodles::fastq;
use noodles::fastq::record::Definition;

use crate::files::open_writer;

pub enum FastqStuff {
    Stuff {
        id: String,
        seq: String,
        qual: String,
    },
    Done,
}

pub struct FastqSender {
    tx: Sender<FastqStuff>,
}

impl FastqSender {
    pub fn new(tx: Sender<FastqStuff>) -> FastqSender {
        FastqSender { tx }
    }

    pub fn send(&self, id: &str, seq: &str, qual: &str) -> Result<()> {
        self.tx
            .send(FastqStuff::Stuff {
                id: id.to_string(),
                seq: seq.to_string(),
                qual: qual.to_string(),
            })
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))
    }

    pub fn done(&self) -> Result<()> {
        self.tx
            .send(FastqStuff::Done)
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))
    }
}

pub struct FastqGatherer {
    tx: Sender<FastqStuff>,
    rx: Receiver<FastqStuff>,
    out_name: String,
}

impl FastqGatherer {
    pub fn new(filename: &str) -> Result<FastqGatherer> {
        let (tx, rx) = channel();
        Ok(FastqGatherer {
            tx,
            rx,
            out_name: filename.to_string(),
        })
    }

    pub fn make_sender(&self) -> FastqSender {
        FastqSender::new(self.tx.clone())
    }

    pub fn send(&self, id: &str, seq: &str, qual: &str) -> Result<()> {
        self.tx
            .send(FastqStuff::Stuff {
                id: id.to_string(),
                seq: seq.to_string(),
                qual: qual.to_string(),
            })
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))
    }

    pub fn done(&self) -> Result<()> {
        self.tx
            .send(FastqStuff::Done)
            .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))
    }

    pub fn gather(&self) -> Result<()> {
        let writer = open_writer(&self.out_name)?;
        let mut out = fastq::Writer::new(writer);
        loop {
            let msg = self
                .rx
                .recv()
                .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
            match msg {
                FastqStuff::Stuff { id, seq, qual } => {
                    let rec = fastq::Record::new(Definition::new(id, ""), seq, qual);
                    out.write_record(&rec)?;
                }
                FastqStuff::Done => {
                    break;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn test_simple_1() -> Result<()> {
        let g = FastqGatherer::new("/tmp/test.fastq")?;
        let s = g.make_sender();
        let t = thread::spawn(move || {
            s.send(&"read1", "ACGT", "FFFF").expect("failed to send");
            s.done().expect("done failed");
        });
        g.gather()?;
        t.join().expect("join failed");
        Ok(())
    }

    #[test]
    fn test_simple_2() -> Result<()> {
        let g = FastqGatherer::new("/tmp/test.fastq")?;
        let s = g.make_sender();
        let t = thread::spawn(move || {
            g.gather().expect("gather failed");
        });
        s.send(&"read1", "ACGT", "FFFF")?;
        s.done()?;
        t.join().expect("join failed");
        Ok(())
    }
}
