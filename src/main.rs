use docopt::Docopt;
use y_tea::discordant::get_disc;
use std::io::{Error, ErrorKind, Result};
use y_tea::clips::get_clips_2;
use y_tea::files::open_writer;
use y_tea::info::BasicInfo;
use y_tea::peaks::peaks;
use y_tea::regions::regions_from_bed;
use y_tea::repeats::prepare_repeats;

const USAGE: &'static str = "
Usage: yTea clips [options] <bamfile> <clips-fastq> <clips-summary>
       yTea peaks [options] <clips-summary> <peaks-summary>
       yTea disc [options] <bamfile> <peaks-summary>
       yTea info <bamfile> <info-file>
       yTea prepare repeats [options] <repeatdb> <reference-genome> <repeat-fasta>

Options:
    -h                      Show this help message.
    -v                      Produce verbose output.
    -R BED                  BED file with repeat regions.
    -C CLASS                Process only the given class of repeat elements (use '_LIST_' to see a list and exit).
    --repeat-match REGEX    Process only repeats with labels matching the given regular expression.
    -W NUM                  Window size. [default: 50]
";

fn main() -> Result<()> {
    let args = Docopt::new(USAGE)
        .and_then(|dopt| dopt.parse())
        .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);

    let verbose = args.get_bool("-v");

    if args.get_bool("prepare") {
        if args.get_bool("repeats") {
            let rpt = args.get_str("<repeatdb>");
            let genome = args.get_str("<reference-genome>");
            let dest_name = args.get_str("<repeat-fasta>");
            let opt_window = if args.get_str("-W").len() == 0 {
                None
            } else {
                let window_str = args.get_str("-W");
                let window: usize = window_str
                    .parse::<usize>()
                    .map_err(|err| Error::new(ErrorKind::Other, err.to_string()))?;
                Some(window)
            };
            let opt_class: Option<String> = if args.get_str("-C").len() > 0 {
                Some(args.get_str("-C").to_string())
            } else {
                None
            };
            let opt_filter: Option<String> = if args.get_str("--repeat-match").len() > 0 {
                Some(args.get_str("--repeat-match").to_string())
            } else {
                None
            };

            prepare_repeats(rpt, genome, dest_name, &opt_window, opt_class, opt_filter)?;
        }
    }
    if args.get_bool("info") {
        let bam = args.get_str("<bamfile>");
        let out_name = args.get_str("<info-file>");
        let info = BasicInfo::from_bam(bam)?;
        let mut out = open_writer(out_name)?;
        write!(out.as_mut(), "{}", serde_json::to_string(&info).unwrap())?;
    }
    if args.get_bool("clips") {
        let bam = args.get_str("<bamfile>");
        let fq = args.get_str("<clips-fastq>");
        let tbl = args.get_str("<clips-summary>");
        let reps = (if args.get_str("-R").len() == 0 {
            None
        } else {
            Some(args.get_str("-R"))
        })
        .map(|path| {
            let res = regions_from_bed(path).expect("failed to read repeats");
            res
        });

        get_clips_2(bam, &Some(fq), &Some(tbl), &reps, verbose)?;
    }
    if args.get_bool("peaks") {
        let clips_path = args.get_str("<clips-summary>");
        let peaks_path = args.get_str("<peaks-summary>");
        peaks(clips_path, peaks_path)?;
    }
    if args.get_bool("disc") {
        let bam = args.get_str("<bamfile>");
        let peaks_path = args.get_str("<peaks-summary>");
        get_disc(bam, peaks_path, &None, false)?
    }
    Ok(())
}
