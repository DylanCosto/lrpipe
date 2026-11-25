//! Counting a contig in pieces and appending the pieces' counters gives the same depth as one
//! counter, in every mode, including overlapping mates split between pieces.

use std::path::Path;

use lrpipe_depth::counter::{ContigCounter, Filter, Mode};
use lrpipe_depth::reader::Alignments;

fn check(file: &str, reference: Option<&str>, piece: i64) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let aln = Alignments::open(&root.join(file), reference.map(|r| root.join(r)).as_deref()).unwrap();
    let mut reader = aln.reader(1).unwrap();
    for mode in [Mode::Cigar, Mode::Fast, Mode::Fragment] {
        for filter in [Filter::default(), Filter { exclude_flags: 0, ..Filter::default() }] {
            for (tid, (name, len)) in aln.contigs().iter().enumerate() {
                if !aln.has_records(tid) {
                    continue;
                }
                let len = *len as usize;
                let mut whole = ContigCounter::new(len, mode, filter.clone());
                let mut pieces: Vec<ContigCounter> = Vec::new();
                reader
                    .for_each(tid, false, |a| {
                        whole.push(a);
                        let k = (a.pos.max(0) / piece) as usize;
                        while pieces.len() <= k {
                            pieces.push(ContigCounter::chunk(len, mode, filter.clone()));
                        }
                        pieces[k].push(a);
                    })
                    .unwrap();
                let mut joined = ContigCounter::chunk(len, mode, filter.clone());
                for p in pieces {
                    joined.append(p);
                }
                let (w, j) = (whole.finish(), joined.finish());
                assert_eq!(w.found, j.found, "{file} {name} {mode:?}");
                assert!(w.depth == j.depth, "{file} {name} {mode:?} -F {}: depth differs", filter.exclude_flags);
            }
        }
    }
}

#[test]
fn illumina_pairs_in_pieces() {
    check("ilmn.bam", None, 7_919);
}

#[test]
fn synthetic_edge_cases_in_pieces() {
    for piece in [1, 7, 100, 333] {
        check("synth.bam", None, piece);
        check("synth.cram", Some("synth.fa"), piece);
    }
}

#[test]
fn long_reads_in_pieces() {
    check("ont.bam", None, 100_000);
    check("hifi.bam", None, 250_000);
}
