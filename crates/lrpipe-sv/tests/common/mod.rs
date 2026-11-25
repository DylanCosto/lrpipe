//! Helpers shared by the integration tests.

use std::{
    fs,
    io::{BufWriter, Write},
    path::Path,
};

/// Rebuilds the fixture contig from the stored window: N up to the window, the window, then N to
/// the contig length. The window's name is `CONTIG:START-END LN:LENGTH` (1-based, inclusive).
/// Each of `extra` (name, length) follows as another contig of arbitrary bases (no N).
pub fn write_reference(window: &Path, out: &Path, extra: &[(&str, usize)]) {
    let text = fs::read_to_string(window).expect("reading the reference window");
    let mut lines = text.lines();
    let header = lines.next().and_then(|l| l.strip_prefix('>')).expect("FASTA header");
    let (name, length) = header.split_once(' ').expect("window name and length");
    let (contig, span) = name.split_once(':').expect("CONTIG:START-END");
    let start: usize = span.split('-').next().and_then(|s| s.parse().ok()).expect("window start");
    let length: usize = length.strip_prefix("LN:").and_then(|s| s.parse().ok()).expect("LN:LENGTH");
    let seq: String = lines.collect();

    let mut f = BufWriter::new(fs::File::create(out).expect("creating the reference"));
    writeln!(f, ">{contig}").unwrap();
    write_n(&mut f, start - 1);
    f.write_all(seq.as_bytes()).unwrap();
    write_n(&mut f, length - (start - 1) - seq.len());
    writeln!(f).unwrap();
    // One line per contig, so each has a single index entry.
    let mut offset = contig.len() + 2;
    let mut fai = format!("{contig}\t{length}\t{offset}\t{length}\t{}\n", length + 1);
    offset += length + 1;
    let mut state: u32 = 1;
    for (name, length) in extra {
        writeln!(f, ">{name}").unwrap();
        offset += name.len() + 2;
        let bases: Vec<u8> = (0..*length)
            .map(|_| {
                state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(state >> 16) as usize % 4]
            })
            .collect();
        f.write_all(&bases).unwrap();
        writeln!(f).unwrap();
        fai.push_str(&format!("{name}\t{length}\t{offset}\t{length}\t{}\n", length + 1));
        offset += length + 1;
    }
    f.flush().unwrap();
    drop(f);
    fs::write(out.with_extension("fa.fai"), fai).expect("writing the reference index");
}

/// Compresses the FASTA `plain` with bgzip into `out`, beside a copy of the plain one's `.fai` (its
/// offsets are into the uncompressed text) and a `.gzi`, as `bgzip -i` and `samtools faidx` make
/// them.
#[allow(dead_code, reason = "not every test that includes this module bgzips a reference")]
pub fn bgzip_reference(plain: &Path, out: &Path) {
    // bgzip's block size.
    const BLOCK: usize = 0xff00;
    let text = fs::read(plain).expect("reading the reference");
    let mut w = pasta::bgzf::io::Writer::new(BufWriter::new(fs::File::create(out).expect("creating the copy")));
    let mut blocks = Vec::new();
    for (i, chunk) in text.chunks(BLOCK).enumerate() {
        if i > 0 {
            blocks.push((w.position(), (i * BLOCK) as u64));
        }
        w.write_all(chunk).unwrap();
        w.flush().unwrap();
    }
    w.finish().expect("compressed").flush().unwrap();
    pasta::bgzf::gzi::fs::write(out.with_added_extension("gzi"), &pasta::bgzf::gzi::Index::from(blocks))
        .expect("writing the .gzi");
    fs::copy(plain.with_added_extension("fai"), out.with_added_extension("fai")).expect("copying the .fai");
}

fn write_n(out: &mut impl Write, mut n: usize) {
    static PAD: [u8; 1 << 16] = [b'N'; 1 << 16];
    while n > 0 {
        let k = n.min(PAD.len());
        out.write_all(&PAD[..k]).unwrap();
        n -= k;
    }
}
