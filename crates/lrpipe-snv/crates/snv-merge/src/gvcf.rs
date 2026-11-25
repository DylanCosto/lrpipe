//! `mergeNonVariant` (`MergeVcf.py`), which runs `gvcfGenerator.mergeCalls`
//! (`preprocess/utils.py`): a contig's merged calls interleaved with its non-variant blocks.
//!
//! The original's merge loses calls and blocks and writes blocks that overlap; this one writes
//! every call and the parts of every block that no call covers (see [`merge_calls`]).

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use snv_call::{Lz4Writer, open_lz4};

use crate::postprocess::{HaploidDecider, Quality, haploid_sample};
use crate::pytext::{TextError, TextLines, parse_int, rstrip};
use crate::{Error, Haploid, RowError};

/// A line read through `lz4 -fdc` (`compressReaderWriter`): a file that cannot be read has no
/// lines, as the original ignores `lz4` failing.
fn lz4_lines(path: &Path) -> impl Iterator<Item = Result<(u64, String), Error>> {
    let reader = open_lz4(path).unwrap_or_else(|e| {
        log::warn!("cannot read {}: {e}; reading it as empty, as the original does", path.display());
        Box::new(io::empty()) as Box<dyn BufRead>
    });
    let path = path.to_owned();
    TextLines::new(reader).zip(1..).map(move |(line, n)| match line {
        Ok(line) => Ok((n, line)),
        Err(TextError::Io(source)) => Err(Error::Io { path: path.clone(), source }),
        Err(TextError::Utf8) => Err(Error::Line { path: path.clone(), line: n, reason: RowError::InvalidUtf8 }),
    })
}

/// `line.strip('\n')`.
fn strip_newlines(line: &str) -> &str {
    line.trim_matches('\n')
}

fn column<'a>(columns: &[&'a str], index: usize) -> Result<&'a str, RowError> {
    columns.get(index).copied().ok_or(RowError::MissingColumn { index })
}

fn int(value: &str) -> Result<i64, RowError> {
    parse_int(value).ok_or_else(|| RowError::InvalidInteger { value: value.to_owned() })
}

/// `ctgName` selects records when it is set and not empty.
fn selected(contig: Option<&str>, name: &str) -> bool {
    contig.is_none_or(|c| c.is_empty() || c == name)
}

/// A variant call as `readCalls(callType='variant')` yields it: the record with a `<NON_REF>`
/// allele, its first and last reference position.
type Call = (String, i64, i64);

/// Uppercase ordinary DNA alleles for readers that require it, preserving symbolic allele IDs.
fn uppercase_dna(allele: &mut String) {
    if allele.as_bytes().iter().all(|base| b"ACGTNacgtn".contains(base)) {
        allele.make_ascii_uppercase();
    }
}

/// Variant records (RefCalls are skipped), with a `<NON_REF>` allele and fields of the
/// appropriate cardinality for its genotype's ploidy.
fn variant_call(line: &str) -> Result<Option<Call>, RowError> {
    let mut columns: Vec<String> = strip_newlines(line).split('\t').map(str::to_owned).collect();
    let fields: Vec<&str> = columns.iter().map(String::as_str).collect();
    let (reference, alternate) = (column(&fields, 3)?, column(&fields, 4)?);
    let n_alt = alternate.split(',').count();
    let start = int(column(&fields, 1)?)?;
    let end = start - 1 + reference.chars().count() as i64;
    if alternate == "." || reference == alternate {
        return Ok(None);
    }
    let last = columns.len() - 1;
    let format = column(&fields, 8)?;
    let genotype = column(&fields, 9)?.split(':').next().unwrap_or(".");
    let ploidy = genotype.split(['/', '|']).count();
    let original_alleles = n_alt + 1;
    let raw_sample = if ploidy == 1 {
        haploid_sample(format, fields[last], genotype, original_alleles, Quality::Keep)
    } else {
        fields[last].to_owned()
    };
    let mut sample: Vec<String> = raw_sample.split(':').map(str::to_owned).collect();
    for (key, value) in format.split(':').zip(&mut sample) {
        if value == "." {
            continue;
        }
        match key {
            "AD" => value.push_str(",0"),
            "AF" => value.push_str(",0.0000"),
            "PL" if ploidy == 1 => value.push_str(",990"),
            "PL" if ploidy == 2 => {
                for _ in 0..=original_alleles {
                    value.push_str(",990");
                }
            }
            "PL" => *value = ".".to_owned(),
            _ => {}
        }
    }
    columns[last] = sample.join(":");
    uppercase_dna(&mut columns[3]);
    columns[4] = columns[4]
        .split(',')
        .map(|allele| {
            let mut allele = allele.to_owned();
            uppercase_dna(&mut allele);
            allele
        })
        .collect::<Vec<_>>()
        .join(",");
    columns[4].push_str(",<NON_REF>");
    Ok(Some((columns.join("\t"), start, end)))
}

/// `int(re.search(r'.*END=(.*)\tGT.*', line).group(1))`: the text between the last `END=`
/// followed by a tab-`GT` and the last tab-`GT` after it, on the line's first line.
fn block_end(line: &str) -> Result<i64, RowError> {
    let first_line = line.split('\n').next().unwrap_or("");
    let captured = first_line
        .rmatch_indices("END=")
        .find_map(|(i, _)| {
            let after = &first_line[i + 4..];
            after.rfind("\tGT").map(|j| &after[..j])
        })
        .ok_or(RowError::MissingColumn { index: 7 })?;
    int(captured)
}

/// `re.sub("END=[0-9]*\t", f"END={end}\t", line)`.
fn replace_end(line: &str, end: i64) -> String {
    let mut out = String::with_capacity(line.len());
    let mut rest = line;
    while let Some(i) = rest.find("END=") {
        let after = &rest[i + 4..];
        let digits = after.bytes().take_while(u8::is_ascii_digit).count();
        if after[digits..].starts_with('\t') {
            out.push_str(&rest[..i]);
            out.push_str(&format!("END={end}\t"));
            rest = &after[digits + 1..];
        } else {
            out.push_str(&rest[..i + 4]);
            rest = after;
        }
    }
    out.push_str(rest);
    out
}

/// A non-variant block as `readCalls(callType='non-variant')` yields it.
#[derive(Debug, Clone)]
struct Block {
    line: String,
    start: i64,
    end: i64,
}

/// `readCalls(callType='non-variant', add_header=True)`. Its header is written to the output
/// when the first record line is read, whatever its contig.
struct NonVariantCalls<'a, I> {
    lines: I,
    path: &'a Path,
    contig: Option<&'a str>,
    header: Vec<String>,
    header_written: bool,
}

impl<I: Iterator<Item = Result<(u64, String), Error>>> NonVariantCalls<'_, I> {
    fn next(&mut self, out: &mut Output) -> Result<Option<Block>, Error> {
        for line in self.lines.by_ref() {
            let (n, line) = line?;
            if line.starts_with('#') {
                if !self.header.contains(&line) {
                    self.header.push(line);
                }
                continue;
            }
            if !self.header.is_empty() && !self.header_written {
                out.write(&format!("{}\n", rstrip(&self.header.concat())))?;
                self.header_written = true;
            }
            let line_err = |reason| Error::Line { path: self.path.to_owned(), line: n, reason };
            let stripped = strip_newlines(&line);
            let columns: Vec<&str> = stripped.split('\t').collect();
            let start = int(column(&columns, 1).map_err(line_err)?).map_err(line_err)?;
            let end = block_end(&line).map_err(line_err)?;
            if !selected(self.contig, columns[0]) {
                continue;
            }
            return Ok(Some(Block { line: stripped.to_owned(), start, end }));
        }
        if !self.header.is_empty() && !self.header_written {
            out.write(&format!("{}\n", rstrip(&self.header.concat())))?;
            self.header_written = true;
        }
        Ok(None)
    }
}

/// The lz4-compressed gVCF being written.
struct Output<'a> {
    path: &'a Path,
    out: Lz4Writer<BufWriter<File>>,
    /// Lines written after the header.
    records: u64,
    decider: &'a HaploidDecider,
}

impl Output<'_> {
    fn write(&mut self, text: &str) -> Result<(), Error> {
        self.out.write_all(text.as_bytes()).map_err(|source| Error::Io { path: self.path.to_owned(), source })
    }

    fn write_line(&mut self, line: &str) -> Result<(), Error> {
        self.records += 1;
        self.write(&format!("{line}\n"))
    }

    /// The part `start..=end` of `block`: as read when that is all of it, else with the new end
    /// (`_writeLeftBlock`) and the new start and its reference base (`_writeRightBlock`).
    fn write_block(
        &mut self,
        block: &Block,
        mut start: i64,
        end: i64,
        reference: &mut ReferenceBases,
    ) -> Result<(), Error> {
        let contig = block.line.split('\t').next().unwrap_or("");
        loop {
            let boundary = self.decider.next_boundary(contig, start).filter(|&boundary| boundary <= end);
            self.write_block_segment(block, start, boundary.map_or(end, |boundary| boundary - 1), reference)?;
            let Some(boundary) = boundary else { return Ok(()) };
            start = boundary;
        }
    }

    fn write_block_segment(
        &mut self,
        block: &Block,
        start: i64,
        end: i64,
        reference: &mut ReferenceBases,
    ) -> Result<(), Error> {
        let line = if end == block.end { block.line.clone() } else { replace_end(&block.line, end) };
        let mut fields = line.split('\t');
        let contig = fields.next().unwrap_or("");
        let haploid = self.decider.decide(contig, start) != Haploid::Off;
        let anchor = fields.nth(2).unwrap_or("");
        if !haploid && start == block.start && matches!(anchor, "A" | "C" | "G" | "T") {
            return self.write_line(&line);
        }
        let mut columns: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if columns.len() < 10 {
            let reason = RowError::MissingColumn { index: columns.len() };
            return Err(Error::Record { row: block.line.clone(), reason });
        }
        if start != block.start || !matches!(anchor.as_bytes(), [b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'])
        {
            columns[3] = reference.base(&columns[0], start)?.to_string();
            columns[1] = start.to_string();
        } else {
            columns[3].make_ascii_uppercase();
        }
        if haploid {
            let genotype = if matches!(columns[9].split(':').next(), Some("0" | "0/0" | "0|0")) { "0" } else { "." };
            columns[9] = haploid_sample(&columns[8], &columns[9], genotype, 2, Quality::FromLikelihoods);
        }
        self.write_line(&columns.join("\t"))
    }
}

/// Reads a gVCF reference block's anchor from the configured FASTA.
struct ReferenceBases {
    path: PathBuf,
    reader: Option<snv_hts::Fasta>,
    contig: String,
    length: i64,
}

impl ReferenceBases {
    fn reader(&mut self) -> Result<&snv_hts::Fasta, Error> {
        if self.reader.is_none() {
            let reader = snv_hts::Fasta::open(&self.path).map_err(|_| Error::Path { path: self.path.clone() })?;
            self.reader = Some(reader);
        }
        Ok(self.reader.as_ref().expect("the reader was just opened"))
    }

    fn length(&mut self, contig: &str) -> Result<i64, Error> {
        if self.contig != contig {
            self.length = self
                .reader()?
                .contig_len(contig)
                .map_err(|source| Error::Io { path: self.path.clone(), source: io::Error::other(source) })?;
            self.contig.clear();
            self.contig.push_str(contig);
        }
        Ok(self.length)
    }

    fn base(&mut self, contig: &str, pos: i64) -> Result<char, Error> {
        let bases = self
            .reader()?
            .fetch(contig, pos - 1, pos - 1)
            .map_err(|source| Error::Io { path: self.path.clone(), source: io::Error::other(source) })?;
        let base = bases.first().ok_or_else(|| Error::Io {
            path: self.path.clone(),
            source: io::Error::new(io::ErrorKind::UnexpectedEof, format!("no base at {contig}:{pos}")),
        })?;
        Ok(char::from(lrpipe_reads::vcf_reference_base(*base)))
    }
}

/// What `mergeCalls` reads and writes.
#[derive(Debug)]
pub(crate) struct GvcfMerge<'a> {
    /// The merged VCF of the contig (`MergeVcf --output_fn`).
    pub(crate) calls: &'a Path,
    /// The sorted non-variant blocks (`--non_var_gvcf_fn`).
    pub(crate) non_variant: &'a Path,
    /// The gVCF to write (`--gvcf_fn`).
    pub(crate) output: &'a Path,
    /// For the reference base of a block cut short by a variant (`--ref_fn`).
    pub(crate) reference: &'a Path,
    /// `--ctgName`.
    pub(crate) contig: Option<&'a str>,
    /// Ploidy of the reference blocks, including PAR boundaries.
    pub(crate) decider: &'a HaploidDecider,
}

/// `mergeCalls`: writes the variant records and the non-variant blocks around them, both in
/// file order, to an lz4-compressed gVCF under the non-variant file's header. Returns the
/// number of lines written after the header.
///
/// Every call is written, at its place, and of each block the parts that no call written before
/// it covers, cut as the original cuts them (a part's new start gets the reference base there).
/// The original's merge instead (Clair3 2.0.3, `preprocess/utils.py`):
/// - writes a second copy of the current block in place of a call that ends before it (one
///   inside an earlier deletion), so the call is lost and the block repeated;
/// - writes the left part of a block that an earlier call already cut with the block's original
///   start, and the block in which its last call lies whole after it, so blocks overlap calls
///   and each other;
/// - writes the blocks after a contig's last call as read, overlapping a deletion that reaches
///   into them;
/// - drops the first block of a contig without calls, which joint genotyping takes for a
///   no-call.
pub(crate) fn merge_calls(merge: &GvcfMerge<'_>) -> Result<u64, Error> {
    let file = File::create(merge.output).map_err(|source| Error::Io { path: merge.output.to_owned(), source })?;
    let mut out = Output {
        path: merge.output,
        out: Lz4Writer::new(BufWriter::with_capacity(1 << 20, file)),
        records: 0,
        decider: merge.decider,
    };
    let mut reference =
        ReferenceBases { path: merge.reference.to_owned(), reader: None, contig: String::new(), length: 0 };

    let calls_path = merge.calls;
    let mut variants = lz4_lines(calls_path).filter_map(|line| {
        let (n, line) = match line {
            Ok(line) => line,
            Err(e) => return Some(Err(e)),
        };
        if line.starts_with('#') {
            return None;
        }
        match variant_call(&line) {
            Ok(Some(call)) if selected(merge.contig, call.0.split('\t').next().unwrap_or("")) => Some(Ok(call)),
            Ok(_) => None,
            Err(reason) => Some(Err(Error::Line { path: calls_path.to_owned(), line: n, reason })),
        }
    });
    let mut blocks = NonVariantCalls {
        lines: lz4_lines(merge.non_variant),
        path: merge.non_variant,
        contig: merge.contig,
        header: Vec::new(),
        header_written: false,
    };

    // The next call, and the first position after every call written so far.
    let mut call = variants.next().transpose()?;
    let mut covered_to = i64::MIN;
    // The first block is read before any call is written, so the header comes first.
    while let Some(block) = blocks.next(&mut out)? {
        let contig = block.line.split('\t').next().unwrap_or("");
        let end = block.end.min(reference.length(contig)?);
        let mut start = block.start.max(covered_to).max(1);
        while start <= end {
            match call.take() {
                // A call that starts before the rest of the block ends (or before the block, in a
                // gap between blocks or inside an earlier deletion): the block up to it, then it.
                Some((line, call_start, call_end)) if call_start <= end => {
                    if call_start > start {
                        out.write_block(&block, start, call_start - 1, &mut reference)?;
                    }
                    out.write_line(&line)?;
                    covered_to = covered_to.max(call_end + 1);
                    start = start.max(covered_to);
                    call = variants.next().transpose()?;
                }
                rest => {
                    call = rest;
                    out.write_block(&block, start, end, &mut reference)?;
                    break;
                }
            }
        }
    }
    while let Some((line, ..)) = call {
        out.write_line(&line)?;
        call = variants.next().transpose()?;
    }
    let records = out.records;
    out.out.finish().map_err(|source| Error::Io { path: merge.output.to_owned(), source })?;
    Ok(records)
}

#[cfg(test)]
mod tests {
    use std::io::Read as _;

    use super::*;

    #[test]
    fn variant_calls_gain_non_ref() {
        let line = "chr1\t10\t.\tAT\tA\t20.5\tPASS\tF\tGT:GQ:DP:AD:AF:PL\t0/1:20:9:4,5:0.5556:20,0,30\n";
        assert_eq!(
            variant_call(line).unwrap(),
            Some((
                "chr1\t10\t.\tAT\tA,<NON_REF>\t20.5\tPASS\tF\tGT:GQ:DP:AD:AF:PL\t0/1:20:9:4,5,0:0.5556,0.0000:20,0,30,990,990,990"
                    .to_owned(),
                10,
                11
            ))
        );
        let multi = "chr1\t10\t.\tA\tC,G\t20\tPASS\tP\tGT:AD:AF:PL\t1/2:0,3,4:0.4,0.5:9,9,9,9,0,9\n";
        assert!(variant_call(multi).unwrap().unwrap().0.ends_with(":9,9,9,9,0,9,990,990,990,990"));
        assert_eq!(variant_call("chr1\t10\t.\tA\t.\t3\tRefCall\tP\tGT:AD:AF:PL\t0/0:3:1.0:990\n").unwrap(), None);
        assert_eq!(variant_call("chr1\t10\n"), Err(RowError::MissingColumn { index: 3 }));
    }

    #[test]
    fn haploid_softmasked_alleles_gain_one_likelihood() {
        let call = "chrX\t10\t.\tAaa\tA\t20.5\tPASS\tF\tGT:GQ:DP:AD:AF:PL\t1:20:9:4,5:0.5556:20,0";
        let (line, _, _) = variant_call(call).unwrap().unwrap();
        assert!(line.contains("\tAAA\tA,<NON_REF>\t"));
        assert!(line.ends_with("\t1:20:9:4,5,0:0.5556,0.0000:20,0,990"));
    }

    /// Each case of `tests/data/gvcf_merge` holds `calls.vcf` and `non_var.gvcf`, and
    /// `expected.gvcf`, the decompressed output of the original's `gvcfGenerator.mergeCalls` on
    /// them (contig `chrT`; `haploid_genotypes` with `haploid_calling=True`) but where its bugs
    /// apply ([`merge_calls`]): the call inside a deletion is kept (`variant_inside_deletion`),
    /// blocks after a call start after it (`adjacent_variants`, `deletion_across_blocks`,
    /// `haploid_genotypes`, `multi_allelic_and_insertion`, `other_contigs`, `snp_inside_block`,
    /// `variant_inside_deletion`), a contig without calls keeps its first block
    /// (`no_variants`, `only_ref_calls`), and a contig without blocks keeps its header (`no_blocks`).
    #[test]
    fn merge_calls_matches_the_original_but_for_its_bugs() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gvcf_merge");
        let out_dir = std::env::temp_dir().join(format!("snv-merge-gvcf-{}", std::process::id()));
        std::fs::create_dir_all(&out_dir).unwrap();
        let mut cases = 0;
        for entry in std::fs::read_dir(&data).unwrap() {
            let case = entry.unwrap().path();
            if !case.is_dir() {
                continue;
            }
            let name = case.file_name().unwrap().to_string_lossy().into_owned();
            let output = out_dir.join(format!("{name}.gvcf"));
            merge_calls(&GvcfMerge {
                calls: &case.join("calls.vcf"),
                non_variant: &case.join("non_var.gvcf"),
                output: &output,
                reference: &data.join("ref.fa"),
                contig: Some("chrT"),
                decider: &HaploidDecider::new(
                    if name == "haploid_genotypes" { Haploid::Precise } else { Haploid::Off },
                    crate::Gender::default(),
                    crate::postprocess::ParRegions::default(),
                ),
            })
            .unwrap();
            let mut got = String::new();
            open_lz4(&output).unwrap().read_to_string(&mut got).unwrap();
            assert_eq!(got, std::fs::read_to_string(case.join("expected.gvcf")).unwrap(), "case {name}");
            cases += 1;
        }
        assert_eq!(cases, 13);
        std::fs::remove_dir_all(&out_dir).unwrap();
    }

    #[test]
    fn block_ends_and_left_cuts() {
        let line = "chr1\t5\t.\tA\t<NON_REF>\t0\t.\tEND=120\tGT:GQ:MIN_DP:PL\t0/0:6:3:0,9,90\n";
        assert_eq!(block_end(line), Ok(120));
        assert_eq!(
            replace_end(strip_newlines(line), 7),
            "chr1\t5\t.\tA\t<NON_REF>\t0\t.\tEND=7\tGT:GQ:MIN_DP:PL\t0/0:6:3:0,9,90"
        );
        assert_eq!(replace_end("END=x\tEND=12\t", 3), "END=x\tEND=3\t");
        assert!(block_end("chr1\t5\t.\n").is_err());
    }
}
