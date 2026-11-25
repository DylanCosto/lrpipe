//! What `lrpipe run --reads` names, told apart in one place, [`classify`]: one sorted, indexed
//! BAM or CRAM, which the stages read as it is; several, which the `merge-reads` stage merges
//! into one BAM first ([`crate::merge_reads`]); or unaligned reads, which the align stage maps
//! with bonefish first: unaligned BAMs (as dorado and PacBio write them), or FASTQ or FASTQ.gz,
//! given once or more.
//!
//! A file's format is told from its contents, not its name. A BAM without reference sequences
//! (`@SQ` lines) is unaligned, and so is one with them whose reads are all unmapped (as the first
//! [`MAPPED_SAMPLE`] show). Any other is aligned, and is never sent through bonefish, though
//! bonefish would take it. bonefish reads no CRAM, so an unaligned CRAM is refused, and aligned
//! and unaligned files together are refused too.
//!
//! An aligned file must also be named as what it is (`.cram` for CRAM, else BAM: every tool
//! tells them apart by name), and is checked by [`crate::checks::aligned`]: its sort order and
//! its index, whose read counts per contig the reference check reads.
//!
//! The read groups must name one sample (`SM`). For a single aligned file, `lrpipe run` accepts
//! other labels for that donor when each is given explicitly with `--sample-alias`. Several
//! files' read groups must also agree on the fields of any shared `@RG` ID ([`combine_headers`]),
//! and several aligned files must have the same `@SQ` lines (names, lengths and order).
//! Their merged header is the first file's with every file's read groups and `@PG` chains. A read
//! in two aligned files would count twice; [`classify`] warns of that from a sample of each
//! file's read names. A file given twice is refused.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use lrpipe_reads::{Decoding, Fields, IndexedReads, MdNm, Reference, Tags};
use pasta::sam::alignment::record::data::field::Tag;
use pasta::{bam, cram, sam};

use crate::plan::{self, BamIndex};

/// What a file of reads is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// A BAM or CRAM of aligned reads.
    Aligned,
    /// A BAM without `@SQ` lines, or of unmapped reads alone: unaligned reads.
    UnalignedBam,
    /// FASTQ, plain or gzipped.
    Fastq,
}

/// The reads `--reads` names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Input {
    /// One sorted, indexed BAM or CRAM, which the stages read.
    Aligned(PathBuf),
    /// Several sorted, indexed BAMs and CRAMs, in the order given, which the `merge-reads` stage
    /// merges into one BAM with this index.
    SeveralAligned(Vec<PathBuf>, BamIndex),
    /// Unaligned reads, each file with its format ([`Format::UnalignedBam`] or [`Format::Fastq`]),
    /// in the order given, which the align stage maps into one BAM.
    Unaligned(Vec<(PathBuf, Format)>),
}

/// What [`classify`] found.
#[derive(Clone, Debug)]
pub struct Classified {
    /// The reads.
    pub input: Input,
    /// Their header, whose read groups say how they were made: the one aligned file's, several
    /// files' merged ([`combine_headers`]), or for unaligned reads the unaligned BAMs' merged
    /// (FASTQ has none: empty).
    pub header: sam::Header,
    /// The aligned files as [`crate::checks::aligned`] found them, in the order given; none for
    /// unaligned reads, whose aligned BAM is checked once the align stage has written it.
    pub aligned: Vec<crate::checks::Aligned>,
    /// What the run warns of.
    pub warnings: Vec<String>,
}

/// Why the reads cannot be taken.
#[derive(Debug, thiserror::Error)]
pub enum InputError {
    /// A file could not be read.
    #[error("reading {}", path.display())]
    Read {
        /// The file.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// A file is none of the formats taken.
    #[error("{}: {what}; --reads takes a BAM or CRAM, aligned or (BAM only) unaligned, or FASTQ (.gz)", path.display())]
    Unknown {
        /// The file.
        path: PathBuf,
        /// What it is, as far as it can be told.
        what: &'static str,
    },
    /// Aligned and unaligned reads together.
    #[error(
        "{} is aligned and {} is not: give either aligned BAMs or CRAMs (several are merged), or unaligned reads \
         (unaligned BAM, FASTQ), which lrpipe aligns",
        aligned.display(),
        unaligned.display()
    )]
    Mixed {
        /// An aligned file.
        aligned: PathBuf,
        /// An unaligned one.
        unaligned: PathBuf,
    },
    /// A file given twice.
    #[error("{} is given twice in --reads", .0.display())]
    Twice(PathBuf),
    /// A file that is not there.
    #[error("--reads {} does not exist", .0.display())]
    Missing(PathBuf),
    /// An aligned file named as the other format.
    #[error(
        "--reads {} is {is}, but named as {named}, and lrpipe's tools take the format from the name: rename it to \
         end in .{ext} (and its index with it)",
        path.display()
    )]
    Misnamed {
        /// The file.
        path: PathBuf,
        /// What it is.
        is: &'static str,
        /// What its name says.
        named: &'static str,
        /// The extension it needs.
        ext: &'static str,
    },
    /// An aligned file fails [`crate::checks::aligned`] (its sort order or its index).
    #[error("{0:#}")]
    Aligned(anyhow::Error),
    /// Unaligned BAM and FASTQ together.
    #[error(
        "{} is an unaligned BAM and {} FASTQ: give unaligned reads of one kind (the BAMs keep the read groups that \
         name the basecall model; FASTQ has none)",
        bam.display(),
        fastq.display()
    )]
    BamAndFastq {
        /// An unaligned BAM.
        bam: PathBuf,
        /// A FASTQ.
        fastq: PathBuf,
    },
    /// A file's header could not be read.
    #[error("reading the header of {}", path.display())]
    Header {
        /// The reads.
        path: PathBuf,
        /// The cause.
        #[source]
        source: io::Error,
    },
    /// Two aligned files' reference sequences differ.
    #[error(
        "{} and {} are aligned to different references ({difference}): every --reads must be aligned to the \
         same one",
        first.display(),
        other.display()
    )]
    ReferenceSequences {
        /// The first file.
        first: PathBuf,
        /// The one that differs from it.
        other: PathBuf,
        /// The first difference.
        difference: String,
    },
    /// Two files have a read group of the same ID with different fields: the reads' `RG` tags
    /// would not say which it is.
    #[error(
        "{} and {} both have read group {id} (@RG ID:{id}), with different fields: give the read groups distinct \
         IDs (`samtools addreplacerg`), or the same fields (`samtools reheader`)",
        first.display(),
        other.display()
    )]
    ReadGroup {
        /// The ID.
        id: String,
        /// The first file with it.
        first: PathBuf,
        /// The one whose fields differ.
        other: PathBuf,
    },
    /// Nearly all the reads sampled from one file are in another: the same reads given twice (a
    /// BAM and its CRAM, a file and a copy, a file and the merge of it).
    #[error(
        "{shared} of the {sampled} reads sampled from {} are in {} too: they hold the same reads (a file and its \
         copy, a BAM and its CRAM, or a merge of them), which would count twice (in depth, allele and SV support); \
         give each read once",
        from.display(),
        other.display()
    )]
    SameReads {
        /// Reads sampled from `from` that `other`'s sample has.
        shared: usize,
        /// Reads sampled from `from`: the smaller sample of the two, which the share is of.
        sampled: usize,
        /// The file of the smaller sample.
        from: PathBuf,
        /// The other file.
        other: PathBuf,
    },
    /// The files' read groups name several samples.
    #[error(
        "the --reads are of several samples (@RG SM: {}): `lrpipe run` calls one; give one sample's files",
        .0.iter().map(|(sample, path)| format!("{sample} in {}", path.display())).collect::<Vec<_>>().join(", ")
    )]
    Samples(Vec<(String, PathBuf)>),
    /// One aligned input's read groups name several samples without a complete alias list.
    #[error(
        "{}: its read groups name several samples (@RG SM: {names}): `lrpipe run` calls one; if these labels are \
         all from the same donor, give --sample-alias for each name other than --sample",
        path.display()
    )]
    SampleAliases {
        /// The aligned input.
        path: PathBuf,
        /// The sample names in its header.
        names: String,
    },
    /// Aliases can only be declared for one aligned input.
    #[error("--sample-alias requires one aligned BAM or CRAM as --reads")]
    SampleAliasInput,
}

/// Records of an aligned-looking BAM (one with `@SQ` lines) [`format`] reads for a mapped one; if
/// the first are all unmapped, it is unaligned. In a coordinate-sorted BAM unmapped reads come
/// last, so there the first record decides.
pub const MAPPED_SAMPLE: usize = 1000;

/// The format of the reads at `path`, by its first bytes, and for a BAM or CRAM its header (and
/// for a BAM with `@SQ` lines its first records).
///
/// # Errors
///
/// If the file cannot be read, or is none of the formats taken (SAM, FASTA, an unaligned CRAM,
/// ...).
pub fn format(path: &Path) -> Result<Format, InputError> {
    let read = |source| InputError::Read { path: path.to_path_buf(), source };
    let unknown = |what| InputError::Unknown { path: path.to_path_buf(), what };
    let mut magic = Vec::with_capacity(4);
    File::open(path).and_then(|f| f.take(4).read_to_end(&mut magic)).map_err(read)?;
    if magic == b"CRAM" {
        let header = File::open(path).and_then(|f| pasta::cram::io::Reader::new(f).read_header()).map_err(read)?;
        if header.reference_sequences().is_empty() {
            return Err(unknown(
                "an unaligned CRAM, which bonefish cannot read (convert it to BAM: samtools view -b -o reads.bam)",
            ));
        }
        if !plan::is_cram(path) {
            return Err(InputError::Misnamed { path: path.to_path_buf(), is: "CRAM", named: "BAM", ext: "cram" });
        }
        return Ok(Format::Aligned);
    }
    // The first lines (of what a gzip or BGZF file holds).
    let file = File::open(path).map_err(read)?;
    let mut text: Box<dyn BufRead> = if magic.starts_with(&[0x1f, 0x8b]) {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut head = Vec::new();
    (&mut text).take(1 << 16).read_to_end(&mut head).map_err(read)?;
    if head.starts_with(b"BAM\x01") {
        let format = bam_format(path).map_err(read)?;
        if format == Format::Aligned && plan::is_cram(path) {
            return Err(InputError::Misnamed { path: path.to_path_buf(), is: "BAM", named: "CRAM", ext: "bam" });
        }
        return Ok(format);
    }
    let mut lines = head.split(|b| *b == b'\n');
    let first = lines.next().unwrap_or_default();
    let third = lines.nth(1);
    match first.first() {
        // A FASTQ record: a name, its bases, and a `+` line (unless the bases wrap).
        Some(b'@') if third.is_some_and(|l| l.starts_with(b"+")) => Ok(Format::Fastq),
        // A SAM header line.
        Some(b'@') if [&b"@HD\t"[..], b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"].iter().any(|t| first.starts_with(t)) => {
            Err(unknown(
                "SAM text, which lrpipe does not read: convert, sort and index it, as with `samtools sort -o \
                 reads.bam READS && samtools index reads.bam`",
            ))
        }
        Some(b'@') => Ok(Format::Fastq),
        Some(b'>') => Err(unknown("FASTA, whose reads have no base qualities, which the callers need")),
        None if head.is_empty() => Err(unknown("empty, with no reads")),
        _ => Err(unknown("not a format lrpipe reads")),
    }
}

/// An aligned or unaligned BAM, by its `@SQ` lines and first records.
fn bam_format(path: &Path) -> io::Result<Format> {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;
    if header.reference_sequences().is_empty() {
        return Ok(Format::UnalignedBam);
    }
    let mut records = 0;
    for record in reader.records().take(MAPPED_SAMPLE) {
        if !record?.flags().is_unmapped() {
            return Ok(Format::Aligned);
        }
        records += 1;
    }
    // A BAM with no reads at all is taken as it says: aligned (to nothing).
    Ok(if records == 0 { Format::Aligned } else { Format::UnalignedBam })
}

/// What `paths` (`--reads`, at least one, in the order given) are: one aligned file, several, or
/// unaligned files of one kind; after checking that every file is there and none is given twice,
/// aligned files pass [`crate::checks::aligned`], and several files are one sample's reads ([`combine_headers`]). Of several
/// aligned files it also samples each one's read names, and warns of reads in two of them. CRAM
/// is decoded with `reference`.
///
/// # Errors
///
/// If a file is missing, cannot be read or is none of the formats taken, a file is given twice,
/// two hold the same reads ([`InputError::SameReads`]), aligned files are given with unaligned ones, unaligned BAMs come with FASTQ, an aligned file is
/// named as the other format or fails [`crate::checks::aligned`], or the files' headers cannot be
/// merged.
pub fn classify(paths: &[PathBuf], reference: &Path) -> Result<Classified, InputError> {
    classify_with_sample_aliases(paths, reference, None)
}

/// [`classify`] with an output sample name and explicit alternative labels for one aligned file.
pub(crate) fn classify_with_sample_aliases(
    paths: &[PathBuf],
    reference: &Path,
    sample_aliases: Option<(&str, &[String])>,
) -> Result<Classified, InputError> {
    let mut seen: Vec<PathBuf> = Vec::new();
    for path in paths {
        if !path.exists() {
            return Err(InputError::Missing(path.clone()));
        }
        let real = std::fs::canonicalize(path).map_err(|source| InputError::Read { path: path.clone(), source })?;
        if seen.contains(&real) {
            return Err(InputError::Twice(path.clone()));
        }
        seen.push(real);
    }
    let formats: Vec<(PathBuf, Format)> =
        paths.iter().map(|p| Ok((p.clone(), format(p)?))).collect::<Result<_, _>>()?;
    let aligned: Vec<&PathBuf> = formats.iter().filter(|(_, f)| *f == Format::Aligned).map(|(p, _)| p).collect();
    let unaligned = formats.iter().find(|(_, f)| *f != Format::Aligned);
    if sample_aliases.is_some() && (aligned.len() != 1 || unaligned.is_some()) {
        return Err(InputError::SampleAliasInput);
    }
    let headers = |files: &[&PathBuf]| -> Result<Vec<(PathBuf, sam::Header)>, InputError> {
        files
            .iter()
            .map(|p| {
                let header = read_header(p).map_err(|source| InputError::Header { path: (*p).clone(), source })?;
                Ok(((*p).clone(), header))
            })
            .collect()
    };
    match (&aligned[..], unaligned) {
        ([], _) => {
            let bam = formats.iter().find(|(_, f)| *f == Format::UnalignedBam);
            let fastq = formats.iter().find(|(_, f)| *f == Format::Fastq);
            if let (Some((bam, _)), Some((fastq, _))) = (bam, fastq) {
                return Err(InputError::BamAndFastq { bam: bam.clone(), fastq: fastq.clone() });
            }
            // bonefish keeps every BAM's read groups and each read's RG tag: they must be one
            // sample's and not clash. Unaligned BAMs' @SQ lines, if any, are not the reference's.
            let bams: Vec<&PathBuf> =
                formats.iter().filter(|(_, f)| *f == Format::UnalignedBam).map(|(p, _)| p).collect();
            let header = combine(&headers(&bams)?, Sequences::Unchecked)?;
            Ok(Classified { input: Input::Unaligned(formats), header, aligned: Vec::new(), warnings: Vec::new() })
        }
        ([first, ..], Some((other, _))) => {
            Err(InputError::Mixed { aligned: (*first).clone(), unaligned: other.clone() })
        }
        (files, None) => {
            let headers = headers(files)?;
            // Each one's sort order and index, before anything reads the index.
            let aligned: Vec<crate::checks::Aligned> = headers
                .iter()
                .map(|(path, header)| crate::checks::aligned(path, header.clone()))
                .collect::<anyhow::Result<_>>()
                .map_err(InputError::Aligned)?;
            let mut warnings: Vec<String> = aligned.iter().flat_map(|a| a.warnings.iter().cloned()).collect();
            if let [(only, header)] = &headers[..] {
                let samples = read_group_samples(&headers);
                if samples.len() > 1 {
                    let names = samples
                        .iter()
                        .map(|(sample, _)| String::from_utf8_lossy(sample))
                        .collect::<Vec<_>>()
                        .join(", ");
                    let Some((sample, _)) = sample_aliases.filter(|(sample, aliases)| {
                        samples.iter().all(|(name, _)| {
                            *name == sample.as_bytes() || aliases.iter().any(|alias| *name == alias.as_bytes())
                        })
                    }) else {
                        return Err(InputError::SampleAliases { path: only.clone(), names });
                    };
                    warnings.push(format!(
                        "{}: treating read-group sample names (@RG SM: {names}) as aliases for --sample {sample}",
                        only.display()
                    ));
                }
                return Ok(Classified {
                    input: Input::Aligned(only.clone()),
                    header: header.clone(),
                    aligned,
                    warnings,
                });
            }
            let header = combine_headers(&headers)?;
            warnings.extend(shared_reads(paths, &header, reference)?);
            let input = Input::SeveralAligned(paths.to_vec(), bam_index(&header));
            Ok(Classified { input, header, aligned, warnings })
        }
    }
}

/// The header of the BAM or CRAM at `path` (CRAM by its name).
///
/// # Errors
///
/// If the file cannot be read, or its header is not BAM's or CRAM's.
pub fn read_header(path: &Path) -> io::Result<sam::Header> {
    let file = File::open(path)?;
    if plan::is_cram(path) {
        cram::io::Reader::new(file).read_header()
    } else {
        bam::io::Reader::new(file).read_header()
    }
}

/// The header of several inputs' reads merged (each file's path and header, in the order
/// given): the first one's, sorted by coordinate, with every input's read groups, in the order
/// they come; every input's `@PG` lines, as `samtools merge` keeps them (a line another file
/// has with other fields renamed, and the lines naming it as their previous program changed to
/// match, so each file's chain stays whole); and their comments, each once.
///
/// # Errors
///
/// If two inputs' `@SQ` lines differ in names, lengths or order, two have a read group of
/// the same ID with different fields, or their read groups name several samples (`SM`).
pub fn combine_headers(headers: &[(PathBuf, sam::Header)]) -> Result<sam::Header, InputError> {
    combine(headers, Sequences::Same)
}

/// Whether [`combine`] requires the same `@SQ` lines.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Sequences {
    /// Aligned reads: the same reference.
    Same,
    /// Unaligned reads: any.
    Unchecked,
}

/// [`combine_headers`], checking the `@SQ` lines as `sequences` says.
fn combine(headers: &[(PathBuf, sam::Header)], sequences: Sequences) -> Result<sam::Header, InputError> {
    use sam::header::record::value::map::{self, header::tag::SORT_ORDER};

    let Some((first_path, first)) = headers.first() else { return Ok(sam::Header::default()) };
    let mut merged = first.clone();
    let hd = merged.header_mut().get_or_insert_with(sam::header::record::value::Map::<map::Header>::default);
    hd.other_fields_mut().insert(SORT_ORDER, "coordinate".into());
    for (input, (path, header)) in headers.iter().enumerate().skip(1) {
        if let Some(difference) = sequence_difference(first, header).filter(|_| sequences == Sequences::Same) {
            return Err(InputError::ReferenceSequences { first: first_path.clone(), other: path.clone(), difference });
        }
        for (id, group) in header.read_groups() {
            match merged.read_groups().get(id) {
                Some(known) if known != group => {
                    let first = headers
                        .iter()
                        .find(|(_, h)| h.read_groups().get(id).is_some())
                        .map_or_else(|| first_path.clone(), |(p, _)| p.clone());
                    return Err(InputError::ReadGroup { id: id.to_string(), first, other: path.clone() });
                }
                Some(_) => {}
                None => {
                    merged.read_groups_mut().insert(id.clone(), group.clone());
                }
            }
        }
        add_programs(&mut merged, header, input + 1);
        for comment in header.comments() {
            if !merged.comments().contains(comment) {
                merged.add_comment(comment.clone());
            }
        }
    }
    let samples = read_group_samples(headers);
    if samples.len() > 1 {
        return Err(InputError::Samples(
            samples
                .into_iter()
                .map(|(name, path)| (String::from_utf8_lossy(name).into_owned(), path.to_path_buf()))
                .collect(),
        ));
    }
    Ok(merged)
}

/// Each sample the read groups name, with the first file naming it.
fn read_group_samples(headers: &[(PathBuf, sam::Header)]) -> Vec<(&[u8], &Path)> {
    let mut samples = Vec::new();
    for (path, header) in headers {
        for group in header.read_groups().values() {
            let Some(sample) = group.other_fields().get(&sam::header::record::value::map::read_group::tag::SAMPLE)
            else {
                continue;
            };
            let sample: &[u8] = sample.as_ref();
            if !samples.iter().any(|(name, _)| *name == sample) {
                samples.push((sample, path.as_path()));
            }
        }
    }
    samples
}

/// Adds the `@PG` lines of `header`, the `number`th input (from 1), to `merged`'s: one `merged` has
/// with the same fields is there already; one whose ID it has with other fields is renamed
/// `ID-number` (or `ID-number.2`, ...), and its file's lines that follow it (`PP`) follow the new
/// name.
fn add_programs(merged: &mut sam::Header, header: &sam::Header, number: usize) {
    use sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;
    let programs = header.programs().as_ref();
    let mut renamed: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    // The line as it would be added: following renamed lines.
    let remapped = |program: &sam::header::record::value::Map<sam::header::record::value::map::Program>,
                    renamed: &std::collections::HashMap<String, String>| {
        let mut program = program.clone();
        let previous = program.other_fields().get(&PREVIOUS_PROGRAM_ID).and_then(|p| renamed.get(&p.to_string()));
        if let Some(previous) = previous.cloned() {
            program.other_fields_mut().insert(PREVIOUS_PROGRAM_ID, previous.into());
        }
        program
    };
    // A line is merged's own when its fields, following renamed lines, are; renaming one can make
    // the lines after it differ, so until none is renamed.
    loop {
        let mut more = false;
        for (id, program) in programs {
            let key = id.to_string();
            if renamed.contains_key(&key) {
                continue;
            }
            let Some(known) = merged.programs().as_ref().get(id) else { continue };
            if *known == remapped(program, &renamed) {
                continue;
            }
            let taken = |candidate: &str| {
                merged.programs().as_ref().keys().any(|k| k == candidate.as_bytes())
                    || programs.keys().any(|k| k == candidate.as_bytes())
                    || renamed.values().any(|v| v == candidate)
            };
            let mut name = format!("{key}-{number}");
            let mut n = 2;
            while taken(&name) {
                name = format!("{key}-{number}.{n}");
                n += 1;
            }
            renamed.insert(key, name);
            more = true;
        }
        if !more {
            break;
        }
    }
    for (id, program) in programs {
        let new_id = renamed.get(&id.to_string());
        if new_id.is_none() && merged.programs().as_ref().contains_key(id) {
            continue;
        }
        let program = remapped(program, &renamed);
        let id = new_id.map_or_else(|| id.clone(), |n| n.as_str().into());
        merged.programs_mut().as_mut().insert(id, program);
    }
}

/// How `other`'s `@SQ` lines first differ from `first`'s, if they do (in name, length or
/// order, or in MD5 checksums where both specify them).
fn sequence_difference(first: &sam::Header, other: &sam::Header) -> Option<String> {
    let lines = |h: &sam::Header| -> Vec<(String, usize)> {
        h.reference_sequences().iter().map(|(name, map)| (name.to_string(), map.length().get())).collect()
    };
    let (a, b) = (lines(first), lines(other));
    let differs = a.iter().zip(&b).position(|(x, y)| x != y);
    match differs {
        Some(i) => Some(format!(
            "@SQ line {} is {} of {} bp in the first and {} of {} bp in the other",
            i + 1,
            a[i].0,
            a[i].1,
            b[i].0,
            b[i].1
        )),
        None if a.len() != b.len() => Some(format!("{} @SQ lines in the first and {} in the other", a.len(), b.len())),
        None => first.reference_sequences().iter().zip(other.reference_sequences()).find_map(|((name, a), (_, b))| {
            let tag = &sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
            let a: &[u8] = a.other_fields().get(tag)?.as_ref();
            let b: &[u8] = b.other_fields().get(tag)?.as_ref();
            (!a.eq_ignore_ascii_case(b)).then(|| format!("{name} has conflicting @SQ M5 checksums"))
        }),
    }
}

/// The merged BAM's index for `header`'s reference sequences: a BAI unless one is too long.
fn bam_index(header: &sam::Header) -> BamIndex {
    let longest = header.reference_sequences().values().map(|r| r.length().get()).max().unwrap_or(1);
    crate::merge_reads::bam_index_for(longest)
}

/// Contigs whose reads [`shared_reads`] samples: the first long ones with reads in two files.
const SAMPLED_CONTIGS: usize = 4;

/// Primary reads whose names [`shared_reads`] takes from each sampled contig of each input.
const SAMPLED_NAMES: usize = 500;

/// Ends a query once enough names are taken.
#[derive(Debug)]
struct Enough;

impl std::fmt::Display for Enough {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("enough names")
    }
}

impl std::error::Error for Enough {}

/// The share of the reads sampled from a file (the smaller sample of two) that, found in another,
/// makes them the same reads given twice, an error; fewer shared reads are a warning.
pub const SAME_READS_SHARE: f64 = 0.9;

/// Warnings of reads in two of `paths`: [`SAMPLED_NAMES`] primary reads' names from the middle on
/// of each of the first [`SAMPLED_CONTIGS`] contigs of 1 Mb or more with reads in two files or
/// more (all contigs, if none is that long), the same contigs and places in each file, looked for
/// in the others' samples. Two files sharing more than [`SAME_READS_SHARE`] of the smaller
/// sample are refused. A few BGZF blocks or CRAM containers of each file are read.
fn shared_reads(paths: &[PathBuf], header: &sam::Header, reference: &Path) -> Result<Vec<String>, InputError> {
    let reference = if paths.iter().any(|p| plan::is_cram(p)) {
        Some(Reference::open(reference).map_err(|e| InputError::Read { path: reference.into(), source: other(e) })?)
    } else {
        None
    };
    let contigs: Vec<(String, usize)> =
        header.reference_sequences().iter().map(|(name, map)| (name.to_string(), map.length().get())).collect();
    let long: Vec<&(String, usize)> = contigs.iter().filter(|(_, len)| *len >= lrpipe_reads::LONG_CONTIG).collect();
    let candidates = if long.is_empty() { contigs.iter().collect() } else { long };
    let mut files = Vec::with_capacity(paths.len());
    for path in paths {
        let mut file = IndexedReads::open(path, reference.as_ref()).map_err(|e| read_error(path, other(e)))?;
        file.set_decoding(Decoding {
            md_nm: MdNm::AsStored,
            fields: Fields::NAME
                | Fields::FLAGS
                | Fields::REFERENCE_SEQUENCE_ID
                | Fields::ALIGNMENT_START
                | Fields::CIGAR,
            tags: Tags::Only(Vec::new()),
        });
        files.push(file);
    }
    // The same contigs in every file: the first with reads in two files or more, where a read
    // can be in two.
    let mut sampled: Vec<(&str, i64)> = Vec::new();
    for (contig, len) in candidates {
        if sampled.len() == SAMPLED_CONTIGS {
            break;
        }
        let len = *len as i64;
        let mut with_reads = 0;
        for (file, path) in files.iter().zip(paths) {
            with_reads += usize::from(file.indexed_bytes(contig, 0, len).map_err(|e| read_error(path, e))? > 0);
        }
        if with_reads >= 2 {
            sampled.push((contig, len));
        }
    }
    let mut samples: Vec<HashSet<Vec<u8>>> = Vec::with_capacity(paths.len());
    for (file, path) in files.iter_mut().zip(paths) {
        let mut names = HashSet::new();
        for &(contig, len) in &sampled {
            let mut taken = 0;
            for (start, end) in [(len / 2, len), (0, len / 2)] {
                let read = file.query(contig, start, end, |record| {
                    let flags = record.flags();
                    if taken == SAMPLED_NAMES {
                        return Err(io::Error::other(Enough));
                    }
                    if !(flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary())
                        && let Some(name) = record.name()
                        && names.insert(name.to_vec())
                    {
                        taken += 1;
                    }
                    Ok(())
                });
                match read {
                    Err(e) if e.get_ref().is_some_and(|e| e.is::<Enough>()) => {}
                    other => other.map_err(|e| read_error(path, e))?,
                }
            }
        }
        samples.push(names);
    }
    let mut warnings = Vec::new();
    for (j, later) in samples.iter().enumerate().skip(1) {
        for (i, earlier) in samples[..j].iter().enumerate() {
            let shared = later.intersection(earlier).count();
            let sampled = [(&paths[j], later.len()), (&paths[i], earlier.len())];
            warnings.extend(shared_verdict(shared, sampled)?);
        }
    }
    Ok(warnings)
}

/// The share of a sample shared with another file's, from which it is said that one of the files
/// may be a merge holding the other: a file's reads are half or so of a merge's of it and as many
/// more, at the same places, so about half its sample is in the merge's.
const MERGE_SHARE: f64 = 0.35;

/// What `shared` reads of two files' samples (each file and its sample's size, the later given
/// first) make: nothing, a warning, or [`InputError::SameReads`]. The share is of the smaller
/// sample, which the messages name.
fn shared_verdict(shared: usize, sampled: [(&PathBuf, usize); 2]) -> Result<Option<String>, InputError> {
    if shared == 0 {
        return Ok(None);
    }
    let [(from, n), (other, _)] = if sampled[1].1 < sampled[0].1 { [sampled[1], sampled[0]] } else { sampled };
    let share = shared as f64 / n.max(1) as f64;
    if share > SAME_READS_SHARE {
        return Err(InputError::SameReads { shared, sampled: n, from: from.clone(), other: other.clone() });
    }
    Ok(Some(format!(
        "{shared} of the {n} reads sampled from {} are in {} too{}: a read given twice counts twice (in depth, \
         allele and SV support); give each read once",
        from.display(),
        other.display(),
        if share >= MERGE_SHARE {
            " (as many as a file and a merge that holds it share: is one of them a merge of the other and more?)"
        } else {
            ""
        }
    )))
}

/// Reading `path` failed with `source`.
fn read_error(path: &Path, source: io::Error) -> InputError {
    InputError::Read { path: path.to_path_buf(), source }
}

/// An [`lrpipe_reads::OpenError`] as the I/O error it is, keeping its message.
fn other(e: lrpipe_reads::OpenError) -> io::Error {
    io::Error::other(e)
}

/// Reads of an unaligned file [`has_modifications`] and [`fastq_read_groups`] look at.
pub const MODIFICATION_SAMPLE: usize = 100;

/// The read groups dorado's FASTQ names in its reads' comments (`RG:Z:RUNID_MODEL`, with a barcode
/// after when demultiplexed), each once with the basecall model it names
/// (`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`), of the first [`MODIFICATION_SAMPLE`] reads; none for
/// FASTQ without them, or whose tags name no versioned model.
///
/// # Errors
///
/// If the file cannot be read.
pub fn fastq_read_groups(path: &Path) -> Result<Vec<(String, String)>, InputError> {
    let mut groups: Vec<(String, String)> = Vec::new();
    fastq_names(path, |name| {
        for tag in name.split(|b| *b == b'\t' || *b == b' ').skip(1) {
            let Some(id) = tag.strip_prefix(b"RG:Z:") else { continue };
            let id = String::from_utf8_lossy(id).into_owned();
            if groups.iter().any(|(g, _)| *g == id) {
                continue;
            }
            // The run's ID, then the model up to the end of its version.
            let Some((_, rest)) = id.split_once('_') else { continue };
            let Some(at) = rest.find("@v") else { continue };
            let end =
                rest[at + 2..].find(|c: char| !(c.is_ascii_digit() || c == '.')).map_or(rest.len(), |e| at + 2 + e);
            let model = rest[..end].trim_end_matches('.').to_owned();
            groups.push((id, model));
        }
        false
    })?;
    Ok(groups)
}

/// Whether the unaligned reads at `path` carry base modifications: an `MM` tag with its `ML`
/// probabilities on any of its first [`MODIFICATION_SAMPLE`] reads. In a BAM they are the
/// reads' tags; in a FASTQ (plain or gzipped, its sequences on one line or wrapped) they are in
/// the records' comments, as dorado writes them and `samtools fastq -T MM,ML` keeps them, which
/// bonefish's `-y` copies to the alignments. Reads without them give no methylation.
///
/// # Errors
///
/// If the file cannot be read.
pub fn has_modifications(path: &Path, format: Format) -> Result<bool, InputError> {
    let read = |source| InputError::Read { path: path.to_path_buf(), source };
    if format != Format::Fastq {
        let mut reader = bam::io::Reader::new(File::open(path).map_err(read)?);
        let header = reader.read_header().map_err(read)?;
        let [mm, ml] = [Tag::new(b'M', b'M'), Tag::new(b'M', b'L')];
        for record in reader.record_bufs(&header).take(MODIFICATION_SAMPLE) {
            let record = record.map_err(read)?;
            if record.data().get(&mm).is_some() && record.data().get(&ml).is_some() {
                return Ok(true);
            }
        }
        return Ok(false);
    }
    let mut found = false;
    fastq_names(path, |name| {
        let (mut mm, mut ml) = (false, false);
        for tag in name.split(|b| *b == b'\t' || *b == b' ').skip(1) {
            mm |= tag.starts_with(b"MM:Z:") || tag.starts_with(b"Mm:Z:");
            ml |= tag.starts_with(b"ML:B:") || tag.starts_with(b"Ml:B:");
        }
        found = mm && ml;
        found
    })?;
    Ok(found)
}

/// Calls `each` with the name line (the name and its comment) of each of the first
/// [`MODIFICATION_SAMPLE`] records of the FASTQ (plain or gzipped, its sequences on one line or
/// wrapped) at `path`, until it returns true.
fn fastq_names(path: &Path, mut each: impl FnMut(&[u8]) -> bool) -> Result<(), InputError> {
    let read = |source| InputError::Read { path: path.to_path_buf(), source };
    let mut magic = Vec::with_capacity(2);
    File::open(path).and_then(|f| f.take(2).read_to_end(&mut magic)).map_err(read)?;
    let file = File::open(path).map_err(read)?;
    let reader: Box<dyn BufRead> = if magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    // A record: `@` and its name and comment; its bases, on one line or more, up to a `+` line;
    // then as many quality characters as bases, on as many lines.
    let mut lines = reader.split(b'\n');
    for _ in 0..MODIFICATION_SAMPLE {
        let Some(name) = lines.next().transpose().map_err(read)? else { break };
        if each(name.strip_suffix(b"\r").unwrap_or(&name)) {
            return Ok(());
        }
        let mut bases = 0;
        loop {
            match lines.next().transpose().map_err(read)? {
                None => return Ok(()),
                Some(line) if line.starts_with(b"+") => break,
                Some(line) => bases += line.strip_suffix(b"\r").unwrap_or(&line).len(),
            }
        }
        let mut qualities = 0;
        while qualities < bases {
            let Some(line) = lines.next().transpose().map_err(read)? else { return Ok(()) };
            qualities += line.strip_suffix(b"\r").unwrap_or(&line).len();
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use pasta::sam::alignment::io::Write as _;
    use pasta::sam::alignment::record::Flags;
    use pasta::sam::alignment::record_buf::data::field::Value;
    use pasta::sam::alignment::record_buf::data::field::value::Array;
    use pasta::sam::alignment::record_buf::{RecordBuf, Sequence};

    fn dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("lrpipe-input-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch directory");
        dir
    }

    /// A BAM with `header` (SAM text) and unmapped reads, one with `tags` per entry (MM/ML or
    /// not), the first `mapped` of them mapped to the header's first reference.
    fn bam(path: &Path, header: &str, reads: &[&[(&[u8; 2], Value)]], mapped: usize) {
        let header: pasta::sam::Header = header.parse().expect("header");
        let mut writer = bam::io::Writer::new(File::create(path).expect("BAM"));
        writer.write_header(&header).expect("header");
        for (i, tags) in reads.iter().enumerate() {
            let mut record = RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(if i < mapped { Flags::empty() } else { Flags::UNMAPPED })
                .set_sequence(Sequence::from(b"ACGT".to_vec()))
                .set_data(tags.iter().map(|(t, v)| (Tag::new(t[0], t[1]), v.clone())).collect());
            if i < mapped {
                use pasta::sam::alignment::record::cigar::{Op, op::Kind};
                record = record
                    .set_reference_sequence_id(0)
                    .set_alignment_start(pasta::core::Position::MIN)
                    .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect());
            }
            writer.write_alignment_record(&header, &record.build()).expect("record");
        }
        writer.try_finish().expect("finished");
    }

    fn gzip(path: &Path, text: &str) {
        let mut out = flate2::write::GzEncoder::new(File::create(path).expect("file"), flate2::Compression::fast());
        out.write_all(text.as_bytes()).expect("written");
        out.finish().expect("gzipped");
    }

    fn cram(path: &Path, header: &str) {
        let header: pasta::sam::Header = header.parse().expect("header");
        let mut writer = pasta::cram::io::Writer::new(File::create(path).expect("CRAM"));
        writer.write_header(&header).expect("header");
        writer.try_finish(&header).expect("finished");
    }

    fn modified() -> [(&'static [u8; 2], Value); 2] {
        [(b"MM", Value::String("C+m?,0;".into())), (b"ML", Value::Array(Array::UInt8(vec![200])))]
    }

    const SQ: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:26\n";

    #[test]
    fn shared_reads_are_counted_of_the_smaller_sample_and_a_merge_is_suggested() {
        let (a, m) = (PathBuf::from("a.bam"), PathBuf::from("merged.bam"));
        // All of a.bam's 76 sampled reads among merged.bam's 152: the share is a.bam's, named.
        let e = shared_verdict(76, [(&m, 152), (&a, 76)]).expect_err("the same reads").to_string();
        assert!(e.starts_with("76 of the 76 reads sampled from a.bam are in merged.bam too"), "{e}");
        // Half of a sample, as a file and a merge holding it share on a whole genome: a warning
        // that says so.
        let w = shared_verdict(250, [(&m, 500), (&a, 500)]).expect("a warning").expect("said");
        assert!(w.starts_with("250 of the 500 reads sampled from merged.bam") && w.contains("a merge"), "{w}");
        let w = shared_verdict(5, [(&m, 500), (&a, 76)]).expect("a warning").expect("said");
        assert!(w.starts_with("5 of the 76 reads sampled from a.bam") && !w.contains("a merge"), "{w}");
        assert_eq!(shared_verdict(0, [(&m, 500), (&a, 76)]).expect("nothing"), None);
    }

    #[test]
    fn formats_are_told_from_the_contents() {
        let d = dir("formats");
        bam(&d.join("aligned.bam"), SQ, &[&[], &[]], 1);
        bam(&d.join("empty.bam"), SQ, &[], 0);
        bam(&d.join("unaligned.bam"), "@HD\tVN:1.6\n@RG\tID:a\tPL:ONT\n", &[&[]], 0);
        // @SQ lines, but only unmapped reads (as some basecallers and `samtools view -f 4` write).
        bam(&d.join("unmapped.bam"), SQ, &[&[], &[]], 0);
        let fastq = "@r1 qs:f:20\tMM:Z:C+m?,0;\tML:B:C,200\nACGT\n+\n!!!!\n";
        std::fs::write(d.join("reads.fq"), fastq).expect("FASTQ");
        gzip(&d.join("reads.fastq.gz"), fastq);
        // A read named as a SAM header line begins is still a read.
        std::fs::write(d.join("tab.fq"), "@HD\tqs:f:20\nACGT\n+\n!!!!\n").expect("FASTQ");
        cram(&d.join("aligned.cram"), SQ);
        cram(&d.join("unaligned.cram"), "@HD\tVN:1.6\n");
        std::fs::write(d.join("reads.sam"), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:26\nr1\t4\t*\n").expect("SAM");
        std::fs::write(d.join("reads.fa"), ">r1\nACGT\n").expect("FASTA");
        // Named as a BAM, but a gzipped FASTQ: the contents decide.
        gzip(&d.join("named.bam"), fastq);
        for (file, want) in [
            ("aligned.bam", Format::Aligned),
            ("empty.bam", Format::Aligned),
            ("unaligned.bam", Format::UnalignedBam),
            ("unmapped.bam", Format::UnalignedBam),
            ("reads.fq", Format::Fastq),
            ("reads.fastq.gz", Format::Fastq),
            ("tab.fq", Format::Fastq),
            ("aligned.cram", Format::Aligned),
            ("named.bam", Format::Fastq),
        ] {
            assert_eq!(format(&d.join(file)).expect(file), want, "{file}");
        }
        for (file, says) in
            [("reads.sam", "SAM text"), ("reads.fa", "FASTA"), ("unaligned.cram", "bonefish cannot read")]
        {
            let e = format(&d.join(file)).expect_err(file).to_string();
            assert!(e.contains(says) && e.contains("--reads takes"), "{e}");
        }
        assert!(matches!(format(&d.join("gone.bam")), Err(InputError::Read { .. })));
        let gone = classify(&[d.join("gone.bam")], &d.join("ref.fa")).expect_err("missing").to_string();
        assert!(gone.contains("gone.bam does not exist"), "{gone}");

        let p = |f: &str| d.join(f);
        let classify = |paths: &[PathBuf]| classify(paths, &d.join("ref.fa")).map(|c| c.input);
        let not_indexed = classify(&[p("aligned.bam")]).expect_err("no index").to_string();
        assert!(
            not_indexed.contains("aligned.bam has no index") && not_indexed.contains("samtools index"),
            "{not_indexed}"
        );
        let bam::Index::Bai(index) = bam::fs::index(p("aligned.bam")).expect("indexed") else { panic!("a BAI") };
        bam::bai::fs::write(p("aligned.bam.bai"), &index).expect("index written");
        assert_eq!(classify(&[p("aligned.bam")]).expect("one"), Input::Aligned(p("aligned.bam")));
        assert_eq!(
            classify(&[p("unaligned.bam"), p("unmapped.bam")]).expect("unaligned"),
            Input::Unaligned(vec![
                (p("unaligned.bam"), Format::UnalignedBam),
                (p("unmapped.bam"), Format::UnalignedBam)
            ])
        );
        let mixed = classify(&[p("unaligned.bam"), p("aligned.bam")]).expect_err("mixed").to_string();
        assert!(mixed.contains("aligned.bam is aligned and") && mixed.contains("unaligned.bam is not"), "{mixed}");
        // Several aligned files are merged (when indexed).
        let several = classify(&[p("aligned.bam"), p("aligned.cram")]).expect_err("an index missing").to_string();
        assert!(several.contains("aligned.cram has no index") && several.contains(".crai"), "{several}");
        let kinds = classify(&[p("unaligned.bam"), p("reads.fq")]).expect_err("BAM and FASTQ").to_string();
        assert!(kinds.contains("unaligned reads of one kind"), "{kinds}");
        let twice = classify(&[p("reads.fq"), d.join("./reads.fq")]).expect_err("twice").to_string();
        assert!(twice.contains("given twice"), "{twice}");
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn read_groups_must_be_one_samples_and_agree() {
        let d = dir("read-groups");
        let rg = |id: &str, sm: &str, pu: &str| format!("@HD\tVN:1.6\n@RG\tID:{id}\tPU:{pu}\tSM:{sm}\n");
        bam(&d.join("a.bam"), &rg("fc1", "HG002", "FC1"), &[&[]], 0);
        bam(&d.join("b.bam"), &rg("fc2", "HG002", "FC2"), &[&[]], 0);
        // dorado's files of one run, with the same read group.
        bam(&d.join("a2.bam"), &rg("fc1", "HG002", "FC1"), &[&[]], 0);
        bam(&d.join("other.bam"), &rg("fc3", "HG003", "FC3"), &[&[]], 0);
        bam(&d.join("clash.bam"), &rg("fc1", "HG002", "FC9"), &[&[]], 0);
        let p = |f: &str| d.join(f);
        let classify = |paths: &[PathBuf]| classify(paths, &d.join("ref.fa"));
        let read_groups = classify(&[p("a.bam"), p("b.bam"), p("a2.bam")]).expect("one sample").header;
        assert_eq!(read_groups.read_groups().len(), 2, "fc1 once, and fc2");
        let samples = classify(&[p("a.bam"), p("other.bam")]).expect_err("two samples").to_string();
        let two = format!(
            "several samples (@RG SM: HG002 in {}, HG003 in {})",
            p("a.bam").display(),
            p("other.bam").display()
        );
        assert!(samples.contains(&two), "{samples}");
        let clash = classify(&[p("a.bam"), p("clash.bam")]).expect_err("one ID, two lines").to_string();
        assert!(clash.contains("both have read group fc1 (@RG ID:fc1), with different fields"), "{clash}");
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn modifications_need_mm_and_ml_on_the_first_reads() {
        let d = dir("modifications");
        // On the second read, after a wrapped first record whose quality starts with `@`.
        let tagged = "@r1 qs:f:20\nACGT\nAC\n+\n@!!!\n!!\n@r2\tqs:f:20\tMM:Z:C+m?,0;\tML:B:C,200\nACGT\n+\n!!!!\n";
        // `MM:Z:` in a sequence or quality line, or a read's name, is not a tag; MM without ML
        // gives no probabilities.
        let untagged = "@MM:Z:x\nMM:Z:\n+\n@MM:Z\n@r2\tMM:Z:C+m?,0;\nACGT\n+\n!!!!\n";
        std::fs::write(d.join("tagged.fq"), tagged).expect("FASTQ");
        gzip(&d.join("tagged.fq.gz"), tagged);
        std::fs::write(d.join("untagged.fq"), untagged).expect("FASTQ");
        bam(&d.join("tagged.bam"), "@HD\tVN:1.6\n", &[&[], &modified()], 0);
        bam(&d.join("untagged.bam"), "@HD\tVN:1.6\n", &[&[], &modified()[..1]], 0);
        for (file, format, want) in [
            ("tagged.fq", Format::Fastq, true),
            ("tagged.fq.gz", Format::Fastq, true),
            ("untagged.fq", Format::Fastq, false),
            ("tagged.bam", Format::UnalignedBam, true),
            ("untagged.bam", Format::UnalignedBam, false),
        ] {
            assert_eq!(has_modifications(&d.join(file), format).expect(file), want, "{file}");
        }
        let _ = std::fs::remove_dir_all(d);
    }

    fn header(text: &str) -> sam::Header {
        format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n{text}").parse().expect("a header")
    }

    /// Every input's `@PG` chain is kept whole: a line of the same ID and fields once, one of the
    /// same ID with other fields renamed, with the line after it following the new name.
    #[test]
    fn every_inputs_programs_are_kept_as_chains() {
        let first =
            header("@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 a\n@PG\tID:samtools\tPN:samtools\tPP:minimap2\n");
        let second = header(
            "@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 b\n@PG\tID:samtools\tPN:samtools\tPP:minimap2\n\
             @PG\tID:dorado\tPN:dorado\n",
        );
        let merged = combine_headers(&[("a.bam".into(), first.clone()), ("b.bam".into(), second)]).expect("merged");
        let lines: Vec<String> = merged
            .programs()
            .as_ref()
            .iter()
            .map(|(id, p)| {
                let pp = p.other_fields().get(&sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID);
                format!("{id}<{}", pp.map(ToString::to_string).unwrap_or_default())
            })
            .collect();
        assert_eq!(lines, ["minimap2<", "samtools<minimap2", "minimap2-2<", "samtools-2<minimap2-2", "dorado<"]);
        // The same file twice: its lines once.
        let twice = combine_headers(&[("a.bam".into(), first.clone()), ("c.bam".into(), first)]).expect("merged");
        assert_eq!(twice.programs().as_ref().len(), 2);
    }
}
