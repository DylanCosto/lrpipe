//! Port of `bam.py`: reading alignments of one sample from one or more
//! indexed BAM/CRAM files, merged by position.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};

use crate::alignment::{AlignmentFile, Aux, Decoding, Fetch, Record};

#[derive(Debug)]
pub enum FetchError {
    SampleNotFound,
    ReferenceNotFound(String),
    Other(anyhow::Error),
}

impl From<anyhow::Error> for FetchError {
    fn from(e: anyhow::Error) -> Self {
        FetchError::Other(e)
    }
}

pub struct SampleBamReader {
    pub source_id: usize,
    path: String,
    reader: AlignmentFile,
    references: HashSet<String>,
    sample_to_group_ids: HashMap<String, HashSet<Vec<u8>>>,
}

impl SampleBamReader {
    /// `threads`: threads decompressing BAM files.
    pub fn open(path: &str, source_id: usize, reference: Option<&str>, threads: usize) -> Result<Self> {
        let reference = reference.map(|r| std::fs::canonicalize(r).unwrap_or_else(|_| Path::new(r).to_path_buf()));
        let mut reader =
            AlignmentFile::open(path, reference.as_deref(), threads).with_context(|| format!("opening {path}"))?;
        // Reads are only read for their alleles here.
        reader.set_decoding(Decoding::Alleles);
        if !reader.load_index()? {
            bail!(
                "The file '{path}' is not indexed. Please create the appropriate BAM/CRAM index with \"samtools index\""
            );
        }
        let header = reader.header().clone();
        let references: HashSet<String> = header.target_names().map(String::from).collect();
        // Empty files are treated as errors (they usually indicate a CRAM whose
        // reference could not be found).
        reader.fetch(Fetch::All)?;
        if !matches!(reader.read(), Ok(Some(_))) {
            bail!(
                "No reads could be retrieved from '{path}'. If this is a CRAM file, possibly the reference could not be found. Try to use --reference=... or check your $REF_PATH/$REF_CACHE settings"
            );
        }
        let mut sample_to_group_ids: HashMap<String, HashSet<Vec<u8>>> = HashMap::new();
        if let Some(rgs) = header.records().get("RG") {
            for rg in rgs {
                let get = |k: &str| rg.iter().find(|(t, _)| t == k).map(|(_, v)| v.clone());
                let id = get("ID").unwrap_or_default();
                match get("SM") {
                    Some(sm) => {
                        sample_to_group_ids.entry(sm).or_default().insert(id.into_bytes());
                    }
                    None => log::warn!(
                        "Read group \"{id}\" does not contain an SM field to assign it to a sample. Use --ignore-read-groups to use these alignments anyway."
                    ),
                }
            }
        }
        Ok(SampleBamReader { source_id, path: path.to_string(), reader, references, sample_to_group_ids })
    }

    pub fn is_cram(&self) -> bool {
        self.reader.is_cram()
    }

    pub fn has_reference(&self, name: &str) -> bool {
        self.references.contains(name)
    }

    pub fn has_sample(&self, sample: &str) -> bool {
        self.sample_to_group_ids.contains_key(sample)
    }

    fn start_fetch(
        &mut self,
        reference: &str,
        region: (i64, Option<i64>),
        sample: Option<&str>,
    ) -> Result<Option<HashSet<Vec<u8>>>, FetchError> {
        if !self.references.contains(reference) {
            return Err(FetchError::ReferenceNotFound(reference.to_string()));
        }
        let groups = match sample {
            None => None,
            Some(s) => Some(self.sample_to_group_ids.get(s).ok_or(FetchError::SampleNotFound)?.clone()),
        };
        let tid = self.reader.header().tid(reference).expect("reference is in the header");
        // pysam uses MAX_POS = 2^31 - 1 as the end of open intervals
        let (start, end) = (region.0, region.1.unwrap_or((1 << 31) - 1));
        self.reader
            .fetch(Fetch::Region(tid, start, end))
            .map_err(|e| FetchError::Other(anyhow!("{e} while fetching {reference} from {}", self.path)))?;
        Ok(groups)
    }

    /// The next alignment (of the selected read groups).
    fn next_record(&mut self, groups: &Option<HashSet<Vec<u8>>>) -> Result<Option<Record>, FetchError> {
        loop {
            let rec = match self.reader.read() {
                Ok(None) => return Ok(None),
                Err(e) => return Err(FetchError::Other(anyhow!("{e} while reading {}", self.path))),
                Ok(Some(rec)) => rec,
            };
            if let Some(g) = groups {
                let keep = match rec.aux(b"RG") {
                    Ok(Aux::String(s)) => g.contains(s.as_bytes()),
                    _ => {
                        return Err(FetchError::Other(anyhow!(
                            "tag 'RG' not present in read {}",
                            String::from_utf8_lossy(rec.qname())
                        )));
                    }
                };
                if !keep {
                    continue;
                }
            }
            return Ok(Some(rec));
        }
    }
}

/// Reads from one or more BAM/CRAM files, merging them by position
/// (`SampleBamReader` / `MultiBamReader`).
pub struct BamReader {
    readers: Vec<SampleBamReader>,
}

impl BamReader {
    /// `threads`: threads decompressing each BAM file.
    pub fn open(paths: &[String], reference: Option<&str>, threads: usize) -> Result<Self> {
        let mut readers = Vec::new();
        for (i, p) in paths.iter().enumerate() {
            readers.push(SampleBamReader::open(p, i, reference, threads)?);
        }
        Ok(BamReader { readers })
    }

    pub fn n_paths(&self) -> usize {
        self.readers.len()
    }

    pub fn has_reference(&self, name: &str) -> bool {
        self.readers.iter().all(|r| r.has_reference(name))
    }

    /// Number of CRAM files (each keeps its own copy of the reference
    /// sequence while decoding)
    pub fn cram_files(&self) -> usize {
        self.readers.iter().filter(|r| r.is_cram()).count()
    }

    /// Decodes every CRAM file with `repository` (see
    /// [`crate::alignment::cram::CramReader::share_reference`]).
    pub fn share_reference(&mut self, repository: &pasta::fasta::Repository) {
        for r in &mut self.readers {
            r.reader.share_reference(repository);
        }
    }

    /// Calls `f(source_id, record)` for all alignments on `reference`
    /// (optionally restricted to `sample`), merged across files by
    /// (reference_start, source_id) as `heapq.merge` does in WhatsHap.
    pub fn for_each(
        &mut self,
        reference: &str,
        sample: Option<&str>,
        f: impl FnMut(usize, &Record) -> Result<()>,
    ) -> Result<(), FetchError> {
        self.for_each_in_regions(reference, &[(0, None)], sample, f)
    }

    /// Like `for_each`, but for each region in turn (alignments overlapping
    /// several regions are reported several times, as in WhatsHap).
    pub fn for_each_in_regions(
        &mut self,
        reference: &str,
        regions: &[(i64, Option<i64>)],
        sample: Option<&str>,
        mut f: impl FnMut(usize, &Record) -> Result<()>,
    ) -> Result<(), FetchError> {
        for &region in regions {
            self.for_each_region(reference, region, sample, &mut f)?;
        }
        Ok(())
    }

    fn for_each_region(
        &mut self,
        reference: &str,
        region: (i64, Option<i64>),
        sample: Option<&str>,
        f: &mut impl FnMut(usize, &Record) -> Result<()>,
    ) -> Result<(), FetchError> {
        let multi = self.readers.len() > 1;
        let mut active: Vec<(usize, Option<HashSet<Vec<u8>>>)> = Vec::new();
        for (k, r) in self.readers.iter_mut().enumerate() {
            if multi && sample.is_some() && !r.has_sample(sample.unwrap()) {
                continue;
            }
            let groups = r.start_fetch(reference, region, sample)?;
            active.push((k, groups));
        }
        if active.is_empty() {
            return Err(FetchError::SampleNotFound);
        }
        let mut current: Vec<Option<Record>> = Vec::new();
        for (k, groups) in &active {
            current.push(self.readers[*k].next_record(groups)?);
        }
        loop {
            // pick the smallest (pos, source_id); ties cannot occur across files
            let mut best: Option<usize> = None;
            for (i, c) in current.iter().enumerate() {
                if let Some(rec) = c {
                    let key = (rec.pos(), self.readers[active[i].0].source_id);
                    match best {
                        None => best = Some(i),
                        Some(b) => {
                            let brec = current[b].as_ref().unwrap();
                            if key < (brec.pos(), self.readers[active[b].0].source_id) {
                                best = Some(i);
                            }
                        }
                    }
                }
            }
            let Some(i) = best else { break };
            let k = active[i].0;
            let sid = self.readers[k].source_id;
            f(sid, current[i].as_ref().unwrap())?;
            current[i] = self.readers[k].next_record(&active[i].1)?;
        }
        Ok(())
    }
}
