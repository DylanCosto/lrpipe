//! FASTA/FASTQ reading as pysam's `FastxFile` does it: htslib's `kseq.h`
//! parser over a plain or gzip-compressed file, and the record formatting
//! of `FastxRecord.__str__`.

use std::fs::File;
use std::io::{BufReader, Read};

use anyhow::{Context, Result, bail};

const BUFSIZE: usize = 16384;

/// `kstream_t`: a buffered byte stream
struct KStream {
    f: Box<dyn Read>,
    buf: Vec<u8>,
    begin: usize,
    end: usize,
    is_eof: bool,
    err: bool,
}

/// Result of `ks_getuntil2`: data was read (length >= 0), or end of file (-1)
enum GetUntil {
    Data,
    Eof,
}

impl KStream {
    fn fill(&mut self) -> Result<()> {
        self.begin = 0;
        let mut n = 0;
        // bgzf_read fills as much as it can
        while n < BUFSIZE {
            match self.f.read(&mut self.buf[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => {
                    self.err = true;
                    self.is_eof = true;
                    return Err(e.into());
                }
            }
        }
        self.end = n;
        if n == 0 {
            self.is_eof = true;
        }
        Ok(())
    }

    fn at_eof(&self) -> bool {
        self.is_eof && self.begin >= self.end
    }

    /// `ks_getc`: the next byte, or None at the end
    fn getc(&mut self) -> Result<Option<u8>> {
        if self.err {
            bail!("error reading stream");
        }
        if self.at_eof() {
            return Ok(None);
        }
        if self.begin >= self.end {
            self.fill()?;
            if self.end == 0 {
                return Ok(None);
            }
        }
        self.begin += 1;
        Ok(Some(self.buf[self.begin - 1]))
    }

    /// `ks_getuntil2`: read up to a delimiter (`line`: '\n', otherwise
    /// whitespace) into `s`, returning the delimiter in `dret` (0 if none).
    fn getuntil(&mut self, line: bool, s: &mut Vec<u8>, dret: &mut u8, append: bool) -> Result<GetUntil> {
        let mut gotany = false;
        *dret = 0;
        if !append {
            s.clear();
        }
        loop {
            if self.err {
                bail!("error reading stream");
            }
            if self.begin >= self.end {
                if !self.is_eof {
                    self.fill()?;
                    if self.end == 0 {
                        break;
                    }
                } else {
                    break;
                }
            }
            let chunk = &self.buf[self.begin..self.end];
            let found = if line {
                chunk.iter().position(|&c| c == b'\n')
            } else {
                // C isspace
                chunk.iter().position(|&c| c == b' ' || (b'\t'..=b'\r').contains(&c))
            };
            let i = found.map_or(self.end, |k| self.begin + k);
            s.extend_from_slice(&self.buf[self.begin..i]);
            gotany = true;
            self.begin = i + 1;
            if i < self.end {
                *dret = self.buf[i];
                break;
            }
        }
        if !gotany && self.at_eof() {
            return Ok(GetUntil::Eof);
        }
        if line && s.len() > 1 && s[s.len() - 1] == b'\r' {
            s.pop();
        }
        Ok(GetUntil::Data)
    }
}

/// A FASTA/FASTQ record (`kseq_t` after `kseq_read`)
#[derive(Default)]
struct KSeq {
    name: Vec<u8>,
    comment: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    last_char: u8,
}

/// A record as pysam's `FastxRecord` holds it
pub struct FastxRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: String,
    pub quality: Option<String>,
}

impl std::fmt::Display for FastxRecord {
    /// `FastxRecord.__str__`
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let comment = self.comment.as_ref().map_or(String::new(), |c| format!(" {c}"));
        match &self.quality {
            None => write!(f, ">{}{comment}\n{}", self.name, self.sequence),
            Some(q) => write!(f, "@{}{comment}\n{}\n+\n{q}", self.name, self.sequence),
        }
    }
}

pub struct FastxFile {
    ks: KStream,
    seq: KSeq,
    path: String,
}

/// pysam's `charptr_to_str`: the C string (up to the first NUL) as UTF-8
fn c_str(bytes: &[u8]) -> Result<String> {
    let end = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
    String::from_utf8(bytes[..end].to_vec()).context("UnicodeDecodeError: 'utf-8' codec can't decode the FASTQ/FASTA")
}

impl FastxFile {
    /// Opens a plain or gzip-compressed FASTA/FASTQ file (`bgzf_open`).
    pub fn open(path: &str) -> Result<FastxFile> {
        let mut file = BufReader::new(File::open(path).with_context(|| format!("file `{path}` not found"))?);
        let gz = std::io::BufRead::fill_buf(&mut file)?.starts_with(b"\x1f\x8b");
        let f: Box<dyn Read> = if gz { Box::new(flate2::bufread::MultiGzDecoder::new(file)) } else { Box::new(file) };
        Ok(FastxFile {
            ks: KStream { f, buf: vec![0; BUFSIZE], begin: 0, end: 0, is_eof: false, err: false },
            seq: KSeq::default(),
            path: path.to_string(),
        })
    }

    /// `kseq_read`: Ok(Some(())) for a record, Ok(None) at the end of file,
    /// Err(-2) for a truncated quality string, Err(-3) for a read error.
    fn kseq_read(&mut self) -> Result<Option<()>, i32> {
        let ks = &mut self.ks;
        let seq = &mut self.seq;
        let getc = |ks: &mut KStream| ks.getc().map_err(|_| -3);
        if seq.last_char == 0 {
            // jump to the next header line
            loop {
                match getc(ks)? {
                    None => return Ok(None),
                    Some(b'>' | b'@') => {
                        seq.last_char = b'>';
                        break;
                    }
                    Some(_) => {}
                }
            }
        }
        seq.comment.clear();
        seq.seq.clear();
        seq.qual.clear();
        let mut c = 0;
        if let GetUntil::Eof = ks.getuntil(false, &mut seq.name, &mut c, false).map_err(|_| -3)? {
            return Ok(None);
        }
        if c != b'\n' {
            ks.getuntil(true, &mut seq.comment, &mut 0, false).map_err(|_| -3)?;
        }
        let mut c = getc(ks)?;
        while let Some(ch) = c {
            if ch == b'>' || ch == b'+' || ch == b'@' {
                break;
            }
            if ch != b'\n' {
                seq.seq.push(ch);
                ks.getuntil(true, &mut seq.seq, &mut 0, true).map_err(|_| -3)?;
            }
            c = getc(ks)?;
        }
        if matches!(c, Some(b'>' | b'@')) {
            seq.last_char = b'>';
        }
        if c != Some(b'+') {
            return Ok(Some(()));
        }
        // skip the rest of the '+' line
        loop {
            match getc(ks)? {
                None => return Err(-2),
                Some(b'\n') => break,
                Some(_) => {}
            }
        }
        loop {
            match ks.getuntil(true, &mut seq.qual, &mut 0, true).map_err(|_| -3)? {
                GetUntil::Eof => break,
                GetUntil::Data if seq.qual.len() >= seq.seq.len() => break,
                GetUntil::Data => {}
            }
        }
        seq.last_char = 0;
        if seq.seq.len() != seq.qual.len() {
            return Err(-2);
        }
        Ok(Some(()))
    }

    /// The next record (`FastxFile.__next__`)
    pub fn next_record(&mut self) -> Result<Option<FastxRecord>> {
        match self.kseq_read() {
            Ok(None) => Ok(None),
            Ok(Some(())) => {
                let s = &self.seq;
                Ok(Some(FastxRecord {
                    name: c_str(&s.name)?,
                    comment: if s.comment.is_empty() { None } else { Some(c_str(&s.comment)?) },
                    sequence: c_str(&s.seq)?,
                    quality: if s.qual.is_empty() { None } else { Some(c_str(&s.qual)?) },
                }))
            }
            Err(-2) => bail!("truncated quality string in {}", self.path),
            Err(_) => bail!("unknown problem parsing {}", self.path),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &[u8]) -> Vec<String> {
        let mut f = FastxFile {
            ks: KStream {
                f: Box::new(std::io::Cursor::new(text.to_vec())),
                buf: vec![0; BUFSIZE],
                begin: 0,
                end: 0,
                is_eof: false,
                err: false,
            },
            seq: KSeq::default(),
            path: String::new(),
        };
        let mut out = Vec::new();
        while let Some(r) = f.next_record().unwrap() {
            out.push(r.to_string());
        }
        out
    }

    #[test]
    fn fastq_and_fasta() {
        assert_eq!(
            parse(b"@r1 c d\nACGT\n+\nIIII\n@r2\nAC\nGT\n+r2\nII\nII\n"),
            vec!["@r1 c d\nACGT\n+\nIIII", "@r2\nACGT\n+\nIIII"]
        );
        assert_eq!(parse(b">s1\r\nAC\r\nGT\r\n>s2 x\n\nA\n"), vec![">s1\nACGT", ">s2 x\nA"]);
        // an empty quality string makes it a FASTA record
        assert_eq!(parse(b"@e\n\n+\n\n"), vec![">e\n"]);
    }
}
