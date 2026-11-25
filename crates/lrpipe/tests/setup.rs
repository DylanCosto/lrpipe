//! `lrpipe setup` against tiny files: served by a local HTTP server the tests start (which can
//! corrupt or cut off a response, and logs each request's Range), or from `file://` URLs.

use std::collections::{HashMap, VecDeque};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Arc, Mutex, PoisonError, RwLock, RwLockReadGuard, RwLockWriteGuard};
use std::time::Duration;

use lrpipe::config::Site;
use lrpipe::plan::Platform;
use lrpipe::setup::fetch::Fetcher;
use lrpipe::setup::{self, Check, Options};

fn temp_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("lrpipe-setup-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir.canonicalize().expect("canonical temp dir")
}

/// Tests that run programs hold this exclusively, the others shared: a setup holds its
/// directory's lock with `flock`, and a child another test forks keeps a copy of that lock until
/// it execs, long enough on a loaded machine for the next setup of that directory to find it busy.
static FORKS: RwLock<()> = RwLock::new(());

fn shared() -> RwLockReadGuard<'static, ()> {
    FORKS.read().unwrap_or_else(PoisonError::into_inner)
}

fn exclusive() -> RwLockWriteGuard<'static, ()> {
    FORKS.write().unwrap_or_else(PoisonError::into_inner)
}

fn fetcher() -> Fetcher {
    Fetcher::new().with_pause(Duration::from_millis(10))
}

/// What the server does to its next response for a path.
#[derive(Clone, Copy, Debug)]
enum Fault {
    /// One byte changed.
    Corrupt,
    /// The headers promise the whole, and the connection closes after this many bytes.
    CutOff(usize),
    /// A Range request is answered with the whole file (200).
    IgnoreRange,
    /// A Range request is answered with 206 and the whole file, as if from byte 0.
    WrongRange,
}

/// A request: its path and the first byte its Range asked for.
type Request = (String, Option<u64>);

#[derive(Default)]
struct Shared {
    files: HashMap<String, Vec<u8>>,
    faults: HashMap<String, VecDeque<Fault>>,
    /// Each file's ETag, where it has one: If-Range is honoured against it.
    etags: HashMap<String, String>,
    log: Vec<Request>,
    /// Each request's path and If-Range.
    if_range: Vec<(String, Option<String>)>,
}

/// An HTTP/1.1 server on 127.0.0.1 serving `files`, with Range requests.
struct Server {
    port: u16,
    shared: Arc<Mutex<Shared>>,
}

impl Server {
    fn start(files: &[(&str, &[u8])]) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().expect("address").port();
        let shared = Arc::new(Mutex::new(Shared {
            files: files.iter().map(|(name, data)| (format!("/{name}"), data.to_vec())).collect(),
            ..Shared::default()
        }));
        let state = Arc::clone(&shared);
        std::thread::spawn(move || {
            for stream in listener.incoming().flatten() {
                let state = Arc::clone(&state);
                std::thread::spawn(move || serve(stream, &state));
            }
        });
        Self { port, shared }
    }

    fn url(&self, name: &str) -> String {
        format!("http://127.0.0.1:{}/{name}", self.port)
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Shared> {
        self.shared.lock().unwrap_or_else(PoisonError::into_inner)
    }

    fn fault(&self, name: &str, fault: Fault) {
        self.lock().faults.entry(format!("/{name}")).or_default().push_back(fault);
    }

    /// The requests so far, and forgets them.
    fn take_log(&self) -> Vec<Request> {
        std::mem::take(&mut self.lock().log)
    }

    /// Serves `data` at `name` from now on, as version `etag`.
    fn publish(&self, name: &str, data: &[u8], etag: &str) {
        let mut shared = self.lock();
        shared.files.insert(format!("/{name}"), data.to_vec());
        shared.etags.insert(format!("/{name}"), format!("\"{etag}\""));
    }

    /// Each request's If-Range so far for `name`.
    fn if_ranges(&self, name: &str) -> Vec<Option<String>> {
        self.lock().if_range.iter().filter(|(p, _)| *p == format!("/{name}")).map(|(_, v)| v.clone()).collect()
    }
}

fn serve(stream: TcpStream, state: &Mutex<Shared>) {
    let mut reader = BufReader::new(stream.try_clone().expect("clone"));
    let mut line = String::new();
    if reader.read_line(&mut line).is_err() {
        return;
    }
    let path = line.split_whitespace().nth(1).unwrap_or_default().to_owned();
    let (mut range, mut if_range) = (None, None);
    loop {
        let mut header = String::new();
        if reader.read_line(&mut header).is_err() || header.trim().is_empty() {
            break;
        }
        let lower = header.to_ascii_lowercase();
        if let Some(value) = lower.strip_prefix("range: bytes=") {
            range = value.trim().trim_end_matches('-').parse::<u64>().ok();
        }
        if lower.starts_with("if-range:") {
            if_range = header.split_once(':').map(|(_, v)| v.trim().to_owned());
        }
    }
    let (data, fault, etag) = {
        let mut shared = state.lock().unwrap_or_else(PoisonError::into_inner);
        shared.log.push((path.clone(), range));
        shared.if_range.push((path.clone(), if_range.clone()));
        let fault = shared.faults.get_mut(&path).and_then(VecDeque::pop_front);
        (shared.files.get(&path).cloned(), fault, shared.etags.get(&path).cloned())
    };
    // If-Range of another version (or of none known): the whole file.
    if if_range.is_some() && if_range != etag {
        range = None;
    }
    if let Some(Fault::IgnoreRange) = fault {
        range = None;
    }
    let etag_header = etag.map_or_else(String::new, |e| format!("ETag: {e}\r\n"));
    let mut out = stream;
    let Some(mut data) = data else {
        let _ = out.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
        return;
    };
    if let Some(Fault::Corrupt) = fault {
        let last = data.len() - 1;
        data[last] ^= 0xff;
    }
    let len = data.len() as u64;
    let (status, body, extra) = match range {
        Some(start) if start >= len => {
            let _ = write!(
                out,
                "HTTP/1.1 416 Range Not Satisfiable\r\nContent-Range: bytes */{len}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
            );
            return;
        }
        Some(_) if matches!(fault, Some(Fault::WrongRange)) => {
            ("206 Partial Content", &data[..], format!("Content-Range: bytes 0-{}/{len}\r\n", len - 1))
        }
        Some(start) => (
            "206 Partial Content",
            &data[start as usize..],
            format!("Content-Range: bytes {start}-{}/{len}\r\n", len - 1),
        ),
        None => ("200 OK", &data[..], String::new()),
    };
    let _ = write!(
        out,
        "HTTP/1.1 {status}\r\nContent-Length: {}\r\n{extra}{etag_header}Connection: close\r\n\r\n",
        body.len()
    );
    let sent = match fault {
        Some(Fault::CutOff(n)) => &body[..n.min(body.len())],
        _ => body,
    };
    let _ = out.write_all(sent);
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).expect("compressed");
    e.finish().expect("compressed")
}

fn sha256(data: &[u8]) -> String {
    use sha2::Digest as _;
    lrpipe::provenance::hex(&sha2::Sha256::digest(data))
}

const FASTA: &[u8] = b">chr1\nACGTNNNNNNACGTAAAAAA\nACGT\n>chr2\nGGGGCCCC\n";
const FAI: &str = "chr1\t24\t6\t20\t21\nchr2\t8\t38\t8\t9\n";
const ADOTTO: &[u8] = b"chr2\t50\t60\tID=b;MOTIFS=A\nchr1\t30\t40\tID=a;MOTIFS=CA\nchr10\t5\t9\tID=c;MOTIFS=G\n";
const STRCHIVE: &[u8] = b"chr1\t7\t12\tID=p;MOTIFS=CAG\n";
const CATALOG: &str = "chr1\t7\t12\tID=p;MOTIFS=CAG\nchr1\t30\t40\tID=a;MOTIFS=CA\nchr10\t5\t9\tID=c;MOTIFS=G\nchr2\t50\t60\tID=b;MOTIFS=A\n";
const TRF: &[u8] = b"chr1\t0\t10\nchr2\t0\t5\n";
const BLACKLIST: &[u8] = b"chr1\t0\t3\tgap\n";

/// The tiny genome's files, by name as served.
fn published() -> Vec<(&'static str, Vec<u8>)> {
    vec![
        ("ref.fa.gz", gzip(FASTA)),
        ("adotto.bed.gz", gzip(ADOTTO)),
        ("strchive.bed", STRCHIVE.to_vec()),
        ("trf.bed", TRF.to_vec()),
        ("blacklist.bed", BLACKLIST.to_vec()),
        ("ont/pileup.pt", b"ont pileup".to_vec()),
        ("ont/full_alignment.pt", b"ont full alignment".to_vec()),
        ("hifi/pileup.pt", b"hifi pileup".to_vec()),
        ("hifi/full_alignment.pt", b"hifi full alignment".to_vec()),
    ]
}

/// A manifest of the tiny genome, each file at `url(name)` with its size and sha256 (or, for
/// those in `unpinned`, none), and the files in `accept` needing their terms accepted.
fn manifest(url: &dyn Fn(&str) -> String, unpinned: &[&str], accept: &[&str]) -> String {
    let files = published();
    let entry = |id: &str, name: &str, path: &str, extra: &str| {
        let data = &files.iter().find(|(n, _)| *n == name).expect("published").1;
        let checksum = if unpinned.contains(&id) {
            "unpinned = \"none published\"".to_owned()
        } else {
            format!("sha256 = \"{}\"", sha256(data))
        };
        let terms = if accept.contains(&id) { "accept = true\nterms = \"research use only\"\n" } else { "" };
        format!(
            "[[file]]\nid = \"{id}\"\nurl = \"{}\"\npath = \"{path}\"\nsize = {}\n{checksum}\nlicence = \"test licence {id}\"\n{terms}{extra}\n",
            url(name),
            data.len()
        )
    };
    let mut text = String::from(
        "format = 1\ngenome = \"tiny\"\nmissing = [\"nothing else\"]\n\n[site]\nreference = \"reference\"\ntandem_repeats = \"trf\"\n\
         catalog = \"catalog\"\nblacklist = \"blacklist\"\nmetadata = \"tiny.mdr\"\n\n[models]\nont = \"m_ont\"\nhifi = \"m_hifi\"\n\n",
    );
    text += &entry("reference", "ref.fa.gz", "tiny.fa", "unpack = \"gunzip\"\n");
    text += &entry("adotto", "adotto.bed.gz", "sources/adotto.bed.gz", "");
    text += &entry("strchive", "strchive.bed", "sources/strchive.bed", "");
    text += &entry("trf", "trf.bed", "trf.bed", "");
    text += &entry("blacklist", "blacklist.bed", "blacklist.bed", "");
    for (platform, model) in [("ont", "m_ont"), ("hifi", "m_hifi")] {
        for file in ["pileup", "full_alignment"] {
            let extra = format!("model = \"{model}\"\nplatform = \"{platform}\"\n");
            text += &entry(
                &format!("{model}/{file}"),
                &format!("{platform}/{file}.pt"),
                &format!("models/{model}/{file}.pt"),
                &extra,
            );
        }
    }
    text += &format!(
        "[[build]]\nid = \"catalog\"\nkind = \"catalog\"\ninputs = [\"adotto\", \"strchive\"]\npath = \"catalog.bed\"\nsha256 = \"{}\"\n",
        sha256(CATALOG.as_bytes())
    );
    text
}

fn serve_tiny() -> Server {
    let files = published();
    let refs: Vec<(&str, &[u8])> = files.iter().map(|(n, d)| (*n, d.as_slice())).collect();
    Server::start(&refs)
}

fn options(dir: &Path, platforms: &[Platform]) -> Options {
    Options { dir: dir.to_path_buf(), platforms: platforms.to_vec(), accept_licences: false, bonefish: None }
}

#[test]
fn setup_fetches_checks_builds_and_writes_a_site_configuration_lrpipe_run_reads() {
    // It copies the directory with cp.
    let _forks = exclusive();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("full");
    let both = [Platform::Ont, Platform::Hifi];
    let site_path = setup::setup(&text, &options(&dir, &both), &fetcher()).expect("set up");
    let dir = fs::canonicalize(&dir).expect("canonical");
    assert_eq!(site_path, dir.join("site.toml"));

    // The site configuration lrpipe run reads, every path absolute and there.
    let site = Site::load(&site_path).expect("lrpipe run's parser takes it");
    assert_eq!(site.reference, dir.join("tiny.fa"));
    assert!(site.tools.bonefish.is_none() && site.align.index.is_none(), "no bonefish index unless asked for");
    assert_eq!(fs::read(&site.reference).expect("unpacked"), FASTA);
    assert!(!dir.join("sources/tiny.fa.gz").exists(), "the download is removed once unpacked");
    assert_eq!(fs::read_to_string(dir.join("tiny.fa.fai")).expect(".fai built"), FAI);
    // The models directory, for lrpipe run to choose from by the reads' header; each platform's
    // default only as a commented-out override.
    assert_eq!(site.snv.models, Some(dir.join("models")));
    assert!(site.snv.ont_model.is_none() && site.snv.hifi_model.is_none());
    assert!(offered(&dir, "ont", "m_ont") && offered(&dir, "hifi", "m_hifi"));
    for model in ["m_ont", "m_hifi"] {
        for file in ["pileup.pt", "full_alignment.pt"] {
            assert!(dir.join("models").join(model).join(file).is_file(), "{model}/{file}");
        }
    }
    assert_eq!(site.sv.tandem_repeats, Some(dir.join("trf.bed")));
    let catalog = site.repeats.expect("a catalogue").catalog;
    assert_eq!(catalog, dir.join("catalog.bed"));
    assert_eq!(fs::read_to_string(&catalog).expect("built"), CATALOG);
    assert_eq!(site.cnv.blacklist, Some(dir.join("blacklist.bed")));
    let mdr = site.cnv.metadata.expect("an .mdr");
    assert!(mdr.is_file(), "made by lrpipe-cnv");
    let text_written = fs::read_to_string(&site_path).expect("site.toml");
    for licence in ["test licence reference", "test licence trf", "test licence m_ont/pileup"] {
        assert!(text_written.contains(licence), "the licences are recorded: {licence}\n{text_written}");
    }
    assert!(text_written.contains("nothing else"), "what is missing is said");

    let verification = setup::verify(&text, &dir, &both).expect("verified");
    assert!(verification.passed(), "{verification}");

    // A rerun fetches nothing and changes nothing, keeping no .previous of what it did not change.
    server.take_log();
    let before = fs::read_to_string(&site_path).expect("site.toml");
    for _ in 0..2 {
        assert_eq!(setup::setup(&text, &options(&dir, &both), &fetcher()).expect("set up again"), site_path);
    }
    assert_eq!(server.take_log(), Vec::<Request>::new(), "a rerun fetches nothing");
    assert_eq!(fs::read_to_string(&site_path).expect("site.toml"), before);
    assert!(!dir.join("site.toml.previous").exists(), "nothing changed to keep");

    // A copy of the directory: its site.toml names the original's files, which --verify says.
    let copy = temp_dir("full-copy");
    let status = Command::new("cp").arg("-a").arg(dir.join(".")).arg(&copy).status().expect("cp runs");
    assert!(status.success());
    let verification = setup::verify(&text, &copy, &both).expect("verified");
    assert!(!verification.passed(), "{verification}");
    let check = |v: &setup::Verification, path: &str| v.files.iter().find(|(p, _)| p == path).map(|(_, c)| c.clone());
    assert!(
        matches!(check(&verification, "site.toml"), Some(Check::Mismatch(why)) if why.contains("not this setup's")),
        "{verification}"
    );
    // Its setup writes the copy's own, the original's kept as .previous.
    setup::setup(&text, &options(&copy, &both), &fetcher()).expect("set up the copy");
    assert!(setup::verify(&text, &copy, &both).expect("verified").passed());
    assert!(copy.join("site.toml.previous").exists());
    let _ = fs::remove_dir_all(copy);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn only_the_chosen_platforms_models_are_fetched() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("platform");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    let site = Site::load(dir.join("site.toml")).expect("site.toml");
    assert_eq!(site.snv.models, Some(dir.join("models")));
    assert!(offered(&dir, "ont", "m_ont") && !offered(&dir, "hifi", "m_hifi"));
    assert!(!dir.join("models/m_hifi").exists());
    assert!(server.take_log().iter().all(|(path, _)| !path.starts_with("/hifi/")));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_corrupted_download_is_fetched_again() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("corrupt");
    server.fault("blacklist.bed", Fault::Corrupt);
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    let requests: Vec<Request> = server.take_log().into_iter().filter(|(p, _)| p == "/blacklist.bed").collect();
    assert_eq!(
        requests,
        [("/blacklist.bed".to_owned(), None), ("/blacklist.bed".to_owned(), None)],
        "fetched whole twice"
    );
    assert_eq!(fs::read(dir.join("blacklist.bed")).expect("fetched"), BLACKLIST);

    // A file corrupted after setup: a rerun sees it changed, checks it and fetches it again.
    fs::write(dir.join("trf.bed"), b"chr1\t0\t10\nchr2\t0\t6\n").expect("corrupted");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up again");
    assert_eq!(server.take_log(), [("/trf.bed".to_owned(), None)]);
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched again"), TRF);

    // One corrupted with its size and time kept: a rerun trusts the record, --verify does not.
    let path = dir.join("blacklist.bed");
    let modified = fs::metadata(&path).and_then(|m| m.modified()).expect("time");
    fs::write(&path, b"chr1\t0\t4\tgap\n").expect("corrupted");
    fs::File::options().write(true).open(&path).and_then(|f| f.set_modified(modified)).expect("time kept");
    let verification = setup::verify(&text, &dir, &[Platform::Ont]).expect("verify");
    assert!(!verification.passed(), "{verification}");
    assert!(
        verification.files.iter().any(|(p, c)| p == "blacklist.bed" && matches!(c, Check::Mismatch(_))),
        "{verification}"
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_cut_off_download_resumes_where_it_stopped() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("resume");
    let trf_len = TRF.len() as u64;
    server.fault("trf.bed", Fault::CutOff(7));
    // A partial file from an earlier setup, too.
    fs::create_dir_all(dir.join("models/m_ont")).expect("dir");
    fs::write(dir.join("models/m_ont/pileup.pt.part"), b"ont p").expect("partial");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    let log = server.take_log();
    let of = |path: &str| -> Vec<Option<u64>> { log.iter().filter(|(p, _)| p == path).map(|(_, r)| *r).collect() };
    assert_eq!(of("/trf.bed"), [None, Some(7)], "resumed from byte 7 of {trf_len}");
    assert_eq!(of("/ont/pileup.pt"), [Some(5)], "resumed from the earlier setup's 5 bytes");
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);
    assert_eq!(fs::read(dir.join("models/m_ont/pileup.pt")).expect("fetched"), b"ont pileup");
    assert!(!dir.join("models/m_ont/pileup.pt.part").exists());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_partial_file_that_is_wrong_is_fetched_again_from_the_start() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("bad-part");
    fs::write(dir.join("trf.bed.part"), b"XXXX").expect("partial");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    let trf: Vec<Option<u64>> =
        server.take_log().into_iter().filter(|(p, _)| p == "/trf.bed").map(|(_, r)| r).collect();
    assert_eq!(trf, [Some(4), None], "resumed, found wrong, fetched whole");
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn files_without_a_checksum_are_fetched_checked_by_size_and_never_verified() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &["trf", "m_ont/pileup"], &[]);
    let dir = temp_dir("unpinned");
    // Fetched and used with no flag, as any other file.
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [None]);
    let site_path = dir.join("site.toml");
    let site = Site::load(&site_path).expect("site.toml");
    assert_eq!(site.sv.tandem_repeats, Some(fs::canonicalize(&dir).expect("canonical").join("trf.bed")));
    let site = fs::read_to_string(&site_path).expect("site.toml");
    assert!(!site.contains("--allow-unpinned") && !site.contains("UNVERIFIED"), "{site}");
    let size_only = site
        .split("# Checked by size only")
        .nth(1)
        .and_then(|s| s.split("\n#\n").next())
        .expect("the files checked by size only");
    assert!(size_only.contains("#   trf.bed\n") && size_only.contains("#   models/m_ont/pileup.pt"), "{site}");
    assert!(!size_only.contains("blacklist.bed"), "only those without a checksum: {site}");

    // --verify passes, saying the two were checked by size only, not verified.
    let verification = setup::verify(&text, &dir, &[Platform::Ont]).expect("verify");
    assert!(verification.passed(), "{verification}");
    let check = |v: &setup::Verification, path: &str| v.files.iter().find(|(p, _)| p == path).map(|(_, c)| c.clone());
    assert!(matches!(check(&verification, "trf.bed"), Some(Check::SizeOnly(_))), "{verification}");
    assert!(matches!(check(&verification, "blacklist.bed"), Some(Check::Verified(_))), "{verification}");
    assert!(format!("{verification}").contains("2 checked by size only"), "{verification}");

    // A rerun fetches nothing.
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up again");
    assert_eq!(server.take_log(), Vec::<Request>::new());

    // One of another size fails --verify, and a rerun fetches it again.
    fs::write(dir.join("trf.bed"), b"chr1\t0\t10\n").expect("truncated");
    let verification = setup::verify(&text, &dir, &[Platform::Ont]).expect("verify");
    assert!(!verification.passed() && matches!(check(&verification, "trf.bed"), Some(Check::Mismatch(_))));
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up again");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [None]);
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched again"), TRF);
    assert!(setup::verify(&text, &dir, &[Platform::Ont]).expect("verify").passed());

    // A missing one fails it.
    fs::remove_file(dir.join("models/m_ont/pileup.pt")).expect("removed");
    let verification = setup::verify(&text, &dir, &[Platform::Ont]).expect("verify");
    assert!(!verification.passed() && check(&verification, "models/m_ont/pileup.pt") == Some(Check::Missing));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn terms_are_accepted_before_a_file_is_fetched() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &["m_ont/pileup"]);
    let dir = temp_dir("terms");
    let err = setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect_err("refused");
    let message = format!("{err:#}");
    assert!(message.contains("research use only") && message.contains("--accept-licences"), "{message}");
    assert_eq!(server.take_log(), Vec::<Request>::new(), "nothing fetched");
    // Not needed without the platform whose model it is.
    setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi needs no terms");
    let accepted = Options { accept_licences: true, ..options(&dir, &[Platform::Ont]) };
    setup::setup(&text, &accepted, &fetcher()).expect("accepted");
    let site = fs::read_to_string(dir.join("site.toml")).expect("site.toml");
    assert!(site.contains("m_ont: test licence m_ont/pileup (its terms accepted"), "{site}");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn file_urls_work_and_pin_writes_sha256s_into_the_manifest() {
    let _forks = shared();
    let publisher = temp_dir("pin-publisher");
    for (name, data) in published() {
        let path = publisher.join(name);
        fs::create_dir_all(path.parent().expect("parent")).expect("dir");
        fs::write(path, data).expect("published");
    }
    let url = |n: &str| format!("file://{}", publisher.join(n).display());
    let ids = ["reference", "trf", "m_ont/pileup", "m_ont/full_alignment", "m_hifi/pileup", "m_hifi/full_alignment"];
    let path = publisher.join("manifest.toml");
    fs::write(&path, format!("# a comment that stays\n{}", manifest(&url, &ids, &[]))).expect("manifest");
    let scratch = temp_dir("pin-scratch");
    let pinned = setup::pin(&path, &scratch, false, &fetcher()).expect("pinned");
    let files = published();
    let data = |name: &str| &files.iter().find(|(n, _)| *n == name).expect("published").1;
    assert_eq!(pinned.len(), 6, "every file without a sha256, every platform's");
    assert!(
        pinned.contains(&("reference".to_owned(), sha256(data("ref.fa.gz")))),
        "the download's, not the unpacked file's"
    );
    assert!(pinned.contains(&("trf".to_owned(), sha256(data("trf.bed")))));
    let text = fs::read_to_string(&path).expect("manifest");
    assert!(text.starts_with("# a comment that stays\n") && !text.contains("unpinned"), "{text}");

    // Now every file is verified, fetched from file:// URLs.
    let dir = temp_dir("pin-setup");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert!(setup::verify(&text, &dir, &[Platform::Ont]).expect("verify").passed());
    assert_eq!(fs::read(dir.join("trf.bed")).expect("copied"), TRF);
    for d in [publisher, scratch, dir] {
        let _ = fs::remove_dir_all(d);
    }
}

#[test]
fn a_bgzipped_reference_gets_its_fai_and_gzi() {
    let _forks = exclusive();
    let publisher = temp_dir("bgzf-publisher");
    // Blocks of 16 bytes, as a BGZF writer flushed at each.
    let mut w = pasta::bgzf::io::Writer::new(Vec::new());
    let mut blocks = Vec::new();
    for (i, chunk) in FASTA.chunks(16).enumerate() {
        if i > 0 {
            blocks.push((w.position(), (i * 16) as u64));
        }
        w.write_all(chunk).expect("written");
        w.flush().expect("flushed");
    }
    let bgzf = w.finish().expect("finished");
    fs::write(publisher.join("ref.fa.gz"), &bgzf).expect("published");
    let text = format!(
        "format = 1\ngenome = \"tiny\"\n[site]\nreference = \"reference\"\n[[file]]\nid = \"reference\"\nurl = \"file://{}\"\n\
         path = \"tiny.fa.gz\"\nsha256 = \"{}\"\nlicence = \"x\"\n",
        publisher.join("ref.fa.gz").display(),
        sha256(&bgzf)
    );
    let dir = temp_dir("bgzf");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(fs::read_to_string(dir.join("tiny.fa.gz.fai")).expect(".fai"), FAI, "offsets in the uncompressed FASTA");
    let gzi = pasta::bgzf::gzi::fs::read(dir.join("tiny.fa.gz.gzi")).expect(".gzi");
    assert_eq!(gzi.as_ref(), blocks.as_slice());
    let site = Site::load(dir.join("site.toml")).expect("site.toml");
    assert!(!offered(&dir, "ont", "m_ont") && site.repeats.is_none() && site.cnv.blacklist.is_none());
    assert!(setup::verify(&text, &dir, &[Platform::Ont]).expect("verify").passed());
    for d in [publisher, dir] {
        let _ = fs::remove_dir_all(d);
    }
}

fn lrpipe(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_lrpipe")).args(args).output().expect("lrpipe runs")
}

#[test]
fn offline_prints_the_urls_checksums_and_commands_and_fetches_nothing() {
    let _forks = exclusive();
    let dir = temp_dir("offline");
    let target = dir.join("ref");
    let out =
        lrpipe(&["setup", "--genome", "grch38", "--platform", "ont", "--dir", &target.to_string_lossy(), "--offline"]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8_lossy(&out.stdout);
    let reference = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/\
                     seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz";
    let dest = target.join("sources/GRCh38_no_alt_analysis_set.fa.gz");
    assert!(text.contains(&format!("curl -fL -C - -o {} {reference}", dest.display())), "{text}");
    assert!(text.contains(&format!("wget -c -O {} {reference}", dest.display())), "{text}");
    assert!(text.contains("md5 a08035b6a6e31780e96a34008ff21bd6"), "{text}");
    assert!(text.contains("r1041_e82_400bps_sup_v500/pileup.pt") && !text.contains("hifi_revio"), "{text}");
    assert!(!text.contains("Research Purposes only"), "ONT's terms need no accepting: {text}");
    assert!(text.contains("--accept-licences") && !text.contains("--allow-unpinned"), "{text}");
    // A model without a checksum is listed with its size and why it has none.
    let model = text.lines().find(|l| l.starts_with("# r1041_e82_400bps_sup_v500/pileup:")).expect("the model");
    assert!(model.contains("8304746 bytes") && model.contains("no checksum: "), "{model}");
    assert!(!model.contains("must be accepted"), "ONT's models are fetched without asking: {model}");
    let trf = text.lines().find(|l| l.starts_with("# tandem_repeats:")).expect("the TRF BED");
    assert!(trf.contains("its terms must be accepted"), "{trf}");
    assert!(text.contains("PacBio Software License Agreement"), "{text}");
    assert!(!target.exists(), "nothing is fetched or made");
    let setup = text.lines().last().expect("the setup to run");
    assert!(setup.starts_with("# lrpipe setup --genome GRCh38 --platform ont --dir "), "{setup}");

    // A manifest of one's own, and a bonefish: the setup to run names them, not a genome.
    let own = dir.join("own.toml");
    fs::write(
        &own,
        lrpipe::setup::manifest::Genome::Grch38.manifest_text().replace("genome = \"GRCh38\"", "genome = \"mine\""),
    )
    .expect("written");
    let out = lrpipe(&[
        "setup",
        "--manifest",
        &own.to_string_lossy(),
        "--platform",
        "ont",
        "--dir",
        &target.to_string_lossy(),
        "--offline",
        "--bonefish",
        "/opt/bonefish",
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8_lossy(&out.stdout);
    let setup = text.lines().last().expect("the setup to run");
    assert!(
        setup.starts_with(&format!("# lrpipe setup --manifest {} --platform ont", own.display()))
            && setup.ends_with(" --bonefish /opt/bonefish")
            && !setup.contains("--genome"),
        "{setup}"
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn the_verify_command_exits_with_1_until_the_setup_is_done() {
    let _forks = exclusive();
    let dir = temp_dir("cli");
    let target = dir.join("ref");
    let target = target.to_string_lossy();
    // The flag that once let files without a checksum through is gone.
    let out = lrpipe(&["setup", "--genome", "GRCh38", "--platform", "hifi", "--dir", &target, "--allow-unpinned"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(!Path::new(&*target).exists(), "nothing fetched or made");

    let server = serve_tiny();
    let manifest_path = dir.join("tiny.toml");
    fs::write(&manifest_path, manifest(&|n| server.url(n), &[], &[])).expect("manifest");
    let m = manifest_path.to_string_lossy();
    let out = lrpipe(&["setup", "--manifest", &m, "--platform", "ont", "--dir", &target, "--verify"]);
    assert_eq!(out.status.code(), Some(1), "nothing there yet");
    assert!(String::from_utf8_lossy(&out.stdout).contains("MISSING"));
    let out = lrpipe(&["setup", "--manifest", &m, "--platform", "ont", "--dir", &target]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    // --verify of a directory set up before needs neither --genome nor --platform.
    let out = lrpipe(&["setup", "--dir", &target, "--verify"]);
    assert!(out.status.success(), "{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    let _ = fs::remove_dir_all(dir);
}

/// `text` with the `size` of file `id` left out.
fn without_size(text: &str, id: &str) -> String {
    let start = text.find(&format!("id = \"{id}\"")).expect("the file");
    let size = start + text[start..].find("size = ").expect("a size");
    let end = size + text[size..].find('\n').expect("a line") + 1;
    format!("{}{}", &text[..size], &text[end..])
}

fn md5(data: &[u8]) -> String {
    use md5::Digest as _;
    lrpipe::provenance::hex(&md5::Md5::digest(data))
}

/// The requests for `name` in `log`, by the first byte asked for.
fn ranges(log: &[Request], name: &str) -> Vec<Option<u64>> {
    log.iter().filter(|(p, _)| *p == format!("/{name}")).map(|(_, r)| *r).collect()
}

#[test]
fn the_resume_budget_is_per_stall_not_in_all() {
    let _forks = shared();
    // Six cut-offs, more than RESUMES, each bringing bytes: the transfer still completes.
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("stalls");
    for _ in 0..6 {
        server.fault("trf.bed", Fault::CutOff(3));
    }
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    let trf = ranges(&server.take_log(), "trf.bed");
    assert_eq!(trf, [None, Some(3), Some(6), Some(9), Some(12), Some(15), Some(18)]);
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_server_that_sends_the_whole_file_on_a_resume_overwrites_the_part() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("whole");
    // The right size of the wrong bytes would pass a size check: only a 200 written from the
    // start (not appended) gives the file.
    fs::write(dir.join("trf.bed.part"), b"XXXXXXXX").expect("partial");
    server.fault("trf.bed", Fault::IgnoreRange);
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [Some(8)], "one request");
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_wrong_content_range_starts_the_file_again() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("wrong-range");
    fs::write(dir.join("trf.bed.part"), &TRF[..5]).expect("partial");
    server.fault("trf.bed", Fault::WrongRange);
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [Some(5), None], "refused, then fetched from the start");
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_part_that_is_already_whole_gets_416_and_is_checked() {
    let _forks = shared();
    let server = serve_tiny();
    // Without a size in the manifest, only the server can say the part is whole.
    let text = without_size(&manifest(&|n| server.url(n), &[], &[]), "trf");
    let dir = temp_dir("416");
    fs::write(dir.join("trf.bed.part"), TRF).expect("whole");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [Some(TRF.len() as u64)], "416, and nothing fetched again");
    assert_eq!(fs::read(dir.join("trf.bed")).expect("there"), TRF);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn an_interrupted_gunzip_is_done_again_from_the_download() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("gunzip");
    // A setup killed while unpacking: the download is there, the FASTA half written.
    fs::create_dir_all(dir.join("sources")).expect("dir");
    fs::write(dir.join("sources/tiny.fa.gz"), gzip(FASTA)).expect("download");
    fs::write(dir.join("tiny.fa.part"), &FASTA[..10]).expect("half");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "ref.fa.gz"), Vec::<Option<u64>>::new(), "not fetched again");
    assert_eq!(fs::read(dir.join("tiny.fa")).expect("unpacked"), FASTA);
    assert!(!dir.join("tiny.fa.part").exists() && !dir.join("sources/tiny.fa.gz").exists());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_resume_of_a_file_without_a_checksum_asks_for_the_same_version() {
    let _forks = shared();
    let server = serve_tiny();
    server.publish("trf.bed", TRF, "v1");
    let text = manifest(&|n| server.url(n), &["trf"], &[]);
    let dir = temp_dir("if-range");
    server.fault("trf.bed", Fault::CutOff(7));
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [None, Some(7)]);
    assert_eq!(server.if_ranges("trf.bed"), [None, Some("\"v1\"".to_owned())]);
    assert_eq!(fs::read(dir.join("trf.bed")).expect("fetched"), TRF);

    // The file changed on the server while a part of the old one waited: the server sends the
    // whole new one, which replaces the part rather than finishing it.
    let new = b"chr3\t0\t10\nchr4\t0\t500\n";
    server.publish("trf.bed", new, "v2");
    let text = manifest(&|n| server.url(n), &["trf"], &[]).replace(&format!("size = {}\n", TRF.len()), "");
    let dir2 = temp_dir("if-range-changed");
    fs::write(dir2.join("trf.bed.part"), &TRF[..7]).expect("old part");
    fs::write(dir2.join("trf.bed.part.validator"), "\"v1\"").expect("its version");
    setup::setup(&text, &options(&dir2, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(fs::read(dir2.join("trf.bed")).expect("fetched"), new);

    // A part of an unknown version of a file without a checksum is not resumed at all.
    let dir3 = temp_dir("if-range-unknown");
    fs::write(dir3.join("trf.bed.part"), b"chr3\t0").expect("a part");
    server.take_log();
    setup::setup(&text, &options(&dir3, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(ranges(&server.take_log(), "trf.bed"), [None]);
    for d in [dir, dir2, dir3] {
        let _ = fs::remove_dir_all(d);
    }
}

#[test]
fn pinning_keeps_a_verified_setup_and_verifies_one_checked_by_size_only() {
    let _forks = shared();
    let server = serve_tiny();
    let files = published();
    let gz = &files.iter().find(|(n, _)| *n == "ref.fa.gz").expect("published").1;
    // The reference with the publisher's md5 only, as GRCh38's and CHM13's are; trf unpinned.
    let text = manifest(&|n| server.url(n), &["trf"], &[])
        .replace(&format!("sha256 = \"{}\"", sha256(gz)), &format!("md5 = \"{}\"", md5(gz)));
    let dir = temp_dir("pin-keeps");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    server.take_log();
    let trf = |v: &setup::Verification| v.files.iter().find(|(p, _)| p == "trf.bed").map(|(_, c)| c.clone());
    let before = setup::verify(&text, &dir, &[Platform::Ont]).expect("verify");
    assert!(matches!(trf(&before), Some(Check::SizeOnly(_))), "{before}");

    // `lrpipe setup --pin` adds the reference's sha256 and trf's.
    let pinned =
        lrpipe::setup::manifest::pin_text(&text, &[("reference".into(), sha256(gz)), ("trf".into(), sha256(TRF))])
            .expect("pinned");
    setup::setup(&pinned, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    assert_eq!(server.take_log(), Vec::<Request>::new(), "nothing fetched again");
    let verification = setup::verify(&pinned, &dir, &[Platform::Ont]).expect("verify");
    assert!(verification.passed(), "{verification}");
    assert!(matches!(trf(&verification), Some(Check::Verified(_))), "{verification}");
    // A checksum the fetched file does not have is noticed, without reading it again.
    let wrong = lrpipe::setup::manifest::pin_text(&pinned, &[("trf".into(), "0".repeat(64))]).expect("pinned");
    setup::setup(&wrong, &options(&dir, &[Platform::Ont]), &fetcher()).expect_err("trf does not match");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn the_verify_command_site_toml_suggests_passes_and_fails_as_it_should() {
    let _forks = exclusive();
    let server = serve_tiny();
    let dir = temp_dir("verify-cmd");
    let target = dir.join("ref");
    let manifest_path = dir.join("tiny.toml");
    fs::write(&manifest_path, manifest(&|n| server.url(n), &["trf"], &[])).expect("manifest");
    let m = manifest_path.to_string_lossy();
    let t = target.to_string_lossy();
    let out = lrpipe(&["setup", "--manifest", &m, "--platform", "ont", "--dir", &t]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let site = fs::read_to_string(target.join("site.toml")).expect("site.toml");
    let command = site.lines().find_map(|l| l.strip_prefix("#   lrpipe ")).expect("the verify command");
    let words: Vec<&str> = command.split_whitespace().collect();
    assert_eq!(words, ["setup", "--verify", "--dir", &*fs::canonicalize(&target).expect("there").to_string_lossy()]);
    let out = lrpipe(&words);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "{stdout}");
    assert!(stdout.contains("size only  trf.bed") && stdout.contains("1 checked by size only:"), "{stdout}");
    // A file without a checksum of another size fails it; so does a wrong or missing one.
    fs::write(target.join("trf.bed"), b"chr1\t0\t10\n").expect("truncated");
    let out = lrpipe(&words);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stdout).contains("MISMATCH   trf.bed"));
    fs::write(target.join("trf.bed"), TRF).expect("restored");
    assert!(lrpipe(&words).status.success());
    fs::write(target.join("blacklist.bed"), b"wrong").expect("corrupted");
    assert_eq!(lrpipe(&words).status.code(), Some(1));
    fs::remove_file(target.join("blacklist.bed")).expect("removed");
    assert_eq!(lrpipe(&words).status.code(), Some(1));
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_catalogue_that_does_not_match_is_not_put_in_place() {
    let _forks = shared();
    let server = serve_tiny();
    let good = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("catalogue");
    setup::setup(&good, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up");
    // A manifest whose catalogue checksum the inputs do not make.
    let bad = good.replace(&sha256(CATALOG.as_bytes()), &"0".repeat(64));
    let err = setup::setup(&bad, &options(&dir, &[Platform::Ont]), &fetcher()).expect_err("refused");
    assert!(format!("{err:#}").contains("not the manifest's"), "{err:#}");
    assert_eq!(fs::read_to_string(dir.join("catalog.bed")).expect("kept"), CATALOG, "the good one stays");
    assert!(!dir.join("catalog.bed.part").exists());
    let fresh = temp_dir("catalogue-fresh");
    setup::setup(&bad, &options(&fresh, &[Platform::Ont]), &fetcher()).expect_err("refused");
    assert!(!fresh.join("catalog.bed").exists(), "a wrong catalogue is never there");
    for d in [dir, fresh] {
        let _ = fs::remove_dir_all(d);
    }
}

#[cfg(unix)]
#[test]
fn an_unpacked_file_is_recorded_before_its_download_goes() {
    let _forks = shared();
    use std::os::unix::fs::PermissionsExt as _;
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("record-first");
    let sources = dir.join("sources");
    fs::create_dir_all(&sources).expect("dir");
    fs::write(sources.join("tiny.fa.gz"), gzip(FASTA)).expect("download");
    fs::write(sources.join("adotto.bed.gz"), gzip(ADOTTO)).expect("download");
    fs::write(sources.join("strchive.bed"), STRCHIVE).expect("download");
    // The download cannot be removed: setup goes on, the unpacked file recorded. (Root, as in
    // the release build's container, removes it anyway; the rerun below holds either way.)
    fs::set_permissions(&sources, fs::Permissions::from_mode(0o555)).expect("read-only");
    let first = setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher());
    fs::set_permissions(&sources, fs::Permissions::from_mode(0o755)).expect("writable");
    first.expect("set up");
    let modified = fs::metadata(dir.join("tiny.fa")).and_then(|m| m.modified()).expect("time");
    match fs::remove_file(sources.join("tiny.fa.gz")) {
        Err(e) if e.kind() != std::io::ErrorKind::NotFound => panic!("removing the download by hand: {e}"),
        _ => {}
    }
    server.take_log();
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up again");
    assert_eq!(server.take_log(), Vec::<Request>::new(), "the recorded FASTA needs no download");
    assert_eq!(
        fs::metadata(dir.join("tiny.fa")).and_then(|m| m.modified()).expect("time"),
        modified,
        "not unpacked again"
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn one_setup_at_a_time() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("lock");
    fs::create_dir_all(dir.join(".lrpipe-setup")).expect("dir");
    let held = fs::File::create(dir.join(".lrpipe-setup/lock")).expect("lock file");
    held.try_lock().expect("locked");
    let err = setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect_err("refused");
    assert!(format!("{err:#}").contains("another lrpipe setup"), "{err:#}");
    assert_eq!(server.take_log(), Vec::<Request>::new());
    drop(held);
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("set up once it is free");
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_transfer_reports_its_progress() {
    let _forks = shared();
    let publisher = temp_dir("progress-publisher");
    let data = vec![b'A'; 300_000];
    fs::write(publisher.join("big"), &data).expect("published");
    let dir = temp_dir("progress");
    let lines = std::cell::RefCell::new(Vec::new());
    let say = |line: &str| lines.borrow_mut().push(line.to_owned());
    let expected = lrpipe::setup::fetch::Expected { size: Some(300_000), ..Default::default() };
    Fetcher::new()
        .with_progress(Duration::ZERO)
        .fetch(&format!("file://{}", publisher.join("big").display()), &dir.join("big"), &expected, &say)
        .expect("fetched");
    let lines = lines.into_inner();
    assert!(lines.iter().any(|l| l == "big: 300.0 kB of 300.0 kB (100%)"), "{lines:?}");
    assert!(lines.len() > 1, "more than once: {lines:?}");
    for d in [publisher, dir] {
        let _ = fs::remove_dir_all(d);
    }
}

#[test]
fn a_rerun_for_another_platform_keeps_the_first_ones_models() {
    let _forks = shared();
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("union");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("ONT");
    setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi");
    assert!(offered(&dir, "ont", "m_ont") && offered(&dir, "hifi", "m_hifi"), "both platforms' models");
    let (_, platforms) = setup::recorded_setup(&dir).expect("recorded");
    assert_eq!(platforms, [Platform::Hifi, Platform::Ont]);
    // One whose models are gone is left out.
    fs::remove_file(dir.join("models/m_ont/pileup.pt")).expect("removed");
    server.take_log();
    setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi again");
    assert!(!offered(&dir, "ont", "m_ont"));
    assert_eq!(server.take_log(), Vec::<Request>::new(), "nothing of the left-out platform is fetched");
    let _ = fs::remove_dir_all(dir);
}

/// Whether DIR/site.toml offers `model` as the commented-out `PLATFORM_model` override.
fn offered(dir: &Path, platform: &str, model: &str) -> bool {
    let text = fs::read_to_string(dir.join("site.toml")).expect("site.toml");
    let line = format!("# {platform}_model = \"{}\"", dir.join("models").join(model).display());
    text.lines().any(|l| l == line)
}

#[test]
fn a_rerun_keeps_a_site_toml_edited_by_hand() {
    // LUX-18: setup's own site.toml is replaced by a rerun's; one edited since is kept, the
    // rerun's going beside it.
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("edited-site");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("ONT");
    setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi");
    assert!(!dir.join("site.toml.new").exists(), "setup's own is replaced");
    let site = dir.join("site.toml");
    let edited = format!("{}\n# my GPU\n", fs::read_to_string(&site).expect("site.toml"));
    fs::write(&site, &edited).expect("edited");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("ONT again");
    assert_eq!(fs::read_to_string(&site).expect("site.toml"), edited, "kept");
    let new = fs::read_to_string(dir.join("site.toml.new")).expect("the rerun's beside it");
    assert!(!new.contains("# my GPU"), "{new}");
    // Moved into place, it is setup's own again.
    fs::rename(dir.join("site.toml.new"), &site).expect("moved");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("ONT once more");
    assert!(!dir.join("site.toml.new").exists());
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn a_site_toml_from_a_setup_that_recorded_no_checksum_is_still_setups_own() {
    // A setup from before lrpipe recorded site.toml's SHA-256: its untouched site.toml is the
    // one its recorded manifest and platforms make, and a rerun replaces it.
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("unrecorded-site");
    setup::setup(&text, &options(&dir, &[Platform::Ont]), &fetcher()).expect("ONT");
    let state = dir.join(".lrpipe-setup/state.json");
    let mut json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&state).expect("state")).expect("json");
    assert!(json.as_object_mut().expect("an object").remove("site_toml").is_some(), "recorded now");
    fs::write(&state, json.to_string()).expect("written");
    setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi");
    assert!(!dir.join("site.toml.new").exists(), "replaced, not kept");
    assert!(offered(&dir, "hifi", "m_hifi"));
    // Edited, it is kept, recorded or not.
    let site = dir.join("site.toml");
    fs::write(&site, format!("{}# mine\n", fs::read_to_string(&site).expect("site.toml"))).expect("edited");
    let mut json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&state).expect("state")).expect("json");
    json.as_object_mut().expect("an object").remove("site_toml");
    fs::write(&state, json.to_string()).expect("written");
    let written = setup::setup(&text, &options(&dir, &[Platform::Hifi]), &fetcher()).expect("HiFi again");
    assert!(dir.join("site.toml.new").exists(), "kept, the rerun's beside it");
    assert_eq!(written.file_name().and_then(|n| n.to_str()), Some("site.toml.new"), "the file it wrote");
    let _ = fs::remove_dir_all(dir);
}

/// A stand-in for bonefish that writes, for `-d INDEX --idx-raw REFERENCE`, the reference into
/// INDEX, and counts its builds in `builds` beside it.
const FAKE_BONEFISH: &str = r#"#!/bin/sh
if [ "$1" = "--version" ]; then echo "bonefish 0.4.0"; exit 0; fi
[ "$1" = "-d" ] && [ "$3" = "--idx-raw" ] || exit 9
echo built >> "$(dirname "$0")/builds"
cp "$4" "$2"
"#;

#[test]
fn bonefish_builds_its_index_once_and_site_toml_names_both() {
    let _forks = exclusive();
    use std::os::unix::fs::PermissionsExt;
    // F6: `lrpipe setup --bonefish PATH` builds bonefish's index of the reference and writes
    // [tools] bonefish and [align] index; a rerun, with or without --bonefish, keeps it, and a new
    // bonefish version builds it again.
    let server = serve_tiny();
    let text = manifest(&|n| server.url(n), &[], &[]);
    let dir = temp_dir("bonefish");
    let tools = temp_dir("bonefish-tools");
    let bonefish = tools.join("bonefish");
    fs::write(&bonefish, FAKE_BONEFISH).expect("written");
    fs::set_permissions(&bonefish, fs::Permissions::from_mode(0o755)).expect("executable");
    let with = Options { bonefish: Some(bonefish.clone()), ..options(&dir, &[Platform::Ont]) };
    let site_path = setup::setup(&text, &with, &fetcher()).expect("set up");
    let dir = fs::canonicalize(&dir).expect("canonical");
    let site = Site::load(&site_path).expect("a site configuration");
    assert_eq!(site.tools.bonefish.as_deref(), Some(bonefish.as_path()));
    let index = dir.join("tiny.fa.bonefish.raw");
    assert_eq!(site.align.index.as_deref(), Some(index.as_path()));
    assert_eq!(fs::read(&index).expect("built"), FASTA, "built from the reference");
    assert!(!dir.join("tiny.fa.bonefish.raw.part").exists());
    let builds = || fs::read_to_string(tools.join("builds")).unwrap_or_default().lines().count();
    assert_eq!(builds(), 1);
    let verification = setup::verify(&text, &dir, &[Platform::Ont]).expect("verified");
    assert!(
        verification.files.iter().any(|(p, c)| p == "tiny.fa.bonefish.raw" && matches!(c, Check::Verified(_))),
        "{verification}"
    );
    assert!(verification.passed(), "{verification}");

    // Again, and without --bonefish: kept, not built again.
    let before = fs::read_to_string(&site_path).expect("site.toml");
    for again in [with.clone(), options(&dir, &[Platform::Ont])] {
        setup::setup(&text, &again, &fetcher()).expect("set up again");
        assert_eq!(fs::read_to_string(&site_path).expect("site.toml"), before);
    }
    assert_eq!(builds(), 1, "the index is built once");

    // A new bonefish version builds it again; a failing one leaves no index behind.
    fs::write(&bonefish, FAKE_BONEFISH.replace("bonefish 0.4.0", "bonefish 0.5.0")).expect("written");
    setup::setup(&text, &with, &fetcher()).expect("set up again");
    assert_eq!(builds(), 2, "another version's index format may differ");
    fs::write(
        &bonefish,
        FAKE_BONEFISH.replace("bonefish 0.4.0", "bonefish 0.6.0").replace("cp \"$4\" \"$2\"", "echo no >&2; exit 3"),
    )
    .expect("written");
    let e = format!("{:#}", setup::setup(&text, &with, &fetcher()).expect_err("bonefish fails"));
    assert!(e.contains("-d failed") && e.contains("no"), "{e}");
    assert!(!dir.join("tiny.fa.bonefish.raw.part").exists());

    // Another program given as bonefish: refused by its --version, before it builds anything (it
    // ran, wrote nothing, and setup failed renaming the index it had not written).
    let other = lrpipe::provenance::resolve(Path::new("true")).expect("the true executable");
    let not = Options { bonefish: Some(other.clone()), ..options(&dir, &[Platform::Ont]) };
    let e = format!("{:#}", setup::setup(&text, &not, &fetcher()).expect_err("not bonefish"));
    assert!(e.contains(&format!("{} is not bonefish: its --version printed", other.display())), "{e}");
    let out = lrpipe(&[
        "setup",
        "--genome",
        "grch38",
        "--platform",
        "ont",
        "--dir",
        &dir.to_string_lossy(),
        "--bonefish",
        &other.to_string_lossy(),
    ]);
    assert_eq!(out.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("is not bonefish"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = fs::remove_dir_all(dir);
    let _ = fs::remove_dir_all(tools);
}
