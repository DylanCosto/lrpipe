# Validating lrpipe-sv on a Linux VM or a Mac

`tools/validate.sh` checks that a machine builds lrpipe-sv and produces exactly Sniffles 2.8.1's
output. It needs only Rust and the repository: the reference output for the bundled test data
is checked in. It writes `validation/validation-HOST-ARCH.txt`, the file to send back, and exits
non-zero if anything fails.

What it does:

1. records the platform: CPU, OS, C library, memory, Rust version, commit;
2. builds the release binary (`cargo build --locked --release`);
3. runs the unit tests and `tests/parity.rs` (`cargo test --locked --release`), which compares
   lrpipe-sv with Sniffles 2.8.1 on `tests/data` for default options at 1, 2 and 4 threads,
   `--reference`, CRAM, `--phase`, `--mosaic`, `--no-qc`, `--symbolic`, `--regions`,
   `--genotype-vcf`, multi-sample calling and bgzipped output;
4. times five runs on the fixture;
5. optionally, compares a full-size run on your own data (see below).

A first run takes a few minutes, mostly compiling. The tests need well under 1 GB of memory and
about 2 GB of disk for the build.

## 1. Get the code onto the machine

The bundled test data (3.8 MB) is part of the git repository, so a clone brings everything.
The GitHub repository is `lrpipe`; this program is `lrpipe-sv` (also `lrpipe sv`).

```bash
git clone https://github.com/DylanCosto/lrpipe.git lrpipe
```

## 2. Install Rust and a linker

lrpipe-sv has no C dependencies, but Rust uses the system linker.

### Linux VM (x86-64 or ARM)

Debian or Ubuntu:

```bash
sudo apt-get update && sudo apt-get install -y build-essential git curl
```

Fedora, RHEL or Rocky: `sudo dnf install -y gcc git curl`.

Then Rust (1.91 or newer):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
rustc -V
```

### Mac (Apple Silicon or Intel)

The Xcode command line tools provide the linker and git:

```bash
xcode-select --install
```

Then Rust, as above:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
rustc -V
```

The script runs with the bash that ships with macOS (3.2); nothing else needs installing.

An Apple Silicon Mac can also cover ARM Linux: run the Linux steps in an arm64 Linux VM (UTM,
Parallels, OrbStack or Lima).

## 3. Run the validation

```bash
cd lrpipe/lrpipe-sv
tools/validate.sh
```

The last lines read `VALIDATION PASSED on aarch64` (or `x86_64`, `arm64`) and name the report.
`THREADS=N` sets the thread count for the timing runs (default: all).

## 4. Optional: a full-size comparison

Point the script at an alignment file and either a Sniffles 2.8.1 VCF of it made with the same
options, or a Sniffles 2.8.1 install to run on the machine:

```bash
INPUT=sample.bam ORACLE=sniffles2.vcf tools/validate.sh
INPUT=sample.cram REFERENCE=ref.fa SNIFFLES=$(command -v sniffles) tools/validate.sh
OPTIONS="--minsupport auto" INPUT=sample.bam ORACLE=sniffles2.vcf tools/validate.sh
```

`REFERENCE` is passed to lrpipe-sv (and to Sniffles2 with `SNIFFLES`); `OPTIONS` adds options to
both. With `ORACLE`, the VCF must have been made with the same reference and options.

A ready-made pair from the build machine: HG003 chr20 (2.5 GB) and the Sniffles 2.8.1 VCF for it
without a reference.

```bash
scp USER@BUILD_HOST:/data2/bioinformatics/lrpipe-sv-bench/hg003.chr20.bam* .
scp USER@BUILD_HOST:/data2/bioinformatics/lrpipe-sv-bench/parity/noref/oracle.vcf hg003.chr20.sniffles2.vcf
INPUT=hg003.chr20.bam ORACLE=hg003.chr20.sniffles2.vcf tools/validate.sh
```

To run Sniffles 2.8.1 on the machine itself, install it with conda (`conda create -n sniffles -c
bioconda -c conda-forge sniffles=2.8.1`) or pip (`pip install sniffles==2.8.1`).

## If something fails

The report holds the failing command's output. For a parity test failure, the message shows the
first differing VCF line, expected and actual.

The one platform-dependent part of the computation is in genotype likelihoods, which call the C
library's `pow` and `log`. Sniffles2 does the same through Python. A last-bit difference between
C libraries (glibc on x86-64 versus ARM, or macOS's libm) could change a rounded QUAL or GQ value.
To tell that apart from a bug, run the full-size comparison with `SNIFFLES=` so that the
reference output is made on the same machine: if Sniffles2 and lrpipe-sv agree there, the
difference is the platform's, not lrpipe-sv's.
