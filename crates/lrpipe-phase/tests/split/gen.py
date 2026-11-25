"""Random inputs for comparing `whatshap split` with `lrpipe-phase split`.

Usage: python gen.py SEED DIR

Writes reads (DIR/reads.fastq[.gz] or DIR/reads.bam), a haplotype list
(DIR/list.txt, possibly compressed) and DIR/cases.txt, one set of extra
command-line arguments per line; the words H1, H2, ... UNTAGGED and HIST are
replaced by per-tool output files (with the extension given after a colon,
e.g. H1:.fastq.gz). The comparison script adds the reads and the list.

FASTQ inputs have comments, multi-line and FASTA records, CRLF line endings
and empty sequences; BAM inputs are written byte by byte (header text with
NUL padding or without a final newline, records without sequence but with a
CIGAR, auxiliary tags) and BGZF-compressed by the comparison script.
"""

import bz2
import gzip
import lzma
import random
import struct
import sys

BASES = "ACGT"


def read_names(rng, n):
    names = []
    for i in range(n):
        name = rng.choice([f"read{i}", f"m64011_190830/{i}/ccs", f"r{i}_x"])
        names.append(name)
    if n > 3 and rng.random() < 0.1:
        names[-1] = names[0]  # duplicate read name
    return names


def write_fastq(rng, path, names):
    crlf = rng.random() < 0.1
    nl = "\r\n" if crlf else "\n"
    parts = []
    for name in names:
        length = rng.choice([0, 1, 5, 50, 300]) if rng.random() < 0.05 else rng.randint(20, 400)
        seq = "".join(rng.choice(BASES) for _ in range(length))
        qual = "".join(chr(33 + rng.randint(0, 40)) for _ in range(length))
        header = name + rng.choice(["", "", " comment here", "\tx=1", " "])
        if rng.random() < 0.05:
            parts.append(">" + header + nl + seq[: length // 2] + nl + seq[length // 2 :] + nl)
        elif rng.random() < 0.1 and length > 2:
            parts.append(
                "@" + header + nl + seq[:2] + nl + seq[2:] + nl + "+" + rng.choice(["", name]) + nl
                + qual[:2] + nl + qual[2:] + nl
            )
        else:
            parts.append("@" + header + nl + seq + nl + "+" + nl + qual + nl)
    data = "".join(parts).encode()
    if path.endswith(".gz"):
        with gzip.open(path, "wb") as f:
            f.write(data)
    else:
        with open(path, "wb") as f:
            f.write(data)


def bam_record(rng, name, refs):
    length = 0 if rng.random() < 0.1 else rng.randint(1, 300)
    seq = [rng.randrange(16) for _ in range(length)]
    mapped = refs and rng.random() < 0.4
    cigar = []
    if mapped and length > 0:
        # the query length of the CIGAR must match the sequence
        clip = rng.randint(0, length - 1)
        cigar = ([(clip, 4)] if clip else []) + [(length - clip, 0)] + ([(5, 5)] if rng.random() < 0.3 else [])
    elif length == 0 and rng.random() < 0.7:
        n = rng.choice([1, 2, 3])
        cigar = [(rng.randint(1, 100), rng.choice([0, 1, 2, 4, 5, 7, 8])) for _ in range(n)]
    tid = rng.randrange(len(refs)) if mapped else -1
    pos = rng.randint(0, 1000) if mapped else -1
    flag = 0 if mapped else 4
    qname = name.encode() + b"\0"
    tags = b""
    if rng.random() < 0.5:
        tags += b"RGZgroup1\0"
    if rng.random() < 0.3:
        tags += b"npi" + struct.pack("<i", rng.randint(0, 99))
    data = qname
    data += b"".join(struct.pack("<I", n << 4 | op) for n, op in cigar)
    data += bytes((seq[i] << 4) | (seq[i + 1] if i + 1 < length else 0) for i in range(0, length, 2))
    data += bytes(rng.randint(0, 40) for _ in range(length)) if rng.random() < 0.9 else b"\xff" * length
    data += tags
    core = struct.pack(
        "<iiBBHHHiiii", tid, pos, len(qname), rng.randint(0, 60), 4680, len(cigar), flag, length, -1, -1, 0
    )
    block = core + data
    return struct.pack("<i", len(block)) + block


def write_bam(rng, path, names):
    refs = [] if rng.random() < 0.5 else [("chr1", 5000), ("chr2", 3000)]
    lines = ["@HD\tVN:1.6\tSO:unknown"]
    lines += [f"@SQ\tSN:{n}\tLN:{l}" for n, l in refs]
    lines.append("@RG\tID:group1\tSM:sample")
    text = "\n".join(lines).encode()
    r = rng.random()
    if r < 0.7:
        text += b"\n"
    elif r < 0.85:
        text += b"\n" + b"\0" * rng.randint(1, 20)
    # else: no final newline
    out = b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(refs))
    for n, l in refs:
        out += struct.pack("<i", len(n) + 1) + n.encode() + b"\0" + struct.pack("<i", l)
    for name in names:
        out += bam_record(rng, name, refs)
    with open(path, "wb") as f:
        f.write(out)


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(seed)
    ploidy = rng.choice([2, 2, 2, 3, 4])
    n = rng.choice([1, 5, 30, 200])
    names = read_names(rng, n)
    is_bam = rng.random() < 0.45
    if is_bam:
        reads = "reads.bam"
        write_bam(rng, f"{out}/reads.raw", names)
    else:
        reads = rng.choice(["reads.fastq", "reads.fastq", "reads.fastq.gz"])
        write_fastq(rng, f"{out}/{reads}", names)

    # haplotype list: a subset of the reads, and some names not among them
    four = rng.random() < 0.6
    rows = []
    listed = rng.sample(names, rng.randint(0, len(names)))
    if rng.random() < 0.3:
        listed.append("not_in_reads")
    bad = rng.randrange(len(listed)) if listed and rng.random() < 0.05 else -1
    for k, name in enumerate(listed):
        hap = rng.choice(["none"] + [f"H{i}" for i in range(1, ploidy + 1)] * 3)
        if k == bad:
            hap = f"H{ploidy + 1}"
        if four:
            rows.append(f"{name}\t{hap}\t{rng.choice(['1000', '2000', '3000', 'none'])}\t{rng.choice(['chr1', 'chr2'])}")
        else:
            rows.append(f"{name}\t{hap}" + rng.choice(["", "", "\textra"]))
    columns = "#readname\thaplotype\tphaseset\tchromosome" if four else "#readname\thaplotype"
    if rng.random() < 0.05:
        columns = "#readname\thaplotype\tphaseset\tchromosome"
    header = columns + "\n" if rng.random() < 0.5 else ""
    nl = "\r\n" if rng.random() < 0.1 else "\n"
    text = (header + nl.join(rows) + (nl if rows and rng.random() < 0.9 else "")).encode()
    compression = rng.choice(["", "", "", ".gz", ".bz2", ".xz"])
    path = f"{out}/list.txt{compression}"
    opener = {"": open, ".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}[compression]
    with opener(path, "wb") as f:
        f.write(text)
    with open(f"{out}/files", "w") as f:
        f.write(f"{reads} list.txt{compression}\n")

    ext = ".bam" if is_bam else rng.choice([".fastq", ".fastq.gz", ".fastq.bz2", ".fastq.xz"])
    hs = " ".join(f"-o H{i}:{ext}" for i in range(1, ploidy + 1))
    cases = [
        hs,
        hs + f" --output-untagged UNTAGGED:{ext} --read-lengths-histogram HIST:.tsv",
        hs + " --add-untagged --discard-unknown-reads",
        hs + f" --only-largest-block --output-untagged UNTAGGED:{ext} --read-lengths-histogram HIST:.tsv.gz",
    ]
    if ploidy == 2:
        cases.append(f"--output-h1 H1:{ext} --output-h2 H2:{ext} --add-untagged")
        cases.append(f"--output-h2 H2:{ext} --output-untagged UNTAGGED:{ext}")
    with open(f"{out}/cases.txt", "w") as f:
        f.write("\n".join(cases) + "\n")


if __name__ == "__main__":
    main()
