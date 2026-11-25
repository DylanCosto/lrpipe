"""Mendelian consistency of three lrpipe output folders (HG002 son, HG003 father, HG004 mother):
at an autosomal locus the son has one allele from each parent; on chrX his one allele is one of
his mother's, on chrY his father's. A call that breaks this is wrong in at least one of the three.

usage: trio.py SON_DIR FATHER_DIR MOTHER_DIR [BED ...]

SNVs and indels (PASS, chr1-22, inside every BED given, such as GIAB's confident regions): a site
missing from a VCF is taken as reference there, which holds where all three are covered.
Calls are matched as written, so normalise the three snv/merge_output.vcf.gz first
(`bcftools norm -f GRCh38.fa -m -any`) for BENCHMARKS.md's indel rate; unnormalised, an indel
written differently in two VCFs counts as an error.
$BCFTOOLS names bcftools (default: on PATH) and $GIAB_TR_TIER1 GIAB's HG002 TR Tier 1 BED.
SVs (PASS, chr1-22): each call is matched in the other two by type, position within 500 bp and
size ratio 0.7 (absent = 0/0); every SV of any of the three is checked once.
Repeats: every catalog locus genotyped in all three, by allele length; exact, within 1 bp, and
within max(1 motif, 5%), split by motif length and by the largest change from the reference."""
import bisect, gzip, os, re, subprocess, sys
from collections import defaultdict

BCF = os.environ.get("BCFTOOLS", "bcftools")
AUTO = {f"chr{i}" for i in range(1, 23)}
TIER1 = os.environ.get("GIAB_TR_TIER1", "HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz")
dirs = sys.argv[1:4]
beds = sys.argv[4:]
ONLY = set(os.environ.get("TRIO_ONLY", "snv,sv,tr").split(","))
names = ["son", "father", "mother"]


def query(path, fmt, include=None):
    cmd = [BCF, "query", "-f", fmt] + (["-i", include] if include else []) + [path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return [line.split("\t") for line in out.splitlines()]


def find(d, kind):
    for n in os.listdir(os.path.join(d, kind)):
        if kind == "snv" and n == "merge_output.vcf.gz" or kind != "snv" and re.match(rf"HG00\d\.{kind}\.vcf(\.gz)?$", n):
            return os.path.join(d, kind, n)
    raise SystemExit(f"no {kind} VCF in {d}")


def read_bed(path):
    iv = defaultdict(list)
    for line in (gzip.open(path, "rt") if path.endswith(".gz") else open(path)):
        c = line.split("\t")
        iv[c[0]].append((int(c[1]), int(c[2])))
    for k in iv:
        iv[k].sort()
    return iv


def inside(iv, chrom, pos):
    a = iv.get(chrom, [])
    i = bisect.bisect_right(a, (pos, float("inf"))) - 1
    return i >= 0 and a[i][0] <= pos < a[i][1]


def mendelian(c, f, m):
    """Son dosage c from one of father's alleles and one of mother's (dosages 0, 1, 2)."""
    al = lambda d: {0} if d == 0 else {1} if d == 2 else {0, 1}
    return any(a + b == c for a in al(f) for b in al(m))


# --- SNVs and small indels -------------------------------------------------------------
regions = [read_bed(b) for b in beds]


def small(path):
    calls = {}
    for chrom, pos, ref, alts, gt in query(path, "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n", 'FILTER="PASS" || FILTER="."'):
        if chrom not in AUTO or not all(inside(r, chrom, int(pos) - 1) for r in regions):
            continue
        g = re.split(r"[/|]", gt)
        if "." in g:
            continue
        for i, a in enumerate(alts.split(","), 1):
            if g.count(str(i)):
                calls[(chrom, int(pos), ref, a)] = g.count(str(i))
    return calls


S = [small(find(d, "snv")) for d in dirs] if "snv" in ONLY else [{}, {}, {}]
res = {"snp": [0, 0], "indel": [0, 0]}
for k in set().union(*S):
    kind = "snp" if len(k[2]) == 1 and len(k[3]) == 1 else "indel"
    res[kind][0] += 1
    res[kind][1] += not mendelian(*(s.get(k, 0) for s in S))
print("== SNVs and indels (Mendelian errors / sites variant in any)")
for kind, (n, bad) in res.items():
    print(f"{kind:6s} {bad:>8,}/{n:<10,} {100 * bad / max(1, n):.3f}%")


# --- SVs -------------------------------------------------------------------------------
def svs(path):
    out = []
    for chrom, pos, t, ln, gt in query(path, "%CHROM\t%POS\t%SVTYPE\t%SVLEN\t[%GT]\n", 'FILTER="PASS"'):
        if chrom not in AUTO or t not in ("INS", "DEL", "DUP", "INV"):
            continue
        g = re.split(r"[/|]", gt)
        if "." in g:
            continue
        out.append((chrom, int(pos), t, abs(int(ln)) if ln not in (".", "") else 0, sum(x != "0" for x in g)))
    return out


def index(calls):
    idx = defaultdict(list)
    for c in calls:
        idx[(c[0], c[2], c[1] // 1000)].append(c)
    return idx


def match(a, idx):
    for w in (-1, 0, 1):
        for b in idx.get((a[0], a[2], a[1] // 1000 + w), []):
            if abs(a[1] - b[1]) <= 500 and (a[2] == "INV" or min(a[3], b[3]) >= 0.7 * max(a[3], b[3])):
                return b
    return None


V = [svs(find(d, "sv")) for d in dirs] if "sv" in ONLY else [[], [], []]
I = [index(v) for v in V]
seen = set()
bins = defaultdict(lambda: [0, 0])
for who, calls in enumerate(V):
    for c in calls:
        hits = [c if j == who else match(c, I[j]) for j in range(3)]
        key = tuple(h for h in hits if h)
        if any(h in seen for h in key):
            continue
        seen.update(key)
        ok = mendelian(*(h[4] if h else 0 for h in hits))
        size = "50-299" if c[3] < 300 else "300-999" if c[3] < 1000 else ">=1000"
        for b in ("all", f"{c[2]} {size}"):
            bins[b][0] += 1
            bins[b][1] += not ok
print(f"\n== SVs (PASS, chr1-22; calls: son {len(V[0]):,}, father {len(V[1]):,}, mother {len(V[2]):,})")
for b in sorted(bins, key=lambda x: (x != "all", x)):
    n, bad = bins[b]
    print(f"{b:16s} {bad:>6,}/{n:<8,} Mendelian errors {100 * bad / max(1, n):.2f}%")


# --- Tandem repeats ---------------------------------------------------------------------
def repeats(path):
    out = {}
    for chrom, pos, trid, end, motifs, gt, al in query(path, "%CHROM\t%POS\t%TRID\t%END\t%MOTIFS\t[%GT]\t[%AL]\n"):
        if "." in gt or al in (".", ""):
            continue
        out[trid] = (chrom, int(pos), int(end) - int(pos), len(motifs.split(",")[0]), [int(x) for x in al.split(",") if x != "."])
    return out


R = [repeats(find(d, "repeats")) for d in dirs] if "tr" in ONLY else [{}, {}, {}]
tols = {"exact": lambda a, b, k: a == b, "1 bp": lambda a, b, k: abs(a - b) <= 1,
        "1 unit/5%": lambda a, b, k: abs(a - b) <= max(k, 0.05 * max(a, b))}


def consistent(chrom, c, f, m, ok, k):
    has = lambda x, alleles: any(ok(x, y, k) for y in alleles)
    if chrom == "chrX":
        return all(has(x, m) for x in c)
    if chrom == "chrY":
        return all(has(x, f) for x in c)
    if len(c) == 1:
        c = c * 2
    return (has(c[0], f) and has(c[1], m)) or (has(c[1], f) and has(c[0], m))


tier1 = read_bed(TIER1)
# chrY: the mother has none, so son and father alone.
loci = [k for k in R[0] if k in R[1] and (k in R[2] or R[0][k][0] == "chrY")]
print(f"\n== Tandem repeats ({len(loci):,} loci genotyped in all three, or son and father on chrY)")
for scope in ("autosomes", "informative", "Tier 1", "chrX", "chrY"):
    counts = {t: [0, 0] for t in tols}
    strata = defaultdict(lambda: [0, 0])
    for k in loci:
        chrom, pos, rl, unit, c = R[0][k]
        f, m = R[1][k][4], R[2][k][4] if k in R[2] else []
        if scope in ("chrX", "chrY"):
            if chrom != scope:
                continue
        elif chrom not in AUTO:
            continue
        if scope == "informative" and all(x == rl for x in c + f + m):
            continue
        if scope == "Tier 1" and not inside(tier1, chrom, pos):
            continue
        for t, ok in tols.items():
            counts[t][0] += 1
            counts[t][1] += consistent(chrom, c, f, m, ok, unit)
        if scope in ("informative", "Tier 1"):
            d = max(abs(x - rl) for x in c + f + m)
            good = consistent(chrom, c, f, m, tols["exact"], unit)
            for s in ("change " + ("0" if d == 0 else "1-4" if d < 5 else "5-49" if d < 50 else ">=50"),
                      "motif " + ("1" if unit == 1 else "2-6" if unit <= 6 else ">6")):
                strata[s][0] += 1
                strata[s][1] += good
    n = counts["exact"][0]
    print(f"{scope:12s} {n:>9,} loci  consistent " + "  ".join(f"{t} {100 * g / max(1, n):.2f}%" for t, (_, g) in counts.items()))
    for s in ("change 0", "change 1-4", "change 5-49", "change >=50", "motif 1", "motif 2-6", "motif >6"):
        if s in strata:
            print(f"    {s:12s} {strata[s][0]:>8,} loci, exact {100 * strata[s][1] / strata[s][0]:.2f}%")
