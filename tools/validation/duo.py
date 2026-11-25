"""Father-son consistency of two lrpipe output folders (HG002 son, HG003 father), a check that
needs no truth set: on an autosome the son carries one of his father's two alleles, and on chrY
his father's one allele. A call that breaks this is wrong in the son or in the father.

usage: duo.py SON_DIR FATHER_DIR [GIAB_SON_VCF GIAB_FATHER_VCF SON_BED FATHER_BED]

SNVs (PASS, chr1-22, inside both GIAB v4.2.1 confident regions when the BEDs are given): a site
missing from one VCF is taken as reference there, so only opposing homozygotes are visible, as
both people are covered. The same test on GIAB's own truth for both people gives its floor.
SVs (PASS, chr1-22): matched by type, position within 500 bp and size ratio 0.7; a son's 1/1
without a father match, or a father's 1/1 without a son match, breaks it.
Repeats: every catalog locus genotyped in both, by allele length; exact, within 1 bp, and within
max(1 motif, 5%). Informative loci are those where either person has an allele unlike the
reference length; Tier 1 loci start inside GIAB's HG002 TR v1.0.1 Tier 1 regions, where HG002's
calls are known to be good, so a break there is mostly the father's error. chrY: the son's allele against the father's.
CNVs: the son's calls also called in the father (50% reciprocal overlap), for reference only.
$BCFTOOLS names bcftools (default: on PATH) and $GIAB_TR_TIER1 GIAB's HG002 TR Tier 1 BED."""
import gzip, os, re, subprocess, sys
from collections import defaultdict

BCF = os.environ.get("BCFTOOLS", "bcftools")
AUTO = {f"chr{i}" for i in range(1, 23)}
TIER1 = os.environ.get("GIAB_TR_TIER1", "HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz")
son, father = sys.argv[1], sys.argv[2]
giab = sys.argv[3:7] if len(sys.argv) >= 7 else None


def query(path, fmt, include=None):
    cmd = [BCF, "query", "-f", fmt] + (["-i", include] if include else []) + [path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return [line.split("\t") for line in out.splitlines()]


def find(d, *names):
    for n in names:
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    raise SystemExit(f"none of {names} in {d}")


def read_bed(path):
    iv = defaultdict(list)
    for line in (gzip.open(path, "rt") if path.endswith(".gz") else open(path)):
        c = line.split("\t")
        iv[c[0]].append((int(c[1]), int(c[2])))
    for k in iv:
        iv[k].sort()
    return iv


def inside(iv, chrom, pos):
    """pos (0-based) inside one of iv[chrom]'s sorted intervals."""
    import bisect
    a = iv.get(chrom, [])
    i = bisect.bisect_right(a, (pos, float("inf"))) - 1
    return i >= 0 and a[i][0] <= pos < a[i][1]


# --- SNVs and small indels -------------------------------------------------------------
def small(path, keep):
    """{(chrom, pos, ref, alt): dosage} for PASS biallelic-split calls on chr1-22 in keep."""
    calls = {}
    for chrom, pos, ref, alts, gt in query(path, "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n", 'FILTER="PASS" || FILTER="."'):
        if chrom not in AUTO or (keep and not keep(chrom, int(pos) - 1)):
            continue
        g = re.split(r"[/|]", gt)
        if "." in g:
            continue
        alts = alts.split(",")
        for i, a in enumerate(alts, 1):
            d = g.count(str(i))
            if d:
                calls[(chrom, int(pos), ref, a)] = d
    return calls


def opposing(s, f):
    """Sites where one is homozygous for an allele the other lacks (absent = reference), split
    into SNPs and indels, over sites variant in either."""
    res = {"snp": [0, 0], "indel": [0, 0]}
    for k in set(s) | set(f):
        kind = "snp" if len(k[2]) == 1 and len(k[3]) == 1 else "indel"
        res[kind][0] += 1
        a, b = s.get(k, 0), f.get(k, 0)
        if (a == 2 and b == 0) or (b == 2 and a == 0):
            res[kind][1] += 1
    return res


keep = None
if giab:
    sb, fb = read_bed(giab[2]), read_bed(giab[3])
    keep = lambda c, p: inside(sb, c, p) and inside(fb, c, p)
snv = lambda d: find(d, "snv/merge_output.vcf.gz", "snv/HG002.vcf.gz", "snv/HG003.vcf.gz")
print("== SNVs and indels (opposing homozygotes / sites variant in either)")
rows = [("lrpipe", snv(son), snv(father))]
if giab:
    rows.append(("GIAB truth", giab[0], giab[1]))
for name, a, b in rows:
    r = opposing(small(a, keep), small(b, keep))
    print(f"{name:11s} SNPs {r['snp'][1]:>7,}/{r['snp'][0]:<10,} {100 * r['snp'][1] / r['snp'][0]:.3f}%   "
          f"indels {r['indel'][1]:>7,}/{r['indel'][0]:<9,} {100 * r['indel'][1] / r['indel'][0]:.3f}%")


# --- SVs -------------------------------------------------------------------------------
def svs(path):
    out = []
    for chrom, pos, t, ln, gt in query(path, "%CHROM\t%POS\t%SVTYPE\t%SVLEN\t[%GT]\n", 'FILTER="PASS"'):
        if chrom not in AUTO or t not in ("INS", "DEL", "DUP", "INV"):
            continue
        g = re.split(r"[/|]", gt)
        if "." in g:
            continue
        ln = abs(int(ln)) if ln not in (".", "") else 0
        out.append((chrom, int(pos), t, ln, sum(x != "0" for x in g)))
    return out


def sv_match(a, idx):
    for b in idx.get((a[0], a[2], a[1] // 1000), []) + idx.get((a[0], a[2], a[1] // 1000 - 1), []) + idx.get((a[0], a[2], a[1] // 1000 + 1), []):
        if abs(a[1] - b[1]) <= 500 and (min(a[3], b[3]) >= 0.7 * max(a[3], b[3]) or a[2] == "INV"):
            return b
    return None


def index(calls):
    idx = defaultdict(list)
    for c in calls:
        idx[(c[0], c[2], c[1] // 1000)].append(c)
    return idx


sv_path = lambda d: find(d, "sv/HG002.sv.vcf.gz", "sv/HG003.sv.vcf.gz", "sv/HG002.sv.vcf", "sv/HG003.sv.vcf")
S, F = svs(sv_path(son)), svs(sv_path(father))
Si, Fi = index(S), index(F)
print("\n== SVs (PASS, chr1-22)")
for who, calls, other in (("son", S, Fi), ("father", F, Si)):
    hom = [c for c in calls if c[4] == 2]
    het = [c for c in calls if c[4] == 1]
    hom_bad = sum(sv_match(c, other) is None for c in hom)
    het_seen = sum(sv_match(c, other) is not None for c in het)
    print(f"{who:6s} {len(calls):,} calls; 1/1 not in the other: {hom_bad:,}/{len(hom):,} ({100 * hom_bad / max(1, len(hom)):.2f}%); "
          f"0/1 also in the other: {het_seen:,}/{len(het):,} ({100 * het_seen / max(1, len(het)):.1f}%)")
for t in ("INS", "DEL"):
    for lo, hi in ((50, 300), (300, 1000), (1000, 10 ** 9)):
        sub = [c for c in S + F if c[2] == t and lo <= c[3] < hi and c[4] == 2]
        bad = sum(sv_match(c, Fi if c in S else Si) is None for c in sub)
        print(f"  {t} {lo}-{hi if hi < 10 ** 9 else ''}: 1/1 not in the other {bad}/{len(sub)} ({100 * bad / max(1, len(sub)):.2f}%)")


# --- Tandem repeats ---------------------------------------------------------------------
def repeats(path):
    out = {}
    for chrom, pos, trid, end, motifs, gt, al in query(path, "%CHROM\t%POS\t%TRID\t%END\t%MOTIFS\t[%GT]\t[%AL]\n"):
        if "." in gt or al in (".", ""):
            continue
        out[trid] = (chrom, int(pos), int(end) - int(pos), len(motifs.split(",")[0]), [int(x) for x in al.split(",") if x != "."])
    return out


rep = lambda d: find(d, "repeats/HG002.repeats.vcf.gz", "repeats/HG003.repeats.vcf.gz")
RS, RF = repeats(rep(son)), repeats(rep(father))
both = [k for k in RS if k in RF]
tols = {"exact": lambda a, b, k: a == b, "1 bp": lambda a, b, k: abs(a - b) <= 1,
        "1 unit/5%": lambda a, b, k: abs(a - b) <= max(k, 0.05 * max(a, b))}
tier1 = read_bed(TIER1) if os.path.exists(TIER1) else {}
print(f"\n== Tandem repeats ({len(both):,} loci genotyped in both; consistent = share an allele)")
for scope in ("autosomes", "informative", "Tier 1", "chrY"):
    counts = {t: [0, 0] for t in tols}
    bins = defaultdict(lambda: [0, 0])
    for k in both:
        chrom, pos, rl, unit, s = RS[k]
        f = RF[k][4]
        if scope == "chrY":
            if chrom != "chrY":
                continue
        elif chrom not in AUTO:
            continue
        if scope == "informative" and all(x == rl for x in s + f):
            continue
        if scope == "Tier 1" and not inside(tier1, chrom, pos):
            continue
        for t, ok in tols.items():
            counts[t][0] += 1
            counts[t][1] += any(ok(a, b, unit) for a in s for b in f)
        if scope in ("informative", "Tier 1"):
            d = max(abs(x - rl) for x in s + f)
            bins["change " + ("0" if d == 0 else "1-4" if d < 5 else "5-49" if d < 50 else ">=50")][0] += 1
            bins["change " + ("0" if d == 0 else "1-4" if d < 5 else "5-49" if d < 50 else ">=50")][1] += any(a == b_ for a in s for b_ in f)
            m = "motif " + ("1" if unit == 1 else "2-6" if unit <= 6 else ">6")
            bins[m][0] += 1
            bins[m][1] += any(a == b_ for a in s for b_ in f)
    n = counts["exact"][0]
    print(f"{scope:12s} {n:>9,} loci  " + "  ".join(f"{t} {100 * g / max(1, n):.2f}%" for t, (_, g) in counts.items()))
    for b in ("change 0", "change 1-4", "change 5-49", "change >=50", "motif 1", "motif 2-6", "motif >6"):
        if b in bins:
            print(f"    {b:12s} {bins[b][0]:>8,} loci, exact {100 * bins[b][1] / bins[b][0]:.2f}%")


# --- CNVs (reference only) --------------------------------------------------------------
def cnvs(path):
    out = []
    for chrom, pos, end, t in query(path, "%CHROM\t%POS\t%INFO/END\t%SVTYPE\n", 'FILTER="PASS"'):
        if chrom in AUTO:
            out.append((chrom, int(pos), int(end), t))
    return out


try:
    cs = cnvs(find(son, "cnv/HG002.vcf.gz"))
    cf = cnvs(find(father, "cnv/HG003.vcf.gz"))
    shared = sum(any(c[0] == d[0] and c[3] == d[3] and min(c[2], d[2]) - max(c[1], d[1]) >= 0.5 * max(c[2] - c[1], d[2] - d[1]) for d in cf) for c in cs)
    print(f"\n== CNVs (PASS, chr1-22): son {len(cs)}, father {len(cf)}, son's also in father {shared}")
except SystemExit as e:
    print(f"\n== CNVs skipped: {e}")
