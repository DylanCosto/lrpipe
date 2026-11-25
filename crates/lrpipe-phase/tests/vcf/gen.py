"""
Random VCFs and record edits for the VCF I/O oracle.

Usage: gen.py SEED DIR

Writes DIR/in.vcf (a random VCF that htslib accepts) and DIR/ops.txt, the
edits and queries that tests/vcf/replay.py (pysam, the oracle) and
examples/vcf_replay.rs (lrpipe-phase) both apply. ops.txt lines:

    mode PS|HP                 phase tag of the output (PhasedVcfWriter tag)
    cmdline TEXT | cmdline -   command line for the output header (- for none)
    fetch CONTIG START STOP|-  region query on the bgzipped, tabix-indexed input
    rec I unphase S            _remove_existing_phasing for sample S of record I
    rec I gt S A,B             call["GT"] = (A, B)
    rec I setps S C A,B        _set_PS: call["PS"] = C; call["GT"] = (A, B); call.phased = True
    rec I ps S .               call["PS"] = None
    rec I hp S VALUE|.         call["HP"] = VALUE (. for None)

The records cover what htslib has to handle when WhatsHap passes a VCF
through: all INFO/FORMAT types and numbers, missing values, missing trailing
sample fields, undeclared contigs/INFO/FORMAT/FILTER, symbolic alleles with
and without END, integers of every BCF width and floats in many spellings.
"""
import random
import string
import sys

seed = int(sys.argv[1])
out_dir = sys.argv[2]
rng = random.Random(seed)


def chance(p):
    return rng.random() < p


def rint():
    c = rng.random()
    if c < 0.5:
        return rng.randint(-5, 100)
    if c < 0.65:
        return rng.randint(-120, 127)
    if c < 0.75:
        return rng.choice([-128, -127, -121, -120, 127, 128, -32768, -32767, -32761, -32760, 32767, 32768])
    if c < 0.9:
        return rng.randint(-32760, 32767)
    return rng.choice([rng.randint(-2147483640, 2147483647), 2147483647, -2147483640, -2147483641, 65536])


def rfloat():
    c = rng.random()
    if c < 0.3:
        return str(round(rng.uniform(-100, 100), rng.randint(0, 6)))
    if c < 0.45:
        return "%.*e" % (rng.randint(0, 8), rng.uniform(-1, 1) * 10 ** rng.randint(-45, 39))
    if c < 0.6:
        return str(rng.randint(-1000, 100000))
    if c < 0.7:
        return repr(rng.random())
    if c < 0.8:
        return "%.*g" % (rng.randint(1, 12), rng.uniform(0, 1) * 10 ** rng.randint(-10, 10))
    if c < 0.9:
        return rng.choice(
            ["0", "-0", "0.0", "1e-45", "3.4e38", "1e39", "-1e39", "inf", "-inf", "nan", "NaN", "Inf",
             ".5", "5.", "+3", "1E5", "1e+05", "0.1", "0.30000000000000004", "123456789", "1.00000001",
             "-0.000001", "99.99999", "1e-7", "2.5e-8", "16777217", "0.333333333"]
        )
    return str(round(rng.uniform(-1, 0), rng.randint(1, 3)))


def rstring(commas=True):
    alphabet = string.ascii_letters + string.digits + "_-+./|"
    if commas:
        alphabet += ","
    s = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 8)))
    return s


def rdesc():
    d = rng.choice(["Depth", "Allele frequency", "Some value, with a comma", "Quality (phred)", "x", "A=B; c"])
    if chance(0.1):
        d += ' with \\"quotes\\"'
    return d


# --- header ---------------------------------------------------------------

n_samples = rng.choice([0, 1, 1, 1, 2, 2, 3, 4])
samples = rng.sample(["HG001", "HG002", "HG003", "NA12878", "s1", "sample_2", "child", "mother"], n_samples)
contigs = rng.sample(["chr1", "chr2", "chr10", "1", "2", "X", "chrM", "chrUn_KI270742v1", "HLA-A*01:01"], rng.randint(1, 4))
tag_mode = rng.choice(["PS", "HP"])

lines = ["##fileformat=VCFv4." + rng.choice(["1", "2", "2", "2", "3"])]
if chance(0.5):
    lines.append("##fileDate=20%02d%02d%02d" % (rng.randint(10, 25), rng.randint(1, 12), rng.randint(1, 28)))
if chance(0.5):
    lines.append("##source=" + rng.choice(["freeBayes v1.3", "GATK", "DeepVariant", "clair3"]))
if chance(0.2):
    lines.append("##phasing=none")
if chance(0.3):
    lines.append("##reference=file:///data/ref.fa")
if chance(0.2):
    lines.append('##bcftools_viewCommand=view -Oz in.vcf; Date=Mon Jan  1 00:00:00 2024')

filters = []
if chance(0.7):
    filters = rng.sample(["q10", "s50", "LowQual", "RefCall"], rng.randint(1, 3))
    for f in filters:
        lines.append('##FILTER=<ID=%s,Description="%s">' % (f, rdesc()))
if chance(0.3):
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')

declared_contigs = [c for c in contigs if chance(0.8)]
for c in declared_contigs:
    extra = ""
    if chance(0.8):
        extra += ",length=%d" % rng.randint(10**5, 10**9)
    if chance(0.2):
        extra += ",assembly=b37"
    lines.append("##contig=<ID=%s%s>" % (c, extra))

if chance(0.3):
    lines.append('##ALT=<ID=DEL,Description="Deletion">')

# INFO fields: (id, number, type)
info_pool = [
    ("DP", "1", "Integer"), ("AF", "A", "Float"), ("MQ", "1", "Float"), ("DB", "0", "Flag"),
    ("AC", "A", "Integer"), ("AN", "1", "Integer"), ("CIGAR", "A", "String"), ("TYPE", ".", "String"),
    ("RO", "1", "Integer"), ("AO", "A", "Integer"), ("QR", "R", "Integer"), ("GLX", "G", "Float"),
    ("XC", "1", "Character"), ("X2", "2", "Integer"), ("XF3", "3", "Float"), ("XS", ".", "Integer"),
    ("SVTYPE", "1", "String"), ("SVLEN", ".", "Integer"), ("END", "1", "Integer"), ("FL", ".", "Float"),
]
infos = rng.sample(info_pool, rng.randint(0, len(info_pool)))
for i, n, t in infos:
    extra = ""
    if chance(0.1):
        extra = ',Source="tool",Version="1.0"'
    lines.append('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s"%s>' % (i, n, t, rdesc(), extra))
# INFO fields used in records without a definition (WhatsHap adds known ones)
# At most one: WhatsHap adds several in Python set order, which depends on the
# hash seed. END counts as used whenever a record has a symbolic allele.
undeclared_infos = [x for x in [("AC", "A", "Integer"), ("SVTYPE", "1", "String"), ("END", "1", "Integer")]
                    if x[0] not in {i[0] for i in infos} and chance(0.15)][:1]
if undeclared_infos and undeclared_infos[0][0] != "END" and "END" not in {i[0] for i in infos}:
    infos.append(("END", "1", "Integer"))
    lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End">')

# FORMAT fields
fmt_pool = [
    ("GQ", "1", rng.choice(["Integer", "Integer", "Float"])),
    ("DP", "1", "Integer"),
    ("AD", rng.choice(["R", "R", "."]), "Integer"),
    ("GL", "G", "Float"),
    ("PL", "G", "Integer"),
    ("PQ", "1", rng.choice(["Float", "Float", "Integer"])),
    ("FT", "1", "String"),
    ("MIN_DP", "1", "Integer"),
    ("VAF", "A", "Float"),
    ("XS", ".", "String"),
    ("XC", "1", "Character"),
    ("XI2", "2", "Integer"),
    ("LEN", "1", "Integer"),
]
formats = []
if n_samples:
    gt_def = ("GT", "1", "String") if chance(0.95) else ("GT", rng.choice([".", "2"]), "String")
    formats = [gt_def] + rng.sample(fmt_pool, rng.randint(0, len(fmt_pool)))
    if chance(0.5):
        formats.append(("PS", "1", "Integer" if chance(0.85) else "Float"))
    if chance(0.3):
        formats.append(("HP", rng.choice([".", ".", "1"]), "String"))
    if chance(0.1):
        formats.append(("HS", ".", "Integer"))
    rng.shuffle(formats)
    for i, n, t in formats:
        lines.append('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">' % (i, n, t, rdesc()))
undeclared_formats = []
if n_samples:
    for x in [("PS", "1", "Integer"), ("PQ", "1", "Float"), ("AD", ".", "Integer"), ("GQ", "1", "Integer")]:
        if x[0] not in {f[0] for f in formats} and chance(0.15):
            undeclared_formats.append(x)

cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
if n_samples:
    cols += ["FORMAT"] + samples
lines.append("\t".join(cols))


# --- records --------------------------------------------------------------

def values(number, typ, n_alt, ploidy=2):
    n_all = n_alt + 1
    if number == "A":
        n = n_alt
    elif number == "R":
        n = n_all
    elif number == "G":
        n = n_all * (n_all + 1) // 2 if ploidy == 2 else n_all
    elif number == ".":
        n = rng.randint(1, 4)
    else:
        n = int(number)
    if number != "1" and chance(0.05):
        n = max(1, n + rng.choice([-1, 1]))
    out = []
    for _ in range(max(n, 1)):
        if chance(0.08):
            out.append(".")
        elif typ == "Integer":
            out.append(str(rint()))
        elif typ == "Float":
            out.append(rfloat())
        elif typ == "Character":
            out.append(rng.choice("ACGTxyz"))
        else:
            out.append(rstring(commas=False))
    if chance(0.05):
        return "."
    return ",".join(out)


def genotype(n_all):
    ploidy = rng.choice([2, 2, 2, 2, 2, 1, 3])
    if chance(0.08):
        return rng.choice([".", "./.", ".|.", "./1", "0|.", "./0"]), 2
    alleles = []
    for _ in range(ploidy):
        alleles.append("." if chance(0.03) else str(rng.randrange(n_all)))
    seps = [rng.choice("/|") for _ in range(ploidy - 1)]
    if chance(0.7):
        seps = [seps[0]] * (ploidy - 1) if seps else []
    gt = alleles[0]
    for s, a in zip(seps, alleles[1:]):
        gt += s + a
    return gt, ploidy


records = []
contig_order = rng.sample(contigs, len(contigs))
for contig in contig_order:
    n = rng.randint(1, 25)
    positions = sorted(rng.randint(1, 2000000) for _ in range(n))
    if chance(0.3) and n > 1:
        positions[1] = positions[0]
    for pos in positions:
        ref_len = rng.choice([1, 1, 1, 1, 2, 3, 5])
        ref = "".join(rng.choice("ACGT") for _ in range(ref_len))
        if chance(0.02):
            ref = "N"
        c = rng.random()
        if c < 0.05:
            alts = []
        elif c < 0.1:
            sym = rng.choice(["<DEL>", "<DEL>", "<DUP>", "<INV>", "<CNV>", "<DEL:ME>", "<INS>", "<*>", "<NON_REF>", "<DUP:TANDEM>"])
            alts = [sym] if chance(0.7) else [sym, "A"]
        elif c < 0.13:
            alts = ["*", rng.choice("ACGT")]
        else:
            alts = []
            for _ in range(rng.choice([1, 1, 1, 1, 2, 3])):
                a = "".join(rng.choice("ACGT") for _ in range(rng.choice([1, 1, 1, 2, 4])))
                if a != ref and a not in alts:
                    alts.append(a)
            if not alts:
                alts = ["T" if ref != "T" else "G"]
        n_alt = len(alts)
        n_all = n_alt + 1
        id_ = "." if chance(0.6) else rng.choice(["rs%d" % rng.randint(1, 10**8), "id1;id2", "x"])
        qual = "." if chance(0.2) else rfloat().lstrip("-")
        if qual in ("inf", "nan", "NaN", "Inf", ""):
            qual = "50"
        if filters and chance(0.5):
            filt = ";".join(rng.sample(filters + ["PASS"], rng.randint(1, 2)))
        else:
            filt = rng.choice([".", "PASS", "PASS"])
        if chance(0.03):
            filt = "undeclaredF"
        info_parts = []
        for i, nn, t in infos + undeclared_infos:
            if i == "END":
                if chance(0.4) or (alts and alts[0].startswith("<") and chance(0.7)):
                    end = pos + ref_len - 1 if chance(0.4) else pos + rng.randint(1, 5000)
                    info_parts.append("END=%d" % end)
                continue
            if not chance(0.5):
                continue
            if t == "Flag":
                info_parts.append(i)
            else:
                info_parts.append("%s=%s" % (i, values(nn, t, max(n_alt, 1))))
        if chance(0.02):
            info_parts.append("UNDECL=1")
        rng.shuffle(info_parts) if chance(0.2) else None
        info = ";".join(info_parts) if info_parts else "."
        fields = [contig, str(pos), id_, ref, ",".join(alts) if alts else ".", qual, filt, info]
        if n_samples:
            fmt_keys = [f for f in formats + undeclared_formats if f[0] == "GT" or chance(0.6)]
            gt_first = [f for f in fmt_keys if f[0] == "GT"]
            rest = [f for f in fmt_keys if f[0] != "GT"]
            rng.shuffle(rest)
            if chance(0.97):
                fmt_keys = gt_first + rest
            else:
                fmt_keys = rest + gt_first
            if chance(0.03):
                fmt_keys = [f for f in fmt_keys if f[0] != "GT"]
            if not fmt_keys:
                fmt_keys = [("GT", "1", "String")]
            fields.append(":".join(f[0] for f in fmt_keys))
            for _ in samples:
                vals = []
                ploidy = 2
                for k, nn, t in fmt_keys:
                    if k == "GT":
                        g, ploidy = genotype(n_all)
                        vals.append(g)
                    elif k == "PS" and t == "Integer":
                        vals.append("." if chance(0.3) else str(rng.choice([rng.randint(1, 100), rng.randint(1, 10**6), pos])))
                    elif k == "HP":
                        vals.append("." if chance(0.4) else "%d-1,%d-2" % (pos, pos))
                    else:
                        vals.append(values(nn, t, max(n_alt, 1), ploidy))
                # (pysam fails if no sample has a value for a GT that is not first)
                if fmt_keys[0][0] == "GT" and chance(0.15) and len(vals) > 1:
                    vals = vals[: rng.randint(1, len(vals) - 1)]
                if fmt_keys[0][0] == "GT" and chance(0.03):
                    vals = ["."]
                fields.append(":".join(vals))
        records.append((contig, pos, n_all, alts, fields))
        lines.append("\t".join(fields))

with open(out_dir + "/in.vcf", "w") as f:
    f.write("\n".join(lines) + "\n")


# --- edits and queries ----------------------------------------------------

ops = ["mode " + tag_mode]
if chance(0.8):
    ops.append("cmdline whatshap phase -o out.vcf --tag=%s in.vcf reads.bam%s" % (tag_mode, ' "quoted"' if chance(0.2) else ""))
else:
    ops.append("cmdline -")
for _ in range(rng.randint(0, 10)):
    c = rng.choice(contigs + ["chrNotThere"])
    near = [p for (cc, p, _, _, _) in records if cc == c]
    if near and chance(0.6):
        # just after a record, to catch records whose end reaches into the region
        start = rng.choice(near) + rng.randint(0, 40000)
    else:
        start = rng.randint(0, 2000000)
    stop = "-" if chance(0.2) else str(start + rng.randint(1, 500000))
    ops.append("fetch %s %d %s" % (c, start, stop))

have_gt = {i for i, (_, _, _, _, fields) in enumerate(records)
           if n_samples and fields[8].split(":")[0] == "GT"}
for i, (contig, pos, n_all, alts, fields) in enumerate(records):
    if i not in have_gt or not chance(0.7):
        continue
    for s in range(n_samples):
        if not chance(0.7):
            continue
        if tag_mode == "PS" and chance(0.6):
            ops.append("rec %d unphase %d" % (i, s))
        if len(alts) != 1:
            continue
        a, b = rng.randrange(2), rng.randrange(2)
        if chance(0.2):
            ops.append("rec %d gt %d %d,%d" % (i, s, a, b))
        if tag_mode == "PS":
            if chance(0.6):
                comp = rng.choice([rng.randint(1, 100), rng.randint(1, 40000), rng.randint(1, 2 * 10**9), pos])
                ops.append("rec %d setps %d %d %d,%d" % (i, s, comp, a, b))
            else:
                ops.append("rec %d ps %d ." % (i, s))
        else:
            if chance(0.6):
                c = rng.randint(1, 10**7)
                ops.append("rec %d hp %d %d-%d,%d-%d" % (i, s, c, a + 1, c, b + 1))
            else:
                ops.append("rec %d hp %d ." % (i, s))

with open(out_dir + "/ops.txt", "w") as f:
    f.write("\n".join(ops) + "\n")
