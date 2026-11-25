"""
The pysam side of the VCF I/O oracle (see tests/vcf/gen.py).

Usage: replay.py DIR > expected.txt

Reads DIR/in.vcf as WhatsHap's VcfReader does and prints what the reads
return, writes DIR/out.vcf through WhatsHap's PhasedVcfWriter with the edits
in DIR/ops.txt applied, and runs the region queries on DIR/in.vcf.gz (created
here, with a tabix index). examples/vcf_replay.rs prints the same lines.
"""
import struct
import sys
import traceback

import pysam
from pysam import VariantFile
from whatshap.vcf import PhasedVcfWriter

d = sys.argv[1]
sys.stdout.reconfigure(line_buffering=True)
ops = [line.rstrip("\n").split(" ", 4) for line in open(d + "/ops.txt")]
mode = next(o[1] for o in ops if o[0] == "mode")
cmdline = next(" ".join(o[1:]) for o in ops if o[0] == "cmdline")
cmdline = None if cmdline == "-" else cmdline
edits = {}
for o in ops:
    if o[0] == "rec":
        edits.setdefault(int(o[1]), []).append(o[2:])


def hexf(x):
    return "%016x" % struct.unpack(">Q", struct.pack(">d", x))[0]


def show_number(v, typ):
    if v is None:
        return "."
    return str(v) if typ == "Integer" else hexf(v)


def show_gt(gt):
    if gt is None:
        return "-"
    return "[" + ",".join("." if a is None else str(a) for a in gt) + "]"


def dump_header(header):
    print("H samples=" + ",".join(header.samples))
    print("H contigs=" + ",".join(header.contigs))
    for name, meta in (("INFO", header.info), ("FORMAT", header.formats)):
        print("H %s=" % name + ",".join("%s:%s:%s" % (k, v.number, v.type) for k, v in meta.items()))
    print("H keys=" + ",".join(r.key for r in header.records))


def dump_sample(rec, header, j):
    call = rec.samples[j]
    out = ["gt=" + show_gt(call.get("GT")), "ph=%d" % call.phased]
    for key in ("GL", "PL", "HP", "PS", "PQ"):
        if key not in header.formats:
            continue
        meta = header.formats[key]
        number, typ = str(meta.number), meta.type
        if typ in ("Integer", "Float") and number == "G":
            v = call.get(key)
            if v is None:
                out.append(key + "=-")
            else:
                if not isinstance(v, tuple):
                    v = (v,)
                out.append(key + "=" + ",".join(show_number(x, typ) for x in v))
        elif typ in ("Integer", "Float") and number == "1":
            try:
                v = call[key]
            except KeyError:
                out.append(key + "=-")
            else:
                if v is None:
                    out.append(key + "=None")
                elif isinstance(v, int):
                    out.append(key + "=i:%d" % v)
                else:
                    out.append(key + "=f:" + hexf(v))
            if key == "PS":
                out.append("hasPS=%d" % ("PS" in call))
                out.append("presentPS=%d" % ("PS" in call and call["PS"] is not None))
        elif typ == "String" and number != "1":
            v = call.get(key)
            if v is None:
                out.append(key + "=-")
            elif isinstance(v, tuple):
                out.append(key + "=" + ",".join("<None>" if x is None else x for x in v))
            else:
                out.append(key + "=" + v)
    return " ".join(out)


def dump_record(prefix, i, rec, header):
    alts = "-" if rec.alts is None else ",".join(rec.alts)
    print(
        "%s %d %s %d %s %s info=%s format=%s"
        % (prefix, i, rec.chrom, rec.start, rec.ref, alts, ",".join(rec.info), ",".join(rec.format))
    )
    for j in range(len(rec.samples)):
        print("  S %d %s" % (j, dump_sample(rec, header, j)))


# 1. Reading (VcfReader)
try:
    with VariantFile(d + "/in.vcf") as f:
        dump_header(f.header)
        for i, rec in enumerate(f):
            dump_record("R", i, rec, f.header)
except Exception:
    traceback.print_exc()
    print("ERROR read")

# 2. Writing (PhasedVcfWriter with the edits applied)
try:
    w = PhasedVcfWriter(d + "/in.vcf", cmdline, d + "/out.vcf", tag=mode)
    header = w._reader.header
    dump_header(header)
    for i, rec in enumerate(w._reader):
        for e in edits.get(i, []):
            op, s = e[0], int(e[1])
            call = rec.samples[s]
            if op == "unphase":
                # PhasedVcfWriter._remove_existing_phasing
                if "GT" not in call:
                    continue
                call.phased = False
                if call["GT"] is not None and all(allele is not None for allele in call["GT"]):
                    call["GT"] = sorted(call["GT"])
            elif op == "gt":
                call["GT"] = tuple(int(a) for a in e[2].split(","))
            elif op == "setps":
                comp, phase = e[2].split(" ")
                call["PS"] = int(comp)
                call["GT"] = tuple(int(a) for a in phase.split(","))
                call.phased = True
            elif op == "ps":
                call["PS"] = None
            elif op == "hp":
                call["HP"] = None if e[2] == "." else e[2]
        if i in edits:
            dump_record("W", i, rec, header)
        w._writer.write(rec)
    w.close()
except Exception:
    traceback.print_exc()
    print("ERROR write")

# 3. Region queries (haplotag)
try:
    pysam.tabix_index(d + "/in.vcf", preset="vcf", force=True, keep_original=True)
except OSError:
    pass  # e.g. a region too long for a tbi index; then fetch() fails
with VariantFile(d + "/in.vcf.gz") as f:
    for o in ops:
        if o[0] != "fetch":
            continue
        contig, start, stop = o[1], int(o[2]), None if o[3] == "-" else int(o[3])
        print("fetch %s %d %s" % (contig, start, o[3]))
        try:
            recs = list(f.fetch(contig, start, stop))
        except ValueError as e:
            print("  ERR %s" % ("invalid contig" if "invalid contig" in str(e) else "other"))
            continue
        except Exception:
            print("  ERR other")
            continue
        for rec in recs:
            print("  F %s %d %s" % (rec.chrom, rec.start, rec.ref))
