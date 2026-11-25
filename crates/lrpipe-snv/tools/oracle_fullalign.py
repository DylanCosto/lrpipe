"""Reference full-alignment tensors from Clair3's own `CreateTensorFullAlignment`.

Two modes, both writing PREFIX.tensors.npy (int8, candidates x 89 x 33 x 8, exactly what
`CreateTensorFullAlignment` returns), PREFIX.info.tsv (`position<TAB>alt_info`, in order) and
PREFIX.seconds (wall time of the `CreateTensorFullAlignment` call alone):

* `pipeline`: runs `CallVariantsFromCffi.main()` with the given arguments (copy them from
  `log/parallel_6_call_var_bam_full_alignment.log`), so the tensors are built in the same process
  state as in a real Clair3 run (torch imported, model loaded; this matters because the C code
  downsamples reads with glibc `rand()`). The run stops right after the tensors are built: no
  model inference and no VCF.

      python tools/oracle_fullalign.py pipeline --clair3-dir $ENV/bin --out OUT/c1 -- \
          --chkpnt_fn ... --bam_fn ... --full_aln_regions ... --phased_vcf_fn ... --platform ont

* `direct`: calls `CreateTensorFullAlignment` in a fresh process with the arguments a
  full-alignment job uses (for synthetic inputs).

      python tools/oracle_fullalign.py direct --clair3-dir $ENV/bin --bam ... --ref ... \
          --ctg-name chr20 --full-aln-regions bed [--phased-vcf vcf.gz] [--no-phasing] --out OUT/x
"""
import argparse
import os
import sys
import time


class _Stop(Exception):
    pass


def _save(out, result, seconds):
    import numpy as np

    tensor, positions, alt_info = result
    tensor = np.asarray(tensor, dtype=np.int8)
    np.save(out + ".tensors.npy", tensor)
    with open(out + ".info.tsv", "w") as fh:
        for pos_info, alt in zip(positions, alt_info):
            fh.write(f"{pos_info}\t{alt}\n")
    with open(out + ".seconds", "w") as fh:
        fh.write(f"{seconds:.6f}\n")
    print(f"{len(positions)} candidates in {seconds:.3f} s", file=sys.stderr)


def pipeline(a, rest):
    sys.path.insert(0, a.clair3_dir)
    import preprocess.CreateTensorFullAlignmentFromCffi as fa

    original = fa.CreateTensorFullAlignment

    def capture(args):
        start = time.perf_counter()
        result = original(args)
        _save(a.out, result, time.perf_counter() - start)
        raise _Stop()

    fa.CreateTensorFullAlignment = capture
    from clair3 import CallVariantsFromCffi

    # The VCF header goes to a scratch file; nothing else is written.
    if "--call_fn" in rest:
        rest[rest.index("--call_fn") + 1] = a.out + ".header.vcf"
    sys.argv = ["CallVariantsFromCffi"] + rest
    try:
        CallVariantsFromCffi.main()
    except _Stop:
        return
    sys.exit("CreateTensorFullAlignment was never called")


def direct(a):
    sys.path.insert(0, a.clair3_dir)
    import shared.param_f as param
    from preprocess.CreateTensorFullAlignmentFromCffi import CreateTensorFullAlignment

    args = argparse.Namespace(
        platform=a.platform, bam_fn=a.bam, ref_fn=a.ref, tensor_can_fn="PIPE", vcf_fn=None,
        ctgName=a.ctg_name, ctgStart=None, ctgEnd=None, bed_fn=None, extend_bed=None,
        full_aln_regions=a.full_aln_regions, phased_vcf_fn=a.phased_vcf,
        no_phasing_for_fa=a.no_phasing, minMQ=a.min_mq, minBQ=param.min_bq,
        enable_long_indel=a.enable_long_indel, enable_dwell_time=False,
    )
    start = time.perf_counter()
    result = CreateTensorFullAlignment(args)
    _save(a.out, result, time.perf_counter() - start)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="mode", required=True)
    pp = sub.add_parser("pipeline")
    pp.add_argument("--clair3-dir", required=True, help="directory holding clair3.py")
    pp.add_argument("--out", required=True)
    pd = sub.add_parser("direct")
    pd.add_argument("--clair3-dir", required=True)
    pd.add_argument("--bam", required=True)
    pd.add_argument("--ref", required=True)
    pd.add_argument("--ctg-name", required=True)
    pd.add_argument("--full-aln-regions", required=True)
    pd.add_argument("--phased-vcf")
    pd.add_argument("--no-phasing", action="store_true")
    pd.add_argument("--min-mq", type=int, default=5)
    pd.add_argument("--enable-long-indel", action="store_true")
    pd.add_argument("--platform", default="ont")
    pd.add_argument("--out", required=True)
    argv = sys.argv[1:]
    rest = []
    if "--" in argv:
        i = argv.index("--")
        argv, rest = argv[:i], argv[i + 1:]
    a = p.parse_args(argv)
    if a.mode == "pipeline":
        pipeline(a, rest)
    else:
        direct(a)


if __name__ == "__main__":
    main()
