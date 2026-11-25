"""Reference inputs for the full-alignment VCF output, from Clair3's own code.

Runs one `CallVariantsFromCffi` full-alignment job exactly as run_clair3 did (the command is
taken from the last column of `parallel_6_call_var_bam_full_alignment.log`): Clair3's
`CreateTensorFullAlignment`, its model, its batching (`predictBatchSize * threads`) and its
`batch_output`. `batch_output` is wrapped to record what it is given, and the job writes its VCF
to PREFIX.clair3.vcf instead of the original `--call_fn`. Then writes

    PREFIX.info.tsv     position string and alt-info string of every candidate, in order
    PREFIX.probs.npy    float32 probabilities, one row per candidate (21 + 3 + 33 + 33)
    PREFIX.replay.vcf   `get_header` + `batch_output` run again on the saved strings and rows
    PREFIX.job.tsv      the job's original VCF, --ref_fn, --cmd_fn and --sampleName (one line)

PREFIX.clair3.vcf and PREFIX.replay.vcf must both equal the job's original VCF. Run it with the
Python of a Clair3 env (its `bin` on PATH, for pigz):

    python tools/oracle_fa_output.py --clair3-dir $ENV/bin \
        --log $OUT/log/parallel_6_call_var_bam_full_alignment.log --chunk 7_13 --out work/c7
"""
import argparse
import shlex
import sys


def job_command(log, chunk):
    """The CallVariantsFromCffi arguments of the job whose --call_fn ends in `_{chunk}.vcf`."""
    with open(log) as fh:
        next(fh)
        for row in fh:
            cmd = shlex.split(row.rstrip("\n").split("\t")[-1])
            call_fn = cmd[cmd.index("--call_fn") + 1]
            if call_fn.endswith(f".{chunk}.vcf"):
                return cmd[cmd.index("CallVariantsFromCffi") + 1:], call_fn
    sys.exit(f"no job for chunk {chunk} in {log}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--clair3-dir", required=True, help="directory holding clair3.py")
    p.add_argument("--log", required=True, help="parallel_6_call_var_bam_full_alignment.log")
    p.add_argument("--chunk", required=True, help="chunk suffix of the job's --call_fn, e.g. 7_13")
    p.add_argument("--out", required=True, help="output prefix")
    a = p.parse_args()

    argv, original_vcf = job_command(a.log, a.chunk)
    argv[argv.index("--call_fn") + 1] = a.out + ".clair3.vcf"
    print(f"job: {original_vcf}", file=sys.stderr)

    sys.path.insert(0, a.clair3_dir)
    import numpy as np
    import clair3.CallVariantsFromCffi as cffi
    from shared.utils import get_header

    positions, alt_infos, rows, configs = [], [], [], []
    original_batch_output = cffi.batch_output

    def recording_batch_output(position, alt_info_list, Y, output_config, output_utilities, args=None):
        positions.extend(position)
        alt_infos.extend(alt_info_list)
        rows.append(np.array(Y, copy=True))
        configs.append(output_config)
        return original_batch_output(position, alt_info_list, Y, output_config, output_utilities, args)

    # call_variants_from_cffi looks `batch_output` up in its module's globals.
    cffi.batch_output = recording_batch_output
    sys.argv = ["CallVariantsFromCffi"] + argv
    cffi.main()

    probs = np.concatenate(rows) if rows else np.zeros((0, 90), dtype=np.float32)
    assert probs.dtype == np.float32 and len(probs) == len(positions) == len(alt_infos)
    np.save(a.out + ".probs.npy", probs)
    with open(a.out + ".info.tsv", "w") as fh:
        for position, alt_info in zip(positions, alt_infos):
            fh.write(f"{position}\t{alt_info}\n")

    def opt(name):
        return argv[argv.index(name) + 1]

    with open(a.out + ".job.tsv", "w") as fh:
        fh.write("\t".join([original_vcf, opt("--ref_fn"), opt("--cmd_fn"), opt("--sampleName")]) + "\n")

    # Replay batch_output on what was saved, read back from disk.
    probs = np.load(a.out + ".probs.npy")
    with open(a.out + ".info.tsv") as fh:
        saved = [row.rstrip("\n").split("\t") for row in fh]
    header = get_header(reference_file_path=opt("--ref_fn"), cmd_fn=opt("--cmd_fn"),
                        sample_name=opt("--sampleName")) + "\n"
    body = original_batch_output([s[0] for s in saved], [s[1] for s in saved], probs, configs[0], None) \
        if saved else ""
    with open(a.out + ".replay.vcf", "w") as fh:
        fh.write(header + body)
    print(f"{len(saved)} candidates", file=sys.stderr)


if __name__ == "__main__":
    main()
