"""Reference full-alignment model inputs and probabilities, from Clair3's own Python code.

Runs one `CallVariantsFromCffi` full-alignment job unchanged (Clair3's own argument parser,
`Run`, `call_variants_from_cffi`, `CreateTensorFullAlignment`, `Clair3_F` and `batch_output`),
taking its arguments from a GNU parallel job log such as `log/parallel_6_call_var_bam_full_alignment.log`.
Only `--call_fn` is redirected to PREFIX.vcf. `_torch_predict` and `batch_output` are wrapped to
record what the job feeds the model and gets back, batch by batch:

    PREFIX.inputs.npy   int8 [n, 89, 33, 8], the model inputs in candidate order
    PREFIX.probs.npy    float32 [n, 90], the probabilities `_torch_predict` returned
    PREFIX.batches.tsv  one line per batch: start, size, seconds spent in `_torch_predict`
    PREFIX.info.tsv     position and alt info per candidate, as `batch_output` received them
    PREFIX.vcf          the job's VCF

Run it with the Python of a Clair3 env:

    python tools/oracle_full_alignment.py --clair3-dir $ENV/bin \
        --log m4-oracle/log/parallel_6_call_var_bam_full_alignment.log --chunk chr20.8_13 --out out/c8
"""
import argparse
import shlex
import sys
import time


def job_args(log, chunk):
    """The CallVariantsFromCffi arguments of the job whose `--full_aln_regions` ends in `chunk`."""
    with open(log) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        command = header.index("Command")
        for line in fh:
            argv = shlex.split(line.rstrip("\n").split("\t", command)[command])
            if argv[argv.index("--full_aln_regions") + 1].endswith("/" + chunk):
                return argv[argv.index("CallVariantsFromCffi") + 1:]
    sys.exit(f"no job for chunk {chunk} in {log}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--clair3-dir", required=True, help="directory holding clair3.py")
    p.add_argument("--log", required=True, help="GNU parallel joblog of the full-alignment jobs")
    p.add_argument("--chunk", required=True, help="candidate BED name, e.g. chr20.8_13")
    p.add_argument("--out", required=True, help="output prefix")
    a = p.parse_args()

    argv = job_args(a.log, a.chunk)
    argv[argv.index("--call_fn") + 1] = a.out + ".vcf"
    sys.path.insert(0, a.clair3_dir)

    import numpy as np
    import clair3.CallVariantsFromCffi as cv

    inputs, probs, batches, info = [], [], [], []
    torch_predict, batch_output = cv._torch_predict, cv.batch_output

    def recording_predict(model, device, X):
        start = time.perf_counter()
        Y = torch_predict(model, device, X)
        batches.append((sum(len(x) for x in inputs), len(X), time.perf_counter() - start))
        inputs.append(X.copy())
        probs.append(Y.copy())
        return Y

    def recording_output(position, alt_info_list, Y, *rest):
        info.extend(zip(position, alt_info_list))
        return batch_output(position, alt_info_list, Y, *rest)

    cv._torch_predict, cv.batch_output = recording_predict, recording_output
    sys.argv = ["CallVariantsFromCffi", *argv]
    cv.main()

    X = np.concatenate(inputs) if inputs else np.zeros((0, 89, 33, 8), dtype=np.int8)
    Y = np.concatenate(probs) if probs else np.zeros((0, 90), dtype=np.float32)
    assert X.dtype == np.int8 and Y.dtype == np.float32, (X.dtype, Y.dtype)
    np.save(a.out + ".inputs.npy", X)
    np.save(a.out + ".probs.npy", Y)
    with open(a.out + ".batches.tsv", "w") as fh:
        for start, size, seconds in batches:
            fh.write(f"{start}\t{size}\t{seconds:.6f}\n")
    with open(a.out + ".info.tsv", "w") as fh:
        for position, alt in info:
            fh.write(f"{position}\t{alt}\n")
    total = sum(s for _, _, s in batches)
    print(f"{a.chunk}: {len(X)} candidates, {len(batches)} batches, {total:.3f} s in _torch_predict",
          file=sys.stderr)


if __name__ == "__main__":
    main()
