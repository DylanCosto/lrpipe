"""Reference outputs for `lrpipe-snv pileup`, from Clair3's own Python code.

Runs `CreateTensorPileup` for one chunk the way `CallVariantsFromCffi --pileup` does, then the
batching, high-depth scaling and `_torch_predict` steps of `call_variants_from_cffi`, and writes
PREFIX.tensors.npy, PREFIX.inputs.npy, PREFIX.info.tsv and PREFIX.probs.npy in the same layout
as the Rust tool. Run it with the Python of a Clair3 env:

    python tools/oracle_pileup.py --clair3-dir $ENV/bin --bam ... --ref ... --ctg-name chr20 \
        --chunk-id 5 --chunk-num 13 --model $ENV/bin/models/.../pileup.pt --out oracle/c5
"""
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--clair3-dir", required=True, help="directory holding clair3.py")
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--ctg-name", required=True)
    p.add_argument("--chunk-id", type=int)
    p.add_argument("--chunk-num", type=int)
    p.add_argument("--ctg-start", type=int)
    p.add_argument("--ctg-end", type=int)
    p.add_argument("--extend-bed")
    p.add_argument("--snp-min-af", type=float, default=0.08)
    p.add_argument("--indel-min-af", type=float, default=0.15)
    p.add_argument("--min-coverage", type=int, default=2)
    p.add_argument("--min-mq", type=int, default=5)
    p.add_argument("--model")
    p.add_argument("--batch-size", type=int, default=200)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    # As in CallVariantsFromCffi.Run, before torch starts any threads.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = str(a.threads)
    sys.path.insert(0, a.clair3_dir)

    import numpy as np
    import torch

    torch.set_num_threads(a.threads)
    torch.set_num_interop_threads(2 if a.threads > 1 else 1)

    import shared.param_p as param
    from preprocess.CreateTensorPileupFromCffi import CreateTensorPileup

    args = argparse.Namespace(
        platform="ont", bam_fn=a.bam, ref_fn=a.ref, tensor_can_fn="PIPE", vcf_fn=None,
        min_af=0.08, snp_min_af=a.snp_min_af, indel_min_af=a.indel_min_af, ctgName=a.ctg_name,
        ctgStart=a.ctg_start, ctgEnd=a.ctg_end, bed_fn=None, gvcf=False, sampleName="SAMPLE",
        samtools="samtools", fast_mode=False, enable_variant_calling_at_sequence_head_and_tail=False,
        minCoverage=a.min_coverage, minMQ=a.min_mq, minBQ=param.min_bq, max_depth=param.max_depth,
        call_snp_only=False, extend_bed=a.extend_bed, temp_file_dir="./", indel_fn=None,
        base_err=param.base_err, gq_bin_size=param.gq_bin_size, bp_resolution=False, zstd=param.zstd,
        test_pos=0, chunk_num=a.chunk_num, chunk_id=a.chunk_id, enable_long_indel=False,
    )
    tensor, positions, alt_info = CreateTensorPileup(args)
    tensor = np.asarray(tensor, dtype=np.int32).reshape(-1, param.no_of_positions, param.channel_size)
    np.save(a.out + ".tensors.npy", tensor)
    with open(a.out + ".info.tsv", "w") as fh:
        for pos_info, alt in zip(positions, alt_info):
            fh.write(f"{pos_info}\t{alt}\n")

    # call_variants_from_cffi, CPU path: batches of predictBatchSize * threads, then scaling.
    inputs = tensor.copy()
    total = len(inputs)
    probs = []
    if a.model:
        from clair3.CallVariantsFromCffi import _load_torch_checkpoint, _torch_predict
        from clair3.model import Clair3_P

        device = torch.device("cpu")
        m = Clair3_P(add_indel_length=False, predict=True, input_channels=param.ont_input_shape[-1])
        m.to(device)
        m.eval()
        _load_torch_checkpoint(m, a.model, device)
    for start in range(0, total, a.batch_size):
        X = inputs[start:start + a.batch_size]
        for alt_idx, alt in enumerate(alt_info[start:start + a.batch_size]):
            depth = int(alt.split('-', maxsplit=1)[0])
            max_depth = param.max_depth_dict["ont"]
            if depth > 0 and depth > max_depth * 1.5:
                scale_factor = depth / max_depth
                X[alt_idx] = X[alt_idx] / scale_factor
        if a.model:
            probs.append(_torch_predict(m, device, X))
    np.save(a.out + ".inputs.npy", inputs)
    if a.model:
        out = np.concatenate(probs) if probs else np.zeros((0, 24), dtype=np.float32)
        np.save(a.out + ".probs.npy", out)
    print(f"{total} candidates", file=sys.stderr)


if __name__ == "__main__":
    main()
