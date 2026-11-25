"""Digest of every Clair3 checkpoint as PyTorch reads it, for crates/snv-infer/tests/checkpoints.rs.

usage: checkpoint_digest.py MODELS_DIR > digests.txt
"""
import sys, glob, hashlib, torch
for f in sorted(glob.glob(sys.argv[1] + "/*/*.pt")):
    sd = torch.load(f, map_location="cpu"); sd = sd.get("state_dict", sd)
    h = hashlib.sha256()
    for k, v in sd.items():
        h.update(k.encode()); h.update(str(tuple(v.shape)).encode()); h.update(v.contiguous().numpy().tobytes())
    print(f.split("/models/")[1], len(sd), h.hexdigest()[:16])
