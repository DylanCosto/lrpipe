"""PyTorch counterpart of `examples/score_npy.rs`: Clair3's own model classes and `_torch_predict`.

usage: score_npy_torch.py CLAIR3_BIN pileup|full_alignment MODEL.pt INPUTS.npy BATCH DEVICE OUT.npy
Run with the Python of the libtorch being compared (e.g. a CUDA build of PyTorch 2.13.0).
"""
import sys
import time

clair3_bin, kind, model_fn, inputs_fn, batch, device_name, out_fn = sys.argv[1:8]
sys.path.insert(0, clair3_bin)
import numpy as np
import torch

torch.set_num_threads(1)
torch.set_num_interop_threads(1)
from clair3.CallVariantsFromCffi import _load_torch_checkpoint, _torch_predict

device = torch.device(device_name)
if kind == "pileup":
    import shared.param_p as param
    from clair3.model import Clair3_P as Model
    m = Model(add_indel_length=False, predict=True, input_channels=param.ont_input_shape[-1])
else:
    import shared.param_f as param
    from clair3.model import Clair3_F as Model
    m = Model(add_indel_length=True, predict=True, input_channels=param.ont_input_shape[-1])
m.to(device)
m.eval()
_load_torch_checkpoint(m, model_fn, device)
X = np.load(inputs_fn)
batch = int(batch)
start = time.time()
probs = [_torch_predict(m, device, X[i:i + batch]) for i in range(0, len(X), batch)]
if device.type == "cuda":
    torch.cuda.synchronize()
print(f"{kind} on {device}: {len(X)} candidates in {time.time() - start:.2f}s", file=sys.stderr)
np.save(out_fn, np.concatenate(probs))
