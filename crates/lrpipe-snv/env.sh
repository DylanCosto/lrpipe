# Environment for this machine's checks and tools: `source env.sh`.
#
# TORCH_ENV is the original's env (Clair3 2.0.3 with PyTorch 2.13.0): it provides the models, and
# the Python, samtools and C tools that run the original for parity. Building lrpipe-snv needs
# nothing from it.
TORCH_ENV=${TORCH_ENV:-/data2/bioinformatics/envs/clair3-t213}
export PATH="$TORCH_ENV/bin:$PATH"
