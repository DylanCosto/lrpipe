"""Read an uncompressed BAM stream on stdin, skip the header, and print the MD5
of all alignment records and their total size. Two BAM files have
byte-identical records (including tag types) exactly when the hashes are equal,
whatever their headers or BGZF block boundaries.

Usage: bgzip -dc file.bam | python3 tests/record_md5.py
"""
import hashlib
import struct
import sys

f = sys.stdin.buffer


def read(n):
    b = f.read(n)
    if len(b) != n:
        sys.exit("truncated BAM header")
    return b


if read(4) != b"BAM\x01":
    sys.exit("not an uncompressed BAM stream")
(l_text,) = struct.unpack("<i", read(4))
read(l_text)
(n_ref,) = struct.unpack("<i", read(4))
for _ in range(n_ref):
    (l_name,) = struct.unpack("<i", read(4))
    read(l_name + 4)
h = hashlib.md5()
total = 0
while chunk := f.read(1 << 22):
    h.update(chunk)
    total += len(chunk)
print(h.hexdigest(), total)
