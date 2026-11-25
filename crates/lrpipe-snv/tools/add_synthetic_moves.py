"""Adds a deterministic synthetic Dorado move table (mv:B:c) to SAM records, for dwell-time
parity tests: the calls mean nothing, but Clair3 and lrpipe-snv must read the tags the same way.

usage: samtools view -h IN.bam REGION | python3 add_synthetic_moves.py [VARIANT] | samtools view -b -o OUT.bam -

Each base gets 1-9 signal blocks (1-25 with VARIANT 2), one base in 61 a stall of 100-259
blocks (past int8), some reads two blocks before their first base, and one read in 17 (7 with
VARIANT 2) no tag at all. The values come from a hash of the read name.
"""
import sys
import zlib

variant = int(sys.argv[1]) if len(sys.argv) > 1 else 1
max_blocks, untagged, stride = (9, 17, '5') if variant == 1 else (25, 7, '6')
for line in sys.stdin:
    if line.startswith('@'):
        sys.stdout.write(line)
        continue
    fields = line.rstrip('\n').split('\t')
    seq = fields[9]
    h = zlib.crc32(fields[0].encode()) or 1
    if seq == '*' or h % untagged == 0:
        sys.stdout.write(line)
        continue
    moves = [stride]
    x = h
    for _ in seq:
        x = (x * 1103515245 + 12345) & 0x7fffffff
        blocks = 1 + (x >> 8) % max_blocks
        if (x >> 4) % 61 == 0:
            blocks = 100 + (x >> 12) % 160
        moves.append('1')
        moves.extend('0' * (blocks - 1))
    if h % 5 == 0:
        moves[1:1] = ['0', '0']
    fields.append('mv:B:c,' + ','.join(moves))
    sys.stdout.write('\t'.join(fields) + '\n')
