#!/usr/bin/env python3
"""Random inputs for tools/ssw_harness.cpp and tools/dbg_harness.cpp.

usage: make_realign_cases.py ssw|dbg N SEED > cases.tsv
"""
import random
import sys

B = 'ACGT'


def mutate(seq, edits):
    x = list(seq)
    for _ in range(edits):
        if not x:
            break
        i = random.randrange(len(x))
        r = random.random()
        if r < 0.4:
            x[i] = random.choice(B)
        elif r < 0.6:
            x.insert(i, ''.join(random.choice(B) for _ in range(random.randint(1, 6))))
        elif r < 0.85:
            del x[i:i + random.randint(1, 6)]
        else:
            x[i] = 'N'
    return ''.join(x)


def ssw_case():
    length = random.choice([20, 40, 80, 150, 250, 400, 700, 1100])
    ref = ''.join(random.choice(B) for _ in range(length))
    if random.random() < 0.1:
        unit = ''.join(random.choice(B) for _ in range(random.randint(1, 6)))
        ref = (unit * length)[:length]
    if random.random() < 0.15:
        query = ''.join(random.choice(B) for _ in range(random.randint(5, 200)))
    else:
        s = random.randint(0, max(0, length - 10))
        query = mutate(ref[s:s + random.randint(10, 200)], random.randint(0, 8))
        if random.random() < 0.2:
            flank = lambda: ''.join(random.choice(B) for _ in range(random.randint(0, 20)))
            query = flank() + query + flank()
    if random.random() < 0.05:
        ref = ref[:length // 2] + 'NNN' + ref[length // 2:]
    return ref + '\t' + (query or 'A')


def dbg_case():
    length = random.choice([30, 60, 120, 200, 400, 800])
    ref = ''.join(random.choice(B) for _ in range(length))
    if random.random() < 0.3:
        unit = ''.join(random.choice(B) for _ in range(random.randint(2, 12)))
        p = random.randrange(length)
        ref = ref[:p] + unit * random.randint(2, 6) + ref[p:]
    haps = [ref] + [mutate(ref, random.randint(1, 4)).replace('N', 'A') for _ in range(random.randint(1, 3))]
    reads, quals = [], []
    for _ in range(random.randint(0, 30)):
        h = random.choice(haps)
        s = random.randint(0, max(0, len(h) - 20))
        read = list(h[s:s + random.randint(20, 150)])
        if read and random.random() < 0.2:
            read[random.randrange(len(read))] = 'N'
        reads.append(''.join(read))
        quals.append(' '.join(str(i) for i in range(len(read)) if random.random() < 0.03))
    return ref + '\t' + ','.join(reads) + '\t' + ','.join(quals)


def main():
    kind, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    random.seed(seed)
    make = ssw_case if kind == 'ssw' else dbg_case
    for _ in range(n):
        print(make())


if __name__ == '__main__':
    main()
