"""One pass over UnnObs.txt, writing <outdir>/<desig>.obs for every designation
in the first column of a sample TSV.

usage: extract2.py <sample.tsv> <outdir>
"""
import os
import sys
import time

SRC = '/home/bmr/src/find_orb/data/UnnObs.txt'

sample, outdir = sys.argv[1], sys.argv[2]
wanted = set()
with open(sample) as f:
    f.readline()
    for line in f:
        if line.startswith('#'):
            continue
        d = line.split('\t')[0].strip()
        if d:
            wanted.add(d)

os.makedirs(outdir, exist_ok=True)
found = {d: [] for d in wanted}
t0 = time.time()
with open(SRC, 'r', errors='replace') as f:
    for line in f:
        if line[5:12] in found:
            found[line[5:12]].append(line.rstrip('\n'))

n = 0
for d, lines in found.items():
    if lines:
        with open(os.path.join(outdir, d + '.obs'), 'w') as g:
            g.write("\n".join(lines) + "\n")
        n += 1
print("extracted %d/%d objects in %.1f s -> %s"
      % (n, len(wanted), time.time() - t0, outdir), file=sys.stderr)
