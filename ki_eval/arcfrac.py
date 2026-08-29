"""Mean and median arc fraction recovered, per stratum, base against ki.

The three-bin summary (full/part/none) is too blunt to see a trade where one
mode converts catastrophic failures into partial recoveries while giving up a
little at the top.  The fraction of the available arc actually fitted is the
continuous version of the same quantity.

usage: arcfrac.py <sample.tsv> <runs.tsv> [runs2.tsv ...]
"""
import collections
import statistics
import sys

sample = sys.argv[1]
meta = {}
for line in open(sample).readlines()[1:]:
    if line.startswith('#'):
        continue
    p = line.rstrip('\n').split('\t')
    if len(p) >= 11:
        meta[p[0]] = dict(s=p[1], arc=float(p[3]), i=float(p[9]))

for runs in sys.argv[2:]:
    R = collections.defaultdict(list)
    for line in open(runs).readlines()[1:]:
        p = line.rstrip('\n').split('\t')
        if len(p) < 13:
            continue
        try:
            R[(p[0], p[1])].append((p[3], float(p[11]), float(p[12])))
        except ValueError:
            R[(p[0], p[1])].append((p[3], None, None))

    def frac(d, mode):
        rs = [r for r in R.get((d, mode), []) if r[0] == 'ok' and r[1] is not None]
        if not rs or meta[d]['arc'] <= 0:
            return None
        return max((t1 - t0) for _, t0, t1 in rs) / meta[d]['arc']

    modes = sorted({k[1] for k in R})
    print("=== %s ===" % runs)
    print("%-8s %5s | %s" % ('', '', ' | '.join("%-22s" % ('%s  mean/med' % m)
                                                for m in modes)))
    tot = {m: [] for m in modes}
    for s in sorted({m['s'] for m in meta.values()}):
        objs = [d for d in meta if meta[d]['s'] == s]
        row = []
        for mode in modes:
            fs = [frac(d, mode) for d in objs]
            fs = [f for f in fs if f is not None]
            tot[mode] += fs
            row.append("%6.3f / %6.3f  n=%2d" % (statistics.mean(fs),
                                                 statistics.median(fs), len(fs))
                       if fs else "        n/a")
        print("%-8s %5d | %s" % (s, len(objs), ' | '.join(row)))
    print()
    for mode in modes:
        fs = tot[mode]
        if not fs:
            continue
        print("TOTAL %-4s n=%d  mean=%.4f median=%.4f  >=0.9: %d  <0.1: %d"
              % (mode, len(fs), statistics.mean(fs), statistics.median(fs),
                 sum(1 for f in fs if f >= .9), sum(1 for f in fs if f < .1)))
    print()
