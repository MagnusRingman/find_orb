"""Does recovering more arc actually produce the right orbit?

Arc recovery on its own is not enough.  A linkage that joins two apparitions
that do not belong together would also 'recover the arc', and would look like
a win on that metric alone while being flatly wrong.  So for every object we
compare the fitted elements against the MPCORB reference, and report agreement
separately for the objects where each mode recovered the full arc.

Note that the reference is not beyond question either: the MPC are
conservative about linkages but do occasionally publish one that turns out to
be wrong, so a handful of objects on which everything disagrees with MPCORB
may be bad linkages rather than bad fits.

usage: elemcheck.py <sample.tsv> <runs.tsv>
"""
import collections
import statistics
import sys

sample, runs = sys.argv[1], sys.argv[2]

meta = {}
for line in open(sample).readlines()[1:]:
    if line.startswith('#'):
        continue
    p = line.rstrip('\n').split('\t')
    if len(p) >= 11:
        meta[p[0]] = dict(s=p[1], arc=float(p[3]), a=float(p[7]),
                          e=float(p[8]), i=float(p[9]))

R = collections.defaultdict(list)
for line in open(runs).readlines()[1:]:
    p = line.rstrip('\n').split('\t')
    if len(p) < 13 or p[3] != 'ok':
        continue
    try:
        R[(p[0], p[1])].append(dict(a=float(p[5]), e=float(p[6]), i=float(p[7]),
                                    t0=float(p[11]), t1=float(p[12])))
    except ValueError:
        pass

modes = sorted({k[1] for k in R})


def best(d, mode):
    rs = R.get((d, mode), [])
    return max(rs, key=lambda r: r['t1'] - r['t0']) if rs else None


def agrees(r, m):
    return (abs(r['a'] - m['a']) / m['a'] < .01 and abs(r['e'] - m['e']) < .02
            and abs(r['i'] - m['i']) < 1.)


print("%-5s | %6s %6s %6s | %s" % ('mode', 'n', 'full', 'agree',
                                   'of the full-arc fits: agree / disagree'))
for mode in modes:
    n = full = agree = full_agree = 0
    da = []
    for d, m in meta.items():
        r = best(d, mode)
        if not r or m['arc'] <= 0:
            continue
        n += 1
        ok = agrees(r, m)
        agree += ok
        if (r['t1'] - r['t0']) / m['arc'] >= .9:
            full += 1
            full_agree += ok
            da.append(abs(r['a'] - m['a']) / m['a'])
    print("%-5s | %6d %6d %6d | %d agree, %d disagree   median |da|/a = %.2e"
          % (mode, n, full, agree, full_agree, full - full_agree,
             statistics.median(da) if da else float('nan')))

if len(modes) > 1:
    print("\nobjects where ki recovered the full arc and base did not:")
    gained = disagreed = 0
    examples = []
    for d, m in meta.items():
        rb, rk = best(d, 'base'), best(d, 'ki')
        if not rb or not rk or m['arc'] <= 0:
            continue
        fb = (rb['t1'] - rb['t0']) / m['arc']
        fk = (rk['t1'] - rk['t0']) / m['arc']
        if fk >= .9 and fb < .9:
            gained += 1
            if not agrees(rk, m):
                disagreed += 1
                if len(examples) < 10:
                    examples.append((d, m, rk))
    print("  %d objects gained;  of those, %d disagree with MPCORB"
          % (gained, disagreed))
    for d, m, r in examples:
        print("    %-8s %-7s ref a=%7.4f e=%6.4f i=%6.2f "
              " got a=%7.4f e=%6.4f i=%6.2f"
              % (d, m['s'], m['a'], m['e'], m['i'], r['a'], r['e'], r['i']))

    print("\nobjects where base recovered the full arc and ki did not:")
    lost = 0
    for d, m in meta.items():
        rb, rk = best(d, 'base'), best(d, 'ki')
        if not rb or not rk or m['arc'] <= 0:
            continue
        if (rb['t1'] - rb['t0']) / m['arc'] >= .9 > (rk['t1'] - rk['t0']) / m['arc']:
            lost += 1
    print("  %d objects lost" % lost)
