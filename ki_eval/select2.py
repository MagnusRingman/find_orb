"""Stratified selection, using the observations actually present in UnnObs.txt
for the strata and MPCORB only for the reference orbit.

Strata, all defined on what the observation file really contains:

  H1  2 apparitions, max gap >= 4yr, nobs <= 20    hardest; Link2 territory
  H2  2 apparitions, max gap 1-4yr
  H3  3 apparitions, max gap >= 4yr
  H4  3-4 apparitions, max gap 1-4yr
  C   >=5 apparitions, nobs >= 60                  control; should be easy

H1..H3 are split by inclination -- lo <3deg, mid 3-15deg, hi >=15deg -- because
the geometry these methods rest on (D_i = q_i x u_i) degenerates as the line of
sight approaches the ecliptic, and that is the axis on which I expect them to
fail first.

usage: select2.py <n_per_stratum> <seed> > sample.tsv
"""
import random
import sys

IDX = '/home/bmr/.claude/jobs/1a505dd0/tmp/unnobs_index.tsv'
ORB = '/home/bmr/.claude/jobs/1a505dd0/tmp/unn.tsv'

n_per = int(sys.argv[1]) if len(sys.argv) > 1 else 20
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20260829

orb = {}
with open(ORB) as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 12:
            orb[p[0]] = p          # desig H a e i nobs nopp arc0 arc1 span rms name


def inc_band(i):
    return 'lo' if i < 3 else ('mid' if i < 15 else 'hi')


def stratum(nobs, gap, nopp):
    yr = gap / 365.25
    if nopp == 2 and yr >= 4 and nobs <= 20:
        return 'H1'
    if nopp == 2 and 1 <= yr < 4:
        return 'H2'
    if nopp == 3 and yr >= 4:
        return 'H3'
    if 3 <= nopp <= 4 and 1 <= yr < 4:
        return 'H4'
    if nopp >= 5 and nobs >= 60:
        return 'C'
    return None


buckets = {}
info = {}
with open(IDX) as f:
    f.readline()
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 7:
            continue
        d = p[0]
        o = orb.get(d)
        if not o:            # no reference orbit -> nothing to score against
            continue
        nobs, arc, gap = int(p[1]), float(p[2]), float(p[3])
        ntrk, nopp = int(p[4]), int(p[5])
        if nobs < 6:         # too few to fit anything at all
            continue
        s = stratum(nobs, gap, nopp)
        if not s:
            continue
        key = s + '-' + inc_band(float(o[4])) if s in ('H1', 'H2', 'H3') else s
        buckets.setdefault(key, []).append(d)
        info[d] = (p, o)

rng = random.Random(seed)
print("desig\tstratum\tnobs\tarc_days\tmax_gap\tn_trk\tn_opp"
      "\tref_a\tref_e\tref_i\tname")
for key in sorted(buckets):
    pool = buckets[key]
    rng.shuffle(pool)
    for d in pool[:n_per]:
        p, o = info[d]
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
              % (d, key, p[1], p[2], p[3], p[4], p[5], o[2], o[3], o[4], o[11]))
    print("# %-8s pool %7d, took %d" % (key, len(pool), min(n_per, len(pool))),
          file=sys.stderr)
