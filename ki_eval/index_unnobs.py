"""Index UnnObs.txt: one row per object describing the data that is actually
there.

MPCORB cannot be used for this.  It is generated later than UnnObs.txt and
knows about linkages the observation file does not: 2013 HA62 is listed with
16 observations over 2013-2025, but UnnObs.txt holds only the eight from 2013.
Selecting hard cases on MPCORB's arc would therefore mis-describe the input fo
is actually given.  So we measure the arc, the gaps and the tracklet count from
the observations themselves, and keep MPCORB only as a reference orbit.

The file is sorted by designation, so objects arrive in contiguous blocks and
we never hold more than one object in memory.  We verify that assumption as we
go rather than trusting it.

usage: index_unnobs.py > unnobs_index.tsv
"""
import sys

SRC = '/home/bmr/src/find_orb/data/UnnObs.txt'
TRACKLET_GAP = 0.3          # days; matches KEPLERIAN_LINK_GAP's default


def mjd_of(line):
    try:
        y = int(line[15:19]); m = int(line[20:22]); d = float(line[23:31])
    except ValueError:
        return None
    a = (14 - m) // 12
    yy = y + 4800 - a
    mm = m + 12 * a - 3
    jdn = (int(d) + (153 * mm + 2) // 5 + 365 * yy + yy // 4
           - yy // 100 + yy // 400 - 32045)
    return jdn - 2400000.5 + (d - int(d))


def emit(desig, mjds, out):
    if not mjds:
        return
    mjds.sort()
    arc = mjds[-1] - mjds[0]
    gaps = [mjds[i + 1] - mjds[i] for i in range(len(mjds) - 1)]
    max_gap = max(gaps) if gaps else 0.
    n_trk = 1 + sum(1 for g in gaps if g > TRACKLET_GAP)
    # "oppositions" as seen in the data: gaps beyond 100 days
    n_opp = 1 + sum(1 for g in gaps if g > 100.)
    out.write("%s\t%d\t%.2f\t%.2f\t%d\t%d\t%.4f\n"
              % (desig, len(mjds), arc, max_gap, n_trk, n_opp, mjds[0]))


out = sys.stdout
out.write("desig\tnobs\tarc_days\tmax_gap\tn_tracklets\tn_opp\tfirst_mjd\n")
seen = set()
cur = None
mjds = []
n_lines = n_obj = n_outoforder = 0
with open(SRC, 'r', errors='replace') as f:
    for line in f:
        n_lines += 1
        d = line[5:12]
        if d != cur:
            if cur is not None:
                emit(cur, mjds, out)
                n_obj += 1
                if cur in seen:
                    n_outoforder += 1
                seen.add(cur)
            cur = d
            mjds = []
        t = mjd_of(line)
        if t is not None:
            mjds.append(t)
if cur is not None:
    emit(cur, mjds, out)
    n_obj += 1

print("lines %d, objects %d, designations seen more than once %d"
      % (n_lines, n_obj, n_outoforder), file=sys.stderr)
