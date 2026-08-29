"""Run fo over a sample, several times per object, with and without the
Keplerian-integral linkage, and record what it found.

The headline metric is arc recovery, not element agreement.  On 1999 EO3 --
33 observations in 1999 and 283 in 2020-21 -- fo drops the 1999 apparition and
fits the recent one alone, yet still lands within 0.06% of MPCORB in a, because
the recent apparition is by itself well observed.  Element agreement therefore
hides exactly the failure we are hunting.  Whether the distant apparition got
linked is visible in the fitted arc.

Each run gets its own working directory and its own config directory, so that
workers cannot tread on each other's elements.txt or debug.txt.

usage: runner.py <sample.tsv> <obsdir> <outdir> <n_runs> <n_workers> [modes]
       modes: comma-separated subset of base,ki  (default both)
"""
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

FO = os.environ.get('FO_BIN', '/home/bmr/src/find_orb/.claude/worktrees/gronchi/fo')

sample, obsdir, outdir, n_runs, n_workers = sys.argv[1:6]
n_runs, n_workers = int(n_runs), int(n_workers)
modes = sys.argv[6].split(',') if len(sys.argv) > 6 else ['base', 'ki']
KI_EXTRA = os.environ.get('KI_ARGS', '').split()


def parse_sof(path):
    """Fixed-width parse driven by the '|' positions in the sof header."""
    try:
        with open(path) as f:
            header = f.readline().rstrip('\n')
            row = f.readline().rstrip('\n')
    except (IOError, IndexError):
        return None
    if not row.strip():
        return None
    bounds, start = [], 0
    for i, c in enumerate(header):
        if c == '|':
            bounds.append((start, i))
            start = i + 1
    bounds.append((start, len(header) + 40))
    names = [header[a:b].strip().rstrip('.^ ').strip() for a, b in bounds]
    vals = {}
    for (a, b), nm in zip(bounds, names):
        vals[nm] = row[a:b].strip()
    return vals


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def sof_date_to_mjd(s):
    """sof dates look like 20200821.60579 (YYYYMMDD.fff)."""
    v = fnum(s)
    if v is None:
        return None
    d = int(v)
    frac = v - d
    y, m, dd = d // 10000, (d // 100) % 100, d % 100
    a = (14 - m) // 12
    yy = y + 4800 - a
    mm = m + 12 * a - 3
    jdn = (dd + (153 * mm + 2) // 5 + 365 * yy + yy // 4
           - yy // 100 + yy // 400 - 32045)
    return jdn - 2400000.5 + frac


def one(job):
    desig, mode, run_idx = job
    tag = '%s_%s_%d' % (desig, mode, run_idx)
    work = os.path.join(outdir, 'w', tag)
    # A config directory per job, not per worker.  fo persists settings into
    # environ.dat in whatever config directory it is given, so a shared one
    # lets a 'ki' run leave KEPLERIAN_LINK=1 behind for a later 'base' run in
    # the same directory -- which silently turns the control into a treatment.
    cfg = os.path.join(outdir, 'w', tag + '.cfg')
    os.makedirs(work, exist_ok=True)
    os.makedirs(cfg, exist_ok=True)
    # '-i' stops fo starting from a stored solution instead of doing IOD.
    args = [FO, os.path.join(os.path.abspath(obsdir), desig + '.obs'),
            '-k', '-i', '-x', os.path.abspath(cfg) + '/']
    # Always state the setting explicitly.  fo persists settings into
    # environ.dat, and reads come from ~/.find_orb whatever -x says, so a
    # single stray interactive run with KEPLERIAN_LINK=1 silently enables it
    # for every later run -- including the controls.  Never rely on absence.
    if mode == 'ki':
        args.append('KEPLERIAN_LINK=1')
        args += KI_EXTRA
    else:
        args.append('KEPLERIAN_LINK=0')
    t0 = time.time()
    try:
        subprocess.run(args, cwd=work, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=600)
    except subprocess.TimeoutExpired:
        shutil.rmtree(work, ignore_errors=True)
        return (desig, mode, run_idx, 'timeout', None)
    dt = time.time() - t0
    v = parse_sof(os.path.join(work, 'sof.txt'))
    shutil.rmtree(work, ignore_errors=True)
    if not v:
        return (desig, mode, run_idx, 'nofit', dt)
    q, e = fnum(v.get('q')), fnum(v.get('e'))
    a = q / (1. - e) if (q is not None and e is not None and e < 1.) else None
    t0m, t1m = sof_date_to_mjd(v.get('Tfirst')), sof_date_to_mjd(v.get('Tlast'))
    return (desig, mode, run_idx, 'ok', dt, a, e, fnum(v.get('i')),
            fnum(v.get('rms')), fnum(v.get('n_u')), fnum(v.get('n_o')),
            t0m, t1m)


def main():
    objs = []
    with open(sample) as f:
        f.readline()
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if p and p[0] and os.path.exists(os.path.join(obsdir, p[0] + '.obs')):
                objs.append(p[0])
    jobs = [(d, m, r) for d in objs for m in modes for r in range(n_runs)]
    os.makedirs(outdir, exist_ok=True)
    out = open(os.path.join(outdir, 'runs.tsv'), 'w')
    out.write("desig\tmode\trun\tstatus\tsecs\ta\te\ti\trms\tn_used\tn_obs"
              "\tfit_t0\tfit_t1\n")
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for r in ex.map(one, jobs, chunksize=1):
            done += 1
            row = list(r) + [''] * (13 - len(r))
            out.write("\t".join('' if x is None else
                                (('%.6f' % x) if isinstance(x, float) else str(x))
                                for x in row) + "\n")
            if done % 50 == 0:
                out.flush()
                print("  %d/%d  %.0fs elapsed" % (done, len(jobs), time.time() - t0),
                      file=sys.stderr)
    out.close()
    shutil.rmtree(os.path.join(outdir, 'w'), ignore_errors=True)
    print("%d jobs (%d objects x %d modes x %d runs) in %.0f s"
          % (len(jobs), len(objs), len(modes), n_runs, time.time() - t0),
          file=sys.stderr)


if __name__ == '__main__':
    main()
