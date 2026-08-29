"""Parse the unnumbered section of MPCORB.DAT into a compact TSV.

MPCORB column layout (1-indexed):
   1-7    packed designation
   9-13   H
  21-25   packed epoch
  27-35   M          38-46  argument of perihelion
  49-57   node       60-68  inclination
  71-79   eccentricity      81-91  mean daily motion
  93-103  semimajor axis
 118-122  number of observations
 124-126  number of oppositions
 128-136  arc: "1982-2026" for multi-opposition, "nnn days" for single
 138-141  rms residual
 167-194  readable designation
"""
import sys

SRC = '/home/bmr/src/find_orb/data/MPCORB.DAT'
OUT = '/home/bmr/.claude/jobs/1a505dd0/tmp/unn.tsv'


def f(line, a, b):
    return line[a - 1:b].strip()


def num(s, default=None):
    try:
        return float(s)
    except ValueError:
        return default


n_read = n_kept = 0
with open(SRC, 'r', errors='replace') as fin, open(OUT, 'w') as fout:
    fout.write("desig\tH\ta\te\ti\tnobs\tnopp\tarc0\tarc1\tarcspan\trms\tname\n")
    for line in fin:
        if len(line) < 194:
            continue
        desig = f(line, 1, 7)
        # packed provisional designations start with I/J/K then two digits
        if not (len(desig) == 7 and desig[0] in 'IJK'
                and desig[1:3].isdigit() and desig[3].isalpha()):
            continue
        n_read += 1
        a = num(f(line, 93, 103))
        e = num(f(line, 71, 79))
        inc = num(f(line, 60, 68))
        H = num(f(line, 9, 13))
        nobs = num(f(line, 118, 122))
        nopp = num(f(line, 124, 126))
        arc = f(line, 128, 136)
        rms = num(f(line, 138, 141))
        name = f(line, 167, 194)
        if None in (a, e, inc, nobs, nopp):
            continue
        arc0 = arc1 = span = ''
        if '-' in arc and 'day' not in arc:
            parts = arc.split('-')
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                arc0, arc1 = parts
                span = str(int(parts[1]) - int(parts[0]))
        elif 'day' in arc:
            arc0 = arc1 = ''
            span = '0'
        fout.write("%s\t%s\t%.4f\t%.4f\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n"
                   % (desig, H if H is not None else '', a, e, inc,
                      int(nobs), int(nopp), arc0, arc1, span,
                      rms if rms is not None else '', name))
        n_kept += 1

print("unnumbered entries read %d, written %d" % (n_read, n_kept), file=sys.stderr)
