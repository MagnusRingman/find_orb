#!/bin/bash
# Side-by-side: find_orb's usual IOD versus the Keplerian-integral linkage.
#
#   usage: ki_sxs.sh <astrometry-file> [extra fo args...]
#          ki_sxs.sh <astrometry-file> --noise      (baseline against itself)
#
# Runs 'fo' twice on the same input -- once as-is, once with KEPLERIAN_LINK=1 --
# and diffs the one-line orbit summaries.
#
# NOTE on interpreting differences: initial_orbit() seeds statistical ranging
# from rand(), so 'fo' is NOT deterministic run to run.  On TestObs.txt, object
# 2007 UQ160 flips between a=5.449,e=0.709 and a=2.989,e=0.204 across two plain
# baseline runs.  Establish that noise floor with --noise before attributing any
# difference to the linkage.
#
# '-k' skips the partial-file unlinking, which asserts at fo.cpp:271 in this
# tree for reasons unrelated to any of this.
#
# Tuning, passed to the KEPLERIAN_LINK run only, via KI_ARGS:
#   KEPLERIAN_LINK_GAP       tracklet break, days (default 0.3)
#   KEPLERIAN_LINK_SPAN      longest baseline, days (default 1000)
#   KEPLERIAN_LINK_RHO_MAX   largest range searched, AU (default 60)
#   KEPLERIAN_LINK_MAX_CHI2  acceptance threshold (default 25)
# e.g. KI_ARGS="KEPLERIAN_LINK_SPAN=600 KEPLERIAN_LINK_MAX_CHI2=10" ki_sxs.sh o.txt

FO=${FO:-./fo}
OUT=${OUT:-$(mktemp -d)}
FILE=$1
shift

if [ -z "$FILE" ]; then
    sed -n '2,25p' "$0"
    exit 1
fi

run() {   # $1 = output file; rest = extra fo args
    local out=$1; shift
    "$FO" "$FILE" -k "$@" 2>/dev/null \
        | sed 's/\x1b\[[0-9;]*m//g' \
        | grep -E "^[0-9]+: " > "$out"
}

if [ "$1" = "--noise" ]; then          # baseline against itself
    shift
    run "$OUT/a.txt" "$@"
    run "$OUT/b.txt" "$@"
    echo "baseline vs baseline (this is the noise floor):"
    diff "$OUT/a.txt" "$OUT/b.txt" | grep -E "^[<>]" || echo "  identical"
    exit 0
fi

run "$OUT/base.txt" "$@"
run "$OUT/ki.txt" KEPLERIAN_LINK=1 $KI_ARGS "$@"

echo "objects: baseline $(wc -l < "$OUT/base.txt"), keplerian $(wc -l < "$OUT/ki.txt")"
if diff -q "$OUT/base.txt" "$OUT/ki.txt" > /dev/null; then
    echo "identical"
else
    echo
    echo "--- differing (< baseline, > KEPLERIAN_LINK=1) ---"
    diff "$OUT/base.txt" "$OUT/ki.txt" | grep -E "^[<>]"
fi
echo
echo "full output in $OUT/base.txt and $OUT/ki.txt"
echo
echo "For firing statistics: fo <file> -k -d1 KEPLERIAN_LINK=1, then grep"
echo "~/.find_orb/debug.txt for 'Keplerian link'."
