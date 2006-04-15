#! /bin/sh

set -e

# Put ../tools on PATH so that alpha-shape and
# surface-distance-transform may be found
#
PATH=../tools:.:$PATH
export PATH

out=/tmp/t-surface-data-term-1-$$
indir=${srcdir}/../test

trap "rm -f $out" 0 1 2 15


for extra in 00 05 15; do
    ./surface-data-term-1 -dt_extra ${extra} ${indir}/poly_5120.obj $out
    diff ${indir}/poly_5120_sdt${extra}.vv $out
done
