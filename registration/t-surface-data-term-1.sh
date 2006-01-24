#! /bin/sh

set -e


out=/tmp/t-surface-data-term-1-$$

trap "rm -f $out" 0 1 2 15


for extra in 00 05 15; do
    ./surface-data-term-1 -dt_extra ${extra} ../test/poly_5120.obj $out
    diff ../test/poly_5120_sdt${extra}.vv $out
done
