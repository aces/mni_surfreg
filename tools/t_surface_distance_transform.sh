#! /bin/sh

set -e


out=/tmp/t_distance_transform-$$

trap "rm -f $out" 0 1 2 15


for extra in 00 05 15; do
    ./surface-distance-transform -extra ${extra} \
	../test/poly_5120.obj ../test/poly_5120_chull.vv $out
    diff ../test/poly_5120_dt${extra}.vv $out
done
