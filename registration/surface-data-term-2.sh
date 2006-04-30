#! /bin/sh
#
# Compute a data value at each term suitable for registering inner
# cortical surface (white/gray interface).
#
# The value assigned to each vertex is the distance transform from a
# set of "seed" vertices.  The distance transform is computed using
# approximate geodesic distances on the surface.
#
# The seed vertices are identified as follows.
# 1. For each vertex, a threshold value is computed.  For alpha above this
#    threshold, the vertex will be on the interior of the alpha-shape of the
#    point set.  A convex hull vertex is never interior, so the threshold
#    is infinity.
# 2. Vertices that have upper threshold greater than alpha=100 
#    are selected as seed points.
#

set -e

PATH=xBINDIRx:/usr/local/gnu/bin:$PATH

dist_trans_extra_vert=5
while [ -n "$1" ]; do
    case "$1" in
      -dt_extra) dist_trans_extra_vert=$2; shift 2;;
      *) break;
    esac
done

if [ $# -ne 2 ]; then
    echo >&2 "usage: $0 [options] cortex.obj output.vv"
    echo >&2 "  -dt_extra N   extra vertices in final distance transform"
    exit 1
fi

obj=$1
out=$2


workdir=/tmp/surfdata-$$
mkdir -m 700 $workdir

trap "rm -rf $workdir" 0 1 2 15

alpha=${workdir}/alpha.vv
ch_seed=${workdir}/generate-seed-1.$$
ch_dt=${workdir}/generate-seed-2.$$
tmp_seed=${workdir}/generate-seed-3.$$
seed=${workdir}/alpha_seed.vv


# Compute alpha-shape intervals in file $alpha.
#
alpha-shape -vv $alpha $obj

# Use alpha shape interval upper value to obtain seeds.
#
awk '{x=(($3 == -1) || ($3 > 100)); print x}' $alpha > $seed

# Compute distance transform from seed vertices.
#
surface-distance-transform -extra $dist_trans_extra_vert $obj $seed $out
