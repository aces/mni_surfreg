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
# 2. A rough distance transform from convex hull vertices is computed in
#    order to mask out points that are "too deep".  The depth threshold is
#    selected to be 35 mm.
# 3. Vertices that are less than 35 mm from a convex hull point 
#    AND have upper threshold greater than alpha=100 
#    are selected as seed points.
#

set -e

PATH=xBINDIRx:$PATH


if [ $# -ne 2 ]; then
    echo >&2 "usage: $0 cortex.obj output.vv"
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

# Let $ch_seed be set of convex hull vertices -- those with threshold infinity.
#     (the value -1 is used as the alpha threshold for such vertices)
# Let $ch_dt be the distance transform of previous.
#
awk '{print ($3 == -1)}' $alpha > $ch_seed
surface-distance-transform $obj $ch_seed $ch_dt

# Use alpha shape interval value and convex hull transform
# to obtain final set of seeds (file $seed, using $tmp_seed as
# temporary storage).
#
cut -f3 -d' ' $alpha | cat -n > $tmp_seed
cat -n $ch_dt | join $tmp_seed - | \
    awk '{x=(($2 == -1) || ($2 > 100)) && ($3 < 35); print x}' > $seed

# Compute distance transform from seed vertices.
#
surface-distance-transform -extra 5 $obj $seed $out
