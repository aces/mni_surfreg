#! /bin/sh

#
# Copyright Alan C. Evans
# Professor of Neurology
# McGill University
#

set -e

indir=${srcdir}/../test
out=/tmp/t_alpha_shape-$$

trap "rm -f $out" 0 1 2 15

./alpha-shape -vv $out ${indir}/poly_5120.obj
diff ${indir}/poly_5120_alpha.vv $out
