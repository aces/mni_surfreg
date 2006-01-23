#! /bin/sh

set -e


out=/tmp/t_alpha_shape-$$

trap "rm -f $out" 0 1 2 15

./alpha-shape -vv $out ../test/poly_5120.obj
diff ../test/poly_5120_alpha.vv $out
