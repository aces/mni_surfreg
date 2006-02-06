#! /bin/sh

set -e


out=/tmp/t-initial-surface-map-$$

trap "rm -f $out" 0 1 2 15

./initial-surface-map ../test/poly_20.obj ../test/poly_20.obj $out
if ! diff -q ../test/poly_20_identity.sm $out; then
    echo "Failed to create identity map."
    exit 1
fi

rm $out

./initial-surface-map -flip-x ../test/poly_20.obj ../test/poly_20.obj $out
if ! diff -q ../test/poly_20_flipx.sm $out; then
    echo "Failed to create x-flipped map."
    exit 1
fi


