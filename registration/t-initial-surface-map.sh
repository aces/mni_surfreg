#! /bin/sh

#
# Copyright Alan C. Evans
# Professor of Neurology
# McGill University
#

set -e

out=/tmp/t-initial-surface-map-$$
indir=${srcdir}/../test

trap "rm -f $out" 0 1 2 15


./initial-surface-map ${indir}/poly_20.obj ${indir}/poly_20.obj $out
if ! diff -q ${indir}/poly_20_identity.sm $out; then
    echo "Failed to create identity map."
    exit 1
fi

rm $out

./initial-surface-map -flip-x ${indir}/poly_20.obj ${indir}/poly_20.obj $out
if ! diff -q ${indir}/poly_20_flipx.sm $out; then
    echo "Failed to create x-flipped map."
    exit 1
fi


