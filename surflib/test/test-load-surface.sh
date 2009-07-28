#! /bin/sh

#
# Copyright Alan C. Evans
# Professor of Neurology
# McGill University
#

set -e

runtest () {
    # f v
    f=$1
    v=$2
    xCREATE_TETRA_PATHx sphere 0 0 0 1 1 1 $f
    ./gen-index $v > sphere.vv
    xBICOBJ2VTK_PATHx -vertex_s index sphere.vv sphere.obj sphere.vtk
    ./t_load_surface_file sphere.obj sphere.vv sphere.vtk
    ./t_load_surface_file sphere.vtk index sphere.vtk
}

runtest    20    12
runtest    80    42
runtest   320   162
runtest  1280   642
runtest  5120  2562
runtest 20480 10242
runtest 81920 40962

rm -f sphere.vv sphere.obj sphere.vtk

exit 0
