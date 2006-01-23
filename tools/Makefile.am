include $(CGAL_MAKEFILE)

INCLUDES = -I$(top_srcdir)

# AM_CXXFLAGS variables defined in CGAL makefile.
AM_CXXFLAGS = $(CGAL_CXXFLAGS) $(LONG_NAME_PROBLEM_CXXFLAGS)
LDADD = ../surflib/libsurface.la $(CGAL_LDFLAGS) 

bin_PROGRAMS = \
	alpha-shape \
	surface-curvature \
	evaluate-between-surfaces \
	label-connected-regions \
	project-label-nearest-vertex \
	surface-smooth-data \
	surface-average \
	surface-distance-transform \
	surface-map-check \
	surface-map-displacement \
	surface-map-stats \
	surface-resample \
	surface-stats

TESTS = t_alpha_shape

# Common include file
EXTRA_DIST = Surface.hpp t_alpha_shape.sh

CLEANFILES = t_alpha_shape

alpha_shape_SOURCES = alpha-shape.cc
surface_curvature_SOURCES = surface-curvature.cc
evaluate_between_surfaces_SOURCES = evaluate-between-surfaces.cc
label_connected_regions_SOURCES = label-connected-regions.cc 
project_label_nearest_vertex_SOURCES = project-label-nearest-vertex.cc 
surface_smooth_data_SOURCES = surface-smooth-data.cc
surface_average_SOURCES = surface-average.cc
surface_distance_transform_SOURCES = surface-distance-transform.cc
surface_map_check_SOURCES = surface-map-check.cc
surface_map_displacement_SOURCES = surface-map-displacement.cc
surface_map_stats_LDADD = ../registration/mesh_area.o $(LDADD)
surface_map_stats_SOURCES = surface-map-stats.cc 
surface_resample_SOURCES = surface-resample.cc
surface_stats_SOURCES = surface-stats.cc

SUFFIXES = .sh
.sh:
	@rm -f $@
	sed < $< > $@-t \
	    -e 's,xCREATE_TETRA_PATHx,$(CREATE_TETRA_PATH),g' \
	    -e 's,xBICOBJ2VTK_PATHx,$(BICOBJ2VTK_PATH),g'
	@chmod +x $@-t
	mv $@-t $@
