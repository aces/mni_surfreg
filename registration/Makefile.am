include $(CGAL_MAKEFILE)

# Both AM_CXXFLAGS variables are from CGAL's makefile.
AM_CXXFLAGS = $(CGAL_CXXFLAGS) $(LONG_NAME_PROBLEM_CXXFLAGS) 
LDADD = ../surflib/libsurface.la $(CGAL_LDFLAGS)

INCLUDES = -I$(top_srcdir) @GSL_CFLAGS@

EXTRA_DIST = surface-data-term-1.sh sphere-register.pl
CLEANFILES = $(generated_scripts)

bin_PROGRAMS = surftracc \
	initial-surface-map \
	refine-surface-map
bin_SCRIPTS = $(generated_scripts)

generated_scripts = surface-data-term-1 sphere-register


surftracc_SOURCES = surftracc.cc surftracc.hpp \
	mesh_area.cc \
	data_term.hpp \
	objective.cc objective.hpp \
	optimize_map.cc \
	smooth_map.cc
surftracc_LDADD = @GSL_LIBS@ $(LDADD)

initial_surface_map_SOURCES = initial-surface-map.cc
refine_surface_map_SOURCES = refine-surface-map.cc



SUFFIXES = .sh .pl
.sh:
	@rm -f $@
	sed < $< > $@-t \
	    -e 's,xBINDIRx,$(bindir),g' \
	    -e 's,xVERSIONx,$(VERSION),'
	chmod +x $@-t
	mv $@-t $@

.pl:
	@rm -f $@
	sed < $< > $@-t \
	    -e 's,xPERLx,$(PERL),g' \
	    -e 's,xINCDIRx,$(pkgdatadir),g' \
	    -e 's,xBINDIRx,$(bindir),g' \
	    -e 's,xVERSIONx,$(VERSION),g'
	chmod +x $@-t
	mv $@-t $@
