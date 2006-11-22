include $(CGAL_MAKEFILE)

noinst_PROGRAMS = test-2 

if HAVE_CPPUNIT
cppunit_programs = cppunit_test cppunit_test_approx
endif

check_PROGRAMS = $(cppunit_programs) t_load_surface_file gen-index
TESTS = $(cppunit_programs) test-load-surface

EXTRA_DIST = octahedron.obj test-load-surface.sh
check_DATA = example1.obj

CLEANFILES = t_SurfaceMap_poly.out t_SurfaceMap.out debug.poly.vtk \
	debug.path.vtk test-load-surface


example1.obj: octahedron.obj
	ln -s $(srcdir)/octahedron.obj example1.obj


# Both AM_CXXFLAGS variables are from CGAL's makefile.
# CGAL_MP_FLOAT_ALLOW_INEXACT is defined to allow dividing 
# an MP_Float by a double (required for CGAL 3.2.1).
AM_CXXFLAGS = $(CGAL_CXXFLAGS) $(LONG_NAME_PROBLEM_CXXFLAGS) \
	-DCGAL_MP_FLOAT_ALLOW_INEXACT
LDADD = ../libsurface.la $(bic_LIBS) $(CGAL_LDFLAGS) $(CPPUNIT_LIBS)

AM_CPPFLAGS = -I$(top_srcdir) $(CPPUNIT_CFLAGS) -DSRCDIR=\"$(srcdir)\"


test_2_SOURCES = test-2.cc

cppunit_test_SOURCES = $(top_srcdir)/test/tester.cc tester.hpp Surface.hpp \
	t_FileReader.cc \
	t_Object_file_inserter.cc \
	t_SurfaceMap.cc \
	t_SurfaceMap2.cc \
	t_SurfaceMap_poly.cc

# The surface stuff uses floating point, so we test using
# approximate equality.
cppunit_test_approx_SOURCES = $(top_srcdir)/test/tester.cc tester_approx.hpp \
	t_rotation.cc \
	t_Facet_coordinate_frame.cc \
	t_Vector_s.cc \
	t_Point_s.cc \
	t_SurfaceLocation.cc \
	t_surface_algebra.cc \
	t_SurfaceFunctions.cc \
	t_ShortestPathQuery.cc

t_load_surface_file_SOURCES = t_load_surface_file.cc Surface.hpp
gen_index_SOURCES = gen-index.cc


SUFFIXES = .sh
.sh:
	@rm -f $@
	sed < $< > $@-t \
	    -e 's,xCREATE_TETRA_PATHx,$(CREATE_TETRA_PATH),g' \
	    -e 's,xBICOBJ2VTK_PATHx,$(BICOBJ2VTK_PATH),g'
	@chmod +x $@-t
	mv $@-t $@


