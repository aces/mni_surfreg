include $(CGAL_MAKEFILE)

SUBDIRS = . test doc

INCLUDES = -I$(top_srcdir)

noinst_LTLIBRARIES = libsurface.la

# AM_CXXFLAGS variables defined in CGAL makefile.
AM_CXXFLAGS = $(CGAL_CXXFLAGS) $(LONG_NAME_PROBLEM_CXXFLAGS)
LDADD = $(noinst_LTLIBRARIES) $(CGAL_LDFLAGS) 

EXTRA_DIST = $(headers)

libsurface_la_SOURCES =				\
	src/Object_file.cc				\
	src/VTKFile.cc
libsurface_la_LIBADD = $(bic_LIBS)

headers_hpp = AlphaShapeUtil.hpp Barycentric.hpp Diagnostic_ostream.hpp	\
	DisjointVertexSets.hpp DistanceTransform.hpp			\
	Facet_coordinate_frame.hpp					\
	FileReader.hpp Histogram.hpp					\
	LinearInterpolation.hpp NearestNeighbourInterpolation.hpp	\
	Object_file.hpp Object_file_inserter.hpp			\
	Point_s.hpp Search.hpp						\
	ShortestPathQuery.hpp Statistic.hpp SurfaceFunctions.hpp	\
	SurfaceMap.hpp SurfaceMapBase.hpp SurfaceMapIO.hpp		\
	SurfaceMap_hash.hpp SurfaceMap_hds.hpp SurfaceMap_poly.hpp	\
	VTKFile.hpp VectorStatistic.hpp					\
	Vector_s.hpp load_surface_file.hpp rotation.hpp			\
	surface_checking.hpp vtk_ostream.hpp

headers_cpp = Facet_coordinate_frame.cpp Point_s.cpp	       \
	Polygons_inserter.cpp Quadmesh_inserter.cpp Search.cpp \
	ShortestPathQuery.cpp SurfaceMap.cpp SurfaceMapIO.cpp  \
	SurfaceMap_hds.cpp SurfaceMap_poly.cpp Vector_s.cpp

headers = $(headers_hpp) $(headers_cpp)
