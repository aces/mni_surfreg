#include "tester_approx.hpp"

#include "Surface.hpp"
#include <surflib/ShortestPathQuery.hpp>
#include "test/tester_meshes.hpp"
#include <surflib/SurfaceFunctions.hpp>

#include <fstream>


using namespace MNI;


class t_ShortestPathQuery : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                 Surface;

    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Facet_iterator         Facet_iterator;

    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Halfedge_iterator      Halfedge_iterator;

    typedef Surface::Vertex_handle          Vertex_handle;
    typedef Surface::Vertex_iterator        Vertex_iterator;

    typedef Surface::Point_3                Point_3;
    typedef Point_s<Surface>                Point_s;

    typedef std::vector<Point_s>::size_type size_type;




    Surface s1;                  // planar
    Surface t0,t1;               // tetrahedron & refinement
    Vertex_handle vert[6];       // handles for s1


    void setUp()
    {
	MakeMesh_1<Surface> mm1;
	s1.delegate( mm1 );

	MakeTetra_0<Surface> mt0;
	t0.delegate( mt0 );

	MakeTetra_1<Surface> mt1;
	t1.delegate( mt1 );

	Vertex_iterator v = s1.vertices_begin();
	for( int i = 0; v != s1.vertices_end(); ++v,++i ) {
	    vert[i] = &*v;
 	}
    }


    void constructors() 
    {
	MNI::ShortestPathQuery<Surface> query(s1,3);
    }


    void find_path()
    {
	using namespace std;

	MNI::ShortestPathQuery<Surface> query(s1,5);

	Point_s a( find_halfedge<Surface>(vert[0],vert[3]), 0.1, 0.1);
	Point_s b( find_halfedge<Surface>(vert[3],vert[1]), 0.3, 0.1);
	
	std::vector<Point_s> path;
	/*
	cout << "query.find_path from: " << a 
	     << ", to: " << b << endl;
	*/
	double dist = query.find_path(a,b,path);
	double dist_true = sqrt(CGAL::to_double(squared_distance(a,b)));

	CPPUNIT_ASSERT( dist >= dist_true );
	CPPUNIT_ASSERT_EQUAL( std::vector<Point_s>::size_type(4), path.size() );

	a = Point_s( find_halfedge<Surface>(vert[0],vert[3]), 0.1, 0.1);
	b = Point_s( find_halfedge<Surface>(vert[0],vert[3]), 0.3, 0.1);
	Point_s c = query.midpoint(a,b);
	Point_3 c_expected =  CGAL::midpoint( a.point(), b.point() );

	CPPUNIT_ASSERT_EQUAL( c_expected, c.point() );

	return;

	cout << "Distance " << dist 
	     << " (" << dist / dist_true << " x actual)" << endl
	     << "Path length " << path.size() << endl;

	std::vector<Point_s>::iterator pi = path.begin();
	for( ; pi != path.end(); ++pi )
	    cout << "    " << pi->point() << endl;
    }


    // Find midpoint on path with source and target on neighbouring
    // facets.
    void assert_good_midpoint( Surface& s,
			       MNI::ShortestPathQuery<Surface>& query,
			       Facet_handle f0, Facet_handle f1 )
    {
	Point_s a( f0->halfedge(), 0.1, 0.1);
	Point_s b( f1->halfedge(), 0.3, 0.1);
	
	std::vector<Point_s> path;
	double dist = query.find_path(a,b,path);

	CPPUNIT_ASSERT( dist > 0 );
	CPPUNIT_ASSERT( path.size() > 2 );
	if (path.size() != 3 )
	    return;

	Point_s c = query.midpoint(a,b);
	Point_s c1 = query.pathpoint_fraction(a,b,0.5);
	Point_s c2 = query.pathpoint_distance(a,b,dist/2.0);
	CPPUNIT_ASSERT_EQUAL( c.point(), c1.point() );
	CPPUNIT_ASSERT_EQUAL( c.point(), c2.point() );

	MNI::write_vtk( "debug.poly.vtk", s, "t_ShortestPathQuery polyhedron" );

	std::ofstream vtk("debug.path.vtk");
	MNI::write_vtk_header( vtk, "t_ShortestPathQuery path" );

	vtk << "POINTS " << path.size()+1 << " double" << std::endl;
	std::vector<Point_s>::iterator pi = path.begin();
 	for( ; pi != path.end(); ++pi ) {
	    Point_3 p = pi->point();
	    vtk << CGAL::to_double(p.x()) << ' '
		<< CGAL::to_double(p.y()) << ' '
		<< CGAL::to_double(p.z()) << std::endl;
	}
	Point_3 p = c.point();
	vtk << CGAL::to_double(p.x()) << ' '
	    << CGAL::to_double(p.y()) << ' '
	    << CGAL::to_double(p.z()) << std::endl;
	
	size_type len = path.size();
	vtk << "LINES 1 " << len+1 << std::endl
	    << len;
	for( size_type i = 0; i < len; ++i )
	    vtk << ' ' << i;
	vtk << std::endl;

	// Points a and b are on neighbouring facets, so point c
	// must be on one of the two facets.
	if ( &* c.halfedge()->facet() == &* a.halfedge()->facet() ) {
	    double d_ac = sqrt(CGAL::to_double(squared_distance(a,c)));
	    CPPUNIT_ASSERT_EQUAL( dist/2.0, d_ac );
	} else {
	    CPPUNIT_ASSERT( &* c.halfedge()->facet() 
			    == &* b.halfedge()->facet() );
	    double d_bc = sqrt(CGAL::to_double(squared_distance(b,c)));
	    CPPUNIT_ASSERT_EQUAL( dist/2.0, d_bc );
	}
    }	


    void find_path_tetra()
    {
	MNI::ShortestPathQuery<Surface> query(t1,5);
	
	Facet_iterator f0 = t1.facets_begin();
	CGAL_For_all(f0,t1.facets_end()) {
	    Facet_iterator f1 = f0;
	    for( ++f1; f1 != t1.facets_end(); ++f1 ) {
		assert_good_midpoint(t1,query,f0,f1);
	    }
	}
    }



    CPPUNIT_TEST_SUITE( t_ShortestPathQuery );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( find_path );
    CPPUNIT_TEST( find_path_tetra );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_ShortestPathQuery );
