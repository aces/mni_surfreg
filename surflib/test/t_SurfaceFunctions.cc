#include "tester_approx.hpp"

#include "Surface.hpp"
#include "test/tester_meshes.hpp"

#include <surflib/SurfaceFunctions.hpp>
#include <surflib/Vector_s.hpp>


using namespace MNI;


class t_SurfaceFunctions : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point_3;
    typedef Surface::Traits::Vector_3       Vector_3;

    typedef Point_s<Surface>                Point_s;
    typedef Vector_s<Surface>               Vector_s;


    Surface s1;
    Vertex_const_handle vert[6];

    void setUp()
    {
	MakeMesh_1<Surface> mm1;
	s1.delegate( mm1 );

	Vertex_const_iterator v = s1.vertices_begin();
	for( int i = 0; v != s1.vertices_end(); ++v,++i ) {
	    vert[i] = &*v;
 	}
    }

    
    void halfedge()
    {
	Halfedge_const_handle h = find_halfedge<Surface>( vert[0],vert[5] );
	CPPUNIT_ASSERT( vert[5] == &* h->vertex() );
	CPPUNIT_ASSERT( vert[0] == &* h->opposite()->vertex() );
    }


    void common_neighbour()
    {
	Vertex_const_iterator n01 
	    = find_common_vertex_neighbour<Surface>( vert[0], vert[1] );
	CPPUNIT_ASSERT( vert[3] == &* n01 );

	Vertex_const_iterator n03
	    = find_common_vertex_neighbour<Surface>( vert[0], vert[3] );
	CPPUNIT_ASSERT( vert[5] == &* n03 );

	Vertex_const_iterator 
	    n12 = find_common_vertex_neighbour<Surface>( vert[1], vert[2] );
	CPPUNIT_ASSERT( vert[4] == &* n12 );

	Vertex_const_iterator 
	    n02 = find_common_vertex_neighbour<Surface>( vert[0], vert[2] );
	CPPUNIT_ASSERT( vert[5] == &* n02 );

	Point_s v0(vert[0]);
	Point_s p = v0.displace( 0.25, 0.5 );
	Point_s q = v0.displace( 0.75, 0.5 );
	Point_s r = v0.displace( 1.25, 0.5 );

	CPPUNIT_ASSERT( !on_adjacent_facets( p,p ));
	CPPUNIT_ASSERT( on_adjacent_facets( p,q ));
	CPPUNIT_ASSERT( !on_adjacent_facets( p,r ));
	CPPUNIT_ASSERT( on_adjacent_facets( q,r ));
    }


    void area()
    {
	Point_s v0(vert[0]);
	Point_s A = v0.displace( 1, 0 );
	Point_s B = v0.displace( 1, 1 );
	Point_s C = v0.displace( 0.5, 1 );

	CPPUNIT_ASSERT_EQUAL( Point_3(2,0,0), A.point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(3,1,0), B.point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(2,1,0), C.point() );

	double area = triangular_region_area( A,B,C );
	CPPUNIT_ASSERT_EQUAL( 0.5, area );
    }


    CPPUNIT_TEST_SUITE( t_SurfaceFunctions );
    CPPUNIT_TEST( halfedge );
    CPPUNIT_TEST( common_neighbour );
    CPPUNIT_TEST( area );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_SurfaceFunctions );
    
