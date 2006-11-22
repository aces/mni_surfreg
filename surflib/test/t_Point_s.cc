#include "tester_approx.hpp"

#include "Surface.hpp"
#include <surflib/Point_s.hpp>

using MNI::Point_s;


class t_Point_s : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;

    typedef Surface::Facet_handle           Facet_handle;

    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Halfedge_const_handle  Halfedge_const_handle;

    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point_3;
    typedef Point_s<Surface>                Point_s;



    /* A surface consisting of a single triangle (a,b,c). */
    Point_3 a,b,c;
    Surface surf;
    Halfedge_handle h_ab;

    
    void setUp()
    {
	a = Point_3(0,0,0);
	b = Point_3(1,0,0);
	c = Point_3(0,1,0);

	h_ab = surf.make_triangle( b,c,a );
    }


    void constructors() 
    {
	Point_s p( h_ab, 0.5,0.5 );
	CPPUNIT_ASSERT_EQUAL( 0.5, p.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.5, p.t2() );

	Vertex_const_handle v = surf.vertices_begin();
	p = Point_s(v);
	CPPUNIT_ASSERT_EQUAL( v->point(), p.point() );
    }


    void accessors()
    {
	Point_s p( h_ab, 0.3, 0.3 );

	CPPUNIT_ASSERT( &* h_ab == &* p.halfedge() );

	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, p.t0(), equality_tolerance );
	CPPUNIT_ASSERT_EQUAL( 0.3, p.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.3, p.t2() );
    }


    void points()
    {
	Point_s p( h_ab, 0.2, 0.1 );

	CPPUNIT_ASSERT_EQUAL( Point_3(0.2,0.1,0), p.point() );
    }

    
    void global_funcs()
    {
	Point_s p( h_ab, 0.2, 0.1 );
	Point_s q( h_ab->next(), 0.2, 0.1 );

	CPPUNIT_ASSERT( is_incident_to_facet(p,h_ab->facet()) );

	Halfedge_const_handle h = find_common_facet(p,q);
	CPPUNIT_ASSERT( h != 0 );
    }


    CPPUNIT_TEST_SUITE( t_Point_s );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( accessors );
    CPPUNIT_TEST( points );
    CPPUNIT_TEST( global_funcs );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_Point_s );
