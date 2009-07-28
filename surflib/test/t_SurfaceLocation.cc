/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "tester.hpp"

#include "Surface.hpp"
#include <surflib/Point_s.hpp>


using MNI::Point_s;



class t_SurfaceLocation : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point_3;

    typedef Point_s<Surface>        Point_s;



    /* A surface consisting of a single triangle (a,b,c). */
    Point_3 a,b,c;
    Surface surf;
    Halfedge_handle h_ca;

    
    void setUp()
    {
	a = Point_3(0,0,0);
	b = Point_3(1,0,0);
	c = Point_3(0,1,0);

	h_ca = surf.make_triangle( a,b,c );
    }


    void constructors() 
    {
	Vertex_const_handle v = surf.vertices_begin();
	Point_s loc1(v);
	CPPUNIT_ASSERT( &*loc1.V0() == &*v );
	CPPUNIT_ASSERT_EQUAL( 0.0, loc1.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.0, loc1.t2() );

	Point_s loc2(h_ca,0,0);
	CPPUNIT_ASSERT_EQUAL( 0.0, loc2.t1() );

	Point_3 mid_ca(0,0.5,0);
	Point_s loc3( h_ca, mid_ca );
	CPPUNIT_ASSERT_EQUAL( 0.5, loc3.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.0, loc3.t2() );
	CPPUNIT_ASSERT_EQUAL( mid_ca, loc3.point() );

	Point_3 q( 0.1, 0.1, 0 );
	loc3 = Point_s( h_ca, q );
	CPPUNIT_ASSERT_EQUAL( q, loc3.point() );
	
	q = Point_3( 0.8, 0.1, 0 );
	loc3 = Point_s( h_ca, q );
	CPPUNIT_ASSERT_EQUAL( q, loc3.point() );
	
	q = Point_3( 0.1, 0.8, 0 );
	loc3 = Point_s( h_ca, q );
	CPPUNIT_ASSERT_EQUAL( q, loc3.point() );
    }


    void point()
    {
	Point_s loc( h_ca, 0.5, 0 );
	CPPUNIT_ASSERT_EQUAL( loc.point(), Point_3(0,0.5,0) );
    }


    CPPUNIT_TEST_SUITE( t_SurfaceLocation );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( point );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_SurfaceLocation );
