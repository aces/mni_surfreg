/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "tester_approx.hpp"

#include "Surface.hpp"
#include <surflib/Facet_coordinate_frame.hpp>

using MNI::Facet_coordinate_frame;



class t_Facet_coordinate_frame : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point;
    typedef Surface::Traits::Vector_3       Vector;

    typedef Facet_coordinate_frame<Surface> Frame;



    /* A surface consisting of a single triangle (a,b,c). */
    Point a,b,c;
    Surface surf;
    Halfedge_handle h_ab;

    
    void setUp()
    {
	a = Point(0,0,0);
	b = Point(1,0,0);
	c = Point(0,1,0);

	h_ab = surf.make_triangle( b,c,a );
    }


    void constructors() 
    {
	Frame f(h_ab);
	CPPUNIT_ASSERT( &* f.halfedge() == &* h_ab );
    }


    void accessors()
    {
	Frame f(h_ab);

	CPPUNIT_ASSERT_EQUAL( a, f.P0() );
	CPPUNIT_ASSERT_EQUAL( b, f.P1() );
	CPPUNIT_ASSERT_EQUAL( c, f.P2() );
    }


    void frame_vectors()
    {
	Frame f(h_ab);

	CPPUNIT_ASSERT_EQUAL( Vector(1,0,0), f.e1() );
	CPPUNIT_ASSERT_EQUAL( Vector(0,1,0), f.e2() );
	CPPUNIT_ASSERT_EQUAL( Vector(0,0,1), f.n() );
    }

    void point_mapping()
    {
	Frame f(h_ab);
	double t1 = 0.5;
	double t2 = 0;
	Point p( 0.5, 0, 0 );

	CPPUNIT_ASSERT_EQUAL( p, f.point(t1,t2) );

	double u,v;
	f.coordinates_of( p, &u, &v );
	CPPUNIT_ASSERT_EQUAL( t1, u );
	CPPUNIT_ASSERT_EQUAL( t2,v );
    }


    void vector_mapping()
    {
	Frame f(h_ab->next());
	double t1 = 0.3;
	double t2 = 0.9;

	Vector vec = (c-b)*t1 + (a-b)*t2;

	CPPUNIT_ASSERT_EQUAL( vec, f.vector(t1,t2) );

	double u,v;
	f.coordinates_of( vec, &u, &v );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( t1, u, equality_tolerance );
	CPPUNIT_ASSERT_DOUBLES_EQUAL( t2, v, equality_tolerance );
    }


    CPPUNIT_TEST_SUITE( t_Facet_coordinate_frame );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( accessors );
    CPPUNIT_TEST( frame_vectors );
    CPPUNIT_TEST( point_mapping );
    CPPUNIT_TEST( vector_mapping );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_Facet_coordinate_frame );
