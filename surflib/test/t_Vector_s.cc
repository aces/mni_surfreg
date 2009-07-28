/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "tester_approx.hpp"

#include "Surface.hpp"
#include <surflib/Vector_s.hpp>

using namespace std;
using MNI::Vector_s;


class t_Vector_s : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point;
    typedef Surface::Traits::Vector_3       Vector_3;
    typedef Vector_s<Surface>               Vector_s;



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
	Vector_s v( h_ab, 1,1 );
	CPPUNIT_ASSERT_EQUAL( 1.0, v.t1() );
	CPPUNIT_ASSERT_EQUAL( 1.0, v.t2() );
    }


    void accessors()
    {
	Vector_s v( h_ab, 0.3, 0.3 );

	CPPUNIT_ASSERT( &* h_ab == &* v.halfedge() );

	//cout << 0.4 - v.t0() << endl;
	CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.4, v.t0(), 1e-15 );
	CPPUNIT_ASSERT_EQUAL( 0.3, v.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.3, v.t2() );
    }


    void coordinates_set_halfedge()
    {
	Vector_s v( h_ab, 0.3, 0.9 );
	double t1,t2;

	v.coordinates_set_halfedge( h_ab, &t1,&t2 );
	CPPUNIT_ASSERT_EQUAL( v.vector(), (b-a)*t1 + (c-a)*t2 );

	v.coordinates_set_halfedge( h_ab->next(), &t1,&t2 );
	CPPUNIT_ASSERT_EQUAL( v.vector(), (c-b)*t1 + (a-b)*t2 );

	v.coordinates_set_halfedge( h_ab->next()->next(), &t1,&t2 );
	CPPUNIT_ASSERT_EQUAL( v.vector(), (a-c)*t1 + (b-c)*t2 );
    }	


    // Input must satisfy: u+v = w
    void vector_algebra_test( const Vector_s& u,
			      const Vector_s& v,
			      const Vector_s& w )
    {
	CPPUNIT_ASSERT_EQUAL( w.vector(), (u+v).vector() );
	CPPUNIT_ASSERT_EQUAL( u.vector(), (w-v).vector() );

    }


    void vectors() 
    {
	Vector_s u( h_ab, 0.2, 0.1 );
	Vector_s v( h_ab, 0.3, 0.3 );
	Vector_s w( h_ab, 0.5, 0.4 );

	CPPUNIT_ASSERT_EQUAL( Vector_3(0.2,0.1,0), u.vector() );
	CPPUNIT_ASSERT_EQUAL( Vector_3(-0.2,-0.1,0), (-u).vector() );
	CPPUNIT_ASSERT_EQUAL( Vector_3(0.1,0.05,0), (u/2).vector() );
	CPPUNIT_ASSERT_EQUAL( Vector_3(0.8,0.4,0), (u*4).vector() );

	vector_algebra_test( u, v, w );
    }

    CPPUNIT_TEST_SUITE( t_Vector_s );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( accessors );
    CPPUNIT_TEST( coordinates_set_halfedge );
    CPPUNIT_TEST( vectors );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_Vector_s );
