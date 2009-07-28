/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "tester_approx.hpp"

#include <cmath>
#include <surflib/rotation.hpp>
#include <CGAL/Simple_cartesian.h>



class t_rotation : public CppUnit::TestCase
{
public:
    typedef CGAL::Simple_cartesian<double>  Kernel;
    typedef Kernel::Point_3                 Point;
    typedef Kernel::Vector_3                Vector;


    // Rotate A about axis by angle_B to obtain B
    Point A,B;
    Vector axis;
    double angle_B;

    // Arbitrarily-chosen point, axis, and angle.
    // Rotate by +angle, then -angle, and get same point.
    Point C;
    Vector axis_C;
    double angle_C;

    
    void setUp()
    {
	A = Point(1,0,0);
	B = Point(0,0,1);
	axis = Vector(0,-1,0);
	angle_B = M_PI_2;

	C = Point(2.3,42.34,-3422.23);
	axis_C = Vector( 0.324, 0.932, -2.342 );
	angle_C = 5.3234;
    }


    void rotate_by_angle()
    {
	Point test_B = MNI::rotate( A, axis, angle_B );
	CPPUNIT_ASSERT_EQUAL( B, test_B );

	// Rotate by +angle, then -angle.
	Point test_C = MNI::rotate( C, axis_C, angle_C );
	test_C = MNI::rotate( test_C, axis_C, -angle_C );
	CPPUNIT_ASSERT_EQUAL( C, test_C );

	// Rotate twice in each direction.
	test_C = MNI::rotate( C, axis_C, angle_C );
	test_C = MNI::rotate( test_C, axis_C, -angle_C );
	test_C = MNI::rotate( test_C, axis_C, -angle_C );
	test_C = MNI::rotate( test_C, axis_C, angle_C );
	CPPUNIT_ASSERT_EQUAL( C, test_C );
    }


    void rotate_by_sin_cos()
    {
	double sin_angle_B = 1;
	double cos_angle_B = 0;

	Point test_B = MNI::rotate( A, axis, cos_angle_B, sin_angle_B );
	CPPUNIT_ASSERT_EQUAL( B, test_B );
    }


    void rotate_with_class()
    {
	MNI::Rotation<Kernel> rotator( axis, angle_B );

	Point test_B = rotator( A );
	CPPUNIT_ASSERT_EQUAL( B, test_B );

	// Rotate twice in each direction.
	MNI::Rotation<Kernel> rotator_C_pos( axis_C, angle_C );
	MNI::Rotation<Kernel> rotator_C_neg( axis_C, -angle_C );

	Point test_C = rotator_C_pos( C );
	test_C = rotator_C_neg( test_C );
	test_C = rotator_C_neg( test_C );
	test_C = rotator_C_pos( test_C );
	CPPUNIT_ASSERT_EQUAL( C, test_C );
    }


    CPPUNIT_TEST_SUITE( t_rotation );
    CPPUNIT_TEST( rotate_by_angle );
    CPPUNIT_TEST( rotate_by_sin_cos );
    CPPUNIT_TEST( rotate_with_class );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_rotation );
