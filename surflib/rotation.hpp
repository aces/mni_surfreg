#ifndef SURFLIB_ROTATION_HPP
#define SURFLIB_ROTATION_HPP


/* Rotation in 3-space.
 */

namespace MNI {


/*! \brief Rotate point about an axis.
 *
 * This function avoids calling \c sin and \c cos functions.
 */
template <class R>
inline
CGAL::Point_3<R> rotate( const CGAL::Point_3<R>& p,
			 const CGAL::Vector_3<R>& axis,
			 double cos_theta,
			 double sin_theta )
{
    CGAL::Vector_3<R> r = p - CGAL::ORIGIN;
    CGAL::Vector_3<R> n( axis );
    
    n = n / sqrt(CGAL::to_double(n*n));

    return CGAL::ORIGIN 
	+ r * cos_theta
	+ n * (n*r)*(1-cos_theta)
	+ CGAL::cross_product(n,r) * sin_theta;
}


/*! \brief Rotate point about an axis.
 *
 * Given point \a P, axis direction \a axis, and an angle \a theta,
 * return the point obtained by rotating \a P in the positive sense
 * about \a axis by angle \a theta.  The positive sense of rotation
 * is given by the usual \em{right hand rule}.
 *
 * This function uses floating-point approximations, so the
 * result is not exact.
 *
 * \param P is the input point
 * \param axis is the axis direction
 * \param theta is rotation angle, in radians
 *
 */
template <class R>
inline
CGAL::Point_3<R> rotate( const CGAL::Point_3<R>& p,
			 const CGAL::Vector_3<R>& axis,
			 double theta )
{
    return rotate( p, axis, std::cos(theta), std::sin(theta) );
}


/*! \brief Rotation functor.
 */
template <class R>
class Rotation
{
    CGAL::Vector_3<R> axis;
    double cos_theta, sin_theta;

    void init()
    {
	double axis_length = std::sqrt(CGAL::to_double(axis*axis));
	CGAL_precondition( axis_length > 0 );
	axis = axis / axis_length;
    }

public:
    //! Construct identity rotation.
    Rotation()
	: axis(1,0,0), cos_theta(1), sin_theta(0)
    {
	init();
    }

    Rotation( const CGAL::Vector_3<R>& a,
	      double c, double s )
	: axis(a), cos_theta(c), sin_theta(s)
    {
	init();
    }

    Rotation( const CGAL::Vector_3<R>& a,
	      double angle )
	: axis(a)
    {
	cos_theta = std::cos(angle);
	sin_theta = std::sin(angle);
	init();
    }

    const CGAL::Point_3<R> operator()( const CGAL::Point_3<R>& p )
    {
	CGAL::Vector_3<R> r = p - CGAL::ORIGIN;

	return CGAL::ORIGIN 
	    + r * cos_theta
	    + axis * (axis*r)*(1-cos_theta)
	    + CGAL::cross_product(axis,r) * sin_theta;
    }
};


/*! \brief Compute rotation of north pole to p.
 *
 * \return Rotation object that rotates the north pole of unit sphere
 *         to point \a p.
 *
 * \pre point \a p is on the unit sphere.
 */
template <class R>
const MNI::Rotation<R> make_pole_rotation( const CGAL::Point_3<R>& p )
{
    // Compute axis and angle of rotation to map north pole to p.  For
    // general p, the axis is given by cross_product( (0,0,1), p ),
    // and the angle cosine is given by dot_product( (0,0,1), p ).
    //
    // If, however, p.x() == p.y() == 0, the transformation sought
    // is either the identity, or a rotation by angle PI.
    //
    CGAL::Vector_3<R> axis( -p.y(), p.x(), 0 );

    double cos_a, sin_a;
    if ( CGAL::to_double(axis*axis) == 0 ) {
	axis = CGAL::Vector_3<R>(1,0,0);
	cos_a = ( p.z() > 0 ? 1 : -1 );
	sin_a = 0;
    } else {
	cos_a = p.z();
	sin_a = std::sqrt( 1 - cos_a*cos_a );
    }

    return MNI::Rotation<R>( axis, cos_a, sin_a );
}


}

#endif
