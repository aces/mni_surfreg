#ifndef MNI_SURFLIB_BARYCENTRIC_HPP
#define MNI_SURFLIB_BARYCENTRIC_HPP


namespace MNI {


/*! \brief Resolve vector along two independent basis vectors.
 *
 * \param v is the vector to resolve
 * \param b1 and \a b2 are the basis vectors
 * \param t1 and \a t2 are the returned coefficients
 *
 * Solve v = t1.b1 + t2.b2, for \a t1 and \a t2, where
 * \a v, \a b1, and \a b2 are vectors (2-vectors or 3-vectors)
 * and \a t1 and \a t2 are double-precision real numbers.
 *
 * \pre \a b1 and \a b2 are linearly independent vectors.
 *
 * Note that this computation is performed using double precision
 * arithmetic, with attendant roundoff problems.  In particular,
 * if \a v corresponds to a point inside triangle with sides given by
 * \a b1 and \a b2, the resulting coordinates may not satisfy the
 * barycentric inequalities for a point inside the triangle.
 */
template <class Vector>
inline void
coordinates_of_vector( const Vector& v, 
		       const Vector& b1, const Vector& b2,
		       double* t1, double* t2 )
{
    /*
     * Solve system v = t1*b1 + t2*b2 by taking dot product
     * with b1 and with b2, resulting in matrix equation:
     *
     *    / a  b \  / t1 \  =  / e \
     *    \ b  c /  \ t2 /     \ f /
     */

    double a = CGAL::to_double( b1*b1 );
    double b = CGAL::to_double( b1*b2 );
    double c = CGAL::to_double( b2*b2 );

    CGAL_precondition( a*c - b*b != 0 );

    double e = CGAL::to_double( v*b1 );
    double f = CGAL::to_double( v*b2 );

    double det = a*c - b*b;

    *t1 = (c*e-b*f)/det;
    *t2 = (-b*e+a*f)/det;
}


/*! Round areal coordinates.
 *
 * Areal coordinates are barycentric coordinates that sum to 1.  Thus
 * only two coordinates are independent, the third is given by
 * t0 = 1 - t1 - t2.
 *
 * Each coordinate must also lie in the range [0,1].  These conditions
 * give six inequalities (though three are redundant) that the two
 * coordinates must satisfy.
 *
 * Coordinates of constructed points generally do not obey these
 * inequalities, due to floating-point roundoff during computation.
 * This function will perturb the two coordinate values, by an amount
 * up to that specified by \a tolerance, so that the inequalities are
 * satisfied.
 *
 * If assertions are enabled and the input is not within tolerance,
 * the program will be aborted.
 */
inline void 
round_areal_coordinates( double& t1,
			 double& t2,
			 double tolerance = 1e-8 )
{
    if ( t1 < 0 ) {
	CGAL_assertion( t1 > -tolerance );
	t1 = 0;
    } else if ( t1 > 1 ) {
	CGAL_assertion( t1 < 1+tolerance );
	t1 = 1;
    }

    if ( t2 < 0 ) {
	CGAL_assertion( t2 > -tolerance );
	t2 = 0;
    } else if ( t2 > 1 ) {
	CGAL_assertion( t2 < 1+tolerance );
	t2 = 1;
    }

    double t0 = 1 - t1 - t2;
    if ( t0 < 0 ) {
	CGAL_assertion( t0 > -tolerance );
	t1 += t0/2.0;
	t2 = 1 - t1;
	t1 = 1 - t2; // in case that t1 is so small that t2 rounds off to 1.
    }
    // Cannot have _t0 > 1, since _t1,_t2 known >= 0.

    // A little paranoia is never misplaced (?)
    CGAL_postcondition( t1 + t2 <= 1 );
    CGAL_postcondition( t1 >= 0 );
    CGAL_postcondition( t2 >= 0 );
}

}

#endif
