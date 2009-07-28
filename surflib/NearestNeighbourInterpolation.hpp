/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef NEARESTNEIGHBOUR_INTERPOLATION_HPP
#define NEARESTNEIGHBOUR_INTERPOLATION_HPP

#include <functional>


namespace MNI {


/*! \brief Compute function value at arbitrary surface point.
 *
 * Adapts a functor that operates on mesh vertices only
 * to a functor that handles arbitrary surface points using
 * nearest-neighbour interpolation.
 *
 * \pre class F must be a unary function of surface verticesn
 * \pre Facets are triangles.
 */

template <class _surface, class Operation>
class NearestNeighbourInterpolation 
    : public std::unary_function< Point_s<_surface>,
				  typename Operation::result_type >
{
public:
    NearestNeighbourInterpolation( const Operation& _f ) : f(_f) {}

    typename Operation::result_type 
    operator() ( const typename _surface::Vertex_const_handle& v )
    {
	return f(v);
    }

    typename Operation::result_type 
    operator() ( const Point_s<_surface>& p )
    {
	typename Operation::result_type val = f(p.V0());

	if ( p.t1() < p.t2() ) {
	    // Choose between V0 and V1
	    if ( p.t1() < p.t0() )
		val = f(p.V1());
	} else {
	    // Choose between V0 and V2
	    if ( p.t2() < p.t0() )  
		val = f(p.V2());
	}
	return val;
    }

protected:
    Operation f;
};


template <class _surface, class Function>
inline NearestNeighbourInterpolation<_surface,Function>
nearest_neighbour_interpolator( const Function& f )
{
    return NearestNeighbourInterpolation<_surface,Function>(f);
}


}

#endif
