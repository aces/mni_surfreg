/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef LINEAR_INTERPOLATION_H   // -*- C++ -*-
#define LINEAR_INTERPOLATION_H

#include <functional>
#include <iostream>


namespace MNI {


/*! \brief Compute function value at arbitrary surface point.
 *
 * Adapts a functor that operates on mesh vertices only
 * to a functor that handles arbitrary surface points by
 * linear interpolation.
 *
 * \pre class F must be a unary function that maps surface
 *      vertices to a double result
 * \pre Facets are triangles.
 */

template <class _surface, class Operation>
class LinearInterpolation 
    : public std::unary_function< Point_s<_surface>,
				  typename Operation::result_type >
{
public:
    LinearInterpolation( const Operation& _f ) : f(_f) {}

    typename Operation::result_type 
    operator() ( const typename _surface::Vertex_const_handle& v )
    {
	return f(v);
    }

    typename Operation::result_type 
    operator() ( const Point_s<_surface>& p )
    {
	typename _surface::Halfedge_const_handle h( p.halfedge() );
	typename Operation::result_type F1 = f(h->vertex());

	h = h->next();
	typename Operation::result_type F2 = f(h->vertex());

	h = h->next();
	typename Operation::result_type F0 = f(h->vertex());

//std::cout << "Inter: t1 = " << p.t1() << " t2 = " << p.t2() << " F0 = " << F0 
//          << " F1 = " << F1 << " F2 = " << F2 << " ret = " 
//          << p.t0()*F0 + p.t1()*F1 + p.t2()*F2 << std::endl;

	return p.t0()*F0 + p.t1()*F1 + p.t2()*F2;
    }

protected:
    Operation f;
};


template <class _surface, class Function>
inline LinearInterpolation<_surface,Function>
linear_interpolator( const Function& f )
{
    return LinearInterpolation<_surface,Function>(f);
}


}

#endif
