/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFACE_SEARCH_H   // -*- C++ -*-
#define SURFACE_SEARCH_H 

#include <surflib/Point_s.hpp>

namespace MNI {


/*! \brief Simple search for vertex that minimizes a function.
 *
 * Walk from a given start vertex until the function no
 * longer decreases.  The vertex returned represents a local
 * minimum in the function \a vf.
 *
 * Function \a vf must be a model of:
 *  double operator() ( typename Vertex::Vertex_const_handle v )
 */
template <class Vertex,class VertexFunctor>
struct VertexWalkSearch
{
    typedef typename Vertex::Vertex_const_handle Vertex_const_handle;

    VertexWalkSearch( VertexFunctor& vf )
	: f(vf) 
    {}

    Vertex_const_handle search( Vertex_const_handle v );

private:
    VertexFunctor f;
};


/* Vertex functor that computes squared distance from vertex 
 * to a fixed point in 3-space.
 */
template <class Vertex>
struct SquaredDistanceToPoint
{
    typedef typename Vertex::Point_3 Point_3;

    SquaredDistanceToPoint( const Point_3& point )
	: p(point) {}

    double operator() ( typename Vertex::Vertex_const_handle v )
    {
	return CGAL::to_double( CGAL::squared_distance( p, v->point() ));
    }

private:
    Point_3 p;
};


/*! \brief Search convex surface for location nearest given point.
 *
 * Search surface containing \a v for the surface point nearest
 * \a p.
 *
 * The precondition that the surface be convex is \em not checked.
 */
template <class Surface_>
const Point_s<Surface_> 
find_nearest_on_convex_surface( typename Surface_::Vertex_const_handle v,
				const typename Surface_::Point_3& p );


/*! \brief Return midpoint shortest path between two points on sphere.
 *
 * \pre Points must lie on sphere with centre at origin.
 * \pre Points must not be antipodal.
 */
template <class Surface_>
const Point_s<Surface_>
find_midpoint_on_sphere( const Point_s<Surface_>& p,
			 const typename Surface_::Point_3& q );


template <class Surface_>
inline
const Point_s<Surface_>
find_midpoint_on_sphere( const Point_s<Surface_>& p,
			 const Point_s<Surface_>& q )
{
    return find_midpoint_on_sphere<Surface_>( p, q.point() );
}


}

#include "Search.cpp"
#endif
