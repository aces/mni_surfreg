/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFACE_FUNCTIONS_H   // -*- C++ -*-
#define SURFACE_FUNCTIONS_H 

#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <surflib/Point_s.hpp>


namespace MNI {


//! Return number of neighbouring vertices.
template <class _surface>
inline int
degree( typename _surface::Vertex_const_handle v )
{
    int d = 0;
    typename _surface::Halfedge_around_vertex_const_circulator h 
	= v->vertex_begin();
    CGAL_For_all(h,v->vertex_begin()) {
	++d;
    }
    return d;
}


//! True if any incident edge is border.
template <class _surface>
inline bool
is_border( typename _surface::Vertex_const_handle v )
{
    typename _surface::Halfedge_around_vertex_const_circulator h 
	= v->vertex_begin();
    CGAL_For_all(h,v->vertex_begin()) {
	if ( h->is_border_edge() )
	    return true;
    }
    return false;
}


//! Find halfedge (u,v).
template <class _surface>
inline typename _surface::Halfedge_const_handle
find_halfedge( typename _surface::Vertex_const_handle u,
	       typename _surface::Vertex_const_handle v )
{
    typename _surface::Halfedge_const_handle h = v->halfedge();
    do {
	if ( &* h->opposite()->vertex() == &* u )
	    return h;
	h = h->next_on_vertex();
    } while ( h != v->halfedge() );

    return 0;
}


//! Find halfedge (u,v).
template <class _surface>
inline typename _surface::Halfedge_handle
find_halfedge( typename _surface::Vertex_handle u,
	       typename _surface::Vertex_handle v )
{
    typename _surface::Halfedge_handle h = v->halfedge();
    do {
	if ( &* h->opposite()->vertex() == &* u )
	    return h;
	h = h->next_on_vertex();
    } while ( h != v->halfedge() );

    return typename _surface::Halfedge_handle(0);
}


//! True if vertices joined by an edge.
template <class _surface>
inline bool
are_neighbours( typename _surface::Vertex_const_handle u,
		typename _surface::Vertex_const_handle v )
{
    return find_halfedge<_surface>(u,v) != 0;
}


//! Find halfedge incident to both.
template <class _surface>
inline typename _surface::Halfedge_handle
find_incident_halfedge( typename _surface::Vertex_handle v,
			typename _surface::Facet_handle f )
{
    typename _surface::Halfedge_around_vertex_circulator h 
	= v->vertex_begin();
    CGAL_For_all(h,v->vertex_begin()) {
	if ( &* f == &* h->facet() )
	    return h;
    }
    return 0;
}


//! Find common neighbour.
template <class _surface>
inline typename _surface::Vertex_const_handle
find_common_vertex_neighbour( typename _surface::Vertex_const_handle u,
			      typename _surface::Vertex_const_handle v )
{
    typename _surface::Halfedge_const_handle hu = u->halfedge();
    do {
	typename _surface::Halfedge_const_handle hv = v->halfedge();
	do {
	    if ( &* hu->opposite()->vertex() == &* hv->opposite()->vertex() )
		return hu->opposite()->vertex();
	    hv = hv->next_on_vertex();
	} while ( hv != v->halfedge() );
	hu = hu->next_on_vertex();
    } while ( hu != u->halfedge() );
    
    return 0;
}


//! Test whether the locations belong to adjacent facets.
/*! This checks for strict adjacency, excluding the possibility
 * that p and q are on the same facet.
 */
template <class _surface>
inline bool
on_adjacent_facets( const Point_s<_surface>& p,
		    const Point_s<_surface>& q )
{
    typename _surface::Halfedge_const_handle h = p.halfedge();
    typename _surface::Facet_const_handle fq = q.halfedge()->facet();

    return ( &* fq == &* h->opposite()->facet() ||
	     &* fq == &* h->next()->opposite()->facet() ||
	     &* fq == &* h->next()->next()->opposite()->facet() );
}


/* This creates a plane passing through three points
 * on the facet.  The points are oriented in a positive
 * sense when seen from the positive side of the plane.
 *
 * Usage:
 *     std::transform( surf.facets_begin(), surf.facets_end(), 
 *                     surf.planes_begin(),
 *                     Plane_equation() );
 */
struct Plane_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) 
    {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3  Plane;
        Plane p( h->vertex()->point(),
		 h->next()->vertex()->point(),
		 h->next()->next()->vertex()->point());
	CGAL_postcondition( !p.is_degenerate() );
	CGAL_postcondition( CGAL::is_valid(p.a()) );
	CGAL_postcondition( CGAL::is_valid(p.b()) );
	CGAL_postcondition( CGAL::is_valid(p.c()) );
	CGAL_postcondition( CGAL::is_valid(p.d()) );
	return p;
    }
};


//! Return facet angle at h->vertex().
template <class _surface>
inline double
facet_angle( typename _surface::Halfedge_const_handle h )
{
    const typename _surface::Point_3& p = h->vertex()->point();
    typename _surface::Traits::Vector_3 v1
	= h->next()->vertex()->point() - p;
    typename _surface::Traits::Vector_3 v2
	= h->opposite()->vertex()->point() - p;

    double v1_len_v2_len = sqrt( CGAL::to_double(v1*v1*v2*v2) );
    double cos_theta = CGAL::to_double( v1*v2 / v1_len_v2_len );
    return acos(cos_theta);
}


/*! Compute area and area gradient for a triangular surface
 * region.  The region is given by three surface points,
 * A, B, and C which are considered to be joined by geodesics.
 *
 * Note!  This code currently computes a heuristic approximation
 * of the area.  We project the three points to a plane (the plane
 * of A_s.coordinate_frame()) and compute the area of the projected
 * triangle.
 *
 * \param h A,B,C are the region vertices
 * \param dAdv return location for rate of change of
 *             area with respect to location of A
 * \return signed area of projected triangle
 */
template <class _surface>
inline double
triangular_region_area( const Point_s<_surface>& A_s,
			const Point_s<_surface>& B_s,
			const Point_s<_surface>& C_s )
{
    typename _surface::Traits::Vector_3 dummy;
    return triangular_region_area( A_s, B_s, C_s, dummy );
}


template <class _surface>
inline double
triangular_region_area( const Point_s<_surface>& A_s,
			const Point_s<_surface>& B_s,
			const Point_s<_surface>& C_s,
			typename _surface::Traits::Vector_3& dAdv )
{
    typedef typename _surface::Point_3            Point_3;
    typedef typename _surface::Traits::Vector_3   Vector_3;

    const Facet_coordinate_frame<_surface>& frame = A_s.coordinate_frame();
    Vector_3 diff_10 = frame.P1() - frame.P0();
    Vector_3 diff_20 = frame.P2() - frame.P0();

    Vector_3 n = CGAL::cross_product( diff_10, diff_20 );
    double n_squared = CGAL::to_double( n*n );
    double n_len = sqrt(n_squared);

    // The vertices of the triangle live at positions A_s.point(),
    // B_s.point(), and C_s.point() in 3-space.  We'll project each of
    // these to the plane of "frame", denoting them A,B, and C,
    // respectively.  Note that A already lives in "frame".

    Point_3 A = A_s.point();
    Point_3 B = B_s.point() + (n*(A - B_s.point())/n_squared)*n;
    Point_3 C = C_s.point() + (n*(A - C_s.point())/n_squared)*n;

    Vector_3 K = CGAL::cross_product( B-A, C-A ) / 2;
    double area = CGAL::to_double( K*n/n_len );

    dAdv = CGAL::cross_product( B-C, n ) / (2*n_len);
    return area;
}


//! Return unsigned area of triangular facet.
template <class _surface>
inline double 
triangular_facet_area( typename _surface::Facet_const_handle f )
{
    return triangular_facet_area<_surface>( f->halfedge() );
}


template <class _surface>
inline double 
triangular_facet_area( typename _surface::Halfedge_const_handle h )
{
    typedef typename _surface::Point_3           Point_3;
    typedef typename _surface::Traits::Vector_3  Vector_3;

    CGAL_precondition( h->next()->next()->next() == h );

    Point_3 v0 = h->vertex()->point();
    h = h->next();
    Point_3 v1 = h->vertex()->point();
    h = h->next();
    Point_3 v2 = h->vertex()->point();
    
    Vector_3 K = CGAL::cross_product(v1 - v0, v2 - v0) / 2;
    return sqrt( CGAL::to_double(K*K) );
}

template <class _surface>
inline double 
triangular_facet_area( const typename _surface::Facet* f )
{
    return triangular_facet_area<_surface>( f->halfedge() );
}


/*! Compute twice the signed area of triangle ABC, and the
 * gradient with respect to the position of vertex A.
 * The sign of the area is determined by the supplied normal.
 *
 * \param h A,B,C are the vertices
 * \param dArea_dA return location for rate of change of
 *             twice the signed area with respect to location of A
 * \return twice signed area of projected triangle
 */
template <class R>
inline typename R::RT
triangle_area2_deriv( const CGAL::Point_3<R>& A,
		      const CGAL::Point_3<R>& B,
		      const CGAL::Point_3<R>& C,
		      const CGAL::Vector_3<R>& normal,
		      CGAL::Vector_3<R>& dArea_dA )
{
    dArea_dA = CGAL::cross_product(B-C,normal);
    return normal * CGAL::cross_product(B-A,C-A);
}


/*! Compute circumcentre in barycentric coordinates.
 *
 * The coordinate frame is <V0,V1,V2>.
 * Result is not rounded.
 */
template <class K>
inline void
circumcentre( const typename K::Point_3& V0,
	      const typename K::Point_3& V1,
	      const typename K::Point_3& V2,
	      double* t1,
	      double* t2 )
{
    typename K::Vector_3 v01(V1-V0);
    typename K::Vector_3 v02(V2-V0);

    double v01_sq = v01*v01;
    double v02_sq = v02*v02;
    double v0102 = v01*v02;

    double D = v01_sq*v02_sq - (v0102*v0102);
    
    *t1 = v02_sq * (v01_sq - v0102) / (2.0 * D);
    *t2 = v01_sq * (v02_sq - v0102) / (2.0 * D);
}


/*! \brief Compute circumcentre.
 *
 * Assumes triangular facet.  Returned barycentric coordinates
 * relative to \a h.
 */
template <class _surface>
inline void
circumcentre( const typename _surface::Halfedge_const_handle h,
	      double* t1, 
	      double* t2 )
{
    circumcentre<typename _surface::Traits>
	( h->opposite()->vertex()->point(),
	  h->vertex()->point(),
	  h->next()->vertex()->point(),
	  t1,t2 );
}


/*! \brief Compute circumcentre.
 *
 * Assumes triangular facet.  Returned barycentric coordinates
 * relative to \a f->halfedge().
 */
template <class _surface>
inline void
circumcentre( const typename _surface::Facet_const_handle f,
	      double* t1, 
	      double* t2 )
{
    circumcentre<_surface>( f->halfedge(), t1, t2 );
}

}

#endif
