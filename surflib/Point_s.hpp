#ifndef POINT_S_H  // -*- C++ -*-
#define POINT_S_H


#include <surflib/Facet_coordinate_frame.hpp>



namespace MNI {


/*! \brief Point on a polyhedral surface.
 *
 */

template <class _Surface>
class Point_s
{
public:
    typedef _Surface                                 Surface;
    typedef Point_s<Surface>                         Self;

    typedef typename Surface::Vertex_const_handle    Vertex_const_handle;
    typedef typename Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef typename Surface::Point_3                Point_3;
    typedef typename Surface::Traits::Vector_3       Vector_3;


    //! Initialize point directly.
    Point_s( Halfedge_const_handle h, double t1, double t2);

    //! Initialize location to a given vertex.
    Point_s( Vertex_const_handle v );

    //! Construct location from point.
    Point_s( Halfedge_const_handle h, const Point_3& p );

    // Default constructor is necessary for use with containers.
    /// \todo set _t1 and _t2 to NaN or some singular value.
    Point_s() { _t1 = _t2 = -999; }


    // Accessors
    const Facet_coordinate_frame<Surface>& coordinate_frame() const 
	{ return frame; }

    Halfedge_const_handle halfedge() const { return frame.halfedge(); }
    double t0() const { return 1-_t1-_t2; }
    double t1() const { return _t1; }
    double t2() const { return _t2; }

    Point_3 point() const;


    // Convenient access to underlying frame.  Maybe this should
    // not be here?
    Vertex_const_handle V0() const  { return coordinate_frame().V0(); }
    Vertex_const_handle V1() const  { return coordinate_frame().V1(); }
    Vertex_const_handle V2() const  { return coordinate_frame().V2(); }

    const Point_3& P0() const { return coordinate_frame().P0(); }
    const Point_3& P1() const { return coordinate_frame().P1(); }
    const Point_3& P2() const { return coordinate_frame().P2(); }



    //! Get coordinates with respect to halfedge V1,V2.
    void coordinates_frame_12( double* t1, double* t2 ) const;

    //! Return point with coordinates relative to halfedge V1,V2.
    Self point_frame_12() const;

    //! Get coordinates with respect to halfedge V2,V0.
    void coordinates_frame_20( double* t1, double* t2 ) const;

    //! Return point with coordinates relative to halfedge V2,V0.
    Self point_frame_20() const;

    //! Get coordinates with respect to any halfedge bounding facet.
    void coordinates_set_halfedge( Halfedge_const_handle h,
				   double* t1, double* t2 ) const;

    
    //! Displace point along surface.
    const Self displace( double du1, double du2 ) const;

protected:
    Facet_coordinate_frame<Surface> frame;
    double _t1,_t2;

    // Test case
    friend class t_Point_s;
};



//! True if p lies on h->facet().
template <class _surface>
bool
is_incident_to_facet( const Point_s<_surface>& p,
		      typename _surface::Facet_const_handle f );


//! Find common facet and return one of its halfedges.
template <class _surface>
typename _surface::Halfedge_const_handle
find_common_facet( const Point_s<_surface>& p,
		   const Point_s<_surface>& q );


template <class _surface>
inline typename _surface::Traits::RT
squared_distance( const Point_s<_surface>& p,
		  const Point_s<_surface>& q )
{
    return CGAL::squared_distance( p.point(),
				   q.point() );
}

}

#include "Point_s.cpp"
#endif
