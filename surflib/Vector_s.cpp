/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

// -*- C++ -*-


namespace MNI {


template <class Surface_>
inline Vector_s<Surface_>::Vector_s( Halfedge_const_handle h, 
				     double t1, double t2)
    : frame(h)
{
    _t1 = t1;
    _t2 = t2;
}


template <class Surface_>
inline void
Vector_s<Surface_>::coordinates_set_halfedge( Halfedge_const_handle g,
					      double* u, double* v ) const
{
    CGAL_precondition( &* g->facet() == &* halfedge()->facet() );

    if ( &* g == &* halfedge() ) {
	*u = t1();  *v = t2();
	return;
    }

    // This should be made more efficient.
    Vector_3 v3( frame.vector(t1(),t2()) );
    Facet_coordinate_frame<Surface_> frame_g(g);
    frame_g.coordinates_of( v3, u,v );
}


template <class _Surface>
inline const typename _Surface::Traits::Vector_3 
Vector_s<_Surface>::vector() const
{
    return frame.vector( t1(),t2() );
}



// -------------------------- //
// --- Algebra on surface --- //
// -------------------------- //



/*! Add surface vector to point.  The result is defined only
 *  when the frames of the vector and the point are on the 
 *  same facet.
 */
template <class Surface_>
inline
const Point_s<Surface_> operator+( const Point_s<Surface_>& p,
				   const Vector_s<Surface_>& v )
{
    // If the frame of the point and of v differ, transform the former
    // to the frame of v.  This is simpler than the other
    // transformation.

    // Reference frame
    typedef Facet_coordinate_frame<Surface_> Frame;
    Frame frame = v.coordinate_frame();

    // Compute barycentric coordinates of new point in (u1,u2)
    double u1,u2;
    p.coordinates_set_halfedge( frame.halfedge(), &u1,&u2 );

    return add_vector( frame, u1,u2, v.t1(),v.t2() );
}


template <class Surface_>
inline
const Point_s<Surface_> operator+( const Vector_s<Surface_>& v,
				   const Point_s<Surface_>& p )
{
    return p + v;
}


template <class Surface_>
inline
const Vector_s<Surface_> operator-( const Point_s<Surface_>& p,
				    const Point_s<Surface_>& q )
{
    double q1,q2;
    q.coordinates_set_halfedge(p.halfedge(), &q1, &q2);
    return Vector_s<Surface_>( p.halfedge(),
			       p.t1() - q1,
			       p.t2() - q2 );
}


template <class _Surface>
inline const Vector_s<_Surface>
Vector_s<_Surface>::operator+( const Self& w ) const
{
    double u,v;
    coordinates_set_halfedge( w.halfedge(), &u, &v );
    return Self( halfedge(), u+w.t1(), v+w.t2() );
}


template <class _Surface>
inline const Vector_s<_Surface>
Vector_s<_Surface>::operator-( const Self& w ) const
{
    double u,v;
    coordinates_set_halfedge( w.halfedge(), &u, &v );
    return Self( halfedge(), u-w.t1(), v-w.t2() );
}


template <class _Surface>
inline const Vector_s<_Surface>
Vector_s<_Surface>::operator-() const
{
    return Self( halfedge(), -t1(), -t2() );
}


template <class _Surface>
inline const Vector_s<_Surface>
Vector_s<_Surface>::operator*(double c) const
{
    return Self( halfedge(), c*t1(), c*t2() );
}


template <class _Surface>
inline const Vector_s<_Surface>
Vector_s<_Surface>::operator/(double c) const
{
    return Self( halfedge(), t1()/c, t2()/c );
}


}

