/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

// -*- C++ -*-

#include <surflib/Barycentric.hpp>


namespace MNI {


// -----------------//
// The constructors //
// -----------------//


/*!
 * \pre t0 >= 0 (t0 = 1 - t1 - t2)
 * \pre \a t1 >= 0
 * \pre \a t2 >= 0
 */
template <class Surface_>
inline
Point_s<Surface_>::Point_s( Halfedge_const_handle h, double t1, double t2)
    : frame(h)
{
    CGAL_precondition( t1 >= 0 );
    CGAL_precondition( t2 >= 0 );
    CGAL_precondition( t1 + t2 <= 1.0 );  // i.e. t0 >= 0
    _t1 = t1;
    _t2 = t2;
}


template <class Surface_>
inline
//Point_s<Surface_>::Point_s( const typename Surface_::Vertex_const_handle v )
Point_s<Surface_>::Point_s( Vertex_const_handle v )
{
    // Need to ensure that the halfedge's facet is interior (not a hole).
    Halfedge_const_handle h = v->halfedge();
    if ( h->is_border() )
	h = h->next_on_vertex();
    frame = Facet_coordinate_frame<Surface>(h->next());
    _t1 = _t2 = 0;
}


/*!
 * \pre \a p lies in \code h->facet()
 */
template <class _surface>
inline
Point_s<_surface>::Point_s( Halfedge_const_handle h,
			    const Point_3& p )
    : frame(h)
{
    // Compute projection of p to the plane of frame(h).
    //
    frame.coordinates_of( p, &_t1, &_t2 );
    MNI::round_areal_coordinates( _t1, _t2 );

    // Check that the surface point agrees with p.
    //
    double mismatch_distance 
	= CGAL::to_double( CGAL::squared_distance( p, 
						   frame.point(_t1,_t2) ));
    if ( mismatch_distance >= 1e-8 )
	std::cerr << "Mismatch by: " << sqrt(mismatch_distance) << std::endl;
	
    CGAL_postcondition( mismatch_distance < 2e-6 );
}



template <class _Surface>
inline typename _Surface::Point_3
Point_s<_Surface>::point() const
{
    return frame.point( t1(),t2() );
}


template <class _Surface>
inline void
Point_s<_Surface>::coordinates_frame_12( double* u, double* v ) const
{
    *u = t2();  *v = t0();
}


template <class _Surface>
inline Point_s<_Surface>
Point_s<_Surface>::point_frame_12() const
{
    double u,v;
    coordinates_frame_12(&u,&v);
    
    /* You might think that no rounding is necessary.
     * However, observe that v = t0() = 1 - t1_ - t2_,
     * and this result might need rounding.
     */
    MNI::round_areal_coordinates(u,v);
    return Self(halfedge()->next(),u,v);
}


template <class _Surface>
inline void
Point_s<_Surface>::coordinates_frame_20( double* u, double* v ) const
{
    *u = t0();  *v = t1();
}


template <class _Surface>
inline Point_s<_Surface>
Point_s<_Surface>::point_frame_20() const
{
    double u,v;
    coordinates_frame_20(&u,&v);
    MNI::round_areal_coordinates(u,v);
    return Self(halfedge()->next()->next(),u,v);
}


/*! \brief Areal coordinates using different base edge.
 *
 * Set \a (u,v) to be the areal coordinates of this point, from
 * a different frame of reference.
 *
 * \pre \a g must be halfedge incident to the facet in which the point lies.
 */
template <class _Surface>
inline void
Point_s<_Surface>::coordinates_set_halfedge( Halfedge_const_handle g,
					     double* u, double* v ) const
{
    CGAL_precondition( &* g->facet() == &* halfedge()->facet() );

    if ( &* g == &* halfedge() ) {
	*u = t1();  *v = t2();
    } else if ( &* g == &* halfedge()->next() ) {
	*u = t2();  *v = t0();
    } else {
	*u = t0();  *v = t1();
    }
}



namespace detail {

/*! \brief Clip the ray p + t*v with the facet containing p.
 *
 * \note \a stepsize is passed by reference and modified by
 * this function.
 */
inline int
clip_ray( double t1, double t2, double dt1, double dt2,
	  double& stepsize )
{
    MNI::round_areal_coordinates(t1,t2);
    double t0 = 1 - t1 - t2;

    double dt0 = - dt1 - dt2;
    int active_constraint = -1;

    /* Figure out how far we can step without violating a constraint.
     * Since we maintain the invariant t0+t1+t2 = 1, we cannot have 
     * ti > 1 without some tj < 0.  So we only check the ti=0
     * constraints.
     */
	    
    if ( dt0 < 0 ) {
	double max_step = -t0/dt0;
	if ( max_step < stepsize ) {
	    stepsize = max_step;
	    active_constraint = 0;
	}
    }
    if ( dt1 < 0 ) {
	double max_step = -t1/dt1;
	if ( max_step < stepsize ) {
	    stepsize = max_step;
	    active_constraint = 1;
	}
    }
    if ( dt2 < 0 ) {
	double max_step = -t2/dt2;
	if ( max_step < stepsize ) {
	    stepsize = max_step;
	    active_constraint = 2;
	}
    }

    return active_constraint;
}

}


/*! \brief Displace point along surface.
 *
 * Return Point_s(frame,u1,u2) + Vector_s(frame,du1,du2).
 */
template <class Surface_>
inline const Point_s<Surface_>
add_vector( Facet_coordinate_frame<Surface_> frame,
	    double u1, double u2,
	    double du1, double du2 )
{
    typedef Facet_coordinate_frame<Surface_>       Frame;
    typedef typename Surface_::Point_3             Point_3;
    typedef typename Surface_::Traits::Vector_3    Vector_3;

    //++add_vector_invocation_count;

    while (1) {
	CGAL_assertion( finite(u1) );
	CGAL_assertion( finite(u2) );
	CGAL_assertion( finite(du1) );
	CGAL_assertion( finite(du2) );

	double step = 1.0;
	int active_constraint =  detail::clip_ray( u1,u2, du1,du2, step );

	if ( active_constraint == -1 )
	    break;

	//++add_vector_cross_count;

	//Point_3 p = frame.point( u1, u2 );
	Point_3 q = frame.point( u1+step*du1, u2+step*du2 );
	Point_3 r = frame.point( u1+du1, u2+du2 );

	if ( active_constraint == 0 ) 
	    frame = Frame(frame.halfedge()->next()->opposite());
	else if ( active_constraint == 1 )
	    frame = Frame(frame.halfedge()->next()->next()->opposite());
	else
	    frame = Frame(frame.halfedge()->opposite());

	Vector_3 rq = r - q;
	Vector_3 v10 = frame.P1() - frame.P0();
	Vector_3 v20 = frame.P2() - frame.P0();

	//frame.coordinates_of( q, &u1,&u2 );
	u1 = CGAL::to_double( (q - frame.P0())*v10 / (v10*v10) );
	u2 = 0;

	double rq_len = sqrt( CGAL::to_double(rq*rq) );
	double v10_len = sqrt( CGAL::to_double(v10*v10) );
	double v20_len = sqrt( CGAL::to_double(v20*v20) );

	double cos_theta = CGAL::to_double( v10*v20 / (v10_len*v20_len) );
	double sin_theta = sqrt(1-cos_theta*cos_theta);

	// Compute components of vector rq in new frame:
	// rq = A1 e1 + A2 e2, 
	// where e1 = v10/v10_len and e2 is perpendicular (also unit length)
	double A1 = CGAL::to_double( (rq*v10) / v10_len);
	double A2_sq = rq_len*rq_len - A1*A1;
	if ( A2_sq < 0 ) {
	    CGAL_assertion( A2_sq > -1e-8 );
	    A2_sq = 0;
	}
	double A2 = sqrt( A2_sq );

	du2 = A2 / (v20_len*sin_theta);
	du1 = (A1 - v20_len*cos_theta*du2) / v10_len;
    }

    u1 += du1;
    u2 += du2;

    MNI::round_areal_coordinates(u1,u2);
    return Point_s<Surface_>( frame.halfedge(), u1, u2 );
}


template <class Surface_>
inline const Point_s<Surface_>
Point_s<Surface_>::displace( double du1, double du2 ) const
{
    return add_vector( coordinate_frame(), t1(),t2(), du1,du2 );
}


template <class _surface>
inline std::ostream& operator<< ( std::ostream& os, 
				  const Point_s<_surface>& p )
{
    os << "Halfedge: " 
       << "[" << p.halfedge()->opposite()->vertex()->point()
       << ", " << p.halfedge()->vertex()->point() 
       << "]";
    //operator<< <_surface>(os, *p.halfedge());
    return os << "  Coordinates: (" << p.t0()
	      << ", " << p.t1()
	      << ", " << p.t2()
	      << ")";
}


template <class _surface>
inline bool
is_incident_to_facet( const Point_s<_surface>& p,
		      typename _surface::Facet_const_handle f_handle )
{
    const typename _surface::Facet* f = &* f_handle;
    if ( &* p.halfedge()->facet() == f )
	return true;

    const int constraint_0 = 0x0101;
    const int constraint_1 = 0x0201;
    const int constraint_2 = 0x0401;

    int constraints = 0;

    if ( p.t0() == 0 ) {
	constraints += constraint_0;
	if ( &* p.halfedge()->next()->opposite()->facet() == f )
	    return true;
    }

    if ( p.t1() == 0 ) {
	constraints += constraint_1;
	if ( &* p.halfedge()->next()->next()->opposite()->facet() == f )
	    return true;
    }

    if ( p.t2() == 0 ) {
	constraints += constraint_2;
	if ( &* p.halfedge()->opposite()->facet() == f )
	    return true;
    }

    if ( (constraints & 0x000f) < 2 )
	return false;


    // At this point, at least two constraints were active.  It is an
    // error for three constraints to be active, so we are at a
    // vertex.  Check all facets incident to vertex.

    typename _surface::Vertex_const_handle v = p.coordinate_frame().V0();

    if ( constraints == 0x0502 )
	v = p.coordinate_frame().V1();
    else if ( constraints == 0x0302 )
	v = p.coordinate_frame().V2();
    else
	CGAL_assertion( constraints == 0x0602 );

    typename _surface::Halfedge_const_handle h = v->halfedge();
    do {
	if ( &* h->facet() == f )
	    return true;
	h = h->next_on_vertex();
    } while ( h != v->halfedge() );
 
    return false;
}
	


template <class _surface>
inline typename _surface::Halfedge_const_handle
find_common_facet( const Point_s<_surface>& p,
		   const Point_s<_surface>& q )
{
    // Test q against all the facets incident to p.

    if ( is_incident_to_facet(q,p.halfedge()->facet()) )
	return p.halfedge();

    const int constraint_0 = 0x0101;
    const int constraint_1 = 0x0201;
    const int constraint_2 = 0x0401;

    int constraints = 0;
    typename _surface::Halfedge_const_handle h;

    if ( p.t0() == 0 ) {
	h = p.halfedge()->next()->opposite();
	constraints += constraint_0;
	if ( is_incident_to_facet(q,h->facet()) )
	    return h;
    }

    if ( p.t1() == 0 ) {
	h = p.halfedge()->next()->next()->opposite();
	constraints += constraint_1;
	if ( is_incident_to_facet(q,h->facet()) )
	    return h;
    }

    if ( p.t2() == 0 ) {
	h = p.halfedge()->opposite();
	constraints += constraint_2;
	if ( is_incident_to_facet(q,h->facet()) )
	    return h;
    }

    if ( constraints & 0x000f < 2 )
	return 0;


    // At this point, at least two constraints were active.  It is an
    // error for three constraints to be active, so we are at a
    // vertex.  Check all facets incident to vertex.

    typename _surface::Vertex_const_handle v = p.coordinate_frame().V0();

    if ( constraints == 0x0502 )
	v = p.coordinate_frame().V1();
    else if ( constraints == 0x0302 )
	v = p.coordinate_frame().V2();
    else
	CGAL_assertion( constraints == 0x0602 );

    h = v->halfedge();
    do {
	if ( is_incident_to_facet(q,h->facet()))
	    return h;
	h = h->next_on_vertex();
    } while ( h != v->halfedge() );
 
    return 0;
}

}

