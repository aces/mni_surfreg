// -*- C++ -*-

namespace MNI {

template <class Vertex, class VertexFunctor>
inline
typename VertexWalkSearch<Vertex,VertexFunctor>::Vertex_const_handle
VertexWalkSearch<Vertex,VertexFunctor>::search( Vertex_const_handle v )
{
    double best_value = f(v);
    Vertex_const_handle best_vertex = v;

    while(1) {
	typename Vertex::Halfedge_around_vertex_const_circulator
	    hc = v->vertex_begin();
	CGAL_For_all( hc, v->vertex_begin() ) {
	    //typename Vertex::Halfedge_const_handle h = hc;
	    double value = f( hc->opposite()->vertex() );
            // This test is safer on RedHat to give good results (Claude).
            double diff = value - best_value;
            if ( 1.0 + diff < 1.0 ) {
		best_value = value;
		best_vertex = hc->opposite()->vertex();
	    }
	}
	if (v == best_vertex)  
	    break;
	v = best_vertex;
    }
    return best_vertex;
}


template <class Surface_>
inline const Point_s<Surface_> 
find_nearest_on_convex_surface( typename Surface_::Vertex_const_handle v,
				const typename Surface_::Point_3& p )
{
    typedef typename Surface_::Vertex  Vertex;

    // 1. Walk to the nearest vertex.
    //
    SquaredDistanceToPoint<Vertex> search_f(p);
    VertexWalkSearch<Vertex, SquaredDistanceToPoint<Vertex> >
	searcher( search_f );

    Point_s<Surface_> p_s = searcher.search( v );

    // 2. Find correction to p_s that takes us closer to p.
    //
    double dt1, dt2;
    p_s.coordinate_frame().coordinates_of( p - p_s.point(),
					   &dt1,&dt2 );
    return p_s.displace( dt1, dt2 );
}


template <class Surface_>
const Point_s<Surface_>
find_midpoint_on_sphere( const Point_s<Surface_>& p,
			 const typename Surface_::Point_3& q )
{
    typedef typename Surface_::Traits::Vector_3 Vector_3;

    // 1. Compute mid_3 = midpoint in 3-space.
    //
    Vector_3 vp( p.point() - CGAL::ORIGIN );
    Vector_3 vq( q - CGAL::ORIGIN );

    vp = vp / std::sqrt( CGAL::to_double(vp*vp) );
    vq = vq / std::sqrt( CGAL::to_double(vq*vq) );

    double cos_theta = vp*vq;

    typename Surface_::Point_3 mid_3 
	= rotate( p.point(),                   // point to rotate
		  CGAL::cross_product(vp,vq),  // axis
		  acos(cos_theta)/2.0          // rotation angle
		  );

    // 2. Find mid_s = target mesh vertex nearest to mid_3.
    //
    return find_nearest_on_convex_surface<Surface_>( p.V0(), mid_3 );
}


}

