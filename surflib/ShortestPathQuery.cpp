// -*- C++ -*-

#include <iostream>
#include <surflib/vtk_ostream.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>


namespace MNI {


/*!  Set up the query structure with 
 *   a number of additional vertices on each edge. 
 */
template <class _surface>
ShortestPathQuery<_surface>::ShortestPathQuery( const _surface& s,
						int n )
    : surface(s),
      graph( s.size_of_vertices() + n*s.size_of_halfedges()/2 + 2 ),
      point_of_vertex( s.size_of_vertices() + n*s.size_of_halfedges()/2 + 2),
      search_length(40,0,800)
{
    vertices_per_edge = n;
    additional_vertex_offset = surface.size_of_vertices();

    CGAL_precondition( vertices_per_edge >= 0 );

    // Set up the maps.
    
    Vertex_const_iterator v = surface.vertices_begin();
    for( int i = 0; v != surface.vertices_end(); ++v,++i ) {
	vertex_index_map[v] = i;
	point_of_vertex[i] = Point_s(v);
    }
    point_of_vertex[search_source()] = point_of_vertex[0];
    point_of_vertex[search_target()] = point_of_vertex[0];

    typename Surface::Edge_const_iterator e = surface.edges_begin();
    for( int i = 0; e != surface.edges_end(); ++e,++i ) {
	Halfedge_const_handle h = e;
	if ( h->is_border() )
	    h = h->opposite();

	halfedge_index_map[h] = i;
	halfedge_index_map[h->opposite()] = i;

	for( int k = 0; k < vertices_per_edge; ++k ) {
	    int index = vertex_index(h,k);
	    CGAL_assertion( index >= additional_vertex_offset );
	    double t1 = double(k+1)/double(1+vertices_per_edge);
	    point_of_vertex[index] = Point_s(h,t1,0);
	}
    }

    // Add weighted edges to search graph.

    e = surface.edges_begin();
    CGAL_For_all(e,surface.edges_end()) {
	Halfedge_const_handle h = e;

	// Insert edges that subdivide surface halfedge h.
	if ( vertices_per_edge == 0 ) {
	    add_edge( vertex_index(h->opposite()->vertex()),
		      vertex_index(h->vertex()) );
	} else {
	    add_edge( vertex_index(h->opposite()->vertex()),
		      vertex_index(h,0) );
	    for( int k = 0; k < vertices_per_edge-1; ++k )
		add_edge( vertex_index(h,k),
			  vertex_index(h,k+1) );
	    add_edge( vertex_index(h,vertices_per_edge-1),
		      vertex_index(h->vertex()) );
	}
    }

    Halfedge_const_iterator h = surface.halfedges_begin();
    CGAL_For_all(h,surface.halfedges_end()) {
	// Insert edges between h and the adjacent halfedge in its facet.
	Halfedge_const_handle g = h->next();
	for( int k = 0; k < vertices_per_edge; ++k )
	    add_edge_set( vertex_index(h,k), g );
    }
}


template <class _surface>
ShortestPathQuery<_surface>::~ShortestPathQuery()
{
    if (debug) {
	std::cerr << "Query statistics: " << visits << std::endl;
	std::cerr << "Histogram of search lengths." << std::endl;
	search_length.dump(std::cerr);
    }
}


template <class _surface>
double
ShortestPathQuery<_surface>::find_path( const Point_s& s, const Point_s& t,
					std::vector<Point_s>& pathlist )
{
    if ( &* s.halfedge()->facet() == &* t.halfedge()->facet() ) {
	// Treat this case special, because we do not add the edge
	// (s,t) to the search graph, and hence would never return the
	// correct answer otherwise.
	pathlist.push_back( s );
	pathlist.push_back( t );
	return sqrt(CGAL::to_double(squared_distance(s,t)));
    }

    vertex_descriptor vs = search_source();
    vertex_descriptor vt = search_target();

    set_terminal_point( vs, s );
    set_terminal_point( vt, t );

    std::vector<vertex_descriptor> pred(boost::num_vertices(graph));
    std::vector<float> dist(boost::num_vertices(graph));

    int count = 0;
    bool found = false;
    try {
	boost::dijkstra_shortest_paths
	    (graph, vt,
	     boost::distance_map(&dist[0]).
	     predecessor_map(&pred[0]).
	     visitor(boost::make_dijkstra_visitor(my_visitor(vs,count))));
    } catch (FoundTarget& ft) {
	found = true;
    }
    CGAL_assertion(found);
    visits.add_sample( count );
    search_length.accumulate(count);

    vertex_descriptor v = vs;
    while (v != vt) {
	CGAL_assertion( pred[v] != v );
	pathlist.push_back( vertex_location(v) );
	v = pred[v];
    }
    pathlist.push_back( vertex_location(v) );

    return dist[vs];
}


template <class _surface>
Point_s<_surface>
ShortestPathQuery<_surface>::pathpoint_distance
( const std::vector<Point_s>& path,
  double distance )
{
    CGAL_precondition( path.size() > 0 );
    if (path.size() == 1 || distance == 0)
	return path[0];
    
    typename std::vector<Point_s>::const_iterator pi = path.begin();
    typename std::vector<Point_s>::const_iterator pi_prev = pi;
    double cumulative_dist = 0;
    double segment_length = 0;

    // cumulative_dist is the path length from path.begin() to pi.
    // pi_prev = pi - 1 (after 1 iteration)
    // segment_length is distance(pi_prev,pi)
    while ( cumulative_dist < distance ) {
	pi_prev = pi++;
	if ( pi == path.end() ) {
	    // This isn't supposed to happen.  
	    // However, due to roundoff, it does.  Sigh.
	    CGAL_assertion( distance - cumulative_dist < 1e-5*distance );
	    return path.back();
	}
	segment_length = sqrt(CGAL::to_double(squared_distance(*pi_prev,*pi)));
	cumulative_dist += segment_length;
    }
    CGAL_assertion( segment_length > 0 );

    // The point sought is (1-a)(pi_prev) + a(pi), where
    // a = (distance - cumulative_dist(pi_prev))/segment_length
    // note that cumulative_dist(pi_prev) = cumulative_dist - segment_length,
    // so a = (distance - cumulative_dist + segment_length)/segment_length
    // 
    // Let b = 1 - a = (cumulative_dist - distance)/segment_length, then
    // point = b(pi_prev) + (1-b)(pi) = (pi) + b(pi_prev - pi)

    /* Should be able to do this on the surface ...
    Vector_s d = (*pi) - (*pi_prev);
    Point_s p = *pi_prev + d*((cumulative_dist-distance)/segment_length);
    return p;
    */

    // Find the common facet containing pi and pi_prev to use as
    // computation frame.
    Halfedge_const_handle h = find_common_facet( *pi_prev, *pi );
    if ( h == 0 ) {
	dump_path( path );
	std::cerr << "pi_prev = " << *pi_prev << std::endl
	          << "pi = " << *pi << std::endl;
	CGAL_assertion_msg(0,"Did not find common facet");
    }
    typename Surface::Traits::Vector_3 d = pi_prev->point() - pi->point();
    Point_3 p = pi->point() + d*((cumulative_dist-distance)/segment_length);
    return Point_s(h,p);
}


template <class _surface>
inline Point_s<_surface>
ShortestPathQuery<_surface>::pathpoint_distance( const Point_s& s, 
						 const Point_s& t,
						 double distance )
{
    CGAL_precondition( distance >= 0 );

    std::vector<Point_s> path;
    double pathlength = find_path(s,t,path);

    CGAL_precondition( distance <= pathlength );

    return pathpoint_distance( path, distance );
}


template <class _surface>
inline Point_s<_surface>
ShortestPathQuery<_surface>::pathpoint_fraction( const Point_s& s, 
						 const Point_s& t,
						 double fraction )
{
    CGAL_precondition( fraction >= 0 );
    CGAL_precondition( fraction <= 1 );

    std::vector<Point_s> path;
    double pathlength = find_path(s,t,path);
    return pathpoint_distance( path, fraction*pathlength );
}


template <class _surface>
inline Point_s<_surface>
ShortestPathQuery<_surface>::midpoint( const Point_s& s, const Point_s& t )
{
    return pathpoint_fraction(s,t,0.5);
}


template <class _surface>
inline int
ShortestPathQuery<_surface>::vertex_index( Vertex_const_handle v ) const
{
    return vertex_index_map[v];
}



template <class _surface>
inline int
ShortestPathQuery<_surface>::vertex_index( Halfedge_const_handle h, 
					   int i ) const
{
    return additional_vertex_offset 
	+ vertices_per_edge * halfedge_index_map[h] + i;
}



template <class _surface>
inline Point_s<_surface>
ShortestPathQuery<_surface>::vertex_location( int index ) const
{
    return point_of_vertex[index];
}



template <class _surface>
inline void 
ShortestPathQuery<_surface>::add_edge( vertex_descriptor u,
				       vertex_descriptor v,
				       float weight )
{
    boost::add_edge(u,v,weight,graph);
}


template <class _surface>
inline void 
ShortestPathQuery<_surface>::add_edge( vertex_descriptor u,
				       vertex_descriptor v )
{
    Point_3 pu( vertex_location_3(u) );
    Point_3 pv( vertex_location_3(v) );
    double sq_dist = CGAL::to_double( CGAL::squared_distance(pu,pv) );
    add_edge(u,v,sqrt(sq_dist));
}



template <class _surface>
inline void 
ShortestPathQuery<_surface>::add_edge_set( vertex_descriptor u,
					   Halfedge_const_handle h )
{
    add_edge( u, vertex_index(h->vertex()) );
    for( int k = 0; k < vertices_per_edge; ++k ) 
	add_edge( u, vertex_index(h,k) );
}


template <class _surface>
void
ShortestPathQuery<_surface>::write_vtk( const char* filename ) const
{
    using namespace boost;

    std::ofstream vtk(filename);
    write_vtk_header( vtk, "ShortestPathQuery edges" );

    vtk << "POINTS " << num_vertices(graph) << " double" << endl;
    graph_traits<Graph>::vertex_iterator vi,vi_end;
    for( tie(vi,vi_end) = vertices(graph); vi != vi_end; ++vi ) {
	Point_3 p( vertex_location_3(*vi) );
	vtk << CGAL::to_double(p.x()) << ' '
	    << CGAL::to_double(p.y()) << ' '
	    << CGAL::to_double(p.z()) << endl;
    }

    vtk << "LINES " << num_edges(graph) << ' ' << 3*num_edges(graph) << endl;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
	graph_traits<Graph>::edge_descriptor e = *ei;
	graph_traits<Graph>::vertex_descriptor u = source(e, graph);
	graph_traits<Graph>::vertex_descriptor v = target(e, graph);
	vtk << "2 " << u << ' ' << v << endl;
    }
}

template <class _surface>
void dump_halfedge( const char* name, 
		    typename _surface::Halfedge_const_handle h )
{
    std::cerr << "    on halfege " << name 
	      << " (" << &* h << ")" << std::endl;
    std::cerr << "    Incident to facets "
	      << &* h->facet() << " and "
	      << &* h->opposite()->facet() << std::endl;
}


template <class _surface>
void dump_vertex( const char* name, 
		  typename _surface::Vertex_const_handle v )
{
    std::cerr << "    at vertex " << name 
	      << " (" << &* v << ")" << std::endl;
    std::cerr << "    Incident to facets ";
    typename _surface::Halfedge_around_vertex_const_circulator h 
	= v->vertex_begin();
    CGAL_For_all(h,v->vertex_begin()) {
	std::cerr << &* h->facet() << " ";
    }
    std::cerr << std::endl;
}


template <class _surface>
void
ShortestPathQuery<_surface>::dump_path( const std::vector<Point_s>& pathlist )
{
    int i = 0;
    typename std::vector<Point_s>::const_iterator pi = pathlist.begin();
    CGAL_For_all(pi,pathlist.end()) {
	std::cerr << "Point " << i++
		  << ": " << *pi 
		  << std::endl;

	const int constraint_0 = 0x0101;
	const int constraint_1 = 0x0201;
	const int constraint_2 = 0x0401;

	int constraints = 0;
	if ( (*pi).t0() == 0 )  constraints += constraint_0;
	if ( (*pi).t1() == 0 )  constraints += constraint_1;
	if ( (*pi).t2() == 0 )  constraints += constraint_2;

	switch(constraints) {
	case 0x0101:  
	    dump_halfedge<_surface>("12",pi->halfedge()->next()); 
	    break;
	case 0x0201:
	    dump_halfedge<_surface>("20",pi->halfedge()->next()->next()); 
	    break;
	case 0x0401:
	    dump_halfedge<_surface>("01",pi->halfedge()); 
	    break;
	case 0x0302:
	    dump_vertex<_surface>("2",pi->V2());
	    break;
	case 0x0502:
	    dump_vertex<_surface>("1",pi->V1());
	    break;
	case 0x0602:
	    dump_vertex<_surface>("0",pi->V0());
	    break;
	default: CGAL_assertion(constraints == 0);
	}
    }
}


template <class _surface>
void
ShortestPathQuery<_surface>::set_terminal_point( vertex_descriptor v, 
						 const Point_s& p )
{
    point_of_vertex[v] = p;
    boost::clear_vertex( v, graph );
    add_edge_set( v, p.halfedge() );
    add_edge_set( v, p.halfedge()->next() );
    add_edge_set( v, p.halfedge()->next()->next() );
}


}
