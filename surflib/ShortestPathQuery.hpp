#ifndef SHORTEST_PATH_QUERY_HPP
#define SHORTEST_PATH_QUERY_HPP

/* A class that supports shortest path queries on a surface.
 */

#include <vector>
#include <CGAL/Unique_hash_map.h>


#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>

#include <surflib/Statistic.hpp>
#include <surflib/Histogram.hpp>
#include <surflib/Point_s.hpp>
#include <surflib/Vector_s.hpp>


namespace MNI {


template <class _surface>
class ShortestPathQuery
{
public:
    typedef _surface Surface;

    typedef typename Surface::Vertex_const_handle      Vertex_const_handle;
    typedef typename Surface::Vertex_const_iterator    Vertex_const_iterator;
    typedef typename Surface::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Surface::Halfedge_const_iterator  Halfedge_const_iterator;

    typedef Point_s<Surface>                           Point_s;
    typedef Vector_s<Surface>                          Vector_s;
    typedef typename Surface::Point_3                  Point_3;
 


    //! Construct query structure for surface s.
    ShortestPathQuery( const _surface& s, 
		       int additional_vertices_per_edge );

    ~ShortestPathQuery();
    

    //! Find path from s to t.  Returns length.
    double find_path( const Point_s& s, const Point_s& t,
		      std::vector<Point_s>& pathlist );

    //! Find point a given distance along path.
    Point_s pathpoint_distance( const std::vector<Point_s>& pathlist,
				double distance );

    //! Find point a given distance along path from s to t.
    Point_s pathpoint_distance( const Point_s& s, const Point_s& t,
				double distance );

    //! Find point a given fraction along path from s to t.
    Point_s pathpoint_fraction( const Point_s& s, const Point_s& t,
				double fraction );

    //! Return midpoint of path from s to t.
    Point_s midpoint( const Point_s& s, const Point_s& t);




    typedef boost::adjacency_list< boost::vecS, boost::vecS, 
				   boost::undirectedS,
				   boost::no_property, 
				   boost::property<boost::edge_weight_t,float> 
    > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor  vertex_descriptor;
    typedef boost::graph_traits<Graph>::edge_descriptor    edge_descriptor;



    //! Index in search graph.
    int vertex_index( Vertex_const_handle v ) const;
    int vertex_index( Halfedge_const_handle h, int i ) const;

    Point_s vertex_location( int index ) const;
    Point_3 vertex_location_3( int index ) const
	{ return vertex_location(index).point(); }


    //! Add one edge to search graph, weighted by Euclidean distance.
    void add_edge( vertex_descriptor u, vertex_descriptor v );

    //! Add edges from u to all vertices on h, Euclidean weight.
    void add_edge_set( vertex_descriptor u, Halfedge_const_handle h );


    //! Add one edge.
    void add_edge( vertex_descriptor u, vertex_descriptor v, float weight );


    //! Write query graph as VTK file.
    void write_vtk( const char* filename ) const;


    //! Enable debugging output.
    static bool debug;


    //! Dump debugging info about the path.
    void dump_path( const std::vector<Point_s>& pathlist );


    const Graph& get_graph() const 
	{ return graph; }

    vertex_descriptor search_source() const 
    {  return boost::num_vertices(graph) - 2; }

    vertex_descriptor search_target() const
    {  return boost::num_vertices(graph) - 1;  }


private:

    //! Connect this point to all the vertices around the facet.
    void set_terminal_point( vertex_descriptor v, const Point_s& p );
    
    class FoundTarget {};

    struct my_visitor 
    {
	typedef boost::on_examine_vertex  event_filter;

	vertex_descriptor t;
	int& visit_count;
	my_visitor( vertex_descriptor target, int& vc )
	    : t(target),visit_count(vc) {}

	void operator()( vertex_descriptor v, const Graph& g )
	{  ++visit_count;  if ( v == t )    throw FoundTarget();  }
    };


private:
    const Surface& surface;
    Graph graph;

    // Maps from surface to graph.
    CGAL::Unique_hash_map<Vertex_const_handle,int> vertex_index_map;
    CGAL::Unique_hash_map<Halfedge_const_handle,int> halfedge_index_map;

    // Map from graph to surface.
    std::vector<Point_s>  point_of_vertex;

    int additional_vertex_offset;
    int vertices_per_edge;

    MNI::Statistic<double> visits;
    MNI::Histogram<int> search_length;
};


template <class _surface> 
bool ShortestPathQuery<_surface>::debug = false;

}


#include "ShortestPathQuery.cpp"

#endif
