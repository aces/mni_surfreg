#ifndef MNI_SURFLIB_DISTANCE_TRANSFORM_HPP   // -*- C++ -*-
#define MNI_SURFLIB_DISTANCE_TRANSFORM_HPP


#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

#include <surflib/ShortestPathQuery.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>


namespace MNI {


/**
 * Compute a distance transform at vertices of a polyhedral surface.
 * 
 * The template type Surface must be a CGAL::Polyhedron_3 subclass.
 * In addition, the Vertex class must have a double field named
 * scalar.
 *
 * After calling this method, the scalar field contains the distance
 * transform value.
 *
 * Template type VertexHandleIterator must be a forward iterator
 * with value_type of Surface::Vertex_handle.  The range
 * [seed_begin,seed_end) is the set of seed vertices; each of
 * which is a vertex of Surface s.
 *
 * The distance transform computation uses the graph of Surface s,
 * augmented with extra_nodes_per_edge extra nodes per edge.
 */

template< class Surface, class VertexHandleIterator >
void distanceTransform( Surface& s, 
			VertexHandleIterator seed_begin,
			VertexHandleIterator seed_end,
			int extra_nodes_per_edge )
{
    typedef typename Surface::Vertex_handle       Vertex_handle;
    typedef typename MNI::ShortestPathQuery<Surface>::Graph     Graph;
    typedef typename MNI::ShortestPathQuery<Surface>::vertex_descriptor  vertex_descriptor;

    std::cout << "Computing query graph ... " << std::flush;
    MNI::ShortestPathQuery<Surface> query( s, extra_nodes_per_edge );
    std::cout << "done.\n";

    const Graph& g(query.get_graph());

    // For debugging.
    // query.write_vtk( "debug.vtk" );

    // Add one source vertex to query graph, with zero-weight
    // edges to each seed vertex.

    vertex_descriptor start = boost::vertex( query.search_source(), g);

    VertexHandleIterator vi = seed_begin;
    CGAL_For_all( vi, seed_end ) {
	vertex_descriptor seed = query.vertex_index(*vi);
	query.add_edge( start, seed, 0 );
    }

    std::vector<float> dist(num_vertices(g));

    std::cout << "Starting shortest-paths computation.\n"
	      << "Graph has " << num_vertices(g)
	      << " vertices and " << num_edges(g)
	      << " edges .... " << std::flush;

    dijkstra_shortest_paths(g, start, boost::distance_map(&dist[0]));

    std::cout << "done.\n";

    // Copy distance information back into vertex scalar of s.

    std::cout << "Storing depths ..." << std::flush;
    typename Surface::Vertex_iterator v = s.vertices_begin();
    for( ; v != s.vertices_end(); ++v ) {
	v->scalar = dist[ query.vertex_index(v) ];
    }
    std::cout << "done.\n";
}




}


#endif
