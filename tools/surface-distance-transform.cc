/* Compute distance transform on surface.
 */

#include <iostream>
#include <iterator>
#include <stdexcept>

#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

#include <surflib/ShortestPathQuery.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>

/* This code depends on having points represented using
 * handles, as we use the handle to identify instances of
 * the same point.  More precisely, the handle is used to
 * map a point of the surface (and of the convex hull) to
 * an integer index.
 */
#include "Surface.hpp"
typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;



using namespace std;


typedef Surface::Vertex_handle       Vertex_handle;
typedef MNI::ShortestPathQuery<Surface>::Graph     Graph;
typedef MNI::ShortestPathQuery<Surface>::vertex_descriptor  vertex_descriptor;



void dump_graph( const char* filename, 
		 MNI::ShortestPathQuery<Surface>& sp )
{
    sp.write_vtk(filename);
}


/* Compute depth in "scalar" element of the surface points.
 * The seed vertices (seed_list) are set to depth 0, the others
 * are assigned the length of the shortest path in a graph augmented
 * with extra_nodes_per_edge.
 */
void compute_distance( Surface& s, 
		       const std::vector<Vertex_handle>& seed_list,
		       int extra_nodes_per_edge )
{
    using namespace boost;

    cout << "Computing query graph ... " << flush;
    MNI::ShortestPathQuery<Surface> query( s, extra_nodes_per_edge );
    cout << "done." << endl;

    const Graph& g(query.get_graph());

    // For debugging.
    //dump_graph( "debug.vtk", query );

    // Add one source vertex to query graph, with zero-weight
    // edges to each seed vertex.

    // HACK!  I know that there were two vertices reserved in graph g
    // for start & terminus.
    //vertex_descriptor start = boost::vertex(boost::num_vertices(g)-1, g);

    vertex_descriptor start = boost::vertex( query.search_source(), g);

    std::vector<Vertex_handle>::const_iterator vi = seed_list.begin();
    CGAL_For_all(vi,seed_list.end()) {
	vertex_descriptor seed = query.vertex_index(*vi);
	//cout << "Seed: " << seed << endl;
	query.add_edge( start, seed, 0 );
    }

    std::vector<float> dist(num_vertices(g));

    cout << "Starting shortest-paths computation." << endl
   	 << "Graph has " << num_vertices(g)
         << " vertices and " << num_edges(g)
         << " edges .... " << flush;

    dijkstra_shortest_paths(g, start, distance_map(&dist[0]));

    cout << "done." << endl;

    // Copy distance information back into vertex scalar of s.

    cout << "Storing depths ..." << flush;
    Surface::Vertex_iterator v = s.vertices_begin();
    for( ; v != s.vertices_end(); ++v ) {
	v->scalar = dist[ query.vertex_index(v) ];
    }
    cout << "done." << endl;
}


int main( int ac, char* av[] )
{
    char* vtk_file = 0;
    int extra_nodes_per_edge = 0;

    // Label value that identifies a seed vertex.  
    // -1 indicates "any nonzero value".
    int seed_label = -1;

    ArgvInfo argTable[] = {
	{ "-vtk", ARGV_STRING, (char*)1, (char*)&vtk_file,
	  "write output in VTK format" },
	{ "-extra", ARGV_INT, (char*)1, (char*)&extra_nodes_per_edge,
	  "number of extra nodes on each edge" },
	{ "-label", ARGV_INT, (char*)1, (char*)&seed_label,
	  "label value that identifies a seed (-1 -> any nonzero value)" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };


    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 4 ) {
	char* s = strrchr( av[0], '/' );
	std::cerr << "usage: " << (s ? s+1 : av[0]) 
		  << " [options] surface_geom surface_seeds distance.vv"
		  << std::endl;
        return 1;
    }

    Surface s;
    MNI::load_surface_with_scalar( s, av[1], av[2] );


    // Generate list of all seed vertices: those with
    // nonzero scalar value.
    //
    std::vector<Vertex_handle> seed_list;

    Surface::Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
	if ( (seed_label == -1 && v->scalar) ||
	     (seed_label == v->scalar) )
	    seed_list.push_back(v);
    }

    cout << seed_list.size() << " seed vertices." << endl;
    compute_distance(s,seed_list,extra_nodes_per_edge);

    if ( vtk_file ) {
	// Write in VTK format.
	// Note: assuming that write_vtk() outputs vertices
	// in order [vertices_begin,vertices_end).
	std::ofstream out( vtk_file );
	MNI::write_vtk( out, s, 
			"Distance transform" );
	out << "POINT_DATA " << s.size_of_vertices() << std::endl
	    << "SCALARS depth double" << std::endl
	    << "LOOKUP_TABLE default" << std::endl;
	for( Surface::Vertex_iterator v = s.vertices_begin();
	     v != s.vertices_end(); ++v )
	    out << v->scalar << std::endl;
    }

    // Dump values to file in vertex order.
    std::ofstream out( av[3] );
    v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
	out << v->scalar << std::endl;
    }

    return 0;
}
