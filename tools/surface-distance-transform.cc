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

#include <surflib/DistanceTransform.hpp>

#include "Surface.hpp"

typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;



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

    try {
	Surface s;
	MNI::load_surface_with_scalar( s, av[1], av[2] );


	// Generate list of all seed vertices: those with
	// nonzero scalar value.
	//
	std::vector<Surface::Vertex_handle> seed_list;

	Surface::Vertex_iterator v = s.vertices_begin();
	CGAL_For_all( v, s.vertices_end() ) {
	    if ( (seed_label == -1 && v->scalar) ||
		 (seed_label == v->scalar) )
		seed_list.push_back(v);
	}

	cout << seed_list.size() << " seed vertices." << endl;
	MNI::distanceTransform( s,
				seed_list.begin(),
				seed_list.end(),
				extra_nodes_per_edge );

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

    } catch ( const std::bad_alloc& e ) {
	std::cerr << "Failed a memory allocation." << "\n"
		  << "No output." << "\n";
	return 2;
    } catch ( const std::exception& e ) {
	std::cerr << "Error: " << e.what() << "\n"
		  << "No output." << "\n";
        return 3;
    } catch ( ... ) {
	std::cerr << "Unknown exception." << "\n"
		  << "No output." << "\n"
		  << "This is likely bug in the code: please report!" << "\n";
        return 4;
    }

    return 0;
}
