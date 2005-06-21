/* Compute mean and standard deviation for a set of co-registered
 * surfaces.
 *
    load the source mesh
    for each (target,map) pair:
	for each vertex v on source/control mesh:
	    accumulate the location smap[v]
    for each vertex v:
	generate mean and RMS values
    write output
*/


#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/circulator.h>

#include <surflib/VectorStatistic.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>
#include <surflib/SurfaceMap_hash.hpp>


typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;
typedef Surface::Halfedge_const_handle      Halfedge_const_handle;
typedef Surface::Vertex_const_handle        Vertex_const_handle;
typedef Surface::Vertex_const_iterator      Vertex_const_iterator;
typedef Surface::Facet_const_handle         Facet_const_handle;
typedef Surface::Facet_const_iterator       Facet_const_iterator;


typedef Surface_kernel::Point_3             Point_3;
typedef Surface_kernel::Vector_3            Vector_3;

typedef MNI::SurfaceMap_hash<Surface>       SurfaceMap;
typedef MNI::VectorStatistic<Point_3>       PointStatistic;



using namespace std;


bool is_map_filename( char* fn )
{
    int len = strlen(fn);
    return (len > 3) && (strcmp(&fn[len-3],".sm") == 0);
}


void accumulate( std::vector<PointStatistic>& points, 
		 const SurfaceMap& smap )
{
    SurfaceMap::ControlMesh::Vertex_const_iterator 
	v = smap.control_mesh().vertices_begin();

    for( int i = 0; v != smap.control_mesh().vertices_end(); ++i,++v ) {
	points[i].add_sample( smap[v].point() );
    }
}


struct point_rms: public unary_function<const PointStatistic&, double> 
{
    double operator()( const PointStatistic& ps ) 
    { return ps.std_dev(); }
};


void usage( char* argv0 )
{
    char* s = strrchr( argv0, '/' );
    std::cerr << "usage: " << (s ? s+1 : argv0) 
	      << " [options] source [map1] target1 [[map2] target2...]"
	      << std::endl
	      << "If mapfile is omitted, an identity map is assumed."
	      << std::endl;
}


int main( int ac, char* av[] )
{
    char* vtk_file = 0;
    char* rms_file = 0;

    static ArgvInfo argTable[] = {
	{ "-vtk", ARGV_STRING, (char*)1, (char*)&vtk_file,
          "write VTK file" },
	{ "-rms", ARGV_STRING, (char*)1, (char*)&rms_file,
          "write rms value file" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || (ac < 3) ) {
	usage( av[0] );
	return 1;
    }

    try {
	Surface source;
	MNI::load_surface_file( source, av[1] );

	std::vector<PointStatistic> meanpoints( source.size_of_vertices() );

	char* smap_fn = 0;
	for( int i = 2; i < ac; ++i ) {
	    if ( is_map_filename(av[i]) ) {
		if ( smap_fn != 0 )
		    throw std::runtime_error( "two consecutive surface map files not permitted" );
		smap_fn = av[i];
	    } else {
		Surface target;
		MNI::load_surface_file( target, av[i] );

		if ( smap_fn ) {
		    SurfaceMap smap( source, target, smap_fn );
		    accumulate( meanpoints, smap );
		} else {
		    SurfaceMap smap( source, target );
		    accumulate( meanpoints, smap );
		}

		smap_fn = 0;
	    }
	}
	if ( smap_fn != 0 )
	    throw std::runtime_error( "last argument is surfacemap" );

	if (vtk_file) {
	    ofstream out( vtk_file );
	    string id( "$Id$" );
	    MNI::write_vtk_header( out, id );

	    // Point data from meanpoints
	    out << "POINTS " << source.size_of_vertices() << " double" << endl;
	    vector<PointStatistic>::iterator i = meanpoints.begin();
	    for( ; i != meanpoints.end(); ++i ) {
		out << (*i).mean() << endl;
	    }
	    MNI::write_vtk_cells( out, source );
	    MNI::write_vtk_point_data_header( out, source );
	    MNI::write_vtk_point_scalars_header( out, "rms", "double" );
	    std::transform( meanpoints.begin(), meanpoints.end(),
			    std::ostream_iterator<double>(out,"\n"),
			    point_rms() );
	}

	if (rms_file) {
	    ofstream out( rms_file );
	    std::transform( meanpoints.begin(), meanpoints.end(),
			    std::ostream_iterator<double>(out,"\n"),
			    point_rms() );
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
