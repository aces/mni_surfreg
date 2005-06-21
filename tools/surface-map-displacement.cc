#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/circulator.h>

#include <ParseArgv.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>
#include <surflib/SurfaceMap.hpp>
#include <surflib/SurfaceFunctions.hpp>

#include <vector>


typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;
typedef Surface::Halfedge_const_handle      Halfedge_const_handle;
typedef Surface::Vertex_const_handle        Vertex_const_handle;
typedef Surface::Vertex_const_iterator      Vertex_const_iterator;
typedef Surface::Facet_const_handle         Facet_const_handle;
typedef Surface::Facet_const_iterator       Facet_const_iterator;


typedef Surface_kernel::Point_3             Point_3;
typedef Surface_kernel::Vector_3            Vector_3;

typedef MNI::SurfaceMap<Surface>            SMap;


using namespace std;



void scale_displacements( std::vector<Vector_3>& displacement, double s )
{
    std::vector<Vector_3>::iterator i = displacement.begin();
    CGAL_For_all(i,displacement.end()) {
	(*i) = s * (*i);
    }
}


void add_displacements( const SMap& map, 
			std::vector<Vector_3>& displacement )
{
    std::vector<Vector_3>::iterator i = displacement.begin();
    Vertex_const_iterator v = map.control_mesh().vertices_begin();

    for( ; v != map.control_mesh().vertices_end(); ++v,++i ) {
	Point_3 p0 = v->point();
	Point_3 p1 = map[v].point();
	*i = (*i) + (p1 - p0);
    }
}


/* Copy source mesh with a scalar or vector at each point.
 */
void write_mesh( const std::string filename, 
		 const Surface& s,
		 const std::vector<Vector_3>& displacement,
		 bool write_scalar,
		 bool write_vector )
{
    std::ofstream out( filename.c_str() );

    MNI::write_vtk( out, s, "Displacement Vectors." );

    out << "POINT_DATA " << s.size_of_vertices() << std::endl;

    if (write_scalar) {
	out << "SCALARS displacement_length double" << std::endl
	    << "LOOKUP_TABLE default" << std::endl;
	std::vector<Vector_3>::const_iterator i = displacement.begin();
	for( ; i != displacement.end(); ++i ) {
	    double sq = CGAL::to_double(*i * *i);
	    out << sqrt(sq) << std::endl;
	}
    }

    if (write_vector) {
	out << "VECTORS displacement double" << std::endl;
	std::ostream_iterator<Vector_3> outiter( out, "\n" );
	std::copy( displacement.begin(), displacement.end(), outiter );
    }
}


void write_vv( const char* filename, 
	       const std::vector<Vector_3>& displacement,
	       bool write_scalar,
	       bool write_vector )
{
    std::ofstream out( filename );

    if (write_scalar && write_vector) {
	cerr << "Cannot write both vectors and scalars to vv file."
	     << endl
	     << "Writing vectors only."
	     << endl;
	write_scalar = false;
    }

    if (write_scalar) {
	std::vector<Vector_3>::const_iterator i = displacement.begin();
	for( ; i != displacement.end(); ++i ) {
	    double sq = CGAL::to_double(*i * *i);
	    out << sqrt(sq) << std::endl;
	}
    }

    if (write_vector) {
	std::ostream_iterator<Vector_3> outiter( out, "\n" );
	std::copy( displacement.begin(), displacement.end(), outiter );
    }
}

	
void usage( char* argv0 )
{
    char* s = strrchr( argv0, '/' );
    std::cerr << "usage: " << (s ? s+1 : argv0) 
	      << " [options] source target [map [map...]]"
	      << std::endl
	      << "If mapfile is omitted, an identity map is assumed."
	      << std::endl;
}


int main( int ac, char* av[] )
{
    int write_scalar = 0;
    int write_vector = 1;
    char* vv_file = 0;
    char* vtk_file = 0;

    static ArgvInfo argTable[] = {
	{ "-scalar", ARGV_CONSTANT, (char*)1, (char*)&write_scalar,
          "write displacement length as scalar" },
	{ "-vector", ARGV_CONSTANT, (char*)1, (char*)&write_vector,
          "write displacement vector" },
	{ "-vv", ARGV_STRING, (char*)1, (char*)&vv_file,
          "write vertex value file" },
	{ "-vtk", ARGV_STRING, (char*)1, (char*)&vtk_file,
          "write VTK file" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || (ac < 3) ) {
	usage( av[0] );
	return 1;
    }

    Surface source;
    Surface target;

    try {
	MNI::load_surface_file( source, av[1] );
	MNI::load_surface_file( target, av[2] );
	
	SMap* map_p;
	if ( ac == 3 ) {
	    map_p = new SMap(source,target);
	    cout << "Using vertex-identity map." << endl;
	} else {
	    cout << "1: " << av[3] << endl;
	    map_p = new SMap(source,target,av[3]);
	}

	double N = 1;
	std::vector<Vector_3> displacement( source.size_of_vertices(),
					    CGAL::NULL_VECTOR );
	add_displacements( *map_p, displacement );

	for( int i = 4; i < ac; ++i ) {
	    ++N;
	    cout << N << ": " << av[i] << endl;
	    delete map_p;
	    map_p = new SMap(source,target,av[i]);
	    add_displacements( *map_p, displacement );
	}

	cout << "Read " << N << " surface map file"
	     << (N == 1 ? "." : "s.") << endl;
	scale_displacements( displacement, 1.0/N );

	if (vv_file)
	    write_vv( vv_file, displacement,
		      write_scalar, write_vector );

	if (vtk_file)
	    write_mesh( vtk_file, source, displacement,
			write_scalar, write_vector );

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
