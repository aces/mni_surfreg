/* Generate initial surface mapping
 * (vertex 0 --> vertex 0, vertex 1 --> vertex 1)
 * from control and target surface files.
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>
                                                                               
#include <surflib/SurfaceMap_hash.hpp>
#include <surflib/load_surface_file.hpp>
                                                                              
typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;
typedef Surface_kernel::Point_3             Point_3;

typedef MNI::SurfaceMap_hash<Surface>       SurfaceMap;


void process( char* control_fn,
	      char* target_fn,
	      char* output_fn )
{
    Surface control;
    MNI::load_surface_file( control, control_fn );

    Surface target;
    MNI::load_surface_file( target, target_fn );

    SurfaceMap sm( control, target );
    std::ofstream os( output_fn );
    MNI::write_surface_map( os, sm );
}


int main( int ac, char* av[] )
{
    if ( ac != 4 ) {
	std::cerr << "usage: " << av[0] << " control target output" 
		  << "\n";
	return 1;
    }

    try {
	process( av[1], av[2], av[3] );
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
