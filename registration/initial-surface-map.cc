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



int main( int ac, char* av[] )
{
    if ( ac != 4 ) {
	std::cerr << "usage: " << av[0] << " control target output" 
		  << std::endl;
	return 1;
    }

    Surface control;
    MNI::load_surface_file( control, av[1] );

    Surface target;
    MNI::load_surface_file( target, av[2] );

    SurfaceMap sm( control, target );
    std::ofstream os( av[3] );
    MNI::write_surface_map( os, sm );

    return 0;
}
