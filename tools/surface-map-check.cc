/* Check that two surface maps are the same (or close).
 * Report differences, e.g. largest difference in t1, t2 values.
 */


#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>

#include <surflib/Statistic.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/SurfaceMap.hpp>

using namespace std;
using MNI::Point_s;


typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;


void compare_locations( MNI::Statistic<double>& discrepancy,
			const Point_s<Surface>& loc1,
			const Point_s<Surface>& loc2 )
{
    double d_sq = CGAL::to_double( CGAL::squared_distance( loc1.point(),
							   loc2.point() ));
    discrepancy.add_sample( sqrt(d_sq) );
}


void compare_maps( const MNI::SurfaceMap<Surface>& sm1,
		   const MNI::SurfaceMap<Surface>& sm2 )
{
    MNI::Statistic<double> discrepancy;

    Surface::Vertex_const_handle v = sm1.control_mesh().vertices_begin();
    for( int i = 0; v != sm1.control_mesh().vertices_end(); ++v,++i ) {
	compare_locations( discrepancy, sm1[v], sm2[v] );
    }

    cout << "Distance discrepancies (3D): " << endl
	 << "  min   : " << discrepancy.min() << endl
	 << "  max   : " << discrepancy.max() << endl
	 << "  mean  : " << discrepancy.mean() << endl
	 << "  stddev: " << discrepancy.std_dev() << endl;
}


int main( int ac, char* av[] )
{
    if ( ac != 5 ) {
	cerr << "usage: " << av[0] << " source target first.sm second.sm"
	     << endl;
	return 1;
    }

    Surface source;
    Surface target;

    int arg_index = 1;
    try {
	MNI::load_surface_file( source, av[arg_index++] );
	MNI::load_surface_file( target, av[arg_index++] );
    } catch ( const std::exception& e ) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch ( ... ) {
        cerr << "Unknown exception." << endl
             << "Most likely, this is a bug in the code.  Please report!"
             << endl;
        return 2;
    }

    MNI::SurfaceMap<Surface> sm1(source,target,av[arg_index++]);
    MNI::SurfaceMap<Surface> sm2(source,target,av[arg_index++]);

    compare_maps(sm1,sm2);

    return 0;
}
