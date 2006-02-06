/* Generate initial surface mapping
 * (vertex 0 --> vertex 0, vertex 1 --> vertex 1)
 * from control and target surface files.
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>

#include <ParseArgv.h>

#include <surflib/SurfaceMap_hash.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/Search.hpp>

                                                                              
typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;
typedef Surface_kernel::Point_3             Point_3;

typedef MNI::SurfaceMap_hash<Surface>       SurfaceMap;


// Booleans that indicate whether to flip across each axis.
static int flip_x = 0;
static int flip_y = 0;
static int flip_z = 0;


Point_3 flipPoint( const Point_3& p )
{
    return Point_3( flip_x ? - p.x() : p.x(), 
		    flip_y ? - p.y() : p.y(), 
		    flip_z ? - p.z() : p.z() );
}



void flip( SurfaceMap& smap )
{
    const SurfaceMap::SourceMesh& source( smap.source_mesh() );
    const SurfaceMap::TargetMesh& target( smap.target_mesh() );

    SurfaceMap::TargetMesh::Vertex_const_iterator 
	vTarget = target.vertices_begin();

    SurfaceMap::SourceMesh::Vertex_const_iterator v = source.vertices_begin();
    for( ; v != source.vertices_end(); ++v ) 
    {
	Point_3 initial = v->point();
	Point_3 final = flipPoint( initial );
	smap[v] = MNI::find_nearest_on_convex_surface<Surface>( vTarget, final );
    }
}

void process( char* control_fn,
	      char* target_fn,
	      char* output_fn )
{
    Surface control;
    MNI::load_surface_file( control, control_fn );

    Surface target;
    MNI::load_surface_file( target, target_fn );

    SurfaceMap sm( control, target );

    if ( flip_x || flip_y || flip_z )
	flip( sm );

    std::ofstream os( output_fn );
    MNI::write_surface_map( os, sm );
}


int main( int ac, char* av[] )
{
    static ArgvInfo argTable[] = {
	{ "-flip-x", ARGV_CONSTANT, (char*)1, (char*)&flip_x,
          "mapping will reflect across plane X=0" },
	{ "-flip-y", ARGV_CONSTANT, (char*)1, (char*)&flip_y,
          "mapping will reflect across plane Y=0" },
	{ "-flip-z", ARGV_CONSTANT, (char*)1, (char*)&flip_z,
          "mapping will reflect across plane Z=0" },

	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 4 ) {
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
