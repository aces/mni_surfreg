#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

#include <ParseArgv.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>
#include <surflib/SurfaceMap.hpp>
#include <surflib/LinearInterpolation.hpp>
#include <surflib/NearestNeighbourInterpolation.hpp>


#include "Surface.hpp"
typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;



using namespace std;



struct ScalarData 
    : public std::unary_function<Surface::Vertex_const_handle, double>
{
    double operator() ( const Surface::Vertex_const_handle& v ) const
    { 
        return v->scalar;
    }
};


template <class Interpolation>
struct vertex_value
    : public std::unary_function<const Surface::Vertex&, double>
{
    vertex_value( MNI::SurfaceMap<Surface>& sm ) 
	: f(f_sd), smap(sm) 
    {}

    double operator()( const Surface::Vertex& v )
    {
	Surface::Vertex_const_handle vh(&v);
	return f( smap[vh] );
    }

    ScalarData f_sd;
    Interpolation f;
    MNI::SurfaceMap<Surface>& smap;
};


template <class Interpolation>
void write_vv_resampled_target( const char* filename,
				MNI::SurfaceMap<Surface>& sm )
{
    std::ofstream out( filename );
    vertex_value<Interpolation> vv(sm);

    const Surface& s = sm.control_mesh();
    std::transform( s.vertices_begin(), s.vertices_end(),
		    std::ostream_iterator<double>(out,"\n"),
		    vv );
}


int main( int ac, char* av[] )
{
    int vv_interpolation_linear = 1;

    static ArgvInfo argTable[] = {
        { "-nearest", ARGV_CONSTANT, (char*)0, 
	  (char*)&vv_interpolation_linear,
          "nearest-neighbour interpolation (default: linear)" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || (ac != 6) ) {
	char* s = strrchr( av[0], '/' );
	std::cerr << "usage: " << (s ? s+1 : av[0]) 
		  << " [options] source_geom target_geom target_attr mapfile"
		  << " resampled_attr"
		  << std::endl;
	return 1;
    }

    Surface source;
    Surface target;

    try {
	MNI::load_surface_file( source, av[1] );
	MNI::load_surface_with_scalar( target, av[2], av[3] );

	MNI::SurfaceMap<Surface> smap(source,target,av[4]);

	if (vv_interpolation_linear)
	    write_vv_resampled_target<MNI::LinearInterpolation
		<Surface,ScalarData> >( av[5], smap );
	else
	    write_vv_resampled_target<MNI::NearestNeighbourInterpolation
	        <Surface,ScalarData> >( av[5], smap );

    } catch ( const std::exception& e ) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch ( ... ) {
        cerr << "Unknown exception." << endl
             << "Most likely, this is a bug in the code.  Please report!"
             << endl;
        return 2;
    }

    return 0;
}
