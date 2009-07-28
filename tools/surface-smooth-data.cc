/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* Average scalar value at each vertex with neighbours.
 */

#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

#include <surflib/VTKFile.hpp>
#include <surflib/vtk_ostream.hpp>

#include "Surface.hpp"
typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;


typedef Surface::Vertex_iterator Vertex_iterator;

using namespace std;



/* Replace the scalar at each vertex v by the weighted sum
 * of scalars in the 1-neighbourhood (i.e. v plus its 1-ring).
 *
 * The scalar at v is weighted by 1/s, and the others are weighted by
 * a/s, where a is specified by the user and s is the sum of all
 * coefficients, so that the weights add up to 1.
 */
void smooth_scalars( Surface& s, double a ) 
{
    std::vector<double> new_scalar( s.size_of_vertices() );

    Vertex_iterator v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++v,++i ) {
	int degree = 0;
	new_scalar[i] = 0;
	Surface::Halfedge_around_vertex_circulator h = v->vertex_begin();
	CGAL_For_all( h, v->vertex_begin() ) {
	    ++degree;
	    new_scalar[i] += h->opposite()->vertex()->scalar;
	}

	new_scalar[i] *= a;
	new_scalar[i] += v->scalar;
	new_scalar[i] /= 1.0 + a*static_cast<double>(degree);
	//cout << i << ": " << new_scalar[i] << endl;
    }

    // Copy scalars back to surface structure
    v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++v,++i ) {
	v->scalar = new_scalar[i];
    }
}


int main( int ac, char* av[] )
{
    if ( ac != 4 && ac != 5 ) {
	cerr << "usage: " << av[0] << " input.vtk output.vtk a [r]" << endl
	     << "  a is the coefficient for neighbours" << endl
	     << "  r is the number of repeats (default 1)" << endl;
	return 1;
    }

    Surface s;
    try {
	MNI::VTKFile vtk( av[1] );

	cout << "File " << av[1] << " header: " 
	     << vtk.get_header() << endl;
	if ( vtk.is_binary() )
	    cout << "File is BINARY" << endl;
	else
	    cout << "File is ASCII" << endl;

	vtk.read_geometry(s);
	MNI::SetScalarValues<Surface> ssv(s);
	vtk.read_scalar_point_attribute(ssv);


	double a = atof( av[3] );
	int r = 1;
	if ( ac == 5 )
	    r = atoi( av[4] );
	if ( r <= 0 ) {
	    cerr << "Number of repeats must be positive.  r = " << r << endl;
	    return 1;
	}

	cout << r << " repetitions of smoothing with a = " << a << endl;

	while( --r >= 0 ) 
	    smooth_scalars(s,a);

	std::ofstream out( av[2] );
	MNI::write_vtk( out, s, std::string("smooth-scalars ") 
			+ av[1] + " " + av[2] + " " + av[3] );
	out << "POINT_DATA " << s.size_of_vertices() << std::endl
	    << "SCALARS data double" << std::endl
	    << "LOOKUP_TABLE default" << std::endl;
	for( Vertex_iterator v = s.vertices_begin();
	     v != s.vertices_end(); ++v )
	    out << v->scalar << std::endl;

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
