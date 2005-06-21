/* Compute the alpha shape, write classification of vertices.
 * See README.alpha-shape for some definitions.
 *
 * The output consists of three values for each vertex.
 * 1. an integer, representing the classification as one of:
 *    EXTERIOR, SINGULAR, REGULAR, or INTERIOR
 * 2. lower threshold; for alpha values below this, the vertex is SINGULAR
 * 3. upper threshold; for alpha values above this, the vertex is INTERIOR
 *
 * For alpha values on the interval [lower,upper], the vertex is REGULAR,
 * i.e., it forms part of the boundary.
 */

#include <iostream>
#include <iterator>
#include <map>

#include <ParseArgv.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Alpha_shape_euclidean_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>


template <class Refs, class Traits>
struct My_vertex 
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
					  typename Traits::Point_3>
{
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> Vertex_base;

    My_vertex( const Point& p ) : Vertex_base(p) {};
    My_vertex() {};

    double low,high;
    int label;
};

struct My_items : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Vertex_wrapper {
	typedef typename Traits::Point_3 Point;
	typedef My_vertex<Refs,Traits> Vertex;
    };
};


typedef CGAL::Cartesian<CGAL::MP_Float>  K;
typedef CGAL::Polyhedron_3<K, My_items>  Surface;

typedef CGAL::Alpha_shape_euclidean_traits_3<K> Gt;
typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;
typedef CGAL::Triangulation_cell_base_3<Gt> Df;
typedef CGAL::Alpha_shape_cell_base_3<Gt, Df>  Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

typedef Surface::Vertex_handle       Vertex_handle;
typedef Surface::Vertex_iterator     Vertex_iterator;

typedef K::Point_3      Point_3;
typedef CGAL::Unique_hash_map<const K::RT*,Vertex_handle> Point_map;


using namespace std;





void classify_vertices( Surface& s, 
			bool classify, Alpha_shape_3::Coord_type alpha )
{
    for( Vertex_iterator v = s.vertices_begin(); v != s.vertices_end(); ++v )
	v->label = -1;

    /* The alpha shape is initialized with a copy of the points of s.
     * In order to make the link from the shape vertices to the
     * vertices of s, we use a map from Point_3 to Vertex, set up
     * using the points of surface s.  For this to work, Point_3 must
     * use Handle/Rep implementation; i.e. Cartesian rather than
     * Simple_cartesian.
     */
    Point_map point_to_vertex;
    Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
        point_to_vertex[& v->point().x()] = v;
    }
    
    cout << "Creating alpha shape ..." << flush;
    Alpha_shape_3 shape( s.points_begin(), s.points_end() );
    cout << "done." << endl;

    if (classify && alpha < 0) {
	cout << "Optimal alpha value for 1 component is " << flush;
	alpha = *shape.find_optimal_alpha(1);
	cout << alpha << endl;
    }

    if ( alpha > 0 ) 
	cout << "Number of solid components: " 
	     << shape.number_of_solid_components(alpha) << endl;


    /* Set classification value in vertices of s.
     */
    Triangulation_3::Finite_vertices_iterator 
	vi = shape.finite_vertices_begin();
    CGAL_For_all( vi, shape.finite_vertices_end() ) {
	CGAL_assertion( point_to_vertex.is_defined( & vi->point().x() ) );
	Surface::Vertex_handle v = point_to_vertex[ & vi->point().x() ];
	CGAL_assertion( v != point_to_vertex.default_value() );

	Alpha_shape_3::Alpha_status* as = vi->get_alpha_status();

	v->low  = CGAL::to_double( as->alpha_mid() );

	if ( as->is_on_chull() )
	    v->high = -1;
	else
	    v->high = CGAL::to_double( as->alpha_max() );

	if (classify)
	    v->label = shape.classify(vi,alpha);
	else
	    v->label = 0;
    }

    for( Vertex_iterator v = s.vertices_begin(); v != s.vertices_end(); ++v )
	CGAL_postcondition( v->label != -1 );
}



void process( char* input_file, 
	      char* vtk_file,
	      char* vv_file,
	      double alpha,
	      int optimal_alpha )
{
    bool classify = optimal_alpha || (alpha > 0);

    Surface s;

    MNI::load_surface_file( s, input_file );
    classify_vertices( s, classify, alpha );

    if ( vtk_file ) {
	// Write in VTK format.
	// Note: assuming that write_vtk() outputs vertices
	// in order [vertices_begin,vertices_end).
	std::ofstream out( vtk_file );
	MNI::write_vtk( out, s, 
			"Vertices classified based on alpha shape." );
	out << "POINT_DATA " << s.size_of_vertices() << std::endl;

	if ( classify ) {
	    out << "SCALARS class int" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	    for( Vertex_iterator v = s.vertices_begin();
		 v != s.vertices_end(); ++v )
		out << v->label << std::endl;
	}

	out << "SCALARS low float" << endl
	    << "LOOKUP_TABLE default" << endl;
	for( Vertex_iterator v = s.vertices_begin();
	     v != s.vertices_end(); ++v )
	    out << v->low << std::endl;

	out << "SCALARS high float" << endl
	    << "LOOKUP_TABLE default" << endl;
	for( Vertex_iterator v = s.vertices_begin();
	     v != s.vertices_end(); ++v )
	    out << v->high << std::endl;
    }

    if ( vv_file ) {
	// Dump (label,low,high) values to file in vertex order.
	std::ofstream out( vv_file );
	for( Vertex_iterator v = s.vertices_begin();
	     v != s.vertices_end(); ++v )
	    out << v->label 
		<< " " << v->low 
		<< " " << v->high
		<< std::endl;
    }
}


void usage( char* argv0 )
{
    char* s = strrchr( argv0, '/' );
    std::cerr << "usage: " << (s ? s+1 : argv0) << " [options] input"
	      << std::endl;
}


int main( int ac, char* av[] )
{
    char* vtk_file = 0;
    char* vv_file = 0;  /* vertex values */

    double alpha = -1;
    int optimal_alpha = 0;

    ArgvInfo argTable[] = {
	{ "-alpha", ARGV_FLOAT, (char*)1, (char*)&alpha,
	  "specify alpha value used to classify each vertex" },
	{ "-opt_alpha", ARGV_CONSTANT, (char*)1, (char*)&optimal_alpha,
	  "compute optimal alpha and classify vertices using it" },
	{ "-vtk", ARGV_STRING, (char*)1, (char*)&vtk_file,
	  "write output in VTK format" },
	{ "-vv", ARGV_STRING, (char*)1, (char*)&vv_file,
	  "write depth values in file, one per line" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };


    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 2 ) {
        usage( av[0] );
        return 1;
    }
    if ( optimal_alpha && (alpha > 0) ) {
	cerr << "Ignoring -alpha when -opt_alpha specified." << endl;
	alpha = -1;
    }
    if ( !optimal_alpha && vtk_file == 0 && vv_file == 0 ) {
	cerr << "Specify one of -opt_alpha, -vtk, or -vv for output.\n";
	return 1;
    }

    try {
	process( av[1], vtk_file, vv_file, alpha, optimal_alpha );
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
