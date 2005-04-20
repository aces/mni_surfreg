#include <CGAL/utility.h>
#include <surflib/load_surface_file.hpp>
#include <surflib/surface_checking.hpp>

#include "Surface.hpp"

typedef Surface_double Surface;



using namespace std;


inline bool approx_equality( double a, double b )
{
    /* The scalars should be equal, in principle.  To accommodate
     * the case that a double is written out as a float, allow
     * a relative error of 1 "units of last place", namely
     * 2^-24 or approximately 1e-7.
     */
    double d = abs(a-b);
    if ( d < 1e-7 * abs(a+b) )    return true;

    cerr << ": " << a << "    " << b << endl
	 << "difference = " << d << endl;
    return false;
}


bool have_equivalent_vertex_list( const Surface& s,
				  const Surface& t )
{
    if ( s.size_of_vertices() != t.size_of_vertices() )
	return false;

    Surface::Vertex_const_iterator vs = s.vertices_begin();
    Surface::Vertex_const_iterator vt = t.vertices_begin();

    CGAL_For_all( vs,s.vertices_end() ) {
	if ( !approx_equality( vs->scalar, vt->scalar ) )
	    return false;
	++vt;
    }
    return true;
}


bool have_equivalent_edge_list( const Surface& s,
				const Surface& t )
{
    if ( s.size_of_halfedges() != t.size_of_halfedges() )
	return false;

    Surface::Halfedge_const_iterator hs = s.halfedges_begin();
    Surface::Halfedge_const_iterator ht = t.halfedges_begin();

    CGAL_For_all( hs,s.halfedges_end() ) {
	if ( !approx_equality( hs->vertex()->scalar,
			       ht->vertex()->scalar )
	     || 
	     !approx_equality( hs->opposite()->vertex()->scalar,
			       ht->opposite()->vertex()->scalar ) )
	    return false;
	++ht;
    }
    return true;
}



int main( int ac, char* av[] ) 
{
    if ( ac != 4 ) {
	cerr << "usage: " << av[0] 
	     << " surf_geom surf_data surf.vtk"
	     << endl;
	return 1;
    }

    try {
	Surface s;
	MNI::load_surface_with_scalar( s, av[1], av[2] );
	MNI::assert_is_sphere_geometry( s );
	MNI::assert_is_sphere_topology( s );

	Surface t;
	MNI::load_vtk_file_with_scalar( t, av[3] );

	if ( have_equivalent_vertex_list(s,t) &&
	     have_equivalent_edge_list(s,t) )
	    return 0;

	return 1;

    } catch ( const std::exception& e ) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    } catch ( ... ) {
        cerr << "Unknown exception." << endl
             << "Most likely, this is a bug in the code.  Please report!"
             << endl;
        return 2;
    }
 
    return 2;
}

