/* Refine control grid of surface map.
 *
 * Inputs: 
 *    target   geometry must be unit sphere
 *    map      input map
 *    control  (new control mesh), must be 1-to-4 subdivision of 
 *             previous control mesh
 *
 */

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polyhedron_3.h>

#include <surflib/Statistic.hpp>
#include <surflib/SurfaceMap_hash.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/Search.hpp>
#include <surflib/rotation.hpp>
#include <surflib/surface_checking.hpp>


typedef CGAL::Simple_cartesian<double>        Surface_kernel;
typedef Surface_kernel::Point_3               Point_3;
typedef Surface_kernel::Vector_3              Vector_3;

typedef CGAL::Polyhedron_3<Surface_kernel>    Surface;

typedef MNI::SurfaceMap_hash<Surface>         SurfaceMap;
typedef MNI::Point_s<SurfaceMap::TargetMesh>  Target_point;


using namespace std;



inline bool is_valid( const Target_point& p )
{
    return p.halfedge() != 0
	&& p.t0() >= 0 && p.t0() <= 1
	&& p.t1() >= 0 && p.t1() <= 1
	&& p.t2() >= 0 && p.t2() <= 1;
}


/* Return point on target surface that is (approximately) midway
 * along the great-circle route from p to q.
 */
Target_point compute_midpoint( const Target_point& p,
			       const Target_point& q )
{
    if (!is_valid(p) || !is_valid(q))
	throw std::domain_error( "compute_midpoint: invalid input" );

    return MNI::find_midpoint_on_sphere( p,q );
}


/* Check that the mesh is a 1-to-4 facet refinement (midpoint refinement)
 * of some initial mesh.  Assume that the initial mesh vertices are 
 * at the beginning of the mesh vertex list.
 *
 * Compute target location for each new control vertex.
 */
void compute_refinement( SurfaceMap& smap )
{
    const Surface& s( smap.control_mesh() );
    CGAL::Inverse_index<Surface::Vertex_const_handle>
	v_index( s.vertices_begin(), s.vertices_end() );

    Surface::size_type size_of_vertices1 = MNI::coarser_size_of_vertices(s);

    // All the existing vertices should have a valid target point.
    //
    Surface::Vertex_const_iterator v = s.vertices_begin();
    for( Surface::size_type i = 0; i < size_of_vertices1; ++i,++v ) {
	if ( !is_valid(smap[v]) )
	    throw std::runtime_error( "expected valid target point" );
    }

    // All new control vertices are assigned new map location.
    // Search neighbours of v for two prexisting vertices (determined
    // by their indices.  Assign the midpoint of the great circle
    // to smap[v].
    //
    CGAL_For_all( v, s.vertices_end() ) {

	if ( is_valid(smap[v]) )
	    throw std::runtime_error( "expected invalid target point" );

	Surface::Vertex_const_handle u1 = 0, u2 = 0;
	Surface::Halfedge_around_vertex_const_circulator 
	    h = v->vertex_begin();

	CGAL_For_all( h, v->vertex_begin() ) {

	    Surface::Vertex_const_handle u = h->opposite()->vertex();
	    if ( v_index[u] < size_of_vertices1 ) {
		if ( u1 == 0 )
		    u1 = u;
		else if ( u2 == 0 )
		    u2 = u;
		else 
		    throw std::runtime_error( "control is not subdivision(2)" );
	    }
	}
	if ( u1 == 0 || u2 == 0 )
	    throw std::runtime_error( "control is not subdivision(3)" );
	
	smap[v] = compute_midpoint( smap[u1], smap[u2] );
    }

}



int main( int ac, char* av[] )
{
    if ( ac != 5 ) {
	std::cerr << "usage: " << av[0] << " map control target output" 
		  << std::endl;
	return 1;
    }

    try {
	Surface control;
	MNI::load_surface_file( control, av[2] );
	cout << "Checking control mesh." << endl;
	MNI::assert_is_sphere_topology( control );
    
	Surface target;
	MNI::load_surface_file( target, av[3] );
	cout << "Checking target mesh." << endl;
	MNI::assert_is_sphere_geometry( target );
	MNI::assert_is_sphere_topology( target );

	SurfaceMap sm( control, target );

	// Set all map locations to invalid surface point
	// -- a sentinel value, for the check in compute_refinement().
	Surface::Vertex_iterator v = control.vertices_begin();
	CGAL_For_all( v, control.vertices_end() )
	    sm[v] = Target_point();

	// Load existing surface map, which covers a subset of
	// the control vertices.

	cout << "Loading existing map." << endl;
	std::ifstream is( av[1] );
	read_surface_map( is, sm, MNI::coarser_size_of_vertices(control) );

	cout << "Refining map." << endl;
	compute_refinement( sm );

	std::ofstream os( av[4] );
	MNI::write_surface_map( os, sm );

    } catch ( const std::exception& e ) {
	cerr << "Error: " << e.what() << endl
	     << "No map output." << endl;
	return 1;
    } catch ( ... ) {
	cerr << "Unknown exception." << endl
	     << "Most likely, this is a bug in the code.  Please report!" 
	     << endl;
	return 2;
    }

    return 0;
}
