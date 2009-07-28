/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* Statistics on area ratios.
 *
 * Inputs:
 *    source   mesh with scalar value
 *    target   mesh with scalar value, geometry must be unit sphere
 *    control  sub-mesh of source
 *    map      from control to target
 */

#include <iostream>
#include <exception>
#include <ParseArgv.h>

/* File utility.h defines template class Quadruple<> which has
 * "template < class U, class V, class W, class X>" methods.
 * Unfortunately, the MNI volume_io library has a "#define X ..."
 * which results in a compile error if utility.h is later included.
 * To avoid this problem, we include it early.
 */
#include <CGAL/utility.h>

#include <surflib/Statistic.hpp>

#include <surflib/load_surface_file.hpp>
#include <surflib/surface_checking.hpp>

// FIXME
#include "registration/surftracc.hpp"



using namespace std;


void process_triangle_area_ratios( const SurfaceMap& smap,
				   char* out_fn,
				   double small_area_ratio = 0 )
{
    ofstream* out = 0;
    if ( out_fn )
	out = new ofstream(out_fn);

    ControlMesh::Facet_const_iterator 
	f = smap.control_mesh().facets_begin();

    MNI::Statistic<double> ratios;
    MNI::Statistic<double> small_ratios;

    CGAL_For_all( f, smap.control_mesh().facets_end() ) {
	double r = normalized_area(smap,f->halfedge());

	ratios.add_sample(r);
	if ( r < small_area_ratio )
	    small_ratios.add_sample(r);

	if (out_fn)
	    (*out) << r << endl;
    }

    if (out_fn)
	delete out;

    cout << "Area ratios:  ";
    ratios.print( std::cout );

    cout << "Ratios < " << small_area_ratio << ":  ";
    small_ratios.print( cout );
}


void process_flipped_triangles( const SurfaceMap& smap,
				char* out_fn )
{
    ofstream* out = 0;
    if ( out_fn )
	out = new ofstream(out_fn);

    ControlMesh::Facet_const_iterator 
	f = smap.control_mesh().facets_begin();

    MNI::Statistic<double> st_area;

    CGAL_For_all( f, smap.control_mesh().facets_end() ) {
	ControlMesh::Halfedge_const_handle h = f->halfedge();

	Point_3 A = smap.get_target(h->vertex()).point();
	Point_3 B = smap.get_target(h->next()->vertex()).point();
	Point_3 C = smap.get_target(h->next()->next()->vertex()).point();

	double area = signed_area(A,B,C);
	if ( area < 0 )
	    st_area.add_sample(-area);

	if (out_fn)
	    (*out) << area << endl;
    }

    cout << "Flipped triangle area statistics:  ";
    st_area.print( std::cout );
    cout << "Total area = " << st_area.sum() << endl;
}


int main( int ac, char* av[] )
{
    int proc_ratios = 0;
    int proc_flipped = 0;

    char* ratio_file = 0;
    char* area_file = 0;

    static ArgvInfo argTable[] = {
	{ "-proc_ratios", ARGV_CONSTANT, (char*)1, (char*)&proc_ratios,
	  "process target / source control triangle ratios." },
	{ "-proc_flipped", ARGV_CONSTANT, (char*)1, (char*)&proc_flipped,
	  "process area of flipped control triangles." },

	{ "-ratio_file", ARGV_STRING, (char*)1, (char*)&ratio_file,
          "write triangle area ratios in file, one per line" },
	{ "-area_file", ARGV_STRING, (char*)1, (char*)&area_file,
          "write signed triangle areas in file, one per line" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 5 ) {
	cerr << "usage: " << av[0] << " source target control in.sm"
	     << endl;
	return 1;
    }

    try {
	Surface source;
	MNI::load_surface_file( source, av[1] );
	MNI::assert_is_sphere_geometry( source );
	MNI::assert_is_sphere_topology( source );

	Surface target;
	MNI::load_surface_file( target, av[2] );
	MNI::assert_is_sphere_geometry( target );
	MNI::assert_is_sphere_topology( target );

	SurfaceMap::ControlMesh control;
	MNI::load_surface_file( control, av[3] );
	MNI::assert_is_sphere_topology( control );

	// Read map from file.
	SurfaceMap smap( source, target, control, av[4] );

	if (proc_ratios)
	    process_triangle_area_ratios( smap, ratio_file );
	if (proc_flipped)
	    process_flipped_triangles( smap, area_file );

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
