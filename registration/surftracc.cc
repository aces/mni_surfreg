/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* Optimize a surface mapping.
 *
 * Inputs:
 *    source   mesh with scalar value
 *    target   mesh with scalar value, geometry must be unit sphere
 *    control  sub-mesh of source
 *    map      from control to target
 *    output   filename of optimized surface map
 */

#include <iostream>
#include <fstream>
#include <exception>
#include <cmath>

#include <ParseArgv.h>

#include "surftracc.hpp"
#include "objective.hpp"

#include <surflib/Diagnostic_ostream.hpp>
#include <surflib/Statistic.hpp>
#include <surflib/load_surface_file.hpp>
#include <surflib/surface_checking.hpp>


// Send diagnostics here
MNI::Diagnostic_ostream diag;



using namespace std;


template <class ComputeCentre>
class NC : public NeighbourhoodCentre
{
public:
    NC( SurfaceMap& sm, ComputeCentre& cc )
	: smap(sm),f(cc) {}

    virtual const Vector_3 operator()( ControlMesh::Vertex_const_handle v )
    {
	return f( smap, v );
    }


protected:
    SurfaceMap& smap;
    ComputeCentre f;
};


/* Returns centroid of 1-ring neighbourhood, equally weighted.
 */
struct ComputeCentre_v1
{
    const Vector_3 operator()( SurfaceMap& smap,
			       ControlMesh::Vertex_const_handle v )
    {
	Vector_3 t_vector( CGAL::NULL_VECTOR );
	double wu_sum = 0;

	ControlMesh::Halfedge_around_vertex_const_circulator h = v->vertex_begin();
	CGAL_For_all( h,v->vertex_begin() ) {
	    t_vector = t_vector +
		(smap.get_target(h->opposite()->vertex()).point() 
		 - CGAL::ORIGIN);
	    wu_sum += 1.0;
	}

	return ( (smap.get_target(v->vertex_begin()->vertex()).point() - CGAL::ORIGIN ) + 
                 2.0 * t_vector / wu_sum ) / 3.0;
    }
};


/* Returns centroid of 1-ring neighbourhood, weighted by edge length.
 */
struct ComputeCentre_v3
{
    const Vector_3 operator()( SurfaceMap& smap,
			       ControlMesh::Vertex_const_handle v )
    {
	Vector_3 t_vector( CGAL::NULL_VECTOR );
	double wu_sum = 0;

	ControlMesh::Halfedge_around_vertex_const_circulator h = v->vertex_begin();
	CGAL_For_all( h,v->vertex_begin() ) {
	    Vector_3 d( CGAL::NULL_VECTOR );
            d = smap.get_target(h->opposite()->vertex()).point() - 
                smap.get_target(v->vertex_begin()->vertex()).point();
	    double mag = std::sqrt(CGAL::to_double(d*d));

	    t_vector = t_vector +
		mag * (smap.get_target(h->opposite()->vertex()).point() - CGAL::ORIGIN);

	    wu_sum += mag;
	}

	return ( (smap.get_target(v->vertex_begin()->vertex()).point() - CGAL::ORIGIN ) + 
                 2.0 * t_vector / wu_sum ) / 3.0;
    }
};

/* Returns centroid of 1-ring neighbourhood, weighted
 * by relative area of incident facets.
  * NOTE: This is very very bad. Don't use it. CL.
 */
struct ComputeCentre_v2 
{
    const Vector_3 operator()( SurfaceMap& smap,
			       ControlMesh::Vertex_const_handle v )
    {
	Vector_3 t_vector( CGAL::NULL_VECTOR );
	double wu_sum = 0;
	
	ControlMesh::Halfedge_const_handle h = v->halfedge();
	do {
	    CGAL_assertion( v == h->vertex() );

	    Vector_3 u = smap.get_target(h->opposite()->vertex()).point()
		- CGAL::ORIGIN;

	    double wu = std::max( 1e-6, normalized_area(smap, h) )
		+ std::max( 1e-6, normalized_area(smap,h->opposite()) );

	    wu_sum += wu;
	    t_vector = t_vector + u * wu;
	    h = h->next()->opposite();
	} while ( h != v->halfedge() );

	CGAL_assertion( wu_sum > 0 );
	return t_vector / wu_sum;
    }
};



struct project_target : public unary_function<const ControlMesh::Vertex&, 
			                      SurfaceMap::Target_point>
{
    const SurfaceMap& smap;

    project_target( const SurfaceMap& sm ) : smap(sm) {}

    const SurfaceMap::Target_point operator()( const ControlMesh::Vertex& v )
    { return smap.get_target( ControlMesh::Vertex_const_handle(&v) ); }
};


/*! Report displacement magnitudes between the point locations
 * of \a t_list and the target locations of \a smap.  Output
 * sent to diag.level(1).
 */
void report_displacement_mag( char* prefix,
			      SurfaceMap& smap,
			      std::vector<SurfaceMap::Target_point>& t_list )
{
    const SurfaceMap::ControlMesh& control( smap.control_mesh() );

    CGAL_precondition( control.size_of_vertices() == t_list.size() );

    MNI::Statistic<double> displacement;
    SurfaceMap::ControlMesh::Vertex_const_iterator 
	v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {
	Vector_3 d = smap[v].point() - t_list[i].point();
	double mag = std::sqrt(CGAL::to_double(d*d));
// cout << " " << i << " " << smap[v].point().x() << " " << smap[v].point().y()
//      << " " << smap[v].point().z() << " to " << t_list[i].point().x() << " " 
//      << t_list[i].point().y() << " " << t_list[i].point().z() << " mag " << mag << endl;
	displacement.add_sample( mag );
    }

    diag.level(1) << prefix;
    displacement.print( diag.level(1) );
}


/* Check whether the current optimization step produced a
 * map with small triangles on target, indicating a danger
 * of non-injectivity.
 *
 * If the step is good, return true.
 *
 * Otherwise, re-set map targets, reduce search radii,
 * and return false.
 */
bool accept_step( int debug,
		  SurfaceMap& smap,
		  const std::vector<SurfaceMap::Target_point>& map_start,
		  std::vector<double>& search_radius,
		  double sr_factor )
{
    ControlMesh& control( smap.control_mesh() );

    std::vector<bool> mark( control.size_of_vertices(), false );
    int count = 0;

    ControlMesh::Vertex_iterator v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {
	if ( is_incident_to_small_facet(smap,v) ) {
	    mark[i] = true;
	    search_radius[i] *= sr_factor;
	    ++count;
	}
    }
    if ( count == 0 ) return true;

    diag.level(1) << "accept_step: small facet count = " << count
		  << endl;

    diag.level(1) << "continuing anyway" << endl;
    return true;
}


int main( int ac, char* av[] )
{
    int debug = 0;

    // Convergence parameters.
    int outer_iter_max = 20;
    int inner_iter_max = 100;
    double abs_tol = 1e-3;

    // Neighbourhood Search parameters.
    double search_radius = 0.5;
    double initial_radius = 0.9;
    int radius_is_absolute = 0;

    // Objective function parameters.
    double neighbourhood_radius = 2.8;
    double penalty_ratio = 0.05;

    // Smoothing parameters.
    double smoothing_neighbour_weight = 1;
    int nc_method = 1;

    // Injectivity ...
    double sr_factor = 1;

    static ArgvInfo argTable[] = {
	{ "-debug", ARGV_INT, 0, (char*)&debug,
	  "set debug level" },

	// --- Convergence parameters --- //

	{ "-outer_iter_max", ARGV_INT, 0, (char*)&outer_iter_max,
	  "iteration limit of outer loop"},
	{ "-inner_iter_max", ARGV_INT, 0, (char*)&inner_iter_max,
	  "iteration limit of inner loop"},
	{ "-abs_tolerance", ARGV_FLOAT, 0, (char*)&abs_tol,
	  "absolute tolerance sought in optimization" },

	// --- Neighbourhood Search parameters --- //

	{ "-search_radius", ARGV_FLOAT, 0, (char*)&search_radius,
          "maximum vertex displacement at each iteration (max. 1)" },
	{ "-initial_radius", ARGV_FLOAT, 0, (char*)&initial_radius,
	  "scale of initial optimizer steps" },
	{ "-absolute_radius", ARGV_CONSTANT, 
	  (char*)1, (char*)&radius_is_absolute,
	  "treat radii as absolute values (do NOT use this option)"},

	// --- Objective Function parameters --- //

	{ "-neighbourhood_radius", ARGV_FLOAT, 0, (char*)&neighbourhood_radius,
          "radius of texture neighbourhood (< 1 if absolute, > 1 if multiplicative factor)" },
	{ "-penalty_ratio", ARGV_FLOAT, 0, (char*)&penalty_ratio,
          "coefficient of distance penalty term" },

	// --- Smoothing parameters --- //

	{ "-smoothing_weight", 
	  ARGV_FLOAT, 0, (char*)&smoothing_neighbour_weight,
          "neighbour weight in smoothing step" },
	{ "-nc_method", ARGV_INT, 0, (char*)&nc_method,
	  "neighbourhood centre computation method (always use 1)" },

	// --- Injectivity --- //

	{ "-sr_factor",
	  ARGV_FLOAT, 0, (char*)&sr_factor,
          "set < 1 to reduce search_radius of bad vertices" },

	{ NULL, ARGV_END, NULL, NULL, NULL }
    };


    // Grab the command line first, since ParseArgv will remove
    // arguments from the argv array.
    std::string command_line(av[0]);
    for( int i = 1; i < ac; ++i ) {
	command_line += " ";
	command_line += av[i];
    }

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 8 ) {
	cerr << "usage: " << av[0] 
	     << " source_geom source_data"
	     << " target_geom target_data"
	     << " control in.sm out.sm"
	     << endl;
        cerr << endl << "Copyright Alan C. Evans" << endl
                     << "Professor of Neurology" << endl
                     << "McGill University" << endl;
	return 1;
    }

#define D(var) diag.level(1) << (#var) << " = " << var << endl;

    diag.maximum_level(debug);
    diag.level(1) << "$Id$" << endl;
    diag.level(1) << "Command line:" << endl
		  << "  " << command_line << endl << endl;

    D(outer_iter_max);
    D(inner_iter_max);
    D(abs_tol);

    D(search_radius);
    D(initial_radius);
    D(radius_is_absolute);

    D(neighbourhood_radius);
    D(penalty_ratio);

    D(smoothing_neighbour_weight);
    D(nc_method);

    D(sr_factor);


    Surface source;
    Surface target;
    SurfaceMap::ControlMesh control;

    try {
	MNI::load_surface_with_scalar( source, av[1], av[2] );
	MNI::assert_is_sphere_geometry( source );
	MNI::assert_is_sphere_topology( source );

	MNI::load_surface_with_scalar( target, av[3], av[4] );
	MNI::assert_is_sphere_geometry( target );
	MNI::assert_is_sphere_topology( target );

	MNI::load_surface_file( control, av[5] );
	MNI::assert_is_sphere_topology( control );

	// Read map from file.
	SurfaceMap smap( source, target, control, av[6] );

	if ( count_small_targets(smap) > 0 ) {
	    cerr << "Input map is non-injective!" << endl;
//	    return 1;
	}

	ComputeCentre_v1 cc_v1;
	ComputeCentre_v2 cc_v2;
	ComputeCentre_v3 cc_v3;
	NeighbourhoodCentre* nc_ptr = 0;

	switch(nc_method) {
	case 1: 
	    nc_ptr = new NC<ComputeCentre_v1>(smap,cc_v1);
	    break;
	case 2: 
	    nc_ptr = new NC<ComputeCentre_v2>(smap,cc_v2);
	    break;
	case 3: 
	    nc_ptr = new NC<ComputeCentre_v3>(smap,cc_v3);
	    break;
	default:
	    cerr << "Illegal -nc_method value: " << nc_method << endl;
	    return 1;
	}


	// Track the number of nonlinear simplex iterations.
	MNI::Statistic<double> simplex_iterations;

        // If neighbourhood_radius < 1, assume it's an absolute length.
        // If it's greater than 1, assume it's a scaling factor so determine
        // a characteristic length from the average edge length on the 
        // control sphere. CL.

	if ( neighbourhood_radius > 1 ) {

            //  Compute an absolute ngh_radius based on the total area of
            //  the unit control sphere. Assume isolateral triangles. Determine
            //  the characteristic edge length. Multiply it by the scaling
            //  factor given by neighbourhood_radius. Project it on the plane
            //  parallel to the tangent plane to the control point. CL.

            int n_cells = smap.control_mesh().size_of_facets();
            double unit_sphere_area = 4.0 * M_PI;  // rad=1.0
            double avg_control_edge_len = std::sqrt( ( unit_sphere_area / n_cells ) *
                                                     ( 4.0 / std::sqrt( 3.0 ) ) );
            neighbourhood_radius *= avg_control_edge_len;
            if( neighbourhood_radius < M_PI / 2.0 ) {
                // the arclength is r*angle, r=1, so neighbourhood_radius is
                // an angle in radians.
                neighbourhood_radius = std::sin( neighbourhood_radius );
            } else {
                neighbourhood_radius = 0.99;
            }

	    diag.level(1) << "Absolute neighbourhood radius = " 
			  << neighbourhood_radius << endl;
	}

	// ---------------------------
	// OUTER LOOP
	// ---------------------------

	vector<double> current_search_radius( control.size_of_vertices(),
					      search_radius );

	// Store starting value of map at each iteration.
	std::vector<SurfaceMap::Target_point> 
	    map_start( control.size_of_vertices() );

	int outer_iter = 1;
	int repeat_count = 0;
	while( outer_iter <= outer_iter_max ) {
	    std::transform( control.vertices_begin(), control.vertices_end(),
			    map_start.begin(), 
			    project_target(smap) );
		       
	    ObjectiveFunction f( smap, 
				 neighbourhood_radius,
				 penalty_ratio );
	    diag.level(1) << "Optimize map iteration #" << outer_iter << endl;
	    optimize_map_gsl( smap, f,
			      !radius_is_absolute,
			      current_search_radius,
			      initial_radius,
			      abs_tol,
			      inner_iter_max,
			      simplex_iterations );

	    diag.level(1) << "Smooth map iteration #" << outer_iter << endl;
	    smooth_map( smap, smoothing_neighbour_weight, *nc_ptr );

	    if ( accept_step( debug, smap, map_start, 
			      current_search_radius, sr_factor ) ) {
		++outer_iter;
		repeat_count = 0;
		report_displacement_mag( (char*)"TOTAL ", smap, map_start );
	    } else {
		if ( repeat_count == 20 ) {
		    cerr << "Repeated current iteration 20 times."
			 << "  Giving up."
			 << endl;
		    break;
		}
		++repeat_count;
		diag.level(1) << "Repeating iteration " << outer_iter
			      << endl;
	    }
	}

	diag.level(1) << "Simplex iteration counts: " << simplex_iterations
		      << endl;

	std::ofstream map_out( av[7] );
	MNI::write_surface_map( map_out, smap );

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
