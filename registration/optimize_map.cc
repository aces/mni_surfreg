#include <iostream>

#include <surflib/Diagnostic_ostream.hpp>
#include <surflib/Statistic.hpp>
#include "surftracc.hpp"
#include "data_term.hpp"
#include "objective.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>


// Send diagnostics here
extern MNI::Diagnostic_ostream diag;


using namespace std;
using MNI::Point_s;



double
radius_to_nearest_neighbour( SurfaceMap& smap,
			     SurfaceMap::ControlMesh::Vertex_const_handle v )
{
    double max_cos_angle = 0;

    Vector_3 vec_v = smap[v].point() - CGAL::ORIGIN;
    ControlMesh::Halfedge_const_handle h = v->halfedge();
    do {
	Vector_3 vec_u = smap[h->opposite()->vertex()].point() - CGAL::ORIGIN;
	double cos_angle = CGAL::to_double(vec_v*vec_u);
	if ( cos_angle > max_cos_angle )
	    max_cos_angle = cos_angle;
	h = h->next()->opposite();
    } while ( h != v->halfedge() );

    double radius_sq = 1 - max_cos_angle*max_cos_angle;
    CGAL_assertion( radius_sq > 0 );
    return std::sqrt( radius_sq );
}


double my_f( const gsl_vector* x, void* params)
{
    ObjectiveFunction& f( *static_cast<ObjectiveFunction*>(params) );
    return f( gsl_vector_get( x, 0 ),
	      gsl_vector_get( x, 1 ) );
}



// TODO investigate GSL error handling

/*! \brief Optimize surface map.
 *
 * If \a radius_is_relative is true, then \a penalty_radius
 * is treated as a fraction of the disance to the nearest neighbour,
 * and \a initial_radius is treated as a fraction of the computed
 * penalty radius.
 */
void optimize_map_gsl( SurfaceMap& smap,
		       ObjectiveFunction& f,
		       bool radius_is_relative,
		       const std::vector<double>& penalty_radius,
		       double initial_radius,
		       double abs_tol,
		       int max_num_iter,
		       MNI::Statistic<double>& simplex_iterations )
{
    SurfaceMap::ControlMesh& control( smap.control_mesh() );

    // Output storage.
    std::vector<Point_s<Surface> > final( control.size_of_vertices() );

    const int dimension = 2;

    // T is the algorithm.
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex;

    // s is the algorithm state.
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc (T, dimension);

    /* Function to minimize */
    gsl_multimin_function minex_func;
    minex_func.f = &my_f;
    minex_func.n = dimension;
    minex_func.params = (void *)&f;

    /* Starting point */
    gsl_vector* x = gsl_vector_alloc(dimension);

    /* Initial vertex size vector */
    gsl_vector* ss = gsl_vector_alloc(dimension);

    // Track the search radii.
    MNI::Statistic<double> search_radius_stat;


    SurfaceMap::ControlMesh::Vertex_iterator v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {
	f.set_vertex(v);

	double penalty_radius_v = penalty_radius[i];
	double initial_radius_v = initial_radius;

	if ( radius_is_relative ) {
	    double r_near = radius_to_nearest_neighbour( smap, v );
	    penalty_radius_v = std::min( penalty_radius_v * r_near, 1.0 );
	    initial_radius_v *= penalty_radius_v;
	}

	search_radius_stat.add_sample( penalty_radius_v );
	diag.level(3) << "Starting value: " << f(0,0) << std::endl
		      << "Search radius set to " << penalty_radius_v
		      << " (simplex scale = " << initial_radius_v << ")"
		      << std::endl;

	f.set_penalty_radius( penalty_radius_v );

	/* Set initial iterate x, and step sizes ss */
	gsl_vector_set_all(x, 0);
	gsl_vector_set_all(ss, initial_radius_v);

	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	int status = GSL_CONTINUE;
	int num_iter = 0;
	for( num_iter = 0;
	     status == GSL_CONTINUE && num_iter < max_num_iter;
	     ++num_iter ) 
	{
	    // nonzero return value is an error code.
	    status = gsl_multimin_fminimizer_iterate(s);
      	    if (status) 
		break;

	    double simplex_size = gsl_multimin_fminimizer_size (s);
	    status = gsl_multimin_test_size (simplex_size, abs_tol);

	    if (status == GSL_SUCCESS) {
		diag.level(10) << "converged to minimum at" << endl;
	    }

	    diag.level(10) << num_iter;
	    for( int i = 0; i < dimension; ++i) {
		diag.level(10) << " " << gsl_vector_get (s->x, i);
	    }
	    diag.level(10) << " f() = " << s->fval 
			  << " size = " << simplex_size
			  << endl;
	}

	double x = gsl_vector_get( s->x, 0 );
	double y = gsl_vector_get( s->x, 1 );
	diag.level(3) << "Ending value: " << f(x,y) << std::endl;
	final[i] = f.target_point(x,y);

	simplex_iterations.add_sample( num_iter );
    }

    if ( radius_is_relative ) {
	diag.level(1) << "Search radius mean: " << search_radius_stat << endl
		      << "    Minimum: " << search_radius_stat.min() << endl
		      << "    Maximum: " << search_radius_stat.max() << endl;
    }

    MNI::Statistic<double> displacement_magnitude;
    v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {
	//smap.get_target(v) = final[i];
	set_new_target( smap, v, final[i], displacement_magnitude );
    }

    diag.level(1) << "OPTIMIZE ";
    displacement_magnitude.print( diag.level(1) );

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}
