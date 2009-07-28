/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <stdexcept>

#include "surftracc.hpp"
#include "data_term.hpp"
#include "objective.hpp"

#include <surflib/Search.hpp>


/*! The constructed object holds no mapping until set_source_vertex()
 *  is called.  The \a neighbourhood_radius must be <= 1.
 */
NeighbourhoodPair::NeighbourhoodPair( SurfaceMap& sm,
				      double neighbourhood_radius )
    : smap(sm), neighbourhood_radius(neighbourhood_radius)
{
    if ( neighbourhood_radius > 1 )
	throw std::domain_error( "NeighbourhoodPair: neighbourhood radius too large" );

    // Compute the neighbourhood stencil.
    //
    CGAL_assertion( neighbourhood_size == 1 + 8 + 16 );

    neighbourhood[0] = Point_3(0,0,1);

    double angle = 0;
    double angle_increment = 2.0*M_PI / 8.0;
    for( int i = 0; i < 8; ++i, angle += angle_increment ) {
	double x = std::cos(angle) * neighbourhood_radius/2.0;
	double y = std::sin(angle) * neighbourhood_radius/2.0;
	double z = std::sqrt( 1.0 - x*x - y*y );
	neighbourhood[1+i] = Point_3(x,y,z);
    }

    angle = 0;
    angle_increment = 2.0*M_PI / 16.0;
    for( int i = 0; i < 16; ++i, angle += angle_increment ) {
	double x = std::cos(angle) * neighbourhood_radius;
	double y = std::sin(angle) * neighbourhood_radius;
	double z = std::sqrt( 1.0 - x*x - y*y );
	neighbourhood[9+i] = Point_3(x,y,z);
    }
}


void NeighbourhoodPair::set_source_vertex( ControlMesh::Vertex_handle v )
{
    rotate_to_surface 
	= MNI::make_pole_rotation( smap.get_source(v)->point() );

    rotate_to_target
	= MNI::make_pole_rotation( smap.get_target(v).point() );

    rotate_to_perturbed_pole = MNI::Rotation<Surface_kernel>();
}


void NeighbourhoodPair::perturb_target_point( const Point_3& perturbed_pole )
{
    rotate_to_perturbed_pole = MNI::make_pole_rotation( perturbed_pole );
}




ObjectiveFunction::ObjectiveFunction( SurfaceMap& sm,
				      double nr,
				      double mr )
    : np(sm,nr),
      smap(sm), 
      interpolate(f_sd),
      model_ratio(mr)
{
}


void ObjectiveFunction::set_penalty_radius( double pr )
{
    penalty_radius = pr;

    if ( penalty_radius > 1 )
	throw std::domain_error( "ObjectiveFunction: penalty radius too large" );

    penalty_radius_sq = penalty_radius*penalty_radius;
}


/*! Compute source samples.
 */
void ObjectiveFunction::set_vertex( ControlMesh::Vertex_handle v_ )
{
    v = v_;

    np.set_source_vertex(v);
    Surface::Vertex_const_handle v_source = smap.get_source(v);

    // In principle, we do not need to assume sample 0 is the central
    // vertex.  We do this now so that we can compare the result with
    // output that uses the previous code, as a check.
    source_sample[0] = v_source->scalar;

    for( int i = 1; i < NeighbourhoodPair::neighbourhood_size; ++i ) {
	MNI::Point_s<Surface>
	    p_s = MNI::find_nearest_on_convex_surface<Surface>
	    ( v_source, np.get_source(i) );
	source_sample[i] = interpolate(p_s);
    }
}



/*! Use to obtain the optimal location of the vertex.
 */
const MNI::Point_s<Surface>
ObjectiveFunction::target_point( double x, double y )
{
    // FIXME: We're making the assumption (again) that point 0
    // will be the central vertex.

    np.perturb_target_point( lift_to_sphere(x,y) );
    return MNI::find_nearest_on_convex_surface<Surface>
	( smap.get_target(v).V0(), np.get_target(0) );
}    


/* This function is only called when the perturbation is
 * less than the penalty radius, which is less than 1.
 */
double ObjectiveFunction::data_term( double x, double y )
{
    np.perturb_target_point( lift_to_sphere(x,y) );

    Point_3 pt0_3 = np.get_target(0);
    MNI::Point_s<Surface> pt0_s = MNI::find_nearest_on_convex_surface<Surface>
        ( smap.get_target(v).V0(), pt0_3 );
    target_sample[0] = interpolate(pt0_s);

    for( int i = 1; i < NeighbourhoodPair::neighbourhood_size; ++i ) {
	Point_3 pt_3 = np.get_target(i);
	MNI::Point_s<Surface> pt_s 
	    = MNI::find_nearest_on_convex_surface<Surface>( pt0_s.V0(), pt_3 );
	target_sample[i] = interpolate(pt_s);
    }

    // Return value in range [0,1], with 0 being optimal match.
    return 1 - correlation_coefficient( source_sample, 
					source_sample+NeighbourhoodPair::neighbourhood_size,
					target_sample );
}


