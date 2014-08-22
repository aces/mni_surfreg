/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
#include <surflib/LinearInterpolation.hpp>
#include <surflib/rotation.hpp>

/* Matched pair of neighbourhoods on source and target unit spheres.
 * The neighbourhoods are parameterized on a unit disc lifted to
 * the z>0 hemisphere.
 *
 */
class NeighbourhoodPair
{
public:

    //! The number of points in the neighbourhood.
    static const int neighbourhood_size = 1+8+16;   //  +32;

    //! Build structure with no mapping.
    NeighbourhoodPair( SurfaceMap& smap, double radius );

    //! Set up mapping between v and smap[v].
    void set_source_vertex( ControlMesh::Vertex_handle v );

    //! Perturb the target of the mapping.
    void perturb_target_point( const Point_3& perturbed_pole );


    //! Return neighbourhood point.
    const Point_3 get_neighbourhood( int i ) 
    { return neighbourhood[i]; }

    //! Return source point.
    const Point_3 get_source( int i ) 
    { return rotate_to_surface(neighbourhood[i]); }

    //! Return target point.
    const Point_3 get_target( int i ) 
    { return rotate_to_target( rotate_to_perturbed_pole
			       (neighbourhood[i]) ); }


private:
    NeighbourhoodPair( const NeighbourhoodPair& );
    NeighbourhoodPair operator=( const NeighbourhoodPair& );

private:
    SurfaceMap& smap;
    double neighbourhood_radius;

    // Store neighbourhood sample locations (on northern hemisphere)
    Point_3 neighbourhood[neighbourhood_size];

    // Store rotations from northern hemisphere neighbourhood to 
    // source location of v, target location smap[v], and from
    // the pole to the perturbed_pole.
    //
    MNI::Rotation<Surface_kernel> rotate_to_surface;
    MNI::Rotation<Surface_kernel> rotate_to_target;
    MNI::Rotation<Surface_kernel> rotate_to_perturbed_pole;
};



/* Objective function for the matching step at control vertex v.  The
 * objective is a traditional two-term function.  The data term
 * computes dissimilarity between two sets of samples using, for
 * example, correlation coefficient.  The "model" term (regularizer)
 * is a penalty on displacement of v.
 *
 * The source and target surface geometry is assumed to be a unit
 * sphere.  The data term is computed on a spherical cap of specified
 * radius centered at v.  At the moment, this neighbourhood mesh is
 * fixed to consist of neighbourhood_size points.
 *
 * The model term is a penalty that increases with the radial
 * displacement of v, tending towards infinity at a finite radius.
 * Thus the search is effectively bounded.
 */
class ObjectiveFunction
{
public:
    //! Construct the search structures.
    ObjectiveFunction( SurfaceMap& smap,
		       double neighbourhood_radius,
		       double penalty_ratio );

    //! Set radius for penalty term.
    void set_penalty_radius( double penalty_radius );

    //! Set the vertex to be optimized.
    void set_vertex( ControlMesh::Vertex_handle v );

    //! Compute objective function.
    double operator() ( double x, double y );

    //! Extract target point.
    const MNI::Point_s<Surface> target_point( double x, double y );

private:
    double data_term( double x, double y );
    double model_term( double x, double y );

    //! Normalize vector.
    static void normalize( Vector_3& x )
    {
	x = x / std::sqrt( CGAL::to_double(x*x) );
    }
	

    //! Lift point (x,y) to northern unit hemisphere.
    static const Point_3 lift_to_sphere( double x, double y )
    {
	return Point_3( x, y, std::sqrt(1-x*x-y*y) );
    }

private:
    NeighbourhoodPair np;

    SurfaceMap& smap;
    ControlMesh::Vertex_handle v;

    ScalarData f_sd;
    MNI::LinearInterpolation<Surface,ScalarData> interpolate;

    double penalty_radius, penalty_radius_sq;
    double model_ratio;  // this is penalty_ratio

    // Store sample values for source & target.
    double source_sample[NeighbourhoodPair::neighbourhood_size];
    double target_sample[NeighbourhoodPair::neighbourhood_size];
};



//! Compute objective function.
inline double ObjectiveFunction::operator() ( double x, double y )
{

    double r_sq = x*x + y*y;
    if ( r_sq >= penalty_radius_sq )
	return std::numeric_limits<double>::max();

    // model_ratio is called penalty_ratio in the main program.

    return data_term(x,y) + model_ratio * model_term(x,y);
}


inline double ObjectiveFunction::model_term( double x, double y )
{
    double r_sq = x*x + y*y;
    return - std::log( 1 - r_sq/penalty_radius_sq );
}


