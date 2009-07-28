/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
#include <functional>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>

#include <surflib/Diagnostic_ostream.hpp>
#include <surflib/Statistic.hpp>
#include <surflib/SurfaceMap_poly.hpp>


/* Put a scalar value in each vertex.
 */
template <class Refs, class Traits>
struct Surftracc_vertex
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
					  typename Traits::Point_3>
{
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_base<Refs, 
					 CGAL::Tag_true, 
					 Point> Vertex_base;

    Surftracc_vertex( const Point& p ) : Vertex_base(p) { scalar = 0; }
    Surftracc_vertex() { scalar = 0; }

    double scalar;
};


struct Polyhedron_items_3 : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Vertex_wrapper {
	typedef typename Traits::Point_3       Point;
	typedef Surftracc_vertex<Refs,Traits>  Vertex;
    };
};


//
// Geometric types
//

typedef CGAL::Simple_cartesian<double>               Surface_kernel;
typedef Surface_kernel::Point_3                      Point_3;
typedef Surface_kernel::Vector_3                     Vector_3;
typedef CGAL::Polyhedron_3<Surface_kernel, 
			   Polyhedron_items_3,
			   CGAL::HalfedgeDS_vector>  Surface;
typedef MNI::SurfaceMap_poly<Surface>                SurfaceMap;
typedef SurfaceMap::ControlMesh                      ControlMesh;



//! Functor to retrieve scalar from vertex.
struct ScalarData
    : public std::unary_function<Surface::Vertex_const_handle, double>
{
    double operator() ( const Surface::Vertex_const_handle& v ) const
    {
        return v->scalar;
    }
};


// Send diagnostics here
extern MNI::Diagnostic_ostream diag;


// Optimize with GSL.  See optimize_map_gsl.cpp.
class ObjectiveFunction;
void optimize_map_gsl( SurfaceMap& smap,
		       ObjectiveFunction& f,
		       bool radius_is_relative,
		       const std::vector<double>& penalty_radius,
		       double initial_radius,
		       double abs_tol,
		       int max_num_iter,
		       MNI::Statistic<double>& simplex_iterations );


// Smooth map.  See smooth_map.cc.
//
class NeighbourhoodCentre
{
public:
    virtual const Vector_3 operator()( ControlMesh::Vertex_const_handle ) = 0;
    virtual ~NeighbourhoodCentre() {}
};

void smooth_map( SurfaceMap& smap, 
		 double neighbourhood_weight,
		 NeighbourhoodCentre& );


// Signed area calculations.  See mesh_area.cpp.

/* Return signed facet area.  The orientation is positive
 * if the vertices are clockwise as seen from the origin.
 */
inline
double signed_area( const Point_3& A,
		    const Point_3& B,
		    const Point_3& C )
{
    Vector_3 K = CGAL::cross_product( B-A, C-A );
    double area = 0.5 * std::sqrt( CGAL::to_double(K*K) );

    if ( K*(A-CGAL::ORIGIN) < 0 )
	area = -area;

    return area;
}

/* Return normalized area of control facet: target over source.
 */
double normalized_area( const SurfaceMap& smap,
			ControlMesh::Halfedge_const_handle h );


/* Return true if any incident control facet has small
 * area ratio.
 */
bool is_incident_to_small_facet( const SurfaceMap& smap,
				 ControlMesh::Vertex_const_handle v,
				 double min_ratio = 1e-3 );

/* Return number of control facets mapped to small regions.
 */
int count_small_targets( const SurfaceMap& smap,
			 double min_ratio = 1e-3 );



/*! Update smap[v] to location p_s, keeping statistics
 * on the displacement magnitude.
 */
inline
void set_new_target( SurfaceMap& smap, 
		     SurfaceMap::ControlMesh::Vertex_const_iterator v,
		     const MNI::Point_s<Surface>& p_s,
		     MNI::Statistic<double>& displacement_mag )
{
    // Sanity check that the target has not moved too much.
    Point_3 p = smap.get_target(v).point();
    Point_3 q = p_s.point();
    double mag = std::sqrt(CGAL::to_double((p-q)*(p-q)));

    diag.level(3) << "Move target from: " << p << std::endl
		  << "              to: " << q << std::endl
		  << "       magnitude: " << mag
		  << std::endl << std::endl;

    CGAL_assertion( CGAL::to_double( (p-CGAL::ORIGIN)*(q-CGAL::ORIGIN)) > 0 );

    displacement_mag.add_sample( mag );
    smap.get_target(v) = p_s;
}


