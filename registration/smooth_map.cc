/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <cmath>

#include "surftracc.hpp"
#include <surflib/Search.hpp>


using namespace std;


/* Smooth the surface mapping by averaging the target location
 * at each control vertex.  Replace the target location at each
 * control vertex v by the weighted sum of v and some neighbourhood
 * centre location:
 *
 *         (1-weight) * v + weight * NeighbourhoodCentre(v)
 *
 * The target surface is assumed to be a sphere.
 * The weighted location is projected back to the sphere.
 */
void smooth_map( SurfaceMap& smap, 
		 double weight,
		 NeighbourhoodCentre& centre )
{
    const SurfaceMap::ControlMesh& control( smap.control_mesh() );
    std::vector<Point_3> t_point( control.size_of_vertices() );

    SurfaceMap::ControlMesh::Vertex_const_iterator v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {

	Vector_3 t_vector = weight * centre(v);
	t_vector = t_vector + (1.0-weight)*(smap.get_target(v).point() - CGAL::ORIGIN);
	t_vector = t_vector / std::sqrt(CGAL::to_double(t_vector*t_vector));

	t_point[i] = CGAL::ORIGIN + t_vector;
    }

    diag.level(3) << "Map Smoothing" << endl;
    MNI::Statistic<double> displacements;

    v = control.vertices_begin();
    for( int i = 0; v != control.vertices_end(); ++i,++v ) {
	set_new_target( smap, v, 
		        MNI::find_nearest_on_convex_surface<Surface>
		        ( smap.get_target(v).V0(), t_point[i] ),
			displacements );
    }

    diag.level(1) << "SMOOTHING ";
    displacements.print( diag.level(1) );
}
