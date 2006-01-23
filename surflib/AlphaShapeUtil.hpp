#ifndef MNI_SURFLIB_ALPHA_SHAPE_UTIL_HPP  // -*- C++ -*-
#define MNI_SURFLIB_ALPHA_SHAPE_UTIL_HPP

#include <iostream>
#include <iterator>
#include <map>

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


namespace MNI {

/**
 * Classify vertices of Surface s according to the alpha shape
 * of its point set.
 *
 * The template type Surface must be a CGAL::Polyhedron_3 subclass
 * that uses a Cartesian kernel (not Simple_cartesian).  In addition, 
 * the Vertex class must have the following three fields: 
 *   int label
 *   double low
 *   double high
 *
 * After calling this method, those fields are modified as follows.
 *
 * 1. If classify is true, then label is set to the vertex classification
 * type an instance of Alpha_shape_3::Classification_type.  This is one
 * of: EXTERIOR, REGULAR (on boundary of alpha complex), or INTERIOR.
 * If classify is false, then label is set to 0.
 *
 * 2. Value low is the lower threshold for this vertex appearing on the
 * boundary of the alpha complex.  For alpha values below this, the vertex 
 * is EXTERIOR.
 *
 * 3. Value high is the upper threshold for this vertex appearing on the
 * boundary of the alpha complex.  For alpha values above this, the vertex
 * is INTERIOR.  If the vertex is on the convex hull, the upper threshold
 * is infinite and high is set to -1.
 *
 * For alpha values on the interval [low,high], the vertex is REGULAR,
 * i.e., it forms part of the boundary.
 *
 * \param s is the surface to classify
 * \param classify set true to have each vertex classified
 * \param alpha if positive, the alpha value to use in classification; 
                if negative, use an alpha value that produces one connected component
 */
template< class Surface >
void classifyVertices( Surface& s, 
		       bool classify, 
		       double alpha )
{
    typedef typename Surface::Vertex_handle       Vertex_handle;
    typedef typename Surface::Vertex_iterator     Vertex_iterator;
    typedef typename Surface::Traits              K;

    typedef typename K::RT RT;

    typedef CGAL::Alpha_shape_euclidean_traits_3<K> Gt;
    typedef CGAL::Alpha_shape_vertex_base_3<Gt> Vb;
    typedef CGAL::Triangulation_cell_base_3<Gt> Df;
    typedef CGAL::Alpha_shape_cell_base_3<Gt, Df>  Fb;
    typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
    typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Triangulation_3;
    typedef CGAL::Alpha_shape_3<Triangulation_3>  Alpha_shape_3;

    typedef CGAL::Unique_hash_map<const RT*,Vertex_handle> Point_map;

    for( Vertex_iterator v = s.vertices_begin(); v != s.vertices_end(); ++v )
	v->label = -1;

    /* The alpha shape is initialized with a copy of the points of s.
     * In order to make the link from an alpha shape vertex to 
     * a polyhedron vertex, we require that the points use a handle/rep
     * implementation and simply use a map from the pointer of the X
     * coordinate to the surface vertex.
     *
     * This means that K must be Cartesian rather than Simple
     * Cartesian.
     */
    Point_map point_to_vertex;
    Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
        point_to_vertex[& v->point().x()] = v;
    }
    
    std::cout << "Creating alpha shape ..." << std::flush;
    Alpha_shape_3 shape( s.points_begin(), s.points_end() );
    std:: cout << "done.\n";

    if (classify && alpha < 0) {
	std::cout << "Optimal alpha value for 1 component is " << std::flush;
	alpha = CGAL::to_double( *shape.find_optimal_alpha(1) );
	std::cout << alpha << "\n";
    }

    if ( alpha > 0 ) 
	std::cout << "Number of solid components: " 
		  << shape.number_of_solid_components(alpha) << "\n";


    /* Set classification value in vertices of s.
     */
    typename Triangulation_3::Finite_vertices_iterator 
	vi = shape.finite_vertices_begin();
    CGAL_For_all( vi, shape.finite_vertices_end() ) {
	CGAL_assertion( point_to_vertex.is_defined( & vi->point().x() ) );
	Vertex_handle v = point_to_vertex[ & vi->point().x() ];
	CGAL_assertion( v != point_to_vertex.default_value() );

	typename Alpha_shape_3::Alpha_status* as = vi->get_alpha_status();

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

}

#endif
