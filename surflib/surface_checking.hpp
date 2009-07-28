/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFLIB_SUBDIVIDED_ICOSAHEDRON_HPP
#define SURFLIB_SUBDIVIDED_ICOSAHEDRON_HPP

#include <cmath>
#include <exception>
#include <CGAL/circulator.h>
#include <CGAL/squared_distance_3.h>

namespace MNI {


/*! \brief Assert that surface geometry is sphere.
 *
 * Check that the 3-space location of each point of surface \a s 
 * lies on a sphere with centre at \a centre and radius \a radius
 * (plus or minus \a tolerance).
 *
 * \throws std::runtime_error if some point is not within tolerance.
 */
template <class Surface_>
void assert_is_sphere_geometry( const Surface_& s, 
				const typename Surface_::Point_3& 
				centre = CGAL::ORIGIN,
				double radius = 1.0,
				double tolerance = 1e-4 )
{
    typename Surface_::Vertex_const_iterator v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
	double d2 = CGAL::to_double(CGAL::squared_distance( v->point(), 
							    centre ));
	if ( std::abs( std::sqrt(d2) - radius ) > tolerance )
	    throw std::runtime_error( "surface geometry is not sphere." );
    }
}


/*! \brief Assert that surface is triangulated sphere.
 *
 * Check that the graph topology of the surface is a valid
 * tesselation of a sphere.  This routine will check that
 * the surface
 * 
 * - satisfies CGAL::Polyhedron_3::is_valid()
 * - has no borders
 * - satisfies Euler equation for genus 0
 * - all facets are triangles
 *
 * \throws std::runtime_error if any check fails.
 */
template <class Surface_>
void assert_is_sphere_topology( const Surface_& s )
{
    if (!s.is_valid())
	throw std::runtime_error( "mesh is not valid" );

    typename Surface_::Halfedge_const_iterator h = s.halfedges_begin();
    CGAL_For_all( h, s.halfedges_end() ) {
	if (h->is_border())
	    throw std::runtime_error( "mesh has a border" );

	if ( h != h->next()->next()->next() )
	    throw std::runtime_error( "mesh has non-triangular facet" );
    }

    if ( s.size_of_vertices() - s.size_of_halfedges()/2 + s.size_of_facets()
	 != 2 )
	throw std::runtime_error( "mesh does not satisfy V-E+F=2" );

    // This should be redundant, since we already know that each facet
    // is a triangle and there are no borders.
    if ( s.size_of_halfedges() != 3*s.size_of_facets() )
	throw std::runtime_error( "mesh does not satisfy 2E = 3F" );
}


/*! \brief Compute size of coarsened mesh.
 *
 * \return size of vertex list for a 4-to-1 merging of facets.
 *
 * \pre Surface \p s must be a 1-to-4 refinement of a previous surface
 * mesh.  
 *
 * \throws std::runtime_error if precondition not met.
 */
template <class Surface_>
typename Surface_::size_type coarser_size_of_vertices( const Surface_& s )
{
    /* Assume that the mesh is genus 0, so Euler's equation is E-V+F=2.
     * Assume that it is maximally-triangulated, so 2E = 3F.
     * Eliminate E, leaving F = 2V-4.
     *
     * Denote the original mesh as "1" and the refined mesh as "2",
     * then we have F_2 = 4F_1, so 
     *    2 V_2 - 4 = 4( 2 V_1 - 4 )   ==>   V_1 = (V_2 + 6) / 4
     */

    // Smallest possible is refinement of tetrahedron, with 9 vertices.
    if ( s.size_of_vertices() < 9 )
	throw std::runtime_error( "mesh is too small" );

    if ( (s.size_of_vertices() + 6) % 4 != 0 )
	throw std::runtime_error( "mesh is not 1-to-4 facet refinement" );

    return (s.size_of_vertices()+6)/4;
}

}


#endif

