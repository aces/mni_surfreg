#include "surftracc.hpp"



/* Return normalized area of control facet: target over source.
 */
double normalized_area( const SurfaceMap& smap,
			ControlMesh::Halfedge_const_handle h )
{
    CGAL_assertion( h->next()->next()->next() == h );

    ControlMesh::Vertex_const_handle vA = h->vertex();
    ControlMesh::Vertex_const_handle vB = h->next()->vertex();
    ControlMesh::Vertex_const_handle vC = h->next()->next()->vertex();

    double s_area = signed_area( smap.get_source(vA)->point(),
				 smap.get_source(vB)->point(),
				 smap.get_source(vC)->point() );

    double t_area = signed_area( smap.get_target(vA).point(),
				 smap.get_target(vB).point(),
				 smap.get_target(vC).point() );
    return t_area / s_area;
}


/* Return true if any incident control facet has small
 * area ratio.
 */
bool is_incident_to_small_facet( const SurfaceMap& smap,
				 ControlMesh::Vertex_const_handle v,
				 double min_ratio )
{
    ControlMesh::Halfedge_const_handle h = v->halfedge();

    // Control mesh is only a HalfedgeDS; we cannot count on having
    // the nice circulators, so we have to do the loop the hard way.
    // We also cannot count on having a consistent orientation for
    // control mesh facets.  Come to think of it, can we even count
    // on having proper triangular facets? ....
    //
    do {
	if ( normalized_area(smap,h) < min_ratio )  return true;
	h = h->next()->opposite();
    } while( h != v->halfedge() );

    return false;
}


/* Return number of control facets mapped to small regions.
 */
int count_small_targets( const SurfaceMap& smap,
			 double min_ratio )
{
    ControlMesh::Halfedge_const_iterator 
	h = smap.control_mesh().halfedges_begin();

    // Since we're iterating over halfedges, we'll triple-count
    // the number of affected facets.

    int count = 0;
    CGAL_For_all( h, smap.control_mesh().halfedges_end() ) {
	if ( normalized_area(smap,h) < min_ratio )
	    ++count;
    }

    return count / 3;  
}


