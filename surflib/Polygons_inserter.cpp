#include <exception>
#include <CGAL/Point_3.h>
#include <CGAL/HalfedgeDS_items_decorator.h>


// Debugging aids
#define CHECK_INDICES 1
#define REPORT_INDICES 0



/*! HalfedgeDS modifier that inserts the BIC polygons structure into the
 * CGAL HalfedgeDS.
 **/
template <class HalfedgeDS>
void MNI::ObjFile_inserter::Polygons<HalfedgeDS>::operator() ( HalfedgeDS& hds )
{
    typedef typename HalfedgeDS::Vertex  Vertex;
    typedef typename Vertex::Point       Point_3;

    int n_vert = poly.n_points;
    int n_face = poly.n_items;

    int n_hedge = 0;
    for( int i = 0; i < n_face; ++i ) {
	n_hedge += GET_OBJECT_SIZE( poly, i );
    }


    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B( hds, true );

    B.begin_surface( n_vert, n_face, n_hedge );

    for( int i = 0; i < n_vert; ++i ) {
	Point_3 p( Point_x(poly.points[i]),
		   Point_y(poly.points[i]),
		   Point_z(poly.points[i]) );
	B.add_vertex( p );
    }

    for( int i = 0; i < n_face; ++i ) {
	int face_size = GET_OBJECT_SIZE( poly, i );
#if REPORT_INDICES
	cout << "Facet " << i << " (size " << face_size << "): ";
#endif

#if CHECK_INDICES
	if ( face_size == 3 ) {
	    int v0 = poly.indices[POINT_INDEX(poly.end_indices,i,0)];
	    int v1 = poly.indices[POINT_INDEX(poly.end_indices,i,1)];
	    int v2 = poly.indices[POINT_INDEX(poly.end_indices,i,2)];

	    CGAL_assertion( v0 != v1 );
	    CGAL_assertion( v1 != v2 );
	    CGAL_assertion( v2 != v0 );
	}
#endif
	B.begin_facet();
	for( int j = 0; j < face_size; ++j ) {
	    int v = poly.indices[POINT_INDEX(poly.end_indices,i,j)];
#if REPORT_INDICES
	    cout << v << " ";
#endif
	    B.add_vertex_to_facet( v );
	}
#if REPORT_INDICES
	cout << "." << endl;
#endif
	B.end_facet();
    }

    B.end_surface();
}


namespace MNI {
namespace ObjFile_inserter {


/*! Insert edges of a BIC polygons structure into the CGAL HalfedgeDS.
 **/
template <class HDS>
void 
insert_polygons_into_hds ( const polygons_struct& poly,
			   HDS& hds )
{
    typedef typename HDS::Halfedge_handle  Halfedge_handle;
    typedef typename HDS::Vertex           Vertex;
    typedef typename HDS::Vertex_handle    Vertex_handle;

    int n_vert = poly.n_points;
    int n_face = poly.n_items;

    int n_hedge = 0;
    for( int i = 0; i < n_face; ++i ) {
	n_hedge += GET_OBJECT_SIZE( poly, i );
    }

    // FIXME: loosen this restriction to allow appending
    // to existing HDS.
    CGAL_assertion( hds.size_of_vertices() == 0 );
    hds.reserve( n_vert, n_hedge, 0 );

    std::vector<Vertex_handle> v_list(n_vert);
    Vertex v_initial;

    // make sure it has a singular value for halfedge
    Halfedge_handle h_sentinel;
    v_initial.set_halfedge(h_sentinel); 

    for( int i = 0; i < n_vert; ++i ) {
	v_list[i] = hds.vertices_push_back( v_initial );
    }

    typename HDS::Halfedge he[2];
    CGAL::HalfedgeDS_items_decorator<HDS> D;

    // Insert all halfedges, linking each to its incident vertex,
    // and linking the incident vertex to one of its halfedges
    // (the last created).
    for( int i = 0; i < n_face; ++i ) {
	int face_size = GET_OBJECT_SIZE( poly, i );

	//cout << "Face " << i << ", size = " << face_size << endl;

	int u = poly.indices[POINT_INDEX(poly.end_indices,i,face_size-1)];
	for( int j = 0; j < face_size; ++j ) {
	    int v = poly.indices[POINT_INDEX(poly.end_indices,i,j)];

	    //cout << "    Vertex " << v << endl;

	    // Generate halfedge pairs.  Do so only once for a given edge.
	    if ( u < v ) {
		Halfedge_handle h = hds.edges_push_back( he[0],he[1] );
		Halfedge_handle v_h = v_list[v]->halfedge();
		if ( v_h == h_sentinel ) {
		    //cout << "    Close tip at " << v << endl;
		    D.close_tip( h, v_list[v] );
		    v_list[v]->set_halfedge( h );
		} else {
		    //cout << "    Insert tip at " << v << endl;
		    CGAL_assertion( v_h->vertex() == v_list[v] );
		    D.insert_tip( h, v_h );
		}
		CGAL_assertion( h->vertex() == v_list[v] );

		// Repeat for h->opposite() and vertex u
		Halfedge_handle g = h->opposite();
		Halfedge_handle u_h = v_list[u]->halfedge();
		if ( u_h == h_sentinel ) {
		    //cout << "    Close tip at " << u << endl;
		    D.close_tip( g, v_list[u] );
		    v_list[u]->set_halfedge( g );
		} else {
		    //cout << "    Insert tip at " << u << endl;
		    CGAL_assertion( u_h->vertex() == v_list[u] );
		    D.insert_tip( g, u_h );
		}
		CGAL_assertion( g->vertex() == v_list[u] );
	    }
	    u = v;
	}
    }	
    
    for( int i = 0; i < n_vert; ++i ) {
	if ( v_list[i]->halfedge() == h_sentinel )
	    throw std::runtime_error
		( "MNI::insert_polygons_into_hds(): isolated vertex." );
    }
}


}
}
