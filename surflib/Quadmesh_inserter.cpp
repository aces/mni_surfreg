/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <CGAL/Point_3.h>


// Debugging aids
#define CHECK_INDICES 1
#define REPORT_INDICES 0



template <class HalfedgeDS>
void MNI::ObjFile_inserter::Quadmesh<HalfedgeDS>::operator() ( HalfedgeDS& hds )
{
    typedef typename HalfedgeDS::Vertex  Vertex;
    typedef typename Vertex::Point       Point_3;

    // A quadmesh has m_points rows and n_points columns of points.

    int m_points = qmesh.m;
    int n_points = qmesh.n;
    int n_vert = m_points * n_points;

    // There are m_patches rows and n_patches columns of patches.
    // (m_patches = m_points-1 (open quadmesh) or m_points (closed).
    // Each patch will turn into a facet of the polyhedron.

    int m_patches, n_patches;
    get_quadmesh_n_objects( &qmesh, &m_patches, &n_patches );
    int n_face = m_patches * n_patches;

    // The number of half-edges is tricky, because there may be
    // 0, 1, or 2 boundary loops, so we can't use 4*n_face.
    // Instead, count the number of edges, and double them.
    // For each row we have one edge for each column of patches, for 
    // a total of (m_points * n_patches) "horizontal" edges.
    // Similarly, there are (n_points * m_patches) "vertical" edges.

    int n_hedge = 2 * ( m_points*n_patches + n_points*m_patches);

    // If we are going to split each patch into two triangles,
    // this adds two new halfedges per patch, and doubles the
    // number of faces.

    if ( split_patches ) {
	n_hedge += 2 * n_face;
	n_face *= 2;
    }
	
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B( hds, true );

    B.begin_surface( n_vert, n_face, n_hedge );

    for( int i = 0; i < n_vert; ++i ) {
	Point_3 p( Point_x(qmesh.points[i]),
		   Point_y(qmesh.points[i]),
		   Point_z(qmesh.points[i]) );
	B.add_vertex( p );
    }

    int facet_index[4];
    for( int i = 0; i < m_patches; ++i ) {
	for( int j = 0; j < n_patches; ++j ) {
	    get_quadmesh_patch_indices( &qmesh, i, j, facet_index );


	    CGAL_assertion( facet_index[0] != facet_index[1] );
	    CGAL_assertion( facet_index[0] != facet_index[2] );
	    CGAL_assertion( facet_index[0] != facet_index[3] );
	    CGAL_assertion( facet_index[1] != facet_index[2] );
	    CGAL_assertion( facet_index[1] != facet_index[3] );
	    CGAL_assertion( facet_index[2] != facet_index[3] );

	    if ( split_patches ) {
		B.begin_facet();
		B.add_vertex_to_facet( facet_index[0] );
		B.add_vertex_to_facet( facet_index[1] );
		B.add_vertex_to_facet( facet_index[2] );
		B.end_facet();
		B.begin_facet();
		B.add_vertex_to_facet( facet_index[2] );
		B.add_vertex_to_facet( facet_index[3] );
		B.add_vertex_to_facet( facet_index[0] );
		B.end_facet();
	    } else {
		B.begin_facet();
		B.add_vertex_to_facet( facet_index[0] );
		B.add_vertex_to_facet( facet_index[1] );
		B.add_vertex_to_facet( facet_index[2] );
		B.add_vertex_to_facet( facet_index[3] );
		B.end_facet();

	    }
	}
    }

    B.end_surface();
}
