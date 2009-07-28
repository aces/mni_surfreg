/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <CGAL/Polyhedron_incremental_builder_3.h>


/* Create a mesh of the following topology
    
                v2
	       /  \
	      /    \
	     /      \
	    /        \
	   /          \
	  v0----------v1
*/

template <class _surface>
class MakeMesh_0 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 3;
	typename _surface::size_type n_facet = 2;
	typename _surface::size_type n_hedge = 6;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, 0 ));
	B.add_vertex( Point_3( 4, 0, 0 ));
	B.add_vertex( Point_3( 2, 2, 0 ));

	add_facet( B, 0,1,2 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};



/* Create a mesh of the following topology 
   (1-to-4 refinement of previous)
    
                v2
	       /  \
	      /    \
	     v5----v4
	    /  \  /  \
	   /    \/    \
	  v0----v3----v1
*/

template <class _surface>
class MakeMesh_1 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 6;
	typename _surface::size_type n_facet = 5;
	typename _surface::size_type n_hedge = 18;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, 0 ));
	B.add_vertex( Point_3( 4, 0, 0 ));
	B.add_vertex( Point_3( 2, 2, 0 ));

	B.add_vertex( Point_3( 2, 0, 0 ));
	B.add_vertex( Point_3( 3, 1, 0 ));
	B.add_vertex( Point_3( 1, 1, 0 ));

	add_facet( B, 0,3,5 );
	add_facet( B, 3,4,5 );
	add_facet( B, 3,1,4 );
	add_facet( B, 4,2,5 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};


/* Create a mesh of the following topology 
   (same as MakeMesh_1; a 1-to-4 refinement of MakeMesh_0)
    
                v2
	       /  \
	      /    \
	     v5----v4
	    /  \  /  \
	   /    \/    \
	  v0----v3----v1

   The difference with MakeMesh_1 is that the mesh is non-planar; vertices
   0,1, and 2 are lifted off the z=0 plane.
*/

template <class _surface>
class MakeMesh_2 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 6;
	typename _surface::size_type n_facet = 5;
	typename _surface::size_type n_hedge = 18;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, -1 ));
	B.add_vertex( Point_3( 4, 0, 1 ));
	B.add_vertex( Point_3( 2, 2, 3 ));

	B.add_vertex( Point_3( 2, 0, 0 ));
	B.add_vertex( Point_3( 3, 1, 0 ));
	B.add_vertex( Point_3( 1, 1, 0 ));

	add_facet( B, 0,3,5 );
	add_facet( B, 3,4,5 );
	add_facet( B, 3,1,4 );
	add_facet( B, 4,2,5 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};



/* Create a mesh of the following topology 
    
                v2
	       /| \
	      / |  \
	     /  v3  \
	    / /   \  \
	   / /     \  \
	  v0----------v1

*/

template <class _surface>
class MakeMesh_3 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 6;
	typename _surface::size_type n_facet = 5;
	typename _surface::size_type n_hedge = 18;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, -1 ));
	B.add_vertex( Point_3( 4, 0, 1 ));
	B.add_vertex( Point_3( 2, 2, 3 ));

	B.add_vertex( Point_3( 2, 0, 0 ));
	B.add_vertex( Point_3( 3, 1, 0 ));
	B.add_vertex( Point_3( 1, 1, 0 ));

	add_facet( B, 0,3,5 );
	add_facet( B, 3,4,5 );
	add_facet( B, 3,1,4 );
	add_facet( B, 4,2,5 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};



/* Create a right tetrahedron.
*/

template <class _surface>
class MakeTetra_0 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 4;
	typename _surface::size_type n_facet = 4;
	typename _surface::size_type n_hedge = 12;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, 0 ));
	B.add_vertex( Point_3( 2, 0, 0 ));
	B.add_vertex( Point_3( 0, 2, 0 ));
	B.add_vertex( Point_3( 0, 0, 2 ));

	add_facet( B, 1,2,3 );
	add_facet( B, 3,2,0 );
	add_facet( B, 3,0,1 );
	add_facet( B, 0,2,1 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};



/* Create a right tetrahedron subdivision.
*/

template <class _surface>
class MakeTetra_1 : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
{
public:
    typedef 
    CGAL::Polyhedron_incremental_builder_3<typename _surface::HalfedgeDS> 
    IncrBuild;

    typedef typename _surface::Point_3 Point_3;

    void add_facet( IncrBuild& B, int v0, int v1, int v2 )
    {
	B.begin_facet();
	B.add_vertex_to_facet( v0 );
	B.add_vertex_to_facet( v1 );
	B.add_vertex_to_facet( v2 );
	B.end_facet();
    }

    void operator() ( typename _surface::HalfedgeDS& hds )
    {
	typename _surface::size_type n_vert = 10;
	typename _surface::size_type n_facet = 16;
	typename _surface::size_type n_hedge = 3*n_facet;

	IncrBuild B( hds, true );
	B.begin_surface( n_vert, n_facet, n_hedge );
	
	B.add_vertex( Point_3( 0, 0, 0 ));
	B.add_vertex( Point_3( 2, 0, 0 ));
	B.add_vertex( Point_3( 0, 2, 0 ));
	B.add_vertex( Point_3( 0, 0, 2 ));

	B.add_vertex( Point_3( 1, 0, 0 ));
	B.add_vertex( Point_3( 1, 1, 0 ));
	B.add_vertex( Point_3( 0, 1, 0 ));

	B.add_vertex( Point_3( 1, 0, 1 ));
	B.add_vertex( Point_3( 0, 1, 1 ));
	B.add_vertex( Point_3( 0, 0, 1 ));

	add_facet( B, 1,5,7 );
	add_facet( B, 5,8,7 );
	add_facet( B, 5,2,8 );
	add_facet( B, 7,8,3 );


	add_facet( B, 2,6,8 );
	add_facet( B, 6,9,8 );
	add_facet( B, 6,0,9 );
	add_facet( B, 8,9,3 );


	add_facet( B, 0,4,9 );
	add_facet( B, 4,7,9 );
	add_facet( B, 4,1,7 );
	add_facet( B, 9,7,3 );


	add_facet( B, 0,6,4 );
	add_facet( B, 4,5,1 );
	add_facet( B, 4,6,5 );
	add_facet( B, 6,2,5 );


	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
    }
};
