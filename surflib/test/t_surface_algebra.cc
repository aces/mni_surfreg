#include "tester_approx.hpp"

#include "Surface.hpp"
#include <surflib/Point_s.hpp>
#include <surflib/Vector_s.hpp>

#include <CGAL/Polyhedron_incremental_builder_3.h>


using MNI::Point_s;
using MNI::Vector_s;



/* Create a mesh of the following topology
       ( all z=0 )

    y=10  p5----p4---p3
          |   / |\   |
          |  /  | \  |
          | /   |  \ |
          |/    |   \|
    y=0   p0---p1---p2
         x=0  x=10  x=20
*/

template <class _surface>
class MakeMesh : public CGAL::Modifier_base<typename _surface::HalfedgeDS>
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
	
	B.add_vertex( Point_3(  0,  0, 0 ));
	B.add_vertex( Point_3( 10,  0, 0 ));
	B.add_vertex( Point_3( 20,  0, 0 ));
	B.add_vertex( Point_3( 20, 10, 0 ));
	B.add_vertex( Point_3( 10, 10, 0 ));
	B.add_vertex( Point_3(  0, 10, 0 ));

	add_facet( B, 0,4,5 );
	add_facet( B, 0,1,4 );
	add_facet( B, 1,2,4 );
	add_facet( B, 2,3,4 );

	B.end_surface();

	CGAL_postcondition( hds.size_of_vertices() == n_vert );
	CGAL_postcondition( hds.size_of_halfedges() == n_hedge );
	//CGAL_postcondition( hds.size_of_facets() == n_facet );
    }
};



class t_surface_algebra : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;

    typedef Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point_3;
    typedef Surface::Traits::Vector_3       Vector_3;

    typedef Point_s<Surface>                Point_s;
    typedef Vector_s<Surface>               Vector_s;


    // Planar mesh from MakeMesh.
    Surface surf;
    Vertex_const_handle vert[6];

    // A tetrahedron.
    Surface t;
    Halfedge_const_handle h20;


    Halfedge_const_handle find_edge( int u, int v )
    {
	Halfedge_const_handle h = vert[v]->halfedge();
	do {
	    if ( &* h->opposite()->vertex() == &* vert[u] ) {
		return h;
	    }
	    h = h->next_on_vertex();
	} while ( h != vert[v]->halfedge() );
	CGAL_postcondition( false );
	return 0;
    }

    void setUp()
    {
	MakeMesh<Surface> mm;
	surf.delegate( mm );

	//std::copy( surf.vertices_begin(), surf.vertices_end(), vert );	
	Vertex_const_iterator v = surf.vertices_begin();
	for( int i = 0; v != surf.vertices_end(); ++v,++i ) {
	    vert[i] = &*v;
 	}

	h20 = t.make_tetrahedron( Point_3(0,0,0),
				  Point_3(1,0,0),
				  Point_3(0,1,0),
				  Point_3(0,0,1) );
    }


    // Test point + vector addition inside facet for points using
    // different facet coordinate frames.
    void inside_facet()
    {
	// All these points are the same: (12, 2, 0)
	Point_s p( find_edge(1,2), 0.2,0.2 );
	Point_s q( find_edge(2,4), 0.2,0.6 );
	Point_s r( find_edge(4,1), 0.6,0.2 );

	CPPUNIT_ASSERT_EQUAL( Point_3(12,2,0), p.point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(12,2,0), q.point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(12,2,0), r.point() );

	// Test add_vector(double,double)
	Point_s z = p.displace( 0.1, 0.2 );
	CPPUNIT_ASSERT_EQUAL( 0.3, z.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.4, z.t2() );

	Vector_s v( find_edge(1,2), 0.2, 0.2 );  // (2,2,0)
	
	CPPUNIT_ASSERT_EQUAL( Vector_3(2,2,0), v.vector() );
	
	// The sum should be (14,4,0)

	CPPUNIT_ASSERT_EQUAL( Point_3(14,4,0), (p+v).point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(14,4,0), (q+v).point() );
	CPPUNIT_ASSERT_EQUAL( Point_3(14,4,0), (r+v).point() );

	z = p + v;
	Vector_s d = z - p;
	CPPUNIT_ASSERT_EQUAL( v.vector(), d.vector() );

	z = q + v;
	d = z - p;
	CPPUNIT_ASSERT_EQUAL( v.vector(), d.vector() );

	z = r + v;
	d = z - r;
	CPPUNIT_ASSERT_EQUAL( v.vector(), d.vector() );
    }


    // Test addition that crosses facet boundaries.
    void across_facet()
    {
	Halfedge_const_handle h12 = find_edge(1,2);

	Point_s p( h12, 0.1,0.1 );  // (11,1,0)

	CPPUNIT_ASSERT_EQUAL( Point_3(19,6,0), 
			      (p + Vector_s(h12,  0.8,0.5)).point() );
	CPPUNIT_ASSERT_EQUAL( Point_3( 7,1,0),
			      (p + Vector_s(h12, -0.4,0 )).point() );
	CPPUNIT_ASSERT_EQUAL( Point_3( 2,6,0), 
			      (p + Vector_s(h12, -0.9,0.5)).point() );
    }

    
    // Starting at a vertex.
    void at_vertex()
    {
	Point_s p(vert[1]);

	// Test add_vector(double,double)
	Point_s z = p.displace( 0.3, 0.4 );
	CPPUNIT_ASSERT_EQUAL( 0.3, z.t1() );
	CPPUNIT_ASSERT_EQUAL( 0.4, z.t2() );

	Halfedge_const_handle h12 = find_edge(1,2);

	p = Point_s( h12, 0,0 );  // (10,0,0)

	z = p.displace(1,0);
	CPPUNIT_ASSERT_EQUAL( Point_3(20,0,0), z.point() );

	z = p.displace(-1,0);
	CPPUNIT_ASSERT_EQUAL( Point_3(0,0,0), z.point() );
    }	


    // Test in 3D using right tetrahedron.
    void tetrahedron()
    {
	//cout << h20->vertex()->point() << endl;
	//cout << h20->opposite()->vertex()->point() << endl;

	Point_s p(h20, 0.2,0.2 );
	CPPUNIT_ASSERT_EQUAL( Point_3(0.2,0.6,0), p.point() );

	Vector_s v_left(h20, 0.4, -0.4);
	CPPUNIT_ASSERT_EQUAL( Vector_3(-0.4,0,0), v_left.vector() );

	CPPUNIT_ASSERT_EQUAL( Point_3(0,0.6,0.2), (p+v_left).point() );

	Vector_s v_down(h20, 1, 0);
	CPPUNIT_ASSERT_EQUAL( Vector_3(0,-1,0), v_down.vector() );

	CPPUNIT_ASSERT_EQUAL( Point_3(0.2,0,0.4), (p+v_down).point() );
    }


    void midpoints()
    {
	// Test the sort of computation used when refining a mesh map.
	// (c.f. registration/SimplifySurface.icc)

	Vertex_const_iterator v = t.vertices_begin();
	for( ; v != t.vertices_end(); ++v ) {
	    Halfedge_const_handle h = v->halfedge();
	    do {
		// Halfedge h is (u,v)
		Vertex_const_handle u = h->opposite()->vertex();
		Point_s u_s(h,0,0);   // surface point u

		Vector_3 d = (v->point() - u->point()) / 2;
		double d1,d2;
		u_s.coordinate_frame().coordinates_of( d, &d1,&d2 );

		Point_s mid = u_s.displace(d1,d2);
		CPPUNIT_ASSERT_EQUAL( u->point() + d,  mid.point() );

		h = h->next_on_vertex();
	    } while ( h != v->halfedge() );
	}
    }	


    CPPUNIT_TEST_SUITE( t_surface_algebra );
    CPPUNIT_TEST( inside_facet );
    CPPUNIT_TEST( across_facet );
    CPPUNIT_TEST( at_vertex );
    CPPUNIT_TEST( tetrahedron );
    CPPUNIT_TEST( midpoints );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_surface_algebra );
