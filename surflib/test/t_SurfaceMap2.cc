#include "tester.hpp"

#include "Surface.hpp"
#include <surflib/SurfaceMap_hds.hpp>

using namespace std;
using MNI::Point_s;


// This is a strict comparison for equality -- must have
// the same halfedge().
template <class _surface>
bool operator==( const Point_s<_surface>& a,
		 const Point_s<_surface>& b )
{
    return &* a.halfedge() == &* b.halfedge()
	&& a.t1() == b.t1()
	&& a.t2() == b.t2();
}



class t_SurfaceMap2 : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef MNI::SurfaceMap_hds<Surface>    SMap;
    typedef SMap::ControlMesh               ControlMesh;

    typedef ControlMesh::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point;
    typedef Point_s<Surface>                Point_s;



    // Points with which to construct surfaces.
    Point a,b,c;
    Surface s,t;
    ControlMesh cmesh;
    
    void setUp()
    {
	a = Point(0,0,0);
	b = Point(1,0,0);
	c = Point(0,1,0);

	s.make_triangle(a,b,c);
	t.make_triangle(a,b,c);

	cmesh.reserve(3,6,2);
	cmesh.vertices_push_back( ControlMesh::Vertex() );
	cmesh.vertices_push_back( ControlMesh::Vertex() );
	cmesh.vertices_push_back( ControlMesh::Vertex() );
    }


    void constructors() 
    {
	SMap map(s,t,cmesh);

	// Check that map is identity
	Surface::Vertex_const_iterator si = s.vertices_begin();
	Surface::Vertex_const_iterator ti = t.vertices_begin();
	SMap::ControlMesh::Vertex_const_iterator 
	    ci = map.control_mesh().vertices_begin();

	for( ; si != s.vertices_end(); ++si,++ti,++ci ) {
	    CPPUNIT_ASSERT( map.get_source(ci) == si );
	    Point_s loc = map.get_target(ci);
	    CPPUNIT_ASSERT( loc.t1() == 0 );
	    CPPUNIT_ASSERT( loc.t2() == 0 );
	    CPPUNIT_ASSERT( &*loc.halfedge()->opposite()->vertex() == &*ti );
	}
    }

    void io() 
    {
	/* Write map, read it back, and check that the maps are the
	 * same.
	 */
	string fn("t_SurfaceMap.out");

	SMap map(s,t,cmesh);
	{
	    std::ofstream os(fn.c_str());
	    MNI::write_surface_map( os, map );
	}

	ControlMesh cmesh2(3,6,2);
	cmesh2.vertices_push_back( ControlMesh::Vertex() );
	cmesh2.vertices_push_back( ControlMesh::Vertex() );
	cmesh2.vertices_push_back( ControlMesh::Vertex() );

	SMap newmap(s,t,cmesh2,fn);

	SMap::ControlMesh::Vertex_const_iterator 
	    v1 = map.control_mesh().vertices_begin();
	SMap::ControlMesh::Vertex_const_iterator 
	    v2 = newmap.control_mesh().vertices_begin();

	for( ; v1 != map.control_mesh().vertices_end(); ++v1,++v2 ) {
	    CPPUNIT_ASSERT( map.get_source(v1) == newmap.get_source(v2) );
	    CPPUNIT_ASSERT_EQUAL( map.get_target(v1), newmap.get_target(v2) );
        }
    }

    CPPUNIT_TEST_SUITE( t_SurfaceMap2 );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( io );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_SurfaceMap2 );
