#include "tester.hpp"

#include "Surface.hpp"
#include <surflib/SurfaceMap_hash.hpp>

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



class t_SurfaceMap : public CppUnit::TestCase
{
public:

    typedef Surface_mpfloat                    Surface;
    typedef Surface::Halfedge_handle        Halfedge_handle;
    typedef Surface::Facet_handle           Facet_handle;
    typedef Surface::Vertex_const_handle    Vertex_const_handle;
    typedef Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Surface::Point_3                Point;
    typedef MNI::Point_s<Surface>           Point_s;


    // Points with which to construct surfaces.
    Point a,b,c;
    Surface s,t;
    
    void setUp()
    {
	a = Point(0,0,0);
	b = Point(1,0,0);
	c = Point(0,1,0);

	s.make_triangle(a,b,c);
	t.make_triangle(a,b,c);
    }


    void constructors() 
    {
	MNI::SurfaceMap_hash<Surface> map(s,t);

	// Check that map is identity
	Vertex_const_iterator vs = s.vertices_begin();
	Vertex_const_iterator vt = t.vertices_begin();

	for( ; vs != s.vertices_end(); ++vs,++vt ) {
	    Point_s loc = map[vs];
	    CPPUNIT_ASSERT( loc.t1() == 0 );
	    CPPUNIT_ASSERT( loc.t2() == 0 );
	    CPPUNIT_ASSERT( &*loc.halfedge()->opposite()->vertex() == &*vt );
	}
    }

    void io() 
    {
	// Write map, read it back.

	MNI::SurfaceMap_hash<Surface> map(s,t);
	std::ofstream os("t_SurfaceMap.out");
	MNI::write_surface_map( os, map );
	
	MNI::SurfaceMap_hash<Surface> newmap(s,t,"t_SurfaceMap.out");

	Vertex_const_iterator v = s.vertices_begin(); 
	for( ; v != s.vertices_end(); ++v ) {
	    CPPUNIT_ASSERT_EQUAL( map[v], newmap[v] );
        }
    }

    CPPUNIT_TEST_SUITE( t_SurfaceMap );
    CPPUNIT_TEST( constructors );
    CPPUNIT_TEST( io );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_SurfaceMap );
