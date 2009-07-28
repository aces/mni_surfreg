/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "tester.hpp"

#include <surflib/Object_file_inserter.hpp>

using namespace std;

struct Traits { typedef int Point_2; };
typedef CGAL_HALFEDGEDS_DEFAULT<Traits> HDS;



class t_Object_file_inserter : public CppUnit::TestCase
{
public:

    typedef HDS::Vertex_iterator        Vertex_iterator;
    typedef HDS::Halfedge_handle        Halfedge_handle;


    HDS mesh;

    
    void check_vertex_pointers() 
    {
	MNI::Object_file infile( "example1.obj" );
	CPPUNIT_ASSERT( infile.is_valid() );
	CPPUNIT_ASSERT_EQUAL( 1, infile.size_of_objects() );
	MNI::ObjFile_inserter::insert_into_hds( infile, 0, true, mesh );

	Vertex_iterator v = mesh.vertices_begin();
	CGAL_For_all( v, mesh.vertices_end() ) {
	    // Cannot assume existence of
	    // Halfedge_around_vertex_circulator, as control mesh may
	    // only be HDS.  Do the loop the old-fashioned way.
	    //
	    Halfedge_handle h = v->halfedge();
	    do {
		CPPUNIT_ASSERT( v == h->vertex() );
		h = h->next()->opposite();
	    } while ( h != v->halfedge() );
	}
    }


    CPPUNIT_TEST_SUITE( t_Object_file_inserter );
    CPPUNIT_TEST( check_vertex_pointers );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_Object_file_inserter );
