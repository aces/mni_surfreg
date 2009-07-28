/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

// -*- C++ -*-

#include <surflib/FileReader.hpp>
#include <CGAL/Inverse_index.h>
#include <algorithm>
#include <fstream>
#include <iomanip>


namespace MNI {

namespace detail {

template <class _Surface>
void 
write_surface_sizes( std::ostream& os, const _Surface& s )
{
    // Cannot count on having facets available.
    // FIXME: do we bother with the halfedge counts?
#if SURFACE_MAP_DESIGN_ONE
    os << s.size_of_vertices() << " "
       << s.size_of_halfedges() << " "
       << s.size_of_facets() << std::endl;
#else
    os << s.size_of_vertices() << std::endl;
#endif
}

}


template <class OutputStream, class SurfaceMap_,
	  class ControlIndexMap, class TargetIndexMap>
void write_surface_map( OutputStream& os,
			const SurfaceMap_& sm,
			const ControlIndexMap& control_index,
			const TargetIndexMap& target_index )
{
    typedef typename SurfaceMap_::ControlMesh          Surface;
    typedef typename Surface::Vertex_const_iterator    Vertex_const_iterator;
    typedef Point_s<typename SurfaceMap_::TargetMesh>  Point_s;

    // First line identifies file type and format number.
    os << "# SurfaceMap 0" << std::endl
       << "# Columns: control vertex, target vertex0, target vertex1, t1, t2"
       << std::endl;

    detail::write_surface_sizes( os, sm.control_mesh() );
    detail::write_surface_sizes( os, sm.target_mesh() );

    Vertex_const_iterator v = sm.control_mesh().vertices_begin(); 
    for( int i = 0; v != sm.control_mesh().vertices_end(); ++v,++i ) {
	Point_s loc = sm[v];
	os << std::setprecision(17)
	   << control_index[v] << " "
	   << target_index[loc.V0()] << " "
	   << target_index[loc.V1()] << " "
	   << loc.t1() << " "
	   << loc.t2() << std::endl;
    }
}


template <class OutputStream, class SurfaceMap_>
void write_surface_map( OutputStream& os, 
			const SurfaceMap_& sm )
{
    typedef typename SurfaceMap_::ControlMesh::Vertex_const_iterator  Cvi;
    typedef typename SurfaceMap_::TargetMesh::Vertex_const_iterator   Tvi;

    CGAL::Inverse_index<Cvi> control_index( sm.control_mesh().vertices_begin(),
					    sm.control_mesh().vertices_end() );
    
    CGAL::Inverse_index<Tvi> target_index( sm.target_mesh().vertices_begin(),
					   sm.target_mesh().vertices_end() );
    
    write_surface_map( os, sm, control_index, target_index );
}




namespace detail {

template <class HDS>
const typename HDS::Vertex_const_handle*
index_to_handle_map( const HDS& s )
{
    typename HDS::Vertex_const_handle* map
	= new typename HDS::Vertex_const_handle[s.size_of_vertices()];

    typename HDS::Vertex_const_iterator v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++v,++i )
	map[i] = typename HDS::Vertex_const_handle(&*v);
    
    return map;
}

template <class _Surface>
const typename _Surface::Vertex** 
index_to_handle_map_v2( const _Surface& s )
{
    const typename _Surface::Vertex** map
	= new const typename _Surface::Vertex*[s.size_of_vertices()];

    typename _Surface::Vertex_const_iterator v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++v,++i )
	map[i] = &*v;
    
    return map;
}

}


template <class SurfaceMap_>
void read_surface_map( std::istream& is, 
		       SurfaceMap_& sm )
{
    read_surface_map( is, sm, sm.control_mesh().size_of_vertices() );
}


template <class SurfaceMap_>
void read_surface_map( std::istream& is, 
		       SurfaceMap_& sm,
		       typename SurfaceMap_::ControlMesh::size_type n_control )
{
    CGAL_precondition( n_control <= sm.control_mesh().size_of_vertices() );

    typedef typename SurfaceMap_::ControlMesh::Vertex_const_handle    Cv;
    typedef typename SurfaceMap_::TargetMesh::Vertex_const_handle     Tv;
    typedef typename SurfaceMap_::TargetMesh::Halfedge_const_handle   Th;

    typedef Point_s<typename SurfaceMap_::TargetMesh>                 Tpoint;

    FileReader fr(is);

    // FIXME: should allow for arbitrary # of comment lines
    //        check the file format version
    fr.expect_line();
    fr.expect_line();

    // FIXME: throw std::runtime_error("message") instead of assertions

    CGAL_assertion( fr.expect_uint() >= n_control );

    // FIXME: Handle old maps by ignoring optional 2 integers?
    //CGAL_assertion( fr.expect_uint() == sm.control_mesh().size_of_halfedges() );
    //CGAL_assertion( fr.expect_uint() == sm.control_mesh().size_of_facets() );

    CGAL_assertion( fr.expect_uint() == sm.target_mesh().size_of_vertices() );
    //CGAL_assertion( fr.expect_uint() == sm.target_mesh().size_of_halfedges() );
    //CGAL_assertion( fr.expect_uint() == sm.target_mesh().size_of_facets() );

    const typename SurfaceMap_::ControlMesh::Vertex_const_handle*
	control_vertex = detail::index_to_handle_map(sm.control_mesh());
    const typename SurfaceMap_::TargetMesh::Vertex_const_handle*
	target_vertex = detail::index_to_handle_map(sm.target_mesh());

    for( unsigned int i = 0; i < n_control; ++i ) {
	Cv vc( control_vertex[fr.expect_int()] );

	// Find halfedge in target mesh
	Tv v0( target_vertex[fr.expect_int()] );
	Tv v1( target_vertex[fr.expect_int()] );
	
	Th h = v1->halfedge();
	while (h->opposite()->vertex() != v0) {
	    h = h->next_on_vertex();
	    // throw runtime_error here
	    CGAL_assertion( h != v1->halfedge() );
	}

	double t1 = fr.expect_double();
	double t2 = fr.expect_double();

	sm[vc] = Tpoint( h,t1,t2 );
    }

    delete[] control_vertex;
    delete[] target_vertex;
}

}
