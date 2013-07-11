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

// #define CLAUDE_DEBUG

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

    unsigned int smap_control_size = fr.expect_uint();
    unsigned int smap_target_size = fr.expect_uint();

    std::cout << "Number of control vertices in input smap file = " 
              << smap_control_size << std::endl;
    std::cout << "Number of vertices in control mesh = "
         << sm.control_mesh().size_of_vertices() << std::endl;

    std::cout << "Number of target vertices in input smap file = " 
              << smap_target_size << std::endl;
    std::cout << "Number of vertices in target mesh = "
         << sm.target_mesh().size_of_vertices() << std::endl;

    if( smap_control_size == sm.control_mesh().size_of_vertices() ) {
      std::cout << "  No refinement in control size of smap." << std::endl;
    } else if( smap_control_size == 4*sm.control_mesh().size_of_vertices()-6 ) {
      std::cout << "  control of smap will be coarsened by 1 level" << std::endl;
      smap_control_size = sm.control_mesh().size_of_vertices();
    } else if( 4*smap_control_size-6 == sm.control_mesh().size_of_vertices() ) {
      std::cout << "  control of smap will be refined by 1 level" << std::endl;
    } else {
      throw std::runtime_error( "Incompatible control mesh sizes." );
      //CGAL_assertion( smap_control_size == sm.control_mesh().size_of_vertices() );
    }

    if( smap_target_size == sm.target_mesh().size_of_vertices() ) {
      std::cout << "  No refinement in target size of smap." << std::endl;
    } else if( smap_target_size == 4*sm.target_mesh().size_of_vertices()-6 ) {
      throw std::runtime_error( "Coarsening of target mesh is not supported." );
    } else if( 4*smap_target_size-6 == sm.target_mesh().size_of_vertices() ) {
      std::cout << "  target of smap will be refined by 1 level" << std::endl;
    } else {
      throw std::runtime_error( "Incompatible target mesh sizes." );
      // CGAL_assertion( smap_target_size == sm.target_mesh().size_of_vertices() );
    }

    const typename SurfaceMap_::ControlMesh::Vertex_const_handle*
	control_vertex = detail::index_to_handle_map(sm.control_mesh());
    const typename SurfaceMap_::TargetMesh::Vertex_const_handle*
	target_vertex = detail::index_to_handle_map(sm.target_mesh());

#if CLAUDE_DEBUG
    typedef typename SurfaceMap_::ControlMesh::Vertex_const_iterator  Cvi;
    typedef typename SurfaceMap_::TargetMesh::Vertex_const_iterator   Tvi;
    CGAL::Inverse_index<Cvi> control_index( sm.control_mesh().vertices_begin(),
					    sm.control_mesh().vertices_end() );
    CGAL::Inverse_index<Tvi> target_index( sm.target_mesh().vertices_begin(),
					   sm.target_mesh().vertices_end() );

    Vertex_const_iterator v = sm.control_mesh().vertices_begin(); 
    for( int i = 0; v != sm.control_mesh().vertices_end(); ++v,++i ) {
      Point_s loc = sm[v];
      os << std::setprecision(17) << control_index[v] << " "
         << target_index[loc.V0()] << " " << target_index[loc.V1()] << " "
         << loc.t1() << " " << loc.t2() << std::endl;
   }
#endif

    for( unsigned int i = 0; i < smap_control_size; ++i ) {
      Cv vc( control_vertex[fr.expect_int()] );

      // Find halfedge in target mesh
      Tv v0( target_vertex[fr.expect_int()] );
      Tv v1( target_vertex[fr.expect_int()] );
      double t1 = fr.expect_double();
      double t2 = fr.expect_double();
      double t0 = 1.0 - t1 - t2;

      if( smap_target_size == sm.target_mesh().size_of_vertices() ) {

	Th h = v1->halfedge();
	while (h->opposite()->vertex() != v0) {
	    h = h->next_on_vertex();
	    // throw runtime_error here
	    CGAL_assertion( h != v1->halfedge() );
	}
	sm[vc] = Tpoint( h,t1,t2 );

      } else if( 4*smap_target_size-6 == sm.target_mesh().size_of_vertices() ) {

        // Find edge in coarse-to-fine mesh. A bit messy.
        //
        //                      v0
        //                     /  \
        //                    /    \
        //                   /      \
        //                  /        \
        //             op2 x----------x op1
        //                / \        / \
        //               /   \      /   \
        //              /     \    /     \
        //             /       \  /       \
        //            /         \/         \
        //           v1----------x---------v2
        //                      op0
        // Point is given on coarse target mesh with edge (v0:v1),
        // linked to v2, with local coords (1-xi-eta,xi,eta). We
        // need to find the (x,y,z) point in the big triangle (v0:v1:v2)
        // and see in which sub-triangle it lies.
        //

	Th h = v1->halfedge();
        Tv v2, op0, op1, op2;

        int done = 0;
	while( !done ) {
          op2 = h->opposite()->vertex();
          Th h2 = op2->halfedge();
	  while( h2->opposite()->vertex() != v0 ) {
	      h2 = h2->next_on_vertex();
              if( h2 == op2->halfedge() ) break;
	  }
          if( h2->opposite()->vertex() == v0 ) {
            op1 = h2->next_on_vertex()->opposite()->vertex();
	    h2 = h2->next_on_vertex();
            op0 = h2->next_on_vertex()->opposite()->vertex();
            // Finally, find v2 from op0 with opposite vertex op1.
	    h2 = op0->halfedge();
	    while (h2->opposite()->vertex() != op1) {
	        h2 = h2->next_on_vertex();
	    }
            v2 = h2->next_on_vertex()->opposite()->vertex();

            // Find in which sub-triangle the point is lying.
            if( t0 > 0.5 ) {
              // use triangle v0:op2:op1
              h = op2->halfedge();
	      while( h->opposite()->vertex() != v0 ) h = h->next_on_vertex();
              t1 *= 2.0;
              t2 *= 2.0;
            } else if( t1 > 0.5 ) {
              // use triangle v1:op0:op2
              h = op0->halfedge();
	      while( h->opposite()->vertex() != v1 ) h = h->next_on_vertex();
              t1 = 2.0*t2;
              t2 = 2.0*t0;
            } else if( t2 > 0.5 ) {
              // use triangle v2:op1:op0
              h = op1->halfedge();
	      while( h->opposite()->vertex() != v2 ) h = h->next_on_vertex();
              t2 = 2.0*t1;
              t1 = 2.0*t0;
            } else {
              // use triangle op0:op1:op2
              h = op1->halfedge();
	      while( h->opposite()->vertex() != op0 ) h = h->next_on_vertex();
              t1 = 1.0 - 2.0*t1;
              t2 = 1.0 - 2.0*t2;
            }
	    sm[vc] = Tpoint( h, t1, t2 );
            done = 1;
          }
          if( !done ) {
	    h = h->next_on_vertex();
            if( h == v1->halfedge() ) break;
          }
        }
      } else {
        // We cannot allow coarsening. Code should have died by now.
	sm[vc] = Tpoint();
      }
    }

    delete[] control_vertex;
    delete[] target_vertex;
}

}
