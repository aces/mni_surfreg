/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFLIB_SURFACEMAP_HASH_HPP
#define SURFLIB_SURFACEMAP_HASH_HPP

#include <CGAL/circulator.h>
#include <CGAL/Unique_hash_map.h>
#include <surflib/Point_s.hpp>
#include <surflib/SurfaceMapBase.hpp>
#include <surflib/SurfaceMapIO.hpp>


namespace MNI {


/*! \brief SurfaceMap implemented with hashing.
 *
 * The control mesh is isomorphic to the source mesh.
 * It is stored implicitly.
 */

template <class SourceMesh_,
	  class TargetMesh_ = SourceMesh_>
class SurfaceMap_hash : public MNI::SurfaceMapBase<SourceMesh_,TargetMesh_>
{
public:

    typedef SourceMesh_ SourceMesh;
    typedef TargetMesh_ TargetMesh;
    typedef SourceMesh_ ControlMesh;

    typedef typename SourceMesh::Vertex_const_handle    Vertex_const_handle;
    typedef Point_s<TargetMesh>                         Point_s;


    /*! Implicitly set control grid to be a copy of the source
     *  surface mesh.
     *
     *  Initialize with mapping source vertex 0 --> target vertex 0, 
     *  vertex 1 --> vertex 1, etc.
     */
    SurfaceMap_hash( const SourceMesh& source_surf,
		     const TargetMesh& target_surf )
	: MNI::SurfaceMapBase<SourceMesh,TargetMesh>( source_surf, 
						      target_surf )
    {
	CGAL_assertion( this->target.size_of_vertices() > 0 );
	typename SourceMesh::Vertex_const_iterator s = this->source.vertices_begin();
	CGAL::Circulator_from_iterator<typename TargetMesh::Vertex_const_iterator> 
	    t( this->target.vertices_begin(), this->target.vertices_end() );

	while( s != this->source.vertices_end() ) {
	    (*this)[s] = Point_s(&*t);
	    ++s;
	    ++t;
	}
    }


    /*! Implicitly set control grid to be a copy of the source
     *  surface mesh.
     *
     *  Read mapping from filename.
     */
    SurfaceMap_hash( const SourceMesh& source_surf,
		     const TargetMesh& target_surf,
		     const std::string& filename )
	: MNI::SurfaceMapBase<SourceMesh,TargetMesh>( source_surf, 
						      target_surf )
    {
	std::ifstream is(filename.c_str());
	MNI::read_surface_map( is, *this );
    }

    
    /*! Map control vertex to source and target locations.
     */
    Vertex_const_handle get_source( Vertex_const_handle v ) const
    { return v; }

    Point_s& get_target( Vertex_const_handle v )
    {  return vertex_map[v];  }

    const Point_s& get_target( Vertex_const_handle v ) const
    {  return vertex_map[v];  }


    //! Convenience notation for get_target().
    Point_s& operator[] ( Vertex_const_handle v )
    {  return get_target(v);  }

    const Point_s& operator[] ( Vertex_const_handle v ) const
    {  return get_target(v);  }


    /*! Access to the polyhedral surface structures.
     */
    const SourceMesh& control_mesh() const
    { return this->source; }


private:
    // For the source only, we store the transformation
    // info at each vertex.
    CGAL::Unique_hash_map<Vertex_const_handle,Point_s> vertex_map;
};


}

#endif

