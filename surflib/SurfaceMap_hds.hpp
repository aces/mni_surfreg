#ifndef SURFLIB_SURFACEMAP_HDS_HPP
#define SURFLIB_SURFACEMAP_HDS_HPP

#include <CGAL/HalfedgeDS_min_items.h>
#include <CGAL/HalfedgeDS_vector.h>

#include <CGAL/circulator.h>

#include <surflib/Point_s.hpp>
#include <surflib/SurfaceMapBase.hpp>
#include <surflib/SurfaceMapIO.hpp>


namespace MNI {

namespace detail {

    /* In order to support vertices at all, the option
     * "Supports_halfedge_vertex" must be Tag_true.  In addition,
     * we need "Supports_vertex_halfedge" to enable the halfedge pointer
     * in each vertex.
     */

    /* The vertex of the control mesh needs the pointer to an incident
     * halfedge, so we use the HalfedgeDS_vertex_base class template with
     * the second template parameter set to Tag_true.
     *
     * We also want a handle for the corresponding source mesh vertex,
     * and a Point_s on the target mesh.  These classes are passed in
     * the Traits structure.
     *
     * Model for CGAL::HalfedgeDSVertex.
     */
    template <class Refs, class Traits>
    struct HalfedgeDS_control_vertex 
	: public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true>
    {
	typedef typename Traits::Source_handle Source_handle;
	typedef typename Traits::Target_point  Target_point;

	Source_handle source_handle;
	Target_point target_point;
    };

    // Model for CGAL::HalfedgeDSItems
    //
    struct HalfedgeDS_control_items 
	: public CGAL::HalfedgeDS_min_items
    {
	template <class Refs, class Traits>
	struct Vertex_wrapper {
	    typedef HalfedgeDS_control_vertex<Refs,Traits> Vertex;
	};

	template <class Refs, class Traits>
	struct Halfedge_wrapper {
	    typedef CGAL::HalfedgeDS_halfedge_base<Refs,
						   CGAL::Tag_false,
						   CGAL::Tag_true,
						   CGAL::Tag_false> Halfedge;
	};
    };

    template <class SourceMesh_,
	      class TargetMesh_>
    struct HalfedgeDS_control_traits
    {
	typedef typename SourceMesh_::Vertex_const_handle   Source_handle;
	typedef Point_s<TargetMesh_>                        Target_point;
    };
}


/*! \brief SurfaceMap implemented using halfedge data structure
 *  for the control mesh.
 *
 */

template <class SourceMesh_,
	  class TargetMesh_ = SourceMesh_>
class SurfaceMap_hds : public MNI::SurfaceMapBase<SourceMesh_,TargetMesh_>
{
public:
    typedef SourceMesh_ SourceMesh;
    typedef TargetMesh_ TargetMesh;

    typedef CGAL::HalfedgeDS_vector
    <detail::HalfedgeDS_control_traits<SourceMesh,TargetMesh>,
     detail::HalfedgeDS_control_items>  ControlMesh;

    typedef typename SourceMesh::Vertex_const_handle    Source_const_handle;
    typedef Point_s<TargetMesh>                         Target_point;

    typedef typename ControlMesh::Vertex_const_handle   Control_const_handle;
    typedef typename ControlMesh::Vertex_handle         Control_handle;


    /*!
     *  Initialize with mapping source vertex 0 --> target vertex 0, 
     *  vertex 1 --> vertex 1, etc.
     */
    SurfaceMap_hds( const SourceMesh& source_surf,
		    const TargetMesh& target_surf,
		    ControlMesh& control_mesh );

    /*!  Read mapping from filename.
     */
    SurfaceMap_hds( const SourceMesh& source_surf,
		    const TargetMesh& target_surf,
		    ControlMesh& control_mesh,
		    const std::string& filename );

    
    /*! Map control vertex to source and target locations.
     */
    Source_const_handle get_source( Control_const_handle v ) const
    { return v->source_handle; }

    Target_point& get_target( Control_const_handle v )
    {  return const_cast<Target_point&>(v->target_point);  }

    const Target_point& get_target( Control_const_handle v ) const
    {  return v->target_point;  }


    //! Convenience notation for get_target().
    Target_point& operator[] ( Control_const_handle v )
    {  return get_target(v);  }

    const Target_point& operator[] ( Control_const_handle v ) const
    {  return get_target(v);  }


    /*! Access to the polyhedral surface structures.
     */
    const ControlMesh& control_mesh() const
    { return control; }

    ControlMesh& control_mesh()
    { return control; }


private:
    //! Initialize to identical-order mapping.
    void initialize_as_identity();


private:
    ControlMesh& control;
};


}

#include "SurfaceMap_hds.cpp"
#endif

