/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFLIB_SURFACEMAP_POLYHEDRON_HPP
#define SURFLIB_SURFACEMAP_POLYHEDRON_HPP

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/circulator.h>

#include <surflib/Point_s.hpp>
#include <surflib/SurfaceMapBase.hpp>
#include <surflib/SurfaceMapIO.hpp>


namespace MNI {

namespace detail {

    /* The vertex stores the target location in the standard "point"
     * member, and the source location in an additional member.
     */

    /* Define a bogus "Point_3" to satisfy code that expects to
     * construct a point from 3 coordinate values.
     */
    template <class TargetMesh_>
    struct Control_Point_3 : public Point_s<TargetMesh_>
    {
	Control_Point_3( double, double, double ) {}
	Control_Point_3() {}
	
	// Bring in necessary Point_s constructors.
	//
	typedef Point_s<TargetMesh_>  Base;
	
	Control_Point_3( const Base& p ) : Base(p) {}

	Control_Point_3( typename TargetMesh_::Halfedge_const_handle h,
			 double t1, double t2 )
	    : Base(h,t1,t2) {}

	Control_Point_3( typename TargetMesh_::Vertex_const_handle v )
	    : Base(v) {}
    };


    template <class SourceMesh_,
	      class TargetMesh_>
    struct Control_PolyhedronTraits_3 
    {
	typedef typename SourceMesh_::Vertex_const_handle   Source_handle;
	typedef Control_Point_3<TargetMesh_>                Target_point;

	// Required types
	typedef Target_point Point_3;
	typedef int Plane_3;
    };


    template <class Refs, class Traits>
    struct Control_vertex
	: public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
					      typename Traits::Target_point>
    {
	typedef typename Traits::Target_point Point;
	typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point>
	Vertex_base;

	Control_vertex( const Point& p ) : Vertex_base(p) {}
	Control_vertex() {}

	typename Traits::Source_handle source_handle;
    };


    struct Control_Polyhedron_items_3 : public CGAL::Polyhedron_items_3
    {
	template <class Refs,class Traits>
	struct Vertex_wrapper {
	    typedef typename Traits::Target_point Point;
	    typedef Control_vertex<Refs,Traits> Vertex;
	};
    };
}


/*! \brief SurfaceMap implemented using polyhedron_3 data structure
 *  for the control mesh.
 */

template <class SourceMesh_,
	  class TargetMesh_ = SourceMesh_>
class SurfaceMap_poly : public MNI::SurfaceMapBase<SourceMesh_,TargetMesh_>
{
public:
    typedef SourceMesh_ SourceMesh;
    typedef TargetMesh_ TargetMesh;

    typedef CGAL::Polyhedron_3
    <detail::Control_PolyhedronTraits_3<SourceMesh,TargetMesh>,
     detail::Control_Polyhedron_items_3,
     CGAL::HalfedgeDS_vector>  ControlMesh;

    typedef typename SourceMesh::Vertex_const_handle    Source_const_handle;
    typedef typename ControlMesh::Traits::Target_point  Target_point;

    typedef typename ControlMesh::Vertex_const_handle   Control_const_handle;
    typedef typename ControlMesh::Vertex_handle         Control_handle;


    /*!
     *  Initialize with mapping source vertex 0 --> target vertex 0, 
     *  vertex 1 --> vertex 1, etc.
     */
    SurfaceMap_poly( const SourceMesh& source_surf,
		     const TargetMesh& target_surf,
		     ControlMesh& control_mesh );

    /*!  Read mapping from filename.
     */
    SurfaceMap_poly( const SourceMesh& source_surf,
		     const TargetMesh& target_surf,
		     ControlMesh& control_mesh,
		     const std::string& filename );

    
    /*! Map control vertex to source and target locations.
     */
    Source_const_handle get_source( Control_const_handle v ) const
    { return v->source_handle; }

    Target_point& get_target( Control_const_handle v )
    {  return const_cast<Target_point&>(v->point());  }

    const Target_point& get_target( Control_const_handle v ) const
    {  return v->point();  }


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

#include "SurfaceMap_poly.cpp"
#endif
