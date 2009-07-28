/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFLIB_SURFACEMAPBASE_HPP
#define SURFLIB_SURFACEMAPBASE_HPP

namespace MNI {

/*! \brief Surface to surface mapping.
 *
 * The mapping associates locations on the source surface to locations
 * on the target surface.  The mapping is stored at a subset of the
 * source vertices, termed the control grid.  The control grid graph
 * topology may differ from that of the source surface, even if there
 * is one control grid for each source grid.
 *
 * The mapping is stored by giving a Point_s located on the target
 * surface, at each vertex of the control grid.  Interpolation is left
 * to the user.
 *
 * A class that models this concept must provide the following
 * types and methods.
 * - class SourceMesh
 * - class TargetMesh
 * - class ControlMesh
 * - Vertex_handle get_source( Vertex_handle )
 * - Point_s get_target( Vertex_handle )
 * - Point_s operator[]( Vertex_handle ), alias for previous
 * - const Surface& source_mesh()
 * - const Surface& target_mesh()
 * - const Surface& control_mesh()
 */


template <class SourceMesh,
	  class TargetMesh = SourceMesh>
class SurfaceMapBase
{
protected:
    //! Constructs undefined surface mapping.
    SurfaceMapBase( const SourceMesh& s,
		    const TargetMesh& t )
	: source(s), target(t)
    {}


private:
    // Disable the copy constructor and assignment operator.
    SurfaceMapBase( const SurfaceMapBase<SourceMesh,TargetMesh>& );
    const SurfaceMapBase& operator=( const SurfaceMapBase<SourceMesh,TargetMesh>& );


public:
    /* Subclass must implement
     * - Vertex_handle get_source( Vertex_handle )
     * - Point_s get_target( Vertex_handle )
     * - Point_s operator[]( Vertex_handle ), alias for previous
     */

    /*! Access to the polyhedral surface structures.
     */
    const SourceMesh& source_mesh() const
    { return source; }

    const TargetMesh& target_mesh() const
    { return target; }

    /* Subclass must implement
     * - const Surface& control_mesh()
     */


protected:
    const SourceMesh& source;
    const TargetMesh& target;
};

}

#endif
