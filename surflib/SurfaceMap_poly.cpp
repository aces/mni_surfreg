
template <class SourceMesh_,class TargetMesh_>
MNI::SurfaceMap_poly<SourceMesh_,TargetMesh_>
::SurfaceMap_poly( const SourceMesh& source_surf,
		   const TargetMesh& target_surf,
		   ControlMesh& control_mesh )
    : MNI::SurfaceMapBase<SourceMesh,TargetMesh>( source_surf, 
						  target_surf ),
      control(control_mesh)
{
    initialize_as_identity();
}


template <class SourceMesh_,class TargetMesh_>
void MNI::SurfaceMap_poly<SourceMesh_,TargetMesh_>::initialize_as_identity()
{
    CGAL_precondition( control.size_of_vertices() > 0 );
    CGAL_precondition( control.size_of_vertices() <= source.size_of_vertices() );
    CGAL_precondition( target.size_of_vertices() > 0 );
    
    typename ControlMesh::Vertex_iterator       vc = control.vertices_begin();
    typename SourceMesh::Vertex_const_iterator 	vs = source.vertices_begin();
    CGAL::Circulator_from_iterator<typename TargetMesh::Vertex_const_iterator> 
	vt( target.vertices_begin(), target.vertices_end() );
    
    CGAL_For_all( vc, control.vertices_end() ) {
	vc->source_handle = vs;
	typename TargetMesh::Vertex_const_handle vt_handle(&*vt);
	vc->point() = Target_point(vt_handle);
	++vs;
	++vt;
    }
}


/*!  Read mapping from filename.
 */
template <class SourceMesh_,class TargetMesh_>
MNI::SurfaceMap_poly<SourceMesh_,TargetMesh_>
::SurfaceMap_poly( const SourceMesh& source_surf,
		   const TargetMesh& target_surf,
		   ControlMesh& control_mesh,
		   const std::string& filename )
    : MNI::SurfaceMapBase<SourceMesh,TargetMesh>( source_surf, 
						  target_surf ),
      control(control_mesh)
{
    initialize_as_identity();
    std::ifstream is(filename.c_str());
    MNI::read_surface_map( is, *this );
}

