#ifndef SURFACEMAP_H    // -*- C++ -*-
#define SURFACEMAP_H

// --------------------------------------------------------------//
// This header is deprecated.  Use SurfaceMap_hash.hpp instead.  //
// --------------------------------------------------------------//

#if MNI_NO_DEPRECATED
#error Header SurfaceMap.h is deprecated.
#endif


#include <surflib/SurfaceMap_hash.hpp>


namespace MNI {


template <class Surface_>
class SurfaceMap : public MNI::SurfaceMap_hash<Surface_>
{
public:

    typedef Surface_ Surface;

    typedef typename Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef typename Surface::Vertex_const_handle    Vertex_const_handle;
    typedef typename Surface::Vertex_const_iterator  Vertex_const_iterator;

    typedef Point_s<Surface>                         Point_s;


    SurfaceMap( const Surface& source_surf, const Surface& target_surf )
	: MNI::SurfaceMap_hash<Surface_>( source_surf, target_surf )
    {}

    SurfaceMap( const Surface& source_surf,
		const Surface& target_surf,
		const std::string& filename )
	: MNI::SurfaceMap_hash<Surface_>( source_surf, 
					  target_surf,
					  filename)
    {}


    void write( const std::string& filename ) const;

    /*! Template parameter VertexIndexMap must map a vertex handle
     *  to the index of the vertex in the list of the appropriate polyhedral
     *  surface.
     */
    template <class ControlIndexMap, class TargetIndexMap>
    void write( const std::string& filename, 
		const ControlIndexMap& control_index,
		const TargetIndexMap& target_index ) const;

    template <class ControlIndexMap, class TargetIndexMap>
    void write( std::ostream& os,
		const ControlIndexMap& control_index,
		const TargetIndexMap& target_index ) const;


    
    // Assignment operator: copy the hash map only
    SurfaceMap<Surface>& operator= ( const SurfaceMap<Surface>& s )
    {
	if ( this == &s )
	    return *this;

	CGAL_precondition( &source == &s.source );
	CGAL_precondition( &target == &s.target );

	//vertex_map = s.vertex_map;
	Vertex_const_iterator v = source.vertices_begin();
	CGAL_For_all( v,source.vertices_end() )
	    (*this)[v] = s[v];

	return *this;
    }
};

}


#include "SurfaceMap.cpp"

#endif
