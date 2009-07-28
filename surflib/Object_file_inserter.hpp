/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef MNI_OBJECT_FILE_INSERTER_H // -*- C++ -*-
#define MNI_OBJECT_FILE_INSERTER_H

#include <stdexcept>
#include <surflib/Object_file.hpp>
#include <CGAL/Polyhedron_incremental_builder_3.h>


// GRR! Both volume_io and gmp insist on defining these two macros
#undef ABS
#undef ALLOC

extern "C" {
    // Define lint so avoid having all kinds of rcsid[] strings compiled
    // in.  This cuts way down on the warning messages.
#  define lint 1
#  include <bicpl.h>
#  undef lint
}


namespace MNI {
namespace ObjFile_inserter {



/*! \internal
 */
template <class HalfedgeDS>
class Polygons : public CGAL::Modifier_base<HalfedgeDS>
{
public:
    const polygons_struct& poly;

    Polygons( const polygons_struct& p ) : poly(p) {};
    void operator() ( HalfedgeDS& hds );
};


/*! \internal
 */
template <class HalfedgeDS>
class Quadmesh : public CGAL::Modifier_base<HalfedgeDS>
{
public:
    const quadmesh_struct& qmesh;
    bool split_patches;

    Quadmesh( const quadmesh_struct& q, bool split )
	: qmesh(q), split_patches(split) {};
    void operator() ( HalfedgeDS& hds );
};


/*!
 * \param index object to insert
 * \param split_quad if true, split quadmesh patches into 
 *                            two triangular facets
 * \param surf the surface data structure in which to insert the object
 */
template <class Polyhedron>
void insert( const Object_file& ofile, int index, bool split_quad, 
	     Polyhedron& surf )
{
    typedef typename Polyhedron::HalfedgeDS HalfedgeDS;

    object_struct* obj_ptr = ofile.object(index);

    if ( ofile.is_polygons(index) ) {
	MNI::ObjFile_inserter::Polygons<HalfedgeDS>
	    modifier( *(get_polygons_ptr(obj_ptr)) );
	surf.delegate( modifier );
    } else if ( ofile.is_quadmesh(index) ) {
	MNI::ObjFile_inserter::Quadmesh<HalfedgeDS>
	    modifier ( *(get_quadmesh_ptr(obj_ptr)), split_quad );
	surf.delegate( modifier );
    } else {
	CGAL_assertion_msg( false, "Unable to handle object type" );
	return;
    }
}



//! Insert edges into hds.  No facet nor point information is used.
template <class HDS>
void insert_polygons_into_hds( const polygons_struct& poly,
			       HDS& hds );


/*!
 * \param index object to insert
 * \param split_quad if true, split quadmesh patches into 
 *                            two triangular facets
 * \param hds the halfedge data structure in which to insert the object
 */
template <class HDS>
void insert_into_hds( const Object_file& ofile, int index, bool split_quad, 
		      HDS& hds )
{
    object_struct* obj_ptr = ofile.object(index);

    if ( ofile.is_polygons(index) ) {
	insert_polygons_into_hds( *(get_polygons_ptr(obj_ptr)), hds );
    } else {
	throw std::runtime_error( "insert_into_hds(): cannot handle object type" );
    }
}


}
}


#include "Polygons_inserter.cpp"
#include "Quadmesh_inserter.cpp"

#endif
