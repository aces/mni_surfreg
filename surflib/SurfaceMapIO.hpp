/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef SURFLIB_SURFACEMAPIO_HPP
#define SURFLIB_SURFACEMAPIO_HPP

/*! \defgroup surfacemap_io SurfaceMap I/O
 *
 * I/O Routines for objects that implement the SurfaceMap concept.
 *
 * \@{
 */

#include <istream>


namespace MNI {

/*! \brief Read surface map from input stream.
 *
 * This function replaces an existing surface mapping, keeping the
 * source, target, and control meshes.  
 *
 * \attention The mapping in the stream must be applicable to the
 * given meshes.  In particular, the stream must have the same control
 * mesh topology as the current surface map, and the same goes for the
 * target mesh.  The stream format of a surface map codes vertices as
 * small integer values.  This function assumes that the control
 * vertex numbers in the stream correspond to the current ordering of
 * the vertices in the surface map control mesh; i.e., control vertex
 * 0 in the stream is the first control vertex in the current surface
 * map, etc.  The same goes for the numbering of target vertices.
 *
 * \bug The stream format of a surface map is subject to change.  
 */
template <class SurfaceMap_>
void read_surface_map( std::istream& is, 
		       SurfaceMap_& sm );


/*! \brief Read partial surface map from input stream.
 *
 * This function overwrites an existing surface mapping, keeping the
 * source, target, and control meshes.  Only the first \a n_control
 * control vertices are updated from the file.
 *
 * \pre \a n_control must be less than or equal to the number of vertices
 *      in the existing control mesh.
 */
template <class SurfaceMap_>
void read_surface_map( std::istream& is, 
		       SurfaceMap_& sm,
		       typename SurfaceMap_::ControlMesh::size_type n_control );




/*! \brief Write surface map to output stream.
 *
 * Outputs the association between control vertices and locations on
 * the target mesh.  Vertices of the meshes are numbered with an
 * integer index value when written to the stream.  Vertices are
 * numbered according to the order in the control and target mesh
 * vertex lists.
 *
 * \bug The stream format of a surface map is subject to change.  
 */
template <class OutputStream, class SurfaceMap_>
void write_surface_map( OutputStream& os, 
			const SurfaceMap_& sm );


/*! \brief Write surface map to output stream with vertex remapping.
 *
 * If the vertex list of the control or target mesh have been reordered,
 * this function can be used to output the map with the correct indices.
 *
 * \arg control_index maps each control mesh vertex handle
 *      to an integer
 * \arg target_index maps each target mesh vertex handle
 *      to an integer
 *
 * The \c ControlIndexMap and \c TargetIndexMap classes must 
 * contain an operator[] that maps a vertex handle to an integer.
 *
 * \bug The stream format of a surface map is subject to change.  
 */
template <class OutputStream, class SurfaceMap_,
	  class ControlIndexMap, class TargetIndexMap>
void write_surface_map( OutputStream& os,
			const SurfaceMap_& sm,
			const ControlIndexMap& control_index,
			const TargetIndexMap& target_index );

}


/* \@} */

// Implementation
#include "SurfaceMapIO.cpp"

#endif
