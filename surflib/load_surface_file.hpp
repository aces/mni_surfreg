#ifndef LOAD_SURFACE_FILE_HPP
#define LOAD_SURFACE_FILE_HPP

#include <cstring>
#include <surflib/Object_file_inserter.hpp>
#include <surflib/VTKFile.hpp>

extern "C" {
#include <bicpl/object_prototypes.h>
}


namespace MNI {


/*! \brief Load a VTK disk file into a surface structure.
 *
 *  Only the geometry is loaded.
 */
template <class Polyhedron>
void load_vtk_file( Polyhedron& s, const char* filename )
{
    VTKFile vtk_in( filename );
    vtk_in.read_geometry(s);
}


/*! \brief Load a VTK disk file into a surface structure.
 *
 *  Geometry and the first scalar is loaded.  The polyhedron
 *  class must have a data member named "scalar".
 */
template <class Polyhedron>
void load_vtk_file_with_scalar( Polyhedron& s, const char* filename )
{
    VTKFile vtk_in( filename );
    vtk_in.read_geometry(s);
    SetScalarValues<Polyhedron> ssv(s);
    vtk_in.read_scalar_point_attribute(ssv);
}


/*! \brief Load a BIC OBJ disk file into a surface structure.
 *
 *  Only the geometry is loaded.
 */
template <class Polyhedron>
void load_bicobj_file( Polyhedron& s, char* filename )
{
    MNI::Object_file infile( filename );
    CGAL_assertion( infile.is_valid() );
    CGAL_assertion( infile.size_of_objects() == 1 );
    MNI::ObjFile_inserter::insert( infile, 0, true, s );
}


/*! \brief Load a BIC OBJ disk file into a surface structure.
 *
 *  The geometry and a scalar field are loaded.  The polyhedron
 *  class must have a data member named "scalar".
 */
template <class Polyhedron>
void load_bicobj_file_with_scalar( Polyhedron& s, 
				   char* obj_filename,
				   char* vv_filename )
{
    load_bicobj_file( s, obj_filename );

    int n;
    double* values;
    if (input_texture_values(vv_filename, &n, &values) != OK )
	throw std::runtime_error( "load_bicobj_file_with_scalar: cannot read vertex value file" );
    if ( static_cast<typename Polyhedron::size_type>(n) != s.size_of_vertices() )
	throw std::runtime_error( "load_bicobj_file_with_scalar: number of values in file != number of vertices" );
	
    typename Polyhedron::Vertex_iterator v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++i,++v ) 
	v->scalar = values[i];

    FREE( values );
}


/*! \brief Load a disk file into a surface structure.
 *
 *  Currently handles BIC OBJ files and VTK files.
 *  Only the geometry is loaded.
 */
template <class Polyhedron>
void load_surface_file( Polyhedron& s, char* filename )
{
    // Guess file format.
    bool is_vtk = true;
    int len = strlen(filename);
    if ( len > 4 && strcmp(&filename[len-4], ".obj") == 0 )	
	is_vtk = false;
    
    if ( is_vtk ) 
	load_vtk_file( s, filename );
    else
	load_bicobj_file( s, filename );
}


/*! \brief Load geometry and scalar from disk file(s) into a surface structure.
 *
 *  Currently handles BIC OBJ files and VTK files.  If the geometry
 *  filename ends in ".obj", the attribute name is interpreted as a
 *  vertex value filename.  Otherwise, it is interpreted as the
 *  attribute name in the VTK file.
 */
template <class Polyhedron>
void load_surface_with_scalar( Polyhedron& s,
			       char* surface_filename,
			       char* attribute_name )
{
    // Guess file format.
    bool is_vtk = true;
    int len = strlen(surface_filename);
    if ( len > 4 && strcmp(&surface_filename[len-4], ".obj") == 0 )	
	is_vtk = false;
    
    if ( is_vtk ) {
	VTKFile vtk_in( surface_filename );
	vtk_in.read_geometry(s);
	SetScalarValues<Polyhedron> ssv(s);
	vtk_in.read_scalar_point_attribute( attribute_name, ssv );
    } else {
	load_bicobj_file_with_scalar( s, surface_filename, attribute_name );
    }
}



}


#endif
