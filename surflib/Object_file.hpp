/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef MNI_OBJFILE_H         // -*- C++ -*-
#define MNI_OBJFILE_H

#include <string>
#include <CGAL/Polyhedron_3.h>


// Defined in <bicpl.h>
struct object_struct;


namespace MNI {


/*! \brief MNI/BIC .obj file of objects.
 *
 * A BIC .obj file contains one or more objects.  The objects are
 * indexed from 0 to size_of_objects() - 1.
 * 
 * Objects may be a mesh of polygons, lines, and other types.  Only
 * polygons and quadmesh types are currently supported.
 */
class Object_file
{
public:

    Object_file( char* filename );
    ~Object_file();


    // Accessor methods

    //! True if file correctly read
    bool is_valid() const;
    int size_of_objects() const;
    object_struct* object( int index ) const;

    bool is_polygons( int index ) const;
    bool is_quadmesh( int index ) const;

private:
    int objects_size;
    object_struct** objects_list;
};


}

#endif
