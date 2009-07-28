/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <ostream>
#include <stdexcept>

#include <surflib/Object_file.hpp>
#include <surflib/Object_file_inserter.hpp>

namespace MNI {



/*! This constructor loads the entire file into memory.
 *
 * \post is_valid() returns false if file was not read correctly
 */
Object_file::Object_file( char* filename )
{
    File_formats format;
    
    if ( input_graphics_file( filename, &format, 
			      &objects_size, &objects_list ) != OK )

	throw std::runtime_error( "cannot read graphics file");
}


Object_file::~Object_file()
{
    delete_object_list( objects_size, objects_list );
}


bool Object_file::is_valid() const
{
    return objects_size > 0;
}



int Object_file::size_of_objects() const
{
    return objects_size;
}


/*! \pre 0 <= index < size_of_objects()
 */
object_struct* Object_file::object( int index ) const
{
    CGAL_assertion( index >= 0 );
    CGAL_assertion( index < objects_size );

    return objects_list[index];
}


/*! \pre 0 <= index < size_of_objects()
 */
bool Object_file::is_polygons( int index ) const
{
    CGAL_assertion( index >= 0 );
    CGAL_assertion( index < objects_size );

    return get_object_type( objects_list[index] ) == POLYGONS;
}


/*! \pre 0 <= index < size_of_objects()
 */
bool Object_file::is_quadmesh( int index ) const
{
    CGAL_assertion( index >= 0 );
    CGAL_assertion( index < objects_size );

    return get_object_type( objects_list[index] ) == QUADMESH;
}



}
