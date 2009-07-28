/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <stdexcept>
#include <surflib/VTKFile.hpp>

using std::ifstream;


namespace MNI {


VTKFile::VTKFile( const char* filename )
    : FileReader( file ), file( filename )
{
    initialize( filename );
}

VTKFile::VTKFile( const std::string filename )
    : FileReader(file), file(filename.c_str())
{
    initialize( filename );
}

void VTKFile::initialize( const std::string& filename )
{
    if ( !file )
	throw std::runtime_error( "failed to read VTK file [" 
				  + filename + "]" );
    file.exceptions( ifstream::failbit | ifstream::badbit );
    read_header();
}

void VTKFile::read_header()
{
    expect_line( "# vtk DataFile Version 2.0" );

    std::getline( file, header );

    std::string s;
    std::getline( file, s );
    if ( s == "ASCII" ) 
	binary_flag = false;
    else if ( s == "BINARY" )
	binary_flag = true;
    else throw FileReaderException( "bad format string [" + s 
				    + "].  Expected ASCII or BINARY." );

    n_point = n_cell = 0;
    state = Header;
}


std::string VTKFile::get_header()
{
    return header;
}


bool VTKFile::is_binary()
{
    return binary_flag;
}


void VTKFile::enter_state_geometry()
{
    if ( state != Header )
	throw FileReaderException( "geometry must follow header" );
    expect_line("DATASET POLYDATA");
    state = Geometry;
}


void VTKFile::enter_state_point_data()
{
    if ( state == PointData )
	return;

    expect_word("POINT_DATA");
    size_t n = expect_int();
    if ( n != n_point )
	throw FileReaderException( "mismatch in number of point data" );

    state = PointData;
}


/*! Move file to beginning of scalar data.
 */
const std::string VTKFile::expect_scalars()
{
    enter_state_point_data();
	
    expect_word("SCALARS");
    std::string data_name = expect_word();
    expect_line();  // FIXME. Check numComp == 1, if specified.
    expect_word("LOOKUP_TABLE");
    expect_line();  // ignore name of lookup table

    return data_name;
}


/*! Move file to beginning of scalar data.
 */
void VTKFile::expect_scalars( const std::string& data_name )
{
    while(1) {
	std::string s = expect_scalars();
	if ( s == data_name )
	    break;
	// else skip over the data
	for( size_t i = 0; i < n_point; ++i ) {
            expect_double();
	}
    }
}


}
