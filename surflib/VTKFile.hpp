#ifndef VTKFILE_H  // -*- C++ -*-
#define VTKFILE_H

#include <string>
#include <fstream>
#include <vector>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <surflib/FileReader.hpp>


namespace MNI {


namespace detail {

template <class HDS>
struct read_vtkfile_geometry 
    : public CGAL::Modifier_base<HDS>, 
      private FileReader
{
    typedef typename HDS::Vertex    Vertex;
    typedef typename Vertex::Point  Point;

    std::ifstream& file;
    std::vector<Point> points;

    void read_point()
    {
	double x = expect_double();
	double y = expect_double();
	double z = expect_double();
	points.push_back( Point(x,y,z) );
    }


    read_vtkfile_geometry( std::ifstream& f ) : FileReader(f),file(f) {}
    void operator() ( HDS& hds )
    {
	std::string s;

	expect_word("POINTS");
	int n_points;
	file >> n_points;
	file >> s;  // datatype.  Ignore and use doubles.

	for( int i = 0; i < n_points; ++i ) {
	    read_point();
	}

	expect_word("POLYGONS");
	int n_polygons, facet_list_size;
	file >> n_polygons;
	file >> facet_list_size;

	CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true );
	B.begin_surface( n_points, n_polygons, (facet_list_size-n_polygons) );
	typename std::vector<Point>::iterator v_it = points.begin();
	while ( v_it != points.end() ) {
	    B.add_vertex( *v_it );
	    ++v_it;
	}

	for( int f = 0; f < n_polygons; ++f ) {
	    int n_vert;
	    file >> n_vert;
	    --facet_list_size;

	    B.begin_facet();
	    while ( n_vert-- > 0 ) {
		int v;
		file >> v;
		--facet_list_size;
		B.add_vertex_to_facet(v);
	    }
	    B.end_facet();
	}
	if ( facet_list_size != 0 )
	    throw FileReaderException( "facet list count is off by " 
				       + facet_list_size );
	B.end_surface();
    }
};

}


class VTKFile : private FileReader
{
public:
    VTKFile( const char* filename );
    VTKFile( const std::string filename );

    std::string get_header();
    bool is_binary();

    /// Read geometry into CGAL Polyhedron_3.
    template <class Polyhedron>
    void read_geometry( Polyhedron& P )
    {
	enter_state_geometry();
	typedef typename Polyhedron::HalfedgeDS HDS;
	detail::read_vtkfile_geometry<HDS> reader(file);
	P.delegate( reader );
	n_point = P.size_of_vertices();
	n_cell = P.size_of_facets();
    }

    //! Return name of scalar data attribute.
    const std::string expect_scalars();

    //! Scan forward to find requested scalar data attribute.
    void expect_scalars( const std::string& data_name );


    template <class Op>
    const std::string read_scalar_point_attribute( Op consume_value )
    {
	std::string data_name = expect_scalars();
	for( size_t i = 0; i < n_point; ++i ) {
	    consume_value( expect_double() );
	}
	return data_name;
    }


    template <class Op>
    void read_scalar_point_attribute( const char* att_name,
				      Op consume_value )
    {
	expect_scalars( att_name );
	for( size_t i = 0; i < n_point; ++i ) {
	    consume_value( expect_double() );
	}
    }


private:
    std::ifstream file;

    // In VTK-speak, "header" is the comment string or title of the data set.
    std::string header;
    bool binary_flag;

    // Sets all of the above.
    void read_header();

    // Set when geometry is read.
    size_t n_point, n_cell;

    enum FileState { Header, Geometry, PointData, CellData };
    FileState state;

    void initialize( const std::string& filename );
    void enter_state_geometry();
    void enter_state_point_data();
    //TODO: void enter_state_cell_data();
};



/*! Functor for setting the scalar value of each vertex.
 *  Useful for reading a VTK file.
 */
template <class Surface_>
class SetScalarValues
{
    Surface_& s;
    typename Surface_::Vertex_iterator v;

public:
    SetScalarValues( Surface_& surf ) : s(surf) 
    {
	v = s.vertices_begin();
    }

    void operator() ( double x ) 
    {
	CGAL_precondition( v != s.vertices_end() );
	v->scalar = x;
	++v;
    }
};


}

#endif
