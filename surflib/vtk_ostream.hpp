#ifndef VTK_OSTREAM_H  // -*- C++ -*-
#define VTK_OSTREAM_H  

#include <string>
#include <iostream>
#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Inverse_index.h>


namespace MNI {


inline 
void write_vtk_header( std::ostream& os, const std::string& message )
{
    os << "# vtk DataFile Version 2.0" << std::endl
       << message << std::endl
       << "ASCII" << std::endl
       << "DATASET POLYDATA" << std::endl;
}


template <class PointIterator>
void write_vtk_points( std::ostream& os,
		       int n_points,
		       PointIterator begin, PointIterator end )
{
    os << "POINTS " << n_points << " double" << std::endl;
    for( PointIterator i = begin; i != end; ++i )
    {
	os << CGAL::to_double(i->x()) << ' '
	   << CGAL::to_double(i->y()) << ' '
	   << CGAL::to_double(i->z()) << std::endl;
    }
}


template <class Polyhedron>
inline
void write_vtk_points( std::ostream& os, const Polyhedron& P )
{
    write_vtk_points( os, P.size_of_vertices(),
		      P.points_begin(), P.points_end() );
}


template <class Polyhedron>
inline
void write_vtk_point_data_header( std::ostream& os, const Polyhedron& P )
{
    os << "POINT_DATA " << P.size_of_vertices() << std::endl;
}


inline
void write_vtk_point_scalars_header( std::ostream& os,
				     const std::string& name,
				     const std::string& type_name )
{
    os << "SCALARS " << name << " " << type_name << std::endl
       << "LOOKUP_TABLE default" << std::endl;
}


template <class Polyhedron, class VertexOp>
void write_vtk_point_scalars( std::ostream& os, const Polyhedron& P,
			      const std::string& name,
			      VertexOp f )
{
    write_vtk_point_scalars_header( os, name, "double" );
    std::transform( P.vertices_begin(), P.vertices_end(),
		    std::ostream_iterator<double>(os,"\n"),
		    f );
}


template <class Polyhedron, class VertexOp>
void write_vtk_point_vectors( std::ostream& os, const Polyhedron& P,
			      const std::string& name,
			      const std::string& type_name,
			      VertexOp f )
{
    os << "VECTORS " << name << " " << type_name << std::endl;
    std::transform( P.vertices_begin(), P.vertices_end(),
		    std::ostream_iterator<typename VertexOp::result_type>(os,"\n"),
		    f );
}


template <class Polyhedron, class VertexOp>
void write_vtk_point_normals( std::ostream& os, const Polyhedron& P,
			      const std::string& name,
			      const std::string& type_name,
			      VertexOp f )
{
    os << "NORMALS " << name << " " << type_name << std::endl;
    std::transform( P.vertices_begin(), P.vertices_end(),
		    std::ostream_iterator<typename VertexOp::result_type>(os,"\n"),
		    f );
}


template <class Polyhedron>
void write_vtk_cells( std::ostream& os, const Polyhedron& P )
{
    typedef typename Polyhedron::Vertex_const_iterator       Vertex_iterator;
    typedef typename Polyhedron::Facet_const_iterator        Facet_iterator;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator
            Halfedge_around_facet_circulator;

    CGAL::Inverse_index<Vertex_iterator> vertex_index( P.vertices_begin(),
						       P.vertices_end() );

    // Write the cell list.
    int cell_list_size = 0;
    for ( Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f) 
    {
        Halfedge_around_facet_circulator h = f->facet_begin();
	cell_list_size += 1 + CGAL::circulator_size(h);
    }

    os << "POLYGONS " << P.size_of_facets() << ' ' 
       << cell_list_size << std::endl;

    for ( Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f) 
    {
        Halfedge_around_facet_circulator h = f->facet_begin();
	os << CGAL::circulator_size(h);
        do {
	    os << ' ' << vertex_index[h->vertex()];
        } while ( ++h != f->facet_begin());
	os << std::endl;
    }
}


template <class Polyhedron>
void write_vtk( std::ostream& os, const Polyhedron& P,
		const std::string& message )
{
    write_vtk_header( os, message );
    write_vtk_points( os, P );
    write_vtk_cells( os, P );
}


template <class Polyhedron>
void write_vtk( const char* filename, const Polyhedron& P,
		const std::string& message )
{
    std::ofstream out(filename);
    write_vtk( out, P, message );
}


}



#endif
