/* Read file containing a label value for each vertex.  Identify
 * connected components of the vertex graph for which each vertex has
 * the same label.  Output new file with the regions labelled from
 * largest to smallest.
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>
#include <surflib/DisjointVertexSets.hpp>

#include "Surface.hpp"
typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;


typedef Surface::Vertex_handle          Vertex_handle;
typedef Surface::Vertex_iterator        Vertex_iterator;


using namespace std;


/* Sort in descending order.
 */
struct sort_by_component_size 
{ 
    MNI::DisjointVertexSets<Surface>& ds;

    sort_by_component_size( MNI::DisjointVertexSets<Surface>& s ) 
	: ds(s) {}

    bool operator()( Vertex_handle u, Vertex_handle v )
    { 
	return ds.component_size(u) > ds.component_size(v); 
    }
};


void compute_components( Surface& s )
{
    MNI::DisjointVertexSets<Surface> dset(s);

    // Join components when both sides of edge have
    // same label.

    Surface::Edge_iterator h = s.edges_begin();
    CGAL_For_all(h,s.edges_end()) {
	Surface::Vertex_handle u = h->vertex();
	Surface::Vertex_handle v = h->opposite()->vertex();
	if ( u->label == v->label )
	    dset.make_union(u,v);
    }

    // Sort resulting set of components by size.

    std::vector<Surface::Vertex_handle> component_vec;
    dset.components(s,std::back_inserter(component_vec));
    std::sort( component_vec.begin(), component_vec.end(),
	       sort_by_component_size(dset) );

    std::cout << "Surface has " << component_vec.size() 
	      << " components."
	      << std::endl;

    // Clear label values.

    Surface::Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v,s.vertices_end() ) {
	v->label = -1;
    }

    // Assign component identifiers to each vertex.

    std::vector<Vertex_handle>::iterator vhi = component_vec.begin();
    for( int i = 0; vhi != component_vec.end(); ++i,++vhi ) {
	(*vhi)->label = i;
    }

    v = s.vertices_begin();
    CGAL_For_all( v,s.vertices_end() ) {
	Vertex_handle v_set = dset.find_set(v);
	if ( v != v_set ) {
	    CGAL_assertion( v->label == -1 );
	    CGAL_assertion( v_set->label >= 0 );
	    v->label = v_set->label;
	}
    }
}



struct vertex_label 
    : public std::unary_function<Surface::Vertex&, int&>
{
    int& operator()( Surface::Vertex& v ) { return v.label; }
};


/* Read label for each vertex of s.
 */
void read_labels( Surface& s, 
		  char* vv_filename )
{
    std::ifstream in( vv_filename );

    Surface::Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v, s.vertices_end() ) {
	in >> v->label;
    }
}


void write_labels( Surface& s, 
		   char* vv_filename )
{
    std::ofstream out( vv_filename );
    transform( s.vertices_begin(), s.vertices_end(),
	       ostream_iterator<int>(out,"\n"), vertex_label() );
}


void usage( char* argv0 )
{
    char* s = strrchr( argv0, '/' );
    std::cerr << "usage: " << (s ? s+1 : argv0) 
	      << " surface in.vv out.vv"
	      << std::endl;
}


int main( int ac, char* av[] )
{
    if ( ac != 4 ) {
        usage( av[0] );
        return 1;
    }

    Surface s;
    MNI::load_surface_file( s, av[1] );

    read_labels( s, av[2] );
    compute_components(s);
    write_labels( s, av[3] );

    return 0;
}
