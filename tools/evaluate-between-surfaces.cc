/* Transfer labels from voxels to surface.
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

#include <surflib/load_surface_file.hpp>


extern "C" {
#   include <volume_io.h>
}


#include "Surface.hpp"
typedef Surface_base< CGAL::Simple_cartesian<double> >  Surface;


typedef Surface::Vertex_iterator Vertex_iterator;
typedef Surface::Point_3  Point_3;

using namespace std;


class QueryableLabelVolume
{
public:
    QueryableLabelVolume( char* fn,
			  int out_value = 0 )
	: outside_value(out_value)
    {
	// We need to get all volumes in x-y-z ordering
	static char* dimensions[] = { MIxspace, MIyspace, MIzspace };

	if ( input_volume( fn, 3, dimensions,
			   MI_ORIGINAL_TYPE, 0, 0, 0,
			   true, &vol, 0 ) != OK ) {
	    throw std::runtime_error( "failed to open volume" );
	}
	if ( get_volume_n_dimensions(vol) != 3 ) {
	    throw std::runtime_error( "expected volume with dimension = 3" );
	}

	// Ensure the labels are byte values.
	BOOLEAN flag;
	if ( get_volume_nc_data_type(vol,&flag) != NC_BYTE )
	    throw std::runtime_error( "expected NC_BYTE volume" );

	get_volume_sizes( vol, sizes );
	get_volume_separations( vol, steps );
    }

    int get( int i, int j, int k )
    {
	if ( i < 0 || i >= sizes[0] ||
	     j < 0 || j >= sizes[1] ||
	     k < 0 || k >= sizes[2] )
	    return outside_value;
	
	return (int) get_volume_real_value( vol, i,j,k,0,0 );
    }

    // Return intersection of query segment, given in world coords.
    int query_segment_world( const Point_3& p, const Point_3& q )
    {
	Real pv[MAX_DIMENSIONS];
	Real qv[MAX_DIMENSIONS];

	convert_world_to_voxel( vol, p.x(),p.y(),p.z(), pv );
	convert_world_to_voxel( vol, q.x(),q.y(),q.z(), qv );

	return query_segment( pv, qv );
    }


    // Return intersection of query segment, given in voxel coords.
    int query_segment( Real p[], Real q[] )
    {
	double dx = (q[0]-p[0]);
	double dy = (q[1]-p[1]);
	double dz = (q[2]-p[2]);

	// Let k be number of steps required so as not to miss any
	// voxels.
	int kx = static_cast<int>(1 + dx/steps[0]);
	int ky = static_cast<int>(1 + dy/steps[1]);
	int kz = static_cast<int>(1 + dz/steps[2]);

	int k = std::max(kx,std::max(ky,kz));
	k = static_cast<int>(k*1.5);

	double x = p[0];
	double y = p[1];
	double z = p[2];

	dx = dx / (k+1);
	dy = dy / (k+1);
	dz = dz / (k+1);

	clear_counts();
	for( int i = 0; i <= k; ++i ) {
	    counts[ get(ROUND(x),ROUND(y),ROUND(z)) ] += 1;
	    x += dx;
	    y += dy;
	    z += dz;
	}

	return largest_count();
    }


private:
    Volume vol;
    int sizes[MAX_DIMENSIONS];
    Real steps[MAX_DIMENSIONS];

    int outside_value;

    int counts[256];


    void clear_counts()
    {
	for( int i = 0; i < 256; ++i )
	    counts[i] = 0;
    }

    int largest_count()
    {
	int i_max = -1;
	for( int i = 1; i < 256; ++i ) {
	    if ( counts[i] > 0 ) {
		if ( i_max == -1 || (counts[i_max] < counts[i]) )
		    i_max = i;
	    }
	}
	return i_max > 0 ? i_max : 0;
    }

};


/* Query for label along segment (u,w), for all w adjacent to v.
 */
int label_to_neighbour( Surface::Vertex_handle u,
			Surface::Vertex_handle v,
			QueryableLabelVolume& v_label )
{
    Surface::Halfedge_around_vertex_circulator
	w = v->vertex_begin();

    CGAL_For_all(w,v->vertex_begin()) {
	int label = v_label.query_segment_world( u->point(),
						 w->opposite()->vertex()->point());
	if (label)
	    return label;
    }
    return 0;
}
    


/* For each vertex v of s1, let w be the corresponding vertex
 * of s2.  Set v->scalar equal to the label encountered most
 * frequently by the line segment (v,w) in the label volume.
 */
void set_vertex_values( Surface& s1,
			Surface& s2,
			QueryableLabelVolume& v_label )
{
    Surface::Vertex_iterator v = s1.vertices_begin();
    Surface::Vertex_iterator w = s2.vertices_begin();

    for( ; v != s1.vertices_end(); ++v,++w ) {
	v->scalar = v_label.query_segment_world( v->point(),
						 w->point() );
	if ( v->scalar == 0 )
	    v->scalar = label_to_neighbour( v, w, v_label );
	
	if ( v->scalar == 0 )
	    v->scalar = label_to_neighbour( w, v, v_label );
	
    }
}


int main( int ac, char* av[] )
{
    if ( ac != 5 ) {
	cerr << "usage: " << av[0] << " surface1 surface2 labels.mnc out.vv" 
	     << endl;
	return 1;
    }


    try {
	Surface s1,s2;
	MNI::load_surface_file( s1, av[1] );
	MNI::load_surface_file( s2, av[2] );

	cout << "Assuming objects " << av[1] 
	     << " and " << av[2] << " are corresponding (without checking)." 
	     << endl;

	QueryableLabelVolume v_label( av[3] );

	set_vertex_values( s1, s2, v_label );

	std::ofstream vv( av[4] );
	Surface::Vertex_iterator v = s1.vertices_begin();
	CGAL_For_all(v,s1.vertices_end())
	    vv << v->scalar << endl;

    } catch ( const std::bad_alloc& e ) {
	std::cerr << "Failed a memory allocation." << "\n"
		  << "No output." << "\n";
	return 2;
    } catch ( const std::exception& e ) {
	std::cerr << "Error: " << e.what() << "\n"
		  << "No output." << "\n";
        return 3;
    } catch ( ... ) {
	std::cerr << "Unknown exception." << "\n"
		  << "No output." << "\n"
		  << "This is likely bug in the code: please report!" << "\n";
        return 4;
    }

    return 0;
}
