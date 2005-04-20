/* Input: MINC label file
 *        cortical surface, s
 *        numerical label, L
 *
 * Output: vertex value file with 1 if the vertex is the nearest
 *         (Euclidean distance) to a voxel with label L
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>

extern "C" {
#   include <volume_io.h>
}

#include <surflib/load_surface_file.hpp>

#include "Surface.hpp"
typedef CGAL::Simple_cartesian<double>  K;
typedef Surface_base< K >               Surface;

typedef K::Point_3 Point_3;
typedef Surface::Vertex_handle Vertex_handle;

// Store a handle and the SQUARED distance between the point and vertex.
typedef std::pair<Vertex_handle,double> Vertex_distance_pair;



using namespace std;


struct vertex_label 
    : public std::unary_function<Surface::Vertex&, int&>
{
    int& operator()( Surface::Vertex& v ) { return v.label; }
};


struct vdp_init 
    : public std::unary_function<const Point_3&,Vertex_distance_pair>
{
    Vertex_handle v;

    vdp_init( Vertex_handle v0 ) : v(v0) {}

    Vertex_distance_pair operator() ( const Point_3& p )
    {
	double d2 = CGAL::to_double( CGAL::squared_distance(p,v->point()) );
	return Vertex_distance_pair( v, d2 );
    }
};


/* Relax the p_near value if v->point() is nearer to p.
 */
inline void relax( const Point_3& p,
		   Vertex_distance_pair& p_near,
		   Vertex_handle v )
{
    double d2 = CGAL::to_double( CGAL::squared_distance(p,v->point()) );
    if ( d2 < p_near.second ) {
	p_near.first = v;
	p_near.second = d2;
    }
}


void scan_surface( Surface& s, 
		   vector<Point_3>& label_points )
{
    vector<Vertex_distance_pair> nearest_vertex;
    std::transform( label_points.begin(), label_points.end(),
		    std::back_inserter(nearest_vertex),
		    vdp_init(s.vertices_begin()) );

    int len = label_points.size();
    cout << "Label list size = " << len << endl
	 << "Distance list size = " << nearest_vertex.size() << endl;
       
    if ( label_points.size() != nearest_vertex.size() )
	throw std::logic_error
	    ( "label_points.size() != nearest_vertex.size()" );

    Surface::Vertex_iterator v = s.vertices_begin();
    CGAL_For_all( v,s.vertices_end() ) {
	for( int i = 0; i < len; ++i ) {
	    relax( label_points[i], 
		   nearest_vertex[i],
		   v );
	}
    }

    v = s.vertices_begin();
    CGAL_For_all( v,s.vertices_end() ) {
	v->label = 0;
    }

    for( int i = 0; i < len; ++i )
	nearest_vertex[i].first->label = 1;
}


/* Scan volume for given label, return vector of
 * world coordinate values.
 */
void get_label_point_list( char* fn,
			   int label,
			   std::vector<Point_3>& points )
{
    // We need to get all volumes in x-y-z ordering
    static char* dimensions[] = { MIxspace, MIyspace, MIzspace };

    Volume vol;
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


    // Volume is acceptable.  Start scan.

    int sizes[MAX_DIMENSIONS];
    get_volume_sizes( vol, sizes );

    for( int i0 = 0; i0 < sizes[0]; ++i0 )
	for( int i1 = 0; i1 < sizes[1]; ++i1 )
	    for( int i2 = 0; i2 < sizes[2]; ++i2 ) {
		double val = get_volume_real_value( vol, i0,i1,i2,0,0 );
		if ( static_cast<int>(val) == label ) {
		    double x,y,z;
		    convert_3D_voxel_to_world( vol, i0,i1,i2, &x,&y,&z );
		    points.push_back(Point_3(x,y,z));
		}
	    }

    delete_volume( vol );
}



int main( int ac, char* av[] )
{
    ArgvInfo argTable[] = {
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 5 ) {
	char* s = strrchr( av[0], '/' );
	std::cerr << "usage: " << (s ? s+1 : av[0]) 
		  << " [options] labels.mnc label surface out.vv"
		  << std::endl;
        return 1;
    }

    try {
	vector<Point_3> label_points;
	get_label_point_list( av[1], atoi(av[2]), label_points );

	Surface s;
	MNI::load_surface_file( s, av[3] );

	scan_surface( s, label_points );

	ofstream out( av[4] );
	transform( s.vertices_begin(), s.vertices_end(),
		   ostream_iterator<int>(out,"\n"),
		   vertex_label() );

    } catch ( const std::exception& e ) {
 	std::cerr << "std::exception: " << e.what() << std::endl;
	return 1;
    } catch ( ... ) {
	std::cerr << "Yikes!  Unknown exception!!" << std::endl;
	return 2;
    }


    return 0;
}
