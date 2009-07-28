/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <cmath>
#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/Polyhedron_3.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/SurfaceFunctions.hpp>


typedef CGAL::Simple_cartesian<double>      Surface_kernel;
typedef CGAL::Polyhedron_3<Surface_kernel>  Surface;
typedef Surface_kernel::Point_3             Point_3;
typedef Surface_kernel::Vector_3            Vector_3;

using namespace std;



/* Assumes triangular facets.
 */
void dump_facet_area( std::ostream& out, const Surface& s )
{
    double total_area = 0;

    Surface::Facet_const_iterator f = s.facets_begin();
    for( int i = 0; f != s.facets_end(); ++i,++f ) {
	Surface::Halfedge_const_handle h = f->halfedge();

	CGAL_assertion( h->next()->next()->next() == h );

	Point_3 v0 = h->vertex()->point();
	h = h->next();
	Point_3 v1 = h->vertex()->point();
	h = h->next();
	Point_3 v2 = h->vertex()->point();
	
	Vector_3 v = 0.5 * CGAL::cross_product(v1 - v0, v2 - v0);
	double v_len = sqrt( CGAL::to_double(v*v) );
	total_area += v_len;
	out << i << ": " << v_len << std::endl;
    }
    std::cerr << "Total Surface Area = " << total_area << std::endl;
}


/* Compute Voronoi cell for tail vertex of \a h.
 * Return area of cell intersecting the triangular facet h->facet().
 * Throw exception if circumcentre of triangle is not in triangle.
 */
inline
double voronoi_area( Surface::Halfedge_const_handle h )
{
    double t1, t2;
    MNI::circumcentre<Surface>( h, &t1,&t2 );
    MNI::Point_s<Surface> c( h, t1,t2);

    double area = MNI::triangular_facet_area<Surface>(h);

    Vector_3 v01( c.P1() - c.P0() );
    Vector_3 v02( c.P2() - c.P0() );

    return area * (2.0*t1*t2 + (v01*v02)*(t2*t2/(v01*v01)
					  +t1*t1/(v02*v02)));
}


/* Assumes triangular facets.
 */
void dump_vertex_area( std::ostream& out, const Surface& s )
{
    double total_area = 0;

    Surface::Vertex_const_iterator v = s.vertices_begin();
    for( int i = 0; v != s.vertices_end(); ++i,++v ) {

	Surface::Halfedge_around_vertex_const_circulator 
	    h = v->vertex_begin();

	double v_area = 0;
	CGAL_For_all( h,v->vertex_begin() ) {
	    v_area += voronoi_area( h->opposite() );
	}

	out << i << ": " << v_area << std::endl;

	total_area += v_area;
    }
    std::cerr << "Total Surface Area = " << total_area << std::endl;
}


/* Assumes triangular facets.
 * Check that circumcentre is located on facet.
 */
void check_circumcentre( const Surface& s )
{
    Surface::Facet_const_iterator f = s.facets_begin();
    for( int i = 0; f != s.facets_end(); ++i,++f ) {

	double t1, t2;
	MNI::circumcentre<Surface>( f, &t1,&t2 );

	MNI::Point_s<Surface> c( f->halfedge(), t1,t2);

#if 0
	cout << "V0 = " << c.P0() << endl
	     << "V1 = " << c.P1() << endl
	     << "V2 = " << c.P2() << endl
	     << "    c = " << c << endl;
#endif

	CGAL_assertion( t1 >= 0 );
	CGAL_assertion( t2 >= 0 );
	CGAL_assertion( t1+t2 <= 1 );

	Vector_3 c_v0 = c.point() - c.P0();
	Vector_3 c_v1 = c.point() - c.P1();
	Vector_3 c_v2 = c.point() - c.P2();

	cout << i << ": " << std::sqrt(CGAL::to_double(c_v0*c_v0)) 
	     << "  " << std::sqrt(CGAL::to_double(c_v1*c_v1)) 
	     << "  " << std::sqrt(CGAL::to_double(c_v2*c_v2)) 
	     << endl;
	    

	CGAL_assertion( abs(c_v0*c_v0 - c_v1*c_v1) < 1e-7 );
	CGAL_assertion( abs(c_v0*c_v0 - c_v2*c_v2) < 1e-7 );
    }
}


void dump_facet_angle( std::ostream& out, const Surface& s )
{
    double min_angle = 8;
    double max_angle = -1;
    int count_acute = 0;
    int count_obtuse = 0;

    Surface::Vertex_const_iterator v = s.vertices_begin();
    CGAL_For_all(v,s.vertices_end()) {
	Surface::Halfedge_around_vertex_const_circulator h = v->vertex_begin();
	CGAL_For_all(h,v->vertex_begin()) {
	    double angle = MNI::facet_angle<Surface>(h);
	    if ( angle < min_angle )
		min_angle = angle;
	    if ( angle > max_angle )
		max_angle = angle;
	    if ( angle < M_PI_2 )
		++count_acute;
	    else
		++count_obtuse;
	    out << angle << std::endl;
	}
    }
    std::cerr << "Minimum angle = " << min_angle << std::endl
	      << "Maximum angle = " << max_angle << std::endl
	      << "There are " << count_acute << " acute angles, and "
	      << count_obtuse << " obtuse." << std::endl;
}


void dump_edge_length( std::ostream& out, const Surface& s )
{
    Surface::Edge_const_iterator e = s.edges_begin();
    for( ; e != s.edges_end(); ++e ) {
	Point_3 p = e->vertex()->point();
	Point_3 q = e->opposite()->vertex()->point();
	out << sqrt( CGAL::to_double( CGAL::squared_distance( p,q )))
	    << std::endl;
    }
}


void usage( char* argv0 )
{
    char* s = strrchr( argv0, '/' );
    std::cerr << "usage: " << (s ? s+1 : argv0) << " [options] file [file...]"
	      << std::endl;
}

	
int main( int ac, char* av[] )
{
    int show_vertex_area = 0;
    int show_edge_length = 0;
    int show_face_area = 0;
    int show_face_angle = 0;
    char* out_filename = 0;

    static ArgvInfo argTable[] = {
	{ "-output", ARGV_STRING, (char*)1, (char*)&out_filename,
	  "output filename [default is stdout]" },
	{ "-vertex_area", ARGV_CONSTANT, (char*)1, (char*)&show_vertex_area,
	  "show vertex areas" },
	{ "-edge_length", ARGV_CONSTANT, (char*)1, (char*)&show_edge_length,
	  "show edge lengths" },
	{ "-face_area", ARGV_CONSTANT, (char*)1, (char*)&show_face_area,
	  "show face areas" },
	{ "-face_angle", ARGV_CONSTANT, (char*)1, (char*)&show_face_angle,
	  "show face angles" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 )) {
	usage( av[0] );
	return 1;
    }

    std::ostream* out = &std::cout;
    if ( out_filename ) {
	out = new std::ofstream(out_filename);
	if (!*out) {
	    cerr << "Cannot open output file " << out_filename << endl;
	    return 1;
	}
    }

    try {
	for( int i = 1; i < ac; ++i ) {
	    
	    Surface s;
	    MNI::load_surface_file( s, av[i] );

	    if ( show_vertex_area )
		dump_vertex_area( *out, s );
	    if ( show_edge_length )
		dump_edge_length( *out, s );
	    if ( show_face_area )
		dump_facet_area( *out, s );
	    if ( show_face_angle )
		dump_facet_angle( *out, s );
	}

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
