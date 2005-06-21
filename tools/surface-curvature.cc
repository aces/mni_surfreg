/* Estimate surface curvature at vertices of triangulated
 * mesh.
 *
 * Ref: Gabriel Taubin
 *      Estimating the Tensor of Curvature of a Surface from a
 *      Polyhedral Approximation.
 */

#include <iostream>
#include <iterator>
#include <algorithm>
#include <ParseArgv.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/double.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>

#include <surflib/load_surface_file.hpp>
#include <surflib/vtk_ostream.hpp>
#include <surflib/SurfaceFunctions.hpp>


using namespace std;
using namespace MNI;


template <class Refs, class Traits>
struct My_vertex
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
					  typename Traits::Point_3>
{
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> 
    Vertex_base;

    My_vertex( const Point& p ) : Vertex_base(p) {}
    My_vertex() {}

    typename Traits::Vector_3 normal;
    double curvature[2];
    double error;
};


template <class Refs, class Traits>
struct My_face
    : public CGAL::HalfedgeDS_face_base<Refs>
{
    double area;
    typename Traits::Vector_3 normal;
};


struct My_polyhedron_items_3 : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Vertex_wrapper {
	typedef typename Traits::Point_3 Point;
	typedef My_vertex<Refs,Traits> Vertex;
    };

    template <class Refs, class Traits>
    struct Face_wrapper {
	typedef My_face<Refs,Traits> Face;
    };
};


typedef CGAL::Simple_cartesian<double>   Surface_kernel;
typedef CGAL::Polyhedron_3<
    Surface_kernel, 
    My_polyhedron_items_3,
    CGAL::HalfedgeDS_list
    > Surface;


typedef Surface_kernel::Vector_3  Vector_3;
typedef Surface::Halfedge_handle  Halfedge_handle;
typedef Surface::Vertex           Vertex;
typedef Surface::Vertex_iterator  Vertex_iterator;




template <class R>
void normalize_vector( CGAL::Vector_3<R>& v )
{
    double l_sq = CGAL::to_double( v*v );
    CGAL_assertion( l_sq > 0 );
    v = v / std::sqrt(l_sq);
}


//! Return vector pointing in direction of halfedge.
const Vector_3 vector_of_halfedge( Halfedge_handle h )
{
    return h->vertex()->point() - h->opposite()->vertex()->point();
}


/* Compute facet normal and area.
 *
 * Usage:
 *     std::for_each( surf.facets_begin(), surf.facets_end(),
 *                    Facet_property() );
 */
struct Facet_property 
{
    template <class Facet>
    void operator()( Facet& f) 
    {
        typename Facet::Halfedge_handle h = f.halfedge();

	Vector_3 u = vector_of_halfedge(h->next());
	Vector_3 v = vector_of_halfedge(h->opposite());
	
	f.normal = CGAL::cross_product(u,v);
	normalize_vector(f.normal);

	f.area = triangular_facet_area<Surface>( &f );
    }
};


/* Compute vertex normal.
 * The facet normal and areas must be valid.
 */
struct Vertex_property
{
    template <class Vertex>
    void operator()( Vertex& v )
    {
	v.normal = CGAL::NULL_VECTOR;

	typename Vertex::Halfedge_around_vertex_circulator 
	    h = v.vertex_begin();
	CGAL_For_all( h, v.vertex_begin() ) {
	    typename Vertex::Facet_handle f = h->facet();
	    v.normal = v.normal + (f->normal * f->area);
	}

	normalize_vector(v.normal);
    }
};


std::ostream& operator<< ( std::ostream& os, const double M[3][3] )
{
    os << "[";
    for( int i = 0; i < 3; ++i ) {
	os << "[";
	for( int j = 0; j < 3; ++j ) {
	    os << M[i][j] << " ";
	}
	os << "];" << std::endl;
    }
    return os << "]" << std::endl;
}


double frobenius_product( const double A[3][3],
			  const double B[3][3] )
{
    double p = 0;
    for( int i = 0; i < 3; ++i )
	for( int j = 0; j < 3; ++j )
	    p +=A[i][j]*B[j][i];
    return p;
}


double frobenius_norm( const double A[3][3] )
{
    return std::sqrt( frobenius_product(A,A) );
}


double taubin_error_unit_sphere( const double A[3][3] )
{
    static double M[3][3] = {{0,0,0},{0,-0.5,0},{0,0,-0.5}};

    return 1 - frobenius_product(A,M) 
	/ (frobenius_norm(A)*frobenius_norm(M));
}


/*! \brief Compute vertex curvature tensor.
 * 
 * \pre Vertex and Facet properties are valid.
 */
struct Vertex_curvature
{
    void clear_matrix()
    {
	for( int i = 0; i < 3; ++i )
	    for( int j = 0; j < 3; ++j )
		M[i][j] = 0;
    }

    void scale_matrix(double s)
    {
	for( int i = 0; i < 3; ++i )
	    for( int j = 0; j < 3; ++j )
		M[i][j] *= s;
    }

    /* Matrix M is the sum of  weight_k * T_k^t * T_k
     * Add one term to M.
     */
    void add_one( double weight, const Vector_3& t )
    {
	for( int i = 0; i < 3; ++i )
	    for( int j = 0; j < 3; ++j )
		M[i][j] += weight * t[i]*t[j];
    }

    /* Housholder transformation of matrix M.
     * Set M to Q^t.M.Q.
     */
    void transform( double Q[3][3] )
    {
	double r[3][3];

#if 0
	cout << "M = " << endl 
	     << M << endl
	     << "Q = " << endl 
	     << Q << endl;
#endif

	for( int i = 0; i < 3; ++i ) {
	    for( int j = 0; j < 3; ++j ) {
		r[i][j] = 0;
		for( int k = 0; k < 3; ++k ) {
		    for( int l = 0; l < 3; ++l ) {
			r[i][j] += Q[k][i]*M[k][l]*Q[l][j];
		    }
		}
	    }
	}

	for( int i = 0; i < 3; ++i )
	    for( int j = 0; j < 3; ++j )
		M[i][j] = r[i][j];
    }

    template <class Vertex>
    void operator()( Vertex& v )
    {
	clear_matrix();
	double weight_sum = 0;

	typename Vertex::Halfedge_around_vertex_circulator 
	    h = v.vertex_begin();
	CGAL_For_all( h, v.vertex_begin() ) {

	    // Let u denote the neighbour of v incident to h->opposite().
	    // Let vu be vector from v to u.
	    // Let t be unit-length projection of vu onto 
            //          plane with normal v.normal.
	    //
	    Vector_3 vu = vector_of_halfedge( h->opposite() );
	    Vector_3 t = vu - v.normal * (v.normal*vu);
	    normalize_vector(t);

	    // Approximate curvature along direction t.
	    double curv = 2.0 * (v.normal*vu) / (vu*vu);

	    // Weight of this term in the matrix is proportional to the sum of
	    // areas of facets incident to h.  (We are assuming here a closed
	    // surface).
	    //
	    double weight = h->facet()->area 
		+ h->opposite()->facet()->area;

	    add_one( weight*curv, t );
	    weight_sum += weight;
	}
	scale_matrix( 1.0 / weight_sum );

	Vector_3 Wplus = Vector_3(1,0,0) + v.normal;
	Vector_3 Wminus = Vector_3(1,0,0) - v.normal;

	Vector_3 W = Wplus;
	if ( Wminus*Wminus > Wplus*Wplus )
	    W = Wminus;
	normalize_vector(W);

	double Q[3][3];
	for( int i = 0; i < 3; ++i )
	    for( int j = 0; j < 3; ++j )
		Q[i][j] = (i==j ? 1 : 0) - 2.0*W[i]*W[j];

	transform(Q);
#if 0
	cout << "Transformed M = " << endl 
	     << M << endl;
#endif

#if 1
	// Set error to the sum of all entries that are
	// supposed to be zero.
	//
	v.error = std::abs(M[0][0]) + std::abs(M[0][1]) + std::abs(M[0][2])
	    + std::abs(M[1][0]) + std::abs(M[2][0]);
#else
	// Set error as in Tabuin's paper, assuming a unit sphere
	// input.
	//
	v.error = taubin_error_unit_sphere(M);
#endif

	// Find the eigenvalues of the lower, right 2x2 submatrix of M
	//
	double half_trace = (M[1][1]+M[2][2])/2.0;
	double determinant = M[1][1]*M[2][2] - M[1][2]*M[2][1];
	double discrim = std::sqrt(half_trace*half_trace - determinant);

	double e1 = half_trace + discrim;
	double e2 = half_trace - discrim;

	// Finally!  We have the curvatures.
	//
	v.curvature[0] = 3*e1 - e2;
	v.curvature[1] = 3*e2 - e1;

	// Sort by absolute value.
	if ( std::abs(v.curvature[0]) > std::abs(v.curvature[1]) )
	    std::swap( v.curvature[0], v.curvature[1] );
    }

private:
    double M[3][3];
};


struct normal: public unary_function<const Vertex&, Vector_3> 
{
    const Vector_3 operator()( const Vertex&v ) { return v.normal; }
};

struct curvature_0: public unary_function<const Vertex&, double> 
{
    double operator()( const Vertex&v ) { return v.curvature[0]; }
};

struct curvature_1: public unary_function<const Vertex&, double> 
{
    double operator()( const Vertex&v ) { return v.curvature[1]; }
};

struct curvature_mean: public unary_function<const Vertex&, double> 
{
    double operator()( const Vertex&v ) 
	{ return (v.curvature[0] + v.curvature[1])/2.0; }
};

struct curvature_gaussian: public unary_function<const Vertex&, double> 
{
    double operator()( const Vertex&v ) 
	{ return v.curvature[0] * v.curvature[1]; }
};

struct curvature_error: public unary_function<const Vertex&, double> 
{
    double operator()( const Vertex&v ) { return v.error; }
};


int main( int ac, char* av[] )
{
    char* vtk_file = 0;
    char* small_file = 0;
    char* large_file = 0;
    char* mean_file = 0;
    char* gauss_file = 0;

    ArgvInfo argTable[] = {
	{ "-small", ARGV_STRING, (char*)1, (char*)&small_file,
	  "write small curvature values in file, one per line" },
	{ "-large", ARGV_STRING, (char*)1, (char*)&large_file,
	  "write large curvature values in file, one per line" },
	{ "-mean", ARGV_STRING, (char*)1, (char*)&mean_file,
	  "write mean curvature values in file, one per line" },
	{ "-gauss", ARGV_STRING, (char*)1, (char*)&gauss_file,
	  "write Gaussian curvature values in file, one per line" },
	{ "-vtk", ARGV_STRING, (char*)1, (char*)&vtk_file,
	  "write output in VTK format (includes all curvatures)" },
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 2 ) {
	cerr << "usage: " << av[0] << " [options] surface" << endl;
	return 1;
    }

    if ( !vtk_file && !small_file && !large_file 
	 && !mean_file && !gauss_file) {
	cerr << "WARNING: need to specify -vtk, -small, etc"
	    " to get output." << endl;
    }

    try {
	Surface s;
	MNI::load_surface_file( s, av[1] );

	std::for_each( s.facets_begin(), s.facets_end(), 
		       Facet_property() );
	std::for_each( s.vertices_begin(), s.vertices_end(),
		       Vertex_property() );
	
	std::for_each( s.vertices_begin(), s.vertices_end(),
		       Vertex_curvature() );


	if ( vtk_file ) {
	    ofstream out( vtk_file );
	    string id( "$Id$ " );
	    write_vtk( out, s, id + string(av[1]) );

	    write_vtk_point_data_header( out, s );

	    write_vtk_point_vectors( out, s, "normals", "double", normal() );

	    write_vtk_point_scalars( out, s, "curv_small", curvature_0() );
	    write_vtk_point_scalars( out, s, "curv_large", curvature_1() );
	    write_vtk_point_scalars( out, s, "mean", curvature_mean() );
	    write_vtk_point_scalars( out, s, "gaussian", curvature_gaussian() );
	    write_vtk_point_scalars( out, s, "curv_error", curvature_error() );
	}

	if ( small_file ) {
	    ofstream out( small_file );
	    transform( s.vertices_begin(), s.vertices_end(),
		       ostream_iterator<double>(out,"\n"), curvature_0() );
	}

	if ( large_file ) {
	    ofstream out( large_file );
	    transform( s.vertices_begin(), s.vertices_end(),
		       ostream_iterator<double>(out,"\n"), curvature_1() );
	}

	if ( mean_file ) {
	    ofstream out( mean_file );
	    transform( s.vertices_begin(), s.vertices_end(),
		       ostream_iterator<double>(out,"\n"), curvature_mean() );
	}

	if ( gauss_file ) {
	    ofstream out( gauss_file );
	    transform( s.vertices_begin(), s.vertices_end(),
		       ostream_iterator<double>(out,"\n"), curvature_gaussian() );
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

