/* Read a surface and a point.  Emit surface location.
 */

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <algorithm>

#include <surflib/Point_s.hpp>
#include <surflib/VTKFile.hpp>

#include <ParseArgv.h>



typedef CGAL::Cartesian<double>     Kernel;
typedef Kernel::Point_3             Point_3;
typedef Kernel::Vector_3            Vector_3;
typedef Kernel::Plane_3             Plane_3;
typedef CGAL::Polyhedron_3<Kernel>  Surface;

typedef MNI::Point_s<Surface>        Location;

using namespace std;




/* This creates a plane passing through three points
 * on the facet.  The points are oriented in a positive
 * sense when seen from the positive side of the plane.
 */
struct Plane_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) 
    {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3  Plane;
        return Plane( h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
    }
};



void read_vtk_file( const std::string filename, Surface& s )
{
    MNI::VTKFile vtk_in( filename );
    vtk_in.read_geometry(s);
}


/* Find a Point_s on surface s corresponding to point p.
 * This version assumes a convex surface.
 *
 * Let f = (V0,V1,V2) be a facet and Pf be its supporting plane.
 * First, locate f for which p is on the positive side of Pf.
 * For each edge = (A,B) of f, let C = A + Pf.orthogonal_vector()
 * and check the orientation of (A,B,C,p).  If the orientation is
 * positive, then p lies outside the triangle.  We update f to
 * be the facet on the other side of edge (A,B).  This new facet 
 * is guaranteed to be closer to p than the original facet.  The
 * distance to facet f forms a decreasing sequence, ending with zero.
 * (Or near zero for non-exact arithmetic...)
 */
void locate_point_convex( const Point_3& p,
			  const Surface& s )
{
    typedef Surface::Vertex_const_handle    V_handle;
    typedef Surface::Halfedge_const_handle  H_handle;
    typedef Surface::Facet_const_handle     F_handle;

    // Look for starting facet.
    F_handle f = s.facets_begin(); 
    int facets_examined = 0;
    while ( f != s.facets_end() ) {
	++facets_examined;
	if ( f->plane().oriented_side(p) != CGAL::ON_NEGATIVE_SIDE )
	    break;
	++f;
    }
    cerr << "Examined " << facets_examined << " facets." << endl;
    CGAL_assertion( f != s.facets_end() );

    // If the surface had no coplanar facets, we would be done.
    // We assume facets may be coplanar, so we start walking over
    // facets until p lies in interior of f.
    bool located_f = false;
    while (!located_f) {
	located_f = true;  // check each edge for falsification.

	Vector_3 normal = f->plane().orthogonal_vector();
	Surface::Halfedge_around_facet_const_circulator h0( f->halfedge() );
	Surface::Halfedge_around_facet_const_circulator h = h0;

	CGAL_For_all( h, h0 ) {
	    const Point_3& A = h->opposite()->vertex()->point();
	    const Point_3& B = h->vertex()->point();
	    Point_3 C = A + normal;
	    
	    if ( CGAL::orientation( A,B,C, p ) == CGAL::POSITIVE ) {
		located_f = false;
		f = h->opposite()->facet();
		break;
	    }
	}
    }

    Surface::Halfedge_const_handle h = f->halfedge();
    cerr << "Found point on facet: " << endl
	 << "    P0 = " << h->opposite()->vertex()->point() << endl
	 << "    P1 = " << h->vertex()->point() << endl
	 << "    P2 = " << h->next()->vertex()->point() << endl;

    MNI::Point_s<Surface> loc( h, p );

    cerr << "Location = " << loc << endl;
}


int main( int ac, char* av[] )
{
    CGAL::set_pretty_mode( std::cout);

    static ArgvInfo argTable[] = {
	{ NULL, ARGV_END, NULL, NULL, NULL }
    };

    if ( ParseArgv( &ac, av, argTable, 0 ) || ac != 5 ) {
	cerr << "usage: " << av[0] << " surface x y z" << endl;
	return 1;
    }

    Surface surf;
    read_vtk_file( av[1], surf );
    std::transform( surf.facets_begin(), surf.facets_end(), surf.planes_begin(),
                    Plane_equation());
	

    Point_3 p( atof(av[2]), atof(av[3]), atof(av[4]) );

    locate_point_convex( p, surf );

    return 0;
}
