/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* This file is a convenience for defining a surface structure
 * that contains two data items at each vertex:
 *
 *     double scalar;
 *     int label;
 *
 * Usage:
 *     #include "Surface.hpp"
 *     typedef Surface_base< CGAL-KERNEL-CLASS >  Surface;
 */

#include <CGAL/Polyhedron_3.h>


namespace detail {

template <class Refs, class Traits>
struct My_vertex 
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
					  typename Traits::Point_3>
{
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> Vertex_base;

    My_vertex( const Point& p ) : Vertex_base(p) { scalar = 0; label = 0; };
    My_vertex() { scalar = 0; label = 0; };

    double scalar;
    int label;
};

struct My_items : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Vertex_wrapper {
	typedef typename Traits::Point_3 Point;
	typedef My_vertex<Refs,Traits> Vertex;
    };
};

}


template <class K>
class Surface_base : public CGAL::Polyhedron_3<K,detail::My_items>
{};

