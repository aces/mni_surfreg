#ifndef SURFACE_H  // -*- C++ -*-
#define SURFACE_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/MP_Float.h>



template <class Refs, class Traits>
struct My_vertex
    : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true,
                                          typename Traits::Point_3>
{
    typedef typename Traits::Point_3 Point;
    typedef CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> 
    Vertex_base;

    My_vertex( const Point& p ) : Vertex_base(p) {};
    My_vertex() {};

    double scalar;
};

struct Polyhedron_items_3 : public CGAL::Polyhedron_items_3
{
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef My_vertex<Refs,Traits> Vertex;
    };
};


//typedef CGAL::Simple_cartesian<CGAL::MP_Float>   Surface_kernel;

typedef CGAL::Polyhedron_3<CGAL::Simple_cartesian<double>,
			   Polyhedron_items_3>  Surface_double;

typedef CGAL::Polyhedron_3<CGAL::Simple_cartesian<CGAL::MP_Float>,
			   Polyhedron_items_3>  Surface_mpfloat;

#endif
