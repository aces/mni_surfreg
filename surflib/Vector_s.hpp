/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef VECTOR_S_H  // -*- C++ -*-
#define VECTOR_S_H


#include <surflib/Point_s.hpp>


namespace MNI {


/*! \brief Vector in tangent space of a surface location.
 *
 */

template <class _Surface>
class Vector_s
{
public:
    typedef _Surface                                 Surface;
    typedef Vector_s<Surface>                        Self;

    typedef typename Surface::Vertex_const_handle    Vertex_const_handle;
    typedef typename Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef typename Surface::Point_3                Point_3;
    typedef typename Surface::Traits::Vector_3       Vector_3;


    //! Direct construction.
    Vector_s( Halfedge_const_handle h, double t1, double t2);


    // Accessors
    const Facet_coordinate_frame<Surface>& coordinate_frame() const 
	{ return frame; }

    Halfedge_const_handle halfedge() const { return frame.halfedge(); }
    double t0() const { return 1-_t1-_t2; }
    double t1() const { return _t1; }
    double t2() const { return _t2; }

    // Get coordinates with respect to given halfedge.
    void coordinates_set_halfedge( Halfedge_const_handle h,
				   double* t1, double* t2 ) const;

    const Vector_3 vector() const;


    // Algebra operators.
    // Compute with vector() if you want lengths and stuff.
    const Self operator+(const Self &w) const;
    const Self operator-(const Self &w) const;
    const Self operator-() const;
    const Self operator/(double c) const;
    const Self operator*(double c) const;

private:
    Facet_coordinate_frame<Surface> frame;
    double _t1,_t2;

    // Test case
    friend class t_Vector_s;
};

}


#include "Vector_s.cpp"
#endif
