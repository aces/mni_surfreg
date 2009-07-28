/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

// -*- C++ -*-

#include <CGAL/basic.h>
#include <surflib/Barycentric.hpp>


namespace MNI {


template <class _Surface>
inline typename _Surface::Traits::Vector_3
Facet_coordinate_frame<_Surface>::e1() const
{
    return P1() - P0();
}

template <class _Surface>
inline typename _Surface::Traits::Vector_3
Facet_coordinate_frame<_Surface>::e2() const
{
    return CGAL::cross_product(n(),e1());
}

template <class _Surface>
inline typename _Surface::Traits::Vector_3
Facet_coordinate_frame<_Surface>::n() const
{
    Vector_3 v10 = P1() - P0();
    Vector_3 v20 = P2() - P0();
    return CGAL::cross_product(v10,v20);
}


template <class _Surface>
inline 
const typename _Surface::Point_3 
Facet_coordinate_frame<_Surface>::point( double t1, double t2 ) const
{
    return P0() + t1*(P1()-P0()) + t2*(P2()-P0());
}


template <class Surface_>
inline 
const typename Surface_::Traits::Vector_3 
Facet_coordinate_frame<Surface_>::vector( double t1, double t2 ) const
{
    return  t1*(P1()-P0()) + t2*(P2()-P0());
}


template <class _Surface>
inline void
Facet_coordinate_frame<_Surface>::coordinates_of( const Point_3& p,
						  double* t1, double* t2 ) const
{
    MNI::coordinates_of_vector( p-P0(), P1() - P0(), P2() - P0(), t1,t2 );
}


template <class _Surface>
inline void
Facet_coordinate_frame<_Surface>::coordinates_of( const Vector_3& v,
						  double* t1, double* t2 ) const
{
    // Solve the (overconstrained) system 
    //    v = t1*(P1-P0) + t2*(P2-P0)
    //
    // in least-squares sense.  Solve the following matrix equation
    // using generalized inverse
    // 
    //     v = t1*v10 + t2*v20     [3 equations, 2 uknowns]

    // TODO Change this to use the ring type up until the very last step.
    
    MNI::coordinates_of_vector( v, P1() - P0(), P2() - P0(), t1,t2 );
}


}
