/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef FACET_COORDINATE_FRAME_H  // -*- C++ -*-
#define FACET_COORDINATE_FRAME_H


namespace MNI {


/*! \brief Coordinate frame on plane tangent to triangular facet.
 *
 * Location on the surface is given by a halfedge together with
 * the areal (normalized barycentric) coordinates of the point
 * in the halfedge's facet.
 *
 *              V2
 *             /  \
 *            /    \
 *           V0----V1
 *               
 * The halfedge, h is (V0,V1).  V0, V1, V2 refer to the polyhedron
 * vertex structures.  The corresponding points in 3-space are P0, P1,
 * and P2.  
 *
 * A point on this plane is specified using areal coordinates
 * (t0,t1,t2), but t0+t1+t2 = 1 so t0 is not stored.  This point has
 * 3-space coordinates t0*P0 + t1*P1 + t2*P2.
 *
 * Areal coordinates are also used for vectors in this plane.
 * Coordinates (t1,t2) coorespond to the 3-space vector 
 * t1*(P1-P0) + t2*(P2-P0).
 *
 * \note Surface point structures are invalidated if the facet
 * handle is invalidated, if a vertex is inserted into the facet,
 * or if the vertex coordinates are changed.
 */

template <class _Surface>
class Facet_coordinate_frame
{
public:
    typedef _Surface Surface;

    typedef typename Surface::Vertex_const_handle    Vertex_const_handle;
    typedef typename Surface::Halfedge_const_handle  Halfedge_const_handle;
    typedef typename Surface::Point_3                Point_3;
    typedef typename Surface::Traits::Vector_3       Vector_3;


    /*! \brief Construct frame.
     *
     * \pre h is not border halfedge
     * \pre h->facet() is triangle.
     */
    Facet_coordinate_frame( Halfedge_const_handle h )
    {
	CGAL_precondition( !h->is_border() );
	CGAL_precondition( h->next()->next()->next() == h );
	this->h = h;
    }

    //! Default constructor returns singular value.
    // Required for use with containers.
    Facet_coordinate_frame() : h(0) {}


    // Accessors
    Halfedge_const_handle halfedge() const { return h; }

    Vertex_const_handle V0() const  { return h->opposite()->vertex(); }
    Vertex_const_handle V1() const  { return h->vertex(); }
    Vertex_const_handle V2() const  { return h->next()->vertex(); }

    const Point_3& P0() const { return V0()->point(); }
    const Point_3& P1() const { return V1()->point(); }
    const Point_3& P2() const { return V2()->point(); }


    //@{
    /*! \brief Coordinate frame in 3-space.
     *
     * The 3-space vectors (e1,e2,n) form a right-handed orthogonal
     * coordinate frame.  The vectors are \em not normalized.
     */
    Vector_3 e1() const;
    Vector_3 e2() const;
    Vector_3 n() const;
    //@}


    //! Compute 3-space point from areal coordinates.
    const Point_3 point(double t1, double t2) const;

    //! Compute 3-space vector from areal coordinates.
    const Vector_3 vector(double t1, double t2) const;

    //! Compute areal coordinates 3-space point.
    void coordinates_of( const Point_3& p, double* t1, double* t2 ) const;

    //! Compute areal coordinates of 3-space vector.
    void coordinates_of( const Vector_3& v, double* t1, double* t2 ) const;


private:
    Halfedge_const_handle h;

    // Test cases.
    friend class t_Facet_coordinate_frame;
};


}


#include "Facet_coordinate_frame.cpp"
#endif

