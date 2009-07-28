/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* Statistics on 3-vectors (or 3-points).
 */

#ifndef SURFLIB_VECTOR_STATISTIC_HPP
#define SURFLIB_VECTOR_STATISTIC_HPP

#include <surflib/Statistic.hpp>


namespace MNI {


/*! 3D Vector statistic.
 */

template <class V>
class VectorStatistic 
{
public:
    void clear() { coord[0].clear(); coord[1].clear(); coord[2].clear(); }

    void add_sample( const V& val )
    {
	coord[0].add_sample( val[0] );
	coord[1].add_sample( val[1] );
	coord[2].add_sample( val[2] );
    }

    int count() const { return coord[0].count(); }
    const V sum() const 
    {
	return V( coord[0].sum(), coord[1].sum(), coord[2].sum() ); 
    }

    const V mean() const 
    { 
	return V( coord[0].mean(), coord[1].mean(), coord[2].mean() );
    }

    double var() const
    {
	return coord[0].var() + coord[1].var() + coord[2].var();
    }

    double std_dev() const { return std::sqrt(var()); }

    double min( int i ) const { return coord[i].min(); }
    double max( int i ) const { return coord[i].max(); }


private:
    Statistic<double> coord[3];
};


}

#endif
