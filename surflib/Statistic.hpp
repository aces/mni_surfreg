#ifndef BIC_STATISTIC_H
#define BIC_STATISTIC_H

#include <iostream>


namespace MNI {


/*! \brief Descriptive statistics of type T.
 *
 * \note T must be convertible to double.
 */


template <class T>
class Statistic {

    T _min, _max, _sum, sum2;
    int _count;

public:

    //! Construct with empty set of samples.
    Statistic() { clear(); }

    //! Clear set of samples.
    void clear() { _min = _max = _sum = sum2 = 0;  _count = 0; }

    //! Add one value to the set of samples.
    void add_sample( T val ) {
	_sum += val;  sum2 += (val*val);  ++_count; 

	if ( _count == 1 ) {
	    _min = _max = val;
	} else {
	    if ( val < _min )  _min = val;
	    if ( val > _max )  _max = val;
	}
    }

    //! Size of sample set.
    int count() const { return _count; }

    //! Total of samples.
    T sum() const { return _sum; }

    //! Mean value of samples.
    T mean() const { return _sum / (T)_count; }

    //! Variance of samples.
    T var() const { T m = mean();  return sum2 / (T)_count - m*m; }

    //! Standard deviation of samples.
    T std_dev() const { return sqrt(var()); }

    //! Minimum value in sample.
    T min() const { return _min; }

    //! Maximum value in sample.
    T max() const { return _max; }

    //! Pretty-print a summary of the statistics.
    template <class OSTREAM>
    void print( OSTREAM& ) const;
};



template <class T>
std::ostream& operator<<( std::ostream& os, const Statistic<T>& stat )
{
    return os << stat.mean() 
	      << " (" << stat.std_dev() << "), N = " << stat.count();
}


template <class T>
template <class OSTREAM>
inline
void Statistic<T>::print( OSTREAM& os ) const
{
    os << count() << " samples in range [" 
       << min() << "," << max() << "]"
       << ", mean = " << mean() << " (" << std_dev() << ")" << std::endl;
}


}
#endif


