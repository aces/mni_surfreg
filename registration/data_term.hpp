/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
// Similarity between two sequences.

/*! \brief Correlation coefficient between two sequences.
 * 
 * \pre (i1,end1] is a sequence of doubles
 * \pre i2 points to a sequence of doubles of same length
 */
template <class InIter>
double correlation_coefficient( InIter i1, InIter end1,
				InIter i2 )
{
    double s1 = 0, // sum of sequence 1
	s2 = 0,    // sum of sequence 2
	s3 = 0,    // sum of squared sequence 1
	s4 = 0,    // sum of squared sequence 2
	s5 = 0;    // sum of cross products
    int n = 0;     // number of samples

    for( ; i1 != end1; ++i1,++i2 ) {
	s1 += *i1;
	s2 += *i2;
	s3 += (*i1)*(*i1);
	s4 += (*i2)*(*i2);
	s5 += (*i1)*(*i2);
	++n;
    }

    double mean_1 = s1 / n;
    double mean_2 = s2 / n;
    double var_1 = s3 / n - mean_1*mean_1;
    double var_2 = s4 / n - mean_2*mean_2;
    double covariance = s5 / n - mean_1*mean_2;

    return covariance / std::sqrt( var_1*var_2 );
}


