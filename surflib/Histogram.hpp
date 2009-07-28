/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* -*- C++ -*- */

#ifndef MNI_HISTOGRAM_H
#define MNI_HISTOGRAM_H

#include <iostream>
#include <stdexcept>


namespace MNI {


/*! \brief Histogram of items of type T.
 *
 */


template <class T>
class Histogram {
public:

//*** Constructors

    //! Divides the interval [min,max] into num_bins.
    Histogram( int num_bins, T min, T max );


//*** Accessors

    int getNumBins( void ) const;
    T getMin( void ) const;
    T getMax( void ) const;
    double getBinWidth( void ) const;


//*** Accumulate data
    
    //! Enter a value into histogram
    void accumulate( T val );

    //! Reset all bin counts
    void reset( void );


//*** Data readout

    //! Number of counts in particular bin
    int getCounts( int bin ) const;

    //! Total of counts contained in first k bins (plus underflow).
    int getCumulativeCounts( int k ) const;

    //! Total number of counts seen
    int getTotalCounts( void ) const;

    //! Number of counts in range [min,max]
    int getValidCounts( void ) const;

    //! Number of counts < min, or > max
    int getUnderflow( void ) const;
    int getOverflow( void ) const;

    //! Return i'th [0...numBin-1] bin# of frequency-sorted sequence
    int getSortedBin( int rank );


    //! Output to stream.
    void dump( std::ostream& os ) const;

//*** Auxiliary

    //! Left edge of bin interval.
    double getBinLow( int bin ) const;

    //! Centre value of bin.
    double getBinCentre( int bin ) const;

    //! Right edge of bin interval.
    double getBinHigh( int bin ) const;

    //! Into which bin would this value fall?
    int findBin( T value ) const;

    //! Return smallest b for which getCumulativeCounts(b) >= counts.
    // FIXME: this is a terrible name
    int findBinForCumulativeCounts( int counts ) const;

private:
    int numBins;
    int totalCounts, underflow, overflow, *counts;
    T min, max;
    double width;

    // May desire bins in frequency-sorted order
    void frequencySort( int lbin, int hbin );
    int* sortedIndex;
    int validSort;
};



template <class T>
Histogram<T>::Histogram( int nBins, T minVal, T maxVal )
{
    numBins = nBins;
    min = minVal;
    max = maxVal;

    if ( maxVal < minVal ) { 
        // user is surely confused ... 
	// do we complain or silently fix?
	min = maxVal;
	max = minVal;
    }

    width = (max-min)/numBins;
    counts = new int[numBins];
    sortedIndex = 0;

    reset();
}


template <class T>
int Histogram<T>::getNumBins( void ) const
{
    return numBins;
}

template <class T>
T Histogram<T>::getMin( void ) const
{
    return min;
}

template <class T>
T Histogram<T>::getMax( void ) const
{
    return max;
}

template <class T>
double Histogram<T>::getBinWidth( void ) const
{
    return width;
}


template <class T>
void Histogram<T>::accumulate( T val )
{
    ++totalCounts;

    if ( val < min )
	++underflow;
    else if ( val >= max )
	++overflow;
    else {
	int bin = findBin(val);
	if (bin < 0 || bin >= numBins)
	    throw std::runtime_error("internal error " __FILE__ );
	++counts[bin];
    }
}

template <class T>
void Histogram<T>::reset( void )
{
    totalCounts = underflow = overflow = 0;

    for ( int i = 0; i < numBins; ++i )
	counts[i] = 0;

    validSort = 0;
}


template <class T>
int Histogram<T>::getCounts( int bin ) const
{
    if ( bin < 0 || bin >= numBins )
	throw std::runtime_error( "Histogram::getCounts: bin out of range" );
    return counts[bin]; 
}

template <class T>
int Histogram<T>::getCumulativeCounts( int bin ) const
{
    if ( bin < 0 || bin >= numBins )
	throw std::runtime_error( "Histogram::getCumulativeCounts: bin out of range" );

    int sum = underflow;
    for( int i = 0; i <= bin; ++i )
	sum += counts[i]; 

    return sum;
}

template <class T>
int Histogram<T>::getTotalCounts( void ) const
{ 
    return totalCounts; 
}

template <class T>
int Histogram<T>::getValidCounts( void ) const
{ 
    return getTotalCounts() - underflow - overflow; 
}

template <class T>
int Histogram<T>::getUnderflow( void ) const
{ 
    return underflow;
}

template <class T>
int Histogram<T>::getOverflow( void ) const
{ 
    return overflow; 
}


template <class T>
void Histogram<T>::dump( std::ostream& os ) const
{ 
    for ( int i = 0; i < numBins; ++i )
	os << getBinCentre(i) << ": " << getCounts(i) << std::endl;
    int u = getUnderflow();
    if (u) 
	os << "Underflow: " << u << std::endl;
    int o = getOverflow();
    if (o)
	os << "Overflow: " << o << std::endl;
}


template <class T>
double Histogram<T>::getBinLow( int i ) const
{
    return min + i * width;
}

template <class T>
double Histogram<T>::getBinCentre( int i ) const
{
    return min + (i + 0.5) * width;
}

template <class T>
double Histogram<T>::getBinHigh( int i ) const
{
    return min + (i + 1) * width;
}


template <class T>
int Histogram<T>::findBin( T x ) const
{
    return static_cast<int>( static_cast<double>(x - min)/width );
}

template <class T>
int Histogram<T>::findBinForCumulativeCounts( int c ) const
{
    int sum = underflow;
    for( int i = 0; i < numBins; ++i ) {
	sum += counts[i]; 
	if ( sum >= c )
	    return i;
    }
    throw std::runtime_error( "Histogram::findBinForCumulativeCounts: counts exceeds total in histogram" );
    return 0; // not reached, but don't want compilers to complain
}

template <class T>
int Histogram<T>::getSortedBin( int rank )
{
    if ( sortedIndex == 0 ) {
	sortedIndex = new int[numBins];
	validSort = 0;
    }

    if ( !validSort ) {
	for (int i = 0; i < numBins; ++i )
	    sortedIndex[i] = i;
	frequencySort( 0, numBins-1 );
	validSort = 1;
    }

    return sortedIndex[rank];
}


// Quicksort array sortedIndex[low_bin ... high_bin] (inclusive)
template <class T>
void Histogram<T>::frequencySort( int low_bin, int high_bin) 
{
    if ( low_bin >= high_bin )
	return;

    // Find pivot value
    int midval = counts[sortedIndex[(high_bin + low_bin) / 2]];

    // Partition array sortedIndex
    int l = low_bin, h = high_bin;
    while ( l <= h ) {
	if ( counts[sortedIndex[l]] <= midval )
	    ++l;
	else if ( counts[sortedIndex[h]] >= midval )
	    --h;
	else {
	    int x = sortedIndex[l];
	    sortedIndex[l] = sortedIndex[h];
	    sortedIndex[h] = x;
	    ++l, --h;
	}
    }

    frequencySort( low_bin, h );
    frequencySort( l, high_bin );
}


}


#endif
