/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#ifndef MNI_SURFLIB_DIAGNOSTIC_OSTREAM_HPP
#define MNI_SURFLIB_DIAGNOSTIC_OSTREAM_HPP

#include <iostream>


namespace MNI {


/*! \brief An output stream that filters according to priority level.
 *
 * This class acts as a filter for an existing IO stream, such as
 * std::cout or std::cerr.  Each <em>put to</em> operation
 * (<tt>operator\<\<</tt>) performed with this class has an associated
 * priority level.  A priority level is simply an integer, for which
 * lower values are more important than higher values.  The Diagnostic_ostream
 * will filter out put to operations which exceed a pre-specified 
 * maximum priority.
 * 
 * A typical use is for runtime-settable levels of output verbosity,
 * e.g. for debugging.
 *
 * \code
 *   Diagnostic_ostream diag(5,std::cout);
 *   diag.level(0) << "This will appear on cout." << endl;
 *   diag.level(5) << "So will this." << endl;
 *   diag.level(6) << "This message priority exceeds the maximum,"
 *                 << " and therefore is suppressed." << endl;
 * \endcode
 *
 * The higher the maximum level, the more verbose is the output.
 *
 */

// Inspired by CGAL's Verbose_ostream.

#define MNI_DO(x)  if ( msg_level <= max_level ) *os << x; return *this


class Diagnostic_ostream
{
    int msg_level, max_level;
    std::ostream* os;

public:
    /*! \param maximum_level Set the maximum level for subsequent
     *         put-to operations.
     *  \param out Reference to the underlying output stream.
     */
    Diagnostic_ostream( int maximum_level = 0,
			std::ostream& out = std::cerr )
	: msg_level(0), max_level(maximum_level), os(&out)
    {}

    //! Set maximum level of output.
    void maximum_level( int l ) { max_level = l; }

    //! Set level for following output.
    Diagnostic_ostream& level( int l ) 
    {
	msg_level = l;
	return *this;
    }

    template <class T>
    Diagnostic_ostream& operator<<( const T& t )
    { MNI_DO(t); }

    Diagnostic_ostream& operator<<( std::ostream& (*f)(std::ostream&) )
    { MNI_DO(f); }

    Diagnostic_ostream& operator<<( std::ios& (*f)(std::ios&) )
    { MNI_DO(f); }

#if needed_not_sure
    Diagnostic_ostream& flush()
    {
	if ( msg_level <= max_level )
	    os->flush();
	return *this;
    }
#endif

};
	

}

#undef MNI_DO

#endif
