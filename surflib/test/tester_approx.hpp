#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>


#include <list>
#include <string>


#include <CGAL/Cartesian.h>

#include <CGAL/circulator.h>
#include <CGAL/Point_3.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/IO/Verbose_ostream.h>


using std::string;

extern CGAL::Verbose_ostream verr;
extern double equality_tolerance;


namespace CppUnit {

    template <class T>
    struct cgal_assertion_traits
    {
	static string toString( const T& x )
	{
	    CppUnit::OStringStream ost;
	    CGAL::set_pretty_mode(ost);
	    ost << "[approx] " << x;
	    return ost.str();
	}
    };

    template<>    
    struct assertion_traits<double>
    {
	static inline bool equal( double a, double b )
	{
	    double diff = a - b;
	    if (fabs(diff) < equality_tolerance)
		return true;
	    std::cout << a << " - " << b << " = " << diff
		 << " >= " << equality_tolerance << std::endl;
	    return false;
	}

      static std::string toString( double x )
      {
          OStringStream ost;
          ost << "[approx] " << x;
          return ost.str();
      }
    };

    template <class R>
    struct assertion_traits<CGAL::Point_3<R> > 
	: public cgal_assertion_traits<CGAL::Point_3<R> >
    {
	static inline bool equal( const CGAL::Point_3<R>& p,
			   const CGAL::Point_3<R>& q )
	{
	    return CGAL::to_double(CGAL::squared_distance(p,q)) 
		< equality_tolerance*equality_tolerance;
	}
    };

    template <class R>
    struct assertion_traits<CGAL::Vector_3<R> > 
	: public cgal_assertion_traits<CGAL::Vector_3<R> >
    {
	static bool equal( const CGAL::Vector_3<R>& u,
			   const CGAL::Vector_3<R>& v )
	{
	    CGAL::Vector_3<R> uv = u - v;
	    double d2 = CGAL::to_double(uv*uv);
	    if ( d2 < equality_tolerance*equality_tolerance )
		return true;
	    std::cerr << "Magnitude difference: " << sqrt(d2) << std::endl;
	    return false;
	}
    };
    
}


