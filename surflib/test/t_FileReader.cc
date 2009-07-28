/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* This comment block is used for testing.  Do not modify.

word: expected
This is the line to expect.
double: 3.123
int: -34
uint: 32

*/


#include "tester.hpp"
#include <surflib/FileReader.hpp>
#include <fstream>


using namespace std;
using MNI::FileReader;



class t_FileReader : public CppUnit::TestCase
{
public:

    void check_read() 
    {
	std::ifstream in( SRCDIR "/t_FileReader.cc" );
	in.exceptions( std::ios::failbit | std::ios::badbit );
	FileReader fr(in);

	// Skip two lines
	fr.expect_line();
	fr.expect_line();

	// From here on, "expect_line()" is used to capture the trailing
	// newline after expect_word or whatever.

	fr.expect_word("word:");
	fr.expect_word("expected");
	fr.expect_line();

	fr.expect_line("This is the line to expect.");

	fr.expect_word("double:");
	CPPUNIT_ASSERT_EQUAL( 3.123, fr.expect_double() );
	fr.expect_line();
	
	fr.expect_word("int:");
	CPPUNIT_ASSERT_EQUAL( -34, fr.expect_int() );
	fr.expect_line();
	
	fr.expect_word("uint:");
	CPPUNIT_ASSERT_EQUAL( 32, fr.expect_int() );
	fr.expect_line();
    }


    CPPUNIT_TEST_SUITE( t_FileReader );
    CPPUNIT_TEST( check_read );
    CPPUNIT_TEST_SUITE_END();
};



CPPUNIT_TEST_SUITE_REGISTRATION( t_FileReader );
