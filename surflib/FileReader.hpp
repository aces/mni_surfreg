#ifndef FILEREADER_H  // -*- C++ -*-
#define FILEREADER_H

#include <string>
#include <istream>
#include <stdexcept>


namespace MNI {


// Original code was written using FileReaderException.
typedef std::runtime_error FileReaderException;


class FileReader
{
    std::istream& in;
    bool skip_comments;
    std::string comment_prefix;

public:
    
    FileReader( std::istream& in_file ) : in( in_file )
    { 
	skip_comments = false; 
    }

    FileReader( std::istream& in_file, 
		const std::string& comment_prefix ) : in( in_file )
    {
	skip_comment_lines( comment_prefix );
    }

    void skip_comment_lines( const std::string& comment_prefix )
    {
	skip_comments = true;
	this->comment_prefix = comment_prefix;
    }
	

    void expect_word( std::string word ) 
    {
	std::string s;
	in >> s;
	if ( s != word )
	    throw FileReaderException( "expected [" + word
				       + "], got [" + s + "]" );
    }

    void expect_line( std::string line ) 
    {
	std::string s = expect_line();
	if ( s != line )
	    throw FileReaderException( "expected [" + line 
				       + "], got [" + s + "]" );
    }

    // FIXME: add error checking.
    std::string expect_word() 
    {
	std::string s;
	in >> s;
	return s;
    }

    // FIXME: add error checking.
    std::string expect_line()
    {
	std::string s;
	do {
	    std::getline( in, s );
	} while ( skip_comments && 
		  s.substr( 0, comment_prefix.length() ) == comment_prefix );
	return s;
    }

    // FIXME: add error checking.
    double expect_double()
    {
	double d;
	in >> d;
	return d;
    }

    // FIXME: add error checking.
    int expect_int()
    {
	int i;
	in >> i;
	return i;
    }


    // FIXME: add error checking.
    unsigned int expect_uint()
    {
	unsigned int i;
	in >> i;
	return i;
    }


};

}

#endif
