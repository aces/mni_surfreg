/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <iostream>

#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/CompilerOutputter.h>


// Used in tester_approx.h
double equality_tolerance = 1e-6;

using namespace std;




int main( int argc, char* argv[] )
{
    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry& registry
	= CppUnit::TestFactoryRegistry::getRegistry();

    runner.addTest( registry.makeTest() );
    runner.setOutputter( CppUnit::CompilerOutputter::defaultOutputter( &runner.result(), std::cout ));

    return ! runner.run( "", false );
}
