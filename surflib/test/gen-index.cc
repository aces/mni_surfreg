/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <iostream>
#include <cstdlib>


int main( int ac, char* av[] )
{
    if ( ac != 2 )
	return 1;

    int n = std::atoi(av[1]);

    for( int i = 1; i <= n; ++i )
	std::cout << i << std::endl;

    return 0;
}
