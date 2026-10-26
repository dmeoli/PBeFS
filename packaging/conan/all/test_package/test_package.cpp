#include "PBeFS/PBeFS.h"

#include <iostream>

using namespace PBeFS_di_unipi_it;

int main()
{
 const EightPuzzle d( EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" ) );
 const auto s = astar( d );
 std::cout << s.cost << '\n';
 return( s.cost == 26 ? 0 : 1 );
 }
