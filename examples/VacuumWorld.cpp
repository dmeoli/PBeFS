/*--------------------------------------------------------------------------*/
/*-------------------------- File VacuumWorld.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The vacuum world with A*: the toy problem of section 3.2.1 of Russell and
 * Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition, 2010),
 * on a row of n cells instead of 2.
 *
 * Usage:   VacuumWorld [<n> [<start> [<dirty cells>]]]
 *
 * Example: VacuumWorld 10 4 1011001101
 *
 * The robot starts in the cell <start> (0 by default), and the string of
 * 0s and 1s tells which cells are dirty (all by default); the actions Left,
 * Right and Suck cost 1 each, and the goal is a clean row. The default is
 * the world of the book, 2 cells both dirty with the robot on the left,
 * whose answer is Suck, Right, Suck. The heuristic, the number of dirty
 * cells, is consistent; with at most 32 cells.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>

#include <cstdlib>

#include <iostream>

#include <stdexcept>

#include <string>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/VacuumWorld.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 try {
  const int n = argc > 1 ? std::atoi( argv[ 1 ] ) : 2;
  const int start = argc > 2 ? std::atoi( argv[ 2 ] ) : 0;
  std::uint32_t dirty = n >= 32 ? ~std::uint32_t( 0 )
                                : ( std::uint32_t( 1 ) << n ) - 1;
  if( argc > 3 ) {
   const std::string d = argv[ 3 ];
   if( int( d.size() ) != n )
    throw( std::invalid_argument( "one 0 or 1 per cell" ) );
   dirty = 0;
   for( int i = 0 ; i < n ; ++i )
    if( d[ i ] == '1' )
     dirty |= std::uint32_t( 1 ) << i;
   }
  const VacuumWorld w( n , start , dirty );
  const auto sol = astar( w );
  for( std::size_t k = 0 ; k < sol.path.size() ; ++k ) {
   const auto s = sol.path[ k ];
   for( int i = 0 ; i < n ; ++i )
    std::cout << ( VacuumWorld::robot( s ) == i ? '[' : ' ' )
              << ( ( VacuumWorld::dirt( s ) >> i ) & 1 ? '*' : '.' )
              << ( VacuumWorld::robot( s ) == i ? ']' : ' ' );
   std::cout << std::endl;
   }
  std::cout << "actions:  " << sol.cost << std::endl
            << "expanded: " << sol.stats.expanded << std::endl
            << "seconds:  " << sol.stats.seconds << std::endl;
  return( 0 );
  }
 catch( std::exception & e ) {
  std::cerr << "VacuumWorld: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File VacuumWorld.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
