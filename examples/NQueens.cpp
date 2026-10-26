/*--------------------------------------------------------------------------*/
/*---------------------------- File NQueens.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The n-queens problem with A*, in the incremental formulation of section
 * 3.2.1 of Russell and Norvig, "Artificial Intelligence: A Modern
 * Approach" (3rd edition, 2010): a queen is added in the leftmost empty
 * column, in a row where no queen already on the board attacks it.
 *
 * Usage:   NQueens [<n>]
 *
 * The default is n = 8, and n is at most 32. Each queen costs 1, and the
 * heuristic, the number of the empty columns, is exact on the states that
 * lead to a solution; A* with the ties broken in favour of the larger g
 * then behaves as a depth-first search, and the counters tell how many
 * dead ends it meets [see NQueens.h].
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

#include <cstdlib>

#include <iostream>

#include <stdexcept>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/NQueens.h"

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
  const int n = argc > 1 ? std::atoi( argv[ 1 ] ) : 8;
  const auto sol = astar( NQueens( n ) );
  if( ! sol.solved() ) {
   std::cout << "no placement of " << n << " queens" << std::endl;
   return( 1 );
   }
  const auto s = sol.path.back();
  for( int r = 0 ; r < n ; ++r ) {
   for( int c = 0 ; c < n ; ++c )
    std::cout << ( s.row[ c ] == r ? " Q" : " ." );
   std::cout << std::endl;
   }
  // the placement checked apart from the Domain: n queens, no two on a
  // row or on a diagonal
  bool valid = s.k == n;
  for( int a = 0 ; valid && ( a < n ) ; ++a )
   for( int b = a + 1 ; b < n ; ++b ) {
    const int d = s.row[ b ] - s.row[ a ];
    if( ( d == 0 ) || ( d == b - a ) || ( d == a - b ) )
     valid = false;
    }
  std::cout << ( valid ? "a valid placement of " : "NOT a placement of " )
            << n << " queens" << std::endl
            << "expanded:  " << sol.stats.expanded << std::endl
            << "generated: " << sol.stats.generated << std::endl
            << "seconds:   " << sol.stats.seconds << std::endl;
  return( valid ? 0 : 1 );
  }
 catch( std::exception & e ) {
  std::cerr << "NQueens: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NQueens.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
