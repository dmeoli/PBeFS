/*--------------------------------------------------------------------------*/
/*--------------------------- File WaterJugs.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The water jugs problem with A*: exercise 3.6(d) of Russell and Norvig,
 * "Artificial Intelligence: A Modern Approach" (3rd edition, 2010), which
 * is the PourProblem of aima-python (github.com/aimacode/aima-python).
 *
 * Usage:   WaterJugs [<goal> <capacity> ...]
 *
 * Example: WaterJugs 4 5 3
 *
 * The jugs are empty at the start, and each move fills a jug from the
 * faucet, empties it on the ground, or pours it into another jug until
 * either the first is empty or the second is full; the goal is a jug
 * holding exactly <goal> units. The default is the exercise of the book:
 * jugs of 12, 8 and 3 gallons, and one gallon to measure. Every move costs
 * 1, and the heuristic is 0 at a goal and 1 elsewhere, which is
 * consistent; with up to 4 jugs of at most 65535 units.
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

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/WaterJugs.h"

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
  int goal = 1;
  std::vector< int > cap = { 12 , 8 , 3 };
  if( argc > 2 ) {
   goal = std::atoi( argv[ 1 ] );
   cap.clear();
   for( int i = 2 ; i < argc ; ++i )
    cap.push_back( std::atoi( argv[ i ] ) );
   }
  const WaterJugs d( cap , goal );
  const auto sol = astar( d );
  if( ! sol.solved() ) {
   std::cout << "the jugs cannot measure " << goal << std::endl;
   return( 1 );
   }
  for( auto s : sol.path ) {
   for( std::size_t j = 0 ; j < d.jugs() ; ++j )
    std::cout << ' ' << WaterJugs::level( s , j ) << '/' << cap[ j ];
   std::cout << std::endl;
   }
  std::cout << "moves:    " << sol.cost << std::endl
            << "expanded: " << sol.stats.expanded << std::endl
            << "seconds:  " << sol.stats.seconds << std::endl;
  return( 0 );
  }
 catch( std::exception & e ) {
  std::cerr << "WaterJugs: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File WaterJugs.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
