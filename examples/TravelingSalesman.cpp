/*--------------------------------------------------------------------------*/
/*----------------------- File TravelingSalesman.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The traveling salesperson problem with A* and the minimum spanning tree
 * heuristic: exercise 3.30 of Russell and Norvig, "Artificial Intelligence:
 * A Modern Approach" (3rd edition, 2010).
 *
 * Usage:   TravelingSalesman [<cities> [<seed>]]
 *
 * The cities (default 10, at most 32) are random points of the unit square
 * drawn from the seed (default 1), and the tour starts and ends at the
 * first one. A state is the current city and the set of the visited ones;
 * the heuristic is the cost of a minimum spanning tree (by the algorithm of
 * Prim) of the unvisited cities together with the current and the first
 * one, which is admissible because the rest of the tour is a path through
 * all of them [see TravelingSalesman.h]. Up to 20 cities, the optimal tour
 * is also computed by the dynamic program of Held and Karp (1962), and the
 * program fails if the two costs differ.
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

#include <cmath>

#include <cstdlib>

#include <iostream>

#include <stdexcept>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/TravelingSalesman.h"

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
  const int n = argc > 1 ? std::atoi( argv[ 1 ] ) : 10;
  const unsigned long seed = argc > 2 ? std::strtoul( argv[ 2 ] , nullptr ,
                                                      10 ) : 1;
  const TravelingSalesman d( TravelingSalesman::random_cities( n , seed ) );
  const auto sol = astar( d );
  if( ! sol.solved() )
   throw( std::logic_error( "no tour found" ) );
  std::cout << "tour:";
  for( auto s : sol.path )
   std::cout << ' ' << TravelingSalesman::city( s );
  std::cout << std::endl << "cost:       " << sol.cost << std::endl
            << "expanded:   " << sol.stats.expanded << std::endl
            << "seconds:    " << sol.stats.seconds << std::endl;
  if( n > 20 )
   return( 0 );
  const double opt = d.held_karp();
  std::cout << "Held-Karp:  " << opt << std::endl;
  if( std::abs( sol.cost - opt ) > 1e-9 * ( 1 + opt ) ) {
   std::cerr << "TravelingSalesman: the tour is not optimal" << std::endl;
   return( 1 );
   }
  return( 0 );
  }
 catch( std::exception & e ) {
  std::cerr << "TravelingSalesman: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File TravelingSalesman.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
