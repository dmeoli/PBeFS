/*--------------------------------------------------------------------------*/
/*--------------------- File MissionariesCannibals.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The missionaries and cannibals problem with A*: exercise 3.9 of Russell
 * and Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition,
 * 2010), with m missionaries, c cannibals and a boat of b seats.
 *
 * Usage:   MissionariesCannibals [<m> <c> <b>]
 *
 * Example: MissionariesCannibals 4 4 2
 *
 * Everybody is on the left bank at the start, with the boat; a crossing
 * carries from 1 to b people, and on neither bank may the missionaries be
 * outnumbered by the cannibals, if any missionary is there. The default is
 * the exercise of the book, with 3 missionaries, 3 cannibals and 2 seats,
 * whose answer is 11 crossings; with 4, 4 and 2 there is none. Every
 * crossing costs 1, and the heuristic is the least number of crossings that
 * would carry the people left over the river with the boat always full,
 * which is admissible.
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

#include "PBeFS/astar/MissionariesCannibals.h"

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
  int m = 3 , c = 3 , b = 2;
  if( argc > 3 ) {
   m = std::atoi( argv[ 1 ] );
   c = std::atoi( argv[ 2 ] );
   b = std::atoi( argv[ 3 ] );
   }
  const MissionariesCannibals d( m , c , b );
  const auto sol = astar( d );
  if( ! sol.solved() ) {
   std::cout << "no crossing of " << m << " missionaries and " << c
             << " cannibals with " << b << " seats" << std::endl;
   return( 1 );
   }
  for( auto s : sol.path )
   std::cout << ' ' << MissionariesCannibals::missionaries( s ) << 'M'
             << MissionariesCannibals::cannibals( s ) << 'C'
             << ( MissionariesCannibals::boat_left( s ) ? " <boat" : "" )
             << "  |  " << m - MissionariesCannibals::missionaries( s )
             << 'M' << c - MissionariesCannibals::cannibals( s ) << 'C'
             << ( MissionariesCannibals::boat_left( s ) ? "" : " boat>" )
             << std::endl;
  std::cout << "crossings: " << sol.cost << std::endl
            << "expanded:  " << sol.stats.expanded << std::endl
            << "seconds:   " << sol.stats.seconds << std::endl;
  return( 0 );
  }
 catch( std::exception & e ) {
  std::cerr << "MissionariesCannibals: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File MissionariesCannibals.cpp -------------------*/
/*--------------------------------------------------------------------------*/
