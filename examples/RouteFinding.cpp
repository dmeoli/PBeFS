/*--------------------------------------------------------------------------*/
/*------------------------- File RouteFinding.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Route finding with A* on the road map of Romania of figure 3.2 of Russell
 * and Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition,
 * 2010).
 *
 * Usage:   RouteFinding [<from> [<to>]]
 *
 * Example: RouteFinding Lugoj Neamt
 *
 * The cities are Arad, Bucharest, Craiova, Drobeta, Eforie, Fagaras,
 * Giurgiu, Hirsova, Iasi, Lugoj, Mehadia, Neamt, Oradea, Pitesti, Rimnicu
 * (Vilcea), Sibiu, Timisoara, Urziceni, Vaslui and Zerind; the default is
 * the trip of the book, from Arad to Bucharest. Towards Bucharest the
 * heuristic is the straight-line distance of figure 3.22 of the book;
 * towards any other city it is the distance between the coordinates of the
 * cities in aima-python (github.com/aimacode/aima-python), which on this
 * map never exceeds the length of a road and is therefore consistent [see
 * Romania.h].
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

#include <iostream>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Romania.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 using City = Romania::City;
 const City from = argc > 1 ? Romania::find( argv[ 1 ] ) : Romania::Arad;
 const City to = argc > 2 ? Romania::find( argv[ 2 ] ) : Romania::Bucharest;
 if( ( from == Romania::NCities ) || ( to == Romania::NCities ) ) {
  std::cerr << "RouteFinding: unknown city" << std::endl;
  return( 2 );
  }

 const auto sol = astar( Romania( from , to ) );
 if( ! sol.solved() ) {
  std::cout << "no route" << std::endl;
  return( 1 );
  }
 double km = 0;
 for( std::size_t k = 0 ; k < sol.path.size() ; ++k ) {
  if( k > 0 )
   for( const Romania::Road & r : Romania::roads() )
    if( ( ( r.a == sol.path[ k - 1 ] ) && ( r.b == sol.path[ k ] ) ) ||
        ( ( r.b == sol.path[ k - 1 ] ) && ( r.a == sol.path[ k ] ) ) )
     km += r.km;
  std::cout << Romania::name( sol.path[ k ] ) << " (" << km << " km)"
            << std::endl;
  }
 std::cout << "expanded: " << sol.stats.expanded << std::endl
           << "seconds:  " << sol.stats.seconds << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File RouteFinding.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
