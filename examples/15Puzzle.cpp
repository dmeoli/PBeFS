/*--------------------------------------------------------------------------*/
/*--------------------------- File 15Puzzle.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Solves an instance of the 15-puzzle, toward the goal of Korf (1985)
 *
 *         1  2  3
 *      4  5  6  7
 *      8  9 10 11
 *     12 13 14 15
 *
 * Usage:   15Puzzle [-n] [-t <threads>] [<16 tiles>]
 *
 * where the 16 tiles are given row by row, 0 for the blank; the default
 * start is the instance 1 of Korf, whose optimal solution has 57 moves
 * (data/15-puzzle/korf100.txt holds all the 100). The search is A* with
 * the Manhattan distance, which gives an optimal solution; with -n the
 * heuristic is the sequence score of Nilsson, which may give a longer one;
 * with -t it is Hash Distributed A* with Abstract Zobrist hashing on the
 * given number of threads.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdlib>

#include <cstring>

#include <iostream>

#include <string>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/SequenceScore.h"

#include "PBeFS/astar/SlidingTile.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// prints the solution of a search

template< class S >
static int report( const S & sol )
{
 if( ! sol.solved() ) {
  std::cout << "no path to the goal" << std::endl;
  return( 1 );
  }
 int k = 0;
 for( auto s : sol.path ) {
  std::cout << "move " << k++ << std::endl;
  FifteenPuzzle::print( std::cout , s );
  std::cout << std::endl;
  }
 std::cout << "moves:     " << sol.cost << std::endl
           << "expanded:  " << sol.stats.expanded << std::endl
           << "seconds:   " << sol.stats.seconds << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 bool nilsson = false;
 unsigned threads = 0;
 std::string board;
 for( int i = 1 ; i < argc ; ++i )
  if( ! std::strcmp( argv[ i ] , "-n" ) )
   nilsson = true;
  else
   if( ! std::strcmp( argv[ i ] , "-t" ) && ( i + 1 < argc ) )
    threads = unsigned( std::atoi( argv[ ++i ] ) );
   else
    board += std::string( argv[ i ] ) + " ";
 if( board.empty() )
  board = "14 13 15 7 11 12 9 5 6 0 2 1 4 8 10 3";

 try {
  const auto start = FifteenPuzzle::parse( board );
  if( nilsson )
   return( report( astar( SequenceScore< 4 , 4 >( start ) ) ) );
  const FifteenPuzzle d( start );
  if( threads == 0 )
   return( report( astar( d ) ) );
  ParallelOptions o;
  o.threads = threads;
  const AbstractZobristDistribution< FifteenPuzzle > az( d );
  return( report( HDAStar< FifteenPuzzle , decltype( az ) >( d , az , o )
                  .solve() ) );
  }
 catch( std::exception & e ) {
  std::cerr << "15Puzzle: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File 15Puzzle.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
