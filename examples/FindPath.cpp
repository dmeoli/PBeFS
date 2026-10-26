/*--------------------------------------------------------------------------*/
/*--------------------------- File FindPath.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Pathfinding with A* on a grid whose cells cost what one pays to leave
 * them, and draws the path on the map.
 *
 * Usage:   FindPath [-8] [-f <map>] [-r <W>x<H>:<p>:<seed>]
 *                   [<sx> <sy> <gx> <gy>]
 *
 * Example: FindPath -r 40x20:0.25:7
 *
 * The map is read from a file in the format of the Moving AI benchmarks
 * (-f), whose free cells cost 1, or it is drawn at random (-r, default
 * 30x15:0.2:1): each cell is a wall with probability p and otherwise costs
 * 1 , ... , 5 to leave, uniformly. The moves are 4, or 8 with -8, in which
 * case a diagonal move costs sqrt( 2 ) times the cell left. The start and
 * the goal default to the bottom-left and the bottom-right cells, which are
 * then kept free. In the drawing, '#' is a wall, a digit the cost of a
 * cell, '*' the path, 'S' the start and 'G' the goal.
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

#include <cstdio>

#include <cstdlib>

#include <cstring>

#include <fstream>

#include <iostream>

#include <string>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Grid.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// solves the grid with Moves moves and draws the path on the map

template< int Moves >
static int solve( GridMap m , int sx , int sy , int gx , int gy )
{
 const Grid< Moves > g( m , sx , sy , gx , gy );
 const auto sol = astar( g );
 if( ! sol.solved() ) {
  std::cout << "no path" << std::endl;
  return( 1 );
  }
 const int w = m.width();
 std::vector< std::string > rows( m.height() , std::string( w , ' ' ) );
 for( int y = 0 ; y < m.height() ; ++y )
  for( int x = 0 ; x < w ; ++x )
   rows[ y ][ x ] = m.cost( x , y ) ? char( '0' + m.cost( x , y ) % 10 ) : '#';
 for( auto s : sol.path )
  rows[ s / w ][ s % w ] = '*';
 rows[ sy ][ sx ] = 'S';
 rows[ gy ][ gx ] = 'G';
 for( const auto & r : rows )
  std::cout << r << std::endl;
 std::cout << "moves:    " << sol.path.size() - 1 << std::endl
           << "cost:     " << sol.cost << std::endl
           << "expanded: " << sol.stats.expanded << std::endl
           << "seconds:  " << sol.stats.seconds << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 try {
  bool eight = false;
  std::string file , random = "30x15:0.2:1";
  std::vector< int > pts;
  for( int i = 1 ; i < argc ; ++i )
   if( ! std::strcmp( argv[ i ] , "-8" ) )
    eight = true;
   else
    if( ! std::strcmp( argv[ i ] , "-f" ) && ( i + 1 < argc ) )
     file = argv[ ++i ];
    else
     if( ! std::strcmp( argv[ i ] , "-r" ) && ( i + 1 < argc ) )
      random = argv[ ++i ];
     else
      pts.push_back( std::atoi( argv[ i ] ) );
  if( ! ( pts.empty() || ( pts.size() == 4 ) ) )
   throw( std::invalid_argument( "the start and the goal need 4 numbers" ) );

  GridMap m;
  if( ! file.empty() ) {
   std::ifstream in( file );
   if( ! in )
    throw( std::invalid_argument( "cannot read " + file ) );
   m = GridMap::read_movingai( in );
   }
  else {
   int w , h;
   double p;
   unsigned long long seed;
   if( std::sscanf( random.c_str() , "%dx%d:%lf:%llu" , & w , & h , & p ,
                    & seed ) != 4 )
    throw( std::invalid_argument( "-r wants <W>x<H>:<p>:<seed>" ) );
   m = GridMap::random( w , h , p , seed , 5 );
   }

  int sx = 0 , sy = m.height() - 1 , gx = m.width() - 1 , gy = sy;
  if( pts.empty() ) {
   m.set_cost( sx , sy , 1 );
   m.set_cost( gx , gy , 1 );
   }
  else {
   sx = pts[ 0 ];
   sy = pts[ 1 ];
   gx = pts[ 2 ];
   gy = pts[ 3 ];
   }
  return( eight ? solve< 8 >( m , sx , sy , gx , gy )
                : solve< 4 >( m , sx , sy , gx , gy ) );
  }
 catch( std::exception & e ) {
  std::cerr << "FindPath: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File FindPath.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
