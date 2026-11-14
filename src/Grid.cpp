/*--------------------------------------------------------------------------*/
/*----------------------------- File Grid.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the class GridMap: the construction of a map, free,
 * random or read from a benchmark of the Moving AI Lab, and the methods
 * that change it; and of the reading of a GridBoard.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <istream>

#include <random>

#include <stdexcept>

#include <string>

#include "PBeFS/astar/Grid.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- CONSTRUCTION OF A MAP --------------------------*/
/*--------------------------------------------------------------------------*/

GridMap::GridMap( int width , int height )
 : f_width( width ) , f_height( height ) ,
   v_cost( std::size_t( width ) * std::size_t( height ) , 1 )
{
 if( ( width < 0 ) || ( height < 0 ) )
  throw( std::invalid_argument( "GridMap::GridMap: negative size" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

GridMap GridMap::random( int width , int height , double p ,
                         std::uint64_t seed , int max_cost )
{
 if( ( max_cost < 1 ) || ( max_cost > 255 ) )
  throw( std::invalid_argument( "GridMap::random: max_cost is 1 to 255" ) );
 GridMap m( width , height );
 std::mt19937_64 gen( seed );
 for( auto & c : m.v_cost )
  if( random_unit( gen ) < p )
   c = 0;
  else
   c = std::uint8_t( max_cost == 1 ? 1 : 1 + random_below(
                                  gen , std::uint64_t( max_cost ) ) );
 m.set_cost( 0 , height - 1 , 1 );
 m.set_cost( width - 1 , height - 1 , 1 );
 return( m );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

GridMap GridMap::read_movingai( std::istream & in )
{
 std::string word;
 int width = -1 , height = -1;
 while( in >> word && ( word != "map" ) ) {
  if( word == "height" )
   in >> height;
  else
   if( word == "width" )
    in >> width;
  }
 if( ( width <= 0 ) || ( height <= 0 ) )
  throw( std::invalid_argument( "GridMap::read_movingai: bad header" ) );
 GridMap m( width , height );
 for( int y = 0 ; y < height ; ++y ) {
  std::string row;
  if( ! ( in >> row ) || ( int( row.size() ) != width ) )
   throw( std::invalid_argument( "GridMap::read_movingai: bad row" ) );
  for( int x = 0 ; x < width ; ++x ) {
   const char c = row[ x ];
   m.set_cost( x , y , ( c == '.' ) || ( c == 'G' ) || ( c == 'S' ) ? 1 : 0 );
   }
  }
 return( m );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

GridBoard GridBoard::read( std::istream & in , bool from_bottom )
{
 int height = -1 , width = -1;
 std::string word;
 if( ! ( in >> height >> width >> word ) || ( height <= 0 ) ||
     ( width <= 0 ) || ( word != "Board:" ) )
  throw( std::invalid_argument( "GridBoard::read: bad header" ) );
 std::getline( in , word );                     // the end of the line
 GridBoard b;
 b.map = GridMap( width , height );
 for( int y = 0 ; y < height ; ++y ) {
  std::string row;
  if( ! std::getline( in , row ) )
   throw( std::invalid_argument( "GridBoard::read: missing row" ) );
  for( int x = 0 ; ( x < width ) && ( x < int( row.size() ) ) ; ++x )
   if( row[ x ] == '#' )
    b.map.set_cost( x , y , 0 );
  }
 std::string cost , moves;
 if( ! ( in >> cost >> moves >> b.sx >> b.sy >> b.gx >> b.gy ) ||
     ( ( cost != "Unit" ) && ( cost != "Life" ) ) ||
     ( ( moves != "Four-way" ) && ( moves != "Eight-way" ) ) )
  throw( std::invalid_argument( "GridBoard::read: bad costs, moves or "
                                "points" ) );
 b.cost = cost == "Life" ? GridCost::life : GridCost::cells;
 b.moves = moves == "Four-way" ? 4 : 8;
 if( from_bottom ) {
  b.sy = height - 1 - b.sy;
  b.gy = height - 1 - b.gy;
  }
 if( ( b.sx < 0 ) || ( b.sy < 0 ) || ( b.gx < 0 ) || ( b.gy < 0 ) ||
     ( b.sx >= width ) || ( b.gx >= width ) || ( b.sy >= height ) ||
     ( b.gy >= height ) )
  throw( std::invalid_argument( "GridBoard::read: start or goal outside" ) );
 if( ! b.map.cost( b.sx , b.sy ) )
  b.map.set_cost( b.sx , b.sy , 1 );
 return( b );
 }

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHANGING AND READING THE MAP ----------------*/
/*--------------------------------------------------------------------------*/

void GridMap::set_cost( int x , int y , int c )
{
 if( ( x < 0 ) || ( y < 0 ) || ( x >= f_width ) || ( y >= f_height ) ||
     ( c < 0 ) || ( c > 255 ) )
  throw( std::invalid_argument( "GridMap::set_cost: bad cell or cost" ) );
 v_cost[ std::size_t( y ) * f_width + x ] = std::uint8_t( c );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

int GridMap::min_cost( void ) const noexcept
{
 int m = 256;
 for( auto c : v_cost )
  if( c && ( c < m ) )
   m = c;
 return( m == 256 ? 1 : m );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File Grid.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
