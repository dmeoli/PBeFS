/*--------------------------------------------------------------------------*/
/*------------------------------ File Grid.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class GridMap, a map of cells with a cost of
 * traversal, and for the class template Grid, which implements the Domain
 * concept [see Domain.h] for pathfinding on a GridMap with 4 or 8 moves.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_Grid
 #define __PBeFS_Grid
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <iosfwd>

#include <stdexcept>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS GridMap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a rectangular map of cells, each with a cost of traversal
/** The cell ( x , y ), with 0 <= x < width() and 0 <= y < height(), has
 * index y width() + x and cost cost( x , y ) in 1 , ... , 255, or 0 if it
 * is blocked. */

class GridMap {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor of a width x height map of free cells of cost 1

 GridMap( int width = 0 , int height = 0 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a random map with blocked cells of probability p
 /** Returns a width x height map in which each cell is blocked with
  * probability \p p, except the bottom-left and the bottom-right corners,
  * which are the start and the goal of the pathfinding instances of Burns
  * et al. (JAIR 2010) and cost 1. A free cell costs 1 if \p max_cost is 1,
  * as in those instances, and otherwise a cost drawn uniformly in 1 , ... ,
  * \p max_cost <= 255, the cells being drawn row by row. */

 [[nodiscard]] static GridMap random( int width , int height , double p ,
                                      std::uint64_t seed ,
                                      int max_cost = 1 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads a map in the format of the benchmarks of Sturtevant (2012)
 /** Reads from \p in a map of the Moving AI Lab benchmarks: the header
  * "type octile", "height H", "width W", "map", then H lines of W
  * characters, '.', 'G' and 'S' being free cells of cost 1, any other one
  * being blocked. Throws std::invalid_argument if the map is malformed. */

 [[nodiscard]] static GridMap read_movingai( std::istream & in );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of columns

 [[nodiscard]] int width( void ) const noexcept { return( f_width ); }

 /// the number of rows

 [[nodiscard]] int height( void ) const noexcept { return( f_height ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost of the cell ( x , y ), 0 if blocked or outside the map

 [[nodiscard]] int cost( int x , int y ) const noexcept {
  if( ( x < 0 ) || ( y < 0 ) || ( x >= f_width ) || ( y >= f_height ) )
   return( 0 );
  return( v_cost[ std::size_t( y ) * f_width + x ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the cost of the cell ( x , y ) to c in 0 , ... , 255

 void set_cost( int x , int y , int c );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the smallest cost of a free cell (1 if there is none)

 [[nodiscard]] int min_cost( void ) const noexcept;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 int f_width;                                    ///< the number of columns
 int f_height;                                   ///< the number of rows
 std::vector< std::uint8_t > v_cost;             ///< the costs, row by row

 };  // end( class( GridMap ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS Grid -------------------------------*/
/*--------------------------------------------------------------------------*/
/// pathfinding on a GridMap with 4 (Moves == 4) or 8 (Moves == 8) moves
/** A state is the index y W + x of a cell. With 4 moves the costs are
 * integers and leaving a cell costs its cost; with 8 moves a diagonal move
 * costs sqrt( 2 ) times the cost of the cell left, and it is allowed only
 * if both the cells it cuts are free. The heuristic is the Manhattan
 * distance (4 moves) or the octile distance (8 moves) times the smallest
 * cost of a free cell.
 *
 * The features of Zobrist hashing are the column and the row of the cell;
 * the abstract features of Abstract Zobrist hashing are the column and the
 * row of the block of block() x block() cells containing it, and the
 * abstract state of AHDA* is that block, which for this domain coincide
 * [see Jinnai and Fukunaga, JAIR 2017, section 4.1]. The nblocks of PBNF
 * are the blocks of nblock_side() x nblock_side() cells, and the perfect
 * hash is the index of the cell. */

template< int Moves = 4 >
class Grid {

 static_assert( ( Moves == 4 ) || ( Moves == 8 ) ,
                "Grid: the moves are either 4 or 8" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint32_t;                     ///< a cell index
 using cost_type = std::conditional_t< Moves == 4 , int , double >;
                                                        ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the map (which must outlive the Grid) and the cells
 /** Constructs the instance of going from ( sx , sy ) to ( gx , gy ) on
  * \p map; \p block is the side of the blocks of AZH and AHDA*, and \p
  * nblock_side that of the nblocks of PBNF. Throws std::invalid_argument
  * if the start or the goal is blocked. */

 Grid( const GridMap & map , int sx , int sy , int gx , int gy ,
       int block = 100 , int nblock_side = 20 )
  : f_map( map ) , f_w( map.width() ) , f_h( map.height() ) ,
    f_start( index( sx , sy ) ) , f_goal( index( gx , gy ) ) ,
    f_gx( gx ) , f_gy( gy ) , f_min( map.min_cost() ) ,
    f_block( block ) , f_nside( nblock_side ) {
  if( ! map.cost( sx , sy ) || ! map.cost( gx , gy ) )
   throw( std::invalid_argument( "Grid::Grid: blocked start or goal" ) );
  if( ( block <= 0 ) || ( nblock_side <= 0 ) )
   throw( std::invalid_argument( "Grid::Grid: blocks of side <= 0" ) );
  f_bw = ( f_w + f_block - 1 ) / f_block;
  f_bh = ( f_h + f_block - 1 ) / f_block;
  f_nw = ( f_w + f_nside - 1 ) / f_nside;
  f_nh = ( f_h + f_nside - 1 ) / f_nside;
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/
 /// the start cell

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 /// true if s is the goal cell

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_goal );
  }

 /// the hash of a cell: its index, mixed

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( mix64( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan (4 moves) or octile (8 moves) distance to the goal

 [[nodiscard]] cost_type heuristic( state_type s ) const noexcept {
  return( distance( int( s % f_w ) , int( s / f_w ) , f_gx , f_gy ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) ) for each free cell t next to s

 template< class F >
 void for_each_successor( state_type s , cost_type , F && f ) const {
  const int x = int( s % f_w ) , y = int( s / f_w );
  const cost_type c = cost_type( f_map.cost( x , y ) );
  static constexpr int dx[ 8 ] = { 0 , -1 , 1 , 0 , -1 , 1 , -1 , 1 };
  static constexpr int dy[ 8 ] = { -1 , 0 , 0 , 1 , -1 , -1 , 1 , 1 };
  for( int k = 0 ; k < Moves ; ++k ) {
   const int nx = x + dx[ k ] , ny = y + dy[ k ];
   if( ! f_map.cost( nx , ny ) )
    continue;
   if( k >= 4 ) {                     // a diagonal: no corner cutting
    if( ! f_map.cost( nx , y ) || ! f_map.cost( x , ny ) )
     continue;
    f( index( nx , ny ) , c * c_sqrt2 , distance( nx , ny , f_gx , f_gy ) );
    }
   else
    f( index( nx , ny ) , c , distance( nx , ny , f_gx , f_gy ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the columns and the rows

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( std::size_t( f_w ) + f_h );
  }

 /// calls f( x ) and f( W + y ) for the cell ( x , y ) = s

 template< class F >
 void for_each_feature( state_type s , F && f ) const {
  f( std::size_t( s % f_w ) );
  f( std::size_t( f_w ) + s / f_w );
  }

 /// the number of abstract features: the columns and rows of blocks

 [[nodiscard]] std::size_t num_abstract_features( void ) const noexcept {
  return( std::size_t( f_bw ) + f_bh );
  }

 /// the abstract feature of a feature: the column or row of its block

 [[nodiscard]] std::size_t abstract_feature( std::size_t i ) const noexcept
 {
  if( i < std::size_t( f_w ) )
   return( i / f_block );
  return( std::size_t( f_bw ) + ( i - f_w ) / f_block );
  }

/*--------------------------- STATE ABSTRACTIONS ---------------------------*/
 /// the abstract state of AHDA*: the block of the cell

 [[nodiscard]] std::uint64_t abstract_state( state_type s ) const noexcept {
  return( std::uint64_t( s / f_w / f_block ) * f_bw + s % f_w / f_block );
  }

 /// the number of nblocks

 [[nodiscard]] std::size_t num_nblocks( void ) const noexcept {
  return( std::size_t( f_nw ) * f_nh );
  }

 /// the nblock of PBNF: the block of nblock_side() x nblock_side() cells

 [[nodiscard]] std::size_t nblock( state_type s ) const noexcept {
  return( std::size_t( s / f_w / f_nside ) * f_nw + s % f_w / f_nside );
  }

 /// calls f( b' ) for each nblock b' next to the nblock b

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const {
  const int bx = int( b % f_nw ) , by = int( b / f_nw );
  for( int dy = -1 ; dy <= 1 ; ++dy )
   for( int dx = -1 ; dx <= 1 ; ++dx ) {
    if( ( ! dx && ! dy ) || ( ( Moves == 4 ) && dx && dy ) )
     continue;
    const int nx = bx + dx , ny = by + dy;
    if( ( nx >= 0 ) && ( ny >= 0 ) && ( nx < f_nw ) && ( ny < f_nh ) )
     f( std::size_t( ny ) * f_nw + nx );
    }
  }

 /// the side of the nblocks

 [[nodiscard]] int nblock_side( void ) const noexcept { return( f_nside ); }

 /// the side of the blocks of AZH and AHDA*

 [[nodiscard]] int block( void ) const noexcept { return( f_block ); }

/*---------------------- OTHER HEURISTICS AND HASHES -----------------------*/
 /// the heuristic of PA*SE: the distance between two cells

 [[nodiscard]] cost_type pairwise_heuristic( state_type s , state_type t )
  const noexcept {
  return( distance( int( s % f_w ) , int( s / f_w ) ,
                    int( t % f_w ) , int( t / f_w ) ) );
  }

 /// the perfect hash of HDA*[P]: the cell index

 [[nodiscard]] std::uint64_t perfect_hash( state_type s ) const noexcept {
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cell index of ( x , y )

 [[nodiscard]] state_type index( int x , int y ) const noexcept {
  return( state_type( y ) * state_type( f_w ) + state_type( x ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 static constexpr double c_sqrt2 = 1.41421356237309504880;   ///< sqrt( 2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the Manhattan or octile distance times the smallest cost

 [[nodiscard]] cost_type distance( int x , int y , int u , int v )
  const noexcept {
  const int ax = std::abs( x - u ) , ay = std::abs( y - v );
  if constexpr( Moves == 4 )
   return( ( ax + ay ) * f_min );
  else
   return( ( std::max( ax , ay ) + ( c_sqrt2 - 1 ) * std::min( ax , ay ) )
           * f_min );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const GridMap & f_map;                           ///< the map
 int f_w;                                         ///< its width
 int f_h;                                         ///< its height
 state_type f_start;                              ///< the start cell
 state_type f_goal;                               ///< the goal cell
 int f_gx;                                        ///< the goal column
 int f_gy;                                        ///< the goal row
 int f_min;                                  ///< the smallest cell cost
 int f_block;                          ///< the side of the AZH blocks
 int f_nside;                          ///< the side of the nblocks
 int f_bw;                             ///< the number of block columns
 int f_bh;                             ///< the number of block rows
 int f_nw;                             ///< the number of nblock columns
 int f_nh;                             ///< the number of nblock rows

 };  // end( class( Grid ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Grid */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Grid.h -----------------------------*/
/*--------------------------------------------------------------------------*/
