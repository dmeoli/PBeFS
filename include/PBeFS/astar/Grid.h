/*--------------------------------------------------------------------------*/
/*------------------------------ File Grid.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class GridMap, a map of cells with a cost of
 * traversal, for the struct GridBoard, a pathfinding instance in the
 * format of the grids of Burns et al. (JAIR 2010), for the class template
 * GridAdjacency, the free cells next to a cell of a GridMap, and for the
 * class template Grid, which implements the Domain concept [see Domain.h]
 * for pathfinding on a GridMap with 4 or 8 moves. The adjacency is that of
 * the set of data structures [see Structures.h]: the one of the library,
 * computed from the map, for the sets own, std and boost (neither the
 * standard library nor Boost has a graph of a grid with blocked cells and
 * diagonal moves, and a graph in compressed rows of the largest maps would
 * take tens of GB), and ADS::GridGraph for the set ads.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
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

#include <optional>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/Structures.h"

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
/** The cell \f$ ( x , y ) \f$, with \f$ 0 \leq x < W \f$ and
 * \f$ 0 \leq y < H \f$ for the width \f$ W \f$ and the height \f$ H \f$
 * of the map, has the index \f$ y W + x \f$ and the cost cost( x , y ) in
 * \f$ \{ 1 , \ldots , 255 \} \f$, or 0 if it is blocked. */

class GridMap {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor of a map of width by height free cells of cost 1

 GridMap( int width = 0 , int height = 0 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a random map with blocked cells of probability p
 /** Returns a map of \f$ \mathrm{width} \times \mathrm{height} \f$
  * cells in which each cell is blocked with
  * probability \p p, except the bottom-left and the bottom-right corners,
  * which are the start and the goal of the pathfinding instances of Burns
  * et al. (JAIR 2010) and cost 1. A free cell costs 1 if \p max_cost is 1,
  * as in those instances, and otherwise a cost drawn uniformly in
  * \f$ \{ 1 , \ldots , c_{\max} \} \f$, \f$ c_{\max} \leq 255 \f$ being
  * \p max_cost, the cells being drawn row by row. */

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
 /// sets the cost of the cell ( x , y ) to c, in
 /// \f$ \{ 0 , \ldots , 255 \} \f$

 void set_cost( int x , int y , int c );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the smallest cost of a free cell (1 if there is none)

 [[nodiscard]] int min_cost( void ) const noexcept;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the costs of the cells, row by row

 [[nodiscard]] const std::vector< std::uint8_t > & costs( void ) const
  noexcept {
  return( v_cost );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 int f_width;                                    ///< the number of columns
 int f_height;                                   ///< the number of rows
 std::vector< std::uint8_t > v_cost;             ///< the costs, row by row

 };  // end( class( GridMap ) )

/*--------------------------------------------------------------------------*/
/*-------------------------------- GridCost --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the cost of leaving a cell of a Grid
/** With cells, leaving the cell \f$ ( x , y ) \f$ costs the cost of the
 * cell in its GridMap; with life, it costs the row \f$ y \f$, whatever the
 * cell, which makes the rows of small \f$ y \f$ cheaper (the first one
 * free) and the shortest path longer than the cheapest one [see Burns et
 * al., JAIR 2010]. A move to a corner costs \f$ \sqrt{2} \f$ times as
 * much in both cases. */

//! [domains-sig-gridboard]
enum class GridCost {
 cells ,                                 ///< the cost of the cell left
 life                                    ///< the row of the cell left
 };

/*--------------------------------------------------------------------------*/
/*---------------------------- STRUCT GridBoard ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a pathfinding instance in the format of the grids of Burns et al.
/** The format of the instances of Burns et al. (JAIR 2010) and of their
 * generators: a line "H W" with the height and the width, a line
 * "Board:", H lines of W characters, '#' for a blocked cell and any other
 * one for a free cell, a line "Unit" or "Life" (the GridCost), a line
 * "Four-way" or "Eight-way" (the moves, a diagonal move being allowed even
 * if it cuts a blocked cell) and the four numbers \f$ x_s \f$,
 * \f$ y_s \f$, \f$ x_g \f$ and \f$ y_g \f$ of the start and the goal.
 * The files of the authors count \f$ y \f$ from the last line of the
 * board up, those of later generators from the first line down: read()
 * takes the convention, and stores the cells with \f$ y \f$ counted from
 * the first line, as GridMap does. */

struct GridBoard {

 GridMap map;                                    ///< the free cells
 int sx = 0;                                     ///< the start column
 int sy = 0;                                     ///< the start row
 int gx = 0;                                     ///< the goal column
 int gy = 0;                                     ///< the goal row
 int moves = 4;                                  ///< 4 or 8
 GridCost cost = GridCost::cells;                ///< unit or life costs

 /// reads a board, \p from_bottom telling whether \f$ y \f$ is counted
 /// from the last line up
 /** Throws std::invalid_argument if the board is malformed or the start
  * or the goal is outside it. A blocked start is made free, which changes
  * no path, since a shortest one never comes back to its start: some
  * instances of the authors have one. */

 [[nodiscard]] static GridBoard read( std::istream & in ,
                                      bool from_bottom = true );

 };  // end( struct( GridBoard ) )
//! [domains-sig-gridboard]

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// \f$ \sqrt{2} \f$, the cost of a diagonal move per unit of cost of a cell

inline constexpr double grid_sqrt2 = 1.41421356237309504880;

/*--------------------------------------------------------------------------*/

#if defined( PBEFS_STRUCTURES_STD ) || defined( PBEFS_STRUCTURES_BOOST )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS GridAdjacency ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the moves from a cell of a GridMap, with 4 or 8 moves, and their costs
/** A move leaves a cell at its cost, times \f$ \sqrt{2} \f$ to a corner,
 * and a move to a corner is allowed only if both the cells it cuts are
 * free, unless the corners may be cut; the moves are tried up, left, right
 * and down, then to the corners. */

template< int Moves >
class GridAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// the cost of a move: an integer with 4 moves, a real with 8
 using cost_type = std::conditional_t< Moves == 4 , int , double >;

 /// constructor, taking the map, which must outlive the GridAdjacency,
 /// and whether a diagonal move may cut a blocked cell

 explicit GridAdjacency( const GridMap & map , bool cut_corners = false )
  : f_map( map ) , f_cut( cut_corners ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , x' , y' , c ) for each free cell t = ( x' , y' ) next to
 /// the free cell s = ( x , y ), c being the cost of the move

 template< class F >
 void for_each( std::uint32_t , int x , int y , F && f ) const {
  static constexpr int dx[ 8 ] = { 0 , -1 , 1 , 0 , -1 , 1 , -1 , 1 };
  static constexpr int dy[ 8 ] = { -1 , 0 , 0 , 1 , -1 , -1 , 1 , 1 };
  const int w = f_map.width();
  const cost_type c = cost_type( f_map.cost( x , y ) );
  for( int k = 0 ; k < Moves ; ++k ) {
   const int nx = x + dx[ k ] , ny = y + dy[ k ];
   if( ! f_map.cost( nx , ny ) )
    continue;
   if( ( k >= 4 ) && ! f_cut &&       // a diagonal: no corner cutting
       ( ! f_map.cost( nx , y ) || ! f_map.cost( x , ny ) ) )
    continue;
   f( std::uint32_t( ny ) * std::uint32_t( w ) + std::uint32_t( nx ) ,
      nx , ny , k >= 4 ? cost_type( c * grid_sqrt2 ) : c );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cell ( x , y ) of the map changed: nothing to do, the moves are
 /// read from the map

 void cell_changed( int , int ) noexcept {}

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const GridMap & f_map;                                         ///< the map
 bool f_cut;                                   ///< corners may be cut

 };  // end( class( GridAdjacency ) )

#endif  /* PBEFS_STRUCTURES_STD and PBEFS_STRUCTURES_BOOST */

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*------------------ THE ADJACENCY OF THE OTHER STRUCTURES -----------------*/
/*--------------------------------------------------------------------------*/

#if defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/GridAdjacency.h"
#endif

/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS GridNBlocks ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the shape of the nblocks of a Grid: rectangles of width times height
/// cells, squares by default

struct GridNBlocks {

 int width;          ///< the columns of an nblock
 int height;         ///< the rows of an nblock

 /// the squares of side cells
 constexpr GridNBlocks( int side = 20 ) noexcept
  : width( side ) , height( side ) {}

 /// the rectangles of w times h cells
 constexpr GridNBlocks( int w , int h ) noexcept : width( w ) , height( h ) {}

 /// the nblocks for cutting a map of w times h cells into cols times rows
 /// rectangles, after the coarse abstraction of the grids of Burns et al.
 /// (JAIR 2010): cols (rows) is raised by one if it does not divide w (h),
 /// and the nblock is \f$ \lfloor w / cols \rfloor \f$ cells wide and
 /// \f$ \lfloor h / rows \rfloor \f$ high, at least 1; the map then has
 /// \f$ \lceil w / \lfloor w / cols \rfloor \rceil \f$ columns of nblocks,
 /// the last one possibly narrower, which may exceed cols by more than one
 /// (and likewise for the rows)
 [[nodiscard]] static constexpr GridNBlocks coarse( int w , int h , int cols ,
                                                     int rows ) noexcept {
  cols = cols > 0 ? cols : 1;
  rows = rows > 0 ? rows : 1;
  cols += ( w % cols ) != 0;
  rows += ( h % rows ) != 0;
  return( GridNBlocks( w / cols > 0 ? w / cols : 1 ,
                       h / rows > 0 ? h / rows : 1 ) );
  }

 };  // end( struct( GridNBlocks ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS Grid -------------------------------*/
/*--------------------------------------------------------------------------*/
/// pathfinding on a GridMap with 4 (Moves == 4) or 8 (Moves == 8) moves
/** A state is the index \f$ y W + x \f$ of a cell. With 4 moves the costs are
 * integers and leaving a cell costs its cost; with 8 moves a diagonal move
 * costs \f$ \sqrt{2} \f$ times the cost of the cell left, and it is allowed
 * only if both the cells it cuts are free, unless the corners may be cut
 * [see GridAdjacency]. The heuristic is the Manhattan distance (4 moves)
 * or the octile distance (8 moves) times the smallest cost of a free cell.
 *
 * With the life costs of GridCost, leaving a cell of the row \f$ y \f$
 * costs \f$ y \f$ (times \f$ \sqrt{2} \f$ to a corner), and the
 * heuristic is the cost of a cheapest path on the board without blocked
 * cells, with the moves to a corner costing as the straight ones when
 * there are 8 moves. Let \f$ V( a , b ) \f$ be the cost of going
 * straight from the row \f$ a \f$ to the row \f$ b \f$, i.e.,
 * \f$ \sum_{k = b + 1}^{a} k \f$ upwards and
 * \f$ \sum_{k = a}^{b - 1} k \f$ downwards. With 4 moves, a path from
 * \f$ ( x , y ) \f$ to \f$ ( u , v ) \f$ whose highest row is
 * \f$ r \leq \min\{ y , v \} \f$ costs at least
 * \f$ V( y , r ) + r | x - u | + V( r , v ) \f$, which is
 * \f$ r | x - u | - r^2 \f$ plus a constant, a concave function of
 * \f$ r \f$: the heuristic is its minimum, at \f$ r = 0 \f$ or at
 * \f$ r = \min\{ y , v \} \f$. With 8 moves, a path of \f$ n \f$
 * moves leaves its \f$ k \f$-th cell from a row no smaller than
 * \f$ \max\{ y - k , v - n + k , 0 \} \f$, and the sum of these bounds
 * does not decrease with \f$ n \f$ (a turn of the rows can be cut out),
 * so the heuristic is that sum for the fewest moves
 * \f$ n = \max\{ | x - u | , | y - v | \} \f$. Both are the exact
 * costs of a relaxation, and therefore consistent.
 *
 * The heuristic of the backward search is the same estimate of the cost
 * from the start, the predecessors of a cell being its free neighbors.
 *
 * The features of Zobrist hashing are the column and the row of the cell;
 * the abstract features of Abstract Zobrist hashing are the column and the
 * row of the block of \f$ b \times b \f$ cells containing it, \f$ b \f$
 * being block(), and the abstract state of AHDA* is that block, which for
 * this domain coincide [see Jinnai and Fukunaga, JAIR 2017, section 4.1].
 * The nblocks of PBNF are the rectangles of nblock_width() times
 * nblock_height() cells [see GridNBlocks], squares by default, and the
 * perfect hash is the index of the cell. The preferred moves are those
 * that lower the heuristic (the distance to the goal, or the cost of the
 * relaxation with the life costs), and the goal features its column and
 * its row [see PreferredDomain and GoalFeatureDomain].
 *
 * The operator of a move by \f$ ( d_x , d_y ) \f$ is \f$ 3 ( d_y + 1 ) +
 * d_x + 1 \f$, whose inverse is 8 minus it [see OperatorDomain]; going back
 * to the parent is never cheaper even when the cost of a move is that of
 * the cell left, or its row with the life costs, since the costs are not
 * negative. The operator selection function of EPEA* computes
 * \f$ \Delta f = c + h( t ) - h( s ) \f$ of each move without passing on
 * the successors of another \f$ \Delta f \f$, and with 8 moves it takes
 * two values within cost_slack() of each other for the same [see
 * PartialExpansionDomain]. */

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
  * nblocks the shape of the nblocks of PBNF (a side, for squares). Throws
  * std::invalid_argument if the start or the goal is blocked. */

 Grid( const GridMap & map , int sx , int sy , int gx , int gy ,
       int block = 100 , GridNBlocks nblocks = {} ,
       GridCost cost = GridCost::cells , bool cut_corners = false )
  : f_map( map ) , f_adj( std::in_place , map , cut_corners ) ,
    f_w( map.width() ) , f_h( map.height() ) , f_start( index( sx , sy ) ) ,
    f_goal( index( gx , gy ) ) , f_sx( sx ) , f_sy( sy ) , f_gx( gx ) ,
    f_gy( gy ) , f_min( map.min_cost() ) , f_block( block ) ,
    f_nsw( nblocks.width ) , f_nsh( nblocks.height ) ,
    f_life( cost == GridCost::life ) , f_cut( cut_corners ) {
  if( ! map.cost( sx , sy ) || ! map.cost( gx , gy ) )
   throw( std::invalid_argument( "Grid::Grid: blocked start or goal" ) );
  if( ( block <= 0 ) || ( f_nsw <= 0 ) || ( f_nsh <= 0 ) )
   throw( std::invalid_argument( "Grid::Grid: blocks of side <= 0" ) );
  f_bw = ( f_w + f_block - 1 ) / f_block;
  f_bh = ( f_h + f_block - 1 ) / f_block;
  f_nw = ( f_w + f_nsw - 1 ) / f_nsw;
  f_nh = ( f_h + f_nsh - 1 ) / f_nsh;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, taking a GridBoard (which must outlive the Grid)
 /** The costs are those of the board, and a diagonal move may cut a
  * blocked cell. Throws std::invalid_argument if the board has not Moves
  * moves, or as the other constructor. */

 explicit Grid( const GridBoard & b , int block = 100 ,
                GridNBlocks nblocks = {} )
  : Grid( b.map , b.sx , b.sy , b.gx , b.gy , block , nblocks ,
          b.cost , true ) {
  if( b.moves != Moves )
   throw( std::invalid_argument( "Grid::Grid: the board has other moves" ) );
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
  return( ADS::splitmix64( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan (4 moves) or octile (8 moves) distance to the goal

 [[nodiscard]] cost_type heuristic( state_type s ) const noexcept {
  return( distance( int( s % f_w ) , int( s / f_w ) , f_gx , f_gy ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) ) for each free cell t next to s, none if s
 /// is blocked (which it may become after a refresh())

 template< class F >
 void for_each_successor( state_type s , cost_type , F && f ) const {
  const int x = int( s % f_w ) , y = int( s / f_w );
  if( ! f_map.cost( x , y ) )
   return;
  f_adj->for_each( s , x , y , [ & ]( state_type t , int nx , int ny ,
                                      cost_type c ) {
   if( f_life )
    c = ( nx != x ) && ( ny != y ) ? cost_type( y * c_sqrt2 )
                                   : cost_type( y );
   f( t , c , distance( nx , ny , f_gx , f_gy ) );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) , op_t ) for each free cell t next to s, but
 /// the one that the inverse of op leads to [see OperatorDomain]

 template< class F >
 void for_each_successor( state_type s , cost_type , std::uint8_t op ,
                          F && f ) const {
  const int x = int( s % f_w ) , y = int( s / f_w );
  if( ! f_map.cost( x , y ) )
   return;
  f_adj->for_each( s , x , y , [ & ]( state_type t , int nx , int ny ,
                                      cost_type c ) {
   const auto o = std::uint8_t( 3 * ( ny - y + 1 ) + nx - x + 1 );
   if( o + op == 8 )
    return;
   if( f_life )
    c = ( nx != x ) && ( ny != y ) ? cost_type( y * c_sqrt2 )
                                   : cost_type( y );
   f( t , c , distance( nx , ny , f_gx , f_gy ) , o );
   } );
  }

/*------------------------- FOR PARTIAL EXPANSION --------------------------*/
 /// the successors of s, but the one that op leads back to, whose
 /// \f$ \Delta f \f$ is delta up to cost_slack( delta ), returning the
 /// least larger one of the others [see PartialExpansionDomain]

 template< class F >
 cost_type for_each_successor_delta( state_type s , cost_type h ,
                                     std::uint8_t op , cost_type delta ,
                                     F && f ) const {
  const double slack = cost_slack( delta );
  cost_type next = infinite_cost< cost_type >();
  for_each_successor( s , h , op , [ & ]( state_type t , cost_type c ,
                                          cost_type ht , std::uint8_t o ) {
   const cost_type df = c + ht - h;
   if( ( double( df ) <= double( delta ) + slack ) &&
       ( double( df ) >= double( delta ) - slack ) )
    f( t , c , ht , o );
   else
    if( ( double( df ) > double( delta ) + slack ) && ( df < next ) )
     next = df;
   } );
  return( next );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the moves to the goal with no obstacle: the Manhattan distance with 4
 /// moves, the larger of the two differences with 8 [see DistanceDomain]

 [[nodiscard]] cost_type distance( state_type s ) const noexcept {
  const int ax = std::abs( int( s % f_w ) - f_gx ) ,
            ay = std::abs( int( s / f_w ) - f_gy );
  return( cost_type( Moves == 4 ? ax + ay : std::max( ax , ay ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c ) for each free cell t next to s, c being the cost of
 /// the move from t to s, i.e., the cost of the cell t (its row with the
 /// life costs), times \f$ \sqrt{2} \f$ for a diagonal [see
 /// PredecessorDomain]; the moves between two cells are the same both
 /// ways
 template< class F >
 void for_each_predecessor( state_type s , F && f ) const {
  const int x = int( s % f_w ) , y = int( s / f_w );
  if( ! f_map.cost( x , y ) )
   return;                                      // a blocked cell
  f_adj->for_each( s , x , y , [ & ]( state_type t , int nx , int ny ,
                                      cost_type ) {
   const cost_type c = cost_type( f_life ? ny : f_map.cost( nx , ny ) );
   f( t , ( nx != x ) && ( ny != y ) ? cost_type( c * c_sqrt2 ) : c );
   } );
  }

/*--------------------- METHODS OF THE BACKWARD SEARCH ---------------------*/
 /// the goal cell

 [[nodiscard]] state_type goal_state( void ) const noexcept {
  return( f_goal );
  }

 /// the estimate of the cost from the start to s, as the heuristic does
 /// from s to the goal

 [[nodiscard]] cost_type backward_heuristic( state_type s ) const noexcept {
  return( distance( f_sx , f_sy , int( s % f_w ) , int( s / f_w ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , hb( t ) ) for each predecessor t of s, as above

 template< class F >
 void for_each_predecessor( state_type s , cost_type , F && f ) const {
  for_each_predecessor( s , [ & ]( state_type t , cost_type c ) {
   f( t , c , backward_heuristic( t ) );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads again the costs of the cells of the map, after some of them
 /// have been changed by GridMap::set_cost(); the heuristic stays the one
 /// of the smallest cost of the map when the Grid was built, which the new
 /// costs of the free cells must not go below for it to stay admissible.
 /// Throws std::invalid_argument if the start or the goal is blocked

 void refresh( void ) {
  check_ends( "Grid::refresh" );
  f_adj.emplace( f_map , f_cut );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the smallest cost of a move: that of the cheapest free cell, 0 with the
 /// life costs (the moves within the row 0)

 [[nodiscard]] cost_type min_edge_cost( void ) const noexcept {
  return( f_life ? cost_type( 0 ) : cost_type( f_min ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells that the cell ( x , y ) has just been changed by
 /// GridMap::set_cost() from the cost old (0 if it was blocked): patches
 /// the moves around it and records the change of each move between two
 /// of the 3 x 3 cells around it, the only ones that the cell can change
 /// (its own, those into it and the diagonals that cut it), for
 /// drain_events() [see EventDomain]; the heuristic must stay admissible as
 /// for refresh(). Throws std::invalid_argument if the start or the goal
 /// is blocked

 void cell_changed( int x , int y , int old ) {
  check_ends( "Grid::cell_changed" );
  const auto before = [ & ]( int a , int b ) {
   return( ( a == x ) && ( b == y ) ? old : f_map.cost( a , b ) );
   };
  const auto after = [ & ]( int a , int b ) { return( f_map.cost( a , b ) );
   };
  for( int ay = y - 1 ; ay <= y + 1 ; ++ay )
   for( int ax = x - 1 ; ax <= x + 1 ; ++ax ) {
    if( ( ax < 0 ) || ( ay < 0 ) || ( ax >= f_w ) || ( ay >= f_h ) )
     continue;
    for( int k = 0 ; k < Moves ; ++k ) {
     const int bx = ax + c_dx[ k ] , by = ay + c_dy[ k ];
     if( ( std::abs( bx - x ) > 1 ) || ( std::abs( by - y ) > 1 ) ||
         ( bx < 0 ) || ( by < 0 ) || ( bx >= f_w ) || ( by >= f_h ) )
      continue;
     const cost_type from = move_cost( before , ax , ay , k );
     const cost_type to = move_cost( after , ax , ay , k );
     if( from == to )
      continue;
     const std::uint64_t e = std::uint64_t( index( ax , ay ) ) * 8 +
                             std::uint64_t( k );
     const auto [ it , fresh ] = f_pending.emplace( e , v_events.size() );
     if( fresh )
      v_events.push_back( { index( ax , ay ) , index( bx , by ) , from ,
                            to } );
     else
      v_events[ it->second ].to = to;
     }
    }
  f_adj->cell_changed( x , y );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( EdgeCostChanged ) for each move whose cost differs from the
 /// one at the last call [see cell_changed()], \f$ \infty \f$ standing
 /// for a move not allowed, and forgets the changes

 template< class F >
 void drain_events( F && f ) {
  for( const auto & e : v_events )
   if( e.from != e.to )
    f( e );
  v_events.clear();
  f_pending.clear();
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the columns and the rows

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( std::size_t( f_w ) + f_h );
  }

 /// calls f( x ) and f( W + y ) for the cell \f$ ( x , y ) \f$ of the
 /// state s

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

 /// the nblock of PBNF: the rectangle of nblock_width() times
 /// nblock_height() cells

 [[nodiscard]] std::size_t nblock( state_type s ) const noexcept {
  return( std::size_t( s / f_w / f_nsh ) * f_nw + s % f_w / f_nsw );
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

 /// the width of an nblock in cells, i.e., its side when the nblocks are
 /// squares

 [[nodiscard]] int nblock_side( void ) const noexcept { return( f_nsw ); }

 /// the width of an nblock in cells

 [[nodiscard]] int nblock_width( void ) const noexcept { return( f_nsw ); }

 /// the height of an nblock in cells

 [[nodiscard]] int nblock_height( void ) const noexcept { return( f_nsh ); }

 /// the side of the blocks of AZH and AHDA*

 [[nodiscard]] int block( void ) const noexcept { return( f_block ); }

/*----------------- FOR THE GREEDY AND THE WIDTH SEARCHES ------------------*/
 /// the preferred successors of s: the moves that lower the heuristic
 /// [see PreferredDomain]

 [[nodiscard]] auto preferences( state_type s ) const noexcept {
  return( [ this , hs = heuristic( s ) ]( state_type t ) {
   return( heuristic( t ) < hs );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( i ) for the features of the goal, its column and its row
 /// [see GoalFeatureDomain]

 template< class F >
 void for_each_goal_feature( F && f ) const {
  f( std::size_t( f_gx ) );
  f( std::size_t( f_w ) + std::size_t( f_gy ) );
  }

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

 /// \f$ \sqrt{2} \f$, the cost of a diagonal move per unit of cost of a
 /// cell
 static constexpr double c_sqrt2 = grid_sqrt2;

 /// the moves up, left, right and down, then to the corners
 static constexpr int c_dx[ 8 ] = { 0 , -1 , 1 , 0 , -1 , 1 , -1 , 1 };
 static constexpr int c_dy[ 8 ] = { -1 , 0 , 0 , 1 , -1 , -1 , 1 , 1 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the cost of the move k from the cell ( x , y ) with the costs of the
 // cells given by cost( x , y ), as for_each_successor() gives it (the
 // corners cut if they may be, the row with the life costs), ∞ if the
 // move is not allowed
 template< class C >
 cost_type move_cost( const C & cost , int x , int y , int k ) const {
  const int nx = x + c_dx[ k ] , ny = y + c_dy[ k ];
  const int c = cost( x , y );
  if( ! c || ! cost( nx , ny ) ||
      ( ( k >= 4 ) && ! f_cut &&
        ( ! cost( nx , y ) || ! cost( x , ny ) ) ) )
   return( infinite_cost< cost_type >() );
  const int m = f_life ? y : c;
  return( k >= 4 ? cost_type( m * c_sqrt2 ) : cost_type( m ) );
  }

 // throws std::invalid_argument if the start or the goal is blocked
 void check_ends( const char * who ) const {
  if( ! f_map.cost( int( f_start % f_w ) , int( f_start / f_w ) ) ||
      ! f_map.cost( f_gx , f_gy ) )
   throw( std::invalid_argument( std::string( who ) +
                                 ": blocked start or goal" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the Manhattan or octile distance times the smallest cost, or the
 // distance with the life costs

 [[nodiscard]] cost_type distance( int x , int y , int u , int v )
  const noexcept {
  if( f_life )
   return( life_distance( x , y , u , v ) );
  const int ax = std::abs( x - u ) , ay = std::abs( y - v );
  if constexpr( Moves == 4 )
   return( ( ax + ay ) * f_min );
  else
   return( ( std::max( ax , ay ) + ( c_sqrt2 - 1 ) * std::min( ax , ay ) )
           * f_min );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 // the cost V( a , b ) of going straight from the row a to the row b with
 // the life costs

 [[nodiscard]] static std::int64_t vertical( std::int64_t a ,
                                             std::int64_t b ) noexcept {
  if( a > b )
   return( ( a * ( a + 1 ) - b * ( b + 1 ) ) / 2 );
  return( ( b * ( b - 1 ) - a * ( a - 1 ) ) / 2 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the cost of a cheapest path from ( x , y ) to ( u , v ) with the life
 // costs, on the board without blocked cells and, with 8 moves, the moves
 // to a corner costing as the straight ones

 [[nodiscard]] cost_type life_distance( int x , int y , int u , int v )
  const noexcept {
  const std::int64_t dx = std::abs( x - u ) , dy = std::abs( y - v );
  if constexpr( Moves == 4 ) {
   const std::int64_t r = std::min( y , v );
   return( cost_type( std::min( vertical( y , 0 ) + vertical( 0 , v ) ,
                                vertical( y , r ) + dx * r +
                                vertical( r , v ) ) ) );
   }
  else {
   // the sum of max{ y - k , v - n + k , 0 } for k = 0 , ... , n - 1: the
   // first term is the largest up to k = c, the second one after
   const std::int64_t n = std::max( dx , dy ) , c = ( y - v + n ) / 2;
   const std::int64_t a = std::min( { c , n - 1 , std::int64_t( y ) } ) + 1;
   const std::int64_t lo = std::max( c + 1 , n - v );
   std::int64_t h = 0;
   if( a > 0 )
    h += a * y - a * ( a - 1 ) / 2;
   if( lo < n )
    h += ( ( v - n + lo ) + ( v - 1 ) ) * ( n - lo ) / 2;
   return( cost_type( h ) );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const GridMap & f_map;                           ///< the map
 std::optional< GridAdjacency< Moves > > f_adj;   ///< its adjacency
 int f_w;                                         ///< its width
 int f_h;                                         ///< its height
 state_type f_start;                              ///< the start cell
 state_type f_goal;                               ///< the goal cell
 int f_sx;                                        ///< the start column
 int f_sy;                                        ///< the start row
 int f_gx;                                        ///< the goal column
 int f_gy;                                        ///< the goal row
 int f_min;                                  ///< the smallest cell cost
 int f_block;                          ///< the side of the AZH blocks
 int f_nsw;                            ///< the columns of an nblock
 int f_nsh;                            ///< the rows of an nblock
 int f_bw;                             ///< the number of block columns
 int f_bh;                             ///< the number of block rows
 int f_nw;                             ///< the number of nblock columns
 int f_nh;                             ///< the number of nblock rows
 bool f_life;                          ///< the life costs
 bool f_cut;                           ///< corners may be cut
 std::vector< EdgeCostChanged< state_type , cost_type > > v_events;
                                       ///< the changes since the last drain
 std::unordered_map< std::uint64_t , std::size_t > f_pending;
                                       ///< the place of a move in v_events

 };  // end( class( Grid ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Grid */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Grid.h -----------------------------*/
/*--------------------------------------------------------------------------*/
