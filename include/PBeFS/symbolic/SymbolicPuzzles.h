/*--------------------------------------------------------------------------*/
/*------------------------- File SymbolicPuzzles.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class templates SymbolicSlidingTile and
 * SymbolicPancake, the sliding-tile puzzle and the pancake problem as
 * SymbolicDomain [see SymbolicDomain.h], with the encodings of their states
 * into bit vectors, their moves as operators that copy variables, and their
 * heuristics as sums of terms: the Manhattan distance and the gap heuristic
 * of Helmert (2010), towards the goal and towards the start.
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

#ifndef __PBeFS_SymbolicPuzzles
 #define __PBeFS_SymbolicPuzzles
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <stdexcept>

#include <vector>

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/symbolic/SymbolicDomain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SymbolicSlidingTile ------------------------*/
/*--------------------------------------------------------------------------*/
/// the sliding-tile puzzle on \f$ W \times H \f$ cells as a SymbolicDomain
/** The variable c is the cell c, its value the tile in it (0 the blank),
 * so that a state takes \f$ W H \lceil \log_2 W H \rceil \f$ bits, 36 for
 * the 8-puzzle and 64 for the 15-puzzle. The move of the blank from the
 * cell b to the next cell c is the operator with the precondition \f$ x_b =
 * 0 \f$ and the effects \f$ x'_b = x_c \f$ and \f$ x'_c = 0 \f$, which
 * reads and writes two variables only; the heuristic is the Manhattan
 * distance, a term per cell, \f$ h_c( t ) \f$ being the distance of the
 * tile t in the cell c from its goal cell (0 for the blank), and the
 * backward heuristic the same distance from the cells of the start. The
 * boards are not checked for parity, so that an unsolvable instance can be
 * given; decode() gives the packed board of SlidingTile [see
 * SlidingTile.h]. */

template< int W , int H >
class SymbolicSlidingTile {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using explicit_type = SlidingTile< W , H >;     ///< the explicit puzzle
 using state_type = typename explicit_type::state_type;  ///< its state
 using cost_type = int;                          ///< the cost type
 using board_type = typename explicit_type::board_type;  ///< a board

 static constexpr int cells = W * H;             ///< the number of cells

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the puzzle from start to goal; throws std::invalid_argument if one of
 /// them is not a permutation of the tiles

 explicit SymbolicSlidingTile( const board_type & start ,
                               const board_type & goal =
                                explicit_type::korf_goal() )
  : f_start( check( start ) ) , f_goal( check( goal ) ) {}

/*-------------------------- METHODS OF THE DOMAIN -------------------------*/

 [[nodiscard]] static std::size_t num_variables( void ) noexcept {
  return( cells );
  }

 [[nodiscard]] static std::size_t range( std::uint32_t ) noexcept {
  return( cells );
  }

 [[nodiscard]] std::vector< std::uint32_t > initial_values( void ) const {
  return( std::vector< std::uint32_t >( f_start.begin() , f_start.end() ) );
  }

 template< class F >
 void for_each_goal( F && f ) const {
  std::vector< SymbolicFact > p;
  for( int c = 0 ; c < cells ; ++c )
   p.push_back( { std::uint32_t( c ) , std::uint32_t( f_goal[ c ] ) } );
  f( p );
  }

 template< class F >
 void for_each_operator( F && f ) const {
  for( int b = 0 ; b < cells ; ++b )
   for( int c : neighbours( b ) ) {
    SymbolicOperator< cost_type > o;
    o.pre = { { std::uint32_t( b ) , 0 } };
    o.effects = { SymbolicEffect::copy_of( std::uint32_t( b ) ,
                                           std::uint32_t( c ) ) ,
                  SymbolicEffect::assign( std::uint32_t( c ) , 0 ) };
    f( o );
    }
  }

 template< class F >
 void for_each_heuristic_term( F && f ) const { manhattan( f_goal , f ); }

 template< class F >
 void for_each_backward_heuristic_term( F && f ) const {
  manhattan( f_start , f );
  }

 /// the packed board of the values x

 [[nodiscard]] static state_type decode( const std::vector< std::uint32_t > &
                                         x ) {
  board_type b{};
  for( int c = 0 ; c < cells ; ++c )
   b[ c ] = static_cast< std::uint8_t >( x[ c ] );
  return( explicit_type::pack( b ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 static const board_type & check( const board_type & b ) {
  std::array< bool , cells > seen{};
  for( auto t : b ) {
   if( ( t >= cells ) || seen[ t ] )
    throw( std::invalid_argument( "SymbolicSlidingTile: not a permutation "
                                  "of the tiles" ) );
   seen[ t ] = true;
   }
  return( b );
  }

 static std::vector< int > neighbours( int b ) {
  std::vector< int > v;
  const int r = b / W , k = b % W;
  if( r > 0 )
   v.push_back( b - W );
  if( k > 0 )
   v.push_back( b - 1 );
  if( k < W - 1 )
   v.push_back( b + 1 );
  if( r < H - 1 )
   v.push_back( b + W );
  return( v );
  }

 template< class F >
 static void manhattan( const board_type & to , F && f ) {
  int cell[ cells ] = {};
  for( int c = 0 ; c < cells ; ++c )
   cell[ to[ c ] ] = c;
  for( int c = 0 ; c < cells ; ++c ) {
   SymbolicTerm< cost_type > t;
   t.vars = { std::uint32_t( c ) };
   t.table.assign( cells , 0 );
   for( int x = 1 ; x < cells ; ++x )
    t.table[ x ] = std::abs( c / W - cell[ x ] / W ) +
                   std::abs( c % W - cell[ x ] % W );
   f( t );
   }
  }

 board_type f_start;                            ///< the start
 board_type f_goal;                             ///< the goal

 };  // end( class( SymbolicSlidingTile ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SymbolicPancake --------------------------*/
/*--------------------------------------------------------------------------*/
/// the pancake problem with N pancakes as a SymbolicDomain
/** The variable p is the position of the pancake p in the stack (0 the
 * top), \f$ \lceil \log_2 N \rceil \f$ bits each, i.e., the inverse of the
 * stack: the flip of the top k pancakes is then the operator whose effect on
 * each pancake p is the map \f$ x'_p = k - 1 - x_p \f$ if \f$ x_p < k \f$,
 * \f$ x'_p = x_p \f$ otherwise, a relation on the bits of p alone, so that
 * the relation of a flip has \f$ O( N^2 \log N ) \f$ nodes. With the stack
 * itself as the state, a flip would copy the position i into the position k
 * - 1 - i, and its relation, which reads the positions in order, would keep
 * up to k values at once [see docs/theory/symbolic-search.md]. The heuristic
 * is the gap heuristic of Helmert (2010): there is a gap above the pancake
 * p + 1 unless the pancakes p and p + 1 are next to each other, one term per
 * pair of pancakes of consecutive sizes, and one for the largest pancake
 * against the plate; the backward one is the same with the sizes renamed by
 * the positions in the start, i.e., the gap heuristic of the problem whose
 * goal is the start. decode() gives the stack of PancakePuzzle [see
 * Pancake.h]. */

template< int N >
class SymbolicPancake {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using explicit_type = PancakePuzzle< N >;          ///< the explicit problem
 using state_type = typename explicit_type::state_type;  ///< a stack
 using cost_type = int;                             ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the problem of sorting start; throws std::invalid_argument if it is
 /// not a permutation

 explicit SymbolicPancake( const state_type & start ) : f_start( start ) {
  std::array< bool , N > seen{};
  for( auto x : start ) {
   if( ( x >= N ) || seen[ x ] )
    throw( std::invalid_argument( "SymbolicPancake: the start is not a "
                                  "permutation" ) );
   seen[ x ] = true;
   }
  }

/*-------------------------- METHODS OF THE DOMAIN -------------------------*/

 [[nodiscard]] static std::size_t num_variables( void ) noexcept {
  return( N );
  }

 [[nodiscard]] static std::size_t range( std::uint32_t ) noexcept {
  return( N );
  }

 [[nodiscard]] std::vector< std::uint32_t > initial_values( void ) const {
  std::vector< std::uint32_t > x( N );
  for( int i = 0 ; i < N ; ++i )
   x[ f_start[ i ] ] = std::uint32_t( i );
  return( x );
  }

 template< class F >
 static void for_each_goal( F && f ) {
  std::vector< SymbolicFact > p;
  for( int i = 0 ; i < N ; ++i )
   p.push_back( { std::uint32_t( i ) , std::uint32_t( i ) } );
  f( p );
  }

 template< class F >
 static void for_each_operator( F && f ) {
  for( int k = 2 ; k <= N ; ++k ) {
   std::vector< std::uint32_t > flip( N );
   for( int i = 0 ; i < N ; ++i )
    flip[ i ] = std::uint32_t( i < k ? k - 1 - i : i );
   SymbolicOperator< cost_type > o;
   for( int p = 0 ; p < N ; ++p )
    o.effects.push_back( SymbolicEffect::map_of( std::uint32_t( p ) ,
                                                 std::uint32_t( p ) ,
                                                 flip ) );
   f( o );
   }
  }

 template< class F >
 static void for_each_heuristic_term( F && f ) {
  std::array< int , N > by_rank;
  for( int i = 0 ; i < N ; ++i )
   by_rank[ i ] = i;
  gaps( by_rank , f );
  }

 template< class F >
 void for_each_backward_heuristic_term( F && f ) const {
  std::array< int , N > by_rank;
  for( int i = 0 ; i < N ; ++i )
   by_rank[ i ] = f_start[ i ];
  gaps( by_rank , f );
  }

 /// the stack of the values x, the positions of the pancakes

 [[nodiscard]] static state_type decode( const std::vector< std::uint32_t > &
                                         x ) {
  state_type s{};
  for( int p = 0 ; p < N ; ++p )
   s[ x[ p ] ] = static_cast< std::uint8_t >( p );
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the gaps of the problem whose goal is the stack by_rank[ 0 ] , ... ,
 // by_rank[ N - 1 ]: a gap above by_rank[ r + 1 ] unless it is next to
 // by_rank[ r ], and one below by_rank[ N - 1 ] unless it is at the bottom
 template< class F >
 static void gaps( const std::array< int , N > & by_rank , F && f ) {
  for( int r = 0 ; r + 1 < N ; ++r ) {
   SymbolicTerm< cost_type > t;
   t.vars = { std::uint32_t( by_rank[ r ] ) ,
              std::uint32_t( by_rank[ r + 1 ] ) };
   t.table.assign( N * N , 0 );
   for( int x = 0 ; x < N ; ++x )
    for( int y = 0 ; y < N ; ++y )
     t.table[ x + N * y ] = std::abs( x - y ) != 1;
   f( t );
   }
  SymbolicTerm< cost_type > t;
  t.vars = { std::uint32_t( by_rank[ N - 1 ] ) };
  t.table.assign( N , 0 );
  for( int x = 0 ; x < N ; ++x )
   t.table[ x ] = x != N - 1;
  f( t );
  }

 state_type f_start;                            ///< the start

 };  // end( class( SymbolicPancake ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SymbolicPuzzles */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SymbolicPuzzles.h ----------------------*/
/*--------------------------------------------------------------------------*/
