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
/** The variable i is the position i of the stack (0 the top), its value
 * the pancake in it, \f$ \lceil \log_2 N \rceil \f$ bits each; the flip of
 * the top k pancakes is the operator with the k effects \f$ x'_i = x_{k - 1
 * - i} \f$, \f$ i < k \f$, and no precondition. The heuristic is the gap
 * heuristic of Helmert (2010), a term per pair of adjacent positions and
 * one for the bottom pancake against the plate, and the backward one the
 * gaps of the stack renamed by the positions of the pancakes in the start,
 * which is the gap heuristic of the problem whose goal is the start. Since
 * a flip is a reversal, its relation keeps \f$ \lfloor k / 2 \rfloor \f$
 * values at once, and its diagram grows as \f$ N^{k / 2} \f$: the encoding
 * suits small stacks [see docs/theory/symbolic-search.md]. decode() gives
 * the stack of PancakePuzzle [see Pancake.h]. */

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
  return( std::vector< std::uint32_t >( f_start.begin() , f_start.end() ) );
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
   SymbolicOperator< cost_type > o;
   for( int i = 0 ; i < k ; ++i )
    o.effects.push_back( SymbolicEffect::copy_of(
     std::uint32_t( i ) , std::uint32_t( k - 1 - i ) ) );
   f( o );
   }
  }

 template< class F >
 static void for_each_heuristic_term( F && f ) {
  std::array< int , N > rank;
  for( int i = 0 ; i < N ; ++i )
   rank[ i ] = i;
  gaps( rank , f );
  }

 template< class F >
 void for_each_backward_heuristic_term( F && f ) const {
  std::array< int , N > rank;
  for( int i = 0 ; i < N ; ++i )
   rank[ f_start[ i ] ] = i;
  gaps( rank , f );
  }

 /// the stack of the values x

 [[nodiscard]] static state_type decode( const std::vector< std::uint32_t > &
                                         x ) {
  state_type s{};
  for( int i = 0 ; i < N ; ++i )
   s[ i ] = static_cast< std::uint8_t >( x[ i ] );
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the gaps of the stack renamed by rank, whose goal is the identity
 template< class F >
 static void gaps( const std::array< int , N > & rank , F && f ) {
  for( int i = 0 ; i + 1 < N ; ++i ) {
   SymbolicTerm< cost_type > t;
   t.vars = { std::uint32_t( i ) , std::uint32_t( i + 1 ) };
   t.table.assign( N * N , 0 );
   for( int x = 0 ; x < N ; ++x )
    for( int y = 0 ; y < N ; ++y )
     t.table[ x + N * y ] = std::abs( rank[ x ] - rank[ y ] ) > 1;
   f( t );
   }
  SymbolicTerm< cost_type > t;
  t.vars = { std::uint32_t( N - 1 ) };
  t.table.assign( N , 0 );
  for( int x = 0 ; x < N ; ++x )
   t.table[ x ] = rank[ x ] != N - 1;
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
