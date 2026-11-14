/*--------------------------------------------------------------------------*/
/*---------------------------- File Pancake.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class templates PancakePuzzle and
 * BurnedPancakePuzzle, which implement the Domain concept [see Domain.h]
 * for the pancake problem of Dweighter (1975) and for its burned variant,
 * with the gap heuristic of Helmert (2010).
 *
 * A stack of N pancakes of distinct sizes \f$ 0 , \ldots , N - 1 \f$ is a
 * permutation p, p[ 0 ] being the top; the move k, \f$ 2 \leq k \leq N
 * \f$, flips the top k pancakes, reversing p[ 0 .. k - 1 ], at cost 1, and
 * the goal is the identity. A gap is a pair of adjacent pancakes, or the
 * bottom one and the plate (size N), whose sizes differ by more than one:
 * a flip changes only the pair at its boundary, between the positions
 * k - 1 and k, and the goal has no gap, so that the number of gaps is a
 * consistent heuristic, computed for a successor in \f$ O( 1 ) \f$ from
 * that of its parent. The variant GAP-x of Holte et al. (AAAI 2016) does
 * not count the gaps between two pancakes one of which is among the x
 * smallest (the plate counts always), which weakens the heuristic in a
 * controlled way.
 *
 * In the burned variant a pancake also has a side burned, and a flip turns
 * every pancake that it moves; the goal is the identity with every burned
 * side down. A pair of adjacent pancakes x above y is good if y = x + 1 and
 * neither is turned, or x = y + 1 and both are turned (the pair is then
 * right after a flip of a block that holds it), and the bottom pancake is
 * good with the plate if it is N - 1, not turned; a flip of k turns the
 * pairs of its block into pairs of the same kind, so that only the pair at
 * its boundary changes, and the number of the pairs that are not good is
 * consistent and zero exactly at the goal.
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

#ifndef __PBeFS_Pancake
 #define __PBeFS_Pancake
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cstddef>

#include <cstdint>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Permutation.h"

#include "PBeFS/astar/PermutationPDB.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS PancakePuzzle ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the pancake problem with N pancakes and the gap heuristic
/** The operator of the flip of k is k, and the flip of k is its own
 * inverse, which the successors of an OperatorDomain skip; the moves are
 * generated from the flip of N down to the flip of 2. The features of
 * Zobrist hashing are the pairs ( position , size ), and the perfect hash
 * the lexicographic rank of the permutation, for N up to 20. */

template< int N >
class PancakePuzzle {

 static_assert( ( N >= 2 ) && ( N < 255 ) , "PancakePuzzle: 2 to 254" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::array< std::uint8_t , N >;   ///< a stack
 using cost_type = int;                               ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the start stack and the x of GAP-x (0: all the
 /// gaps); throws std::invalid_argument if start is not a permutation

 explicit PancakePuzzle( const state_type & start , unsigned ignore = 0 )
  : f_start( start ) , f_ignore( ignore ) {
  std::array< bool , N > seen{};
  for( std::uint8_t x : start ) {
   if( ( x >= N ) || seen[ x ] )
    throw( std::invalid_argument( "PancakePuzzle: the start is not a "
                                  "permutation" ) );
   seen[ x ] = true;
   }
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 [[nodiscard]] static bool is_goal( const state_type & s ) noexcept {
  for( int i = 0 ; i < N ; ++i )
   if( s[ i ] != i )
    return( false );
  return( true );
  }

 [[nodiscard]] static std::uint64_t hash( const state_type & s ) noexcept {
  return( bytes_hash( s.data() , N ) );
  }

 [[nodiscard]] cost_type heuristic( const state_type & s ) const noexcept {
  int h = 0;
  for( int i = 0 ; i + 1 < N ; ++i )
   h += gap( s[ i ] , s[ i + 1 ] );
  return( h + gap( s[ N - 1 ] , N ) );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f )
  const {
  for_each_successor( s , h , no_op , [ & ]( const state_type & t ,
                                             cost_type c , cost_type ht ,
                                             std::uint8_t ) {
   f( t , c , ht );
   } );
  }

 /// the successors but the flip op, which undoes the move that made s
 template< class F >
 void for_each_successor( const state_type & s , cost_type h ,
                          std::uint8_t op , F && f ) const {
  for( int k = N ; k >= 2 ; --k ) {
   if( k == op )
    continue;
   const int below = k < N ? s[ k ] : N;
   const int ht = h - gap( s[ k - 1 ] , below ) + gap( s[ 0 ] , below );
   state_type t = s;
   std::reverse( t.begin() , t.begin() + k );
   f( t , 1 , ht , std::uint8_t( k ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] static constexpr std::size_t num_features( void ) noexcept {
  return( std::size_t( N ) * N );
  }

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const {
  for( int i = 0 ; i < N ; ++i )
   f( std::size_t( i ) * N + s[ i ] );
  }

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const {
  for( int i = 0 ; i < N ; ++i )
   if( s[ i ] != t[ i ] ) {
    f( std::size_t( i ) * N + s[ i ] );
    f( std::size_t( i ) * N + t[ i ] );
    }
  }

/*------------------------------ OTHER HASHES ------------------------------*/
 /// the perfect hash: the lexicographic rank of the stack

 [[nodiscard]] std::uint64_t perfect_hash( const state_type & s ) const
  noexcept requires( N <= 20 ) {
  static const LexRanking r( N , N );
  return( r.rank( s.data() ) );
  }

/*------------------------- THE PERMUTATION PUZZLE -------------------------*/
 /// one group of N positions [see
 /// PermutationDomain]

 [[nodiscard]] static std::vector< PermutationGroup > permutation_groups(
  void ) {
  return( std::vector< PermutationGroup >{
   PermutationGroup{ unsigned( N ) , 1u , false } } );
  }

 /// the flips of 2 to N

 [[nodiscard]] static std::vector< PermutationMove > permutation_moves( void )
 {
  std::vector< PermutationMove > v;
  for( int k = 2 ; k <= N ; ++k ) {
   PermutationMove m;
   m.from.assign( 1 , std::vector< std::uint8_t >( N ) );
   m.turn.assign( 1 , std::vector< std::uint8_t >( N , 0 ) );
   for( int i = 0 ; i < N ; ++i )
    m.from[ 0 ][ i ] = std::uint8_t( i < k ? k - 1 - i : i );
   v.push_back( std::move( m ) );
   }
  return( v );
  }

 /// the pancake at each position

 static void permutation_of( const state_type & s , unsigned ,
                             std::uint8_t * obj , std::uint8_t * ori )
  noexcept {
  for( int i = 0 ; i < N ; ++i ) {
   obj[ i ] = std::uint8_t( s[ i ] );
   ori[ i ] = 0;
   }
  }

 /// the goal, the identity

 template< class F >
 static void for_each_goal_permutation( unsigned , F && f ) {
  std::uint8_t obj[ N ] , ori[ N ];
  for( int i = 0 ; i < N ; ++i ) {
   obj[ i ] = std::uint8_t( i );
   ori[ i ] = 0;
   }
  f( obj , ori );
  }

/*-------------------------------- HELPERS ---------------------------------*/
 /// the x of GAP-x

 [[nodiscard]] unsigned ignored( void ) const noexcept { return( f_ignore ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the stack of the sizes in text, from the top, separated by blanks;
 /// throws std::invalid_argument if they are not N

 [[nodiscard]] static state_type parse( const std::string & text ) {
  std::istringstream in( text );
  state_type s{};
  for( auto & x : s ) {
   int v;
   if( ! ( in >> v ) || ( v < 0 ) || ( v >= N ) )
    throw( std::invalid_argument( "PancakePuzzle::parse: bad stack" ) );
   x = std::uint8_t( v );
   }
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a uniformly random stack, the same for a seed on every platform

 [[nodiscard]] static state_type random( std::uint64_t seed ) {
  std::mt19937_64 gen( seed );
  state_type s{};
  for( int i = 0 ; i < N ; ++i )
   s[ i ] = std::uint8_t( i );
  for( int i = N - 1 ; i > 0 ; --i )
   std::swap( s[ i ] , s[ random_below( gen , std::uint64_t( i + 1 ) ) ] );
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// 1 if x above y is a gap that counts
 [[nodiscard]] int gap( int x , int y ) const noexcept {
  if( y == N )
   return( x != N - 1 );
  if( ( unsigned( std::min( x , y ) ) < f_ignore ) )
   return( 0 );
  return( ( x - y > 1 ) || ( y - x > 1 ) );
  }

 state_type f_start;                            ///< the start
 unsigned f_ignore;                             ///< the x of GAP-x

 };  // end( class( PancakePuzzle ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BurnedPancakePuzzle ------------------------*/
/*--------------------------------------------------------------------------*/
/// the burned pancake problem with N pancakes and its gap heuristic
/** A pancake is a byte, its size in the low 7 bits and the bit 7 set if
 * its burned side is up; the operators are as in PancakePuzzle, from the
 * flip of N down to the flip of 1, which only turns the top pancake. In
 * text a pancake is its size plus one, negative if turned. */

template< int N >
class BurnedPancakePuzzle {

 static_assert( ( N >= 1 ) && ( N <= 127 ) , "BurnedPancakePuzzle: 1 to "
                "127" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::array< std::uint8_t , N >;   ///< a stack
 using cost_type = int;                               ///< the cost type

 static constexpr std::uint8_t turned = 0x80;  ///< the bit of a turned one

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the start stack; throws std::invalid_argument if
 /// its sizes are not a permutation

 explicit BurnedPancakePuzzle( const state_type & start )
  : f_start( start ) {
  std::array< bool , N > seen{};
  for( std::uint8_t x : start ) {
   const int v = x & 0x7F;
   if( ( v >= N ) || seen[ v ] )
    throw( std::invalid_argument( "BurnedPancakePuzzle: the start is not "
                                  "a permutation" ) );
   seen[ v ] = true;
   }
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 [[nodiscard]] static bool is_goal( const state_type & s ) noexcept {
  for( int i = 0 ; i < N ; ++i )
   if( s[ i ] != i )
    return( false );
  return( true );
  }

 [[nodiscard]] static std::uint64_t hash( const state_type & s ) noexcept {
  return( bytes_hash( s.data() , N ) );
  }

 [[nodiscard]] static cost_type heuristic( const state_type & s ) noexcept {
  int h = 0;
  for( int i = 0 ; i + 1 < N ; ++i )
   h += bad( s[ i ] , s[ i + 1 ] );
  return( h + bad( s[ N - 1 ] , N ) );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f )
  const {
  for_each_successor( s , h , no_op , [ & ]( const state_type & t ,
                                             cost_type c , cost_type ht ,
                                             std::uint8_t ) {
   f( t , c , ht );
   } );
  }

 /// the successors but the flip op, which undoes the move that made s
 template< class F >
 void for_each_successor( const state_type & s , cost_type h ,
                          std::uint8_t op , F && f ) const {
  for( int k = N ; k >= 1 ; --k ) {
   if( k == op )
    continue;
   const int below = k < N ? s[ k ] : N;
   const int ht = h - bad( s[ k - 1 ] , below ) +
                  bad( std::uint8_t( s[ 0 ] ^ turned ) , below );
   state_type t = s;
   std::reverse( t.begin() , t.begin() + k );
   for( int i = 0 ; i < k ; ++i )
    t[ i ] ^= turned;
   f( t , 1 , ht , std::uint8_t( k ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] static constexpr std::size_t num_features( void ) noexcept {
  return( std::size_t( N ) * N * 2 );
  }

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const {
  for( int i = 0 ; i < N ; ++i )
   f( feature( i , s[ i ] ) );
  }

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const {
  for( int i = 0 ; i < N ; ++i )
   if( s[ i ] != t[ i ] ) {
    f( feature( i , s[ i ] ) );
    f( feature( i , t[ i ] ) );
    }
  }

/*------------------------- THE PERMUTATION PUZZLE -------------------------*/
 /// one group of N positions, the orientation being the burned side [see
 /// PermutationDomain]

 [[nodiscard]] static std::vector< PermutationGroup > permutation_groups(
  void ) {
  return( std::vector< PermutationGroup >{
   PermutationGroup{ unsigned( N ) , 2u , false } } );
  }

 /// the flips of 1 to N

 [[nodiscard]] static std::vector< PermutationMove > permutation_moves( void )
 {
  std::vector< PermutationMove > v;
  for( int k = 1 ; k <= N ; ++k ) {
   PermutationMove m;
   m.from.assign( 1 , std::vector< std::uint8_t >( N ) );
   m.turn.assign( 1 , std::vector< std::uint8_t >( N , 0 ) );
   for( int i = 0 ; i < N ; ++i )
    m.from[ 0 ][ i ] = std::uint8_t( i < k ? k - 1 - i : i );
   for( int i = 0 ; i < k ; ++i )
    m.turn[ 0 ][ i ] = 1;
   v.push_back( std::move( m ) );
   }
  return( v );
  }

 /// the pancake at each position and whether it is turned

 static void permutation_of( const state_type & s , unsigned ,
                             std::uint8_t * obj , std::uint8_t * ori )
  noexcept {
  for( int i = 0 ; i < N ; ++i ) {
   obj[ i ] = std::uint8_t( s[ i ] & 0x7F );
   ori[ i ] = std::uint8_t( s[ i ] >> 7 );
   }
  }

 /// the goal, the identity

 template< class F >
 static void for_each_goal_permutation( unsigned , F && f ) {
  std::uint8_t obj[ N ] , ori[ N ];
  for( int i = 0 ; i < N ; ++i ) {
   obj[ i ] = std::uint8_t( i );
   ori[ i ] = 0;
   }
  f( obj , ori );
  }

/*-------------------------------- HELPERS ---------------------------------*/
 /// the stack in text, from the top: N sizes from 1 to N, negative if the
 /// burned side is up; throws std::invalid_argument if it is not that

 [[nodiscard]] static state_type parse( const std::string & text ) {
  std::istringstream in( text );
  state_type s{};
  for( auto & x : s ) {
   int v;
   if( ! ( in >> v ) || ( v == 0 ) || ( v > N ) || ( v < -N ) )
    throw( std::invalid_argument( "BurnedPancakePuzzle::parse: bad stack" ) );
   x = std::uint8_t( v > 0 ? v - 1 : ( -v - 1 ) | turned );
   }
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a uniformly random stack, each pancake turned with probability 1/2

 [[nodiscard]] static state_type random( std::uint64_t seed ) {
  std::mt19937_64 gen( seed );
  state_type s{};
  for( int i = 0 ; i < N ; ++i )
   s[ i ] = std::uint8_t( i );
  for( int i = N - 1 ; i > 0 ; --i )
   std::swap( s[ i ] , s[ random_below( gen , std::uint64_t( i + 1 ) ) ] );
  for( auto & x : s )
   if( random_below( gen , 2 ) )
    x |= turned;
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the feature of the pancake x at the position i
 [[nodiscard]] static std::size_t feature( int i , std::uint8_t x ) noexcept {
  return( ( std::size_t( i ) * N + ( x & 0x7F ) ) * 2 + ( x >> 7 ) );
  }

 /// 1 unless x above y (N for the plate) is a good pair
 [[nodiscard]] static int bad( std::uint8_t x , int y ) noexcept {
  const int a = x & 0x7F;
  const bool ta = x & turned;
  if( y == N )
   return( ( a != N - 1 ) || ta );
  const int b = y & 0x7F;
  const bool tb = y & turned;
  if( ta != tb )
   return( 1 );
  return( ta ? a != b + 1 : b != a + 1 );
  }

 state_type f_start;                            ///< the start

 };  // end( class( BurnedPancakePuzzle ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Pancake */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Pancake.h ----------------------------*/
/*--------------------------------------------------------------------------*/
