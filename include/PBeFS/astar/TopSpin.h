/*--------------------------------------------------------------------------*/
/*----------------------------- File TopSpin.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template TopSpin, which implements the Domain
 * concept [see Domain.h] for the TopSpin puzzle of N tokens on a ring with
 * a turnstile of K, the permutation puzzle of Felner, Zahavi, Schaeffer and
 * Holte (IJCAI 2005) and of the pattern databases of Holte et al. (2006).
 *
 * A state is the ring, the token at each of the N positions; the move i
 * reverses the K tokens at the positions i , ... , i + K - 1 (modulo N), at
 * cost 1, and the goals are the N rotations of the identity, i.e., the
 * tokens in increasing order around the ring wherever the token 0 is. The
 * heuristic counts the breakpoints, the pairs of adjacent positions whose
 * tokens are not consecutive in either order: a move changes only the two
 * pairs at the ends of the turnstile, since the reversal keeps the pairs
 * inside it, so that half the breakpoints, rounded up, is consistent.
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

#ifndef __PBeFS_TopSpin
 #define __PBeFS_TopSpin
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
/*----------------------------- CLASS TopSpin ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the TopSpin puzzle of N tokens and a turnstile of K
/** The operator of the move i is i, its own inverse, which the successors
 * of an OperatorDomain skip; those of a tree search also skip the move j
 * after i when the two turnstiles are disjoint, and hence the moves
 * commute, and j < i. The features of Zobrist hashing are the pairs
 * ( position , token ), and the perfect hash the lexicographic rank of the
 * ring, for N up to 20. */

template< int N , int K = 4 >
class TopSpin {

 static_assert( ( K >= 2 ) && ( K <= N ) && ( N < 255 ) ,
                "TopSpin: 2 <= K <= N < 255" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::array< std::uint8_t , N >;   ///< a ring
 using cost_type = int;                               ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the start ring; throws std::invalid_argument if it
 /// is not a permutation

 explicit TopSpin( const state_type & start ) : f_start( start ) {
  std::array< bool , N > seen{};
  for( std::uint8_t x : start ) {
   if( ( x >= N ) || seen[ x ] )
    throw( std::invalid_argument( "TopSpin: the start is not a "
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
   if( s[ ( i + 1 ) % N ] != ( s[ i ] + 1 ) % N )
    return( false );
  return( true );
  }

 [[nodiscard]] static std::uint64_t hash( const state_type & s ) noexcept {
  return( bytes_hash( s.data() , N ) );
  }

 [[nodiscard]] static cost_type heuristic( const state_type & s ) noexcept {
  int b = 0;
  for( int i = 0 ; i < N ; ++i ) {
   const int x = s[ i ] , y = s[ ( i + 1 ) % N ];
   b += ( y != ( x + 1 ) % N ) && ( x != ( y + 1 ) % N );
   }
  return( ( b + 1 ) / 2 );
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

 /// the successors but the move op, which undoes the move that made s
 template< class F >
 void for_each_successor( const state_type & s , cost_type ,
                          std::uint8_t op , F && f ) const {
  for( int i = 0 ; i < N ; ++i )
   if( i != op ) {
    const state_type t = apply( s , i );
    f( t , 1 , heuristic( t ) , std::uint8_t( i ) );
    }
  }

 /// the successors but the move op and the moves j < op that commute with
 /// it
 template< class F >
 void for_each_tree_successor( const state_type & s , cost_type ,
                               std::uint8_t op , F && f ) const {
  for( int i = 0 ; i < N ; ++i ) {
   if( ( i == op ) || ( ( op != no_op ) && ( i < op ) &&
                        commute( i , op ) ) )
    continue;
   const state_type t = apply( s , i );
   f( t , 1 , heuristic( t ) , std::uint8_t( i ) );
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
 /// the perfect hash: the lexicographic rank of the ring

 [[nodiscard]] std::uint64_t perfect_hash( const state_type & s ) const
  noexcept requires( N <= 20 ) {
  static const LexRanking r( N , N );
  return( r.rank( s.data() ) );
  }

/*------------------------- THE PERMUTATION PUZZLE -------------------------*/
 /// one group of N positions [see PermutationDomain]

 [[nodiscard]] static std::vector< PermutationGroup > permutation_groups(
  void ) {
  return( std::vector< PermutationGroup >{
   PermutationGroup{ unsigned( N ) , 1u , false } } );
  }

 /// the N moves, in the order of their operators

 [[nodiscard]] static std::vector< PermutationMove > permutation_moves(
  void ) {
  std::vector< PermutationMove > v;
  for( int i = 0 ; i < N ; ++i ) {
   PermutationMove m;
   m.from.assign( 1 , std::vector< std::uint8_t >( N ) );
   for( int j = 0 ; j < N ; ++j )
    m.from[ 0 ][ j ] = std::uint8_t( j );
   for( int x = 0 ; x < K ; ++x )
    m.from[ 0 ][ ( i + x ) % N ] = std::uint8_t( ( i + K - 1 - x ) % N );
   v.push_back( std::move( m ) );
   }
  return( v );
  }

 /// the token at each position

 static void permutation_of( const state_type & s , unsigned ,
                             std::uint8_t * obj , std::uint8_t * ori )
  noexcept {
  for( int i = 0 ; i < N ; ++i ) {
   obj[ i ] = s[ i ];
   ori[ i ] = 0;
   }
  }

 /// the goals, the N rotations of the identity

 template< class F >
 static void for_each_goal_permutation( unsigned , F && f ) {
  std::uint8_t obj[ N ] , ori[ N ];
  for( int r = 0 ; r < N ; ++r ) {
   for( int i = 0 ; i < N ; ++i ) {
    obj[ i ] = std::uint8_t( ( i + r ) % N );
    ori[ i ] = 0;
    }
   f( obj , ori );
   }
  }

/*-------------------------------- HELPERS ---------------------------------*/
 /// the ring s after the move i

 [[nodiscard]] static state_type apply( const state_type & s , int i )
  noexcept {
  state_type t = s;
  for( int x = 0 ; x < K / 2 ; ++x )
   std::swap( t[ ( i + x ) % N ] , t[ ( i + K - 1 - x ) % N ] );
  return( t );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the turnstiles of the moves i and j are disjoint

 [[nodiscard]] static constexpr bool commute( int i , int j ) noexcept {
  const int d = ( ( j - i ) % N + N ) % N;
  return( ( d >= K ) && ( N - d >= K ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the ring in text, the tokens from the position 0, separated by blanks

 [[nodiscard]] static state_type parse( const std::string & text ) {
  std::istringstream in( text );
  state_type s{};
  for( auto & x : s ) {
   int v;
   if( ! ( in >> v ) || ( v < 0 ) || ( v >= N ) )
    throw( std::invalid_argument( "TopSpin::parse: bad ring" ) );
   x = std::uint8_t( v );
   }
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the identity after n random moves, the same for a seed on every
 /// platform (a uniform ring if n is 0)

 [[nodiscard]] static state_type random( std::uint64_t seed , int n = 0 ) {
  std::mt19937_64 gen( seed );
  state_type s{};
  for( int i = 0 ; i < N ; ++i )
   s[ i ] = std::uint8_t( i );
  if( n == 0 )
   for( int i = N - 1 ; i > 0 ; --i )
    std::swap( s[ i ] , s[ random_below( gen , std::uint64_t( i + 1 ) ) ] );
  else
   for( int k = 0 ; k < n ; ++k )
    s = apply( s , int( random_below( gen , N ) ) );
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 state_type f_start;                            ///< the start

 };  // end( class( TopSpin ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_TopSpin */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File TopSpin.h ---------------------------*/
/*--------------------------------------------------------------------------*/
