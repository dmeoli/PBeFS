/*--------------------------------------------------------------------------*/
/*------------------------ File TravelingSalesman.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class TravelingSalesman, which implements the Domain
 * concept [see Domain.h] for the traveling salesperson problem with the
 * minimum spanning tree heuristic of exercise 3.30 of Russell and Norvig,
 * "Artificial Intelligence: A Modern Approach" (3rd edition, 2010).
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

#ifndef __PBeFS_TravelingSalesman
 #define __PBeFS_TravelingSalesman
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <limits>

#include <random>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS TravelingSalesman -------------------------*/
/*--------------------------------------------------------------------------*/
/// the traveling salesperson problem on at most 32 cities
/** The tour starts and ends at the city 0. A state is the current city and
 * the set of the visited ones, packed as the visited cities in the bits
 * \f$ 0 , \ldots , 31 \f$, the current city in the bits
 * \f$ 32 , \ldots , 36 \f$ and, in the bit 37, whether the tour has gone
 * back to the city 0; the word is also its perfect hash. The heuristic is
 * the cost of a minimum spanning tree (by the algorithm of Prim) of the
 * unvisited cities together with the current and the first one, which is
 * admissible because the rest of the tour is a path through all of them.
 *
 * The features of Zobrist hashing are the visit of the city \f$ k \f$
 * (numbered \f$ k < n \f$), the current city \f$ c \f$ (numbered
 * \f$ n + c \f$) and the closing of the tour (numbered \f$ 2 n \f$). */

class TravelingSalesman {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint64_t;               ///< a packed state
 using cost_type = double;                       ///< the cost type
 using matrix_type = std::vector< std::vector< double > >;
                                                 ///< the distances

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the \f$ n \times n \f$ distances,
 /// \f$ 2 \leq n \leq 32 \f$

 explicit TravelingSalesman( matrix_type dist )
  : v_d( std::move( dist ) ) , f_n( int( v_d.size() ) ) {
  if( ( f_n < 2 ) || ( f_n > 32 ) )
   throw( std::invalid_argument(
                   "TravelingSalesman::TravelingSalesman: 2 to 32 cities" ) );
  f_all = ( std::uint64_t( 1 ) << f_n ) - 1;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Euclidean distances of n random points of the unit square, drawn
 /// from the seed

 [[nodiscard]] static matrix_type random_cities( int n ,
                                                 std::uint64_t seed ) {
  if( n < 0 )
   throw( std::invalid_argument(
                "TravelingSalesman::random_cities: negative n" ) );
  std::mt19937_64 gen( seed );
  std::vector< double > x( n ) , y( n );
  for( int c = 0 ; c < n ; ++c ) {
   x[ c ] = random_unit( gen );
   y[ c ] = random_unit( gen );
   }
  matrix_type d( n , std::vector< double >( n ) );
  for( int a = 0 ; a < n ; ++a )
   for( int b = 0 ; b < n ; ++b )
    d[ a ][ b ] = std::hypot( x[ a ] - x[ b ] , y[ a ] - y[ b ] );
  return( d );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the current city of the state s

 [[nodiscard]] static int city( std::uint64_t s ) noexcept {
  return( int( ( s >> 32 ) & 31 ) );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] std::uint64_t initial_state( void ) const noexcept {
  return( 1 );
  }

 [[nodiscard]] bool is_goal( std::uint64_t s ) const noexcept {
  return( ( s >> 37 ) & 1 );
  }

 [[nodiscard]] std::uint64_t hash( std::uint64_t s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost of a minimum spanning tree of the cities still to be touched

 [[nodiscard]] double heuristic( std::uint64_t s ) const {
  if( is_goal( s ) )
   return( 0 );
  int u[ 32 ];
  std::size_t m = 0;
  for( int c = 0 ; c < f_n ; ++c )
   if( ! ( s & ( std::uint64_t( 1 ) << c ) ) || ( c == 0 ) ||
       ( c == city( s ) ) )
    u[ m++ ] = c;
  double key[ 32 ];
  bool in[ 32 ];
  for( std::size_t k = 0 ; k < m ; ++k ) {
   key[ k ] = std::numeric_limits< double >::infinity();
   in[ k ] = false;
   }
  key[ 0 ] = 0;
  double tot = 0;
  for( std::size_t it = 0 ; it < m ; ++it ) {
   std::size_t b = 0;
   double best = std::numeric_limits< double >::infinity();
   for( std::size_t k = 0 ; k < m ; ++k )
    if( ! in[ k ] && ( key[ k ] < best ) )
     best = key[ b = k ];
   in[ b ] = true;
   tot += best;
   for( std::size_t k = 0 ; k < m ; ++k )
    if( ! in[ k ] )
     key[ k ] = std::min( key[ k ] , v_d[ u[ b ] ][ u[ k ] ] );
   }
  return( tot );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) ) for each unvisited city, or for the way
 /// back to the city 0 once all are visited

 template< class F >
 void for_each_successor( std::uint64_t s , double , F && f ) const {
  const int c = city( s );
  if( ( s & f_all ) == f_all ) {                  // back to the first city
   const std::uint64_t t = ( s & f_all ) | ( std::uint64_t( 1 ) << 37 );
   f( t , v_d[ c ][ 0 ] , 0.0 );
   return;
   }
  for( int k = 1 ; k < f_n ; ++k )
   if( ! ( s & ( std::uint64_t( 1 ) << k ) ) ) {
    const std::uint64_t t = ( ( s & f_all ) | ( std::uint64_t( 1 ) << k ) ) |
                            ( std::uint64_t( k ) << 32 );
    f( t , v_d[ c ][ k ] , heuristic( t ) );
    }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the visited cities, the current one, the end

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( 2 * std::size_t( f_n ) + 1 );
  }

 /// calls f( k ) for each visited city \f$ k \f$, f( n + c ) for the
 /// current city \f$ c \f$ and f( 2 n ) if the tour is closed

 template< class F >
 void for_each_feature( std::uint64_t s , F && f ) const {
  for( std::uint64_t v = s & f_all ; v ; v &= v - 1 )
   f( std::size_t( std::countr_zero( v ) ) );
  f( std::size_t( f_n + city( s ) ) );
  if( is_goal( s ) )
   f( 2 * std::size_t( f_n ) );
  }

/*------------------------------ OTHER HASHES ------------------------------*/
 /// the perfect hash: the packed state

 [[nodiscard]] std::uint64_t perfect_hash( std::uint64_t s ) const noexcept {
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost of an optimal tour by the dynamic program of Held and Karp
 /** With \f$ d_{jk} \f$ the distance from the city \f$ j \f$ to the city
  * \f$ k \f$ and \f$ C( S , j ) \f$ the least cost of a path that starts
  * at the city 0, visits exactly the cities of the set \f$ S \ni 0 \f$ and
  * ends at \f$ j \in S \f$, the dynamic program of Held and Karp (Journal
  * of the Society for Industrial and Applied Mathematics, 1962) rests on
  * the recurrence
  * \f[
  *  C( S \cup \{ k \} , k ) = \min_{j \in S} \{ C( S , j ) + d_{jk} \} ,
  *  \qquad k \notin S ,
  * \f]
  * and the cost of an optimal tour is
  * \f$ \min_{j \neq 0} \{ C( V , j ) + d_{j0} \} \f$, \f$ V \f$ being
  * the set of all the cities:
  *
  * \code{.unparsed}
  *  1  procedure HELDKARP(d, n)
  *  2     C({0}, 0) ← 0;  C(S, j) ← ∞ for each other pair
  *  3     for each set S ∋ 0, by increasing value of its bits, do
  *  4        for each j ∈ S with C(S, j) < ∞ do
  *  5           for each k ∉ S do
  *  6              C(S ∪ {k}, k) ← min(C(S ∪ {k}, k), C(S, j) + d(j, k))
  *  7     return min over j ≠ 0 of C(V, j) + d(j, 0)
  * \endcode
  *
  * Line 3 visits each set after all its subsets, since adding a city
  * raises the value of the bits; the table has \f$ 2^n n \f$ entries and
  * lines 3 to 6 take \f$ O( 2^n n^2 ) \f$ time, and therefore the method
  * is meant for checking the searches on small instances (the example
  * TravelingSalesman uses it up to 20 cities). */

 [[nodiscard]] double held_karp( void ) const {
  const int n = f_n;
  const double inf = std::numeric_limits< double >::infinity();
  std::vector< double > c( ( std::size_t( 1 ) << n ) * std::size_t( n ) ,
                           inf );
  const auto at = [ & ]( std::uint64_t s , int j ) -> double & {
   return( c[ std::size_t( s ) * std::size_t( n ) + std::size_t( j ) ] );
   };
  at( 1 , 0 ) = 0;
  for( std::uint64_t s = 1 ; s <= f_all ; s += 2 )
   for( int j = 0 ; j < n ; ++j )
    if( ( s & ( std::uint64_t( 1 ) << j ) ) && ( at( s , j ) < inf ) )
     for( int k = 1 ; k < n ; ++k )
      if( ! ( s & ( std::uint64_t( 1 ) << k ) ) ) {
       double & t = at( s | ( std::uint64_t( 1 ) << k ) , k );
       t = std::min( t , at( s , j ) + v_d[ j ][ k ] );
       }
  double best = inf;
  for( int j = 1 ; j < n ; ++j )
   best = std::min( best , at( f_all , j ) + v_d[ j ][ 0 ] );
  return( best );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of cities

 [[nodiscard]] int size( void ) const noexcept { return( f_n ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 matrix_type v_d;                               ///< the distances
 int f_n;                                       ///< the number of cities
 std::uint64_t f_all;                           ///< the set of all cities

 };  // end( class( TravelingSalesman ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_TravelingSalesman */

/*--------------------------------------------------------------------------*/
/*---------------------- End File TravelingSalesman.h ----------------------*/
/*--------------------------------------------------------------------------*/
