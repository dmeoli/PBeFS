/*--------------------------------------------------------------------------*/
/*--------------------------- File Permutation.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the rankings of the k-permutations of n items, i.e., the
 * bijections between the sequences \f$ a = ( a_0 , \ldots , a_{k-1} ) \f$
 * of k distinct values of \f$ [ 0 , n ) \f$ and the integers of
 * \f$ [ 0 , N ) \f$, \f$ N = n! / ( n - k )! \f$, which index the tables
 * of the pattern databases of the permutation domains [see
 * PermutationPDB.h]: the k items of a pattern are tracked, and
 * \f$ a_i \f$ is the location of the i-th of them.
 *
 * The three rankings of the file differ in the order and in the cost:
 *
 * - LexRanking gives the lexicographic rank, whose digits are
 *   \f$ d_i = a_i - | \{ j < i : a_j < a_i \} | \in [ 0 , n - i ) \f$ and
 *   \f$ \mathrm{rank}( a ) = \sum_i d_i \prod_{j > i} ( n - j ) \f$; with
 *   a word of n bits for the values already seen, each digit is a
 *   population count and each inverse digit the selection of a bit, so that
 *   both directions take \f$ O( k ) \f$ operations for \f$ n \leq 64 \f$.
 *   Two consecutive ranks differ only in the last item, so that close
 *   abstract states fall in close entries.
 *
 * - TreeRanking gives the same lexicographic rank for any \f$ n \leq 256
 *   \f$ by a complete binary tree of counters over the values, after Bonet
 *   (2008): the number of the values below \f$ a_i \f$ already seen is a
 *   walk from a leaf to the root, and the \f$ d \f$-th value not yet seen
 *   a walk from the root to a leaf, \f$ O( k \log n ) \f$ in all.
 *
 * - MyrvoldRuskeyRanking is the linear-time ranking of Myrvold and Ruskey
 *   (2001) restricted to k items: the unranking runs the first k steps of
 *   their unrank1 on the identity, \f$ \pi[ i - 1 ] \leftrightarrow
 *   \pi[ r \bmod i ] \f$ for \f$ i = n , n - 1 , \ldots , n - k + 1 \f$
 *   with \f$ r \gets \lfloor r / i \rfloor \f$ after each one, and reads
 *   \f$ a_j = \pi[ n - 1 - j ] \f$; the ranking replays the same swaps with
 *   the inverse of \f$ \pi \f$ at hand, the digit of the step \f$ j \f$
 *   being \f$ \pi^{-1}[ a_j ] \f$. Both take \f$ O( n ) \f$ operations,
 *   the identity included, and the order of the ranks has no locality.
 *
 * The ranks are of 64 bits, and the constructors throw if \f$ N \f$ does
 * not fit (for the full permutations, \f$ n \leq 20 \f$).
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

#ifndef __PBeFS_Permutation
 #define __PBeFS_Permutation
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <bit>

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <cstring>

#include <limits>

#include <stdexcept>

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the number \f$ n! / ( n - k )! \f$ of the k-permutations of n items, or
/// 0 if it does not fit in 64 bits; 0 also if k > n

[[nodiscard]] constexpr std::uint64_t k_permutations( unsigned n ,
                                                      unsigned k ) noexcept
{
 if( k > n )
  return( 0 );
 std::uint64_t r = 1;
 for( unsigned i = 0 ; i < k ; ++i ) {
  const std::uint64_t f = n - i;
  if( r > std::numeric_limits< std::uint64_t >::max() / f )
   return( 0 );
  r *= f;
  }
 return( r );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a hash of the n bytes at b, by words of 8 bytes mixed as in ADS (the
/// states of the permutation puzzles)

[[nodiscard]] inline std::uint64_t bytes_hash( const std::uint8_t * b ,
                                               std::size_t n ) noexcept
{
 std::uint64_t h = 0x9E3779B97F4A7C15ULL ^ n;
 std::size_t i = 0;
 for( ; i + 8 <= n ; i += 8 ) {
  std::uint64_t w;
  std::memcpy( & w , b + i , 8 );
  h = ADS::splitmix64( h ^ w );
  }
 if( i < n ) {
  std::uint64_t w = 0;
  std::memcpy( & w , b + i , n - i );
  h = ADS::splitmix64( h ^ w );
  }
 return( h );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a ranking of the k-permutations of n items
/** A type R is a Ranking if R( n , k ) builds the ranking of the
 * k-permutations of n items, size() is their number N, rank( a ) the rank
 * in [ 0 , N ) of the k distinct values a[ 0 ] , ... , a[ k - 1 ] of
 * [ 0 , n ), and unrank( r , a ) its inverse, which writes them in a; both
 * are const and use no state but the stack, so that many threads may call
 * them at once. */

//! [heuristics-sig-ranking]
template< class R >
concept Ranking = std::constructible_from< R , unsigned , unsigned > &&
 requires( const R & r , const std::uint8_t * a , std::uint8_t * b ,
           std::uint64_t x ) {
  { r.size() } -> std::convertible_to< std::uint64_t >;
  { r.rank( a ) } -> std::convertible_to< std::uint64_t >;
  r.unrank( x , b );
  };
//! [heuristics-sig-ranking]

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS LexRanking ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the lexicographic ranking of the k-permutations of \f$ n \leq 64 \f$
/// items, by population counts and selections of bits

class LexRanking {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 static constexpr const char * name = "lex";       ///< the name

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the ranking of the k-permutations of n items; throws
 /// std::invalid_argument if n > 64, k > n or N needs more than 64 bits

 LexRanking( unsigned n , unsigned k ) : f_n( n ) , f_k( k ) ,
                                         f_size( k_permutations( n , k ) ) {
  if( ( n > 64 ) || ( k > n ) || ( f_size == 0 ) )
   throw( std::invalid_argument( "LexRanking::LexRanking: n > 64, k > n "
                                 "or more than 2^64 ranks" ) );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the items

 [[nodiscard]] unsigned n( void ) const noexcept { return( f_n ); }

 /// the length of the sequences

 [[nodiscard]] unsigned k( void ) const noexcept { return( f_k ); }

 /// the number of the ranks, \f$ n! / ( n - k )! \f$

 [[nodiscard]] std::uint64_t size( void ) const noexcept { return( f_size ); }

/*-------------------------------- RANKING ---------------------------------*/
 /// the rank of the k distinct values of a

 [[nodiscard]] std::uint64_t rank( const std::uint8_t * a ) const noexcept {
  std::uint64_t r = 0 , seen = 0;
  for( unsigned i = 0 ; i < f_k ; ++i ) {
   const std::uint64_t below = seen & ( ( std::uint64_t( 1 ) << a[ i ] ) -
                                        1u );
   r = r * std::uint64_t( f_n - i ) +
       std::uint64_t( a[ i ] - unsigned( std::popcount( below ) ) );
   seen |= std::uint64_t( 1 ) << a[ i ];
   }
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in a the k values of rank r < size()

 void unrank( std::uint64_t r , std::uint8_t * a ) const noexcept {
  std::uint8_t d[ 64 ];
  for( unsigned i = f_k ; i-- > 0 ; ) {
   d[ i ] = std::uint8_t( r % std::uint64_t( f_n - i ) );
   r /= std::uint64_t( f_n - i );
   }
  std::uint64_t free = f_n == 64 ? ~std::uint64_t( 0 )
                                 : ( std::uint64_t( 1 ) << f_n ) - 1u;
  for( unsigned i = 0 ; i < f_k ; ++i ) {
   std::uint64_t x = free;
   for( unsigned j = d[ i ] ; j > 0 ; --j )
    x &= x - 1u;
   a[ i ] = std::uint8_t( std::countr_zero( x ) );
   free &= ~( std::uint64_t( 1 ) << a[ i ] );
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 unsigned f_n;                                  ///< the items
 unsigned f_k;                                  ///< the length
 std::uint64_t f_size;                          ///< the ranks

 };  // end( class( LexRanking ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS TreeRanking ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the lexicographic ranking of the k-permutations of \f$ n \leq 256 \f$
/// items, by a binary tree of counters
/** The tree has \f$ 2^L \geq n \f$ leaves, one per value, and each node
 * counts the values of its leaves already seen (by rank()) or not yet
 * taken (by unrank()); the left sibling of each node on the path of a leaf
 * gives the values below it. The counters are a local array of 511 bytes
 * at most. */

class TreeRanking {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 static constexpr const char * name = "tree";      ///< the name

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the ranking of the k-permutations of n items; throws
 /// std::invalid_argument if n > 256, k > n or N needs more than 64 bits

 TreeRanking( unsigned n , unsigned k ) : f_n( n ) , f_k( k ) ,
                                          f_size( k_permutations( n , k ) ) {
  if( ( n > 256 ) || ( k > n ) || ( f_size == 0 ) )
   throw( std::invalid_argument( "TreeRanking::TreeRanking: n > 256, k > "
                                 "n or more than 2^64 ranks" ) );
  f_levels = 0;
  while( ( 1u << f_levels ) < n )
   ++f_levels;
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the items

 [[nodiscard]] unsigned n( void ) const noexcept { return( f_n ); }

 /// the length of the sequences

 [[nodiscard]] unsigned k( void ) const noexcept { return( f_k ); }

 /// the number of the ranks, \f$ n! / ( n - k )! \f$

 [[nodiscard]] std::uint64_t size( void ) const noexcept { return( f_size ); }

/*-------------------------------- RANKING ---------------------------------*/
 /// the rank of the k distinct values of a
 /** The node \f$ v \geq 1 \f$ of the tree has the children \f$ 2 v \f$ and
  * \f$ 2 v + 1 \f$, and the leaf of the value x is \f$ 2^L + x \f$; going
  * up from a right child adds the count of its left sibling. */

 [[nodiscard]] std::uint64_t rank( const std::uint8_t * a ) const noexcept {
  std::uint8_t c[ 512 ] = {};
  const unsigned leaves = 1u << f_levels;
  std::uint64_t r = 0;
  for( unsigned i = 0 ; i < f_k ; ++i ) {
   unsigned v = leaves + a[ i ] , below = 0;
   while( v > 1 ) {
    if( v & 1u )
     below += c[ v - 1 ];
    ++c[ v ];
    v >>= 1;
    }
   r = r * std::uint64_t( f_n - i ) + std::uint64_t( a[ i ] - below );
   }
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in a the k values of rank r < size()
 /** Each node counts the values of its leaves not yet taken; the digit
  * d is found by going down to the left child if it has more than d free
  * values, and else to the right one with d decreased by them. */

 void unrank( std::uint64_t r , std::uint8_t * a ) const noexcept {
  std::uint8_t d[ 256 ];
  for( unsigned i = f_k ; i-- > 0 ; ) {
   d[ i ] = std::uint8_t( r % std::uint64_t( f_n - i ) );
   r /= std::uint64_t( f_n - i );
   }
  std::uint16_t c[ 512 ];
  const unsigned leaves = 1u << f_levels;
  for( unsigned v = 0 ; v < leaves ; ++v )
   c[ leaves + v ] = v < f_n ? 1 : 0;
  for( unsigned v = leaves ; v-- > 1 ; )
   c[ v ] = std::uint16_t( c[ 2 * v ] + c[ 2 * v + 1 ] );
  for( unsigned i = 0 ; i < f_k ; ++i ) {
   unsigned v = 1 , x = d[ i ];
   while( v < leaves ) {
    --c[ v ];
    v *= 2;
    if( x >= c[ v ] ) {
     x -= c[ v ];
     ++v;
     }
    }
   c[ v ] = 0;
   a[ i ] = std::uint8_t( v - leaves );
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 unsigned f_n;                                  ///< the items
 unsigned f_k;                                  ///< the length
 std::uint64_t f_size;                          ///< the ranks
 unsigned f_levels;                             ///< L, the depth

 };  // end( class( TreeRanking ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MyrvoldRuskeyRanking -----------------------*/
/*--------------------------------------------------------------------------*/
/// the ranking of Myrvold and Ruskey (2001) of the k-permutations of
/// \f$ n \leq 256 \f$ items, in linear time
/** The rank is \f$ \sum_j e_j \prod_{l < j} ( n - l ) \f$, the digit
 * \f$ e_j \in [ 0 , n - j ) \f$ being the position swapped with
 * \f$ n - 1 - j \f$ at the step \f$ j \f$ of the unranking [see the file
 * header]; for k = n it is the rank1 of their paper. */

class MyrvoldRuskeyRanking {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 static constexpr const char * name = "mr";        ///< the name

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the ranking of the k-permutations of n items; throws
 /// std::invalid_argument if n > 256, k > n or N needs more than 64 bits

 MyrvoldRuskeyRanking( unsigned n , unsigned k )
  : f_n( n ) , f_k( k ) , f_size( k_permutations( n , k ) ) {
  if( ( n > 256 ) || ( k > n ) || ( f_size == 0 ) )
   throw( std::invalid_argument( "MyrvoldRuskeyRanking::MyrvoldRuskeyRanking"
                                 ": n > 256, k > n or more than 2^64 ranks" )
          );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the items

 [[nodiscard]] unsigned n( void ) const noexcept { return( f_n ); }

 /// the length of the sequences

 [[nodiscard]] unsigned k( void ) const noexcept { return( f_k ); }

 /// the number of the ranks, \f$ n! / ( n - k )! \f$

 [[nodiscard]] std::uint64_t size( void ) const noexcept { return( f_size ); }

/*-------------------------------- RANKING ---------------------------------*/
 /// the rank of the k distinct values of a

 [[nodiscard]] std::uint64_t rank( const std::uint8_t * a ) const noexcept {
  std::uint8_t p[ 256 ] , q[ 256 ];              // pi and its inverse
  for( unsigned i = 0 ; i < f_n ; ++i )
   p[ i ] = q[ i ] = std::uint8_t( i );
  std::uint64_t r = 0 , m = 1;
  for( unsigned j = 0 ; j < f_k ; ++j ) {
   const unsigned top = f_n - 1 - j , e = q[ a[ j ] ];
   r += m * e;
   m *= std::uint64_t( f_n - j );
   const std::uint8_t x = p[ top ];               // swap p[ top ] , p[ e ]
   p[ top ] = p[ e ];
   p[ e ] = x;
   q[ p[ top ] ] = std::uint8_t( top );
   q[ x ] = std::uint8_t( e );
   }
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in a the k values of rank r < size()

 void unrank( std::uint64_t r , std::uint8_t * a ) const noexcept {
  std::uint8_t p[ 256 ];
  for( unsigned i = 0 ; i < f_n ; ++i )
   p[ i ] = std::uint8_t( i );
  for( unsigned j = 0 ; j < f_k ; ++j ) {
   const unsigned top = f_n - 1 - j;
   const unsigned e = unsigned( r % std::uint64_t( f_n - j ) );
   r /= std::uint64_t( f_n - j );
   const std::uint8_t x = p[ top ];
   p[ top ] = p[ e ];
   p[ e ] = x;
   a[ j ] = p[ top ];
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 unsigned f_n;                                  ///< the items
 unsigned f_k;                                  ///< the length
 std::uint64_t f_size;                          ///< the ranks

 };  // end( class( MyrvoldRuskeyRanking ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Permutation */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Permutation.h -------------------------*/
/*--------------------------------------------------------------------------*/
