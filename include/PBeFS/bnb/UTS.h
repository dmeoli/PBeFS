/*--------------------------------------------------------------------------*/
/*------------------------------- File UTS.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for UTS, the Unbalanced Tree Search benchmark of Olivier et
 * al. (2007) as a problem of the Branch-and-Bound (B&B) in the mode
 * enumerate [see BBMode], whose tree is implicit and whose count of nodes
 * is the answer. It is an example of the B&B layer, outside the engine.
 *
 * A node holds the 20 bytes of the state of a random generator and its
 * depth: the state of the child \f$ i \f$ is the SHA-1 digest of the state
 * of its parent followed by the four bytes of \f$ i \f$ (most significant
 * first), the state of the root the digest of 16 zero bytes followed by
 * the four bytes of the seed, and the last four bytes of a state, with the
 * sign bit cleared, are the random number \f$ r \f$ of the node, read as
 * the probability \f$ u = r / 2^{31} \f$. In a binomial tree the root has
 * \f$ b_0 \f$ children, and any other node \f$ m \f$ children if
 * \f$ u < q \f$ and none otherwise; in a geometric tree the number of
 * children of a node of depth \f$ d \f$ is
 * \f$ \lfloor \log( 1 - u ) / \log( 1 - p ) \rfloor \f$ with
 * \f$ p = 1 / ( 1 + b_d ) \f$, i.e., geometric of mean \f$ b_d \f$, where
 * \f$ b_d = b_0 \f$ at the root and otherwise
 * \f$ b_0 ( 1 - d / g ) \f$ (linear),
 * \f$ b_0 d^{- \log b_0 / \log g } \f$ (exponential decrease),
 * \f$ b_0^{\sin( 2 \pi d / g )} \f$ up to \f$ d = 5 g \f$ and 0 below
 * (cyclic), or \f$ b_0 \f$ up to \f$ d = g \f$ and 0 below (fixed), with
 * \f$ g \f$ the depth parameter; a non-root node has at most 100 children,
 * and the root of a binomial tree at most \f$ \lceil b_0 \rceil \f$. The
 * three trees of about four million nodes of the benchmark are T1
 * (geometric, fixed, \f$ b_0 = 4 \f$, \f$ g = 10 \f$, seed 19), T2
 * (geometric, cyclic, \f$ b_0 = 6 \f$, \f$ g = 16 \f$, seed 502) and T3
 * (binomial, \f$ b_0 = 2000 \f$, \f$ m = 8 \f$, \f$ q = 0.124875 \f$, seed
 * 42), given by UTS::t1(), t2() and t3().
 *
 * Every node is feasible, of value 0, so that BBStats::solutions counts the
 * nodes of the tree; the children come one at a time [see BBCursor], and
 * child( n , k ) gives the child k of a node at once, so that the searches
 * that send a node as the path that leads to it rebuild it without an
 * evaluation. The nodes are trivially copyable, and therefore travel over
 * MPI as they are.
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

#ifndef __PBeFS_UTS
 #define __PBeFS_UTS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cmath>

#include <cstdint>

#include <cstring>

#include <numbers>

#include <optional>

#include "PBeFS/bnb/BranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS SHA1 -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the SHA-1 digest of a short message, as in FIPS 180-4

struct SHA1 {

 /// the digest of the n bytes at m
 static std::array< std::uint8_t , 20 > digest( const std::uint8_t * m ,
                                                std::size_t n ) {
  std::uint32_t h[ 5 ] = { 0x67452301u , 0xEFCDAB89u , 0x98BADCFEu ,
                           0x10325476u , 0xC3D2E1F0u };
  std::uint8_t b[ 64 ];
  const std::uint64_t bits = std::uint64_t( n ) * 8;
  std::size_t at = 0;
  bool one = false , done = false;
  while( ! done ) {
   std::size_t k = 0;
   for( ; ( k < 64 ) && ( at < n ) ; ++k , ++at )
    b[ k ] = m[ at ];
   if( ( k < 64 ) && ! one ) {
    b[ k++ ] = 0x80;
    one = true;
    }
   if( k <= 56 && one ) {
    for( ; k < 56 ; ++k )
     b[ k ] = 0;
    for( int j = 7 ; j >= 0 ; --j )
     b[ k++ ] = std::uint8_t( bits >> ( 8 * j ) );
    done = true;
    }
   else
    for( ; k < 64 ; ++k )
     b[ k ] = 0;
   block( h , b );
   }
  std::array< std::uint8_t , 20 > d;
  for( int i = 0 ; i < 5 ; ++i )
   for( int j = 0 ; j < 4 ; ++j )
    d[ 4 * i + j ] = std::uint8_t( h[ i ] >> ( 24 - 8 * j ) );
  return( d );
  }

 private:

 static std::uint32_t rol( std::uint32_t x , int s ) {
  return( ( x << s ) | ( x >> ( 32 - s ) ) );
  }

 // the compression of a block of 64 bytes into h
 static void block( std::uint32_t * h , const std::uint8_t * b ) {
  std::uint32_t w[ 80 ];
  for( int t = 0 ; t < 16 ; ++t )
   w[ t ] = ( std::uint32_t( b[ 4 * t ] ) << 24 ) |
            ( std::uint32_t( b[ 4 * t + 1 ] ) << 16 ) |
            ( std::uint32_t( b[ 4 * t + 2 ] ) << 8 ) |
            std::uint32_t( b[ 4 * t + 3 ] );
  for( int t = 16 ; t < 80 ; ++t )
   w[ t ] = rol( w[ t - 3 ] ^ w[ t - 8 ] ^ w[ t - 14 ] ^ w[ t - 16 ] , 1 );
  std::uint32_t a = h[ 0 ] , bb = h[ 1 ] , c = h[ 2 ] , d = h[ 3 ] ,
                e = h[ 4 ];
  for( int t = 0 ; t < 80 ; ++t ) {
   std::uint32_t f , k;
   if( t < 20 ) {
    f = ( bb & c ) | ( ~bb & d );
    k = 0x5A827999u;
    }
   else
    if( t < 40 ) {
     f = bb ^ c ^ d;
     k = 0x6ED9EBA1u;
     }
    else
     if( t < 60 ) {
      f = ( bb & c ) | ( bb & d ) | ( c & d );
      k = 0x8F1BBCDCu;
      }
     else {
      f = bb ^ c ^ d;
      k = 0xCA62C1D6u;
      }
   const std::uint32_t x = rol( a , 5 ) + f + e + k + w[ t ];
   e = d;
   d = c;
   c = rol( bb , 30 );
   bb = a;
   a = x;
   }
  h[ 0 ] += a;
  h[ 1 ] += bb;
  h[ 2 ] += c;
  h[ 3 ] += d;
  h[ 4 ] += e;
  }

 };  // end( struct( SHA1 ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS UTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the unbalanced tree search, as a B&B that counts the nodes

struct UTS {

 using value_type = long;          ///< the type of the values

 /// a node: the state of its generator and its depth
 struct node_type {
  std::array< std::uint8_t , 20 > state;  ///< the state
  unsigned depth = 0;                     ///< the depth
  };

 static constexpr Sense sense = Sense::minimize;

 /// the shape of a geometric tree
 enum class Shape { linear = 0 , cyclic = 2 , fixed = 3 , expdec = 1 };

 bool binomial = false;            ///< binomial, or else geometric
 double b0 = 4;                    ///< the children of the root, or mean
 int gen_mx = 6;                   ///< the depth parameter g
 Shape shape = Shape::linear;      ///< the shape of a geometric tree
 double q = 0;                     ///< binomial: the probability of m
 int m = 0;                        ///< binomial: the children of a node
 int seed = 0;                     ///< the seed of the root

 /// the tree T1 of the benchmark (4130071 nodes)
 static UTS t1( void ) {
  UTS u;
  u.b0 = 4;
  u.gen_mx = 10;
  u.shape = Shape::fixed;
  u.seed = 19;
  return( u );
  }

 /// the tree T2 of the benchmark (4117769 nodes)
 static UTS t2( void ) {
  UTS u;
  u.b0 = 6;
  u.gen_mx = 16;
  u.shape = Shape::cyclic;
  u.seed = 502;
  return( u );
  }

 /// the tree T3 of the benchmark (4112897 nodes)
 static UTS t3( void ) {
  UTS u;
  u.binomial = true;
  u.b0 = 2000;
  u.m = 8;
  u.q = 0.124875;
  u.seed = 42;
  return( u );
  }

 /// a geometric tree of the given shape
 static UTS geometric( Shape s , double b , int g , int r ) {
  UTS u;
  u.shape = s;
  u.b0 = b;
  u.gen_mx = g;
  u.seed = r;
  return( u );
  }

 /// a binomial tree
 static UTS binomial_tree( double b , int mm , double qq , int r ) {
  UTS u;
  u.binomial = true;
  u.b0 = b;
  u.m = mm;
  u.q = qq;
  u.seed = r;
  return( u );
  }

 node_type root( void ) const {
  std::uint8_t b[ 20 ] = {};
  put( b + 16 , std::uint32_t( seed ) );
  return( node_type{ SHA1::digest( b , 20 ) , 0 } );
  }

 /// the number of children of n
 int arity( const node_type & n ) const {
  const std::uint32_t r = ( std::uint32_t( n.state[ 16 ] ) << 24 ) |
                          ( std::uint32_t( n.state[ 17 ] ) << 16 ) |
                          ( std::uint32_t( n.state[ 18 ] ) << 8 ) |
                          std::uint32_t( n.state[ 19 ] );
  const double u = double( r & 0x7fffffffu ) / 2147483648.0;
  if( binomial ) {
   if( n.depth == 0 )
    return( int( std::floor( b0 ) ) );
   return( u < q ? m : 0 );
   }
  double b = b0;
  const double d = double( n.depth );
  if( n.depth > 0 )
   switch( shape ) {
    case( Shape::expdec ):
     b = b0 * std::pow( d , - std::log( b0 ) / std::log( double( gen_mx ) ) );
     break;
    case( Shape::cyclic ):
     b = n.depth > unsigned( 5 * gen_mx ) ? 0 :
         std::pow( b0 , std::sin( 2 * std::numbers::pi * d /
                                  double( gen_mx ) ) );
     break;
    case( Shape::fixed ):
     b = n.depth < unsigned( gen_mx ) ? b0 : 0;
     break;
    default:
     b = b0 * ( 1 - d / double( gen_mx ) );
    }
  const double p = 1 / ( 1 + b );
  const int k = int( std::floor( std::log( 1 - u ) / std::log( 1 - p ) ) );
  return( std::min( k , 100 ) );
  }

 BasicEvaluation< long > evaluate( node_type & ) const {
  BasicEvaluation< long > e;
  e.feasible = 0;
  return( e );
  }

 /// the child k of n, if n has one
 std::optional< BBChild< node_type , long > > child( const node_type & n ,
                                                     std::size_t k ) const {
  if( long( k ) >= arity( n ) )
   return( std::nullopt );
  std::uint8_t b[ 24 ];
  std::memcpy( b , n.state.data() , 20 );
  put( b + 20 , std::uint32_t( k ) );
  return( BBChild< node_type , long >{
           node_type{ SHA1::digest( b , 24 ) , n.depth + 1 } , 0 } );
  }

 /// the children of n, one at a time
 struct Cursor {
  const UTS * u;
  node_type n;
  int count;
  int k = 0;
  std::optional< BBChild< node_type , long > > next( void ) {
   if( k >= count )
    return( std::nullopt );
   return( u->child( n , std::size_t( k++ ) ) );
   }
  };

 Cursor children( const node_type & n ) const {
  return( Cursor{ this , n , arity( n ) } );
  }

 /// the number of nodes of the subtree of n, by a plain recursion
 long count( const node_type & n ) const {
  long c = 1;
  const int a = arity( n );
  for( int k = 0 ; k < a ; ++k )
   c += count( child( n , std::size_t( k ) )->node );
  return( c );
  }

 private:

 static void put( std::uint8_t * b , std::uint32_t v ) {
  b[ 0 ] = std::uint8_t( v >> 24 );
  b[ 1 ] = std::uint8_t( v >> 16 );
  b[ 2 ] = std::uint8_t( v >> 8 );
  b[ 3 ] = std::uint8_t( v );
  }

 };  // end( struct( UTS ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_UTS */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File UTS.h -----------------------------*/
/*--------------------------------------------------------------------------*/
