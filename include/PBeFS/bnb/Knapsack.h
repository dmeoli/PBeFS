/*--------------------------------------------------------------------------*/
/*---------------------------- File Knapsack.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the struct template BasicKnapsack, a 0-1 knapsack
 * problem that implements the BranchAndBound concept [see
 * BranchAndBound.h] with values of a given type, and Knapsack, the one with
 * values of type double, with the
 * random instances of the classes of Pisinger and the dynamic program that
 * gives their optimal value; the tests and the benchmark of the
 * Branch-and-Bound share it.
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

#ifndef __PBeFS_Knapsack
 #define __PBeFS_Knapsack
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <cmath>

#include <cstring>

#include <limits>

#include <random>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/bnb/BranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BasicKnapsack ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a 0-1 knapsack, its items sorted by decreasing ratio of profit to weight
/** A node fixes the first \f$ k \f$ items, its bound is the one of the
 * fractional relaxation of Dantzig, the greedy completion of the node is
 * its feasible value, and it is branched on the item \f$ k \f$; a node
 * over the capacity is infeasible. The children get the bound of their
 * parent, and the item is taken only if it fits; with explicit_children,
 * they get their own bound, and the item is taken anyway, giving a child
 * that is found infeasible when it is evaluated if it does not fit. With
 * an integer type V the bound is rounded down, which is still a bound since
 * the profits are integer, so that a node whose bound rounds down to the
 * incumbent is pruned. The nodes travel between the processes of MPI by
 * serialize() and deserialize(). */

template< class V = double >
struct BasicKnapsack {

 using value_type = V;            ///< the type of the values

 struct node_type {
  std::size_t k = 0;              ///< the items fixed
  long weight = 0;                ///< the weight of those taken
  long profit = 0;                ///< their profit
  std::vector< char > x;          ///< the items taken
  };

 static constexpr Sense sense = Sense::maximize;

 std::vector< long > p;           ///< the profits
 std::vector< long > w;           ///< the weights
 long cap = 0;                    ///< the capacity
 bool explicit_children = false;  ///< own bounds, infeasible children

 node_type root( void ) const {
  node_type n;
  n.x.assign( p.size() , 0 );
  return( n );
  }

 // the value of the greedy completion of n
 long greedy( const node_type & n ) const {
  long r = cap - n.weight , v = n.profit;
  for( std::size_t i = n.k ; i < p.size() ; ++i )
   if( w[ i ] <= r ) {
    r -= w[ i ];
    v += p[ i ];
    }
  return( v );
  }

 BasicEvaluation< V > evaluate( node_type & n ) const {
  BasicEvaluation< V > e;
  if( n.weight > cap ) {
   e.infeasible = true;
   return( e );
   }
  long r = cap - n.weight;
  double b = double( n.profit );
  std::size_t j = n.k;
  for( ; ( j < p.size() ) && ( w[ j ] <= r ) ; ++j ) {
   r -= w[ j ];
   b += double( p[ j ] );
   }
  if( j < p.size() )
   b += double( p[ j ] ) * double( r ) / double( w[ j ] );
  if constexpr( std::is_integral_v< V > )
   e.bound = V( std::floor( b ) );
  else
   e.bound = V( b );
  e.feasible = V( greedy( n ) );
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  if( n.k == p.size() )
   return;
  // the item k left out first and taken last, so that with ties a depth
  // first search goes on with the item taken, as the greedy order does
  const V none = BBValueTraits< V >::infinity();
  const auto bound = [ & ]( node_type & x ) {
   if( ! explicit_children )
    return( none );
   const BasicEvaluation< V > e = evaluate( x );
   return( e.infeasible ? none : e.bound );
   };
  node_type c = n;
  ++c.k;
  const bool fits = explicit_children || ( n.weight + w[ n.k ] <= cap );
  node_type t = c;
  const V bc = bound( c );
  f( std::move( c ) , bc );
  if( fits ) {
   t.weight += w[ n.k ];
   t.profit += p[ n.k ];
   t.x[ n.k ] = 1;
   const V bt = bound( t );
   f( std::move( t ) , bt );
   }
  }

 std::uint64_t hash( const node_type & n ) const {
  return( ADS::splitmix64( ( std::uint64_t( n.k ) << 40 ) ^
                 ( std::uint64_t( n.weight ) << 20 ) ^
                 std::uint64_t( n.profit ) ) );
  }

 V estimate( const node_type & n ) const {
  return( V( greedy( n ) ) );
  }

 // the bytes of a node, for MPI
 void serialize( const node_type & n , std::vector< char > & out ) const {
  const auto put = [ & ]( const auto & v ) {
   const auto * b = reinterpret_cast< const char * >( & v );
   out.insert( out.end() , b , b + sizeof( v ) );
   };
  put( n.k );
  put( n.weight );
  put( n.profit );
  out.insert( out.end() , n.x.begin() , n.x.end() );
  }

 node_type deserialize( const char * & in ) const {
  node_type n;
  std::memcpy( & n.k , in , sizeof( n.k ) );
  in += sizeof( n.k );
  std::memcpy( & n.weight , in , sizeof( n.weight ) );
  in += sizeof( n.weight );
  std::memcpy( & n.profit , in , sizeof( n.profit ) );
  in += sizeof( n.profit );
  n.x.assign( in , in + p.size() );
  in += p.size();
  return( n );
  }

 // the optimal value, by the dynamic program on the capacity
 long optimum( void ) const {
  std::vector< long > best( std::size_t( cap ) + 1 , 0 );
  for( std::size_t i = 0 ; i < p.size() ; ++i )
   for( long c = cap ; c >= w[ i ] ; --c )
    best[ c ] = std::max( best[ c ] , best[ c - w[ i ] ] + p[ i ] );
  return( best[ cap ] );
  }

 // a random instance of the class of Pisinger cls (0 uncorrelated, 1
 // weakly correlated, 2 strongly correlated), of n items of weight at most
 // r, the capacity half the total weight
 static BasicKnapsack random( int cls , std::size_t n , long r ,
                         std::uint64_t seed ) {
  std::mt19937_64 gen( seed );
  std::vector< std::pair< long , long > > it( n );
  long tw = 0;
  for( auto & [ pi , wi ] : it ) {
   wi = 1 + long( random_below( gen , std::uint64_t( r ) ) );
   if( cls == 0 )
    pi = 1 + long( random_below( gen , std::uint64_t( r ) ) );
   else
    if( cls == 1 )
     pi = std::max( 1L , wi - r / 10 + long( random_below(
                                   gen , std::uint64_t( r / 5 + 1 ) ) ) );
    else
     pi = wi + r / 10;
   tw += wi;
   }
  std::sort( it.begin() , it.end() , []( const auto & a , const auto & b ) {
   return( a.first * b.second > b.first * a.second );
   } );
  BasicKnapsack kp;
  for( const auto & [ pi , wi ] : it ) {
   kp.p.push_back( pi );
   kp.w.push_back( wi );
   }
  kp.cap = tw / 2;
  return( kp );
  }

 // true if n is a node of value v: its greedy completion fits and gives v
 bool check( const node_type & n , V v ) const {
  long wt = 0 , pt = 0;
  for( std::size_t i = 0 ; i < n.k ; ++i )
   if( n.x[ i ] ) {
    wt += w[ i ];
    pt += p[ i ];
    }
  return( ( wt == n.weight ) && ( pt == n.profit ) && ( wt <= cap ) &&
          ( V( greedy( n ) ) == v ) );
  }

 };  // end( struct( BasicKnapsack ) )

/*--------------------------------------------------------------------------*/
/// the 0-1 knapsack with values of type double

using Knapsack = BasicKnapsack<>;

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Knapsack */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Knapsack.h ---------------------------*/
/*--------------------------------------------------------------------------*/
