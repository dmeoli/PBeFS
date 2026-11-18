/*--------------------------------------------------------------------------*/
/*------------------------------- File QAP.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for QAP, the quadratic assignment problem of Koopmans and
 * Beckmann as a problem of the Branch-and-Bound (B&B), with the bound of
 * Gilmore (1962) and Lawler (1963), and a reader of the instances of
 * QAPLIB. It is an example of the B&B layer, outside the engine.
 *
 * Given the flows \f$ F = ( f_{ij} ) \f$ between \f$ n \f$ facilities and
 * the distances \f$ D = ( d_{kl} ) \f$ between \f$ n \f$ locations, the
 * problem is \f$ \min_\pi \sum_{i,j} f_{ij} d_{\pi(i)\pi(j)} \f$ over the
 * permutations \f$ \pi \f$. The facilities are assigned in the order of
 * their decreasing total flow, and a node assigns the first \f$ d \f$ of
 * them; its children assign the next facility to each free location, in
 * increasing order, with no estimate (the least long, so that in the
 * search a child gets the bound of its parent until it is evaluated), and
 * child( n , k ) gives the child k at once, without evaluating n, so that
 * the work stealing of intervals applies [see IntervalBranchAndBound].
 * For the set \f$ U \f$ of the facilities left and \f$ L \f$ of the
 * locations left, the bound is the cost of the pairs of facilities
 * assigned plus the value of the linear assignment of \f$ U \f$ to
 * \f$ L \f$ whose cost of \f$ i \to l \f$ is \f$ f_{ii} d_{ll} \f$ plus
 * the cost of the flows of \f$ i \f$ with the facilities assigned, plus
 * the least scalar product of the flows \f$ f_{ij} \f$,
 * \f$ j \in U \setminus \{ i \} \f$, with the distances \f$ d_{lk} \f$,
 * \f$ k \in L \setminus \{ l \} \f$ (the first sorted up, the second
 * down); it is solved by the Hungarian method in \f$ O( |U|^3 ) \f$. Each
 * ordered pair of facilities is counted once, by its first facility, so
 * that the bound holds for asymmetric matrices too. A complete assignment
 * is feasible, of its cost.
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

#ifndef __PBeFS_QAP
 #define __PBeFS_QAP
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstdint>

#include <fstream>

#include <functional>

#include <limits>

#include <numeric>

#include <optional>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/bnb/BranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the least cost of an assignment of the rows to the columns of the
/// square matrix c of size m (by rows), by the Hungarian method with the
/// potentials, in \f$ O( m^3 ) \f$

inline long linear_assignment( const std::vector< long > & c , std::size_t m )
{
 if( m == 0 )
  return( 0 );
 const long inf = std::numeric_limits< long >::max() / 4;
 std::vector< long > u( m + 1 , 0 ) , v( m + 1 , 0 ) , way( m + 1 , 0 );
 std::vector< std::size_t > row( m + 1 , 0 );   // the row of each column
 for( std::size_t i = 1 ; i <= m ; ++i ) {
  row[ 0 ] = i;
  std::size_t j0 = 0;
  std::vector< long > mins( m + 1 , inf );
  std::vector< char > used( m + 1 , 0 );
  do {
   used[ j0 ] = 1;
   const std::size_t i0 = row[ j0 ];
   long delta = inf;
   std::size_t j1 = 0;
   for( std::size_t j = 1 ; j <= m ; ++j )
    if( ! used[ j ] ) {
     const long cur = c[ ( i0 - 1 ) * m + j - 1 ] - u[ i0 ] - v[ j ];
     if( cur < mins[ j ] ) {
      mins[ j ] = cur;
      way[ j ] = long( j0 );
      }
     if( mins[ j ] < delta ) {
      delta = mins[ j ];
      j1 = j;
      }
     }
   for( std::size_t j = 0 ; j <= m ; ++j )
    if( used[ j ] ) {
     u[ row[ j ] ] += delta;
     v[ j ] -= delta;
     }
    else
     mins[ j ] -= delta;
   j0 = j1;
   } while( row[ j0 ] != 0 );
  do {
   const std::size_t j1 = std::size_t( way[ j0 ] );
   row[ j0 ] = row[ j1 ];
   j0 = j1;
   } while( j0 );
  }
 long cost = 0;
 for( std::size_t j = 1 ; j <= m ; ++j )
  cost += c[ ( row[ j ] - 1 ) * m + j - 1 ];
 return( cost );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS QAP --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the quadratic assignment problem, by the Gilmore-Lawler bound

struct QAP {

 using value_type = long;          ///< the type of the values

 /// a node: the locations of the first facilities of the order
 struct node_type {
  std::vector< std::int16_t > loc;  ///< the location of each facility, or -1
  std::uint16_t d = 0;              ///< the facilities assigned
  };

 static constexpr Sense sense = Sense::minimize;

 std::size_t n = 0;                ///< the facilities and the locations
 std::vector< long > f;            ///< the flows, by rows
 std::vector< long > dist;         ///< the distances, by rows
 std::vector< std::size_t > order; ///< the order of the facilities

 /// the flow from i to j
 [[nodiscard]] long flow( std::size_t i , std::size_t j ) const {
  return( f[ i * n + j ] );
  }

 /// the distance from k to l
 [[nodiscard]] long distance( std::size_t k , std::size_t l ) const {
  return( dist[ k * n + l ] );
  }

 /// the instance of a file of QAPLIB: n, then the matrices A and B, the
 /// cost being \f$ \sum a_{ij} b_{\pi(i)\pi(j)} \f$; throws
 /// std::runtime_error if it cannot be read
 static QAP qaplib( const std::string & file ) {
  std::ifstream in( file );
  QAP q;
  if( ! ( in >> q.n ) )
   throw( std::runtime_error( "QAP::qaplib: cannot read " + file ) );
  q.f.resize( q.n * q.n );
  q.dist.resize( q.n * q.n );
  for( auto & x : q.f )
   in >> x;
  for( auto & x : q.dist )
   in >> x;
  if( ! in )
   throw( std::runtime_error( "QAP::qaplib: bad matrices in " + file ) );
  q.prepare();
  return( q );
  }

 /// the order of the facilities, by decreasing total flow (the lower first
 /// on ties); to be called after f is set
 void prepare( void ) {
  order.resize( n );
  std::iota( order.begin() , order.end() , std::size_t( 0 ) );
  std::vector< long > t( n , 0 );
  for( std::size_t i = 0 ; i < n ; ++i )
   for( std::size_t j = 0 ; j < n ; ++j )
    t[ i ] += flow( i , j ) + flow( j , i );
  std::stable_sort( order.begin() , order.end() ,
                    [ & ]( std::size_t a , std::size_t b ) {
                     return( t[ a ] > t[ b ] );
                     } );
  }

 node_type root( void ) const {
  return( node_type{ std::vector< std::int16_t >( n , -1 ) , 0 } );
  }

 /// the cost of the pairs of the facilities assigned
 [[nodiscard]] long fixed_cost( const node_type & x ) const {
  long c = 0;
  for( std::size_t a = 0 ; a < x.d ; ++a )
   for( std::size_t b = 0 ; b < x.d ; ++b ) {
    const std::size_t i = order[ a ] , j = order[ b ];
    c += flow( i , j ) * distance( std::size_t( x.loc[ i ] ) ,
                                   std::size_t( x.loc[ j ] ) );
    }
  return( c );
  }

 /// the Gilmore-Lawler bound of x, its cost if it is complete
 [[nodiscard]] long glb( const node_type & x ) const {
  const long fixed = fixed_cost( x );
  if( x.d == n )
   return( fixed );
  std::vector< std::size_t > U( order.begin() + x.d , order.end() ) , L;
  std::vector< char > used( n , 0 );
  for( std::size_t a = 0 ; a < x.d ; ++a )
   used[ std::size_t( x.loc[ order[ a ] ] ) ] = 1;
  for( std::size_t l = 0 ; l < n ; ++l )
   if( ! used[ l ] )
    L.push_back( l );
  const std::size_t m = U.size();
  std::vector< long > c( m * m );
  std::vector< long > fu , dl;
  for( std::size_t a = 0 ; a < m ; ++a ) {
   const std::size_t i = U[ a ];
   fu.clear();
   for( const std::size_t j : U )
    if( j != i )
     fu.push_back( flow( i , j ) );
   std::sort( fu.begin() , fu.end() );
   for( std::size_t b = 0 ; b < m ; ++b ) {
    const std::size_t l = L[ b ];
    long lin = flow( i , i ) * distance( l , l );
    for( std::size_t k = 0 ; k < x.d ; ++k ) {
     const std::size_t j = order[ k ];
     const std::size_t lj = std::size_t( x.loc[ j ] );
     lin += flow( i , j ) * distance( l , lj ) +
            flow( j , i ) * distance( lj , l );
     }
    dl.clear();
    for( const std::size_t k : L )
     if( k != l )
      dl.push_back( distance( l , k ) );
    std::sort( dl.begin() , dl.end() , std::greater<>() );
    long q = 0;
    for( std::size_t t = 0 ; t < fu.size() ; ++t )
     q += fu[ t ] * dl[ t ];
    c[ a * m + b ] = lin + q;
    }
   }
  return( fixed + linear_assignment( c , m ) );
  }

 BasicEvaluation< long > evaluate( node_type & x ) const {
  BasicEvaluation< long > e;
  e.bound = glb( x );
  if( x.d == n )
   e.feasible = e.bound;
  return( e );
  }

 template< class F >
 void branch( const node_type & x , F && sink ) const {
  if( x.d == n )
   return;
  std::vector< char > used( n , 0 );
  for( std::size_t a = 0 ; a < x.d ; ++a )
   used[ std::size_t( x.loc[ order[ a ] ] ) ] = 1;
  const std::size_t i = order[ x.d ];
  for( std::size_t l = 0 ; l < n ; ++l )
   if( ! used[ l ] ) {
    node_type c = x;
    c.loc[ i ] = std::int16_t( l );
    ++c.d;
    sink( std::move( c ) , std::numeric_limits< long >::min() );
    }
  }

 /// the child k of x, the next facility at its k-th free location, without
 /// the evaluation of x
 std::optional< BBChild< node_type , long > > child( const node_type & x ,
                                                     std::size_t k ) const {
  if( x.d == n )
   return( std::nullopt );
  std::vector< char > used( n , 0 );
  for( std::size_t a = 0 ; a < x.d ; ++a )
   used[ std::size_t( x.loc[ order[ a ] ] ) ] = 1;
  for( std::size_t l = 0 ; l < n ; ++l )
   if( ! used[ l ] && ( k-- == 0 ) ) {
    node_type c = x;
    c.loc[ order[ x.d ] ] = std::int16_t( l );
    ++c.d;
    const long none = std::numeric_limits< long >::min();
    return( BBChild< node_type , long >{ std::move( c ) , none } );
    }
  return( std::nullopt );
  }

 /// the number of children of a node of depth d
 [[nodiscard]] std::size_t arity( unsigned d ) const {
  return( d < n ? n - d : 0 );
  }

 /// the cost of the permutation of a complete node
 [[nodiscard]] long cost( const node_type & x ) const {
  long c = 0;
  for( std::size_t i = 0 ; i < n ; ++i )
   for( std::size_t j = 0 ; j < n ; ++j )
    c += flow( i , j ) * distance( std::size_t( x.loc[ i ] ) ,
                                   std::size_t( x.loc[ j ] ) );
  return( c );
  }

 };  // end( struct( QAP ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_QAP */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File QAP.h -----------------------------*/
/*--------------------------------------------------------------------------*/
