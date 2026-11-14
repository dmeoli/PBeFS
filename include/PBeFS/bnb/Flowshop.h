/*--------------------------------------------------------------------------*/
/*---------------------------- File Flowshop.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for Flowshop, the permutation flowshop of the least makespan
 * (\f$ F_m | prmu | C_{\max} \f$) as a problem of the Branch-and-Bound
 * (B&B), with the one-machine and the two-machine (Johnson) bounds of
 * Lageweg, Lenstra and Rinnooy Kan (1978), the branching that fixes the
 * next job of the head or of the tail of the sequence, and the instances of
 * Taillard (1993), generated from their seeds. It is an example of the B&B
 * layer, outside the engine.
 *
 * A node fixes a head \f$ \sigma_1 \f$ and a tail \f$ \sigma_2 \f$ of the
 * sequence, the unscheduled jobs \f$ U \f$ being in between. Let
 * \f$ f_k \f$ be the completion time of \f$ \sigma_1 \f$ on the machine
 * \f$ k \f$ (the least head \f$ \min_j \sum_{l < k} p_{lj} \f$ if
 * \f$ \sigma_1 \f$ is empty), \f$ b_k \f$ the time from the start of
 * \f$ \sigma_2 \f$ on \f$ k \f$ to the end (the least tail if it is
 * empty), and \f$ r_k = \sum_{j \in U} p_{kj} \f$. The one-machine bound is
 * \f$ \max_k ( t_k + b_k ) \f$ with \f$ t_0 = f_0 + r_0 \f$ and
 * \f$ t_k = \max \{ t_{k-1} , f_k + r_k \} \f$, and the two-machine bound
 * is the largest, over the pairs \f$ k < l \f$ of machines, of the
 * makespan of the two-machine flowshop with time lags
 * \f$ \sum_{k < i < l} p_{ij} \f$ on \f$ U \f$, started at \f$ f_k \f$ and
 * \f$ f_l \f$ and followed by \f$ b_k \f$ and \f$ b_l \f$, which the
 * order of Johnson for the pair (computed once on all the jobs) solves
 * exactly. A complete sequence is feasible, of its makespan.
 *
 * The children of a node append a job of \f$ U \f$ to \f$ \sigma_1 \f$
 * (forward), or prepend it to \f$ \sigma_2 \f$ (backward); the direction is
 * fixed, alternates with the depth, or (dynamic) is the one with the fewer
 * children whose bound is below the incumbent, the larger sum of their
 * bounds on ties, which needs the bounds of all the children of both
 * directions. The bound of each child is computed when the node is
 * branched, and is its estimate; its evaluation does not compute it again.
 * The children come by increasing job, or, with sort_children, by
 * increasing bound (the job on ties). make_worker( i ) gives an evaluator
 * that keeps the incumbent of cutoff( z ), which the dynamic branching
 * needs; child( n , k ) gives the child k of a node without evaluating it,
 * for the searches that rebuild a node from a path or an interval, but not
 * with the dynamic branching, whose tree depends on the incumbent.
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

#ifndef __PBeFS_Flowshop
 #define __PBeFS_Flowshop
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <climits>

#include <cstdint>

#include <numeric>

#include <optional>

#include <stdexcept>

#include <vector>

#include "PBeFS/bnb/BranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Flowshop ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the permutation flowshop of the least makespan

struct Flowshop {

 using value_type = int;          ///< the type of the values

 /// the bound of a node
 enum class Bound { one_machine , two_machine };

 /// the direction of the children of a node
 enum class Branching { forward , backward , alternate , dynamic };

 /// a node: the sequence, whose first l1 jobs are the head and whose last
 /// l2 the tail, the times of the head and of the tail on each machine,
 /// and the bound computed when it was generated (-1 if none)
 struct node_type {
  std::vector< std::uint16_t > perm;   ///< the jobs
  std::uint16_t l1 = 0;                ///< the length of the head
  std::uint16_t l2 = 0;                ///< ... and of the tail
  std::vector< int > front;            ///< the head on each machine
  std::vector< int > back;             ///< the tail on each machine
  int lb = -1;                         ///< its bound, if known
  };

 static constexpr Sense sense = Sense::minimize;

 int n = 0;                            ///< the jobs
 int m = 0;                            ///< the machines
 std::vector< int > p;                 ///< p[ k * n + j ], machine k
 Bound bound = Bound::two_machine;     ///< the bound
 Branching branching = Branching::forward;  ///< the direction
 bool sort_children = false;           ///< the children by their bound

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the processing time of the job j on the machine k

 [[nodiscard]] int time( int k , int j ) const {
  return( p[ std::size_t( k ) * std::size_t( n ) + std::size_t( j ) ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the instance taNNN of Taillard (1 to 120), generated from its seed by
 /// the generator of the paper (the uniform integers of [1, 99] drawn
 /// machine by machine, job by job)

 static Flowshop taillard( int id ) {
  static const long seeds[ 120 ] = {
   873654221 , 379008056 , 1866992158 , 216771124 , 495070989 ,
   402959317 , 1369363414 , 2021925980 , 573109518 , 88325120 ,
   587595453 , 1401007982 , 873136276 , 268827376 , 1634173168 ,
   691823909 , 73807235 , 1273398721 , 2065119309 , 1672900551 ,
   479340445 , 268827376 , 1958948863 , 918272953 , 555010963 ,
   2010851491 , 1519833303 , 1748670931 , 1923497586 , 1829909967 ,
   1328042058 , 200382020 , 496319842 , 1203030903 , 1730708564 ,
   450926852 , 1303135678 , 1273398721 , 587288402 , 248421594 ,
   1958948863 , 575633267 , 655816003 , 1977864101 , 93805469 ,
   1803345551 , 49612559 , 1899802599 , 2013025619 , 578962478 ,
   1539989115 , 691823909 , 655816003 , 1315102446 , 1949668355 ,
   1923497586 , 1805594913 , 1861070898 , 715643788 , 464843328 ,
   896678084 , 1179439976 , 1122278347 , 416756875 , 267829958 ,
   1835213917 , 1328833962 , 1418570761 , 161033112 , 304212574 ,
   1539989115 , 655816003 , 960914243 , 1915696806 , 2013025619 ,
   1168140026 , 1923497586 , 167698528 , 1528387973 , 993794175 ,
   450926852 , 1462772409 , 1021685265 , 83696007 , 508154254 ,
   1861070898 , 26482542 , 444956424 , 2115448041 , 118254244 ,
   471503978 , 1215892992 , 135346136 , 1602504050 , 160037322 ,
   551454346 , 519485142 , 383947510 , 1968171878 , 540872513 ,
   2013025619 , 475051709 , 914834335 , 810642687 , 1019331795 ,
   2056065863 , 1342855162 , 1325809384 , 1988803007 , 765656702 ,
   1368624604 , 450181436 , 1927888393 , 1759567256 , 606425239 ,
   19268348 , 1298201670 , 2041736264 , 379756761 , 28837162 };
  if( ( id < 1 ) || ( id > 120 ) )
   throw( std::invalid_argument( "Flowshop::taillard: no instance" ) );
  Flowshop f;
  f.n = id > 110 ? 500 : id > 90 ? 200 : id > 60 ? 100 : id > 30 ? 50 : 20;
  const int g = ( id - 1 ) / 10;       // the group of ten
  static const int machines[ 12 ] = { 5 , 10 , 20 , 5 , 10 , 20 , 5 , 10 ,
                                      20 , 10 , 20 , 20 };
  f.m = machines[ g ];
  long seed = seeds[ id - 1 ];
  f.p.resize( std::size_t( f.n ) * std::size_t( f.m ) );
  for( int k = 0 ; k < f.m ; ++k )
   for( int j = 0 ; j < f.n ; ++j ) {
    const long a = 16807 , b = 127773 , c = 2836 , mm = 2147483647;
    const long q = seed / b;
    seed = a * ( seed % b ) - q * c;
    if( seed < 0 )
     seed += mm;
    const float u = float( seed ) / float( mm );
    f.p[ std::size_t( k ) * std::size_t( f.n ) + std::size_t( j ) ] =
     1 + int( u * 99 );
    }
  f.prepare();
  return( f );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the optimal makespans of the instances ta001 to ta030

 static int taillard_optimum( int id ) {
  static const int opt[ 30 ] = {
   1278 , 1359 , 1081 , 1293 , 1235 , 1195 , 1234 , 1206 , 1230 , 1108 ,
   1582 , 1659 , 1496 , 1377 , 1419 , 1397 , 1484 , 1538 , 1593 , 1591 ,
   2297 , 2099 , 2326 , 2223 , 2291 , 2226 , 2273 , 2200 , 2237 , 2178 };
  return( ( id >= 1 ) && ( id <= 30 ) ? opt[ id - 1 ] : -1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// computes what the bounds need once: the least heads and tails, and the
 /// order of Johnson and the lags of each pair of machines; to be called
 /// after p is set

 void prepare( void ) {
  head.assign( std::size_t( m ) , INT_MAX );
  tail.assign( std::size_t( m ) , INT_MAX );
  for( int j = 0 ; j < n ; ++j ) {
   int h = 0;
   for( int k = 0 ; k < m ; ++k ) {
    head[ std::size_t( k ) ] = std::min( head[ std::size_t( k ) ] , h );
    h += time( k , j );
    }
   int t = 0;
   for( int k = m - 1 ; k >= 0 ; --k ) {
    tail[ std::size_t( k ) ] = std::min( tail[ std::size_t( k ) ] , t );
    t += time( k , j );
    }
   }
  pairs.clear();
  for( int a = 0 ; a + 1 < m ; ++a )
   for( int b = a + 1 ; b < m ; ++b ) {
    Pair x{ a , b , std::vector< int >( std::size_t( n ) , 0 ) , {} };
    for( int j = 0 ; j < n ; ++j )
     for( int k = a + 1 ; k < b ; ++k )
      x.lag[ std::size_t( j ) ] += time( k , j );
    x.order.resize( std::size_t( n ) );
    std::iota( x.order.begin() , x.order.end() , 0 );
    const auto t1 = [ & ]( int j ) {
     return( time( a , j ) + x.lag[ std::size_t( j ) ] ); };
    const auto t2 = [ & ]( int j ) {
     return( time( b , j ) + x.lag[ std::size_t( j ) ] ); };
    std::stable_sort( x.order.begin() , x.order.end() ,
                      [ & ]( int i , int j ) {
                       const bool fi = t1( i ) < t2( i ) ,
                                  fj = t1( j ) < t2( j );
                       if( fi != fj )
                        return( fi );
                       if( fi )
                        return( t1( i ) < t1( j ) );
                       return( t2( i ) > t2( j ) );
                       } );
    pairs.push_back( std::move( x ) );
    }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 node_type root( void ) const {
  node_type r;
  r.perm.resize( std::size_t( n ) );
  std::iota( r.perm.begin() , r.perm.end() , std::uint16_t( 0 ) );
  r.front.assign( std::size_t( m ) , 0 );
  r.back.assign( std::size_t( m ) , 0 );
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bound of a node (its makespan if it is complete)

 [[nodiscard]] int lower_bound( const node_type & x ) const {
  if( x.l1 + x.l2 == n )
   return( makespan( x ) );
  std::vector< int > f( x.front ) , b( x.back );
  if( x.l1 == 0 )
   f = head;
  if( x.l2 == 0 )
   b = tail;
  if( bound == Bound::one_machine ) {
   std::vector< int > r( std::size_t( m ) , 0 );
   for( int i = x.l1 ; i < n - x.l2 ; ++i )
    for( int k = 0 ; k < m ; ++k )
     r[ std::size_t( k ) ] += time( k , x.perm[ std::size_t( i ) ] );
   int t = f[ 0 ] + r[ 0 ];
   int lb = t + b[ 0 ];
   for( int k = 1 ; k < m ; ++k ) {
    t = std::max( t , f[ std::size_t( k ) ] + r[ std::size_t( k ) ] );
    lb = std::max( lb , t + b[ std::size_t( k ) ] );
    }
   return( lb );
   }
  std::vector< char > in( std::size_t( n ) , 0 );
  for( int i = x.l1 ; i < n - x.l2 ; ++i )
   in[ x.perm[ std::size_t( i ) ] ] = 1;
  int lb = 0;
  for( const Pair & q : pairs ) {
   int t0 = f[ std::size_t( q.a ) ] , t1 = f[ std::size_t( q.b ) ];
   for( const int j : q.order )
    if( in[ std::size_t( j ) ] ) {
     t0 += time( q.a , j );
     t1 = std::max( t1 , t0 + q.lag[ std::size_t( j ) ] ) + time( q.b , j );
     }
   lb = std::max( { lb , t1 + b[ std::size_t( q.b ) ] ,
                    t0 + b[ std::size_t( q.a ) ] } );
   }
  return( lb );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the makespan of a complete node

 [[nodiscard]] int makespan( const node_type & x ) const {
  int c = 0;
  for( int k = 0 ; k < m ; ++k )
   c = std::max( c , x.front[ std::size_t( k ) ] +
                     x.back[ std::size_t( k ) ] );
  return( c );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the makespan of a sequence of all the jobs

 [[nodiscard]] int makespan( const std::vector< std::uint16_t > & s ) const {
  std::vector< int > c( std::size_t( m ) , 0 );
  for( const auto j : s ) {
   c[ 0 ] += time( 0 , j );
   for( int k = 1 ; k < m ; ++k )
    c[ std::size_t( k ) ] = std::max( c[ std::size_t( k ) ] ,
                                      c[ std::size_t( k - 1 ) ] ) +
                            time( k , j );
   }
  return( c.back() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the child of x that schedules the job at the position i of x.perm
 /// (an unscheduled one) next to the head (forward) or to the tail

 [[nodiscard]] node_type make_child( const node_type & x , int i ,
                                     bool forward ) const {
  node_type c = x;
  c.lb = -1;
  if( forward ) {
   std::swap( c.perm[ c.l1 ] , c.perm[ std::size_t( i ) ] );
   const int j = c.perm[ c.l1++ ];
   c.front[ 0 ] += time( 0 , j );
   for( int k = 1 ; k < m ; ++k )
    c.front[ std::size_t( k ) ] = std::max( c.front[ std::size_t( k ) ] ,
                                            c.front[ std::size_t( k - 1 ) ] ) +
                                  time( k , j );
   }
  else {
   const std::size_t e = std::size_t( n - c.l2 - 1 );
   std::swap( c.perm[ e ] , c.perm[ std::size_t( i ) ] );
   const int j = c.perm[ e ];
   ++c.l2;
   c.back[ std::size_t( m - 1 ) ] += time( m - 1 , j );
   for( int k = m - 2 ; k >= 0 ; --k )
    c.back[ std::size_t( k ) ] = std::max( c.back[ std::size_t( k ) ] ,
                                           c.back[ std::size_t( k + 1 ) ] ) +
                                 time( k , j );
   }
  return( c );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the positions in x.perm of the unscheduled jobs, by increasing job

 [[nodiscard]] std::vector< int > free_positions( const node_type & x ) const {
  std::vector< int > v;
  for( int i = x.l1 ; i < n - x.l2 ; ++i )
   v.push_back( i );
  std::sort( v.begin() , v.end() , [ & ]( int a , int b ) {
   return( x.perm[ std::size_t( a ) ] < x.perm[ std::size_t( b ) ] );
   } );
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the children of x in the direction forward, with their bounds, in the
 /// order of the children

 [[nodiscard]] std::vector< node_type > children_of( const node_type & x ,
                                                     bool forward ) const {
  std::vector< node_type > v;
  for( const int i : free_positions( x ) ) {
   v.push_back( make_child( x , i , forward ) );
   v.back().lb = lower_bound( v.back() );
   }
  if( sort_children )
   std::stable_sort( v.begin() , v.end() ,
                     []( const node_type & a , const node_type & b ) {
                      return( a.lb < b.lb );
                      } );
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the direction of the children of x with the incumbent ub, true for
 /// forward, and the children in it, with their bounds

 [[nodiscard]] std::pair< bool , std::vector< node_type > > split(
                                                     const node_type & x ,
                                                     int ub ) const {
  const unsigned d = unsigned( x.l1 + x.l2 );
  switch( branching ) {
   case( Branching::forward ):
    return( std::pair( true , children_of( x , true ) ) );
   case( Branching::backward ):
    return( std::pair( false , children_of( x , false ) ) );
   case( Branching::alternate ):
    return( std::pair( d % 2 == 0 , children_of( x , d % 2 == 0 ) ) );
   default:
    break;
   }
  auto f = children_of( x , true );
  auto b = children_of( x , false );
  const auto score = [ & ]( const std::vector< node_type > & v ) {
   long alive = 0 , sum = 0;
   for( const auto & c : v ) {
    alive += c.lb < ub;
    sum += c.lb;
    }
   return( std::pair( alive , - sum ) );
   };
  if( score( b ) < score( f ) )
   return( std::pair( false , std::move( b ) ) );
  return( std::pair( true , std::move( f ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the evaluator of a worker, which keeps the incumbent of cutoff( z )

 struct Worker {

  const Flowshop * f;
  int ub = INT_MAX;

  void cutoff( int z ) { ub = z; }

  BasicEvaluation< int > evaluate( node_type & x ) const {
   BasicEvaluation< int > e;
   e.bound = x.lb >= 0 ? x.lb : f->lower_bound( x );
   x.lb = e.bound;
   if( x.l1 + x.l2 == f->n )
    e.feasible = e.bound;
   return( e );
   }

  template< class F >
  void branch( const node_type & x , F && sink ) const {
   for( auto & c : f->split( x , ub ).second ) {
    const int lb = c.lb;
    sink( std::move( c ) , lb );
    }
   }

  std::optional< BBChild< node_type , int > > child( const node_type & x ,
                                                     std::size_t k ) const {
   return( f->child( x , k ) );
   }

  };  // end( struct( Worker ) )

 Worker make_worker( unsigned ) const { return( Worker{ this } ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the child k of x, for a branching that does not depend on the
 /// incumbent

 [[nodiscard]] std::optional< BBChild< node_type , int > > child(
                                                     const node_type & x ,
                                                     std::size_t k ) const {
  if( branching == Branching::dynamic )
   throw( std::logic_error( "Flowshop::child: dynamic branching" ) );
  if( k >= std::size_t( n - x.l1 - x.l2 ) )
   return( std::nullopt );
  const unsigned d = unsigned( x.l1 + x.l2 );
  const bool fw = ( branching == Branching::forward ) ||
                  ( ( branching == Branching::alternate ) && ( d % 2 == 0 ) );
  if( sort_children ) {
   auto v = children_of( x , fw );
   const int lb = v[ k ].lb;
   return( BBChild< node_type , int >{ std::move( v[ k ] ) , lb } );
   }
  node_type c = make_child( x , free_positions( x )[ k ] , fw );
  c.lb = lower_bound( c );
  const int lb = c.lb;
  return( BBChild< node_type , int >{ std::move( c ) , lb } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of children of a node of depth d

 [[nodiscard]] std::size_t arity( unsigned d ) const {
  return( d < unsigned( n ) ? std::size_t( n ) - d : 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the sequence of a complete node

 [[nodiscard]] std::vector< std::uint16_t > sequence(
                                               const node_type & x ) const {
  return( x.perm );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 private:

 /// a pair of machines: the lags of the jobs between them and their order
 /// of Johnson
 struct Pair {
  int a , b;
  std::vector< int > lag;
  std::vector< int > order;
  };

 std::vector< int > head;          ///< the least heads
 std::vector< int > tail;          ///< the least tails
 std::vector< Pair > pairs;        ///< the pairs of machines

 };  // end( struct( Flowshop ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Flowshop */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Flowshop.h ---------------------------*/
/*--------------------------------------------------------------------------*/
