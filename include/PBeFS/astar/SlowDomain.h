/*--------------------------------------------------------------------------*/
/*--------------------------- File SlowDomain.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SlowDomain, a Domain [see Domain.h]
 * whose expansions are made artificially slow, to compare the searches when
 * the expansion dominates their cost, and for the class template
 * SlowEdgeDomain, an EdgeDomain whose edges are made artificially slow to
 * evaluate, cheap or expensive, and possibly missing, to compare the
 * searches by edges with those by states.
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

#ifndef __PBeFS_SlowDomain
 #define __PBeFS_SlowDomain
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <utility>

#include "PBeFS/astar/Domain.h"

#include "ADS/dictionary/RandomHashing.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS SlowDomain ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a Domain whose expansions cost a given amount of extra work
/** Wraps the Domain D and adds to each expansion work iterations of the
 * computation \f$ j = 11 i \bmod 9999943 \f$, \f$ i = j \f$, which keeps a
 * core busy with no memory access, as in the experiment of Jinnai and
 * Fukunaga (JAIR 2017, section 3.1.1) on the effect of slow expansions. It
 * is meant to compare the searches when the expansion dominates, as in
 * planning with an expensive heuristic or in motion planning. */

template< class D >
class SlowDomain {

 public:

 using state_type = state_t< D >;                  ///< the state type
 using cost_type = cost_t< D >;                    ///< the cost type

 /// constructor, taking the Domain (which must outlive this) and the work

 SlowDomain( const D & d , unsigned long work ) : f_d( d ) , f_work( work ) {}

 [[nodiscard]] state_type initial_state( void ) const {
  return( f_d.initial_state() );
  }

 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( f_d.is_goal( s ) );
  }

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d.hash( s ) );
  }

 [[nodiscard]] cost_type heuristic( const state_type & s ) const {
  return( f_d.heuristic( s ) );
  }

 /// the successors of D, after work iterations of busy work

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f ) const
 {
  work();
  f_d.for_each_successor( s , h , std::forward< F >( f ) );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type h ,
                          std::uint8_t op , F && f ) const
  requires OperatorDomain< D > {
  work();
  f_d.for_each_successor( s , h , op , std::forward< F >( f ) );
  }

 // the other methods of D, when D has them

 [[nodiscard]] std::size_t num_features( void ) const
  requires FeatureDomain< D > { return( f_d.num_features() ); }

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const
  requires FeatureDomain< D > {
  f_d.for_each_feature( s , std::forward< F >( f ) );
  }

 [[nodiscard]] std::size_t num_abstract_features( void ) const
  requires AbstractFeatureDomain< D > {
  return( f_d.num_abstract_features() );
  }

 [[nodiscard]] std::size_t abstract_feature( std::size_t i ) const
  requires AbstractFeatureDomain< D > { return( f_d.abstract_feature( i ) ); }

 [[nodiscard]] std::uint64_t abstract_state( const state_type & s ) const
  requires AbstractStateDomain< D > { return( f_d.abstract_state( s ) ); }

 [[nodiscard]] std::uint64_t perfect_hash( const state_type & s ) const
  requires PerfectHashDomain< D > { return( f_d.perfect_hash( s ) ); }

 [[nodiscard]] std::size_t num_nblocks( void ) const
  requires NBlockDomain< D > { return( f_d.num_nblocks() ); }

 [[nodiscard]] std::size_t nblock( const state_type & s ) const
  requires NBlockDomain< D > { return( f_d.nblock( s ) ); }

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const
  requires NBlockDomain< D > {
  f_d.for_each_nblock_successor( b , std::forward< F >( f ) );
  }

 [[nodiscard]] cost_type pairwise_heuristic( const state_type & s ,
                                             const state_type & t ) const
  requires PairwiseDomain< D > { return( f_d.pairwise_heuristic( s , t ) ); }

 private:

 // the busy work of an expansion
 void work( void ) const {
  std::uint64_t i = 7;
  for( unsigned long k = 0 ; k < f_work ; ++k )
   i = ( 11 * i ) % 9999943;
  c_sink += i;
  }

 const D & f_d;                                    ///< the Domain
 unsigned long f_work;                             ///< iterations per call
 inline static thread_local std::uint64_t c_sink = 0;   ///< keeps the work

 };  // end( class( SlowDomain ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS EdgeWork ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the work of the edges of a SlowEdgeDomain

//! [shared-memory-sig-edgework]
struct EdgeWork {

 unsigned long cheap = 0;        ///< iterations of a cheap edge
 unsigned long expensive = 0;    ///< iterations of an expensive edge
 double expensive_share = 0;     ///< the share of the expensive edges
 double missing_share = 0;       ///< the share of the edges that, once
                                 ///< evaluated, turn out not to exist
 std::uint64_t seed = 1;         ///< the seed of the two choices

 };  // end( struct( EdgeWork ) )
//! [shared-memory-sig-edgework]

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SlowEdgeDomain --------------------------*/
/*--------------------------------------------------------------------------*/
/// an EdgeDomain whose edges cost a given amount of work to evaluate
/** Wraps the Domain D: the actions of a state s are the indices
 * \f$ 0 , \dots , k - 1 \f$ of its successors in the order in which D
 * gives them, and evaluating the edge \f$ ( s , a ) \f$ costs the work of
 * SlowDomain [see there] for cheap or for expensive iterations. An edge is
 * expensive, or missing, if a hash of the hash of s, of a and of the seed,
 * mapped to \f$ [ 0 , 1 ) \f$, falls below expensive_share, or below
 * missing_share for an independent hash, so that a given edge is always of
 * the same kind; a missing edge is reported by lazy_edge with its cost, and
 * found missing only by evaluate_edge, as a collision check would do. The
 * expansion of a state, for_each_successor, evaluates all its edges, and
 * then costs what the search by edges pays for them; the graph is that of
 * D without the missing edges, on which the heuristics of D stay admissible
 * and consistent. With cheap and expensive work in the ratio of 1 to 30 on
 * half of the edges, this is the regime of the 2D grids of GePA*SE
 * (section "Evaluation" of Mukherjee and Likhachev, SoCS 2023). */

template< class D >
class SlowEdgeDomain {

 public:

 using state_type = state_t< D >;                  ///< the state type
 using cost_type = cost_t< D >;                    ///< the cost type

 /// constructor, taking the Domain (which must outlive this) and the work

 SlowEdgeDomain( const D & d , const EdgeWork & w ) : f_d( d ) , f_w( w ) {}

 [[nodiscard]] state_type initial_state( void ) const {
  return( f_d.initial_state() );
  }

 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( f_d.is_goal( s ) );
  }

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d.hash( s ) );
  }

 [[nodiscard]] cost_type heuristic( const state_type & s ) const {
  return( f_d.heuristic( s ) );
  }

 /// the existing successors of D, each after the work of its edge

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f ) const
 {
  std::size_t a = 0;
  f_d.for_each_successor( s , h , [ & ]( const state_type & t , cost_type c ,
                                         cost_type ht ) {
   work( s , a );
   if( ! missing( s , a ) )
    f( t , c , ht );
   ++a;
   } );
  }

 /// the indices of the successors of s in D

 template< class F >
 void for_each_action( const state_type & s , F && f ) const {
  std::size_t a = 0;
  f_d.for_each_successor( s , f_d.heuristic( s ) ,
                          [ & ]( const state_type & , cost_type ,
                                 cost_type ) { f( a++ ); } );
  }

 /// the successor a of s in D, after the work of the edge, if not missing

 template< class F >
 void evaluate_edge( const state_type & s , cost_type h , std::size_t a ,
                     F && f ) const {
  work( s , a );
  if( ! missing( s , a ) )
   nth( s , h , a , f );
  }

 /// the successor a of s in D, with no work, missing or not

 template< class F >
 void lazy_edge( const state_type & s , cost_type h , std::size_t a ,
                 F && f ) const { nth( s , h , a , f ); }

 /// true if the edge ( s , a ) is expensive

 [[nodiscard]] bool expensive_edge( const state_type & s ,
                                    std::size_t a ) const {
  return( draw( s , a , 0 ) < f_w.expensive_share );
  }

 /// true if the edge ( s , a ) does not exist, once evaluated

 [[nodiscard]] bool missing( const state_type & s , std::size_t a ) const {
  return( ( f_w.missing_share > 0 ) &&
          ( draw( s , a , 1 ) < f_w.missing_share ) );
  }

 [[nodiscard]] cost_type pairwise_heuristic( const state_type & s ,
                                             const state_type & t ) const
  requires PairwiseDomain< D > { return( f_d.pairwise_heuristic( s , t ) ); }

 private:

 // a number in [ 0 , 1 ) fixed by s, a, the seed and the choice k
 double draw( const state_type & s , std::size_t a , unsigned k ) const {
  const std::uint64_t x = ADS::splitmix64( f_d.hash( s ) ^
                          ADS::splitmix64( ( std::uint64_t( a ) << 1 | k ) ^
                                           f_w.seed ) );
  return( double( x >> 11 ) * 0x1.0p-53 );
  }

 // the successor a of s in D
 template< class F >
 void nth( const state_type & s , cost_type h , std::size_t a ,
           F & f ) const {
  std::size_t i = 0;
  f_d.for_each_successor( s , h , [ & ]( const state_type & t , cost_type c ,
                                         cost_type ht ) {
   if( i++ == a )
    f( t , c , ht );
   } );
  }

 // the busy work of the edge ( s , a ), as that of SlowDomain
 void work( const state_type & s , std::size_t a ) const {
  const unsigned long n = expensive_edge( s , a ) ? f_w.expensive
                                                  : f_w.cheap;
  std::uint64_t i = 7;
  for( unsigned long k = 0 ; k < n ; ++k )
   i = ( 11 * i ) % 9999943;
  c_sink += i;
  }

 const D & f_d;                                    ///< the Domain
 EdgeWork f_w;                                     ///< the work
 inline static thread_local std::uint64_t c_sink = 0;   ///< keeps the work

 };  // end( class( SlowEdgeDomain ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SlowDomain */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SlowDomain.h --------------------------*/
/*--------------------------------------------------------------------------*/
