/*--------------------------------------------------------------------------*/
/*--------------------------- File SlowDomain.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SlowDomain, a Domain [see Domain.h]
 * whose expansions are made artificially slow, to compare the searches when
 * the expansion dominates their cost.
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
  std::uint64_t i = 7;
  for( unsigned long k = 0 ; k < f_work ; ++k )
   i = ( 11 * i ) % 9999943;
  c_sink += i;
  f_d.for_each_successor( s , h , std::forward< F >( f ) );
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

 const D & f_d;                                    ///< the Domain
 unsigned long f_work;                             ///< iterations per call
 inline static thread_local std::uint64_t c_sink = 0;   ///< keeps the work

 };  // end( class( SlowDomain ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SlowDomain */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SlowDomain.h --------------------------*/
/*--------------------------------------------------------------------------*/
