/*--------------------------------------------------------------------------*/
/*------------------------ File ExplicitAdjacency.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ExplicitAdjacency of the set of data
 * structures of ADS [see Structures.h], the arcs of a weighted graph, on
 * the GraphBuilder of ADS, frozen into its ArrayGraph before the first
 * visit; ExplicitGraph.h includes it in place of its own.
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

#ifndef __PBeFS_ads_ExplicitAdjacency
 #define __PBeFS_ads_ExplicitAdjacency
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <mutex>

#include <optional>

#include "ADS/graph/ArrayGraph.h"

#include "ADS/graph/GraphBuilder.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ExplicitAdjacency -------------------------*/
/*--------------------------------------------------------------------------*/
/// the arcs of a weighted graph, on the GraphBuilder and ArrayGraph of ADS
/** The interface is that of the ExplicitAdjacency of ExplicitGraph.h. The
 * arcs are collected by the GraphBuilder; the first for_each() after an
 * add() freezes them into an ArrayGraph, in compressed rows that keep the
 * order of insertion, under a lock taken only then, so that the searches
 * of many threads may start together. A loop ( u , u ) is not kept, as the
 * GraphBuilder wants two distinct nodes and a loop never shortens a path.
 * An add() must not run with a for_each(). */

template< class Cost >
class ExplicitAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor of n nodes and no arc

 explicit ExplicitAdjacency( std::size_t n )
  : f_n( n ) , f_b( int( n ) , true ) {}

 /// copy constructor: the arcs, to be frozen again

 ExplicitAdjacency( const ExplicitAdjacency & o )
  : f_n( o.f_n ) , f_b( o.f_b ) {}

 /// copy assignment: the arcs, to be frozen again

 ExplicitAdjacency & operator=( const ExplicitAdjacency & o ) {
  if( this != & o ) {
   f_n = o.f_n;
   f_b = o.f_b;
   f_g.reset();
   f_frozen.store( false , std::memory_order_relaxed );
   }
  return( *this );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the arc ( u , v ) of cost c

 void add( std::uint32_t u , std::uint32_t v , Cost c ) {
  if( u == v )
   return;
  f_b.addEdge( int( u ) , int( v ) , c );
  f_frozen.store( false , std::memory_order_relaxed );
  }

 /// calls f( v , c ) for each arc ( u , v ) of cost c

 template< class F >
 void for_each( std::uint32_t u , F && f ) const {
  if( ! f_frozen.load( std::memory_order_acquire ) )
   freeze();
  f_g->forEachEdge( int( u ) , [ & ]( int v , const Cost & c ) {
   f( std::uint32_t( v ) , c );
   } );
  }

 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept { return( f_n ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // builds the ArrayGraph of the arcs, once

 void freeze( void ) const {
  std::lock_guard< std::mutex > lock( f_lock );
  if( f_frozen.load( std::memory_order_relaxed ) )
   return;
  f_g.emplace( f_b.toArrayGraph() );
  f_frozen.store( true , std::memory_order_release );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 std::size_t f_n;                                ///< the number of nodes
 ADS::GraphBuilder< Cost > f_b;                  ///< the arcs, as added
 mutable std::optional< ADS::ArrayGraph< Cost > > f_g;
                                                 ///< the arcs, frozen
 mutable std::atomic< bool > f_frozen{ false };  ///< whether f_g is current
 mutable std::mutex f_lock;                      ///< the lock of freeze()

 };  // end( class( ExplicitAdjacency ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_ExplicitAdjacency */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ExplicitAdjacency.h ----------------------*/
/*--------------------------------------------------------------------------*/
