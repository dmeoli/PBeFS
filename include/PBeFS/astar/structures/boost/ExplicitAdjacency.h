/*--------------------------------------------------------------------------*/
/*------------------------ File ExplicitAdjacency.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ExplicitAdjacency of the set of data
 * structures of Boost [see Structures.h], the arcs of a weighted graph, in
 * the boost::compressed_sparse_row_graph built from them before the first
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

#ifndef __PBeFS_boost_ExplicitAdjacency
 #define __PBeFS_boost_ExplicitAdjacency
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <mutex>

#include <optional>

#include <utility>

#include <vector>

#include <boost/graph/compressed_sparse_row_graph.hpp>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ExplicitAdjacency -------------------------*/
/*--------------------------------------------------------------------------*/
/// the arcs of a weighted graph, in a compressed_sparse_row_graph of Boost
/** The interface is that of the ExplicitAdjacency of ExplicitGraph.h. The
 * arcs are collected in a vector; the first for_each() after an add()
 * builds from them the graph in compressed rows, under a lock taken only
 * then, so that the searches of many threads may start together. An add()
 * must not run with a for_each(). */

template< class Cost >
class ExplicitAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor of n nodes and no arc

 explicit ExplicitAdjacency( std::size_t n ) : f_n( n ) {}

 /// copy constructor: the arcs, to be built again

 ExplicitAdjacency( const ExplicitAdjacency & o )
  : f_n( o.f_n ) , v_arc( o.v_arc ) , v_cost( o.v_cost ) {}

 /// copy assignment: the arcs, to be built again

 ExplicitAdjacency & operator=( const ExplicitAdjacency & o ) {
  if( this != & o ) {
   f_n = o.f_n;
   v_arc = o.v_arc;
   v_cost = o.v_cost;
   f_g.reset();
   f_built.store( false , std::memory_order_relaxed );
   }
  return( *this );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the arc ( u , v ) of cost c

 void add( std::uint32_t u , std::uint32_t v , Cost c ) {
  v_arc.emplace_back( u , v );
  v_cost.push_back( c );
  f_built.store( false , std::memory_order_relaxed );
  }

 /// calls f( v , c ) for each arc ( u , v ) of cost c

 template< class F >
 void for_each( std::uint32_t u , F && f ) const {
  if( ! f_built.load( std::memory_order_acquire ) )
   build();
  const Csr & g = *f_g;
  for( auto [ e , end ] = boost::out_edges( u , g ) ; e != end ; ++e )
   f( std::uint32_t( boost::target( *e , g ) ) , g[ *e ] );
  }

 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept { return( f_n ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the graph in compressed rows, with the cost of each arc
 using Csr = boost::compressed_sparse_row_graph< boost::directedS ,
  boost::no_property , Cost , boost::no_property , std::uint32_t ,
  std::uint32_t >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // builds the graph of the arcs, once

 void build( void ) const {
  std::lock_guard< std::mutex > lock( f_lock );
  if( f_built.load( std::memory_order_relaxed ) )
   return;
  f_g.emplace( boost::edges_are_unsorted_multi_pass , v_arc.begin() ,
               v_arc.end() , v_cost.begin() , std::uint32_t( f_n ) );
  f_built.store( true , std::memory_order_release );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 std::size_t f_n;                                ///< the number of nodes
 std::vector< std::pair< std::uint32_t , std::uint32_t > > v_arc;
                                                 ///< the arcs, as added
 std::vector< Cost > v_cost;                     ///< their costs
 mutable std::optional< Csr > f_g;               ///< the arcs, built
 mutable std::atomic< bool > f_built{ false };   ///< whether f_g is current
 mutable std::mutex f_lock;                      ///< the lock of build()

 };  // end( class( ExplicitAdjacency ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_boost_ExplicitAdjacency */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ExplicitAdjacency.h ----------------------*/
/*--------------------------------------------------------------------------*/
