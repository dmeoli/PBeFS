/*--------------------------------------------------------------------------*/
/*------------------------ File ConcurrentClosed.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ConcurrentClosed of the sets of data
 * structures of the standard library and of Boost [see Structures.h]: the
 * NodeTable of the set split into 256 shards by the hash of the state,
 * each with its own lock, which is also the lock of the nodes of its
 * states, so that update() takes the lock of the shard of the state.
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

#ifndef __PBeFS_std_ConcurrentClosed
 #define __PBeFS_std_ConcurrentClosed
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <memory>

#include <mutex>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/core/SearchUtils.h"


/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ConcurrentClosed -------------------------*/
/*--------------------------------------------------------------------------*/
/// the closed list shared by the workers, in 256 shards each with its lock
/** The interface is the one described in ConcurrentClosed.h. */

template< Domain D , class Lock >
class ConcurrentClosed {

 public:

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type

 /// the nodes that a worker allocates
 struct Local {
  NodePool< NodeT > pool;                         ///< the nodes

  /// the number of the nodes in the table
  [[nodiscard]] std::size_t stored( void ) const { return( pool.size() ); }
  };

 /// a state generated, with its g, its h and its hash
 struct Successor {
  State state;                                    ///< the state
  Cost g;                                         ///< its g
  Cost h;                                         ///< its h
  std::uint64_t hash;                             ///< its hash
  std::uint8_t op = no_op;                        ///< the operator of t
  };

 /// constructor, taking the Domain (which must outlive the table)

 explicit ConcurrentClosed( const D & domain ) : f_domain( domain ) {}

 /// the lock of the node n, that of the shard of its state

 Lock & lock( const NodeT * n ) {
  return( shard( f_domain.hash( n->state ) ).lock );
  }

 /// nothing, since the table of a shard may grow under its lock

 void prefetch( std::uint64_t ) const noexcept {}

 /// records the path to t of cost g through parent and op [see
 /// ConcurrentClosed.h]

 NodeT * update( const State & t , std::uint64_t hash , Cost g , Cost h ,
                 NodeT * parent , std::uint8_t op , Local & local ,
                 std::uint64_t & duplicates ) {
  Shard & sd = shard( hash );
  std::lock_guard< Lock > guard( sd.lock );
  auto [ x , fresh ] = sd.table.find_or_insert( t , hash , [ & ] {
   NodeT * y = local.pool.allocate();
   *y = NodeT{ t , g , h , parent , false , false , op };
   return( y );
   } );
  if( fresh )
   return( x );
  if( g < x->g ) {
   x->g = g;
   x->parent = parent;
   x->op = op;
   x->closed = false;
   return( x );
   }
  ++duplicates;
  return( nullptr );
  }

 private:

 static constexpr unsigned c_shards = 256;         ///< shards of the table

 // a shard of the table and its lock
 struct alignas( cache_line ) Shard {
  Lock lock;
  NodeTable< NodeT > table;
  };

 // the shard of a hash value
 Shard & shard( std::uint64_t h ) { return( f_shards[ h >> 56 ] ); }

 const D & f_domain;                               ///< the Domain
 std::unique_ptr< Shard[] > f_shards{ new Shard[ c_shards ] };  ///< shards

 };  // end( class( ConcurrentClosed ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_std_ConcurrentClosed */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConcurrentClosed.h -----------------------*/
/*--------------------------------------------------------------------------*/
