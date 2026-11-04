/*--------------------------------------------------------------------------*/
/*------------------------ File ConcurrentClosed.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ConcurrentClosed of the set of data
 * structures of ADS [see Structures.h]: the nodes are kept in the
 * ConcurrentFlatTable of ADS as a set of pointers compared by their
 * states, a word of 8 bytes for each, which a lookup reads and an
 * insertion fills by a compare-and-swap with no lock. The fields that
 * change are protected by 256 locks, the lock of a node chosen by its
 * address, and update() takes it only when the state is already in the
 * table and g is smaller than its g read without the lock: since the g of
 * a node only decreases, a g read without the lock is never smaller than
 * the one under it, and a duplicate, the most frequent case, takes no lock.
 * The g of a node is read without the lock, and written under it, as an
 * std::atomic_ref. A worker fills a spare node with the state before
 * inserting it, and keeps it for the next update() if the state was in
 * the table, so that the node is allocated once per state inserted.
 * prefetch() fetches the line of the tags and that of the word of the
 * place of a hash, as the ConcurrentFlatTable does.
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

#ifndef __PBeFS_ads_ConcurrentClosed
 #define __PBeFS_ads_ConcurrentClosed
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <memory>

#include <mutex>

#include <variant>

#include "ADS/concurrent/ConcurrentFlatTable.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

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
/// the closed list shared by the workers, a lock-free set of the nodes
/** The interface is the one described in ConcurrentClosed.h. */

template< Domain D , class Lock >
class ConcurrentClosed {

 public:

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type

 static_assert( alignof( Cost ) >=
                std::atomic_ref< Cost >::required_alignment );

 /// the nodes that a worker allocates, and its spare node
 struct Local {
  NodePool< NodeT > pool;                         ///< the nodes
  NodeT * spare = nullptr;                        ///< not in the table

  /// the number of the nodes in the table
  [[nodiscard]] std::size_t stored( void ) const {
   return( pool.size() - ( spare ? 1 : 0 ) );
   }
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

 explicit ConcurrentClosed( const D & domain )
  : f_domain( domain ) , f_table( 0 , ByState{ &domain } ) {}

 /// the lock of the node n, chosen by its address

 Lock & lock( const NodeT * n ) {
  const auto a = std::uint64_t( reinterpret_cast< std::uintptr_t >( n ) );
  return( f_stripes[ ( ( a >> 4 ) * 0x9E3779B97F4A7C15ULL ) >> 56 ].lock );
  }

 /// the tags and the word that the lookup of the hash reads first,
 /// fetched with no wait

 void prefetch( std::uint64_t hash ) const { f_table.prefetch( hash ); }

 /// records the path to t of cost g through parent and op [see
 /// ConcurrentClosed.h]

 NodeT * update( const State & t , std::uint64_t hash , Cost g , Cost h ,
                 NodeT * parent , std::uint8_t op , Local & local ,
                 std::uint64_t & duplicates ) {
  if( ! local.spare )
   local.spare = local.pool.allocate();
  *local.spare = NodeT{ t , g , h , parent , false , false , op };
  const auto [ x , fresh ] = f_table.findOrInsertKey( local.spare , hash );
  if( fresh ) {
   local.spare = nullptr;
   return( x );
   }
  if( g < std::atomic_ref< Cost >( x->g ).load(
           std::memory_order_relaxed ) ) {
   std::lock_guard< Lock > guard( lock( x ) );
   if( g < x->g ) {
    std::atomic_ref< Cost >( x->g ).store( g , std::memory_order_relaxed );
    x->parent = parent;
    x->op = op;
    x->closed = false;
    return( x );
    }
   }
  ++duplicates;
  return( nullptr );
  }

 private:

 // the hash of the state of a node, well spread by the Domain
 struct ByState {
  using is_avalanching = void;
  const D * d;
  std::uint64_t operator()( const NodeT * n ) const {
   return( d->hash( n->state ) );
   }
  };

 // the equality of the states of two nodes
 struct SameState {
  bool operator()( const NodeT * a , const NodeT * b ) const {
   return( a->state == b->state );
   }
  };

 static constexpr unsigned c_stripes = 256;        ///< locks of the nodes

 // a lock on a line of its own
 struct alignas( cache_line ) Stripe {
  Lock lock;
  };

 using Table = ADS::ConcurrentFlatTable< NodeT * , std::monostate , ByState ,
                                         SameState >;
 static_assert( Table::packed );

 const D & f_domain;                               ///< the Domain
 Table f_table;                                    ///< the nodes
 std::unique_ptr< Stripe[] > f_stripes{ new Stripe[ c_stripes ] }; ///< locks

 };  // end( class( ConcurrentClosed ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_ConcurrentClosed */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConcurrentClosed.h -----------------------*/
/*--------------------------------------------------------------------------*/
