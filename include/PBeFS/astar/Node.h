/*--------------------------------------------------------------------------*/
/*------------------------------ File Node.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template Node, a node of the search graph, for
 * the class template NodePool, the arena the nodes are allocated from, and
 * for the infinite cost and the path to a goal.
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

#ifndef __PBeFS_Node
 #define __PBeFS_Node
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <limits>

#include <memory>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Structures.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS Node -------------------------------*/
/*--------------------------------------------------------------------------*/
/// a node of the search graph
/** There is one canonical Node per generated state. The field closed tells
 * whether the Node is in CLOSED, i.e., it has been expanded with its current
 * g; a Node reached again through a cheaper path is moved back to OPEN by
 * clearing it [see line 11 of Algorithm 1 of Fukunaga et al., 2017], while
 * expanded stays set, so that the re-expansions can be counted. The field
 * op is used by the OperatorDomain only, and fits in the padding. */

template< class State , class Cost >
struct Node {

 State state;     ///< the state
 Cost g;          ///< the cost of the best path found from the root
 Cost h;          ///< the heuristic estimate of the cost to a goal
 Node * parent;   ///< the predecessor on that path (nullptr at the root)
 bool closed;     ///< true if expanded with the current g
 bool expanded;   ///< true if ever expanded
 std::uint8_t op = no_op;  ///< the operator that generated the state

 /// the estimate \f$ f = g + h \f$ of the cost of a solution through the
 /// Node
 [[nodiscard]] Cost f( void ) const noexcept { return( g + h ); }

 };  // end( struct( Node ) )

#ifdef PBEFS_OWN_STRUCTURES

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS NodePool -----------------------------*/
/*--------------------------------------------------------------------------*/
/// an arena of nodes allocated in large chunks
/** Nodes are never freed one by one during the search, so allocating them
 * in blocks removes the cost of the allocator and keeps them contiguous in
 * memory; their addresses never change, which the parent pointers rely on.
 * The chunks of \f$ 2^{16} \f$ nodes take a few MB each and go on huge pages
 * [see HugePageAllocator in SearchUtils.h]. */

template< class NodeT >
class NodePool {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor, taking the number of nodes of a chunk

 explicit NodePool( std::size_t chunk = std::size_t( 1 ) << 16 )
  : f_chunk( chunk ) , f_used( 0 ) , f_size( 0 ) {}

 /// destructor: destroys the nodes and frees the chunks

 ~NodePool() {
  for( NodeT * c : v_chunks ) {
   std::destroy_n( c , f_chunk );
   HugePageAllocator< NodeT >().deallocate( c , f_chunk );
   }
  }

 NodePool( const NodePool & ) = delete;
 NodePool & operator=( const NodePool & ) = delete;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns a new (uninitialized) node

 [[nodiscard]] NodeT * allocate( void ) {
  if( v_chunks.empty() || ( f_used == f_chunk ) ) {
   NodeT * c = HugePageAllocator< NodeT >().allocate( f_chunk );
   std::uninitialized_default_construct_n( c , f_chunk );
   v_chunks.push_back( c );
   f_used = 0;
   }
  ++f_size;
  return( v_chunks.back() + f_used++ );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of nodes allocated so far

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_size );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the memory taken by the chunks, in bytes

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( v_chunks.size() * f_chunk * sizeof( NodeT ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 std::size_t f_chunk;                            ///< nodes per chunk
 std::size_t f_used;                   ///< nodes used in the last chunk
 std::size_t f_size;                             ///< nodes allocated
 std::vector< NodeT * > v_chunks;                          ///< the chunks

 };  // end( class( NodePool ) )

#endif  /* PBEFS_OWN_STRUCTURES */

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the largest value of a cost type, used as "no incumbent yet"

template< class Cost >
[[nodiscard]] constexpr Cost infinite_cost( void ) noexcept
{
 if constexpr( std::numeric_limits< Cost >::has_infinity )
  return( std::numeric_limits< Cost >::infinity() );
 else
  return( std::numeric_limits< Cost >::max() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the states from the root to goal, following the parent pointers

template< class NodeT >
[[nodiscard]] auto trace_path( const NodeT * goal )
{
 std::vector< decltype( goal->state ) > path;
 for( const NodeT * n = goal ; n ; n = n->parent )
  path.push_back( n->state );
 std::reverse( path.begin() , path.end() );
 return( path );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

// the structures of another set, in place of those of the library
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/NodePool.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/NodePool.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/NodePool.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Node */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Node.h -----------------------------*/
/*--------------------------------------------------------------------------*/
