/*--------------------------------------------------------------------------*/
/*------------------------------ File Node.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template Node, a node of the search graph, for
 * the class template NodePool, the arena the nodes are allocated from,
 * taken from the set of data structures chosen [see Structures.h], and for
 * the infinite cost and the path to a goal. A NodePool< NodeT > has a
 * constructor taking the number of nodes of a chunk, allocate(), which
 * returns a new uninitialized node, size() and bytes(); the nodes are never
 * freed one by one, and their addresses never change, which the parent
 * pointers rely on.
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

// the structures of the chosen set [see Structures.h]
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
