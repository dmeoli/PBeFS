/*--------------------------------------------------------------------------*/
/*---------------------------- File OpenList.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the priority queues of the OPEN list, taken from the set
 * of data structures chosen [see Structures.h]: the class templates
 * HeapOpenList, a heap for any cost type, and BucketOpenList, a two-level
 * bucket queue for integer costs, and their entries.
 *
 * Both queues order the entries by increasing f and break the ties in
 * favour of the larger g, which reaches the goal earlier inside the last f
 * layer [see Burns et al., SoCS 2012], LIFO among the entries with the same
 * f and g. Neither supports decrease-key: a node reached through a cheaper
 * path is pushed again, and an entry is stale when the node is closed or its
 * g no longer matches the g of the entry. The caller discards the stale
 * entries as they surface. A queue over NodeT has push( node ), which pushes
 * node with its current f and g, push( node , g ), which pushes it with the
 * given g (and \f$ f = g + h \f$) without reading node->g, top(), the
 * OpenEntry with the least f and the largest g among those, top_f(), its f,
 * pop(), which removes and returns it, empty() and size(), the number of
 * entries, the stale ones included.
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

#ifndef __PBeFS_OpenList
 #define __PBeFS_OpenList
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cstddef>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the kind of OPEN list a search uses
enum class OpenKind {
 automatic ,  ///< buckets for integer costs, heap otherwise
 heap ,       ///< HeapOpenList
 buckets      ///< BucketOpenList (integer costs only)
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an entry of the OPEN list: the node and the g it had when pushed

template< class NodeT >
struct OpenEntry {

 NodeT * node;               ///< the node
 decltype( NodeT::g ) g;     ///< the g of the node when it was pushed

 /// true if the node has been closed or improved since it was pushed
 [[nodiscard]] bool stale( void ) const noexcept {
  return( node->closed || ( node->g != g ) );
  }

 };  // end( struct( OpenEntry ) )


/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

// the structures of the chosen set [see Structures.h]
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/OpenList.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/OpenList.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/OpenList.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_OpenList */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File OpenList.h ---------------------------*/
/*--------------------------------------------------------------------------*/
