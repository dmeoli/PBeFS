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
 * favor of the larger g, which reaches the goal earlier inside the last f
 * layer [see Burns et al., SoCS 2012], LIFO among the entries with the same
 * f and g in BucketOpenList, in no fixed order in the heaps. Neither
 * supports decrease-key: a node reached through a cheaper path is pushed
 * again, and an entry is stale when the node is closed or its g no longer
 * matches the g of the entry. The caller discards the stale
 * entries as they surface. A queue over NodeT has push( node ), which pushes
 * node with its current f and g, push( node , g ), which pushes it with the
 * given g (and \f$ f = g + h \f$) without reading node->g, push( node ,
 * g , key ), which pushes it with the given g and a key of the same type
 * that takes the place of f, such as the \f$ g + w h \f$ of the weighted
 * searches [see Weight], top(), the OpenEntry with the least f (or key) and
 * the largest g among those, top_f(), its f (or key), pop(), which removes
 * and returns it, empty() and size(), the number of entries, the stale
 * ones included.
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

#ifndef __PBeFS_OpenList
 #define __PBeFS_OpenList
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cstddef>

#include <limits>

#include <type_traits>

#include <utility>

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

//! [data-structures-sig-openkind]
/// the kind of OPEN list a search uses
enum class OpenKind {
 automatic ,  ///< buckets for integer costs, heap otherwise
 heap ,       ///< HeapOpenList
 buckets ,    ///< BucketOpenList (integer costs only)
 heap_low_g   ///< HeapOpenList with the ties to the smaller g [see
              ///< LowGOpenList]
 };
//! [data-structures-sig-openkind]

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
/*--------------------------- CLASS LowGOpenList ---------------------------*/
/*--------------------------------------------------------------------------*/
/// an OPEN list Base whose ties of f (or of the key) go to the smaller g
/** Base is a HeapOpenList, which breaks the ties in favor of the larger g;
 * LowGOpenList pushes each entry with the g mirrored, \f$ M - g \f$ with
 * integer costs (\f$ M \f$ half the largest cost) and \f$ - g \f$
 * otherwise, both exact, and mirrors it back in the entries it returns, so
 * that the ties go to the smaller g and the stale test is unchanged. With a
 * weight \f$ w > 1 \f$, the entries of equal key \f$ g + w h \f$ and smaller
 * g are those of smaller f: on grids this order usually expands fewer nodes
 * and finds cheaper paths than the larger g, which is better on the
 * sliding-tile puzzles. A BucketOpenList cannot be the Base, since it keeps
 * its layers dense in g. */

template< class Base >
class LowGOpenList {

 using NodeT = std::remove_pointer_t< decltype( std::declval< Base & >()
                                                .top().node ) >;
 using Cost = decltype( NodeT::g );

 public:

 /// pushes a node with its current f and g
 void push( NodeT * node ) { push( node , node->g ); }

 /// pushes a node with the given g (and f = g + h)
 void push( NodeT * node , Cost g ) { push( node , g , Cost( g + node->h ) ); }

 /// pushes a node with the given g and key
 void push( NodeT * node , Cost g , Cost key ) {
  f_base.push( node , mirror( g ) , key );
  }

 /// the entry with the least f (or key) and the smallest g among those
 [[nodiscard]] OpenEntry< NodeT > top( void ) {
  const OpenEntry< NodeT > e = f_base.top();
  return( OpenEntry< NodeT >{ e.node , mirror( e.g ) } );
  }

 /// the f (or key) of top()
 [[nodiscard]] Cost top_f( void ) { return( f_base.top_f() ); }

 /// removes and returns top()
 OpenEntry< NodeT > pop( void ) {
  const OpenEntry< NodeT > e = f_base.pop();
  return( OpenEntry< NodeT >{ e.node , mirror( e.g ) } );
  }

 /// true if there is no entry
 [[nodiscard]] bool empty( void ) const { return( f_base.empty() ); }

 /// the number of entries, the stale ones included
 [[nodiscard]] std::size_t size( void ) const { return( f_base.size() ); }

 private:

 /// g mirrored, an involution
 static Cost mirror( Cost g ) noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( std::numeric_limits< Cost >::max() / 2 - g ) );
  else
   return( - g );
  }

 Base f_base;                                          ///< the OPEN list

 };  // end( class( LowGOpenList ) )


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
