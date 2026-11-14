/*--------------------------------------------------------------------------*/
/*---------------------------- File OpenList.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the priority queues of the OPEN list of the set of data
 * structures of ADS [see Structures.h], HeapOpenList and BucketOpenList, on
 * the DaryHeap and the BucketQueue of ADS; OpenList.h includes it in place
 * of its own, with the same order of the entries but for the ties of f and
 * g.
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

#ifndef __PBeFS_ads_OpenList
 #define __PBeFS_ads_OpenList
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <type_traits>

#include "ADS/priority/BucketQueue.h"

#include "ADS/priority/DaryHeap.h"

#include "PBeFS/astar/OpenList.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS HeapOpenList ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list on the DaryHeap of ADS
template< class NodeT , unsigned Arity = 4 >
class HeapOpenList {

 using Cost = decltype( NodeT::g );

 public:

 /// pushes a node with its current f and g

 void push( NodeT * node ) { push( node , node->g ); }

 /// pushes a node with the given g (and f = g + h), not reading node->g

 void push( NodeT * node , Cost g ) {
  f_heap.push( Item{ g + node->h , g , node } );
  }

 /// pushes a node with the given g and the given key in place of f [see
 /// Weight]

 void push( NodeT * node , Cost g , Cost key ) {
  f_heap.push( Item{ key , g , node } );
  }

 /// the entry with the smallest f (and the largest g among those)

 [[nodiscard]] OpenEntry< NodeT > top( void ) const {
  return( OpenEntry< NodeT >{ f_heap.top().node , f_heap.top().g } );
  }

 /// the f of top(), or its key if it was pushed with one

 [[nodiscard]] Cost top_f( void ) const { return( f_heap.top().f ); }

 /// removes and returns top()

 OpenEntry< NodeT > pop( void ) {
  const Item x = f_heap.pop();
  return( OpenEntry< NodeT >{ x.node , x.g } );
  }

 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept {
  return( f_heap.empty() );
  }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_heap.size() );
  }

 private:

 /// an item of the heap
 //! [data-structures-heap-order]
 struct Item {
  Cost f;
  Cost g;
  NodeT * node;
  };

 /// the order of the heap: smaller f first, then larger g
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.g > b.g ) ) );
   }
  };
 //! [data-structures-heap-order]

 ADS::DaryHeap< Item , Better , Arity > f_heap;         ///< the heap

 };  // end( class( HeapOpenList ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BucketOpenList --------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list on the BucketQueue of ADS, for integer costs
template< class NodeT >
class BucketOpenList {

 using Cost = decltype( NodeT::g );

 static_assert( std::is_integral_v< Cost > ,
                "BucketOpenList requires an integer cost type" );

 public:

 /// pushes a node with its current f and g

 void push( NodeT * node ) { push( node , node->g ); }

 /// pushes a node with the given g (and f = g + h), not reading node->g

 void push( NodeT * node , Cost g ) { f_queue.push( node , g + node->h , g ); }

 /// pushes a node with the given g and the given key in place of f [see
 /// Weight]

 void push( NodeT * node , Cost g , Cost key ) {
  f_queue.push( node , key , g );
  }

 /// the entry with the smallest f (and the largest g among those)

 [[nodiscard]] OpenEntry< NodeT > top( void ) {
  return( OpenEntry< NodeT >{ f_queue.top() , f_queue.topG() } );
  }

 /// the f of top(), or its key if it was pushed with one

 [[nodiscard]] Cost top_f( void ) { return( f_queue.topF() ); }

 /// removes and returns top()

 OpenEntry< NodeT > pop( void ) {
  const auto e = f_queue.popEntry();
  return( OpenEntry< NodeT >{ e.item , e.g } );
  }

 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept {
  return( f_queue.empty() );
  }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_queue.size() );
  }

 private:

 ADS::BucketQueue< NodeT * , Cost > f_queue;               ///< the queue

 };  // end( class( BucketOpenList ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_OpenList */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File OpenList.h ---------------------------*/
/*--------------------------------------------------------------------------*/
