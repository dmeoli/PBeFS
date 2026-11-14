/*--------------------------------------------------------------------------*/
/*---------------------------- File OpenList.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the priority queues of the OPEN list of the set of data
 * structures of Boost [see Structures.h], HeapOpenList and BucketOpenList,
 * on boost::heap::d_ary_heap; OpenList.h includes it in place of its own,
 * with the same order of the entries but for the ties of f and g.
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

#ifndef __PBeFS_boost_OpenList
 #define __PBeFS_boost_OpenList
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <boost/heap/d_ary_heap.hpp>

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
/// the OPEN list on boost::heap::d_ary_heap, a d-ary heap
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
  const Item x = f_heap.top();
  f_heap.pop();
  return( OpenEntry< NodeT >{ x.node , x.g } );
  }

 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept { return( f_heap.empty() ); }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_heap.size() );
  }

 private:

 /// an item of the heap
 struct Item {
  Cost f;
  Cost g;
  NodeT * node;
  };

 /// the order of boost::heap (a max-heap): a below b if its f is larger,
 /// or if its g is smaller for equal f
 struct Below {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f > b.f ) || ( ( a.f == b.f ) && ( a.g < b.g ) ) );
   }
  };

 boost::heap::d_ary_heap< Item , boost::heap::arity< Arity > ,
                          boost::heap::compare< Below > > f_heap;

 };  // end( class( HeapOpenList ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BucketOpenList --------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list on boost::heap::d_ary_heap, a 4-ary heap
/** Boost has no bucket queue, and hence the heap is used for integer
 * costs too. */
template< class NodeT >
class BucketOpenList {

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
  const Item x = f_heap.top();
  f_heap.pop();
  return( OpenEntry< NodeT >{ x.node , x.g } );
  }

 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept { return( f_heap.empty() ); }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_heap.size() );
  }

 private:

 /// an item of the heap
 struct Item {
  Cost f;
  Cost g;
  NodeT * node;
  };

 /// the order of boost::heap (a max-heap): a below b if its f is larger,
 /// or if its g is smaller for equal f
 struct Below {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f > b.f ) || ( ( a.f == b.f ) && ( a.g < b.g ) ) );
   }
  };

 boost::heap::d_ary_heap< Item , boost::heap::arity< 4 > ,
                          boost::heap::compare< Below > > f_heap;

 };  // end( class( BucketOpenList ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_boost_OpenList */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File OpenList.h ---------------------------*/
/*--------------------------------------------------------------------------*/
