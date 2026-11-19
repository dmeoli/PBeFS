/*--------------------------------------------------------------------------*/
/*--------------------------- File MOStructures.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the data structures of the multi-objective searches of
 * the set of data structures of the standard library [see Structures.h]:
 * MOHeap, the OPEN list, on std::priority_queue, and LexSet, an ordered set
 * of vectors, on std::set; MOStructures.h includes it in place of its own.
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

#ifndef __PBeFS_std_MOStructures
 #define __PBeFS_std_MOStructures
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <iterator>

#include <queue>

#include <set>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS MOHeap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list of the multi-objective searches, on std::priority_queue
/** The interface is the one described in MOStructures.h. */

template< class Item , class Better >
class MOHeap {

 // the order of std::priority_queue, whose top is the largest item
 struct Worse {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( Better()( b , a ) );
   }
  };

 public:

 void push( const Item & x ) { f_heap.push( x ); }

 [[nodiscard]] const Item & top( void ) const { return( f_heap.top() ); }

 Item pop( void ) {
  Item x = f_heap.top();
  f_heap.pop();
  return( x );
  }

 [[nodiscard]] bool empty( void ) const noexcept {
  return( f_heap.empty() );
  }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_heap.size() );
  }

 private:

 std::priority_queue< Item , std::vector< Item > , Worse > f_heap;
                                                        ///< the heap

 };  // end( class( MOHeap ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LexSet ------------------------------*/
/*--------------------------------------------------------------------------*/
/// an ordered set of vectors, on std::set
/** The interface is the one described in MOStructures.h; the keys up to a
 * vector are visited in increasing order. */

template< class Key >
class LexSet {

 public:

 void insert( const Key & k ) { f_set.insert( k ); }

 void erase( const Key & k ) { f_set.erase( k ); }

 [[nodiscard]] std::size_t size( void ) const { return( f_set.size() ); }

 template< class P >
 [[nodiscard]] bool any_upto( const Key & b , P && p ) const {
  for( auto it = f_set.begin() , e = f_set.upper_bound( b ) ; it != e ;
       ++it )
   if( p( *it ) )
    return( true );
  return( false );
  }

 [[nodiscard]] const Key * floor( const Key & b ) const {
  auto it = f_set.upper_bound( b );
  return( it == f_set.begin() ? nullptr : & *std::prev( it ) );
  }

 template< class P >
 void erase_from( const Key & b , P && p , bool run ) {
  for( auto it = f_set.lower_bound( b ) ; it != f_set.end() ; )
   if( p( *it ) )
    it = f_set.erase( it );
   else
    if( run )
     return;
    else
     ++it;
  }

 private:

 std::set< Key > f_set;                           ///< the set

 };  // end( class( LexSet ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_std_MOStructures */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MOStructures.h -------------------------*/
/*--------------------------------------------------------------------------*/
