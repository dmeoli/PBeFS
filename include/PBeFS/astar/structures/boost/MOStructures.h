/*--------------------------------------------------------------------------*/
/*--------------------------- File MOStructures.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the data structures of the multi-objective searches of
 * the set of data structures of Boost [see Structures.h]: MOHeap, the OPEN
 * list, on boost::heap::d_ary_heap, and LexSet, an ordered set of vectors,
 * on boost::container::set with an AVL tree; MOStructures.h includes it in
 * place of its own.
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

#ifndef __PBeFS_boost_MOStructures
 #define __PBeFS_boost_MOStructures
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <functional>

#include <iterator>

#include <boost/container/set.hpp>

#include <boost/heap/d_ary_heap.hpp>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS MOHeap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list of the multi-objective searches, on
/// boost::heap::d_ary_heap
/** The interface is the one described in MOStructures.h. */

template< class Item , class Better >
class MOHeap {

 // the order of boost::heap (a max-heap): a below b if b is better
 struct Below {
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

 boost::heap::d_ary_heap< Item , boost::heap::arity< 4 > ,
                          boost::heap::compare< Below > > f_heap;

 };  // end( class( MOHeap ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LexSet ------------------------------*/
/*--------------------------------------------------------------------------*/
/// an ordered set of vectors, on boost::container::set with an AVL tree
/** The interface is the one described in MOStructures.h; the keys up to a
 * vector are visited in increasing order. */

template< class Key >
class LexSet {

 using Set = boost::container::set< Key , std::less< Key > , void ,
              typename boost::container::tree_assoc_options<
               boost::container::tree_type<
                boost::container::avl_tree > >::type >;

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

 Set f_set;                                       ///< the set

 };  // end( class( LexSet ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_boost_MOStructures */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MOStructures.h -------------------------*/
/*--------------------------------------------------------------------------*/
