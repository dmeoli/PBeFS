/*--------------------------------------------------------------------------*/
/*--------------------------- File MOStructures.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the data structures of the multi-objective searches of
 * the set of data structures of ADS [see Structures.h]: MOHeap, the OPEN
 * list, on the DaryHeap of ADS, and LexSet, an ordered set of vectors, on
 * the AVLTree of ADS; MOStructures.h includes it in place of its own.
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

#ifndef __PBeFS_ads_MOStructures
 #define __PBeFS_ads_MOStructures
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <vector>

#include "ADS/priority/DaryHeap.h"

#include "ADS/tree/AVLTree.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS MOHeap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list of the multi-objective searches, on the DaryHeap of ADS
/** The interface is the one described in MOStructures.h. */

template< class Item , class Better >
class MOHeap {

 public:

 void push( const Item & x ) { f_heap.push( x ); }

 [[nodiscard]] const Item & top( void ) const { return( f_heap.top() ); }

 Item pop( void ) { return( f_heap.pop() ); }

 [[nodiscard]] bool empty( void ) const noexcept {
  return( f_heap.empty() );
  }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_heap.size() );
  }

 private:

 ADS::DaryHeap< Item , Better , 4 > f_heap;             ///< the heap

 };  // end( class( MOHeap ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LexSet ------------------------------*/
/*--------------------------------------------------------------------------*/
/// an ordered set of vectors, on the AVLTree of ADS
/** The interface is the one described in MOStructures.h. The keys up to a
 * vector are visited by the descent of Algorithm 2 of Ren et al. (SoCS
 * 2022): the root of a subtree first, then its left subtree, and its right
 * subtree only if the root precedes the vector. */

template< class Key >
class LexSet {

 using Tree = ADS::AVLTree< Key , char >;
 using Shape = typename Tree::Shape;

 public:

 void insert( const Key & k ) { f_tree.insert( k , 0 ); }

 void erase( const Key & k ) { f_tree.remove( k ); }

 [[nodiscard]] std::size_t size( void ) const { return( f_tree.size() ); }

 template< class P >
 [[nodiscard]] bool any_upto( const Key & b , P && p ) const {
  return( upto( Shape::root( f_tree ) , b , p ) );
  }

 [[nodiscard]] const Key * floor( const Key & b ) const {
  const Key * best = nullptr;
  for( auto n = Shape::root( f_tree ) ; n ; ) {
   const Key & k = Shape::read( n ).first;
   if( b < k )
    n = Shape::left( n );
   else {
    best = & k;
    n = Shape::right( n );
    }
   }
  return( best );
  }

 template< class P >
 void erase_from( const Key & b , P && p , bool run ) {
  v_gone.clear();
  from( Shape::root( f_tree ) , b , p , run );
  for( const Key & k : v_gone )
   f_tree.remove( k );
  }

 private:

 // true if p( k ) for a key k <= b of the subtree of n
 template< class N , class P >
 static bool upto( const N * n , const Key & b , P & p ) {
  if( ! n )
   return( false );
  const Key & k = Shape::read( n ).first;
  if( p( k ) )
   return( true );
  if( b < k )
   return( upto( Shape::left( n ) , b , p ) );
  return( upto( Shape::left( n ) , b , p ) ||
          upto( Shape::right( n ) , b , p ) );
  }

 // the keys k >= b of the subtree of n with p( k ) into v_gone, in order;
 // with run, up to the first one without; returns false if it stopped
 template< class N , class P >
 bool from( const N * n , const Key & b , P & p , bool run ) {
  if( ! n )
   return( true );
  const Key & k = Shape::read( n ).first;
  if( k < b )
   return( from( Shape::right( n ) , b , p , run ) );
  if( ! from( Shape::left( n ) , b , p , run ) )
   return( false );
  if( p( k ) )
   v_gone.push_back( k );
  else
   if( run )
    return( false );
  return( from( Shape::right( n ) , b , p , run ) );
  }

 Tree f_tree;                                     ///< the tree
 std::vector< Key > v_gone;                       ///< the keys to remove

 };  // end( class( LexSet ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_MOStructures */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MOStructures.h -------------------------*/
/*--------------------------------------------------------------------------*/
