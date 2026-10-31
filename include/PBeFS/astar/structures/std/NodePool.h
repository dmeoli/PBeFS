/*--------------------------------------------------------------------------*/
/*---------------------------- File NodePool.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodePool of the set of data structures
 * of the standard library [see Structures.h], the arena the nodes of A* are
 * allocated from, on the std::deque; Node.h includes it in place of its own.
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

#ifndef __PBeFS_std_NodePool
 #define __PBeFS_std_NodePool
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <deque>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS NodePool -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the arena of the nodes, on the std::deque
/** The interface is the one described in Node.h. */
template< class NodeT >
class NodePool {

 public:

 /// constructor; the chunk is ignored, std::deque has its own blocks

 explicit NodePool( std::size_t = 0 ) {}

 /// returns a new (value-initialized) node, whose address never changes

 [[nodiscard]] NodeT * allocate( void ) {
  return( & v_nodes.emplace_back() );
  }

 /// the number of nodes allocated so far

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_nodes.size() );
  }

 /// the memory taken by the nodes, in bytes (the blocks of the deque apart)

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( v_nodes.size() * sizeof( NodeT ) );
  }

 private:

 std::deque< NodeT > v_nodes;                              ///< the nodes

 };  // end( class( NodePool ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_std_NodePool */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodePool.h ---------------------------*/
/*--------------------------------------------------------------------------*/
