/*--------------------------------------------------------------------------*/
/*---------------------------- File NodePool.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodePool of the set of data structures
 * of Boost [see Structures.h], the arena the nodes of A* are allocated from,
 * on the object_pool of Boost.Pool; Node.h includes it in place of its own.
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

#ifndef __PBeFS_boost_NodePool
 #define __PBeFS_boost_NodePool
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <boost/pool/object_pool.hpp>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS NodePool -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the arena of the nodes, on the object_pool of Boost.Pool
/** The interface is the one described in Node.h. */
template< class NodeT >
class NodePool {

 public:

 /// constructor, taking the number of nodes of the first chunk

 explicit NodePool( std::size_t chunk = std::size_t( 1 ) << 16 )
  : f_pool( chunk ) {}

 /// returns a new (value-initialized) node, whose address never changes

 [[nodiscard]] NodeT * allocate( void ) {
  ++f_size;
  return( f_pool.construct() );
  }

 /// the number of nodes allocated so far

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_size );
  }

 /// the memory taken by the nodes, in bytes (the chunks of the pool apart)

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_size * sizeof( NodeT ) );
  }

 private:

 boost::object_pool< NodeT > f_pool;                       ///< the nodes
 std::size_t f_size = 0;                             ///< nodes allocated

 };  // end( class( NodePool ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_boost_NodePool */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodePool.h ---------------------------*/
/*--------------------------------------------------------------------------*/
