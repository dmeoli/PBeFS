/*--------------------------------------------------------------------------*/
/*---------------------------- File NodePool.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodePool of the set of data structures
 * of ADS [see Structures.h], the arena the nodes of A* are allocated from,
 * on the Arena of ADS, on huge pages; Node.h includes it in place of its
 * own.
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

#ifndef __PBeFS_ads_NodePool
 #define __PBeFS_ads_NodePool
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include "ADS/memory/Arena.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS NodePool -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the arena of the nodes, on the Arena of ADS, on huge pages
/** The interface is the one described in Node.h. */
template< class NodeT >
class NodePool : public ADS::Arena< NodeT , ADS::HugePageAllocator< NodeT > > {

 public:

 /// constructor, taking the number of nodes of a chunk

 explicit NodePool( std::size_t chunk = std::size_t( 1 ) << 16 )
  : ADS::Arena< NodeT , ADS::HugePageAllocator< NodeT > >( chunk ) {}

 /// returns a new (value-initialized) node, whose address never changes

 [[nodiscard]] NodeT * allocate( void ) { return( this->make() ); }

 };  // end( class( NodePool ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_NodePool */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodePool.h ---------------------------*/
/*--------------------------------------------------------------------------*/
