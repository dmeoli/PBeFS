/*--------------------------------------------------------------------------*/
/*---------------------------- File NodeTable.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodeTable of the set of data structures
 * of ADS [see Structures.h], the closed list of A*, on the PointerTable of
 * ADS, on huge pages; NodeTable.h includes it in place of its own.
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

#ifndef __PBeFS_ads_NodeTable
 #define __PBeFS_ads_NodeTable
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include "ADS/dictionary/PointerTable.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS NodeTable -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the closed list, on the PointerTable of ADS, on huge pages
/** The interface is the one described in NodeTable.h. */
template< class NodeT >
class NodeTable
 : public ADS::PointerTable< NodeT , ADS::StateEqual ,
                            ADS::HugePageAllocator< std::byte > > {

 public:

 /// the number of slots that holds n nodes without growing

 [[nodiscard]] static std::size_t slots_for( std::size_t n ) noexcept {
  return( n / 7 * 10 + 16 );
  }

 /// constructor, taking the initial number of slots

 explicit NodeTable( std::size_t capacity = 1024 )
  : ADS::PointerTable< NodeT , ADS::StateEqual ,
                       ADS::HugePageAllocator< std::byte > >( capacity ) {}

 };  // end( class( NodeTable ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_NodeTable */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodeTable.h --------------------------*/
/*--------------------------------------------------------------------------*/
