/*--------------------------------------------------------------------------*/
/*---------------------------- File NodeTable.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodeTable, the closed list of the
 * searches: a map from the generated states to their canonical Node, taken
 * from the set of data structures chosen [see Structures.h]. A
 * NodeTable< NodeT > has a constructor taking the initial number of slots,
 * the static slots_for( n ), the number of slots that holds n nodes without
 * growing, find_or_insert( state , hash , make ), which looks up state,
 * whose hash is hash, and inserts the node returned by make() if it is
 * absent, returning the node of state and true if it has just been
 * inserted, find( state , hash ), which returns nullptr if state has never
 * been generated, prefetch( hash ), which brings in the cache the first
 * slot probed for hash, so that a search can overlap the cache misses of
 * the lookups of all the successors of a node, size() and bytes().
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

#ifndef __PBeFS_NodeTable
 #define __PBeFS_NodeTable
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/

// the structures of the chosen set [see Structures.h]
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/NodeTable.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/NodeTable.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/NodeTable.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_NodeTable */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodeTable.h --------------------------*/
/*--------------------------------------------------------------------------*/
