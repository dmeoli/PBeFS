/*--------------------------------------------------------------------------*/
/*--------------------------- File Structures.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the choice of the set of data structures of A*: the arena
 * of the nodes (NodePool), the closed list (NodeTable), the OPEN lists
 * (HeapOpenList and BucketOpenList) and the adjacency of the explicit
 * graphs and of the grids. The set is that of ADS unless one of the macros
 * PBEFS_STRUCTURES_STD and PBEFS_STRUCTURES_BOOST is defined, which CMake
 * does by the option PBEFS_STRUCTURES (ads, std or boost) and the makefile
 * by PBeFSSTRUCT; then Node.h, NodeTable.h, OpenList.h, ExplicitGraph.h and
 * Grid.h take the structures of the standard library or of Boost from
 * structures/std/ or structures/boost/ in place of those of ADS in
 * structures/ads/, with the same interface. The searches are the same with
 * every set, and so is the order of their expansions but for the ties of f
 * and g in OPEN, which the sets break differently; the standard library and
 * Boost are there to measure the searches against them.
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

#ifndef __PBeFS_Structures
 #define __PBeFS_Structures
                      /* self-identification: #endif at the end of the file */

#if ( defined( PBEFS_STRUCTURES_STD ) + defined( PBEFS_STRUCTURES_BOOST ) + \
      defined( PBEFS_STRUCTURES_ADS ) ) > 1
 #error "PBeFS: more than one set of data structures is chosen"
#endif

#if ! ( defined( PBEFS_STRUCTURES_STD ) || defined( PBEFS_STRUCTURES_BOOST ) )
 #ifndef PBEFS_STRUCTURES_ADS
  /// defined when the data structures are those of ADS, the default
  #define PBEFS_STRUCTURES_ADS
 #endif
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Structures */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Structures.h --------------------------*/
/*--------------------------------------------------------------------------*/
