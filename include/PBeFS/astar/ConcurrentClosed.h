/*--------------------------------------------------------------------------*/
/*------------------------ File ConcurrentClosed.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ConcurrentClosed, the closed list
 * that the workers of SPAStar and MQAStar share, taken from the set of
 * data structures chosen [see Structures.h]. A ConcurrentClosed< D , Lock >
 * maps the states of the Domain D to their canonical Node and protects the
 * fields of the nodes that change during the search (g, parent, op,
 * closed and expanded) with locks of type Lock. It has a constructor
 * taking the Domain, the type Local of the nodes that a worker allocates,
 * whose stored() is their number in the table, lock( n ), the lock of the
 * node n, under which a worker reads and writes those fields, the type
 * Successor of a state generated with its g, h, hash and operator [see
 * OperatorDomain], prefetch( hash ), which brings in the cache what an
 * update() of a state of that hash reads first, so that a worker overlaps
 * the misses of all the successors of a node, and update( t , hash , g ,
 * h , parent , op , local , duplicates ), which records the path to the
 * state t, of hash hash, of cost g through parent and the operator op: it
 * returns the node of t if t is new or g is smaller than its g, which then
 * is set to g with the parent and the operator and the node is no longer
 * closed, and nullptr otherwise, adding 1 to duplicates. The fields of a
 * node returned may be read only under its lock, but for h, which never
 * changes.
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

#ifndef __PBeFS_ConcurrentClosed
 #define __PBeFS_ConcurrentClosed
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/

// the structures of the chosen set [see Structures.h]
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/ConcurrentClosed.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/ConcurrentClosed.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/ConcurrentClosed.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ConcurrentClosed */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ConcurrentClosed.h -----------------------*/
/*--------------------------------------------------------------------------*/
