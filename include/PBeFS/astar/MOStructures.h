/*--------------------------------------------------------------------------*/
/*--------------------------- File MOStructures.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the data structures of the multi-objective searches,
 * those of the chosen set [see Structures.h]:
 *
 * - MOHeap< Item , Better >, a priority queue of plain items, the best by
 *   Better on top, with push( x ), top(), pop() (which returns the item
 *   removed), empty() and size();
 *
 * - LexSet< Key >, an ordered set of vectors (std::array) by the
 *   lexicographic order, with insert( k ), erase( k ), size(),
 *   any_upto( b , p ), true if p( k ) holds for a key \f$ k \leq_{lex} b
 *   \f$, which stops at the first one, floor( b ), the largest key
 *   \f$ \leq_{lex} b \f$ (nullptr if none), and erase_from( b , p , run ),
 *   which removes the keys \f$ k \geq_{lex} b \f$ with p( k ), all of them,
 *   or, with run, those up to the first key without.
 *
 * The set of ADS has the DaryHeap and the AVLTree of ADS, whose keys up to
 * b are visited from the root down as in the check of Ren et al. (SoCS
 * 2022); that of the standard library has std::priority_queue and
 * std::set, and that of Boost boost::heap::d_ary_heap and
 * boost::container::set on an AVL tree, whose keys are visited in order.
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

#ifndef __PBeFS_MOStructures
 #define __PBeFS_MOStructures
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/

// the structures of the chosen set [see Structures.h]
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/MOStructures.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/MOStructures.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/MOStructures.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MOStructures */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MOStructures.h -------------------------*/
/*--------------------------------------------------------------------------*/
