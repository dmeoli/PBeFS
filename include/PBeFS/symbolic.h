/*--------------------------------------------------------------------------*/
/*---------------------------- File symbolic.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file that includes the symbolic searches of PBeFS, which search
 * sets of states represented by binary decision diagrams of Sylvan and
 * need the library PBeFS::Symbolic, built only if Sylvan is found:
 *
 * - Bdd and BddRuntime, the diagrams and the scope of Sylvan [see Bdd.h];
 *
 * - the concept SymbolicDomain and its refinements [see SymbolicDomain.h];
 *
 * - SymbolicModel, a SymbolicDomain compiled into diagrams [see
 *   SymbolicModel.h];
 *
 * - SymbolicSearch, the breadth-first search, Dijkstra's algorithm and
 *   BDDA*, forward, backward and bidirectional [see SymbolicSearch.h];
 *
 * - the domains SymbolicSlidingTile, SymbolicPancake [see
 *   SymbolicPuzzles.h] and SymbolicSAS, a SAS+ task [see SymbolicSAS.h].
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

#ifndef __PBeFS_symbolic
 #define __PBeFS_symbolic
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/symbolic/Bdd.h"

#include "PBeFS/symbolic/SymbolicDomain.h"

#include "PBeFS/symbolic/SymbolicModel.h"

#include "PBeFS/symbolic/SymbolicPuzzles.h"

#include "PBeFS/symbolic/SymbolicSAS.h"

#include "PBeFS/symbolic/SymbolicSearch.h"

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_symbolic */

/*--------------------------------------------------------------------------*/
/*------------------------- End File symbolic.h ----------------------------*/
/*--------------------------------------------------------------------------*/
