/*--------------------------------------------------------------------------*/
/*------------------------------- File bnb.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file that includes the tree searches of PBeFS, i.e., the
 * best-first searches with an incumbent (Branch-and-Bound) of a problem
 * satisfying the concept BranchAndBound, save the parts that need an
 * external library: FarmBranchAndBound.h (FastFlow) and MPIBranchAndBound.h
 * (MPI). They are the sequential BranchAndBoundSearch [see
 * BranchAndBound.h], SharedBranchAndBound and DistributedBranchAndBound
 * on a runtime of the core [see ParallelBranchAndBound.h and core.h],
 * DeterministicBranchAndBound, by synchronous rounds whose outcome does not
 * depend on the timing of the workers [see DeterministicBranchAndBound.h],
 * SelfSplitBranchAndBound, whose workers search apart the subtrees of a
 * first sampling of the tree [see SelfSplitBranchAndBound.h], and
 * StealingBranchAndBound, whose workers keep their children and steal a
 * node when they have none [see StealingBranchAndBound.h], all of which can
 * also be chosen by a value [see BranchAndBoundScheme.h]. Knapsack is a 0-1
 * knapsack problem to try them on [see Knapsack.h].
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

#ifndef __PBeFS_bnb
 #define __PBeFS_bnb
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/core.h"

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "PBeFS/bnb/DeterministicBranchAndBound.h"

#include "PBeFS/bnb/SelfSplitBranchAndBound.h"

#include "PBeFS/bnb/StealingBranchAndBound.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

#include "PBeFS/bnb/Knapsack.h"

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_bnb */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File bnb.h -----------------------------*/
/*--------------------------------------------------------------------------*/
