/*--------------------------------------------------------------------------*/
/*------------------------------- File bnb.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file that includes the tree searches of PBeFS, i.e., the
 * best-first searches with an incumbent (Branch-and-Bound) of a problem
 * satisfying the concept BranchAndBound, save those that need MPI. They
 * are the sequential BranchAndBoundSearch [see BranchAndBound.h],
 * SharedBranchAndBound and DistributedBranchAndBound on an executor of
 * ADS [see ParallelBranchAndBound.h],
 * FarmBranchAndBound, master-worker on a farm with feedback [see
 * FarmBranchAndBound.h],
 * DeterministicBranchAndBound, by synchronous rounds whose outcome does not
 * depend on the timing of the workers [see DeterministicBranchAndBound.h],
 * SelfSplitBranchAndBound, whose workers search apart the subtrees of a
 * first sampling of the tree [see SelfSplitBranchAndBound.h], and
 * StealingBranchAndBound, whose workers keep their children and steal a
 * node when they have none [see StealingBranchAndBound.h], and
 * MultiQueueBranchAndBound, whose workers share the nodes in the MultiQueue
 * of ADS [see MultiQueueBranchAndBound.h], SkeletonBranchAndBound, the
 * depth-first searches of the skeletons [see SkeletonBranchAndBound.h],
 * all of which can also be chosen by a value [see BranchAndBoundScheme.h]
 * and started by a race [see RacingBranchAndBound.h], and
 * IntervalBranchAndBound, by the work stealing of intervals of leaves [see
 * IntervalBranchAndBound.h]; any of them can run in segments with
 * checkpoints [see BBCheckpoint.h], and the master-hub-worker search over
 * MPI is MPIHubBranchAndBound.h. Knapsack, UTS, MaxClique,
 * MaxCommonSubgraph, SubgraphIsomorphism, QAP and Flowshop are problems to
 * try them on.
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

#include "PBeFS/bnb/FarmBranchAndBound.h"

#include "PBeFS/bnb/SelfSplitBranchAndBound.h"

#include "PBeFS/bnb/StealingBranchAndBound.h"

#include "PBeFS/bnb/MultiQueueBranchAndBound.h"

#include "PBeFS/bnb/SkeletonBranchAndBound.h"

#include "PBeFS/bnb/RacingBranchAndBound.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

#include "PBeFS/bnb/IntervalBranchAndBound.h"

#include "PBeFS/bnb/BBCheckpoint.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/bnb/UTS.h"

#include "PBeFS/bnb/MaxClique.h"

#include "PBeFS/bnb/MaxCommonSubgraph.h"

#include "PBeFS/bnb/SubgraphIsomorphism.h"

#include "PBeFS/bnb/QAP.h"

#include "PBeFS/bnb/Flowshop.h"

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_bnb */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File bnb.h -----------------------------*/
/*--------------------------------------------------------------------------*/
