/*--------------------------------------------------------------------------*/
/*------------------------------ File astar.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file that includes the graph searches of PBeFS save those that
 * need MPI, HDAStarMPI.h and ParallelIDAStarMPI.h.
 *
 * The searches look for a path of least cost in a Domain [see Domain.h],
 * i.e., an implicit graph given by an initial state, a goal test and an
 * expansion procedure, with:
 *
 * - AStar, the sequential A* [see AStar.h];
 *
 * - IDAStar, the iterative-deepening A*, which holds only the current path
 *   [see IDAStar.h];
 *
 * - ParallelIDAStar, the iterations of IDA* below a frontier of units
 *   built by A*, by uniform-cost search or by depth, on the workers of an
 *   executor of ADS, synchronous or asynchronous, with work stealing and
 *   the splitting of the stacks, and ParallelIDAStarMPI on the processes
 *   of MPI [see ParallelIDAStar.h and ParallelIDAStarMPI.h];
 *
 * - DDDAStar, the A* with delayed duplicate detection, whose OPEN and
 *   CLOSED lists are files on the disk [see DDDAStar.h];
 *
 * - SPAStar, Simple Parallel A*, with one OPEN list shared by the workers
 *   [see SPAStar.h];
 *
 * - MQAStar, MultiQueue A*, as SPAStar with a relaxed OPEN list, the
 *   MultiQueue of ADS [see MQAStar.h];
 *
 * - HDAStar, Hash Distributed A*, in which a distribution function
 *   [see Distribution.h] assigns every state to the worker that owns it,
 *   on an executor of ADS, and HDAStarMPI on the processes
 *   of MPI [see HDAStarMPI.h];
 *
 * - SafePBNF, Parallel Best-NBlock-First search [see SafePBNF.h];
 *
 * - PASE, Parallel A* for Slow Expansions [see PASE.h];
 *
 * - FarmAStar, a master-worker A* on a farm with feedback [see
 *   FarmAStar.h].
 *
 * The domains are SlidingTile (the puzzle of \f$ W H - 1 \f$ tiles),
 * SequenceScore (the same with the heuristic of Nilsson), Grid (pathfinding
 * on a GridMap with 4 or 8 moves, or on a GridBoard with unit or life
 * costs), SequenceAlignment (the optimal alignment of several sequences of
 * amino acids), ExplicitGraph (a weighted graph given by adjacency lists)
 * and Romania (the road map of Romania as one of them), and
 * the problems of the textbook of Russell and Norvig that the examples
 * solve, i.e., NQueens, VacuumWorld, MissionariesCannibals, WaterJugs and
 * TravelingSalesman; SlowDomain slows down the expansions of any of them.
 * SlidingTilePDB is SlidingTile with the additive pattern databases of
 * PatternDatabase.h in place of the Manhattan distance. The permutation
 * puzzles PancakePuzzle and BurnedPancakePuzzle (with the gap heuristic),
 * TopSpin and RubiksCube, and the Towers of Hanoi with four pegs
 * (TowersOfHanoi), take the pattern databases of any abstraction of
 * AbstractionPDB.h and PermutationPDB.h [see Permutation.h for the
 * rankings], and DimacsGraph.h reads the road graphs of DIMACS into an
 * ExplicitGraph.
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

#ifndef __PBeFS_astar
 #define __PBeFS_astar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/core.h"

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/AnytimeAStar.h"

#include "PBeFS/astar/AbstractionPDB.h"

#include "PBeFS/astar/DDDAStar.h"

#include "PBeFS/astar/DimacsGraph.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/ExplicitGraph.h"

#include "PBeFS/astar/FarmAStar.h"

#include "PBeFS/astar/FocalSearch.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/ImprovedOptimisticSearch.h"

#include "PBeFS/astar/LPAStar.h"

#include "PBeFS/astar/MissionariesCannibals.h"

#include "PBeFS/astar/NQueens.h"

#include "PBeFS/astar/OptimisticSearch.h"

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/ParallelIDAStar.h"

#include "PBeFS/astar/MQAStar.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/Permutation.h"

#include "PBeFS/astar/PermutationPDB.h"

#include "PBeFS/astar/PotentialSearch.h"

#include "PBeFS/astar/Romania.h"

#include "PBeFS/astar/RubiksCube.h"

#include "PBeFS/astar/SPAStar.h"

#include "PBeFS/astar/SafePBNF.h"

#include "PBeFS/astar/SequenceAlignment.h"

#include "PBeFS/astar/SequenceScore.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/SlowDomain.h"

#include "PBeFS/astar/TopSpin.h"

#include "PBeFS/astar/TowersOfHanoi.h"

#include "PBeFS/astar/TravelingSalesman.h"

#include "PBeFS/astar/VacuumWorld.h"

#include "PBeFS/astar/WaterJugs.h"

#include "PBeFS/astar/Weight.h"

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_astar */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File astar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
