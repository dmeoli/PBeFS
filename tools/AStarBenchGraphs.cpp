/*--------------------------------------------------------------------------*/
/*----------------------- File AStarBenchGraphs.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The searches of AStarBench on the alignments, the road graphs and Rubik's
 * cube [see AStarBenchSolve.h].
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AStarBenchSolve.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

PBEFS_BENCH_SOLVE( SequenceAlignment<> );
PBEFS_BENCH_SOLVE( ExplicitGraph< std::int64_t > );
PBEFS_BENCH_SOLVE( RubiksCube );
PBEFS_BENCH_SOLVE( PDBDomain< RubiksCube , RubikPDBHeuristic > );

/*--------------------------------------------------------------------------*/
/*--------------------- End File AStarBenchGraphs.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
