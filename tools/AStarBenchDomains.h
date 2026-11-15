/*--------------------------------------------------------------------------*/
/*------------------------ File AStarBenchDomains.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The Domains of AStarBench, one PBEFS_BENCH_DOMAIN( D ) per line, which
 * the includer defines, with PermutationPDB< P > the permutation puzzle P
 * with its pattern databases: AStarBenchPart.cpp numbers them from 1 in
 * this order, and the build compiles it for each of them and each of its
 * three variants, counting the lines that start with PBEFS_BENCH_DOMAIN
 * [see AStarBenchPart.cpp].
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- DOMAINS ---------------------------------*/
/*--------------------------------------------------------------------------*/

// the sliding-tile puzzles
PBEFS_BENCH_DOMAIN( EightPuzzle )
PBEFS_BENCH_DOMAIN( FifteenPuzzle )
PBEFS_BENCH_DOMAIN( TwentyFourPuzzle )

// the grids, the boards and the 2D navigation
PBEFS_BENCH_DOMAIN( Grid< 4 > )
PBEFS_BENCH_DOMAIN( Grid< 8 > )
PBEFS_BENCH_DOMAIN( RobotNav2D )

// the alignment of sequences, the graphs and the Rubik's cube
PBEFS_BENCH_DOMAIN( SequenceAlignment<> )
PBEFS_BENCH_DOMAIN( ExplicitGraph< std::int64_t > )
PBEFS_BENCH_DOMAIN( RubiksCube )
PBEFS_BENCH_DOMAIN( PDBDomain< RubiksCube , RubikPDBHeuristic > )

// the permutation puzzles, alone and with their pattern databases
PBEFS_BENCH_DOMAIN( PancakePuzzle< 12 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< PancakePuzzle< 12 > > )
PBEFS_BENCH_DOMAIN( PancakePuzzle< 16 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< PancakePuzzle< 16 > > )
PBEFS_BENCH_DOMAIN( PancakePuzzle< 20 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< PancakePuzzle< 20 > > )
PBEFS_BENCH_DOMAIN( PancakePuzzle< 32 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< PancakePuzzle< 32 > > )
PBEFS_BENCH_DOMAIN( PancakePuzzle< 64 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< PancakePuzzle< 64 > > )
PBEFS_BENCH_DOMAIN( BurnedPancakePuzzle< 8 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< BurnedPancakePuzzle< 8 > > )
PBEFS_BENCH_DOMAIN( BurnedPancakePuzzle< 10 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< BurnedPancakePuzzle< 10 > > )
PBEFS_BENCH_DOMAIN( BurnedPancakePuzzle< 12 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< BurnedPancakePuzzle< 12 > > )
PBEFS_BENCH_DOMAIN( BurnedPancakePuzzle< 16 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< BurnedPancakePuzzle< 16 > > )
PBEFS_BENCH_DOMAIN( TopSpin< 12 , 4 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< TopSpin< 12 , 4 > > )
PBEFS_BENCH_DOMAIN( TopSpin< 14 , 4 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< TopSpin< 14 , 4 > > )
PBEFS_BENCH_DOMAIN( TopSpin< 16 , 4 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< TopSpin< 16 , 4 > > )
PBEFS_BENCH_DOMAIN( TopSpin< 17 , 4 > )
PBEFS_BENCH_DOMAIN( PermutationPDB< TopSpin< 17 , 4 > > )

// the Towers of Hanoi with their pattern databases
PBEFS_BENCH_DOMAIN( PDBDomain< TowersOfHanoi< 12 > , HanoiPDB > )
PBEFS_BENCH_DOMAIN( PDBDomain< TowersOfHanoi< 16 > , HanoiPDB > )
PBEFS_BENCH_DOMAIN( PDBDomain< TowersOfHanoi< 20 > , HanoiPDB > )

// the problems of the examples
PBEFS_BENCH_DOMAIN( NQueens )
PBEFS_BENCH_DOMAIN( VacuumWorld )
PBEFS_BENCH_DOMAIN( MissionariesCannibals )
PBEFS_BENCH_DOMAIN( WaterJugs )
PBEFS_BENCH_DOMAIN( TravelingSalesman )
PBEFS_BENCH_DOMAIN( Romania )

/*--------------------------------------------------------------------------*/
/*---------------------- End File AStarBenchDomains.h ----------------------*/
/*--------------------------------------------------------------------------*/
