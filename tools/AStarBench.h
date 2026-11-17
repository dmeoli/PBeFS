/*--------------------------------------------------------------------------*/
/*--------------------------- File AStarBench.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The parts of AStarBench that its translation units share: the options of
 * the command line, the output of a run and the runs of a Domain. The
 * searches are in AStarBenchSolve.h, compiled only by AStarBenchPart.cpp,
 * which instantiates them for one of the Domains of AStarBenchDomains.h
 * [see AStarBench.cpp for the options].
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

#ifndef __PBeFS_AStarBench
 #define __PBeFS_AStarBench
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstdint>

#include <cstdio>

#include <cstdlib>

#include <fstream>

#include <iostream>

#include <limits>

#include <map>

#include <memory>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/DDDAStar.h"

#include "PBeFS/astar/DimacsGraph.h"

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/PermutationPDB.h"

#include "PBeFS/astar/RubiksCube.h"

#include "PBeFS/astar/TopSpin.h"

#include "PBeFS/astar/TowersOfHanoi.h"

#include "PBeFS/astar/MissionariesCannibals.h"

#include "PBeFS/astar/NQueens.h"

#include "PBeFS/astar/Romania.h"

#include "PBeFS/astar/SequenceAlignment.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/RobotNav2D.h"

#include "PBeFS/astar/SlowDomain.h"

#include "PBeFS/astar/TravelingSalesman.h"

#include "PBeFS/astar/VacuumWorld.h"

#include "PBeFS/astar/WaterJugs.h"

#ifdef PBEFS_BENCH_MPI
 #include "PBeFS/astar/HDAStarMPI.h"
 #include "PBeFS/astar/ParallelIDAStarMPI.h"
#else
 #include "PBeFS/astar/AnytimeAStar.h"
 #include "PBeFS/astar/FocalSearch.h"
 #include "PBeFS/astar/ParallelIDAStar.h"
 #include "PBeFS/astar/HDAStar.h"
 #include "PBeFS/astar/ImprovedOptimisticSearch.h"
 #include "PBeFS/astar/LPAStar.h"
 #include "PBeFS/astar/OptimisticSearch.h"
 #include "PBeFS/astar/PotentialSearch.h"
 #include "PBeFS/astar/SPAStar.h"
 #include "PBeFS/astar/MQAStar.h"
 #include "PBeFS/astar/PASE.h"
 #include "PBeFS/astar/EPASE.h"
 #include "PBeFS/astar/MPLP.h"
 #include "PBeFS/astar/PWAStar.h"
 #include "PBeFS/astar/SafePBNF.h"
 #ifdef PBEFS_HAS_OPENMP
  #include "ADS/parallel/OpenMPExecutor.h"
 #endif
 #include "PBeFS/astar/FarmAStar.h"
 #include "PBeFS/astar/BiBestFirst.h"
 #include "PBeFS/astar/NBS.h"
 #include "PBeFS/astar/DVCBS.h"
 #include "PBeFS/astar/GBFHS.h"
 #include "PBeFS/astar/ParallelBidirectional.h"
 #include "PBeFS/astar/PEMSearch.h"
 #include "PBeFS/astar/IDDAStar.h"
 #include "PBeFS/astar/BeamSearch.h"
 #include "PBeFS/astar/EnforcedHillClimbing.h"
 #include "PBeFS/astar/FringeSearch.h"
 #include "PBeFS/astar/GreedySearch.h"
 #include "PBeFS/astar/HPAStar.h"
 #include "PBeFS/astar/IteratedSearch.h"
 #include "PBeFS/astar/ParallelGreedy.h"
 #include "PBeFS/astar/PartialExpansion.h"
 #include "PBeFS/astar/RBFS.h"
 #include "PBeFS/astar/RealTimeSearch.h"
 #include "PBeFS/astar/WidthSearch.h"
 #ifdef PBEFS_HAS_SYLVAN
  #include "PBeFS/symbolic.h"
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  #include "ADS/parallel/FastFlowExecutor.h"
 #endif
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------- TYPES AND GLOBAL VARIABLES -----------------------*/
/*--------------------------------------------------------------------------*/

// the options of the command line, by name
inline std::map< std::string , std::string > opt = {
 { "domain" , "tiles15" } , { "file" , "" } , { "ids" , "" } ,
 { "random" , "" } , { "points" , "" } , { "algo" , "astar" } ,
 { "runtime" , "thread" } , { "dist" , "Z" } , { "threads" , "1" } ,
 { "batch" , "64" } , { "flush" , "64" } , { "open" , "buckets" } ,
 { "block" , "" } , { "limit" , "0" } , { "repeat" , "1" } ,
 { "minexp" , "32" } , { "nblock" , "" } , { "eps" , "1" } ,
 { "slow" , "0" } , { "expect" , "0" } , { "size" , "" } ,
 { "parts" , "64" } , { "dir" , "" } , { "origin" , "bottom" } ,
 { "matrix" , "" } , { "metric" , "euclidean" } , { "gap" , "0" } ,
 { "plane" , "auto" } , { "weight" , "" } , { "dd" , "0" } ,
 { "anytime" , "0" } , { "schedule" , "" } , { "optimism" , "0" } ,
 { "focal" , "d" } , { "paths" , "0" } , { "frontier" , "0" } ,
 { "rule" , "" } , { "units" , "steal" } , { "async" , "0" } ,
 { "bounds" , "1" } , { "split" , "0" } , { "table" , "0" } ,
 { "det" , "0" } , { "pdb" , "" } , { "cache" , "" } , { "reverse" , "0" } ,
 { "edge" , "-1" } , { "fraction" , "0.5" } , { "side" , "auto" } ,
 { "allbounds" , "0" } , { "eager" , "0" } , { "disk" , "0" } ,
 { "slots" , "1048576" } ,
 { "eval" , "eager" } ,
 { "queues" , "best:0" } , { "tie" , "fifo" } ,
 { "boost" , "0" } , { "epsilon" , "0.2" } , { "reopen" , "0" } ,
 { "seed" , "1" } , { "pref" , "ignore" } , { "width" , "2" } ,
 { "sge" , "0" } , { "rounds" , "0" } ,
 { "cutoff" , "0" } , { "beam" , "100" } , { "lookahead" , "1" } ,
 { "trials" , "1" } , { "cluster" , "10" } , { "scan" , "64" } ,
 { "gen" , "0" } , { "dispatch" , "1" } , { "step" , "0.5" } ,
 { "time" , "0" } , { "edges" , "" } , { "scale" , "5" } ,
 { "factors" , "0" } , { "diag" , "1" } ,
 { "symdir" , "f" } , { "images" , "seq" } , { "lace" , "1" } ,
 { "bddmem" , "2048" } , { "perimeter" , "0" } , { "list" , "0" } ,
 { "symstats" , "0" } };

inline int rank = 0;                   // the MPI rank (0 without MPI)

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the indices in a list "a-b,c,..." (1-based), or all of [ 1 , n ]

inline std::vector< int > parse_ids( const std::string & s , int n )
{
 std::vector< int > v;
 if( s.empty() ) {
  for( int i = 1 ; i <= n ; ++i )
   v.push_back( i );
  return( v );
  }
 std::istringstream in( s );
 std::string part;
 while( std::getline( in , part , ',' ) ) {
  const auto dash = part.find( '-' );
  const int a = std::stoi( part.substr( 0 , dash ) );
  const int b = dash == std::string::npos
                ? a : std::stoi( part.substr( dash + 1 ) );
  for( int i = a ; i <= b ; ++i )
   v.push_back( i );
  }
 return( v );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the integer option o, or d if it is not given

inline int option( const std::string & o , int d )
{
 return( opt[ o ].empty() ? d : std::stoi( opt[ o ] ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the contents of the file of --file

inline std::string read_file( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
 std::ostringstream os;
 os << in.rdbuf();
 return( os.str() );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the weight of --weight, or that of the search when it is not given:
// infinite (GBFS) for the searches of GreedyOptions, 1 for the others

inline std::string weight_option( void )
{
 if( ! opt[ "weight" ].empty() )
  return( opt[ "weight" ] );
 const std::string & a = opt[ "algo" ];
 return( ( a == "gbfs" ) || ( a == "ehc" ) || ( a == "pgbfs" ) ||
         ( a == "hdgbfs" ) ? "inf" : "1" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// prints the line of a run

template< class S >
inline void print( const std::string & instance , const S & s )
{
 if( rank != 0 )
  return;
 const auto & st = s.stats;
 const char * status[] = { "solved" , "unsolvable" , "node_limit" ,
                           "disk_limit" };
 std::cout.precision( 12 );
 std::cout << opt[ "domain" ] << ',' << instance << ',' << opt[ "algo" ]
           << ',' << opt[ "runtime" ] << ',' << opt[ "dist" ] << ','
           << opt[ "threads" ] << ',' << opt[ "batch" ] << ','
           << opt[ "flush" ] << ',' << opt[ "open" ] << ','
           << status[ int( s.status ) ] << ',' << s.cost << ','
           << st.expanded << ',' << st.reexpanded << ',' << st.generated
           << ',' << st.duplicates << ',' << st.sent << ',' << st.messages
           << ',' << st.stored << ',' << st.load_balance() << ','
           << st.communication_overhead() << ',' << st.seconds << ','
           << weight_option() << ',' << s.bound << ','
           << s.improvements.size() << ',' << st.evaluated << ','
           << st.iterations << std::endl;
 if( opt[ "paths" ] == "1" )
  for( const auto & x : s.improvements )
   std::cerr << "path," << instance << ',' << x.cost << ',' << x.bound
             << ',' << x.expanded << ',' << x.seconds << std::endl;
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d, slowed down by --slow iterations per expansion if not 0, or
// with the edges of --edges: defined in AStarBenchSolve.h, and instantiated
// for each Domain of AStarBenchDomains.h by AStarBenchPart.cpp, so that no
// translation unit compiles every search for every Domain

template< class D >
Solution< state_t< D > , cost_t< D > > solve_maybe_slow( const D & d );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Domain d --repeat times, the line labelled by instance

template< class D >
void run_domain( const std::string & instance , const D & d )
{
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int r = 0 ; r < repeat ; ++r )
  print( instance , solve_maybe_slow( d ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the fields of --size separated by ':', at least one

inline std::vector< std::string > size_fields( void )
{
 std::vector< std::string > v;
 std::istringstream in( opt[ "size" ] );
 std::string part;
 while( std::getline( in , part , ':' ) )
  v.push_back( part );
 if( v.empty() || v[ 0 ].empty() )
  throw( std::invalid_argument( "--size is needed for " + opt[ "domain" ] ) );
 return( v );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the objects of the tables of --pdb, k at a time out of n: the first k,
// the next k, ... (and the last ones, fewer)

inline std::vector< std::vector< std::uint8_t > > pdb_patterns( unsigned n )
{
 std::vector< std::vector< std::uint8_t > > v;
 if( opt[ "pdb" ].empty() )
  return( v );
 const unsigned k = unsigned( std::stoul( opt[ "pdb" ] ) );
 if( ( k == 0 ) || ( k > n ) )
  throw( std::invalid_argument( "--pdb wants 1 to the size of the puzzle" ) );
 for( unsigned a = 0 ; a < n ; a += k ) {
  v.emplace_back();
  for( unsigned x = a ; x < std::min( n , a + k ) ; ++x )
   v.back().push_back( std::uint8_t( x ) );
  }
 return( v );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the permutation puzzle d, with the maximum of the tables of --pdb
// if given

template< class D >
inline void run_permutation( const std::string & instance , const D & d ,
                             unsigned n )
{
 const auto pats = pdb_patterns( n );
 if( pats.empty() ) {
  run_domain( instance , d );
  return;
  }
 using PDB = PermutationPDB< LexRanking , 8 >;
 std::vector< PDB > v;
 for( const auto & p : pats )
  v.emplace_back( d , 0 , p , unsigned( std::stoul( opt[ "threads" ] ) ) ,
                  opt[ "cache" ] );
 const PermutationHeuristic< D > h( d , std::move( v ) );
 run_domain( instance , PDBDomain< D , PermutationHeuristic< D > >( d , h ) );
 }

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AStarBench */

/*--------------------------------------------------------------------------*/
/*------------------------- End File AStarBench.h --------------------------*/
/*--------------------------------------------------------------------------*/
