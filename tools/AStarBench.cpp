/*--------------------------------------------------------------------------*/
/*-------------------------- File AStarBench.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the searches of PBeFS: solves instances of the
 * 15-puzzle, of the 8-puzzle, of the 24-puzzle, of grid pathfinding, of
 * multiple sequence alignment or of the problems of the examples with a
 * chosen search, runtime, work distribution and number of workers, and
 * prints one line of comma-separated values per run.
 *
 * Compiled with PBEFS_BENCH_MPI defined, it is the benchmark of HDAStarMPI,
 * to be started under mpiexec; only the process of rank 0 prints.
 *
 *     AStarBench [options]
 *
 *     --domain  tiles15 | tiles8 | tiles24 | grid4 | grid8 | board | msa |
 *               queens | vacuum | missionaries | jugs | tsp | romania |
 *               pancake | burned | topspin | rubik | hanoi | road
 *                                                          (tiles15)
 *     --file    <instances or map>     korf100.txt, instances50.txt, a
 *                                      Moving AI map, a board, the
 *                                      sequences, the cities, korf10.txt,
 *                                      toh20.txt, or a road graph (the
 *                                      name of its .gr without .gr)
 *                                      (below)
 *     --ids     <a-b,c,...>            instances of the file (all)
 *     --random  <W>x<H>:<p>:<seed>[:<c>]  a random grid instead of a map,
 *                                      its free cells costing 1 to c (1),
 *               or <seed>[:<n>[:<total>]]  random sequences (below)
 *     --size    <size>                 the size of an example (below)
 *     --points  <sx>,<sy>,<gx>,<gy>    start and goal (bottom corners)
 *     --origin  bottom | top           the row 0 of a board   (bottom)
 *     --matrix  <scores>               the substitution matrix of msa
 *                                                          (PAM250)
 *     --gap     <penalty>              a linear gap penalty of msa (none),
 *               or the x of GAP-x of the pancakes          (0)
 *     --plane   <d> | 1/<k> | auto       the thickness of the hyperplanes
 *                                      of --dist HWD              (auto)
 *     --metric  euclidean | manhattan  the distances of tsp   (euclidean)
 *     --algo    astar | ida | pida | ddd | spa | mqa | hda | pbnf | pase |
 *               farm | awa | ara | optimistic | ios | focal | dps | lpa
 *                                                          (astar)
 *     --weight  <weight w of the heuristic, or bound>     (1)
 *     --dd      0 | 1   duplicate dropping                 (0)
 *     --anytime 0 | 1   anytime search                     (0)
 *     --schedule <w0>,<w1>,...   the weights of ARA* (w, 1) and of
 *               ARPBNF with --algo pbnf               (none: wPBNF)
 *     --optimism <weight of optimistic and ios>           (2 w - 1)
 *     --focal   d | h | opt | ees  the order of focal      (d)
 *     --paths   0 | 1   print the paths of an anytime search on stderr (0)
 *     --frontier <units of pida>                           (128 per worker)
 *     --rule    astar | breadth | depth  (the frontier of pida)   (astar)
 *     --units   steal | shared           (the units of pida)      (steal)
 *     --async   1 for the asynchronous scheme of pida      (0)
 *     --bounds  0 to search every unit at every bound      (1)
 *     --split   <least slack of a node given away by pida> (0: none)
 *     --table   <levels of pida looked up in the frontier> (0)
 *     --det     1 for the deterministic path of pida       (0)
 *     --pdb     the tables of rubik (korf), pancake, burned and topspin
 *               (<k>: the first k objects, then the next k, ...), hanoi
 *               (<k>+<k>+...: the ranges of disks from the smallest)
 *     --cache   <directory of the tables>                  (none)
 *     --reverse 1 to search hanoi from the goal to the start  (0)
 *     --parts   <files of a bucket of ddd>                 (64)
 *     --dir     <directory of the files of ddd>     (the temporary one)
 *     --runtime thread | openmp | fastflow                 (thread)
 *     --dist    Z | AZ | A | P | H | R | HWD               (Z)
 *     --threads <p>                                        (1)
 *     --batch   <nodes per message>                        (64)
 *     --flush   <expansions between flushes>               (64)
 *     --open    buckets | heap | lowg (the heap with the ties of f to the
 *               smaller g)                                 (buckets)
 *     --block   <side of the blocks of AZH and AHDA*>   (100, msa 4)
 *     --minexp  <expansions before PBNF switches nblock>  (32)
 *     --nblock  <side of the nblocks of PBNF>        (20, msa 16)
 *     --eps     <inflation of the independence of PA*SE>   (1)
 *     --slow    <iterations of busy work per expansion>    (0)
 *     --limit   <nodes per worker>                         (none)
 *     --expect  <nodes expected per worker>               (unknown)
 *     --repeat  <runs of each instance>                    (1)
 *     --header  print the header line
 *
 * A file of the tiles domains holds one board per line, or a single
 * instance given by the cells of the tiles [see
 * SlidingTile::read_positions()]. The domain board reads a grid in the
 * format of Burns et al. (JAIR 2010) [see GridBoard], with its moves and
 * costs, the row 0 of its start and goal being the last line (bottom) or
 * the first one (top). The domain msa aligns the sequences of a FASTA file,
 * or the random ones of --random <seed>[:<n>[:<total>]] [see
 * SequenceAlignment::random_sequences()], scored by PAM250 or by the 625
 * scores of --matrix [see SequenceAlignment::read_scores()], with the gap
 * scores of the matrix or the linear gap penalty of --gap; --dist HWD is
 * its hyperplane work distribution [see HyperplaneDistribution], of the
 * thickness that --plane gives or that the rule of Kobayashi et al. gives
 * for --threads.
 *
 * The problems of the examples are sized by --size, and the instance of
 * the line is the size, or the size and the id after a slash: queens
 * <n>; vacuum <n>, with all the cells dirty and the robot in each cell
 * of --ids (the middle one by default, the cells counted from 0);
 * missionaries <m>:<c>:<b>; jugs <goal>:<capacity>:...; tsp <n>, with the
 * random cities of each seed of --ids (1 by default), or with the first
 * <n> cities of --file (all of them if --size is not given) [see
 * TravelingSalesman::read_cities()]; romania <from>:<to> (Arad:Bucharest
 * by default). The permutation puzzles are sized by --size too, among the
 * sizes compiled: pancake 12, 16, 20, 32, 64 and burned 8, 10, 12, 16 (a
 * random stack for each seed of --ids), topspin 12, 14, 16, 17 with a
 * turnstile of 4 (a random walk of 1000 moves for each seed); rubik takes
 * the scrambles of --file, or random walks of --size moves for each seed of
 * --ids, and hanoi the lines of --file (12, 16 or 20 disks, by --size, 20
 * by default). A road graph is read from the files <file>.gr and <file>.co,
 * and searched between the nodes of --points <s>,<t>, or between random
 * nodes for each seed of --ids, with the heuristic of the coordinates
 * (--open heap suits its large costs). The line ends with the weight, the
 * lower bound on the optimal cost and the number of paths found.
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

#include <cstdint>

#include <cstdio>

#include <cstdlib>

#include <fstream>

#include <iostream>

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
 #include "PBeFS/astar/SafePBNF.h"
 #ifdef PBEFS_HAS_OPENMP
  #include "ADS/parallel/OpenMPExecutor.h"
 #endif
 #include "PBeFS/astar/FarmAStar.h"
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
static std::map< std::string , std::string > opt = {
 { "domain" , "tiles15" } , { "file" , "" } , { "ids" , "" } ,
 { "random" , "" } , { "points" , "" } , { "algo" , "astar" } ,
 { "runtime" , "thread" } , { "dist" , "Z" } , { "threads" , "1" } ,
 { "batch" , "64" } , { "flush" , "64" } , { "open" , "buckets" } ,
 { "block" , "" } , { "limit" , "0" } , { "repeat" , "1" } ,
 { "minexp" , "32" } , { "nblock" , "" } , { "eps" , "1" } ,
 { "slow" , "0" } , { "expect" , "0" } , { "size" , "" } ,
 { "parts" , "64" } , { "dir" , "" } , { "origin" , "bottom" } ,
 { "matrix" , "" } , { "metric" , "euclidean" } , { "gap" , "0" } ,
 { "plane" , "auto" } , { "weight" , "1" } , { "dd" , "0" } ,
 { "anytime" , "0" } , { "schedule" , "" } , { "optimism" , "0" } ,
 { "focal" , "d" } , { "paths" , "0" } , { "frontier" , "0" } ,
 { "rule" , "astar" } , { "units" , "steal" } , { "async" , "0" } ,
 { "bounds" , "1" } , { "split" , "0" } , { "table" , "0" } ,
 { "det" , "0" } , { "pdb" , "" } , { "cache" , "" } , { "reverse" , "0" } };

static int rank = 0;                   // the MPI rank (0 without MPI)

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the indices in a list "a-b,c,..." (1-based), or all of [ 1 , n ]

static std::vector< int > parse_ids( const std::string & s , int n )
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

static int option( const std::string & o , int d )
{
 return( opt[ o ].empty() ? d : std::stoi( opt[ o ] ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the contents of the file of --file

static std::string read_file( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
 std::ostringstream os;
 os << in.rdbuf();
 return( os.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel options from the command line

static ParallelOptions parallel_options( void )
{
 ParallelOptions o;
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoul( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.node_limit = std::stoul( opt[ "limit" ] );
 o.expected_nodes = std::stoul( opt[ "expect" ] );
 o.weight = std::stod( opt[ "weight" ] );
 o.drop_duplicates = opt[ "dd" ] == "1";
 o.anytime = opt[ "anytime" ] == "1";
 return( o );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of ParallelIDAStar from the command line

static PIDAOptions pida_options( void )
{
 PIDAOptions o;
 static_cast< ParallelOptions & >( o ) = parallel_options();
 o.frontier = std::stoul( opt[ "frontier" ] );
 const std::string & r = opt[ "rule" ];
 if( ( r != "astar" ) && ( r != "breadth" ) && ( r != "depth" ) )
  throw( std::invalid_argument( "unknown frontier rule " + r ) );
 o.rule = r == "astar" ? FrontierRule::astar : r == "breadth"
                       ? FrontierRule::breadth : FrontierRule::depth;
 o.schedule = opt[ "units" ] == "shared" ? UnitSchedule::shared
                                            : UnitSchedule::steal;
 o.asynchronous = opt[ "async" ] != "0";
 o.unit_bounds = opt[ "bounds" ] != "0";
 o.split = std::stod( opt[ "split" ] );
 o.table_depth = unsigned( std::stoul( opt[ "table" ] ) );
 o.deterministic = opt[ "det" ] != "0";
 return( o );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// prints the line of a run

template< class S >
static void print( const std::string & instance , const S & s )
{
 if( rank != 0 )
  return;
 const auto & st = s.stats;
 const char * status[] = { "solved" , "unsolvable" , "node_limit" };
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
           << opt[ "weight" ] << ',' << s.bound << ','
           << s.improvements.size() << std::endl;
 if( opt[ "paths" ] == "1" )
  for( const auto & x : s.improvements )
   std::cerr << "path," << instance << ',' << x.cost << ',' << x.bound
             << ',' << x.expanded << ',' << x.seconds << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with HDA* and the distribution Dist on the chosen runtime

template< class D , class Dist >
static auto run_hda( const D & d , const Dist & dist ,
                     const ParallelOptions & o )
{
 #ifdef PBEFS_BENCH_MPI
  return( HDAStarMPI< D , Dist >( d , dist , o ).solve() );
 #else
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( HDAStar< D , Dist , ADS::ThreadPool >( d , dist , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( HDAStar< D , Dist , ADS::OpenMPExecutor >( d , dist , o ).solve() );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( HDAStar< D , Dist , ADS::FastFlowExecutor >( d , dist , o )
           .solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with the search, runtime and distribution of the command line

template< class D >
static Solution< state_t< D > , cost_t< D > > solve( const D & d )
{
 const std::string & algo = opt[ "algo" ];
 const ParallelOptions o = parallel_options();

 #ifdef PBEFS_BENCH_MPI
 if( algo == "pida" ) {
  if constexpr( std::is_trivially_copyable_v< state_t< D > > )
   return( ParallelIDAStarMPI< D >( d , pida_options() ).solve() );
  else
   throw( std::invalid_argument( "pida over MPI needs states that are "
                                 "trivially copyable" ) );
  }
 #else
 if( algo == "ida" )
  return( idastar( d , o ) );
 if( algo == "pida" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( ParallelIDAStar< D , ADS::ThreadPool >( d , pida_options() )
           .solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( ParallelIDAStar< D , ADS::OpenMPExecutor >( d , pida_options() )
           .solve() );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( ParallelIDAStar< D , ADS::FastFlowExecutor >( d , pida_options() )
           .solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "ddd" ) {
  if constexpr( DiskDomain< D > ) {
   DDDOptions x;
   x.partitions = std::stoul( opt[ "parts" ] );
   x.directory = opt[ "dir" ];
   return( ddd_astar( d , o , x ) );
   }
  else
   throw( std::invalid_argument( "ddd needs integer costs and states that"
                                 " are trivially copyable" ) );
  }
 const OpenKind ok = opt[ "open" ] == "heap" ? OpenKind::heap
                     : opt[ "open" ] == "lowg" ? OpenKind::heap_low_g
                                               : OpenKind::buckets;
 if( algo == "astar" )
  return( astar( d , ok , o ) );
 if( ( algo == "awa" ) || ( algo == "ara" ) ) {
  AnytimeOptions x;
  static_cast< SearchOptions & >( x ) = o;
  if( algo == "ara" ) {
   std::istringstream in( opt[ "schedule" ] );
   std::string w;
   while( std::getline( in , w , ',' ) )
    x.schedule.push_back( std::stod( w ) );
   if( x.schedule.empty() )
    x.schedule = { o.weight };
   }
  return( anytime_astar( d , x , ok ) );
  }
 if( ( algo == "optimistic" ) || ( algo == "ios" ) ) {
  OptimisticOptions x;
  static_cast< SearchOptions & >( x ) = o;
  x.optimism = std::stod( opt[ "optimism" ] );
  if( algo == "ios" )
   return( ImprovedOptimisticSearch< D >( d , x ).solve() );
  return( optimistic_search( d , x , ok ) );
  }
 if( algo == "focal" ) {
  FocalOptions x;
  static_cast< SearchOptions & >( x ) = o;
  const std::string & f = opt[ "focal" ];
  x.order = f == "h" ? FocalOrder::heuristic
            : f == "opt" ? FocalOrder::optimistic
            : f == "ees" ? FocalOrder::explicit_estimation
            : FocalOrder::distance;
  return( FocalSearch< D >( d , x ).solve() );
  }
 if( algo == "dps" )
  return( PotentialSearch< D >( d , o ).solve() );
 if( algo == "lpa" ) {
  if constexpr( PredecessorDomain< D > )
   return( LPAStar< D >( d , o ).solve() );
  throw( std::invalid_argument( "no predecessors for this domain" ) );
  }
 if( algo == "spa" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( SPAStar< D , ADS::ThreadPool >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( SPAStar< D , ADS::OpenMPExecutor >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "mqa" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( MQAStar< D , ADS::ThreadPool >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( MQAStar< D , ADS::OpenMPExecutor >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "pbnf" ) {
  if constexpr( NBlockDomain< D > ) {
   PBNFOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.min_expansions = unsigned( std::stoul( opt[ "minexp" ] ) );
   std::istringstream in( opt[ "schedule" ] );
   std::string sw;
   while( std::getline( in , sw , ',' ) )
    po.schedule.push_back( std::stod( sw ) );
   const std::string & r = opt[ "runtime" ];
   if( r == "thread" )
    return( SafePBNF< D , ADS::ThreadPool >( d , po ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( SafePBNF< D , ADS::OpenMPExecutor >( d , po ).solve() );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   }
  throw( std::invalid_argument( "no nblocks for this domain" ) );
  }
 if( algo == "farm" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( FarmAStar< D , ADS::ThreadPool >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( FarmAStar< D , ADS::OpenMPExecutor >( d , o ).solve() );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( FarmAStar< D , ADS::FastFlowExecutor >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "pase" ) {
  if constexpr( PairwiseDomain< D > ) {
   PASEOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.eps = std::stod( opt[ "eps" ] );
   const std::string & r = opt[ "runtime" ];
   if( r == "thread" )
    return( PASE< D , ADS::ThreadPool >( d , po ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( PASE< D , ADS::OpenMPExecutor >( d , po ).solve() );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   }
  throw( std::invalid_argument( "no pairwise heuristic for this domain" ) );
  }
 #endif
 if( algo != "hda" )
  throw( std::invalid_argument( "unknown search " + algo ) );

 const std::string & dist = opt[ "dist" ];
 if constexpr( FeatureDomain< D > )
  if( dist == "Z" )
   return( run_hda( d , ZobristDistribution< D >( d ) , o ) );
 if constexpr( AbstractFeatureDomain< D > )
  if( dist == "AZ" )
   return( run_hda( d , AbstractZobristDistribution< D >( d ) , o ) );
 if constexpr( AbstractStateDomain< D > )
  if( dist == "A" )
   return( run_hda( d , AbstractStateDistribution< D >( d ) , o ) );
 if( dist == "H" )
  return( run_hda( d , StateHashDistribution< D >( d ) , o ) );
 if constexpr( LevelDomain< D > )
  if( dist == "HWD" ) {
   auto [ planes , split ] = HyperplaneDistribution< D >::thickness(
                              d , o.threads );
   const std::string & t = opt[ "plane" ];
   if( t.rfind( "1/" , 0 ) == 0 ) {
    planes = 1;
    split = unsigned( std::stoul( t.substr( 2 ) ) );
    }
   else
    if( t != "auto" ) {
     planes = unsigned( std::stoul( t ) );
     split = 1;
     }
   return( run_hda( d , HyperplaneDistribution< D >( d , 1 , planes ,
                                                       split ) , o ) );
   }
 if constexpr( PerfectHashDomain< D > )
  if( dist == "P" )
   return( run_hda( d , PerfectHashDistribution< D >( d ) , o ) );
 #ifndef PBEFS_BENCH_MPI
 if( dist == "R" )
  return( run_hda( d , RandomDistribution< D >( d ) , o ) );
 #endif
 throw( std::invalid_argument( "distribution " + dist +
                               " not available for this domain" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d, slowed down by --slow iterations per expansion if not 0

template< class D >
static Solution< state_t< D > , cost_t< D > > solve_maybe_slow( const D & d )
{
 const unsigned long work = std::stoul( opt[ "slow" ] );
 if( work == 0 )
  return( solve( d ) );
 return( solve( SlowDomain< D >( d , work ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Domain d --repeat times, the line labelled by instance

template< class D >
static void run_domain( const std::string & instance , const D & d )
{
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int r = 0 ; r < repeat ; ++r )
  print( instance , solve_maybe_slow( d ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the sliding-tile instances of the file

template< class D >
static void run_tiles( void )
{
 const std::string text = read_file();
 if( text.find( "starting positions" ) != std::string::npos ) {
  std::istringstream is( text );
  const auto [ start , goal ] = D::read_positions( is );
  run_domain( opt[ "file" ] , D( start , goal ) );
  return;
  }
 std::istringstream in( text );
 std::vector< typename D::board_type > boards;
 std::string line;
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  boards.push_back( D::parse( line ) );
  }
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int id : parse_ids( opt[ "ids" ] , int( boards.size() ) ) ) {
  const D d( boards.at( std::size_t( id - 1 ) ) );
  for( int r = 0 ; r < repeat ; ++r )
   print( std::to_string( id ) , solve_maybe_slow( d ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a grid instance, from a Moving AI map or random

template< int Moves >
static void run_grid( void )
{
 GridMap m;
 std::string name;
 if( ! opt[ "random" ].empty() ) {
  int w , h , c = 1;
  double p;
  unsigned long long seed;
  if( std::sscanf( opt[ "random" ].c_str() , "%dx%d:%lf:%llu:%d" , & w , & h ,
                   & p , & seed , & c ) < 4 )
   throw( std::invalid_argument( "--random wants WxH:p:seed[:c]" ) );
  m = GridMap::random( w , h , p , seed , c );
  name = opt[ "random" ];
  }
 else {
  std::ifstream in( opt[ "file" ] );
  if( ! in )
   throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
  m = GridMap::read_movingai( in );
  name = opt[ "file" ];
  }
 int sx = 0 , sy = m.height() - 1 , gx = m.width() - 1 , gy = m.height() - 1;
 if( ! opt[ "points" ].empty() &&
     ( std::sscanf( opt[ "points" ].c_str() , "%d,%d,%d,%d" , & sx , & sy ,
                    & gx , & gy ) != 4 ) )
  throw( std::invalid_argument( "--points wants sx,sy,gx,gy" ) );
 const Grid< Moves > g( m , sx , sy , gx , gy , option( "block" , 100 ) ,
                        option( "nblock" , 20 ) );
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int r = 0 ; r < repeat ; ++r )
  print( name , solve_maybe_slow( g ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a grid in the format of Burns et al., with its moves

static void run_board( void )
{
 const std::string & o = opt[ "origin" ];
 if( ( o != "bottom" ) && ( o != "top" ) )
  throw( std::invalid_argument( "--origin wants bottom or top" ) );
 std::istringstream in( read_file() );
 const GridBoard b = GridBoard::read( in , o == "bottom" );
 const int block = option( "block" , 100 ) , nside = option( "nblock" , 20 );
 if( b.moves == 4 )
  run_domain( opt[ "file" ] , Grid< 4 >( b , block , nside ) );
 else
  run_domain( opt[ "file" ] , Grid< 8 >( b , block , nside ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the alignment of the sequences of the file

static void run_msa( void )
{
 std::vector< std::string > seqs;
 std::string name = opt[ "file" ];
 if( ! opt[ "random" ].empty() ) {
  unsigned long long seed;
  int n = 0 , total = 5000;
  if( std::sscanf( opt[ "random" ].c_str() , "%llu:%d:%d" , & seed , & n ,
                   & total ) < 1 )
   throw( std::invalid_argument( "--random wants seed[:n[:total]]" ) );
  seqs = SequenceAlignment<>::random_sequences( seed , n , total );
  name = opt[ "random" ];
  }
 else {
  std::istringstream in( read_file() );
  seqs = SequenceAlignment<>::read_fasta( in );
  }
 ScoreMatrix scores = pam250;
 if( ! opt[ "matrix" ].empty() ) {
  std::ifstream ms( opt[ "matrix" ] );
  if( ! ms )
   throw( std::invalid_argument( "cannot read " + opt[ "matrix" ] ) );
  scores = SequenceAlignment<>::read_scores( ms );
  }
 run_domain( name , SequenceAlignment<>(
              seqs , scores , option( "block" , 4 ) ,
              option( "nblock" , 16 ) , std::stoi( opt[ "gap" ] ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the fields of --size separated by ':', at least one

static std::vector< std::string > size_fields( void )
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

static std::vector< std::vector< std::uint8_t > > pdb_patterns( unsigned n )
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
static void run_permutation( const std::string & instance , const D & d ,
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the pancakes of N, burned or not, a random stack for each seed

template< int N , bool Burned >
static void run_pancake( void )
{
 const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 );
 for( int seed : ids ) {
  const std::string name = std::to_string( N ) + "/" + std::to_string( seed );
  if constexpr( Burned ) {
   const BurnedPancakePuzzle< N > d( BurnedPancakePuzzle< N >::random(
                                      std::uint64_t( seed ) ) );
   run_permutation( name , d , N );
   }
  else {
   const PancakePuzzle< N > d( PancakePuzzle< N >::random(
                                std::uint64_t( seed ) ) ,
                               unsigned( std::stoul( opt[ "gap" ] ) ) );
   run_permutation( name , d , N );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs TopSpin( N , 4 ), a random walk of 1000 moves for each seed

template< int N >
static void run_topspin( void )
{
 const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 );
 for( int seed : ids ) {
  const TopSpin< N , 4 > d( TopSpin< N , 4 >::random( std::uint64_t( seed ) ,
                                                      1000 ) );
  run_permutation( std::to_string( N ) + "/" + std::to_string( seed ) , d ,
                   N );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs Rubik's cube: the scrambles of --file, or random walks

static void run_rubik( void )
{
 std::vector< std::pair< std::string , RubiksCube::state_type > > v;
 if( ! opt[ "file" ].empty() ) {
  std::ifstream in( opt[ "file" ] );
  if( ! in )
   throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
  std::vector< RubiksCube::state_type > all;
  for( std::string line ; std::getline( in , line ) ; )
   if( ! line.empty() && ( line[ 0 ] != '#' ) )
    all.push_back( RubiksCube::scramble( line ) );
  for( int id : parse_ids( opt[ "ids" ] , int( all.size() ) ) )
   v.emplace_back( std::to_string( id ) , all.at( std::size_t( id - 1 ) ) );
  }
 else {
  const int moves = std::stoi( size_fields()[ 0 ] );
  const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                        : parse_ids( opt[ "ids" ] , 0 );
  for( int seed : ids )
   v.emplace_back( std::to_string( moves ) + "/" + std::to_string( seed ) ,
                   RubiksCube::random_walk( moves , std::uint64_t( seed ) ) );
  }
 const RubiksCube base;
 std::unique_ptr< RubikPDBHeuristic > h;
 if( opt[ "pdb" ] == "korf" )
  h = std::make_unique< RubikPDBHeuristic >(
       base , unsigned( std::stoul( opt[ "threads" ] ) ) , opt[ "cache" ] );
 else
  if( ! opt[ "pdb" ].empty() )
   throw( std::invalid_argument( "rubik knows --pdb korf only" ) );
 for( const auto & [ name , s ] : v ) {
  const RubiksCube d( s );
  if( h )
   run_domain( name , PDBDomain< RubiksCube , RubikPDBHeuristic >( d , *h ) );
  else
   run_domain( name , d );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Towers of Hanoi of D disks on the lines of --file, with the
// additive tables of --pdb (the ranges of disks, k+k+...)

template< int D >
static void run_hanoi( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
 std::vector< std::pair< std::uint64_t , std::uint64_t > > all;
 for( std::string line ; std::getline( in , line ) ; ) {
  if( line.empty() || ( line[ 0 ] != '(' ) )
   continue;
  std::istringstream is( line );
  const auto s = TowersOfHanoi< D >::read( is );
  all.emplace_back( s , TowersOfHanoi< D >::read( is ) );
  }
 std::vector< unsigned > ks;
 {
  std::istringstream is( opt[ "pdb" ].empty() ? std::to_string( D / 2 ) +
                         "+" + std::to_string( D - D / 2 ) : opt[ "pdb" ] );
  for( std::string part ; std::getline( is , part , '+' ) ; )
   ks.push_back( unsigned( std::stoul( part ) ) );
 }
 for( int id : parse_ids( opt[ "ids" ] , int( all.size() ) ) ) {
  auto [ s , g ] = all.at( std::size_t( id - 1 ) );
  if( opt[ "reverse" ] != "0" )
   std::swap( s , g );
  const TowersOfHanoi< D > d( s , g );
  const HanoiPDB h( ks , g , unsigned( std::stoul( opt[ "threads" ] ) ) ,
                    opt[ "cache" ] );
  run_domain( std::to_string( id ) , PDBDomain< TowersOfHanoi< D > ,
                                               HanoiPDB >( d , h ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a road graph of DIMACS, between --points or random pairs

static void run_road( void )
{
 std::ifstream gr( opt[ "file" ] + ".gr" ) , co( opt[ "file" ] + ".co" );
 if( ! gr || ! co )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] +
                                ".gr and .co" ) );
 auto g = read_dimacs_graph< std::int64_t >( gr );
 const auto c = read_dimacs_coordinates( co , g.size() );
 const double kappa = dimacs_scale( g , c );
 std::vector< std::pair< std::uint32_t , std::uint32_t > > pairs;
 if( ! opt[ "points" ].empty() ) {
  unsigned long s , t;
  if( std::sscanf( opt[ "points" ].c_str() , "%lu,%lu" , & s , & t ) != 2 )
   throw( std::invalid_argument( "--points wants s,t for a road graph" ) );
  pairs.emplace_back( std::uint32_t( s - 1 ) , std::uint32_t( t - 1 ) );
  }
 else
  for( int seed : opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 ) ) {
   std::mt19937_64 gen{ std::uint64_t( seed ) };
   const auto s = std::uint32_t( random_below( gen , g.size() ) );
   pairs.emplace_back( s , std::uint32_t( random_below( gen , g.size() ) ) );
   }
 for( const auto & [ s , t ] : pairs ) {
  g.set_source( s );
  g.set_target( t );
  set_dimacs_heuristic( g , c , kappa );
  run_domain( std::to_string( s + 1 ) + ":" + std::to_string( t + 1 ) , g );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the pancakes, the burned pancakes or TopSpin of the size --size

static void run_permutation_sized( void )
{
 const std::string & d = opt[ "domain" ];
 const int n = std::stoi( size_fields()[ 0 ] );
 if( d == "pancake" )
  switch( n ) {
   case 12: run_pancake< 12 , false >(); return;
   case 16: run_pancake< 16 , false >(); return;
   case 20: run_pancake< 20 , false >(); return;
   case 32: run_pancake< 32 , false >(); return;
   case 64: run_pancake< 64 , false >(); return;
   default: throw( std::invalid_argument( "pancake: 12, 16, 20, 32, 64" ) );
   }
 if( d == "burned" )
  switch( n ) {
   case 8: run_pancake< 8 , true >(); return;
   case 10: run_pancake< 10 , true >(); return;
   case 12: run_pancake< 12 , true >(); return;
   case 16: run_pancake< 16 , true >(); return;
   default: throw( std::invalid_argument( "burned: 8, 10, 12, 16" ) );
   }
 switch( n ) {
  case 12: run_topspin< 12 >(); return;
  case 14: run_topspin< 14 >(); return;
  case 16: run_topspin< 16 >(); return;
  case 17: run_topspin< 17 >(); return;
  default: throw( std::invalid_argument( "topspin: 12, 14, 16, 17" ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Towers of Hanoi of --size disks (20 by default)

static void run_hanoi_sized( void )
{
 switch( opt[ "size" ].empty() ? 20 : std::stoi( opt[ "size" ] ) ) {
  case 12: run_hanoi< 12 >(); return;
  case 16: run_hanoi< 16 >(); return;
  case 20: run_hanoi< 20 >(); return;
  default: throw( std::invalid_argument( "hanoi: 12, 16, 20 disks" ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a problem of the examples, sized by --size

static void run_example( void )
{
 const std::string & dom = opt[ "domain" ];
 if( ( dom == "tsp" ) && ! opt[ "file" ].empty() ) {
  const std::string & m = opt[ "metric" ];
  if( ( m != "euclidean" ) && ( m != "manhattan" ) )
   throw( std::invalid_argument( "--metric wants euclidean or manhattan" ) );
  std::istringstream in( read_file() );
  run_domain( opt[ "file" ] + "/" + opt[ "size" ] , TravelingSalesman(
               TravelingSalesman::read_cities(
                in , m == "euclidean" ? TravelingSalesman::Metric::euclidean
                                     : TravelingSalesman::Metric::manhattan ,
                option( "size" , 0 ) ) ) );
  return;
  }
 if( dom == "romania" ) {
  if( opt[ "size" ].empty() )
   opt[ "size" ] = "Arad:Bucharest";
  const auto f = size_fields();
  const auto from = Romania::find( f[ 0 ] );
  const auto to = Romania::find( f.size() > 1 ? f[ 1 ] : "Bucharest" );
  if( ( from == Romania::NCities ) || ( to == Romania::NCities ) )
   throw( std::invalid_argument( "unknown city in " + opt[ "size" ] ) );
  run_domain( opt[ "size" ] , Romania( from , to ) );
  return;
  }
 const auto f = size_fields();
 const int n = std::stoi( f[ 0 ] );
 if( dom == "queens" )
  run_domain( opt[ "size" ] , NQueens( n ) );
 else
  if( dom == "vacuum" ) {
   const std::uint32_t dirty = n >= 32 ? ~std::uint32_t( 0 )
                                       : ( std::uint32_t( 1 ) << n ) - 1;
   const auto ids = opt[ "ids" ].empty()
                    ? std::vector< int >{ n / 2 } : parse_ids( opt[ "ids" ] ,
                                                               0 );
   for( int start : ids )
    run_domain( opt[ "size" ] + "/" + std::to_string( start ) ,
                VacuumWorld( n , start , dirty ) );
   }
  else
   if( dom == "missionaries" ) {
    if( f.size() != 3 )
     throw( std::invalid_argument( "--size wants <m>:<c>:<b>" ) );
    run_domain( opt[ "size" ] , MissionariesCannibals(
                 n , std::stoi( f[ 1 ] ) , std::stoi( f[ 2 ] ) ) );
    }
   else
    if( dom == "jugs" ) {
     std::vector< int > cap;
     for( std::size_t i = 1 ; i < f.size() ; ++i )
      cap.push_back( std::stoi( f[ i ] ) );
     run_domain( opt[ "size" ] , WaterJugs( cap , n ) );
     }
    else
     if( dom == "tsp" ) {
      const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                           : parse_ids( opt[ "ids" ] , 0 );
      for( int seed : ids )
       run_domain( opt[ "size" ] + "/" + std::to_string( seed ) ,
                   TravelingSalesman( TravelingSalesman::random_cities(
                                       n , std::uint64_t( seed ) ) ) );
      }
     else
      throw( std::invalid_argument( "unknown domain " + dom ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 #ifdef PBEFS_BENCH_MPI
  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  opt[ "algo" ] = "hda";
  opt[ "runtime" ] = "mpi";
 #endif

 int ret = 0;
 try {
  bool header = false;
  for( int i = 1 ; i < argc ; ++i ) {
   std::string a = argv[ i ];
   if( a == "--header" ) {
    header = true;
    continue;
    }
   if( ( a.rfind( "--" , 0 ) != 0 ) || ! opt.count( a.substr( 2 ) ) ||
       ( i + 1 >= argc ) )
    throw( std::invalid_argument( "bad option " + a ) );
   opt[ a.substr( 2 ) ] = argv[ ++i ];
   }
  #ifdef PBEFS_BENCH_MPI
   int np;
   MPI_Comm_size( MPI_COMM_WORLD , & np );
   opt[ "threads" ] = std::to_string( np );
  #endif
  if( header && ( rank == 0 ) )
   std::cout << "domain,instance,algo,runtime,dist,threads,batch,flush,open,"
                "status,cost,expanded,reexpanded,generated,duplicates,sent,"
                "messages,stored,LB,CO,seconds,weight,bound,paths"
             << std::endl;

  const std::string & d = opt[ "domain" ];
  if( d == "tiles15" )
   run_tiles< FifteenPuzzle >();
  else
   if( d == "tiles8" )
    run_tiles< EightPuzzle >();
   else
    if( d == "tiles24" )
     run_tiles< TwentyFourPuzzle >();
    else
     if( d == "grid4" )
      run_grid< 4 >();
     else
      if( d == "grid8" )
       run_grid< 8 >();
      else
       if( d == "board" )
        run_board();
       else
        if( d == "msa" )
         run_msa();
        else
         if( ( d == "pancake" ) || ( d == "burned" ) || ( d == "topspin" ) )
          run_permutation_sized();
         else
          if( d == "rubik" )
           run_rubik();
          else
           if( d == "hanoi" )
            run_hanoi_sized();
           else
            if( d == "road" )
             run_road();
            else
             run_example();
  }
 catch( std::exception & e ) {
  std::cerr << "AStarBench: " << e.what() << std::endl;
  ret = 2;
  }

 #ifdef PBEFS_BENCH_MPI
  MPI_Finalize();
 #endif
 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File AStarBench.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
