/*--------------------------------------------------------------------------*/
/*-------------------------- File AStarBench.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the searches of PBeFS: solves instances of the
 * 15-puzzle, of the 8-puzzle, of grid pathfinding or of the problems of
 * the examples with a chosen search, runtime, work distribution and number
 * of workers, and prints one line of comma-separated values per run.
 *
 * Compiled with PBEFS_BENCH_MPI defined, it is the benchmark of HDAStarMPI,
 * to be started under mpiexec; only the process of rank 0 prints.
 *
 *     AStarBench [options]
 *
 *     --domain  tiles15 | tiles8 | grid4 | grid8 | queens | vacuum |
 *               missionaries | jugs | tsp | romania           (tiles15)
 *     --file    <instances or map>     korf100.txt, or a Moving AI map
 *     --ids     <a-b,c,...>            instances of the file (all)
 *     --random  <W>x<H>:<p>:<seed>[:<c>]  a random grid instead of a map,
 *                                      its free cells costing 1 to c (1)
 *     --size    <size>                 the size of an example (below)
 *     --points  <sx>,<sy>,<gx>,<gy>    start and goal (bottom corners)
 *     --algo    astar | spa | hda | pbnf | pase | farm     (astar)
 *     --runtime thread | openmp | fastflow                 (thread)
 *     --dist    Z | AZ | A | P | H | R                     (Z)
 *     --threads <p>                                        (1)
 *     --batch   <nodes per message>                        (64)
 *     --flush   <expansions between flushes>               (64)
 *     --open    buckets | heap                             (buckets)
 *     --block   <side of the grid blocks of AZH and AHDA*> (100)
 *     --minexp  <expansions before PBNF switches nblock>  (32)
 *     --nblock  <side of the grid nblocks of PBNF>         (20)
 *     --eps     <inflation of the independence of PA*SE>   (1)
 *     --slow    <iterations of busy work per expansion>    (0)
 *     --limit   <nodes per worker>                         (none)
 *     --expect  <nodes expected per worker>               (unknown)
 *     --repeat  <runs of each instance>                    (1)
 *     --header  print the header line
 *
 * The problems of the examples are sized by --size, and the instance of
 * the line is the size, or the size and the id after a slash: queens
 * <n>; vacuum <n>, with all the cells dirty and the robot in each cell
 * of --ids (the middle one by default, the cells counted from 0);
 * missionaries <m>:<c>:<b>; jugs <goal>:<capacity>:...; tsp <n>, with the
 * random cities of each seed of --ids (1 by default); romania
 * <from>:<to> (Arad:Bucharest by default).
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
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

#include <sstream>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/MissionariesCannibals.h"

#include "PBeFS/astar/NQueens.h"

#include "PBeFS/astar/Romania.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/SlowDomain.h"

#include "PBeFS/astar/TravelingSalesman.h"

#include "PBeFS/astar/VacuumWorld.h"

#include "PBeFS/astar/WaterJugs.h"

#ifdef PBEFS_BENCH_MPI
 #include "PBeFS/astar/HDAStarMPI.h"
#else
 #include "PBeFS/astar/HDAStar.h"
 #include "PBeFS/astar/SPAStar.h"
 #include "PBeFS/astar/PASE.h"
 #include "PBeFS/astar/SafePBNF.h"
 #ifdef PBEFS_HAS_OPENMP
  #include "PBeFS/core/OpenMPRuntime.h"
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  #include "PBeFS/astar/FarmAStar.h"
  #include "PBeFS/core/FastFlowRuntime.h"
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
 { "block" , "100" } , { "limit" , "0" } , { "repeat" , "1" } ,
 { "minexp" , "32" } , { "nblock" , "20" } , { "eps" , "1" } ,
 { "slow" , "0" } , { "expect" , "0" } , { "size" , "" } };

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
// the parallel options from the command line

static ParallelOptions parallel_options( void )
{
 ParallelOptions o;
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoul( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.node_limit = std::stoul( opt[ "limit" ] );
 o.expected_nodes = std::stoul( opt[ "expect" ] );
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
           << st.communication_overhead() << ',' << st.seconds << std::endl;
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
   return( HDAStar< D , Dist , ThreadRuntime >( d , dist , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( HDAStar< D , Dist , OpenMPRuntime >( d , dist , o ).solve() );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( HDAStar< D , Dist , FastFlowRuntime >( d , dist , o ).solve() );
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

 #ifndef PBEFS_BENCH_MPI
 if( algo == "astar" )
  return( astar( d , opt[ "open" ] == "heap" ? OpenKind::heap
                                             : OpenKind::buckets , o ) );
 if( algo == "spa" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( SPAStar< D , ThreadRuntime >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( SPAStar< D , OpenMPRuntime >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "pbnf" ) {
  if constexpr( NBlockDomain< D > ) {
   PBNFOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.min_expansions = unsigned( std::stoul( opt[ "minexp" ] ) );
   const std::string & r = opt[ "runtime" ];
   if( r == "thread" )
    return( SafePBNF< D , ThreadRuntime >( d , po ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( SafePBNF< D , OpenMPRuntime >( d , po ).solve() );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   }
  throw( std::invalid_argument( "no nblocks for this domain" ) );
  }
 #ifdef PBEFS_HAS_FASTFLOW
 if( algo == "farm" )
  return( FarmAStar< D >( d , o ).solve() );
 #endif
 if( algo == "pase" ) {
  if constexpr( PairwiseDomain< D > ) {
   PASEOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.eps = std::stod( opt[ "eps" ] );
   const std::string & r = opt[ "runtime" ];
   if( r == "thread" )
    return( PASE< D , ThreadRuntime >( d , po ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( PASE< D , OpenMPRuntime >( d , po ).solve() );
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
// runs the sliding-tile instances of the file

template< class D >
static void run_tiles( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
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
 const Grid< Moves > g( m , sx , sy , gx , gy , std::stoi( opt[ "block" ] ) ,
                        std::stoi( opt[ "nblock" ] ) );
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int r = 0 ; r < repeat ; ++r )
  print( name , solve_maybe_slow( g ) );
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
// runs the Domain d --repeat times, the line labelled by instance

template< class D >
static void run_domain( const std::string & instance , const D & d )
{
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int r = 0 ; r < repeat ; ++r )
  print( instance , solve_maybe_slow( d ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a problem of the examples, sized by --size

static void run_example( void )
{
 const std::string & dom = opt[ "domain" ];
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
                "messages,stored,LB,CO,seconds" << std::endl;

  const std::string & d = opt[ "domain" ];
  if( d == "tiles15" )
   run_tiles< FifteenPuzzle >();
  else
   if( d == "tiles8" )
    run_tiles< EightPuzzle >();
   else
    if( d == "grid4" )
     run_grid< 4 >();
    else
     if( d == "grid8" )
      run_grid< 8 >();
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
