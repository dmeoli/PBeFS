/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndBoundBench.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the Branch-and-Bound (B&B) of PBeFS: solves random 0-1
 * knapsacks [see Knapsack.h] or random traveling salesperson problems with
 * a chosen scheme, runtime, selection and number of workers, and prints
 * one line of comma-separated values per run, with the columns of
 * AStarBench.
 *
 * Compiled with PBEFS_BENCH_MPI defined, it is the benchmark of
 * MPIBranchAndBoundSearch, to be started under mpiexec; only the process
 * of rank 0 prints.
 *
 *     BranchAndBoundBench [options]
 *
 *     --domain  knapsack | tsp                             (knapsack)
 *     --size    <class>:<n>:<r> (knapsack) | <n> (tsp)
 *     --ids     <a-b,c,...>            the seeds of the instances (1)
 *     --algo    sequential | shared | distributed | deterministic |
 *               selfsplit | stealing | farm               (sequential)
 *     --runtime thread | openmp | fastflow                 (thread)
 *     --select  best_first | depth_first | breadth_first |
 *               best_estimate | plunging                   (best_first)
 *     --threads <p>                                        (1)
 *     --batch   <nodes per message>                        (16)
 *     --flush   <evaluations between flushes>              (1)
 *     --round   <nodes per round of the deterministic B&B> (0, 1 each)
 *     --split   <open nodes per worker of SelfSplit>       (0, i.e., 4)
 *     --limit   <evaluated nodes>                          (none)
 *     --slow    <iterations of busy work per evaluation>   (0)
 *     --repeat  <runs of each instance>                    (1)
 *     --header  print the header line
 *
 * The knapsack of class 0, 1 or 2 of Pisinger (uncorrelated, weakly and
 * strongly correlated) has n items of weight at most r and the capacity
 * half the total weight; the traveling salesperson problem has n random
 * cities of the unit square [see TravelingSalesman.h], and its B&B is the
 * tree of the partial tours of the Domain, each bounded by its cost plus
 * the minimum spanning tree heuristic. The columns are those of
 * AStarBench, so that the same scripts read both: dist is the selection,
 * open is empty, status is solved for an optimal value, cost the optimal
 * value, expanded the nodes evaluated, reexpanded those branched,
 * duplicates the pruned ones, stored the most open nodes of a queue, LB the
 * largest evaluations of a worker over their mean, and CO the nodes sent
 * to another worker over those generated. With --slow, each evaluation
 * does the busy work of SlowDomain [see SlowDomain.h] first, as a bound
 * computed by an expensive relaxation would.
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

#include <algorithm>

#include <cstdint>

#include <cstdlib>

#include <iostream>

#include <map>

#include <numeric>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/TravelingSalesman.h"

#include "PBeFS/bnb/Knapsack.h"

#ifdef PBEFS_BENCH_MPI
 #include "PBeFS/bnb/MPIBranchAndBound.h"
#else
 #include "PBeFS/bnb/BranchAndBoundScheme.h"
 #ifdef PBEFS_HAS_OPENMP
  #include "PBeFS/core/OpenMPRuntime.h"
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  #include "PBeFS/bnb/FarmBranchAndBound.h"
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
 { "domain" , "knapsack" } , { "size" , "" } , { "ids" , "1" } ,
 { "algo" , "sequential" } , { "runtime" , "thread" } ,
 { "select" , "best_first" } , { "threads" , "1" } , { "batch" , "16" } ,
 { "flush" , "1" } , { "round" , "0" } , { "split" , "0" } ,
 { "limit" , "0" } , { "repeat" , "1" } , { "slow" , "0" } };

static int rank = 0;                   // the MPI rank (0 without MPI)

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS TourTree -----------------------------*/
/*--------------------------------------------------------------------------*/
// the B&B of a Domain whose search space is a tree, here the partial tours
// of the TravelingSalesman: a node is a state with its cost g and its
// heuristic h, its bound is g + h, a goal is feasible with value g, and its
// children are the successors of the state, each estimated by its own
// g + h; the nodes are trivially copyable, and travel between processes as
// their bytes

struct TourTree {

 struct node_type {
  std::uint64_t s;                 // the state
  double g;                        // the cost of the partial tour
  double h;                        // the heuristic of s
  };

 static constexpr Sense sense = Sense::minimize;

 const TravelingSalesman & d;      // the Domain

 node_type root( void ) const {
  const auto s = d.initial_state();
  return( node_type{ s , 0 , d.heuristic( s ) } );
  }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e;
  e.bound = n.g + n.h;
  if( d.is_goal( n.s ) )
   e.feasible = n.g;
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  d.for_each_successor( n.s , n.h , [ & ]( std::uint64_t t , double c ,
                                           double ht ) {
   f( node_type{ t , n.g + c , ht } , n.g + c + ht );
   } );
  }

 std::uint64_t hash( const node_type & n ) const {
  return( mix64( n.s ) );
  }

 double estimate( const node_type & n ) const {
  return( n.g + n.h );
  }

 };  // end( struct( TourTree ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SlowProblem ----------------------------*/
/*--------------------------------------------------------------------------*/
// the B&B problem Q whose evaluations cost work iterations of busy work
// more, the computation of SlowDomain

template< class Q >
struct SlowProblem {

 using node_type = typename Q::node_type;

 static constexpr Sense sense = Q::sense;

 const Q & q;                      // the problem
 unsigned long work;               // the iterations per evaluation

 node_type root( void ) const { return( q.root() ); }

 Evaluation evaluate( node_type & n ) const {
  std::uint64_t i = 7;
  for( unsigned long k = 0 ; k < work ; ++k )
   i = ( 11 * i ) % 9999943;
  c_sink += i;
  return( q.evaluate( n ) );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  q.branch( n , std::forward< F >( f ) );
  }

 std::uint64_t hash( const node_type & n ) const { return( q.hash( n ) ); }

 double estimate( const node_type & n ) const {
  return( q.estimate( n ) );
  }

 void serialize( const node_type & n , std::vector< char > & out ) const
  requires requires( const Q & x , const node_type & m ,
                     std::vector< char > & o ) { x.serialize( m , o ); } {
  q.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const
  requires requires( const Q & x , const char * & i ) {
                     x.deserialize( i ); } {
  return( q.deserialize( in ) );
  }

 inline static thread_local std::uint64_t c_sink = 0;  // keeps the work

 };  // end( struct( SlowProblem ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the indices in a list "a-b,c,..."

static std::vector< int > parse_ids( const std::string & s )
{
 std::vector< int > v;
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
// the fields of --size separated by ':'

static std::vector< std::string > size_fields( void )
{
 std::vector< std::string > v;
 std::istringstream in( opt[ "size" ] );
 std::string part;
 while( std::getline( in , part , ':' ) )
  v.push_back( part );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options from the command line

static ParallelBBOptions bb_options( void )
{
 ParallelBBOptions o;
 const std::string & sel = opt[ "select" ];
 if( sel == "best_first" )
  o.selection = Selection::best_first;
 else
  if( sel == "depth_first" )
   o.selection = Selection::depth_first;
  else
   if( sel == "breadth_first" )
    o.selection = Selection::breadth_first;
   else
    if( sel == "best_estimate" )
     o.selection = Selection::best_estimate;
    else
     if( sel == "plunging" )
      o.selection = Selection::plunging;
     else
      throw( std::invalid_argument( "unknown selection " + sel ) );
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoul( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.round = std::stoul( opt[ "round" ] );
 o.split = std::stoul( opt[ "split" ] );
 o.node_limit = std::stoul( opt[ "limit" ] );
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
 const char * status[] = { "solved" , "unsolvable" , "node_limit" ,
                           "time_limit" };
 double lb = 1;
 if( ! st.evaluated_by.empty() ) {
  const double sum = std::accumulate( st.evaluated_by.begin() ,
                                      st.evaluated_by.end() , 0.0 );
  if( sum > 0 )
   lb = double( *std::max_element( st.evaluated_by.begin() ,
                                   st.evaluated_by.end() ) ) *
        double( st.evaluated_by.size() ) / sum;
  }
 const double co = st.generated ? double( st.sent ) / double( st.generated )
                                : 0;
 std::cout.precision( 12 );
 std::cout << opt[ "domain" ] << ',' << instance << ',' << opt[ "algo" ]
           << ',' << opt[ "runtime" ] << ',' << opt[ "select" ] << ','
           << opt[ "threads" ] << ',' << opt[ "batch" ] << ','
           << opt[ "flush" ] << ",," << status[ int( s.status ) ] << ','
           << s.value << ',' << st.evaluated << ',' << st.branched << ','
           << st.generated << ','
           << st.pruned_generated + st.pruned_extracted + st.pruned_evaluated
           << ',' << st.sent << ',' << st.messages << ',' << st.max_open
           << ',' << lb << ',' << co << ',' << st.seconds << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves q with the scheme, runtime and selection of the command line

template< class Q >
static BBSolution< typename Q::node_type > solve( const Q & q )
{
 const ParallelBBOptions o = bb_options();
 #ifdef PBEFS_BENCH_MPI
  return( MPIBranchAndBoundSearch< Q >( q , o ).solve() );
 #else
  const std::string & algo = opt[ "algo" ];
  #ifdef PBEFS_HAS_FASTFLOW
  if( algo == "farm" )
   return( FarmBranchAndBound< Q >( q , o ).solve() );
  #endif
  const BBScheme s = bb_scheme( algo );
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( parallel_branch_and_bound< Q , ThreadRuntime >( q , s , o ) );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( parallel_branch_and_bound< Q , OpenMPRuntime >( q , s , o ) );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( parallel_branch_and_bound< Q , FastFlowRuntime >( q , s , o ) );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs q --repeat times, the line labelled by instance

template< class Q >
static void run( const std::string & instance , const Q & q )
{
 const int repeat = std::stoi( opt[ "repeat" ] );
 const unsigned long work = std::stoul( opt[ "slow" ] );
 for( int r = 0 ; r < repeat ; ++r )
  if( work == 0 )
   print( instance , solve( q ) );
  else
   print( instance , solve( SlowProblem< Q >{ q , work } ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 #ifdef PBEFS_BENCH_MPI
  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  opt[ "algo" ] = "mpi";
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

  const auto f = size_fields();
  const std::string & d = opt[ "domain" ];
  for( int seed : parse_ids( opt[ "ids" ] ) ) {
   const std::string instance = opt[ "size" ] + "/" + std::to_string( seed );
   if( d == "knapsack" ) {
    if( f.size() != 3 )
     throw( std::invalid_argument( "--size wants <class>:<n>:<r>" ) );
    run( instance , Knapsack::random( std::stoi( f[ 0 ] ) ,
                                      std::stoul( f[ 1 ] ) ,
                                      std::stol( f[ 2 ] ) ,
                                      std::uint64_t( seed ) ) );
    }
   else
    if( d == "tsp" ) {
     if( f.size() != 1 )
      throw( std::invalid_argument( "--size wants <n>" ) );
     const TravelingSalesman t( TravelingSalesman::random_cities(
                                 std::stoi( f[ 0 ] ) ,
                                 std::uint64_t( seed ) ) );
     run( instance , TourTree{ t } );
     }
    else
     throw( std::invalid_argument( "unknown domain " + d ) );
   }
  }
 catch( std::exception & e ) {
  std::cerr << "BranchAndBoundBench: " << e.what() << std::endl;
  ret = 2;
  }

 #ifdef PBEFS_BENCH_MPI
  MPI_Finalize();
 #endif
 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File BranchAndBoundBench.cpp --------------------*/
/*--------------------------------------------------------------------------*/
