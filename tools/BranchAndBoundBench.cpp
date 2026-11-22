/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndBoundBench.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the Branch-and-Bound (B&B) of PBeFS: solves random 0-1
 * knapsacks [see Knapsack.h], random traveling salesperson problems, the
 * trees of UTS [see UTS.h], the maximum cliques of DIMACS graphs [see
 * MaxClique.h], the quadratic assignments of QAPLIB [see QAP.h], the
 * flowshops of Taillard [see Flowshop.h], and the subgraph isomorphisms
 * and the common subgraphs of random graphs [see SubgraphIsomorphism.h and
 * MaxCommonSubgraph.h], with a chosen scheme, runtime, selection and number
 * of workers, and prints one line of comma-separated values per run, with
 * the columns of AStarBench.
 *
 * Compiled with PBEFS_BENCH_MPI defined, it is the benchmark of
 * MPIBranchAndBoundSearch (--algo mpi) and of MPIHubBranchAndBoundSearch
 * (--algo hub), to be started under mpiexec; only the process of rank 0
 * prints.
 *
 *     BranchAndBoundBench [options]
 *
 *     --domain  knapsack | tsp | uts | clique | qap | flowshop | sip |
 *               mcs                                        (knapsack)
 *     --size    <class>:<n>:<r> (knapsack) | <n> (tsp) | t1 | t2 | t3 |
 *               geo:<shape>:<b>:<g>:<seed> | bin:<b>:<m>:<q>:<seed> (uts)
 *               | <np>:<pp>:<nt>:<pt> (sip, mcs: sizes and percents of
 *               the edges of the two random graphs)
 *     --file    <file of a DIMACS graph or of QAPLIB> (clique, qap)
 *     --values  double | long     the values of the knapsack (double)
 *     --ids     <a-b,c,...>  the seeds, or the instances of Taillard (1)
 *     --algo    sequential | shared | distributed | deterministic |
 *               selfsplit | stealing | multiqueue | depthbounded |
 *               budget | stackstealing | ordered | farm | ivm
 *               (sequential), or mpi | hub with MPI
 *     --runtime thread | openmp | fastflow                 (thread)
 *     --select  best_first | depth_first | breadth_first |
 *               best_estimate | plunging                   (best_first)
 *     --threads <p>                                        (1)
 *     --batch   <nodes per message>                        (16)
 *     --flush   <evaluations between flushes>              (1)
 *     --round   <nodes per round of the deterministic B&B> (0, 1 each)
 *     --split   <open nodes per worker of SelfSplit>       (0, i.e., 4)
 *     --dive    0 | 1: the dive that ends the SelfSplit sampling     (1)
 *     --share   0 | 1: the incumbent shared by the SelfSplit workers (0)
 *     --spawn   <spawn depth of depthbounded and ordered>  (2)
 *     --budget  <backtracks of budget>                     (10000)
 *     --stealall 0 | 1: a steal takes all of a level      (0)
 *     --discrepancy 0 | 1: ordered by discrepancy          (1)
 *     --indexed 0 | 1: the tasks as paths                  (0)
 *     --race    <selection,...>: the settings of a race    (none)
 *     --racenodes <evaluations of each racer>              (0, 100 p)
 *     --mode    optimize | decide | enumerate              (optimize)
 *     --target  <value sought by decide>                   (none)
 *     --incumbent <value> | opt: the initial incumbent     (none)
 *     --prune   0 | 1: prune_level                         (0)
 *     --bound   1 | 2: the bound of the flowshop           (2)
 *     --branch  forward | backward | alternate | dynamic   (forward)
 *     --sort    0 | 1: the children of the flowshop sorted (0)
 *     --checkpoint <file>: segments with checkpoints       (none)
 *     --every   <evaluations of a segment>                 (100000)
 *     --cluster <processes of a cluster of hub>            (4)
 *     --unit    <evaluations of a unit of work of hub>     (32)
 *     --limit   <evaluated nodes>                          (none)
 *     --slow    <iterations of busy work per evaluation>   (0)
 *     --repeat  <runs of each instance>                    (1)
 *     --header  print the header line
 *
 * The knapsack of class 0, 1 or 2 of Pisinger (uncorrelated, weakly and
 * strongly correlated) has n items of weight at most r and the capacity
 * half the total weight, and with --values long its values are integer,
 * its bounds rounded down [see BasicKnapsack]; the traveling salesperson
 * problem has n random
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
 * computed by an expensive relaxation would. In the column algo, SelfSplit
 * with --share 1 is selfsplit-shared, and with --dive 0 selfsplit-nodive.
 * With --mode enumerate, cost is the number of the feasible nodes (the
 * nodes of a tree of UTS, the embeddings of a subgraph), and the last
 * column, recomputed, holds the evaluations that rebuild the nodes of the
 * tasks given as paths, or the ancestors of the intervals taken.
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

#include <algorithm>

#include <cstdint>

#include <cstdio>

#include <cstdlib>

#include <iostream>

#include <map>

#include <numeric>

#include <optional>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/TravelingSalesman.h"

#include "PBeFS/bnb/Flowshop.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/bnb/MaxClique.h"

#include "PBeFS/bnb/MaxCommonSubgraph.h"

#include "PBeFS/bnb/QAP.h"

#include "PBeFS/bnb/SubgraphIsomorphism.h"

#include "PBeFS/bnb/UTS.h"

#ifdef PBEFS_BENCH_MPI
 #include "PBeFS/bnb/MPIHubBranchAndBound.h"
#else
 #include "PBeFS/bnb/BBCheckpoint.h"
 #include "PBeFS/bnb/BranchAndBoundScheme.h"
 #include "PBeFS/bnb/IntervalBranchAndBound.h"
 #ifdef PBEFS_HAS_OPENMP
  #include "ADS/parallel/OpenMPExecutor.h"
 #endif
 #include "PBeFS/bnb/FarmBranchAndBound.h"
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
 { "domain" , "knapsack" } , { "size" , "" } , { "ids" , "1" } ,
 { "algo" , "sequential" } , { "runtime" , "thread" } ,
 { "select" , "best_first" } , { "threads" , "1" } , { "batch" , "16" } ,
 { "flush" , "1" } , { "round" , "0" } , { "split" , "0" } ,
 { "dive" , "1" } , { "share" , "0" } ,
 { "limit" , "0" } , { "repeat" , "1" } , { "slow" , "0" } ,
 { "values" , "double" } , { "file" , "" } , { "spawn" , "2" } ,
 { "budget" , "10000" } , { "stealall" , "0" } , { "discrepancy" , "1" } ,
 { "indexed" , "0" } , { "race" , "" } , { "racenodes" , "0" } ,
 { "mode" , "optimize" } , { "target" , "" } , { "incumbent" , "" } ,
 { "prune" , "0" } , { "bound" , "2" } , { "branch" , "forward" } ,
 { "sort" , "0" } , { "checkpoint" , "" } , { "every" , "100000" } ,
 { "cluster" , "4" } , { "unit" , "32" } };

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
  return( ADS::splitmix64( n.s ) );
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
 using value_type = bb_value_t< Q >;

 static constexpr Sense sense = Q::sense;

 const Q & q;                      // the problem
 unsigned long work;               // the iterations per evaluation

 node_type root( void ) const { return( q.root() ); }

 BasicEvaluation< value_type > evaluate( node_type & n ) const {
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

 value_type estimate( const node_type & n ) const {
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
// the options from the command line, for values of type V

static Selection selection( const std::string & sel )
{
 if( sel == "best_first" )
  return( Selection::best_first );
 if( sel == "depth_first" )
  return( Selection::depth_first );
 if( sel == "breadth_first" )
  return( Selection::breadth_first );
 if( sel == "best_estimate" )
  return( Selection::best_estimate );
 if( sel == "plunging" )
  return( Selection::plunging );
 throw( std::invalid_argument( "unknown selection " + sel ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options from the command line, for values of type V; opt is the
// optimum of the instance, if known

template< class V >
static BasicParallelBBOptions< V > bb_options( std::optional< V > best =
                                                std::nullopt )
{
 BasicParallelBBOptions< V > o;
 o.selection = selection( opt[ "select" ] );
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoul( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.round = std::stoul( opt[ "round" ] );
 o.split = std::stoul( opt[ "split" ] );
 o.split_dive = opt[ "dive" ] != "0";
 o.split_shared = opt[ "share" ] != "0";
 o.node_limit = std::stoul( opt[ "limit" ] );
 o.spawn_depth = unsigned( std::stoul( opt[ "spawn" ] ) );
 o.budget = std::stoull( opt[ "budget" ] );
 o.steal_all = opt[ "stealall" ] != "0";
 o.discrepancy = opt[ "discrepancy" ] != "0";
 o.indexed = opt[ "indexed" ] != "0";
 o.prune_level = opt[ "prune" ] != "0";
 o.race_nodes = std::stoul( opt[ "racenodes" ] );
 if( ! opt[ "race" ].empty() ) {
  std::istringstream in( opt[ "race" ] );
  for( std::string part ; std::getline( in , part , ',' ) ; ) {
   BasicBBOptions< V > r;
   r.selection = selection( part );
   o.race.push_back( r );
   }
  }
 const std::string & m = opt[ "mode" ];
 o.mode = m == "decide" ? BBMode::decide :
          m == "enumerate" ? BBMode::enumerate : BBMode::optimize;
 if( ( m != "decide" ) && ( m != "enumerate" ) && ( m != "optimize" ) )
  throw( std::invalid_argument( "unknown mode " + m ) );
 if( ! opt[ "target" ].empty() )
  o.target = V( std::stod( opt[ "target" ] ) );
 if( opt[ "incumbent" ] == "opt" ) {
  if( ! best )
   throw( std::invalid_argument( "no known optimum" ) );
  o.incumbent = *best;
  }
 else
  if( ! opt[ "incumbent" ].empty() )
   o.incumbent = V( std::stod( opt[ "incumbent" ] ) );
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
 // the variants of SelfSplit by the name of the search
 std::string algo = opt[ "algo" ];
 if( algo == "selfsplit" ) {
  if( opt[ "share" ] != "0" )
   algo += "-shared";
  if( opt[ "dive" ] == "0" )
   algo += "-nodive";
  }
 std::cout.precision( 12 );
 std::cout << opt[ "domain" ] << ',' << instance << ',' << algo
           << ',' << opt[ "runtime" ] << ',' << opt[ "select" ] << ','
           << opt[ "threads" ] << ',' << opt[ "batch" ] << ','
           << opt[ "flush" ] << ",," << status[ int( s.status ) ] << ',';
 if( opt[ "mode" ] == "enumerate" )
  std::cout << st.solutions;
 else
  std::cout << s.value;
 std::cout << ',' << st.evaluated << ',' << st.branched << ','
           << st.generated << ','
           << st.pruned_generated + st.pruned_extracted + st.pruned_evaluated
           << ',' << st.sent << ',' << st.messages << ',' << st.max_open
           << ',' << lb << ',' << co << ',' << st.seconds << ','
           << st.recomputed << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves q with the scheme, runtime and selection of the command line, the
// optimum of the instance being best, if known

template< class Q >
static BBSolution< typename Q::node_type , bb_value_t< Q > > solve(
                                     const Q & q ,
                                     std::optional< bb_value_t< Q > > best )
{
 const auto o = bb_options< bb_value_t< Q > >( best );
 const std::string & algo = opt[ "algo" ];
 #ifdef PBEFS_BENCH_MPI
  if constexpr( MPIBranchAndBound< Q > ) {
   if( algo == "hub" ) {
    BBHubParams h;
    h.cluster = unsigned( std::stoul( opt[ "cluster" ] ) );
    h.unit_work = std::stoul( opt[ "unit" ] );
    return( MPIHubBranchAndBoundSearch< Q >( q , o , h ).solve() );
    }
   return( MPIBranchAndBoundSearch< Q >( q , o ).solve() );
   }
  else
   throw( std::invalid_argument( "the nodes of this domain do not travel "
                                 "over MPI" ) );
 #else
  const std::string & r = opt[ "runtime" ];
  if( algo == "farm" ) {
   if( r == "thread" )
    return( FarmBranchAndBound< Q , ADS::ThreadPool >( q , o ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( FarmBranchAndBound< Q , ADS::OpenMPExecutor >( q , o )
            .solve() );
   #endif
   #ifdef PBEFS_HAS_FASTFLOW
   if( r == "fastflow" )
    return( FarmBranchAndBound< Q , ADS::FastFlowExecutor >( q , o )
            .solve() );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   }
  if( algo == "ivm" ) {
   if constexpr( RadixBranchAndBound< Q > ) {
    if( r == "thread" )
     return( IntervalBranchAndBound< Q , ADS::ThreadPool >( q , o )
             .solve() );
    #ifdef PBEFS_HAS_OPENMP
    if( r == "openmp" )
     return( IntervalBranchAndBound< Q , ADS::OpenMPExecutor >( q , o )
             .solve() );
    #endif
    #ifdef PBEFS_HAS_FASTFLOW
    if( r == "fastflow" )
     return( IntervalBranchAndBound< Q , ADS::FastFlowExecutor >( q , o )
             .solve() );
    #endif
    throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
    }
   else
    throw( std::invalid_argument( "ivm is for flowshop and qap" ) );
   }
  const BBScheme s = bb_scheme( algo );
  const auto by = [ & ]( const Q & x , const auto & so , auto && starts ) {
   if( r == "thread" )
    return( parallel_branch_and_bound< Q , ADS::ThreadPool >( x , s , so ,
             std::move( starts ) ) );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( parallel_branch_and_bound< Q , ADS::OpenMPExecutor >( x , s ,
             so , std::move( starts ) ) );
   #endif
   #ifdef PBEFS_HAS_FASTFLOW
   if( r == "fastflow" )
    return( parallel_branch_and_bound< Q , ADS::FastFlowExecutor >( x , s ,
             so , std::move( starts ) ) );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   };
  if( ! opt[ "checkpoint" ].empty() ) {
   if constexpr( SerializableBranchAndBound< Q > )
    return( checkpointed_branch_and_bound( q , o , by ,
                                           opt[ "checkpoint" ] ,
                                           std::stoul( opt[ "every" ] ) ) );
   else
    throw( std::invalid_argument( "no checkpoint for these nodes" ) );
   }
  return( by( q , o , std::optional< std::vector< BBLeaf<
               typename Q::node_type , bb_value_t< Q > > > >() ) );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs q --repeat times, the line labeled by instance, the optimum of
// the instance being best, if known

template< class Q >
static void run( const std::string & instance , const Q & q ,
                 std::optional< bb_value_t< Q > > best = std::nullopt )
{
 const int repeat = std::stoi( opt[ "repeat" ] );
 const unsigned long work = std::stoul( opt[ "slow" ] );
 for( int r = 0 ; r < repeat ; ++r )
  if( work == 0 )
   print( instance , solve( q , best ) );
  else
   if constexpr( requires( const Q & x , typename Q::node_type & n ) {
                  x.branch( n , detail::bb_sink< typename Q::node_type ,
                                                 bb_value_t< Q > >{} );
                  x.hash( n ); x.estimate( n ); } )
    print( instance , solve( SlowProblem< Q >{ q , work } , best ) );
   else
    throw( std::invalid_argument( "--slow is for knapsack and tsp" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack of the fields f of --size and of the seed, with values of
// type V

template< class V >
static BasicKnapsack< V > random_knapsack(
                              const std::vector< std::string > & f , int seed )
{
 return( BasicKnapsack< V >::random( std::stoi( f[ 0 ] ) ,
                                     std::stoul( f[ 1 ] ) ,
                                     std::stol( f[ 2 ] ) ,
                                     std::uint64_t( seed ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the tree of UTS of the fields f of --size: t1, t2, t3, or
// geo:<shape>:<b>:<g>:<seed> (the shape 0 linear, 1 exponential, 2 cyclic,
// 3 fixed) or bin:<b>:<m>:<q>:<seed>

static UTS uts_tree( const std::vector< std::string > & f )
{
 if( f.size() == 1 ) {
  if( f[ 0 ] == "t1" )
   return( UTS::t1() );
  if( f[ 0 ] == "t2" )
   return( UTS::t2() );
  if( f[ 0 ] == "t3" )
   return( UTS::t3() );
  }
 if( ( f.size() == 5 ) && ( f[ 0 ] == "geo" ) )
  return( UTS::geometric( UTS::Shape( std::stoi( f[ 1 ] ) ) ,
                          std::stod( f[ 2 ] ) , std::stoi( f[ 3 ] ) ,
                          std::stoi( f[ 4 ] ) ) );
 if( ( f.size() == 5 ) && ( f[ 0 ] == "bin" ) )
  return( UTS::binomial_tree( std::stod( f[ 1 ] ) , std::stoi( f[ 2 ] ) ,
                              std::stod( f[ 3 ] ) , std::stoi( f[ 4 ] ) ) );
 throw( std::invalid_argument( "--size of uts: t1, t2, t3, geo:... or "
                               "bin:..." ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a random graph of n vertices, each edge there with the given percent

static BitGraph random_graph( std::size_t n , std::size_t pct ,
                              std::uint64_t seed )
{
 std::mt19937_64 gen( seed );
 BitGraph g( n );
 for( std::size_t i = 0 ; i < n ; ++i )
  for( std::size_t j = i + 1 ; j < n ; ++j )
   if( random_below( gen , 100 ) < pct )
    g.add_edge( i , j );
 return( g );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 #ifdef PBEFS_BENCH_MPI
  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  opt[ "algo" ] = "mpi";        // or hub
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
                "messages,stored,LB,CO,seconds,recomputed" << std::endl;

  const auto f = size_fields();
  const std::string & d = opt[ "domain" ];
  for( int seed : parse_ids( opt[ "ids" ] ) ) {
   const std::string instance = opt[ "size" ] + "/" + std::to_string( seed );
   if( d == "knapsack" ) {
    if( f.size() != 3 )
     throw( std::invalid_argument( "--size wants <class>:<n>:<r>" ) );
    if( opt[ "values" ] == "double" )
     run( instance , random_knapsack< double >( f , seed ) );
    else
     if( opt[ "values" ] == "long" )
      run( instance , random_knapsack< long >( f , seed ) );
     else
      throw( std::invalid_argument( "unknown values " + opt[ "values" ] ) );
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
     if( d == "uts" )
      run( opt[ "size" ] , uts_tree( f ) );
     else
      if( d == "clique" )
       run( opt[ "file" ] , MaxClique( BitGraph::dimacs( opt[ "file" ] ) ) );
      else
       if( d == "qap" )
        run( opt[ "file" ] , QAP::qaplib( opt[ "file" ] ) );
       else
        if( d == "flowshop" ) {
         Flowshop fs = Flowshop::taillard( seed );
         fs.bound = opt[ "bound" ] == "1" ? Flowshop::Bound::one_machine
                                          : Flowshop::Bound::two_machine;
         const std::string & b = opt[ "branch" ];
         fs.branching = b == "backward" ? Flowshop::Branching::backward :
                        b == "alternate" ? Flowshop::Branching::alternate :
                        b == "dynamic" ? Flowshop::Branching::dynamic :
                        Flowshop::Branching::forward;
         fs.sort_children = opt[ "sort" ] != "0";
         const int best = Flowshop::taillard_optimum( seed );
         char name[ 8 ];
         std::snprintf( name , sizeof( name ) , "ta%03d" , seed );
         run( name , fs , best > 0 ? std::optional< int >( best )
                                   : std::nullopt );
         }
        else
         if( ( d == "sip" ) || ( d == "mcs" ) ) {
          if( f.size() != 4 )
           throw( std::invalid_argument(
                   "--size wants <np>:<pp>:<nt>:<pt>" ) );
          const BitGraph a = random_graph( std::stoul( f[ 0 ] ) ,
                                           std::stoul( f[ 1 ] ) , 2 * seed );
          const BitGraph b = random_graph( std::stoul( f[ 2 ] ) ,
                                           std::stoul( f[ 3 ] ) ,
                                           2 * seed + 1 );
          if( d == "sip" )
           run( instance , SubgraphIsomorphism{ a , b , false } );
          else
           run( instance , MaxCommonSubgraph{ a , b } );
          }
         else
          throw( std::invalid_argument( "unknown domain " + d ) );
   if( ( d == "uts" ) || ( d == "clique" ) || ( d == "qap" ) )
    break;                              // a single instance
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
