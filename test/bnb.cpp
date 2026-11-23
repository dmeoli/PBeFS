/*--------------------------------------------------------------------------*/
/*------------------------------ File bnb.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the skeletons of the
 * Branch-and-Bound (B&B), its modes, the racing ramp-up, the checkpoints,
 * the work stealing of intervals, the master-hub-worker search over MPI,
 * an evaluator that throws and the example domains. The references do not
 * come from the library: the dynamic program of the knapsack, a plain
 * recursion over the trees of UTS, the enumeration of all the permutations
 * of the small flowshops and quadratic assignments, the enumeration of all
 * the injections of the small subgraph isomorphisms, the algorithm of Bron
 * and Kerbosch for the small cliques, and the published optima of the
 * instances of QAPLIB, of Taillard and of the DIMACS graph brock200_1.
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

#include <atomic>

#include <chrono>

#include <condition_variable>

#include <cstdio>

#include <cstdlib>

#include <filesystem>

#include <functional>

#include <iostream>

#include <mutex>

#include <numeric>

#include <random>

#include <stdexcept>

#include <string>

#include <thread>

#include <vector>

#include "PBeFS/bnb/BBCheckpoint.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

#include "PBeFS/bnb/Flowshop.h"

#include "PBeFS/bnb/IntervalBranchAndBound.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/bnb/MaxClique.h"

#include "PBeFS/bnb/MaxCommonSubgraph.h"

#include "PBeFS/bnb/QAP.h"

#include "PBeFS/bnb/SkeletonBranchAndBound.h"

#include "PBeFS/bnb/SubgraphIsomorphism.h"

#include "PBeFS/bnb/UTS.h"

#include "test.h"

#ifdef PBEFS_HAS_OPENMP
 #include "ADS/parallel/OpenMPExecutor.h"
#endif

#ifdef PBEFS_HAS_FASTFLOW
 #include "ADS/parallel/FastFlowExecutor.h"
#endif

#ifdef PBEFS_HAS_MPI
 #include "PBeFS/bnb/MPIHubBranchAndBound.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/*--------------------------------------------------------------------------*/
// the name of an executor

template< class Executor >
std::string ex_name( void )
{
 if constexpr( std::is_same_v< Executor , ADS::ThreadPool > )
  return( "thread" );
#ifdef PBEFS_HAS_OPENMP
 else
  if constexpr( std::is_same_v< Executor , ADS::OpenMPExecutor > )
   return( "openmp" );
#endif
 else
  return( "fastflow" );
 }

/*--------------------------------------------------------------------------*/
// random knapsacks, small enough for every scheme

std::vector< Knapsack > small_knapsacks( std::size_t n = 30 )
{
 std::vector< Knapsack > v;
 for( int cls = 0 ; cls < 3 ; ++cls )
  for( std::uint64_t s = 1 ; s <= 2 ; ++s )
   v.push_back( Knapsack::random( cls , n , 1000 , 100 * cls + s ) );
 return( v );
 }

/*--------------------------------------------------------------------------*/
// the names of the rules of the skeletons

const char * rule_name( BBRule r )
{
 static const char * names[] = { "depthbounded" , "budget" ,
                                 "stackstealing" , "ordered" };
 return( names[ int( r ) ] );
 }

/*--------------------------------------------------------------------------*/
// the number of the feasible nodes of the tree of a knapsack, by a plain
// recursion on its evaluate() and branch()

template< class K >
long count_tree( const K & kp , typename K::node_type n )
{
 const auto e = kp.evaluate( n );
 if( e.infeasible )
  return( 0 );
 long c = e.feasible ? 1 : 0;
 std::vector< typename K::node_type > kids;
 kp.branch( n , [ & ]( auto && x , auto ) { kids.push_back( x ); } );
 for( auto & k : kids )
  c += count_tree( kp , k );
 return( c );
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*------------------------------- THE CASES --------------------------------*/
/*--------------------------------------------------------------------------*/
// the four skeletons, with nodes and with paths as tasks, on the knapsack
// (the optimum, the node of the incumbent, a node limit and the restart
// from its open leaves, the leaves of a whole search) and on UTS (the
// count of its nodes), with 1 to 8 workers

template< class Executor >
void test_bb_skeleton( const std::string & )
{
 const std::string ex = ex_name< Executor >();
 const UTS u = UTS::geometric( UTS::Shape::fixed , 4 , 6 , 19 );
 const long nodes = u.count( u.root() );
 for( const BBRule r : { BBRule::depth_bounded , BBRule::budget ,
                         BBRule::stack_stealing , BBRule::ordered } )
  for( const bool ix : { false , true } )
   for( const unsigned p : { 1u , 2u , 3u , 8u } ) {
    const std::string tag = std::string( rule_name( r ) ) +
                            ( ix ? " by paths" : "" ) + " on " + ex +
                            " with " + std::to_string( p ) + " workers";
    {
     BasicParallelBBOptions< long > o;
     o.threads = p;
     o.mode = BBMode::enumerate;
     o.indexed = ix;
     o.budget = 20;
     o.steal_all = ( p == 3 );
     const auto s = SkeletonBranchAndBound< UTS , Executor >( u , r ,
                                                             o ).solve();
     check( ( long( s.stats.solutions ) == nodes ) &&
            ( long( s.stats.evaluated ) == nodes ) ,
            "UTS: " + std::to_string( s.stats.solutions ) + " nodes of " +
            std::to_string( nodes ) + ", " + tag );
     }
    for( const Knapsack & kp : small_knapsacks() ) {
     const long opt = kp.optimum();
     ParallelBBOptions o;
     o.threads = p;
     o.indexed = ix;
     o.budget = 10;
     o.spawn_depth = 3;
     o.discrepancy = ( p != 2 );
     o.steal_all = ( p == 8 );
     const auto s = SkeletonBranchAndBound< Knapsack , Executor >( kp , r ,
                                                                  o ).solve();
     check( ( s.status == BBStatus::optimal ) && ( long( s.value ) == opt ) &&
            s.node && kp.check( *s.node , s.value ) ,
            "knapsack: optimum " + std::to_string( opt ) + ", " + tag );
     ParallelBBOptions l = o;
     l.node_limit = 15;
     l.keep_frontier = true;
     const auto a = SkeletonBranchAndBound< Knapsack , Executor >( kp , r ,
                                                                  l ).solve();
     ParallelBBOptions m = o;
     if( a.feasible() )
      m.incumbent = a.value;
     std::vector< BBLeaf< Knapsack::node_type > > open;
     for( const auto & x : a.frontier )
      if( x.kind == BBLeafKind::open )
       open.push_back( x );
     const auto b = SkeletonBranchAndBound< Knapsack , Executor >( kp , r ,
                                                         m ).solve( open );
     const double v = std::max( a.feasible() ? a.value : 0. ,
                                b.feasible() ? b.value : 0. );
     check( ( a.status != BBStatus::node_limit ) || ( long( v ) == opt ) ,
            "knapsack: the restart from the open leaves, " + tag );
     check( a.bound >= opt - 1e-6 ,
            "knapsack: the bound of a stopped search, " + tag );
     }
    }

 // prune_level on the cliques, whose children come by their colors
 std::mt19937_64 gen( 7 );
 BitGraph g( 60 );
 for( std::size_t i = 0 ; i < 60 ; ++i )
  for( std::size_t j = i + 1 ; j < 60 ; ++j )
   if( random_below( gen , 100 ) < 50 )
    g.add_edge( i , j );
 const MaxClique mc( g );
 int best = -1;
 std::uint64_t gen0 = 0;
 for( const bool pl : { false , true } )
  for( const BBRule r : { BBRule::stack_stealing , BBRule::budget } )
   for( const unsigned p : { 1u , 4u } ) {
    BasicParallelBBOptions< int > o;
    o.threads = p;
    o.prune_level = pl;
    const auto s = SkeletonBranchAndBound< MaxClique , Executor >( mc , r ,
                                                                  o ).solve();
    if( best < 0 ) {
     best = s.value;
     gen0 = s.stats.generated;
     }
    check( ( s.value == best ) && mc.check( *s.node ) &&
           ( int( s.node->clique.size() ) == s.value ) ,
           "clique of a random graph by " + std::string( rule_name( r ) ) +
           ( pl ? " with prune_level" : "" ) + " on " + ex );
    if( pl && ( p == 1 ) && ( r == BBRule::stack_stealing ) )
     check( s.stats.generated < gen0 ,
            "prune_level builds fewer children (" +
            std::to_string( s.stats.generated ) + " < " +
            std::to_string( gen0 ) + ")" );
    }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the modes decide and enumerate with every scheme

void test_bb_modes( const std::string & )
{
 for( const Knapsack & kp : small_knapsacks( 14 ) ) {
  const long opt = kp.optimum();
  const long tree = count_tree( kp , kp.root() );
  for( std::size_t s = 0 ; s < 11 ; ++s )
   for( const unsigned p : { 1u , 3u } ) {
    const BBScheme sc = BBScheme( s );
    const std::string tag = std::string( bb_scheme_name( sc ) ) + " with " +
                            std::to_string( p ) + " workers";
    ParallelBBOptions o;
    o.threads = p;
    o.mode = BBMode::decide;
    o.target = double( opt - 50 );
    auto d = parallel_branch_and_bound( kp , sc , o );
    check( ( d.status == BBStatus::optimal ) && ( d.value >= opt - 50 ) &&
           d.node && kp.check( *d.node , d.value ) ,
           "decide: a value within the target, " + tag );
    o.target = double( opt + 1 );
    d = parallel_branch_and_bound( kp , sc , o );
    check( d.status == BBStatus::infeasible ,
           "decide: no value beyond the optimum, " + tag );
    o.target.reset();
    d = parallel_branch_and_bound( kp , sc , o );
    check( ( d.status == BBStatus::optimal ) && d.feasible() ,
           "decide: any feasible value, " + tag );
    o.mode = BBMode::enumerate;
    const auto e = parallel_branch_and_bound( kp , sc , o );
    check( ( long( e.stats.solutions ) == tree ) &&
           ( long( e.value ) == opt ) ,
           "enumerate: " + std::to_string( e.stats.solutions ) + " of " +
           std::to_string( tree ) + " feasible nodes, " + tag );
    }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the racing ramp-up before every scheme

void test_bb_race( const std::string & )
{
 for( const Knapsack & kp : small_knapsacks() ) {
  const long opt = kp.optimum();
  for( std::size_t s = 0 ; s < 11 ; ++s )
   for( const std::size_t budget : { std::size_t( 5 ) , std::size_t( 50 ) ,
                                     std::size_t( 100000 ) } ) {
    const BBScheme sc = BBScheme( s );
    ParallelBBOptions o;
    o.threads = 3;
    o.race.resize( 3 );
    o.race[ 1 ].selection = Selection::depth_first;
    o.race[ 2 ].selection = Selection::plunging;
    o.race_nodes = budget;
    const auto r = parallel_branch_and_bound( kp , sc , o );
    check( ( r.status == BBStatus::optimal ) && ( long( r.value ) == opt ) &&
           r.node && kp.check( *r.node , r.value ) &&
           ( r.bound >= opt - 1e-6 ) && ( r.bound <= opt + 1e-6 ) ,
           "race of " + std::to_string( budget ) + " nodes, then " +
           std::string( bb_scheme_name( sc ) ) );
    }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the checkpoints: a search in segments gives the optimum; one stopped by
// a limit is resumed from its file; a truncated or corrupted file throws

void test_bb_checkpoint( const std::string & )
{
 // a name of its own, since other runs of the tester may share the folder
 const std::string file = ( std::filesystem::temp_directory_path() /
                            ( "pbefs-checkpoint-" + std::to_string(
                               std::random_device()() ) + ".bin" ) ).string();
 for( const Knapsack & kp : small_knapsacks() ) {
  const long opt = kp.optimum();
  for( const BBScheme sc : { BBScheme::sequential , BBScheme::shared ,
                             BBScheme::stack_stealing ,
                             BBScheme::deterministic } ) {
   const auto search = [ sc ]( const Knapsack & q , const auto & so ,
                               auto && starts ) {
    return( parallel_branch_and_bound( q , sc , so ,
                                       std::move( starts ) ) );
    };
   const std::string tag = bb_scheme_name( sc ).data();
   std::remove( file.c_str() );
   ParallelBBOptions o;
   o.threads = 3;
   const auto s = checkpointed_branch_and_bound( kp , o , search , file ,
                                                 7 );
   check( ( s.status == BBStatus::optimal ) && ( long( s.value ) == opt ) &&
          s.node && kp.check( *s.node , s.value ) &&
          ( s.bound <= opt + 1e-6 ) && ( s.bound >= opt - 1e-6 ) ,
          "segments of 7 evaluations, " + tag );
   const auto c = bb_read_checkpoint( kp , file );
   check( c && c->open.empty() , "the last checkpoint is empty, " + tag );

   // stopped after 20 evaluations, then resumed
   std::remove( file.c_str() );
   o.node_limit = 20;
   const auto a = checkpointed_branch_and_bound( kp , o , search , file ,
                                                 6 );
   const auto mid = bb_read_checkpoint( kp , file );
   check( mid && ( ( a.status != BBStatus::node_limit ) ||
                   ! mid->open.empty() ) &&
          ( a.bound >= opt - 1e-6 ) ,
          "the checkpoint of a run stopped by the limit, " + tag );
   o.node_limit = 0;
   const auto b = checkpointed_branch_and_bound( kp , o , search , file ,
                                                 6 );
   check( ( b.status == BBStatus::optimal ) && ( long( b.value ) == opt ) ,
          "the run resumed from the checkpoint, " + tag );
   check( b.stats.evaluated >= a.stats.evaluated ,
          "the counters go on from the checkpoint, " + tag );
   }
  }

 // truncated and corrupted files: bb_read_checkpoint() throws, and reads
 // nothing past their end (which the sanitizers would see)
 const auto ks = small_knapsacks();
 const Knapsack & kp = ks.front();
 using V = bb_value_t< Knapsack >;
 const auto bytes = [ & ]( void ) {
  std::ifstream in( file , std::ios::binary );
  return( std::string( ( std::istreambuf_iterator< char >( in ) ) ,
                       std::istreambuf_iterator< char >() ) );
  };
 const auto put = [ & ]( const std::string & b ) {
  std::ofstream( file , std::ios::binary | std::ios::trunc ) << b;
  };
 const auto throws = [ & ]( void ) {
  try {
   static_cast< void >( bb_read_checkpoint( kp , file ) );
   }
  catch( const std::runtime_error & ) {
   return( true );
   }
  return( false );
  };
 BBCheckpoint< Knapsack::node_type , V > c;
 for( unsigned d = 0 ; d < 3 ; ++d )
  c.open.push_back( { kp.root() , V( 10 + d ) , d , BBLeafKind::open } );
 bb_write_checkpoint( kp , c , file );
 const std::string whole = bytes();
 const auto back = bb_read_checkpoint( kp , file );
 check( back && ( back->open.size() == 3 ) && ! back->incumbent &&
        ! back->node , "a checkpoint of 3 open nodes read back" );
 c.incumbent = V( 5 );
 c.node = kp.root();
 bb_write_checkpoint( kp , c , file );
 const std::string with_node = bytes();
 int missed = 0;
 for( const std::string * w : { & whole , & with_node } )
  for( std::size_t l = 0 ; l < w->size() ; ++l ) {
   put( w->substr( 0 , l ) );
   missed += ! throws();
   }
 check( missed == 0 , "every truncation of two checkpoints of " +
        std::to_string( whole.size() ) + " and " +
        std::to_string( with_node.size() ) + " bytes throws" );

 // the size of the values, the flag of the incumbent and the number of
 // open nodes (after the magic, the size, the flag, the incumbent, the flag
 // of its node, the closed bound and 12 counters) changed
 const std::size_t at_n = 8 + 4 + 1 + sizeof( V ) + 1 + sizeof( V ) + 96;
 const auto with = [ & ]( std::size_t at , const auto x ) {
  std::string b = whole;
  std::memcpy( b.data() + at , & x , sizeof( x ) );
  put( b );
  return( throws() );
  };
 const bool all = with( 8 , std::uint32_t( sizeof( V ) + 1 ) ) &&
                  with( 12 , std::uint8_t( 7 ) ) &&
                  with( 13 + sizeof( V ) , std::uint8_t( 1 ) ) &&
                  with( at_n , std::uint64_t( 1 ) << 61 ) &&
                  with( at_n , ~ std::uint64_t( 0 ) ) &&
                  with( at_n , std::uint64_t( 4 ) ) &&
                  with( at_n , std::uint64_t( 2 ) );
 check( all , "a corrupted size of the values, flag or number of open nodes "
        "throws" );
 std::remove( file.c_str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the work stealing of intervals on small flowshops and quadratic
// assignments, against the enumeration of all the permutations

namespace {

Flowshop random_flowshop( int n , int m , std::uint64_t seed )
{
 std::mt19937_64 gen( seed );
 Flowshop f;
 f.n = n;
 f.m = m;
 f.p.resize( std::size_t( n * m ) );
 for( auto & x : f.p )
  x = 1 + int( random_below( gen , 99 ) );
 f.prepare();
 return( f );
 }

int brute_flowshop( const Flowshop & f )
{
 std::vector< std::uint16_t > s( std::size_t( f.n ) );
 std::iota( s.begin() , s.end() , std::uint16_t( 0 ) );
 int best = INT_MAX;
 do
  best = std::min( best , f.makespan( s ) );
 while( std::next_permutation( s.begin() , s.end() ) );
 return( best );
 }

QAP random_qap( std::size_t n , std::uint64_t seed , bool asym )
{
 std::mt19937_64 gen( seed );
 QAP q;
 q.n = n;
 q.f.assign( n * n , 0 );
 q.dist.assign( n * n , 0 );
 for( std::size_t i = 0 ; i < n ; ++i )
  for( std::size_t j = 0 ; j < n ; ++j )
   if( i != j ) {
    q.f[ i * n + j ] = long( random_below( gen , 10 ) );
    q.dist[ i * n + j ] = long( 1 + random_below( gen , 20 ) );
    if( ! asym && ( j < i ) ) {
     q.f[ i * n + j ] = q.f[ j * n + i ];
     q.dist[ i * n + j ] = q.dist[ j * n + i ];
     }
    }
 q.prepare();
 return( q );
 }

long brute_qap( const QAP & q )
{
 std::vector< std::size_t > s( q.n );
 std::iota( s.begin() , s.end() , std::size_t( 0 ) );
 long best = std::numeric_limits< long >::max();
 do {
  long c = 0;
  for( std::size_t i = 0 ; i < q.n ; ++i )
   for( std::size_t j = 0 ; j < q.n ; ++j )
    c += q.flow( i , j ) * q.distance( s[ i ] , s[ j ] );
  best = std::min( best , c );
  } while( std::next_permutation( s.begin() , s.end() ) );
 return( best );
 }

}  // end( namespace )

template< class Executor >
void test_bb_interval( const std::string & )
{
 const std::string ex = ex_name< Executor >();
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  Flowshop f = random_flowshop( 7 , 3 + int( seed % 2 ) , seed );
  f.branching = seed % 2 ? Flowshop::Branching::forward
                         : Flowshop::Branching::alternate;
  f.bound = seed < 3 ? Flowshop::Bound::one_machine
                     : Flowshop::Bound::two_machine;
  f.sort_children = ( seed == 4 );
  const int opt = brute_flowshop( f );
  std::uint64_t ev1 = 0;
  for( const unsigned p : { 1u , 2u , 3u , 8u } ) {
   BasicParallelBBOptions< int > o;
   o.threads = p;
   const auto s = IntervalBranchAndBound< Flowshop , Executor >( f ,
                                                                o ).solve();
   check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
          s.node && ( f.makespan( f.sequence( *s.node ) ) == opt ) ,
          "intervals on a flowshop 7x" + std::to_string( f.m ) + ", " + ex +
          " with " + std::to_string( p ) + " workers" );
   // with the optimum given, the nodes evaluated are those of the
   // sequential depth-first search, whatever the steals
   o.incumbent = opt;
   const auto t = IntervalBranchAndBound< Flowshop , Executor >( f ,
                                                                o ).solve();
   if( p == 1 )
    ev1 = t.stats.evaluated;
   check( t.stats.evaluated == ev1 ,
          "intervals: the same nodes evaluated with the optimum given (" +
          std::to_string( t.stats.evaluated ) + " and " +
          std::to_string( ev1 ) + "), " + std::to_string( p ) + " workers" );
   // a node limit: the leaves left cover the rest of the tree
   o.incumbent.reset();
   o.node_limit = 30;
   o.keep_frontier = true;
   o.keep_open_only = true;
   const auto a = IntervalBranchAndBound< Flowshop , Executor >( f ,
                                                                o ).solve();
   if( a.status == BBStatus::node_limit ) {
    BasicParallelBBOptions< int > r;
    r.threads = p;
    if( a.feasible() )
     r.incumbent = a.value;
    const auto b = SkeletonBranchAndBound< Flowshop , Executor >( f ,
                    BBRule::stack_stealing , r ).solve( a.frontier );
    const int v = std::min( a.feasible() ? a.value : INT_MAX ,
                            b.feasible() ? b.value : INT_MAX );
    check( v == opt , "intervals: the restart from the leaves left by a "
                      "limit, " + std::to_string( p ) + " workers" );
    }
   }
  }
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const QAP q = random_qap( 7 , seed , seed == 3 );
  const long opt = brute_qap( q );
  for( const unsigned p : { 1u , 4u } ) {
   BasicParallelBBOptions< long > o;
   o.threads = p;
   const auto s = IntervalBranchAndBound< QAP , Executor >( q , o ).solve();
   check( ( s.value == opt ) && s.node && ( q.cost( *s.node ) == opt ) ,
          "intervals on a quadratic assignment of 7" +
          std::string( seed == 3 ? " (asymmetric)" : "" ) + ", " + ex +
          " with " + std::to_string( p ) + " workers" );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the example domains against their references

namespace {

// the largest clique of g, by Bron and Kerbosch with a pivot
std::size_t bron_kerbosch( const BitGraph & g )
{
 std::size_t best = 0;
 std::function< void( std::size_t , BitSet , BitSet ) > rec =
  [ & ]( std::size_t r , BitSet p , BitSet x ) {
   if( p.empty() && x.empty() ) {
    best = std::max( best , r );
    return;
    }
   BitSet u = p;
   for( std::size_t j = 0 ; j < u.w.size() ; ++j )
    u.w[ j ] |= x.w[ j ];
   const std::size_t piv = u.first();
   BitSet c = p;
   c.remove( g.adj[ piv ] );
   while( ! c.empty() ) {
    const std::size_t v = c.first();
    c.reset( v );
    BitSet p2 = p , x2 = x;
    p2 &= g.adj[ v ];
    x2 &= g.adj[ v ];
    rec( r + 1 , p2 , x2 );
    p.reset( v );
    x.set( v );
    }
   };
 BitSet all( g.n );
 for( std::size_t v = 0 ; v < g.n ; ++v )
  all.set( v );
 rec( 0 , all , BitSet( g.n ) );
 return( best );
 }

BitGraph random_graph( std::size_t n , unsigned pct , std::uint64_t seed )
{
 std::mt19937_64 gen( seed );
 BitGraph g( n );
 for( std::size_t i = 0 ; i < n ; ++i )
  for( std::size_t j = i + 1 ; j < n ; ++j )
   if( random_below( gen , 100 ) < pct )
    g.add_edge( i , j );
 return( g );
 }

// the number of the injections of the pattern into the target that map
// edges to edges (and non-edges to non-edges if induced)
long count_embeddings( const BitGraph & pa , const BitGraph & ta ,
                       bool induced )
{
 std::vector< std::size_t > f( pa.n );
 std::vector< char > used( ta.n , 0 );
 std::function< long( std::size_t ) > rec = [ & ]( std::size_t k ) -> long {
  if( k == pa.n )
   return( 1 );
  long c = 0;
  for( std::size_t t = 0 ; t < ta.n ; ++t ) {
   if( used[ t ] )
    continue;
   bool ok = true;
   for( std::size_t q = 0 ; ok && ( q < k ) ; ++q ) {
    if( pa.edge( k , q ) && ! ta.edge( t , f[ q ] ) )
     ok = false;
    if( induced && ! pa.edge( k , q ) && ta.edge( t , f[ q ] ) )
     ok = false;
    }
   if( ! ok )
    continue;
   used[ t ] = 1;
   f[ k ] = t;
   c += rec( k + 1 );
   used[ t ] = 0;
   }
  return( c );
  };
 return( rec( 0 ) );
 }

}  // end( namespace )

void test_bb_domains( const std::string & data )
{
 // the cliques: random graphs against Bron and Kerbosch, brock200_1; in a
 // light run [see light_run()] the densest graph and the search of
 // brock200_1 are left out, the file is only read
 for( std::uint64_t s = 1 ; s <= ( light_run() ? 3 : 4 ) ; ++s ) {
  const BitGraph g = random_graph( 40 + 10 * s , 30 + 15 * unsigned( s ) ,
                                   s );
  const MaxClique mc( g );
  const auto r = parallel_branch_and_bound( mc , BBScheme::shared ,
                                            BasicParallelBBOptions< int >() );
  check( ( std::size_t( r.value ) == bron_kerbosch( g ) ) &&
         mc.check( *r.node ) , "max clique of a random graph " +
                               std::to_string( s ) );
  }
 if( light_run() ) {
  const BitGraph g = BitGraph::dimacs( data + "/dimacs/brock200_1.clq" );
  std::size_t deg = 0;
  for( std::size_t v = 0 ; v < g.n ; ++v )
   deg += g.degree( v );
  check( ( g.n == 200 ) && ( deg == 2 * 14834 ) ,
         "brock200_1 has 200 vertices and 14834 edges" );
  }
 else {
  const MaxClique mc( BitGraph::dimacs( data + "/dimacs/brock200_1.clq" ) );
  BasicParallelBBOptions< int > o;
  o.threads = 4;
  o.prune_level = true;
  const auto r = SkeletonBranchAndBound< MaxClique >( mc ,
                  BBRule::stack_stealing , o ).solve();
  check( ( r.value == 21 ) && mc.check( *r.node ) ,
         "brock200_1 has a clique of 21 vertices (" +
         std::to_string( r.value ) + ")" );
  }

 // the common subgraphs: against the cliques of the modular product
 for( std::uint64_t s = 1 ; s <= 4 ; ++s ) {
  MaxCommonSubgraph mcs;
  mcs.g = random_graph( 7 + s % 2 , 40 , s );
  mcs.h = random_graph( 8 , 50 , 10 + s );
  const MaxClique mp( mcs.modular_product() );
  const auto c = branch_and_bound( mp );
  for( const BBScheme sc : { BBScheme::sequential ,
                             BBScheme::stack_stealing ,
                             BBScheme::distributed } ) {
   BasicParallelBBOptions< int > o;
   o.threads = 3;
   const auto r = parallel_branch_and_bound( mcs , sc , o );
   check( ( r.value == c.value ) && mcs.check( *r.node ) &&
          ( int( r.node->match.size() ) == r.value ) ,
          "common subgraph " + std::to_string( s ) + " of " +
          std::to_string( r.value ) + " vertices, as the modular product, " +
          std::string( bb_scheme_name( sc ) ) );
   }
  }

 // the subgraph isomorphisms: decide and count against the enumeration
 for( std::uint64_t s = 1 ; s <= 6 ; ++s ) {
  SubgraphIsomorphism si;
  si.target = random_graph( 9 , 45 , 20 + s );
  si.pattern = random_graph( 4 + s % 2 , 50 , 40 + s );
  si.induced = ( s % 3 == 0 );
  const long cnt = count_embeddings( si.pattern , si.target , si.induced );
  for( const BBScheme sc : { BBScheme::sequential , BBScheme::shared ,
                             BBScheme::depth_bounded ,
                             BBScheme::stack_stealing } ) {
   BasicParallelBBOptions< int > o;
   o.threads = 3;
   o.mode = BBMode::enumerate;
   const auto e = parallel_branch_and_bound( si , sc , o );
   check( long( e.stats.solutions ) == cnt ,
          "subgraph isomorphism " + std::to_string( s ) + ": " +
          std::to_string( e.stats.solutions ) + " of " +
          std::to_string( cnt ) + " embeddings, " +
          std::string( bb_scheme_name( sc ) ) );
   o.mode = BBMode::decide;
   const auto d = parallel_branch_and_bound( si , sc , o );
   check( ( d.status == ( cnt ? BBStatus::optimal : BBStatus::infeasible ) ) &&
          ( ! cnt || ( d.node && si.check( *d.node ) ) ) ,
          "subgraph isomorphism " + std::to_string( s ) + " decided, " +
          std::string( bb_scheme_name( sc ) ) );
   }
  }

 // QAPLIB, with their published optima; nug12, of 34116 nodes, not in a
 // light run
 const std::pair< const char * , long > qaps[] = {
  { "chr12a" , 9552 } , { "had12" , 1652 } , { "nug12" , 578 } };
 for( const auto & [ name , opt ] : qaps ) {
  if( light_run() && ( std::string( name ) == "nug12" ) )
   continue;
  const QAP q = QAP::qaplib( data + "/qaplib/" + name + ".dat" );
  BasicParallelBBOptions< long > o;
  o.threads = 4;
  o.incumbent = opt + 1;
  const auto r = SkeletonBranchAndBound< QAP >( q , BBRule::stack_stealing ,
                                                o ).solve();
  check( ( r.value == opt ) && r.node && ( q.cost( *r.node ) == opt ) ,
         std::string( name ) + " has the optimum " + std::to_string( opt ) +
         " (" + std::to_string( r.value ) + ", " +
         std::to_string( r.stats.evaluated ) + " nodes)" );
  }

 // the flowshops: small ones against the enumeration, every bound and
 // branching, and ta001 with the two-machine bound
 for( std::uint64_t s = 1 ; s <= 4 ; ++s ) {
  Flowshop f = random_flowshop( 7 , 4 , 50 + s );
  const int opt = brute_flowshop( f );
  for( const auto bd : { Flowshop::Bound::one_machine ,
                         Flowshop::Bound::two_machine } )
   for( const auto br : { Flowshop::Branching::forward ,
                          Flowshop::Branching::backward ,
                          Flowshop::Branching::alternate ,
                          Flowshop::Branching::dynamic } ) {
    f.bound = bd;
    f.branching = br;
    const auto r = parallel_branch_and_bound( f , BBScheme::shared ,
                                        BasicParallelBBOptions< int >() );
    check( ( r.value == opt ) && r.node &&
           ( f.makespan( f.sequence( *r.node ) ) == opt ) ,
           "flowshop " + std::to_string( s ) + ", bound " +
           std::to_string( int( bd ) + 1 ) + ", branching " +
           std::to_string( int( br ) ) );
    }
  }
 {
  Flowshop f = Flowshop::taillard( 1 );
  f.branching = Flowshop::Branching::dynamic;
  BasicParallelBBOptions< int > o;
  o.threads = 4;
  o.incumbent = Flowshop::taillard_optimum( 1 ) + 1;
  const auto r = SkeletonBranchAndBound< Flowshop >( f ,
                  BBRule::stack_stealing , o ).solve();
  check( ( r.value == 1278 ) && r.node &&
         ( f.makespan( f.sequence( *r.node ) ) == 1278 ) ,
         "ta001 has the makespan 1278 (" + std::to_string( r.value ) +
         ", " + std::to_string( r.stats.evaluated ) + " nodes)" );
  }

 // UTS: the sizes of the trees against the plain recursion
 for( const UTS & u : { UTS::geometric( UTS::Shape::linear , 4 , 12 , 34 ) ,
                        UTS::geometric( UTS::Shape::cyclic , 3 , 5 , 50 ) ,
                        UTS::binomial_tree( 200 , 4 , 0.24 , 42 ) } ) {
  const long n = u.count( u.root() );
  for( const BBScheme sc : { BBScheme::sequential , BBScheme::stealing ,
                             BBScheme::budget , BBScheme::ordered } ) {
   BasicParallelBBOptions< long > o;
   o.threads = 4;
   o.mode = BBMode::enumerate;
   o.indexed = true;
   const auto r = parallel_branch_and_bound( u , sc , o );
   check( long( r.stats.solutions ) == n ,
          "UTS of " + std::to_string( n ) + " nodes (" +
          std::to_string( r.stats.solutions ) + "), " +
          std::string( bb_scheme_name( sc ) ) );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// an evaluator that throws: every scheme, the farm, the race and the
// intervals must stop all their workers and throw the exception of the
// evaluator again, which a watchdog checks they do in bounded time

namespace {

// the problem q whose evaluators throw std::domain_error at the call k of
// evaluate( n ), or of branch( n , sink ) and child( n , k ) with
// in_branch, counted from 1 over all the workers, and with every at each
// later call too

template< class Q >
struct Throwing {

 using node_type = typename Q::node_type;
 using value_type = bb_value_t< Q >;

 static constexpr Sense sense = Q::sense;

 using Inner = decltype( detail::bb_worker( std::declval< const Q & >() ,
                                            0u ) );

 struct Worker {
  Inner w;
  const Throwing * t;

  void cutoff( value_type z )
   requires requires( Inner & x ) { x.cutoff( z ); } {
   w.cutoff( z );
   }

  BasicEvaluation< value_type > evaluate( node_type & n ) {
   t->tick( false );
   return( w.evaluate( n ) );
   }

  template< class F >
  void branch( const node_type & n , F && f ) {
   t->tick( true );
   w.branch( n , std::forward< F >( f ) );
   }

  auto child( const node_type & n , std::size_t k )
   requires requires( Inner & x ) { x.child( n , k ); } {
   t->tick( true );
   return( w.child( n , k ) );
   }
  };

 const Q & q;
 std::atomic< std::uint64_t > & calls;
 std::uint64_t k;
 bool in_branch = false;
 bool every = false;

 node_type root( void ) const { return( q.root() ); }

 Worker make_worker( unsigned i ) const {
  return( Worker{ detail::bb_worker( q , i ) , this } );
  }

 std::size_t arity( unsigned d ) const
  requires requires( const Q & x ) { x.arity( d ); } {
  return( q.arity( d ) );
  }

 std::uint64_t hash( const node_type & n ) const
  requires HashedBranchAndBound< Q > {
  return( q.hash( n ) );
  }

 void serialize( const node_type & n , std::vector< char > & out ) const
  requires requires( const Q & x ) { x.serialize( n , out ); } {
  q.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const
  requires requires( const Q & x ) { x.deserialize( in ); } {
  return( q.deserialize( in ) );
  }

 void tick( bool branching ) const {
  if( branching != in_branch )
   return;
  const std::uint64_t c = ++calls;
  if( ( c == k ) || ( every && ( c > k ) ) )
   throw( std::domain_error( "Throwing::tick: call " +
                             std::to_string( c ) ) );
  }

 };  // end( struct( Throwing ) )

// the ways of throwing: at the call k of evaluate or branch, once or at
// every later call too
struct Throw {
 std::uint64_t k;
 bool in_branch;
 bool every;
 };

const std::vector< Throw > throws = { { 1 , false , false } ,
                                      { 10 , false , false } ,
                                      { 40 , false , false } ,
                                      { 10 , true , false } ,
                                      { 10 , false , true } };

std::string throw_name( const Throw & t )
{
 return( std::string( t.in_branch ? "branch" : "evaluate" ) + " call " +
         std::to_string( t.k ) + ( t.every ? " and later" : "" ) );
 }

// runs f, ending the tester with a failure if f has not returned after the
// given seconds, so that a search that hangs fails its case
template< class F >
void within( int seconds , const std::string & what , F && f )
{
 std::mutex m;
 std::condition_variable cv;
 bool done = false;
 std::thread dog( [ & ] {
  std::unique_lock< std::mutex > lock( m );
  if( ! cv.wait_for( lock , std::chrono::seconds( seconds ) ,
                     [ & ] { return( done ); } ) ) {
   std::cout << "FAILED: " << what << ", not ended after " << seconds
             << " s" << std::endl;
   std::_Exit( EXIT_FAILURE );
   }
  } );
 struct End {
  std::mutex & m;
  std::condition_variable & cv;
  bool & done;
  std::thread & dog;
  ~End() {
   {
    std::lock_guard< std::mutex > guard( m );
    done = true;
    }
   cv.notify_one();
   dog.join();
   }
  } end{ m , cv , done , dog };
 f();
 }

// true if solve() throws the std::domain_error of the evaluator, within
// the watchdog
template< class Solve >
bool throws_domain_error( Solve && solve , const std::string & what )
{
 bool ok = false;
 within( 120 , what , [ & ] {
  try {
   ( void ) solve();
   }
  catch( const std::domain_error & ) {
   ok = true;
   }
  } );
 return( ok );
 }

}  // end( namespace )

template< class Executor >
void test_bb_throw( const std::string & )
{
 const std::string ex = ex_name< Executor >();
 const Knapsack kp = small_knapsacks().at( 2 );
 const Flowshop f = random_flowshop( 7 , 4 , 3 );
 for( const Throw & t : throws ) {
  const std::string tag = throw_name( t ) + " on " + ex;
  std::atomic< std::uint64_t > calls{ 0 };
  const Throwing< Knapsack > q{ kp , calls , t.k , t.in_branch , t.every };
  ParallelBBOptions o;
  o.threads = 4;
  for( std::size_t s = 0 ; s < 11 ; ++s ) {
   const BBScheme sc = BBScheme( s );
   const std::string what = std::string( bb_scheme_name( sc ) ) + ", " +
                            tag;
   calls = 0;
   check( throws_domain_error( [ & ] {
           return( parallel_branch_and_bound< Throwing< Knapsack > ,
                                              Executor >( q , sc , o ) );
           } , what ) , "the exception of " + what );
   }

  calls = 0;
  check( throws_domain_error( [ & ] {
          return( FarmBranchAndBound< Throwing< Knapsack > , Executor >(
                   q , o ).solve() );
          } , "farm, " + tag ) , "the exception of the farm, " + tag );

  ParallelBBOptions r = o;
  r.race.resize( 3 );
  r.race[ 1 ].selection = Selection::depth_first;
  r.race[ 2 ].selection = Selection::plunging;
  r.race_nodes = 20;
  calls = 0;
  check( throws_domain_error( [ & ] {
          return( parallel_branch_and_bound< Throwing< Knapsack > ,
                                             Executor >( q ,
                                              BBScheme::shared , r ) );
          } , "race, " + tag ) , "the exception of the race, " + tag );

  std::atomic< std::uint64_t > fcalls{ 0 };
  const Throwing< Flowshop > qf{ f , fcalls , t.k , t.in_branch ,
                                 t.every };
  BasicParallelBBOptions< int > of;
  of.threads = 4;
  check( throws_domain_error( [ & ] {
          return( IntervalBranchAndBound< Throwing< Flowshop > , Executor >(
                   qf , of ).solve() );
          } , "intervals, " + tag ) ,
         "the exception of the intervals, " + tag );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the master-hub-worker search over MPI

#ifdef PBEFS_HAS_MPI
void test_bb_hub_mpi( const std::string & )
{
 int np;
 MPI_Comm_size( MPI_COMM_WORLD , & np );
 // with need_work 6 a hub whose cluster sums more than 6 open nodes may
 // has no worker above 6, and tells the master, which tries the next hub
 for( const std::size_t need : { std::size_t( 2 ) , std::size_t( 6 ) } )
  for( const unsigned cluster : { 1u , 2u , 4u } )
   for( const Knapsack & kp : small_knapsacks() ) {
    const long opt = kp.optimum();
    BBHubParams h;
    h.cluster = cluster;
    h.unit_work = 3;
    h.need_work = need;
    h.elite = cluster;
    ParallelBBOptions o;
    const auto s = MPIHubBranchAndBoundSearch< Knapsack >( kp , o ,
                                                           h ).solve();
    check( ( s.status == BBStatus::optimal ) && ( long( s.value ) == opt ) &&
           s.node && kp.check( *s.node , s.value ) ,
           "hubs of " + std::to_string( cluster ) + ", need_work " +
           std::to_string( need ) + " on " + std::to_string( np ) +
           " processes: optimum " + std::to_string( opt ) );
    o.node_limit = 20;
    const auto l = MPIHubBranchAndBoundSearch< Knapsack >( kp , o ,
                                                           h ).solve();
    check( l.bound >= opt - 1e-6 , "hubs: the bound of a stopped search" );
    }
 const UTS u = UTS::geometric( UTS::Shape::fixed , 4 , 6 , 19 );
 const long n = u.count( u.root() );
 BasicParallelBBOptions< long > o;
 o.mode = BBMode::enumerate;
 BBHubParams h;
 h.cluster = 2;
 const auto s = MPIHubBranchAndBoundSearch< UTS >( u , o , h ).solve();
 check( long( s.stats.solutions ) == n ,
        "hubs: the " + std::to_string( n ) + " nodes of a tree of UTS (" +
        std::to_string( s.stats.solutions ) + ")" );
 }

// an evaluator that throws over MPI, in both searches: the process whose
// evaluator throws throws its exception again, and the others
// std::runtime_error, all within the watchdog

void test_bb_throw_mpi( const std::string & )
{
 int np;
 MPI_Comm_size( MPI_COMM_WORLD , & np );
 const Knapsack kp = small_knapsacks().at( 2 );
 for( const Throw & t : throws ) {
  if( t.k > 10 )                       // the calls of each process
   continue;
  std::atomic< std::uint64_t > calls{ 0 };
  const Throwing< Knapsack > q{ kp , calls , t.k , t.in_branch , t.every };
  for( const bool hub : { false , true } ) {
   const std::string what = std::string( hub ? "hubs" : "MPI B&B" ) +
                            " on " + std::to_string( np ) +
                            " processes, " + throw_name( t );
   calls = 0;
   int mine[ 2 ] = { 0 , 1 };       // the evaluator's, none
   within( 120 , what , [ & ] {
    try {
     ParallelBBOptions o;
     o.batch = 4;
     if( hub ) {
      BBHubParams h;
      h.cluster = 2;
      h.unit_work = 3;
      ( void ) MPIHubBranchAndBoundSearch< Throwing< Knapsack > >( q , o ,
                                                                   h ).solve();
      }
     else
      ( void ) MPIBranchAndBoundSearch< Throwing< Knapsack > >( q ,
                                                                o ).solve();
     }
    catch( const std::domain_error & ) {
     mine[ 0 ] = 1;
     mine[ 1 ] = 0;
     }
    catch( const std::runtime_error & ) {
     mine[ 1 ] = 0;
     }
    } );
   int any[ 2 ];
   MPI_Allreduce( mine , any , 2 , MPI_INT , MPI_MAX , MPI_COMM_WORLD );
   check( ( any[ 0 ] == 1 ) && ( any[ 1 ] == 0 ) ,
          "the exception of " + what );
   }
  }
 }
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- THE CASES OF THE TESTER ------------------------*/
/*--------------------------------------------------------------------------*/

void test_bb_skeleton_thread( const std::string & d ) {
 test_bb_skeleton< ADS::ThreadPool >( d ); }

void test_bb_interval_thread( const std::string & d ) {
 test_bb_interval< ADS::ThreadPool >( d ); }

void test_bb_throw_thread( const std::string & d ) {
 test_bb_throw< ADS::ThreadPool >( d ); }

#ifdef PBEFS_HAS_OPENMP
void test_bb_skeleton_openmp( const std::string & d ) {
 test_bb_skeleton< ADS::OpenMPExecutor >( d ); }

void test_bb_interval_openmp( const std::string & d ) {
 test_bb_interval< ADS::OpenMPExecutor >( d ); }

void test_bb_throw_openmp( const std::string & d ) {
 test_bb_throw< ADS::OpenMPExecutor >( d ); }
#endif

#ifdef PBEFS_HAS_FASTFLOW
void test_bb_skeleton_fastflow( const std::string & d ) {
 test_bb_skeleton< ADS::FastFlowExecutor >( d ); }

void test_bb_interval_fastflow( const std::string & d ) {
 test_bb_interval< ADS::FastFlowExecutor >( d ); }

void test_bb_throw_fastflow( const std::string & d ) {
 test_bb_throw< ADS::FastFlowExecutor >( d ); }
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- End File bnb.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
