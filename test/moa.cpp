/*--------------------------------------------------------------------------*/
/*------------------------------ File moa.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the multi-objective
 * searches. The references do not come from the library: the Pareto fronts
 * of the small graphs and grids by the enumeration of all their simple
 * paths, and the traces of the examples of Hernández et al. (AIJ 2023,
 * figure 1 and tables 1 and 2) and of Hernández et al. (IJCAI 2023,
 * figure 1 and tables 1 and 2), whose fronts and numbers of extractions
 * are those of the papers; on the larger instances every search, with every
 * order of the objectives, has to find the front of the others.
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

#include <array>

#include <fstream>

#include <functional>

#include <iostream>

#include <numeric>

#include <random>

#include <sstream>

#include <string>

#include <vector>

#include "PBeFS/astar/MOAStar.h"

#include "PBeFS/astar/MultiObjectiveGraph.h"

#include "PBeFS/astar/MultiObjectiveGrid.h"

#include "PBeFS/astar/ParallelMOAStar.h"

#ifdef PBEFS_HAS_OPENMP
 #include "ADS/parallel/OpenMPExecutor.h"
#endif

#ifdef PBEFS_HAS_FASTFLOW
 #include "ADS/parallel/FastFlowExecutor.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

void check( bool ok , const std::string & what );    // [see test.cpp]

namespace {

template< std::size_t K >
using Graph = MultiObjectiveGraph< long , K >;

template< std::size_t K >
using Front = std::vector< std::array< long , K > >;

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
// the Pareto front of d by the enumeration of its simple paths

template< class D >
auto brute_front( const D & d )
{
 using State = mo_state_t< D >;
 using Vec = mo_cost_t< D >;
 std::vector< Vec > all;
 std::vector< State > path{ d.initial_state() };
 std::function< void( const Vec & ) > dfs = [ & ]( const Vec & g ) {
  const State s = path.back();
  if( d.is_goal( s ) ) {
   all.push_back( g );
   return;
   }
  d.for_each_successor( s , d.heuristic( s ) , [ & ]( const State & t ,
                                                      const Vec & c ,
                                                      const Vec & ) {
   if( std::find( path.begin() , path.end() , t ) != path.end() )
    return;
   path.push_back( t );
   dfs( add_costs( g , c ) );
   path.pop_back();
   } );
  };
 dfs( Vec{} );
 return( pareto_front( all ) );
 }

/*--------------------------------------------------------------------------*/
// true if the costs of sol are front and its paths are paths of d of their
// cost (any of the parallel arcs of a step)

template< class D , class S >
bool good_front( const D & d , const S & sol ,
                 const std::vector< mo_cost_t< D > > & front )
{
 using Vec = mo_cost_t< D >;
 if( sol.costs() != front )
  return( false );
 if( sol.status != ( front.empty() ? Status::unsolvable : Status::solved ) )
  return( false );
 for( const auto & p : sol.front ) {
  if( p.states.empty() || ! ( p.states.front() == d.initial_state() ) ||
      ! d.is_goal( p.states.back() ) )
   return( false );
  std::vector< Vec > sums{ Vec{} };
  for( std::size_t i = 1 ; i < p.states.size() ; ++i ) {
   std::vector< Vec > next;
   d.for_each_successor( p.states[ i - 1 ] ,
                         d.heuristic( p.states[ i - 1 ] ) ,
                         [ & ]( const auto & t , const Vec & c ,
                                const Vec & ) {
    if( t == p.states[ i ] )
     for( const Vec & s : sums )
      next.push_back( add_costs( s , c ) );
    } );
   std::sort( next.begin() , next.end() );
   next.erase( std::unique( next.begin() , next.end() ) , next.end() );
   sums = std::move( next );
   }
  if( std::find( sums.begin() , sums.end() , p.cost ) == sums.end() )
   return( false );
  }
 return( true );
 }

/*--------------------------------------------------------------------------*/
// a random graph of n nodes and m arcs with costs in [ 0 , c ], parallel
// arcs, loops and arcs of cost 0 included, searched from 0 to n - 1

template< std::size_t K >
Graph< K > random_graph( std::uint32_t n , std::size_t m , long c ,
                         std::mt19937_64 & gen , bool ideal )
{
 Graph< K > g( n , 0 , n - 1 );
 for( std::size_t a = 0 ; a < m ; ++a ) {
  std::array< long , K > w;
  for( auto & x : w )
   x = long( random_below( gen , std::uint64_t( c + 1 ) ) );
  g.add_arc( std::uint32_t( random_below( gen , n ) ) ,
             std::uint32_t( random_below( gen , n ) ) , w );
  }
 if( ideal )
  g.set_ideal_heuristic();
 return( g );
 }

/*--------------------------------------------------------------------------*/
// the sequential searches on d, in every order of the objectives, lazy
// and not: each one has to find front

template< class D >
void check_sequential( const D & d , const std::vector< mo_cost_t< D > > &
                       front , const std::string & what )
{
 constexpr std::size_t K = objectives_of< D >;
 std::array< unsigned , K > p;
 std::iota( p.begin() , p.end() , 0u );
 bool ok = true;
 std::string bad;
 do
  for( int a = 0 ; a < 4 ; ++a )
   for( bool lazy : { false , true } ) {
    if( ( a == 0 ) && ( K != 2 ) )
     continue;
    MOSearchOptions o;
    o.lazy = lazy;
    o.order.assign( p.begin() , p.end() );
    const auto s = mo_search( d , MOAlgorithm( a ) , o );
    if( ! good_front( d , s , front ) ) {
     ok = false;
     bad = " (algorithm " + std::to_string( a ) + ( lazy ? ", lazy" : "" ) +
           ")";
     }
    }
 while( std::next_permutation( p.begin() , p.end() ) );
 check( ok , "the sequential searches in every order, " + what + bad );
 }

/*--------------------------------------------------------------------------*/
// the parallel searches on d on the executor ex: each one, with 1 to 8
// workers, has to find front

template< class Executor , class D >
void check_parallel( Executor & ex , const D & d ,
                     const std::vector< mo_cost_t< D > > & front ,
                     const std::string & what )
{
 constexpr std::size_t K = objectives_of< D >;
 bool ok = true;
 std::string bad;
 for( unsigned p : { 1u , 2u , 3u , 4u , 8u } )
  for( bool lazy : { false , true } ) {
   ParallelMOOptions o;
   o.threads = p;
   o.lazy = lazy;
   for( bool share : { true , false } ) {
    o.share_bounds = share;
    MOAStarPar< D , ArrayFrontier , Executor > a( d , o );
    a.set_executor( ex );
    if( ! good_front( d , a.solve() , front ) ) {
     ok = false;
     bad = " (LTMOA*-par, " + std::to_string( p ) + " workers)";
     }
    }
   o.share_bounds = true;
   MOAStarPar< D , ArrayFrontier , Executor , MOBucketQueue > q( d , o );
   q.set_executor( ex );
   if( ! good_front( d , q.solve() , front ) ) {
    ok = false;
    bad = " (LTMOA*-par on buckets, " + std::to_string( p ) + " workers)";
    }
   MOAStarPar< D , EMOAFrontier , Executor > e( d , o );
   e.set_executor( ex );
   if( ! good_front( d , e.solve() , front ) ) {
    ok = false;
    bad = " (EMOA*-par, " + std::to_string( p ) + " workers)";
    }
   for( std::size_t batch : { std::size_t( 1 ) , std::size_t( 64 ) } ) {
    o.batch = batch;
    o.flush_period = unsigned( batch );
    o.window = batch == 1 ? -1 : 0;
    o.acknowledge = batch != 1;
    HDMOAStar< D , Executor > h( d , o );
    h.set_executor( ex );
    if( ! good_front( d , h.solve() , front ) ) {
     ok = false;
     bad = " (HDMOA*, " + std::to_string( p ) + " workers)";
     }
    if( batch == 1 ) {                  // the window without the counts
     o.window = 0;
     o.acknowledge = false;
     HDMOAStar< D , Executor > w( d , o );
     w.set_executor( ex );
     if( ! good_front( d , w.solve() , front ) ) {
      ok = false;
      bad = " (HDMOA* without acknowledge, " + std::to_string( p ) +
            " workers)";
      }
     }
    }
   }
 check( ok , "the parallel searches on " + ex_name< Executor >() + ", " +
        what + " (" + std::to_string( K ) + " objectives)" + bad );
 }

/*--------------------------------------------------------------------------*/
// the example of figure 1 of Hernández et al. (AIJ 2023), its 9 arcs on a
// path to the goal recovered from tables 1 and 2 (the tenth arc leaves no
// trace in them), with the perfect heuristic of the paper

Graph< 2 > boa_example( void )
{
 // the states s_start, s_1, s_2, s_3, s_4, s_goal
 Graph< 2 > g( 6 , 0 , 5 );
 g.add_arc( 0 , 1 , { 1 , 1 } );
 g.add_arc( 0 , 2 , { 1 , 5 } );
 g.add_arc( 0 , 3 , { 1 , 1 } );
 g.add_arc( 1 , 2 , { 1 , 2 } );
 g.add_arc( 1 , 4 , { 4 , 6 } );
 g.add_arc( 1 , 5 , { 7 , 5 } );
 g.add_arc( 2 , 5 , { 2 , 4 } );
 g.add_arc( 3 , 2 , { 2 , 1 } );
 g.add_arc( 3 , 5 , { 5 , 7 } );
 g.set_ideal_heuristic();
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the example of figure 1 of Hernández et al. (IJCAI 2023), the arcs and
// the heuristic recovered from tables 1 and 2

Graph< 3 > ltmoa_example( void )
{
 // the states s_start, s_1, s_2, s_goal
 Graph< 3 > g( 4 , 0 , 3 );
 g.add_arc( 0 , 1 , { 1 , 1 , 1 } );
 g.add_arc( 0 , 2 , { 1 , 3 , 2 } );
 g.add_arc( 1 , 2 , { 3 , 1 , 3 } );
 g.add_arc( 1 , 3 , { 1 , 3 , 4 } );
 g.add_arc( 2 , 3 , { 1 , 1 , 1 } );
 g.add_arc( 2 , 1 , { 1 , 1 , 1 } );
 g.set_heuristic( 0 , { 2 , 3 , 3 } );
 g.set_heuristic( 1 , { 1 , 2 , 4 } );
 g.set_heuristic( 2 , { 1 , 1 , 1 } );
 return( g );
 }

/*--------------------------------------------------------------------------*/
// the hand-checked examples: fronts and extractions of the papers

void check_examples( void )
{
 const auto b = boa_example();
 check( brute_front( b ) == Front< 2 >{ { 3 , 9 } , { 4 , 7 } , { 5 , 6 } } ,
        "the front of the example of BOA*" );
 auto s = MOAStar< Graph< 2 > , ScalarFrontier >( b ).solve();
 check( good_front( b , s , brute_front( b ) ) &&
        ( s.stats.extracted == 10 ) && ( s.stats.generated == 11 ) &&
        ( s.stats.inserted == 10 ) ,
        "BOA* on its example: 10 iterations, 2 of the 11 successors "
        "pruned when generated (tables 1 and 2)" );

 const auto l = ltmoa_example();
 check( brute_front( l ) == Front< 3 >{ { 2 , 4 , 3 } , { 5 , 3 , 5 } } ,
        "the front of the example of LTMOA*" );
 auto t = MOAStar< Graph< 3 > , ArrayFrontier >( l ).solve();
 check( good_front( l , t , brute_front( l ) ) &&
        ( t.stats.extracted == 7 ) ,
        "LTMOA* on its example: 7 iterations (table 1)" );
 MOSearchOptions lazy;
 lazy.lazy = true;
 t = MOAStar< Graph< 3 > , ArrayFrontier >( l , lazy ).solve();
 check( good_front( l , t , brute_front( l ) ) &&
        ( t.stats.extracted == 9 ) ,
        "LazyLTMOA* on its example: 9 iterations (table 2)" );
 t = NAMOAStarDR< Graph< 3 > >( l ).solve();
 check( good_front( l , t , brute_front( l ) ) &&
        ( t.stats.extracted + t.stats.duplicates == 7 ) &&
        ( t.stats.duplicates == 1 ) ,
        "NAMOA*dr on the example of LTMOA*: 7 iterations, one label "
        "removed from OPEN (section 4)" );
 t = MOAStar< Graph< 3 > , TreeFrontier >( l ).solve();
 check( good_front( l , t , brute_front( l ) ) &&
        ( t.stats.extracted == 7 ) ,
        "EMOA* on the example of LTMOA*: the extractions of LTMOA*" );
 }

/*--------------------------------------------------------------------------*/
// an index of vectors against a scan of all the vectors added, on random
// vectors of K components added and checked in any order, or mostly by
// increasing first component

template< std::size_t K , class Index >
void check_index( const std::string & name )
{
 std::mt19937_64 gen( 23 + K );
 bool ok = true;
 for( int r = 0 ; r < 20 ; ++r ) {
  Index b;
  std::vector< std::array< long , K > > all;
  const long c = 1 + r % 7;
  for( int k = 0 ; k < 600 ; ++k ) {
   std::array< long , K > v;
   for( auto & x : v )
    x = long( random_below( gen , std::uint64_t( c * 10 ) ) );
   if( r % 2 )                          // mostly by increasing v_1
    v[ 0 ] = k / 4 + long( random_below( gen , 3 ) );
   std::uint64_t n = 0;
   const bool d = b.dominated( v , n );
   const bool e = std::any_of( all.begin() , all.end() , [ & ]( auto & u ) {
    return( weakly_dominates( u , v ) );
    } );
   ok = ok && ( d == e );
   if( ( k % 3 ) || ! d ) {
    b.insert( v , n );
    all.push_back( v );
    }
   }
  ok = ok && ( b.size() == all.size() );
  }
 check( ok , name + " against a scan, " + std::to_string( K ) +
        " components" );
 }

/*--------------------------------------------------------------------------*/
// the front of d by LTMOA*, the reference of the instances too large for
// the enumeration

template< class D >
auto ltmoa_front( const D & d )
{
 return( MOAStar< D , ArrayFrontier >( d ).solve().costs() );
 }

/*--------------------------------------------------------------------------*/
// the edge cases: the start a goal, a goal not reachable, the node limit,
// the orders and the objectives refused

void check_edge_cases( void )
{
 Graph< 2 > g( 3 , 1 , 1 );
 g.add_arc( 1 , 0 , { 1 , 1 } );
 g.set_ideal_heuristic();
 for( int a = 0 ; a < 4 ; ++a ) {
  const auto s = mo_search( g , MOAlgorithm( a ) );
  check( good_front( g , s , Front< 2 >{ { 0 , 0 } } ) ,
         "the start is the goal, algorithm " + std::to_string( a ) );
  }
 g.set_target( 2 );
 g.set_ideal_heuristic();
 for( int a = 0 ; a < 4 ; ++a ) {
  const auto s = mo_search( g , MOAlgorithm( a ) );
  check( ( s.status == Status::unsolvable ) && s.front.empty() &&
         ( s.stats.extracted == 0 ) ,
         "a goal not reachable, a dead end at the start, algorithm " +
         std::to_string( a ) );
  }
 g.set_heuristic( 1 , { 0 , 0 } );
 g.set_heuristic( 0 , { 0 , 0 } );
 for( int a = 0 ; a < 4 ; ++a ) {
  const auto s = mo_search( g , MOAlgorithm( a ) );
  check( ( s.status == Status::unsolvable ) && s.front.empty() ,
         "a goal not reachable, h = 0, algorithm " + std::to_string( a ) );
  }

 std::mt19937_64 gen( 7 );
 const auto r = random_graph< 3 >( 30 , 200 , 9 , gen , true );
 const GridMap map( 10 , 10 );
 MultiObjectiveGrid< long , 3 > grid( map , 4 , 0 , 0 , 9 , 9 , 1 , 10 , 8 );
 const auto all = MOAStar< decltype( grid ) >( grid ).solve();
 const auto front = all.costs();
 MOSearchOptions lim;
 lim.node_limit = all.stats.stored / 3;
 for( int a = 1 ; a < 4 ; ++a ) {
  const auto s = mo_search( grid , MOAlgorithm( a ) , lim );
  bool ok = ( s.status == Status::node_limit ) && ( front.size() > 3 );
  for( const auto & v : s.costs() )
   ok = ok && ( std::find( front.begin() , front.end() , v ) !=
                front.end() );
  check( ok , "the node limit, only Pareto-optimal paths, algorithm " +
         std::to_string( a ) );
  }

 bool threw = false;
 try {
  MOSearchOptions o;
  o.order = { 0 , 0 };
  static_cast< void >( mo_search( g , MOAlgorithm::boa , o ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "an order that is not a permutation is refused" );
 threw = false;
 try {
  static_cast< void >( mo_search( r , MOAlgorithm::boa ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "BOA* on three objectives is refused" );
 }

/*--------------------------------------------------------------------------*/
// small random graphs of K objectives, with the ideal heuristic and with
// h = 0, against the enumeration of their paths: the sequential searches
// if ex is nullptr, the parallel ones on *ex otherwise

template< std::size_t K , class Executor >
void check_random_graphs( Executor * ex , int count )
{
 std::mt19937_64 gen( 11 + K );
 int bad = 0 , nonempty = 0;
 for( int k = 0 ; k < count ; ++k ) {
  const std::uint32_t n = std::uint32_t( 4 + random_below( gen , 6 ) );
  const std::size_t m = 2 * n + random_below( gen , 3 * n );
  const bool ideal = k % 4 != 3;
  const auto g = random_graph< K >( n , m , 1 + long( k % 6 ) , gen ,
                                    ideal );
  const auto front = brute_front( g );
  nonempty += ! front.empty();
  const std::string what = "random graph " + std::to_string( k );
  if( ex )
   check_parallel( *ex , g , front , what );
  else
   check_sequential( g , front , what );
  bad += ! good_front( g , MOAStar< Graph< K > >( g ).solve() , front );
  // the buckets of f_1 extract the labels in the order of the heap
  const auto a = MOAStar< Graph< K > , ArrayFrontier >( g ).solve();
  const auto b = MOAStar< Graph< K > , ArrayFrontier ,
                          MOBucketQueue >( g ).solve();
  bad += ! good_front( g , b , front ) ||
         ( a.stats.extracted != b.stats.extracted ) ||
         ( a.stats.expanded != b.stats.expanded ) ||
         ( a.stats.checks != b.stats.checks );
  }
 check( ( bad == 0 ) && ( nonempty > count / 2 ) ,
        std::to_string( count ) + " random graphs of " + std::to_string( K ) +
        " objectives against their enumeration" );
 }

/*--------------------------------------------------------------------------*/
// the grids: the small ones against the enumeration, those of Ren et al.
// (SoCS 2022, 10 x 10, costs in [ 1 , 10 ]) and with obstacles across the
// searches

template< std::size_t K , class Executor >
void check_grids( Executor * ex )
{
 const GridMap empty5( 5 , 5 ) , empty4( 4 , 3 );
 const GridMap walls = GridMap::random( 5 , 5 , 0.2 , 3 );
 const auto run = [ & ]( auto & d , const std::string & what , bool brute ) {
  const auto front = brute ? brute_front( d ) : ltmoa_front( d );
  check( ! front.empty() , what + ": a goal is reachable" );
  using G = std::remove_cvref_t< decltype( d ) >;
  const auto a = MOAStar< G , EMOAFrontier >( d ).solve();
  const auto b = MOAStar< G , EMOAFrontier , MOBucketQueue >( d ).solve();
  check( ( a.costs() == b.costs() ) &&
         ( a.stats.extracted == b.stats.extracted ) &&
         ( a.stats.expanded == b.stats.expanded ) ,
         what + ": the buckets of f_1 extract as the heap" );
  if( ex )
   check_parallel( *ex , d , front , what );
  else
   check_sequential( d , front , what );
  };
 MultiObjectiveGrid< long , K > a( empty5 , 4 , 0 , 0 , 4 , 4 , 1 , 10 , 5 );
 run( a , "the empty grid of 5 x 5 cells, 4 moves" , true );
 a.set_ideal_heuristic();
 run( a , "the empty grid of 5 x 5 cells, 4 moves, ideal h" , true );
 MultiObjectiveGrid< long , K > b( empty4 , 8 , 0 , 0 , 3 , 2 , 0 , 3 , 6 );
 run( b , "the empty grid of 4 x 3 cells, 8 moves, costs from 0" , true );
 MultiObjectiveGrid< long , K > c( walls , 8 , 0 , 4 , 4 , 4 , 1 , 10 , 7 );
 run( c , "a grid of 5 x 5 cells with obstacles, 8 moves" , true );
 const GridMap empty10( 10 , 10 );
 MultiObjectiveGrid< long , K > e( empty10 , 4 , 0 , 0 , 9 , 9 , 1 , 10 ,
                                   8 );
 e.set_ideal_heuristic();
 run( e , "the empty grid of 10 x 10 cells of Ren et al." , false );
 }

/*--------------------------------------------------------------------------*/
// the road graph of the sample of DIMACS with K objectives: the distances,
// a second file of times written from them, and the objectives computed
// from the graph, on some pairs of nodes

template< std::size_t K , class Executor >
void check_dimacs( Executor * ex , const std::string & data )
{
 std::ifstream in( data + "/dimacs/sample.gr" );
 std::stringstream d , t;
 d << in.rdbuf();
 // the second file: the same arcs, a time of each derived from its length
 for( std::string line ; std::getline( d , line ) ; ) {
  if( line.empty() || ( line[ 0 ] != 'a' ) ) {
   t << line << "\n";
   continue;
   }
  std::istringstream is( line );
  char a;
  long u , v , w;
  is >> a >> u >> v >> w;
  t << "a " << u << " " << v << " " << ( w * 7 ) % 97 + w / 3 << "\n";
  }
 std::vector< ExtraCost > extra;
 const ExtraCost all[] = { ExtraCost::degree , ExtraCost::random ,
                           ExtraCost::avg_degree , ExtraCost::arcs };
 for( std::size_t i = 2 ; i < K ; ++i )
  extra.push_back( all[ i - 2 ] );
 std::mt19937_64 gen( 13 );
 for( int k = 0 ; k < 3 ; ++k ) {
  d.clear();
  d.seekg( 0 );
  t.clear();
  t.seekg( 0 );
  const auto s = std::uint32_t( 1 + random_below( gen , 900 ) );
  const auto g0 = std::uint32_t( 1 + random_below( gen , 900 ) );
  auto g = read_dimacs_objectives< long , K >( { & d , & t } , extra , s ,
                                               g0 , 5 );
  g.set_ideal_heuristic();
  const auto front = ltmoa_front( g );
  const std::string what = "the DIMACS sample, pair " + std::to_string( k );
  check( ( g.size() == 900 ) && ( g.arcs() == 3106 ) && ! front.empty() ,
         what + ": read with a front" );
  if( ex )
   check_parallel( *ex , g , front , what );
  else
   check_sequential( g , front , what );
  }
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
// the sequential searches

void test_moa( const std::string & data )
{
 check_examples();
 check_edge_cases();
 check_index< 3 , LayerIndex< long , 3 > >( "LayerIndex" );
 check_index< 4 , LayerIndex< long , 4 > >( "LayerIndex" );
 check_index< 3 , BlockIndex< long , 3 > >( "BlockIndex" );
 check_index< 5 , BlockIndex< long , 5 > >( "BlockIndex" );
 ADS::ThreadPool * none = nullptr;
 check_random_graphs< 2 >( none , 60 );
 check_random_graphs< 3 >( none , 60 );
 check_random_graphs< 4 >( none , 30 );
 check_grids< 2 >( none );
 check_grids< 3 >( none );
 check_dimacs< 2 >( none , data );
 check_dimacs< 3 >( none , data );
 check_dimacs< 4 >( none , data );

 // the readers refuse the files that are not of one graph
 std::istringstream a( "p sp 2 1\na 1 2 3\n" ) , b( "p sp 2 1\na 2 1 3\n" );
 bool threw = false;
 try {
  static_cast< void >( read_dimacs_objectives< long , 2 >( { & a , & b } ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "two files of different arcs are refused" );
 std::istringstream c( "p sp 2 1\na 1 2 3\n" );
 threw = false;
 try {
  static_cast< void >( read_dimacs_objectives< long , 3 >( { & c } ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "objectives that are not K are refused" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel searches on the executor Executor

template< class Executor >
void test_moa_parallel( const std::string & data )
{
 Executor ex( 8 );
 check_random_graphs< 2 >( & ex , 12 );
 check_random_graphs< 3 >( & ex , 12 );
 check_random_graphs< 4 >( & ex , 6 );
 check_grids< 3 >( & ex );
 check_dimacs< 2 >( & ex , data );
 check_dimacs< 3 >( & ex , data );

 // the node limit
 const GridMap map( 10 , 10 );
 MultiObjectiveGrid< long , 3 > r( map , 4 , 0 , 0 , 9 , 9 , 1 , 10 , 8 );
 const auto all = MOAStar< decltype( r ) >( r ).solve();
 const auto front = all.costs();
 ParallelMOOptions o;
 o.threads = 3;
 o.node_limit = all.stats.stored / 6;
 MOAStarPar< decltype( r ) , ArrayFrontier , Executor > a( r , o );
 a.set_executor( ex );
 HDMOAStar< decltype( r ) , Executor > h( r , o );
 h.set_executor( ex );
 // MOAStarPar finds only Pareto-optimal paths, HDMOAStar, which corrects
 // its labels, paths that are weakly dominated by some of them
 const auto sa = a.solve() , sh = h.solve();
 bool ok = sa.status == Status::node_limit;
 for( const auto & v : sa.costs() )
  ok = ok && ( std::find( front.begin() , front.end() , v ) != front.end() );
 check( ok , "the node limit of LTMOA*-par on " + ex_name< Executor >() +
        ", only Pareto-optimal paths" );
 ok = sh.status == Status::node_limit;
 for( const auto & v : sh.costs() )
  ok = ok && std::any_of( front.begin() , front.end() , [ & ]( auto & u ) {
   return( weakly_dominates( u , v ) );
   } );
 check( ok , "the node limit of HDMOA* on " + ex_name< Executor >() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void test_moa_thread( const std::string & data )
{
 test_moa_parallel< ADS::ThreadPool >( data );
 }

#ifdef PBEFS_HAS_OPENMP
void test_moa_openmp( const std::string & data )
{
 test_moa_parallel< ADS::OpenMPExecutor >( data );
 }
#endif

#ifdef PBEFS_HAS_FASTFLOW
void test_moa_fastflow( const std::string & data )
{
 test_moa_parallel< ADS::FastFlowExecutor >( data );
 }
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- End File moa.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
