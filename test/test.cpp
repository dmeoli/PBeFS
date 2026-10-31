/*--------------------------------------------------------------------------*/
/*----------------------------- File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of PBeFS: each case, selected by the first argument, solves
 * instances of known optimal cost with one family of searches and checks
 * both the cost and the path, which must start at the initial state, end at
 * a goal, be made of edges of the Domain and cost what the Solution says.
 *
 * The references do not come from the library: a breadth-first search of
 * the whole state space of the 8-puzzle, the optimal costs of the 100
 * instances of Korf (1985), Dijkstra's algorithm on the grids, and the
 * distances on the map of Romania of Russell and Norvig.
 *
 *     test astar | grid | graph | limits | reopen
 *     test pdb | pdb-search
 *     test hda-thread | hda-openmp | hda-fastflow
 *     test farm-thread | farm-openmp | farm-fastflow
 *     test spa-thread | spa-openmp
 *     test pbnf-thread | pbnf-openmp
 *     test pase-thread
 *     test bb | bb-shared-thread | bb-dist-thread | bb-steal-thread
 *     test bb-shared-openmp | bb-dist-openmp | bb-steal-openmp
 *     test bb-shared-fastflow | bb-dist-fastflow | bb-steal-fastflow
 *     test bb-farm-thread | bb-farm-openmp | bb-farm-fastflow
 *     test bb-det-thread | bb-det-openmp | bb-det-fastflow
 *     test bb-split-thread | bb-split-openmp | bb-split-fastflow
 *     test bb-scheme | bb-cutoff | bb-values
 *     test mpi | bb-mpi             (under mpiexec)
 *     test registry <tests.txt>
 *
 * The cases are listed in test/tests.txt, which CMake and the makefile
 * read; the case registry checks that the list and the cases of the tester
 * are the same.
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

#include <atomic>

#include <cmath>

#include <compare>

#include <cstring>

#include <deque>

#include <filesystem>

#include <fstream>

#include <initializer_list>

#include <iostream>

#include <limits>

#include <map>

#include <mutex>

#include <numeric>

#include <optional>

#include <queue>

#include <random>

#include <set>

#include <sstream>

#include <stdexcept>

#include <string>

#include <string_view>

#include <tuple>

#include <type_traits>

#include <unordered_map>

#include <utility>

#include <vector>

#include "PBeFSConfig.h"

#include "PBeFS/astar/AStar.h"

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "PBeFS/bnb/DeterministicBranchAndBound.h"

#include "PBeFS/bnb/SelfSplitBranchAndBound.h"

#include "PBeFS/bnb/FarmBranchAndBound.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/ExplicitGraph.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/SPAStar.h"

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/FarmAStar.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/SafePBNF.h"

#include "PBeFS/astar/SlowDomain.h"

#include "PBeFS/astar/SlidingTile.h"

#ifdef PBEFS_HAS_OPENMP
 #include "ADS/parallel/OpenMPExecutor.h"
#endif

#ifdef PBEFS_HAS_FASTFLOW
 #include "ADS/parallel/FastFlowExecutor.h"
#endif

#ifdef PBEFS_HAS_MPI
 #include "PBeFS/astar/HDAStarMPI.h"
 #include "PBeFS/bnb/MPIBranchAndBound.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- GLOBAL VARIABLES ----------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;              // the checks failed so far
static bool quiet = false;            // print only the failures

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// records the outcome of a check

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 else
  if( ! quiet )
   std::cout << "ok: " << what << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the name of the executor in the names of the checks

template< class Executor >
static std::string executor_name( void )
{
 if constexpr( std::is_same_v< Executor , ADS::ThreadPool > )
  return( "thread" );
#ifdef PBEFS_HAS_OPENMP
 else
  if constexpr( std::is_same_v< Executor , ADS::OpenMPExecutor > )
   return( "openmp" );
#endif
 return( "fastflow" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the costs are equal, up to a relative 1e-9 for real costs

template< class Cost >
static bool same_cost( Cost a , Cost b )
{
 if constexpr( std::is_floating_point_v< Cost > )
  return( std::abs( a - b ) <= 1e-9 * std::max( Cost( 1 ) , std::abs( b ) ) );
 else
  return( a == b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that sol is solved with cost expected and a valid path

template< class D , class S >
static void check_solution( const D & d , const S & sol ,
                            cost_t< D > expected , const std::string & what )
{
 bool ok = sol.solved() && same_cost( sol.cost , expected ) &&
           ! sol.path.empty() && ( sol.path.front() == d.initial_state() ) &&
           d.is_goal( sol.path.back() );
 cost_t< D > sum = 0;
 for( std::size_t i = 1 ; ok && ( i < sol.path.size() ) ; ++i ) {
  bool edge = false;
  cost_t< D > best = infinite_cost< cost_t< D > >();
  d.for_each_successor( sol.path[ i - 1 ] ,
                        d.heuristic( sol.path[ i - 1 ] ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > c ,
                               cost_t< D > ) {
   if( ( t == sol.path[ i ] ) && ( c < best ) ) {
    edge = true;
    best = c;
    }
   } );
  ok = edge;
  sum += best;
  }
 ok = ok && same_cost( sum , sol.cost );
 std::ostringstream os;
 os << what << " (cost " << sol.cost << ", expected " << expected
    << ", expanded " << sol.stats.expanded << ")";
 check( ok , os.str() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS Blind -------------------------------*/
/*--------------------------------------------------------------------------*/
// a Domain with the heuristic of D replaced by 0: A* becomes Dijkstra

template< class D >
struct Blind {

 using state_type = state_t< D >;
 using cost_type = cost_t< D >;

 const D & d;

 state_type initial_state( void ) const { return( d.initial_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 cost_type heuristic( const state_type & ) const { return( 0 ); }

 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_type & t , cost_type c ,
                               cost_type ) { f( t , c , cost_type( 0 ) ); } );
  }

 };

/*--------------------------------------------------------------------------*/
/*------------------------------- INSTANCES --------------------------------*/
/*--------------------------------------------------------------------------*/
// the instances of Korf (1985) with their optimal costs

struct KorfInstance {
 int id;
 FifteenPuzzle::board_type board;
 int cost;
 };

static std::vector< KorfInstance > korf( const std::string & dir )
{
 std::ifstream in( dir + "/15-puzzle/korf100.txt" );
 std::vector< KorfInstance > v;
 std::string line;
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  std::istringstream is( line );
  KorfInstance k;
  k.id = int( v.size() ) + 1;
  for( auto & t : k.board ) {
   int x;
   is >> x;
   t = std::uint8_t( x );
   }
  is >> k.cost;
  v.push_back( k );
  }
 if( v.size() != 100 )
  throw( std::runtime_error( "cannot read korf100.txt from " + dir ) );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the easiest instances of Korf, by number (cost 41 to 45)

static const std::vector< int > easy_korf = { 55 , 16 , 42 , 79 , 71 , 12 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the map of Romania of Russell and Norvig, from a city to Bucharest

enum City { Arad , Bucharest , Craiova , Drobeta , Eforie , Fagaras ,
            Giurgiu , Hirsova , Iasi , Lugoj , Mehadia , Neamt , Oradea ,
            Pitesti , RimnicuVilcea , Sibiu , Timisoara , Urziceni , Vaslui ,
            Zerind , n_cities };

static ExplicitGraph< int > romania( City from )
{
 ExplicitGraph< int > g( n_cities , from , Bucharest );
 const int e[][ 3 ] = {
  { Arad , Sibiu , 140 } , { Arad , Zerind , 75 } ,
  { Arad , Timisoara , 118 } ,
  { Bucharest , Giurgiu , 90 } , { Bucharest , Urziceni , 85 } ,
  { Bucharest , Fagaras , 211 } , { Bucharest , Pitesti , 101 } ,
  { Craiova , Drobeta , 120 } , { Craiova , RimnicuVilcea , 146 } ,
  { Craiova , Pitesti , 138 } , { Drobeta , Mehadia , 75 } ,
  { Eforie , Hirsova , 86 } , { Fagaras , Sibiu , 99 } ,
  { Hirsova , Urziceni , 98 } , { Iasi , Vaslui , 92 } ,
  { Iasi , Neamt , 87 } ,
  { Lugoj , Timisoara , 111 } , { Lugoj , Mehadia , 70 } ,
  { Oradea , Zerind , 71 } , { Oradea , Sibiu , 151 } ,
  { Pitesti , RimnicuVilcea , 97 } , { RimnicuVilcea , Sibiu , 80 } ,
  { Urziceni , Vaslui , 142 } };
 for( const auto & x : e )
  g.add_edge( std::uint32_t( x[ 0 ] ) , std::uint32_t( x[ 1 ] ) , x[ 2 ] );
 // the straight-line distances to Bucharest
 const int sld[ n_cities ] = { 366 , 0 , 160 , 242 , 161 , 176 , 77 , 151 ,
                               226 , 244 , 241 , 234 , 380 , 100 , 193 , 253 ,
                               329 , 80 , 199 , 374 };
 for( int c = 0 ; c < n_cities ; ++c )
  g.set_heuristic( std::uint32_t( c ) , sld[ c ] );
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the distances to Bucharest, by Dijkstra on the map

static int romania_distance( City from )
{
 auto g = romania( from );
 return( astar( Blind< ExplicitGraph< int > >{ g } ).cost );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 8-puzzle with the blank moved from the goal through the given cells,
// each next to the one before

static EightPuzzle::board_type walk8( std::initializer_list< int > cells )
{
 auto b = EightPuzzle::korf_goal();
 int z = 0;
 for( int c : cells ) {
  std::swap( b[ z ] , b[ c ] );
  z = c;
  }
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the limit cases of a search of the 8-puzzle, run by solve( d , o ) for
// the domain d and the options o: the initial state that is the goal (cost
// 0, a path of one state), a state one move away, and a node limit of 5 on
// a state farther than that

template< class O , class Solve >
static void check_astar_edges( const O & base , Solve && solve ,
                               const std::string & name )
{
 const EightPuzzle g( EightPuzzle::korf_goal() );
 const auto s0 = solve( g , base );
 check( s0.solved() && ( s0.cost == 0 ) && ( s0.path.size() == 1 ) ,
        name + " initial state that is the goal" );

 const EightPuzzle one( walk8( { 1 } ) );
 const auto s1 = solve( one , base );
 check( s1.solved() && ( s1.cost == 1 ) && ( s1.path.size() == 2 ) ,
        name + " one move" );

 const EightPuzzle far( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                 0 , 3 , 4 , 7 , 6 } ) );
 O o = base;
 o.node_limit = 5;
 const auto sl = solve( far , o );
 check( ( astar( far ).cost >= 6 ) && ( sl.status == Status::node_limit ) &&
        sl.path.empty() , name + " node limit of 5" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the distance from the goal of every state of the 8-puzzle, by a
// breadth-first search from the goal (the moves are reversible)

static std::unordered_map< std::uint64_t , int > eight_distances(
 const EightPuzzle & d )
{
 std::unordered_map< std::uint64_t , int > dist;
 std::uint64_t goal = EightPuzzle::pack( EightPuzzle::korf_goal() );
 std::queue< std::uint64_t > q;
 dist[ goal ] = 0;
 q.push( goal );
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  d.for_each_successor( s , 0 , [ & ]( std::uint64_t t , int , int ) {
   if( dist.emplace( t , dist[ s ] + 1 ).second )
    q.push( t );
   } );
  }
 return( dist );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- CASES ----------------------------------*/
/*--------------------------------------------------------------------------*/
// sequential A*, both OPEN lists, on the 8-puzzle and the 15-puzzle

static void test_astar( const std::string & data )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 check( dist.size() == 181440 , "the 8-puzzle has 9! / 2 states" );
 std::size_t k = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 907 ) != 0 )                     // 200 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto sb = astar( d , OpenKind::buckets );
  const auto sh = astar( d , OpenKind::heap );
  if( ! ( sb.solved() && ( sb.cost == dst ) && sh.solved() &&
          ( sh.cost == dst ) ) ) {
   check_solution( d , sb , dst , "8-puzzle, buckets" );
   check_solution( d , sh , dst , "8-puzzle, heap" );
   }
  }
 check( failures == 0 , "A* on 200 states of the 8-puzzle" );

 const auto inst = korf( data );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  check_solution( d , astar( d , OpenKind::buckets ) , ki.cost ,
                  "A* buckets, Korf #" + std::to_string( id ) );
  check_solution( d , astar( d , OpenKind::heap ) , ki.cost ,
                  "A* heap, Korf #" + std::to_string( id ) );
  }

 for( OpenKind k : { OpenKind::buckets , OpenKind::heap } )
  check_astar_edges( SearchOptions() , [ k ]( const EightPuzzle & d ,
                                              const SearchOptions & o ) {
   return( astar( d , k , o ) );
   } , std::string( "A* " ) + ( k == OpenKind::heap ? "heap" : "buckets" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// both OPEN lists against an ordered set on random pushes and pops: the
// least f moving back below the f of the last pop, as in a worker of HDA*,
// and the values of g below and above BucketOpenList::DenseG, the first of
// these making the layers heaps while they hold entries

template< class Open >
static void check_open( const std::string & name )
{
 using N = Node< unsigned int , int >;
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed )
  for( int top_g : { 200 , 2000 } ) {
   std::mt19937_64 rng( seed );
   std::deque< N > nodes;
   Open open;
   std::multiset< std::pair< int , int > > ref;  // ( f , - g )
   bool ok = true;
   int last_f = 0;
   for( int i = 0 ; ( i < 20000 ) && ok ; ++i ) {
    if( ( rng() % 3 != 0 ) || ref.empty() ) {
     // a g up to top_g, grown along the run, and an f near the last pop
     const int g = static_cast< int >( rng() % ( 1 + ( top_g * i ) / 20000 ) );
     const int f = std::max( g , last_f - 20 +
                                 static_cast< int >( rng() % 60 ) );
     nodes.push_back( N{ 0u , g , f - g , nullptr , false , false } );
     open.push( &nodes.back() , g );
     ref.emplace( f , -g );
     }
    else {
     const auto e = open.pop();
     const auto it = ref.begin();
     ok = ( e.g == -it->second ) && ( e.g + e.node->h == it->first );
     last_f = it->first;
     ref.erase( it );
     }
    ok = ok && ( open.size() == ref.size() );
    }
   for( ; ok && ! ref.empty() ; ref.erase( ref.begin() ) ) {
    const auto e = open.pop();
    ok = ( e.g == -ref.begin()->second ) &&
         ( e.g + e.node->h == ref.begin()->first );
    }
   check( ok && open.empty() , name + ", seed " + std::to_string( seed ) +
                               ", g up to " + std::to_string( top_g ) );
   }
 }

static void test_open( void )
{
 using N = Node< unsigned int , int >;
 check_open< HeapOpenList< N > >( "the heap against an ordered set" );
 check_open< BucketOpenList< N > >( "the buckets against an ordered set" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cost of a cheapest path on the map m from ( sx , sy ) to ( gx , gy )
// with 4 or 8 moves, -1 if there is none: Dijkstra on the map itself, with
// the rules of Grid written again (leaving a cell costs its cost, times
// sqrt( 2 ) to a corner, whose two cut cells must be free), so that it does
// not depend on the adjacency of the set of data structures

static double grid_dijkstra( const GridMap & m , int moves , int sx , int sy ,
                             int gx , int gy )
{
 static constexpr int dx[ 8 ] = { 0 , -1 , 1 , 0 , -1 , 1 , -1 , 1 };
 static constexpr int dy[ 8 ] = { -1 , 0 , 0 , 1 , -1 , -1 , 1 , 1 };
 const int w = m.width();
 using Entry = std::pair< double , int >;
 std::vector< double > d( std::size_t( w ) * m.height() ,
                          std::numeric_limits< double >::infinity() );
 std::priority_queue< Entry , std::vector< Entry > ,
                      std::greater< Entry > > q;
 d[ sy * w + sx ] = 0;
 q.emplace( 0.0 , sy * w + sx );
 while( ! q.empty() ) {
  const auto [ du , u ] = q.top();
  q.pop();
  if( du > d[ u ] )
   continue;
  const int x = u % w , y = u / w;
  if( ( x == gx ) && ( y == gy ) )
   return( du );
  const double c = m.cost( x , y );
  for( int k = 0 ; k < moves ; ++k ) {
   const int nx = x + dx[ k ] , ny = y + dy[ k ];
   if( ! m.cost( nx , ny ) || ( ( k >= 4 ) && ( ! m.cost( nx , y ) ||
                                                ! m.cost( x , ny ) ) ) )
    continue;
   const double nd = du + ( k >= 4 ? c * std::sqrt( 2.0 ) : c );
   if( nd < d[ ny * w + nx ] ) {
    d[ ny * w + nx ] = nd;
    q.emplace( nd , ny * w + nx );
    }
   }
  }
 return( -1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* on random grids against Dijkstra, with 4 and 8 moves, on maps of unit
// costs and of costs 1 to 5

static void test_grid( void )
{
 for( std::uint64_t seed = 1 ; seed <= 6 ; ++seed ) {
  const auto m = GridMap::random( 300 , 200 , 0.25 , seed );
  const Grid< 4 > g4( m , 0 , 199 , 299 , 199 );
  const auto ref4 = astar( Blind< Grid< 4 > >{ g4 } );
  if( ref4.solved() )
   check_solution( g4 , astar( g4 ) , ref4.cost ,
                   "grid 4, seed " + std::to_string( seed ) );
  else
   check( ! astar( g4 ).solved() , "grid 4 unsolvable, seed " +
                                    std::to_string( seed ) );
  const Grid< 8 > g8( m , 0 , 199 , 299 , 199 );
  const auto ref8 = astar( Blind< Grid< 8 > >{ g8 } );
  if( ref8.solved() )
   check_solution( g8 , astar( g8 ) , ref8.cost ,
                   "grid 8, seed " + std::to_string( seed ) );
  }

 // the weighted map of the old FindPath example: costs 1 to 5, 9 = wall
 GridMap w( 20 , 20 );
 for( int y = 0 ; y < 20 ; ++y )
  for( int x = 0 ; x < 20 ; ++x )
   w.set_cost( x , y , 1 + ( 7 * x + 3 * y ) % 5 );
 for( int y = 2 ; y < 18 ; ++y )
  w.set_cost( 10 , y , 0 );
 const Grid< 4 > gw( w , 3 , 6 , 17 , 15 );
 check_solution( gw , astar( gw ) , astar( Blind< Grid< 4 > >{ gw } ).cost ,
                 "weighted grid 4" );
 const double dw4 = grid_dijkstra( w , 4 , 3 , 6 , 17 , 15 );
 check_solution( gw , astar( gw ) , int( std::lround( dw4 ) ) ,
                 "weighted grid 4 against Dijkstra on the map" );
 const Grid< 8 > gw8( w , 3 , 6 , 17 , 15 );
 const double dw8 = grid_dijkstra( w , 8 , 3 , 6 , 17 , 15 );
 check_solution( gw8 , astar( gw8 ) , dw8 ,
                 "weighted grid 8 against Dijkstra on the map" );

 // random maps of costs 1 to 5, whose moves cost differently both ways
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const auto mc = GridMap::random( 120 , 80 , 0.25 , seed , 5 );
  const std::string t = ", costs 1 to 5, seed " + std::to_string( seed );
  const double r4 = grid_dijkstra( mc , 4 , 0 , 79 , 119 , 79 );
  const Grid< 4 > c4( mc , 0 , 79 , 119 , 79 );
  if( r4 >= 0 )
   check_solution( c4 , astar( c4 ) , int( std::lround( r4 ) ) ,
                   "grid 4" + t );
  else
   check( ! astar( c4 ).solved() , "grid 4 unsolvable" + t );
  const double r8 = grid_dijkstra( mc , 8 , 0 , 79 , 119 , 79 );
  const Grid< 8 > c8( mc , 0 , 79 , 119 , 79 );
  if( r8 >= 0 )
   check_solution( c8 , astar( c8 ) , r8 , "grid 8" + t );
  else
   check( ! astar( c8 ).solved() , "grid 8 unsolvable" + t );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* on the map of Romania

static void test_graph( void )
{
 auto g = romania( Arad );
 const auto s = astar( g );
 check_solution( g , s , 418 , "Romania, Arad to Bucharest" );
 check( ( s.path.size() == 5 ) && ( s.path[ 1 ] == Sibiu ) &&
        ( s.path[ 2 ] == RimnicuVilcea ) && ( s.path[ 3 ] == Pitesti ) ,
        "Romania, through Sibiu, Rimnicu Vilcea and Pitesti" );
 for( int c = 0 ; c < n_cities ; ++c ) {
  auto gc = romania( City( c ) );
  check_solution( gc , astar( gc ) , romania_distance( City( c ) ) ,
                  "Romania, city " + std::to_string( c ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the node limit and the unsolvable instances

static void test_limits( const std::string & data )
{
 const auto inst = korf( data );
 const FifteenPuzzle d( inst[ 0 ].board );
 SearchOptions o;
 o.node_limit = 10000;
 const auto s = astar( d , OpenKind::automatic , o );
 check( ( s.status == Status::node_limit ) && ( s.stats.stored >= 10000 ) &&
        s.path.empty() , "A* stops at the node limit" );

 GridMap m( 50 , 50 );
 for( int y = 0 ; y < 50 ; ++y )
  m.set_cost( 25 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 49 , 49 );
 const auto u = astar( g );
 check( ( u.status == Status::unsolvable ) && ( u.stats.expanded == 25 * 50 ) ,
        "A* proves the wall unsolvable after the 1250 cells it reaches" );

 ParallelOptions po;
 po.threads = 4;
 const StateHashDistribution< Grid< 4 > > hd( g );
 const auto hu = HDAStar< Grid< 4 > , decltype( hd ) >( g , hd , po ).solve();
 check( hu.status == Status::unsolvable , "HDA* proves the wall unsolvable" );
 const auto su = SPAStar< Grid< 4 > >( g , po ).solve();
 check( su.status == Status::unsolvable , "SPA* proves the wall unsolvable" );

 po.node_limit = 5000;
 const ZobristDistribution< FifteenPuzzle > z( d );
 const auto hl = HDAStar< FifteenPuzzle , decltype( z ) >( d , z , po )
                  .solve();
 check( hl.status == Status::node_limit , "HDA* stops at the node limit" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the board of D reached from the goal of Korf by n random moves of the
// blank, none of which undoes the one before

template< class D >
static typename D::board_type random_walk( int n , std::uint64_t seed )
{
 auto b = D::korf_goal();
 std::mt19937_64 gen( seed );
 int blank = 0 , last = -1;
 for( int i = 0 ; i < n ; ++i ) {
  int next[ 4 ] , k = 0;
  const int r = blank / D::width , c = blank % D::width;
  if( r > 0 )
   next[ k++ ] = blank - D::width;
  if( c > 0 )
   next[ k++ ] = blank - 1;
  if( c < D::width - 1 )
   next[ k++ ] = blank + 1;
  if( r < D::height - 1 )
   next[ k++ ] = blank + D::width;
  int to;
  do
   to = next[ gen() % unsigned( k ) ];
  while( to == last );
  std::swap( b[ blank ] , b[ to ] );
  last = blank;
  blank = to;
  }
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// IDA* against the exact distances of the 8-puzzle, the optima of Korf on
// the 15-puzzle, A* on the grids and the map of Romania (where the parent
// is pruned by its state, with no operator), and the node limit

static void test_idastar( const std::string & data )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 std::size_t k = 0;
 int wrong = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 1814 ) != 0 )                    // 100 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto r = idastar( d );
  if( ! ( r.solved() && ( r.cost == dst ) ) ) {
   ++wrong;
   check_solution( d , r , dst , "IDA* 8-puzzle" );
   }
  }
 check( wrong == 0 , "IDA* on 100 states of the 8-puzzle" );

 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] , easy_korf[ 2 ] } ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto r = idastar( d );
  check_solution( d , r , ki.cost , "IDA*, Korf #" + std::to_string( id ) );
  check( r.stats.stored == std::uint64_t( ki.cost ) + 1 , "IDA* holds only "
         "the path, Korf #" + std::to_string( id ) );
  }

 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 7 , 5 , 0.2 , seed , 3 );
  const Grid< 4 > g4( m , 0 , 4 , 6 , 4 );
  const auto a4 = astar( g4 );
  if( a4.solved() )
   check_solution( g4 , idastar( g4 ) , a4.cost , "IDA* grid 4, seed " +
                   std::to_string( seed ) );
  const Grid< 8 > g8( m , 0 , 4 , 6 , 4 );
  const auto a8 = astar( g8 );
  if( a8.solved() )
   check_solution( g8 , idastar( g8 ) , a8.cost , "IDA* grid 8, seed " +
                   std::to_string( seed ) );
  }

 auto g = romania( Arad );
 check_solution( g , idastar( g ) , 418 , "IDA* Romania, Arad to Bucharest" );

 ExplicitGraph< int > dag( 4 , 0 , 3 );         // the target not reachable
 dag.add_edge( 0 , 1 , 1 , true );
 dag.add_edge( 1 , 2 , 1 , true );
 dag.add_edge( 0 , 2 , 1 , true );
 check( idastar( dag ).status == Status::unsolvable ,
        "IDA* proves an acyclic graph unsolvable" );

 SearchOptions o;
 o.node_limit = 1000;
 const FifteenPuzzle d( inst[ 0 ].board );
 const auto l = idastar( d , o );
 check( ( l.status == Status::node_limit ) && ( l.stats.expanded == 1000 ) &&
        l.path.empty() , "IDA* stops at the expansion limit" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 24-puzzle: the integer of 128 bits, the packed board, A* and IDA*
// with the Manhattan distance and with a pattern database, the parallel
// searches, and the 50 instances of the data

static void test_tiles24( const std::string & data )
{
 const Word128 one( 1 );
 const Word128 w = one << 70;
 check( ( w.lo == 0 ) && ( w.hi == 64 ) && ( ( w >> 70 ) == one ) &&
        ( countr_zero( w ) == 70 ) && ( countr_zero( Word128() ) == 128 ) &&
        ( ( Word128( ~std::uint64_t( 0 ) ) << 4 ).hi == 15 ) &&
        ( ( Word128( 0 , 1 ) >> 1 ).lo == ( std::uint64_t( 1 ) << 63 ) ) &&
        ( ( ~Word128() & w ) == w ) , "the shifts and the bits of Word128" );

 using D = TwentyFourPuzzle;
 static_assert( std::is_same_v< D::state_type , Word128 > &&
                std::is_trivially_copyable_v< Word128 > &&
                ! PerfectHashDomain< D > && PerfectHashDomain< FifteenPuzzle >
                && OperatorDomain< D > && FeatureDomain< D > );
 bool cells = true;
 for( int c = 0 ; c < D::cells ; ++c ) {
  auto b = D::korf_goal();
  std::swap( b[ 0 ] , b[ c ] );
  const auto s = D::pack( b );
  cells = cells && ( D::blank( s ) == c ) && ( D::unpack( s ) == b );
  for( int i = 0 ; i < D::cells ; ++i )
   cells = cells && ( D::tile( s , i ) == b[ i ] );
  }
 check( cells , "the blank and the tiles of every cell of the 24-puzzle" );

 PDBOptions po;
 Partition triples;
 for( int t = 1 ; t < D::cells ; t += 3 )
  triples.push_back( { t , t + 1 , t + 2 } );
 const AdditivePDB< 5 , 5 , NibbleEntries > pdb( triples , D::korf_goal() ,
                                                 po );
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const auto b = random_walk< D >( 40 , seed );
  const D d( b );
  const std::string tag = ", walk " + std::to_string( seed );
  const auto a = astar( d );
  check( a.solved() , "A* on the 24-puzzle" + tag );
  check_solution( d , a , a.cost , "A* on the 24-puzzle" + tag );
  check_solution( d , idastar( d ) , a.cost , "IDA* on the 24-puzzle" + tag );
  const SlidingTilePDB< 5 , 5 , NibbleEntries > dp( b , pdb );
  check( dp.heuristic( dp.initial_state() ) >= d.heuristic(
          d.initial_state() ) , "the 3-tile patterns dominate Manhattan" +
         tag );
  check_solution( dp , idastar( dp ) , a.cost , "IDA* with the 3-tile "
                  "patterns" + tag );
  ParallelOptions o;
  o.threads = 3;
  const ZobristDistribution< D > z( d );
  check_solution( d , HDAStar< D , decltype( z ) >( d , z , o ).solve() ,
                  a.cost , "HDA*[Z] on the 24-puzzle" + tag );
  check_solution( d , SPAStar< D >( d , o ).solve() , a.cost ,
                  "SPA* on the 24-puzzle" + tag );
  }

 std::ifstream in( data + "/24-puzzle/instances50.txt" );
 std::string line;
 int n = 0 , bad = 0;
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  try {
   const D d( D::parse( line ) );
   if( d.heuristic( d.initial_state() ) <= 0 )
    ++bad;
   }
  catch( std::invalid_argument & ) {
   ++bad;
   }
  ++n;
  }
 check( ( n == 50 ) && ( bad == 0 ) , "the 50 instances of the 24-puzzle "
        "are boards from which the goal is reachable" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// an admissible but inconsistent heuristic, which forces a reopening: A*
// must still find the optimal cost

static void test_reopen( void )
{
 // s = 0, goal = 4; the cheap path 0 - 1 - 2 - 3 - 4 costs 4, while the
 // heuristic of 1 lures A* into 0 - 2 first
 ExplicitGraph< int > g( 5 , 0 , 4 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 2 , 1 , true );
 g.add_edge( 0 , 2 , 3 , true );
 g.add_edge( 2 , 3 , 1 , true );
 g.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };   // admissible: h* = 4 3 2 1 0
 for( std::uint32_t u = 0 ; u < 5 ; ++u )
  g.set_heuristic( u , h[ u ] );
 const auto s = astar( g );
 check_solution( g , s , 4 , "A*, inconsistent heuristic" );
 check( s.stats.reexpanded > 0 , "A* reopens a closed node" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the distances of the placements of the tiles pat of the 8-puzzle to
// goal, by a breadth-first search of the abstraction written apart from
// the library: on the placements alone (the blank ignored), or by a 0-1
// search on the placements with the cell of the blank, keeping the least
// distance over the cells of the blank (the blank minimized)

static std::map< std::vector< int > , int > abstract_distances(
 const std::vector< int > & pat , const EightPuzzle::board_type & goal ,
 bool minimized )
{
 const auto next = []( int a , int b ) {
  return( std::abs( a / 3 - b / 3 ) + std::abs( a % 3 - b % 3 ) == 1 );
  };
 std::vector< int > g;
 for( int t : pat )
  for( int c = 0 ; c < 9 ; ++c )
   if( goal[ c ] == t )
    g.push_back( c );
 std::map< std::pair< std::vector< int > , int > , int > dist;
 std::deque< std::pair< std::vector< int > , int > > q;
 for( int b = 0 ; b < 9 ; ++b )
  if( minimized && ( std::find( g.begin() , g.end() , b ) == g.end() ) ) {
   dist[ { g , b } ] = 0;
   q.push_back( { g , b } );
   }
 if( ! minimized ) {
  dist[ { g , -1 } ] = 0;
  q.push_back( { g , -1 } );
  }
 while( ! q.empty() ) {
  const auto u = q.front();
  q.pop_front();
  const int d = dist[ u ];
  const auto & [ p , b ] = u;
  for( std::size_t i = 0 ; i < p.size() ; ++i )
   for( int c = 0 ; c < 9 ; ++c ) {
    if( ! next( p[ i ] , c ) ||
        ( std::find( p.begin() , p.end() , c ) != p.end() ) ||
        ( minimized && ( c != b ) ) )
     continue;
    auto v = p;
    v[ i ] = c;
    const std::pair< std::vector< int > , int > w{ v , minimized ? p[ i ]
                                                                 : -1 };
    if( dist.emplace( w , d + 1 ).second )
     q.push_back( w );
    else
     if( dist[ w ] > d + 1 ) {
      dist[ w ] = d + 1;
      q.push_back( w );
      }
    }
  if( minimized )                        // the blank moves for free
   for( int c = 0 ; c < 9 ; ++c ) {
    if( ! next( b , c ) || ( std::find( p.begin() , p.end() , c ) !=
                             p.end() ) )
     continue;
    const std::pair< std::vector< int > , int > w{ p , c };
    const auto it = dist.find( w );
    if( ( it == dist.end() ) || ( it->second > d ) ) {
     dist[ w ] = d;
     q.push_front( w );
     }
    }
  }
 std::map< std::vector< int > , int > h;
 for( const auto & [ u , d ] : dist ) {
  const auto it = h.find( u.first );
  if( ( it == h.end() ) || ( it->second > d ) )
   h[ u.first ] = d;
  }
 return( h );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cells of the tiles of the state s of the 8-puzzle

static std::array< std::uint8_t , 9 > cells8( std::uint64_t s )
{
 std::array< std::uint8_t , 9 > q{};
 SlidingTilePDB< 3 , 3 , ByteEntries >::cells_of( s , q.data() );
 return( q );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the heuristic of the additive tables p checked on every state of the
// 8-puzzle (dist, the distances to the goal of Korf): 0 at the goal, at
// most the distance, consistent along every move if so said, and, if not
// null, the values of the exact tables e no smaller; returns the number of
// the moves along which it is not consistent

template< class P , class E = P >
static std::size_t check_heuristic8(
 const P & p , const std::unordered_map< std::uint64_t , int > & dist ,
 bool consistent , const std::string & what , const E * e = nullptr )
{
 const auto h = [ & ]( std::uint64_t s ) {
  const auto q = cells8( s );
  if constexpr( requires { p.exact( q.data() ); } )
   return( p.exact( q.data() ) );
  else
   return( p.heuristic( q.data() ) );
  };
 const EightPuzzle any( EightPuzzle::korf_goal() );
 bool admissible = true , dominated = true;
 std::size_t bad = 0;
 for( const auto & [ s , d ] : dist ) {
  const int hs = h( s );
  admissible = admissible && ( hs <= d );
  if( e )
   dominated = dominated && ( hs <= e->heuristic( cells8( s ).data() ) );
  any.for_each_successor( s , 0 , [ & ]( std::uint64_t t , int , int ) {
   if( std::abs( hs - h( t ) ) > 1 )
    ++bad;
   } );
  }
 const std::uint64_t goal = EightPuzzle::pack( EightPuzzle::korf_goal() );
 check( h( goal ) == 0 , what + ": h = 0 at the goal" );
 check( admissible , what + ": admissible on all the 9! / 2 states" );
 if( e )
  check( dominated , what + ": at most the exact tables" );
 if( consistent )
  check( bad == 0 , what + ": consistent along every move" );
 return( bad );
 }

/*--------------------------------------------------------------------------*/
// the pattern databases: the values of the tables against a search of the
// abstraction written apart, the admissibility, the consistency and the
// encodings on the whole 8-puzzle, the reflection, the 15-puzzle against
// the optimal costs of Korf, and the files

static void test_pdb( const std::string & data )
{
 // the tables against the abstract search, with two goals, both rules
 auto centre = EightPuzzle::korf_goal();          // the blank in the middle
 std::swap( centre[ 0 ] , centre[ 4 ] );
 for( const auto & goal : { EightPuzzle::korf_goal() , centre } )
  for( const auto & tiles : std::vector< std::vector< int > >{
        { 1 } , { 1 , 2 , 3 } , { 5 , 6 , 7 , 8 } , { 2 , 4 , 6 , 8 } ,
        { 1 , 2 , 3 , 4 , 5 , 6 , 7 } } )
   for( BlankRule rule : { BlankRule::ignored , BlankRule::minimized } ) {
    const TilePattern< 3 , 3 > pat( tiles , goal );
    using Table = PatternTable< 3 , 3 , ByteEntries >;
    const auto h1 = Table::distances( pat , rule , 1 );
    const auto h3 = Table::distances( pat , rule , 3 );
    const auto ref = abstract_distances( tiles , goal ,
                                         rule == BlankRule::minimized );
    bool same = ( h1 == h3 ) && ( ref.size() == h1.size() );
    std::uint8_t p[ 9 ];
    for( std::uint64_t r = 0 ; same && ( r < h1.size() ) ; ++r ) {
     pat.unrank( r , p );
     const std::vector< int > v( p , p + tiles.size() );
     same = ( pat.rank( p ) == r ) && ( ref.at( v ) == h1[ r ] );
     }
    std::string name = "pattern {";
    for( int t : tiles )
     name += " " + std::to_string( t );
    check( same , name + " }, blank " +
                  ( rule == BlankRule::ignored ? "ignored" : "minimized" ) +
                  ( goal == centre ? ", blank in the middle" : "" ) +
                  ": the abstract distances, on 1 and 3 workers" );
    }

 // the heuristic on the whole 8-puzzle, two partitions
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 const auto goal8 = EightPuzzle::korf_goal();
 for( const Partition & part : { Partition{ { 1 , 2 , 3 , 4 } ,
                                            { 5 , 6 , 7 , 8 } } ,
                                 Partition{ { 1 , 2 , 3 } , { 4 , 5 , 6 } ,
                                            { 7 , 8 } } } ) {
  const std::string tag = part.size() == 2 ? " 4-4" : " 3-3-2";
  for( bool reflect : { false , true } ) {
   PDBOptions o;
   o.reflect = reflect;
   o.threads = 2;
   const std::string rt = tag + ( reflect ? " reflected" : "" );
   const AdditivePDB< 3 , 3 , ByteEntries > byte( part , goal8 , o );
   check( byte.reflected() == reflect , "the reflection" + rt );
   check_heuristic8( byte , dist , true , "byte" + rt );
   const AdditivePDB< 3 , 3 , NibbleEntries > nib( part , goal8 , o );
   check_heuristic8( nib , dist , true , "nibble" + rt , & byte );
   bool clamped = false;
   for( const auto & t : nib.tables() )
    clamped = clamped || t.entries().clamped();
   check( ! clamped , "nibble" + rt + ": no excess beyond 15" );
   o.compress = 4;
   const AdditivePDB< 3 , 3 , MinNibbleEntries > mn( part , goal8 , o );
   check_heuristic8( mn , dist , false , "min 4" + rt , & byte );
   o.compress = 1;
   o.blank = BlankRule::minimized;
   const AdditivePDB< 3 , 3 , ByteEntries > mb( part , goal8 , o );
   const auto bad = check_heuristic8( mb , dist , false ,
                                      "minimized" + rt );
   std::cout << "  minimized" << rt << ": inconsistent along " << bad
             << " moves" << std::endl;
   o.blank = BlankRule::ignored;
   bool dom = true;
   for( const auto & [ s , d ] : dist ) {
    const auto q = cells8( s );
    dom = dom && ( byte.heuristic( q.data() ) <= mb.heuristic( q.data() ) ) &&
          ( byte.heuristic( q.data() ) == nib.heuristic( q.data() ) );
    }
   check( dom , "nibble = byte <= minimized" + rt );
   }

  // the residues: without the reflection, the value of every state from
  // the residues alone, and that of every successor from its parent
  PDBOptions o;
  o.reflect = false;
  const AdditivePDB< 3 , 3 , ByteEntries > byte( part , goal8 , o );
  const AdditivePDB< 3 , 3 , Mod3Entries > m3( part , goal8 , o );
  const AdditivePDB< 3 , 3 , ParityBitEntries > bit( part , goal8 , o );
  check_heuristic8( m3 , dist , true , "mod 3" + tag , & byte );
  check_heuristic8( bit , dist , true , "parity bit" + tag , & byte );
  const EightPuzzle any( goal8 );
  bool exact3 = true , exact4 = true , step = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   const int hb = byte.heuristic( q.data() );
   exact3 = exact3 && ( m3.exact( q.data() ) == hb );
   exact4 = exact4 && ( bit.exact( q.data() ) == hb );
   const SlidingTilePDB< 3 , 3 , Mod3Entries > d3( EightPuzzle::unpack( s ) ,
                                                   m3 );
   const SlidingTilePDB< 3 , 3 , ParityBitEntries > d4(
    EightPuzzle::unpack( s ) , bit );
   d3.for_each_successor( s , hb , [ & ]( std::uint64_t t , int , int ht ) {
    step = step && ( ht == byte.heuristic( cells8( t ).data() ) );
    } );
   d4.for_each_successor( s , hb , [ & ]( std::uint64_t t , int , int ht ) {
    step = step && ( ht == byte.heuristic( cells8( t ).data() ) );
    } );
   }
  check( exact3 && exact4 , "mod 3 and parity bit" + tag +
                            ": the exact values from the residues alone" );
  check( step , "mod 3 and parity bit" + tag +
                ": the value of every successor from its parent" );
  bool thrown = false;
  try {
   o.reflect = true;
   const AdditivePDB< 3 , 3 , Mod3Entries > r3( part , goal8 , o );
   const SlidingTilePDB< 3 , 3 , Mod3Entries > d( goal8 , r3 );
   }
  catch( std::invalid_argument & ) {
   thrown = true;
   }
  check( thrown , "mod 3" + tag + ": the domain refuses the reflection" );
  }

 // the reflection: an automorphism of the puzzle that fixes the goal
 {
  const AdditivePDB< 3 , 3 , ByteEntries > p( Partition{ { 1 , 2 , 3 } ,
                                                         { 4 , 5 , 6 } ,
                                                         { 7 , 8 } } );
  const std::uint64_t goal = EightPuzzle::pack( goal8 );
  bool fixed = true , same = true , lookup = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   std::array< std::uint8_t , 9 > r{};
   p.reflection( q.data() , r.data() );
   EightPuzzle::board_type b{};
   for( int t = 0 ; t < 9 ; ++t )
    b[ r[ t ] ] = std::uint8_t( t );
   const std::uint64_t phi = EightPuzzle::pack( b );
   same = same && ( dist.at( phi ) == d );
   lookup = lookup && ( p.heuristic( q.data() ) ==
                        std::max( p.sum( q.data() ) , p.sum( r.data() ) ) );
   if( s == goal )
    fixed = phi == goal;
   }
  check( fixed && same && lookup , "the reflection fixes the goal, keeps "
         "the distances and gives the maximum of the two sums" );
  bool none = true;
  auto off = goal8;                     // the blank off the diagonal
  std::swap( off[ 0 ] , off[ 1 ] );
  none = ! AdditivePDB< 3 , 3 , ByteEntries >( Partition{ { 2 , 3 } } , off )
          .reflected();
  check( none , "no reflection when the blank of the goal is off the "
                "diagonal" );
  }

 // bad partitions
 bool thrown = false;
 try {
  const AdditivePDB< 3 , 3 , ByteEntries > p( Partition{ { 1 , 2 } ,
                                                         { 2 , 3 } } );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "a tile in two patterns is refused" );
 thrown = false;
 try {
  const TilePattern< 3 , 3 > p( { 0 , 1 } );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "the blank in a pattern is refused" );

 // the 15-puzzle, 5-5-5: admissible on the instances of Korf; the tables
 // in bytes are built once, in a cache, and encoded from there
 namespace fs = std::filesystem;
 const fs::path dir = fs::temp_directory_path() /
                      ( "pbefs-pdb-test-" + std::to_string(
                         std::random_device()() ) );
 fs::create_directories( dir );
 const auto inst = korf( data );
 for( BlankRule rule : { BlankRule::ignored , BlankRule::minimized } ) {
  PDBOptions o;
  o.blank = rule;
  o.cache = ( dir / "cache15" ).string();
  const AdditivePDB< 4 , 4 , ByteEntries > p( partition_15_5_5_5() ,
                                              FifteenPuzzle::korf_goal() ,
                                              o );
  const AdditivePDB< 4 , 4 , NibbleEntries > n( partition_15_5_5_5() ,
                                                FifteenPuzzle::korf_goal() ,
                                                o );
  o.compress = 8;
  const AdditivePDB< 4 , 4 , MinNibbleEntries > m( partition_15_5_5_5() ,
                                                   FifteenPuzzle::korf_goal() ,
                                                   o );
  bool ok = true;
  double sum = 0 , md = 0;
  for( const auto & ki : inst ) {
   const SlidingTilePDB< 4 , 4 , ByteEntries > d( ki.board , p );
   const FifteenPuzzle f( ki.board );
   std::array< std::uint8_t , 16 > q{};
   d.cells_of( d.initial_state() , q.data() );
   const int h = d.heuristic( d.initial_state() );
   ok = ok && ( h <= ki.cost ) && ( h >= f.heuristic( f.initial_state() ) )
        && ( n.heuristic( q.data() ) == h ) &&
        ( m.heuristic( q.data() ) <= h );
   sum += h;
   md += f.heuristic( f.initial_state() );
   }
  const std::string r = rule == BlankRule::ignored ? "ignored" : "minimized";
  check( ok , "5-5-5, blank " + r + ": Manhattan <= min 8 <= nibble = byte "
              "<= the optimum on the 100 instances of Korf" );
  std::cout << "  5-5-5, blank " << r << ": mean h " << sum / 100
            << " (Manhattan " << md / 100 << ")" << std::endl;
  std::array< std::uint8_t , 16 > g{};
  SlidingTilePDB< 4 , 4 , ByteEntries >::cells_of(
   FifteenPuzzle::pack( FifteenPuzzle::korf_goal() ) , g.data() );
  check( ( p.heuristic( g.data() ) == 0 ) && ( m.heuristic( g.data() ) == 0 ) ,
         "5-5-5, blank " + r + ": h = 0 at the goal" );
  }

 // the files: round trip, corruption, cache
 const TilePattern< 3 , 3 > pat( { 1 , 2 , 3 , 4 } );
 const auto h = PatternTable< 3 , 3 , ByteEntries >::distances(
                 pat , BlankRule::ignored );
 const auto round = [ & ]< class Enc >( const std::string & what ,
                                        unsigned c ) {
  using T = PatternTable< 3 , 3 , Enc >;
  const T t = T::encode( pat , BlankRule::ignored , h , c );
  const std::string f = ( dir / ( what + ".pdb" ) ).string();
  t.save( f );
  const T u = T::load( f , pat , BlankRule::ignored , c );
  check( u.from_file() && ! t.from_file() &&
         ( u.entries().payload() == t.entries().payload() ) ,
         what + ": saved and loaded" );
  };
 round.operator()< ByteEntries >( "byte" , 1 );
 round.operator()< NibbleEntries >( "nibble" , 1 );
 round.operator()< Mod3Entries >( "mod3" , 1 );
 round.operator()< ParityBitEntries >( "bit" , 1 );
 round.operator()< MinNibbleEntries >( "min" , 4 );
 round.operator()< SparseEntries >( "sparse" , 1 );
 {
  PDBOptions o;
  const Partition half{ { 1 , 2 , 3 , 4 } , { 5 , 6 , 7 , 8 } };
  const AdditivePDB< 3 , 3 , NibbleEntries > n( half , goal8 , o );
  const AdditivePDB< 3 , 3 , SparseEntries > sp( half , goal8 , o );
  bool same = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   same = same && ( n.heuristic( q.data() ) == sp.heuristic( q.data() ) );
   }
  check( same , "sparse = nibble on all the states of the 8-puzzle" );
  }

 using T = PatternTable< 3 , 3 , NibbleEntries >;
 const std::string f = ( dir / "nibble.pdb" ).string();
 const auto rejected = [ & ]( const std::string & what , auto && change ) {
  std::ifstream in( f , std::ios::binary );
  std::vector< char > b( ( std::istreambuf_iterator< char >( in ) ) ,
                         std::istreambuf_iterator< char >() );
  change( b );
  const std::string g = ( dir / "changed.pdb" ).string();
  std::ofstream( g , std::ios::binary ).write( b.data() ,
                                               std::streamsize( b.size() ) );
  bool bad = false;
  try {
   static_cast< void >( T::load( g , pat , BlankRule::ignored ) );
   }
  catch( std::runtime_error & ) {
   bad = true;
   }
  check( bad , "a file " + what + " is rejected" );
  };
 rejected( "with a byte changed" , []( auto & b ) {
  b[ b.size() / 2 ] ^= 1;
  } );
 rejected( "truncated" , []( auto & b ) { b.resize( b.size() - 9 ); } );
 rejected( "empty" , []( auto & b ) { b.clear(); } );
 rejected( "of another version" , []( auto & b ) { b[ 8 ] = 2; } );
 bool other = false;
 try {
  static_cast< void >( T::load( f , TilePattern< 3 , 3 >( { 1 , 2 , 3 , 5 } ) ,
                                BlankRule::ignored ) );
  }
 catch( std::runtime_error & ) {
  other = true;
  }
 check( other , "a file of another pattern is rejected" );

 const std::string cache = ( dir / "cache" ).string();
 const auto first = T::cached( cache , pat , BlankRule::ignored );
 const auto second = T::cached( cache , pat , BlankRule::ignored );
 check( ! first.from_file() && second.from_file() &&
        ( first.entries().payload() == second.entries().payload() ) ,
        "the cache builds a table once, then loads it" );
 const std::string own = ( fs::path( cache ) /
                           T::file_name( pat , BlankRule::ignored ) ).string();
 {
  std::fstream io( own , std::ios::binary | std::ios::in | std::ios::out );
  io.seekp( 60 );
  io.put( 'x' );
  }
 const auto third = T::cached( cache , pat , BlankRule::ignored );
 check( ! third.from_file() &&
        ( third.entries().payload() == first.entries().payload() ) ,
        "the cache encodes again a corrupted table" );
 fs::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A*, HDA* and SPA* with the pattern databases: the optimal costs, and
// fewer expansions than with the Manhattan distance

static void test_pdb_search( const std::string & data )
{
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 const Partition part{ { 1 , 2 , 3 , 4 } , { 5 , 6 , 7 , 8 } };
 const auto goal8 = EightPuzzle::korf_goal();
 PDBOptions o;
 const AdditivePDB< 3 , 3 , NibbleEntries > nib( part , goal8 , o );
 o.compress = 4;
 const AdditivePDB< 3 , 3 , MinNibbleEntries > mn( part , goal8 , o );
 o.compress = 1;
 o.blank = BlankRule::minimized;
 const AdditivePDB< 3 , 3 , NibbleEntries > mz( part , goal8 , o );
 o.blank = BlankRule::ignored;
 o.reflect = false;
 const AdditivePDB< 3 , 3 , ParityBitEntries > bit( part , goal8 , o );
 std::size_t k = 0 , em = 0 , ep = 0;
 for( const auto & [ s , d ] : dist ) {
  if( ( k++ % 907 ) != 0 )                     // 200 states, spread out
   continue;
  const auto b = EightPuzzle::unpack( s );
  const EightPuzzle m( b );
  const auto sm = astar( m );
  const SlidingTilePDB< 3 , 3 , NibbleEntries > dn( b , nib );
  const auto sn = astar( dn );
  const SlidingTilePDB< 3 , 3 , MinNibbleEntries > dmn( b , mn );
  const SlidingTilePDB< 3 , 3 , NibbleEntries > dmz( b , mz );
  const SlidingTilePDB< 3 , 3 , ParityBitEntries > db( b , bit );
  const auto s1 = astar( dmn );
  const auto s2 = astar( dmz );
  const auto s3 = astar( db );
  if( ! ( sn.solved() && ( sn.cost == d ) && ( s1.cost == d ) &&
          ( s2.cost == d ) && ( s3.cost == d ) ) ) {
   check_solution( dn , sn , d , "A* nibble, 8-puzzle" );
   check_solution( dmn , s1 , d , "A* min 4, 8-puzzle" );
   check_solution( dmz , s2 , d , "A* minimized, 8-puzzle" );
   check_solution( db , s3 , d , "A* parity bit, 8-puzzle" );
   }
  em += sm.stats.expanded;
  ep += sn.stats.expanded;
  }
 check( failures == 0 , "A* with the tables on 200 states of the 8-puzzle" );
 check( ep < em , "A* with 4-4 expands fewer nodes than with Manhattan (" +
                  std::to_string( ep ) + " against " + std::to_string( em ) +
                  ")" );

 const auto inst = korf( data );
 const AdditivePDB< 4 , 4 , NibbleEntries > p( partition_15_5_5_5() );
 PDBOptions nr;
 nr.reflect = false;
 const AdditivePDB< 4 , 4 , ParityBitEntries > pb( partition_15_5_5_5() ,
                                                  FifteenPuzzle::korf_goal() ,
                                                  nr );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle m( ki.board );
  const SlidingTilePDB< 4 , 4 , NibbleEntries > d( ki.board , p );
  const SlidingTilePDB< 4 , 4 , ParityBitEntries > db( ki.board , pb );
  const auto sm = astar( m );
  const auto sp = astar( d );
  const auto sh = astar( d , OpenKind::heap );
  const auto sb = astar( db );
  const std::string tag = ", Korf #" + std::to_string( id );
  check_solution( d , sp , ki.cost , "A* 5-5-5" + tag );
  check_solution( d , sh , ki.cost , "A* heap 5-5-5" + tag );
  check_solution( db , sb , ki.cost , "A* 5-5-5 parity bit" + tag );
  check( sp.stats.expanded < sm.stats.expanded , "A* 5-5-5 expands fewer "
         "nodes than Manhattan" + tag + " (" +
         std::to_string( sp.stats.expanded ) + " against " +
         std::to_string( sm.stats.expanded ) + ")" );
  }

 using D = SlidingTilePDB< 4 , 4 , NibbleEntries >;
 for( unsigned t : { 1u , 3u , 8u } ) {
  ParallelOptions po;
  po.threads = t;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const D d( ki.board , p );
   const std::string tag = " p " + std::to_string( t ) + ", Korf #" +
                           std::to_string( id );
   const ZobristDistribution< D > z( d );
   check_solution( d , HDAStar< D , decltype( z ) , ADS::ThreadPool >(
                        d , z , po ).solve() , ki.cost , "HDA*[Z] 5-5-5" +
                                                         tag );
   check_solution( d , SPAStar< D , ADS::ThreadPool >( d , po ).solve() ,
                   ki.cost , "SPA* 5-5-5" + tag );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// HDA* with every distribution, 1 to 8 workers, on the given runtime

template< class Executor >
static void test_hda( const std::string & data )
{
 const auto inst = korf( data );
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  ParallelOptions o;
  o.threads = p;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   const std::string tag = " p " + std::to_string( p ) + ", Korf #" +
                           std::to_string( id );
   const ZobristDistribution< FifteenPuzzle > z( d );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( z ) , Executor >(
                        d , z , o ).solve() , ki.cost , "HDA*[Z]" + tag );
   const AbstractZobristDistribution< FifteenPuzzle > az( d );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( az ) , Executor >(
                        d , az , o ).solve() , ki.cost , "HDA*[AZ]" + tag );
   const AbstractStateDistribution< FifteenPuzzle > a( d );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( a ) , Executor >(
                        d , a , o ).solve() , ki.cost , "HDA*[A]" + tag );
   const PerfectHashDistribution< FifteenPuzzle > ph( d );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( ph ) , Executor >(
                        d , ph , o ).solve() , ki.cost , "HDA*[P]" + tag );
   }

  // the random owner on the 8-puzzle, whose 9! / 2 states bound the p
  // copies that it may store of each state, whatever the timing
  std::size_t k = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 9007 ) != 0 )                   // 21 states, spread out
    continue;
   const EightPuzzle e( EightPuzzle::unpack( s ) );
   const RandomDistribution< EightPuzzle > r( e );
   const auto sol = HDAStar< EightPuzzle , decltype( r ) , Executor >(
                     e , r , o ).solve();
   if( ! sol.solved() || ( sol.cost != dst ) )
    check_solution( e , sol , dst , "HDA*[R] 8-puzzle, p " +
                                    std::to_string( p ) );
   }
  check( failures == 0 , "HDA*[R] on 21 states of the 8-puzzle, p " +
                         std::to_string( p ) );

  const auto m = GridMap::random( 300 , 200 , 0.25 , 2 );
  const Grid< 8 > g( m , 0 , 199 , 299 , 199 , 20 );
  const auto ref = astar( Blind< Grid< 8 > >{ g } );
  const AbstractZobristDistribution< Grid< 8 > > gz( g );
  check_solution( g , HDAStar< Grid< 8 > , decltype( gz ) , Executor >(
                       g , gz , o ).solve() , ref.cost ,
                  "HDA*[AZ] grid 8, p " + std::to_string( p ) );

  auto rg = romania( Arad );
  const ZobristDistribution< ExplicitGraph< int > > rz( rg );
  check_solution( rg , HDAStar< ExplicitGraph< int > , decltype( rz ) ,
                                Executor >( rg , rz , o ).solve() , 418 ,
                  "HDA*[Z] Romania, p " + std::to_string( p ) );
  }

 const auto hda = []( const EightPuzzle & d , const ParallelOptions & o ) {
  const ZobristDistribution< EightPuzzle > z( d );
  return( HDAStar< EightPuzzle , decltype( z ) , Executor >( d , z , o )
          .solve() );
  };
 ParallelOptions o;
 o.threads = 8;
 check_astar_edges( o , hda , "HDA* 8 workers" );
 o.threads = 0;
 o.batch = 0;
 o.flush_period = 0;
 check_astar_edges( o , hda , "HDA* parameters 0" );
 o = ParallelOptions();
 o.threads = 3;
 o.batch = 1;
 o.flush_period = 1;
 check_astar_edges( o , hda , "HDA* batches of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// SPA* with 1 to 8 workers, on the given runtime

template< class Executor >
static void test_spa( const std::string & data )
{
 const auto inst = korf( data );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  ParallelOptions o;
  o.threads = p;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   check_solution( d , SPAStar< FifteenPuzzle , Executor >( d , o ).solve() ,
                   ki.cost , "SPA* p " + std::to_string( p ) + ", Korf #" +
                             std::to_string( id ) );
   }
  auto rg = romania( Arad );
  check_solution( rg , SPAStar< ExplicitGraph< int > , Executor >( rg , o )
                        .solve() , 418 ,
                  "SPA* Romania, p " + std::to_string( p ) );
  }

 const auto spa = []( const EightPuzzle & d , const ParallelOptions & o ) {
  return( SPAStar< EightPuzzle , Executor >( d , o ).solve() );
  };
 ParallelOptions o;
 o.threads = 8;
 check_astar_edges( o , spa , "SPA* 8 workers" );
 o.threads = 0;
 check_astar_edges( o , spa , "SPA* 0 workers" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// SafePBNF with 1 to 8 workers, on the given runtime

template< class Executor >
static void test_pbnf( const std::string & data )
{
 const auto inst = korf( data );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  PBNFOptions o;
  o.threads = p;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   check_solution( d , SafePBNF< FifteenPuzzle , Executor >( d , o ).solve() ,
                   ki.cost , "SafePBNF p " + std::to_string( p ) +
                             ", Korf #" + std::to_string( id ) );
   }
  for( std::uint64_t seed : { 2u , 3u } ) {
   const auto m = GridMap::random( 300 , 200 , 0.25 , seed );
   const Grid< 4 > g4( m , 0 , 199 , 299 , 199 , 100 , 10 );
   const auto r4 = astar( Blind< Grid< 4 > >{ g4 } );
   if( r4.solved() )
    check_solution( g4 , SafePBNF< Grid< 4 > , Executor >( g4 , o ).solve() ,
                    r4.cost , "SafePBNF grid 4, p " + std::to_string( p ) );
   const Grid< 8 > g8( m , 0 , 199 , 299 , 199 , 100 , 10 );
   const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
   if( r8.solved() )
    check_solution( g8 , SafePBNF< Grid< 8 > , Executor >( g8 , o ).solve() ,
                    r8.cost , "SafePBNF grid 8, p " + std::to_string( p ) );
   }
  }

 GridMap m( 60 , 60 );
 for( int y = 0 ; y < 60 ; ++y )
  m.set_cost( 30 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 59 , 59 , 100 , 7 );
 PBNFOptions o;
 o.threads = 4;
 check( SafePBNF< Grid< 4 > , Executor >( g , o ).solve().status ==
        Status::unsolvable , "SafePBNF proves the wall unsolvable" );

 const auto pbnf = []( const EightPuzzle & d , const PBNFOptions & x ) {
  return( SafePBNF< EightPuzzle , Executor >( d , x ).solve() );
  };
 o.threads = 8;
 check_astar_edges( o , pbnf , "SafePBNF 8 workers" );
 o.threads = 0;
 check_astar_edges( o , pbnf , "SafePBNF 0 workers" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// PA*SE with 1 to 8 workers: every state is expanded at most once

template< class Executor >
static void test_pase( const std::string & data )
{
 const auto inst = korf( data );
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  PASEOptions o;
  o.threads = p;
  std::size_t k = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 9007 ) != 0 )                   // 21 states, spread out
    continue;
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const auto sol = PASE< EightPuzzle , Executor >( d , o ).solve();
   if( ! sol.solved() || ( sol.cost != dst ) )
    check_solution( d , sol , dst , "PA*SE 8-puzzle, p " +
                                    std::to_string( p ) );
   }
  check( failures == 0 , "PA*SE on 21 states of the 8-puzzle, p " +
                         std::to_string( p ) );

  const auto m = GridMap::random( 200 , 120 , 0.25 , 3 );
  const Grid< 8 > g8( m , 0 , 119 , 199 , 119 );
  const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
  const SlowDomain< Grid< 8 > > slow( g8 , 2000 );
  if( r8.solved() ) {
   check_solution( g8 , PASE< Grid< 8 > , Executor >( g8 , o ).solve() ,
                   r8.cost , "PA*SE grid 8, p " + std::to_string( p ) );
   const auto ss = PASE< SlowDomain< Grid< 8 > > , Executor >( slow , o )
                    .solve();
   check_solution( slow , ss , r8.cost , "PA*SE slow grid 8, p " +
                                          std::to_string( p ) );
   check( ss.stats.reexpanded == 0 , "PA*SE expands each state once" );
   }

  o.eps = 1.5;
  const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto se = PASE< FifteenPuzzle , Executor >( d , o ).solve();
  check( se.solved() && ( se.cost >= ki.cost ) &&
         ( se.cost <= 1.5 * ki.cost ) ,
         "PA*SE eps 1.5 within the bound, p " + std::to_string( p ) );
  }

 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 8 > g( m , 0 , 0 , 39 , 39 );
 PASEOptions o;
 o.threads = 4;
 check( PASE< Grid< 8 > , Executor >( g , o ).solve().status ==
        Status::unsolvable , "PA*SE proves the wall unsolvable" );

 const auto pase = []( const EightPuzzle & d , const PASEOptions & x ) {
  return( PASE< EightPuzzle , Executor >( d , x ).solve() );
  };
 o.threads = 8;
 check_astar_edges( o , pase , "PA*SE 8 workers" );
 o.threads = 0;
 check_astar_edges( o , pase , "PA*SE 0 workers" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the master-worker FarmAStar, 1 to 8 workers, 1 or 4 nodes per task

template< class Executor >
static void test_farm( const std::string & data )
{
 const auto inst = korf( data );
 for( unsigned p : { 1u , 2u , 3u , 8u } )
  for( std::size_t b : { std::size_t( 1 ) , std::size_t( 4 ) } ) {
   ParallelOptions o;
   o.threads = p;
   o.batch = b;
   const std::string tag = " p " + std::to_string( p ) + " batch " +
                           std::to_string( b );
   const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
   const FifteenPuzzle d( ki.board );
   check_solution( d , FarmAStar< FifteenPuzzle , Executor >( d , o )
                       .solve() , ki.cost , "FarmAStar Korf #" +
                             std::to_string( easy_korf[ 0 ] ) + tag );
   const auto m = GridMap::random( 200 , 120 , 0.25 , 3 );
   const Grid< 8 > g8( m , 0 , 119 , 199 , 119 );
   const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
   if( r8.solved() )
    check_solution( g8 , FarmAStar< Grid< 8 > , Executor >( g8 , o )
                         .solve() , r8.cost , "FarmAStar grid 8" + tag );
   auto rg = romania( Arad );
   check_solution( rg , FarmAStar< ExplicitGraph< int > , Executor >( rg ,
                                                                      o )
                        .solve() , 418 , "FarmAStar Romania" + tag );
   }
 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 39 , 39 );
 ParallelOptions o;
 o.threads = 3;
 check( FarmAStar< Grid< 4 > , Executor >( g , o ).solve().status ==
        Status::unsolvable , "FarmAStar proves the wall unsolvable" );
 o.node_limit = 300;
 check( FarmAStar< Grid< 4 > , Executor >( g , o ).solve().status ==
        Status::node_limit , "FarmAStar stops at the node limit" );

 const auto farm = []( const EightPuzzle & d , const ParallelOptions & x ) {
  return( FarmAStar< EightPuzzle , Executor >( d , x ).solve() );
  };
 o = ParallelOptions();
 o.threads = 8;
 check_astar_edges( o , farm , "FarmAStar 8 workers" );
 o.threads = 0;
 o.batch = 0;
 check_astar_edges( o , farm , "FarmAStar parameters 0" );
 }


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#ifdef PBEFS_HAS_MPI

// HDA* over MPI with every deterministic distribution

static void test_mpi( const std::string & data )
{
 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 2 ] } ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const std::string tag = ", Korf #" + std::to_string( id );
  ParallelOptions o;
  o.batch = 32;
  const ZobristDistribution< FifteenPuzzle > z( d );
  check_solution( d , HDAStarMPI< FifteenPuzzle , decltype( z ) >( d , z , o )
                       .solve() , ki.cost , "MPI HDA*[Z]" + tag );
  const AbstractZobristDistribution< FifteenPuzzle > az( d );
  check_solution( d , HDAStarMPI< FifteenPuzzle , decltype( az ) >( d , az ,
                       o ).solve() , ki.cost , "MPI HDA*[AZ]" + tag );
  const AbstractStateDistribution< FifteenPuzzle > a( d );
  check_solution( d , HDAStarMPI< FifteenPuzzle , decltype( a ) >( d , a , o )
                       .solve() , ki.cost , "MPI HDA*[A]" + tag );
  const PerfectHashDistribution< FifteenPuzzle > ph( d );
  check_solution( d , HDAStarMPI< FifteenPuzzle , decltype( ph ) >( d , ph ,
                       o ).solve() , ki.cost , "MPI HDA*[P]" + tag );
  }
 const TwentyFourPuzzle d24( random_walk< TwentyFourPuzzle >( 40 , 1 ) );
 const ZobristDistribution< TwentyFourPuzzle > z24( d24 );
 check_solution( d24 , HDAStarMPI< TwentyFourPuzzle , decltype( z24 ) >(
                        d24 , z24 ).solve() , astar( d24 ).cost ,
                 "MPI HDA*[Z] on the 24-puzzle" );

 auto rg = romania( Arad );
 const StateHashDistribution< ExplicitGraph< int > > rh( rg );
 check_solution( rg , HDAStarMPI< ExplicitGraph< int > , decltype( rh ) >(
                       rg , rh ).solve() , 418 , "MPI HDA*[H] Romania" );

 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 39 , 39 );
 const StateHashDistribution< Grid< 4 > > gh( g );
 check( HDAStarMPI< Grid< 4 > , decltype( gh ) >( g , gh ).solve().status ==
        Status::unsolvable , "MPI HDA* proves the wall unsolvable" );

 ParallelOptions lo;
 lo.node_limit = 2000;
 const FifteenPuzzle d( inst[ 0 ].board );
 const ZobristDistribution< FifteenPuzzle > z( d );
 check( HDAStarMPI< FifteenPuzzle , decltype( z ) >( d , z , lo ).solve()
         .status == Status::node_limit , "MPI HDA* stops at the node limit" );

 const auto hda = []( const EightPuzzle & e , const ParallelOptions & x ) {
  const ZobristDistribution< EightPuzzle > ez( e );
  return( HDAStarMPI< EightPuzzle , decltype( ez ) >( e , ez , x ).solve() );
  };
 ParallelOptions o;
 check_astar_edges( o , hda , "MPI HDA*" );
 o.batch = 0;
 o.flush_period = 0;
 check_astar_edges( o , hda , "MPI HDA* parameters 0" );
 }

#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack as a minimization of the opposite profit, which exercises
// the other sense

struct MinKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::minimize;

 const Knapsack & kp;

 node_type root( void ) const { return( kp.root() ); }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e = kp.evaluate( n );
  if( e.infeasible )
   return( e );
  e.bound = - e.bound;
  e.feasible = - *e.feasible;
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  kp.branch( n , [ & ]( node_type && c , double e ) {
   f( std::move( c ) , - e );
   } );
  }

 void serialize( const node_type & n , std::vector< char > & out ) const {
  kp.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const {
  return( kp.deserialize( in ) );
  }

 };  // end( struct( MinKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack whose workers take the incumbent by cutoff( z ), and find
// infeasible a node whose bound cannot exceed z, as a relaxation stopped
// by its cutoff does; a z above the optimum is recorded in wrong

struct CutoffKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::maximize;

 struct Worker {
  const Knapsack & kp;
  std::atomic< bool > & wrong;
  double opt;
  double z = - std::numeric_limits< double >::infinity();

  void cutoff( double v ) {
   if( v > opt )
    wrong = true;
   z = v;
   }

  Evaluation evaluate( node_type & n ) const {
   Evaluation e = kp.evaluate( n );
   if( ( ! e.infeasible ) && ( e.bound <= z ) ) {
    e = Evaluation();
    e.infeasible = true;
    }
   return( e );
   }

  template< class F >
  void branch( const node_type & n , F && f ) const {
   kp.branch( n , std::forward< F >( f ) );
   }
  };

 const Knapsack & kp;
 std::atomic< bool > & wrong;
 double opt;                      // the optimum of kp

 node_type root( void ) const { return( kp.root() ); }

 Worker make_worker( unsigned ) const {
  return( Worker{ kp , wrong , opt } );
  }

 };  // end( struct( CutoffKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack whose worker i counts its evaluations in count[ i ]

struct CountKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::maximize;

 struct Worker {
  const Knapsack & kp;
  std::atomic< std::uint64_t > & count;

  Evaluation evaluate( node_type & n ) const {
   ++count;
   return( kp.evaluate( n ) );
   }

  template< class F >
  void branch( const node_type & n , F && f ) const {
   kp.branch( n , std::forward< F >( f ) );
   }
  };

 const Knapsack & kp;
 std::vector< std::atomic< std::uint64_t > > & count;

 node_type root( void ) const { return( kp.root() ); }

 Worker make_worker( unsigned i ) const {
  return( Worker{ kp , count.at( i ) } );
  }

 };  // end( struct( CountKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the instances of the knapsack of the tests: 2 of each class of Pisinger

static std::vector< Knapsack > knapsacks( void )
{
 std::vector< Knapsack > v;
 for( int cls : { 0 , 1 , 2 } )
  for( std::uint64_t seed : { 1 , 2 } )
   v.push_back( Knapsack::random( cls , cls == 2 ? 30 : 50 , 1000 ,
                                  seed ) );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the leaves of the B&B on the knapsack kp, solved by solve( q , o , seeds )
// for the problem q (Knapsack or MinKnapsack), the options o and the nodes
// to start from (none for the root): the leaves of a search carried out to
// the end bound the optimum and cover the tree, so that a search started
// from them solves the knapsack of a smaller capacity, in both senses; a
// search started from the open leaves of one stopped by a node limit, with
// its incumbent, ends it; and a search started from no node ends at once.
// With explicit children, the leaves of the smaller knapsack, the
// infeasible ones among them, cover the tree of the larger one, which a
// search started from them solves again. With MPI, each process has its
// own leaves, which may be none.

using Leaves = std::vector< BBLeaf< Knapsack::node_type > >;

template< class O , class Solve >
static void check_frontier_of( const Knapsack & kp , O o , Solve && solve ,
                               const std::string & name , bool local )
{
 constexpr double inf = std::numeric_limits< double >::infinity();
 const double opt = double( kp.optimum() );
 Knapsack kq = kp;
 kq.cap = kp.cap * 3 / 4;
 const double optq = double( kq.optimum() );
 o.keep_frontier = true;

 const auto s = solve( kp , o , std::optional< Leaves >() );
 bool ok = ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           ( local || ! s.frontier.empty() );
 Leaves seeds;
 for( const auto & l : s.frontier ) {
  ok = ok && ( l.kind != BBLeafKind::open ) && ( l.bound <= opt + 1e-9 );
  seeds.push_back( l );
  seeds.back().bound = inf;            // unknown, in the maximization
  }
 check( ok , name + " leaves of the whole tree" );
 const auto r = solve( kq , o , std::optional< Leaves >( seeds ) );
 check( ( r.status == BBStatus::optimal ) && ( r.value == optq ) &&
        r.node && kq.check( *r.node , r.value ) ,
        name + " from the leaves, a smaller capacity" );
 if( kp.explicit_children ) {
  seeds = r.frontier;
  ok = true;
  for( auto & x : seeds ) {
   ok = ok && ( ( x.kind == BBLeafKind::infeasible ) == ( x.bound == - inf ) );
   x.bound = inf;
   }
  const auto b = solve( kp , o , std::optional< Leaves >( seeds ) );
  check( ok && ( b.status == BBStatus::optimal ) && ( b.value == opt ) &&
         b.node && kp.check( *b.node , b.value ) ,
         name + " from the leaves, the capacity back" );
  }

 const auto m = solve( MinKnapsack{ kp } , o , std::optional< Leaves >() );
 ok = m.status == BBStatus::optimal;
 seeds.clear();
 for( const auto & l : m.frontier ) {
  ok = ok && ( l.bound >= - opt - 1e-9 );
  seeds.push_back( l );
  seeds.back().bound = - inf;          // unknown, in the minimization
  }
 const auto mr = solve( MinKnapsack{ kq } , o ,
                        std::optional< Leaves >( seeds ) );
 check( ok && ( mr.status == BBStatus::optimal ) && ( mr.value == - optq ) &&
        mr.node && kq.check( *mr.node , - mr.value ) ,
        name + " from the leaves, a smaller capacity, minimization" );

 O ol = o;
 ol.node_limit = 5;
 const auto l = solve( kp , ol , std::optional< Leaves >() );
 seeds.clear();
 for( const auto & x : l.frontier )
  if( x.kind == BBLeafKind::open )
   seeds.push_back( x );
 O oc = o;
 if( l.feasible() )
  oc.incumbent = l.value;
 const auto c = solve( kp , oc , std::optional< Leaves >( seeds ) );
 check( ( l.status == BBStatus::node_limit ) &&
        ( local || ! seeds.empty() ) && ( c.status == BBStatus::optimal ) &&
        ( c.value == opt ) && ( c.gap <= 1e-9 ) ,
        name + " from the open leaves of a node limit" );

 const auto e = solve( kp , o , std::optional< Leaves >( Leaves() ) );
 check( ( e.status == BBStatus::infeasible ) && ( e.stats.evaluated == 0 ) &&
        e.frontier.empty() , name + " from no node" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the same, with the children of the knapsack getting the bound of their
// parent and with explicit children, which prunes some of them when they
// are generated and gives infeasible leaves

template< class O , class Solve >
static void check_frontier( const Knapsack & kp , const O & o ,
                            Solve && solve , const std::string & name ,
                            bool local = false )
{
 Knapsack kb = kp;
 kb.explicit_children = true;
 check_frontier_of( kp , o , solve , name , local );
 check_frontier_of( kb , o , solve , name + " explicit children" , local );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a tree of a root and two children without children, both of bound 0:
// the first with the feasible value 1, the second with none, so that
// neither is pruned and both are leaves, and the optimum is 1

struct Fork {

 using node_type = int;

 static constexpr Sense sense = Sense::minimize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e;
  if( n == 1 )
   e.feasible = 1;
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  if( n == 0 ) {
   f( 1 , 0.0 );
   f( 2 , 0.0 );
   }
  }

 };  // end( struct( Fork ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a root that is infeasible

struct NoRoot {

 using node_type = int;

 static constexpr Sense sense = Sense::minimize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & ) const {
  Evaluation e;
  e.infeasible = true;
  return( e );
  }

 template< class F >
 void branch( const node_type & , F && ) const {}

 };  // end( struct( NoRoot ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a root solved without branching: its bound is its feasible value, 5

struct Leaf {

 using node_type = int;

 static constexpr Sense sense = Sense::maximize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & ) const {
  Evaluation e;
  e.bound = 5;
  e.feasible = 5;
  return( e );
  }

 template< class F >
 void branch( const node_type & , F && ) const {}

 };  // end( struct( Leaf ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the limit cases of a B&B solved by run( q , o ) for the problem q and the
// options o: a root that is infeasible, a root solved without branching,
// the tree of Fork, and, on the knapsack kp, an initial incumbent better
// than the optimum (which nothing improves), an absolute tolerance, a time
// limit that is reached at once, a node limit of 1, and plunging with the
// dive going always and never

template< class O , class Run >
static void check_bb_edges( const Knapsack & kp , const O & base , Run && run ,
                            const std::string & name )
{
 const double opt = double( kp.optimum() );

 const auto n = run( NoRoot{} , base );
 check( ( n.status == BBStatus::infeasible ) && ! n.node && ! n.feasible() ,
        name + " infeasible root" );

 const auto l = run( Leaf{} , base );
 check( ( l.status == BBStatus::optimal ) && ( l.value == 5 ) &&
        ( l.bound == 5 ) && ( l.gap == 0 ) && ( l.stats.evaluated == 1 ) ,
        name + " root solved without branching" );

 const auto f = run( Fork{} , base );
 check( ( f.status == BBStatus::optimal ) && ( f.value == 1 ) &&
        ( f.bound == 1 ) , name + " fork" );

 O o = base;
 o.incumbent = opt + 1;
 const auto z = run( kp , o );
 check( ( z.status == BBStatus::optimal ) && ( z.value == opt + 1 ) &&
        ! z.node && ( z.bound >= opt + 1 - 1e-9 ) ,
        name + " initial incumbent better than the optimum" );

 o = base;
 o.abs_gap = 0.05 * opt;
 const auto a = run( kp , o );
 check( ( a.status == BBStatus::optimal ) &&
        ( a.value >= opt - o.abs_gap - 1e-9 ) && ( a.bound >= opt - 1e-9 ) &&
        ( a.gap <= o.abs_gap + 1e-9 ) , name + " absolute tolerance" );

 o = base;
 o.time_limit = 1e-12;
 const auto t = run( kp , o );
 check( ( t.status == BBStatus::time_limit ) && ( t.bound >= opt - 1e-9 ) &&
        ( ! t.feasible() || ( t.value <= opt ) ) , name + " time limit" );

 o = base;
 o.node_limit = 1;
 const auto one = run( kp , o );
 check( ( one.status == BBStatus::node_limit ) &&
        ( one.stats.evaluated >= 1 ) && ( one.bound >= opt - 1e-9 ) ,
        name + " node limit of 1" );

 for( double g : { 0.0 , 1.0 } ) {
  o = base;
  o.selection = Selection::plunging;
  o.plunge_gap = g;
  const auto p = run( kp , o );
  check( ( p.status == BBStatus::optimal ) && ( p.value == opt ) ,
         name + " plunging with plunge_gap " + std::to_string( g ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the sequential B&B on the knapsack against the dynamic program, with
// every selection, in both senses, with a relative tolerance, with an
// initial incumbent and with a node limit, its leaves, and its limit cases

static void test_bb( const std::string & )
{
 const auto solve = []( const auto & q , const BBOptions & o ,
                        std::optional< Leaves > seeds ) {
  BranchAndBoundSearch< std::decay_t< decltype( q ) > > b( q , o );
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   BBOptions o;
   o.selection = sel;
   const auto s = branch_and_bound( kp , o );
   const std::string t = tag + " selection " +
                         std::to_string( int( sel ) );
   check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
          s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
          "B&B optimum" + t );
   const auto m = branch_and_bound( MinKnapsack{ kp } , o );
   check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
          m.node && kp.check( *m.node , - m.value ) ,
          "B&B minimization" + t );
   check_frontier( kp , o , solve , "B&B" + t );
   }

  BBOptions o;
  o.rel_gap = 0.05;
  const auto r = branch_and_bound( kp , o );
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) &&
         ( r.gap <= 0.05 * r.value + 1e-9 ) ,
         "B&B within a relative gap of 5%" + tag );

  o = BBOptions();
  o.incumbent = opt;
  double seen = 0;
  o.on_incumbent = [ & ]( double v ) { seen = v; };
  const auto z = branch_and_bound( kp , o );
  check( ( z.status == BBStatus::optimal ) && ( z.value == opt ) &&
         ! z.node && ( seen == 0 ) , "B&B with the optimum given" + tag );

  o = BBOptions();
  o.node_limit = 3;
  const auto l = branch_and_bound( kp , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) &&
         ( l.stats.evaluated == 3 ) , "B&B node limit" + tag );
  }

 BBOptions o;
 o.keep_frontier = true;
 const auto f = branch_and_bound( Fork{} , o );
 std::set< std::pair< int , int > > leaves;
 for( const auto & x : f.frontier )
  leaves.insert( { x.node , int( x.kind ) } );
 check( ( f.status == BBStatus::optimal ) && ( f.value == 1 ) &&
        ( f.bound == 1 ) &&
        ( leaves == std::set< std::pair< int , int > >{
                     { 1 , int( BBLeafKind::feasible ) } ,
                     { 2 , int( BBLeafKind::pruned ) } } ) ,
        "B&B leaves without children" );

 check_bb_edges( knapsacks().front() , BBOptions() ,
                 []( const auto & q , const BBOptions & x ) {
                  return( branch_and_bound( q , x ) );
                  } , "B&B" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a parallel B&B of the problem q: K = 0 the shared queue, K = 1 the
// distributed queues, K = 2 the farm, K = 3 the work stealing

template< int K , class Executor , class Q >
static auto bb_parallel( const Q & q , const ParallelBBOptions & o ,
                         std::optional< std::vector< BBLeaf<
                          typename Q::node_type > > > seeds = std::nullopt )
{
 const auto solve = [ & ]( const auto & b ) {
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 if constexpr( K == 0 )
  return( solve( SharedBranchAndBound< Q , Executor >( q , o ) ) );
 else
  if constexpr( K == 1 )
   return( solve( DistributedBranchAndBound< Q , Executor >( q , o ) ) );
  else
   if constexpr( K == 3 )
    return( solve( StealingBranchAndBound< Q , Executor >( q , o ) ) );
   else
    return( solve( FarmBranchAndBound< Q , Executor >( q , o ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a parallel B&B on the knapsack against the dynamic program, with 1 to 8
// workers and every selection, in both senses (the minimization has no
// hash( n ), and hence random owners), with a relative tolerance and with
// a node limit, its leaves, and its limit cases with 8 workers (more than
// the nodes of the small trees) and with the parameters at 0, which mean 1

template< int K , class Executor >
static void test_bb_parallel( const std::string & )
{
 const std::string name = std::string( K == 0 ? "shared" :
                                       K == 1 ? "distributed" :
                                       K == 2 ? "farm" : "stealing" ) +
                          " B&B on " + executor_name< Executor >();
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( unsigned p : { 1u , 2u , 3u , 8u } )
   for( Selection sel : { Selection::best_first , Selection::depth_first ,
                          Selection::breadth_first ,
                          Selection::best_estimate , Selection::plunging } ) {
    ParallelBBOptions o;
    o.threads = p;
    o.selection = sel;
    const auto s = bb_parallel< K , Executor >( kp , o );
    check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
           name + tag + " p " + std::to_string( p ) + " selection " +
           std::to_string( int( sel ) ) );
    }

  ParallelBBOptions o;
  o.threads = 4;
  const auto m = bb_parallel< K , Executor >( MinKnapsack{ kp } , o );
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         name + " minimization" + tag );

  o.rel_gap = 0.05;
  const auto r = bb_parallel< K , Executor >( kp , o );
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) ,
         name + " within a relative gap of 5%" + tag );

  o = ParallelBBOptions();
  o.threads = 4;
  o.node_limit = 3;
  const auto l = bb_parallel< K , Executor >( kp , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) ,
         name + " node limit" + tag );

  for( Selection sel : { Selection::best_first , Selection::plunging } ) {
   o = ParallelBBOptions();
   o.threads = 3;
   o.selection = sel;
   check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                                std::optional< Leaves > seeds ) {
    return( bb_parallel< K , Executor >( q , x , std::move( seeds ) ) );
    } , name + tag + " selection " + std::to_string( int( sel ) ) );
   }
  }

 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( bb_parallel< K , Executor >( q , x ) );
  };
 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.batch = 0;
 o.flush_period = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.batch = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " batches of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the deterministic B&B on the knapsack against the dynamic program, with
// every selection: with 4 nodes per round its outcome (the value, the node,
// the bound and the counters) is the same with 1 to 8 workers and at every
// repetition, also when a node limit stops it, and with one node per worker
// it is the optimum in both senses; its leaves and its limit cases

template< class Executor >
static void test_bb_det( const std::string & )
{
 const std::string name = "deterministic B&B on " +
                          executor_name< Executor >();
 using Out = std::tuple< int , double , double , std::uint64_t ,
                         std::uint64_t , std::uint64_t , std::uint64_t ,
                         std::uint64_t , std::vector< char > >;
 const auto outcome = []( const auto & s ) {
  return( Out{ int( s.status ) , s.value , s.bound , s.stats.evaluated ,
               s.stats.generated , s.stats.pruned_generated ,
               s.stats.pruned_extracted , s.stats.pruned_evaluated ,
               s.node ? s.node->x : std::vector< char >() } );
  };
 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( DeterministicBranchAndBound< std::decay_t< decltype( q ) > ,
                                       Executor >( q , x ).solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( std::size_t limit : { std::size_t( 0 ) , std::size_t( 7 ) } ) {
    ParallelBBOptions o;
    o.selection = sel;
    o.round = 4;
    o.node_limit = limit;
    std::optional< Out > ref;
    bool ok = true , same = true;
    for( unsigned p : { 1u , 2u , 3u , 8u } )
     for( int rep = 0 ; rep < 2 ; ++rep ) {
      o.threads = p;
      const auto s = run( kp , o );
      if( limit )
       ok = ok && ( s.status == BBStatus::node_limit ) &&
            ( s.stats.evaluated <= limit ) && ( s.bound >= opt - 1e-9 );
      else
       ok = ok && ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
            s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 );
      const Out out = outcome( s );
      if( ! ref )
       ref = out;
      else
       same = same && ( out == *ref );
      }
    const std::string t = name + tag + " selection " +
                          std::to_string( int( sel ) ) +
                          ( limit ? " node limit" : "" );
    check( ok , t );
    check( same , t + ", the same with 1 to 8 workers" );
    }

  ParallelBBOptions o;
  o.threads = 3;
  const auto m = run( MinKnapsack{ kp } , o );
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         name + " minimization, one node per worker" + tag );
  check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                               std::optional< Leaves > seeds ) {
   DeterministicBranchAndBound< std::decay_t< decltype( q ) > , Executor >
    b( q , x );
   return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
   } , name + tag );
  }

 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.round = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.round = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " rounds of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the SelfSplit B&B on the knapsack against the dynamic program, with 1 to
// 8 workers, every selection and a split of 1 and 4 nodes per worker: its
// outcome is the same at every repetition, also when a node limit stops it,
// each worker evaluates its nodes by its own evaluator, as many as the
// counters say, and the values given to on_incumbent improve each time and
// end at the optimum; in both senses, its leaves and its limit cases

template< class Executor >
static void test_bb_split( const std::string & )
{
 const std::string name = "SelfSplit B&B on " + executor_name< Executor >();
 using Out = std::tuple< int , double , double , std::uint64_t ,
                         std::uint64_t , std::vector< std::uint64_t > ,
                         std::vector< char > >;
 const auto outcome = []( const auto & s ) {
  return( Out{ int( s.status ) , s.value , s.bound , s.stats.evaluated ,
               s.stats.generated , s.stats.evaluated_by ,
               s.node ? s.node->x : std::vector< char >() } );
  };
 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( SelfSplitBranchAndBound< std::decay_t< decltype( q ) > ,
                                   Executor >( q , x ).solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( std::size_t split : { std::size_t( 1 ) , std::size_t( 4 ) } )
    for( std::size_t limit : { std::size_t( 0 ) , std::size_t( 7 ) } ) {
     bool ok = true , same = true , own = true;
     for( unsigned p : { 1u , 2u , 3u , 8u } ) {
      ParallelBBOptions o;
      o.selection = sel;
      o.split = split;
      o.node_limit = limit;
      o.threads = p;
      std::optional< Out > ref;
      for( int rep = 0 ; rep < 2 ; ++rep ) {
       std::vector< std::atomic< std::uint64_t > > count( p );
       const auto s = run( CountKnapsack{ kp , count } , o );
       if( limit )
        ok = ok && ( s.status == BBStatus::node_limit ) &&
             ( s.stats.evaluated <= limit ) && ( s.bound >= opt - 1e-9 );
       else
        ok = ok && ( s.status == BBStatus::optimal ) &&
             ( s.value == opt ) && s.node && kp.check( *s.node , s.value ) &&
             ( s.gap <= 1e-9 );
       own = own && ( s.stats.evaluated_by.size() == p );
       for( unsigned j = 0 ; own && ( j < p ) ; ++j )
        own = s.stats.evaluated_by[ j ] == count[ j ].load();
       const Out out = outcome( s );
       if( ! ref )
        ref = out;
       else
        same = same && ( out == *ref );
       }
      }
     const std::string t = name + tag + " selection " +
                           std::to_string( int( sel ) ) + " split " +
                           std::to_string( split ) +
                           ( limit ? " node limit" : "" );
     check( ok , t );
     check( same , t + ", the same at every repetition" );
     check( own , t + ", the evaluations of each worker" );
     }

  ParallelBBOptions o;
  o.threads = 3;
  std::vector< double > seen;
  std::mutex m;
  o.on_incumbent = [ & ]( double v ) {
   std::lock_guard< std::mutex > guard( m );
   seen.push_back( v );
   };
  const auto s = run( kp , o );
  bool up = ( ! seen.empty() ) && ( seen.back() == opt ) &&
            ( s.value == opt );
  for( std::size_t j = 1 ; j < seen.size() ; ++j )
   up = up && ( seen[ j ] > seen[ j - 1 ] );
  check( up , name + " incumbents reported improving" + tag );

  o = ParallelBBOptions();
  o.threads = 3;
  const auto mn = run( MinKnapsack{ kp } , o );
  check( ( mn.status == BBStatus::optimal ) && ( mn.value == - opt ) &&
         mn.node && kp.check( *mn.node , - mn.value ) ,
         name + " minimization" + tag );
  check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                               std::optional< Leaves > seeds ) {
   SelfSplitBranchAndBound< std::decay_t< decltype( q ) > , Executor >
    b( q , x );
   return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
   } , name + tag );
  }

 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.split = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.split = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " split of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the B&B chosen by its scheme: each name gives back its scheme, an unknown
// name is refused, and every scheme solves the knapsack in both senses, from
// the root and from the leaves of a search, with 1 and 3 workers

static void test_bb_scheme( const std::string & )
{
 bool names = true;
 for( BBScheme s : { BBScheme::sequential , BBScheme::shared ,
                     BBScheme::distributed , BBScheme::deterministic ,
                     BBScheme::selfsplit , BBScheme::stealing } )
  names = names && ( bb_scheme( bb_scheme_name( s ) ) == s );
 bool refused = false;
 try {
  ( void ) bb_scheme( "nonesuch" );
  }
 catch( const std::invalid_argument & ) {
  refused = true;
  }
 check( names && refused , "B&B schemes by name" );

 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( std::string_view name : { "sequential" , "shared" , "distributed" ,
                                 "deterministic" , "selfsplit" ,
                                 "stealing" } )
   for( unsigned p : { 1u , 3u } ) {
    ParallelBBOptions o;
    o.threads = p;
    const BBScheme s = bb_scheme( name );
    const auto r = parallel_branch_and_bound( kp , s , o );
    const auto m = parallel_branch_and_bound( MinKnapsack{ kp } , s , o );
    const std::string t = "B&B scheme " + std::string( name ) + tag +
                          " p " + std::to_string( p );
    check( ( r.status == BBStatus::optimal ) && ( r.value == opt ) &&
           r.node && kp.check( *r.node , r.value ) &&
           ( m.status == BBStatus::optimal ) && ( m.value == - opt ) , t );
    check_frontier( kp , o , [ & ]( const auto & q ,
                                    const ParallelBBOptions & x ,
                                    std::optional< Leaves > seeds ) {
     return( parallel_branch_and_bound( q , s , x , std::move( seeds ) ) );
     } , t );
    }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a rational number of an exact arithmetic, as a B&B value: num / den with
// den > 0 and no common factor, or an infinity, den = 0 and num = +-1,
// whose products fit a long long for the values of the tests; it is
// trivially copyable and 16 bytes wide, so its incumbent goes under a lock
// where a std::atomic of that width is not lock-free (e.g., on x86-64)

struct Rational {

 long long num = 0;                             // the numerator
 long long den = 1;                             // the denominator

 Rational( void ) = default;

 Rational( long long n , long long d = 1 ) : num( n ) , den( d ) {
  if( den < 0 ) {
   num = - num;
   den = - den;
   }
  const long long g = std::gcd( num , den );
  if( g > 1 ) {
   num /= g;
   den /= g;
   }
  }

 static Rational infinite( int sign ) {
  Rational r;
  r.num = sign;
  r.den = 0;
  return( r );
  }

 bool operator==( const Rational & o ) const {
  return( ( num == o.num ) && ( den == o.den ) );
  }

 auto operator<=>( const Rational & o ) const {
  if( ! den || ! o.den )                       // an infinity on either side
   return( ( den ? 0 : num ) <=> ( o.den ? 0 : o.num ) );
  return( ( num * o.den ) <=> ( o.num * den ) );
  }

 Rational operator-( void ) const {
  Rational r = *this;
  r.num = - r.num;
  return( r );
  }

 Rational operator-( const Rational & o ) const {
  return( Rational( num * o.den - o.num * den , den * o.den ) );
  }

 };  // end( struct( Rational ) )

// the traits of Rational: an exact tolerance, the relative one taken as a
// fraction with denominator 10^6, and the test of a dive in long double

template<>
struct PBeFS_di_unipi_it::BBValueTraits< Rational > {
 static Rational infinity( void ) { return( Rational::infinite( 1 ) ); }
 static Rational none( void ) { return( infinity() ); }
 static bool finite( const Rational & v ) { return( v.den != 0 ); }
 static Rational tolerance( const Rational & z , const Rational & a ,
                            double r ) {
  const Rational t( ( z.num < 0 ? - z.num : z.num ) *
                    (long long)( r * 1e6 ) , z.den * 1000000 );
  return( std::max( a , t ) );
  }
 static bool within( const Rational & b , const Rational & l ,
                     const Rational & z , double g ) {
  if( ! finite( l ) )
   return( false );
  using L = long double;
  const auto x = []( const Rational & v ) { return( L( v.num ) / v.den ); };
  return( x( b ) <= x( l ) + L( g ) * ( x( z ) - x( l ) ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack with exact values: the fractional bound of Dantzig as a
// Rational, and the profits as integers

struct ExactKnapsack {

 using node_type = Knapsack::node_type;
 using value_type = Rational;

 static constexpr Sense sense = Sense::maximize;

 const Knapsack & kp;

 node_type root( void ) const { return( kp.root() ); }

 BasicEvaluation< Rational > evaluate( node_type & n ) const {
  BasicEvaluation< Rational > e;
  if( n.weight > kp.cap ) {
   e.infeasible = true;
   return( e );
   }
  long r = kp.cap - n.weight , b = n.profit;
  std::size_t j = n.k;
  for( ; ( j < kp.p.size() ) && ( kp.w[ j ] <= r ) ; ++j ) {
   r -= kp.w[ j ];
   b += kp.p[ j ];
   }
  e.bound = j < kp.p.size() ? Rational( b * kp.w[ j ] + kp.p[ j ] * r ,
                                        kp.w[ j ] )
                            : Rational( b );
  e.feasible = Rational( kp.greedy( n ) );
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  kp.branch( n , [ & ]( node_type && c , double ) {
   f( std::move( c ) , BBValueTraits< Rational >::infinity() );
   } );
  }

 std::uint64_t hash( const node_type & n ) const { return( kp.hash( n ) ); }

 void serialize( const node_type & n , std::vector< char > & out ) const {
  kp.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const {
  return( kp.deserialize( in ) );
  }

 };  // end( struct( ExactKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack kp with values of type long

static BasicKnapsack< long > integer_knapsack( const Knapsack & kp )
{
 BasicKnapsack< long > ki;
 ki.p = kp.p;
 ki.w = kp.w;
 ki.cap = kp.cap;
 ki.explicit_children = kp.explicit_children;
 return( ki );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the B&B with values other than double, by every scheme of the core with
// 1 and 3 workers: the knapsack with integer values (bounds rounded down),
// in both selections that order by the bound and with a node limit and an
// initial incumbent, and with the exact values of Rational (by the farm
// too); the optima are those of the dynamic program,
// and the integer bounds evaluate no more nodes than the double ones in
// the sequential best-first search, summed over the instances

static void test_bb_values( const std::string & )
{
 static_assert( detail::bb_lock_free< long >::value );
 std::uint64_t ev_double = 0 , ev_long = 0;
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const long opt = kp.optimum();
  const std::string tag = " knapsack " + std::to_string( i++ );
  const BasicKnapsack< long > ki = integer_knapsack( kp );
  const ExactKnapsack ke{ kp };

  ev_double += branch_and_bound( kp ).stats.evaluated;
  const auto s = branch_and_bound( ki );
  ev_long += s.stats.evaluated;
  check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
         s.node && ki.check( *s.node , s.value ) && ( s.gap == 0 ) &&
         ( s.bound == opt ) , "B&B long" + tag );

  for( Selection sel : { Selection::depth_first , Selection::plunging ,
                         Selection::best_estimate } ) {
   BasicBBOptions< long > o;
   o.selection = sel;
   const auto t = branch_and_bound( ki , o );
   check( ( t.status == BBStatus::optimal ) && ( t.value == opt ) ,
          "B&B long" + tag + " selection " + std::to_string( int( sel ) ) );
   }

  BasicBBOptions< long > o;
  o.node_limit = 3;
  const auto l = branch_and_bound( ki , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt ) &&
         ( l.stats.evaluated == 3 ) , "B&B long node limit" + tag );

  o = BasicBBOptions< long >();
  o.incumbent = opt;
  long seen = 0;
  o.on_incumbent = [ & ]( long v ) { seen = v; };
  const auto z = branch_and_bound( ki , o );
  check( ( z.status == BBStatus::optimal ) && ( z.value == opt ) &&
         ! z.node && ( seen == 0 ) , "B&B long with the optimum given" +
         tag );

  const auto x = branch_and_bound( ke );
  check( ( x.status == BBStatus::optimal ) && ( x.value == Rational( opt ) ) &&
         x.node && kp.check( *x.node , double( opt ) ) &&
         ( x.gap == Rational( 0 ) ) , "B&B Rational" + tag );

  for( std::string_view name : { "shared" , "distributed" ,
                                 "deterministic" , "selfsplit" ,
                                 "stealing" } )
   for( unsigned p : { 1u , 3u } ) {
    const BBScheme sch = bb_scheme( name );
    const std::string t = std::string( name ) + tag + " p " +
                          std::to_string( p );
    BasicParallelBBOptions< long > ol;
    ol.threads = p;
    const auto r = parallel_branch_and_bound( ki , sch , ol );
    check( ( r.status == BBStatus::optimal ) && ( r.value == opt ) &&
           r.node && ki.check( *r.node , r.value ) ,
           "B&B long scheme " + t );
    BasicParallelBBOptions< Rational > oq;
    oq.threads = p;
    const auto q = parallel_branch_and_bound( ke , sch , oq );
    check( ( q.status == BBStatus::optimal ) &&
           ( q.value == Rational( opt ) ) && q.node &&
           kp.check( *q.node , double( opt ) ) ,
           "B&B Rational scheme " + t );
    const auto fl = FarmBranchAndBound< BasicKnapsack< long > >( ki , ol )
                    .solve();
    const auto fq = FarmBranchAndBound< ExactKnapsack >( ke , oq ).solve();
    check( ( fl.status == BBStatus::optimal ) && ( fl.value == opt ) &&
           ( fq.status == BBStatus::optimal ) &&
           ( fq.value == Rational( opt ) ) , "B&B farm " + t );
    }
  }
 check( ev_long <= ev_double , "B&B long evaluates no more nodes (" +
        std::to_string( ev_long ) + " against " +
        std::to_string( ev_double ) + ")" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the incumbent given to the workers by cutoff( z ): the sequential, the
// shared, the distributed, the stealing and the deterministic B&B find the
// optimum with workers that fence the nodes by it, which they do (some
// evaluated node is infeasible, which none is in the knapsack alone), z
// never exceeds the optimum, and the outcome of the deterministic B&B is
// the same with 1 to 8 workers

static void test_bb_cutoff( const std::string & )
{
 std::atomic< bool > wrong = false;
 bool fenced = false;
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  const CutoffKnapsack ck{ kp , wrong , opt };
  const auto good = [ & ]( const auto & s ) {
   fenced = fenced || ( s.stats.infeasible > 0 );
   return( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           s.node && kp.check( *s.node , s.value ) );
   };
  BBOptions so;
  check( good( branch_and_bound( ck , so ) ) , "B&B by cutoff" + tag );
  ParallelBBOptions o;
  o.threads = 4;
  check( good( SharedBranchAndBound< CutoffKnapsack >( ck , o ).solve() ) ,
         "shared B&B by cutoff" + tag );
  check( good( DistributedBranchAndBound< CutoffKnapsack >( ck , o )
               .solve() ) , "distributed B&B by cutoff" + tag );
  check( good( StealingBranchAndBound< CutoffKnapsack >( ck , o ).solve() ) ,
         "stealing B&B by cutoff" + tag );
  o.round = 4;
  std::optional< std::tuple< double , double , std::uint64_t ,
                             std::uint64_t > > ref;
  bool ok = true , same = true;
  for( unsigned p : { 1u , 2u , 3u , 8u } ) {
   o.threads = p;
   const auto s = DeterministicBranchAndBound< CutoffKnapsack >( ck , o )
                   .solve();
   ok = ok && good( s );
   const auto out = std::make_tuple( s.value , s.bound , s.stats.evaluated ,
                                     s.stats.infeasible );
   if( ! ref )
    ref = out;
   else
    same = same && ( out == *ref );
   }
  check( ok , "deterministic B&B by cutoff" + tag );
  check( same , "deterministic B&B by cutoff" + tag +
                ", the same with 1 to 8 workers" );
  }
 check( fenced , "B&B by cutoff: some node fenced by the incumbent" );
 check( ! wrong , "B&B by cutoff: no incumbent above the optimum" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef PBEFS_HAS_MPI

// the B&B over MPI on the knapsack against the dynamic program, with every
// selection, in both senses (the minimization has no hash( n ), and hence
// random owners), with a relative tolerance and with a node limit, its
// leaves, its limit cases, and with integer and rational values, which
// travel as their bytes

static void test_bb_mpi( const std::string & )
{
 const auto solve = []( const auto & q , const ParallelBBOptions & o ,
                        std::optional< Leaves > seeds ) {
  MPIBranchAndBoundSearch< std::decay_t< decltype( q ) > > b( q , o );
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   ParallelBBOptions o;
   o.selection = sel;
   o.batch = 4;
   const auto s = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
   check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
          s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
          "MPI B&B" + tag + " selection " + std::to_string( int( sel ) ) );
   }

  ParallelBBOptions o;
  const MinKnapsack mk{ kp };
  const auto m = MPIBranchAndBoundSearch< MinKnapsack >( mk , o ).solve();
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         "MPI B&B minimization" + tag );

  o.rel_gap = 0.05;
  const auto r = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) ,
         "MPI B&B within a relative gap of 5%" + tag );

  o = ParallelBBOptions();
  o.node_limit = 3;
  const auto l = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) ,
         "MPI B&B node limit" + tag );

  o = ParallelBBOptions();
  o.batch = 4;
  check_frontier( kp , o , solve , "MPI B&B" + tag , true );

  const BasicKnapsack< long > ki = integer_knapsack( kp );
  BasicParallelBBOptions< long > ol;
  ol.batch = 4;
  const auto si = MPIBranchAndBoundSearch< BasicKnapsack< long > >( ki , ol )
                  .solve();
  check( ( si.status == BBStatus::optimal ) && ( si.value == long( opt ) ) &&
         si.node && ki.check( *si.node , si.value ) && ( si.gap == 0 ) ,
         "MPI B&B long" + tag );
  const ExactKnapsack ke{ kp };
  BasicParallelBBOptions< Rational > oq;
  oq.batch = 4;
  const auto sq = MPIBranchAndBoundSearch< ExactKnapsack >( ke , oq ).solve();
  check( ( sq.status == BBStatus::optimal ) &&
         ( sq.value == Rational( (long long)( opt ) ) ) && sq.node &&
         kp.check( *sq.node , opt ) , "MPI B&B Rational" + tag );
  }

 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( MPIBranchAndBoundSearch< std::decay_t< decltype( q ) > >( q , x )
          .solve() );
  };
 ParallelBBOptions o;
 check_bb_edges( knapsacks().front() , o , run , "MPI B&B" );
 o.batch = 0;
 o.flush_period = 0;
 check_bb_edges( knapsacks().front() , o , run , "MPI B&B parameters 0" );
 }

#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the tester, each with the runtime it needs ("-" if none):
// test/tests.txt lists the same ones [see test_registry()]

struct Case {
 const char * name;                          // the name on the command line
 const char * needs;                         // the runtime it needs
 void ( * run )( const std::string & );      // the test, given the data
 };

static const std::vector< Case > cases = {
 { "astar" , "-" , test_astar } ,
 { "open" , "-" , []( const std::string & ) { test_open(); } } ,
 { "grid" , "-" , []( const std::string & ) { test_grid(); } } ,
 { "idastar" , "-" , test_idastar } ,
 { "tiles24" , "-" , test_tiles24 } ,
 { "graph" , "-" , []( const std::string & ) { test_graph(); } } ,
 { "limits" , "-" , test_limits } ,
 { "reopen" , "-" , []( const std::string & ) { test_reopen(); } } ,
 { "pdb" , "-" , test_pdb } ,
 { "pdb-search" , "-" , test_pdb_search } ,
 { "hda-thread" , "-" , test_hda< ADS::ThreadPool > } ,
 { "spa-thread" , "-" , test_spa< ADS::ThreadPool > } ,
 { "pbnf-thread" , "-" , test_pbnf< ADS::ThreadPool > } ,
 { "pase-thread" , "-" , test_pase< ADS::ThreadPool > } ,
 { "bb" , "-" , test_bb } ,
 { "bb-shared-thread" , "-" , test_bb_parallel< 0 , ADS::ThreadPool > } ,
 { "bb-dist-thread" , "-" , test_bb_parallel< 1 , ADS::ThreadPool > } ,
 { "bb-steal-thread" , "-" , test_bb_parallel< 3 , ADS::ThreadPool > } ,
 { "bb-det-thread" , "-" , test_bb_det< ADS::ThreadPool > } ,
 { "bb-split-thread" , "-" , test_bb_split< ADS::ThreadPool > } ,
 { "bb-farm-thread" , "-" , test_bb_parallel< 2 , ADS::ThreadPool > } ,
 { "farm-thread" , "-" , test_farm< ADS::ThreadPool > } ,
 { "bb-scheme" , "-" , test_bb_scheme } ,
 { "bb-cutoff" , "-" , test_bb_cutoff } ,
 { "bb-values" , "-" , test_bb_values } ,
 #ifdef PBEFS_HAS_OPENMP
  { "hda-openmp" , "openmp" , test_hda< ADS::OpenMPExecutor > } ,
  { "spa-openmp" , "openmp" , test_spa< ADS::OpenMPExecutor > } ,
  { "pbnf-openmp" , "openmp" , test_pbnf< ADS::OpenMPExecutor > } ,
  { "bb-shared-openmp" , "openmp" ,
    test_bb_parallel< 0 , ADS::OpenMPExecutor > } ,
  { "bb-dist-openmp" , "openmp" ,
    test_bb_parallel< 1 , ADS::OpenMPExecutor > } ,
  { "bb-steal-openmp" , "openmp" ,
    test_bb_parallel< 3 , ADS::OpenMPExecutor > } ,
  { "bb-det-openmp" , "openmp" , test_bb_det< ADS::OpenMPExecutor > } ,
  { "bb-split-openmp" , "openmp" , test_bb_split< ADS::OpenMPExecutor > } ,
  { "bb-farm-openmp" , "openmp" ,
    test_bb_parallel< 2 , ADS::OpenMPExecutor > } ,
  { "farm-openmp" , "openmp" , test_farm< ADS::OpenMPExecutor > } ,
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  { "hda-fastflow" , "fastflow" , test_hda< ADS::FastFlowExecutor > } ,
  { "farm-fastflow" , "fastflow" , test_farm< ADS::FastFlowExecutor > } ,
  { "bb-shared-fastflow" , "fastflow" ,
    test_bb_parallel< 0 , ADS::FastFlowExecutor > } ,
  { "bb-dist-fastflow" , "fastflow" ,
    test_bb_parallel< 1 , ADS::FastFlowExecutor > } ,
  { "bb-steal-fastflow" , "fastflow" ,
    test_bb_parallel< 3 , ADS::FastFlowExecutor > } ,
  { "bb-farm-fastflow" , "fastflow" ,
    test_bb_parallel< 2 , ADS::FastFlowExecutor > } ,
  { "bb-det-fastflow" , "fastflow" , test_bb_det< ADS::FastFlowExecutor > } ,
  { "bb-split-fastflow" , "fastflow" ,
    test_bb_split< ADS::FastFlowExecutor > } ,
 #endif
 #ifdef PBEFS_HAS_MPI
  { "mpi" , "mpi" , test_mpi } ,
  { "bb-mpi" , "mpi" , test_bb_mpi } ,
 #endif
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the list (test/tests.txt) whose runtime is built are those
// of the tester, and vice versa

static void test_registry( const std::string & list )
{
 const auto built = []( const std::string & needs ) {
  if( needs == "-" )
   return( true );
  #ifdef PBEFS_HAS_OPENMP
   if( needs == "openmp" )
    return( true );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
   if( needs == "fastflow" )
    return( true );
  #endif
  #ifdef PBEFS_HAS_MPI
   if( needs == "mpi" )
    return( true );
  #endif
  return( false );
  };

 std::ifstream in( list );
 check( bool( in ) , "the list " + list + " is read" );
 std::set< std::string > listed;
 for( std::string line ; std::getline( in , line ) ; ) {
  std::istringstream is( line );
  std::string name , needs;
  if( ( is >> name >> needs ) && ( name[ 0 ] != '#' ) && built( needs ) )
   listed.insert( name );
  }
 std::set< std::string > tester;
 for( const Case & c : cases )
  tester.insert( c.name );
 for( const auto & n : listed )
  check( tester.count( n ) == 1 , "the case " + n + " is in the tester" );
 for( const auto & n : tester )
  check( listed.count( n ) == 1 , "the case " + n + " is in the list" );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 if( argc < 2 ) {
  std::cerr << "usage: " << argv[ 0 ] << " <case> [<data dir>]" << std::endl;
  return( 2 );
  }
 const std::string what = argv[ 1 ];
 const std::string data = argc > 2 ? argv[ 2 ] : "../data";
 const auto c = std::find_if( cases.begin() , cases.end() ,
                              [ & ]( const Case & x ) {
                               return( what == x.name );
                               } );
 if( ( what != "registry" ) && ( c == cases.end() ) ) {
  std::cerr << "unknown case " << what << std::endl;
  return( 2 );
  }
 [[maybe_unused]] const bool mpi = ( c != cases.end() ) &&
                                   ( std::string( c->needs ) == "mpi" );

 #ifdef PBEFS_HAS_MPI
  int rank = 0;
  if( mpi ) {
   MPI_Init( & argc , & argv );
   MPI_Comm_rank( MPI_COMM_WORLD , & rank );
   quiet = rank != 0;
   }
 #endif

 try {
  if( what == "registry" )
   test_registry( data );
  else
   c->run( data );
  }
 catch( std::exception & e ) {
  std::cout << "FAILED: exception " << e.what() << std::endl;
  ++failures;
  }

 #ifdef PBEFS_HAS_MPI
  if( mpi ) {
   int all = 0;
   MPI_Allreduce( & failures , & all , 1 , MPI_INT , MPI_SUM ,
                  MPI_COMM_WORLD );
   failures = all;
   MPI_Finalize();
   }
 #endif

 if( ! quiet )
  std::cout << ( failures ? "FAILED: " : "passed: " ) << what << " ("
            << failures << " failures)" << std::endl;
 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File test.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
