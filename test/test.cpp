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
 * instances of Korf (1985), Dijkstra's algorithm on the grids and on the
 * boards, the dynamic program on the lattice of the alignments, the one of
 * Held and Karp on the cities, and the distances on the map of Romania of
 * Russell and Norvig.
 *
 *     test astar | grid | graph | limits | reopen
 *     test msa | board | files
 *     test pdb | pdb-search
 *     test hda-thread | hda-openmp | hda-fastflow
 *     test farm-thread | farm-openmp | farm-fastflow
 *     test spa-thread | spa-openmp
 *     test mqa-thread | mqa-openmp
 *     test pbnf-thread | pbnf-openmp
 *     test pase-thread
 *     test pida-thread | pida-openmp | pida-fastflow
 *     test permutation | abstraction | pancake | topspin | hanoi | dimacs
 *     test rubik | rubik-pdb
 *     test bb | bb-shared-thread | bb-dist-thread | bb-steal-thread
 *     test bb-mq-thread | bb-mq-openmp
 *     test bb-shared-openmp | bb-dist-openmp | bb-steal-openmp
 *     test bb-shared-fastflow | bb-dist-fastflow | bb-steal-fastflow
 *     test bb-farm-thread | bb-farm-openmp | bb-farm-fastflow
 *     test bb-det-thread | bb-det-openmp | bb-det-fastflow
 *     test bb-split-thread | bb-split-openmp | bb-split-fastflow
 *     test bb-scheme | bb-cutoff | bb-values
 *     test bb-skel-thread | bb-skel-openmp | bb-skel-fastflow
 *     test bb-ivm-thread | bb-ivm-openmp | bb-ivm-fastflow
 *     test bb-modes | bb-race | bb-checkpoint | bb-domains
 *     test mpi | bb-mpi | pida-mpi | bb-hub-mpi   (under mpiexec)
 *     test registry <tests.txt>
 *
 * The cases are listed in test/tests.txt, which CMake and the makefile
 * read; the case registry checks that the list and the cases of the tester
 * are the same.
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

#include <cmath>

#include <cctype>

#include <compare>

#include <cstring>

#include <deque>

#include <filesystem>

#include <fstream>

#include <functional>

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

#include <thread>

#include <tuple>

#include <type_traits>

#include <unordered_map>

#include <utility>

#include <vector>

#include "PBeFSConfig.h"

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/AnytimeAStar.h"

#include "PBeFS/astar/FocalSearch.h"

#include "PBeFS/astar/ImprovedOptimisticSearch.h"

#include "PBeFS/astar/LPAStar.h"

#include "PBeFS/astar/OptimisticSearch.h"

#include "PBeFS/astar/PotentialSearch.h"

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "PBeFS/bnb/DeterministicBranchAndBound.h"

#include "PBeFS/bnb/SelfSplitBranchAndBound.h"

#include "PBeFS/bnb/FarmBranchAndBound.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/astar/AbstractionPDB.h"

#include "PBeFS/astar/DimacsGraph.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/ExplicitGraph.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/DDDAStar.h"

#include "PBeFS/astar/SPAStar.h"
#include "PBeFS/astar/MQAStar.h"

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/ParallelIDAStar.h"

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/FarmAStar.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/Permutation.h"

#include "PBeFS/astar/PermutationPDB.h"

#include "PBeFS/astar/RubiksCube.h"

#include "PBeFS/astar/SafePBNF.h"

#include "PBeFS/astar/SequenceAlignment.h"

#include "PBeFS/astar/SlowDomain.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/TopSpin.h"

#include "PBeFS/astar/TowersOfHanoi.h"

#include "PBeFS/astar/TravelingSalesman.h"

#ifdef PBEFS_HAS_OPENMP
 #include "ADS/parallel/OpenMPExecutor.h"
#endif

#ifdef PBEFS_HAS_FASTFLOW
 #include "ADS/parallel/FastFlowExecutor.h"
#endif

#ifdef PBEFS_HAS_MPI
 #include "PBeFS/astar/HDAStarMPI.h"
 #include "PBeFS/astar/ParallelIDAStarMPI.h"
 #include "PBeFS/bnb/MPIBranchAndBound.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- GLOBAL VARIABLES ----------------------------*/
/*--------------------------------------------------------------------------*/

int failures = 0;                     // the checks failed so far
bool quiet = false;                   // print only the failures

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// records the outcome of a check

void check( bool ok , const std::string & what )
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// 16 states of the 8-puzzle spread out, with their distances, in the order
// of their packed boards, so that every platform takes the same ones

static std::vector< std::pair< std::uint64_t , int > > eight_sample(
 const std::unordered_map< std::uint64_t , int > & dist )
{
 std::vector< std::pair< std::uint64_t , int > > all( dist.begin() ,
                                                      dist.end() );
 std::sort( all.begin() , all.end() );
 std::vector< std::pair< std::uint64_t , int > > v;
 for( std::size_t i = 0 ; i < all.size() ; i += 12007 )
  v.push_back( all[ i ] );
 return( v );
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
 const auto mu = MQAStar< Grid< 4 > >( g , po ).solve();
 check( mu.status == Status::unsolvable ,
        "MultiQueue A* proves the wall unsolvable" );

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
// A* with delayed duplicate detection on the disk against the exact
// distances of the 8-puzzle (with one file per bucket, a prime number of
// them and the default), A* on the 15-puzzle, the grids and Romania; the
// directory of its files removed at the end, the unsolvable wall, the node
// limit, and the errors of an inconsistent heuristic and of no partition

static void test_ddd( const std::string & data )
{
 ADS::TempDirectory dir;                       // where the runs go
 DDDOptions o;
 o.directory = dir.path();

 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 for( std::size_t parts :
      { std::size_t( 1 ) , std::size_t( 7 ) , std::size_t( 64 ) } ) {
  o.partitions = parts;
  std::size_t k = 0;
  int wrong = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 3628 ) != 0 )                   // 50 states, spread out
    continue;
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const auto r = ddd_astar( d , {} , o );
   if( ! ( r.solved() && ( r.cost == dst ) ) ) {
    ++wrong;
    check_solution( d , r , dst , "DDD-A* 8-puzzle" );
    }
   }
  check( wrong == 0 , "DDD-A* on 50 states of the 8-puzzle, " +
                       std::to_string( parts ) + " partitions" );
  }
 check( std::filesystem::is_empty( dir.path() ) ,
        "DDD-A* removes the directory of its files" );

 o.partitions = 16;
 o.buffer_bytes = 4096;                        // some records per write
 const auto inst = korf( data );
 {
  const int id = easy_korf[ 0 ];               // 2 * 10^5 expansions
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto r = ddd_astar( d , {} , o );
  check_solution( d , r , ki.cost , "DDD-A*, Korf #" + std::to_string( id ) );
  check( r.stats.expanded <= r.stats.stored , "DDD-A* expands a state at"
                                             " most once, Korf #" +
                                               std::to_string( id ) );
  }

 o = DDDOptions{};
 o.directory = dir.path();
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 7 , 5 , 0.2 , seed , 3 );
  const Grid< 4 > g4( m , 0 , 4 , 6 , 4 );
  const auto a4 = astar( g4 );
  const auto r4 = ddd_astar( g4 , {} , o );
  if( a4.solved() )
   check_solution( g4 , r4 , a4.cost ,
                   "DDD-A* grid 4, seed " + std::to_string( seed ) );
  else
   check( r4.status == Status::unsolvable ,
          "DDD-A* unsolvable grid, seed " + std::to_string( seed ) );
  }
 auto g = romania( Arad );
 check_solution( g , ddd_astar( g , {} , o ) , 418 ,
                 "DDD-A* Romania, Arad to Bucharest" );

 GridMap w( 50 , 50 );
 for( int y = 0 ; y < 50 ; ++y )
  w.set_cost( 25 , y , 0 );
 const Grid< 4 > wall( w , 0 , 0 , 49 , 49 );
 const auto u = ddd_astar( wall , {} , o );
 check( ( u.status == Status::unsolvable ) &&
          ( u.stats.expanded == 25 * 50 ) && ( u.stats.stored == 25 * 50 ) ,
        "DDD-A* proves the wall unsolvable"
        " after the 1250 cells it reaches" );
 check( u.stats.generated + 1 == u.stats.stored + u.stats.duplicates ,
        "DDD-A* stores each record of the wall once or drops it as a copy" );

 SearchOptions lim;
 lim.node_limit = 10000;
 const FifteenPuzzle hard( inst[ 0 ].board );
 const auto l = ddd_astar( hard , lim , o );
 check( ( l.status == Status::node_limit ) && ( l.stats.stored >= 10000 ) &&
          l.path.empty() ,
        "DDD-A* stops at the node limit" );

 ExplicitGraph< int > bad( 5 , 0 , 4 );         // that of test_reopen
 bad.add_edge( 0 , 1 , 1 , true );
 bad.add_edge( 1 , 2 , 1 , true );
 bad.add_edge( 0 , 2 , 3 , true );
 bad.add_edge( 2 , 3 , 1 , true );
 bad.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };
 for( std::uint32_t v = 0 ; v < 5 ; ++v )
  bad.set_heuristic( v , h[ v ] );
 bool threw = false;
 try {
  ( void ) ddd_astar( bad , {} , o );
  }
 catch( const std::logic_error & ) {
  threw = true;
  }
 check( threw , "DDD-A* refuses an inconsistent heuristic" );
 threw = false;
 try {
  DDDOptions none = o;
  none.partitions = 0;
  ( void ) ddd_astar( g , {} , none );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "DDD-A* refuses no partition" );
 check( std::filesystem::is_empty( dir.path() ) ,
        "DDD-A* leaves no file, also when it stops or throws" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the hash of the states of a Domain for the containers of the standard
// library

template< class D >
struct StateHash {
 const D * d;
 std::size_t operator()( const state_t< D > & s ) const {
  return( std::size_t( d->hash( s ) ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cost g*( s ) of a cheapest path from the initial state to each state
// that it reaches, by Dijkstra's algorithm on the edges of the Domain

template< class D >
static std::unordered_map< state_t< D > , cost_t< D > , StateHash< D > >
 all_distances( const D & d )
{
 using S = state_t< D >;
 using C = cost_t< D >;
 std::unordered_map< S , C , StateHash< D > > g( 1024 , StateHash< D >{ & d }
                                                 );
 using Entry = std::pair< C , S >;
 const auto later = []( const Entry & a , const Entry & b ) {
  return( a.first > b.first );
  };
 std::priority_queue< Entry , std::vector< Entry > , decltype( later ) >
  q( later );
 g[ d.initial_state() ] = 0;
 q.emplace( C( 0 ) , d.initial_state() );
 while( ! q.empty() ) {
  const C gs = q.top().first;
  const S s = q.top().second;
  q.pop();
  if( gs > g[ s ] )
   continue;
  d.for_each_successor( s , d.heuristic( s ) , [ & ]( const S & t , C c ,
                                                       C ) {
   const auto it = g.find( t );
   if( ( it == g.end() ) || ( gs + c < it->second ) ) {
    g[ t ] = gs + c;
    q.emplace( gs + c , t );
    }
   } );
  }
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks the heuristic of d on the states of g: h( t ) >= h( s ) - c on
// every edge ( s , t ) of cost c, h = 0 at the goals, and A* (sol) solved
// with the optimal cost, expanding no state twice (but for 10% of them with
// real costs), at least the states of f* = g* + h smaller than the optimal
// cost C* and at most those of f* not larger: the bounds that every
// tie-breaking of A* meets with a consistent heuristic

template< class D , class G , class S >
static void check_consistent_astar( const D & d , const G & g , const S & sol ,
                                    const std::string & what )
{
 using C = cost_t< D >;
 const auto eps = []( C x ) {
  if constexpr( std::is_floating_point_v< C > )
   return( 1e-9 * std::max( C( 1 ) , std::abs( x ) ) );
  else
   return( C( 0 ) );
  };
 C opt = infinite_cost< C >();
 std::size_t bad = 0;
 for( const auto & [ s , gs ] : g ) {
  const C hs = d.heuristic( s );
  if( d.is_goal( s ) ) {
   opt = std::min( opt , gs );
   if( hs != 0 )
    ++bad;
   }
  d.for_each_successor( s , hs , [ & ]( const state_t< D > & , C c ,
                                        C ht ) {
   if( ht + c < hs - eps( hs ) )
    ++bad;
   } );
  }
 check( bad == 0 , what + ": the heuristic is consistent and 0 at the goal"
        " on the " + std::to_string( g.size() ) + " reachable states" );
 check_solution( d , sol , opt , what + ": A* against Dijkstra" );
 std::size_t lt = 0 , le = 0;
 for( const auto & [ s , gs ] : g ) {
  const C f = gs + d.heuristic( s );
  lt += f < opt - eps( opt );
  le += f <= opt + eps( opt );
  }
 // with real costs, two paths of the same cost summed in another order may
 // differ in the last bit, and the cheaper one found later reopens a state
 const auto r = sol.stats.reexpanded , e = sol.stats.expanded - r;
 check( ( std::is_floating_point_v< C > ? 10 * r <= e : r == 0 ) &&
        ( lt <= e ) && ( e <= le ) ,
        what + ": A* expands " + std::to_string( e ) + " states (and " +
        std::to_string( r ) + " again), between " + std::to_string( lt ) +
        " (f* < C*) and " + std::to_string( le ) + " (f* <= C*)" );
 }

/*--------------------------------------------------------------------------*/
// the cost of an optimal alignment by the dynamic program on the whole
// lattice, written from the definition of the sum-of-pairs cost with no
// part of SequenceAlignment but its pair costs

static int alignment_dp( const std::vector< std::string > & seq ,
                         const SequenceAlignment<> & d )
{
 const int n = int( seq.size() );
 std::vector< std::size_t > w( n );
 std::size_t size = 1;
 for( int i = n ; i-- > 0 ; ) {
  w[ i ] = size;
  size *= seq[ i ].size() + 1;
  }
 std::vector< int > best( size , std::numeric_limits< int >::max() );
 best[ 0 ] = 0;
 std::vector< std::size_t > p( n );
 for( std::size_t k = 0 ; k < size ; ++k ) {
  std::size_t r = k;
  for( int i = 0 ; i < n ; ++i ) {
   p[ i ] = r / w[ i ];
   r %= w[ i ];
   }
  for( std::uint32_t m = 1 ; m < ( 1u << n ) ; ++m ) {
   bool ok = true;
   std::size_t from = k;
   for( int i = 0 ; i < n ; ++i )
    if( m & ( 1u << i ) ) {
     ok = ok && ( p[ i ] > 0 );
     from -= w[ i ];
     }
   if( ! ok || ( best[ from ] == std::numeric_limits< int >::max() ) )
    continue;
   int c = 0;
   for( int i = 0 ; i < n ; ++i )
    for( int j = i + 1 ; j < n ; ++j ) {
     const int a = m & ( 1u << i ) ? SequenceAlignment<>::symbol(
                                      seq[ i ][ p[ i ] - 1 ] ) : amino_gap;
     const int b = m & ( 1u << j ) ? SequenceAlignment<>::symbol(
                                      seq[ j ][ p[ j ] - 1 ] ) : amino_gap;
     c += d.pair_cost( a , b );
     }
   best[ k ] = std::min( best[ k ] , best[ from ] + c );
   }
  }
 return( best[ size - 1 ] );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the sequences of a FASTA file of the data

static std::vector< std::string > fasta( const std::string & file )
{
 std::ifstream in( file );
 if( ! in )
  throw( std::runtime_error( "cannot read " + file ) );
 return( SequenceAlignment<>::read_fasta( in ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// multiple sequence alignment: the reader of FASTA files, the costs of
// PAM250, A* against the dynamic program on the lattice, Dijkstra and the
// bounds on its expansions, the alignment of the path, IDA* and DDD-A*,
// the features of the Zobrist distributions and the graph of the nblocks

static void test_msa( const std::string & data )
{
 using D = SequenceAlignment<>;
 static_assert( IncrementalFeatureDomain< D > && AbstractFeatureDomain< D >
                && AbstractStateDomain< D > && NBlockDomain< D > &&
                DiskDomain< D > && ! PairwiseDomain< D > );

 std::istringstream fa( "junk\n>a b\nac-g\n d\n>c\n\nWW.\n" );
 const auto rf = D::read_fasta( fa );
 check( ( rf.size() == 2 ) && ( rf[ 0 ] == "ACGD" ) && ( rf[ 1 ] == "WW" ) ,
        "FASTA: names, lines, case, gaps and blank lines" );
 bool threw = false;
 try {
  std::istringstream none( "ACGT\n" );
  ( void ) D::read_fasta( none );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "FASTA without sequences refused" );
 int refused = 0;
 for( const auto & bad : std::vector< std::vector< std::string > >{
       { "ACD" } , { "AC" , "" } , { "AC" , "AOC" } , { "AC" , "A*" } ,
       std::vector< std::string >( 11 , "A" ) } )
  try {
   const D x( bad );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 5 , "a single sequence, an empty one, an unknown symbol,"
                       " a gap and 11 sequences are refused" );

 const D two( { "W" , "W" } );
 const auto sym = []( char a ) { return( D::symbol( a ) ); };
 check( ( two.pair_cost( sym( 'W' ) , sym( 'W' ) ) == 0 ) &&
        ( two.pair_cost( sym( 'C' ) , sym( 'C' ) ) == 5 ) &&
        ( two.pair_cost( sym( 'A' ) , amino_gap ) == 25 ) &&
        ( two.pair_cost( amino_gap , amino_gap ) == 0 ) &&
        ( two.pair_cost( sym( 'W' ) , sym( 'C' ) ) == 25 ) ,
        "the costs 17 - PAM250, 25 against a gap, 0 for two gaps" );
 bool symmetric = true;
 for( int a = 0 ; a < amino_count ; ++a )
  for( int b = 0 ; b < amino_count ; ++b )
   symmetric = symmetric && ( pam250[ a ][ b ] == pam250[ b ][ a ] );
 check( symmetric , "PAM250 is symmetric" );
 // a pair aligned with the linear gap penalty 8: the cost against the
 // largest score sum( S ) - 8 g, by the recurrence of Needleman and Wunsch
 {
  const std::string u = "HEAGAWGHEE" , v = "PAWHEAE";
  std::vector< std::vector< int > > best( u.size() + 1 ,
                                          std::vector< int >( v.size() + 1 ) );
  for( std::size_t i = 0 ; i <= u.size() ; ++i )
   for( std::size_t j = 0 ; j <= v.size() ; ++j )
    if( i + j == 0 )
     best[ i ][ j ] = 0;
    else {
     int b = std::numeric_limits< int >::min();
     if( i )
      b = std::max( b , best[ i - 1 ][ j ] - 8 );
     if( j )
      b = std::max( b , best[ i ][ j - 1 ] - 8 );
     if( i && j )
      b = std::max( b , best[ i - 1 ][ j - 1 ] +
                        pam250[ sym( u[ i - 1 ] ) ][ sym( v[ j - 1 ] ) ] );
     best[ i ][ j ] = b;
     }
  const D p( { u , v } , pam250 , 8 , 16 , 8 );
  const int score = best[ u.size() ][ v.size() ];
  check_solution( p , astar( p ) , 17 * int( u.size() + v.size() ) -
                  2 * score , "the gap penalty 8: twice the cost is 17 ( L1 "
                  "+ L2 ) less twice the score " + std::to_string( score ) );
  const D q( { u , v } );
  int score16 = 0;                     // the same with the penalty 16.5
  std::vector< std::vector< double > > b2( u.size() + 1 ,
                                           std::vector< double >(
                                            v.size() + 1 ) );
  for( std::size_t i = 0 ; i <= u.size() ; ++i )
   for( std::size_t j = 0 ; j <= v.size() ; ++j )
    if( i + j > 0 ) {
     double b = -1e18;
     if( i )
      b = std::max( b , b2[ i - 1 ][ j ] - 16.5 );
     if( j )
      b = std::max( b , b2[ i ][ j - 1 ] - 16.5 );
     if( i && j )
      b = std::max( b , b2[ i - 1 ][ j - 1 ] +
                        pam250[ sym( u[ i - 1 ] ) ][ sym( v[ j - 1 ] ) ] );
     b2[ i ][ j ] = b;
     }
  score16 = int( std::lround( 2 * b2[ u.size() ][ v.size() ] ) );
  check_solution( q , astar( q ) , ( 17 * int( u.size() + v.size() ) -
                                     score16 ) / 2 ,
                  "the gap costs of PAM250: the score with the penalty 16.5" );
  }
 std::ostringstream os;
 for( const auto & row : pam250 )
  for( int x : row )
   os << x << ' ';
 std::istringstream is( os.str() );
 check( D::read_scores( is ) == pam250 , "a matrix read back" );

 // the costs of the instances, against the dynamic program on the lattice
 // and the costs that the code of HDA* of Jinnai and Fukunaga finds
 struct Inst { const char * file; int cost; };
 const Inst inst[] = { { "test0" , 1741 } , { "test1" , 2665 } ,
                       { "test3" , 2022 } , { "test4" , 2312 } ,
                       { "test5" , 3990 } };
 for( const auto & [ name , cost ] : inst ) {
  const auto seq = fasta( data + "/msa/easy/" + name + ".tfa" );
  const D d( seq );
  const std::string what = std::string( "MSA " ) + name;
  const auto sol = astar( d );
  if( seq.size() <= 4 ) {      // test5: 483840 points of 127 moves each
   check( alignment_dp( seq , d ) == cost , what + ": the dynamic program "
          "on the lattice costs " + std::to_string( cost ) );
   check_consistent_astar( d , all_distances( d ) , sol , what );
   }
  else
   check_solution( d , sol , cost , what + ", the cost of the reference "
                   "code" );
  const auto rows = d.alignment( sol.path );
  bool same = rows.size() == seq.size();
  int columns = 0;
  for( std::size_t c = 0 ; same && ( c < rows[ 0 ].size() ) ; ++c ) {
   for( std::size_t i = 0 ; i < rows.size() ; ++i )
    for( std::size_t j = i + 1 ; j < rows.size() ; ++j )
     columns += d.pair_cost( rows[ i ][ c ] == '-' ? amino_gap
                                                   : sym( rows[ i ][ c ] ) ,
                             rows[ j ][ c ] == '-' ? amino_gap
                                                   : sym( rows[ j ][ c ] ) );
   }
  for( std::size_t i = 0 ; same && ( i < rows.size() ) ; ++i ) {
   std::string r = rows[ i ];
   r.erase( std::remove( r.begin() , r.end() , '-' ) , r.end() );
   same = ( r == seq[ i ] ) && ( rows[ i ].size() == rows[ 0 ].size() );
   }
  check( same && ( columns == cost ) , what + ": the rows of the alignment "
         "are the sequences and their columns cost " +
         std::to_string( columns ) );
  }

 // the instance of 7 sequences, IDA* and DDD-A*
 {
  const auto seq = fasta( data + "/msa/easy/test2.tfa" );
  const D d( seq );
  check_solution( d , astar( d ) , 6055 , "MSA test2, 7 sequences" );
  const auto seq0 = fasta( data + "/msa/easy/test0.tfa" );
  const D d0( seq0 );
  check_solution( d0 , idastar( d0 ) , 1741 , "IDA* MSA test0" );
  ADS::TempDirectory dir;
  DDDOptions o;
  o.directory = dir.path();
  check_solution( d0 , ddd_astar( d0 , {} , o ) , 1741 , "DDD-A* MSA test0" );
  }

 // the instances of BAliBASE read and their heuristics computed
 int read = 0;
 for( const auto & f : std::filesystem::directory_iterator(
                        data + "/msa/balibase" ) ) {
  const auto seq = fasta( f.path().string() );
  const D d( seq );
  read += ( seq.size() >= 4 ) && ( d.heuristic( d.initial_state() ) > 0 );
  }
 check( read == 15 , "the 15 instances of BAliBASE, 4 or 5 sequences each" );

 // the Zobrist keys of the successors from the key of the parent, and the
 // graph of the nblocks: symmetric, and holding every move between boxes
 const auto seq = fasta( data + "/msa/easy/test3.tfa" );
 const D d( seq , pam250 , 3 , 4 );
 const ZobristDistribution< D > z( d );
 const AbstractZobristDistribution< D > az( d );
 static_assert( IncrementalDistribution< decltype( z ) , D > &&
                IncrementalDistribution< decltype( az ) , D > );
 const auto g = all_distances( d );
 std::set< std::pair< std::size_t , std::size_t > > nb;
 for( std::size_t b = 0 ; b < d.num_nblocks() ; ++b )
  d.for_each_nblock_successor( b , [ & ]( std::size_t c ) {
   nb.emplace( b , c );
   } );
 int bad = 0 , moves = 0 , kz = 0 , kaz = 0 , ka = 0 , visited = 0;
 const AbstractStateDistribution< D > as( d );
 for( const auto & sg : g ) {
  const auto & s = sg.first;
  if( ( visited++ % 5 ) != 0 )
   continue;
  for( unsigned p : { 1u , 3u , 64u } ) {
   const auto oz = successor_owners< D >( z , s , z.key( s ) , p );
   const auto oaz = successor_owners< D >( az , s , az.key( s ) , p );
   d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
    if( ( oz( t ) != std::pair( z.owner( t , p ) , z.key( t ) ) ) ||
        ( oaz( t ) != std::pair( az.owner( t , p ) , az.key( t ) ) ) )
     ++bad;
    } );
   }
  d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
   ++moves;
   kz += z.owner( t , 64 ) == z.owner( s , 64 );
   kaz += az.owner( t , 64 ) == az.owner( s , 64 );
   ka += as.owner( t , 64 ) == as.owner( s , 64 );
   const auto a = d.nblock( s ) , b = d.nblock( t );
   if( ( a != b ) && ( ! nb.count( { a , b } ) || ! nb.count( { b , a } ) ) )
    ++bad;
   } );
  }
 bool sym_nb = true;
 for( const auto & [ a , b ] : nb )
  sym_nb = sym_nb && nb.count( { b , a } ) && ( a != b ) &&
           ( a < d.num_nblocks() ) && ( b < d.num_nblocks() );
 check( ( bad == 0 ) && sym_nb , "MSA: the keys of the successors from the "
        "key of the parent, and the nblocks, symmetric and next to each "
        "other along every move" );
 // the hyperplanes: the levels, and the bound of Remark 1 of Kobayashi,
 // Kishimoto and Watanabe on the workers that the successors of the states
 // of a worker go to, ceil( n / d ) + 1 for d >= 1 (n / d + 1 there, which
 // the thickness 5 exceeds) and ( n + 1 ) / d for d < 1
 static_assert( LevelDomain< D > );
 bool levels = d.max_level() == 19 + 19 + 20 + 15;
 for( const auto & sg : g )
  d.for_each_successor( sg.first , 0 , [ & ]( const auto & t , int , int ) {
   const auto & s = sg.first;
   levels = levels && ( d.level( t ) > d.level( s ) ) &&
            ( d.level( t ) <= d.level( s ) + 4 );
   } );
 check( levels , "MSA: every move raises the level by 1 to 4, up to 73" );
 const unsigned n = unsigned( seq.size() );
 for( const auto [ planes , split ] : { std::pair( 1u , 1u ) ,
                                        std::pair( 2u , 1u ) ,
                                        std::pair( 5u , 1u ) ,
                                        std::pair( 1u , 2u ) ,
                                        std::pair( 1u , 4u ) } ) {
  const HyperplaneDistribution< D > h( d , 3 , planes , split );
  std::vector< std::set< unsigned > > to( 64 );
  for( const auto & sg : g )
   d.for_each_successor( sg.first , 0 , [ & ]( const auto & t , int ,
                                               int ) {
    to[ h.owner( sg.first , 64 ) ].insert( h.owner( t , 64 ) );
    } );
  std::size_t most = 0;
  for( const auto & w : to )
   most = std::max( most , w.size() );
  const std::size_t bound = split == 1 ? ( n + planes - 1 ) / planes + 1
                                        : ( n + 1 ) * split;
  check( most <= bound , "MSA, 64 workers: the successors of the states of "
         "a worker reach at most " + std::to_string( most ) + " workers "
         "with the planes of thickness " + std::to_string( planes ) + "/" +
         std::to_string( split ) + ", against the bound " +
         std::to_string( bound ) );
  }
 bool refused_hwd = false;
 try {
  const HyperplaneDistribution< D > bad( d , 1 , 2 , 2 );
  }
 catch( const std::invalid_argument & ) {
  refused_hwd = true;
  }
 const auto th = HyperplaneDistribution< D >::thickness( d , 64 );
 check( refused_hwd && ( th == std::pair( 1u , 27u ) ) &&
        ( HyperplaneDistribution< D >::thickness( d , 1 ) ==
          std::pair( 1u , 1u ) ) , "the thickness 1/27 of the rule of "
        "Kobayashi et al. for 73 residues and 64 workers, and 2/2 refused" );

 // the random instances of Jinnai and Fukunaga
 bool rnd = true;
 for( std::uint64_t seed = 1 ; seed <= 20 ; ++seed ) {
  const auto r = D::random_sequences( seed );
  const int k = int( r.size() );
  rnd = rnd && ( k >= 4 ) && ( k <= 9 ) &&
        ( D::random_sequences( seed ) == r );
  for( const auto & x : r ) {
   rnd = rnd && ( 10 * k * int( x.size() ) > 9 * 5000 ) &&
         ( 10 * k * int( x.size() ) < 11 * 5000 );
   for( char c : x )
    rnd = rnd && ( D::symbol( c ) >= 0 ) && ( D::symbol( c ) < 20 );
   }
  }
 const auto small = D::random_sequences( 7 , 3 , 30 );
 const D ds( small );
 check( rnd && ( small.size() == 3 ) && ( small[ 0 ].size() == 10 ) &&
        ( alignment_dp( small , ds ) == astar( ds ).cost ) ,
        "random sequences: 4 to 9 of lengths within 10% of 5000 / n, of "
        "the 20 amino acids, the same for the same seed" );
 check( ( ka > 4 * kz ) && ( kaz > 4 * kz ) , "MSA, 64 workers: of " +
        std::to_string( moves ) + " moves, the boxes of side 3 keep " +
        std::to_string( ka ) + " on the owner of their parent with AHDA*, "
        + std::to_string( kaz ) + " with AZH, and Zobrist hashing " +
        std::to_string( kz ) );
 }

/*--------------------------------------------------------------------------*/
// the cost of a cheapest path of a GridBoard by Dijkstra's algorithm,
// written from the rules of the format: leaving the row y costs 1 (unit)
// or y (life), times sqrt( 2 ) to a corner, and a move to a corner never
// looks at the cells it cuts

static double board_dijkstra( const GridBoard & b )
{
 const GridMap & m = b.map;
 const int w = m.width();
 using Entry = std::pair< double , int >;
 std::vector< double > d( std::size_t( w ) * m.height() ,
                          std::numeric_limits< double >::infinity() );
 std::priority_queue< Entry , std::vector< Entry > ,
                      std::greater< Entry > > q;
 d[ b.sy * w + b.sx ] = 0;
 q.emplace( 0.0 , b.sy * w + b.sx );
 while( ! q.empty() ) {
  const auto [ du , u ] = q.top();
  q.pop();
  if( du > d[ u ] )
   continue;
  const int x = u % w , y = u / w;
  if( ( x == b.gx ) && ( y == b.gy ) )
   return( du );
  const double c = b.cost == GridCost::life ? y : 1;
  for( int dy = -1 ; dy <= 1 ; ++dy )
   for( int dx = -1 ; dx <= 1 ; ++dx ) {
    if( ( ! dx && ! dy ) || ( ( b.moves == 4 ) && dx && dy ) ||
        ! m.cost( x + dx , y + dy ) )
     continue;
    const double nd = du + ( dx && dy ? c * std::sqrt( 2.0 ) : c );
    const int v = ( y + dy ) * w + x + dx;
    if( nd < d[ v ] ) {
     d[ v ] = nd;
     q.emplace( nd , v );
     }
    }
  }
 return( -1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a random GridBoard of the given moves and costs, its start and goal in
// the first and the last third of the board

static GridBoard random_board( int w , int h , double p , std::uint64_t seed ,
                               int moves , GridCost cost )
{
 GridBoard b;
 b.map = GridMap::random( w , h , p , seed );
 std::mt19937_64 gen( seed );
 b.sx = int( random_below( gen , std::uint64_t( w / 3 ) ) );
 b.sy = int( random_below( gen , std::uint64_t( h / 3 ) ) );
 b.gx = w - 1 - int( random_below( gen , std::uint64_t( w / 3 ) ) );
 b.gy = h - 1 - int( random_below( gen , std::uint64_t( h / 3 ) ) );
 b.map.set_cost( b.sx , b.sy , 1 );
 b.map.set_cost( b.gx , b.gy , 1 );
 b.moves = moves;
 b.cost = cost;
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the grids in the format of Burns et al.: the reader and its errors, the
// life costs and the heuristics, exact on the boards without blocked cells,
// A* against Dijkstra and the bounds on its expansions, the corners cut,
// the boards of the data with both origins of the rows

static void test_board( const std::string & data )
{
 std::istringstream tiny( "3 4\nBoard:\n  # \n#\n    \nLife\nEight-way\n"
                          "0 0 3 2\n" );
 const GridBoard t = GridBoard::read( tiny );
 check( ( t.map.width() == 4 ) && ( t.map.height() == 3 ) &&
        ( t.map.cost( 2 , 0 ) == 0 ) && ( t.map.cost( 0 , 1 ) == 0 ) &&
        ( t.map.cost( 1 , 1 ) == 1 ) && ( t.sx == 0 ) && ( t.sy == 2 ) &&
        ( t.gx == 3 ) && ( t.gy == 0 ) && ( t.moves == 8 ) &&
        ( t.cost == GridCost::life ) , "a board read with y from the bottom"
        ", short rows free" );
 std::istringstream top( "2 2\nBoard:\n# \n  \nUnit\nFour-way\n0 0 1 1\n" );
 const GridBoard tb = GridBoard::read( top , false );
 check( ( tb.sy == 0 ) && ( tb.gy == 1 ) && ( tb.map.cost( 0 , 0 ) == 1 ) &&
        ( tb.moves == 4 ) && ( tb.cost == GridCost::cells ) ,
        "a board read with y from the top, the blocked start freed" );
 int refused = 0;
 for( const char * bad : { "2 2\nBoard\n  \n  \nUnit\nFour-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nOther\nFour-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nUnit\nSix-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nUnit\nFour-way\n0 0 2 1\n" ,
                           "2 2\nBoard:\n  \n" } )
  try {
   std::istringstream in( bad );
   ( void ) GridBoard::read( in );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 5 , "boards without Board:, with other costs or moves,"
                       " a goal outside or missing rows are refused" );

 // the life heuristics: the exact cost on a board without blocked cells,
 // the moves to a corner costing as the straight ones (relaxed by hand)
 int wrong = 0;
 for( int moves : { 4 , 8 } )
  for( int gy = 0 ; gy < 7 ; gy += 3 )
   for( int gx = 0 ; gx < 9 ; gx += 4 ) {
    GridBoard e;
    e.map = GridMap( 9 , 7 );
    e.sx = 0;
    e.sy = 0;
    e.gx = gx;
    e.gy = gy;
    e.moves = moves;
    e.cost = GridCost::life;
    std::vector< std::int64_t > dist( 63 , std::int64_t( 1 ) << 40 );
    dist[ gy * 9 + gx ] = 0;
    for( int it = 0 ; it < 63 ; ++it )        // Bellman-Ford to the goal
     for( int u = 0 ; u < 63 ; ++u )
      for( int dy = -1 ; dy <= 1 ; ++dy )
       for( int dx = -1 ; dx <= 1 ; ++dx ) {
        const int x = u % 9 + dx , y = u / 9 + dy;
        if( ( ! dx && ! dy ) || ( ( moves == 4 ) && dx && dy ) ||
            ( x < 0 ) || ( y < 0 ) || ( x >= 9 ) || ( y >= 7 ) )
         continue;
        dist[ u ] = std::min( dist[ u ] ,
                              std::int64_t( u / 9 ) + dist[ y * 9 + x ] );
        }
    for( int u = 0 ; u < 63 ; ++u ) {
     double h;
     if( moves == 4 )
      h = Grid< 4 >( e ).heuristic( std::uint32_t( u ) );
     else
      h = Grid< 8 >( e ).heuristic( std::uint32_t( u ) );
     wrong += std::abs( h - double( dist[ u ] ) ) > 1e-9;
     }
    }
 check( wrong == 0 , "the life heuristics are the exact costs of the "
                     "relaxation on boards without blocked cells" );

 // random boards of both moves and costs, against Dijkstra
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed )
  for( GridCost cost : { GridCost::cells , GridCost::life } ) {
   const std::string tag = std::string( cost == GridCost::life ? " life"
                                                               : " unit" ) +
                           ", seed " + std::to_string( seed );
   const auto b4 = random_board( 60 , 40 , 0.3 , seed , 4 , cost );
   const Grid< 4 > g4( b4 , 10 , 8 );
   const double r4 = board_dijkstra( b4 );
   if( r4 >= 0 ) {
    check_solution( g4 , astar( g4 ) , int( std::lround( r4 ) ) ,
                    "board 4" + tag + " against Dijkstra on the map" );
    check_consistent_astar( g4 , all_distances( g4 ) , astar( g4 ) ,
                            "board 4" + tag );
    }
   else
    check( ! astar( g4 ).solved() , "board 4 unsolvable" + tag );
   const auto b8 = random_board( 60 , 40 , 0.3 , seed , 8 , cost );
   const Grid< 8 > g8( b8 , 10 , 8 );
   const double r8 = board_dijkstra( b8 );
   if( r8 >= 0 ) {
    check_solution( g8 , astar( g8 ) , r8 , "board 8" + tag +
                    " against Dijkstra on the map" );
    check_consistent_astar( g8 , all_distances( g8 ) , astar( g8 ) ,
                            "board 8" + tag );
    }
   else
    check( ! astar( g8 ).solved() , "board 8 unsolvable" + tag );
   }

 // a corner cut between two blocked cells
 GridBoard c;
 c.map = GridMap( 2 , 2 );
 c.map.set_cost( 1 , 0 , 0 );
 c.map.set_cost( 0 , 1 , 0 );
 c.gx = 1;
 c.gy = 1;
 c.moves = 8;
 const Grid< 8 > cut( c );
 check_solution( cut , astar( cut ) , grid_sqrt2 , "a corner cut between "
                 "two blocked cells" );
 const Grid< 8 > nocut( c.map , 0 , 0 , 1 , 1 );
 check( ! astar( nocut ).solved() , "no corner cut by default" );
 bool threw = false;
 try {
  const Grid< 4 > wrong4( c );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "a board of 8 moves refused by Grid< 4 >" );

 // the boards of the data, their rows from the bottom or from the top,
 // against Dijkstra and the costs found by the code of their authors
 struct File { const char * name; bool bottom; int cost; };
 const File files[] = { { "ethan_unit_four-way_1" , true , 23 } ,
                        { "jordan_unit_four-way_100x100_1" , true , 143 } ,
                        { "jordan_unit_four-way_500x500_1" , true , 763 } ,
                        { "jordan_life_four-way_500x500_1" , true , 362873 } ,
                        { "grid_500-500-0.7_0" , false , 755 } };
 for( const auto & [ name , bottom , cost ] : files ) {
  std::ifstream in( data + "/grid/" + name );
  const GridBoard b = GridBoard::read( in , bottom );
  const Grid< 4 > g( b );
  check( std::lround( board_dijkstra( b ) ) == cost , std::string( "board " )
         + name + ": Dijkstra costs " + std::to_string( cost ) );
  check_solution( g , astar( g ) , cost , std::string( "board " ) + name );
  }
 }

/*--------------------------------------------------------------------------*/
// the other files of the data: the 15-puzzle by the cells of the tiles,
// the 24-puzzle sets, and the cities of the traveling salesperson

static void test_files( const std::string & data )
{
 const auto inst = korf( data );
 for( int id : { 1 , 12 , 55 } ) {
  std::ifstream in( data + "/15-puzzle/korf_" + std::to_string( id ) );
  const auto [ start , goal ] = FifteenPuzzle::read_positions( in );
  check( ( start == inst[ id - 1 ].board ) &&
         ( goal == FifteenPuzzle::korf_goal() ) , "korf_" +
         std::to_string( id ) + " is the instance " + std::to_string( id ) +
         " of Korf" );
  if( id != 1 ) {
   const FifteenPuzzle d( start , goal );
   check_solution( d , astar( d ) , inst[ id - 1 ].cost , "A* on korf_" +
                   std::to_string( id ) );
   }
  }
 std::istringstream small( "3 3\nstarting positions for each tile:\n"
                           "1 0 2 3 4 5 6 7 8\ngoal positions:\n"
                           "0 1 2 3 4 5 6 7 8\n" );
 const auto [ s8 , g8 ] = EightPuzzle::read_positions( small );
 check( ( s8[ 0 ] == 1 ) && ( s8[ 1 ] == 0 ) && ( g8 == EightPuzzle::
        korf_goal() ) , "an 8-puzzle by the cells of its tiles" );
 int refused = 0;
 for( const char * bad : { "4 4\nstarting positions for each tile:\n" ,
                           "3 3\nstarting positions for each tile:\n"
                           "0 0 1 2 3 4 5 6 7\n" } )
  try {
   std::istringstream in( bad );
   ( void ) EightPuzzle::read_positions( in );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 2 , "a board of another size and a repeated cell are "
                       "refused" );

 for( const auto & [ file , count ] : { std::pair( "instances100.txt" , 100 ) ,
                                        std::pair( "fromkorf8.txt" , 8 ) } ) {
  std::ifstream in( data + "/24-puzzle/" + file );
  std::string line;
  int n = 0 , bad = 0;
  while( std::getline( in , line ) ) {
   if( line.empty() || ( line[ 0 ] == '#' ) )
    continue;
   try {
    const TwentyFourPuzzle d( TwentyFourPuzzle::parse( line ) );
    bad += d.heuristic( d.initial_state() ) <= 0;
    }
   catch( std::invalid_argument & ) {
    ++bad;
    }
   ++n;
   }
  check( ( n == count ) && ( bad == 0 ) , std::string( "the instances of " )
         + file + " are boards from which the goal is reachable" );
  }

 for( auto metric : { TravelingSalesman::Metric::euclidean ,
                      TravelingSalesman::Metric::manhattan } ) {
  const std::string tag = metric == TravelingSalesman::Metric::euclidean
                          ? "Euclidean" : "Manhattan";
  for( int k : { 6 , 9 , 12 } ) {
   std::ifstream in( data + "/tsp/tsp_40_0" );
   const auto m = TravelingSalesman::read_cities( in , metric , k );
   bool integral = true;
   for( const auto & row : m )
    for( double x : row )
     integral = integral && ( x == std::floor( x ) );
   const TravelingSalesman d( m );
   check_solution( d , astar( d ) , d.held_karp() , "A* on the first " +
                   std::to_string( k ) + " cities of tsp_40_0, " + tag );
   if( metric == TravelingSalesman::Metric::manhattan )
    check( integral , "the Manhattan distances are integers" );
   }
  }
 // the optimal tours of the first 6, 8, 10 and 12 cities of tsp_40_0 with
 // the Manhattan distances, those of the code of HDA* of their authors
 const std::pair< int , double > tours[] = { { 6 , 28975 } , { 8 , 31254 } ,
                                             { 10 , 32734 } ,
                                             { 12 , 33861 } };
 for( const auto & [ k , cost ] : tours ) {
  std::ifstream tf( data + "/tsp/tsp_40_0" );
  const TravelingSalesman d( TravelingSalesman::read_cities(
                       tf , TravelingSalesman::Metric::manhattan , k ) );
  check_solution( d , astar( d ) , cost , "the tour of " + std::to_string( k )
                  + " cities of tsp_40_0 of the reference code" );
  }
 std::ifstream in( data + "/tsp/tsp_40_0" );
 const auto all = TravelingSalesman::read_cities( in );
 std::istringstream few( "3\n0 0\n1 1\n" );
 bool threw = false;
 try {
  ( void ) TravelingSalesman::read_cities( few );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( ( all.size() == 40 ) && threw , "40 cities read, and a missing one "
                                        "refused" );
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
// the owners and the keys that the Zobrist distributions give to the
// successors of the states of random walks from the key of their parent,
// against those computed from scratch

template< class D >
static void check_successor_owners( void )
{
 int bad = 0;
 for( std::uint64_t seed = 1 ; seed <= 20 ; ++seed ) {
  const D d( random_walk< D >( 30 , seed ) );
  const ZobristDistribution< D > z( d );
  const AbstractZobristDistribution< D > az( d );
  static_assert( IncrementalDistribution< decltype( z ) , D > &&
                 IncrementalDistribution< decltype( az ) , D > );
  const auto s = d.initial_state();
  for( unsigned p : { 1u , 2u , 3u , 7u , 64u } ) {
   const auto oz = successor_owners< D >( z , s , z.key( s ) , p );
   const auto oaz = successor_owners< D >( az , s , az.key( s ) , p );
   d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
    if( ( oz( t ) != std::pair( z.owner( t , p ) , z.key( t ) ) ) ||
        ( oaz( t ) != std::pair( az.owner( t , p ) , az.key( t ) ) ) )
     ++bad;
    } );
   }
  }
 check( bad == 0 , "the owners and the keys of the successors from the "
                   "key of the parent, " + std::to_string( D::cells ) +
                   " cells" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// HDA* with every distribution, 1 to 8 workers, on the given runtime

template< class Executor >
static void test_hda( const std::string & data )
{
 check_successor_owners< EightPuzzle >();
 check_successor_owners< FifteenPuzzle >();
 check_successor_owners< TwentyFourPuzzle >();
 const auto inst = korf( data );
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 // with more workers than cores, a worker without a core may hold the path
 // to the goal while the others expand above the optimum, which no
 // incumbent prunes yet: the 15-puzzle only with a core for each worker
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  ParallelOptions o;
  o.threads = p;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   if( p > cores )
    break;
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

  using M = SequenceAlignment<>;
  const M md( fasta( data + "/msa/easy/test3.tfa" ) , pam250 , 4 , 6 );
  const std::string mt = " MSA test3, p " + std::to_string( p );
  const ZobristDistribution< M > mz( md );
  check_solution( md , HDAStar< M , decltype( mz ) , Executor >(
                        md , mz , o ).solve() , 2022 , "HDA*[Z]" + mt );
  const AbstractZobristDistribution< M > maz( md );
  check_solution( md , HDAStar< M , decltype( maz ) , Executor >(
                        md , maz , o ).solve() , 2022 , "HDA*[AZ]" + mt );
  const AbstractStateDistribution< M > ma( md );
  check_solution( md , HDAStar< M , decltype( ma ) , Executor >(
                        md , ma , o ).solve() , 2022 , "HDA*[A]" + mt );
  for( const auto [ planes , split ] : { std::pair( 1u , 1u ) ,
                                         std::pair( 3u , 1u ) ,
                                         std::pair( 1u , 4u ) ,
                                         HyperplaneDistribution< M >::
                                          thickness( md , p ) } ) {
   const HyperplaneDistribution< M > mh( md , 1 , planes , split );
   check_solution( md , HDAStar< M , decltype( mh ) , Executor >(
                         md , mh , o ).solve() , 2022 ,
                   "HDA*[Hyperplane] d " + std::to_string( planes ) + "/" +
                   std::to_string( split ) + mt );
   }

  const auto lb = random_board( 60 , 40 , 0.3 , 2 , 8 , GridCost::life );
  const Grid< 8 > lg( lb , 10 , 8 );
  const AbstractZobristDistribution< Grid< 8 > > lz( lg );
  check_solution( lg , HDAStar< Grid< 8 > , decltype( lz ) , Executor >(
                        lg , lz , o ).solve() , board_dijkstra( lb ) ,
                  "HDA*[AZ] life board 8, p " + std::to_string( p ) );
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
  const SequenceAlignment<> md( fasta( data + "/msa/easy/test0.tfa" ) );
  check_solution( md , SPAStar< SequenceAlignment<> , Executor >( md , o )
                        .solve() , 1741 ,
                  "SPA* MSA test0, p " + std::to_string( p ) );
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
// MultiQueue A* with 1 to 8 workers, on the given runtime

template< class Executor >
static void test_mqa( const std::string & data )
{
 const auto inst = korf( data );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  ParallelOptions o;
  o.threads = p;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   check_solution( d , MQAStar< FifteenPuzzle , Executor >( d , o ).solve() ,
                   ki.cost , "MQA* p " + std::to_string( p ) + ", Korf #" +
                             std::to_string( id ) );
   }
  auto rg = romania( Arad );
  check_solution( rg , MQAStar< ExplicitGraph< int > , Executor >( rg , o )
                        .solve() , 418 ,
                  "MQA* Romania, p " + std::to_string( p ) );
  const SequenceAlignment<> md( fasta( data + "/msa/easy/test0.tfa" ) );
  check_solution( md , MQAStar< SequenceAlignment<> , Executor >( md , o )
                        .solve() , 1741 ,
                  "MQA* MSA test0, p " + std::to_string( p ) );
  }

 const auto mqa = []( const EightPuzzle & d , const ParallelOptions & o ) {
  return( MQAStar< EightPuzzle , Executor >( d , o ).solve() );
  };
 ParallelOptions o;
 o.threads = 8;
 check_astar_edges( o , mqa , "MQA* 8 workers" );
 o.threads = 0;
 check_astar_edges( o , mqa , "MQA* 0 workers" );
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
  using M = SequenceAlignment<>;
  const M md( fasta( data + "/msa/easy/test3.tfa" ) , pam250 , 4 , 5 );
  check_solution( md , SafePBNF< M , Executor >( md , o ).solve() , 2022 ,
                  "SafePBNF MSA test3, p " + std::to_string( p ) );
  const auto lb = random_board( 60 , 40 , 0.3 , 3 , 4 , GridCost::life );
  const Grid< 4 > lg( lb , 10 , 8 );
  check_solution( lg , SafePBNF< Grid< 4 > , Executor >( lg , o ).solve() ,
                  int( std::lround( board_dijkstra( lb ) ) ) ,
                  "SafePBNF life board 4, p " + std::to_string( p ) );
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
  const auto lb = random_board( 60 , 40 , 0.3 , 1 , 8 , GridCost::life );
  const Grid< 8 > lg( lb );
  check_solution( lg , PASE< Grid< 8 > , Executor >( lg , o ).solve() ,
                  board_dijkstra( lb ) , "PA*SE life board 8, p " +
                                         std::to_string( p ) );

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
   const SequenceAlignment<> md( fasta( data + "/msa/easy/test0.tfa" ) );
   check_solution( md , FarmAStar< SequenceAlignment<> , Executor >( md , o )
                        .solve() , 1741 , "FarmAStar MSA test0" + tag );
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

/*--------------------------------------------------------------------------*/
/*------------------ BOUNDED-SUBOPTIMAL AND ANYTIME SEARCH -----------------*/
/*--------------------------------------------------------------------------*/
// the cost of the path of sol on d, the cheapest edge between two states,
// and whether the path is one from the initial state to a goal

template< class D , class S >
static std::pair< bool , cost_t< D > > path_cost( const D & d , const S & sol )
{
 bool ok = ! sol.path.empty() && ( sol.path.front() == d.initial_state() ) &&
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
 return( std::make_pair( ok , sum ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that sol is solved with a valid path of cost within [ C* , w C* ],
// whose bound is within [ 0 , C* ], and, if it has improvements, that
// they decrease, stay within [ C* , ... ] with bounds below C*, and end
// with the cost of the path

template< class D , class S >
static void check_bounded( const D & d , const S & sol , cost_t< D > cstar ,
                           double w , const std::string & what )
{
 const auto [ valid , sum ] = path_cost( d , sol );
 const double c = double( cstar );
 const double tol = 1e-9 * std::max( 1.0 , c );
 bool ok = sol.solved() && valid && same_cost( sum , sol.cost ) &&
           ( double( sol.cost ) >= c - tol ) &&
           ( double( sol.cost ) <= w * c + tol ) &&
           ( double( sol.bound ) <= c + tol ) &&
           ( double( sol.bound ) >= 0 );
 for( std::size_t i = 0 ; ok && ( i < sol.improvements.size() ) ; ++i ) {
  const auto & x = sol.improvements[ i ];
  ok = ( double( x.cost ) >= c - tol ) && ( double( x.bound ) <= c + tol ) &&
       ( ( i == 0 ) || ( x.cost < sol.improvements[ i - 1 ].cost ) );
  }
 if( ok && ! sol.improvements.empty() )
  ok = same_cost( sol.improvements.back().cost , sol.cost );
 std::ostringstream os;
 os << what << " (w " << w << ", cost " << sol.cost << ", C* " << cstar
    << ", bound " << sol.bound << ", expanded " << sol.stats.expanded
    << ", paths " << sol.improvements.size() << ")";
 check( ok , os.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the instances of the sequential bounded searches, each with its optimal
// cost: states of the 8-puzzle, easy instances of Korf, grids of costs 1 to
// 5 with 4 and 8 moves, Romania and the inconsistent heuristic of
// test_reopen(); solve( d , w , what ) checks one of them

template< class Solve >
static void bounded_instances( const std::string & data , Solve && solve ,
                               bool consistent_only = false ,
                               std::size_t korfs = 3 )
{
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 for( const auto & [ s , dst ] : eight_sample( dist ) )
  solve( EightPuzzle( EightPuzzle::unpack( s ) ) , dst , "8-puzzle" );
 // states whose first path the second search of IOS reaches anew, giving
 // its nodes other parents
 for( const std::uint64_t s : { std::uint64_t( 1380455286 ) ,
                                std::uint64_t( 338851941 ) ,
                                std::uint64_t( 4336407651 ) } )
  solve( EightPuzzle( EightPuzzle::unpack( s ) ) , dist.at( s ) ,
         "8-puzzle" );
 const auto inst = korf( data );
 for( std::size_t i = 0 ; i < korfs ; ++i ) {
  const int id = easy_korf[ i ];
  const auto & ki = inst[ id - 1 ];
  solve( FifteenPuzzle( ki.board ) , ki.cost ,
         "Korf #" + std::to_string( id ) );
  }
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto mc = GridMap::random( 120 , 80 , 0.25 , seed , 5 );
  const std::string t = " seed " + std::to_string( seed );
  const double r4 = grid_dijkstra( mc , 4 , 0 , 79 , 119 , 79 );
  if( r4 >= 0 )
   solve( Grid< 4 >( mc , 0 , 79 , 119 , 79 ) , int( std::lround( r4 ) ) ,
          "grid 4" + t );
  const double r8 = grid_dijkstra( mc , 8 , 0 , 79 , 119 , 79 );
  if( r8 >= 0 )
   solve( Grid< 8 >( mc , 0 , 79 , 119 , 79 ) , r8 , "grid 8" + t );
  }
 solve( romania( Arad ) , 418 , "Romania" );
 if( consistent_only )
  return;
 ExplicitGraph< int > g( 5 , 0 , 4 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 2 , 1 , true );
 g.add_edge( 0 , 2 , 3 , true );
 g.add_edge( 2 , 3 , 1 , true );
 g.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };
 for( std::uint32_t u = 0 ; u < 5 ; ++u )
  g.set_heuristic( u , h[ u ] );
 solve( g , 4 , "inconsistent heuristic" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the edge cases of a sequential search: the initial state a goal, an
// unsolvable grid, and the node limit

template< class Solve >
static void bounded_edges( Solve && solve , const std::string & name )
{
 const EightPuzzle g( EightPuzzle::korf_goal() );
 const auto s0 = solve( g , std::size_t( 0 ) );
 check( s0.solved() && ( s0.cost == 0 ) && ( s0.path.size() == 1 ) ,
        name + " initial state that is the goal" );
 GridMap m( 30 , 30 );
 for( int y = 0 ; y < 30 ; ++y )
  m.set_cost( 15 , y , 0 );
 const Grid< 4 > wall( m , 0 , 0 , 29 , 29 );
 check( solve( wall , std::size_t( 0 ) ).status == Status::unsolvable ,
        name + " proves the wall unsolvable" );
 const EightPuzzle far( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                 0 , 3 , 4 , 7 , 6 } ) );
 check( solve( far , std::size_t( 5 ) ).status == Status::node_limit ,
        name + " node limit of 5" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the class Weight: its fractions, keys, rules and bounds

static void check_weight( void )
{
 const Weight< int > a( 1.5 );
 check( ( a.key( 2 , 3 ) == 13 ) && ( a.scale( 10 ) == 20 ) &&
        ( a.lower_bound( 10 ) == 7 ) && ( a.times( 10 ) == 15 ) &&
        ( a.value() == 1.5 ) , "Weight 1.5 is 3 / 2 with integer costs" );
 const Weight< int > b( 1.1 );
 check( ( b.key( 1 , 1 ) == 21 ) && ( b.times( 10 ) == 11 ) ,
        "Weight 1.1 is 11 / 10 with integer costs" );
 const Weight< double > c( 2 );
 check( ( c.key( 1 , 1 ) == 3 ) && ( c.lower_bound( 5 ) == 2.5 ) &&
        c.prunes( 2 , 1 , 6 ) && ! c.prunes( 2 , 1 , 6.5 ) ,
        "Weight 2 with real costs" );
 const Weight< int > one( 1 , true , false );
 check( one.unit() && ! one.anytime() && ( one.key( 3 , 4 ) == 7 ) &&
        ( one.lower_bound( 9 ) == 9 ) && one.prunes( 3 , 4 , 7 ) &&
        ! one.prunes( 3 , 4 , 8 ) , "Weight 1 is A*" );
 const Weight< int > dd( 2 , false , true );
 const auto hu = []( int h ) { return( [ h ] { return( h ); } ); };
 check( dd.improves( 10 , 8 , true , hu( 5 ) , 4 ) &&
        ! dd.improves( 10 , 9 , true , hu( 6 ) , 4 ) &&
        dd.improves( 10 , 9 , false , hu( 6 ) , 4 ) &&
        ! dd.improves( 10 , 10 , false , hu( 0 ) , 0 ) ,
        "Weight::improves() with duplicate dropping" );
 const Weight< int > at( 3 , true );
 check( at.anytime() && at.live( 1000 , 5 ) && at.prunes( 3 , 2 , 5 ) &&
        ! at.prunes( 3 , 1 , 5 ) && ( at.proven( 7 ) == 7 ) ,
        "Weight 3, anytime" );
 bool thrown = false;
 try { Weight< int > x( 0.5 ); } catch( std::invalid_argument & ) {
  thrown = true; }
 check( thrown , "Weight below 1 throws" );
 thrown = false;
 try { Weight< int > x( 2 , true , true ); } catch( std::invalid_argument & ) {
  thrown = true; }
 check( thrown , "Weight anytime with duplicate dropping throws" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// weighted A* (with and without duplicate dropping, both OPEN lists) and
// weighted IDA*

static void test_weighted( const std::string & data )
{
 check_weight();
 for( double w : { 1.0 , 1.5 , 2.0 , 5.0 } ) {
  const std::string tw = " w " + std::to_string( w ).substr( 0 , 4 ) + " ";
  for( bool dd : { false , true } )
   bounded_instances( data , [ & ]( const auto & d , auto cstar ,
                                    const std::string & what ) {
    SearchOptions o;
    o.weight = w;
    o.drop_duplicates = dd;
    check_bounded( d , astar( d , OpenKind::automatic , o ) , cstar , w ,
                   std::string( dd ? "wA*dd" : "wA*" ) + tw + what );
    if( dd )
     check_bounded( d , astar( d , OpenKind::heap , o ) , cstar , w ,
                    "wA*dd heap" + tw + what );
    else
     check_bounded( d , astar( d , OpenKind::heap_low_g , o ) , cstar , w ,
                    "wA* low g" + tw + what );
    } );
  bounded_instances( data , [ & ]( const auto & d , auto cstar ,
                                   const std::string & what ) {
   if( ( what.rfind( "grid" , 0 ) == 0 ) ||
       ( ( what.rfind( "Korf" , 0 ) == 0 ) && ( w < 2 ) ) )
    return;              // the many paths of a grid, Korf below 2: too slow
   SearchOptions o;
   o.weight = w;
   check_bounded( d , idastar( d , o ) , cstar , w , "wIDA*" + tw + what );
   } );
  }

 // the weight 1 is A*, with or without duplicate dropping
 const auto inst = korf( data );
 const FifteenPuzzle d( inst[ easy_korf[ 0 ] - 1 ].board );
 SearchOptions o;
 o.drop_duplicates = true;
 const auto x = astar( d ) , y = astar( d , OpenKind::automatic , o );
 check( ( x.stats.expanded == y.stats.expanded ) && ( x.cost == y.cost ) &&
        ( x.bound == x.cost ) , "the weight 1 is A*, with duplicate dropping"
                                " too" );
 // a larger weight expands fewer nodes on Korf #55, as it does in general
 o = SearchOptions();
 o.weight = 3;
 check( astar( d , OpenKind::automatic , o ).stats.expanded <
        x.stats.expanded , "wA* with w 3 expands fewer nodes than A*" );

 for( double w : { 1.0 , 2.0 } )
  bounded_edges( [ w ]( const auto & e , std::size_t limit ) {
   SearchOptions p;
   p.weight = w;
   p.node_limit = limit;
   return( astar( e , OpenKind::automatic , p ) );
   } , "wA* w " + std::to_string( int( w ) ) );
 bool thrown = false;
 try {
  SearchOptions p;
  p.weight = 2;
  p.anytime = true;
  const EightPuzzle e( EightPuzzle::korf_goal() );
  static_cast< void >( AStar< EightPuzzle >( e , p ) );
  }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "AStar refuses the anytime search" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// AWA* and ARA*: the optimal path at the end, the paths found on the way

static void test_anytime( const std::string & data )
{
 for( double w : { 1.5 , 3.0 } ) {
  const std::string tw = " w " + std::to_string( w ).substr( 0 , 3 ) + " ";
  bounded_instances( data , [ & ]( const auto & d , auto cstar ,
                                   const std::string & what ) {
   AnytimeOptions o;
   o.weight = w;
   const auto s = anytime_astar( d , o );
   check_bounded( d , s , cstar , 1 , "AWA*" + tw + what );
   check( ! s.improvements.empty() &&
          ( double( s.improvements.front().cost ) <=
            w * double( cstar ) + 1e-9 * double( cstar ) ) &&
          same_cost( s.bound , s.cost ) ,
          "AWA* first path within w C*, bound C* at the end," + tw + what );
   o.schedule = { w , 1 + ( w - 1 ) / 2 , 1 + ( w - 1 ) / 4 };
   const auto r = anytime_astar( d , o , OpenKind::heap );
   check_bounded( d , r , cstar , 1 , "ARA*" + tw + what );
   check( ! r.improvements.empty() &&
          ( double( r.improvements.front().cost ) <=
            w * double( cstar ) + 1e-9 * double( cstar ) ) ,
          "ARA* first path within w C*," + tw + what );
   } );
  }

 for( bool ara : { false , true } )
  bounded_edges( [ ara ]( const auto & e , std::size_t limit ) {
   AnytimeOptions p;
   p.weight = 2;
   if( ara )
    p.schedule = { 2 , 1.5 };
   p.node_limit = limit;
   return( anytime_astar( e , p ) );
   } , ara ? "ARA*" : "AWA*" );

 // with the node limit, the best path found is kept
 const auto inst = korf( data );
 const FifteenPuzzle d( inst[ easy_korf[ 0 ] - 1 ].board );
 AnytimeOptions o;
 o.weight = 5;
 const auto full = anytime_astar( d , o );
 o.node_limit = full.stats.stored / 2;
 const auto part = anytime_astar( d , o );
 check( ( part.status == Status::node_limit ) &&
        ( full.improvements.size() > 0 ) &&
        ( part.path.empty() || path_cost( d , part ).first ) ,
        "AWA* at the node limit keeps the best path found" );
 bool thrown = false;
 try {
  o.schedule = { 2 , 3 };
  static_cast< void >( AnytimeAStar< FifteenPuzzle >( d , o ) );
  }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "ARA* refuses a schedule that increases" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the optimistic search, IOS, the focal searches, EES and DPS

static void test_focal( const std::string & data )
{
 for( double w : { 1.2 , 2.0 , 4.0 } ) {
  const std::string tw = " w " + std::to_string( w ).substr( 0 , 3 ) + " ";
  bounded_instances( data , [ & ]( const auto & d , auto cstar ,
                                   const std::string & what ) {
   OptimisticOptions o;
   o.weight = w;
   check_bounded( d , optimistic_search( d , o ) , cstar , w ,
                  "optimistic" + tw + what );
   using D = std::decay_t< decltype( d ) >;
   FocalOptions f;
   f.weight = w;
   for( FocalOrder k : { FocalOrder::distance , FocalOrder::heuristic ,
                         FocalOrder::optimistic ,
                         FocalOrder::explicit_estimation } ) {
    f.order = k;
    const char * name[] = { "A*eps d" , "A*eps h" , "focal 2w-1" , "EES" };
    check_bounded( d , FocalSearch< D >( d , f ).solve() , cstar , w ,
                   name[ int( k ) ] + tw + what );
    }
   SearchOptions p;
   p.weight = w;
   check_bounded( d , PotentialSearch< D >( d , p ).solve() , cstar , w ,
                  "DPS" + tw + what );
   } , false , 1 );
  // IOS needs a consistent heuristic
  bounded_instances( data , [ & ]( const auto & d , auto cstar ,
                                   const std::string & what ) {
   using D = std::decay_t< decltype( d ) >;
   OptimisticOptions o;
   o.weight = w;
   check_bounded( d , ImprovedOptimisticSearch< D >( d , o ).solve() ,
                  cstar , w , "IOS" + tw + what );
   } , true , 1 );
  }

 bounded_edges( []( const auto & e , std::size_t limit ) {
  OptimisticOptions p;
  p.weight = 2;
  p.node_limit = limit;
  return( optimistic_search( e , p ) );
  } , "optimistic" );
 bounded_edges( []( const auto & e , std::size_t limit ) {
  OptimisticOptions p;
  p.weight = 2;
  p.node_limit = limit;
  return( ImprovedOptimisticSearch< std::decay_t< decltype( e ) > >( e , p )
          .solve() );
  } , "IOS" );
 for( FocalOrder k : { FocalOrder::distance ,
                       FocalOrder::explicit_estimation } )
  bounded_edges( [ k ]( const auto & e , std::size_t limit ) {
   FocalOptions p;
   p.weight = 2;
   p.order = k;
   p.node_limit = limit;
   return( FocalSearch< std::decay_t< decltype( e ) > >( e , p ).solve() );
   } , k == FocalOrder::distance ? "A*eps" : "EES" );
 bounded_edges( []( const auto & e , std::size_t limit ) {
  SearchOptions p;
  p.weight = 2;
  p.node_limit = limit;
  return( PotentialSearch< std::decay_t< decltype( e ) > >( e , p )
          .solve() );
  } , "DPS" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// LPA* on grids whose costs change, against A* from scratch after each
// change, on the 8-puzzle and on directed graphs

template< int Moves >
static void check_lpa_grid( std::uint64_t seed )
{
 using G = Grid< Moves >;
 auto m = GridMap::random( 60 , 40 , 0.2 , seed , 5 );
 G g( m , 0 , 39 , 59 , 39 );
 LPAStar< G > lpa( g );
 const std::string t = " grid " + std::to_string( Moves ) + ", seed " +
                       std::to_string( seed );
 std::uint64_t replan = 0 , scratch = 0;
 std::mt19937_64 rng( seed );
 for( int round = 0 ; round <= 10 ; ++round ) {
  if( round > 0 ) {
   std::vector< std::pair< int , int > > changed;
   for( int k = 0 ; k < 15 ; ++k ) {
    const int x = int( rng() % 60 ) , y = int( rng() % 40 );
    if( ( ( x == 0 ) && ( y == 39 ) ) || ( ( x == 59 ) && ( y == 39 ) ) )
     continue;
    m.set_cost( x , y , int( rng() % 6 ) );   // 0, a wall, to 5
    changed.emplace_back( x , y );
    }
   g.refresh();
   for( auto [ x , y ] : changed )
    for( int dy = -1 ; dy <= 1 ; ++dy )
     for( int dx = -1 ; dx <= 1 ; ++dx )
      if( ( x + dx >= 0 ) && ( x + dx < 60 ) && ( y + dy >= 0 ) &&
          ( y + dy < 40 ) )
       lpa.update( g.index( x + dx , y + dy ) );
   }
  const auto s = lpa.solve();
  const auto a = astar( g );
  replan += s.stats.expanded;
  scratch += a.stats.expanded;
  if( a.solved() )
   check_solution( g , s , a.cost , "LPA*" + t + ", round " +
                                    std::to_string( round ) );
  else
   check( s.status == Status::unsolvable , "LPA* unsolvable" + t +
                                           ", round " +
                                           std::to_string( round ) );
  }
 check( replan < 2 * scratch , "LPA* replanning on" + t + " (" +
        std::to_string( replan ) + " expansions, A* from scratch " +
        std::to_string( scratch ) + ")" );
 }

static void test_lpa( const std::string & )
{
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  check_lpa_grid< 4 >( seed );
  check_lpa_grid< 8 >( seed );
  }

 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 for( const auto & [ s , dst ] : eight_sample( dist ) ) {
  const EightPuzzle e( EightPuzzle::unpack( s ) );
  LPAStar< EightPuzzle > lpa( e );
  const auto x = lpa.solve();
  check_solution( e , x , dst , "LPA* 8-puzzle" );
  const auto y = lpa.solve();
  check( ( y.cost == x.cost ) && ( y.stats.expanded == 0 ) ,
         "LPA* with no change expands nothing" );
  }

 // a directed graph: 0 -> 1 -> 3 costs 2, 0 -> 2 -> 3 costs 4
 ExplicitGraph< int > g( 4 , 0 , 3 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 3 , 1 , true );
 g.add_edge( 0 , 2 , 2 , true );
 g.add_edge( 2 , 3 , 2 , true );
 g.add_edge( 3 , 0 , 1 , true );
 LPAStar< ExplicitGraph< int > > lg( g );
 check_solution( g , lg.solve() , 2 , "LPA* directed graph" );
 auto rg = romania( Arad );
 check_solution( rg , LPAStar< ExplicitGraph< int > >( rg ).solve() , 418 ,
                 "LPA* Romania" );
 GridMap m( 30 , 30 );
 for( int y = 0 ; y < 30 ; ++y )
  m.set_cost( 15 , y , 0 );
 const Grid< 4 > wall( m , 0 , 0 , 29 , 29 );
 check( LPAStar< Grid< 4 > >( wall ).solve().status == Status::unsolvable ,
        "LPA* proves the wall unsolvable" );
 const EightPuzzle goal( EightPuzzle::korf_goal() );
 const auto z = LPAStar< EightPuzzle >( goal ).solve();
 check( z.solved() && ( z.cost == 0 ) && ( z.path.size() == 1 ) ,
        "LPA* initial state that is the goal" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel searches with a weight, duplicate dropping and anytime, on
// the given runtime: within the bound, and optimal when anytime; SPA*,
// MQA* and SafePBNF only if Shared, as in their own cases

template< class Executor , bool Shared = true >
static void test_bounded( const std::string & data )
{
 const auto inst = korf( data );
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 const auto m = GridMap::random( 200 , 120 , 0.25 , 4 , 5 );
 const double r8 = grid_dijkstra( m , 8 , 0 , 119 , 199 , 119 );
 const Grid< 8 > g8( m , 0 , 119 , 199 , 119 , 20 , 10 );
 auto rg = romania( Arad );
 const auto mode = []( int k ) {
  return( std::string( k == 0 ? "" : ( k == 1 ? "dd " : "anytime " ) ) );
  };
 for( unsigned p : { 1u , 3u , 8u } ) {
  if( p > cores )
   break;
  for( double w : { 1.5 , 3.0 } )
   for( int k = 0 ; k < 3 ; ++k ) {
    PBNFOptions o;
    o.threads = p;
    o.weight = w;
    o.drop_duplicates = k == 1;
    o.anytime = k == 2;
    const double wb = k == 2 ? 1 : w;
    const std::string tag = " " + mode( k ) + "p " + std::to_string( p ) +
                            ", ";
    for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
     const auto & ki = inst[ id - 1 ];
     const FifteenPuzzle d( ki.board );
     const std::string t = tag + "Korf #" + std::to_string( id );
     const ZobristDistribution< FifteenPuzzle > z( d );
     check_bounded( d , HDAStar< FifteenPuzzle , decltype( z ) , Executor >(
                         d , z , o ).solve() , ki.cost , wb , "wHDA*" + t );
     check_bounded( d , FarmAStar< FifteenPuzzle , Executor >( d , o )
                         .solve() , ki.cost , wb , "wFarmA*" + t );
     if constexpr( Shared ) {
      check_bounded( d , SPAStar< FifteenPuzzle , Executor >( d , o )
                          .solve() , ki.cost , wb , "wSPA*" + t );
      check_bounded( d , MQAStar< FifteenPuzzle , Executor >( d , o )
                          .solve() , ki.cost , wb , "wMQA*" + t );
      check_bounded( d , SafePBNF< FifteenPuzzle , Executor >( d , o )
                          .solve() , ki.cost , wb , "wPBNF" + t );
      }
     }
    const AbstractZobristDistribution< Grid< 8 > > gz( g8 );
    check_bounded( g8 , HDAStar< Grid< 8 > , decltype( gz ) , Executor >(
                         g8 , gz , o ).solve() , r8 , wb , "wHDA*" + tag +
                         "grid 8" );
    const ZobristDistribution< ExplicitGraph< int > > rz( rg );
    check_bounded( rg , HDAStar< ExplicitGraph< int > , decltype( rz ) ,
                                 Executor >( rg , rz , o ).solve() , 418 , wb ,
                   "wHDA*" + tag + "Romania" );
    check_bounded( g8 , FarmAStar< Grid< 8 > , Executor >( g8 , o ).solve() ,
                   r8 , wb , "wFarmA*" + tag + "grid 8" );
    if constexpr( Shared )
     if( k == 0 ) {                         // ARPBNF, optimal at the end
      PBNFOptions a = o;
      a.schedule = { w , 1 + ( w - 1 ) / 2 , 1 + ( w - 1 ) / 4 };
      const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
      const FifteenPuzzle d( ki.board );
      check_bounded( d , SafePBNF< FifteenPuzzle , Executor >( d , a )
                          .solve() , ki.cost , 1 , "ARPBNF" + tag +
                          "Korf #" + std::to_string( easy_korf[ 0 ] ) );
      check_bounded( g8 , SafePBNF< Grid< 8 > , Executor >( g8 , a )
                           .solve() , r8 , 1 , "ARPBNF" + tag + "grid 8" );
      }
    if constexpr( Shared ) {
     check_bounded( g8 , SPAStar< Grid< 8 > , Executor >( g8 , o ).solve() ,
                    r8 , wb , "wSPA*" + tag + "grid 8" );
     check_bounded( g8 , MQAStar< Grid< 8 > , Executor >( g8 , o ).solve() ,
                    r8 , wb , "wMQA*" + tag + "grid 8" );
     check_bounded( g8 , SafePBNF< Grid< 8 > , Executor >( g8 , o ).solve() ,
                    r8 , wb , "wPBNF" + tag + "grid 8" );
     }
    }
  }
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
 // the weight, the duplicate dropping and the anytime search
 for( int k = 0 ; k < 3 ; ++k )
  for( double w : { 1.5 , 3.0 } ) {
   ParallelOptions x;
   x.weight = w;
   x.drop_duplicates = k == 1;
   x.anytime = k == 2;
   const std::string tag = std::string( k == 0 ? "" : k == 1 ? "dd " :
                                        "anytime " ) + "w " +
                           std::to_string( w ).substr( 0 , 3 );
   const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
   const FifteenPuzzle d15( ki.board );
   const ZobristDistribution< FifteenPuzzle > z15( d15 );
   check_bounded( d15 , HDAStarMPI< FifteenPuzzle , decltype( z15 ) >(
                         d15 , z15 , x ).solve() , ki.cost , k == 2 ? 1 : w ,
                  "MPI wHDA* " + tag + ", Korf #" +
                  std::to_string( easy_korf[ 0 ] ) );
   check_bounded( rg , HDAStarMPI< ExplicitGraph< int > , decltype( rh ) >(
                        rg , rh , x ).solve() , 418 , k == 2 ? 1 : w ,
                  "MPI wHDA* " + tag + ", Romania" );
   }

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

  // the optimum read from outside prunes as the one given, but it is not a
  // value of the search; a worse one leaves the search its own optimum
  o = BBOptions();
  o.shared_incumbent = [ & ] { return( opt ); };
  const auto x = branch_and_bound( kp , o );
  check( ( x.status == BBStatus::infeasible ) && ! x.node &&
         ( x.bound <= opt + 1e-9 ) &&
         ( x.stats.evaluated == z.stats.evaluated ) ,
         "B&B with the optimum shared" + tag );
  o.shared_incumbent = [ & ] { return( opt / 2 ); };
  const auto y = branch_and_bound( kp , o );
  check( ( y.status == BBStatus::optimal ) && ( y.value == opt ) &&
         y.node && kp.check( *y.node , y.value ) ,
         "B&B with a worse value shared" + tag );

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
// distributed queues, K = 2 the farm, K = 3 the work stealing, K = 4 the
// MultiQueue

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
    if constexpr( K == 4 )
     return( solve( MultiQueueBranchAndBound< Q , Executor >( q , o ) ) );
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
                                       K == 2 ? "farm" :
                                       K == 3 ? "stealing" : "multiqueue" ) +
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
// 8 workers, every selection, a split of 1 and 4 nodes per worker, and the
// sampling with and without its dive: its outcome is the same at every
// repetition, also when a node limit stops it, each worker evaluates its
// nodes by its own evaluator, as many as the counters say, and the values
// given to on_incumbent improve each time and end at the optimum; in both
// senses, its leaves and its limit cases; with the incumbent shared, the
// optimum with every selection and number of workers, and its limit cases;
// with 1 worker, the nodes of the sequential search

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
    for( std::size_t limit : { std::size_t( 0 ) , std::size_t( 7 ) } )
     for( bool dive : { true , false } ) {
      bool ok = true , same = true , own = true;
      for( unsigned p : { 1u , 2u , 3u , 8u } ) {
       ParallelBBOptions o;
       o.selection = sel;
       o.split = split;
       o.split_dive = dive;
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
                            ( limit ? " node limit" : "" ) +
                            ( dive ? "" : " without the dive" );
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

  bool alone = true;
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   ParallelBBOptions a;
   a.selection = sel;
   a.threads = 1;
   BBOptions b;
   b.selection = sel;
   const auto x = run( kp , a );
   const auto y = BranchAndBoundSearch< Knapsack >( kp , b ).solve();
   alone = alone && ( x.value == y.value ) &&
           ( x.stats.evaluated == y.stats.evaluated ) &&
           ( x.stats.generated == y.stats.generated );
   }
  check( alone , name + " with 1 worker evaluates the nodes of the "
                        "sequential search" + tag );

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

  bool shared = true;
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( unsigned p : { 1u , 2u , 3u , 8u } ) {
    o = ParallelBBOptions();
    o.selection = sel;
    o.threads = p;
    o.split_shared = true;
    const auto s = run( kp , o );
    const auto m = run( MinKnapsack{ kp } , o );
    shared = shared && ( s.status == BBStatus::optimal ) &&
             ( s.value == opt ) && s.node && kp.check( *s.node , s.value ) &&
             ( s.gap <= 1e-9 ) && ( m.status == BBStatus::optimal ) &&
             ( m.value == - opt ) && m.node &&
             kp.check( *m.node , - m.value );
    }
  check( shared , name + " incumbent shared" + tag );
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
 o.split_shared = true;
 check_bb_edges( knapsacks().front() , o , run , name + " incumbent shared" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.split_dive = false;
 check_bb_edges( knapsacks().front() , o , run , name + " without the dive" );
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
                     BBScheme::selfsplit , BBScheme::stealing ,
                     BBScheme::multiqueue } )
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
                                 "stealing" , "multiqueue" } )
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
    const auto ml = parallel_branch_and_bound( ki , BBScheme::multiqueue ,
                                               ol );
    bool refused = false;
    try {
     ( void ) parallel_branch_and_bound( ke , BBScheme::multiqueue , oq );
     }
    catch( const std::invalid_argument & ) {
     refused = true;
     }
    check( ( ml.status == BBStatus::optimal ) && ( ml.value == opt ) &&
           refused , "B&B long multiqueue, Rational refused " + t );
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

/*--------------------------------------------------------------------------*/
/*---------------- PERMUTATION PUZZLES AND PARALLEL IDA* -------------------*/
/*--------------------------------------------------------------------------*/
// the distance to a goal of every state reachable from the goals of d, by
// a breadth-first search (the moves of these domains are reversible)

template< class D >
static std::map< state_t< D > , int > distances_from(
 const D & d , const std::vector< state_t< D > > & goals )
{
 std::map< state_t< D > , int > dist;
 std::queue< state_t< D > > q;
 for( const auto & g : goals )
  if( dist.emplace( g , 0 ).second )
   q.push( g );
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  const int ds = dist[ s ];
  d.for_each_successor( s , cost_t< D >( 0 ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > ,
                               cost_t< D > ) {
   if( dist.emplace( t , ds + 1 ).second )
    q.push( t );
   } );
  }
 return( dist );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if h is admissible and consistent on every edge of the states of
// dist, and 0 at the goals

template< class D , class H >
static bool admissible_consistent(
 const D & d , const std::map< state_t< D > , int > & dist , H && h )
{
 for( const auto & [ s , ds ] : dist ) {
  const int hs = int( h( s ) );
  if( ( hs > ds ) || ( ( ds == 0 ) && ( hs != 0 ) ) )
   return( false );
  bool ok = true;
  d.for_each_successor( s , cost_t< D >( 0 ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > c ,
                               cost_t< D > ) {
   ok = ok && ( hs <= int( c ) + int( h( t ) ) );
   } );
  if( ! ok )
   return( false );
  }
 return( true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the heuristic that for_each_successor() gives each successor is
// h of it, on the states of dist (the incremental heuristics)

template< class D >
static bool successor_h( const D & d ,
                         const std::map< state_t< D > , int > & dist )
{
 bool ok = true;
 for( const auto & [ s , ds ] : dist )
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > ,
                               cost_t< D > ht ) {
   ok = ok && ( ht == d.heuristic( t ) );
   } );
 return( ok );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the rankings: bijections onto [ 0 , N ), lexicographic for LexRanking
// and TreeRanking (the same ranks), and the factorials

template< Ranking R >
static bool ranking_ok( unsigned n , unsigned k , bool lex )
{
 const R r( n , k );
 std::vector< std::uint8_t > a( k ) , b( k ) , prev( k );
 for( std::uint64_t x = 0 ; x < r.size() ; ++x ) {
  r.unrank( x , a.data() );
  std::uint64_t seen = 0;
  for( std::uint8_t v : a ) {
   if( ( v >= n ) || ( ( seen >> v ) & 1u ) )
    return( false );
   seen |= std::uint64_t( 1 ) << v;
   }
  if( r.rank( a.data() ) != x )
   return( false );
  if( lex && ( x > 0 ) && ! std::lexicographical_compare(
                              prev.begin() , prev.end() , a.begin() ,
                              a.end() ) )
   return( false );
  prev = a;
  }
 return( true );
 }

static void test_permutation( const std::string & )
{
 check( ( k_permutations( 16 , 6 ) == 5765760 ) &&
        ( k_permutations( 20 , 20 ) == 2432902008176640000ULL ) &&
        ( k_permutations( 21 , 21 ) == 0 ) && ( k_permutations( 3 , 4 ) == 0 ),
        "the numbers of k-permutations" );
 for( unsigned n : { 1u , 4u , 7u } )
  for( unsigned k = 1 ; k <= n ; ++k ) {
   const std::string t = " n " + std::to_string( n ) + " k " +
                         std::to_string( k );
   check( ranking_ok< LexRanking >( n , k , true ) , "lex ranking" + t );
   check( ranking_ok< TreeRanking >( n , k , true ) , "tree ranking" + t );
   check( ranking_ok< MyrvoldRuskeyRanking >( n , k , false ) ,
          "Myrvold-Ruskey ranking" + t );
   }
 // the same lexicographic ranks on 64 items, and the large n of the tree
 const LexRanking l( 64 , 5 );
 const TreeRanking t( 64 , 5 );
 const MyrvoldRuskeyRanking m( 200 , 4 );
 const TreeRanking t2( 200 , 4 );
 std::mt19937_64 gen( 7 );
 bool same = true , back = true;
 for( int i = 0 ; i < 1000 ; ++i ) {
  std::uint8_t a[ 5 ] , b[ 5 ];
  l.unrank( random_below( gen , l.size() ) , a );
  same = same && ( t.rank( a ) == l.rank( a ) );
  const std::uint64_t x = random_below( gen , m.size() );
  m.unrank( x , b );
  back = back && ( m.rank( b ) == x );
  const std::uint64_t y = random_below( gen , t2.size() );
  t2.unrank( y , b );
  back = back && ( t2.rank( b ) == y );
  }
 check( same , "lex and tree ranks agree on 64 items" );
 check( back , "Myrvold-Ruskey and tree on 200 items" );
 bool threw = false;
 try {
  LexRanking bad( 65 , 2 );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "lex ranking refuses 65 items" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the construction of the tables of abstract spaces: unit and integer
// costs against Dijkstra, the unreachable states, the files and the cache

// n states, the goal 0, the edges ( r , r + 1 ) of cost 1 + r % c and
// ( r , 2 r ) of cost c + 1, and the states from n on unreachable
struct ChordSpace {
 std::uint64_t n , extra;
 unsigned c;
 bool chords;
 std::uint64_t size( void ) const { return( n + extra ); }
 template< class F > void for_each_goal( F && f ) const { f( 0 ); }
 template< class F >
 void for_each_predecessor( std::uint64_t r , F && f ) const {
  if( r >= n )
   return;
  if( r + 1 < n )
   f( r + 1 , unsigned( 1 + r % c ) );
  if( r > 0 )
   f( r - 1 , unsigned( 1 + ( r - 1 ) % c ) );
  if( chords && ( 2 * r < n ) && ( r > 0 ) )
   f( 2 * r , c + 1 );
  if( chords && ( r % 2 == 0 ) && ( r > 0 ) )
   f( r / 2 , c + 1 );
  }
 };

static void test_abstraction( const std::string & )
{
 for( unsigned c : { 1u , 3u } ) {
  const ChordSpace sp{ 100000 , 7 , c , true };    // 7 blocks
  const auto h = abstract_distances( sp , 3 );
  std::vector< long > ref( sp.n , -1 );        // Dijkstra
  std::priority_queue< std::pair< long , std::uint64_t > ,
                       std::vector< std::pair< long , std::uint64_t > > ,
                       std::greater<> > q;
  ref[ 0 ] = 0;
  q.emplace( 0 , 0 );
  while( ! q.empty() ) {
   const long dv = q.top().first;
   const std::uint64_t v = q.top().second;
   q.pop();
   if( dv != ref[ v ] )
    continue;
   sp.for_each_predecessor( v , [ & ]( std::uint64_t u , unsigned w ) {
    if( ( ref[ u ] < 0 ) || ( dv + long( w ) < ref[ u ] ) ) {
     ref[ u ] = dv + long( w );
     q.emplace( ref[ u ] , u );
     }
    } );
   }
  bool ok = true;
  for( std::uint64_t r = 0 ; r < sp.n ; ++r )
   ok = ok && ( long( h[ r ] ) == ref[ r ] );
  for( std::uint64_t r = sp.n ; r < sp.size() ; ++r )
   ok = ok && ( h[ r ] == 255 );
  check( ok , "abstract distances against Dijkstra, costs 1 to " +
         std::to_string( c + 1 ) );
  }
 bool threw = false;
 try {
  static_cast< void >( abstract_distances( ChordSpace{ 300 , 0 , 1 , false } ,
                                           2 ) );
  }
 catch( std::logic_error & ) {
  threw = true;
  }
 check( threw , "a distance beyond 254 is refused" );

 const ChordSpace small{ 200 , 5 , 1 , false };
 const auto h = abstract_distances< ADS::ThreadPool >( small , 2 );
 const AbstractionTable< 8 > t8( h );
 const AbstractionTable< 4 > t4( h );
 bool ok = true;
 for( std::uint64_t r = 0 ; r < small.size() ; ++r )
  ok = ok && ( t8[ r ] == h[ r ] ) &&
       ( t4[ r ] == std::min< unsigned >( h[ r ] , 15 ) );
 check( ok && ( 2 * t4.bytes() <= t8.bytes() + 8 ) ,
        "the tables of 8 and 4 bits" );
 const auto dir = std::filesystem::temp_directory_path() /
                  ( "pbefs-abs-" + std::to_string( std::random_device()() ) );
 std::filesystem::create_directories( dir );
 const std::string f = ( dir / "t.pdb" ).string();
 t4.save( f , 42 );
 const auto back = AbstractionTable< 4 >::load( f , 42 , small.size() );
 ok = true;
 for( std::uint64_t r = 0 ; r < small.size() ; ++r )
  ok = ok && ( back[ r ] == t4[ r ] );
 check( ok , "a table saved and loaded" );
 threw = false;
 try {
  static_cast< void >( AbstractionTable< 4 >::load( f , 43 , small.size() ) );
  }
 catch( std::runtime_error & ) {
  threw = true;
  }
 check( threw , "a table is not loaded with another key" );
 int built = 0;
 const auto make = [ & ] {
  ++built;
  return( abstract_distances( small , 1 ) );
  };
 static_cast< void >( AbstractionTable< 8 >::cached( dir.string() , "c.pdb" ,
                                                     7 , small.size() ,
                                                     make ) );
 const auto c2 = AbstractionTable< 8 >::cached( dir.string() , "c.pdb" , 7 ,
                                                small.size() , make );
 check( ( built == 1 ) && ( c2[ 150 ] == 150 ) ,
        "a table built once, then read from the cache" );
 std::filesystem::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the pancakes: the gap heuristics against the exact distances of every
// stack of 7, the searches against them, and the pattern databases

static void test_pancake( const std::string & )
{
 using P7 = PancakePuzzle< 7 >;
 const P7 any( P7::random( 1 ) );
 P7::state_type id{};
 for( int i = 0 ; i < 7 ; ++i )
  id[ i ] = std::uint8_t( i );
 const auto dist = distances_from( any , { id } );
 check( dist.size() == 5040 , "the 7! stacks of 7 pancakes" );
 check( admissible_consistent( any , dist , [ & ]( const P7::state_type & s ) {
         return( any.heuristic( s ) ); } ) && successor_h( any , dist ) ,
        "the gap heuristic is consistent, and incremental" );
 const P7 gap2( P7::random( 1 ) , 2 );
 check( admissible_consistent( gap2 , dist , [ & ]( const P7::state_type & s )
         { return( gap2.heuristic( s ) ); } ) && successor_h( gap2 , dist ) ,
        "GAP-2 is consistent, and incremental" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 97 ) != 0 )
   continue;
  const P7 d( s );
  wrong += ( astar( d ).cost != ds ) || ( idastar( d ).cost != ds );
  }
 check( wrong == 0 , "A* and IDA* on 52 stacks of 7" );

 using B5 = BurnedPancakePuzzle< 5 >;
 const B5 bany( B5::random( 1 ) );
 B5::state_type bid{};
 for( int i = 0 ; i < 5 ; ++i )
  bid[ i ] = std::uint8_t( i );
 const auto bdist = distances_from( bany , { bid } );
 check( bdist.size() == 3840 , "the 5! 2^5 stacks of 5 burned pancakes" );
 check( admissible_consistent( bany , bdist , []( const B5::state_type & s ) {
         return( B5::heuristic( s ) ); } ) && successor_h( bany , bdist ) ,
        "the burned gap heuristic is consistent, and incremental" );
 check( B5::parse( "-3 1 -2 5 4" ) == B5::state_type{ 0x82 , 0 , 0x81 , 4 , 3 }
        , "a burned stack in text" );

 // the tables of 4 of 7 pancakes and of 3 of 5 burned ones, and their
 // maximum, within the exact distances; the files of the cache
 PermutationPDB< LexRanking , 8 > p( any , 0 , { 0 , 1 , 2 , 3 } , 2 );
 PermutationPDB< MyrvoldRuskeyRanking , 8 > q( any , 0 , { 3 , 4 , 5 , 6 } ,
                                              2 );
 check( p.table().size() == 840 , "the table of 4 of 7 pancakes" );
 bool ok = true;
 for( const auto & [ s , ds ] : dist )
  ok = ok && ( int( p.value( any , s ) ) <= ds ) &&
       ( int( q.value( any , s ) ) <= ds );
 check( ok , "the tables of 4 pancakes are admissible" );
 const PermutationHeuristic< P7 > ph( any , { p } );
 const PDBDomain< P7 , PermutationHeuristic< P7 > > pd( any , ph );
 check( admissible_consistent( any , dist , [ & ]( const P7::state_type & s ) {
         return( pd.heuristic( s ) ); } ) ,
        "the maximum of the gap and of a table is consistent" );
 PermutationPDB< TreeRanking , 4 > b( bany , 0 , { 0 , 1 , 2 } , 1 );
 ok = b.table().size() == 60 * 8;
 for( const auto & [ s , ds ] : bdist )
  ok = ok && ( int( b.value( bany , s ) ) <= ds );
 check( ok , "the table of 3 burned pancakes is admissible" );

 const auto dir = std::filesystem::temp_directory_path() /
                  ( "pbefs-perm-" + std::to_string( std::random_device()() ) );
 PermutationPDB< LexRanking , 8 > c1( any , 0 , { 0 , 1 , 2 } , 1 ,
                                      dir.string() );
 PermutationPDB< LexRanking , 8 > c2( any , 0 , { 0 , 1 , 2 } , 1 ,
                                      dir.string() );
 ok = std::distance( std::filesystem::directory_iterator( dir ) ,
                     std::filesystem::directory_iterator() ) == 1;
 for( std::uint64_t r = 0 ; r < c1.table().size() ; ++r )
  ok = ok && ( c1.table()[ r ] == c2.table()[ r ] );
 check( ok , "a table of pancakes in the cache" );
 std::filesystem::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// TopSpin: the breakpoints and a table against the exact distances of
// TopSpin( 8 , 4 ), and the searches

static void test_topspin( const std::string & )
{
 using T = TopSpin< 8 , 4 >;
 const T any( T::random( 1 , 10 ) );
 std::vector< T::state_type > goals;
 T::for_each_goal_permutation( 0 , [ & ]( const std::uint8_t * o ,
                                          const std::uint8_t * ) {
  T::state_type g;
  std::copy( o , o + 8 , g.begin() );
  goals.push_back( g );
  } );
 const auto dist = distances_from( any , goals );
 check( dist.size() == 40320 , "the 8! rings of TopSpin( 8 , 4 )" );
 check( admissible_consistent( any , dist , []( const T::state_type & s ) {
         return( T::heuristic( s ) ); } ) ,
        "the breakpoints of TopSpin are consistent" );
 PermutationPDB< LexRanking , 8 > p( any , 0 , { 0 , 1 , 2 , 3 , 4 } , 2 );
 const PermutationHeuristic< T > ph( any , { p } );
 const PDBDomain< T , PermutationHeuristic< T > > pd( any , ph );
 check( admissible_consistent( any , dist , [ & ]( const T::state_type & s ) {
         return( pd.heuristic( s ) ); } ) ,
        "a table of 5 tokens of TopSpin is consistent" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 2017 ) != 0 )
   continue;
  const T d( s );
  const PDBDomain< T , PermutationHeuristic< T > > e( d , ph );
  wrong += ( idastar( d ).cost != ds ) || ( idastar( e ).cost != ds ) ||
           ( astar( e ).cost != ds );
  }
 check( wrong == 0 , "IDA* and A* on 20 rings of TopSpin( 8 , 4 )" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// Rubik's cube: the size of the tree of table 1 of Korf (1997), the states
// at each depth, the instances of Korf, and small tables

static void test_rubik( const std::string & data )
{
 const RubiksCube c;
 std::vector< std::uint64_t > tree( 6 , 0 );
 std::function< void( const RubiksCube::state_type & , std::uint8_t , int )
                > down = [ & ]( const RubiksCube::state_type & s ,
                                std::uint8_t op , int depth ) {
  if( depth == 5 )
   return;
  c.for_each_tree_successor( s , 0 , op , [ & ](
   const RubiksCube::state_type & t , int , int , std::uint8_t o ) {
   ++tree[ std::size_t( depth + 1 ) ];
   down( t , o , depth + 1 );
   } );
  };
 down( c.initial_state() , no_op , 0 );
 check( tree == std::vector< std::uint64_t >{ 0 , 18 , 243 , 3240 , 43254 ,
                                              577368 } ,
        "the tree of the cube, table 1 of Korf (1997), to depth 5" );
 std::set< RubiksCube::state_type > seen{ c.initial_state() } ,
                                    front{ c.initial_state() };
 std::vector< std::size_t > layer;
 for( int depth = 1 ; depth <= 4 ; ++depth ) {
  std::set< RubiksCube::state_type > next;
  for( const auto & s : front )
   c.for_each_successor( s , 0 , [ & ]( const RubiksCube::state_type & t ,
                                        int , int ) {
    if( seen.insert( t ).second )
     next.insert( t );
    } );
  layer.push_back( next.size() );
  front.swap( next );
  }
 check( layer == std::vector< std::size_t >{ 18 , 243 , 3240 , 43239 } ,
        "the states of the cube at the depths 1 to 4" );

 // every move has order 4 (a quarter turn) or 2, and its inverse
 bool ok = true;
 const auto s0 = RubiksCube::random_walk( 30 , 3 );
 for( int m = 0 ; m < 18 ; ++m ) {
  auto s = s0;
  for( int i = 0 ; i < 4 ; ++i )
   s = RubiksCube::apply( s , m );
  ok = ok && ( s == s0 ) &&
       ( RubiksCube::apply( RubiksCube::apply( s0 , m ) ,
                            3 * ( m / 3 ) + 2 - m % 3 ) == s0 );
  }
 check( ok , "the moves of the cube and their inverses" );

 std::ifstream in( data + "/rubik/korf10.txt" );
 int n = 0;
 for( std::string line ; std::getline( in , line ) ; )
  if( ! line.empty() ) {
   ++n;
   ok = ok && ( RubiksCube::parse_moves( line ).size() == 100 );
   }
 check( ok && ( n == 10 ) , "the 10 instances of Korf, of 100 moves" );

 // small tables: 4 corners and 4 + 4 edges, against IDA* without them
 std::vector< PermutationPDB< LexRanking , 4 > > v;
 v.emplace_back( c , 0 , std::vector< std::uint8_t >{ 0 , 1 , 2 , 3 } , 0 );
 v.emplace_back( c , 1 , std::vector< std::uint8_t >{ 0 , 1 , 2 , 3 } , 0 );
 v.emplace_back( c , 1 , std::vector< std::uint8_t >{ 4 , 5 , 6 , 7 } , 0 );
 check( ( v[ 0 ].table().size() == 1680 * 81 ) &&
        ( v[ 1 ].table().size() == 11880 * 16 ) ,
        "the sizes of the tables of 4 corners and of 4 edges" );
 const PermutationHeuristic< RubiksCube , LexRanking , 4 > h( c , v );
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const RubiksCube d( RubiksCube::random_walk( 5 , seed ) );
  const PDBDomain< RubiksCube , decltype( h ) > e( d , h );
  const auto ref = idastar( d );
  check_solution( e , idastar( e ) , ref.cost , "IDA* with tables, walk " +
                  std::to_string( seed ) );
  PIDAOptions o;
  o.threads = 3;
  check_solution( e , ParallelIDAStar< decltype( e ) >( e , o ).solve() ,
                  ref.cost , "PIDA* with tables, walk " +
                  std::to_string( seed ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the tables of Korf (1997) for Rubik's cube: their sizes, largest values
// and means, against those of his paper (88179840 entries of mean 8.764
// and largest value 11; 42577920 of mean 7.668 and largest 10)

static void test_rubik_pdb( const std::string & )
{
 #if ! defined( NDEBUG ) || ( defined( ADS_CONTRACTS ) && ADS_CONTRACTS )
  check( true , "the tables of Korf only in an optimized build without the "
         "checks of ADS" );
  return;
 #endif
 const RubiksCube c;
 const RubikPDBHeuristic h( c );
 const std::uint64_t size[] = { 88179840 , 42577920 , 42577920 };
 const unsigned top[] = { 11 , 10 , 10 };
 const double mean[] = { 8.764 , 7.668 , 7.668 };
 for( std::size_t i = 0 ; i < 3 ; ++i ) {
  const auto hist = h.pdbs()[ i ].table().histogram();
  double sum = 0;
  unsigned largest = 0;
  for( std::size_t x = 0 ; x < hist.size() ; ++x )
   if( hist[ x ] ) {
    sum += double( x ) * double( hist[ x ] );
    largest = unsigned( x );
    }
  const double m = sum / double( size[ i ] );
  std::ostringstream os;
  os << "the table " << i << " of Korf: " << h.pdbs()[ i ].table().size()
     << " entries, largest " << largest << ", mean " << m;
  check( ( h.pdbs()[ i ].table().size() == size[ i ] ) &&
         ( largest == top[ i ] ) && ( std::abs( m - mean[ i ] ) < 5e-4 ) ,
         os.str() );
  }
 for( std::uint64_t seed = 1 ; seed <= 2 ; ++seed ) {
  const RubiksCube d( RubiksCube::random_walk( 10 , seed ) );
  const PDBDomain< RubiksCube , RubikPDBHeuristic > e( d , h );
  const auto r = idastar( e );
  PIDAOptions o;
  o.threads = 4;
  check_solution( e , ParallelIDAStar< decltype( e ) >( e , o ).solve() ,
                  r.cost , "PIDA* with the tables of Korf, walk of 10, seed "
                  + std::to_string( seed ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the Towers of Hanoi with four pegs: the instances of Siag et al., the
// tables against the exact distances of 8 disks, and the searches

static void test_hanoi( const std::string & data )
{
 std::ifstream in( data + "/hanoi/toh20.txt" );
 int n = 0;
 bool ok = true;
 for( std::string line ; std::getline( in , line ) ; ) {
  if( line.empty() || ( line[ 0 ] != '(' ) )
   continue;
  std::istringstream is( line );
  const auto s = TowersOfHanoi< 20 >::read( is );
  const auto g = TowersOfHanoi< 20 >::read( is );
  while( ! line.empty() &&
         std::isspace( static_cast< unsigned char >( line.back() ) ) )
   line.pop_back();
  ok = ok && ( TowersOfHanoi< 20 >::write( s ) + " " +
               TowersOfHanoi< 20 >::write( g ) == line );
  ++n;
  }
 check( ok && ( n == 20 ) , "the 20 instances of Siag et al., read and "
        "written back" );

 using H8 = TowersOfHanoi< 8 >;
 std::mt19937_64 gen( 5 );
 std::uint64_t goal = 0;
 for( int k = 0 ; k < 8 ; ++k )
  goal |= random_below( gen , 4 ) << ( 2 * k );
 const H8 any( 0 , goal );
 const auto dist = distances_from( any , { goal } );
 check( dist.size() == 65536 , "the 4^8 states of 8 disks" );
 const HanoiPDB p( { 3 , 5 } , goal , 2 );
 const HanoiPDB q( { 8 } , goal , 2 );
 ok = true;
 for( const auto & [ s , ds ] : dist )
  ok = ok && ( int( q( s ) ) == ds );
 check( ok , "the table of all the disks is the exact distance" );
 check( admissible_consistent( any , dist , [ & ]( std::uint64_t s ) {
         return( p( s ) ); } ) , "the tables of 3 + 5 disks are additive "
        "and consistent" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 4099 ) != 0 )
   continue;
  const H8 d( s , goal );
  const PDBDomain< H8 , HanoiPDB > e( d , p );
  PIDAOptions o;
  o.threads = 2;
  o.frontier = 256;
  o.table_depth = 8;
  wrong += ( astar( e ).cost != ds ) ||
           ( ParallelIDAStar< decltype( e ) >( e , o ).solve().cost != ds );
  }
 check( wrong == 0 , "A* and PIDA* on 16 instances of 8 disks" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the road graphs of DIMACS: a graph and its coordinates in the format of
// the challenge, the consistency of the heuristic, A* against Dijkstra

static void test_dimacs( const std::string & data )
{
 std::ifstream gr( data + "/dimacs/sample.gr" ) ,
               co( data + "/dimacs/sample.co" );
 check( gr && co , "the sample graph is read" );
 auto g = read_dimacs_graph< std::int64_t >( gr , 1 , 2 );
 const auto c = read_dimacs_coordinates( co , g.size() );
 const double kappa = dimacs_scale( g , c );
 check( kappa > 0 , "the factor of the heuristic is positive" );
 std::mt19937_64 gen( 3 );
 for( int i = 0 ; i < 5 ; ++i ) {
  g.set_source( std::uint32_t( random_below( gen , g.size() ) ) );
  g.set_target( std::uint32_t( random_below( gen , g.size() ) ) );
  set_dimacs_heuristic( g , c , kappa );
  bool ok = g.heuristic( g.target() ) == 0;
  for( std::uint32_t u = 0 ; u < g.size() ; ++u )
   g.for_each_successor( u , 0 , [ & ]( std::uint32_t , std::int64_t w ,
                                        std::int64_t hv ) {
    ok = ok && ( g.heuristic( u ) <= w + hv );
    } );
  check( ok , "the heuristic of the coordinates is consistent, pair " +
         std::to_string( i ) );
  const auto ref = astar( Blind< ExplicitGraph< std::int64_t > >{ g } ,
                          OpenKind::heap );
  const auto a = astar( g , OpenKind::heap );
  if( ref.solved() )
   check_solution( g , a , ref.cost , "A* on the sample, pair " +
                   std::to_string( i ) );
  else
   check( ! a.solved() , "A* on the sample, pair " + std::to_string( i ) );
  }
 std::istringstream bad( "p sp 2 1\na 1 3 5\n" );
 bool threw = false;
 try {
  static_cast< void >( read_dimacs_graph( bad ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "an arc out of range is refused" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// ParallelIDAStar with every frontier rule, schedule and scheme, on the
// given runtime: the 8-puzzle against its exact distances, the 15-puzzle,
// the pancakes, the map of Romania, the node limit, the unsolvable graphs,
// IDA* reproduced with one unit and one worker, the deterministic path

template< class Executor >
static void test_pida( const std::string & data )
{
 const auto inst = korf( data );
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 // the variants: rule, schedule, unit bounds, scheme, split, table depth
 struct Variant {
  FrontierRule rule;
  UnitSchedule schedule;
  bool unit_bounds , asynchronous;
  double split;
  unsigned table;
  const char * name;
  };
 const Variant variants[] = {
  { FrontierRule::astar , UnitSchedule::steal , true , false , 0 , 0 ,
    "astar steal" } ,
  { FrontierRule::breadth , UnitSchedule::shared , true , false , 0 , 0 ,
    "breadth shared" } ,
  { FrontierRule::depth , UnitSchedule::shared , false , false , 0 , 0 ,
    "depth classical" } ,
  { FrontierRule::astar , UnitSchedule::steal , true , true , 0 , 0 ,
    "astar asynchronous" } ,
  { FrontierRule::depth , UnitSchedule::steal , true , false , 2 , 0 ,
    "depth split" } ,
  { FrontierRule::astar , UnitSchedule::shared , true , true , 2 , 3 ,
    "astar asynchronous split table" } ,
  { FrontierRule::breadth , UnitSchedule::steal , true , false , 0 , 4 ,
    "breadth table" } };
 const auto options = []( const Variant & v , unsigned p ) {
  PIDAOptions o;
  o.threads = p;
  o.rule = v.rule;
  o.schedule = v.schedule;
  o.unit_bounds = v.unit_bounds;
  o.asynchronous = v.asynchronous;
  o.split = v.split;
  o.table_depth = v.table;
  return( o );
  };
 for( const Variant & v : variants )
  for( unsigned p : { 1u , 2u , 4u } ) {
   const PIDAOptions o = options( v , p );
   const std::string tag = std::string( " " ) + v.name + ", p " +
                           std::to_string( p );
   std::size_t k = 0;
   int wrong = 0;
   for( const auto & [ s , dst ] : dist ) {
    if( ( k++ % 18143 ) != 0 )                  // 10 states
     continue;
    const EightPuzzle e( EightPuzzle::unpack( s ) );
    const auto r = ParallelIDAStar< EightPuzzle , Executor >( e , o )
                   .solve();
    if( ! r.solved() || ( r.cost != dst ) ) {
     ++wrong;
     check_solution( e , r , dst , "PIDA* 8-puzzle" + tag );
     }
    }
   check( wrong == 0 , "PIDA* on 10 states of the 8-puzzle" + tag );
   const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
   const FifteenPuzzle d( ki.board );
   check_solution( d , ParallelIDAStar< FifteenPuzzle , Executor >( d , o )
                        .solve() , ki.cost , "PIDA* Korf #" +
                   std::to_string( easy_korf[ 0 ] ) + tag );
   const PancakePuzzle< 12 > pc( PancakePuzzle< 12 >::random( p ) );
   check_solution( pc , ParallelIDAStar< PancakePuzzle< 12 > , Executor >(
                         pc , o ).solve() , astar( pc ).cost ,
                   "PIDA* pancakes of 12" + tag );
   if( v.table == 0 ) {
    auto g = romania( Arad );
    check_solution( g , ParallelIDAStar< ExplicitGraph< int > , Executor >(
                         g , o ).solve() , 418 , "PIDA* Romania" + tag );
    }
   }

 // one unit and one worker: the expansions of IDA*
 for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
  const FifteenPuzzle d( inst[ id - 1 ].board );
  PIDAOptions o;
  o.frontier = 1;
  const auto r = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  const auto i = idastar( d );
  check( ( r.cost == i.cost ) && ( r.stats.expanded == i.stats.expanded ) &&
         ( r.stats.generated == i.stats.generated ) && ( r.path == i.path ) ,
         "PIDA* with one unit is IDA*, Korf #" + std::to_string( id ) );
  }

 // the deterministic path, the same with 1 and 3 workers
 {
  const FifteenPuzzle d( inst[ easy_korf[ 2 ] - 1 ].board );
  PIDAOptions o;
  o.deterministic = true;
  o.frontier = 64;
  const auto a = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  o.threads = 3;
  const auto b = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  check( a.solved() && ( a.path == b.path ) && ( a.cost == b.cost ) ,
         "PIDA* deterministic: the same path with 1 and 3 workers" );
  bool threw = false;
  o.asynchronous = true;
  try {
   ParallelIDAStar< FifteenPuzzle , Executor > x( d , o );
   }
  catch( std::invalid_argument & ) {
   threw = true;
   }
  check( threw , "PIDA* refuses deterministic and asynchronous" );
  }

 // the node limit, a goal in the frontier, the unsolvable graphs
 {
  const FifteenPuzzle d( inst[ 0 ].board );
  PIDAOptions o;
  o.threads = 3;
  o.node_limit = 1000;
  const auto r = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  check( ( r.status == Status::node_limit ) && r.path.empty() ,
         "PIDA* stops at the expansion limit" );
  const FifteenPuzzle e( random_walk< FifteenPuzzle >( 6 , 2 ) );
  PIDAOptions f;
  f.frontier = 100000;
  check_solution( e , ParallelIDAStar< FifteenPuzzle , Executor >( e , f )
                       .solve() , idastar( e ).cost ,
                  "PIDA* whose frontier reaches the goal" );
  ExplicitGraph< int > dag( 5 , 0 , 4 );
  dag.add_edge( 0 , 1 , 1 , true );
  dag.add_edge( 1 , 2 , 1 , true );
  dag.add_edge( 0 , 3 , 2 , true );
  for( FrontierRule rule : { FrontierRule::astar , FrontierRule::depth } )
   for( bool async : { false , true } ) {
    PIDAOptions u;
    u.threads = 2;
    u.frontier = 2;
    u.rule = rule;
    u.asynchronous = async;
    check( ParallelIDAStar< ExplicitGraph< int > , Executor >( dag , u )
           .solve().status == Status::unsolvable ,
           "PIDA* proves an acyclic graph unsolvable" );
    }
  }
 }

#ifdef PBEFS_HAS_MPI

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// ParallelIDAStarMPI on the 15-puzzle, the pancakes, the Towers of Hanoi
// and Romania, with every frontier rule

static void test_pida_mpi( const std::string & data )
{
 const auto inst = korf( data );
 for( FrontierRule rule : { FrontierRule::astar , FrontierRule::breadth ,
                            FrontierRule::depth } ) {
  PIDAOptions o;
  o.rule = rule;
  const std::string tag = rule == FrontierRule::astar ? " astar"
                          : rule == FrontierRule::breadth ? " breadth"
                                                          : " depth";
  for( int id : { easy_korf[ 0 ] , easy_korf[ 2 ] } ) {
   const FifteenPuzzle d( inst[ id - 1 ].board );
   check_solution( d , ParallelIDAStarMPI< FifteenPuzzle >( d , o ).solve() ,
                   inst[ id - 1 ].cost , "MPI PIDA* Korf #" +
                   std::to_string( id ) + tag );
   }
  const PancakePuzzle< 14 > pc( PancakePuzzle< 14 >::random( 3 ) );
  check_solution( pc , ParallelIDAStarMPI< PancakePuzzle< 14 > >( pc , o )
                        .solve() , astar( pc ).cost , "MPI PIDA* pancakes"
                  + tag );
  auto g = romania( Arad );
  check_solution( g , ParallelIDAStarMPI< ExplicitGraph< int > >( g , o )
                       .solve() , 418 , "MPI PIDA* Romania" + tag );
  }
 PIDAOptions l;
 l.node_limit = 500;
 const FifteenPuzzle d( inst[ 0 ].board );
 check( ParallelIDAStarMPI< FifteenPuzzle >( d , l ).solve().status ==
        Status::node_limit , "MPI PIDA* stops at the expansion limit" );
 ExplicitGraph< int > dag( 4 , 0 , 3 );
 dag.add_edge( 0 , 1 , 1 , true );
 dag.add_edge( 1 , 2 , 1 , true );
 check( ParallelIDAStarMPI< ExplicitGraph< int > >( dag ).solve().status ==
        Status::unsolvable , "MPI PIDA* proves an acyclic graph unsolvable" );
 }

#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the skeletons, the modes, the ramp-up, the checkpoints, the
// intervals, the hubs and the domains of the B&B [see bnb.cpp]

void test_bb_skeleton_thread( const std::string & );
void test_bb_interval_thread( const std::string & );
void test_bb_modes( const std::string & );
void test_bb_race( const std::string & );
void test_bb_checkpoint( const std::string & );
void test_bb_domains( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 void test_bb_skeleton_openmp( const std::string & );
 void test_bb_interval_openmp( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 void test_bb_skeleton_fastflow( const std::string & );
 void test_bb_interval_fastflow( const std::string & );
#endif
#ifdef PBEFS_HAS_MPI
 void test_bb_hub_mpi( const std::string & );
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
 { "ddd" , "-" , test_ddd } ,
 { "tiles24" , "-" , test_tiles24 } ,
 { "graph" , "-" , []( const std::string & ) { test_graph(); } } ,
 { "limits" , "-" , test_limits } ,
 { "reopen" , "-" , []( const std::string & ) { test_reopen(); } } ,
 { "msa" , "-" , test_msa } ,
 { "board" , "-" , test_board } ,
 { "files" , "-" , test_files } ,
 { "pdb" , "-" , test_pdb } ,
 { "pdb-search" , "-" , test_pdb_search } ,
 { "weighted" , "-" , test_weighted } ,
 { "anytime" , "-" , test_anytime } ,
 { "focal" , "-" , test_focal } ,
 { "lpa" , "-" , test_lpa } ,
 { "bounded-thread" , "-" , test_bounded< ADS::ThreadPool > } ,
 { "hda-thread" , "-" , test_hda< ADS::ThreadPool > } ,
 { "spa-thread" , "-" , test_spa< ADS::ThreadPool > } ,
 { "mqa-thread" , "-" , test_mqa< ADS::ThreadPool > } ,
 { "pbnf-thread" , "-" , test_pbnf< ADS::ThreadPool > } ,
 { "pase-thread" , "-" , test_pase< ADS::ThreadPool > } ,
 { "pida-thread" , "-" , test_pida< ADS::ThreadPool > } ,
 { "permutation" , "-" , test_permutation } ,
 { "abstraction" , "-" , test_abstraction } ,
 { "pancake" , "-" , test_pancake } ,
 { "topspin" , "-" , test_topspin } ,
 { "rubik" , "-" , test_rubik } ,
 { "rubik-pdb" , "-" , test_rubik_pdb } ,
 { "hanoi" , "-" , test_hanoi } ,
 { "dimacs" , "-" , test_dimacs } ,
 { "bb" , "-" , test_bb } ,
 { "bb-shared-thread" , "-" , test_bb_parallel< 0 , ADS::ThreadPool > } ,
 { "bb-dist-thread" , "-" , test_bb_parallel< 1 , ADS::ThreadPool > } ,
 { "bb-steal-thread" , "-" , test_bb_parallel< 3 , ADS::ThreadPool > } ,
 { "bb-mq-thread" , "-" , test_bb_parallel< 4 , ADS::ThreadPool > } ,
 { "bb-det-thread" , "-" , test_bb_det< ADS::ThreadPool > } ,
 { "bb-split-thread" , "-" , test_bb_split< ADS::ThreadPool > } ,
 { "bb-farm-thread" , "-" , test_bb_parallel< 2 , ADS::ThreadPool > } ,
 { "farm-thread" , "-" , test_farm< ADS::ThreadPool > } ,
 { "bb-scheme" , "-" , test_bb_scheme } ,
 { "bb-cutoff" , "-" , test_bb_cutoff } ,
 { "bb-values" , "-" , test_bb_values } ,
 { "bb-skel-thread" , "-" , test_bb_skeleton_thread } ,
 { "bb-ivm-thread" , "-" , test_bb_interval_thread } ,
 { "bb-modes" , "-" , test_bb_modes } ,
 { "bb-race" , "-" , test_bb_race } ,
 { "bb-checkpoint" , "-" , test_bb_checkpoint } ,
 { "bb-domains" , "-" , test_bb_domains } ,
 #ifdef PBEFS_HAS_OPENMP
  { "hda-openmp" , "openmp" , test_hda< ADS::OpenMPExecutor > } ,
  { "spa-openmp" , "openmp" , test_spa< ADS::OpenMPExecutor > } ,
  { "mqa-openmp" , "openmp" , test_mqa< ADS::OpenMPExecutor > } ,
  { "pbnf-openmp" , "openmp" , test_pbnf< ADS::OpenMPExecutor > } ,
  { "bb-shared-openmp" , "openmp" ,
    test_bb_parallel< 0 , ADS::OpenMPExecutor > } ,
  { "bb-dist-openmp" , "openmp" ,
    test_bb_parallel< 1 , ADS::OpenMPExecutor > } ,
  { "bb-steal-openmp" , "openmp" ,
    test_bb_parallel< 3 , ADS::OpenMPExecutor > } ,
  { "bb-mq-openmp" , "openmp" ,
    test_bb_parallel< 4 , ADS::OpenMPExecutor > } ,
  { "bb-det-openmp" , "openmp" , test_bb_det< ADS::OpenMPExecutor > } ,
  { "bb-split-openmp" , "openmp" , test_bb_split< ADS::OpenMPExecutor > } ,
  { "bb-skel-openmp" , "openmp" , test_bb_skeleton_openmp } ,
  { "bb-ivm-openmp" , "openmp" , test_bb_interval_openmp } ,
  { "bb-farm-openmp" , "openmp" ,
    test_bb_parallel< 2 , ADS::OpenMPExecutor > } ,
  { "farm-openmp" , "openmp" , test_farm< ADS::OpenMPExecutor > } ,
  { "bounded-openmp" , "openmp" , test_bounded< ADS::OpenMPExecutor > } ,
  { "pida-openmp" , "openmp" , test_pida< ADS::OpenMPExecutor > } ,
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  { "hda-fastflow" , "fastflow" , test_hda< ADS::FastFlowExecutor > } ,
  { "farm-fastflow" , "fastflow" , test_farm< ADS::FastFlowExecutor > } ,
  { "bounded-fastflow" , "fastflow" ,
    test_bounded< ADS::FastFlowExecutor , false > } ,
  { "pida-fastflow" , "fastflow" , test_pida< ADS::FastFlowExecutor > } ,
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
  { "bb-skel-fastflow" , "fastflow" , test_bb_skeleton_fastflow } ,
  { "bb-ivm-fastflow" , "fastflow" , test_bb_interval_fastflow } ,
 #endif
 #ifdef PBEFS_HAS_MPI
  { "mpi" , "mpi" , test_mpi } ,
  { "bb-mpi" , "mpi" , test_bb_mpi } ,
  { "pida-mpi" , "mpi" , test_pida_mpi } ,
  { "bb-hub-mpi" , "mpi" , test_bb_hub_mpi } ,
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
