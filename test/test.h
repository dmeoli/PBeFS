/*--------------------------------------------------------------------------*/
/*------------------------------ File test.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The parts of the tester of PBeFS that its translation units share
 * [see test.cpp]: the includes, the checks, the instances and their
 * references, and the declarations of the cases, each defined in the
 * file of its family.
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

#ifndef __PBeFS_test
 #define __PBeFS_test
                      /* self-identification: #endif at the end of the file */

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

#include "PBeFS/astar/IDDAStar.h"

#include "PBeFS/astar/BiBestFirst.h"

#include "PBeFS/astar/NBS.h"

#include "PBeFS/astar/DVCBS.h"

#include "PBeFS/astar/GBFHS.h"

#include "PBeFS/astar/ParallelBidirectional.h"

#include "PBeFS/astar/PEMSearch.h"

#include "PBeFS/astar/SPAStar.h"
#include "PBeFS/astar/MQAStar.h"

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/ParallelIDAStar.h"

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/EPASE.h"

#include "PBeFS/astar/MPLP.h"

#include "PBeFS/astar/PWAStar.h"

#include "PBeFS/astar/RobotNav2D.h"

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

#include "PBeFS/astar/GreedySearch.h"

#include "PBeFS/astar/EnforcedHillClimbing.h"

#include "PBeFS/astar/IteratedSearch.h"

#include "PBeFS/astar/WidthSearch.h"

#include "PBeFS/astar/ParallelGreedy.h"

#include "PBeFS/astar/PartialExpansion.h"

#include "PBeFS/astar/FringeSearch.h"

#include "PBeFS/astar/RBFS.h"

#include "PBeFS/astar/BeamSearch.h"

#include "PBeFS/astar/RealTimeSearch.h"

#include "PBeFS/astar/HPAStar.h"

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

extern int failures;                  // the checks failed so far
extern bool quiet;                    // print only the failures

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// records the outcome of a check

void check( bool ok , const std::string & what );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the name of the executor in the names of the checks

template< class Executor >
inline std::string executor_name( void )
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
inline bool same_cost( Cost a , Cost b )
{
 if constexpr( std::is_floating_point_v< Cost > )
  return( std::abs( a - b ) <= 1e-9 * std::max( Cost( 1 ) , std::abs( b ) ) );
 else
  return( a == b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that sol is solved with cost expected and a valid path

template< class D , class S >
inline void check_solution( const D & d , const S & sol ,
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

inline std::vector< KorfInstance > korf( const std::string & dir )
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

inline const std::vector< int > easy_korf = { 55 , 16 , 42 , 79 , 71 , 12 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the map of Romania of Russell and Norvig, from a city to Bucharest

enum City { Arad , Bucharest , Craiova , Drobeta , Eforie , Fagaras ,
            Giurgiu , Hirsova , Iasi , Lugoj , Mehadia , Neamt , Oradea ,
            Pitesti , RimnicuVilcea , Sibiu , Timisoara , Urziceni , Vaslui ,
            Zerind , n_cities };

inline ExplicitGraph< int > romania( City from )
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

inline int romania_distance( City from )
{
 auto g = romania( from );
 return( astar( Blind< ExplicitGraph< int > >{ g } ).cost );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 8-puzzle with the blank moved from the goal through the given cells,
// each next to the one before

inline EightPuzzle::board_type walk8( std::initializer_list< int > cells )
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
inline void check_astar_edges( const O & base , Solve && solve ,
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

inline std::unordered_map< std::uint64_t , int > eight_distances(
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

inline std::vector< std::pair< std::uint64_t , int > > eight_sample(
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
/*-------------------------- SHARED BY THE CASES ---------------------------*/
/*--------------------------------------------------------------------------*/
// the cost of a cheapest path on the map m from ( sx , sy ) to ( gx , gy )
// with 4 or 8 moves, -1 if there is none: Dijkstra on the map itself, with
// the rules of Grid written again (leaving a cell costs its cost, times
// sqrt( 2 ) to a corner, whose two cut cells must be free), so that it does
// not depend on the adjacency of the set of data structures

inline double grid_dijkstra( const GridMap & m , int moves , int sx , int sy ,
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
// the board of D reached from the goal of Korf by n random moves of the
// blank, none of which undoes the one before

template< class D >
inline typename D::board_type random_walk( int n , std::uint64_t seed )
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
inline std::unordered_map< state_t< D > , cost_t< D > , StateHash< D > >
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
inline void check_consistent_astar( const D & d , const G & g , const S & sol ,
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the sequences of a FASTA file of the data

inline std::vector< std::string > fasta( const std::string & file )
{
 std::ifstream in( file );
 if( ! in )
  throw( std::runtime_error( "cannot read " + file ) );
 return( SequenceAlignment<>::read_fasta( in ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cost of a cheapest path of a GridBoard by Dijkstra's algorithm,
// written from the rules of the format: leaving the row y costs 1 (unit)
// or y (life), times sqrt( 2 ) to a corner, and a move to a corner never
// looks at the cells it cuts

inline double board_dijkstra( const GridBoard & b )
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

inline GridBoard random_board( int w , int h , double p , std::uint64_t seed ,
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
// the distances of the placements of the tiles pat of the 8-puzzle to
// goal, by a breadth-first search of the abstraction written apart from
// the library: on the placements alone (the blank ignored), or by a 0-1
// search on the placements with the cell of the blank, keeping the least
// distance over the cells of the blank (the blank minimized)

inline std::map< std::vector< int > , int > abstract_distances(
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

inline std::array< std::uint8_t , 9 > cells8( std::uint64_t s )
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
inline std::size_t check_heuristic8(
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that sol is solved with a valid path whose cost is sol.cost, and
// at most bound if given

template< class D , class S >
inline void check_path( const D & d , const S & sol ,
                        const std::string & what ,
                        double bound = infinite_cost< double >() )
{
 bool ok = sol.solved() && ! sol.path.empty() &&
           ( sol.path.front() == d.initial_state() ) &&
           d.is_goal( sol.path.back() );
 cost_t< D > sum = 0;
 for( std::size_t i = 1 ; ok && ( i < sol.path.size() ) ; ++i ) {
  bool edge = false;
  cost_t< D > best = infinite_cost< cost_t< D > >();
  d.for_each_successor( sol.path[ i - 1 ] , d.heuristic( sol.path[ i - 1 ] ) ,
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
 ok = ok && same_cost( sum , sol.cost ) &&
      ( double( sol.cost ) <= bound * ( 1 + 1e-9 ) );
 std::ostringstream os;
 os << what << " (cost " << sol.cost << ", expanded " << sol.stats.expanded
    << ")";
 check( ok , os.str() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- THE CASES --------------------------------*/
/*--------------------------------------------------------------------------*/
// each defined in the file of its family, the templates instantiated
// there for the executors of the table of test.cpp

// [see search.cpp]
void test_astar( const std::string & );
void test_open( void );
void test_grid( void );
void test_graph( void );
void test_limits( const std::string & );
void test_idastar( const std::string & );
void test_tiles24( const std::string & );
void test_ddd( const std::string & );

// [see files.cpp]
void test_msa( const std::string & );
void test_board( const std::string & );
void test_files( const std::string & );
void test_reopen( void );

// [see pdb.cpp]
void test_pdb( const std::string & );
void test_pdb_search( const std::string & );

// [see parallel.cpp]
template< class Executor >
void test_hda( const std::string & );
template< class Executor >
void test_spa( const std::string & );
template< class Executor >
void test_mqa( const std::string & );
template< class Executor >
void test_pbnf( const std::string & );
template< class Executor >
void test_pase( const std::string & );

// [see edge.cpp]
void test_edge( const std::string & );
template< class Executor >
void test_edge_parallel( const std::string & );
template< class Executor >
void test_farm( const std::string & );

// [see bounded.cpp]
void test_weighted( const std::string & );
void test_anytime( const std::string & );
void test_focal( const std::string & );
void test_lpa( const std::string & );
template< class Executor , bool Shared = true >
void test_bounded( const std::string & );
void test_mpi( const std::string & );

// [see bb.cpp]
void test_bb( const std::string & );
template< int K , class Executor >
void test_bb_parallel( const std::string & );
template< class Executor >
void test_bb_det( const std::string & );
template< class Executor >
void test_bb_split( const std::string & );
void test_bb_scheme( const std::string & );
void test_bb_values( const std::string & );
void test_bb_cutoff( const std::string & );
void test_bb_mpi( const std::string & );

// [see permutation.cpp]
void test_permutation( const std::string & );
void test_abstraction( const std::string & );
void test_pancake( const std::string & );
void test_topspin( const std::string & );
void test_rubik( const std::string & );
void test_rubik_pdb( const std::string & );
void test_hanoi( const std::string & );
void test_dimacs( const std::string & );
template< class Executor >
void test_pida( const std::string & );
void test_pida_mpi( const std::string & );

// [see bidirectional.cpp]
void test_bidirectional( const std::string & );

// [see pbidirectional.cpp]
template< class Executor >
void test_bidirectional_parallel( const std::string & );

// [see pem.cpp]
template< class Executor >
void test_pem( const std::string & );
void test_idd( const std::string & );

// [see greedy.cpp]
void test_gbfs( const std::string & );
void test_gbfs_queues( const std::string & );
void test_ehc( const std::string & );
void test_iterated( const std::string & );
void test_width( const std::string & );
void test_classic( const std::string & );
void test_realtime( const std::string & );
void test_hpa( const std::string & );

// [see pgreedy.cpp]
void test_pgbfs_rounds( const std::string & );
template< class Executor >
void test_pgbfs( const std::string & );

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_test */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test.h -----------------------------*/
/*--------------------------------------------------------------------------*/
