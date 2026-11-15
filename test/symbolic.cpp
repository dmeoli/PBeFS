/*--------------------------------------------------------------------------*/
/*--------------------------- File symbolic.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the symbolic
 * searches, compiled only if Sylvan is found. The references do not come
 * from the symbolic searches: a breadth-first search of the whole state
 * space of the 8-puzzle and of the pancakes, the optimal costs of the
 * instances of Korf (1985), and A* of the library on the same instances,
 * on the 8-puzzle whose moves cost the tile moved and on small SAS+ tasks.
 * Every path is checked on the explicit domain: it must start at the
 * initial state, end at a goal and be made of its moves, whose costs add up
 * to the one reported. The counters of the searches, which do not depend
 * on the workers since the diagrams are canonical, are checked to be the
 * same with 1 to 4 workers of Lace and of each executor.
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

#include <cmath>

#include <fstream>

#include <map>

#include <queue>

#include <set>

#include <sstream>

#include <string>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/symbolic.h"

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

/*--------------------------------------------------------------------------*/
/*---------------------------- EXPLICIT DOMAINS ----------------------------*/
/*--------------------------------------------------------------------------*/
// the 8-puzzle whose move of the tile t costs t + shift (heavy tiles, with
// moves of cost 0 if shift is -1), symbolic and explicit

struct HeavyEight {

 using cost_type = int;
 using explicit_type = EightPuzzle;
 using state_type = EightPuzzle::state_type;
 using board_type = EightPuzzle::board_type;

 board_type start;
 int shift = 0;

 // the symbolic side
 static std::size_t num_variables( void ) { return( 9 ); }
 static std::size_t range( std::uint32_t ) { return( 9 ); }
 std::vector< std::uint32_t > initial_values( void ) const {
  return( std::vector< std::uint32_t >( start.begin() , start.end() ) );
  }
 template< class F >
 void for_each_goal( F && f ) const {
  std::vector< SymbolicFact > p;
  for( std::uint32_t c = 0 ; c < 9 ; ++c )
   p.push_back( { c , c } );
  f( p );
  }
 template< class F >
 void for_each_operator( F && f ) const {
  for( int b = 0 ; b < 9 ; ++b )
   for( int c : { b - 3 , b - 1 , b + 1 , b + 3 } ) {
    if( ( c < 0 ) || ( c > 8 ) || ( ( c == b - 1 ) && ( b % 3 == 0 ) ) ||
        ( ( c == b + 1 ) && ( b % 3 == 2 ) ) )
     continue;
    for( std::uint32_t t = 1 ; t < 9 ; ++t ) {
     SymbolicOperator< int > o;
     o.cost = int( t ) + shift;
     o.pre = { { std::uint32_t( b ) , 0 } , { std::uint32_t( c ) , t } };
     o.effects = { SymbolicEffect::assign( std::uint32_t( b ) , t ) ,
                   SymbolicEffect::assign( std::uint32_t( c ) , 0 ) };
     f( o );
     }
    }
  }
 static state_type decode( const std::vector< std::uint32_t > & x ) {
  board_type b{};
  for( int c = 0 ; c < 9 ; ++c )
   b[ c ] = std::uint8_t( x[ c ] );
  return( EightPuzzle::pack( b ) );
  }

 // the explicit side, with h = 0
 state_type initial_state( void ) const { return( EightPuzzle::pack( start ) );
  }
 static bool is_goal( state_type s ) {
  return( s == EightPuzzle::pack( EightPuzzle::korf_goal() ) );
  }
 static int heuristic( state_type ) { return( 0 ); }
 static std::uint64_t hash( state_type s ) { return( ADS::splitmix64( s ) ); }
 template< class F >
 void for_each_successor( state_type s , int , F && f ) const {
  const auto b = EightPuzzle::unpack( s );
  int z = 0;
  while( b[ z ] != 0 )
   ++z;
  for( int c : { z - 3 , z - 1 , z + 1 , z + 3 } ) {
   if( ( c < 0 ) || ( c > 8 ) || ( ( c == z - 1 ) && ( z % 3 == 0 ) ) ||
       ( ( c == z + 1 ) && ( z % 3 == 2 ) ) )
    continue;
   auto t = b;
   std::swap( t[ z ] , t[ c ] );
   f( EightPuzzle::pack( t ) , int( b[ c ] ) + shift , 0 );
   }
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a SAS+ task with the members that SymbolicSAS reads

struct Task {
 struct Fact { int var , val; };
 struct Effect {
  int var , pre , post;
  std::vector< Fact > cond = {};
  };
 struct Operator {
  std::vector< Fact > prevail;
  std::vector< Effect > effects;
  int cost;
  };
 std::vector< int > range , init;
 std::vector< Fact > goal;
 std::vector< Operator > ops;
 bool use_costs = true;
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a robot in rooms carrying balls: the variable 0 is the room of the
// robot, 1 .. balls the place of each ball (a room, or rooms for the
// hand), the last one the hand (0 free, 1 full); a move costs move, a
// pick pick and a drop drop, and the goal puts the balls in the last room

Task robot_task( int rooms , int balls , int move , int pick , int drop )
{
 Task t;
 t.range.push_back( rooms );
 t.init.push_back( 0 );
 for( int b = 0 ; b < balls ; ++b ) {
  t.range.push_back( rooms + 1 );
  t.init.push_back( b % rooms );
  t.goal.push_back( { 1 + b , rooms - 1 } );
  }
 const int hand = 1 + balls;
 t.range.push_back( 2 );
 t.init.push_back( 0 );
 for( int r = 0 ; r < rooms ; ++r )
  for( int s = 0 ; s < rooms ; ++s )
   if( s != r )
    t.ops.push_back( { {} , { { 0 , r , s , {} } } , move } );
 for( int b = 0 ; b < balls ; ++b )
  for( int r = 0 ; r < rooms ; ++r ) {
   t.ops.push_back( { { { 0 , r } } ,
                      { { 1 + b , r , rooms , {} } , { hand , 0 , 1 , {} } } ,
                      pick } );
   t.ops.push_back( { { { 0 , r } } ,
                      { { 1 + b , rooms , r , {} } , { hand , 1 , 0 , {} } } ,
                      drop } );
   }
 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the briefcase of Pednault: the variable 0 is the room of the case, 1 + 2 o
// the room of the object o and 2 + 2 o whether it is in the case; a move of
// the case (cost 2) carries the objects in it, by conditional effects, and
// putting an object in or taking it out costs 1; the goal puts the object o
// in the room ( o + 1 ) mod rooms

Task briefcase_task( int rooms , int objects )
{
 Task t;
 t.range.push_back( rooms );
 t.init.push_back( 0 );
 for( int o = 0 ; o < objects ; ++o ) {
  t.range.push_back( rooms );
  t.init.push_back( o % rooms );
  t.range.push_back( 2 );
  t.init.push_back( 0 );
  t.goal.push_back( { 1 + 2 * o , ( o + 1 ) % rooms } );
  }
 for( int r = 0 ; r < rooms ; ++r )
  for( int q = 0 ; q < rooms ; ++q ) {
   if( q == r )
    continue;
   Task::Operator m{ {} , { { 0 , r , q , {} } } , 2 };
   for( int o = 0 ; o < objects ; ++o )
    m.effects.push_back( { 1 + 2 * o , -1 , q , { { 2 + 2 * o , 1 } } } );
   t.ops.push_back( m );
   }
 for( int o = 0 ; o < objects ; ++o )
  for( int r = 0 ; r < rooms ; ++r ) {
   t.ops.push_back( { { { 0 , r } , { 1 + 2 * o , r } } ,
                      { { 2 + 2 * o , 0 , 1 , {} } } , 1 } );
   t.ops.push_back( { {} , { { 2 + 2 * o , 1 , 0 , {} } } , 1 } );
   }
 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// k robots, each with a ball and a hand of its own in rooms of its own, as
// in robot_task(): the robot i has the variables 3 i (its room), 3 i + 1
// (its ball) and 3 i + 2 (its hand), and no operator touches two robots

Task robots_task( int k , int rooms )
{
 Task t;
 for( int i = 0 ; i < k ; ++i ) {
  const int R = 3 * i , B = 3 * i + 1 , H = 3 * i + 2;
  t.range.insert( t.range.end() , { rooms , rooms + 1 , 2 } );
  t.init.insert( t.init.end() , { 0 , ( i + 1 ) % rooms , 0 } );
  t.goal.push_back( { B , rooms - 1 } );
  for( int r = 0 ; r < rooms ; ++r ) {
   if( r + 1 < rooms ) {
    t.ops.push_back( { {} , { { R , r , r + 1 , {} } } , 1 + i } );
    t.ops.push_back( { {} , { { R , r + 1 , r , {} } } , 1 + i } );
    }
   t.ops.push_back( { { { R , r } } ,
                      { { B , r , rooms , {} } , { H , 0 , 1 , {} } } , 1 } );
   t.ops.push_back( { { { R , r } } ,
                      { { B , rooms , r , {} } , { H , 1 , 0 , {} } } , 1 } );
   }
  }
 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a Task as an explicit Domain, with h = 0

struct ExplicitTask {

 using state_type = std::vector< int >;
 using cost_type = int;

 const Task & t;

 state_type initial_state( void ) const { return( t.init ); }
 bool is_goal( const state_type & s ) const {
  for( const auto & g : t.goal )
   if( s[ std::size_t( g.var ) ] != g.val )
    return( false );
  return( true );
  }
 static int heuristic( const state_type & ) { return( 0 ); }
 static std::uint64_t hash( const state_type & s ) {
  std::uint64_t h = 0;
  for( int x : s )
   h = ADS::splitmix64( h ^ std::uint64_t( x ) );
  return( h );
  }
 template< class F >
 void for_each_successor( const state_type & s , int , F && f ) const {
  for( const auto & o : t.ops ) {
   bool ok = true;
   for( const auto & p : o.prevail )
    ok = ok && ( s[ std::size_t( p.var ) ] == p.val );
   for( const auto & e : o.effects )
    ok = ok && ( ( e.pre < 0 ) || ( s[ std::size_t( e.var ) ] == e.pre ) );
   if( ! ok )
    continue;
   auto u = s;
   for( const auto & e : o.effects ) {
    bool holds = true;
    for( const auto & c : e.cond )
     holds = holds && ( s[ std::size_t( c.var ) ] == c.val );
    if( holds )
     u[ std::size_t( e.var ) ] = e.post;
    }
   f( u , t.use_costs ? o.cost : 1 , 0 );
   }
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a SymbolicSlidingTile whose heuristic is twice the Manhattan distance,
// inconsistent

template< int W , int H >
struct DoubledTile : SymbolicSlidingTile< W , H > {
 using SymbolicSlidingTile< W , H >::SymbolicSlidingTile;
 template< class F >
 void for_each_heuristic_term( F && f ) const {
  SymbolicSlidingTile< W , H >::for_each_heuristic_term(
   [ & ]( SymbolicTerm< int > t ) {
    for( auto & x : t.table )
     x *= 2;
    f( t );
    } );
  }
 };

/*--------------------------------------------------------------------------*/
/*-------------------------------- HELPERS ---------------------------------*/
/*--------------------------------------------------------------------------*/
// the names of the directions

const char * dir_name( SymbolicDirection d )
{
 return( d == SymbolicDirection::forward ? "forward" :
         d == SymbolicDirection::backward ? "backward" : "bidirectional" );
 }

const SymbolicDirection all_dirs[] = { SymbolicDirection::forward ,
                                       SymbolicDirection::backward ,
                                       SymbolicDirection::bidirectional };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the number of states at each distance from the initial state of the
// explicit domain d, by a breadth-first search (every move counted 1)

template< class E >
std::vector< double > bfs_layers( const E & d )
{
 struct Hash {
  const E * d;
  std::size_t operator()( const typename E::state_type & s ) const {
   return( std::size_t( d->hash( s ) ) );
   }
  };
 using Map = std::unordered_map< typename E::state_type , int , Hash >;
 Map dist( 1024 , Hash{ & d } );
 std::queue< typename E::state_type > q;
 dist.emplace( d.initial_state() , 0 );
 q.push( d.initial_state() );
 std::vector< double > layers;
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  const int ds = dist[ s ];
  if( std::size_t( ds ) >= layers.size() )
   layers.push_back( 0 );
  ++layers[ std::size_t( ds ) ];
  d.for_each_successor( s , 0 , [ & ]( const auto & t , auto , auto ) {
   if( dist.emplace( t , ds + 1 ).second )
    q.push( t );
   } );
  }
 return( layers );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that sol is solved with cost expected and that its path, decoded
// by the symbolic domain sd, is a path of the explicit domain e from its
// initial state to a goal of that cost

template< class S , class E , class Sol >
void check_path( const S & sd , const E & e , const Sol & sol , int expected ,
                 const std::string & what )
{
 bool ok = sol.solved() && ( sol.cost == expected ) &&
           ( sol.bound == expected ) && ! sol.path.empty();
 std::vector< typename E::state_type > p;
 if( ok )
  for( const auto & x : sol.path )
   p.push_back( sd.decode( x ) );
 ok = ok && ( p.front() == e.initial_state() ) && e.is_goal( p.back() );
 long sum = 0;
 for( std::size_t i = 1 ; ok && ( i < p.size() ) ; ++i ) {
  long best = -1;
  e.for_each_successor( p[ i - 1 ] , 0 , [ & ]( const auto & t , auto c ,
                                                auto ) {
   if( ( t == p[ i ] ) && ( ( best < 0 ) || ( long( c ) < best ) ) )
    best = long( c );
   } );
  ok = best >= 0;
  sum += best;
  }
 ok = ok && ( sum == expected );
 std::ostringstream os;
 os << what << " (cost " << sol.cost << ", expected " << expected
    << ", steps " << sol.stats.steps << ", images " << sol.stats.images
    << ", largest set " << sol.stats.max_set_nodes << " nodes)";
 check( ok , os.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the instances of Korf (1985) with their optimal costs

std::vector< std::pair< FifteenPuzzle::board_type , int > > korf(
 const std::string & dir )
{
 std::ifstream in( dir + "/15-puzzle/korf100.txt" );
 std::vector< std::pair< FifteenPuzzle::board_type , int > > v;
 for( std::string line ; std::getline( in , line ) ; ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  std::istringstream is( line );
  FifteenPuzzle::board_type b{};
  for( auto & t : b ) {
   int x;
   is >> x;
   t = std::uint8_t( x );
   }
  int c;
  is >> c;
  v.emplace_back( b , c );
  }
 if( v.size() != 100 )
  throw( std::runtime_error( "cannot read korf100.txt from " + dir ) );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// 16 states of the 8-puzzle by walks of the blank from the goal, and the
// first of the farthest ones (at 31 moves), by an explicit breadth-first
// search from the goal

std::vector< EightPuzzle::board_type > eight_boards( void )
{
 std::vector< EightPuzzle::board_type > v;
 for( std::uint64_t seed = 1 ; seed <= 16 ; ++seed ) {
  auto b = EightPuzzle::korf_goal();
  std::uint64_t x = seed;
  int z = 0;
  for( int k = 0 ; k < int( 4 * seed ) ; ++k ) {
   x = ADS::splitmix64( x );
   int c = z;
   switch( x % 4 ) {
    case( 0 ): c = z >= 3 ? z - 3 : z; break;
    case( 1 ): c = z % 3 ? z - 1 : z; break;
    case( 2 ): c = z % 3 < 2 ? z + 1 : z; break;
    default: c = z < 6 ? z + 3 : z;
    }
   std::swap( b[ z ] , b[ c ] );
   z = c;
   }
  v.push_back( b );
  }
 const EightPuzzle g( EightPuzzle::korf_goal() );
 std::unordered_map< std::uint64_t , int > dist;
 std::queue< std::uint64_t > q;
 dist[ g.initial_state() ] = 0;
 q.push( g.initial_state() );
 std::uint64_t far = g.initial_state();
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  const int ds = dist[ s ];
  if( ( ds > dist[ far ] ) || ( ( ds == dist[ far ] ) && ( s < far ) ) )
   far = s;
  g.for_each_successor( s , 0 , [ & ]( std::uint64_t t , int , int ) {
   if( dist.emplace( t , ds + 1 ).second )
    q.push( t );
   } );
  }
 v.push_back( EightPuzzle::unpack( far ) );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the counters that must not depend on the workers

bool same_counts( const SymbolicStats & a , const SymbolicStats & b )
{
 return( ( a.steps == b.steps ) && ( a.images == b.images ) &&
         ( a.expanded == b.expanded ) &&
         ( a.expanded_nodes == b.expanded_nodes ) &&
         ( a.max_set_nodes == b.max_set_nodes ) &&
         ( a.max_image_nodes == b.max_image_nodes ) &&
         ( a.path_images == b.path_images ) &&
         ( a.relation_nodes == b.relation_nodes ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- CASES ---------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace )

/*--------------------------------------------------------------------------*/
// the symbolic breadth-first search: the whole state space of the 8-puzzle
// and of the pancakes by layers, and the number of moves forward, backward
// and bidirectionally, against the explicit breadth-first search

void test_sym_bfs( const std::string & )
{
 BddRuntime rt( 2 );

 const auto goal = EightPuzzle::korf_goal();
 const SymbolicSlidingTile< 3 , 3 > g8( goal );
 const auto ref8 = bfs_layers( EightPuzzle( goal ) );
 for( bool backward : { false , true } ) {
  const auto r = symbolic_reachable( g8 , backward );
  check( ( r.states == 181440 ) && ( r.layers == ref8 ) &&
         ( r.layers.size() == 32 ) , std::string( "8-puzzle, the 181440 "
         "states in 32 layers as by the explicit search, " ) +
         ( backward ? "backward" : "forward" ) + " (largest layer " +
         std::to_string( r.max_set_nodes ) + " nodes, all of them " +
         std::to_string( r.reached_nodes ) + ")" );
  }

 PancakePuzzle< 7 >::state_type id7{};
 for( int i = 0 ; i < 7 ; ++i )
  id7[ i ] = std::uint8_t( i );
 const auto rp = symbolic_reachable( SymbolicPancake< 7 >( id7 ) );
 check( ( rp.states == 5040 ) &&
        ( rp.layers == bfs_layers( PancakePuzzle< 7 >( id7 ) ) ) ,
        "7 pancakes, the 5040 stacks by layers as by the explicit search" );

 // a board of the other parity: the half of the space it reaches, and no
 // path to the goal in any direction
 auto odd = goal;
 std::swap( odd[ 1 ] , odd[ 2 ] );
 const SymbolicSlidingTile< 3 , 3 > u8( odd );
 check( symbolic_reachable( u8 ).states == 181440 ,
        "8-puzzle of odd parity, 181440 states reachable" );
 for( auto dir : all_dirs ) {
  SymbolicOptions o;
  o.direction = dir;
  check( symbolic_bfs( u8 , o ).status == Status::unsolvable ,
         std::string( "8-puzzle of odd parity, unsolvable, " ) +
         dir_name( dir ) );
  }

 for( auto dir : all_dirs ) {
  SymbolicOptions o;
  o.direction = dir;
  for( const auto & b : eight_boards() ) {
   const EightPuzzle e( b );
   const int c = astar( e ).cost;
   check_path( SymbolicSlidingTile< 3 , 3 >( b ) , e ,
               symbolic_bfs( SymbolicSlidingTile< 3 , 3 >( b ) , o ) , c ,
               std::string( "8-puzzle BFS, " ) + dir_name( dir ) );
   }
  for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
   const auto s = PancakePuzzle< 8 >::random( seed );
   const PancakePuzzle< 8 > e( s );
   check_path( SymbolicPancake< 8 >( s ) , e ,
               symbolic_bfs( SymbolicPancake< 8 >( s ) , o ) ,
               astar( e ).cost , std::string( "8 pancakes BFS, " ) +
               dir_name( dir ) );
   }
  // the start that is the goal, and one move away
  const auto r0 = symbolic_bfs( g8 , o );
  check( r0.solved() && ( r0.cost == 0 ) && ( r0.path.size() == 1 ) ,
         std::string( "8-puzzle, the start is the goal, " ) +
         dir_name( dir ) );
  auto one = goal;
  std::swap( one[ 0 ] , one[ 1 ] );
  check_path( SymbolicSlidingTile< 3 , 3 >( one ) , EightPuzzle( one ) ,
              symbolic_bfs( SymbolicSlidingTile< 3 , 3 >( one ) , o ) , 1 ,
              std::string( "8-puzzle, one move, " ) + dir_name( dir ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// symbolic Dijkstra's algorithm on the 8-puzzle whose moves cost the tile
// moved (one partition per cost) and on the same with moves of cost 0,
// against A* with h = 0, every direction, merged and unmerged partitions

void test_sym_dijkstra( const std::string & )
{
 BddRuntime rt( 2 );
 for( int shift : { 0 , -1 } )
  for( auto dir : all_dirs )
   for( std::size_t merge : { std::size_t( 100000 ) , std::size_t( 0 ) } ) {
    SymbolicOptions o;
    o.direction = dir;
    o.max_relation_nodes = merge;
    const auto boards = eight_boards();
    for( std::size_t i = 0 ; i < boards.size() ; i += 3 ) {
     const HeavyEight h{ boards[ i ] , shift };
     const int c = astar( h ).cost;
     std::ostringstream name;
     name << "8-puzzle of heavy tiles" << ( shift ? " and moves of cost 0"
          : "" ) << ", Dijkstra " << dir_name( dir ) << ( merge ? ", merged"
          : ", one operator per partition" ) << ", board " << i;
     check_path( h , h , symbolic_search( h , o ) , c , name.str() );
     }
    }
 const HeavyEight h{ eight_boards()[ 0 ] , 0 };
 SymbolicSearch< HeavyEight > merged( h ) , apart( h , [] {
  SymbolicOptions o;
  o.max_relation_nodes = 0;
  return( o );
  }() );
 check( ( merged.model().partitions().size() == 8 ) &&
        ( apart.model().partitions().size() == 192 ) &&
        ( apart.model().num_operators() == 192 ) ,
        "8-puzzle of heavy tiles, 192 operators in 8 partitions by cost" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// BDDA* with the Manhattan distance on the easiest instances of Korf and
// with the gap heuristic on the pancakes, every direction, against A*; the
// heuristics checked consistent, and an inconsistent one found so

void test_sym_astar( const std::string & data )
{
 BddRuntime rt( 2 , std::size_t( 4 ) << 30 );

 const SymbolicSlidingTile< 3 , 3 > t8( eight_boards().back() );
 const SymbolicModel< SymbolicSlidingTile< 3 , 3 > > m8( t8 );
 check( m8.consistent( false ) && m8.consistent( true ) ,
        "8-puzzle, the Manhattan distance consistent both ways" );
 const DoubledTile< 3 , 3 > x8( eight_boards().back() );
 const SymbolicModel< DoubledTile< 3 , 3 > > mx( x8 );
 check( ! mx.consistent( false ) ,
        "8-puzzle, twice the Manhattan distance found inconsistent" );
 const SymbolicPancake< 8 > p8( PancakePuzzle< 8 >::random( 7 ) );
 const SymbolicModel< SymbolicPancake< 8 > > mp( p8 );
 check( mp.consistent( false ) && mp.consistent( true ) ,
        "8 pancakes, the gap heuristic consistent both ways" );

 const auto inst = korf( data );
 for( auto dir : all_dirs )
  for( int id : { 55 , 16 , 42 , 79 , 71 , 12 } ) {
   const auto & [ b , c ] = inst[ std::size_t( id - 1 ) ];
   const FifteenPuzzle e( b );
   const int a = astar( e ).cost;
   SymbolicOptions o;
   o.direction = dir;
   check_path( SymbolicSlidingTile< 4 , 4 >( b ) , e ,
               symbolic_search( SymbolicSlidingTile< 4 , 4 >( b ) , o ) , c ,
               std::string( "15-puzzle BDDA*, " ) + dir_name( dir ) +
               ", Korf #" + std::to_string( id ) + ", A* " +
               std::to_string( a ) );
   }

 for( auto dir : all_dirs )
  for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
   const auto s = PancakePuzzle< 9 >::random( seed );
   const PancakePuzzle< 9 > e( s );
   SymbolicOptions o;
   o.direction = dir;
   check_path( SymbolicPancake< 9 >( s ) , e ,
               symbolic_search( SymbolicPancake< 9 >( s ) , o ) ,
               astar( e ).cost , std::string( "9 pancakes BDDA*, " ) +
               dir_name( dir ) + ", seed " + std::to_string( seed ) );
   }

 // the 8-puzzle without and with the heuristic: the same costs, fewer
 // states expanded with it
 for( auto dir : all_dirs ) {
  const auto b = eight_boards().back();
  SymbolicOptions o;
  o.direction = dir;
  const auto withh = symbolic_search( SymbolicSlidingTile< 3 , 3 >( b ) , o );
  o.heuristic = false;
  const auto blind = symbolic_search( SymbolicSlidingTile< 3 , 3 >( b ) , o );
  check( withh.solved() && blind.solved() && ( withh.cost == 31 ) &&
         ( blind.cost == 31 ) &&
         ( withh.stats.expanded < blind.stats.expanded ) ,
         std::string( "8-puzzle, BDDA* expands fewer states than "
                      "Dijkstra, " ) + dir_name( dir ) + " (" +
         std::to_string( withh.stats.expanded ) + " against " +
         std::to_string( blind.stats.expanded ) + ")" );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// SAS+ tasks through SymbolicSAS: a robot carrying balls between rooms,
// with costs, with moves of cost 0, without costs and unsolvable, against
// A* with h = 0 on the explicit task

void test_sym_sas( const std::string & )
{
 BddRuntime rt( 2 );
 struct Case { int rooms , balls , move , pick , drop; bool costs; };
 const Case cases[] = { { 3 , 3 , 2 , 1 , 0 , true } ,
                        { 4 , 3 , 1 , 1 , 1 , false } ,
                        { 3 , 4 , 3 , 0 , 2 , true } };
 for( const auto & k : cases ) {
  Task t = robot_task( k.rooms , k.balls , k.move , k.pick , k.drop );
  t.use_costs = k.costs;
  const ExplicitTask e{ t };
  const int c = astar( e ).cost;
  const SymbolicSAS< Task > sd( t );
  for( auto dir : all_dirs ) {
   SymbolicOptions o;
   o.direction = dir;
   const auto sol = symbolic_search( sd , o );
   bool ok = sol.solved() && ( sol.cost == c );
   // the path, a sequence of states of the task
   long sum = 0;
   for( std::size_t i = 1 ; ok && ( i < sol.path.size() ) ; ++i ) {
    long best = -1;
    std::vector< int > s( sol.path[ i - 1 ].begin() ,
                          sol.path[ i - 1 ].end() ) ,
                       u( sol.path[ i ].begin() , sol.path[ i ].end() );
    e.for_each_successor( s , 0 , [ & ]( const std::vector< int > & v ,
                                         int w , int ) {
     if( ( v == u ) && ( ( best < 0 ) || ( w < best ) ) )
      best = w;
     } );
    ok = best >= 0;
    sum += best;
    }
   ok = ok && ! sol.path.empty() && ( sum == c ) &&
        ( std::vector< int >( sol.path.front().begin() ,
                              sol.path.front().end() ) == t.init ) &&
        e.is_goal( std::vector< int >( sol.path.back().begin() ,
                                       sol.path.back().end() ) );
   check( ok , "robot with " + std::to_string( k.balls ) + " balls in " +
          std::to_string( k.rooms ) + " rooms, " + dir_name( dir ) +
          " (cost " + std::to_string( sol.cost ) + ", expected " +
          std::to_string( c ) + ")" );
   }
  }
 // a goal no operator reaches: a ball in the hand that cannot be filled
 Task t = robot_task( 3 , 2 , 1 , 1 , 1 );
 t.ops.erase( std::remove_if( t.ops.begin() , t.ops.end() ,
                              []( const Task::Operator & o ) {
                               return( o.effects.size() == 2 &&
                                       o.effects[ 1 ].post == 1 );
                               } ) , t.ops.end() );
 t.goal = { { 1 , 3 } };
 for( auto dir : all_dirs ) {
  SymbolicOptions o;
  o.direction = dir;
  check( symbolic_search( SymbolicSAS< Task >( t ) , o ).status ==
         Status::unsolvable , std::string( "robot that cannot pick, "
                                           "unsolvable, " ) +
         dir_name( dir ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the pattern databases and the perimeter: BDDA* with them alone (without
// the terms of the domain) on the 8-puzzle, on its heavy tiles with moves
// of cost 0, on the pancakes and on the robot, every direction, against A*;
// the heuristics checked consistent and admissible by consistent()

void test_sym_pdb( const std::string & )
{
 BddRuntime rt( 2 );
 struct Setting {
  std::vector< std::vector< std::uint32_t > > patterns;
  std::size_t perimeter;
  };
 const std::vector< Setting > settings = {
  { { { 0 , 1 , 2 , 3 , 4 } , { 4 , 5 , 6 , 7 , 8 } } , 0 } ,
  { {} , 6 } ,
  { { { 0 , 1 , 2 , 3 , 4 } , { 4 , 5 , 6 , 7 , 8 } } , 6 } };
 const auto name = []( const Setting & z ) {
  return( std::to_string( z.patterns.size() ) + " patterns, perimeter " +
          std::to_string( z.perimeter ) );
  };
 const auto options = []( const Setting & z , SymbolicDirection dir ) {
  SymbolicOptions o;
  o.direction = dir;
  o.terms = false;
  o.patterns = z.patterns;
  o.perimeter = z.perimeter;
  return( o );
  };
 const auto certified = []( const auto & search ) {
  return( search.model().consistent( false ) &&
          search.model().consistent( true ) &&
          search.model().has_heuristic( false ) &&
          search.model().has_heuristic( true ) );
  };

 const auto far = eight_boards().back();
 const SymbolicSlidingTile< 3 , 3 > t8( far );
 SymbolicOptions blind;
 blind.heuristic = false;
 const auto ref8 = symbolic_search( t8 , blind );
 for( const auto & z : settings ) {
  const SymbolicSearch< SymbolicSlidingTile< 3 , 3 > > c8(
   t8 , options( z , SymbolicDirection::bidirectional ) );
  check( certified( c8 ) , "8-puzzle, " + name( z ) +
         ", consistent and admissible both ways" );
  for( auto dir : all_dirs ) {
   const auto sol = symbolic_search( t8 , options( z , dir ) );
   check_path( t8 , EightPuzzle( far ) , sol , 31 , "8-puzzle, " +
               name( z ) + ", " + dir_name( dir ) + ", " +
               std::to_string( sol.stats.expanded ) + " states expanded "
               "against " + std::to_string( ref8.stats.expanded ) +
               " of Dijkstra" );
   check( sol.stats.expanded <= ref8.stats.expanded ,
          "8-puzzle, " + name( z ) + ", no more states than Dijkstra, " +
          dir_name( dir ) );
   }
  }

 const auto boards = eight_boards();
 for( int shift : { 0 , -1 } )
  for( const auto & z : settings )
   for( std::size_t i = 0 ; i < boards.size() ; i += 4 ) {
    const HeavyEight h{ boards[ i ] , shift };
    const int c = astar( h ).cost;
    const SymbolicSearch< HeavyEight > ch(
     h , options( z , SymbolicDirection::bidirectional ) );
    check( certified( ch ) , "heavy tiles" + std::string( shift ?
           " with moves of cost 0" : "" ) + ", " + name( z ) +
           ", consistent and admissible, board " + std::to_string( i ) );
    for( auto dir : all_dirs )
     check_path( h , h , symbolic_search( h , options( z , dir ) ) , c ,
                 "heavy tiles" + std::string( shift ? " with moves of "
                 "cost 0" : "" ) + ", " + name( z ) + ", " +
                 dir_name( dir ) + ", board " + std::to_string( i ) );
    }

 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto st = PancakePuzzle< 8 >::random( seed );
  const SymbolicPancake< 8 > p8( st );
  const PancakePuzzle< 8 > e( st );
  const int c = astar( e ).cost;
  const Setting z{ { { 0 , 1 , 2 , 3 } , { 4 , 5 , 6 , 7 } } , 3 };
  const SymbolicSearch< SymbolicPancake< 8 > > cp(
   p8 , options( z , SymbolicDirection::bidirectional ) );
  check( certified( cp ) , "8 pancakes, " + name( z ) + ", consistent and "
         "admissible, seed " + std::to_string( seed ) );
  for( auto dir : all_dirs )
   check_path( p8 , e , symbolic_search( p8 , options( z , dir ) ) , c ,
               "8 pancakes, " + name( z ) + ", " + dir_name( dir ) +
               ", seed " + std::to_string( seed ) );
  }

 Task t = robot_task( 3 , 3 , 2 , 1 , 0 );
 const ExplicitTask e{ t };
 const int c = astar( e ).cost;
 const SymbolicSAS< Task > sd( t );
 for( const Setting & z : { Setting{ { { 0 , 1 , 4 } , { 0 , 2 , 3 } } , 0 } ,
                            Setting{ { { 0 , 1 , 4 } } , 4 } } ) {
  const SymbolicSearch< SymbolicSAS< Task > > cs(
   sd , options( z , SymbolicDirection::bidirectional ) );
  check( certified( cs ) , "robot, " + name( z ) + ", consistent and "
         "admissible" );
  for( auto dir : all_dirs ) {
   const auto sol = symbolic_search( sd , options( z , dir ) );
   check( sol.solved() && ( sol.cost == c ) && ! sol.path.empty() ,
          "robot, " + name( z ) + ", " + dir_name( dir ) + " (cost " +
          std::to_string( sol.cost ) + ", expected " + std::to_string( c ) +
          ", " + std::to_string( sol.stats.abstraction_images ) +
          " products of the abstractions)" );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the variants of the model and of the open list: conditional effects (the
// briefcase), additive pattern databases (independent robots, where their
// sum is exact), the order of the variables (given and searched) and List
// BDDA*, every direction, against A* and against the plain search

void test_sym_variants( const std::string & )
{
 BddRuntime rt( 2 );

 // a SAS+ task searched in every direction against A*
 const auto against = [ & ]( const Task & t , SymbolicOptions o ,
                             const std::string & what ) {
  const ExplicitTask e{ t };
  const int c = astar( e ).cost;
  const SymbolicSAS< Task > sd( t );
  for( auto dir : all_dirs ) {
   o.direction = dir;
   const auto sol = symbolic_search( sd , o );
   bool ok = sol.solved() && ( sol.cost == c ) && ! sol.path.empty();
   long sum = 0;
   for( std::size_t i = 1 ; ok && ( i < sol.path.size() ) ; ++i ) {
    long best = -1;
    const std::vector< int > x( sol.path[ i - 1 ].begin() ,
                                sol.path[ i - 1 ].end() ) ,
                             y( sol.path[ i ].begin() , sol.path[ i ].end() );
    e.for_each_successor( x , 0 , [ & ]( const std::vector< int > & v ,
                                         int w , int ) {
     if( ( v == y ) && ( ( best < 0 ) || ( w < best ) ) )
      best = w;
     } );
    ok = best >= 0;
    sum += best;
    }
   check( ok && ( sum == c ) , what + ", " + dir_name( dir ) + " (cost " +
          std::to_string( sol.cost ) + ", expected " + std::to_string( c ) +
          ")" );
   }
  };

 for( int objects : { 2 , 3 } ) {
  const Task t = briefcase_task( 3 , objects );
  against( t , SymbolicOptions() , "briefcase with " +
           std::to_string( objects ) + " objects, conditional effects" );
  SymbolicOptions o;
  o.max_relation_nodes = 0;
  o.patterns = { { 0 , 1 , 2 } , { 0 , 3 , 4 } };
  against( t , o , "briefcase with " + std::to_string( objects ) +
           " objects, projected conditional effects" );
  const SymbolicSAS< Task > sd( t );
  o.direction = SymbolicDirection::bidirectional;
  const SymbolicSearch< SymbolicSAS< Task > > s( sd , o );
  check( s.model().consistent( false ) && s.model().consistent( true ) ,
         "briefcase, the pattern databases consistent and admissible" );
  }

 {
  const Task t = robots_task( 3 , 4 );
  SymbolicOptions o;
  o.patterns = { { 0 , 1 , 2 } , { 3 , 4 , 5 } , { 6 , 7 , 8 } };
  o.additive = true;
  against( t , o , "3 robots, additive pattern databases" );
  const SymbolicSAS< Task > sd( t );
  o.direction = SymbolicDirection::bidirectional;
  const SymbolicSearch< SymbolicSAS< Task > > s( sd , o );
  const ExplicitTask e{ t };
  const int c = astar( e ).cost;
  const auto & H = s.model().levels( false );
  const auto at = std::find_if( H.begin() , H.end() , [ & ]( const auto & x ) {
   return( ! x.second.disjoint( s.model().initial() ) );
   } );
  check( s.model().consistent( false ) && s.model().consistent( true ) &&
         ( at != H.end() ) && ( at->first == c ) ,
         "3 robots, the sum of the pattern databases is exact at the start" );
  bool threw = false;
  try {
   SymbolicOptions bad;
   bad.patterns = { { 1 } , { 2 } };
   bad.additive = true;
   const SymbolicSearch< SymbolicSAS< Task > > x( sd , bad );
   }
  catch( const std::invalid_argument & ) {
   threw = true;
   }
  check( threw , "additive patterns written by one operator are refused" );
  }

 // the order of the variables: reversed, and searched
 {
  const auto b = eight_boards().back();
  const SymbolicSlidingTile< 3 , 3 > t8( b );
  const auto plain = symbolic_search( t8 );
  SymbolicOptions o;
  o.layout.order = { 8 , 7 , 6 , 5 , 4 , 3 , 2 , 1 , 0 };
  const auto rev = symbolic_search( t8 , o );
  check_path( t8 , EightPuzzle( b ) , rev , plain.cost , "8-puzzle, the "
              "variables reversed (" + std::to_string( rev.stats.
              relation_nodes ) + " nodes of the relation against " +
              std::to_string( plain.stats.relation_nodes ) + ")" );
  const Task t = robot_task( 3 , 3 , 2 , 1 , 0 );
  SymbolicOptions q;
  q.layout.optimize = true;
  against( t , q , "robot, the order of the variables searched" );
  const auto p9 = PancakePuzzle< 9 >::random( 2 );
  const SymbolicPancake< 9 > sp( p9 );
  const auto sol = symbolic_search( sp , q );
  check_path( sp , PancakePuzzle< 9 >( p9 ) , sol ,
              astar( PancakePuzzle< 9 >( p9 ) ).cost ,
              "9 pancakes, the order of the variables searched" );
  bool threw = false;
  try {
   SymbolicOptions bad;
   bad.layout.order = { 0 , 1 , 1 , 3 , 4 , 5 , 6 , 7 , 8 };
   const SymbolicSearch< SymbolicSlidingTile< 3 , 3 > > x( t8 , bad );
   }
  catch( const std::invalid_argument & ) {
   threw = true;
   }
  check( threw , "an order that is not a permutation is refused" );
  }

 // List BDDA*: the same sets expanded, fewer conjunctions with H_k
 for( auto dir : all_dirs ) {
  const auto b = eight_boards().back();
  const SymbolicSlidingTile< 3 , 3 > t8( b );
  SymbolicOptions o;
  o.direction = dir;
  const auto matrix = symbolic_search( t8 , o );
  o.list = true;
  const auto list = symbolic_search( t8 , o );
  check_path( t8 , EightPuzzle( b ) , list , matrix.cost ,
              std::string( "8-puzzle, List BDDA*, " ) + dir_name( dir ) +
              " (" + std::to_string( list.stats.splits ) + " conjunctions "
              "with H_k against " + std::to_string( matrix.stats.splits ) +
              ")" );
  check( list.stats.expanded == matrix.stats.expanded ,
         std::string( "8-puzzle, List BDDA* expands the states of BDDA*, " ) +
         dir_name( dir ) );
  const auto p9 = PancakePuzzle< 9 >::random( 4 );
  const SymbolicPancake< 9 > sp( p9 );
  check_path( sp , PancakePuzzle< 9 >( p9 ) , symbolic_search( sp , o ) ,
              astar( PancakePuzzle< 9 >( p9 ) ).cost ,
              std::string( "9 pancakes, List BDDA*, " ) + dir_name( dir ) );
  const HeavyEight h{ b , -1 };
  check_path( h , h , symbolic_search( h , o ) , astar( h ).cost ,
              std::string( "heavy tiles with moves of cost 0, List BDDA*, " )
              + dir_name( dir ) );
  }
 }

namespace {

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the costs of the k cheapest paths of the explicit domain e from its
// initial state to a goal, states repeated, by a uniform-cost search on the
// paths themselves

template< class E >
std::vector< long > k_cheapest( const E & e , std::size_t k )
{
 using S = typename E::state_type;
 using Item = std::pair< long , std::vector< S > >;
 const auto later = []( const Item & a , const Item & b ) {
  return( a.first > b.first );
  };
 std::priority_queue< Item , std::vector< Item > , decltype( later ) >
  q( later );
 q.push( { 0 , { e.initial_state() } } );
 std::vector< long > costs;
 while( ! q.empty() && ( costs.size() < k ) ) {
  const Item it = q.top();
  q.pop();
  if( e.is_goal( it.second.back() ) )
   costs.push_back( it.first );
  e.for_each_successor( it.second.back() , 0 , [ & ]( const S & t , auto c ,
                                                      auto ) {
   auto p = it.second;
   p.push_back( t );
   q.push( { it.first + long( c ) , std::move( p ) } );
   } );
  }
 return( costs );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks the k cheapest plans of sd against those of the explicit e: the
// same costs, by nondecreasing cost, each a path of e of its cost, and no
// two equal

template< class SD , class E >
void check_top_k( const SD & sd , const E & e , std::size_t k ,
                  const std::string & what )
{
 SymbolicSearch< SD > s( sd );
 const auto plans = s.top_k( k );
 const auto ref = k_cheapest( e , k );
 bool ok = plans.size() == ref.size();
 std::set< std::vector< std::vector< std::uint32_t > > > seen;
 for( std::size_t i = 0 ; ok && ( i < plans.size() ) ; ++i ) {
  ok = ( plans[ i ].cost == ref[ i ] ) &&
       seen.insert( plans[ i ].path ).second;
  long sum = 0;
  for( std::size_t j = 1 ; ok && ( j < plans[ i ].path.size() ) ; ++j ) {
   long best = -1;
   const auto x = sd.decode( plans[ i ].path[ j - 1 ] );
   const auto y = sd.decode( plans[ i ].path[ j ] );
   e.for_each_successor( x , 0 , [ & ]( const auto & v , auto c , auto ) {
    if( ( v == y ) && ( ( best < 0 ) || ( long( c ) < best ) ) )
     best = long( c );
    } );
   ok = best >= 0;
   sum += best;
   }
  ok = ok && ( sum == plans[ i ].cost ) &&
       ( sd.decode( plans[ i ].path.front() ) == e.initial_state() ) &&
       e.is_goal( sd.decode( plans[ i ].path.back() ) );
  }
 check( ok , what + ", the " + std::to_string( k ) + " cheapest plans (" +
        std::to_string( plans.size() ) + " found, the last of cost " +
        ( plans.empty() ? std::string( "-" ) :
          std::to_string( plans.back().cost ) ) + ")" );
 }

}  // end( namespace )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the k cheapest plans, states repeated, against the explicit enumeration
// of the paths, on the 8-puzzle and on its heavy tiles; a domain with moves
// of cost 0 is refused

void test_sym_topk( const std::string & )
{
 BddRuntime rt( 2 );
 const auto boards = eight_boards();
 {
  check_top_k( SymbolicSlidingTile< 3 , 3 >( boards[ 0 ] ) ,
               EightPuzzle( boards[ 0 ] ) , 15 , "8-puzzle, board 0" );
  const HeavyEight h{ boards[ 0 ] , 0 };
  check_top_k( h , h , 10 , "heavy tiles, board 0" );
  }
 auto g = EightPuzzle::korf_goal();
 check_top_k( SymbolicSlidingTile< 3 , 3 >( g ) , EightPuzzle( g ) , 10 ,
              "8-puzzle, the start is the goal" );
 bool threw = false;
 try {
  const HeavyEight h{ boards[ 0 ] , -1 };
  SymbolicSearch< HeavyEight > s( h );
  ( void ) s.top_k( 3 );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "the cheapest plans refused with moves of cost 0" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel images: the partitions as tasks of Lace and on 1 to 4
// workers of the executor, with 1 to 4 workers of Lace, give the costs and
// every counter of the sequential images

template< class Executor >
static void test_sym_parallel( const std::string & ex )
{
 struct Run {
  SymbolicStats stats;
  int cost;
  };
 const auto runs = [ & ]( SymbolicImages im , unsigned threads ) {
  std::vector< Run > v;
  for( auto dir : all_dirs ) {
   SymbolicOptions o;
   o.direction = dir;
   o.images = im;
   o.threads = threads;
   o.max_relation_nodes = 0;
   const auto b = eight_boards().back();
   const auto s1 = symbolic_search< SymbolicSlidingTile< 3 , 3 > ,
                                    Executor >(
    SymbolicSlidingTile< 3 , 3 >( b ) , o );
   v.push_back( { s1.stats , s1.cost } );
   const auto p = PancakePuzzle< 9 >::random( 3 );
   const auto s2 = symbolic_search< SymbolicPancake< 9 > , Executor >(
    SymbolicPancake< 9 >( p ) , o );
   v.push_back( { s2.stats , s2.cost } );
   o.max_relation_nodes = 100000;
   const HeavyEight h{ b , -1 };
   const auto s3 = symbolic_search< HeavyEight , Executor >( h , o );
   v.push_back( { s3.stats , s3.cost } );
   const auto r = symbolic_reachable< SymbolicSlidingTile< 3 , 3 > ,
                                      Executor >(
    SymbolicSlidingTile< 3 , 3 >( b ) , false , o );
   SymbolicStats st;
   st.images = r.images;
   st.max_set_nodes = r.max_set_nodes;
   st.expanded = r.states;
   v.push_back( { st , int( r.layers.size() ) } );
   }
  return( v );
  };
 std::vector< Run > ref;
 {
  BddRuntime rt( 1 );
  ref = runs( SymbolicImages::sequential , 1 );
  }
 const auto same = [ & ]( const std::vector< Run > & v ) {
  bool ok = v.size() == ref.size();
  for( std::size_t i = 0 ; ok && ( i < v.size() ) ; ++i )
   ok = ( v[ i ].cost == ref[ i ].cost ) &&
        same_counts( v[ i ].stats , ref[ i ].stats );
  return( ok );
  };
 for( unsigned lace : { 1u , 2u , 4u } ) {
  BddRuntime rt( lace );
  check( same( runs( SymbolicImages::sequential , 1 ) ) ,
         "sequential images, " + std::to_string( lace ) +
         " workers of Lace, the counters of 1" );
  check( same( runs( SymbolicImages::lace , 1 ) ) ,
         "the partitions as tasks, " + std::to_string( lace ) +
         " workers of Lace" );
  for( unsigned p : { 1u , 2u , 4u } ) {
   const auto v = runs( SymbolicImages::executor , p );
   bool split = true;
   for( const auto & r : v )
    if( ! r.stats.images_by.empty() ) {
     std::uint64_t sum = 0;
     for( auto x : r.stats.images_by )
      sum += x;
     split = split && ( sum == r.stats.images ) &&
             ( r.stats.images_by.size() == p );
     }
   check( same( v ) && split , "the partitions on " + std::to_string( p ) +
          " workers of " + ex + ", " + std::to_string( lace ) +
          " workers of Lace" );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void test_sym_par_thread( const std::string & )
{
 test_sym_parallel< ADS::ThreadPool >( "thread" );
 }

#ifdef PBEFS_HAS_OPENMP
void test_sym_par_openmp( const std::string & )
{
 test_sym_parallel< ADS::OpenMPExecutor >( "openmp" );
 }
#endif

#ifdef PBEFS_HAS_FASTFLOW
void test_sym_par_fastflow( const std::string & )
{
 test_sym_parallel< ADS::FastFlowExecutor >( "fastflow" );
 }
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File symbolic.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
