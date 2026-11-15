/*--------------------------------------------------------------------------*/
/*---------------------------- File manual.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] that check the
 * exercises and the examples of use of the manual: a shifted and a
 * weighted heuristic, the example domains, the options of SafePBNF and the
 * examples between their markers, which the chapters include.
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

#include "test.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- CASES ----------------------------------*/
/*--------------------------------------------------------------------------*/
// a Domain with the heuristic of D replaced by W h + Delta on every state,
// goals included: weighted A* with Delta = 0, a shifted heuristic with
// W = 1

template< class D , int W , int Delta >
struct Reweighted {

 using state_type = state_t< D >;
 using cost_type = cost_t< D >;

 const D & d;

 state_type initial_state( void ) const { return( d.initial_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 cost_type heuristic( const state_type & s ) const {
  return( W * d.heuristic( s ) + Delta );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_type & t , cost_type c ,
                               cost_type ht ) {
                         f( t , c , W * ht + Delta );
                         } );
  }

 };
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* with the Manhattan distance shifted by a constant on every state of
// the 8-puzzle: the same costs and the same expansions, although the
// shifted heuristic is not admissible

void test_shift( const std::string & )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 std::size_t k = 0;
 bool same = true;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 907 ) != 0 )                     // 200 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto a = astar( Reweighted< EightPuzzle , 1 , 0 >{ d } );
  const auto b = astar( Reweighted< EightPuzzle , 1 , 5 >{ d } );
  same = same && a.solved() && b.solved() && ( a.cost == dst ) &&
         ( b.cost == dst ) && ( a.stats.expanded == b.stats.expanded );
  }
 check( same , "A* with h + 5 on every state: the optimal costs and the "
               "same expansions on 200 states of the 8-puzzle" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// weighted A* by a wrapper of the domain: with twice the Manhattan
// distance, a cost between C* and 2 C* on 200 states of the 8-puzzle

void test_weighted_wrapper( const std::string & )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 std::size_t k = 0;
 bool bounded = true;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 907 ) != 0 )
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const Reweighted< EightPuzzle , 2 , 0 > w{ d };
  const auto r = astar( w );
  bounded = bounded && r.solved() && ( r.cost >= dst ) &&
            ( r.cost <= 2 * dst );
  if( r.solved() )
   check_solution( w , r , r.cost , "weighted A*, 8-puzzle" );
  }
 check( bounded , "weighted A* with w = 2: C* <= cost <= 2 C* on 200 "
                  "states of the 8-puzzle" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// what a breadth-first scan of the states reachable from the initial state
// of a Domain finds: the states, the edges, the goals, the edges on which
// h( u ) > c( u , v ) + h( v ), and the successors whose h passed by the
// expansion differs from heuristic()

struct DomainScan {
 std::size_t states = 0;
 std::size_t edges = 0;
 std::size_t goals = 0;
 std::size_t inconsistent = 0;
 std::size_t wrong_h = 0;
 };

template< class D >
static DomainScan scan_domain( const D & d )
{
 using S = state_t< D >;
 using C = cost_t< D >;
 struct Hash {
  const D * d;
  std::size_t operator()( const S & s ) const {
   return( std::size_t( d->hash( s ) ) );
   }
  };
 std::unordered_map< S , char , Hash > seen( 64 , Hash{ & d } );
 std::queue< S > q;
 DomainScan r;
 seen.emplace( d.initial_state() , 0 );
 q.push( d.initial_state() );
 while( ! q.empty() ) {
  const S s = q.front();
  q.pop();
  const C hs = d.heuristic( s );
  if( d.is_goal( s ) )
   ++r.goals;
  d.for_each_successor( s , hs , [ & ]( const S & t , C c , C ht ) {
   ++r.edges;
   if( std::abs( double( ht ) - double( d.heuristic( t ) ) ) > 1e-9 )
    ++r.wrong_h;
   if( double( hs ) > double( c ) + double( ht ) + 1e-9 )
    ++r.inconsistent;
   if( seen.emplace( t , 0 ).second )
    q.push( t );
   } );
  }
 r.states = seen.size();
 return( r );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the domains of the examples: their state spaces against the properties
// of their heuristics, their branching factors and their known answers

void test_domains( void )
{
 // the 8-puzzle: 9! / 2 states, 8 / 3 moves per state, Manhattan consistent
 const EightPuzzle e( EightPuzzle::korf_goal() );
 const auto se = scan_domain( e );
 check( ( se.states == 181440 ) && ( 3 * se.edges == 8 * se.states ) &&
        ! se.inconsistent && ! se.wrong_h ,
        "8-puzzle: 9!/2 states, 8/3 moves per state, consistent h" );

 // the sequence score overestimates one move away from the goal
 const SequenceScore< 3 , 3 > sq( EightPuzzle::korf_goal() );
 int worst = 0;
 e.for_each_successor( e.initial_state() , 0 ,
                       [ & ]( std::uint64_t t , int , int ) {
  worst = std::max( worst , sq.heuristic( t ) );
  } );
 check( worst > 1 , "sequence score above h* = 1 one move from the goal" );

 // a grid of 5 x 5 free cells: 80 moves with 4, 144 with 8
 const GridMap open5( 5 , 5 );
 const auto s4 = scan_domain( Grid< 4 >( open5 , 0 , 4 , 4 , 4 ) );
 const auto s8 = scan_domain( Grid< 8 >( open5 , 0 , 4 , 4 , 4 ) );
 check( ( s4.edges == 80 ) && ( s8.edges == 144 ) && ! s4.inconsistent &&
        ! s8.inconsistent , "5 x 5 grid: 80 and 144 moves, consistent h" );
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 60 , 40 , 0.25 , seed , 5 );
  check( ! scan_domain( Grid< 4 >( m , 0 , 39 , 59 , 39 ) ).inconsistent &&
         ! scan_domain( Grid< 8 >( m , 0 , 39 , 59 , 39 ) ).inconsistent ,
         "grids of costs 1 to 5: Manhattan and octile consistent, seed " +
         std::to_string( seed ) );
  }

 // Romania: both heuristics consistent, towards every city
 int bad = 0;
 for( std::uint32_t c = 0 ; c < Romania::NCities ; ++c )
  bad += scan_domain( Romania( Romania::Arad ,
                               Romania::City( c ) ) ).inconsistent > 0;
 check( bad == 0 , "Romania: consistent heuristic towards every city" );

 // the queens: 2057 partial placements and 92 solutions for 8
 const auto q8 = scan_domain( NQueens( 8 ) );
 check( ( q8.states == 2057 ) && ( q8.goals == 92 ) && ! q8.inconsistent ,
        "8 queens: 2057 states, 92 solutions, consistent h" );
 check( astar( NQueens( 3 ) ).status == Status::unsolvable ,
        "3 queens: unsolvable" );

 // the vacuum world: n 2^n states, consistent h, the answer of the book
 const auto sv = scan_domain( VacuumWorld( 10 , 4 , 0x3FF ) );
 check( ( sv.states == 10 * 1024 ) && ! sv.inconsistent ,
        "vacuum world of 10 cells: 10 * 2^10 states, consistent h" );
 check( astar( VacuumWorld( 2 , 0 , 3 ) ).cost == 3 ,
        "vacuum world of the book: Suck, Right, Suck" );

 // the water jugs
 const WaterJugs j( { 5 , 3 } , 4 );
 check( ( astar( j ).cost == 6 ) && ! scan_domain( j ).inconsistent ,
        "jugs of 5 and 3: 4 units in 6 moves, consistent h" );

 // the missionaries: admissible (h <= the distance to the goal, by a
 // breadth-first search from the goal, the crossings being reversible)
 // but not consistent
 const MissionariesCannibals mc( 3 , 3 , 2 );
 const std::uint64_t goal = std::uint64_t( 1 ) << 32;  // 0 , 0 , right
 std::unordered_map< std::uint64_t , int > dist;
 std::queue< std::uint64_t > q;
 dist[ goal ] = 0;
 q.push( goal );
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  mc.for_each_successor( s , 0 , [ & ]( std::uint64_t t , int , int ) {
   if( dist.emplace( t , dist[ s ] + 1 ).second )
    q.push( t );
   } );
  }
 bool admissible = true;
 for( const auto & [ s , d ] : dist )
  admissible = admissible && ( mc.heuristic( s ) <= d );
 check( admissible && ( dist[ mc.initial_state() ] == 11 ) &&
        ( astar( mc ).cost == 11 ) && scan_domain( mc ).inconsistent ,
        "missionaries 3 3 2: 11 crossings, admissible, not consistent" );
 check( astar( MissionariesCannibals( 4 , 4 , 2 ) ).status ==
        Status::unsolvable , "missionaries 4 4 2: unsolvable" );

 // the traveling salesperson: A* against Held and Karp, consistent h, and
 // ( n - 1 ) 2^( n - 2 ) + 2 states
 for( int n = 5 ; n <= 9 ; ++n ) {
  const TravelingSalesman t( TravelingSalesman::random_cities(
                              n , std::uint64_t( n ) ) );
  const auto st = scan_domain( t );
  check( ( std::abs( astar( t ).cost - t.held_karp() ) < 1e-9 ) &&
         ! st.inconsistent && ( st.states == std::size_t(
                                 ( n - 1 ) << ( n - 2 ) ) + 2 ) ,
         "TSP of " + std::to_string( n ) + " cities: Held and Karp, "
         "consistent h, all the states" );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// SafePBNF with min_expansions 0 (taken as 1), 1, 7 and 100000, on the
// 15-puzzle, on a grid with 8 moves and on the vacuum world

template< class Executor >
void test_pbnf_options( const std::string & data )
{
 const auto inst = korf( data );
 const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
 const FifteenPuzzle d( ki.board );
 const auto m = GridMap::random( 300 , 200 , 0.25 , 2 );
 const Grid< 8 > g8( m , 0 , 199 , 299 , 199 , 100 , 10 );
 const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
 const VacuumWorld v( 12 , 5 , 0xFFFu );
 const auto rv = astar( v );
 for( unsigned e : { 0u , 1u , 7u , 100000u } )
  for( unsigned p : { 1u , 3u , 8u } ) {
   PBNFOptions o;
   o.threads = p;
   o.min_expansions = e;
   const std::string tag = "SafePBNF min_expansions " + std::to_string( e ) +
                           ", p " + std::to_string( p );
   check_solution( d , SafePBNF< FifteenPuzzle , Executor >( d , o ).solve() ,
                   ki.cost , tag + ", Korf #" +
                             std::to_string( easy_korf[ 0 ] ) );
   if( r8.solved() )
    check_solution( g8 , SafePBNF< Grid< 8 > , Executor >( g8 , o ).solve() ,
                    r8.cost , tag + ", grid 8" );
   check_solution( v , SafePBNF< VacuumWorld , Executor >( v , o ).solve() ,
                   rv.cost , tag + ", vacuum world" );
   }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the Domain of the first example of use of the manual: the states 0 ... n
// of a line, the moves + 1 of cost 1 and + 3 of cost 2, and the exact
// distance to n as the heuristic

//! [introduction-usage-domain]
struct Line {

 using state_type = int;
 using cost_type = int;

 int n;                                        // the goal

 int initial_state( void ) const { return( 0 ); }
 bool is_goal( int s ) const { return( s == n ); }
 int heuristic( int s ) const { return( ( n - s ) / 3 * 2 + ( n - s ) % 3 ); }
 std::uint64_t hash( int s ) const {
  return( std::uint64_t( s ) * 0x9E3779B97F4A7C15ULL );
  }

 template< class F >
 void for_each_successor( int s , int , F && f ) const {
  if( s + 1 <= n )
   f( s + 1 , 1 , heuristic( s + 1 ) );
  if( s + 3 <= n )
   f( s + 3 , 2 , heuristic( s + 3 ) );
  }

 };

static_assert( Domain< Line > );
//! [introduction-usage-domain]
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the examples of use of the manual, each between its markers, and the
// checks of their results after them

void test_manual( const std::string & data )
{
 const auto inst = korf( data );
 const auto & ki = inst[ easy_korf[ 0 ] - 1 ];  // the easiest, cost 41

 {
  //! [introduction-usage]
  const auto sol = astar( Line{ 10 } );    // the bucket queue: int costs
  // sol.status, sol.cost, sol.path from 0 to 10, and the counters sol.stats
  //! [introduction-usage]
  const Line line{ 10 };
  check( sol.solved() && ( sol.cost == 7 ) && ( sol.path.front() == 0 ) &&
         ( sol.path.back() == 10 ) && ( sol.stats.expanded >= 1 ) &&
         ( astar( Blind< Line >{ line } ).cost == 7 ) ,
         "manual: the line from 0 to 10 in 3 moves of 3 and one of 1" );
  }

 {
  //! [best-first-usage]
  const EightPuzzle d( EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" ) );
  SearchOptions o;
  o.expected_nodes = 1 << 12;                // CLOSED sized at the start
  const auto s = astar( d , OpenKind::automatic , o );   // the buckets
  const auto t = astar( d , OpenKind::heap , o );        // the heap
  // the path goes from the initial state to a goal, one move per cost unit
  const bool path = s.solved() && ( s.status == Status::solved ) &&
                    ( s.path.front() == d.initial_state() ) &&
                    d.is_goal( s.path.back() ) &&
                    ( s.path.size() == std::size_t( s.cost ) + 1 );
  const auto & st = s.stats;      // expanded, reexpanded, generated,
                                  // duplicates and stored
  o.node_limit = 5;               // stop after storing 5 nodes
  const auto l = astar( d , OpenKind::automatic , o );
  // l.status == Status::node_limit, and l.path is empty
  //! [best-first-usage]
  const auto dist = eight_distances( d );
  check( path && t.solved() && ( t.cost == s.cost ) &&
         ( s.cost == dist.at( d.initial_state() ) ) ,
         "manual: A* on the 8-puzzle, both OPEN lists, cost " +
         std::to_string( s.cost ) );
  check( st.reexpanded == 0 , "manual: no reexpansion, h consistent" );
  check( st.generated >= st.stored - 1 + st.duplicates ,
         "manual: generated >= stored - 1 + duplicates" );
  check( ( l.status == Status::node_limit ) && l.path.empty() ,
         "manual: A* stops at the node limit of 5" );
  }

 {
  //! [data-structures-usage]
  const FifteenPuzzle d( ki.board );
  SearchOptions o;
  o.expected_nodes = 1 << 20;                // CLOSED sized at the start
  const auto b = astar( d , OpenKind::buckets , o );     // the buckets
  const auto h = astar( d , OpenKind::heap , o );        // the heap
  const auto a = astar( d , OpenKind::automatic , o );   // here the buckets
  // the same choice at compile time, by the OPEN list of AStar
  using NodeT = Node< state_t< FifteenPuzzle > , cost_t< FifteenPuzzle > >;
  const auto e = AStar< FifteenPuzzle , HeapOpenList< NodeT > >( d , o )
                 .solve();
  //! [data-structures-usage]
  check_solution( d , b , ki.cost , "manual: A* with the buckets" );
  check_solution( d , h , ki.cost , "manual: A* with the heap" );
  check_solution( d , a , ki.cost , "manual: A* with the automatic OPEN" );
  check_solution( d , e , ki.cost , "manual: AStar with HeapOpenList" );
  check( ( a.stats.expanded == b.stats.expanded ) &&
         ( e.stats.expanded == h.stats.expanded ) ,
         "manual: automatic is buckets, AStar with HeapOpenList is heap" );
  }

 {
  ADS::TempDirectory dir;                       // where the files go
  //! [memory-bounded-usage]
  const FifteenPuzzle d( ki.board );
  const auto i = idastar( d );                  // memory: the path only
  DDDOptions x;
  x.partitions = 16;              // files per bucket (default 64)
  x.buffer_bytes = 1 << 12;       // buffer of each file written (2^16)
  x.directory = dir.path();       // if empty, the temporary of the system
  const auto e = ddd_astar( d , SearchOptions{} , x );
  // the node limit bounds the expansions of IDA*, the states closed by
  // DDD-A*
  SearchOptions lim;
  lim.node_limit = 1000;
  const auto l = idastar( d , lim );
  //! [memory-bounded-usage]
  check_solution( d , i , ki.cost , "manual: IDA*" );
  check( i.stats.stored == std::uint64_t( ki.cost ) + 1 ,
         "manual: IDA* stores the path only" );
  check_solution( d , e , ki.cost , "manual: DDD-A*" );
  check( e.stats.expanded <= e.stats.stored ,
         "manual: DDD-A* expands a state at most once" );
  check( std::filesystem::is_empty( dir.path() ) ,
         "manual: DDD-A* removes its files" );
  check( l.status == Status::node_limit , "manual: IDA* at the node limit" );
  }

 {
  //! [heuristics-usage]
  PDBOptions o;
  o.blank = BlankRule::ignored;   // consistent tables (the default)
  o.reflect = true;               // also the reflection about the diagonal
  o.threads = 2;                  // workers of the construction
  // a directory in o.cache: the tables built once, then read from there
  const Partition part{ { 1 , 2 , 3 , 4 } , { 5 , 6 , 7 , 8 } };
  const AdditivePDB< 3 , 3 > pdb( part , EightPuzzle::korf_goal() , o );
  const auto board = EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" );
  const SlidingTilePDB< 3 , 3 > d( board , pdb );
  const auto s = astar( d );      // s.cost and s.stats.expanded
  // the residues modulo 3 need the tables without the reflection
  o.reflect = false;
  const AdditivePDB< 3 , 3 , Mod3Entries > pdb3( part ,
                                                 EightPuzzle::korf_goal() ,
                                                 o );
  const SlidingTilePDB< 3 , 3 , Mod3Entries > d3( board , pdb3 );
  const auto s3 = astar( d3 );
  //! [heuristics-usage]
  const EightPuzzle m( board );
  const auto sm = astar( m );
  std::uint8_t q[ 9 ];
  SlidingTilePDB< 3 , 3 >::cells_of( d.initial_state() , q );
  check( pdb.reflected() , "manual: the 4-4 tables with the reflection" );
  check( pdb.heuristic( q ) >= m.heuristic( m.initial_state() ) ,
         "manual: the tables at least the Manhattan distance" );
  check( s.solved() && s3.solved() && sm.solved() && ( s.cost == sm.cost ) &&
         ( s3.cost == sm.cost ) ,
         "manual: Manhattan, nibbles and mod 3 give cost " +
         std::to_string( sm.cost ) );
  check( ( s.stats.reexpanded == 0 ) && ( s3.stats.reexpanded == 0 ) ,
         "manual: no reexpansion with the tables, blank ignored" );
  }

 {
  using Executor = ADS::ThreadPool;             // the default of the searches
  //! [shared-memory-usage]
  const FifteenPuzzle d( ki.board );
  ParallelOptions o;
  o.threads = 4;          // the only option of SPA*, MQA* also reads o.seed
  SPAStar< FifteenPuzzle > spa( d , o );
  const auto s = spa.solve();     // on a pool of o.threads threads of its own
  Executor ex( 4 );               // one executor for several searches,
                                  // which must outlive them
  MQAStar< FifteenPuzzle > mq( d , o );
  mq.set_executor( ex );
  const auto m = mq.solve();
  const auto expanded = s.stats.expanded;       // by all the workers
  const auto reexpanded = s.stats.reexpanded;
  const double lb = s.stats.load_balance();     // 1 if perfectly balanced
  //! [shared-memory-usage]
  const auto ref = AStar< FifteenPuzzle >( d ).solve();
  check( ( s.status == Status::solved ) && ( m.status == Status::solved ) &&
         ( s.cost == m.cost ) && ( s.cost == ref.cost ) &&
         ( ref.cost == ki.cost ) , "manual: SPA* and MQA* on 4 workers" );
  check( ( reexpanded <= expanded ) && ( lb >= 1 ) ,
         "manual: the counters of SPA*" );
  }

 {
  //! [hash-distributed-usage]
  const EightPuzzle d( EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" ) );
  ParallelOptions o;
  o.threads = 4;
  o.batch = 64;           // nodes per message (the default)
  o.flush_period = 64;    // expansions between two flushes of an outbox
  // the seed changes the owners of the states, not the cost
  const ZobristDistribution< EightPuzzle > z( d , 7 );
  const auto a = HDAStar< EightPuzzle , decltype( z ) >( d , z , o )
                 .solve();
  const AbstractZobristDistribution< EightPuzzle > az( d );
  const auto b = HDAStar< EightPuzzle , decltype( az ) >( d , az , o )
                 .solve();
  const double co = a.stats.communication_overhead();  // sent / generated
  const double lb = a.stats.load_balance();
  const auto sent = a.stats.sent;           // nodes sent to another worker
  const auto messages = a.stats.messages;   // the batches they went in
  //! [hash-distributed-usage]
  check( ( a.status == Status::solved ) && ( b.status == Status::solved ) &&
         ( a.cost == 26 ) && ( b.cost == 26 ) ,
         "manual: HDA* with both Zobrist distributions" );
  check( ( std::abs( co - double( sent ) / double( a.stats.generated ) ) <
           1e-12 ) && ( messages <= sent ) && ( lb >= 1 ) ,
         "manual: the counters of HDA*" );
  }

 {
  //! [safepbnf-usage]
  const auto m = GridMap::random( 300 , 200 , 0.25 , 2 );
  // from the bottom-left to the bottom-right corner, with nblocks of
  // 10 x 10 cells (100 is the side of the blocks of the other searches)
  const Grid< 4 > g( m , 0 , 199 , 299 , 199 , 100 , 10 );
  PBNFOptions o;
  o.threads = 4;
  o.min_expansions = 32;   // expansions in an nblock before a switch
  const auto r = SafePBNF< Grid< 4 > >( g , o ).solve();
  // r.status, r.cost, r.path and r.stats as for A*
  //! [safepbnf-usage]
  const auto ref = astar( Blind< Grid< 4 > >{ g } );
  if( ref.solved() )
   check_solution( g , r , ref.cost , "manual: SafePBNF on a grid" );
  else
   check( r.status == Status::unsolvable ,
          "manual: SafePBNF on an unsolvable grid" );
  check( r.stats.expanded > 0 , "manual: SafePBNF expands" );
  }

 {
  //! [branch-and-bound-usage]
  // a random 0-1 knapsack of 50 items, of class 0
  const Knapsack kp = Knapsack::random( 0 , 50 , 1000 , 1 );
  BBOptions o;
  o.selection = Selection::best_first;     // the least bound first
  const auto s = branch_and_bound( kp , o );
  // s.status, s.value and the incumbent node s.node
  ParallelBBOptions po;
  po.threads = 4;
  const BBScheme sch = bb_scheme( "stealing" );   // a name, e.g., from a file
  const auto r = parallel_branch_and_bound( kp , sch , po );
  bool refused = false;
  try {
   ( void ) bb_scheme( "nonesuch" );
   }
  catch( const std::invalid_argument & ) {     // not the name of a scheme
   refused = true;
   }
  //! [branch-and-bound-usage]
  check( ( s.status == BBStatus::optimal ) &&
         ( s.value == double( kp.optimum() ) ) && s.node &&
         kp.check( *s.node , s.value ) , "manual: the sequential B&B" );
  check( ( r.status == BBStatus::optimal ) && ( r.value == s.value ) ,
         "manual: the B&B by work stealing, chosen by name" );
  check( refused , "manual: no scheme of an unknown name" );
  }

 {
  //! [domains-usage]
  // the 8-puzzle towards the goal with the blank in the cell 0
  const SlidingTile< 3 , 3 > p( { 7 , 2 , 4 , 5 , 0 , 6 , 8 , 3 , 1 } );
  const auto sp = astar( p );                   // the buckets: int costs
  const auto sh = astar( p , OpenKind::heap );
  // a random grid of 300 x 200 cells, each blocked with probability 1/4,
  // from the bottom-left to the bottom-right corner
  const auto m = GridMap::random( 300 , 200 , 0.25 , 1 );
  const Grid< 4 > g( m , 0 , 199 , 299 , 199 );    // 4 moves, int costs
  const auto sg = astar( g );
  const Grid< 8 > g8( m , 0 , 199 , 299 , 199 );   // 8 moves, double costs
  const auto s8 = astar( g8 );
  //! [domains-usage]
  const auto dist = eight_distances( p );
  check( sp.solved() && ( sp.cost == 26 ) && p.is_goal( sp.path.back() ) &&
         ( sh.cost == 26 ) && ( dist.at( p.initial_state() ) == 26 ) ,
         "manual: the 8-puzzle in 26 moves, both OPEN lists" );
  const auto ref = astar( Blind< Grid< 4 > >{ g } );
  check( ( sg.solved() == ref.solved() ) &&
         ( ! sg.solved() || ( sg.cost == ref.cost ) ) ,
         "manual: A* on the grid with 4 moves against Dijkstra" );
  check( ! ( sg.solved() && s8.solved() ) || ( s8.cost <= sg.cost ) ,
         "manual: 8 moves cost no more than 4" );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_pbnf_options< ADS::ThreadPool >( const std::string & );

/*--------------------------------------------------------------------------*/
/*-------------------------- End File manual.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
