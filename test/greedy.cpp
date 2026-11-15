/*--------------------------------------------------------------------------*/
/*---------------------------- File greedy.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the greedy and the
 * other searches.
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
/*------------------- THE GREEDY AND THE OTHER SEARCHES --------------------*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the two evaluations of the greedy searches

static const GreedyEvaluation evaluations[] = { GreedyEvaluation::eager ,
                                                GreedyEvaluation::lazy };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a Domain whose states in dead are dead ends, h = infinite_cost()

template< class D >
struct DeadEnds {

 using state_type = state_t< D >;
 using cost_type = cost_t< D >;

 const D & d;
 std::set< state_type > dead;

 state_type initial_state( void ) const { return( d.initial_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 cost_type heuristic( const state_type & s ) const {
  return( dead.count( s ) ? infinite_cost< cost_type >() : d.heuristic( s ) );
  }
 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_type & t , cost_type c ,
                               cost_type ) { f( t , c , heuristic( t ) ); } );
  }

 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a Domain with the successors of D but no heuristic computed when they
// are generated (a DeferredDomain), counting the evaluations

template< class D >
struct Deferred {

 using state_type = state_t< D >;
 using cost_type = cost_t< D >;

 const D & d;
 mutable std::uint64_t evaluations = 0;

 state_type initial_state( void ) const { return( d.initial_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 cost_type heuristic( const state_type & s ) const {
  ++evaluations;
  return( d.heuristic( s ) );
  }
 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  d.for_each_successor( s , d.heuristic( s ) , std::forward< F >( f ) );
  }
 template< class F >
 void for_each_successor_state( const state_type & s , F && f ) const {
  d.for_each_successor( s , 0 , [ & ]( const state_type & t , cost_type c ,
                                       cost_type ) { f( t , c ); } );
  }

 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// 200 states of the 8-puzzle with their distances to the goal

static std::vector< std::pair< std::uint64_t , int > > eight_sample( void )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 std::vector< std::pair< std::uint64_t , int > > v;
 std::size_t k = 0;
 for( const auto & [ s , dst ] : dist )
  if( ( k++ % 907 ) == 0 )
   v.emplace_back( s , dst );
 std::sort( v.begin() , v.end() );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// GBFS and the search by g + w h, eager and lazy: valid paths, the bound
// of w, the bound of the costs, the dead ends, the limits

void test_gbfs( const std::string & data )
{
 const auto sample = eight_sample();
 const auto inst = korf( data );
 for( GreedyEvaluation ev : evaluations ) {
  const std::string e = ev == GreedyEvaluation::eager ? "eager" : "lazy";
  GreedyOptions o;
  o.evaluation = ev;
  int bad = 0;
  // the bound w C* holds with eager evaluation (with lazy evaluation the
  // key of a successor has the h of its parent, and no bound is known)
  const bool eager = ev == GreedyEvaluation::eager;
  for( const auto & [ s , dst ] : sample ) {
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const auto g = gbfs( d , o );
   GreedyOptions w = o;
   w.weight = 1;
   w.reopen = true;
   const auto a = gbfs( d , w );
   w.weight = 2;
   const auto b = gbfs( d , w );
   if( ! ( g.solved() && a.solved() && b.solved() &&
           ( ( ! eager ) || ( ( a.cost == dst ) && ( b.cost <= 2 * dst ) ) )
           ) )
    ++bad;
   }
  check( bad == 0 , "GBFS " + e + ( eager ? ", w = 1 optimal and w = 2 "
                                             "within 2 C*" : ", w = 1 , 2" ) +
                    " on 200 states of the 8-puzzle" );
  for( int id : easy_korf ) {
   const FifteenPuzzle d( inst[ id - 1 ].board );
   check_path( d , gbfs( d , o ) , "GBFS " + e + ", Korf #" +
                                   std::to_string( id ) );
   GreedyOptions w = o;
   w.weight = 3;
   w.reopen = true;
   check_path( d , gbfs( d , w ) , "wA* w = 3 " + e + ", Korf #" +
                                   std::to_string( id ) ,
               eager ? 3.0 * inst[ id - 1 ].cost
                     : infinite_cost< double >() );
   }
  check_astar_edges( o , []( const EightPuzzle & d , const GreedyOptions & x )
                     { return( gbfs( d , x ) ); } , "GBFS " + e );
  for( City c : { Arad , Oradea , Neamt } ) {
   const auto g = romania( c );
   check_path( g , gbfs( g , o ) , "GBFS " + e + " on Romania" );
   }
  // the bound: g + h >= C* prunes the optimal path, C* + 1 does not
  const FifteenPuzzle d( inst[ easy_korf[ 0 ] - 1 ].board );
  const int cs = inst[ easy_korf[ 0 ] - 1 ].cost;
  GreedyOptions w = o;
  w.weight = 1;
  w.reopen = true;
  w.prune_by_f = true;
  w.bound = cs;
  check( gbfs( d , w ).status == Status::unsolvable ,
         "GBFS " + e + ", the bound C* prunes every path" );
  w.bound = cs + 1;
  const auto x = gbfs( d , w );
  check( x.solved() && ( x.cost == cs ) , "GBFS " + e + ", the bound C* + 1"
                                          " keeps the optimal one" );
  }

 // dead ends: the states next to the start but one are dead, the path
 // avoids them
 const EightPuzzle p( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 } ) );
 DeadEnds< EightPuzzle > de{ p , {} };
 int kept = 0;
 p.for_each_successor( p.initial_state() , 0 ,
                       [ & ]( std::uint64_t t , int , int ) {
  if( kept++ > 0 )
   de.dead.insert( t );
  } );
 for( GreedyEvaluation ev : evaluations ) {
  GreedyOptions o;
  o.evaluation = ev;
  const auto s = gbfs( de , o );
  bool avoids = s.solved();
  for( auto t : s.path )
   avoids = avoids && ! de.dead.count( t );
  check( avoids , "GBFS avoids the dead ends" );
  }

 // the deferred evaluation evaluates only the states selected
 {
  const FifteenPuzzle d( inst[ easy_korf[ 1 ] - 1 ].board );
  Deferred< FifteenPuzzle > df{ d };
  GreedyOptions o;
  o.evaluation = GreedyEvaluation::lazy;
  const auto s = gbfs( df , o );
  const std::uint64_t ev = df.evaluations;
  check( s.stats.evaluated == ev ,
         "lazy GBFS on a DeferredDomain evaluates what it counts (" +
         std::to_string( ev ) + ")" );
  check_path( df , s , "lazy GBFS on a DeferredDomain" );
  check( s.stats.evaluated < s.stats.generated ,
         "lazy GBFS evaluates fewer states than it generates" );
  }

 // a graph with no path to the goal
 ExplicitGraph< int > g( 4 , 0 , 3 );
 g.add_edge( 0 , 1 , 1 );
 g.add_edge( 1 , 2 , 1 );
 for( GreedyEvaluation ev : evaluations ) {
  GreedyOptions o;
  o.evaluation = ev;
  const auto s = gbfs( g , o );
  check( ( s.status == Status::unsolvable ) && ( s.stats.expanded == 3 ) ,
         "GBFS on a graph with no path" );
  }

 // goal at generation saves the last expansion
 {
  const FifteenPuzzle d( inst[ easy_korf[ 2 ] - 1 ].board );
  GreedyOptions o;
  const auto a = gbfs( d , o );
  o.goal_at_generation = true;
  const auto b = gbfs( d , o );
  check_path( d , b , "GBFS with the goal test at the generation" );
  check( b.stats.expanded < a.stats.expanded ,
         "the goal test at the generation saves expansions" );
  }

 check( [] {
  try {
   GreedyOptions o;
   o.queues = { QueueSpec{ QueueKind::best , { 2 } , false } };
   const EightPuzzle d( walk8( { 1 } ) );
   static_cast< void >( gbfs( d , o ) );
   return( false );
   }
  catch( std::invalid_argument & ) {
   return( true );
   }
  }() , "a queue on a heuristic out of range is refused" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the queues of OPEN: alternation, preferred successors and boosting,
// type-based, epsilon-greedy, Pareto, the breaking of the ties

void test_gbfs_queues( const std::string & data )
{
 const auto inst = korf( data );
 struct Config {
  const char * name;
  std::vector< QueueSpec > queues;
  int boost;
  TieBreak tie;
  };
 const std::vector< Config > configs = {
  { "best FIFO" , { { QueueKind::best , { 0 } , false } } , 0 ,
    TieBreak::fifo } ,
  { "best LIFO" , { { QueueKind::best , { 0 } , false } } , 0 ,
    TieBreak::lifo } ,
  { "best low g" , { { QueueKind::best , { 0 } , false } } , 0 ,
    TieBreak::low_g } ,
  { "best high g" , { { QueueKind::best , { 0 } , false } } , 0 ,
    TieBreak::high_g } ,
  { "best random" , { { QueueKind::best , { 0 } , false } } , 0 ,
    TieBreak::random } ,
  { "lexicographic h0, h1" , { { QueueKind::best , { 0 , 1 } , false } } ,
    0 , TieBreak::fifo } ,
  { "alternation h0, h1" , { { QueueKind::best , { 0 } , false } ,
                             { QueueKind::best , { 1 } , false } } , 0 ,
    TieBreak::fifo } ,
  { "preferred, boost 1000" , { { QueueKind::best , { 0 } , false } ,
                                { QueueKind::best , { 0 } , true } } ,
    1000 , TieBreak::fifo } ,
  { "type-based" , { { QueueKind::best , { 0 } , false } ,
                     { QueueKind::type , { 0 } , false } } , 0 ,
    TieBreak::fifo } ,
  { "epsilon-greedy" , { { QueueKind::epsilon , { 0 } , false } } , 0 ,
    TieBreak::fifo } ,
  { "Pareto h0, h1" , { { QueueKind::pareto , { 0 , 1 } , false } } , 0 ,
    TieBreak::fifo } ,
  { "Pareto and preferred" , { { QueueKind::pareto , { 0 , 1 } , false } ,
                               { QueueKind::best , { 0 } , true } } , 1000 ,
    TieBreak::fifo } };
 for( const Config & c : configs )
  for( GreedyEvaluation ev : evaluations ) {
   const std::string n = std::string( c.name ) +
                         ( ev == GreedyEvaluation::eager ? ", eager"
                                                           : ", lazy" );
   GreedyOptions o;
   o.evaluation = ev;
   o.queues = c.queues;
   o.boost = c.boost;
   o.tie = c.tie;
   o.preferred_first = c.boost > 0;
   for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] , easy_korf[ 2 ] } ) {
    const FifteenPuzzle d( inst[ id - 1 ].board );
    const auto a = gbfs( d , o );
    if( id == easy_korf[ 0 ] ) {
     check_path( d , a , n + ", Korf #" + std::to_string( id ) );
     const auto b = gbfs( d , o );
     check( ( a.stats.expanded == b.stats.expanded ) && ( a.cost == b.cost ) ,
            n + ", the same seed gives the same search" );
     }
    else
     if( ! a.solved() )
      check_path( d , a , n + ", Korf #" + std::to_string( id ) );
    }
   bool one = true;                       // the grid has one heuristic
   for( const QueueSpec & q : c.queues )
    for( std::size_t i : q.heuristics )
     one = one && ( i == 0 );
   if( ! one )
    continue;
   GridMap m = GridMap::random( 40 , 40 , 0.25 , 7 , 5 );
   const Grid< 8 > g( m , 0 , 39 , 39 , 39 );
   check_path( g , gbfs( g , o ) , n + ", a random grid of costs 1 to 5" );
   }
 // the epsilon queue with epsilon 1 is a random walk on OPEN, with 0 the
 // best queue
 {
  const FifteenPuzzle d( inst[ easy_korf[ 3 ] - 1 ].board );
  GreedyOptions o;
  const auto a = gbfs( d , o );
  o.queues = { QueueSpec{ QueueKind::epsilon , { 0 } , false } };
  o.epsilon = 0;
  const auto b = gbfs( d , o );
  check( a.stats.expanded == b.stats.expanded ,
         "epsilon-greedy with epsilon 0 is GBFS" );
  o.epsilon = 1;                        // a random walk on OPEN
  const EightPuzzle e( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 } ) );
  check_path( e , gbfs( e , o ) , "epsilon-greedy with epsilon 1" );
  }
 // a queue of the preferred successors only: on the 8-puzzle the
 // preferred moves lower the Manhattan distance
 {
  const EightPuzzle d( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 } ) );
  int good = 0 , all = 0;
  const auto pref = d.preferences( d.initial_state() );
  d.for_each_successor( d.initial_state() , d.heuristic( d.initial_state() )
                        , [ & ]( std::uint64_t t , int , int ht ) {
   ++all;
   good += ( bool( pref( t ) ) ==
             ( ht < d.heuristic( d.initial_state() ) ) ) ? 1 : 0;
   } );
  check( ( all > 0 ) && ( good == all ) ,
         "the preferred moves of the puzzle lower the Manhattan distance" );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// enforced hill-climbing

void test_ehc( const std::string & data )
{
 const auto inst = korf( data );
 const auto sample = eight_sample();
 for( PreferredUsage u : { PreferredUsage::ignore , PreferredUsage::prune ,
                           PreferredUsage::first } ) {
  const std::string n = u == PreferredUsage::ignore ? "EHC" :
                        u == PreferredUsage::prune ? "EHC, preferred only"
                                                   : "EHC, preferred first";
  EHCOptions o;
  o.preferred = u;
  o.fallback = true;
  int bad = 0;
  for( const auto & [ s , dst ] : sample ) {
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const auto x = ehc( d , o );
   if( ! x.solved() )
    ++bad;
   else
    if( dst > 0 )
     check_path( d , x , n + " on the 8-puzzle" );
   }
  check( bad == 0 , n + " with the fallback solves 200 states of the "
                        "8-puzzle" );
  for( int id : { 55 , 42 } ) {
   const FifteenPuzzle d( inst[ id - 1 ].board );
   check_path( d , ehc( d , o ) , n + ", Korf #" + std::to_string( id ) );
   }
  }
 // without the fallback, the plateau of a graph without a better state
 ExplicitGraph< int > g( 4 , 0 , 3 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 2 , 1 , true );
 g.add_edge( 2 , 3 , 1 , true );
 g.set_heuristic( 0 , 1 );
 g.set_heuristic( 1 , 5 );
 g.set_heuristic( 2 , 5 );
 const auto x = ehc( g );
 check_path( g , x , "EHC across a plateau of h = 5" );
 ExplicitGraph< int > dead( 3 , 0 , 2 );
 dead.add_edge( 0 , 1 , 1 , true );
 dead.set_heuristic( 0 , 2 );
 dead.set_heuristic( 1 , 3 );
 check( ehc( dead ).status == Status::unsolvable , "EHC fails at a dead end"
        );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the iterated searches: LAMA reaches the optimal cost, the phases improve
// and share the evaluations

void test_iterated( const std::string & data )
{
 const auto inst = korf( data );
 const auto sample = eight_sample();
 int bad = 0;
 for( std::size_t i = 0 ; i < sample.size() ; i += 10 ) {
  const EightPuzzle d( EightPuzzle::unpack( sample[ i ].first ) );
  IteratedSearch< EightPuzzle > it( d , lama_phases( 2 , true ) );
  const auto x = it.solve();
  bool ok = x.solved() && ( x.cost == sample[ i ].second );
  int last = std::numeric_limits< int >::max();
  for( const auto & p : it.phases() )
   if( p.cost != infinite_cost< int >() ) {
    ok = ok && ( p.cost < last );
    last = p.cost;
    }
  bad += ok ? 0 : 1;
  }
 check( bad == 0 , "LAMA reaches the optimal cost on 20 states of the "
                   "8-puzzle, each phase better than the one before" );
 const FifteenPuzzle d( inst[ 54 ].board );
 auto o = lama_phases( 2 , true );
 o.max_phases = 3;
 IteratedSearch< FifteenPuzzle > a( d , o );
 const auto x = a.solve();
 check_path( d , x , "LAMA with 3 phases, Korf #55" );
 o.share_cache = false;
 IteratedSearch< FifteenPuzzle > b( d , o );
 const auto y = b.solve();
 check( y.cost == x.cost , "the cache does not change the paths" );
 check( x.stats.evaluated < y.stats.evaluated ,
        "the phases share the evaluations (" +
        std::to_string( x.stats.evaluated ) + " instead of " +
        std::to_string( y.stats.evaluated ) + ")" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the width-based searches

void test_width( const std::string & data )
{
 // the novelty table against all the tuples recorded by brute force
 {
  std::mt19937_64 rng( 3 );
  for( unsigned k = 1 ; k <= 3 ; ++k ) {
   NoveltyTable nt( 12 , k );
   std::set< std::vector< std::size_t > > seen;
   int bad = 0;
   for( int it = 0 ; it < 300 ; ++it ) {
    std::vector< std::size_t > F;
    for( std::size_t i = 0 ; i < 12 ; ++i )
     if( random_below( rng , 3 ) == 0 )
      F.push_back( i );
    unsigned w = k + 1;
    for( unsigned m = 1 ; m <= k ; ++m ) {
     if( F.size() < m )
      break;
     std::vector< bool > pick( F.size() , false );
     std::fill( pick.begin() , pick.begin() + m , true );
     do {
      std::vector< std::size_t > t;
      for( std::size_t i = 0 ; i < F.size() ; ++i )
       if( pick[ i ] )
        t.push_back( F[ i ] );
      if( seen.insert( t ).second )
       w = std::min( w , m );
      } while( std::prev_permutation( pick.begin() , pick.end() ) );
     }
    if( nt.record( F , F ) != w )
     ++bad;
    }
   check( bad == 0 , "the novelty of width " + std::to_string( k ) +
                     " against brute force" );
   }
  }
 // IW(2) on a grid of 4 moves never prunes (a cell is a pair of features)
 // and is a breadth-first search, optimal on unit costs; SIW and BFWS
 // find a path
 for( std::uint64_t seed : { 1 , 2 , 3 } ) {
  GridMap m = GridMap::random( 50 , 50 , 0.2 , seed );
  const Grid< 4 > g( m , 0 , 49 , 49 , 49 );
  const auto a = astar( g );
  if( ! a.solved() )
   continue;
  WidthOptions o;
  o.width = 2;
  const auto w2 = iw( g , o );
  check_solution( g , w2 , a.cost , "IW(2) on a grid is breadth-first" );
  o.width = 1;
  const auto w1 = iw( g , o );
  check( ( ! w1.solved() ) || ( w1.cost >= a.cost ) ,
         "IW(1) on a grid prunes" );
  o.width = 2;
  check_path( g , siw( g , o ) , "SIW on a random grid" );
  check_path( g , bfws( g , o ) , "BFWS on a random grid" );
  o.use_heuristic = true;
  check_path( g , bfws( g , o ) , "BFWS on h on a random grid" );
  }
 // BFWS on the puzzles; SIW on the 8-puzzle with a large width
 const auto inst = korf( data );
 for( int id : { 55 , 42 } ) {
  const FifteenPuzzle d( inst[ id - 1 ].board );
  WidthOptions o;
  check_path( d , bfws( d , o ) , "BFWS on #g, Korf #" +
                                  std::to_string( id ) );
  o.use_heuristic = true;
  check_path( d , bfws( d , o ) , "BFWS on h, Korf #" + std::to_string( id ) );
  o.prune = true;
  o.width = 2;
  const auto p = bfws( d , o );
  check( ( ! p.solved() ) || ( p.cost >= inst[ id - 1 ].cost ) ,
         "k-BFWS prunes" );
  }
 const EightPuzzle e( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                               0 , 3 , 4 , 7 , 6 } ) );
 WidthOptions o;
 o.width = 8;
 SerializedWidth< EightPuzzle > s( e , o );
 check_path( e , s.solve() , "SIW on the 8-puzzle, width up to 8" );
 check( s.calls() >= 1 , "SIW calls IW" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// PEA* (with and without skip_known) and EPEA* on a grid: the cost of A*,
// and at most old selections, those of PEA* before it left out the known
// successors, compared f up to cost_slack() and skipped the parent

template< class G >
static void check_partial( const G & g , std::uint64_t old ,
                           const std::string & what )
{
 const auto a = astar( g );
 PartialOptions keep;
 keep.skip_known = false;
 const auto p = pea_star( g );
 const auto q = pea_star( g , keep );
 const auto e = epea_star( g );
 check_solution( g , p , a.cost , "PEA*" + what );
 check_solution( g , q , a.cost , "PEA* keeping the known successors" +
                                  what );
 check_solution( g , e , a.cost , "EPEA*" + what );
 check( ( p.stats.expanded <= q.stats.expanded ) &&
        ( q.stats.expanded <= old ) && ( e.stats.expanded <= old ) &&
        ( e.stats.generated < p.stats.generated ) ,
        "PEA* and EPEA* select at most " + std::to_string( old ) +
        " nodes, skip_known no more, EPEA* generates less" + what );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// partial expansion, Fringe Search, RBFS and beam search against A*

void test_classic( const std::string & data )
{
 const auto sample = eight_sample();
 int bad = 0;
 for( const auto & [ s , dst ] : sample ) {
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  PartialOptions c;
  c.cutoff = 2;
  BeamOptions bo;
  bo.width = std::numeric_limits< std::size_t >::max();
  for( int x : { pea_star( d ).cost , pea_star( d , c ).cost ,
                 epea_star( d ).cost , fringe_search( d ).cost ,
                 rbfs( d ).cost , beam_search( d , bo ).cost } )
   bad += ( x == dst ) ? 0 : 1;
  }
 check( bad == 0 , "PEA*, PEA* (C = 2), EPEA*, Fringe, RBFS and the beam "
                   "of infinite width optimal on 200 states of the "
                   "8-puzzle" );
 const auto inst = korf( data );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto a = astar( d );
  const auto p = pea_star( d );
  const auto e = epea_star( d );
  const std::string n = ", Korf #" + std::to_string( id );
  check_solution( d , p , ki.cost , "PEA*" + n );
  check_solution( d , e , ki.cost , "EPEA*" + n );
  check( ( e.stats.generated < a.stats.generated ) &&
         ( e.stats.stored < a.stats.stored ) &&
         ( p.stats.stored < a.stats.stored ) ,
         "EPEA* generates and stores less than A*, PEA* stores less" + n );
  check_solution( d , fringe_search( d ) , ki.cost , "Fringe" + n );
  if( id == easy_korf[ 0 ] )
   check_solution( d , rbfs( d ) , ki.cost , "RBFS" + n );
  BeamOptions bo;
  bo.width = 2000;
  check_path( d , beam_search( d , bo ) , "beam of 2000" + n );
  }
 for( std::uint64_t seed : { 1 , 2 } ) {
  GridMap m = GridMap::random( 50 , 50 , 0.25 , seed , 4 );
  const Grid< 8 > g( m , 0 , 49 , 49 , 49 );
  const auto a = astar( g );
  if( ! a.solved() )
   continue;
  check_solution( g , pea_star( g ) , a.cost , "PEA* on a grid of 8 moves" );
  check_solution( g , epea_star( g ) , a.cost ,
                  "EPEA* on a grid of 8 moves, cells of cost 1 to 4" );
  check_solution( g , fringe_search( g ) , a.cost ,
                  "Fringe on a grid of 8 moves" );
  const Grid< 4 > g4( m , 0 , 49 , 49 , 49 );
  const auto a4 = astar( g4 );
  if( a4.solved() )
   check_solution( g4 , epea_star( g4 ) , a4.cost ,
                   "EPEA* on a grid of 4 moves, cells of cost 1 to 4" );
  }
 static constexpr std::uint64_t old4[] = { 466 , 750 , 574 };
 static constexpr std::uint64_t old8[] = { 752 , 2167 , 1495 };
 for( int k = 0 ; k < 3 ; ++k ) {
  const auto m = GridMap::random( 60 , 60 , 0.25 , std::uint64_t( k + 1 ) );
  const std::string n = " on the random map " + std::to_string( k + 1 ) +
                        " of 60 x 60 cells, ";
  check_partial( Grid< 4 >( m , 0 , 59 , 59 , 59 ) , old4[ k ] ,
                 n + "4 moves" );
  check_partial( Grid< 8 >( m , 0 , 59 , 59 , 59 ) , old8[ k ] ,
                 n + "8 moves" );
  }
 for( std::uint64_t seed : { 1 , 2 } ) {
  GridMap m = GridMap::random( 8 , 8 , 0.2 , seed );
  const Grid< 4 > g4( m , 0 , 7 , 7 , 7 );
  const auto a = astar( g4 );
  if( a.solved() )
   check_solution( g4 , rbfs( g4 ) , a.cost , "RBFS on a grid of 4 moves" );
  }
 // the limits
 const EightPuzzle d( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                               0 , 3 , 4 , 7 , 6 } ) );
 SearchOptions l;
 l.node_limit = 5;
 PartialOptions pl;
 pl.node_limit = 5;
 check( ( pea_star( d , pl ).status == Status::node_limit ) &&
        ( epea_star( d , l ).status == Status::node_limit ) &&
        ( fringe_search( d , l ).status == Status::node_limit ) &&
        ( rbfs( d , l ).status == Status::node_limit ) ,
        "the node limit of PEA*, EPEA*, Fringe and RBFS" );
 ExplicitGraph< int > g( 4 , 0 , 3 );
 g.add_edge( 0 , 1 , 1 );
 g.add_edge( 1 , 2 , 1 );
 check( ( pea_star( g ).status == Status::unsolvable ) &&
        ( fringe_search( g ).status == Status::unsolvable ) &&
        ( beam_search( g ).status == Status::unsolvable ) ,
        "PEA*, Fringe and beam on a graph with no path" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// LRTA* and LSS-LRTA*: the trials converge to an optimal path, and H stays
// admissible

void test_realtime( const std::string & )
{
 for( std::uint64_t seed : { 1 , 2 , 3 } ) {
  GridMap m = GridMap::random( 30 , 30 , 0.2 , seed );
  const Grid< 4 > g( m , 0 , 29 , 29 , 29 );
  const auto a = astar( g );
  if( ! a.solved() )
   continue;
  for( std::size_t L : { 1 , 8 , 64 } ) {
   RealTimeOptions o;
   o.lookahead = L;
   o.trials = 2000;
   o.until_convergence = true;
   RealTimeSearch< Grid< 4 > > r( g , o );
   const auto x = r.solve();
   const std::string n = ( L == 1 ? std::string( "LRTA*" )
                                  : "LSS-LRTA*, L = " + std::to_string( L ) )
                         + ", seed " + std::to_string( seed );
   check_solution( g , x , a.cost , n + ", converged to an optimal path" );
   // H admissible on the path of the last trial
   bool adm = true;
   int rest = x.cost;
   for( std::size_t i = 0 ; i + 1 < x.path.size() ; ++i ) {
    adm = adm && ( r.learned( x.path[ i ] ) <= rest );
    rest -= 1;
    }
   check( adm , n + ", H admissible on the optimal path" );
   check( r.trials().front() >= a.cost , n + ", the first trial is not "
                                             "cheaper than optimal" );
   }
  }
 // 8 moves and cells of cost 1 to 3, whose moves cost what the cell left
 // costs: H against h*, computed by Dijkstra backward from the goal
 for( std::uint64_t seed : { 1 , 2 , 3 } ) {
  GridMap m = GridMap::random( 40 , 40 , 0.25 , seed , 3 );
  const Grid< 8 > g( m , 0 , 39 , 39 , 39 );
  const auto a = astar( g );
  if( ! a.solved() )
   continue;
  const std::size_t nc = 40 * 40;
  std::vector< std::vector< std::pair< std::uint32_t , double > > > back( nc );
  for( std::uint32_t u = 0 ; u < nc ; ++u )
   if( m.cost( int( u % 40 ) , int( u / 40 ) ) )
    g.for_each_successor( u , 0.0 , [ & ]( std::uint32_t v , double c ,
                                           double ) {
     back[ v ].emplace_back( u , c );
     } );
  std::vector< double > hs( nc , infinite_cost< double >() );
  using QI = std::pair< double , std::uint32_t >;
  std::priority_queue< QI , std::vector< QI > , std::greater< QI > > q;
  hs[ g.index( 39 , 39 ) ] = 0;
  q.push( { 0.0 , g.index( 39 , 39 ) } );
  while( ! q.empty() ) {
   const double d = q.top().first;
   const std::uint32_t v = q.top().second;
   q.pop();
   if( d != hs[ v ] )
    continue;
   for( const auto & [ u , c ] : back[ v ] )
    if( d + c < hs[ u ] ) {
     hs[ u ] = d + c;
     q.push( { d + c , u } );
     }
   }
  for( std::size_t L : { 1 , 8 , 32 } ) {
   RealTimeOptions o;
   o.lookahead = L;
   o.trials = 3000;
   o.until_convergence = true;
   RealTimeSearch< Grid< 8 > > r( g , o );
   const auto x = r.solve();
   const std::string n = ( L == 1 ? std::string( "LRTA*" )
                                  : "LSS-LRTA*, L = " + std::to_string( L ) )
                         + ", 8 moves, costs 1 to 3, seed " +
                         std::to_string( seed );
   check_solution( g , x , a.cost , n + ", converged to an optimal path" );
   int over = 0;
   for( std::uint32_t u = 0 ; u < nc ; ++u )
    if( ( hs[ u ] != infinite_cost< double >() ) &&
        ( r.learned( u ) > hs[ u ] * ( 1 + 1e-9 ) + 1e-9 ) )
     ++over;
   check( over == 0 , n + ", H admissible on every state (" +
                      std::to_string( over ) + " above h*)" );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// HPA*: valid paths, not cheaper than A*, close to it on open maps

void test_hpa( const std::string & )
{
 for( std::uint64_t seed : { 1 , 2 , 3 , 4 } ) {
  GridMap m = GridMap::random( 64 , 64 , 0.15 , seed );
  for( int C : { 8 , 16 } ) {
   HPAOptions o;
   o.cluster = C;
   HPAStar< 8 > h( m , o );
   HPAStar< 4 > h4( m , o );
   std::mt19937_64 rng( seed );
   int found = 0 , queries = 0;
   double worst = 1;
   for( int q = 0 ; q < 20 ; ++q ) {
    const int sx = int( random_below( rng , 64 ) ) ,
              sy = int( random_below( rng , 64 ) ) ,
              gx = int( random_below( rng , 64 ) ) ,
              gy = int( random_below( rng , 64 ) );
    if( ! m.cost( sx , sy ) || ! m.cost( gx , gy ) )
     continue;
    const Grid< 8 > g( m , sx , sy , gx , gy );
    const auto a = astar( g );
    if( ! a.solved() )
     continue;
    ++queries;
    const auto x = h.solve( sx , sy , gx , gy );
    if( ! x.solved() )
     continue;
    ++found;
    check_path( g , x , "HPA*, 8 moves, C = " + std::to_string( C ) ,
                infinite_cost< double >() );
    check( x.cost >= a.cost * ( 1 - 1e-9 ) , "HPA* not cheaper than A*" );
    worst = std::max( worst , x.cost / std::max( a.cost , 1.0 ) );
    const Grid< 4 > g4( m , sx , sy , gx , gy );
    const auto y = h4.solve( sx , sy , gx , gy );
    if( y.solved() )
     check( y.cost >= astar( g4 ).cost , "HPA*, 4 moves, not cheaper" );
    }
   std::ostringstream os;
   os << "HPA* finds " << found << " of " << queries << " queries, C = " << C
      << ", the worst ratio to A* " << worst;
   check( ( found * 10 >= queries * 9 ) && ( worst < 1.5 ) , os.str() );
   }
  }
 // the start in the goal and in the same cluster
 GridMap m( 20 , 20 );
 HPAStar< 4 > h( m );
 const auto x = h.solve( 3 , 3 , 3 , 3 );
 check( x.solved() && ( x.cost == 0 ) && ( x.path.size() == 1 ) ,
        "HPA* from the goal" );
 const auto y = h.solve( 1 , 1 , 4 , 5 );
 check( y.solved() && ( y.cost == 7 ) , "HPA* inside one cluster, optimal" );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File greedy.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
