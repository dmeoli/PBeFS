/*--------------------------------------------------------------------------*/
/*---------------------------- File pgreedy.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the parallel greedy
 * searches.
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
// the bench transition system of a finite Domain: the states that GBFS
// expands with some breaking of the ties [Heusner, Keller and Helmert,
// SoCS 2017, theorem 1], from definitions 3 to 6 of Shimoda and Fukunaga
// (AAAI 2025); the goals have no successors and are progress states

template< class D >
static std::set< state_t< D > > bench_transition_system( const D & d )
{
 using S = state_t< D >;
 using C = cost_t< D >;
 const C inf = infinite_cost< C >();
 // the states, by a breadth-first search, the goals not expanded
 std::map< S , std::size_t > id;
 std::vector< S > st;
 std::vector< std::vector< std::size_t > > succ , pred;
 std::deque< S > q;
 id[ d.initial_state() ] = 0;
 st.push_back( d.initial_state() );
 q.push_back( d.initial_state() );
 succ.emplace_back();
 while( ! q.empty() ) {
  const S s = q.front();
  q.pop_front();
  if( d.is_goal( s ) )
   continue;
  const std::size_t i = id[ s ];
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const S & t , C , C ) {
   auto [ it , fresh ] = id.emplace( t , st.size() );
   if( fresh ) {
    st.push_back( t );
    succ.emplace_back();
    q.push_back( t );
    }
   if( it->second != i )
    succ[ i ].push_back( it->second );
   } );
  }
 const std::size_t n = st.size();
 pred.assign( n , {} );
 for( std::size_t i = 0 ; i < n ; ++i )
  for( std::size_t j : succ[ i ] )
   pred[ j ].push_back( i );
 std::vector< C > h( n ) , hwm( n , inf );
 for( std::size_t i = 0 ; i < n ; ++i )
  h[ i ] = d.heuristic( st[ i ] );
 // the high-water marks, by a bottleneck Dijkstra from the goals
 using QI = std::pair< C , std::size_t >;
 std::priority_queue< QI , std::vector< QI > , std::greater< QI > > pq;
 for( std::size_t i = 0 ; i < n ; ++i )
  if( d.is_goal( st[ i ] ) ) {
   hwm[ i ] = h[ i ];
   pq.push( { hwm[ i ] , i } );
   }
 while( ! pq.empty() ) {
  const auto [ v , j ] = pq.top();
  pq.pop();
  if( v != hwm[ j ] )
   continue;
  for( std::size_t i : pred[ j ] ) {
   const C c = std::max( h[ i ] , v );
   if( c < hwm[ i ] ) {
    hwm[ i ] = c;
    pq.push( { c , i } );
    }
   }
  }
 const auto hwm_succ = [ & ]( std::size_t i ) {
  C m = inf;
  for( std::size_t j : succ[ i ] )
   m = std::min( m , hwm[ j ] );
  return( m );
  };
 const auto progress = [ & ]( std::size_t i ) {
  return( d.is_goal( st[ i ] ) ||
          ( ( hwm[ i ] != inf ) && ( hwm[ i ] > hwm_succ( i ) ) ) );
  };
 // the benches, from that of the initial state
 std::set< S > out;
 std::set< std::size_t > done;
 std::deque< std::size_t > todo{ 0 };
 while( ! todo.empty() ) {
  const std::size_t s = todo.front();
  todo.pop_front();
  if( ! done.insert( s ).second )
   continue;
  const C level = hwm_succ( s );
  out.insert( st[ s ] );
  std::set< std::size_t > inner;
  std::deque< std::size_t > w{ s };
  std::set< std::size_t > exits;
  while( ! w.empty() ) {
   const std::size_t u = w.front();
   w.pop_front();
   for( std::size_t v : succ[ u ] ) {
    if( v == s )
     continue;
    if( progress( v ) ) {
     if( h[ v ] == level )
      exits.insert( v );
     continue;
     }
    if( ( h[ v ] <= level ) && inner.insert( v ).second )
     w.push_back( v );
    }
   }
  for( std::size_t v : inner )
   out.insert( st[ v ] );
  for( std::size_t v : exits ) {
   out.insert( st[ v ] );
   if( ! d.is_goal( st[ v ] ) )
    todo.push_back( v );
   }
  }
 return( out );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the states committed by OBAT, in order, are an order in which
// GBFS may expand them [theorem 2 of Shimoda and Fukunaga, AAAI 2025]

template< class D >
static bool gbfs_order( const D & d , const std::vector< state_t< D > > & c )
{
 using S = state_t< D >;
 std::multimap< cost_t< D > , S > open;
 std::set< S > closed{ d.initial_state() };
 open.emplace( d.heuristic( d.initial_state() ) , d.initial_state() );
 for( const S & s : c ) {
  if( open.empty() )
   return( false );
  const auto least = open.begin()->first;
  auto it = open.end();
  for( auto x = open.equal_range( least ).first ;
       x != open.equal_range( least ).second ; ++x )
   if( x->second == s )
    it = x;
  if( it == open.end() )
   return( false );
  open.erase( it );
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const S & t , cost_t< D > , cost_t< D > ht ) {
   if( ( ht != infinite_cost< cost_t< D > >() ) &&
       closed.insert( t ).second )
    open.emplace( ht , t );
   } );
  }
 return( true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a random directed graph of n nodes with integer h, goal n - 1 without
// successors and h( goal ) = 0

static ExplicitGraph< int > random_bench_graph( std::uint32_t n ,
                                                std::uint64_t seed )
{
 std::mt19937_64 rng( seed );
 ExplicitGraph< int > g( n , 0 , n - 1 );
 for( std::uint32_t u = 0 ; u + 1 < n ; ++u ) {
  const int k = 1 + int( random_below( rng , 3 ) );
  for( int e = 0 ; e < k ; ++e ) {
   const auto v = std::uint32_t( random_below( rng , n ) );
   if( v != u )
    g.add_edge( u , v , 1 , true );
   }
  g.set_heuristic( u , 1 + int( random_below( rng , 6 ) ) );
  }
 g.set_heuristic( n - 1 , 0 );
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the names of the rules

static const char * rule_name( GreedyRule r )
{
 static const char * n[] = { "KPGBFS" , "PUHF" , "PUHF2" , "PUHF3" ,
                             "PUHF4" , "SPUHF" , "OBAT" };
 return( n[ int( r ) ] );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the guarantees of the parallel greedy searches, in synchronous rounds:
// the BTS, OBAT against GBFS, the pathologies of KPGBFS and PUHF

void test_pgbfs_rounds( const std::string & data )
{
 const GreedyRule rules[] = { GreedyRule::kpgbfs , GreedyRule::puhf ,
                              GreedyRule::puhf2 , GreedyRule::puhf3 ,
                              GreedyRule::puhf4 , GreedyRule::spuhf ,
                              GreedyRule::obat };
 // one worker is GBFS with the first path kept, FIFO and LIFO (PUHF4
 // prefers the certain nodes among those of equal h, another breaking of
 // the ties)
 const auto inst = korf( data );
 for( int id : easy_korf ) {
  const FifteenPuzzle d( inst[ id - 1 ].board );
  for( TieBreak t : { TieBreak::fifo , TieBreak::lifo } ) {
   GreedyOptions g;
   g.improve_open = false;
   g.tie = t;
   const auto x = gbfs( d , g );
   for( GreedyRule r : rules )
    for( bool sge : { false , true } ) {
     if( r == GreedyRule::puhf4 )
      continue;
     ParallelGreedyOptions o;
     o.rule = r;
     o.rounds = true;
     o.sge = sge;
     o.tie = t;
     const auto y = ParallelGreedy< FifteenPuzzle >( d , o ).solve();
     if( ( y.stats.expanded != x.stats.expanded ) || ( y.cost != x.cost ) )
      check( false , std::string( rule_name( r ) ) + ( sge ? " with SGE" : "" )
             + ", one worker is GBFS, Korf #" + std::to_string( id ) );
     }
   }
  }
 check( failures == 0 , "every rule but PUHF4 with one worker is GBFS on 6"
                        " instances of Korf" );

 // the BTS of random graphs, against GBFS with every breaking of the ties
 int checked = 0 , outside = 0 , gbfs_out = 0 , obat_bad = 0 , def_bad = 0;
 for( std::uint64_t seed = 1 ; seed <= 60 ; ++seed ) {
  const auto g = random_bench_graph( 40 , seed );
  if( ! astar( g ).solved() )
   continue;
  ++checked;
  const auto bts = bench_transition_system( g );
  // GBFS, FIFO and LIFO, as a parallel search of one worker (which is
  // GBFS, see above) that records the states it expands
  for( TieBreak t : { TieBreak::fifo , TieBreak::lifo } ) {
   ParallelGreedyOptions p;
   p.rounds = true;
   p.tie = t;
   ParallelGreedy< ExplicitGraph< int > > y( g , p );
   y.set_record( true );
   static_cast< void >( y.solve() );
   for( auto u : y.trace().expanded )
    gbfs_out += bts.count( u ) ? 0 : 1;
   }
  for( GreedyRule r : rules )
   for( unsigned k : { 2u , 3u , 5u } )
    for( bool sge : { false , true } ) {
     ParallelGreedyOptions o;
     o.rule = r;
     o.rounds = true;
     o.threads = k;
     o.sge = sge;
     ParallelGreedy< ExplicitGraph< int > > x( g , o );
     x.set_record( true );
     const auto s = x.solve();
     if( ! s.solved() ) {
      check_path( g , s , std::string( rule_name( r ) ) +
                          " on a random graph" );
      continue;
      }
     if( r != GreedyRule::kpgbfs )
      for( auto u : x.trace().expanded )
       outside += bts.count( u ) ? 0 : 1;
     if( r == GreedyRule::obat ) {
      obat_bad += gbfs_order( g , x.trace().committed ) ? 0 : 1;
      std::map< int , int > per_h;
      for( auto u : x.trace().deferred )
       if( ++per_h[ g.heuristic( u ) ] > int( k ) )
        ++def_bad;
      if( x.trace().deferred.size() > k * s.path.size() )
       ++def_bad;
      }
     }
  }
 check( checked >= 20 , "random graphs with a path: " +
                        std::to_string( checked ) );
 check( gbfs_out == 0 , "GBFS expands only states of the BTS (the BTS of "
                        "the tester is consistent with theorem 1)" );
 check( outside == 0 , "PUHF to PUHF4, SPUHF and OBAT, with and without "
                       "SGE, 2 to 5 workers, expand only states of the BTS" );
 check( obat_bad == 0 , "OBAT commits the states in an order of GBFS "
                        "(theorem 2)" );
 check( def_bad == 0 , "OBAT leaves at most k states of each h in DEFERRED,"
                       " and at most k | p | (theorem 3)" );

 // the same on the 8-puzzle, whose space is enumerated
 {
  const EightPuzzle d( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 , 3 , 0 , 1 , 4 , 5 ,
                                8 } ) );
  const auto bts = bench_transition_system( d );
  int out = 0 , bad = 0;
  for( GreedyRule r : rules )
   for( unsigned k : { 2u , 4u , 8u } ) {
    ParallelGreedyOptions o;
    o.rule = r;
    o.rounds = true;
    o.threads = k;
    ParallelGreedy< EightPuzzle > x( d , o );
    x.set_record( true );
    const auto s = x.solve();
    check_path( d , s , std::string( rule_name( r ) ) + ", " +
                        std::to_string( k ) + " workers, the 8-puzzle" );
    if( r != GreedyRule::kpgbfs )
     for( auto u : x.trace().expanded )
      out += bts.count( u ) ? 0 : 1;
    if( r == GreedyRule::obat )
     bad += gbfs_order( d , x.trace().committed ) ? 0 : 1;
    }
  check( out == 0 , "the constrained rules expand only states of the BTS "
                    "of the 8-puzzle (" + std::to_string( bts.size() ) +
                    " states)" );
  check( bad == 0 , "OBAT commits in an order of GBFS on the 8-puzzle" );
  }

 // KPGBFS leaves the BTS, PUHF does not: s0 -> a (h 1), b (h 2), a -> goal,
 // b -> a chain of h 2 that leads nowhere
 {
  ExplicitGraph< int > g( 8 , 0 , 3 );
  g.set_heuristic( 0 , 3 );
  g.set_heuristic( 1 , 1 );
  g.set_heuristic( 2 , 2 );
  g.set_heuristic( 3 , 0 );
  g.add_edge( 0 , 1 , 1 , true );
  g.add_edge( 0 , 2 , 1 , true );
  g.add_edge( 1 , 3 , 1 , true );
  g.add_edge( 2 , 4 , 1 , true );
  for( std::uint32_t u = 4 ; u < 8 ; ++u ) {
   g.set_heuristic( u , 2 );
   if( u < 7 )
    g.add_edge( u , u + 1 , 1 , true );
   }
  const auto bts = bench_transition_system( g );
  ParallelGreedyOptions o;
  o.rounds = true;
  o.threads = 2;
  ParallelGreedy< ExplicitGraph< int > > k( g , o );
  k.set_record( true );
  static_cast< void >( k.solve() );
  int out = 0;
  for( auto u : k.trace().expanded )
   out += bts.count( u ) ? 0 : 1;
  check( out > 0 , "KPGBFS with 2 workers expands a state outside the BTS" );
  o.rule = GreedyRule::puhf;
  ParallelGreedy< ExplicitGraph< int > > p( g , o );
  p.set_record( true );
  static_cast< void >( p.solve() );
  out = 0;
  for( auto u : p.trace().expanded )
   out += bts.count( u ) ? 0 : 1;
  check( out == 0 , "PUHF with 2 workers does not" );
  }

 // PUHF explores two benches at once, OBAT one [figure 1 of Shimoda and
 // Fukunaga, AAAI 2025]: s0 (h 10) -> A, B (h 5), each the head of a line
 // of x states of h 4 to its own goal
 for( std::uint32_t x : { 10u , 100u , 1000u } ) {
  const std::uint32_t n = 3 + 2 * x + 1;
  ExplicitGraph< int > g( n , 0 , n - 1 );
  g.set_heuristic( 0 , 10 );
  g.set_heuristic( 1 , 5 );
  g.set_heuristic( 2 , 5 );
  g.add_edge( 0 , 1 , 1 , true );
  g.add_edge( 0 , 2 , 1 , true );
  for( std::uint32_t side = 0 ; side < 2 ; ++side ) {
   std::uint32_t prev = 1 + side;
   for( std::uint32_t i = 0 ; i < x ; ++i ) {
    const std::uint32_t v = 3 + side * x + i;
    g.set_heuristic( v , 4 );
    g.add_edge( prev , v , 1 , true );
    prev = v;
    }
   g.add_edge( prev , n - 1 , 1 , true );
   }
  // GBFS expands x + 3 states with every breaking of the ties (s0, A or B,
  // its line, the goal), and the path has x + 3 states
  const std::uint64_t ngbfs = x + 3 , path = x + 3;
  for( unsigned k : { 2u , 4u } ) {
   ParallelGreedyOptions o;
   o.rounds = true;
   o.threads = k;
   o.rule = GreedyRule::puhf;
   const auto p = ParallelGreedy< ExplicitGraph< int > >( g , o ).solve();
   o.rule = GreedyRule::obat;
   const auto b = ParallelGreedy< ExplicitGraph< int > >( g , o ).solve();
   std::ostringstream os;
   os << "x = " << x << ", " << k << " workers: GBFS " << ngbfs << ", PUHF "
      << p.stats.expanded << ", OBAT " << b.stats.expanded << " <= "
      << ngbfs + k * path;
   check( ( b.stats.expanded <= ngbfs + k * path ) &&
          ( p.stats.expanded >= ngbfs + x - 1 ) , os.str() );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel greedy searches on an executor: valid paths with every rule,
// the BTS under real threads, HDGBFS

template< class Executor >
void test_pgbfs( const std::string & data )
{
 const std::string ex = executor_name< Executor >();
 const GreedyRule rules[] = { GreedyRule::kpgbfs , GreedyRule::puhf ,
                              GreedyRule::puhf2 , GreedyRule::puhf3 ,
                              GreedyRule::puhf4 , GreedyRule::spuhf ,
                              GreedyRule::obat };
 const auto inst = korf( data );
 Executor pool( 8 );
 for( GreedyRule r : rules )
  for( bool sge : { false , true } )
   for( unsigned k : { 1u , 2u , 4u , 8u } ) {
    ParallelGreedyOptions o;
    o.rule = r;
    o.sge = sge;
    o.threads = k;
    const std::string n = std::string( rule_name( r ) ) +
                          ( sge ? " with SGE, " : ", " ) +
                          std::to_string( k ) + " workers on " + ex;
    for( int id : { 55 , 42 } ) {
     const FifteenPuzzle d( inst[ id - 1 ].board );
     ParallelGreedy< FifteenPuzzle , Executor > s( d , o );
     s.set_executor( pool );
     check_path( d , s.solve() , n + ", Korf #" + std::to_string( id ) );
     }
    GridMap m = GridMap::random( 60 , 60 , 0.25 , 5 , 3 );
    const Grid< 8 > g( m , 0 , 59 , 59 , 59 );
    ParallelGreedy< Grid< 8 > , Executor > s( g , o );
    s.set_executor( pool );
    check_path( g , s.solve() , n + ", a grid" );
    }
 // the BTS and the order of OBAT under real threads, on the 8-puzzle
 {
  const EightPuzzle d( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 , 3 , 0 , 1 , 4 , 5 ,
                                8 } ) );
  const auto bts = bench_transition_system( d );
  int out = 0 , bad = 0;
  for( int rep = 0 ; rep < 5 ; ++rep )
   for( GreedyRule r : rules ) {
    if( r == GreedyRule::kpgbfs )
     continue;
    for( bool sge : { false , true } ) {
     ParallelGreedyOptions o;
     o.rule = r;
     o.sge = sge;
     o.threads = 4;
     ParallelGreedy< EightPuzzle , Executor > x( d , o );
     x.set_executor( pool );
     x.set_record( true );
     static_cast< void >( x.solve() );
     for( auto u : x.trace().expanded )
      out += bts.count( u ) ? 0 : 1;
     if( r == GreedyRule::obat )
      bad += gbfs_order( d , x.trace().committed ) ? 0 : 1;
     }
    }
  check( out == 0 , "the constrained rules on " + ex + " expand only states"
                    " of the BTS of the 8-puzzle" );
  check( bad == 0 , "OBAT on " + ex + " commits in an order of GBFS" );
  }
 // HDGBFS
 for( unsigned k : { 1u , 2u , 4u } ) {
  ParallelGreedyOptions o;
  o.threads = k;
  const FifteenPuzzle d( inst[ 54 ].board );
  ZobristDistribution< FifteenPuzzle > z( d );
  HDGBFS< FifteenPuzzle , ZobristDistribution< FifteenPuzzle > , Executor >
   s( d , z , o );
  s.set_executor( pool );
  check_path( d , s.solve() , "HDGBFS, " + std::to_string( k ) +
                              " workers on " + ex + ", Korf #55" );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_pgbfs< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_pgbfs< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_pgbfs< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File pgreedy.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
