/*--------------------------------------------------------------------------*/
/*------------------------- File bidirectional.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the bidirectional
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
// a BidirectionalDomain that records the states expanded in each
// direction, i.e., those whose successors or predecessors are asked

template< class D >
struct Spy {

 using state_type = state_t< D >;
 using cost_type = cost_t< D >;

 const D & d;
 mutable std::vector< state_type > seen[ 2 ];

 state_type initial_state( void ) const { return( d.initial_state() ); }
 state_type goal_state( void ) const { return( d.goal_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 cost_type heuristic( const state_type & s ) const {
  return( d.heuristic( s ) );
  }
 cost_type backward_heuristic( const state_type & s ) const {
  return( d.backward_heuristic( s ) );
  }
 cost_type min_edge_cost( void ) const { return( d.min_edge_cost() ); }

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f )
  const {
  seen[ 0 ].push_back( s );
  d.for_each_successor( s , h , std::forward< F >( f ) );
  }

 template< class F >
 void for_each_predecessor( const state_type & s , cost_type h , F && f )
  const {
  seen[ 1 ].push_back( s );
  d.for_each_predecessor( s , h , std::forward< F >( f ) );
  }

 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the distances from s of every state of the 8-puzzle, by a breadth-first
// search (the moves are reversible)

static std::unordered_map< std::uint64_t , int > eight_from(
 std::uint64_t s )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 std::unordered_map< std::uint64_t , int > dist;
 std::queue< std::uint64_t > q;
 dist[ s ] = 0;
 q.push( s );
 while( ! q.empty() ) {
  const auto u = q.front();
  q.pop();
  any.for_each_successor( u , 0 , [ & ]( std::uint64_t t , int , int ) {
   if( dist.emplace( t , dist[ u ] + 1 ).second )
    q.push( t );
   } );
  }
 return( dist );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// every in-memory bidirectional search on d, against the cost expected
// (unsolvable if negative), with the options o

template< class D >
static void check_bidirectional( const D & d , cost_t< D > expected ,
                                 const BidirectionalOptions & o ,
                                 const std::string & tag )
{
 const auto one = [ & ]( const std::string & name , const auto & s ) {
  if( expected < 0 )
   check( s.status == Status::unsolvable , name + " unsolvable, " + tag );
  else
   check_solution( d , s , expected , name + ", " + tag );
  };
 one( "bidirectional Dijkstra" , bidijkstra( d , o ) );
 one( "MM" , mm( d , o ) );
 one( "MMe" , mme( d , o ) );
 one( "fMM(0.5)" , fmm( d , o ) );
 BidirectionalOptions q = o;
 q.fraction = 0.25;
 one( "fMM(0.25)" , fmm( d , q ) );
 if( o.consistent )
  one( "BAE*" , bae( d , o ) );
 one( "NBS" , nbs( d , o ) );
 one( "DVCBS" , dvcbs( d , o ) );
 if constexpr( std::is_integral_v< cost_t< D > > ) {
  one( "GBFHS" , gbfhs( d , o ) );
  one( "bidirectional BFS" , bibfs( d , o ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the variants of the options of the bidirectional searches: the bounds,
// the eager NBS and DVCBS, no epsilon, the choice of the direction

static std::vector< std::pair< BidirectionalOptions , std::string > >
 bidirectional_variants( void )
{
 std::vector< std::pair< BidirectionalOptions , std::string > > v;
 for( int k = 0 ; k < 6 ; ++k ) {
  BidirectionalOptions o;
  std::string n = "default";
  switch( k ) {
   case 1: o.all_bounds = true; n = "all bounds"; break;
   case 2: o.eager = true; n = "eager"; break;
   case 3: o.epsilon = 0; n = "no epsilon"; break;
   case 4: o.side = SideChoice::cardinality; n = "cardinality"; break;
   case 5: o.side = SideChoice::alternate; o.fraction = 0.75;
           n = "alternate"; break;
   default: break;
   }
  v.emplace_back( o , n );
  }
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the sequential bidirectional searches: the predecessors of the domains
// against their successors, every search on states of the 8-puzzle, on
// the 15-puzzle, on random grids with 4 and 8 moves and weighted cells, on
// Romania and on directed graphs, the meeting in the middle of MM and
// fMM, no expansion above C* for MM, NBS, DVCBS and GBFHS, the limit
// cases and the reverse A*

void test_bidirectional( const std::string & data )
{
 static_assert( BidirectionalDomain< FifteenPuzzle > &&
                BackwardOperatorDomain< FifteenPuzzle > &&
                BidirectionalDomain< Grid< 8 > > &&
                BidirectionalDomain< ExplicitGraph< double > > &&
                MinCostDomain< Grid< 4 > > &&
                BidirectionalDomain< ReverseDomain< FifteenPuzzle > > &&
                OperatorDomain< ReverseDomain< FifteenPuzzle > > );

 // the edges seen backward are those seen forward, with their costs
 bool same = true;
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  auto m = GridMap::random( 11 , 7 , 0.25 , seed , 5 );
  m.set_cost( 0 , 0 , 1 );
  m.set_cost( 10 , 6 , 1 );
  const Grid< 8 > g( m , 0 , 0 , 10 , 6 );
  for( std::uint32_t s = 0 ; s < 77 ; ++s ) {
   if( ! m.cost( int( s % 11 ) , int( s / 11 ) ) )
    continue;
   g.for_each_successor( s , 0.0 , [ & ]( std::uint32_t t , double c ,
                                          double ) {
    bool f = false;
    g.for_each_predecessor( t , 0.0 , [ & ]( std::uint32_t u , double b ,
                                             double ) {
     f = f || ( ( u == s ) && ( b == c ) );
     } );
    same = same && f;
    } );
   }
  }
 check( same , "the predecessors of the grids are the reversed moves" );

 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 const auto variants = bidirectional_variants();
 std::size_t k = 0;
 int wrong = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 9071 ) != 0 )                    // 20 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  for( const auto & [ o , n ] : variants ) {
   const int before = failures;
   check_bidirectional( d , dst , o , "8-puzzle, " + n );
   wrong += failures - before;
   }
  }
 check( wrong == 0 , "every bidirectional search on 20 states of the "
                     "8-puzzle and 6 variants" );

 // MM and fMM( p ) meet where they should, and no search expands a node
 // of f above C* (with the true distances, which g bounds)
 bool middle = true , below = true;
 k = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 18143 ) != 0 )                   // 10 states
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto from = eight_from( s );
  const auto to = eight_from( EightPuzzle::pack( EightPuzzle::korf_goal() ) );
  for( double p : { 0.5 , 0.25 , 0.8 } ) {
   const Spy< EightPuzzle > spy{ d , {} };
   BidirectionalOptions o;
   o.fraction = p;
   const auto r = p == 0.5 ? mm( spy , o ) : fmm( spy , o );
   for( auto x : spy.seen[ 0 ] )
    middle = middle && ( from.at( x ) <= p * r.cost + 1e-9 );
   for( auto x : spy.seen[ 1 ] )
    middle = middle && ( to.at( x ) <= ( 1 - p ) * r.cost + 1e-9 );
   }
  for( int a = 0 ; a < 4 ; ++a ) {
   const Spy< EightPuzzle > spy{ d , {} };
   const auto r = a == 0 ? mm( spy ) : a == 1 ? nbs( spy ) :
                  a == 2 ? dvcbs( spy ) : gbfhs( spy );
   for( auto x : spy.seen[ 0 ] )
    below = below && ( from.at( x ) + d.heuristic( x ) <= r.cost );
   for( auto x : spy.seen[ 1 ] )
    below = below && ( to.at( x ) + d.backward_heuristic( x ) <= r.cost );
   }
  }
 check( middle , "MM and fMM(p) expand no node beyond p C*" );
 check( below , "MM, NBS, DVCBS and GBFHS expand no node of f above C*" );

 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 2 ] } ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  check_bidirectional( d , ki.cost , BidirectionalOptions() ,
                       "Korf #" + std::to_string( id ) );
  const auto r = astar( ReverseDomain< FifteenPuzzle >( d ) );
  check( r.solved() && ( r.cost == ki.cost ) &&
         ( r.path.front() == d.goal_state() ) &&
         ( r.path.back() == d.initial_state() ) ,
         "reverse A*, Korf #" + std::to_string( id ) );
  }

 for( std::uint64_t seed = 1 ; seed <= 6 ; ++seed ) {
  const auto m = GridMap::random( 15 , 10 , 0.3 , seed , seed % 2 ? 1 : 4 );
  for( int t = 0 ; t < 3 ; ++t ) {
   const int sx = int( ( seed * 7 + t * 3 ) % 15 ) ,
             sy = int( ( seed + t ) % 10 ) ,
             gx = int( ( seed * 3 + 11 ) % 15 ) ,
             gy = int( ( seed * 5 + t * 7 ) % 10 );
   if( ! m.cost( sx , sy ) || ! m.cost( gx , gy ) )
    continue;
   const std::string tag = "grid, seed " + std::to_string( seed ) +
                           ", pair " + std::to_string( t );
   const Grid< 4 > g4( m , sx , sy , gx , gy );
   const auto a4 = astar( g4 );
   check_bidirectional( g4 , a4.solved() ? a4.cost : -1 ,
                        BidirectionalOptions() , "4 moves, " + tag );
   const Grid< 8 > g8( m , sx , sy , gx , gy );
   const auto a8 = astar( g8 );
   for( const auto & [ o , n ] : variants )
    check_bidirectional( g8 , a8.solved() ? a8.cost : -1.0 , o ,
                         "8 moves, " + n + ", " + tag );
   }
  }

 auto rom = romania( Arad );
 check_bidirectional( rom , 418 , BidirectionalOptions() ,
                      "Romania, Arad to Bucharest" );

 // directed graphs with arcs of cost 0 to 9, without and with a heuristic
 // admissible but not consistent in the variant that allows it
 std::mt19937_64 gen( 3 );
 for( int t = 0 ; t < 30 ; ++t ) {
  const std::uint32_t n = 2 + std::uint32_t( random_below( gen , 20 ) );
  ExplicitGraph< int > g( n , std::uint32_t( random_below( gen , n ) ) ,
                          std::uint32_t( random_below( gen , n ) ) );
  const int arcs = int( random_below( gen , 3 * n ) );
  for( int i = 0 ; i < arcs ; ++i )
   g.add_edge( std::uint32_t( random_below( gen , n ) ) ,
               std::uint32_t( random_below( gen , n ) ) ,
               int( random_below( gen , 10 ) ) , random_below( gen , 3 ) > 0 );
  const auto a = astar( g );
  const std::string tag = "graph " + std::to_string( t );
  check_bidirectional( g , a.solved() ? a.cost : -1 , BidirectionalOptions() ,
                       tag );
  }
 ExplicitGraph< int > bad( 5 , 0 , 4 );        // that of test_reopen
 bad.add_edge( 0 , 1 , 1 , true );
 bad.add_edge( 1 , 2 , 1 , true );
 bad.add_edge( 0 , 2 , 3 , true );
 bad.add_edge( 2 , 3 , 1 , true );
 bad.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };
 for( std::uint32_t v = 0 ; v < 5 ; ++v )
  bad.set_heuristic( v , h[ v ] );
 BidirectionalOptions adm;
 adm.consistent = false;
 check_bidirectional( bad , 4 , adm , "an admissible inconsistent heuristic" );
 bool threw = false;
 try {
  ( void ) bae( bad , adm );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "BAE* refuses heuristics not said consistent" );

 // the limit cases: start and goal the same, one move, the node limit
 const EightPuzzle g0( EightPuzzle::korf_goal() );
 check_bidirectional( g0 , 0 , BidirectionalOptions() , "start = goal" );
 const EightPuzzle g1( walk8( { 1 } ) );
 check_bidirectional( g1 , 1 , BidirectionalOptions() , "one move" );
 BidirectionalOptions lim;
 lim.node_limit = 50;
 const FifteenPuzzle hard( inst[ 0 ].board );
 check( ( mm( hard , lim ).status == Status::node_limit ) &&
        ( bae( hard , lim ).status == Status::node_limit ) &&
        ( nbs( hard , lim ).status == Status::node_limit ) &&
        ( dvcbs( hard , lim ).status == Status::node_limit ) &&
        ( gbfhs( hard , lim ).status == Status::node_limit ) ,
        "the bidirectional searches stop at the node limit" );
 ExplicitGraph< int > cut( 4 , 0 , 3 );        // the goal not reachable
 cut.add_edge( 0 , 1 , 1 , true );
 cut.add_edge( 1 , 2 , 1 , true );
 cut.add_edge( 3 , 2 , 1 , true );
 check_bidirectional( cut , -1 , BidirectionalOptions() ,
                      "the goal not reachable" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the two directions in parallel on the Executor: PNBA* and the parallel
// Dijkstra's algorithm, MM and BAE* on the 8-puzzle, the 15-puzzle, the
// grids and the directed graphs, with an Executor of its own and given

template< class Executor >
void test_bidirectional_parallel( const std::string & data )
{
 const std::string rt = executor_name< Executor >();
 const auto all = [ & ]( const auto & d , auto expected ,
                         const std::string & tag , Executor * ex ) {
  using DD = std::decay_t< decltype( d ) >;
  const auto run = [ & ]( auto search , const std::string & name ) {
   if( ex )
    search.set_executor( *ex );
   const auto s = search.solve();
   if( expected < 0 )
    check( s.status == Status::unsolvable , name + " unsolvable, " + tag );
   else
    check_solution( d , s , expected , name + ", " + tag );
   };
  run( ParallelBidirectional< DD , ParallelBiKind::nba , Executor >( d ) ,
       "PNBA* " + rt );
  run( ParallelBidirectional< DD , ParallelBiKind::dijkstra , Executor >( d ) ,
       "parallel bidirectional Dijkstra " + rt );
  run( ParallelBidirectional< DD , ParallelBiKind::mm , Executor >( d ) ,
       "parallel MM " + rt );
  run( ParallelBidirectional< DD , ParallelBiKind::bae , Executor >( d ) ,
       "parallel BAE* " + rt );
  };
 Executor ex( 2 );
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 std::size_t k = 0;
 int wrong = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 3628 ) != 0 )                    // 50 states, spread out
   continue;
  const int before = failures;
  all( EightPuzzle( EightPuzzle::unpack( s ) ) , dst , "8-puzzle" ,
       k % 2 ? &ex : nullptr );
  wrong += failures - before;
  }
 check( wrong == 0 , "the parallel bidirectional searches on 50 states of "
                     "the 8-puzzle, " + rt );
 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 3 ] } )
  all( FifteenPuzzle( inst[ id - 1 ].board ) , inst[ id - 1 ].cost ,
       "Korf #" + std::to_string( id ) , &ex );
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const auto m = GridMap::random( 20 , 12 , 0.3 , seed , 3 );
  if( ! m.cost( 0 , 0 ) || ! m.cost( 19 , 11 ) )
   continue;
  const Grid< 8 > g( m , 0 , 0 , 19 , 11 );
  const auto a = astar( g );
  all( g , a.solved() ? a.cost : -1.0 , "grid, seed " +
       std::to_string( seed ) , &ex );
  }
 std::mt19937_64 gen( 5 );
 for( int t = 0 ; t < 20 ; ++t ) {
  const std::uint32_t n = 2 + std::uint32_t( random_below( gen , 20 ) );
  ExplicitGraph< int > g( n , std::uint32_t( random_below( gen , n ) ) ,
                          std::uint32_t( random_below( gen , n ) ) );
  const int arcs = int( random_below( gen , 3 * n ) );
  for( int i = 0 ; i < arcs ; ++i )
   g.add_edge( std::uint32_t( random_below( gen , n ) ) ,
               std::uint32_t( random_below( gen , n ) ) ,
               1 + int( random_below( gen , 9 ) ) ,
               random_below( gen , 3 ) > 0 );
  const auto a = astar( g );
  all( g , a.solved() ? a.cost : -1 , "graph " + std::to_string( t ) , &ex );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel external-memory searches with 1, 2 and 4 workers: PEM-A*,
// PEM-rA*, PEMM and PEM-BAE* on the 8-puzzle, the grids and Romania, and on
// random walks of the 15-puzzle with 2 workers, the directory removed at
// the end, the cap on the disk, the node limit and the inconsistent
// heuristic refused

template< class Executor >
void test_pem( const std::string & data )
{
 const std::string rt = executor_name< Executor >();
 ADS::TempDirectory dir;
 PEMOptions po;
 po.directory = dir.path();
 po.partitions = 5;
 po.buffer_bytes = 1024;
 const auto all = [ & ]( const auto & d , auto expected ,
                         const std::string & tag , unsigned p ) {
  using DD = std::decay_t< decltype( d ) >;
  const auto one = [ & ]( const std::string & name , const auto & s ) {
   if( expected < 0 )
    check( s.status == Status::unsolvable , name + " unsolvable, " + tag );
   else
    check_solution( d , s , expected , name + ", " + tag );
   };
  const std::string w = " " + rt + " p " + std::to_string( p );
  one( "PEM-A*" + w , PEMSearch< DD , PEMKind::astar , Executor >(
                        d , p , {} , po ).solve() );
  one( "PEMM" + w , PEMSearch< DD , PEMKind::mm , Executor >(
                      d , p , {} , po ).solve() );
  one( "PEM-BAE*" + w , PEMSearch< DD , PEMKind::bae , Executor >(
                          d , p , {} , po ).solve() );
  auto r = PEMSearch< ReverseDomain< DD > , PEMKind::astar , Executor >(
            ReverseDomain< DD >( d ) , p , {} , po ).solve();
  std::reverse( r.path.begin() , r.path.end() );
  one( "PEM-rA*" + w , r );
  };
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 for( unsigned p : { 1u , 2u , 4u } ) {
  std::size_t k = 0;
  int wrong = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 30241 ) != 0 )                  // 6 states
    continue;
   const int before = failures;
   all( EightPuzzle( EightPuzzle::unpack( s ) ) , dst , "8-puzzle" , p );
   wrong += failures - before;
   }
  check( wrong == 0 , "the PEM searches on 6 states of the 8-puzzle, " +
                      rt + " p " + std::to_string( p ) );
  for( std::uint64_t seed = 1 ; seed <= 2 ; ++seed ) {
   const auto m = GridMap::random( 14 , 9 , 0.25 , seed , 3 );
   if( ! m.cost( 0 , 8 ) || ! m.cost( 13 , 0 ) )
    continue;
   const Grid< 4 > g( m , 0 , 8 , 13 , 0 );
   const auto a = astar( g );
   all( g , a.solved() ? a.cost : -1 , "grid, seed " +
        std::to_string( seed ) , p );
   }
  all( romania( Arad ) , 418 , "Romania" , p );
  }
 for( std::uint64_t seed = 1 ; seed <= 2 ; ++seed ) {
  const FifteenPuzzle w( random_walk< FifteenPuzzle >( 60 , seed ) );
  all( w , astar( w ).cost , "15-puzzle, walk " + std::to_string( seed ) ,
       2 );
  }
 const auto inst = korf( data );
 const FifteenPuzzle d( inst[ easy_korf[ 0 ] - 1 ].board );
 check( std::filesystem::is_empty( dir.path() ) ,
        "the PEM searches remove their files, " + rt );

 PEMOptions cap = po;
 cap.disk_limit = 4096;
 check( PEMSearch< FifteenPuzzle , PEMKind::bae , Executor >(
         d , 2 , {} , cap ).solve().status == Status::disk_limit ,
        "PEM-BAE* stops at the cap on the disk, " + rt );
 SearchOptions lim;
 lim.node_limit = 1000;
 check( PEMSearch< FifteenPuzzle , PEMKind::mm , Executor >(
         d , 2 , lim , po ).solve().status == Status::node_limit ,
        "PEMM stops at the node limit, " + rt );
 ExplicitGraph< int > bad( 5 , 0 , 4 );
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
  ( void ) PEMSearch< ExplicitGraph< int > , PEMKind::astar , Executor >(
            bad , 2 , {} , po ).solve();
  }
 catch( const std::logic_error & ) {
  threw = true;
  }
 check( threw , "PEM-A* refuses an inconsistent heuristic, " + rt );
 check( std::filesystem::is_empty( dir.path() ) ,
        "the PEM searches leave no file when they stop or throw, " + rt );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A*-IDD against the exact distances of the 8-puzzle (with blocks small
// enough that most of the closed list is on the disk), on the 15-puzzle,
// the grids and Romania, the reopening of an inconsistent heuristic, the
// limits on the table and on the disk, and the files removed

void test_idd( const std::string & data )
{
 ADS::TempDirectory dir;
 IDDOptions o;
 o.directory = dir.path();
 o.slots = 1 << 19;                        // 458752 closed states
 o.block = 16;
 o.partitions = 7;
 o.open_block = 8;
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 std::size_t k = 0;
 int wrong = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 3628 ) != 0 )                    // 50 states
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto r = idd_astar( d , {} , o );
  if( ! ( r.solved() && ( r.cost == dst ) ) ) {
   ++wrong;
   check_solution( d , r , dst , "A*-IDD 8-puzzle" );
   }
  }
 check( wrong == 0 , "A*-IDD on 50 states of the 8-puzzle" );
 const auto inst = korf( data );
 const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
 const FifteenPuzzle d( ki.board );
 const auto r = idd_astar( d , {} , o );
 check_solution( d , r , ki.cost , "A*-IDD, Korf #55" );
 check( ( r.stats.expanded == r.stats.stored ) && ( r.stats.messages > 0 ) ,
        "A*-IDD expands each state once and reads the closed file" );
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 9 , 7 , 0.2 , seed , 3 );
  const Grid< 4 > g( m , 0 , 6 , 8 , 6 );
  const auto a = astar( g );
  const auto x = idd_astar( g , {} , o );
  if( a.solved() )
   check_solution( g , x , a.cost , "A*-IDD grid, seed " +
                   std::to_string( seed ) );
  else
   check( x.status == Status::unsolvable , "A*-IDD unsolvable grid, seed " +
          std::to_string( seed ) );
  }
 check_solution( romania( Arad ) , idd_astar( romania( Arad ) , {} , o ) ,
                 418 , "A*-IDD Romania" );
 ExplicitGraph< int > bad( 5 , 0 , 4 );
 bad.add_edge( 0 , 1 , 1 , true );
 bad.add_edge( 1 , 2 , 1 , true );
 bad.add_edge( 0 , 2 , 3 , true );
 bad.add_edge( 2 , 3 , 1 , true );
 bad.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };
 for( std::uint32_t v = 0 ; v < 5 ; ++v )
  bad.set_heuristic( v , h[ v ] );
 const auto b = idd_astar( bad , {} , o );
 check_solution( bad , b , 4 , "A*-IDD, inconsistent heuristic" );
 check( b.stats.reexpanded > 0 , "A*-IDD reopens a closed state" );
 IDDOptions small = o;
 small.slots = 1024;
 small.block = 64;
 const FifteenPuzzle hard( inst[ 0 ].board );
 check( idd_astar( hard , {} , small ).status == Status::node_limit ,
        "A*-IDD stops when its table is full" );
 IDDOptions cap = o;
 cap.disk_limit = 4096;
 check( idd_astar( hard , {} , cap ).status == Status::disk_limit ,
        "A*-IDD stops at the cap on the disk" );
 check( std::filesystem::is_empty( dir.path() ) ,
        "A*-IDD removes its files" );
 bool threw = false;
 try {
  IDDOptions none = o;
  none.slots = 8;
  ( void ) idd_astar( d , {} , none );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "A*-IDD refuses a table smaller than 8 blocks" );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_bidirectional_parallel< ADS::ThreadPool >
  ( const std::string & );
template void test_pem< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_bidirectional_parallel< ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_pem< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_bidirectional_parallel< ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_pem< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File bidirectional.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
