/*--------------------------------------------------------------------------*/
/*---------------------------- File bounded.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the
 * bounded-suboptimal and anytime searches and HDA* over MPI, the planning
 * tasks among its instances [see incremental.cpp for LPA*].
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

void test_weighted( const std::string & data )
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

void test_anytime( const std::string & data )
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

void test_focal( const std::string & data )
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
// the parallel searches with a weight, duplicate dropping and anytime, on
// the given runtime: within the bound, and optimal when anytime; SPA*,
// MQA* and SafePBNF only if Shared, as in their own cases; in a light run
// [see light_run()] with 1 and 3 workers, on the walk of fifteen_cases()
// and on a grid of 100 by 60 cells

template< class Executor , bool Shared>
void test_bounded( const std::string & data )
{
 const auto fifteen = fifteen_cases( korf( data ) , 2 );
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 const int gw = light_run() ? 100 : 200 , gh = light_run() ? 60 : 120;
 const auto m = GridMap::random( gw , gh , 0.25 , 4 , 5 );
 const double r8 = grid_dijkstra( m , 8 , 0 , gh - 1 , gw - 1 , gh - 1 );
 const Grid< 8 > g8( m , 0 , gh - 1 , gw - 1 , gh - 1 , 20 , 10 );
 auto rg = romania( Arad );
 const auto mode = []( int k ) {
  return( std::string( k == 0 ? "" : ( k == 1 ? "dd " : "anytime " ) ) );
  };
 for( unsigned p : case_workers( { 1 , 3 , 8 } , { 1 , 3 } ) ) {
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
    for( const auto & ki : fifteen ) {
     const FifteenPuzzle d( ki.board );
     const std::string t = tag + ki.name;
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
      const auto & ki = fifteen.front();
      const FifteenPuzzle d( ki.board );
      check_bounded( d , SafePBNF< FifteenPuzzle , Executor >( d , a )
                          .solve() , ki.cost , 1 , "ARPBNF" + tag +
                          ki.name );
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

void test_mpi( const std::string & data )
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

 const auto pd = planning( data , "blocks-6-2.sas" ,
                           PlanningHeuristic::mas , 512 );
 const ZobristDistribution< PlanningD > pz( pd );
 check_plan( pd , HDAStarMPI< PlanningD , decltype( pz ) >( pd , pz ).solve() ,
             20 , "MPI HDA*[Z] on blocks-6-2" );
 const AbstractZobristDistribution< PlanningD > paz( pd );
 check_plan( pd , HDAStarMPI< PlanningD , decltype( paz ) >( pd , paz )
                   .solve() , 20 , "MPI HDA*[AZ] on blocks-6-2" );

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

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_bounded< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_bounded< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_bounded< ADS::FastFlowExecutor , false >
   ( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File bounded.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
