/*--------------------------------------------------------------------------*/
/*--------------------------- File parallel.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the parallel
 * searches: HDA*, SPA*, MQA*, PBNF, PA*SE, the searches with expensive
 * edges and the farm.
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
void test_hda( const std::string & data )
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
void test_spa( const std::string & data )
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
void test_mqa( const std::string & data )
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
void test_pbnf( const std::string & data )
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
void test_pase( const std::string & data )
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

  // wPA*SE: the weight in the keys, eps raised to it, and the scan of all
  // of OPEN; with integer costs a weight of 1.3 is taken as 333 / 256
  for( double w : { 1.3 , 2.0 , 3.0 } )
   for( std::size_t k : { std::size_t( 64 ) , std::size_t( 0 ) } ) {
    PASEOptions ow;
    ow.threads = p;
    ow.weight = w;
    ow.candidates = k;
    const auto sw = PASE< FifteenPuzzle , Executor >( d , ow ).solve();
    check( sw.solved() && ( sw.cost >= ki.cost ) &&
           ( sw.cost <= w * ki.cost ) && ( sw.stats.reexpanded == 0 ) ,
           "wPA*SE w " + std::to_string( w ) + " scan " +
           std::to_string( k ) + " within the bound, p " +
           std::to_string( p ) );
    if( r8.solved() ) {
     const auto gw = PASE< Grid< 8 > , Executor >( g8 , ow ).solve();
     check( gw.solved() && ( gw.cost <= w * r8.cost + 1e-9 ) ,
            "wPA*SE grid 8 w " + std::to_string( w ) + ", p " +
            std::to_string( p ) );
     }
    }
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
// the 2D navigation of the edge-based searches on a map of the data, its
// k-th pair of points, with the expensive diagonal motions if repeat > 1

static std::optional< std::pair< NavMap , std::array< int , 4 > > > nav_map(
                       const std::string & data , const std::string & name ,
                       int scale , std::size_t k )
{
 const std::string dir = data + "/robot-nav-2d/" + name + "/";
 std::ifstream in( dir + name + ".map" ) , fs( dir + "nav2d_starts.txt" ) ,
               fg( dir + "nav2d_goals.txt" );
 if( ! in || ! fs || ! fg )
  return( std::nullopt );
 auto m = NavMap::read_movingai( in , scale );
 const auto pts = RobotNav2D::read_points( fs , fg , scale );
 if( k >= pts.size() )
  return( std::nullopt );
 return( std::pair( std::move( m ) , pts[ k ] ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a map of 240 by 160 cells with two walls, the gap of the first above the
// line from ( 30 , 80 ) to ( 210 , 80 ) and that of the second below it,
// for the navigation of the parallel tests, small enough for the sanitizers

static NavMap nav_walls( void )
{
 NavMap m( 240 , 160 );
 for( int y = 0 ; y < 160 ; ++y ) {
  if( ( y < 40 ) || ( y > 90 ) )
   for( int x = 76 ; x < 84 ; ++x )
    m.set_blocked( x , y , true );
  if( ( y < 96 ) || ( y > 150 ) )
   for( int x = 156 ; x < 164 ; ++x )
    m.set_blocked( x , y , true );
  }
 return( m );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the edge domains and the sequential edge searches: the edges of
// SlowEdgeDomain and RobotNav2D against their expansions, eA*, PwA* and the
// lazy search with one worker against A*, on several domains

void test_edge( const std::string & data )
{
 // the actions and the evaluated edges give the successors of an expansion
 const auto same_edges = []( const auto & d , const auto & s ,
                             const std::string & what ) {
  using DD = std::remove_cvref_t< decltype( d ) >;
  std::vector< std::pair< state_t< DD > , cost_t< DD > > > a , b;
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const auto & t , auto c , auto ) {
   a.emplace_back( t , c );
   } );
  d.for_each_action( s , [ & ]( std::size_t x ) {
   d.evaluate_edge( s , d.heuristic( s ) , x ,
                    [ & ]( const auto & t , auto c , auto ) {
    b.emplace_back( t , c );
    } );
   } );
  check( a == b , what + ": the edges are the successors" );
  };
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 EdgeWork ew;
 ew.cheap = 10;
 ew.expensive = 100;
 ew.expensive_share = 0.5;
 ew.missing_share = 0.15;
 std::size_t k = 0 , bad = 0 , n = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 9007 ) != 0 )                      // 21 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const SlowEdgeDomain< EightPuzzle > e( d , ew );
  same_edges( e , e.initial_state() , "SlowEdgeDomain 8-puzzle" );
  const auto r = astar( e );                     // the graph with missing
  for( double w : { 1.0 , 2.0 } ) {
   PWAStarOptions po;
   po.weight = w;
   EPASEOptions eo;
   eo.weight = w;
   MPLPOptions mo;
   mo.weight = w;
   const auto a = PWAStar< decltype( e ) >( e , po ).solve();
   const auto b = EPASE< decltype( e ) >( e , eo ).solve();
   eo.generalized = true;
   const auto c = EPASE< decltype( e ) >( e , eo ).solve();
   const auto m = MPLP< decltype( e ) >( e , mo ).solve();
   for( const auto * x : { & a , & b , & c , & m } ) {
    ++n;
    if( x->solved() != r.solved() )
     ++bad;
    else
     if( r.solved() && ( ( x->cost < r.cost ) || ( x->cost > w * r.cost ) ) )
      ++bad;
     else
      if( r.solved() )
       check_solution( e , *x , x->cost , "edge search 8-puzzle path" );
    }
   if( ( w == 1 ) && r.solved() )
    check( ( a.stats.evaluated <= 4 * a.stats.expanded ) &&
           ( b.stats.reexpanded == 0 ) && ( c.stats.reexpanded == 0 ) ,
           "the edges of the 8-puzzle evaluated once" );
   }
  }
 check( bad == 0 , "PwA*, eA*, GeA* and the lazy search on 21 states of the "
                   "8-puzzle with missing edges, within the bound (" +
                   std::to_string( n ) + " runs)" );

 // the lazy edges of SlowEdgeDomain are its edges, missing ones included
 {
  const EightPuzzle d( walk8( { 1 , 2 , 5 , 4 , 3 } ) );
  const SlowEdgeDomain< EightPuzzle > e( d , ew );
  std::size_t lazy = 0 , real = 0;
  e.for_each_action( e.initial_state() , [ & ]( std::size_t x ) {
   e.lazy_edge( e.initial_state() , 0 , x , [ & ]( auto , int , int ) {
    ++lazy;
    } );
   e.evaluate_edge( e.initial_state() , 0 , x , [ & ]( auto , int , int ) {
    ++real;
    } );
   } );
  check( ( lazy >= real ) && ( lazy > 0 ) , "the lazy edges include the "
                                            "evaluated ones" );
  }

 // a random grid with 8 moves and missing edges: real costs
 const auto gm = GridMap::random( 60 , 40 , 0.2 , 5 );
 const Grid< 8 > g8( gm , 0 , 39 , 59 , 39 );
 const SlowEdgeDomain< Grid< 8 > > ge( g8 , ew );
 same_edges( ge , ge.initial_state() , "SlowEdgeDomain grid 8" );
 const auto rg = astar( ge );
 if( rg.solved() ) {
  check_solution( ge , PWAStar< decltype( ge ) >( ge ).solve() , rg.cost ,
                  "PwA* grid 8 with missing edges" );
  check_solution( ge , EPASE< decltype( ge ) >( ge ).solve() , rg.cost ,
                  "eA* grid 8 with missing edges" );
  check_solution( ge , MPLP< decltype( ge ) >( ge ).solve() , rg.cost ,
                  "lazy search grid 8 with missing edges" );
  }

 // the navigation: a map with a wall and a gap, and a map of the data
 NavMap wall( 200 , 120 );
 for( int y = 0 ; y < 120 ; ++y )
  if( ( y < 40 ) || ( y > 80 ) )
   for( int x = 96 ; x < 104 ; ++x )
    wall.set_blocked( x , y , true );
 NavOptions no;
 no.repeat = 30;
 const RobotNav2D nw( wall , 30 , 60 , 170 , 60 , no );
 same_edges( nw , nw.initial_state() , "RobotNav2D" );
 check( nw.pairwise_heuristic( RobotNav2D::pack( 0 , 0 ) ,
                               RobotNav2D::pack( 3 , 4 ) ) == 5 ,
        "RobotNav2D pairwise heuristic" );
 const auto rn = astar( nw );
 check( rn.solved() , "RobotNav2D through the gap of the wall" );
 if( rn.solved() ) {
  check_solution( nw , EPASE< RobotNav2D >( nw ).solve() , rn.cost ,
                  "eA* RobotNav2D" );
  check_solution( nw , MPLP< RobotNav2D >( nw ).solve() , rn.cost ,
                  "lazy search RobotNav2D" );
  }
 NavMap closed( 200 , 120 );
 for( int y = 0 ; y < 120 ; ++y )
  for( int x = 96 ; x < 104 ; ++x )
   closed.set_blocked( x , y , true );
 const RobotNav2D nc( closed , 30 , 60 , 170 , 60 , no );
 check( ( EPASE< RobotNav2D >( nc ).solve().status == Status::unsolvable ) &&
        ( MPLP< RobotNav2D >( nc ).solve().status == Status::unsolvable ) &&
        ( PWAStar< RobotNav2D >( nc ).solve().status ==
          Status::unsolvable ) , "the edge searches prove the closed wall "
                                 "unsolvable" );
 if( const auto mp = nav_map( data , "den520d" , 5 , 0 ) ) {
  const auto & [ m , q ] = *mp;
  for( bool factors : { false , true } ) {
   NavOptions o;
   o.factors = factors;
   const RobotNav2D d( m , q[ 0 ] , q[ 1 ] , q[ 2 ] , q[ 3 ] , o );
   const auto r = astar( d );
   check( r.solved() , "RobotNav2D den520d #1 solved" );
   if( ! r.solved() )
    continue;
   const std::string tag = factors ? " with factors" : "";
   check_solution( d , PWAStar< RobotNav2D >( d ).solve() , r.cost ,
                   "PwA* den520d #1" + tag );
   const auto e = EPASE< RobotNav2D >( d ).solve();
   check_solution( d , e , r.cost , "eA* den520d #1" + tag );
   check( e.stats.evaluated <= 8 * r.stats.expanded ,
          "eA* evaluates no more edges than A* den520d #1" + tag );
   PWAStarOptions pw;
   pw.weight = 50;
   EPASEOptions ew50;
   ew50.weight = 50;
   const auto a50 = PWAStar< RobotNav2D >( d , pw ).solve();
   const auto e50 = EPASE< RobotNav2D >( d , ew50 ).solve();
   check( a50.solved() && e50.solved() && ( a50.cost <= 50 * r.cost ) &&
          ( e50.cost <= 50 * r.cost ) , "w 50 den520d #1" + tag );
   }
  }
 else
  check( false , "cannot read the map den520d from " + data );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel edge searches, 1 to 8 workers: ePA*SE, w-ePA*SE and
// GePA*SE with and without the dispatcher, A-ePA*SE and A-GePA*SE, PwA*
// and MPLP, on the 8-puzzle, a grid and the navigation

template< class Executor >
void test_edge_parallel( const std::string & )
{
 const std::string rt = " on " + executor_name< Executor >();
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 EdgeWork ew;
 ew.cheap = 5;
 ew.expensive = 100;
 ew.expensive_share = 0.5;
 ew.missing_share = 0.1;
 const auto gm = GridMap::random( 60 , 40 , 0.2 , 5 );
 const Grid< 8 > g8( gm , 0 , 39 , 59 , 39 );
 const SlowEdgeDomain< Grid< 8 > > ge( g8 , ew );
 const auto rg = astar( ge );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  const std::string tag = ", p " + std::to_string( p ) + rt;
  std::size_t k = 0 , bad = 0 , n = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 30011 ) != 0 )                    // 7 states, spread out
    continue;
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const SlowEdgeDomain< EightPuzzle > e( d , ew );
   const auto r = astar( e );
   for( double w : { 1.0 , 2.5 } )
    for( int v = 0 ; v < 4 ; ++v ) {
     EPASEOptions o;
     o.threads = p;
     o.weight = w;
     o.generalized = v & 1;
     o.dispatcher = ! ( v & 2 );
     const auto x = EPASE< decltype( e ) , Executor >( e , o ).solve();
     ++n;
     if( ( x.solved() != r.solved() ) || ( x.stats.reexpanded != 0 ) ||
         ( r.solved() && ( ( x.cost < r.cost ) ||
                           ( x.cost > w * r.cost ) ) ) )
      ++bad;
     }
   PWAStarOptions po;
   po.threads = p;
   const auto a = PWAStar< decltype( e ) , Executor >( e , po ).solve();
   const auto a1 = PWAStar< decltype( e ) >( e ).solve();
   if( ( a.cost != a1.cost ) || ( a.stats.expanded != a1.stats.expanded ) ||
       ( a.stats.evaluated != a1.stats.evaluated ) )
    ++bad;
   MPLPOptions mo;
   mo.threads = p;
   const auto m = MPLP< decltype( e ) , Executor >( e , mo ).solve();
   if( ( m.solved() != r.solved() ) || ( r.solved() && ( m.cost != r.cost ) ) )
    ++bad;
   n += 2;
   }
  check( bad == 0 , "ePA*SE, w-ePA*SE, GePA*SE, PwA* and MPLP on 7 states "
                    "of the 8-puzzle (" + std::to_string( n ) + " runs)" +
                    tag );

  if( rg.solved() )
   for( int v = 0 ; v < 4 ; ++v ) {
    EPASEOptions o;
    o.threads = p;
    o.generalized = v & 1;
    o.dispatcher = ! ( v & 2 );
    o.candidates = v == 3 ? 0 : 64;
    check_solution( ge , EPASE< decltype( ge ) , Executor >( ge , o )
                         .solve() , rg.cost , "ePA*SE grid 8, variant " +
                                              std::to_string( v ) + tag );
    o.anytime = true;
    o.weight = 3;
    o.weight_step = 1;
    const auto x = EPASE< decltype( ge ) , Executor >( ge , o ).solve();
    check_solution( ge , x , rg.cost , "A-ePA*SE grid 8, variant " +
                                       std::to_string( v ) + tag );
    double last = infinite_cost< double >();
    bool mono = true;
    for( const auto & y : x.improvements ) {
     mono = mono && ( y.cost <= last ) && ( y.bound <= rg.cost + 1e-9 );
     last = y.cost;
     }
    check( mono && ( x.improvements.size() == 3 ) &&
           ( x.stats.iterations == 3 ) &&
           ( x.stats.evaluated <= 8 * x.stats.stored ) ,
           "A-ePA*SE grid 8 improves within the bounds, each edge "
           "evaluated once, variant " + std::to_string( v ) + tag );
    MPLPOptions mo;
    mo.threads = p;
    mo.weight = 2;
    const auto m = MPLP< decltype( ge ) , Executor >( ge , mo ).solve();
    check( m.solved() && ( m.cost <= 2 * rg.cost + 1e-9 ) ,
           "MPLP w 2 grid 8, variant " + std::to_string( v ) + tag );
    }

  {
   const NavMap m = nav_walls();
   NavOptions no;
   no.factors = true;
   no.repeat = 3;
   const RobotNav2D d( m , 30 , 80 , 210 , 80 , no );
   const auto r = astar( d );
   if( r.solved() )
    for( double w : { 1.0 , 50.0 } ) {
     EPASEOptions o;
     o.threads = p;
     o.weight = w;
     const auto e = EPASE< RobotNav2D , Executor >( d , o ).solve();
     o.generalized = true;
     const auto g = EPASE< RobotNav2D , Executor >( d , o ).solve();
     PASEOptions pa;
     pa.threads = p;
     pa.weight = w;
     const auto ps = PASE< RobotNav2D , Executor >( d , pa ).solve();
     MPLPOptions mo;
     mo.threads = p;
     mo.weight = w;
     const auto ml = MPLP< RobotNav2D , Executor >( d , mo ).solve();
     const std::string tw = " w " + std::to_string( int( w ) ) + tag;
     if( w == 1 ) {
      check_solution( d , e , r.cost , "ePA*SE walls" + tw );
      check_solution( d , g , r.cost , "GePA*SE walls" + tw );
      check_solution( d , ps , r.cost , "PA*SE walls" + tw );
      check_solution( d , ml , r.cost , "MPLP walls" + tw );
      }
     else
      check( e.solved() && g.solved() && ps.solved() && ml.solved() &&
             ( e.cost <= w * r.cost ) && ( g.cost <= w * r.cost ) &&
             ( ps.cost <= w * r.cost ) && ( ml.cost <= w * r.cost ) ,
             "ePA*SE, GePA*SE, wPA*SE and MPLP walls" + tw );
     }
   }
  }

 // the edge cases: no worker, the goal at the start, one move, the limit
 const auto edges = []( const EightPuzzle & d ) {
  return( SlowEdgeDomain< EightPuzzle >( d , EdgeWork{} ) );
  };
 const auto epase = [ & ]( const EightPuzzle & d , const EPASEOptions & x ) {
  return( EPASE< SlowEdgeDomain< EightPuzzle > , Executor >( edges( d ) , x )
          .solve() );
  };
 const auto pwa = [ & ]( const EightPuzzle & d , const PWAStarOptions & x ) {
  return( PWAStar< SlowEdgeDomain< EightPuzzle > , Executor >( edges( d ) ,
                                                               x ).solve() );
  };
 const auto mplp = [ & ]( const EightPuzzle & d , const MPLPOptions & x ) {
  return( MPLP< SlowEdgeDomain< EightPuzzle > , Executor >( edges( d ) , x )
          .solve() );
  };
 for( unsigned p : { 0u , 4u } ) {
  EPASEOptions eo;
  eo.threads = p;
  check_astar_edges( eo , epase , "ePA*SE " + std::to_string( p ) +
                                  " workers" + rt );
  eo.dispatcher = false;
  check_astar_edges( eo , epase , "ePA*SE without dispatcher " +
                                  std::to_string( p ) + " workers" + rt );
  eo.anytime = true;
  eo.weight = 2;
  check_astar_edges( eo , epase , "A-ePA*SE " + std::to_string( p ) +
                                  " workers" + rt );
  PWAStarOptions po;
  po.threads = p;
  check_astar_edges( po , pwa , "PwA* " + std::to_string( p ) +
                                " workers" + rt );
  MPLPOptions mo;
  mo.threads = p;
  check_astar_edges( mo , mplp , "MPLP " + std::to_string( p ) +
                                 " workers" + rt );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the master-worker FarmAStar, 1 to 8 workers, 1 or 4 nodes per task

template< class Executor >
void test_farm( const std::string & data )
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
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_hda< ADS::ThreadPool >( const std::string & );
template void test_spa< ADS::ThreadPool >( const std::string & );
template void test_mqa< ADS::ThreadPool >( const std::string & );
template void test_pbnf< ADS::ThreadPool >( const std::string & );
template void test_pase< ADS::ThreadPool >( const std::string & );
template void test_edge_parallel< ADS::ThreadPool >( const std::string & );
template void test_farm< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_hda< ADS::OpenMPExecutor >( const std::string & );
 template void test_spa< ADS::OpenMPExecutor >( const std::string & );
 template void test_mqa< ADS::OpenMPExecutor >( const std::string & );
 template void test_pbnf< ADS::OpenMPExecutor >( const std::string & );
 template void test_pase< ADS::OpenMPExecutor >( const std::string & );
 template void test_edge_parallel< ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_farm< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_hda< ADS::FastFlowExecutor >( const std::string & );
 template void test_pase< ADS::FastFlowExecutor >( const std::string & );
 template void test_edge_parallel< ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_farm< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File parallel.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
