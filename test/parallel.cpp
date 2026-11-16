/*--------------------------------------------------------------------------*/
/*--------------------------- File parallel.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the parallel
 * searches: HDA*, SPA*, MQA*, PBNF and PA*SE.
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
  const FeatureModel< D > m( d , 200 );
  const auto oz = make_distribution( m , "OZ:1" );
  static_assert( IncrementalDistribution< decltype( z ) , D > &&
                 IncrementalDistribution< decltype( az ) , D > &&
                 IncrementalDistribution< decltype( oz ) , D > );
  const auto s = d.initial_state();
  for( unsigned p : { 1u , 2u , 3u , 7u , 64u } ) {
   const auto sz = successor_owners< D >( z , s , z.key( s ) , p );
   const auto saz = successor_owners< D >( az , s , az.key( s ) , p );
   const auto soz = successor_owners< D >( oz , s , oz.key( s ) , p );
   d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
    if( ( sz( t ) != std::pair( z.owner( t , p ) , z.key( t ) ) ) ||
        ( saz( t ) != std::pair( az.owner( t , p ) , az.key( t ) ) ) ||
        ( soz( t ) != std::pair( oz.owner( t , p ) , oz.key( t ) ) ) )
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// HDA* and SPA* with one worker on the random maps of 60 x 60 cells with 8
// moves, whose real costs reached in two orders differ in the last bits:
// no node is reopened by the rounding, and the cost is that of A*

template< class Executor >
static void check_rounding( void )
{
 for( std::uint64_t k = 1 ; k <= 3 ; ++k ) {
  const auto m = GridMap::random( 60 , 60 , 0.25 , k );
  const Grid< 8 > g( m , 0 , 59 , 59 , 59 );
  const auto a = astar( g );
  if( ! a.solved() )
   continue;
  const std::string n = " on the random map " + std::to_string( k ) +
                        " of 60 x 60 cells, 8 moves";
  ParallelOptions o;
  o.threads = 1;
  const AbstractZobristDistribution< Grid< 8 > > z( g );
  const auto h = HDAStar< Grid< 8 > , decltype( z ) , Executor >( g , z , o )
                  .solve();
  const auto s = SPAStar< Grid< 8 > , Executor >( g , o ).solve();
  check_solution( g , h , a.cost , "HDA* with one worker" + n );
  check_solution( g , s , a.cost , "SPA* with one worker" + n );
  check( ( a.stats.reexpanded == 0 ) && ( h.stats.reexpanded == 0 ) &&
         ( s.stats.reexpanded == 0 ) ,
         "A*, HDA* and SPA* reopen no node, the rounding of g aside" + n );
  }
 }

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
   if( id != easy_korf[ 0 ] )
    continue;
   // the synchronous sends, the locked inboxes (PRA*, APRA*) and the
   // lending of the nodes
   using Open = BucketOpenList< Node< FifteenPuzzle::state_type , int > >;
   using Locked = LockedChannelsOf< Executor >;
   HDAOptions so( o );
   so.sync_send = true;
   const auto ss = HDAStar< FifteenPuzzle , decltype( z ) , Executor >(
                    d , z , so ).solve();
   check_solution( d , ss , ki.cost , "HDA*[Z] synchronous sends" + tag );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( z ) , Executor ,
                                Open , Locked >( d , z , o ).solve() ,
                   ki.cost , "HDA*[Z] locked inboxes" + tag );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( z ) , Executor ,
                                Open , Locked >( d , z , so ).solve() ,
                   ki.cost , "PRA*" + tag );
   check_solution( d , HDAStar< FifteenPuzzle , decltype( a ) , Executor ,
                                Open , Locked >( d , a , so ).solve() ,
                   ki.cost , "APRA*" + tag );
   HDAOptions lo( o );
   lo.outsource = 2;
   check_solution( d , HDAStar< FifteenPuzzle , decltype( z ) , Executor >(
                        d , z , lo ).solve() , ki.cost ,
                   "HDA*[Z] lending at a gap of 2" + tag );
   // the distributions built on a sample, and two hand-made ones
   const FeatureModel< FifteenPuzzle > m( d , 4096 );
   std::vector< TableZobristDistribution< FifteenPuzzle > > tz;
   for( const auto & n : auto_distribution_names() )
    tz.push_back( make_distribution( m , n ) );
   tz.push_back( auto_distribution( m , p ) );
   const auto blocks = FifteenPuzzle::tile_abstraction( TileRegions::block ,
                                                        { 1 , 2 , 3 , 4 } );
   tz.emplace_back( d , FeatureAbstraction::of( blocks ) , 1 ,
                    "tile:block:1-4" );
   for( const auto & x : tz )
    check_solution( d , HDAStar< FifteenPuzzle , TableZobristDistribution<
                                  FifteenPuzzle > , Executor >( d , x , o )
                         .solve() , ki.cost , "HDA*[" + x.name() + "]" + tag );
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
  HDAOptions lo( o );
  lo.outsource = 3;
  lo.sync_send = true;
  check_solution( g , HDAStar< Grid< 8 > , decltype( gz ) , Executor >(
                       g , gz , lo ).solve() , ref.cost ,
                  "HDA*[AZ] grid 8, synchronous sends and lending, p " +
                  std::to_string( p ) );

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
  for( const auto & [ planes , split ] : { std::pair( 1u , 1u ) ,
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
 const auto pra = []( const EightPuzzle & d , const HDAOptions & x ) {
  const ZobristDistribution< EightPuzzle > z( d );
  return( HDAStar< EightPuzzle , decltype( z ) , Executor ,
                   BucketOpenList< Node< EightPuzzle::state_type , int > > ,
                   LockedChannelsOf< Executor > >( d , z , x ).solve() );
  };
 HDAOptions x;
 x.threads = 8;
 x.sync_send = true;
 x.outsource = 1;
 check_astar_edges( x , pra , "PRA* with lending, 8 workers" );
 check_rounding< Executor >();
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
   PBNFOptions u = o;
   u.safe = false;
   check_solution( d , SafePBNF< FifteenPuzzle , Executor >( d , u ).solve() ,
                   ki.cost , "PBNF p " + std::to_string( p ) +
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
   PBNFOptions u = o;
   u.safe = false;
   if( r4.solved() )
    check_solution( g4 , SafePBNF< Grid< 4 > , Executor >( g4 , u ).solve() ,
                    r4.cost , "PBNF grid 4, p " + std::to_string( p ) );
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
 o.safe = false;
 check_astar_edges( o , pbnf , "PBNF 0 workers" );
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

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_hda< ADS::ThreadPool >( const std::string & );
template void test_spa< ADS::ThreadPool >( const std::string & );
template void test_mqa< ADS::ThreadPool >( const std::string & );
template void test_pbnf< ADS::ThreadPool >( const std::string & );
template void test_pase< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_hda< ADS::OpenMPExecutor >( const std::string & );
 template void test_spa< ADS::OpenMPExecutor >( const std::string & );
 template void test_mqa< ADS::OpenMPExecutor >( const std::string & );
 template void test_pbnf< ADS::OpenMPExecutor >( const std::string & );
 template void test_pase< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_hda< ADS::FastFlowExecutor >( const std::string & );
 template void test_pbnf< ADS::FastFlowExecutor >( const std::string & );
 template void test_pase< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File parallel.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
