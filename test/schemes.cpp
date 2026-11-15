/*--------------------------------------------------------------------------*/
/*---------------------------- File schemes.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the searches of the
 * parallel schemes of Burns et al. and of Jinnai and Fukunaga: the
 * breadth-first and cost-bounded depth-first baselines, the batch weighted
 * A*, PSDD and its variants, TDS and the portfolio of A* searches.
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
// the breadth-first search and the cost-bounded depth-first search

void test_bfs( const std::string & data )
{
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 std::size_t k = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 9007 ) != 0 )                    // 21 states, spread out
   continue;
  const EightPuzzle e( EightPuzzle::unpack( s ) );
  const auto sol = breadth_first( e );
  if( ! sol.solved() || ( sol.cost != dst ) )
   check_solution( e , sol , dst , "BFS 8-puzzle" );
  }
 check( failures == 0 , "BFS on 21 states of the 8-puzzle, unit costs" );

 // the path of fewest edges, Arad-Sibiu-Fagaras-Bucharest, not the
 // cheapest one, of cost 418
 auto rg = romania( Arad );
 check_solution( rg , breadth_first( rg ) , 450 ,
                 "BFS Romania, the path of fewest edges" );

 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const std::string tag = ", Korf #" + std::to_string( id );
  check_solution( d , cost_bounded_dfs( d , ki.cost ) , ki.cost ,
                  "cost-bounded DFS at the optimum" + tag );
  check( cost_bounded_dfs( d , ki.cost - 1 ).status == Status::unsolvable ,
         "cost-bounded DFS below the optimum" + tag );
  }
 const auto r = cost_bounded_dfs( rg , 450 );
 check( r.solved() && ( r.cost <= 450 ) && ( r.cost >= 418 ) ,
        "cost-bounded DFS Romania within 450" );
 check( cost_bounded_dfs( rg , 417 ).status == Status::unsolvable ,
        "cost-bounded DFS Romania below the optimum" );
 const EightPuzzle g( EightPuzzle::korf_goal() );
 check( cost_bounded_dfs( g , 0 ).solved() &&
        ( cost_bounded_dfs( EightPuzzle( walk8( { 1 } ) ) , 0 ).status ==
          Status::unsolvable ) , "cost-bounded DFS with the bound 0" );

 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 4 > w( m , 0 , 0 , 39 , 39 );
 check( breadth_first( w ).status == Status::unsolvable ,
        "BFS proves the wall unsolvable" );
 check_astar_edges( SearchOptions() , []( const EightPuzzle & d ,
                                          const SearchOptions & o ) {
  return( breadth_first( d , o ) );
  } , "BFS" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 15-puzzle with the heuristic computed in batches: 0 with the
// successors, the Manhattan distance by heuristic_batch(), which counts
// the states and the calls

struct BatchTiles {

 using state_type = state_t< FifteenPuzzle >;
 using cost_type = int;

 explicit BatchTiles( const FifteenPuzzle & x ) : d( x ) {}

 const FifteenPuzzle & d;
 mutable std::uint64_t states = 0;              // the states evaluated
 mutable std::uint64_t calls = 0;               // the batches

 state_type initial_state( void ) const { return( d.initial_state() ); }
 bool is_goal( const state_type & s ) const { return( d.is_goal( s ) ); }
 std::uint64_t hash( const state_type & s ) const { return( d.hash( s ) ); }
 int heuristic( const state_type & s ) const { return( d.heuristic( s ) ); }

 template< class F >
 void for_each_successor( const state_type & s , int , F && f ) const {
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_type & t , int c , int ) {
   f( t , c , 0 );
   } );
  }

 void heuristic_batch( std::span< const state_type > s ,
                       std::span< int > h ) const {
  ++calls;
  states += s.size();
  for( std::size_t i = 0 ; i < s.size() ; ++i )
   h[ i ] = d.heuristic( s[ i ] );
  }

 };
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the batch weighted A*: A* with batches of 1, optimal with the weight 1,
// within the weight otherwise, and the heuristic in batches

void test_batch( const std::string & data )
{
 static_assert( BatchHeuristicDomain< BatchTiles > &&
                ! BatchHeuristicDomain< FifteenPuzzle > );
 const auto inst = korf( data );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const std::string tag = ", Korf #" + std::to_string( id );
  const auto ref = astar( d );
  BatchOptions o;
  const auto b1 = batch_astar( d , o );
  check_solution( d , b1 , ki.cost , "batch A* K 1" + tag );
  check( b1.stats.expanded == ref.stats.expanded ,
         "batch A* K 1 expands as A*" + tag );
  for( std::size_t k : { std::size_t( 8 ) , std::size_t( 64 ) ,
                         std::size_t( 512 ) } ) {
   o.batch = k;
   check_solution( d , batch_astar( d , o ) , ki.cost ,
                   "batch A* K " + std::to_string( k ) + tag );
   }
  o.batch = 64;
  o.weight = 2;
  const auto w = batch_astar( d , o );
  check( w.solved() && ( w.cost <= 2 * ki.cost ) ,
         "batch A* K 64 w 2 within the weight" + tag );
  check_solution( d , w , w.cost , "batch A* K 64 w 2 path" + tag );

  // the heuristic in batches: once per state stored, the root included
  const BatchTiles bt( d );
  o.weight = 1;
  o.batch = 1;
  const auto t1 = batch_astar( bt , o );
  check( t1.solved() && ( t1.cost == ki.cost ) &&
         ( t1.stats.expanded == ref.stats.expanded ) &&
         ( bt.states == t1.stats.stored ) ,
         "batch A* K 1 with heuristic_batch" + tag );
  o.batch = 64;
  bt.states = bt.calls = 0;
  const auto t64 = batch_astar( bt , o );
  check( t64.solved() && ( t64.cost == ki.cost ) &&
         ( bt.states == t64.stats.stored ) &&
         ( bt.calls * 32 <= t64.stats.expanded ) ,
         "batch A* K 64 with heuristic_batch, " +
         std::to_string( bt.calls ) + " batches" + tag );
  }

 for( std::uint64_t seed : { 2u , 3u } ) {
  const auto m = GridMap::random( 300 , 200 , 0.25 , seed , 10 );
  const Grid< 8 > g8( m , 0 , 199 , 299 , 199 );
  const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
  if( ! r8.solved() )
   continue;
  BatchOptions o;
  o.batch = 100;
  check_solution( g8 , batch_astar( g8 , o ) , r8.cost ,
                  "batch A* K 100 grid 8, costs 1 to 10" );
  }

 BatchOptions o;
 o.batch = 0;
 o.weight = 0;
 check_astar_edges( o , []( const EightPuzzle & d , const BatchOptions & x ) {
  return( batch_astar( d , x ) );
  } , "batch A* parameters 0" );
 o.batch = 16;
 check_astar_edges( o , []( const EightPuzzle & d , const BatchOptions & x ) {
  return( batch_astar( d , x ) );
  } , "batch A* K 16" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// PSDD, BFPSDD, IDPSDD and DynamicBoundedPSDD with 1 to 8 workers, on the
// given runtime, also with costs other than 1

template< class Executor >
void test_psdd( const std::string & data )
{
 constexpr GridNBlocks nb = GridNBlocks::coarse( 2000 , 1200 , 50 , 50 );
 constexpr GridNBlocks nc = GridNBlocks::coarse( 2001 , 7 , 50 , 50 );
 static_assert( ( nb.width == 40 ) && ( nb.height == 24 ) &&
                ( nc.width == 39 ) && ( nc.height == 1 ) );
 const auto inst = korf( data );
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  PSDDOptions o;
  o.threads = p;
  const std::string tp = " p " + std::to_string( p );
  for( int id : { easy_korf[ 0 ] , easy_korf[ 5 ] } ) {
   if( ( id != easy_korf[ 0 ] ) && ( ( p > cores ) || ( p == 8 ) ) )
    break;
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   const std::string tag = tp + ", Korf #" + std::to_string( id );
   using P = PSDD< FifteenPuzzle , Executor >;
   check_solution( d , P( d , o ).solve( ki.cost ) , ki.cost ,
                   "PSDD bound at the optimum" + tag );
   int pruned = 0;
   check( ( P( d , o ).solve( ki.cost - 1 , & pruned ).status ==
            Status::unsolvable ) && ( pruned == ki.cost ) ,
          "PSDD bound below the optimum, least f pruned" + tag );
   check_solution( d , BFPSDD< FifteenPuzzle , Executor >( d , o ).solve() ,
                   ki.cost , "BFPSDD" + tag );
   if( id == easy_korf[ 0 ] )
    check_solution( d , IDPSDD< FifteenPuzzle , Executor >( d , o )
                        .solve() , ki.cost , "IDPSDD" + tag );
   }

  // the grids of the test of the batch A*, solvable, with 1 and 3 workers
  for( std::uint64_t seed : { 2u , 3u } ) {
   if( ( p != 1 ) && ( p != 3 ) )
    break;
   const auto m = GridMap::random( 300 , 200 , 0.25 , seed , 10 );
   const Grid< 4 > g4( m , 0 , 199 , 299 , 199 , 100 , 10 );
   const auto r4 = astar( Blind< Grid< 4 > >{ g4 } );
   const Grid< 8 > g8( m , 0 , 199 , 299 , 199 , 100 , 10 );
   const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
   const std::string tag = tp + ", costs 1 to 10, seed " +
                           std::to_string( seed );
   check( r4.solved() && r8.solved() , "the grids are solvable" + tag );
   check_solution( g4 , PSDD< Grid< 4 > , Executor >( g4 , o )
                        .solve( r4.cost ) , r4.cost ,
                   "PSDD grid 4 bound at the optimum" + tag );
   check_solution( g4 , BFPSDD< Grid< 4 > , Executor >( g4 , o ).solve() ,
                   r4.cost , "BFPSDD grid 4" + tag );
   check_solution( g4 , DynamicBoundedPSDD< Grid< 4 > , Executor >(
                         g4 , 1.5 , o ).solve() , r4.cost ,
                   "DynamicBoundedPSDD w 1.5 grid 4" + tag );
   check_solution( g8 , PSDD< Grid< 8 > , Executor >( g8 , o )
                        .solve( r8.cost ) , r8.cost ,
                   "PSDD grid 8 bound at the optimum" + tag );
   check_solution( g8 , BFPSDD< Grid< 8 > , Executor >( g8 , o ).solve() ,
                   r8.cost , "BFPSDD grid 8" + tag );
   }
  const auto ms = GridMap::random( 40 , 30 , 0.2 , 4 , 3 );
  const Grid< 8 > gs( ms , 0 , 29 , 39 , 29 , 100 , 6 );
  const auto rs = astar( Blind< Grid< 8 > >{ gs } );
  check( rs.solved() , "the grid 40 x 30 is solvable" + tp );
  PSDDOptions cr = o;
  cr.growth = 2;
  check_solution( gs , IDPSDD< Grid< 8 > , Executor >( gs , cr ).solve() ,
                  rs.cost , "IDPSDD-CR grid 8 40 x 30, costs 1 to 3" + tp );
  const Grid< 8 > gr( ms , 0 , 29 , 39 , 29 , 100 ,
                      GridNBlocks::coarse( 40 , 30 , 7 , 4 ) );
  check_solution( gr , BFPSDD< Grid< 8 > , Executor >( gr , o ).solve() ,
                  rs.cost , "BFPSDD grid 8 40 x 30, nblocks of 5 x 6 cells" +
                  tp );
  check_solution( gs , IDPSDD< Grid< 8 > , Executor >( gs , o ).solve() ,
                  rs.cost , "IDPSDD grid 8 40 x 30, costs 1 to 3" + tp );
  check_solution( gs , PSDD< Grid< 8 > , Executor >( gs , o ).solve() ,
                  rs.cost , "PSDD unbounded grid 8 40 x 30, costs 1 to 3" +
                  tp );
  }

 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 std::size_t k = 0;
 PSDDOptions o;
 o.threads = 3;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 30011 ) != 0 )                   // 7 states, spread out
   continue;
  const EightPuzzle e( EightPuzzle::unpack( s ) );
  check_solution( e , DynamicBoundedPSDD< EightPuzzle , Executor >( e , 1.5 ,
                                                                 o ).solve() ,
                  dst , "DynamicBoundedPSDD w 1.5 8-puzzle" );
  }

 GridMap m( 60 , 60 );
 for( int y = 0 ; y < 60 ; ++y )
  m.set_cost( 30 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 59 , 59 , 100 , 7 );
 o.threads = 4;
 check( ( PSDD< Grid< 4 > , Executor >( g , o ).solve().status ==
          Status::unsolvable ) &&
        ( BFPSDD< Grid< 4 > , Executor >( g , o ).solve().status ==
          Status::unsolvable ) &&
        ( IDPSDD< Grid< 4 > , Executor >( g , o ).solve().status ==
          Status::unsolvable ) , "PSDD, BFPSDD and IDPSDD prove the wall "
                                 "unsolvable" );

 const auto psdd = []( const EightPuzzle & d , const PSDDOptions & x ) {
  return( PSDD< EightPuzzle , Executor >( d , x ).solve() );
  };
 const auto bfpsdd = []( const EightPuzzle & d , const PSDDOptions & x ) {
  return( BFPSDD< EightPuzzle , Executor >( d , x ).solve() );
  };
 o.threads = 8;
 check_astar_edges( o , psdd , "PSDD 8 workers" );
 check_astar_edges( o , bfpsdd , "BFPSDD 8 workers" );
 o.threads = 0;
 o.multiplier = 0;
 o.min_expansions = 0;
 check_astar_edges( o , bfpsdd , "BFPSDD parameters 0" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// TDS with 1 to 8 workers, on the given runtime

template< class Executor >
void test_tds( const std::string & data )
{
 const auto inst = korf( data );
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  ParallelOptions o;
  o.threads = p;
  const std::string tp = " p " + std::to_string( p );
  for( int id : { easy_korf[ 0 ] , easy_korf[ 5 ] } ) {
   if( p > cores )
    break;
   const auto & ki = inst[ id - 1 ];
   const FifteenPuzzle d( ki.board );
   const std::string tag = tp + ", Korf #" + std::to_string( id );
   const ZobristDistribution< FifteenPuzzle > z( d );
   const auto sz = TDS< FifteenPuzzle , decltype( z ) , Executor >( d , z , o )
                   .solve();
   check_solution( d , sz , ki.cost , "TDS[Z]" + tag );
   if( id != easy_korf[ 0 ] )
    continue;
   const AbstractZobristDistribution< FifteenPuzzle > az( d );
   check_solution( d , TDS< FifteenPuzzle , decltype( az ) , Executor >(
                        d , az , o ).solve() , ki.cost , "TDS[AZ]" + tag );
   if( p == 1 )
    check( sz.stats.expanded ==
           TDS< FifteenPuzzle , decltype( z ) , Executor >( d , z , o )
           .solve().stats.expanded ,
           "TDS with one worker expands the same at every run" + tag );
   }
  std::size_t k = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 30011 ) != 0 )                  // 7 states, spread out
    continue;
   const EightPuzzle e( EightPuzzle::unpack( s ) );
   const ZobristDistribution< EightPuzzle > z( e );
   const auto sol = TDS< EightPuzzle , decltype( z ) , Executor >( e , z , o )
                    .solve();
   if( ! sol.solved() || ( sol.cost != dst ) )
    check_solution( e , sol , dst , "TDS 8-puzzle" + tp );
   }
  check( failures == 0 , "TDS on 7 states of the 8-puzzle" + tp );
  const auto m = GridMap::random( 30 , 20 , 0.2 , 4 , 3 );
  const Grid< 8 > g8( m , 0 , 19 , 29 , 19 , 10 );
  const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
  const AbstractZobristDistribution< Grid< 8 > > gz( g8 );
  if( r8.solved() ) {
   check_solution( g8 , TDS< Grid< 8 > , decltype( gz ) , Executor >(
                         g8 , gz , o ).solve() , r8.cost ,
                   "TDS[AZ] grid 8 30 x 20, costs 1 to 3" + tp );
   TDSOptions cr( o );
   cr.growth = 2;
   check_solution( g8 , TDS< Grid< 8 > , decltype( gz ) , Executor >(
                         g8 , gz , cr ).solve() , r8.cost ,
                   "TDS-CR grid 8 30 x 20" + tp );
   cr.optimal = false;
   const auto f = TDS< Grid< 8 > , decltype( gz ) , Executor >( g8 , gz , cr )
                  .solve();
   check_solution( g8 , f , f.cost , "TDS-CR first path, grid 8" + tp );
   check( f.cost >= r8.cost , "TDS-CR first path not below the optimum" + tp );
   }
  auto rg = romania( Arad );
  const ZobristDistribution< ExplicitGraph< int > > rz( rg );
  check_solution( rg , TDS< ExplicitGraph< int > , decltype( rz ) ,
                            Executor >( rg , rz , o ).solve() , 418 ,
                  "TDS[Z] Romania" + tp );
  }
 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 39 , 39 );
 const StateHashDistribution< Grid< 4 > > gh( g );
 ParallelOptions o;
 o.threads = 3;
 check( TDS< Grid< 4 > , decltype( gh ) , Executor >( g , gh , o ).solve()
        .status == Status::unsolvable , "TDS proves the wall unsolvable" );

 const auto tds = []( const EightPuzzle & e , const ParallelOptions & x ) {
  const ZobristDistribution< EightPuzzle > z( e );
  return( TDS< EightPuzzle , decltype( z ) , Executor >( e , z , x )
          .solve() );
  };
 o.threads = 8;
 check_astar_edges( o , tds , "TDS 8 workers" );
 o.threads = 0;
 o.batch = 0;
 o.flush_period = 0;
 check_astar_edges( o , tds , "TDS parameters 0" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the portfolio of A* searches with 1 to 8 workers, on the given runtime

template< class Executor >
void test_multi( const std::string & data )
{
 const auto inst = korf( data );
 const auto & ki = inst[ easy_korf[ 0 ] - 1 ];
 const FifteenPuzzle d( ki.board );
 ParallelOptions o;
 const auto one = MultiAStar< FifteenPuzzle , Executor >( d , o ).solve();
 check_solution( d , one , ki.cost , "MultiAStar p 1" );
 check( one.stats.expanded == astar( d ).stats.expanded ,
        "MultiAStar with one worker expands as A*" );
 for( unsigned p : { 2u , 3u , 8u } ) {
  o.threads = p;
  const std::string tp = " p " + std::to_string( p );
  const auto s = MultiAStar< FifteenPuzzle , Executor >( d , o ).solve();
  check_solution( d , s , ki.cost , "MultiAStar Korf #" +
                                    std::to_string( easy_korf[ 0 ] ) + tp );
  check( s.stats.expanded_by.size() == p ,
         "MultiAStar counts each worker" + tp );
  const auto m = GridMap::random( 200 , 120 , 0.25 , 3 , 10 );
  const Grid< 8 > g8( m , 0 , 119 , 199 , 119 );
  const auto r8 = astar( Blind< Grid< 8 > >{ g8 } );
  if( r8.solved() )
   check_solution( g8 , MultiAStar< Grid< 8 > , Executor >( g8 , o )
                        .solve() , r8.cost , "MultiAStar grid 8" + tp );
  auto rg = romania( Arad );
  check_solution( rg , MultiAStar< ExplicitGraph< int > , Executor >( rg ,
                                                                       o )
                       .solve() , 418 , "MultiAStar Romania" + tp );
  }
 GridMap m( 40 , 40 );
 for( int y = 0 ; y < 40 ; ++y )
  m.set_cost( 20 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 39 , 39 );
 o.threads = 3;
 check( MultiAStar< Grid< 4 > , Executor >( g , o ).solve().status ==
        Status::unsolvable , "MultiAStar proves the wall unsolvable" );

 const auto multi = []( const EightPuzzle & e , const ParallelOptions & x ) {
  return( MultiAStar< EightPuzzle , Executor >( e , x ).solve() );
  };
 o.threads = 8;
 check_astar_edges( o , multi , "MultiAStar 8 workers" );
 o.threads = 0;
 check_astar_edges( o , multi , "MultiAStar 0 workers" );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_psdd< ADS::ThreadPool >( const std::string & );
template void test_tds< ADS::ThreadPool >( const std::string & );
template void test_multi< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_psdd< ADS::OpenMPExecutor >( const std::string & );
 template void test_tds< ADS::OpenMPExecutor >( const std::string & );
 template void test_multi< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_psdd< ADS::FastFlowExecutor >( const std::string & );
 template void test_tds< ADS::FastFlowExecutor >( const std::string & );
 template void test_multi< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- End File schemes.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
