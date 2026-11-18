/*--------------------------------------------------------------------------*/
/*------------------------ File pbidirectional.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the parallel
 * bidirectional searches.
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

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_bidirectional_parallel< ADS::ThreadPool >
  ( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_bidirectional_parallel< ADS::OpenMPExecutor >
   ( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_bidirectional_parallel< ADS::FastFlowExecutor >
   ( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*---------------------- End File pbidirectional.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
