/*--------------------------------------------------------------------------*/
/*------------------------------ File pem.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the parallel external
 * MM and A* with immediate duplicate detection.
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
// the parallel external-memory searches with 1, 2 and 4 workers: PEM-A*,
// PEM-rA*, PEMM and PEM-BAE* on the 8-puzzle, the grids and Romania, and on
// random walks of the 15-puzzle with 2 workers, the directory removed at
// the end, the cap on the disk, the node limit and the inconsistent
// heuristic refused

// the PEM searches run on 3 states of the 8-puzzle and on 2 walks of 40
// moves of the 15-puzzle, and in a light run [see light_run()] with 1 and 4
// workers only, on 2 states and on 1 walk of 30 moves

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
 const std::size_t stride = light_run() ? 90720 : 60482;
 for( unsigned p : case_workers( { 1 , 2 , 4 } , { 1 , 4 } ) ) {
  std::size_t k = 0;
  int wrong = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % stride ) != 0 )
    continue;
   const int before = failures;
   all( EightPuzzle( EightPuzzle::unpack( s ) ) , dst , "8-puzzle" , p );
   wrong += failures - before;
   }
  check( wrong == 0 , "the PEM searches on " +
                      std::to_string( ( dist.size() + stride - 1 ) /
                                      stride ) +
                      " states of the 8-puzzle, " + rt + " p " +
                      std::to_string( p ) );
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
 for( std::uint64_t seed = 1 ; seed <= ( light_run() ? 1u : 2u ) ;
      ++seed ) {
  const FifteenPuzzle w( random_walk< FifteenPuzzle >( light_run() ? 30
                                                                 : 40 ,
                                                      seed ) );
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

template void test_pem< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_pem< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_pem< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- End File pem.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
