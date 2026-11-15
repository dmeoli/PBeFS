/*--------------------------------------------------------------------------*/
/*---------------------------- File search.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the sequential
 * searches of the 8-puzzle, the 15-puzzle, the 24-puzzle, the grids and the
 * graphs: A* with both OPEN lists, IDA*, the limits and DDD.
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
// sequential A*, both OPEN lists, on the 8-puzzle and the 15-puzzle

void test_astar( const std::string & data )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 check( dist.size() == 181440 , "the 8-puzzle has 9! / 2 states" );
 std::size_t k = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 907 ) != 0 )                     // 200 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto sb = astar( d , OpenKind::buckets );
  const auto sh = astar( d , OpenKind::heap );
  if( ! ( sb.solved() && ( sb.cost == dst ) && sh.solved() &&
          ( sh.cost == dst ) ) ) {
   check_solution( d , sb , dst , "8-puzzle, buckets" );
   check_solution( d , sh , dst , "8-puzzle, heap" );
   }
  }
 check( failures == 0 , "A* on 200 states of the 8-puzzle" );

 const auto inst = korf( data );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  check_solution( d , astar( d , OpenKind::buckets ) , ki.cost ,
                  "A* buckets, Korf #" + std::to_string( id ) );
  check_solution( d , astar( d , OpenKind::heap ) , ki.cost ,
                  "A* heap, Korf #" + std::to_string( id ) );
  }

 for( OpenKind k : { OpenKind::buckets , OpenKind::heap } )
  check_astar_edges( SearchOptions() , [ k ]( const EightPuzzle & d ,
                                              const SearchOptions & o ) {
   return( astar( d , k , o ) );
   } , std::string( "A* " ) + ( k == OpenKind::heap ? "heap" : "buckets" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// both OPEN lists against an ordered set on random pushes and pops: the
// least f moving back below the f of the last pop, as in a worker of HDA*,
// and the values of g below and above BucketOpenList::DenseG, the first of
// these making the layers heaps while they hold entries

template< class Open >
static void check_open( const std::string & name )
{
 using N = Node< unsigned int , int >;
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed )
  for( int top_g : { 200 , 2000 } ) {
   std::mt19937_64 rng( seed );
   std::deque< N > nodes;
   Open open;
   std::multiset< std::pair< int , int > > ref;  // ( f , - g )
   bool ok = true;
   int last_f = 0;
   for( int i = 0 ; ( i < 20000 ) && ok ; ++i ) {
    if( ( rng() % 3 != 0 ) || ref.empty() ) {
     // a g up to top_g, grown along the run, and an f near the last pop
     const int g = static_cast< int >( rng() % ( 1 + ( top_g * i ) / 20000 ) );
     const int f = std::max( g , last_f - 20 +
                                 static_cast< int >( rng() % 60 ) );
     nodes.push_back( N{ 0u , g , f - g , nullptr , false , false } );
     open.push( &nodes.back() , g );
     ref.emplace( f , -g );
     }
    else {
     const auto e = open.pop();
     const auto it = ref.begin();
     ok = ( e.g == -it->second ) && ( e.g + e.node->h == it->first );
     last_f = it->first;
     ref.erase( it );
     }
    ok = ok && ( open.size() == ref.size() );
    }
   for( ; ok && ! ref.empty() ; ref.erase( ref.begin() ) ) {
    const auto e = open.pop();
    ok = ( e.g == -ref.begin()->second ) &&
         ( e.g + e.node->h == ref.begin()->first );
    }
   check( ok && open.empty() , name + ", seed " + std::to_string( seed ) +
                               ", g up to " + std::to_string( top_g ) );
   }
 }

void test_open( void )
{
 using N = Node< unsigned int , int >;
 check_open< HeapOpenList< N > >( "the heap against an ordered set" );
 check_open< BucketOpenList< N > >( "the buckets against an ordered set" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* on random grids against Dijkstra, with 4 and 8 moves, on maps of unit
// costs and of costs 1 to 5

void test_grid( void )
{
 for( std::uint64_t seed = 1 ; seed <= 6 ; ++seed ) {
  const auto m = GridMap::random( 300 , 200 , 0.25 , seed );
  const Grid< 4 > g4( m , 0 , 199 , 299 , 199 );
  const auto ref4 = astar( Blind< Grid< 4 > >{ g4 } );
  if( ref4.solved() )
   check_solution( g4 , astar( g4 ) , ref4.cost ,
                   "grid 4, seed " + std::to_string( seed ) );
  else
   check( ! astar( g4 ).solved() , "grid 4 unsolvable, seed " +
                                    std::to_string( seed ) );
  const Grid< 8 > g8( m , 0 , 199 , 299 , 199 );
  const auto ref8 = astar( Blind< Grid< 8 > >{ g8 } );
  if( ref8.solved() )
   check_solution( g8 , astar( g8 ) , ref8.cost ,
                   "grid 8, seed " + std::to_string( seed ) );
  }

 // the weighted map of the old FindPath example: costs 1 to 5, 9 = wall
 GridMap w( 20 , 20 );
 for( int y = 0 ; y < 20 ; ++y )
  for( int x = 0 ; x < 20 ; ++x )
   w.set_cost( x , y , 1 + ( 7 * x + 3 * y ) % 5 );
 for( int y = 2 ; y < 18 ; ++y )
  w.set_cost( 10 , y , 0 );
 const Grid< 4 > gw( w , 3 , 6 , 17 , 15 );
 check_solution( gw , astar( gw ) , astar( Blind< Grid< 4 > >{ gw } ).cost ,
                 "weighted grid 4" );
 const double dw4 = grid_dijkstra( w , 4 , 3 , 6 , 17 , 15 );
 check_solution( gw , astar( gw ) , int( std::lround( dw4 ) ) ,
                 "weighted grid 4 against Dijkstra on the map" );
 const Grid< 8 > gw8( w , 3 , 6 , 17 , 15 );
 const double dw8 = grid_dijkstra( w , 8 , 3 , 6 , 17 , 15 );
 check_solution( gw8 , astar( gw8 ) , dw8 ,
                 "weighted grid 8 against Dijkstra on the map" );

 // random maps of costs 1 to 5, whose moves cost differently both ways
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const auto mc = GridMap::random( 120 , 80 , 0.25 , seed , 5 );
  const std::string t = ", costs 1 to 5, seed " + std::to_string( seed );
  const double r4 = grid_dijkstra( mc , 4 , 0 , 79 , 119 , 79 );
  const Grid< 4 > c4( mc , 0 , 79 , 119 , 79 );
  if( r4 >= 0 )
   check_solution( c4 , astar( c4 ) , int( std::lround( r4 ) ) ,
                   "grid 4" + t );
  else
   check( ! astar( c4 ).solved() , "grid 4 unsolvable" + t );
  const double r8 = grid_dijkstra( mc , 8 , 0 , 79 , 119 , 79 );
  const Grid< 8 > c8( mc , 0 , 79 , 119 , 79 );
  if( r8 >= 0 )
   check_solution( c8 , astar( c8 ) , r8 , "grid 8" + t );
  else
   check( ! astar( c8 ).solved() , "grid 8 unsolvable" + t );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* on the map of Romania

void test_graph( void )
{
 auto g = romania( Arad );
 const auto s = astar( g );
 check_solution( g , s , 418 , "Romania, Arad to Bucharest" );
 check( ( s.path.size() == 5 ) && ( s.path[ 1 ] == Sibiu ) &&
        ( s.path[ 2 ] == RimnicuVilcea ) && ( s.path[ 3 ] == Pitesti ) ,
        "Romania, through Sibiu, Rimnicu Vilcea and Pitesti" );
 for( int c = 0 ; c < n_cities ; ++c ) {
  auto gc = romania( City( c ) );
  check_solution( gc , astar( gc ) , romania_distance( City( c ) ) ,
                  "Romania, city " + std::to_string( c ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the node limit and the unsolvable instances

void test_limits( const std::string & data )
{
 const auto inst = korf( data );
 const FifteenPuzzle d( inst[ 0 ].board );
 SearchOptions o;
 o.node_limit = 10000;
 const auto s = astar( d , OpenKind::automatic , o );
 check( ( s.status == Status::node_limit ) && ( s.stats.stored >= 10000 ) &&
        s.path.empty() , "A* stops at the node limit" );

 GridMap m( 50 , 50 );
 for( int y = 0 ; y < 50 ; ++y )
  m.set_cost( 25 , y , 0 );
 const Grid< 4 > g( m , 0 , 0 , 49 , 49 );
 const auto u = astar( g );
 check( ( u.status == Status::unsolvable ) && ( u.stats.expanded == 25 * 50 ) ,
        "A* proves the wall unsolvable after the 1250 cells it reaches" );

 ParallelOptions po;
 po.threads = 4;
 const StateHashDistribution< Grid< 4 > > hd( g );
 const auto hu = HDAStar< Grid< 4 > , decltype( hd ) >( g , hd , po ).solve();
 check( hu.status == Status::unsolvable , "HDA* proves the wall unsolvable" );
 const auto su = SPAStar< Grid< 4 > >( g , po ).solve();
 check( su.status == Status::unsolvable , "SPA* proves the wall unsolvable" );
 const auto mu = MQAStar< Grid< 4 > >( g , po ).solve();
 check( mu.status == Status::unsolvable ,
        "MultiQueue A* proves the wall unsolvable" );

 po.node_limit = 5000;
 const ZobristDistribution< FifteenPuzzle > z( d );
 const auto hl = HDAStar< FifteenPuzzle , decltype( z ) >( d , z , po )
                  .solve();
 check( hl.status == Status::node_limit , "HDA* stops at the node limit" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// IDA* against the exact distances of the 8-puzzle, the optima of Korf on
// the 15-puzzle, A* on the grids and the map of Romania (where the parent
// is pruned by its state, with no operator), and the node limit

void test_idastar( const std::string & data )
{
 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 std::size_t k = 0;
 int wrong = 0;
 for( const auto & [ s , dst ] : dist ) {
  if( ( k++ % 1814 ) != 0 )                    // 100 states, spread out
   continue;
  const EightPuzzle d( EightPuzzle::unpack( s ) );
  const auto r = idastar( d );
  if( ! ( r.solved() && ( r.cost == dst ) ) ) {
   ++wrong;
   check_solution( d , r , dst , "IDA* 8-puzzle" );
   }
  }
 check( wrong == 0 , "IDA* on 100 states of the 8-puzzle" );

 const auto inst = korf( data );
 for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] , easy_korf[ 2 ] } ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto r = idastar( d );
  check_solution( d , r , ki.cost , "IDA*, Korf #" + std::to_string( id ) );
  check( r.stats.stored == std::uint64_t( ki.cost ) + 1 , "IDA* holds only "
         "the path, Korf #" + std::to_string( id ) );
  }

 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 7 , 5 , 0.2 , seed , 3 );
  const Grid< 4 > g4( m , 0 , 4 , 6 , 4 );
  const auto a4 = astar( g4 );
  if( a4.solved() )
   check_solution( g4 , idastar( g4 ) , a4.cost , "IDA* grid 4, seed " +
                   std::to_string( seed ) );
  const Grid< 8 > g8( m , 0 , 4 , 6 , 4 );
  const auto a8 = astar( g8 );
  if( a8.solved() )
   check_solution( g8 , idastar( g8 ) , a8.cost , "IDA* grid 8, seed " +
                   std::to_string( seed ) );
  }

 auto g = romania( Arad );
 check_solution( g , idastar( g ) , 418 , "IDA* Romania, Arad to Bucharest" );

 ExplicitGraph< int > dag( 4 , 0 , 3 );         // the target not reachable
 dag.add_edge( 0 , 1 , 1 , true );
 dag.add_edge( 1 , 2 , 1 , true );
 dag.add_edge( 0 , 2 , 1 , true );
 check( idastar( dag ).status == Status::unsolvable ,
        "IDA* proves an acyclic graph unsolvable" );

 SearchOptions o;
 o.node_limit = 1000;
 const FifteenPuzzle d( inst[ 0 ].board );
 const auto l = idastar( d , o );
 check( ( l.status == Status::node_limit ) && ( l.stats.expanded == 1000 ) &&
        l.path.empty() , "IDA* stops at the expansion limit" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 24-puzzle: the integer of 128 bits, the packed board, A* and IDA*
// with the Manhattan distance and with a pattern database, the parallel
// searches, and the 50 instances of the data

void test_tiles24( const std::string & data )
{
 const Word128 one( 1 );
 const Word128 w = one << 70;
 check( ( w.lo == 0 ) && ( w.hi == 64 ) && ( ( w >> 70 ) == one ) &&
        ( countr_zero( w ) == 70 ) && ( countr_zero( Word128() ) == 128 ) &&
        ( ( Word128( ~std::uint64_t( 0 ) ) << 4 ).hi == 15 ) &&
        ( ( Word128( 0 , 1 ) >> 1 ).lo == ( std::uint64_t( 1 ) << 63 ) ) &&
        ( ( ~Word128() & w ) == w ) , "the shifts and the bits of Word128" );

 using D = TwentyFourPuzzle;
 static_assert( std::is_same_v< D::state_type , Word128 > &&
                std::is_trivially_copyable_v< Word128 > &&
                ! PerfectHashDomain< D > && PerfectHashDomain< FifteenPuzzle >
                && OperatorDomain< D > && FeatureDomain< D > );
 bool cells = true;
 for( int c = 0 ; c < D::cells ; ++c ) {
  auto b = D::korf_goal();
  std::swap( b[ 0 ] , b[ c ] );
  const auto s = D::pack( b );
  cells = cells && ( D::blank( s ) == c ) && ( D::unpack( s ) == b );
  for( int i = 0 ; i < D::cells ; ++i )
   cells = cells && ( D::tile( s , i ) == b[ i ] );
  }
 check( cells , "the blank and the tiles of every cell of the 24-puzzle" );

 PDBOptions po;
 Partition triples;
 for( int t = 1 ; t < D::cells ; t += 3 )
  triples.push_back( { t , t + 1 , t + 2 } );
 const AdditivePDB< 5 , 5 , NibbleEntries > pdb( triples , D::korf_goal() ,
                                                 po );
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const auto b = random_walk< D >( 40 , seed );
  const D d( b );
  const std::string tag = ", walk " + std::to_string( seed );
  const auto a = astar( d );
  check( a.solved() , "A* on the 24-puzzle" + tag );
  check_solution( d , a , a.cost , "A* on the 24-puzzle" + tag );
  check_solution( d , idastar( d ) , a.cost , "IDA* on the 24-puzzle" + tag );
  const SlidingTilePDB< 5 , 5 , NibbleEntries > dp( b , pdb );
  check( dp.heuristic( dp.initial_state() ) >= d.heuristic(
          d.initial_state() ) , "the 3-tile patterns dominate Manhattan" +
         tag );
  check_solution( dp , idastar( dp ) , a.cost , "IDA* with the 3-tile "
                  "patterns" + tag );
  ParallelOptions o;
  o.threads = 3;
  const ZobristDistribution< D > z( d );
  check_solution( d , HDAStar< D , decltype( z ) >( d , z , o ).solve() ,
                  a.cost , "HDA*[Z] on the 24-puzzle" + tag );
  check_solution( d , SPAStar< D >( d , o ).solve() , a.cost ,
                  "SPA* on the 24-puzzle" + tag );
  }

 std::ifstream in( data + "/24-puzzle/instances50.txt" );
 std::string line;
 int n = 0 , bad = 0;
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  try {
   const D d( D::parse( line ) );
   if( d.heuristic( d.initial_state() ) <= 0 )
    ++bad;
   }
  catch( std::invalid_argument & ) {
   ++bad;
   }
  ++n;
  }
 check( ( n == 50 ) && ( bad == 0 ) , "the 50 instances of the 24-puzzle "
        "are boards from which the goal is reachable" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* with delayed duplicate detection on the disk against the exact
// distances of the 8-puzzle (with one file per bucket, a prime number of
// them and the default), A* on the 15-puzzle, the grids and Romania; the
// directory of its files removed at the end, the unsolvable wall, the node
// limit, and the errors of an inconsistent heuristic and of no partition

void test_ddd( const std::string & data )
{
 ADS::TempDirectory dir;                       // where the runs go
 DDDOptions o;
 o.directory = dir.path();

 const EightPuzzle any( EightPuzzle::korf_goal() );
 const auto dist = eight_distances( any );
 for( std::size_t parts :
      { std::size_t( 1 ) , std::size_t( 7 ) , std::size_t( 64 ) } ) {
  o.partitions = parts;
  std::size_t k = 0;
  int wrong = 0;
  for( const auto & [ s , dst ] : dist ) {
   if( ( k++ % 3628 ) != 0 )                   // 50 states, spread out
    continue;
   const EightPuzzle d( EightPuzzle::unpack( s ) );
   const auto r = ddd_astar( d , {} , o );
   if( ! ( r.solved() && ( r.cost == dst ) ) ) {
    ++wrong;
    check_solution( d , r , dst , "DDD-A* 8-puzzle" );
    }
   }
  check( wrong == 0 , "DDD-A* on 50 states of the 8-puzzle, " +
                       std::to_string( parts ) + " partitions" );
  }
 check( std::filesystem::is_empty( dir.path() ) ,
        "DDD-A* removes the directory of its files" );

 o.partitions = 16;
 o.buffer_bytes = 4096;                        // some records per write
 const auto inst = korf( data );
 {
  const int id = easy_korf[ 0 ];               // 2 * 10^5 expansions
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle d( ki.board );
  const auto r = ddd_astar( d , {} , o );
  check_solution( d , r , ki.cost , "DDD-A*, Korf #" + std::to_string( id ) );
  check( r.stats.expanded <= r.stats.stored , "DDD-A* expands a state at"
                                             " most once, Korf #" +
                                               std::to_string( id ) );
  }

 o = DDDOptions{};
 o.directory = dir.path();
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  const auto m = GridMap::random( 7 , 5 , 0.2 , seed , 3 );
  const Grid< 4 > g4( m , 0 , 4 , 6 , 4 );
  const auto a4 = astar( g4 );
  const auto r4 = ddd_astar( g4 , {} , o );
  if( a4.solved() )
   check_solution( g4 , r4 , a4.cost ,
                   "DDD-A* grid 4, seed " + std::to_string( seed ) );
  else
   check( r4.status == Status::unsolvable ,
          "DDD-A* unsolvable grid, seed " + std::to_string( seed ) );
  }
 auto g = romania( Arad );
 check_solution( g , ddd_astar( g , {} , o ) , 418 ,
                 "DDD-A* Romania, Arad to Bucharest" );

 GridMap w( 50 , 50 );
 for( int y = 0 ; y < 50 ; ++y )
  w.set_cost( 25 , y , 0 );
 const Grid< 4 > wall( w , 0 , 0 , 49 , 49 );
 const auto u = ddd_astar( wall , {} , o );
 check( ( u.status == Status::unsolvable ) &&
          ( u.stats.expanded == 25 * 50 ) && ( u.stats.stored == 25 * 50 ) ,
        "DDD-A* proves the wall unsolvable"
        " after the 1250 cells it reaches" );
 check( u.stats.generated + 1 == u.stats.stored + u.stats.duplicates ,
        "DDD-A* stores each record of the wall once or drops it as a copy" );

 SearchOptions lim;
 lim.node_limit = 10000;
 const FifteenPuzzle hard( inst[ 0 ].board );
 const auto l = ddd_astar( hard , lim , o );
 check( ( l.status == Status::node_limit ) && ( l.stats.stored >= 10000 ) &&
          l.path.empty() ,
        "DDD-A* stops at the node limit" );

 ExplicitGraph< int > bad( 5 , 0 , 4 );         // that of test_reopen
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
  ( void ) ddd_astar( bad , {} , o );
  }
 catch( const std::logic_error & ) {
  threw = true;
  }
 check( threw , "DDD-A* refuses an inconsistent heuristic" );
 threw = false;
 try {
  DDDOptions none = o;
  none.partitions = 0;
  ( void ) ddd_astar( g , {} , none );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "DDD-A* refuses no partition" );
 check( std::filesystem::is_empty( dir.path() ) ,
        "DDD-A* leaves no file, also when it stops or throws" );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File search.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
