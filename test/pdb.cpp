/*--------------------------------------------------------------------------*/
/*------------------------------ File pdb.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the pattern databases
 * and the searches with them.
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
// the pattern databases: the values of the tables against a search of the
// abstraction written apart, the admissibility, the consistency and the
// encodings on the whole 8-puzzle, the reflection, the 15-puzzle against
// the optimal costs of Korf, and the files

void test_pdb( const std::string & data )
{
 // the tables against the abstract search, with two goals, both rules
 auto centre = EightPuzzle::korf_goal();          // the blank in the middle
 std::swap( centre[ 0 ] , centre[ 4 ] );
 for( const auto & goal : { EightPuzzle::korf_goal() , centre } )
  for( const auto & tiles : std::vector< std::vector< int > >{
        { 1 } , { 1 , 2 , 3 } , { 5 , 6 , 7 , 8 } , { 2 , 4 , 6 , 8 } ,
        { 1 , 2 , 3 , 4 , 5 , 6 , 7 } } )
   for( BlankRule rule : { BlankRule::ignored , BlankRule::minimized } ) {
    const TilePattern< 3 , 3 > pat( tiles , goal );
    using Table = PatternTable< 3 , 3 , ByteEntries >;
    const auto h1 = Table::distances( pat , rule , 1 );
    const auto h3 = Table::distances( pat , rule , 3 );
    const auto ref = abstract_distances( tiles , goal ,
                                         rule == BlankRule::minimized );
    bool same = ( h1 == h3 ) && ( ref.size() == h1.size() );
    std::uint8_t p[ 9 ];
    for( std::uint64_t r = 0 ; same && ( r < h1.size() ) ; ++r ) {
     pat.unrank( r , p );
     const std::vector< int > v( p , p + tiles.size() );
     same = ( pat.rank( p ) == r ) && ( ref.at( v ) == h1[ r ] );
     }
    std::string name = "pattern {";
    for( int t : tiles )
     name += " " + std::to_string( t );
    check( same , name + " }, blank " +
                  ( rule == BlankRule::ignored ? "ignored" : "minimized" ) +
                  ( goal == centre ? ", blank in the middle" : "" ) +
                  ": the abstract distances, on 1 and 3 workers" );
    }

 // the heuristic on the whole 8-puzzle, two partitions
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 const auto goal8 = EightPuzzle::korf_goal();
 for( const Partition & part : { Partition{ { 1 , 2 , 3 , 4 } ,
                                            { 5 , 6 , 7 , 8 } } ,
                                 Partition{ { 1 , 2 , 3 } , { 4 , 5 , 6 } ,
                                            { 7 , 8 } } } ) {
  const std::string tag = part.size() == 2 ? " 4-4" : " 3-3-2";
  for( bool reflect : { false , true } ) {
   PDBOptions o;
   o.reflect = reflect;
   o.threads = 2;
   const std::string rt = tag + ( reflect ? " reflected" : "" );
   const AdditivePDB< 3 , 3 , ByteEntries > byte( part , goal8 , o );
   check( byte.reflected() == reflect , "the reflection" + rt );
   check_heuristic8( byte , dist , true , "byte" + rt );
   const AdditivePDB< 3 , 3 , NibbleEntries > nib( part , goal8 , o );
   check_heuristic8( nib , dist , true , "nibble" + rt , & byte );
   bool clamped = false;
   for( const auto & t : nib.tables() )
    clamped = clamped || t.entries().clamped();
   check( ! clamped , "nibble" + rt + ": no excess beyond 15" );
   o.compress = 4;
   const AdditivePDB< 3 , 3 , MinNibbleEntries > mn( part , goal8 , o );
   check_heuristic8( mn , dist , false , "min 4" + rt , & byte );
   o.compress = 1;
   o.blank = BlankRule::minimized;
   const AdditivePDB< 3 , 3 , ByteEntries > mb( part , goal8 , o );
   const auto bad = check_heuristic8( mb , dist , false ,
                                      "minimized" + rt );
   std::cout << "  minimized" << rt << ": inconsistent along " << bad
             << " moves" << std::endl;
   o.blank = BlankRule::ignored;
   bool dom = true;
   for( const auto & [ s , d ] : dist ) {
    const auto q = cells8( s );
    dom = dom && ( byte.heuristic( q.data() ) <= mb.heuristic( q.data() ) ) &&
          ( byte.heuristic( q.data() ) == nib.heuristic( q.data() ) );
    }
   check( dom , "nibble = byte <= minimized" + rt );
   }

  // the residues: without the reflection, the value of every state from
  // the residues alone, and that of every successor from its parent
  PDBOptions o;
  o.reflect = false;
  const AdditivePDB< 3 , 3 , ByteEntries > byte( part , goal8 , o );
  const AdditivePDB< 3 , 3 , Mod3Entries > m3( part , goal8 , o );
  const AdditivePDB< 3 , 3 , ParityBitEntries > bit( part , goal8 , o );
  check_heuristic8( m3 , dist , true , "mod 3" + tag , & byte );
  check_heuristic8( bit , dist , true , "parity bit" + tag , & byte );
  const EightPuzzle any( goal8 );
  bool exact3 = true , exact4 = true , step = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   const int hb = byte.heuristic( q.data() );
   exact3 = exact3 && ( m3.exact( q.data() ) == hb );
   exact4 = exact4 && ( bit.exact( q.data() ) == hb );
   const SlidingTilePDB< 3 , 3 , Mod3Entries > d3( EightPuzzle::unpack( s ) ,
                                                   m3 );
   const SlidingTilePDB< 3 , 3 , ParityBitEntries > d4(
    EightPuzzle::unpack( s ) , bit );
   d3.for_each_successor( s , hb , [ & ]( std::uint64_t t , int , int ht ) {
    step = step && ( ht == byte.heuristic( cells8( t ).data() ) );
    } );
   d4.for_each_successor( s , hb , [ & ]( std::uint64_t t , int , int ht ) {
    step = step && ( ht == byte.heuristic( cells8( t ).data() ) );
    } );
   }
  check( exact3 && exact4 , "mod 3 and parity bit" + tag +
                            ": the exact values from the residues alone" );
  check( step , "mod 3 and parity bit" + tag +
                ": the value of every successor from its parent" );
  bool thrown = false;
  try {
   o.reflect = true;
   const AdditivePDB< 3 , 3 , Mod3Entries > r3( part , goal8 , o );
   const SlidingTilePDB< 3 , 3 , Mod3Entries > d( goal8 , r3 );
   }
  catch( std::invalid_argument & ) {
   thrown = true;
   }
  check( thrown , "mod 3" + tag + ": the domain refuses the reflection" );
  }

 // the reflection: an automorphism of the puzzle that fixes the goal
 {
  const AdditivePDB< 3 , 3 , ByteEntries > p( Partition{ { 1 , 2 , 3 } ,
                                                         { 4 , 5 , 6 } ,
                                                         { 7 , 8 } } );
  const std::uint64_t goal = EightPuzzle::pack( goal8 );
  bool fixed = true , same = true , lookup = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   std::array< std::uint8_t , 9 > r{};
   p.reflection( q.data() , r.data() );
   EightPuzzle::board_type b{};
   for( int t = 0 ; t < 9 ; ++t )
    b[ r[ t ] ] = std::uint8_t( t );
   const std::uint64_t phi = EightPuzzle::pack( b );
   same = same && ( dist.at( phi ) == d );
   lookup = lookup && ( p.heuristic( q.data() ) ==
                        std::max( p.sum( q.data() ) , p.sum( r.data() ) ) );
   if( s == goal )
    fixed = phi == goal;
   }
  check( fixed && same && lookup , "the reflection fixes the goal, keeps "
         "the distances and gives the maximum of the two sums" );
  bool none = true;
  auto off = goal8;                     // the blank off the diagonal
  std::swap( off[ 0 ] , off[ 1 ] );
  none = ! AdditivePDB< 3 , 3 , ByteEntries >( Partition{ { 2 , 3 } } , off )
          .reflected();
  check( none , "no reflection when the blank of the goal is off the "
                "diagonal" );
  }

 // bad partitions
 bool thrown = false;
 try {
  const AdditivePDB< 3 , 3 , ByteEntries > p( Partition{ { 1 , 2 } ,
                                                         { 2 , 3 } } );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "a tile in two patterns is refused" );
 thrown = false;
 try {
  const TilePattern< 3 , 3 > p( { 0 , 1 } );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "the blank in a pattern is refused" );

 // the 15-puzzle, 5-5-5: admissible on the instances of Korf; the tables
 // in bytes are built once, in a cache, and encoded from there
 namespace fs = std::filesystem;
 const fs::path dir = fs::temp_directory_path() /
                      ( "pbefs-pdb-test-" + std::to_string(
                         std::random_device()() ) );
 fs::create_directories( dir );
 const auto inst = korf( data );
 for( BlankRule rule : { BlankRule::ignored , BlankRule::minimized } ) {
  PDBOptions o;
  o.blank = rule;
  o.cache = ( dir / "cache15" ).string();
  const AdditivePDB< 4 , 4 , ByteEntries > p( partition_15_5_5_5() ,
                                              FifteenPuzzle::korf_goal() ,
                                              o );
  const AdditivePDB< 4 , 4 , NibbleEntries > n( partition_15_5_5_5() ,
                                                FifteenPuzzle::korf_goal() ,
                                                o );
  o.compress = 8;
  const AdditivePDB< 4 , 4 , MinNibbleEntries > m( partition_15_5_5_5() ,
                                                   FifteenPuzzle::korf_goal() ,
                                                   o );
  bool ok = true;
  double sum = 0 , md = 0;
  for( const auto & ki : inst ) {
   const SlidingTilePDB< 4 , 4 , ByteEntries > d( ki.board , p );
   const FifteenPuzzle f( ki.board );
   std::array< std::uint8_t , 16 > q{};
   d.cells_of( d.initial_state() , q.data() );
   const int h = d.heuristic( d.initial_state() );
   ok = ok && ( h <= ki.cost ) && ( h >= f.heuristic( f.initial_state() ) )
        && ( n.heuristic( q.data() ) == h ) &&
        ( m.heuristic( q.data() ) <= h );
   sum += h;
   md += f.heuristic( f.initial_state() );
   }
  const std::string r = rule == BlankRule::ignored ? "ignored" : "minimized";
  check( ok , "5-5-5, blank " + r + ": Manhattan <= min 8 <= nibble = byte "
              "<= the optimum on the 100 instances of Korf" );
  std::cout << "  5-5-5, blank " << r << ": mean h " << sum / 100
            << " (Manhattan " << md / 100 << ")" << std::endl;
  std::array< std::uint8_t , 16 > g{};
  SlidingTilePDB< 4 , 4 , ByteEntries >::cells_of(
   FifteenPuzzle::pack( FifteenPuzzle::korf_goal() ) , g.data() );
  check( ( p.heuristic( g.data() ) == 0 ) && ( m.heuristic( g.data() ) == 0 ) ,
         "5-5-5, blank " + r + ": h = 0 at the goal" );
  }

 // the files: round trip, corruption, cache
 const TilePattern< 3 , 3 > pat( { 1 , 2 , 3 , 4 } );
 const auto h = PatternTable< 3 , 3 , ByteEntries >::distances(
                 pat , BlankRule::ignored );
 const auto round = [ & ]< class Enc >( const std::string & what ,
                                        unsigned c ) {
  using T = PatternTable< 3 , 3 , Enc >;
  const T t = T::encode( pat , BlankRule::ignored , h , c );
  const std::string f = ( dir / ( what + ".pdb" ) ).string();
  t.save( f );
  const T u = T::load( f , pat , BlankRule::ignored , c );
  check( u.from_file() && ! t.from_file() &&
         ( u.entries().payload() == t.entries().payload() ) ,
         what + ": saved and loaded" );
  };
 round.operator()< ByteEntries >( "byte" , 1 );
 round.operator()< NibbleEntries >( "nibble" , 1 );
 round.operator()< Mod3Entries >( "mod3" , 1 );
 round.operator()< ParityBitEntries >( "bit" , 1 );
 round.operator()< MinNibbleEntries >( "min" , 4 );
 round.operator()< SparseEntries >( "sparse" , 1 );
 {
  PDBOptions o;
  const Partition half{ { 1 , 2 , 3 , 4 } , { 5 , 6 , 7 , 8 } };
  const AdditivePDB< 3 , 3 , NibbleEntries > n( half , goal8 , o );
  const AdditivePDB< 3 , 3 , SparseEntries > sp( half , goal8 , o );
  bool same = true;
  for( const auto & [ s , d ] : dist ) {
   const auto q = cells8( s );
   same = same && ( n.heuristic( q.data() ) == sp.heuristic( q.data() ) );
   }
  check( same , "sparse = nibble on all the states of the 8-puzzle" );
  }

 using T = PatternTable< 3 , 3 , NibbleEntries >;
 const std::string f = ( dir / "nibble.pdb" ).string();
 const auto rejected = [ & ]( const std::string & what , auto && change ) {
  std::ifstream in( f , std::ios::binary );
  std::vector< char > b( ( std::istreambuf_iterator< char >( in ) ) ,
                         std::istreambuf_iterator< char >() );
  change( b );
  const std::string g = ( dir / "changed.pdb" ).string();
  std::ofstream( g , std::ios::binary ).write( b.data() ,
                                               std::streamsize( b.size() ) );
  bool bad = false;
  try {
   static_cast< void >( T::load( g , pat , BlankRule::ignored ) );
   }
  catch( std::runtime_error & ) {
   bad = true;
   }
  check( bad , "a file " + what + " is rejected" );
  };
 rejected( "with a byte changed" , []( auto & b ) {
  b[ b.size() / 2 ] ^= 1;
  } );
 rejected( "truncated" , []( auto & b ) { b.resize( b.size() - 9 ); } );
 rejected( "empty" , []( auto & b ) { b.clear(); } );
 rejected( "of another version" , []( auto & b ) { b[ 8 ] = 2; } );
 bool other = false;
 try {
  static_cast< void >( T::load( f , TilePattern< 3 , 3 >( { 1 , 2 , 3 , 5 } ) ,
                                BlankRule::ignored ) );
  }
 catch( std::runtime_error & ) {
  other = true;
  }
 check( other , "a file of another pattern is rejected" );

 const std::string cache = ( dir / "cache" ).string();
 const auto first = T::cached( cache , pat , BlankRule::ignored );
 const auto second = T::cached( cache , pat , BlankRule::ignored );
 check( ! first.from_file() && second.from_file() &&
        ( first.entries().payload() == second.entries().payload() ) ,
        "the cache builds a table once, then loads it" );
 const std::string own = ( fs::path( cache ) /
                           T::file_name( pat , BlankRule::ignored ) ).string();
 {
  std::fstream io( own , std::ios::binary | std::ios::in | std::ios::out );
  io.seekp( 60 );
  io.put( 'x' );
  }
 const auto third = T::cached( cache , pat , BlankRule::ignored );
 check( ! third.from_file() &&
        ( third.entries().payload() == first.entries().payload() ) ,
        "the cache encodes again a corrupted table" );
 fs::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A*, HDA* and SPA* with the pattern databases: the optimal costs, and
// fewer expansions than with the Manhattan distance

void test_pdb_search( const std::string & data )
{
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 const Partition part{ { 1 , 2 , 3 , 4 } , { 5 , 6 , 7 , 8 } };
 const auto goal8 = EightPuzzle::korf_goal();
 PDBOptions o;
 const AdditivePDB< 3 , 3 , NibbleEntries > nib( part , goal8 , o );
 o.compress = 4;
 const AdditivePDB< 3 , 3 , MinNibbleEntries > mn( part , goal8 , o );
 o.compress = 1;
 o.blank = BlankRule::minimized;
 const AdditivePDB< 3 , 3 , NibbleEntries > mz( part , goal8 , o );
 o.blank = BlankRule::ignored;
 o.reflect = false;
 const AdditivePDB< 3 , 3 , ParityBitEntries > bit( part , goal8 , o );
 std::size_t k = 0 , em = 0 , ep = 0;
 for( const auto & [ s , d ] : dist ) {
  if( ( k++ % 907 ) != 0 )                     // 200 states, spread out
   continue;
  const auto b = EightPuzzle::unpack( s );
  const EightPuzzle m( b );
  const auto sm = astar( m );
  const SlidingTilePDB< 3 , 3 , NibbleEntries > dn( b , nib );
  const auto sn = astar( dn );
  const SlidingTilePDB< 3 , 3 , MinNibbleEntries > dmn( b , mn );
  const SlidingTilePDB< 3 , 3 , NibbleEntries > dmz( b , mz );
  const SlidingTilePDB< 3 , 3 , ParityBitEntries > db( b , bit );
  const auto s1 = astar( dmn );
  const auto s2 = astar( dmz );
  const auto s3 = astar( db );
  if( ! ( sn.solved() && ( sn.cost == d ) && ( s1.cost == d ) &&
          ( s2.cost == d ) && ( s3.cost == d ) ) ) {
   check_solution( dn , sn , d , "A* nibble, 8-puzzle" );
   check_solution( dmn , s1 , d , "A* min 4, 8-puzzle" );
   check_solution( dmz , s2 , d , "A* minimized, 8-puzzle" );
   check_solution( db , s3 , d , "A* parity bit, 8-puzzle" );
   }
  em += sm.stats.expanded;
  ep += sn.stats.expanded;
  }
 check( failures == 0 , "A* with the tables on 200 states of the 8-puzzle" );
 check( ep < em , "A* with 4-4 expands fewer nodes than with Manhattan (" +
                  std::to_string( ep ) + " against " + std::to_string( em ) +
                  ")" );

 const auto inst = korf( data );
 const AdditivePDB< 4 , 4 , NibbleEntries > p( partition_15_5_5_5() );
 PDBOptions nr;
 nr.reflect = false;
 const AdditivePDB< 4 , 4 , ParityBitEntries > pb( partition_15_5_5_5() ,
                                                  FifteenPuzzle::korf_goal() ,
                                                  nr );
 for( int id : easy_korf ) {
  const auto & ki = inst[ id - 1 ];
  const FifteenPuzzle m( ki.board );
  const SlidingTilePDB< 4 , 4 , NibbleEntries > d( ki.board , p );
  const SlidingTilePDB< 4 , 4 , ParityBitEntries > db( ki.board , pb );
  const auto sm = astar( m );
  const auto sp = astar( d );
  const auto sh = astar( d , OpenKind::heap );
  const auto sb = astar( db );
  const std::string tag = ", Korf #" + std::to_string( id );
  check_solution( d , sp , ki.cost , "A* 5-5-5" + tag );
  check_solution( d , sh , ki.cost , "A* heap 5-5-5" + tag );
  check_solution( db , sb , ki.cost , "A* 5-5-5 parity bit" + tag );
  check( sp.stats.expanded < sm.stats.expanded , "A* 5-5-5 expands fewer "
         "nodes than Manhattan" + tag + " (" +
         std::to_string( sp.stats.expanded ) + " against " +
         std::to_string( sm.stats.expanded ) + ")" );
  }

 using D = SlidingTilePDB< 4 , 4 , NibbleEntries >;
 for( unsigned t : { 1u , 3u , 8u } ) {
  ParallelOptions po;
  po.threads = t;
  for( int id : { easy_korf[ 0 ] , easy_korf[ 1 ] } ) {
   const auto & ki = inst[ id - 1 ];
   const D d( ki.board , p );
   const std::string tag = " p " + std::to_string( t ) + ", Korf #" +
                           std::to_string( id );
   const ZobristDistribution< D > z( d );
   check_solution( d , HDAStar< D , decltype( z ) , ADS::ThreadPool >(
                        d , z , po ).solve() , ki.cost , "HDA*[Z] 5-5-5" +
                                                         tag );
   check_solution( d , SPAStar< D , ADS::ThreadPool >( d , po ).solve() ,
                   ki.cost , "SPA* 5-5-5" + tag );
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File pdb.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
