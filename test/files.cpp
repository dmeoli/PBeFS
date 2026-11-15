/*--------------------------------------------------------------------------*/
/*----------------------------- File files.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the files of the
 * data: the alignment of sequences, the boards, the other files and the
 * reopenings.
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
// the cost of an optimal alignment by the dynamic program on the whole
// lattice, written from the definition of the sum-of-pairs cost with no
// part of SequenceAlignment but its pair costs

static int alignment_dp( const std::vector< std::string > & seq ,
                         const SequenceAlignment<> & d )
{
 const int n = int( seq.size() );
 std::vector< std::size_t > w( n );
 std::size_t size = 1;
 for( int i = n ; i-- > 0 ; ) {
  w[ i ] = size;
  size *= seq[ i ].size() + 1;
  }
 std::vector< int > best( size , std::numeric_limits< int >::max() );
 best[ 0 ] = 0;
 std::vector< std::size_t > p( n );
 for( std::size_t k = 0 ; k < size ; ++k ) {
  std::size_t r = k;
  for( int i = 0 ; i < n ; ++i ) {
   p[ i ] = r / w[ i ];
   r %= w[ i ];
   }
  for( std::uint32_t m = 1 ; m < ( 1u << n ) ; ++m ) {
   bool ok = true;
   std::size_t from = k;
   for( int i = 0 ; i < n ; ++i )
    if( m & ( 1u << i ) ) {
     ok = ok && ( p[ i ] > 0 );
     from -= w[ i ];
     }
   if( ! ok || ( best[ from ] == std::numeric_limits< int >::max() ) )
    continue;
   int c = 0;
   for( int i = 0 ; i < n ; ++i )
    for( int j = i + 1 ; j < n ; ++j ) {
     const int a = m & ( 1u << i ) ? SequenceAlignment<>::symbol(
                                      seq[ i ][ p[ i ] - 1 ] ) : amino_gap;
     const int b = m & ( 1u << j ) ? SequenceAlignment<>::symbol(
                                      seq[ j ][ p[ j ] - 1 ] ) : amino_gap;
     c += d.pair_cost( a , b );
     }
   best[ k ] = std::min( best[ k ] , best[ from ] + c );
   }
  }
 return( best[ size - 1 ] );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// multiple sequence alignment: the reader of FASTA files, the costs of
// PAM250, A* against the dynamic program on the lattice, Dijkstra and the
// bounds on its expansions, the alignment of the path, IDA* and DDD-A*,
// the features of the Zobrist distributions and the graph of the nblocks

void test_msa( const std::string & data )
{
 using D = SequenceAlignment<>;
 static_assert( IncrementalFeatureDomain< D > && AbstractFeatureDomain< D >
                && AbstractStateDomain< D > && NBlockDomain< D > &&
                DiskDomain< D > && ! PairwiseDomain< D > );

 std::istringstream fa( "junk\n>a b\nac-g\n d\n>c\n\nWW.\n" );
 const auto rf = D::read_fasta( fa );
 check( ( rf.size() == 2 ) && ( rf[ 0 ] == "ACGD" ) && ( rf[ 1 ] == "WW" ) ,
        "FASTA: names, lines, case, gaps and blank lines" );
 bool threw = false;
 try {
  std::istringstream none( "ACGT\n" );
  ( void ) D::read_fasta( none );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "FASTA without sequences refused" );
 int refused = 0;
 for( const auto & bad : std::vector< std::vector< std::string > >{
       { "ACD" } , { "AC" , "" } , { "AC" , "AOC" } , { "AC" , "A*" } ,
       std::vector< std::string >( 11 , "A" ) } )
  try {
   const D x( bad );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 5 , "a single sequence, an empty one, an unknown symbol,"
                       " a gap and 11 sequences are refused" );

 const D two( { "W" , "W" } );
 const auto sym = []( char a ) { return( D::symbol( a ) ); };
 check( ( two.pair_cost( sym( 'W' ) , sym( 'W' ) ) == 0 ) &&
        ( two.pair_cost( sym( 'C' ) , sym( 'C' ) ) == 5 ) &&
        ( two.pair_cost( sym( 'A' ) , amino_gap ) == 25 ) &&
        ( two.pair_cost( amino_gap , amino_gap ) == 0 ) &&
        ( two.pair_cost( sym( 'W' ) , sym( 'C' ) ) == 25 ) ,
        "the costs 17 - PAM250, 25 against a gap, 0 for two gaps" );
 bool symmetric = true;
 for( int a = 0 ; a < amino_count ; ++a )
  for( int b = 0 ; b < amino_count ; ++b )
   symmetric = symmetric && ( pam250[ a ][ b ] == pam250[ b ][ a ] );
 check( symmetric , "PAM250 is symmetric" );
 // a pair aligned with the linear gap penalty 8: the cost against the
 // largest score sum( S ) - 8 g, by the recurrence of Needleman and Wunsch
 {
  const std::string u = "HEAGAWGHEE" , v = "PAWHEAE";
  std::vector< std::vector< int > > best( u.size() + 1 ,
                                          std::vector< int >( v.size() + 1 ) );
  for( std::size_t i = 0 ; i <= u.size() ; ++i )
   for( std::size_t j = 0 ; j <= v.size() ; ++j )
    if( i + j == 0 )
     best[ i ][ j ] = 0;
    else {
     int b = std::numeric_limits< int >::min();
     if( i )
      b = std::max( b , best[ i - 1 ][ j ] - 8 );
     if( j )
      b = std::max( b , best[ i ][ j - 1 ] - 8 );
     if( i && j )
      b = std::max( b , best[ i - 1 ][ j - 1 ] +
                        pam250[ sym( u[ i - 1 ] ) ][ sym( v[ j - 1 ] ) ] );
     best[ i ][ j ] = b;
     }
  const D p( { u , v } , pam250 , 8 , 16 , 8 );
  const int score = best[ u.size() ][ v.size() ];
  check_solution( p , astar( p ) , 17 * int( u.size() + v.size() ) -
                  2 * score , "the gap penalty 8: twice the cost is 17 ( L1 "
                  "+ L2 ) less twice the score " + std::to_string( score ) );
  const D q( { u , v } );
  int score16 = 0;                     // the same with the penalty 16.5
  std::vector< std::vector< double > > b2( u.size() + 1 ,
                                           std::vector< double >(
                                            v.size() + 1 ) );
  for( std::size_t i = 0 ; i <= u.size() ; ++i )
   for( std::size_t j = 0 ; j <= v.size() ; ++j )
    if( i + j > 0 ) {
     double b = -1e18;
     if( i )
      b = std::max( b , b2[ i - 1 ][ j ] - 16.5 );
     if( j )
      b = std::max( b , b2[ i ][ j - 1 ] - 16.5 );
     if( i && j )
      b = std::max( b , b2[ i - 1 ][ j - 1 ] +
                        pam250[ sym( u[ i - 1 ] ) ][ sym( v[ j - 1 ] ) ] );
     b2[ i ][ j ] = b;
     }
  score16 = int( std::lround( 2 * b2[ u.size() ][ v.size() ] ) );
  check_solution( q , astar( q ) , ( 17 * int( u.size() + v.size() ) -
                                     score16 ) / 2 ,
                  "the gap costs of PAM250: the score with the penalty 16.5" );
  }
 std::ostringstream os;
 for( const auto & row : pam250 )
  for( int x : row )
   os << x << ' ';
 std::istringstream is( os.str() );
 check( D::read_scores( is ) == pam250 , "a matrix read back" );

 // the costs of the instances, against the dynamic program on the lattice
 // and the costs that the code of HDA* of Jinnai and Fukunaga finds
 struct Inst { const char * file; int cost; };
 const Inst inst[] = { { "test0" , 1741 } , { "test1" , 2665 } ,
                       { "test3" , 2022 } , { "test4" , 2312 } ,
                       { "test5" , 3990 } };
 for( const auto & [ name , cost ] : inst ) {
  const auto seq = fasta( data + "/msa/easy/" + name + ".tfa" );
  const D d( seq );
  const std::string what = std::string( "MSA " ) + name;
  const auto sol = astar( d );
  if( seq.size() <= 4 ) {      // test5: 483840 points of 127 moves each
   check( alignment_dp( seq , d ) == cost , what + ": the dynamic program "
          "on the lattice costs " + std::to_string( cost ) );
   check_consistent_astar( d , all_distances( d ) , sol , what );
   }
  else
   check_solution( d , sol , cost , what + ", the cost of the reference "
                   "code" );
  const auto rows = d.alignment( sol.path );
  bool same = rows.size() == seq.size();
  int columns = 0;
  for( std::size_t c = 0 ; same && ( c < rows[ 0 ].size() ) ; ++c ) {
   for( std::size_t i = 0 ; i < rows.size() ; ++i )
    for( std::size_t j = i + 1 ; j < rows.size() ; ++j )
     columns += d.pair_cost( rows[ i ][ c ] == '-' ? amino_gap
                                                   : sym( rows[ i ][ c ] ) ,
                             rows[ j ][ c ] == '-' ? amino_gap
                                                   : sym( rows[ j ][ c ] ) );
   }
  for( std::size_t i = 0 ; same && ( i < rows.size() ) ; ++i ) {
   std::string r = rows[ i ];
   r.erase( std::remove( r.begin() , r.end() , '-' ) , r.end() );
   same = ( r == seq[ i ] ) && ( rows[ i ].size() == rows[ 0 ].size() );
   }
  check( same && ( columns == cost ) , what + ": the rows of the alignment "
         "are the sequences and their columns cost " +
         std::to_string( columns ) );
  }

 // the instance of 7 sequences, IDA* and DDD-A*
 {
  const auto seq = fasta( data + "/msa/easy/test2.tfa" );
  const D d( seq );
  check_solution( d , astar( d ) , 6055 , "MSA test2, 7 sequences" );
  const auto seq0 = fasta( data + "/msa/easy/test0.tfa" );
  const D d0( seq0 );
  check_solution( d0 , idastar( d0 ) , 1741 , "IDA* MSA test0" );
  ADS::TempDirectory dir;
  DDDOptions o;
  o.directory = dir.path();
  check_solution( d0 , ddd_astar( d0 , {} , o ) , 1741 , "DDD-A* MSA test0" );
  }

 // the instances of BAliBASE read and their heuristics computed
 int read = 0;
 for( const auto & f : std::filesystem::directory_iterator(
                        data + "/msa/balibase" ) ) {
  const auto seq = fasta( f.path().string() );
  const D d( seq );
  read += ( seq.size() >= 4 ) && ( d.heuristic( d.initial_state() ) > 0 );
  }
 check( read == 15 , "the 15 instances of BAliBASE, 4 or 5 sequences each" );

 // the Zobrist keys of the successors from the key of the parent, and the
 // graph of the nblocks: symmetric, and holding every move between boxes
 const auto seq = fasta( data + "/msa/easy/test3.tfa" );
 const D d( seq , pam250 , 3 , 4 );
 const ZobristDistribution< D > z( d );
 const AbstractZobristDistribution< D > az( d );
 static_assert( IncrementalDistribution< decltype( z ) , D > &&
                IncrementalDistribution< decltype( az ) , D > );
 const auto g = all_distances( d );
 std::set< std::pair< std::size_t , std::size_t > > nb;
 for( std::size_t b = 0 ; b < d.num_nblocks() ; ++b )
  d.for_each_nblock_successor( b , [ & ]( std::size_t c ) {
   nb.emplace( b , c );
   } );
 int bad = 0 , moves = 0 , kz = 0 , kaz = 0 , ka = 0 , visited = 0;
 const AbstractStateDistribution< D > as( d );
 for( const auto & sg : g ) {
  const auto & s = sg.first;
  if( ( visited++ % 5 ) != 0 )
   continue;
  for( unsigned p : { 1u , 3u , 64u } ) {
   const auto oz = successor_owners< D >( z , s , z.key( s ) , p );
   const auto oaz = successor_owners< D >( az , s , az.key( s ) , p );
   d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
    if( ( oz( t ) != std::pair( z.owner( t , p ) , z.key( t ) ) ) ||
        ( oaz( t ) != std::pair( az.owner( t , p ) , az.key( t ) ) ) )
     ++bad;
    } );
   }
  d.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
   ++moves;
   kz += z.owner( t , 64 ) == z.owner( s , 64 );
   kaz += az.owner( t , 64 ) == az.owner( s , 64 );
   ka += as.owner( t , 64 ) == as.owner( s , 64 );
   const auto a = d.nblock( s ) , b = d.nblock( t );
   if( ( a != b ) && ( ! nb.count( { a , b } ) || ! nb.count( { b , a } ) ) )
    ++bad;
   } );
  }
 bool sym_nb = true;
 for( const auto & [ a , b ] : nb )
  sym_nb = sym_nb && nb.count( { b , a } ) && ( a != b ) &&
           ( a < d.num_nblocks() ) && ( b < d.num_nblocks() );
 check( ( bad == 0 ) && sym_nb , "MSA: the keys of the successors from the "
        "key of the parent, and the nblocks, symmetric and next to each "
        "other along every move" );
 // the hyperplanes: the levels, and the bound of Remark 1 of Kobayashi,
 // Kishimoto and Watanabe on the workers that the successors of the states
 // of a worker go to, ceil( n / d ) + 1 for d >= 1 (n / d + 1 there, which
 // the thickness 5 exceeds) and ( n + 1 ) / d for d < 1
 static_assert( LevelDomain< D > );
 bool levels = d.max_level() == 19 + 19 + 20 + 15;
 for( const auto & sg : g )
  d.for_each_successor( sg.first , 0 , [ & ]( const auto & t , int , int ) {
   const auto & s = sg.first;
   levels = levels && ( d.level( t ) > d.level( s ) ) &&
            ( d.level( t ) <= d.level( s ) + 4 );
   } );
 check( levels , "MSA: every move raises the level by 1 to 4, up to 73" );
 const unsigned n = unsigned( seq.size() );
 for( const auto [ planes , split ] : { std::pair( 1u , 1u ) ,
                                        std::pair( 2u , 1u ) ,
                                        std::pair( 5u , 1u ) ,
                                        std::pair( 1u , 2u ) ,
                                        std::pair( 1u , 4u ) } ) {
  const HyperplaneDistribution< D > h( d , 3 , planes , split );
  std::vector< std::set< unsigned > > to( 64 );
  for( const auto & sg : g )
   d.for_each_successor( sg.first , 0 , [ & ]( const auto & t , int ,
                                               int ) {
    to[ h.owner( sg.first , 64 ) ].insert( h.owner( t , 64 ) );
    } );
  std::size_t most = 0;
  for( const auto & w : to )
   most = std::max( most , w.size() );
  const std::size_t bound = split == 1 ? ( n + planes - 1 ) / planes + 1
                                        : ( n + 1 ) * split;
  check( most <= bound , "MSA, 64 workers: the successors of the states of "
         "a worker reach at most " + std::to_string( most ) + " workers "
         "with the planes of thickness " + std::to_string( planes ) + "/" +
         std::to_string( split ) + ", against the bound " +
         std::to_string( bound ) );
  }
 bool refused_hwd = false;
 try {
  const HyperplaneDistribution< D > bad( d , 1 , 2 , 2 );
  }
 catch( const std::invalid_argument & ) {
  refused_hwd = true;
  }
 const auto th = HyperplaneDistribution< D >::thickness( d , 64 );
 check( refused_hwd && ( th == std::pair( 1u , 27u ) ) &&
        ( HyperplaneDistribution< D >::thickness( d , 1 ) ==
          std::pair( 1u , 1u ) ) , "the thickness 1/27 of the rule of "
        "Kobayashi et al. for 73 residues and 64 workers, and 2/2 refused" );

 // the random instances of Jinnai and Fukunaga
 bool rnd = true;
 for( std::uint64_t seed = 1 ; seed <= 20 ; ++seed ) {
  const auto r = D::random_sequences( seed );
  const int k = int( r.size() );
  rnd = rnd && ( k >= 4 ) && ( k <= 9 ) &&
        ( D::random_sequences( seed ) == r );
  for( const auto & x : r ) {
   rnd = rnd && ( 10 * k * int( x.size() ) > 9 * 5000 ) &&
         ( 10 * k * int( x.size() ) < 11 * 5000 );
   for( char c : x )
    rnd = rnd && ( D::symbol( c ) >= 0 ) && ( D::symbol( c ) < 20 );
   }
  }
 const auto small = D::random_sequences( 7 , 3 , 30 );
 const D ds( small );
 check( rnd && ( small.size() == 3 ) && ( small[ 0 ].size() == 10 ) &&
        ( alignment_dp( small , ds ) == astar( ds ).cost ) ,
        "random sequences: 4 to 9 of lengths within 10% of 5000 / n, of "
        "the 20 amino acids, the same for the same seed" );
 check( ( ka > 4 * kz ) && ( kaz > 4 * kz ) , "MSA, 64 workers: of " +
        std::to_string( moves ) + " moves, the boxes of side 3 keep " +
        std::to_string( ka ) + " on the owner of their parent with AHDA*, "
        + std::to_string( kaz ) + " with AZH, and Zobrist hashing " +
        std::to_string( kz ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the grids in the format of Burns et al.: the reader and its errors, the
// life costs and the heuristics, exact on the boards without blocked cells,
// A* against Dijkstra and the bounds on its expansions, the corners cut,
// the boards of the data with both origins of the rows

void test_board( const std::string & data )
{
 std::istringstream tiny( "3 4\nBoard:\n  # \n#\n    \nLife\nEight-way\n"
                          "0 0 3 2\n" );
 const GridBoard t = GridBoard::read( tiny );
 check( ( t.map.width() == 4 ) && ( t.map.height() == 3 ) &&
        ( t.map.cost( 2 , 0 ) == 0 ) && ( t.map.cost( 0 , 1 ) == 0 ) &&
        ( t.map.cost( 1 , 1 ) == 1 ) && ( t.sx == 0 ) && ( t.sy == 2 ) &&
        ( t.gx == 3 ) && ( t.gy == 0 ) && ( t.moves == 8 ) &&
        ( t.cost == GridCost::life ) , "a board read with y from the bottom"
        ", short rows free" );
 std::istringstream top( "2 2\nBoard:\n# \n  \nUnit\nFour-way\n0 0 1 1\n" );
 const GridBoard tb = GridBoard::read( top , false );
 check( ( tb.sy == 0 ) && ( tb.gy == 1 ) && ( tb.map.cost( 0 , 0 ) == 1 ) &&
        ( tb.moves == 4 ) && ( tb.cost == GridCost::cells ) ,
        "a board read with y from the top, the blocked start freed" );
 int refused = 0;
 for( const char * bad : { "2 2\nBoard\n  \n  \nUnit\nFour-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nOther\nFour-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nUnit\nSix-way\n0 0 1 1\n" ,
                           "2 2\nBoard:\n  \n  \nUnit\nFour-way\n0 0 2 1\n" ,
                           "2 2\nBoard:\n  \n" } )
  try {
   std::istringstream in( bad );
   ( void ) GridBoard::read( in );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 5 , "boards without Board:, with other costs or moves,"
                       " a goal outside or missing rows are refused" );

 // the life heuristics: the exact cost on a board without blocked cells,
 // the moves to a corner costing as the straight ones (relaxed by hand)
 int wrong = 0;
 for( int moves : { 4 , 8 } )
  for( int gy = 0 ; gy < 7 ; gy += 3 )
   for( int gx = 0 ; gx < 9 ; gx += 4 ) {
    GridBoard e;
    e.map = GridMap( 9 , 7 );
    e.sx = 0;
    e.sy = 0;
    e.gx = gx;
    e.gy = gy;
    e.moves = moves;
    e.cost = GridCost::life;
    std::vector< std::int64_t > dist( 63 , std::int64_t( 1 ) << 40 );
    dist[ gy * 9 + gx ] = 0;
    for( int it = 0 ; it < 63 ; ++it )        // Bellman-Ford to the goal
     for( int u = 0 ; u < 63 ; ++u )
      for( int dy = -1 ; dy <= 1 ; ++dy )
       for( int dx = -1 ; dx <= 1 ; ++dx ) {
        const int x = u % 9 + dx , y = u / 9 + dy;
        if( ( ! dx && ! dy ) || ( ( moves == 4 ) && dx && dy ) ||
            ( x < 0 ) || ( y < 0 ) || ( x >= 9 ) || ( y >= 7 ) )
         continue;
        dist[ u ] = std::min( dist[ u ] ,
                              std::int64_t( u / 9 ) + dist[ y * 9 + x ] );
        }
    for( int u = 0 ; u < 63 ; ++u ) {
     double h;
     if( moves == 4 )
      h = Grid< 4 >( e ).heuristic( std::uint32_t( u ) );
     else
      h = Grid< 8 >( e ).heuristic( std::uint32_t( u ) );
     wrong += std::abs( h - double( dist[ u ] ) ) > 1e-9;
     }
    }
 check( wrong == 0 , "the life heuristics are the exact costs of the "
                     "relaxation on boards without blocked cells" );

 // random boards of both moves and costs, against Dijkstra
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed )
  for( GridCost cost : { GridCost::cells , GridCost::life } ) {
   const std::string tag = std::string( cost == GridCost::life ? " life"
                                                               : " unit" ) +
                           ", seed " + std::to_string( seed );
   const auto b4 = random_board( 60 , 40 , 0.3 , seed , 4 , cost );
   const Grid< 4 > g4( b4 , 10 , 8 );
   const double r4 = board_dijkstra( b4 );
   if( r4 >= 0 ) {
    check_solution( g4 , astar( g4 ) , int( std::lround( r4 ) ) ,
                    "board 4" + tag + " against Dijkstra on the map" );
    check_consistent_astar( g4 , all_distances( g4 ) , astar( g4 ) ,
                            "board 4" + tag );
    }
   else
    check( ! astar( g4 ).solved() , "board 4 unsolvable" + tag );
   const auto b8 = random_board( 60 , 40 , 0.3 , seed , 8 , cost );
   const Grid< 8 > g8( b8 , 10 , 8 );
   const double r8 = board_dijkstra( b8 );
   if( r8 >= 0 ) {
    check_solution( g8 , astar( g8 ) , r8 , "board 8" + tag +
                    " against Dijkstra on the map" );
    check_consistent_astar( g8 , all_distances( g8 ) , astar( g8 ) ,
                            "board 8" + tag );
    }
   else
    check( ! astar( g8 ).solved() , "board 8 unsolvable" + tag );
   }

 // a corner cut between two blocked cells
 GridBoard c;
 c.map = GridMap( 2 , 2 );
 c.map.set_cost( 1 , 0 , 0 );
 c.map.set_cost( 0 , 1 , 0 );
 c.gx = 1;
 c.gy = 1;
 c.moves = 8;
 const Grid< 8 > cut( c );
 check_solution( cut , astar( cut ) , grid_sqrt2 , "a corner cut between "
                 "two blocked cells" );
 const Grid< 8 > nocut( c.map , 0 , 0 , 1 , 1 );
 check( ! astar( nocut ).solved() , "no corner cut by default" );
 bool threw = false;
 try {
  const Grid< 4 > wrong4( c );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "a board of 8 moves refused by Grid< 4 >" );

 // the boards of the data, their rows from the bottom or from the top,
 // against Dijkstra and the costs found by the code of their authors
 struct File { const char * name; bool bottom; int cost; };
 const File files[] = { { "ethan_unit_four-way_1" , true , 23 } ,
                        { "jordan_unit_four-way_100x100_1" , true , 143 } ,
                        { "jordan_unit_four-way_500x500_1" , true , 763 } ,
                        { "jordan_life_four-way_500x500_1" , true , 362873 } ,
                        { "grid_500-500-0.7_0" , false , 755 } };
 for( const auto & [ name , bottom , cost ] : files ) {
  std::ifstream in( data + "/grid/" + name );
  const GridBoard b = GridBoard::read( in , bottom );
  const Grid< 4 > g( b );
  check( std::lround( board_dijkstra( b ) ) == cost , std::string( "board " )
         + name + ": Dijkstra costs " + std::to_string( cost ) );
  check_solution( g , astar( g ) , cost , std::string( "board " ) + name );
  }
 }

/*--------------------------------------------------------------------------*/
// the other files of the data: the 15-puzzle by the cells of the tiles,
// the 24-puzzle sets, and the cities of the traveling salesperson

void test_files( const std::string & data )
{
 const auto inst = korf( data );
 for( int id : { 1 , 12 , 55 } ) {
  std::ifstream in( data + "/15-puzzle/korf_" + std::to_string( id ) );
  const auto [ start , goal ] = FifteenPuzzle::read_positions( in );
  check( ( start == inst[ id - 1 ].board ) &&
         ( goal == FifteenPuzzle::korf_goal() ) , "korf_" +
         std::to_string( id ) + " is the instance " + std::to_string( id ) +
         " of Korf" );
  if( id != 1 ) {
   const FifteenPuzzle d( start , goal );
   check_solution( d , astar( d ) , inst[ id - 1 ].cost , "A* on korf_" +
                   std::to_string( id ) );
   }
  }
 std::istringstream small( "3 3\nstarting positions for each tile:\n"
                           "1 0 2 3 4 5 6 7 8\ngoal positions:\n"
                           "0 1 2 3 4 5 6 7 8\n" );
 const auto [ s8 , g8 ] = EightPuzzle::read_positions( small );
 check( ( s8[ 0 ] == 1 ) && ( s8[ 1 ] == 0 ) && ( g8 == EightPuzzle::
        korf_goal() ) , "an 8-puzzle by the cells of its tiles" );
 int refused = 0;
 for( const char * bad : { "4 4\nstarting positions for each tile:\n" ,
                           "3 3\nstarting positions for each tile:\n"
                           "0 0 1 2 3 4 5 6 7\n" } )
  try {
   std::istringstream in( bad );
   ( void ) EightPuzzle::read_positions( in );
   }
  catch( const std::invalid_argument & ) {
   ++refused;
   }
 check( refused == 2 , "a board of another size and a repeated cell are "
                       "refused" );

 for( const auto & [ file , count ] : { std::pair( "instances100.txt" , 100 ) ,
                                        std::pair( "fromkorf8.txt" , 8 ) } ) {
  std::ifstream in( data + "/24-puzzle/" + file );
  std::string line;
  int n = 0 , bad = 0;
  while( std::getline( in , line ) ) {
   if( line.empty() || ( line[ 0 ] == '#' ) )
    continue;
   try {
    const TwentyFourPuzzle d( TwentyFourPuzzle::parse( line ) );
    bad += d.heuristic( d.initial_state() ) <= 0;
    }
   catch( std::invalid_argument & ) {
    ++bad;
    }
   ++n;
   }
  check( ( n == count ) && ( bad == 0 ) , std::string( "the instances of " )
         + file + " are boards from which the goal is reachable" );
  }

 for( auto metric : { TravelingSalesman::Metric::euclidean ,
                      TravelingSalesman::Metric::manhattan } ) {
  const std::string tag = metric == TravelingSalesman::Metric::euclidean
                          ? "Euclidean" : "Manhattan";
  for( int k : { 6 , 9 , 12 } ) {
   std::ifstream in( data + "/tsp/tsp_40_0" );
   const auto m = TravelingSalesman::read_cities( in , metric , k );
   bool integral = true;
   for( const auto & row : m )
    for( double x : row )
     integral = integral && ( x == std::floor( x ) );
   const TravelingSalesman d( m );
   check_solution( d , astar( d ) , d.held_karp() , "A* on the first " +
                   std::to_string( k ) + " cities of tsp_40_0, " + tag );
   if( metric == TravelingSalesman::Metric::manhattan )
    check( integral , "the Manhattan distances are integers" );
   }
  }
 // the optimal tours of the first 6, 8, 10 and 12 cities of tsp_40_0 with
 // the Manhattan distances, those of the code of HDA* of their authors
 const std::pair< int , double > tours[] = { { 6 , 28975 } , { 8 , 31254 } ,
                                             { 10 , 32734 } ,
                                             { 12 , 33861 } };
 for( const auto & [ k , cost ] : tours ) {
  std::ifstream tf( data + "/tsp/tsp_40_0" );
  const TravelingSalesman d( TravelingSalesman::read_cities(
                       tf , TravelingSalesman::Metric::manhattan , k ) );
  check_solution( d , astar( d ) , cost , "the tour of " + std::to_string( k )
                  + " cities of tsp_40_0 of the reference code" );
  }
 std::ifstream in( data + "/tsp/tsp_40_0" );
 const auto all = TravelingSalesman::read_cities( in );
 std::istringstream few( "3\n0 0\n1 1\n" );
 bool threw = false;
 try {
  ( void ) TravelingSalesman::read_cities( few );
  }
 catch( const std::invalid_argument & ) {
  threw = true;
  }
 check( ( all.size() == 40 ) && threw , "40 cities read, and a missing one "
                                        "refused" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// an admissible but inconsistent heuristic, which forces a reopening: A*
// must still find the optimal cost

void test_reopen( void )
{
 // s = 0, goal = 4; the cheap path 0 - 1 - 2 - 3 - 4 costs 4, while the
 // heuristic of 1 lures A* into 0 - 2 first
 ExplicitGraph< int > g( 5 , 0 , 4 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 2 , 1 , true );
 g.add_edge( 0 , 2 , 3 , true );
 g.add_edge( 2 , 3 , 1 , true );
 g.add_edge( 3 , 4 , 1 , true );
 const int h[ 5 ] = { 0 , 3 , 0 , 1 , 0 };   // admissible: h* = 4 3 2 1 0
 for( std::uint32_t u = 0 ; u < 5 ; ++u )
  g.set_heuristic( u , h[ u ] );
 const auto s = astar( g );
 check_solution( g , s , 4 , "A*, inconsistent heuristic" );
 check( s.stats.reexpanded > 0 , "A* reopens a closed node" );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File files.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
