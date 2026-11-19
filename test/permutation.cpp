/*--------------------------------------------------------------------------*/
/*-------------------------- File permutation.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the permutation
 * puzzles, the Rubik's cube, the Towers of Hanoi and the parallel IDA*.
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
// the distance to a goal of every state reachable from the goals of d, by
// a breadth-first search (the moves of these domains are reversible)

template< class D >
static std::map< state_t< D > , int > distances_from(
 const D & d , const std::vector< state_t< D > > & goals )
{
 std::map< state_t< D > , int > dist;
 std::queue< state_t< D > > q;
 for( const auto & g : goals )
  if( dist.emplace( g , 0 ).second )
   q.push( g );
 while( ! q.empty() ) {
  const auto s = q.front();
  q.pop();
  const int ds = dist[ s ];
  d.for_each_successor( s , cost_t< D >( 0 ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > ,
                               cost_t< D > ) {
   if( dist.emplace( t , ds + 1 ).second )
    q.push( t );
   } );
  }
 return( dist );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if h is admissible and consistent on every edge of the states of
// dist, and 0 at the goals

template< class D , class H >
static bool admissible_consistent(
 const D & d , const std::map< state_t< D > , int > & dist , H && h )
{
 for( const auto & [ s , ds ] : dist ) {
  const int hs = int( h( s ) );
  if( ( hs > ds ) || ( ( ds == 0 ) && ( hs != 0 ) ) )
   return( false );
  bool ok = true;
  d.for_each_successor( s , cost_t< D >( 0 ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > c ,
                               cost_t< D > ) {
   ok = ok && ( hs <= int( c ) + int( h( t ) ) );
   } );
  if( ! ok )
   return( false );
  }
 return( true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the heuristic that for_each_successor() gives each successor is
// h of it, on the states of dist (the incremental heuristics)

template< class D >
static bool successor_h( const D & d ,
                         const std::map< state_t< D > , int > & dist )
{
 bool ok = true;
 for( const auto & [ s , ds ] : dist )
  d.for_each_successor( s , d.heuristic( s ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > ,
                               cost_t< D > ht ) {
   ok = ok && ( ht == d.heuristic( t ) );
   } );
 return( ok );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the rankings: bijections onto [ 0 , N ), lexicographic for LexRanking
// and TreeRanking (the same ranks), and the factorials

template< Ranking R >
static bool ranking_ok( unsigned n , unsigned k , bool lex )
{
 const R r( n , k );
 std::vector< std::uint8_t > a( k ) , b( k ) , prev( k );
 for( std::uint64_t x = 0 ; x < r.size() ; ++x ) {
  r.unrank( x , a.data() );
  std::uint64_t seen = 0;
  for( std::uint8_t v : a ) {
   if( ( v >= n ) || ( ( seen >> v ) & 1u ) )
    return( false );
   seen |= std::uint64_t( 1 ) << v;
   }
  if( r.rank( a.data() ) != x )
   return( false );
  if( lex && ( x > 0 ) && ! std::lexicographical_compare(
                              prev.begin() , prev.end() , a.begin() ,
                              a.end() ) )
   return( false );
  prev = a;
  }
 return( true );
 }

void test_permutation( const std::string & )
{
 check( ( k_permutations( 16 , 6 ) == 5765760 ) &&
        ( k_permutations( 20 , 20 ) == 2432902008176640000ULL ) &&
        ( k_permutations( 21 , 21 ) == 0 ) && ( k_permutations( 3 , 4 ) == 0 ),
        "the numbers of k-permutations" );
 for( unsigned n : { 1u , 4u , 7u } )
  for( unsigned k = 1 ; k <= n ; ++k ) {
   const std::string t = " n " + std::to_string( n ) + " k " +
                         std::to_string( k );
   check( ranking_ok< LexRanking >( n , k , true ) , "lex ranking" + t );
   check( ranking_ok< TreeRanking >( n , k , true ) , "tree ranking" + t );
   check( ranking_ok< MyrvoldRuskeyRanking >( n , k , false ) ,
          "Myrvold-Ruskey ranking" + t );
   }
 // the same lexicographic ranks on 64 items, and the large n of the tree
 const LexRanking l( 64 , 5 );
 const TreeRanking t( 64 , 5 );
 const MyrvoldRuskeyRanking m( 200 , 4 );
 const TreeRanking t2( 200 , 4 );
 std::mt19937_64 gen( 7 );
 bool same = true , back = true;
 for( int i = 0 ; i < 1000 ; ++i ) {
  std::uint8_t a[ 5 ] , b[ 5 ];
  l.unrank( random_below( gen , l.size() ) , a );
  same = same && ( t.rank( a ) == l.rank( a ) );
  const std::uint64_t x = random_below( gen , m.size() );
  m.unrank( x , b );
  back = back && ( m.rank( b ) == x );
  const std::uint64_t y = random_below( gen , t2.size() );
  t2.unrank( y , b );
  back = back && ( t2.rank( b ) == y );
  }
 check( same , "lex and tree ranks agree on 64 items" );
 check( back , "Myrvold-Ruskey and tree on 200 items" );
 bool threw = false;
 try {
  LexRanking bad( 65 , 2 );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "lex ranking refuses 65 items" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the construction of the tables of abstract spaces: unit and integer
// costs against Dijkstra, the unreachable states, the files and the cache

// n states, the goal 0, the edges ( r , r + 1 ) of cost 1 + r % c and
// ( r , 2 r ) of cost c + 1, and the states from n on unreachable
struct ChordSpace {
 std::uint64_t n , extra;
 unsigned c;
 bool chords;
 std::uint64_t size( void ) const { return( n + extra ); }
 template< class F > void for_each_goal( F && f ) const { f( 0 ); }
 template< class F >
 void for_each_predecessor( std::uint64_t r , F && f ) const {
  if( r >= n )
   return;
  if( r + 1 < n )
   f( r + 1 , unsigned( 1 + r % c ) );
  if( r > 0 )
   f( r - 1 , unsigned( 1 + ( r - 1 ) % c ) );
  if( chords && ( 2 * r < n ) && ( r > 0 ) )
   f( 2 * r , c + 1 );
  if( chords && ( r % 2 == 0 ) && ( r > 0 ) )
   f( r / 2 , c + 1 );
  }
 };

void test_abstraction( const std::string & )
{
 for( unsigned c : { 1u , 3u } ) {
  const ChordSpace sp{ 100000 , 7 , c , true };    // 7 blocks
  const auto h = abstract_distances( sp , 3 );
  std::vector< long > ref( sp.n , -1 );        // Dijkstra
  std::priority_queue< std::pair< long , std::uint64_t > ,
                       std::vector< std::pair< long , std::uint64_t > > ,
                       std::greater<> > q;
  ref[ 0 ] = 0;
  q.emplace( 0 , 0 );
  while( ! q.empty() ) {
   const long dv = q.top().first;
   const std::uint64_t v = q.top().second;
   q.pop();
   if( dv != ref[ v ] )
    continue;
   sp.for_each_predecessor( v , [ & ]( std::uint64_t u , unsigned w ) {
    if( ( ref[ u ] < 0 ) || ( dv + long( w ) < ref[ u ] ) ) {
     ref[ u ] = dv + long( w );
     q.emplace( ref[ u ] , u );
     }
    } );
   }
  bool ok = true;
  for( std::uint64_t r = 0 ; r < sp.n ; ++r )
   ok = ok && ( long( h[ r ] ) == ref[ r ] );
  for( std::uint64_t r = sp.n ; r < sp.size() ; ++r )
   ok = ok && ( h[ r ] == 255 );
  check( ok , "abstract distances against Dijkstra, costs 1 to " +
         std::to_string( c + 1 ) );
  }
 bool threw = false;
 try {
  static_cast< void >( abstract_distances( ChordSpace{ 300 , 0 , 1 , false } ,
                                           2 ) );
  }
 catch( std::logic_error & ) {
  threw = true;
  }
 check( threw , "a distance beyond 254 is refused" );

 const ChordSpace small{ 200 , 5 , 1 , false };
 const auto h = abstract_distances< ADS::ThreadPool >( small , 2 );
 const AbstractionTable< 8 > t8( h );
 const AbstractionTable< 4 > t4( h );
 bool ok = true;
 for( std::uint64_t r = 0 ; r < small.size() ; ++r )
  ok = ok && ( t8[ r ] == h[ r ] ) &&
       ( t4[ r ] == std::min< unsigned >( h[ r ] , 15 ) );
 check( ok && ( 2 * t4.bytes() <= t8.bytes() + 8 ) ,
        "the tables of 8 and 4 bits" );
 const auto dir = std::filesystem::temp_directory_path() /
                  ( "pbefs-abs-" + std::to_string( std::random_device()() ) );
 std::filesystem::create_directories( dir );
 const std::string f = ( dir / "t.pdb" ).string();
 t4.save( f , 42 );
 const auto back = AbstractionTable< 4 >::load( f , 42 , small.size() );
 ok = true;
 for( std::uint64_t r = 0 ; r < small.size() ; ++r )
  ok = ok && ( back[ r ] == t4[ r ] );
 check( ok , "a table saved and loaded" );
 threw = false;
 try {
  static_cast< void >( AbstractionTable< 4 >::load( f , 43 , small.size() ) );
  }
 catch( std::runtime_error & ) {
  threw = true;
  }
 check( threw , "a table is not loaded with another key" );
 int built = 0;
 const auto make = [ & ] {
  ++built;
  return( abstract_distances( small , 1 ) );
  };
 static_cast< void >( AbstractionTable< 8 >::cached( dir.string() , "c.pdb" ,
                                                     7 , small.size() ,
                                                     make ) );
 const auto c2 = AbstractionTable< 8 >::cached( dir.string() , "c.pdb" , 7 ,
                                                small.size() , make );
 check( ( built == 1 ) && ( c2[ 150 ] == 150 ) ,
        "a table built once, then read from the cache" );
 std::filesystem::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the pancakes: the gap heuristics against the exact distances of every
// stack of 7, the searches against them, and the pattern databases

void test_pancake( const std::string & )
{
 using P7 = PancakePuzzle< 7 >;
 const P7 any( P7::random( 1 ) );
 P7::state_type id{};
 for( int i = 0 ; i < 7 ; ++i )
  id[ i ] = std::uint8_t( i );
 const auto dist = distances_from( any , { id } );
 check( dist.size() == 5040 , "the 7! stacks of 7 pancakes" );
 check( admissible_consistent( any , dist , [ & ]( const P7::state_type & s ) {
         return( any.heuristic( s ) ); } ) && successor_h( any , dist ) ,
        "the gap heuristic is consistent, and incremental" );
 const P7 gap2( P7::random( 1 ) , 2 );
 check( admissible_consistent( gap2 , dist , [ & ]( const P7::state_type & s )
         { return( gap2.heuristic( s ) ); } ) && successor_h( gap2 , dist ) ,
        "GAP-2 is consistent, and incremental" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 97 ) != 0 )
   continue;
  const P7 d( s );
  wrong += ( astar( d ).cost != ds ) || ( idastar( d ).cost != ds );
  }
 check( wrong == 0 , "A* and IDA* on 52 stacks of 7" );

 using B5 = BurnedPancakePuzzle< 5 >;
 const B5 bany( B5::random( 1 ) );
 B5::state_type bid{};
 for( int i = 0 ; i < 5 ; ++i )
  bid[ i ] = std::uint8_t( i );
 const auto bdist = distances_from( bany , { bid } );
 check( bdist.size() == 3840 , "the 5! 2^5 stacks of 5 burned pancakes" );
 check( admissible_consistent( bany , bdist , []( const B5::state_type & s ) {
         return( B5::heuristic( s ) ); } ) && successor_h( bany , bdist ) ,
        "the burned gap heuristic is consistent, and incremental" );
 check( B5::parse( "-3 1 -2 5 4" ) == B5::state_type{ 0x82 , 0 , 0x81 , 4 , 3 }
        , "a burned stack in text" );

 // the tables of 4 of 7 pancakes and of 3 of 5 burned ones, and their
 // maximum, within the exact distances; the files of the cache
 PermutationPDB< LexRanking , 8 > p( any , 0 , { 0 , 1 , 2 , 3 } , 2 );
 PermutationPDB< MyrvoldRuskeyRanking , 8 > q( any , 0 , { 3 , 4 , 5 , 6 } ,
                                              2 );
 check( p.table().size() == 840 , "the table of 4 of 7 pancakes" );
 bool ok = true;
 for( const auto & [ s , ds ] : dist )
  ok = ok && ( int( p.value( any , s ) ) <= ds ) &&
       ( int( q.value( any , s ) ) <= ds );
 check( ok , "the tables of 4 pancakes are admissible" );
 const PermutationHeuristic< P7 > ph( any , { p } );
 const PDBDomain< P7 , PermutationHeuristic< P7 > > pd( any , ph );
 check( admissible_consistent( any , dist , [ & ]( const P7::state_type & s ) {
         return( pd.heuristic( s ) ); } ) ,
        "the maximum of the gap and of a table is consistent" );
 PermutationPDB< TreeRanking , 4 > b( bany , 0 , { 0 , 1 , 2 } , 1 );
 ok = b.table().size() == 60 * 8;
 for( const auto & [ s , ds ] : bdist )
  ok = ok && ( int( b.value( bany , s ) ) <= ds );
 check( ok , "the table of 3 burned pancakes is admissible" );

 const auto dir = std::filesystem::temp_directory_path() /
                  ( "pbefs-perm-" + std::to_string( std::random_device()() ) );
 PermutationPDB< LexRanking , 8 > c1( any , 0 , { 0 , 1 , 2 } , 1 ,
                                      dir.string() );
 PermutationPDB< LexRanking , 8 > c2( any , 0 , { 0 , 1 , 2 } , 1 ,
                                      dir.string() );
 ok = std::distance( std::filesystem::directory_iterator( dir ) ,
                     std::filesystem::directory_iterator() ) == 1;
 for( std::uint64_t r = 0 ; r < c1.table().size() ; ++r )
  ok = ok && ( c1.table()[ r ] == c2.table()[ r ] );
 check( ok , "a table of pancakes in the cache" );
 std::filesystem::remove_all( dir );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// TopSpin: the breakpoints and a table against the exact distances of
// TopSpin( 8 , 4 ), and the searches

void test_topspin( const std::string & )
{
 using T = TopSpin< 8 , 4 >;
 const T any( T::random( 1 , 10 ) );
 std::vector< T::state_type > goals;
 T::for_each_goal_permutation( 0 , [ & ]( const std::uint8_t * o ,
                                          const std::uint8_t * ) {
  T::state_type g;
  std::copy( o , o + 8 , g.begin() );
  goals.push_back( g );
  } );
 const auto dist = distances_from( any , goals );
 check( dist.size() == 40320 , "the 8! rings of TopSpin( 8 , 4 )" );
 check( admissible_consistent( any , dist , []( const T::state_type & s ) {
         return( T::heuristic( s ) ); } ) ,
        "the breakpoints of TopSpin are consistent" );
 PermutationPDB< LexRanking , 8 > p( any , 0 , { 0 , 1 , 2 , 3 , 4 } , 2 );
 const PermutationHeuristic< T > ph( any , { p } );
 const PDBDomain< T , PermutationHeuristic< T > > pd( any , ph );
 check( admissible_consistent( any , dist , [ & ]( const T::state_type & s ) {
         return( pd.heuristic( s ) ); } ) ,
        "a table of 5 tokens of TopSpin is consistent" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 2017 ) != 0 )
   continue;
  const T d( s );
  const PDBDomain< T , PermutationHeuristic< T > > e( d , ph );
  wrong += ( idastar( d ).cost != ds ) || ( idastar( e ).cost != ds ) ||
           ( astar( e ).cost != ds );
  }
 check( wrong == 0 , "IDA* and A* on 20 rings of TopSpin( 8 , 4 )" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// Rubik's cube: the size of the tree of table 1 of Korf (1997), the states
// at each depth, the instances of Korf, and small tables

void test_rubik( const std::string & data )
{
 const RubiksCube c;
 std::vector< std::uint64_t > tree( 6 , 0 );
 std::function< void( const RubiksCube::state_type & , std::uint8_t , int )
                > down = [ & ]( const RubiksCube::state_type & s ,
                                std::uint8_t op , int depth ) {
  if( depth == 5 )
   return;
  c.for_each_tree_successor( s , 0 , op , [ & ](
   const RubiksCube::state_type & t , int , int , std::uint8_t o ) {
   ++tree[ std::size_t( depth + 1 ) ];
   down( t , o , depth + 1 );
   } );
  };
 down( c.initial_state() , no_op , 0 );
 check( tree == std::vector< std::uint64_t >{ 0 , 18 , 243 , 3240 , 43254 ,
                                              577368 } ,
        "the tree of the cube, table 1 of Korf (1997), to depth 5" );
 std::set< RubiksCube::state_type > seen{ c.initial_state() } ,
                                    front{ c.initial_state() };
 std::vector< std::size_t > layer;
 for( int depth = 1 ; depth <= 4 ; ++depth ) {
  std::set< RubiksCube::state_type > next;
  for( const auto & s : front )
   c.for_each_successor( s , 0 , [ & ]( const RubiksCube::state_type & t ,
                                        int , int ) {
    if( seen.insert( t ).second )
     next.insert( t );
    } );
  layer.push_back( next.size() );
  front.swap( next );
  }
 check( layer == std::vector< std::size_t >{ 18 , 243 , 3240 , 43239 } ,
        "the states of the cube at the depths 1 to 4" );

 // every move has order 4 (a quarter turn) or 2, and its inverse
 bool ok = true;
 const auto s0 = RubiksCube::random_walk( 30 , 3 );
 for( int m = 0 ; m < 18 ; ++m ) {
  auto s = s0;
  for( int i = 0 ; i < 4 ; ++i )
   s = RubiksCube::apply( s , m );
  ok = ok && ( s == s0 ) &&
       ( RubiksCube::apply( RubiksCube::apply( s0 , m ) ,
                            3 * ( m / 3 ) + 2 - m % 3 ) == s0 );
  }
 check( ok , "the moves of the cube and their inverses" );

 std::ifstream in( data + "/rubik/korf10.txt" );
 int n = 0;
 for( std::string line ; std::getline( in , line ) ; )
  if( ! line.empty() ) {
   ++n;
   ok = ok && ( RubiksCube::parse_moves( line ).size() == 100 );
   }
 check( ok && ( n == 10 ) , "the 10 instances of Korf, of 100 moves" );

 // small tables: 4 corners and 4 + 4 edges, against IDA* without them
 std::vector< PermutationPDB< LexRanking , 4 > > v;
 v.emplace_back( c , 0 , std::vector< std::uint8_t >{ 0 , 1 , 2 , 3 } , 0 );
 v.emplace_back( c , 1 , std::vector< std::uint8_t >{ 0 , 1 , 2 , 3 } , 0 );
 v.emplace_back( c , 1 , std::vector< std::uint8_t >{ 4 , 5 , 6 , 7 } , 0 );
 check( ( v[ 0 ].table().size() == 1680 * 81 ) &&
        ( v[ 1 ].table().size() == 11880 * 16 ) ,
        "the sizes of the tables of 4 corners and of 4 edges" );
 const PermutationHeuristic< RubiksCube , LexRanking , 4 > h( c , v );
 for( std::uint64_t seed = 1 ; seed <= 4 ; ++seed ) {
  const RubiksCube d( RubiksCube::random_walk( 5 , seed ) );
  const PDBDomain< RubiksCube , decltype( h ) > e( d , h );
  const auto ref = idastar( d );
  check_solution( e , idastar( e ) , ref.cost , "IDA* with tables, walk " +
                  std::to_string( seed ) );
  PIDAOptions o;
  o.threads = 3;
  check_solution( e , ParallelIDAStar< decltype( e ) >( e , o ).solve() ,
                  ref.cost , "PIDA* with tables, walk " +
                  std::to_string( seed ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the tables of Korf (1997) for Rubik's cube: their sizes, largest values
// and means, against those of his paper (88179840 entries of mean 8.764
// and largest value 11; 42577920 of mean 7.668 and largest 10)

void test_rubik_pdb( const std::string & )
{
 #if ! defined( NDEBUG ) || ( defined( ADS_CONTRACTS ) && ADS_CONTRACTS )
  check( true , "the tables of Korf only in an optimized build without the "
         "checks of ADS" );
  return;
 #endif
 const RubiksCube c;
 const RubikPDBHeuristic h( c );
 const std::uint64_t size[] = { 88179840 , 42577920 , 42577920 };
 const unsigned top[] = { 11 , 10 , 10 };
 const double mean[] = { 8.764 , 7.668 , 7.668 };
 for( std::size_t i = 0 ; i < 3 ; ++i ) {
  const auto hist = h.pdbs()[ i ].table().histogram();
  double sum = 0;
  unsigned largest = 0;
  for( std::size_t x = 0 ; x < hist.size() ; ++x )
   if( hist[ x ] ) {
    sum += double( x ) * double( hist[ x ] );
    largest = unsigned( x );
    }
  const double m = sum / double( size[ i ] );
  std::ostringstream os;
  os << "the table " << i << " of Korf: " << h.pdbs()[ i ].table().size()
     << " entries, largest " << largest << ", mean " << m;
  check( ( h.pdbs()[ i ].table().size() == size[ i ] ) &&
         ( largest == top[ i ] ) && ( std::abs( m - mean[ i ] ) < 5e-4 ) ,
         os.str() );
  }
 for( std::uint64_t seed = 1 ; seed <= 2 ; ++seed ) {
  const RubiksCube d( RubiksCube::random_walk( 10 , seed ) );
  const PDBDomain< RubiksCube , RubikPDBHeuristic > e( d , h );
  const auto r = idastar( e );
  PIDAOptions o;
  o.threads = 4;
  check_solution( e , ParallelIDAStar< decltype( e ) >( e , o ).solve() ,
                  r.cost , "PIDA* with the tables of Korf, walk of 10, seed "
                  + std::to_string( seed ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the Towers of Hanoi with four pegs: the instances of Siag et al., the
// tables against the exact distances of 8 disks, and the searches

void test_hanoi( const std::string & data )
{
 std::ifstream in( data + "/hanoi/toh20.txt" );
 int n = 0;
 bool ok = true;
 for( std::string line ; std::getline( in , line ) ; ) {
  if( line.empty() || ( line[ 0 ] != '(' ) )
   continue;
  std::istringstream is( line );
  const auto s = TowersOfHanoi< 20 >::read( is );
  const auto g = TowersOfHanoi< 20 >::read( is );
  while( ! line.empty() &&
         std::isspace( static_cast< unsigned char >( line.back() ) ) )
   line.pop_back();
  ok = ok && ( TowersOfHanoi< 20 >::write( s ) + " " +
               TowersOfHanoi< 20 >::write( g ) == line );
  ++n;
  }
 check( ok && ( n == 20 ) , "the 20 instances of Siag et al., read and "
        "written back" );

 using H8 = TowersOfHanoi< 8 >;
 std::mt19937_64 gen( 5 );
 std::uint64_t goal = 0;
 for( int k = 0 ; k < 8 ; ++k )
  goal |= random_below( gen , 4 ) << ( 2 * k );
 const H8 any( 0 , goal );
 const auto dist = distances_from( any , { goal } );
 check( dist.size() == 65536 , "the 4^8 states of 8 disks" );
 const HanoiPDB p( { 3 , 5 } , goal , 2 );
 const HanoiPDB q( { 8 } , goal , 2 );
 ok = true;
 for( const auto & [ s , ds ] : dist )
  ok = ok && ( int( q( s ) ) == ds );
 check( ok , "the table of all the disks is the exact distance" );
 check( admissible_consistent( any , dist , [ & ]( std::uint64_t s ) {
         return( p( s ) ); } ) , "the tables of 3 + 5 disks are additive "
        "and consistent" );
 int k = 0 , wrong = 0;
 for( const auto & [ s , ds ] : dist ) {
  if( ( k++ % 4099 ) != 0 )
   continue;
  const H8 d( s , goal );
  const PDBDomain< H8 , HanoiPDB > e( d , p );
  PIDAOptions o;
  o.threads = 2;
  o.frontier = 256;
  o.table_depth = 8;
  wrong += ( astar( e ).cost != ds ) ||
           ( ParallelIDAStar< decltype( e ) >( e , o ).solve().cost != ds );
  }
 check( wrong == 0 , "A* and PIDA* on 16 instances of 8 disks" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the road graphs of DIMACS: a graph and its coordinates in the format of
// the challenge, the consistency of the heuristic, A* against Dijkstra

void test_dimacs( const std::string & data )
{
 std::ifstream gr( data + "/dimacs/sample.gr" ) ,
               co( data + "/dimacs/sample.co" );
 check( gr && co , "the sample graph is read" );
 auto g = read_dimacs_graph< std::int64_t >( gr , 1 , 2 );
 const auto c = read_dimacs_coordinates( co , g.size() );
 const double kappa = dimacs_scale( g , c );
 check( kappa > 0 , "the factor of the heuristic is positive" );
 std::mt19937_64 gen( 3 );
 for( int i = 0 ; i < 5 ; ++i ) {
  g.set_source( std::uint32_t( random_below( gen , g.size() ) ) );
  g.set_target( std::uint32_t( random_below( gen , g.size() ) ) );
  set_dimacs_heuristic( g , c , kappa );
  bool ok = g.heuristic( g.target() ) == 0;
  for( std::uint32_t u = 0 ; u < g.size() ; ++u )
   g.for_each_successor( u , 0 , [ & ]( std::uint32_t , std::int64_t w ,
                                        std::int64_t hv ) {
    ok = ok && ( g.heuristic( u ) <= w + hv );
    } );
  check( ok , "the heuristic of the coordinates is consistent, pair " +
         std::to_string( i ) );
  const auto ref = astar( Blind< ExplicitGraph< std::int64_t > >{ g } ,
                          OpenKind::heap );
  const auto a = astar( g , OpenKind::heap );
  if( ref.solved() )
   check_solution( g , a , ref.cost , "A* on the sample, pair " +
                   std::to_string( i ) );
  else
   check( ! a.solved() , "A* on the sample, pair " + std::to_string( i ) );
  }
 std::istringstream bad( "p sp 2 1\na 1 3 5\n" );
 bool threw = false;
 try {
  static_cast< void >( read_dimacs_graph( bad ) );
  }
 catch( std::invalid_argument & ) {
  threw = true;
  }
 check( threw , "an arc out of range is refused" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// ParallelIDAStar with every frontier rule, schedule and scheme, on the
// given runtime: the 8-puzzle against its exact distances, the 15-puzzle,
// the pancakes, the map of Romania, the node limit, the unsolvable graphs,
// IDA* reproduced with one unit and one worker, the deterministic path

template< class Executor >
void test_pida( const std::string & data )
{
 const auto inst = korf( data );
 // the 15-puzzle of the cases: the instance 55 of Korf, or in a light run a
 // walk of 50 moves
 struct Fifteen {
  FifteenPuzzle::board_type board;
  int cost;
  std::string name;
  };
 Fifteen fifteen{ inst[ easy_korf[ 0 ] - 1 ].board ,
                  inst[ easy_korf[ 0 ] - 1 ].cost ,
                  "Korf #" + std::to_string( easy_korf[ 0 ] ) };
 if( light_run() ) {
  fifteen.board = random_walk< FifteenPuzzle >( 50 , 3 );
  fifteen.cost = astar( FifteenPuzzle( fifteen.board ) ).cost;
  fifteen.name = "walk";
  }
 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 // the variants: rule, schedule, unit bounds, scheme, split, table depth
 struct Variant {
  FrontierRule rule;
  UnitSchedule schedule;
  bool unit_bounds , asynchronous;
  double split;
  unsigned table;
  const char * name;
  };
 const Variant variants[] = {
  { FrontierRule::astar , UnitSchedule::steal , true , false , 0 , 0 ,
    "astar steal" } ,
  { FrontierRule::breadth , UnitSchedule::shared , true , false , 0 , 0 ,
    "breadth shared" } ,
  { FrontierRule::depth , UnitSchedule::shared , false , false , 0 , 0 ,
    "depth classical" } ,
  { FrontierRule::astar , UnitSchedule::steal , true , true , 0 , 0 ,
    "astar asynchronous" } ,
  { FrontierRule::depth , UnitSchedule::steal , true , false , 2 , 0 ,
    "depth split" } ,
  { FrontierRule::astar , UnitSchedule::shared , true , true , 2 , 3 ,
    "astar asynchronous split table" } ,
  { FrontierRule::breadth , UnitSchedule::steal , true , false , 0 , 4 ,
    "breadth table" } };
 const auto options = []( const Variant & v , unsigned p ) {
  PIDAOptions o;
  o.threads = p;
  o.rule = v.rule;
  o.schedule = v.schedule;
  o.unit_bounds = v.unit_bounds;
  o.asynchronous = v.asynchronous;
  o.split = v.split;
  o.table_depth = v.table;
  return( o );
  };
 for( const Variant & v : variants )
  for( unsigned p : light_run() ? std::vector< unsigned >{ 1 , 4 }
                                 : std::vector< unsigned >{ 1 , 2 , 4 } ) {
   const PIDAOptions o = options( v , p );
   const std::string tag = std::string( " " ) + v.name + ", p " +
                           std::to_string( p );
   std::size_t k = 0;
   int wrong = 0;
   for( const auto & [ s , dst ] : dist ) {
    if( ( k++ % 18143 ) != 0 )                  // 10 states
     continue;
    const EightPuzzle e( EightPuzzle::unpack( s ) );
    const auto r = ParallelIDAStar< EightPuzzle , Executor >( e , o )
                   .solve();
    if( ! r.solved() || ( r.cost != dst ) ) {
     ++wrong;
     check_solution( e , r , dst , "PIDA* 8-puzzle" + tag );
     }
    }
   check( wrong == 0 , "PIDA* on 10 states of the 8-puzzle" + tag );
   const FifteenPuzzle d( fifteen.board );
   check_solution( d , ParallelIDAStar< FifteenPuzzle , Executor >( d , o )
                        .solve() , fifteen.cost , "PIDA* " + fifteen.name +
                   tag );
   const PancakePuzzle< 12 > pc( PancakePuzzle< 12 >::random( p ) );
   check_solution( pc , ParallelIDAStar< PancakePuzzle< 12 > , Executor >(
                         pc , o ).solve() , astar( pc ).cost ,
                   "PIDA* pancakes of 12" + tag );
   if( v.table == 0 ) {
    auto g = romania( Arad );
    check_solution( g , ParallelIDAStar< ExplicitGraph< int > , Executor >(
                         g , o ).solve() , 418 , "PIDA* Romania" + tag );
    }
   }

 // one unit and one worker: the expansions of IDA*
 for( int id : light_run() ? std::vector< int >{ 0 }
                           : std::vector< int >{ easy_korf[ 0 ] ,
                                                 easy_korf[ 1 ] } ) {
  const FifteenPuzzle d( id ? inst[ id - 1 ].board : fifteen.board );
  PIDAOptions o;
  o.frontier = 1;
  const auto r = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  const auto i = idastar( d );
  check( ( r.cost == i.cost ) && ( r.stats.expanded == i.stats.expanded ) &&
         ( r.stats.generated == i.stats.generated ) && ( r.path == i.path ) ,
         "PIDA* with one unit is IDA*, " +
         ( id ? "Korf #" + std::to_string( id ) : fifteen.name ) );
  }

 // the deterministic path, the same with 1 and 3 workers
 {
  const FifteenPuzzle d( inst[ easy_korf[ 2 ] - 1 ].board );
  PIDAOptions o;
  o.deterministic = true;
  o.frontier = 64;
  const auto a = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  o.threads = 3;
  const auto b = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  check( a.solved() && ( a.path == b.path ) && ( a.cost == b.cost ) ,
         "PIDA* deterministic: the same path with 1 and 3 workers" );
  bool threw = false;
  o.asynchronous = true;
  try {
   ParallelIDAStar< FifteenPuzzle , Executor > x( d , o );
   }
  catch( std::invalid_argument & ) {
   threw = true;
   }
  check( threw , "PIDA* refuses deterministic and asynchronous" );
  }

 // the node limit, a goal in the frontier, the unsolvable graphs
 {
  const FifteenPuzzle d( inst[ 0 ].board );
  PIDAOptions o;
  o.threads = 3;
  o.node_limit = 1000;
  const auto r = ParallelIDAStar< FifteenPuzzle , Executor >( d , o ).solve();
  check( ( r.status == Status::node_limit ) && r.path.empty() ,
         "PIDA* stops at the expansion limit" );
  const FifteenPuzzle e( random_walk< FifteenPuzzle >( 6 , 2 ) );
  PIDAOptions f;
  f.frontier = 100000;
  check_solution( e , ParallelIDAStar< FifteenPuzzle , Executor >( e , f )
                       .solve() , idastar( e ).cost ,
                  "PIDA* whose frontier reaches the goal" );
  ExplicitGraph< int > dag( 5 , 0 , 4 );
  dag.add_edge( 0 , 1 , 1 , true );
  dag.add_edge( 1 , 2 , 1 , true );
  dag.add_edge( 0 , 3 , 2 , true );
  for( FrontierRule rule : { FrontierRule::astar , FrontierRule::depth } )
   for( bool async : { false , true } ) {
    PIDAOptions u;
    u.threads = 2;
    u.frontier = 2;
    u.rule = rule;
    u.asynchronous = async;
    check( ParallelIDAStar< ExplicitGraph< int > , Executor >( dag , u )
           .solve().status == Status::unsolvable ,
           "PIDA* proves an acyclic graph unsolvable" );
    }
  }
 }

#ifdef PBEFS_HAS_MPI

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// ParallelIDAStarMPI on the 15-puzzle, the pancakes, the Towers of Hanoi
// and Romania, with every frontier rule

void test_pida_mpi( const std::string & data )
{
 const auto inst = korf( data );
 for( FrontierRule rule : { FrontierRule::astar , FrontierRule::breadth ,
                            FrontierRule::depth } ) {
  PIDAOptions o;
  o.rule = rule;
  const std::string tag = rule == FrontierRule::astar ? " astar"
                          : rule == FrontierRule::breadth ? " breadth"
                                                          : " depth";
  for( int id : { easy_korf[ 0 ] , easy_korf[ 2 ] } ) {
   const FifteenPuzzle d( inst[ id - 1 ].board );
   check_solution( d , ParallelIDAStarMPI< FifteenPuzzle >( d , o ).solve() ,
                   inst[ id - 1 ].cost , "MPI PIDA* Korf #" +
                   std::to_string( id ) + tag );
   }
  const PancakePuzzle< 14 > pc( PancakePuzzle< 14 >::random( 3 ) );
  check_solution( pc , ParallelIDAStarMPI< PancakePuzzle< 14 > >( pc , o )
                        .solve() , astar( pc ).cost , "MPI PIDA* pancakes"
                  + tag );
  auto g = romania( Arad );
  check_solution( g , ParallelIDAStarMPI< ExplicitGraph< int > >( g , o )
                       .solve() , 418 , "MPI PIDA* Romania" + tag );
  }
 PIDAOptions l;
 l.node_limit = 500;
 const FifteenPuzzle d( inst[ 0 ].board );
 check( ParallelIDAStarMPI< FifteenPuzzle >( d , l ).solve().status ==
        Status::node_limit , "MPI PIDA* stops at the expansion limit" );
 ExplicitGraph< int > dag( 4 , 0 , 3 );
 dag.add_edge( 0 , 1 , 1 , true );
 dag.add_edge( 1 , 2 , 1 , true );
 check( ParallelIDAStarMPI< ExplicitGraph< int > >( dag ).solve().status ==
        Status::unsolvable , "MPI PIDA* proves an acyclic graph unsolvable" );
 }

#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_pida< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_pida< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_pida< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------ End File permutation.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
