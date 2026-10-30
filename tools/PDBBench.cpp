/*--------------------------------------------------------------------------*/
/*--------------------------- File PDBBench.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the pattern databases of the sliding-tile puzzles [see
 * PatternDatabase.h]. In the mode table it builds (or loads from the cache
 * directory) the tables of a partition with an encoding and prints their
 * size, the time to have them, the mean heuristic and the mean Manhattan
 * distance on random states, and the time of a lookup; in the mode search
 * it solves instances of the 15-puzzle or of the 24-puzzle by A*, IDA* or
 * HDA* with the Manhattan distance or with the tables, one line of
 * comma-separated values per run.
 *
 *     PDBBench [options]
 *
 *     --mode     search | table                            (search)
 *     --puzzle   15 | 24                                   (15)
 *     --file     <instances>              korf100.txt, instances50.txt
 *     --ids      <a-b,c,...>              instances of the file (all)
 *     --heur     manhattan | 555 | 663 | 78 | 6666         (78)
 *     --enc      byte | nibble | mod3 | bit | min | sparse (nibble)
 *     --compress <ranks per entry of min>                  (1)
 *     --blank    ignored | minimized                       (ignored)
 *     --reflect  0 | 1                                     (1)
 *     --cache    <directory of the tables>                 (none)
 *     --build    <workers of the construction, 0 all>     (0)
 *     --algo     astar | ida | hda                         (astar)
 *     --threads  <workers of HDA*>                         (1)
 *     --limit    <nodes, expansions for IDA*, of a run>    (none)
 *     --samples  <random states of the mode table>        (1000000)
 *     --repeat   <runs of each instance>                   (1)
 *     --header   print the header line
 *
 * The partition 6666 is that of the 24-puzzle, the others those of the
 * 15-puzzle; on the 24-puzzle the search of Korf and Felner is IDA*, since
 * A* and HDA* store every state they generate; mod3 and bit need
 * --reflect 0, sparse a build with the structures of ADS. In the mode
 * table the lookup of mod3 and bit is the change of the heuristic along a
 * move (AdditivePDB::delta()), that of the others the whole heuristic.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>

#include <cstdint>

#include <fstream>

#include <iostream>

#include <map>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/SlidingTile.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------- TYPES AND GLOBAL VARIABLES -----------------------*/
/*--------------------------------------------------------------------------*/

// the options of the command line, by name
static std::map< std::string , std::string > opt = {
 { "mode" , "search" } , { "puzzle" , "15" } , { "file" , "" } ,
 { "ids" , "" } ,
 { "heur" , "78" } , { "enc" , "nibble" } , { "compress" , "1" } ,
 { "blank" , "ignored" } , { "reflect" , "1" } , { "cache" , "" } ,
 { "build" , "0" } , { "algo" , "astar" } , { "threads" , "1" } ,
 { "samples" , "1000000" } , { "repeat" , "1" } , { "limit" , "0" } };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the indices in a list "a-b,c,..." (1-based), or all of [ 1 , n ]

static std::vector< int > parse_ids( const std::string & s , int n )
{
 std::vector< int > v;
 if( s.empty() ) {
  for( int i = 1 ; i <= n ; ++i )
   v.push_back( i );
  return( v );
  }
 std::istringstream in( s );
 std::string part;
 while( std::getline( in , part , ',' ) ) {
  const auto dash = part.find( '-' );
  const int a = std::stoi( part.substr( 0 , dash ) );
  const int b = dash == std::string::npos
                ? a : std::stoi( part.substr( dash + 1 ) );
  for( int i = a ; i <= b ; ++i )
   v.push_back( i );
  }
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the seconds elapsed since t

static double since( std::chrono::steady_clock::time_point t )
{
 return( std::chrono::duration< double >( std::chrono::steady_clock::now() -
                                          t ).count() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the partition of --heur

static Partition partition( void )
{
 const std::string & h = opt[ "heur" ];
 if( h == "555" )
  return( partition_15_5_5_5() );
 if( h == "663" )
  return( partition_15_6_6_3() );
 if( h == "78" )
  return( partition_15_7_8() );
 if( h == "6666" )
  return( partition_24_6_6_6_6() );
 throw( std::invalid_argument( "unknown partition " + h ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of the construction from the command line

static PDBOptions pdb_options( void )
{
 PDBOptions o;
 if( opt[ "blank" ] == "minimized" )
  o.blank = BlankRule::minimized;
 else
  if( opt[ "blank" ] != "ignored" )
   throw( std::invalid_argument( "unknown rule " + opt[ "blank" ] ) );
 o.reflect = opt[ "reflect" ] != "0";
 o.compress = unsigned( std::stoul( opt[ "compress" ] ) );
 o.threads = unsigned( std::stoul( opt[ "build" ] ) );
 o.cache = opt[ "cache" ];
 return( o );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// prints the line of a run

template< class S >
static void print( int id , const S & s )
{
 const auto & st = s.stats;
 const char * status[] = { "solved" , "unsolvable" , "node_limit" };
 std::cout.precision( 12 );
 std::cout << id << ',' << opt[ "heur" ] << ',' << opt[ "enc" ] << ','
           << opt[ "compress" ] << ',' << opt[ "blank" ] << ','
           << opt[ "reflect" ] << ',' << opt[ "algo" ] << ','
           << opt[ "threads" ] << ',' << status[ int( s.status ) ] << ','
           << s.cost << ',' << st.expanded << ',' << st.reexpanded << ','
           << st.generated << ',' << st.seconds << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d by --algo, --repeat times

template< class D >
static void solve( int id , const D & d )
{
 SearchOptions so;
 so.node_limit = std::stoul( opt[ "limit" ] );
 for( int r = std::stoi( opt[ "repeat" ] ) ; r > 0 ; --r )
  if( opt[ "algo" ] == "astar" )
   print( id , AStar< D >( d , so ).solve() );
  else
   if( opt[ "algo" ] == "ida" )
    print( id , IDAStar< D >( d , so ).solve() );
   else
    if( opt[ "algo" ] == "hda" ) {
     ParallelOptions po;
     static_cast< SearchOptions & >( po ) = so;
     po.threads = unsigned( std::stoul( opt[ "threads" ] ) );
     const ZobristDistribution< D > z( d );
     print( id , HDAStar< D , ZobristDistribution< D > , ThreadRuntime >(
                  d , z , po ).solve() );
     }
    else
     throw( std::invalid_argument( "unknown search " + opt[ "algo" ] ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the instances of the puzzle D of --file

template< class D >
static std::vector< typename D::board_type > instances( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
 std::vector< typename D::board_type > boards;
 std::string line;
 while( std::getline( in , line ) )
  if( ! line.empty() && ( line[ 0 ] != '#' ) )
   boards.push_back( D::parse( line ) );
 return( boards );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the mode search on the board W x H with the tables encoded by Enc

template< int W , int H , class Enc >
static void search( void )
{
 using D = SlidingTile< W , H >;
 const auto boards = instances< D >();
 if( opt[ "heur" ] == "manhattan" ) {
  for( int id : parse_ids( opt[ "ids" ] , int( boards.size() ) ) )
   solve( id , D( boards.at( std::size_t( id - 1 ) ) ) );
  return;
  }
 const auto pdb = AdditivePDB< W , H , Enc >::template make< ThreadRuntime >(
                   partition() , D::korf_goal() , pdb_options() );
 for( int id : parse_ids( opt[ "ids" ] , int( boards.size() ) ) )
  solve( id , SlidingTilePDB< W , H , Enc >( boards.at( std::size_t( id - 1 )
                                                       ) , pdb ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the mode table on the board W x H with the tables encoded by Enc

template< int W , int H , class Enc >
static void table( void )
{
 constexpr int n = W * H;
 using PDB = AdditivePDB< W , H , Enc >;
 const auto t0 = std::chrono::steady_clock::now();
 const PDB pdb = PDB::template make< ThreadRuntime >(
                  partition() , TilePattern< W , H >::korf_goal() ,
                  pdb_options() );
 const double build = since( t0 );
 std::uint64_t entries = 0;
 bool loaded = true;
 for( const auto & t : pdb.tables() ) {
  entries += t.pattern().entries();
  loaded = loaded && t.from_file();
  }

 // random states, the cells of their tiles
 const std::size_t m = std::stoul( opt[ "samples" ] );
 std::vector< std::array< std::uint8_t , n > > q( m );
 std::mt19937_64 rng( 1 );
 std::array< std::uint8_t , n > perm{};
 for( int i = 0 ; i < n ; ++i )
  perm[ i ] = std::uint8_t( i );
 for( auto & x : q ) {
  std::shuffle( perm.begin() , perm.end() , rng );
  x = perm;
  }
 const auto manhattan = [ & ]( const std::array< std::uint8_t , n > & x ) {
  int s = 0;
  for( int t = 1 ; t < n ; ++t )
   s += TilePattern< W , H >::distance( x[ t ] , t );
  return( s );
  };

 double mh = 0 , md = 0;
 long long sink = 0;
 const auto t1 = std::chrono::steady_clock::now();
 if constexpr( Enc::relative ) {
  // the change along the move of a tile next to the blank into it
  for( const auto & x : q ) {
   const int b = x[ 0 ];
   const int c = TilePattern< W , H >::moves( b ).cell[ 0 ];
   int tile = 0;
   for( int t = 1 ; t < n ; ++t )
    if( x[ t ] == c )
     tile = t;
   sink += pdb.delta( x.data() , tile , b );
   }
  }
 else
  for( const auto & x : q )
   sink += pdb.heuristic( x.data() );
 const double look = since( t1 );
 const std::size_t me = Enc::relative ? std::min< std::size_t >( m , 1000 )
                                      : m;
 for( std::size_t i = 0 ; i < me ; ++i ) {
  if constexpr( Enc::relative )
   mh += pdb.exact( q[ i ].data() );
  else
   mh += pdb.heuristic( q[ i ].data() );
  md += manhattan( q[ i ] );
  }
 std::cout.precision( 6 );
 std::cout << opt[ "heur" ] << ',' << opt[ "enc" ] << ','
           << opt[ "compress" ] << ',' << opt[ "blank" ] << ','
           << ( pdb.reflected() ? 1 : 0 ) << ',' << entries << ','
           << pdb.bytes() << ',' << 8.0 * double( pdb.bytes() ) /
                                    double( entries ) << ','
           << ( loaded ? "loaded" : "built" ) << ',' << build << ','
           << mh / double( me ) << ',' << md / double( me ) << ','
           << 1e9 * look / double( m ) << ',' << ( sink & 1 ) << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the mode with the encoding Enc

template< class Enc >
static void run( void )
{
 if( opt[ "mode" ] == "search" ) {
  if( opt[ "puzzle" ] == "24" )
   search< 5 , 5 , Enc >();
  else
   if( opt[ "puzzle" ] == "15" )
    search< 4 , 4 , Enc >();
   else
    throw( std::invalid_argument( "unknown puzzle " + opt[ "puzzle" ] ) );
  }
 else
  if( opt[ "mode" ] != "table" )
   throw( std::invalid_argument( "unknown mode " + opt[ "mode" ] ) );
  else
   if( opt[ "heur" ] == "6666" )
    table< 5 , 5 , Enc >();
   else
    table< 4 , 4 , Enc >();
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 try {
  bool header = false;
  for( int i = 1 ; i < argc ; ++i ) {
   std::string a = argv[ i ];
   if( a == "--header" ) {
    header = true;
    continue;
    }
   if( ( a.rfind( "--" , 0 ) != 0 ) || ! opt.count( a.substr( 2 ) ) ||
       ( i + 1 >= argc ) )
    throw( std::invalid_argument( "bad option " + a ) );
   opt[ a.substr( 2 ) ] = argv[ ++i ];
   }
  if( header )
   std::cout << ( opt[ "mode" ] == "table"
                  ? "heur,enc,compress,blank,reflect,entries,bytes,bits,"
                    "source,seconds,mean_h,mean_manhattan,ns_lookup,sink"
                  : "instance,heur,enc,compress,blank,reflect,algo,threads,"
                    "status,cost,expanded,reexpanded,generated,seconds" )
             << std::endl;

  const std::string & e = opt[ "enc" ];
  if( e == "byte" )
   run< ByteEntries >();
  else
   if( e == "nibble" )
    run< NibbleEntries >();
   else
    if( e == "mod3" )
     run< Mod3Entries >();
    else
     if( e == "bit" )
      run< ParityBitEntries >();
     else
      if( e == "min" )
       run< MinNibbleEntries >();
      else
       #ifdef PBEFS_STRUCTURES_ADS
        if( e == "sparse" )
         run< SparseEntries >();
        else
       #endif
         throw( std::invalid_argument( "unknown encoding " + e ) );
  }
 catch( std::exception & ex ) {
  std::cerr << "PDBBench: " << ex.what() << std::endl;
  return( 1 );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File PDBBench.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
