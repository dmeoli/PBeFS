/*--------------------------------------------------------------------------*/
/*-------------------------- File TraceSearch.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Records the trace of a search [see TraceRecorder.h]: it runs one search
 * on one instance with a TraceRecorder and writes the events, one JSON
 * object per line after a header line, with a summary of the search on the
 * standard error.
 *
 *     TraceSearch [options]
 *
 *     --domain   8puzzle | 15puzzle | grid | queens | rubik | hanoi |
 *                pancake | topspin | knapsack               (8puzzle)
 *     --heur     h1 | h2 (8puzzle), h2 | 555 (15puzzle), none | korf
 *                (rubik), none | pdb (hanoi)                (h2, none)
 *     --algo     astar | idastar | spa | hda | pbnf | pase | farm for the
 *                Domains, bb | bb-shared | bb-dist for knapsack
 *                                                           (astar, bb)
 *     --threads  <workers of the parallel searches>         (1)
 *     --board    <tiles of the start, the blank 0, row by row>
 *     --instance <line of data/15-puzzle/korf100.txt or of
 *                data/rubik/korf10.txt>                     (79, 1)
 *     --walk     <moves of a random walk from the goal, in place of
 *                --board or --instance>                     (none)
 *     --size     <n of queens, items of knapsack, side of grid,
 *                disks of hanoi, 8 or 10, pancakes, 10 or 14,
 *                tokens of topspin, 12>          (8, 40, 64, 8, 10, 12)
 *     --seed     <seed of --walk, of grid and of knapsack>  (1)
 *     --data     <directory of the data>                    (data)
 *     --cache    <directory of the tables of korf>          (none)
 *     --cap      <events recorded, 0 none>                  (1000000)
 *     --every    <one event recorded every k>               (1)
 *     --events   <kinds recorded, comma-separated, or all>  (all)
 *     --limit    <nodes stored, expansions of IDA*, evaluations of B&B,
 *                0 none>                                    (0)
 *     --out      <file of the trace, - for the standard output>  (-)
 *
 * The 8-puzzle starts by default from the board of figure 3.28 of Russell
 * and Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition,
 * 2010), whose optimal path costs 26, with the goal of Korf (the blank in
 * the top-left corner); h1 is the number of the misplaced tiles, h2 the
 * Manhattan distance, 555 the additive tables of the partition 5-5-5. The
 * grid is a random map of side --size with a quarter of the cells blocked
 * and 8 moves, from the bottom-left to the bottom-right corner; the queens
 * are placed column by column; the cube with --heur korf takes the tables
 * of Korf (1997), 82 megabytes, built in minutes or read from --cache, and
 * an instance of korf10.txt asks for billions of expansions, which --limit
 * bounds; the Towers of Hanoi have four pegs, all the disks on the first
 * at the start (or --walk random moves from the goal) and on the last at
 * the goal, and --heur pdb takes the additive tables of the smaller and
 * the larger half of the disks; the knapsack is a random instance of the
 * weakly correlated class of Pisinger with weights up to 1000; the
 * pancakes are a random stack of --seed with the gap heuristic, and TopSpin
 * a ring of 12 tokens with a turnstile of 4, --walk random moves from the
 * identity, with the heuristic of its breakpoints. The ids of the trace
 * are the states themselves where they are a word (the packed boards of
 * the puzzles, the cells of the grid, two bits per disk of the pegs of
 * Hanoi), so that a reader decodes them, else their hashes in the Domain,
 * and for the knapsack the number of the items fixed (in the top byte) and
 * those taken among them (in the bits); a last line, after the events,
 * gives the result: for a Domain its status, cost, path (by the ids) and
 * counters, and the instance (the start of a puzzle, the map of the grid,
 * the cubies of the cube, the stack of pancakes, the ring of TopSpin), and
 * for the cube, the pancakes and TopSpin, whose ids are hashes, also the
 * states of the path (and for the cube the moves between them); for the
 * knapsack the value found, the optimum, the counters and the instance
 * (profits, weights and capacity). The kinds of --events are those that
 * the trace names (expand, generate, duplicate, push, pop, message,
 * threshold, solution, evaluate, branch, prune). PBNF takes the
 * sliding-tile puzzles and the grid, PA*SE the sliding-tile puzzles.
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

#include <cstdint>

#include <fstream>

#include <iostream>

#include <map>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/FarmAStar.h"

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/IDAStar.h"

#include "PBeFS/astar/NQueens.h"

#include "PBeFS/astar/Pancake.h"

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/PermutationPDB.h"

#include "PBeFS/astar/RubiksCube.h"

#include "PBeFS/astar/SafePBNF.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/SPAStar.h"

#include "PBeFS/astar/TopSpin.h"

#include "PBeFS/astar/TowersOfHanoi.h"

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/bnb/Knapsack.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "PBeFS/core/TraceRecorder.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------- TYPES AND GLOBAL VARIABLES -----------------------*/
/*--------------------------------------------------------------------------*/

// the options of the command line, by name
// the instance of the run, as JSON, for the result line

static std::string instance_json = "null";

static std::map< std::string , std::string > opt = {
 { "domain" , "8puzzle" } , { "heur" , "" } , { "algo" , "" } ,
 { "threads" , "1" } , { "board" , "" } , { "instance" , "" } ,
 { "walk" , "" } , { "size" , "" } , { "seed" , "1" } ,
 { "data" , "data" } , { "cache" , "" } , { "cap" , "1000000" } ,
 { "every" , "1" } , { "events" , "all" } , { "limit" , "0" } ,
 { "out" , "-" } };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the 8-puzzle with h1, the number of the tiles out of their goal cell

class MisplacedTiles {

 public:

 using state_type = EightPuzzle::state_type;
 using cost_type = int;

 explicit MisplacedTiles( const EightPuzzle::board_type & start )
  : f_base( start ) ,
    f_goal( EightPuzzle::pack( EightPuzzle::korf_goal() ) ) {}

 state_type initial_state( void ) const { return( f_base.initial_state() ); }

 bool is_goal( state_type s ) const { return( f_base.is_goal( s ) ); }

 std::uint64_t hash( state_type s ) const { return( f_base.hash( s ) ); }

 int heuristic( state_type s ) const {
  int h = 0;
  for( int c = 0 ; c < 9 ; ++c )
   if( const int t = EightPuzzle::tile( s , c ) ;
       t && ( t != EightPuzzle::tile( f_goal , c ) ) )
    ++h;
  return( h );
  }

 template< class F >
 void for_each_successor( state_type s , int , F && f ) const {
  f_base.for_each_successor( s , 0 , [ & ]( state_type t , int c , int ) {
   f( t , c , heuristic( t ) );
   } );
  }

 private:

 EightPuzzle f_base;                             // the puzzle
 state_type f_goal;                              // the goal, packed

 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the id of a state in a trace: the state if it is a word, else its hash
// in the Domain

template< class D >
struct HashId {
 const D * d;
 std::uint64_t operator()( const state_t< D > & s ) const {
  if constexpr( std::is_integral_v< state_t< D > > )
   return( std::uint64_t( s ) );
  else
   return( d->hash( s ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the id of a node of the knapsack: the number k of the items fixed in
// the top byte and the items taken among them in the bits 0 to k - 1, so
// that a reader finds the depth, the items and the parent of a node (the
// same id with the bit k - 1 cleared and k - 1 fixed); 56 items at most

struct KnapsackId {
 std::uint64_t operator()( const Knapsack::node_type & n ) const {
  std::uint64_t id = std::uint64_t( n.k ) << 56;
  for( std::size_t i = 0 ; i < n.k ; ++i )
   if( n.x[ i ] )
    id |= std::uint64_t( 1 ) << i;
  return( id );
  }
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the value of the option name as an unsigned number

static std::uint64_t number( const std::string & name )
{
 return( std::stoull( opt[ name ] ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the recorder of the options, ids by id

template< class Id >
static TraceRecorder< Id > recorder( Id id )
{
 std::uint32_t mask = 0;
 if( opt[ "events" ] == "all" )
  mask = trace_all;
 else {
  std::istringstream in( opt[ "events" ] );
  for( std::string k ; std::getline( in , k , ',' ) ; ) {
   bool known = false;
   for( unsigned e = 0 ; e <= unsigned( TraceEvent::prune ) ; ++e )
    if( k == trace_name( TraceEvent( e ) ) ) {
     mask |= trace_bit( TraceEvent( e ) );
     known = true;
     }
   if( ! known )
    throw( std::invalid_argument( "unknown kind of event " + k ) );
   }
  }
 return( TraceRecorder< Id >( id , number( "cap" ) , mask ,
                              number( "every" ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// writes the trace of r to --out and the summary to the standard error

template< class R >
static void write( const R & r , const std::string & summary ,
                   const std::string & result = "" )
{
 const auto put = [ & ]( std::ostream & out ) {
  r.write( out );
  if( ! result.empty() )
   out << result << "\n";
  };
 if( opt[ "out" ] == "-" )
  put( std::cout );
 else {
  std::ofstream out( opt[ "out" ] );
  if( ! out )
   throw( std::runtime_error( "cannot write " + opt[ "out" ] ) );
  put( out );
  }
 std::cerr << summary << ", " << r.size() << " events recorded, "
           << r.dropped() << " dropped" << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cubies of a cube as a JSON list, 3 c + o at the corner positions and
// 2 e + o at the edge positions [see RubiksCube.h]

static std::string cube_json( const RubiksCube::state_type & s )
{
 std::ostringstream os;
 os << "[";
 for( std::size_t i = 0 ; i < s.size() ; ++i )
  os << ( i ? "," : "" ) << int( s[ i ] );
 os << "]";
 return( os.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the path of the cube, whose ids are hashes, as its states ("states") and
// as the moves between them in the notation of Singmaster ("moves")

static void cube_path( std::ostream & res ,
                       const std::vector< RubiksCube::state_type > & path )
{
 static const char * const turn[ 3 ] = { "" , "2" , "'" };
 res << ",\"states\":[";
 for( std::size_t i = 0 ; i < path.size() ; ++i )
  res << ( i ? "," : "" ) << cube_json( path[ i ] );
 res << "],\"moves\":[";
 for( std::size_t i = 1 ; i < path.size() ; ++i )
  for( int m = 0 ; m < RubiksCube::moves ; ++m )
   if( RubiksCube::apply( path[ i - 1 ] , m ) == path[ i ] ) {
    res << ( i > 1 ? ",\"" : "\"" ) << "URFDLB"[ m / 3 ] << turn[ m % 3 ]
        << "\"";
    break;
    }
 res << "]";
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the path of a Domain whose states are arrays of small numbers named by
// their hashes (the pancakes, TopSpin), as its states ("states")

template< class S >
static void array_path( std::ostream & res , const std::vector< S > & path )
{
 res << ",\"states\":[";
 for( std::size_t i = 0 ; i < path.size() ; ++i ) {
  res << ( i ? ",[" : "[" );
  for( std::size_t j = 0 ; j < path[ i ].size() ; ++j )
   res << ( j ? "," : "" ) << int( path[ i ][ j ] );
  res << "]";
  }
 res << "]";
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs --algo on the Domain d with --threads workers and writes its trace

template< class D >
static void run( const D & d )
{
 using R = TraceRecorder< HashId< D > >;
 using Open = default_open_t< D >;
 R rec = recorder( HashId< D >{ & d } );
 const std::string & a = opt[ "algo" ].empty() ? "astar" : opt[ "algo" ];
 ParallelOptions o;
 o.threads = unsigned( number( "threads" ) );
 o.node_limit = std::size_t( number( "limit" ) );
 const auto solve = [ & ]( auto && s ) {
  s.set_observer( rec );
  return( s.solve() );
  };
 Solution< state_t< D > , cost_t< D > > sol;
 if( a == "astar" )
  sol = solve( AStar< D , Open , R >( d , o ) );
 else
  if( a == "idastar" )
   sol = solve( IDAStar< D , R >( d , o ) );
  else
   if( a == "spa" )
    sol = solve( SPAStar< D , ADS::ThreadPool , Open , R >( d , o ) );
   else
    if( a == "hda" ) {
     const StateHashDistribution< D > z( d );
     sol = solve( HDAStar< D , StateHashDistribution< D > , ADS::ThreadPool ,
                           Open , R >( d , z , o ) );
     }
    else
     if( a == "farm" )
      sol = solve( FarmAStar< D , ADS::ThreadPool , Open , R >( d , o ) );
     else
      if( a == "pbnf" ) {
       if constexpr( NBlockDomain< D > ) {
        PBNFOptions po;
        static_cast< ParallelOptions & >( po ) = o;
        sol = solve( SafePBNF< D , ADS::ThreadPool , Open , R >( d , po ) );
        }
       else
        throw( std::invalid_argument( "PBNF takes the puzzles and the "
                                      "grid" ) );
       }
      else
       if( a == "pase" ) {
        if constexpr( PairwiseDomain< D > ) {
         PASEOptions po;
         static_cast< ParallelOptions & >( po ) = o;
         sol = solve( PASE< D , ADS::ThreadPool , R >( d , po ) );
         }
        else
         throw( std::invalid_argument( "PA*SE takes the puzzles" ) );
        }
       else
        throw( std::invalid_argument( "unknown --algo " + a ) );
 std::ostringstream os;
 os << a << ": " << ( sol.solved() ? "solved" : "not solved" ) << ", cost "
    << sol.cost << ", expanded " << sol.stats.expanded << ", generated "
    << sol.stats.generated;
 const HashId< D > id{ & d };
 const SearchStats & st = sol.stats;
 std::ostringstream res;
 res << "{\"result\":\"" << ( sol.solved() ? "solved"
     : sol.status == Status::unsolvable ? "unsolvable" : "limit" )
     << "\",\"cost\":" << sol.cost << ",\"path\":[";
 for( std::size_t i = 0 ; i < sol.path.size() ; ++i )
  res << ( i ? ",\"" : "\"" ) << std::hex << id( sol.path[ i ] ) << std::dec
      << "\"";
 res << "],\"stats\":{\"expanded\":" << st.expanded << ",\"reexpanded\":"
     << st.reexpanded << ",\"generated\":" << st.generated
     << ",\"duplicates\":" << st.duplicates << ",\"sent\":" << st.sent
     << ",\"messages\":" << st.messages << ",\"stored\":" << st.stored
     << ",\"iterations\":" << st.iterations << ",\"seconds\":"
     << st.seconds << ",\"expanded_by\":[";
 for( std::size_t i = 0 ; i < st.expanded_by.size() ; ++i )
  res << ( i ? "," : "" ) << st.expanded_by[ i ];
 res << "]}";
 if constexpr( std::is_same_v< state_t< D > , RubiksCube::state_type > )
  cube_path( res , sol.path );
 else
  if constexpr( ! std::is_integral_v< state_t< D > > &&
                requires( const state_t< D > & x ) {
                 x.size(); int( x[ 0 ] ); } )
   array_path( res , sol.path );
 res << ",\"instance\":" << instance_json << "}";
 write( rec , os.str() , res.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the start of a sliding-tile puzzle as the instance of the result

template< class B >
static void tiles_instance( const B & b )
{
 std::ostringstream os;
 os << "{\"board\":[";
 for( std::size_t i = 0 ; i < b.size() ; ++i )
  os << ( i ? "," : "" ) << int( b[ i ] );
 os << "]}";
 instance_json = os.str();
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the start of a sliding-tile puzzle: --board, or --walk random moves of
// the blank from the goal, or the board b

template< class P >
static typename P::board_type tiles_start( typename P::board_type b )
{
 if( ! opt[ "board" ].empty() ) {
  std::istringstream in( opt[ "board" ] );
  for( auto & t : b ) {
   int x;
   if( ! ( in >> x ) )
    throw( std::invalid_argument( "--board wants the tiles" ) );
   t = std::uint8_t( x );
   }
  }
 else
  if( ! opt[ "walk" ].empty() ) {
   const P g( P::korf_goal() );
   std::mt19937_64 gen( number( "seed" ) );
   auto s = g.initial_state();
   for( std::uint64_t i = number( "walk" ) ; i > 0 ; --i ) {
    std::vector< typename P::state_type > next;
    g.for_each_successor( s , 0 , [ & ]( const auto & t , int , int ) {
     next.push_back( t );
     } );
    s = next[ random_below( gen , next.size() ) ];
    }
   b = P::unpack( s );
   }
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the line --instance (counted from 1) of a file of the data, without the
// comments

static std::string data_line( const std::string & file , std::size_t def )
{
 std::ifstream in( opt[ "data" ] + "/" + file );
 if( ! in )
  throw( std::runtime_error( "cannot read " + opt[ "data" ] + "/" + file ) );
 const std::size_t k = opt[ "instance" ].empty() ? def
                       : std::size_t( number( "instance" ) );
 std::size_t i = 0;
 for( std::string line ; std::getline( in , line ) ; )
  if( ! line.empty() && ( line[ 0 ] != '#' ) && ( ++i == k ) )
   return( line );
 throw( std::invalid_argument( "no instance " + std::to_string( k ) +
                               " in " + file ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the Towers of Hanoi of D disks with four pegs, from all the disks on peg
// 0 (or --walk random moves from the goal) to all of them on peg 3

template< int D >
static void run_hanoi( void )
{
 const std::uint64_t goal = TowersOfHanoi< D >::mask();
 std::uint64_t s = 0;
 if( ! opt[ "walk" ].empty() ) {
  const TowersOfHanoi< D > g( goal , goal );
  std::mt19937_64 gen( number( "seed" ) );
  s = goal;
  for( std::uint64_t i = number( "walk" ) ; i > 0 ; --i ) {
   std::vector< std::uint64_t > next;
   g.for_each_successor( s , 0 , [ & ]( std::uint64_t t , int , int ) {
    next.push_back( t );
    } );
   s = next[ random_below( gen , next.size() ) ];
   }
  }
 instance_json = "{\"disks\":" + std::to_string( D ) + ",\"start\":\"" +
                 TowersOfHanoi< D >::write( s ) + "\",\"goal\":\"" +
                 TowersOfHanoi< D >::write( goal ) + "\"}";
 const TowersOfHanoi< D > d( s , goal );
 const std::string & h = opt[ "heur" ];
 if( h.empty() || ( h == "none" ) )
  run( d );
 else
  if( h == "pdb" ) {
   const HanoiPDB pdb( { unsigned( D / 2 ) , unsigned( D - D / 2 ) } , goal ,
                       0 , opt[ "cache" ] );
   run( PDBDomain< TowersOfHanoi< D > , HanoiPDB >( d , pdb ) );
   }
  else
   throw( std::invalid_argument( "the towers of Hanoi take none and pdb" ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the pancakes: N of them in a uniformly random stack of --seed, with the
// gap heuristic of the Domain

template< int N >
static void run_pancake( void )
{
 const auto s = PancakePuzzle< N >::random( number( "seed" ) );
 std::ostringstream os;
 os << "{\"pancakes\":" << N << ",\"start\":[";
 for( int i = 0 ; i < N ; ++i )
  os << ( i ? "," : "" ) << int( s[ i ] );
 os << "]}";
 instance_json = os.str();
 run( PancakePuzzle< N >( s ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// TopSpin: N tokens and a turnstile of 4, --walk random moves from the
// identity (a uniform ring if none), with the heuristic of the Domain

template< int N >
static void run_topspin( void )
{
 const int walk = opt[ "walk" ].empty() ? 0 : int( number( "walk" ) );
 const auto s = TopSpin< N , 4 >::random( number( "seed" ) , walk );
 std::ostringstream os;
 os << "{\"tokens\":" << N << ",\"turnstile\":4,\"start\":[";
 for( int i = 0 ; i < N ; ++i )
  os << ( i ? "," : "" ) << int( s[ i ] );
 os << "]" << ( opt[ "walk" ].empty() ? "" : ",\"walk\":" + opt[ "walk" ] )
    << ",\"seed\":" << opt[ "seed" ] << "}";
 instance_json = os.str();
 run( TopSpin< N , 4 >( s ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the Domain of --domain and --heur, run

static void run_domain( void )
{
 const std::string & dom = opt[ "domain" ];
 const std::string & h = opt[ "heur" ];
 if( dom == "8puzzle" ) {
  const auto b = tiles_start< EightPuzzle >( { 7 , 2 , 4 , 5 , 0 , 6 , 8 ,
                                               3 , 1 } );
  tiles_instance( b );
  if( h == "h1" )
   run( MisplacedTiles( b ) );
  else
   if( h.empty() || ( h == "h2" ) )
    run( EightPuzzle( b ) );
   else
    throw( std::invalid_argument( "the 8-puzzle takes h1 and h2" ) );
  return;
  }
 if( dom == "15puzzle" ) {
  FifteenPuzzle::board_type b{};
  if( opt[ "board" ].empty() && opt[ "walk" ].empty() ) {
   std::istringstream in( data_line( "15-puzzle/korf100.txt" , 79 ) );
   for( auto & t : b ) {
    int x;
    in >> x;
    t = std::uint8_t( x );
    }
   }
  b = tiles_start< FifteenPuzzle >( b );
  tiles_instance( b );
  if( h.empty() || ( h == "h2" ) )
   run( FifteenPuzzle( b ) );
  else
   if( h == "555" ) {
    const AdditivePDB< 4 , 4 > pdb( partition_15_5_5_5() ,
                                    FifteenPuzzle::korf_goal() );
    run( SlidingTilePDB< 4 , 4 >( b , pdb ) );
    }
   else
    throw( std::invalid_argument( "the 15-puzzle takes h2 and 555" ) );
  return;
  }
 if( dom == "grid" ) {
  const int n = opt[ "size" ].empty() ? 64 : int( number( "size" ) );
  const auto m = GridMap::random( n , n , 0.25 , number( "seed" ) );
  std::ostringstream os;
  os << "{\"width\":" << n << ",\"height\":" << n << ",\"start\":"
     << n * ( n - 1 ) << ",\"goal\":" << n * n - 1 << ",\"rows\":[";
  for( int y = 0 ; y < n ; ++y ) {
   os << ( y ? ",\"" : "\"" );
   for( int x = 0 ; x < n ; ++x )
    os << ( m.cost( x , y ) ? '.' : '@' );
   os << "\"";
   }
  os << "]}";
  instance_json = os.str();
  run( Grid< 8 >( m , 0 , n - 1 , n - 1 , n - 1 ) );
  return;
  }
 if( dom == "queens" ) {
  run( NQueens( opt[ "size" ].empty() ? 8 : int( number( "size" ) ) ) );
  return;
  }
 if( dom == "rubik" ) {
  const RubiksCube::state_type s = opt[ "walk" ].empty()
   ? RubiksCube::scramble( data_line( "rubik/korf10.txt" , 1 ) )
   : RubiksCube::random_walk( int( number( "walk" ) ) , number( "seed" ) );
  const RubiksCube d( s );
  instance_json = "{\"cubies\":" + cube_json( s ) + ( opt[ "walk" ].empty()
   ? ",\"korf\":" + ( opt[ "instance" ].empty() ? "1" : opt[ "instance" ] )
   : ",\"walk\":" + opt[ "walk" ] + ",\"seed\":" + opt[ "seed" ] ) + "}";
  if( h.empty() || ( h == "none" ) )
   run( d );
  else
   if( h == "korf" ) {
    const RubikPDBHeuristic pdb( d , 0 , opt[ "cache" ] );
    run( PDBDomain< RubiksCube , RubikPDBHeuristic >( d , pdb ) );
    }
   else
    throw( std::invalid_argument( "the cube takes none and korf" ) );
  return;
  }
 if( dom == "pancake" ) {
  if( opt[ "size" ].empty() || ( opt[ "size" ] == "10" ) )
   run_pancake< 10 >();
  else
   if( opt[ "size" ] == "14" )
    run_pancake< 14 >();
   else
    throw( std::invalid_argument( "the pancakes take 10 or 14" ) );
  return;
  }
 if( dom == "topspin" ) {
  if( opt[ "size" ].empty() || ( opt[ "size" ] == "12" ) )
   run_topspin< 12 >();
  else
   throw( std::invalid_argument( "TopSpin takes 12 tokens" ) );
  return;
  }
 if( dom == "hanoi" ) {
  if( opt[ "size" ].empty() || ( opt[ "size" ] == "8" ) )
   run_hanoi< 8 >();
  else
   if( opt[ "size" ] == "10" )
    run_hanoi< 10 >();
   else
    throw( std::invalid_argument( "the towers of Hanoi take 8 or 10 "
                                  "disks" ) );
  return;
  }
 throw( std::invalid_argument( "unknown --domain " + dom ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack, by --algo

static void run_knapsack( void )
{
 using R = TraceRecorder< KnapsackId >;
 const auto kp = Knapsack::random( 1 , opt[ "size" ].empty() ? 40
                                    : std::size_t( number( "size" ) ) ,
                                   1000 , number( "seed" ) );
 R rec = recorder( KnapsackId{} );
 const std::string & a = opt[ "algo" ].empty() ? "bb" : opt[ "algo" ];
 ParallelBBOptions o;
 o.threads = unsigned( number( "threads" ) );
 o.node_limit = std::size_t( number( "limit" ) );
 const auto solve = [ & ]( auto && s ) {
  s.set_observer( rec );
  return( s.solve() );
  };
 BBSolution< Knapsack::node_type , double > sol;
 if( a == "bb" )
  sol = solve( BranchAndBoundSearch< Knapsack , R >( kp , o ) );
 else
  if( a == "bb-shared" )
   sol = solve( SharedBranchAndBound< Knapsack , ADS::ThreadPool , R >( kp ,
                                                                       o ) );
  else
   if( a == "bb-dist" )
    sol = solve( DistributedBranchAndBound< Knapsack , ADS::ThreadPool ,
                                            R >( kp , o ) );
   else
    throw( std::invalid_argument( "the knapsack takes bb, bb-shared and "
                                  "bb-dist" ) );
 std::ostringstream os;
 os << a << ": value " << sol.value << " (optimum " << kp.optimum()
    << "), evaluated " << sol.stats.evaluated << ", branched "
    << sol.stats.branched;
 std::ostringstream res;
 res << "{\"result\":\"" << ( sol.value == double( kp.optimum() ) ? "optimal"
                                                                : "limit" )
     << "\",\"value\":" << sol.value << ",\"optimum\":" << kp.optimum()
     << ",\"stats\":{\"evaluated\":" << sol.stats.evaluated
     << ",\"branched\":" << sol.stats.branched << "},\"instance\":{\"p\":[";
 for( std::size_t i = 0 ; i < kp.p.size() ; ++i )
  res << ( i ? "," : "" ) << kp.p[ i ];
 res << "],\"w\":[";
 for( std::size_t i = 0 ; i < kp.w.size() ; ++i )
  res << ( i ? "," : "" ) << kp.w[ i ];
 res << "],\"cap\":" << kp.cap << ",\"seed\":" << opt[ "seed" ] << "}}";
 write( rec , os.str() , res.str() );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 try {
  for( int i = 1 ; i < argc ; ++i ) {
   const std::string a = argv[ i ];
   if( ( a.rfind( "--" , 0 ) != 0 ) || ! opt.count( a.substr( 2 ) ) ||
       ( i + 1 >= argc ) )
    throw( std::invalid_argument( "bad option " + a ) );
   opt[ a.substr( 2 ) ] = argv[ ++i ];
   }
  if( opt[ "domain" ] == "knapsack" )
   run_knapsack();
  else
   run_domain();
  }
 catch( const std::exception & e ) {
  std::cerr << "TraceSearch: " << e.what() << std::endl;
  return( 1 );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File TraceSearch.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
