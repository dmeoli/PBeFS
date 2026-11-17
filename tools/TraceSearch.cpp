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
 *     --domain   8puzzle | 15puzzle | grid | queens | rubik | knapsack
 *                                                           (8puzzle)
 *     --heur     h1 | h2 (8puzzle), h2 | 555 (15puzzle), none | korf
 *                (rubik)                                    (h2, none)
 *     --algo     astar | idastar | spa | hda | pbnf | pase | farm for the
 *                Domains, bb | bb-shared | bb-dist for knapsack
 *                                                           (astar, bb)
 *     --threads  <workers of the parallel searches>         (1)
 *     --board    <tiles of the start, the blank 0, row by row>
 *     --instance <line of data/15-puzzle/korf100.txt or of
 *                data/rubik/korf10.txt>                     (79, 1)
 *     --walk     <moves of a random walk from the goal, in place of
 *                --board or --instance>                     (none)
 *     --size     <n of queens, items of knapsack, side of grid>
 *                                                           (8, 40, 64)
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
 * bounds; the knapsack is a random instance of the weakly correlated
 * class of Pisinger with weights up to 1000. The ids of the trace are the
 * states themselves where they are a word (the packed boards of the
 * puzzles, the cells of the grid), so that a reader decodes them, else
 * their hashes in the Domain, and for the knapsack the items fixed and the
 * profit of the node; a last line, after the events, gives the result of
 * a Domain: its status, cost, path (by the ids) and counters, and the
 * instance (the start of a puzzle, the map of the grid); the kinds of --events are those that
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

#include "PBeFS/astar/PASE.h"

#include "PBeFS/astar/PatternDatabase.h"

#include "PBeFS/astar/PermutationPDB.h"

#include "PBeFS/astar/RubiksCube.h"

#include "PBeFS/astar/SafePBNF.h"

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/SPAStar.h"

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
// the id of a node of the knapsack: the items fixed and the profit

struct KnapsackId {
 std::uint64_t operator()( const Knapsack::node_type & n ) const {
  return( ( std::uint64_t( n.k ) << 40 ) ^ std::uint64_t( n.profit ) );
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
 res << "]},\"instance\":" << instance_json << "}";
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
 write( rec , os.str() );
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
