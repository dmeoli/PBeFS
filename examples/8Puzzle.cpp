/*--------------------------------------------------------------------------*/
/*---------------------------- File 8Puzzle.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The 8-puzzle with A*, after section 3.6 of Russell and Norvig,
 * "Artificial Intelligence: A Modern Approach" (3rd edition, 2010).
 *
 * Usage:   8Puzzle [-h h1|h2|seq] [-g aima|nilsson] [<start>|<level>]
 *          8Puzzle [-g aima|nilsson] -t [<samples>]
 *
 * The first form solves one instance and prints its path. The goal is the
 * one of figure 3.28 of the book (the blank in the top-left corner, then
 * 1 , ... , 8 row by row) or, with -g nilsson, the "snail" goal of Nilsson
 * (1 2 3 / 8 _ 4 / 7 6 5). The start is given as the nine tiles row by row,
 * 0 being the blank, or as a level, i.e., easy, medium, hard or worst: the
 * first state found at distance 5, 10 or 20, or at the largest distance,
 * by a breadth-first search from the goal. Without a start, it is the
 * instance of figure 3.28, whose optimal solution has 26 moves. The
 * heuristic is the number of misplaced tiles (h1), the Manhattan distance
 * (h2, the default) or the sequence score of Nilsson (seq), which is not
 * admissible and may therefore give a longer path.
 *
 * The second form repeats the experiment of figure 3.29 of the book: for
 * each even depth from 2 to 24 it draws <samples> (default 100) random
 * states at that distance from the goal, solves them with h1 and with h2,
 * and prints the average number of nodes generated and the effective
 * branching factor b*, i.e., the solution of N + 1 = 1 + b* + ... +
 * ( b* )^d. The figure is for A* without a closed list, while this is A*
 * on the graph, and the number of the samples is not the same; hence, the
 * counts are close to those of the book but not equal.
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

#include <cctype>

#include <cstring>

#include <iomanip>

#include <iostream>

#include <random>

#include <string>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/SequenceScore.h"

#include "PBeFS/astar/SlidingTile.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

using Board = EightPuzzle::board_type;
using State = EightPuzzle::state_type;

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MisplacedTiles --------------------------*/
/*--------------------------------------------------------------------------*/
/// the 8-puzzle with h1, the number of the tiles out of their goal cell
/** A Domain written on top of another one: everything but the heuristic
 * comes from EightPuzzle. */

class MisplacedTiles {

 public:

 using state_type = State;                       ///< a packed board
 using cost_type = int;                          ///< the cost type

 MisplacedTiles( const Board & start , const Board & goal )
  : f_base( start , goal ) , f_goal( EightPuzzle::pack( goal ) ) {}

 [[nodiscard]] State initial_state( void ) const noexcept {
  return( f_base.initial_state() );
  }

 [[nodiscard]] bool is_goal( State s ) const noexcept {
  return( f_base.is_goal( s ) );
  }

 [[nodiscard]] std::uint64_t hash( State s ) const noexcept {
  return( f_base.hash( s ) );
  }

 [[nodiscard]] int heuristic( State s ) const noexcept {
  int h = 0;
  for( int c = 0 ; c < 9 ; ++c )
   if( const int t = EightPuzzle::tile( s , c ) ;
       t && ( t != EightPuzzle::tile( f_goal , c ) ) )
    ++h;
  return( h );
  }

 template< class F >
 void for_each_successor( State s , int , F && f ) const {
  f_base.for_each_successor( s , 0 , [ & ]( State t , int c , int ) {
   f( t , c , heuristic( t ) );
   } );
  }

 private:

 EightPuzzle f_base;                             ///< the puzzle
 State f_goal;                                   ///< the goal, packed

 };  // end( class( MisplacedTiles ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the states at distance 0 , 1 , ... from goal, layer by layer, in the
// order of a breadth-first search from it

static std::vector< std::vector< State > > layers( const Board & goal )
{
 const EightPuzzle d( goal , goal );
 std::unordered_map< State , int > seen;
 std::vector< std::vector< State > > lay( 1 , { d.initial_state() } );
 seen.emplace( d.initial_state() , 0 );
 while( ! lay.back().empty() ) {
  std::vector< State > next;
  for( State s : lay.back() )
   d.for_each_successor( s , 0 , [ & ]( State t , int , int ) {
    if( seen.emplace( t , int( lay.size() ) ).second )
     next.push_back( t );
    } );
  lay.push_back( std::move( next ) );
  }
 lay.pop_back();
 return( lay );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the board of a packed state

static Board unpack( State s )
{
 Board b{};
 for( int c = 0 ; c < 9 ; ++c )
  b[ c ] = std::uint8_t( EightPuzzle::tile( s , c ) );
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// b* such that 1 + b* + ... + ( b* )^d = n + 1, by bisection

static double branching( double n , int d )
{
 double lo = 1 , hi = 10;
 for( int it = 0 ; it < 100 ; ++it ) {
  const double b = ( lo + hi ) / 2;
  double sum = 1 , p = 1;
  for( int k = 1 ; k <= d ; ++k )
   sum += p *= b;
  ( sum < n + 1 ? lo : hi ) = b;
  }
 return( ( lo + hi ) / 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the experiment of figure 3.29: samples random states per even depth

static int table( const Board & goal , int samples )
{
 const auto lay = layers( goal );
 std::mt19937_64 gen( 1 );
 std::cout << " d    nodes h1    nodes h2     b* h1   b* h2" << std::endl;
 for( int d = 2 ; ( d <= 24 ) && ( d < int( lay.size() ) ) ; d += 2 ) {
  double n1 = 0 , n2 = 0;
  for( int k = 0 ; k < samples ; ++k ) {
   const Board start = unpack( lay[ d ][ random_below( gen ,
                                                        lay[ d ].size() ) ] );
   n1 += double( astar( MisplacedTiles( start , goal ) ).stats.generated );
   n2 += double( astar( EightPuzzle( start , goal ) ).stats.generated );
   }
  n1 /= samples;
  n2 /= samples;
  std::cout << std::setw( 2 ) << d << std::fixed << std::setprecision( 0 )
            << std::setw( 12 ) << n1 << std::setw( 12 ) << n2
            << std::setprecision( 2 ) << std::setw( 10 )
            << branching( n1 , d ) << std::setw( 8 ) << branching( n2 , d )
            << std::endl;
  }
 return( 0 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// prints the path found by a search and its counters

template< class S >
static int report( const S & sol )
{
 if( ! sol.solved() ) {
  std::cout << "no path to the goal" << std::endl;
  return( 1 );
  }
 int k = 0;
 for( State s : sol.path ) {
  std::cout << "move " << k++ << std::endl;
  EightPuzzle::print( std::cout , s );
  std::cout << std::endl;
  }
 std::cout << "moves:     " << sol.cost << std::endl
           << "expanded:  " << sol.stats.expanded << std::endl
           << "generated: " << sol.stats.generated << std::endl
           << "seconds:   " << sol.stats.seconds << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 std::string heur = "h2" , goal_name = "aima" , start_arg;
 int samples = 0;
 try {
  for( int i = 1 ; i < argc ; ++i )
   if( ! std::strcmp( argv[ i ] , "-h" ) && ( i + 1 < argc ) )
    heur = argv[ ++i ];
   else
    if( ! std::strcmp( argv[ i ] , "-g" ) && ( i + 1 < argc ) )
     goal_name = argv[ ++i ];
    else
     if( ! std::strcmp( argv[ i ] , "-t" ) )
      samples = ( i + 1 < argc ) ? std::stoi( argv[ ++i ] ) : 100;
     else
      start_arg = argv[ i ];

  if( ( goal_name != "aima" ) && ( goal_name != "nilsson" ) )
   throw( std::invalid_argument( "the goal is aima or nilsson" ) );
  const Board goal = goal_name == "aima"
                      ? Board{ 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 }
                      : Board{ 1 , 2 , 3 , 8 , 0 , 4 , 7 , 6 , 5 };
  if( samples > 0 )
   return( table( goal , samples ) );

  Board start{ 7 , 2 , 4 , 5 , 0 , 6 , 8 , 3 , 1 };   // figure 3.28
  const char * levels[] = { "easy" , "medium" , "hard" , "worst" };
  const int depth[] = { 5 , 10 , 20 , -1 };
  bool level = false;
  for( int l = 0 ; l < 4 ; ++l )
   if( start_arg == levels[ l ] ) {
    const auto lay = layers( goal );
    start = unpack( lay[ depth[ l ] < 0 ? lay.size() - 1
                                        : std::size_t( depth[ l ] ) ][ 0 ] );
    level = true;
    }
  if( ! ( level || start_arg.empty() ) ) {
   int k = 0;
   for( char c : start_arg )
    if( std::isdigit( static_cast< unsigned char >( c ) ) && ( k < 9 ) )
     start[ k++ ] = std::uint8_t( c - '0' );
   if( k != 9 )
    throw( std::invalid_argument( "the start needs nine tiles" ) );
   }

  if( heur == "h1" )
   return( report( astar( MisplacedTiles( start , goal ) ) ) );
  if( heur == "h2" )
   return( report( astar( EightPuzzle( start , goal ) ) ) );
  if( heur == "seq" )
   return( report( astar( SequenceScore< 3 , 3 >( start , goal ) ) ) );
  throw( std::invalid_argument( "the heuristic is h1, h2 or seq" ) );
  }
 catch( std::exception & e ) {
  std::cerr << "8Puzzle: " << e.what() << std::endl;
  return( 2 );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File 8Puzzle.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
