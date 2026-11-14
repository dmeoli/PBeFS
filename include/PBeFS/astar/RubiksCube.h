/*--------------------------------------------------------------------------*/
/*--------------------------- File RubiksCube.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class RubiksCube, which implements the Domain concept
 * [see Domain.h] for Rubik's cube in the half-turn metric, and for the
 * pattern databases of Korf (1997): the 8 corners and two sets of 6 edges.
 *
 * A state is the cubie and the orientation at each of the 8 corner and 12
 * edge positions, in the order of Kociemba: the corners URF, UFL, ULB,
 * UBR, DFR, DLF, DBL, DRB and the edges UR, UF, UL, UB, DR, DF, DL, DB, FR,
 * FL, BL, BR. A move turns a face (U, R, F, D, L, B) by a quarter turn
 * clockwise, a half turn or a quarter turn counterclockwise, at cost 1,
 * and is the permutation of the positions with the turns of the cubies of
 * the cubie model of Kociemba: the cubie that comes to the position i is
 * the one at \f$ \mathrm{from}[ i ] \f$, its orientation grown by
 * \f$ \mathrm{turn}[ i ] \f$ modulo 3 (corners) or 2 (edges), so that the
 * cube is a permutation puzzle [see PermutationPDB.h] of two groups whose
 * orientations sum to 0.
 *
 * The operator of the turn t (0, 1, 2) of the face f is 3 f + t. The
 * successors of an OperatorDomain skip the turns of the face just turned,
 * which lead back to the parent or to a sibling, and those of a tree
 * search also the turns of the face opposite to the one just turned when
 * it comes first in the order U, R, F, D, L, B, since two turns of
 * opposite faces commute: the branching factor is then 13.34847 after the
 * first move, and the tree of depth d has as many nodes as in table 1 of
 * Korf (1997).
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_RubiksCube
 #define __PBeFS_RubiksCube
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cstddef>

#include <cstdint>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Permutation.h"

#include "PBeFS/astar/PermutationPDB.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS RubiksCube ----------------------------*/
/*--------------------------------------------------------------------------*/
/// Rubik's cube in the half-turn metric
/** The state holds 3 c + o at the corner positions 0 to 7 and 2 e + o at
 * the edge positions 8 to 19, c and e being the cubies and o their
 * orientations. The heuristic is 0: the pattern databases are given by
 * RubikPDBHeuristic and PDBDomain [see PermutationPDB.h]. The features of
 * Zobrist hashing are the pairs ( position , value ). */

class RubiksCube {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::array< std::uint8_t , 20 >;  ///< a cube
 using cost_type = int;                               ///< the cost type

 static constexpr int moves = 18;                    ///< the moves

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the start cube

 explicit RubiksCube( const state_type & start = solved() )
  : f_start( start ) {}

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 [[nodiscard]] static bool is_goal( const state_type & s ) noexcept {
  return( s == solved() );
  }

 [[nodiscard]] static std::uint64_t hash( const state_type & s ) noexcept {
  return( bytes_hash( s.data() , s.size() ) );
  }

 [[nodiscard]] static cost_type heuristic( const state_type & ) noexcept {
  return( 0 );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f )
  const {
  for_each_successor( s , h , no_op , [ & ]( const state_type & t ,
                                             cost_type c , cost_type ht ,
                                             std::uint8_t ) {
   f( t , c , ht );
   } );
  }

 /// the successors but the turns of the face of op
 template< class F >
 void for_each_successor( const state_type & s , cost_type ,
                          std::uint8_t op , F && f ) const {
  const int last = op == no_op ? -1 : op / 3;
  for( int m = 0 ; m < moves ; ++m )
   if( m / 3 != last )
    f( apply( s , m ) , 1 , 0 , std::uint8_t( m ) );
  }

 /// the successors but the turns of the face of op and, if the face of op
 /// comes second in its pair of opposite faces, those of the first
 template< class F >
 void for_each_tree_successor( const state_type & s , cost_type ,
                               std::uint8_t op , F && f ) const {
  const int last = op == no_op ? -1 : op / 3;
  for( int m = 0 ; m < moves ; ++m ) {
   const int face = m / 3;
   if( ( face == last ) || ( ( last >= 3 ) && ( face == last - 3 ) ) )
    continue;
   f( apply( s , m ) , 1 , 0 , std::uint8_t( m ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] static constexpr std::size_t num_features( void ) noexcept {
  return( 20 * 24 );
  }

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const {
  for( std::size_t i = 0 ; i < 20 ; ++i )
   f( i * 24 + s[ i ] );
  }

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const {
  for( std::size_t i = 0 ; i < 20 ; ++i )
   if( s[ i ] != t[ i ] ) {
    f( i * 24 + s[ i ] );
    f( i * 24 + t[ i ] );
    }
  }

/*------------------------- THE PERMUTATION PUZZLE -------------------------*/
 /// the corners (8, orientations modulo 3) and the edges (12, modulo 2),
 /// whose orientations sum to 0

 [[nodiscard]] static std::vector< PermutationGroup > permutation_groups(
  void ) {
  return( std::vector< PermutationGroup >{
   PermutationGroup{ 8 , 3 , true } , PermutationGroup{ 12 , 2 , true } } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the 18 moves, in the order of their operators

 [[nodiscard]] static std::vector< PermutationMove > permutation_moves(
  void ) {
  std::vector< PermutationMove > v;
  for( const Table & t : tables() ) {
   PermutationMove m;
   m.from = { std::vector< std::uint8_t >( t.from.begin() ,
                                           t.from.begin() + 8 ) ,
              std::vector< std::uint8_t >( t.from.begin() + 8 ,
                                           t.from.end() ) };
   m.turn = { std::vector< std::uint8_t >( t.turn.begin() ,
                                           t.turn.begin() + 8 ) ,
              std::vector< std::uint8_t >( t.turn.begin() + 8 ,
                                           t.turn.end() ) };
   for( auto & x : m.from[ 1 ] )
    x = std::uint8_t( x - 8 );
   v.push_back( std::move( m ) );
   }
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cubies of the group g (0 the corners, 1 the edges) in s, with
 /// their orientations

 static void permutation_of( const state_type & s , unsigned g ,
                             std::uint8_t * obj , std::uint8_t * ori )
  noexcept {
  if( g == 0 )
   for( int i = 0 ; i < 8 ; ++i ) {
    obj[ i ] = std::uint8_t( s[ i ] / 3 );
    ori[ i ] = std::uint8_t( s[ i ] % 3 );
    }
  else
   for( int i = 0 ; i < 12 ; ++i ) {
    obj[ i ] = std::uint8_t( s[ 8 + i ] / 2 );
    ori[ i ] = std::uint8_t( s[ 8 + i ] % 2 );
    }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the goal, the solved cube

 template< class F >
 static void for_each_goal_permutation( unsigned g , F && f ) {
  std::uint8_t obj[ 12 ] , ori[ 12 ];
  for( int i = 0 ; i < 12 ; ++i ) {
   obj[ i ] = std::uint8_t( i );
   ori[ i ] = 0;
   }
  ( void ) g;
  f( obj , ori );
  }

/*-------------------------------- HELPERS ---------------------------------*/
 /// the solved cube

 [[nodiscard]] static constexpr state_type solved( void ) noexcept {
  state_type s{};
  for( int i = 0 ; i < 8 ; ++i )
   s[ i ] = std::uint8_t( 3 * i );
  for( int i = 0 ; i < 12 ; ++i )
   s[ 8 + i ] = std::uint8_t( 2 * i );
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cube s after the move m (operator 3 f + t)

 [[nodiscard]] static state_type apply( const state_type & s , int m )
  noexcept {
  const Table & t = tables()[ std::size_t( m ) ];
  state_type r;
  for( int i = 0 ; i < 8 ; ++i ) {
   const int x = s[ t.from[ i ] ];
   r[ i ] = std::uint8_t( x - x % 3 + ( x % 3 + t.turn[ i ] ) % 3 );
   }
  for( int i = 8 ; i < 20 ; ++i )
   r[ i ] = std::uint8_t( s[ t.from[ i ] ] ^ t.turn[ i ] );
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the operators of the moves in text, in the notation of Singmaster: a
 /// face (U, R, F, D, L, B) followed by nothing (a quarter turn
 /// clockwise), 2 (a half turn), or - or ' (a quarter turn
 /// counterclockwise); throws std::invalid_argument on a bad move

 [[nodiscard]] static std::vector< std::uint8_t > parse_moves(
  const std::string & text ) {
  static const std::string faces = "URFDLB";
  std::vector< std::uint8_t > v;
  for( std::size_t i = 0 ; i < text.size() ; ) {
   const char c = text[ i ];
   if( ( c == ' ' ) || ( c == '\t' ) || ( c == '\r' ) || ( c == '\n' ) ) {
    ++i;
    continue;
    }
   const auto f = faces.find( c );
   if( f == std::string::npos )
    throw( std::invalid_argument( "RubiksCube::parse_moves: bad face" ) );
   int t = 0;
   if( ( i + 1 < text.size() ) && ( text[ i + 1 ] == '2' ) ) {
    t = 1;
    ++i;
    }
   else
    if( ( i + 1 < text.size() ) && ( ( text[ i + 1 ] == '-' ) ||
                                      ( text[ i + 1 ] == '\'' ) ) ) {
     t = 2;
     ++i;
     }
   v.push_back( std::uint8_t( 3 * f + t ) );
   ++i;
   }
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the solved cube after the moves of text [see parse_moves()]

 [[nodiscard]] static state_type scramble( const std::string & text ) {
  state_type s = solved();
  for( std::uint8_t m : parse_moves( text ) )
   s = apply( s , m );
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the solved cube after n random moves, each of a face other than the
 /// last one, the same for a seed on every platform

 [[nodiscard]] static state_type random_walk( int n , std::uint64_t seed ) {
  std::mt19937_64 gen( seed );
  state_type s = solved();
  int last = -1;
  for( int i = 0 ; i < n ; ++i ) {
   int m;
   do
    m = int( random_below( gen , moves ) );
   while( m / 3 == last );
   last = m / 3;
   s = apply( s , m );
   }
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a move on the 20 positions
 struct Table {
  std::array< std::uint8_t , 20 > from;
  std::array< std::uint8_t , 20 > turn;
  };

 /// the 18 moves, from the six quarter turns of Kociemba composed
 [[nodiscard]] static const std::array< Table , moves > & tables( void ) {
  static const std::array< Table , moves > t = [] {
   // the quarter turns U, R, F, D, L, B: the corners, then the edges
   static constexpr std::uint8_t cp[ 6 ][ 8 ] = {
    { 3 , 0 , 1 , 2 , 4 , 5 , 6 , 7 } , { 4 , 1 , 2 , 0 , 7 , 5 , 6 , 3 } ,
    { 1 , 5 , 2 , 3 , 0 , 4 , 6 , 7 } , { 0 , 1 , 2 , 3 , 5 , 6 , 7 , 4 } ,
    { 0 , 2 , 6 , 3 , 4 , 1 , 5 , 7 } , { 0 , 1 , 3 , 7 , 4 , 5 , 2 , 6 } };
   static constexpr std::uint8_t co[ 6 ][ 8 ] = {
    { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } , { 2 , 0 , 0 , 1 , 1 , 0 , 0 , 2 } ,
    { 1 , 2 , 0 , 0 , 2 , 1 , 0 , 0 } , { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
    { 0 , 1 , 2 , 0 , 0 , 2 , 1 , 0 } , { 0 , 0 , 1 , 2 , 0 , 0 , 2 , 1 } };
   static constexpr std::uint8_t ep[ 6 ][ 12 ] = {
    { 3 , 0 , 1 , 2 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 } ,
    { 8 , 1 , 2 , 3 , 11 , 5 , 6 , 7 , 4 , 9 , 10 , 0 } ,
    { 0 , 9 , 2 , 3 , 4 , 8 , 6 , 7 , 1 , 5 , 10 , 11 } ,
    { 0 , 1 , 2 , 3 , 5 , 6 , 7 , 4 , 8 , 9 , 10 , 11 } ,
    { 0 , 1 , 10 , 3 , 4 , 5 , 9 , 7 , 8 , 2 , 6 , 11 } ,
    { 0 , 1 , 2 , 11 , 4 , 5 , 6 , 10 , 8 , 9 , 3 , 7 } };
   static constexpr std::uint8_t eo[ 6 ][ 12 ] = {
    { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
    { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
    { 0 , 1 , 0 , 0 , 0 , 1 , 0 , 0 , 1 , 1 , 0 , 0 } ,
    { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
    { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
    { 0 , 0 , 0 , 1 , 0 , 0 , 0 , 1 , 0 , 0 , 1 , 1 } };
   std::array< Table , moves > r{};
   for( int f = 0 ; f < 6 ; ++f ) {
    Table q;
    for( int i = 0 ; i < 8 ; ++i ) {
     q.from[ i ] = cp[ f ][ i ];
     q.turn[ i ] = co[ f ][ i ];
     }
    for( int i = 0 ; i < 12 ; ++i ) {
     q.from[ 8 + i ] = std::uint8_t( 8 + ep[ f ][ i ] );
     q.turn[ 8 + i ] = eo[ f ][ i ];
     }
    // a then b: from[ i ] = a.from[ b.from[ i ] ], turn[ i ] =
    // a.turn[ b.from[ i ] ] + b.turn[ i ]
    Table x = q;
    for( int t = 0 ; t < 3 ; ++t ) {
     r[ std::size_t( 3 * f + t ) ] = x;
     Table y;
     for( int i = 0 ; i < 20 ; ++i ) {
      y.from[ i ] = x.from[ q.from[ i ] ];
      y.turn[ i ] = std::uint8_t( ( x.turn[ q.from[ i ] ] + q.turn[ i ] ) %
                                  ( i < 8 ? 3 : 2 ) );
      }
     x = y;
     }
    }
   return( r );
   }();
  return( t );
  }

 state_type f_start;                            ///< the start

 };  // end( class( RubiksCube ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS RubikPDBHeuristic ------------------------*/
/*--------------------------------------------------------------------------*/
/// the heuristic of Korf (1997): the maximum of the tables of the 8 corners
/// and of two sets of 6 edges, 4 bits per entry
/** The tables have \f$ 8! \cdot 3^7 = 88179840 \f$ and \f$ 12! / 6! \cdot
 * 2^6 = 42577920 \f$ entries, 82 megabytes in all; the edges are split
 * into UR, UF, UL, UB, FR, FL and DR, DF, DL, DB, BL, BR unless given, two
 * sets that the half turn about the axis of L and R maps into each other,
 * so that both tables have the mean 7.668 of Korf (1997) (the split into
 * the first and the last six edges gives 7.619 and 7.626). */

class RubikPDBHeuristic
 : public PermutationHeuristic< RubiksCube , LexRanking , 4 > {

 public:

 /// the tables of the cube d (which must outlive the heuristic), built on
 /// threads workers (all the cores if 0) or read from the cache directory
 explicit RubikPDBHeuristic( const RubiksCube & d , unsigned threads = 0 ,
                             const std::string & cache = "" ,
                             const std::vector< std::uint8_t > & edges_a = {
                              0 , 1 , 2 , 3 , 8 , 9 } ,
                             const std::vector< std::uint8_t > & edges_b = {
                              4 , 5 , 6 , 7 , 10 , 11 } )
  : PermutationHeuristic< RubiksCube , LexRanking , 4 >( d , tables(
     d , threads , cache , edges_a , edges_b ) ) {}

 private:

 // the three tables, built in place
 [[nodiscard]] static std::vector< PDB > tables(
  const RubiksCube & d , unsigned threads , const std::string & cache ,
  const std::vector< std::uint8_t > & a ,
  const std::vector< std::uint8_t > & b ) {
  std::vector< PDB > v;
  v.reserve( 3 );
  v.emplace_back( d , 0 , std::vector< std::uint8_t >{ 0 , 1 , 2 , 3 , 4 , 5 ,
                                                       6 , 7 } ,
                  threads , cache );
  v.emplace_back( d , 1 , a , threads , cache );
  v.emplace_back( d , 1 , b , threads , cache );
  return( v );
  }

 };  // end( class( RubikPDBHeuristic ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_RubiksCube */

/*--------------------------------------------------------------------------*/
/*------------------------- End File RubiksCube.h --------------------------*/
/*--------------------------------------------------------------------------*/
