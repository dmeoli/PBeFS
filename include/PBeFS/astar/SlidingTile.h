/*--------------------------------------------------------------------------*/
/*--------------------------- File SlidingTile.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SlidingTile, which implements the
 * Domain concept [see Domain.h] for the sliding-tile puzzles with the
 * Manhattan distance heuristic. The puzzle with \f$ W H - 1 \f$ tiles is
 * played on a board of \f$ W \times H \f$ cells.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_SlidingTile
 #define __PBeFS_SlidingTile
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <ostream>

#include <sstream>

#include <stdexcept>

#include <string>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SlidingTile ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the sliding-tile puzzle on a board of at most 16 cells
/** The board has \f$ W \times H \f$ cells, with \f$ W H \leq 16 \f$. A
 * state is packed into one 64-bit word, 4 bits per cell, the cell
 * \f$ i = r W + k \f$ holding the tile in row \f$ r \f$ and column \f$ k \f$
 * (0 being the blank): the word is at the same time the state, its equality
 * test and the key of its hash, as advocated by Burns et al. (SoCS 2012).
 * The blank is found as the zero nibble of the word with a few bitwise
 * operations, and the Manhattan distance is updated incrementally from a
 * table of the change caused by each move of each tile.
 *
 * The feature i of Zobrist hashing is "tile t is in cell c", with
 * \f$ i = t W H + c \f$, for the tiles \f$ t = 1 , \ldots , W H - 1 \f$. The
 * abstract feature of Abstract Zobrist hashing keeps of the cell only
 * whether it lies in the upper or in the lower half of the board, which is
 * the hand-crafted projection of Jinnai and Fukunaga (JAIR 2017, fig. 8b)
 * for the 15-puzzle. The state abstraction of AHDA* looks at the cells of
 * tiles 1, 2 and 3, and the nblocks of PBNF at those of the blank and of
 * tiles 1 and 2 [see Burns et al., JAIR 2010]. */

template< int W , int H >
class SlidingTile {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr int cells = W * H;       ///< the number of cells

 static_assert( ( W >= 2 ) && ( H >= 2 ) && ( cells <= 16 ) ,
                "SlidingTile: the packed state holds at most 16 cells" );

 using state_type = std::uint64_t;         ///< a packed board
 using cost_type = int;                    ///< the cost type
 using board_type = std::array< std::uint8_t , cells >;  ///< a board

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// the goal of Korf (1985): the blank in the first cell, then 1 , 2 , ...

 [[nodiscard]] static board_type korf_goal( void ) {
  board_type b{};
  for( int i = 0 ; i < cells ; ++i )
   b[ i ] = static_cast< std::uint8_t >( i );
  return( b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, taking the start and the goal boards
 /** Constructs the puzzle of going from \p start to \p goal; throws
  * std::invalid_argument if either is not a permutation of 0 , ... ,
  * W H - 1, or if the goal cannot be reached from the start. */

 explicit SlidingTile( const board_type & start ,
                       const board_type & goal = korf_goal() )
  : f_start( pack( check( start ) ) ) , f_goal( pack( check( goal ) ) ) {
  if( parity( start ) != parity( goal ) )
   throw( std::invalid_argument( "SlidingTile::SlidingTile: the goal is "
                                 "not reachable from the start" ) );
  for( int c = 0 ; c < cells ; ++c )
   v_goal_cell[ goal[ c ] ] = c;
  for( int t = 1 ; t < cells ; ++t )
   for( int c = 0 ; c < cells ; ++c )
    v_md[ t ][ c ] = distance( c , v_goal_cell[ t ] );
  for( int b = 0 ; b < cells ; ++b ) {
   const int r = b / W , k = b % W;
   Moves & o = v_moves[ b ];
   if( r > 0 )
    o.cell[ o.n++ ] = b - W;
   if( k > 0 )
    o.cell[ o.n++ ] = b - 1;
   if( k < W - 1 )
    o.cell[ o.n++ ] = b + 1;
   if( r < H - 1 )
    o.cell[ o.n++ ] = b + W;
   }
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/
 /// the start board

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if s is the goal board

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_goal );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the hash of a board: its packed word, mixed

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( mix64( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan distance of s from the goal

 [[nodiscard]] cost_type heuristic( state_type s ) const noexcept {
  int h = 0;
  for( int c = 0 ; c < cells ; ++c )
   h += v_md[ tile( s , c ) ][ c ];
  return( h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h( t ) ) for each successor t of s
 /** Calls \p f( t , 1 , h( t ) ) for each successor t of \p s, h( t ) being
  * computed from \p h = h( s ) by the change of the moved tile. */

 template< class F >
 void for_each_successor( state_type s , cost_type h , F && f ) const {
  const int b = blank( s );
  const Moves & o = v_moves[ b ];
  for( int i = 0 ; i < o.n ; ++i ) {
   const int c = o.cell[ i ];
   const int t = tile( s , c );
   const state_type child = ( s & ~( state_type( 0xF ) << ( 4 * c ) ) )
                            | ( state_type( t ) << ( 4 * b ) );
   f( child , 1 , h + v_md[ t ][ b ] - v_md[ t ][ c ] );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// as above, but without the successor that op, the cell where the blank
 /// was, leads back to; the operator of a successor is the cell of the
 /// blank in s

 template< class F >
 void for_each_successor( state_type s , cost_type h , std::uint8_t op ,
                          F && f ) const {
  const int b = blank( s );
  const Moves & o = v_moves[ b ];
  for( int i = 0 ; i < o.n ; ++i ) {
   const int c = o.cell[ i ];
   if( c == op )
    continue;
   const int t = tile( s , c );
   const state_type child = ( s & ~( state_type( 0xF ) << ( 4 * c ) ) )
                            | ( state_type( t ) << ( 4 * b ) );
   f( child , 1 , h + v_md[ t ][ b ] - v_md[ t ][ c ] , std::uint8_t( b ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: one per pair ( tile , cell )

 [[nodiscard]] static constexpr std::size_t num_features( void ) noexcept {
  return( std::size_t( cells ) * cells );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t W H + c ) for each tile t != 0 in each cell c of s

 template< class F >
 void for_each_feature( state_type s , F && f ) const {
  for( int c = 0 ; c < cells ; ++c )
   if( const int t = tile( s , c ) )
    f( std::size_t( t * cells + c ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of abstract features: one per pair ( tile , half board )

 [[nodiscard]] static constexpr std::size_t num_abstract_features( void )
  noexcept {
  return( std::size_t( cells ) * 2 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the abstract feature of a feature: the tile and the half of its cell

 [[nodiscard]] static constexpr std::size_t abstract_feature( std::size_t i )
  noexcept {
  const std::size_t t = i / cells , c = i % cells;
  return( t * 2 + ( c / W >= H / 2 ? 1 : 0 ) );
  }

/*--------------------------- STATE ABSTRACTIONS ---------------------------*/
 /// the abstract state of AHDA*: the cells of tiles 1, 2 and 3

 [[nodiscard]] std::uint64_t abstract_state( state_type s ) const noexcept {
  int p[ 4 ] = { 0 , 0 , 0 , 0 };
  for( int c = 0 ; c < cells ; ++c )
   if( const int t = tile( s , c ) ; ( t >= 1 ) && ( t <= 3 ) )
    p[ t ] = c;
  return( ( std::uint64_t( p[ 1 ] ) * cells + p[ 2 ] ) * cells + p[ 3 ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of nblocks: one per triple of cells, not all of them used

 [[nodiscard]] static constexpr std::size_t num_nblocks( void ) noexcept {
  return( std::size_t( cells ) * cells * cells );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the nblock of PBNF: the cells of the blank and of tiles 1 and 2

 [[nodiscard]] std::size_t nblock( state_type s ) const noexcept {
  int p[ 3 ] = { 0 , 0 , 0 };
  for( int c = 0 ; c < cells ; ++c )
   if( const int t = tile( s , c ) ; t <= 2 )
    p[ t ] = c;
  return( nblock_id( p[ 0 ] , p[ 1 ] , p[ 2 ] ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( b' ) for each abstract successor b' of the nblock b

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const {
  const int p2 = int( b % cells );
  const int p1 = int( b / cells % cells );
  const int p0 = int( b / cells / cells );
  if( ( p0 == p1 ) || ( p0 == p2 ) || ( p1 == p2 ) )
   return;                                    // not an abstract state
  const Moves & o = v_moves[ p0 ];
  for( int i = 0 ; i < o.n ; ++i ) {
   const int c = o.cell[ i ];
   if( c == p1 )
    f( nblock_id( c , p0 , p2 ) );
   else
    if( c == p2 )
     f( nblock_id( c , p1 , p0 ) );
    else
     f( nblock_id( c , p1 , p2 ) );
   }
  }

/*---------------------- OTHER HEURISTICS AND HASHES -----------------------*/
 /// the heuristic of PA*SE: the sum of the displacements of the tiles

 [[nodiscard]] cost_type pairwise_heuristic( state_type s , state_type t )
  const noexcept {
  int ps[ cells ] , pt[ cells ];
  for( int c = 0 ; c < cells ; ++c ) {
   ps[ tile( s , c ) ] = c;
   pt[ tile( t , c ) ] = c;
   }
  int h = 0;
  for( int k = 1 ; k < cells ; ++k )
   h += distance( ps[ k ] , pt[ k ] );
  return( h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the lexicographic rank of the permutation of the cells of the tiles
 /** It is the perfect hash of Korf and Schultze (2005), by which Burns et
  * al. (JAIR 2010) distributed the states of HDA*: as Jinnai and Fukunaga
  * (JAIR 2017) show, modulo the number of processors it depends on very
  * few tiles, and it is kept as a poor work distribution to compare to. */

 [[nodiscard]] std::uint64_t perfect_hash( state_type s ) const noexcept {
  int pos[ cells ];
  for( int c = 0 ; c < cells ; ++c )
   pos[ tile( s , c ) ] = c;
  std::uint64_t rank = 0;
  std::uint32_t used = 0;
  for( int k = 0 ; k < cells ; ++k ) {
   const std::uint32_t below = used & ( ( 1u << pos[ k ] ) - 1u );
   rank = rank * std::uint64_t( cells - k ) +
          std::uint64_t( pos[ k ] - std::popcount( below ) );
   used |= 1u << pos[ k ];
   }
  return( rank );
  }

/*----------------------------- BOARD HELPERS ------------------------------*/
 /// the tile in the cell c of s

 [[nodiscard]] static int tile( state_type s , int c ) noexcept {
  return( int( ( s >> ( 4 * c ) ) & 0xF ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cell of the blank in s, i.e., of the zero nibble of the word

 [[nodiscard]] static int blank( state_type s ) noexcept {
  std::uint64_t t = s | ( s >> 1 );
  t |= t >> 2;
  t = ~t & 0x1111111111111111ULL & c_cell_mask;
  return( std::countr_zero( t ) / 4 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// packs a board into a word

 [[nodiscard]] static state_type pack( const board_type & b ) noexcept {
  state_type s = 0;
  for( int c = 0 ; c < cells ; ++c )
   s |= state_type( b[ c ] ) << ( 4 * c );
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// unpacks a word into a board

 [[nodiscard]] static board_type unpack( state_type s ) noexcept {
  board_type b{};
  for( int c = 0 ; c < cells ; ++c )
   b[ c ] = static_cast< std::uint8_t >( tile( s , c ) );
  return( b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads a board: W H numbers, the tiles cell by cell, 0 for the blank

 [[nodiscard]] static board_type parse( const std::string & text ) {
  std::istringstream in( text );
  board_type b{};
  for( int c = 0 ; c < cells ; ++c ) {
   int t = -1;
   if( ! ( in >> t ) || ( t < 0 ) || ( t >= cells ) )
    throw( std::invalid_argument( "SlidingTile::parse: bad board" ) );
   b[ c ] = static_cast< std::uint8_t >( t );
   }
  return( check( b ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// prints a board, one row per line

 static void print( std::ostream & os , state_type s ) {
  for( int r = 0 ; r < H ; ++r ) {
   for( int k = 0 ; k < W ; ++k ) {
    const int t = tile( s , r * W + k );
    if( t == 0 )
     os << "   ";
    else
     os << ( t < 10 ? "  " : " " ) << t;
    }
   os << '\n';
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the cells adjacent to a cell, in the order up, left, right, down
 struct Moves {
  int n = 0;
  int cell[ 4 ] = { 0 , 0 , 0 , 0 };
  };

 /// the bits of the cells in the word
 static constexpr std::uint64_t c_cell_mask =
  ~std::uint64_t( 0 ) >> ( 64 - 4 * cells );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the Manhattan distance between two cells

 static int distance( int a , int b ) noexcept {
  return( std::abs( a / W - b / W ) + std::abs( a % W - b % W ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the index of the nblock with the blank in p0 and tiles 1 , 2 in p1 , p2

 static std::size_t nblock_id( int p0 , int p1 , int p2 ) noexcept {
  return( ( std::size_t( p0 ) * cells + p1 ) * cells + p2 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // throws unless b is a permutation of 0 , ... , cells - 1

 static const board_type & check( const board_type & b ) {
  std::uint32_t seen = 0;
  for( std::uint8_t t : b ) {
   if( ( t >= cells ) || ( ( seen >> t ) & 1u ) )
    throw( std::invalid_argument( "SlidingTile::check: not a permutation" )
           );
   seen |= 1u << t;
   }
  return( b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the invariant of solvability: the parity of the permutation of the tiles
 // plus, on boards of even width, the row of the blank

 static int parity( const board_type & b ) noexcept {
  int inv = 0 , blank_row = 0;
  for( int i = 0 ; i < cells ; ++i ) {
   if( b[ i ] == 0 ) {
    blank_row = i / W;
    continue;
    }
   for( int j = i + 1 ; j < cells ; ++j )
    if( b[ j ] && ( b[ j ] < b[ i ] ) )
     ++inv;
   }
  return( ( W % 2 == 1 ? inv : inv + blank_row ) % 2 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 state_type f_start;                             ///< the start board
 state_type f_goal;                              ///< the goal board
 int v_goal_cell[ cells ] = {};                  ///< the goal cell of a tile
 int v_md[ cells ][ cells ] = {};    ///< distance of tile t in c from goal
 Moves v_moves[ cells ];                    ///< the cells next to each cell

 };  // end( class( SlidingTile ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

using EightPuzzle = SlidingTile< 3 , 3 >;       ///< the 8-puzzle
using FifteenPuzzle = SlidingTile< 4 , 4 >;     ///< the 15-puzzle

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SlidingTile */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SlidingTile.h -------------------------*/
/*--------------------------------------------------------------------------*/
