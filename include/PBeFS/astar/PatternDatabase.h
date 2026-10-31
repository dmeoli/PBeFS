/*--------------------------------------------------------------------------*/
/*------------------------ File PatternDatabase.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the pattern databases of the sliding-tile puzzles, i.e.,
 * the heuristics that store the exact distances of an abstraction of the
 * puzzle. The class templates are TilePattern (a subset of the tiles and
 * the perfect ranking of their placements), PatternTable (the distances of
 * a pattern, built by a backward breadth-first search and stored by one of
 * the encodings), AdditivePDB (the sum of the tables of a partition of the
 * tiles, with the reflection about the main diagonal) and SlidingTilePDB,
 * the domain SlidingTile with the additive heuristic in place of the
 * Manhattan distance. The theory (the abstraction, the admissibility of the
 * additive sum and of the reflection, the consistency, the encodings) is on
 * the page "Pattern databases" of the documentation.
 *
 * A pattern \f$ \pi \f$ is a set of \f$ k \f$ tiles of a board of
 * \f$ n = W H \f$ cells. A placement of \f$ \pi \f$ is the vector
 * \f$ p = ( p_1 , \ldots , p_k ) \f$ of the distinct cells of its tiles,
 * and its rank is the lexicographic rank of the \f$ k \f$-permutation
 * \f$ p \f$ of \f$ n \f$ items, a bijection onto
 * \f$ [ 0 , n! / ( n - k )! ) \f$. With the blank ignored, a move of the
 * abstraction takes one tile of \f$ \pi \f$ to an adjacent cell that no
 * tile of \f$ \pi \f$ holds, at cost 1; with the blank minimized, the
 * blank is part of the abstract state, the moves of the blank into a cell
 * that no tile of \f$ \pi \f$ holds cost 0 (they move a tile out of the
 * pattern), and the table keeps the least distance over the cells of the
 * blank. The table is filled by a breadth-first search backward from the
 * goal placement, by layers, each layer scanned in parallel by blocks:
 *
 * \code{.unparsed}
 *  1  procedure BUILD-IGNORED(π)               ▷ the blank ignored
 *  2     for each r ∈ [0, N) do h[r] ← ∞          ▷ N = n! / (n − k)!
 *  3     h[rank(goal(π))] ← 0;  d ← 0
 *  4     loop
 *  5        new ← 0
 *  6        for each r with h[r] = d do        ▷ by blocks, in parallel
 *  7           p ← unrank(r)
 *  8           for each tile i of π, each cell c next to p_i, c ∉ p do
 *  9              r' ← rank(p with p_i ← c)
 * 10              if h[r'] = ∞ then          ▷ by compare-and-swap
 * 11                 h[r'] ← d + 1;  new ← new + 1
 * 12        if new = 0 then return h
 * 13        d ← d + 1
 *
 *  1  procedure BUILD-MINIMIZED(π)             ▷ the blank minimized
 *  2     S[r, b] ← unseen for each r and each cell b ∉ unrank(r)
 *  3     S[rank(goal(π)), b] ← current for each b ∉ goal(π)
 *  4     h[r] ← ∞ for each r;  h[rank(goal(π))] ← 0;  d ← 0
 *  5     loop
 *  6        for each (r, b) with S[r, b] = current do    ▷ in parallel
 *  7           p ← unrank(r);  R ← region(p, b)  ▷ the cells b reaches
 *  8           if b ≠ min R then continue      ▷ one cell per region
 *  9           for each c ∈ R, each tile i of π with p_i next to c do
 * 10              p' ← p with p_i ← c;  r' ← rank(p')
 * 11              for each b' ∈ region(p', old p_i) do
 * 12                 if S[r', b'] = unseen then S[r', b'] ← next
 * 13              if h[r'] = ∞ then h[r'] ← d + 1
 * 14        if no state is next then return h
 * 15        S ← S with current → old and next → current;  d ← d + 1
 * \endcode
 *
 * The states \f$ ( r , b ) \f$ of BUILD-MINIMIZED take 2 bits each, as in
 * the two-bit breadth-first search of Korf (AAAI 2008); the cell of the
 * blank is coded by its index among the \f$ n - k \f$ cells free of the
 * pattern. The table of the distances takes 1 byte per placement during
 * the construction, then it is encoded by one of the following classes,
 * all of which store the excess \f$ e( p ) = ( h( p ) - m( p ) ) / 2 \f$
 * over the Manhattan distance \f$ m( p ) \f$ of the tiles of the pattern,
 * which is a nonnegative integer since a move changes \f$ m \f$ by one unit
 * at the cost of one unit, save ByteEntries (the distance itself) and
 * Mod3Entries (the distance modulo 3):
 *
 * - ByteEntries, 8 bits per entry;
 * - NibbleEntries, 4 bits per entry (the excess saturated at 15);
 * - Mod3Entries, \f$ 8 / 5 = 1.6 \f$ bits per entry, five residues in a
 *   byte, after Breyer and Korf (AAAI 2010): the exact value of a
 *   successor follows from that of its parent, since the two differ by at
 *   most one;
 * - MinNibbleEntries, 4 bits per group of \f$ c \f$ consecutive ranks
 *   (\f$ c \f$ a power of 2), the least excess of the group, the lossy
 *   compression of Felner, Korf, Meshulam and Holte (JAIR 2007);
 * - SparseEntries, a bit per entry that tells whether its excess is
 *   nonzero, with rank and select, and the nonzero excesses in 4 bits
 *   each.
 *
 * The lookup of a state \f$ s \f$, whose tile \f$ t \f$ is in the cell
 * \f$ q_t \f$, sums the tables of the partition \f$ \pi_1 , \ldots ,
 * \pi_m \f$ and, on a square board whose goal has the blank on the main
 * diagonal, takes the maximum with the same sum on the reflection
 * \f$ \phi( s ) \f$, which moves the tile \f$ t \f$ of the cell
 * \f$ q_t \f$ to the transposed cell \f$ T( q_t ) \f$ and renames it
 * \f$ \sigma( t ) \f$, the tile that the goal holds in the transposed goal
 * cell of \f$ t \f$:
 *
 * \code{.unparsed}
 *  1  function H(s)
 *  2     σ₁ ← Σ_j VALUE_j(q restricted to π_j)
 *  3     if not reflected then return σ₁
 *  4     for each tile t do q'[σ(t)] ← T(q_t)
 *  5     σ₂ ← Σ_j VALUE_j(q' restricted to π_j)
 *  6     return max(σ₁, σ₂)
 *  7  function VALUE_j(p)                    ▷ r = rank(p), m = m(p)
 *  8     byte: e[r];  nibble, sparse: m + 2 e[r];  min: m + 2 e[r / c]
 *  9  function CHILD-H(h(s), s, t)        ▷ mod 3, tile x of π_j moved
 * 10     δ ← (res_j(t) − res_j(s) + 1) mod 3 − 1
 * 11     return h(s) + δ
 * 12  function EXACT_j(p)       ▷ mod 3: the value from the residues alone
 * 13     v ← 0;  r ← rank(p)
 * 14     while r ≠ rank(goal(π_j)) do
 * 15        r ← a neighbour r' of r with res_j(r') = (res_j(r) − 1) mod 3
 * 16        v ← v + 1
 * 17     return v
 * \endcode
 *
 * A table can be saved to a file and loaded back (a versioned header, the
 * pattern, the goal, the encoding and a checksum of the whole file), and
 * a cache directory keeps the tables built, so that the large ones (those
 * of the 24-puzzle take minutes) are built once.
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

#ifndef __PBeFS_PatternDatabase
 #define __PBeFS_PatternDatabase
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <atomic>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <cstdio>

#include <cstdlib>

#include <cstring>

#include <filesystem>

#include <fstream>

#include <iterator>

#include <stdexcept>

#include <string>

#include <thread>

#include <utility>

#include <vector>

#include "PBeFS/astar/SlidingTile.h"

#include "PBeFS/astar/Structures.h"

#include "PBeFS/core/SearchUtils.h"

#include "ADS/succinct/PackedArray.h"

#include "ADS/succinct/RankSelect.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// how the abstraction of a pattern treats the blank

enum class BlankRule : std::uint8_t {
 ignored ,   ///< the blank is not in the abstract state: consistent tables
 minimized   ///< the blank is, and the table keeps the least distance
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a partition of the tiles into disjoint patterns

using Partition = std::vector< std::vector< int > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the options of the construction of an AdditivePDB

struct PDBOptions {
 BlankRule blank = BlankRule::ignored;     ///< the rule of the blank
 bool reflect = true;          ///< also the reflection, where it is valid
 unsigned compress = 1;        ///< ranks per entry of MinNibbleEntries
 unsigned threads = 0;         ///< workers of the construction (0: all)
 std::string cache;            ///< directory of the tables ("": none)
 };

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS TilePattern ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a pattern of the sliding-tile puzzle and the ranking of its placements
/** The board has W x H cells, at most 32, the cell r W + k being in row r
 * and column k as in SlidingTile. A placement of the k tiles of the
 * pattern is the array of their cells, in the order of the tiles given to
 * the constructor, and its rank is that of the k-permutation in
 * lexicographic order: the digit \f$ d_i \f$ of the tile \f$ i \f$ is the
 * number of the cells below \f$ p_i \f$ not held by the tiles before it,
 * and the rank is \f$ \sum_i d_i \prod_{j > i} ( n - j ) \f$, which takes
 * \f$ k \f$ population counts, the inverse \f$ k \f$ divisions and
 * \f$ k \f$ selections of a bit. Two consecutive ranks differ only in the
 * cell of the last tile, among the cells free of the others. */

template< int W , int H >
class TilePattern {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr int cells = W * H;       ///< the number of cells

 static_assert( ( W >= 2 ) && ( H >= 2 ) && ( cells <= 32 ) ,
                "TilePattern: at most 32 cells" );

 using board_type = std::array< std::uint8_t , cells >;  ///< a board

 /// the cells adjacent to a cell, in the order up, left, right, down
 struct Moves {
  int n = 0;                                ///< their number
  int cell[ 4 ] = { 0 , 0 , 0 , 0 };        ///< the cells
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// the goal of Korf (1985): the blank in the first cell, then 1 , 2 , ...

 [[nodiscard]] static board_type korf_goal( void ) {
  board_type b{};
  for( int i = 0 ; i < cells ; ++i )
   b[ i ] = static_cast< std::uint8_t >( i );
  return( b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, taking the tiles of the pattern and the goal board
 /** Throws std::invalid_argument if \p tiles is empty, has a tile twice or
  * a tile out of 1 , ... , W H - 1, or if \p goal is not a permutation of
  * the cells. */

 explicit TilePattern( const std::vector< int > & tiles ,
                       const board_type & goal = korf_goal() )
  : v_tiles( tiles ) , v_goal( goal ) {
  std::uint64_t seen = 0;
  for( std::uint8_t t : goal ) {
   if( ( t >= cells ) || ( ( seen >> t ) & 1u ) )
    throw( std::invalid_argument( "TilePattern::TilePattern: the goal is "
                                  "not a permutation" ) );
   seen |= std::uint64_t( 1 ) << t;
   }
  if( tiles.empty() || ( int( tiles.size() ) >= cells ) )
   throw( std::invalid_argument( "TilePattern::TilePattern: a pattern has "
                                 "1 to W H - 1 tiles" ) );
  seen = 0;
  for( int t : tiles ) {
   if( ( t <= 0 ) || ( t >= cells ) || ( ( seen >> t ) & 1u ) )
    throw( std::invalid_argument( "TilePattern::TilePattern: bad tile " +
                                  std::to_string( t ) ) );
   seen |= std::uint64_t( 1 ) << t;
   }
  f_k = int( tiles.size() );
  board_type cell{};
  for( int c = 0 ; c < cells ; ++c )
   cell[ goal[ c ] % cells ] = std::uint8_t( c );
  for( std::size_t i = 0 ; i < v_goal_cell.size() ; ++i )
   if( i < tiles.size() )
    v_goal_cell[ i ] = cell[ std::size_t( tiles[ i ] ) % cells ];
  f_entries = 1;
  for( int i = 0 ; i < f_k ; ++i )
   f_entries *= std::uint64_t( cells - i );
  f_goal_rank = rank( v_goal_cell.data() );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number k of the tiles of the pattern

 [[nodiscard]] int size( void ) const noexcept { return( f_k ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the tiles of the pattern

 [[nodiscard]] const std::vector< int > & tiles( void ) const noexcept {
  return( v_tiles );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the goal board

 [[nodiscard]] const board_type & goal( void ) const noexcept {
  return( v_goal );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the goal cell of the i-th tile of the pattern

 [[nodiscard]] int goal_cell( int i ) const noexcept {
  return( v_goal_cell[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of the placements, n! / ( n - k )!

 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_entries );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the rank of the goal placement

 [[nodiscard]] std::uint64_t goal_rank( void ) const noexcept {
  return( f_goal_rank );
  }

/*------------------------- RANKING THE PLACEMENTS -------------------------*/
 /// the rank of the placement p, p[ i ] the cell of the i-th tile

 [[nodiscard]] std::uint64_t rank( const std::uint8_t * p ) const noexcept {
  std::uint64_t r = 0;
  std::uint32_t used = 0;
  for( int i = 0 ; i < f_k ; ++i ) {
   const std::uint32_t below = used & ( ( std::uint32_t( 1 ) << p[ i ] ) -
                                        1u );
   r = r * std::uint64_t( cells - i ) +
       std::uint64_t( p[ i ] - std::popcount( below ) );
   used |= std::uint32_t( 1 ) << p[ i ];
   }
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in p the placement of rank r < entries()

 void unrank( std::uint64_t r , std::uint8_t * p ) const noexcept {
  int d[ cells ];
  for( int i = f_k - 1 ; i >= 0 ; --i ) {
   d[ i ] = int( r % std::uint64_t( cells - i ) );
   r /= std::uint64_t( cells - i );
   }
  std::uint32_t free = c_full;
  for( int i = 0 ; i < f_k ; ++i ) {
   p[ i ] = std::uint8_t( select( free , d[ i ] ) );
   free &= ~( std::uint32_t( 1 ) << p[ i ] );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan distance of the tiles of the placement p from the goal

 [[nodiscard]] int manhattan( const std::uint8_t * p ) const noexcept {
  int m = 0;
  for( int i = 0 ; i < f_k ; ++i )
   m += distance( p[ i ] , v_goal_cell[ i ] );
  return( m );
  }

/*----------------------------- BOARD HELPERS ------------------------------*/
 /// the cells adjacent to the cell c

 [[nodiscard]] static const Moves & moves( int c ) noexcept {
  return( c_moves[ std::size_t( c ) ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan distance between two cells

 [[nodiscard]] static int distance( int a , int b ) noexcept {
  return( std::abs( a / W - b / W ) + std::abs( a % W - b % W ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the mask of all the cells

 static constexpr std::uint32_t c_full = cells == 32 ? ~std::uint32_t( 0 )
  : ( std::uint32_t( 1 ) << cells ) - 1u;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the place of the j-th set bit of x (from 0), which must exist

 [[nodiscard]] static int select( std::uint32_t x , int j ) noexcept {
  for( ; j > 0 ; --j )
   x &= x - 1u;
  return( std::countr_zero( x ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cells that a blank in the cell b reaches through the cells of
 /// free, i.e., the region of b in the mask free (b in free)

 [[nodiscard]] static std::uint32_t region( std::uint32_t free , int b )
  noexcept {
  std::uint32_t r = std::uint32_t( 1 ) << b , last = 0;
  while( r != last ) {
   last = r;
   std::uint32_t g = r | ( r >> W ) | ( r << W ) |
                     ( ( r & ~c_right ) << 1 ) | ( ( r & ~c_left ) >> 1 );
   r = g & free;
   }
  return( r );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the cells of the first column
 static constexpr std::uint32_t c_left = [] {
  std::uint32_t m = 0;
  for( int r = 0 ; r < H ; ++r )
   m |= std::uint32_t( 1 ) << ( r * W );
  return( m );
  }();

 /// the cells of the last column
 static constexpr std::uint32_t c_right = c_left << ( W - 1 );

 /// the moves of each cell
 static constexpr std::array< Moves , cells > c_moves = [] {
  std::array< Moves , cells > m{};
  for( int b = 0 ; b < cells ; ++b ) {
   const int r = b / W , k = b % W;
   Moves & o = m[ std::size_t( b ) ];
   if( r > 0 )
    o.cell[ o.n++ ] = b - W;
   if( k > 0 )
    o.cell[ o.n++ ] = b - 1;
   if( k < W - 1 )
    o.cell[ o.n++ ] = b + 1;
   if( r < H - 1 )
    o.cell[ o.n++ ] = b + W;
   }
  return( m );
  }();

 std::vector< int > v_tiles;                     ///< the tiles
 board_type v_goal;                              ///< the goal board
 std::array< std::uint8_t , cells > v_goal_cell{};  ///< goal cell of tile i
 int f_k = 0;                                    ///< the number of tiles
 std::uint64_t f_entries = 0;                    ///< the placements
 std::uint64_t f_goal_rank = 0;                  ///< rank of the goal

 };  // end( class( TilePattern ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- THE ENCODINGS ------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// the first n bytes of the words w, the least significant byte first,
/// i.e., the bytes of a packed array as they are written in a file
[[nodiscard]] inline std::vector< std::uint8_t > packed_bytes(
 const std::vector< std::uint64_t > & w , std::size_t n )
{
 std::vector< std::uint8_t > b( n );
 for( std::size_t i = 0 ; i < n ; ++i )
  b[ i ] = std::uint8_t( w[ i >> 3 ] >> ( 8 * ( i & 7 ) ) );
 return( b );
 }

/// the field i of B bits of the bytes d, the least significant bits first

template< unsigned B >
[[nodiscard]] inline unsigned packed_field( const std::vector< std::uint8_t >
                                            & d , std::uint64_t i ) noexcept
{
 const std::uint64_t bit = i * B;
 return( ( d[ bit >> 3 ] >> ( bit & 7 ) ) & ( ( 1u << B ) - 1 ) );
 }

}  // end( namespace detail )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table, one byte each: the distance itself
/** Every encoding has the code and the name written in the files, tells
 * whether its value needs that of a neighbour (relative), and is built by
 * encode( h , m , c ) from the distances h, the Manhattan distance m( r )
 * of the placement of rank r and the parameter c, or by decode() from the
 * bytes of payload(). The entries are held by the packed arrays of ADS,
 * whose fields lie in the bytes as they do in the files. */

class ByteEntries {

 public:

 static constexpr std::uint32_t code = 1;           ///< the code in a file
 static constexpr const char * name = "byte";       ///< the name
 static constexpr bool relative = false;            ///< absolute values

 /// the encoding of the distances h

 template< class M >
 [[nodiscard]] static ByteEntries encode( const std::vector< std::uint8_t >
                                          & h , M && , unsigned = 1 ) {
  ADS::Sequential ex;
  ByteEntries e;
  e.f_data = Array( ex , h.size() , [ & ]( std::size_t r ) {
   return( h[ r ] );
   } );
  return( e );
  }

 /// the encoding of n entries read from a file

 [[nodiscard]] static ByteEntries decode( std::vector< std::uint8_t > data ,
                                          std::uint64_t n , unsigned = 1 ) {
  if( data.size() != n )
   throw( std::runtime_error( "ByteEntries::decode: bad size" ) );
  return( encode( data , 0 ) );
  }

 /// the value of the rank r, whose placement has Manhattan distance m

 [[nodiscard]] int value( std::uint64_t r , int ) const noexcept {
  return( int( f_data.get( r ) ) );
  }

 /// the number of the entries
 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_data.size() );
  }

 /// the bytes of the entries
 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_data.size() );
  }

 /// the bytes written in a file
 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  return( detail::packed_bytes( f_data.words() , bytes() ) );
  }

 private:

 using Array = ADS::PackedArray< 8 >;            ///< a byte an entry

 Array f_data;                                   ///< the entries

 };  // end( class( ByteEntries ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table, 4 bits each: the excess over Manhattan
/** The excess e = ( h - m ) / 2 is stored, saturated at 15; a saturated
 * entry gives a smaller value, still admissible, and clamped() counts them
 * at the encoding. */

class NibbleEntries {

 public:

 static constexpr std::uint32_t code = 2;           ///< the code in a file
 static constexpr const char * name = "nibble";     ///< the name
 static constexpr bool relative = false;            ///< absolute values

 /// the encoding of the distances h

 template< class M >
 [[nodiscard]] static NibbleEntries encode( const std::vector< std::uint8_t >
                                            & h , M && m , unsigned = 1 ) {
  ADS::Sequential ex;
  NibbleEntries e;
  e.f_data = Array( ex , h.size() , [ & ]( std::size_t r ) {
   return( std::uint64_t( excess( h[ r ] , m( r ) ) ) );
   } );
  return( e );
  }

 /// the encoding of n entries read from a file

 [[nodiscard]] static NibbleEntries decode( std::vector< std::uint8_t > data ,
                                            std::uint64_t n , unsigned = 1 ) {
  if( data.size() != ( n + 1 ) / 2 )
   throw( std::runtime_error( "NibbleEntries::decode: bad size" ) );
  ADS::Sequential ex;
  NibbleEntries e;
  e.f_data = Array( ex , n , [ & ]( std::size_t r ) {
   return( detail::packed_field< 4 >( data , r ) );
   } );
  return( e );
  }

 /// the value of the rank r, whose placement has Manhattan distance m

 [[nodiscard]] int value( std::uint64_t r , int m ) const noexcept {
  return( m + 2 * int( f_data.get( r ) ) );
  }

 /// the entries saturated at the encoding
 [[nodiscard]] std::uint64_t clamped( void ) const noexcept {
  return( f_data.clamped() );
  }

 /// the number of the entries
 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_data.size() );
  }

 /// the bytes of the entries
 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( ( f_data.size() + 1 ) / 2 );
  }

 /// the bytes written in a file
 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  return( detail::packed_bytes( f_data.words() , bytes() ) );
  }

 /// the excess ( h - m ) / 2, checked to be a nonnegative integer

 [[nodiscard]] static int excess( int h , int m ) {
  if( ( h < m ) || ( ( h - m ) % 2 ) )
   throw( std::logic_error( "NibbleEntries::excess: h - m is not even "
                            "and nonnegative" ) );
  return( ( h - m ) / 2 );
  }

 private:

 using Array = ADS::PackedArray< 4 , true >;     ///< saturated at 15

 Array f_data;                                   ///< two entries a byte

 };  // end( class( NibbleEntries ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table, the distance modulo 3, five in a byte
/** The byte holds \f$ \sum_{j < 5} x_j 3^j \f$ for the residues
 * \f$ x_j \f$ of five consecutive ranks, i.e., \f$ 8 / 5 = 1.6 \f$ bits
 * per entry (Breyer and Korf, AAAI 2010). The value of a placement is not
 * stored: it is that of a neighbour plus the difference of the residues
 * [see AdditivePDB::delta()], or it is found from the residues alone by a
 * descent to the goal [see PatternTable::exact()]; both need a table whose
 * neighbours differ by at most one, i.e., with the blank ignored. */

class Mod3Entries {

 public:

 static constexpr std::uint32_t code = 3;           ///< the code in a file
 static constexpr const char * name = "mod3";       ///< the name
 static constexpr bool relative = true;             ///< relative values
 static constexpr int modulus = 3;                  ///< values modulo 3

 /// the encoding of the distances h

 template< class M >
 [[nodiscard]] static Mod3Entries encode( const std::vector< std::uint8_t >
                                          & h , M && , unsigned = 1 ) {
  ADS::Sequential ex;
  Mod3Entries e;
  e.f_data = ADS::TritArray( ex , h.size() , [ & ]( std::size_t r ) {
   return( unsigned( h[ r ] % 3 ) );
   } );
  return( e );
  }

 /// the encoding of n entries read from a file

 [[nodiscard]] static Mod3Entries decode( std::vector< std::uint8_t > data ,
                                          std::uint64_t n , unsigned = 1 ) {
  if( data.size() != ( n + 4 ) / 5 )
   throw( std::runtime_error( "Mod3Entries::decode: bad size" ) );
  for( std::uint8_t b : data )
   if( b >= 243 )
    throw( std::runtime_error( "Mod3Entries::decode: bad byte" ) );
  static constexpr unsigned power[ 5 ] = { 1 , 3 , 9 , 27 , 81 };
  ADS::Sequential ex;
  Mod3Entries e;
  e.f_data = ADS::TritArray( ex , n , [ & ]( std::size_t r ) {
   return( data[ r / 5 ] / power[ r % 5 ] % 3 );
   } );
  return( e );
  }

 /// the residue modulo 3 of the distance of the rank r, whose placement
 /// has Manhattan distance m

 [[nodiscard]] int residue( std::uint64_t r , int ) const noexcept {
  return( int( f_data.get( r ) ) );
  }

 /// the number of the entries
 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_data.size() );
  }

 /// the bytes of the entries
 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_data.bytes() );
  }

 /// the bytes written in a file
 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  return( f_data.data() );
  }

 private:

 ADS::TritArray f_data;                          ///< five entries a byte

 };  // end( class( Mod3Entries ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table, the second bit of the distance, eight in
/// a byte
/** With the blank ignored, every move of the abstraction costs 1 and
 * changes the Manhattan distance m of the pattern by one unit, so that the
 * distance h has the parity of m and two neighbours differ by exactly
 * one; hence the bit \f$ \lfloor h / 2 \rfloor \bmod 2 \f$ and the parity
 * of m give h modulo 4, which tells \f$ +1 \f$ from \f$ -1 \f$. This is
 * the 1-bit encoding of Clausecker and Reinefeld (SoCS 2019); like
 * Mod3Entries, its values are relative [see Mod3Entries]. */

class ParityBitEntries {

 public:

 static constexpr std::uint32_t code = 6;           ///< the code in a file
 static constexpr const char * name = "bit";        ///< the name
 static constexpr bool relative = true;             ///< relative values
 static constexpr int modulus = 4;                  ///< values modulo 4

 /// the encoding of the distances h

 template< class M >
 [[nodiscard]] static ParityBitEntries encode(
  const std::vector< std::uint8_t > & h , M && m , unsigned = 1 ) {
  ADS::Sequential ex;
  ParityBitEntries e;
  e.f_data = Array( ex , h.size() , [ & ]( std::size_t r ) {
   static_cast< void >( NibbleEntries::excess( h[ r ] , m( r ) ) );
   return( std::uint64_t( ( h[ r ] >> 1 ) & 1 ) );
   } );
  return( e );
  }

 /// the encoding of n entries read from a file

 [[nodiscard]] static ParityBitEntries decode(
  std::vector< std::uint8_t > data , std::uint64_t n , unsigned = 1 ) {
  if( data.size() != ( n + 7 ) / 8 )
   throw( std::runtime_error( "ParityBitEntries::decode: bad size" ) );
  ADS::Sequential ex;
  ParityBitEntries e;
  e.f_data = Array( ex , n , [ & ]( std::size_t r ) {
   return( detail::packed_field< 1 >( data , r ) );
   } );
  return( e );
  }

 /// the residue modulo 4 of the distance of the rank r, whose placement
 /// has Manhattan distance m

 [[nodiscard]] int residue( std::uint64_t r , int m ) const noexcept {
  return( 2 * int( f_data.get( r ) ) + ( m & 1 ) );
  }

 /// the number of the entries
 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_data.size() );
  }

 /// the bytes of the entries
 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( ( f_data.size() + 7 ) / 8 );
  }

 /// the bytes written in a file
 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  return( detail::packed_bytes( f_data.words() , bytes() ) );
  }

 private:

 using Array = ADS::PackedArray< 1 >;            ///< a bit an entry

 Array f_data;                                   ///< eight entries a byte

 };  // end( class( ParityBitEntries ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table, the least excess of each group of ranks
/** The ranks are grouped by c consecutive ones, c a power of 2, and the
 * group stores the least excess over Manhattan of its placements in 4
 * bits, saturated at 15, i.e., \f$ 4 / c \f$ bits per entry: the lossy
 * compression of Felner, Korf, Meshulam and Holte (JAIR 2007), which here
 * keeps the Manhattan distance of each placement. The value is at most the
 * exact one, hence admissible, but a table so compressed is not consistent
 * in general. Consecutive ranks differ in the cell of the last tile of the
 * pattern only, so that a group holds placements that are close. */

class MinNibbleEntries {

 public:

 static constexpr std::uint32_t code = 4;           ///< the code in a file
 static constexpr const char * name = "min";        ///< the name
 static constexpr bool relative = false;            ///< absolute values

 /// the encoding of the distances h, by groups of c ranks

 template< class M >
 [[nodiscard]] static MinNibbleEntries encode(
  const std::vector< std::uint8_t > & h , M && m , unsigned c = 1 ) {
  ADS::Sequential ex;
  MinNibbleEntries e;
  e.set_group( h.size() , c );
  e.f_data = Array( ex , e.groups() , [ & ]( std::size_t g ) {
   const std::uint64_t a = std::uint64_t( g ) << e.f_shift;
   const std::uint64_t b = std::min( e.f_n , a + ( 1u << e.f_shift ) );
   int x = 15;
   for( std::uint64_t r = a ; r < b ; ++r )
    x = std::min( x , NibbleEntries::excess( h[ r ] , m( r ) ) );
   return( std::uint64_t( x ) );
   } );
  return( e );
  }

 /// the encoding of n entries read from a file

 [[nodiscard]] static MinNibbleEntries decode(
  std::vector< std::uint8_t > data , std::uint64_t n , unsigned c = 1 ) {
  MinNibbleEntries e;
  e.set_group( n , c );
  if( data.size() != ( e.groups() + 1 ) / 2 )
   throw( std::runtime_error( "MinNibbleEntries::decode: bad size" ) );
  ADS::Sequential ex;
  e.f_data = Array( ex , e.groups() , [ & ]( std::size_t g ) {
   return( detail::packed_field< 4 >( data , g ) );
   } );
  return( e );
  }

 /// the value of the rank r, whose placement has Manhattan distance m

 [[nodiscard]] int value( std::uint64_t r , int m ) const noexcept {
  return( m + 2 * int( f_data.get( r >> f_shift ) ) );
  }

 /// the number of the entries
 [[nodiscard]] std::uint64_t entries( void ) const noexcept { return( f_n ); }

 /// the ranks per group
 [[nodiscard]] unsigned group( void ) const noexcept {
  return( 1u << f_shift );
  }

 /// the bytes of the entries
 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( ( groups() + 1 ) / 2 );
  }

 /// the bytes written in a file
 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  return( detail::packed_bytes( f_data.words() , bytes() ) );
  }

 private:

 using Array = ADS::PackedArray< 4 >;            ///< 4 bits a group

 // n entries in groups of c, a power of 2
 void set_group( std::uint64_t n , unsigned c ) {
  if( ( c == 0 ) || ! std::has_single_bit( c ) )
   throw( std::invalid_argument( "MinNibbleEntries::set_group: the group "
                                 "is a power of 2" ) );
  f_n = n;
  f_shift = std::countr_zero( c );
  }

 // the number of the groups
 [[nodiscard]] std::uint64_t groups( void ) const noexcept {
  return( ( f_n + ( std::uint64_t( 1 ) << f_shift ) - 1 ) >> f_shift );
  }

 Array f_data;                                   ///< two groups a byte
 std::uint64_t f_n = 0;                          ///< the entries
 int f_shift = 0;                                ///< log2 of the group

 };  // end( class( MinNibbleEntries ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the entries of a pattern table whose excess is nonzero, by rank
/** A vector of bits B tells which entries have a nonzero excess over the
 * Manhattan distance [see NibbleEntries], and these excesses are stored in
 * 4 bits each, in the order of their ranks, saturated at 15. The excess
 * of the rank r is 0 if \f$ B[ r ] = 0 \f$, and otherwise the excess number
 * \f$ rank_1( r ) \f$ of the list, on the RankSelect of ADS. It takes
 * \f$ n ( 1 + \rho ) + 4 z \f$ bits for n entries, z of them nonzero, the
 * directory of rank and select taking \f$ \rho n \f$ bits (about 0.4 n),
 * and hence it is smaller than NibbleEntries when fewer than about 64% of
 * the entries are nonzero, at the price of a rank per lookup. */

class SparseEntries {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 static constexpr std::uint32_t code = 5;           ///< the code in a file
 static constexpr const char * name = "sparse";     ///< the name
 static constexpr bool relative = false;            ///< absolute values

/*------------------------------ CONSTRUCTION ------------------------------*/
 /// the encoding of the distances h, m( r ) the Manhattan distance of r

 template< class M >
 [[nodiscard]] static SparseEntries encode( const std::vector< std::uint8_t >
                                            & h , M && m , unsigned = 1 ) {
  ADS::BitVector b( h.size() );
  Excesses x;
  for( std::uint64_t r = 0 ; r < h.size() ; ++r ) {
   const int mr = m( r );
   if( ( h[ r ] < mr ) || ( ( h[ r ] - mr ) % 2 ) )
    throw( std::logic_error( "SparseEntries::encode: h - m is not even "
                             "and nonnegative" ) );
   const int e = ( h[ r ] - mr ) / 2;
   if( e == 0 )
    continue;
   b.set( r );
   x.push_back( Excesses::Word( e ) );
   }
  SparseEntries s;
  s.f_bits = ADS::RankSelect( std::move( b ) );
  s.f_excess = std::move( x );
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the encoding of n entries read from a file: the words of B, then the
 /// excesses

 [[nodiscard]] static SparseEntries decode( std::vector< std::uint8_t > data ,
                                            std::uint64_t n , unsigned = 1 ) {
  const std::uint64_t words = ( n + 63 ) / 64;
  if( data.size() < 8 * words )
   throw( std::runtime_error( "SparseEntries::decode: bad size" ) );
  ADS::BitVector b( n );
  std::uint64_t z = 0;
  for( std::uint64_t w = 0 ; w < words ; ++w ) {
   std::uint64_t x = 0;
   for( int i = 0 ; i < 8 ; ++i )
    x |= std::uint64_t( data[ 8 * w + std::uint64_t( i ) ] ) << ( 8 * i );
   if( ( w + 1 == words ) && ( n % 64 ) && ( x >> ( n % 64 ) ) )
    throw( std::runtime_error( "SparseEntries::decode: bits past the "
                               "end" ) );
   z += std::uint64_t( std::popcount( x ) );
   for( ; x ; x &= x - 1 )
    b.set( 64 * w + std::uint64_t( std::countr_zero( x ) ) );
   }
  if( data.size() != 8 * words + ( z + 1 ) / 2 )
   throw( std::runtime_error( "SparseEntries::decode: bad size" ) );
  const std::vector< std::uint8_t > tail( data.begin() +
                                         std::ptrdiff_t( 8 * words ) ,
                                         data.end() );
  ADS::Sequential ex;
  SparseEntries s;
  s.f_bits = ADS::RankSelect( std::move( b ) );
  s.f_excess = Excesses( ex , z , [ & ]( std::size_t i ) {
   return( detail::packed_field< 4 >( tail , i ) );
   } );
  return( s );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the rank r, whose placement has Manhattan distance m

 [[nodiscard]] int value( std::uint64_t r , int m ) const {
  if( ! f_bits[ r ] )
   return( m );
  return( m + 2 * int( f_excess.get( f_bits.rank1( r ) ) ) );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the entries

 [[nodiscard]] std::uint64_t entries( void ) const noexcept {
  return( f_bits.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the entries whose excess is nonzero

 [[nodiscard]] std::uint64_t nonzero( void ) const noexcept {
  return( f_bits.ones() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes of B, of its directory and of the excesses

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( ( f_bits.size() + 63 ) / 64 * 8 + f_bits.directoryBits() / 8 +
          ( f_excess.size() + 1 ) / 2 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes written in a file

 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  std::vector< std::uint8_t > d;
  const auto & w = f_bits.bits().words();
  const std::size_t words = ( f_bits.size() + 63 ) / 64;
  const auto x = detail::packed_bytes( f_excess.words() ,
                                      ( f_excess.size() + 1 ) / 2 );
  d.reserve( 8 * words + x.size() );
  for( std::size_t i = 0 ; i < words ; ++i )
   for( int k = 0 ; k < 8 ; ++k )
    d.push_back( std::uint8_t( w[ i ] >> ( 8 * k ) ) );
  d.insert( d.end() , x.begin() , x.end() );
  return( d );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Excesses = ADS::PackedArray< 4 , true >;  ///< saturated at 15

 ADS::RankSelect f_bits;            ///< B, with its rank and select
 Excesses f_excess;                 ///< the nonzero excesses

 };  // end( class( SparseEntries ) )

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS PatternTable --------------------------*/
/*--------------------------------------------------------------------------*/
/// the distances of a pattern of the sliding-tile puzzle, encoded by Enc
/** The table maps the rank of a placement of the tiles of the pattern to
 * its distance from the goal placement in the abstraction given by the
 * rule of the blank [see BlankRule], encoded by Enc (ByteEntries,
 * NibbleEntries, Mod3Entries, ParityBitEntries, MinNibbleEntries or
 * SparseEntries). It is built by build(), read from a file by load(), and
 * taken from a cache directory, or built and saved there, by cached(). */

template< int W , int H , class Enc = NibbleEntries >
class PatternTable {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Pattern = TilePattern< W , H >;           ///< the pattern
 using board_type = typename Pattern::board_type;   ///< a board

 static constexpr int cells = W * H;             ///< the cells

 /// the version of the format of the files
 static constexpr std::uint32_t file_version = 1;

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the pattern, its rule, the parameter of Enc and the
 /// entries already encoded

 PatternTable( const Pattern & p , BlankRule rule , unsigned compress ,
               Enc entries , bool from_file = false )
  : f_pattern( p ) , f_rule( rule ) , f_compress( compress ) ,
    f_entries( std::move( entries ) ) , f_from_file( from_file ) {
  if constexpr( Enc::relative )
   if( rule != BlankRule::ignored )
    throw( std::invalid_argument( "PatternTable::PatternTable: the "
                                  "residues need the blank ignored" ) );
  if( f_entries.entries() != p.entries() )
   throw( std::invalid_argument( "PatternTable::PatternTable: entries and "
                                 "pattern differ" ) );
  for( int i = 0 ; i < p.size() ; ++i )
   for( int c = 0 ; c < cells ; ++c )
    v_md[ std::size_t( i ) ][ std::size_t( c ) ] =
     std::uint8_t( Pattern::distance( c , p.goal_cell( i ) ) );
  }

/*------------------------------ CONSTRUCTION ------------------------------*/
 /// the exact distances of the placements of p, one byte each
 /** Runs the breadth-first search of BUILD-IGNORED or BUILD-MINIMIZED [see
  * the file header] on \p threads workers of Executor (all the cores if 0).
  */

 template< ADS::WorkerExecutor Executor = ADS::ThreadPool >
 [[nodiscard]] static std::vector< std::uint8_t > distances(
  const Pattern & p , BlankRule rule , unsigned threads = 0 ) {
  if( threads == 0 )
   threads = std::max( 1u , std::thread::hardware_concurrency() );
  auto h = rule == BlankRule::ignored
           ? distances_ignored< Executor >( p , threads )
           : distances_minimized< Executor >( p , threads );
  for( std::uint8_t x : h )
   if( x == c_unseen )
    throw( std::logic_error( "PatternTable::distances: a placement is "
                             "not reached" ) );
  return( h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the table of p, from the exact distances h

 [[nodiscard]] static PatternTable encode( const Pattern & p ,
                                           BlankRule rule ,
                                           const std::vector< std::uint8_t >
                                           & h , unsigned compress = 1 ) {
  std::vector< std::uint8_t > q( std::size_t( p.size() ) );
  auto m = [ & ]( std::uint64_t r ) {
   p.unrank( r , q.data() );
   return( p.manhattan( q.data() ) );
   };
  return( PatternTable( p , rule , compress ,
                        Enc::encode( h , m , compress ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the table of p, built on threads workers of Executor

 template< ADS::WorkerExecutor Executor = ADS::ThreadPool >
 [[nodiscard]] static PatternTable build( const Pattern & p ,
                                          BlankRule rule ,
                                          unsigned compress = 1 ,
                                          unsigned threads = 0 ) {
  return( encode( p , rule , distances< Executor >( p , rule , threads ) ,
                  compress ) );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the placement q, q[ i ] the cell of the i-th tile

 [[nodiscard]] int value( const std::uint8_t * q ) const noexcept
  requires( ! Enc::relative ) {
  return( f_entries.value( f_pattern.rank( q ) , manhattan( q ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Manhattan distance of the tiles of the placement q, by a table

 [[nodiscard]] int manhattan( const std::uint8_t * q ) const noexcept {
  int m = 0;
  for( int i = 0 ; i < f_pattern.size() ; ++i )
   m += v_md[ std::size_t( i ) ][ q[ i ] ];
  return( m );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the residue of the value of the placement q, modulo Enc::modulus

 [[nodiscard]] int residue( const std::uint8_t * q ) const noexcept
  requires( Enc::relative ) {
  return( f_entries.residue( f_pattern.rank( q ) , manhattan( q ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the value of the placement q from the residues alone
 /** Walks from \p q to the goal placement, each step to a neighbour whose
  * residue is one less modulo 3, which has a value one less (EXACT of the
  * file header); the number of steps is the value. With the blank
  * ignored, the abstraction is an undirected graph with unit costs, and
  * hence such a neighbour exists until the goal is reached. */

 [[nodiscard]] int exact( const std::uint8_t * q ) const
  requires( Enc::relative ) {
  const int k = f_pattern.size();
  std::uint8_t p[ cells ];
  std::copy( q , q + k , p );
  constexpr int M = Enc::modulus;
  std::uint64_t r = f_pattern.rank( p );
  int v = 0;
  while( r != f_pattern.goal_rank() ) {
   const int want = ( residue( p ) + M - 1 ) % M;
   std::uint32_t occ = 0;
   for( int i = 0 ; i < k ; ++i )
    occ |= std::uint32_t( 1 ) << p[ i ];
   bool moved = false;
   for( int i = 0 ; ( i < k ) && ! moved ; ++i ) {
    const auto & mv = Pattern::moves( p[ i ] );
    for( int j = 0 ; j < mv.n ; ++j ) {
     const int c = mv.cell[ j ];
     if( ( occ >> c ) & 1u )
      continue;
     const std::uint8_t old = p[ i ];
     p[ i ] = std::uint8_t( c );
     if( residue( p ) == want ) {
      r = f_pattern.rank( p );
      moved = true;
      break;
      }
     p[ i ] = old;
     }
    }
   if( ! moved )
    throw( std::logic_error( "PatternTable::exact: no neighbour one step "
                             "closer to the goal" ) );
   ++v;
   }
  return( v );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the pattern

 [[nodiscard]] const Pattern & pattern( void ) const noexcept {
  return( f_pattern );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the rule of the blank

 [[nodiscard]] BlankRule rule( void ) const noexcept { return( f_rule ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the encoded entries

 [[nodiscard]] const Enc & entries( void ) const noexcept {
  return( f_entries );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes of the entries

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_entries.bytes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the table was read from a file

 [[nodiscard]] bool from_file( void ) const noexcept {
  return( f_from_file );
  }

/*------------------------------ PERSISTENCE -------------------------------*/
 /// writes the table to the file path
 /** The file holds, in little-endian order, the magic "PBeFSPDB", the
  * version of the format, W, H, the rule of the blank, the code of the
  * encoding and its parameter, the number k of the tiles and the tiles,
  * the goal board, the number of the entries, the length of the payload
  * and the payload, and at the end a checksum of all that precedes it. It
  * is written to a temporary file in the same directory, then renamed, so
  * that a reader never sees it half written. */

 void save( const std::string & path ) const {
  std::vector< std::uint8_t > b = header( f_pattern , f_rule ,
                                          f_compress );
  const auto pay = f_entries.payload();
  put( b , f_entries.entries() , 8 );
  put( b , pay.size() , 8 );
  b.insert( b.end() , pay.begin() , pay.end() );
  put( b , checksum( b.data() , b.size() ) , 8 );
  const std::string tmp = path + ".tmp" +
                          std::to_string( std::hash< std::thread::id >()(
                                           std::this_thread::get_id() ) );
  {
   std::ofstream out( tmp , std::ios::binary | std::ios::trunc );
   out.write( reinterpret_cast< const char * >( b.data() ) ,
              std::streamsize( b.size() ) );
   if( ! out )
    throw( std::runtime_error( "PatternTable::save: cannot write " + tmp )
           );
   }
  std::filesystem::rename( tmp , path );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads the table of p from the file path
 /** Throws std::runtime_error if the file cannot be read, if its magic,
  * version, checksum or lengths are wrong, or if it holds the table of
  * another pattern, goal, rule or encoding. */

 [[nodiscard]] static PatternTable load( const std::string & path ,
                                         const Pattern & p , BlankRule rule ,
                                         unsigned compress = 1 ) {
  std::ifstream in( path , std::ios::binary );
  if( ! in )
   throw( std::runtime_error( "PatternTable::load: cannot read " + path ) );
  std::vector< std::uint8_t > b( ( std::istreambuf_iterator< char >( in ) ) ,
                                 std::istreambuf_iterator< char >() );
  const std::vector< std::uint8_t > head = header( p , rule , compress );
  const std::size_t fixed = head.size() + 16;
  if( b.size() < fixed + 8 )
   throw( std::runtime_error( "PatternTable::load: " + path +
                              " is truncated" ) );
  if( ! std::equal( head.begin() , head.begin() + 12 , b.begin() ) )
   throw( std::runtime_error( "PatternTable::load: " + path +
                              " is not a table of this version" ) );
  if( get( b.data() + b.size() - 8 , 8 ) !=
      checksum( b.data() , b.size() - 8 ) )
   throw( std::runtime_error( "PatternTable::load: " + path +
                              " is corrupted (checksum)" ) );
  if( ! std::equal( head.begin() , head.end() , b.begin() ) )
   throw( std::runtime_error( "PatternTable::load: " + path +
                              " holds another table" ) );
  const std::uint64_t n = get( b.data() + head.size() , 8 );
  const std::uint64_t len = get( b.data() + head.size() + 8 , 8 );
  if( ( n != p.entries() ) || ( len != b.size() - fixed - 8 ) )
   throw( std::runtime_error( "PatternTable::load: " + path +
                              " has wrong lengths" ) );
  std::vector< std::uint8_t > pay( b.begin() + std::ptrdiff_t( fixed ) ,
                                   b.end() - 8 );
  return( PatternTable( p , rule , compress ,
                        Enc::decode( std::move( pay ) , n , compress ) ,
                        true ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the name of the file of the table of p in a cache directory

 [[nodiscard]] static std::string file_name( const Pattern & p ,
                                             BlankRule rule ,
                                             unsigned compress = 1 ) {
  std::string s = "pdb-" + std::to_string( W ) + "x" + std::to_string( H ) +
                  ( rule == BlankRule::ignored ? "-ignored-" : "-minimized-" )
                  + Enc::name + "-" + std::to_string( compress ) + "-t";
  for( std::size_t i = 0 ; i < p.tiles().size() ; ++i )
   s += ( i ? "." : "" ) + std::to_string( p.tiles()[ i ] );
  const auto & g = p.goal();
  char hex[ 17 ];
  std::snprintf( hex , sizeof( hex ) , "%016llx" , static_cast<
                 unsigned long long >( checksum( g.data() , g.size() ) ) );
  return( s + "-g" + hex + ".pdb" );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the table of p from the cache directory dir, or built and saved there
 /** Looks for the file of the table in \p dir; if it is not there, or it
  * cannot be loaded, looks for the table of the same pattern in bytes
  * (ByteEntries) and encodes it, else builds it by build() on \p threads
  * workers of Executor; the tables built or encoded are saved in \p dir,
  * which is created if needed. */

 template< ADS::WorkerExecutor Executor = ADS::ThreadPool >
 [[nodiscard]] static PatternTable cached( const std::string & dir ,
                                           const Pattern & p ,
                                           BlankRule rule ,
                                           unsigned compress = 1 ,
                                           unsigned threads = 0 ) {
  namespace fs = std::filesystem;
  const fs::path own = fs::path( dir ) / file_name( p , rule , compress );
  if( fs::exists( own ) )
   try {
    return( load( own.string() , p , rule , compress ) );
    }
   catch( std::runtime_error & ) {}            // built again below
  fs::create_directories( dir );
  using Bytes = PatternTable< W , H , ByteEntries >;
  const fs::path raw = fs::path( dir ) / Bytes::file_name( p , rule );
  std::vector< std::uint8_t > h;
  if( fs::exists( raw ) )
   try {
    h = Bytes::load( raw.string() , p , rule ).entries().payload();
    }
   catch( std::runtime_error & ) {}
  if( h.empty() ) {
   h = distances< Executor >( p , rule , threads );
   Bytes( p , rule , 1 , ByteEntries::encode( h , 0 ) ).save(
    raw.string() );
   }
  PatternTable t = encode( p , rule , h , compress );
  if( own != raw )
   t.save( own.string() );
  return( t );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the checksum of the n bytes at b, by words of 8 bytes

 [[nodiscard]] static std::uint64_t checksum( const std::uint8_t * b ,
                                              std::size_t n ) noexcept {
  std::uint64_t h = ADS::splitmix64( 0x5042654653504442ULL ^ n );
  std::size_t i = 0;
  for( ; i + 8 <= n ; i += 8 )
   h = ADS::splitmix64( h ^ get( b + i , 8 ) ) + i;
  if( i < n )
   h = ADS::splitmix64( h ^ get( b + i , int( n - i ) ) );
  return( h );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the value of the entries not reached yet
 static constexpr std::uint8_t c_unseen = 0xFF;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the first bytes of a file: magic, version, sizes, rule, encoding,
 // pattern and goal

 [[nodiscard]] static std::vector< std::uint8_t > header( const Pattern & p ,
                                                          BlankRule rule ,
                                                          unsigned compress
                                                          ) {
  std::vector< std::uint8_t > b = { 'P' , 'B' , 'e' , 'F' , 'S' , 'P' , 'D' ,
                                    'B' };
  put( b , file_version , 4 );
  put( b , std::uint64_t( W ) , 4 );
  put( b , std::uint64_t( H ) , 4 );
  put( b , std::uint64_t( rule ) , 4 );
  put( b , Enc::code , 4 );
  put( b , compress , 4 );
  put( b , std::uint64_t( p.size() ) , 4 );
  for( int t : p.tiles() )
   b.push_back( std::uint8_t( t ) );
  b.insert( b.end() , p.goal().begin() , p.goal().end() );
  return( b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // appends the k low bytes of x, the lowest first

 static void put( std::vector< std::uint8_t > & b , std::uint64_t x ,
                  int k ) {
  for( int i = 0 ; i < k ; ++i )
   b.push_back( std::uint8_t( x >> ( 8 * i ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the integer of the k bytes at b, the lowest first

 [[nodiscard]] static std::uint64_t get( const std::uint8_t * b , int k )
  noexcept {
  std::uint64_t x = 0;
  for( int i = 0 ; i < k ; ++i )
   x |= std::uint64_t( b[ i ] ) << ( 8 * i );
  return( x );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // calls f( a , b ) for the blocks [ a , b ) of [ 0 , n ) on p workers of
 // Executor, each worker taking the next block when it ends one

 template< class Executor , class F >
 static void blocks( std::uint64_t n , unsigned p , F && f ) {
  constexpr std::uint64_t block = 1 << 14;
  std::atomic< std::uint64_t > next{ 0 };
  run_workers< Executor >( nullptr , p , [ & ]( unsigned ) {
   for( ; ; ) {
    const std::uint64_t a = next.fetch_add( block ,
                                            std::memory_order_relaxed );
    if( a >= n )
     break;
    f( a , std::min( n , a + block ) );
    }
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // BUILD-IGNORED: the distances with the blank ignored

 template< class Executor >
 [[nodiscard]] static std::vector< std::uint8_t > distances_ignored(
  const Pattern & pat , unsigned threads ) {
  const int k = pat.size();
  std::vector< std::uint8_t > h( pat.entries() , c_unseen );
  h[ pat.goal_rank() ] = 0;
  for( int d = 0 ; ; ++d ) {
   if( d + 1 >= c_unseen )
    throw( std::logic_error( "PatternTable::distances_ignored: a distance "
                             "beyond 254" ) );
   std::atomic< std::uint64_t > fresh{ 0 };
   blocks< Executor >( h.size() , threads , [ & ]( std::uint64_t a ,
                                                  std::uint64_t b ) {
    std::uint64_t mine = 0;
    std::uint8_t p[ cells ];
    for( std::uint64_t r = a ; r < b ; ++r ) {
     if( std::atomic_ref< std::uint8_t >( h[ r ] ).load(
          std::memory_order_relaxed ) != d )
      continue;
     pat.unrank( r , p );
     std::uint32_t occ = 0;
     for( int i = 0 ; i < k ; ++i )
      occ |= std::uint32_t( 1 ) << p[ i ];
     for( int i = 0 ; i < k ; ++i ) {
      const std::uint8_t old = p[ i ];
      const auto & mv = Pattern::moves( old );
      for( int j = 0 ; j < mv.n ; ++j ) {
       const int c = mv.cell[ j ];
       if( ( occ >> c ) & 1u )
        continue;
       p[ i ] = std::uint8_t( c );
       std::uint8_t e = c_unseen;
       if( std::atomic_ref< std::uint8_t >( h[ pat.rank( p ) ] )
            .compare_exchange_strong( e , std::uint8_t( d + 1 ) ,
                                      std::memory_order_relaxed ) )
        ++mine;
       }
      p[ i ] = old;
      }
     }
    fresh.fetch_add( mine , std::memory_order_relaxed );
    } );
   if( fresh.load() == 0 )
    return( h );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // BUILD-MINIMIZED: the distances with the blank minimized; the state of
 // ( r , b ) takes 2 bits (0 unseen, 1 current, 2 next, 3 old) at the
 // index r ( n - k ) + f, f being the index of b among the free cells

 template< class Executor >
 [[nodiscard]] static std::vector< std::uint8_t > distances_minimized(
  const Pattern & pat , unsigned threads ) {
  const int k = pat.size();
  const std::uint64_t nf = std::uint64_t( cells - k );
  const std::uint64_t states = pat.entries() * nf;
  std::vector< std::uint64_t > s( ( states + 31 ) / 32 , 0 );
  std::vector< std::uint8_t > h( pat.entries() , c_unseen );

  // the index of the free cell b, and the mark of the state i as next
  const auto index = [ nf ]( std::uint64_t r , std::uint32_t occ , int b ) {
   const std::uint32_t below = occ & ( ( std::uint32_t( 1 ) << b ) - 1u );
   return( r * nf + std::uint64_t( b - std::popcount( below ) ) );
   };
  const auto mark = [ & ]( std::uint64_t i ) {
   std::atomic_ref< std::uint64_t > w( s[ i / 32 ] );
   const int sh = int( 2 * ( i % 32 ) );
   std::uint64_t old = w.load( std::memory_order_relaxed );
   for( ; ; ) {
    if( ( old >> sh ) & 3u )
     return( false );
    if( w.compare_exchange_weak( old , old | ( std::uint64_t( 2 ) << sh ) ,
                                 std::memory_order_relaxed ) )
     return( true );
    }
   };

  std::uint8_t g[ cells ];
  pat.unrank( pat.goal_rank() , g );
  std::uint32_t gocc = 0;
  for( int i = 0 ; i < k ; ++i )
   gocc |= std::uint32_t( 1 ) << g[ i ];
  for( int b = 0 ; b < cells ; ++b )
   if( ! ( ( gocc >> b ) & 1u ) ) {
    const std::uint64_t i = index( pat.goal_rank() , gocc , b );
    s[ i / 32 ] |= std::uint64_t( 1 ) << ( 2 * ( i % 32 ) );
    }
  h[ pat.goal_rank() ] = 0;

  constexpr std::uint64_t lo = 0x5555555555555555ULL;
  for( int d = 0 ; ; ++d ) {
   if( d + 1 >= c_unseen )
    throw( std::logic_error( "PatternTable::distances_minimized: a "
                             "distance beyond 254" ) );
   blocks< Executor >( s.size() , threads , [ & ]( std::uint64_t a ,
                                                  std::uint64_t e ) {
    std::uint8_t p[ cells ];
    std::int8_t at[ cells ];
    for( std::uint64_t wi = a ; wi < e ; ++wi ) {
     const std::uint64_t w = std::atomic_ref< std::uint64_t >( s[ wi ] )
                              .load( std::memory_order_relaxed );
     std::uint64_t cur = w & lo & ~( w >> 1 );        // the fields = 1
     while( cur ) {
      const std::uint64_t i = wi * 32 +
                              std::uint64_t( std::countr_zero( cur ) / 2 );
      cur &= cur - 1;
      if( i >= states )
       break;
      const std::uint64_t r = i / nf;
      pat.unrank( r , p );
      std::uint32_t occ = 0;
      std::fill( at , at + cells , std::int8_t( -1 ) );
      for( int t = 0 ; t < k ; ++t ) {
       occ |= std::uint32_t( 1 ) << p[ t ];
       at[ p[ t ] ] = std::int8_t( t );
       }
      const std::uint32_t free = Pattern::c_full & ~occ;
      const int b = Pattern::select( free , int( i % nf ) );
      const std::uint32_t reg = Pattern::region( free , b );
      if( std::countr_zero( reg ) != b )
       continue;                              // not the least of its region
      for( std::uint32_t rr = reg ; rr ; rr &= rr - 1u ) {
       const int c = std::countr_zero( rr );
       const auto & mv = Pattern::moves( c );
       for( int j = 0 ; j < mv.n ; ++j ) {
        const int a2 = mv.cell[ j ];
        const int t = at[ a2 ];
        if( t < 0 )
         continue;
        p[ t ] = std::uint8_t( c );
        const std::uint64_t r2 = pat.rank( p );
        p[ t ] = std::uint8_t( a2 );
        const std::uint32_t occ2 = ( occ & ~( std::uint32_t( 1 ) << a2 ) ) |
                                   ( std::uint32_t( 1 ) << c );
        const std::uint32_t reg2 = Pattern::region( Pattern::c_full & ~occ2 ,
                                                    a2 );
        for( std::uint32_t q = reg2 ; q ; q &= q - 1u )
         mark( index( r2 , occ2 , std::countr_zero( q ) ) );
        std::uint8_t un = c_unseen;
        std::atomic_ref< std::uint8_t >( h[ r2 ] ).compare_exchange_strong(
         un , std::uint8_t( d + 1 ) , std::memory_order_relaxed );
        }
       }
      }
     }
    } );

   // current -> old, next -> current
   std::atomic< std::uint64_t > fresh{ 0 };
   blocks< Executor >( s.size() , threads , [ & ]( std::uint64_t a ,
                                                  std::uint64_t e ) {
    std::uint64_t mine = 0;
    for( std::uint64_t wi = a ; wi < e ; ++wi ) {
     std::uint64_t w = s[ wi ];
     const std::uint64_t cur = w & lo & ~( w >> 1 );
     const std::uint64_t nxt = ( w >> 1 ) & lo & ~w;
     w |= cur << 1;
     w = ( w & ~( nxt << 1 ) ) | nxt;
     s[ wi ] = w;
     mine += std::uint64_t( std::popcount( nxt ) );
     }
    fresh.fetch_add( mine , std::memory_order_relaxed );
    } );
   if( fresh.load() == 0 )
    return( h );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 Pattern f_pattern;                              ///< the pattern
 BlankRule f_rule;                               ///< the rule of the blank
 unsigned f_compress;                            ///< the parameter of Enc
 Enc f_entries;                                  ///< the entries
 bool f_from_file;                               ///< read from a file
 std::array< std::array< std::uint8_t , cells > , cells > v_md{};
                              ///< distance of the i-th tile in c from goal

 };  // end( class( PatternTable ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS AdditivePDB ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the additive heuristic of disjoint pattern tables
/** The tables are those of a partition of (some of) the tiles into disjoint
 * patterns, and the heuristic of a state is the sum of the values of its
 * placements in the tables. With the blank ignored or minimized, a move
 * of the puzzle moves one tile, hence it is a move of the abstraction of
 * at most one pattern, which counts it once; the sum of the abstract
 * distances is thus at most the true one (admissible), and with the
 * blank ignored it changes by at most one along a move (consistent). On a
 * square board whose goal has the blank on the main diagonal, the
 * heuristic is also the maximum with the same sum on the reflection of the
 * state about the diagonal, which is an automorphism of the puzzle fixing
 * the goal [see the file header]. */

template< int W , int H , class Enc = NibbleEntries >
class AdditivePDB {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Table = PatternTable< W , H , Enc >;      ///< a table
 using Pattern = TilePattern< W , H >;           ///< a pattern
 using board_type = typename Pattern::board_type;   ///< a board

 static constexpr int cells = W * H;             ///< the cells

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the tables, all with the same goal

 explicit AdditivePDB( std::vector< Table > tables , bool reflect = true )
  : v_tables( std::move( tables ) ) {
  if( v_tables.empty() )
   throw( std::invalid_argument( "AdditivePDB::AdditivePDB: no table" ) );
  v_goal = v_tables.front().pattern().goal();
  std::fill( v_owner.begin() , v_owner.end() , std::int8_t( -1 ) );
  for( std::size_t j = 0 ; j < v_tables.size() ; ++j ) {
   const Pattern & p = v_tables[ j ].pattern();
   if( p.goal() != v_goal )
    throw( std::invalid_argument( "AdditivePDB::AdditivePDB: the tables "
                                  "have different goals" ) );
   for( int i = 0 ; i < p.size() ; ++i ) {
    const int t = p.tiles()[ std::size_t( i ) ];
    if( v_owner[ std::size_t( t ) ] >= 0 )
     throw( std::invalid_argument( "AdditivePDB::AdditivePDB: the tile " +
                                   std::to_string( t ) +
                                   " is in two patterns" ) );
    v_owner[ std::size_t( t ) ] = std::int8_t( j );
    v_place[ std::size_t( t ) ] = std::int8_t( i );
    }
   }
  int blank = 0;
  for( int c = 0 ; c < cells ; ++c )
   if( v_goal[ c ] == 0 )
    blank = c;
  f_reflect = reflect && ( W == H ) && ( blank / W == blank % W );
  if( f_reflect ) {
   board_type cell{};
   for( int c = 0 ; c < cells ; ++c ) {
    cell[ v_goal[ c ] ] = std::uint8_t( c );
    v_transpose[ c ] = std::uint8_t( ( c % W ) * W + c / W );
    }
   for( int t = 0 ; t < cells ; ++t )
    v_sigma[ t ] = v_goal[ v_transpose[ cell[ t ] ] ];
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, building the tables of a partition for the goal
 /** Builds (or takes from the cache directory of \p o) the table of each
  * pattern of \p partition, each one with all the workers of \p o, on the
  * C++ threads; make() does the same on another runtime. */

 AdditivePDB( const Partition & partition ,
              const board_type & goal = Pattern::korf_goal() ,
              const PDBOptions & o = {} )
  : AdditivePDB( make< ADS::ThreadPool >( partition , goal , o ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the heuristic of a partition, its tables built on Executor

 template< class Executor >
 [[nodiscard]] static AdditivePDB make(
  const Partition & partition ,
  const board_type & goal = Pattern::korf_goal() ,
  const PDBOptions & o = {} ) {
  std::vector< Table > t;
  for( const auto & tiles : partition ) {
   const Pattern p( tiles , goal );
   if( o.cache.empty() )
    t.push_back( Table::template build< Executor >( p , o.blank ,
                                                   o.compress , o.threads ) );
   else
    t.push_back( Table::template cached< Executor >( o.cache , p , o.blank ,
                                                    o.compress ,
                                                    o.threads ) );
   }
  return( AdditivePDB( std::move( t ) , o.reflect ) );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the heuristic of the state whose tile t is in the cell q[ t ]

 [[nodiscard]] int heuristic( const std::uint8_t * q ) const noexcept
  requires( ! Enc::relative ) {
  const int a = sum( q );
  if( ! f_reflect )
   return( a );
  std::uint8_t r[ cells ];
  reflection( q , r );
  return( std::max( a , sum( r ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the sum of the values of the tables for the cells q
 /** The sum of the tables only, without the reflection: for the tests and
  * the benchmarks. */

 [[nodiscard]] int sum( const std::uint8_t * q ) const noexcept
  requires( ! Enc::relative ) {
  int s = 0;
  std::uint8_t p[ cells ];
  for( const Table & t : v_tables ) {
   gather( t , q , p );
   s += t.value( p );
   }
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the heuristic of the cells q, from the residues alone [see
 /// PatternTable::exact()]

 [[nodiscard]] int exact( const std::uint8_t * q ) const
  requires( Enc::relative ) {
  const auto add = [ & ]( const std::uint8_t * x ) {
   int s = 0;
   std::uint8_t p[ cells ];
   for( const Table & t : v_tables ) {
    gather( t , x , p );
    s += t.exact( p );
    }
   return( s );
   };
  const int a = add( q );
  if( ! f_reflect )
   return( a );
  std::uint8_t r[ cells ];
  reflection( q , r );
  return( std::max( a , add( r ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the change of the sum when the tile t moves from the cells q to the
 /// cell c, from the residues of the table of t (CHILD-H of the file
 /// header); 0 for a tile in no pattern

 [[nodiscard]] int delta( const std::uint8_t * q , int t , int c ) const
  noexcept requires( Enc::relative ) {
  const int j = v_owner[ std::size_t( t ) ];
  if( j < 0 )
   return( 0 );
  const Table & tb = v_tables[ std::size_t( j ) ];
  std::uint8_t p[ cells ];
  gather( tb , q , p );
  const int before = tb.residue( p );
  p[ std::size_t( v_place[ std::size_t( t ) ] ) ] = std::uint8_t( c );
  constexpr int M = Enc::modulus;
  return( ( tb.residue( p ) - before + 1 + M ) % M - 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in r the cells of the reflection of the state of cells q

 void reflection( const std::uint8_t * q , std::uint8_t * r ) const noexcept
 {
  for( int t = 0 ; t < cells ; ++t )
   r[ v_sigma[ t ] ] = v_transpose[ q[ t ] ];
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the goal board

 [[nodiscard]] const board_type & goal( void ) const noexcept {
  return( v_goal );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the heuristic also looks up the reflection

 [[nodiscard]] bool reflected( void ) const noexcept { return( f_reflect ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the tables

 [[nodiscard]] const std::vector< Table > & tables( void ) const noexcept {
  return( v_tables );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes of the entries of all the tables

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  std::size_t b = 0;
  for( const Table & t : v_tables )
   b += t.bytes();
  return( b );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // writes in p the placement of the tiles of t, from the cells q
 static void gather( const Table & t , const std::uint8_t * q ,
                     std::uint8_t * p ) noexcept {
  const auto & ts = t.pattern().tiles();
  for( std::size_t i = 0 ; i < ts.size() ; ++i )
   p[ i ] = q[ ts[ i ] ];
  }

 std::vector< Table > v_tables;                  ///< the tables
 board_type v_goal{};                            ///< the goal board
 bool f_reflect = false;                         ///< the reflection used
 std::array< std::int8_t , cells > v_owner{};    ///< table of a tile, or -1
 std::array< std::int8_t , cells > v_place{};    ///< its place there
 std::array< std::uint8_t , cells > v_transpose{};  ///< transposed cell
 std::array< std::uint8_t , cells > v_sigma{};   ///< tile of the reflection

 };  // end( class( AdditivePDB ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SlidingTilePDB --------------------------*/
/*--------------------------------------------------------------------------*/
/// the sliding-tile puzzle with an additive pattern database heuristic
/** The domain SlidingTile, with all its features, hashes and abstractions,
 * whose heuristic is that of an AdditivePDB (which must outlive it) in
 * place of the Manhattan distance; the goal is that of the tables. With
 * Mod3Entries the heuristic of a successor is that of its parent plus the
 * change of the table of the tile moved, which the residues give, and
 * only the initial state is looked up from the residues alone; this needs
 * the sum without the reflection, since from the maximum of two sums the
 * parts cannot be told apart. */

template< int W , int H , class Enc = NibbleEntries >
class SlidingTilePDB : public SlidingTile< W , H > {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Base = SlidingTile< W , H >;                  ///< the puzzle
 using PDB = AdditivePDB< W , H , Enc >;             ///< the heuristic
 using state_type = typename Base::state_type;       ///< a packed board
 using cost_type = typename Base::cost_type;         ///< the cost type
 using board_type = typename Base::board_type;       ///< a board

 static constexpr int cells = Base::cells;           ///< the cells

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the start board and the heuristic

 SlidingTilePDB( const board_type & start , const PDB & pdb )
  : Base( start , pdb.goal() ) , f_pdb( pdb ) {
  if constexpr( Enc::relative )
   if( pdb.reflected() )
    throw( std::invalid_argument( "SlidingTilePDB::SlidingTilePDB: the "
                                  "residues need the sum without the "
                                  "reflection" ) );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/
 /// the heuristic of s, looked up in the tables

 [[nodiscard]] cost_type heuristic( state_type s ) const {
  std::uint8_t q[ cells ];
  cells_of( s , q );
  if constexpr( Enc::relative )
   return( f_pdb.exact( q ) );
  else
   return( f_pdb.heuristic( q ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h( t ) ) for each successor t of s

 template< class F >
 void for_each_successor( state_type s , cost_type h , F && f ) const {
  std::uint8_t q[ cells ];
  cells_of( s , q );
  const int b = Base::blank( s );
  Base::for_each_successor( s , h , [ & ]( state_type t , cost_type c ,
                                           cost_type ) {
   f( t , c , child_h( q , h , t , b ) );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// as above, but without the successor that op leads back to [see
 /// SlidingTile]

 template< class F >
 void for_each_successor( state_type s , cost_type h , std::uint8_t op ,
                          F && f ) const {
  std::uint8_t q[ cells ];
  cells_of( s , q );
  const int b = Base::blank( s );
  Base::for_each_successor( s , h , op , [ & ]( state_type t , cost_type c ,
                                                cost_type , std::uint8_t o ) {
   f( t , c , child_h( q , h , t , b ) , o );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the heuristic

 [[nodiscard]] const PDB & pdb( void ) const noexcept { return( f_pdb ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in q the cell of each tile of s

 static void cells_of( state_type s , std::uint8_t * q ) noexcept {
  for( int c = 0 ; c < cells ; ++c )
   q[ Base::tile( s , c ) ] = std::uint8_t( c );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the heuristic of the successor t of the state of cells q and heuristic
 // h, the blank of which was in b: the tile in b moved there from the
 // blank of t
 cost_type child_h( std::uint8_t * q , cost_type h , state_type t , int b )
  const noexcept {
  const int x = Base::tile( t , b );
  if constexpr( Enc::relative )
   return( h + f_pdb.delta( q , x , b ) );
  else {
   const int from = q[ x ];
   q[ x ] = std::uint8_t( b );
   q[ 0 ] = std::uint8_t( from );
   const cost_type v = f_pdb.heuristic( q );
   q[ x ] = std::uint8_t( from );
   q[ 0 ] = std::uint8_t( b );
   return( v );
   }
  }

 const PDB & f_pdb;                              ///< the heuristic

 };  // end( class( SlidingTilePDB ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the partition 7-8 of the 15-puzzle of Korf and Felner (AIJ 2002)

[[nodiscard]] inline Partition partition_15_7_8( void ) {
 return( Partition{ { 1 , 2 , 3 , 4 , 5 , 6 , 7 } ,
                    { 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 } } );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the partition 5-5-5 of the 15-puzzle of Felner, Korf and Hanan (JAIR
/// 2004)

[[nodiscard]] inline Partition partition_15_5_5_5( void ) {
 return( Partition{ { 1 , 2 , 4 , 5 , 8 } , { 3 , 6 , 7 , 10 , 11 } ,
                    { 9 , 12 , 13 , 14 , 15 } } );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the partition 6-6-3 of the 15-puzzle of Felner, Korf and Hanan (JAIR
/// 2004)

[[nodiscard]] inline Partition partition_15_6_6_3( void ) {
 return( Partition{ { 1 , 2 , 4 , 5 , 8 , 9 } ,
                    { 3 , 6 , 7 , 10 , 11 , 15 } , { 12 , 13 , 14 } } );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the partition 6-6-6-6 of the 24-puzzle of Korf and Felner (AIJ 2002)

[[nodiscard]] inline Partition partition_24_6_6_6_6( void ) {
 return( Partition{ { 1 , 5 , 6 , 10 , 11 , 12 } ,
                    { 2 , 3 , 4 , 7 , 8 , 9 } ,
                    { 13 , 14 , 18 , 19 , 23 , 24 } ,
                    { 15 , 16 , 17 , 20 , 21 , 22 } } );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PatternDatabase */

/*--------------------------------------------------------------------------*/
/*---------------------- End File PatternDatabase.h ------------------------*/
/*--------------------------------------------------------------------------*/
