/*--------------------------------------------------------------------------*/
/*------------------------ File SequenceAlignment.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SequenceAlignment, which implements
 * the Domain concept [see Domain.h] for the optimal alignment of several
 * protein sequences as a shortest path in an \f$ n \f$-dimensional lattice,
 * with the sum-of-pairs cost of the PAM250 matrix and the heuristic of the
 * optimal pairwise alignments.
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

#ifndef __PBeFS_SequenceAlignment
 #define __PBeFS_SequenceAlignment
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cctype>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <cstring>

#include <istream>

#include <limits>

#include <random>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
//! [msa-instances-sig-symbols]
/// the symbols of a substitution matrix: the 20 amino acids, the
/// ambiguity codes B, J, Z and X, and the gap '*', in this order

inline constexpr char amino_symbols[] = "ARNDCQEGHILKMFPSTWYVBJZX*";

/// the number of symbols, the gap included

inline constexpr int amino_count = 25;

/// the index of the gap among the symbols

inline constexpr int amino_gap = 24;

/// a substitution matrix: the score of each pair of symbols

using ScoreMatrix = std::array< std::array< int , amino_count > ,
                                amino_count >;
//! [msa-instances-sig-symbols]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the PAM250 matrix of Dayhoff, Schwartz and Orcutt (1978)
/** The scores of the mutation data matrix at 250 PAMs of "A model of
 * evolutionary change in proteins" (Atlas of Protein Sequence and
 * Structure 5, supplement 3, 1978), with the ambiguity codes B, J, Z and X
 * and the score -8 of a residue against a gap. */

// clang-format off
inline constexpr ScoreMatrix pam250 = { {
 {  2 , -2 ,  0 ,  0 , -2 ,  0 ,  0 ,  1 , -1 , -1 , -2 , -1 , -1 , -3 ,  1 ,
    1 ,  1 , -6 , -3 ,  0 ,  0 , -1 ,  0 , -1 , -8 } ,
 { -2 ,  6 ,  0 , -1 , -4 ,  1 , -1 , -3 ,  2 , -2 , -3 ,  3 ,  0 , -4 ,  0 ,
    0 , -1 ,  2 , -4 , -2 , -1 , -3 ,  0 , -1 , -8 } ,
 {  0 ,  0 ,  2 ,  2 , -4 ,  1 ,  1 ,  0 ,  2 , -2 , -3 ,  1 , -2 , -3 ,  0 ,
    1 ,  0 , -4 , -2 , -2 ,  2 , -3 ,  1 , -1 , -8 } ,
 {  0 , -1 ,  2 ,  4 , -5 ,  2 ,  3 ,  1 ,  1 , -2 , -4 ,  0 , -3 , -6 , -1 ,
    0 ,  0 , -7 , -4 , -2 ,  3 , -3 ,  3 , -1 , -8 } ,
 { -2 , -4 , -4 , -5 , 12 , -5 , -5 , -3 , -3 , -2 , -6 , -5 , -5 , -4 , -3 ,
    0 , -2 , -8 ,  0 , -2 , -4 , -5 , -5 , -1 , -8 } ,
 {  0 ,  1 ,  1 ,  2 , -5 ,  4 ,  2 , -1 ,  3 , -2 , -2 ,  1 , -1 , -5 ,  0 ,
   -1 , -1 , -5 , -4 , -2 ,  1 , -2 ,  3 , -1 , -8 } ,
 {  0 , -1 ,  1 ,  3 , -5 ,  2 ,  4 ,  0 ,  1 , -2 , -3 ,  0 , -2 , -5 , -1 ,
    0 ,  0 , -7 , -4 , -2 ,  3 , -3 ,  3 , -1 , -8 } ,
 {  1 , -3 ,  0 ,  1 , -3 , -1 ,  0 ,  5 , -2 , -3 , -4 , -2 , -3 , -5 ,  0 ,
    1 ,  0 , -7 , -5 , -1 ,  0 , -4 ,  0 , -1 , -8 } ,
 { -1 ,  2 ,  2 ,  1 , -3 ,  3 ,  1 , -2 ,  6 , -2 , -2 ,  0 , -2 , -2 ,  0 ,
   -1 , -1 , -3 ,  0 , -2 ,  1 , -2 ,  2 , -1 , -8 } ,
 { -1 , -2 , -2 , -2 , -2 , -2 , -2 , -3 , -2 ,  5 ,  2 , -2 ,  2 ,  1 , -2 ,
   -1 ,  0 , -5 , -1 ,  4 , -2 ,  3 , -2 , -1 , -8 } ,
 { -2 , -3 , -3 , -4 , -6 , -2 , -3 , -4 , -2 ,  2 ,  6 , -3 ,  4 ,  2 , -3 ,
   -3 , -2 , -2 , -1 ,  2 , -3 ,  5 , -3 , -1 , -8 } ,
 { -1 ,  3 ,  1 ,  0 , -5 ,  1 ,  0 , -2 ,  0 , -2 , -3 ,  5 ,  0 , -5 , -1 ,
    0 ,  0 , -3 , -4 , -2 ,  1 , -3 ,  0 , -1 , -8 } ,
 { -1 ,  0 , -2 , -3 , -5 , -1 , -2 , -3 , -2 ,  2 ,  4 ,  0 ,  6 ,  0 , -2 ,
   -2 , -1 , -4 , -2 ,  2 , -2 ,  3 , -2 , -1 , -8 } ,
 { -3 , -4 , -3 , -6 , -4 , -5 , -5 , -5 , -2 ,  1 ,  2 , -5 ,  0 ,  9 , -5 ,
   -3 , -3 ,  0 ,  7 , -1 , -4 ,  2 , -5 , -1 , -8 } ,
 {  1 ,  0 ,  0 , -1 , -3 ,  0 , -1 ,  0 ,  0 , -2 , -3 , -1 , -2 , -5 ,  6 ,
    1 ,  0 , -6 , -5 , -1 , -1 , -2 ,  0 , -1 , -8 } ,
 {  1 ,  0 ,  1 ,  0 ,  0 , -1 ,  0 ,  1 , -1 , -1 , -3 ,  0 , -2 , -3 ,  1 ,
    2 ,  1 , -2 , -3 , -1 ,  0 , -2 ,  0 , -1 , -8 } ,
 {  1 , -1 ,  0 ,  0 , -2 , -1 ,  0 ,  0 , -1 ,  0 , -2 ,  0 , -1 , -3 ,  0 ,
    1 ,  3 , -5 , -3 ,  0 ,  0 , -1 , -1 , -1 , -8 } ,
 { -6 ,  2 , -4 , -7 , -8 , -5 , -7 , -7 , -3 , -5 , -2 , -3 , -4 ,  0 , -6 ,
   -2 , -5 , 17 ,  0 , -6 , -5 , -3 , -6 , -1 , -8 } ,
 { -3 , -4 , -2 , -4 ,  0 , -4 , -4 , -5 ,  0 , -1 , -1 , -4 , -2 ,  7 , -5 ,
   -3 , -3 ,  0 , 10 , -2 , -3 , -1 , -4 , -1 , -8 } ,
 {  0 , -2 , -2 , -2 , -2 , -2 , -2 , -1 , -2 ,  4 ,  2 , -2 ,  2 , -1 , -1 ,
   -1 ,  0 , -6 , -2 ,  4 , -2 ,  2 , -2 , -1 , -8 } ,
 {  0 , -1 ,  2 ,  3 , -4 ,  1 ,  3 ,  0 ,  1 , -2 , -3 ,  1 , -2 , -4 , -1 ,
    0 ,  0 , -5 , -3 , -2 ,  3 , -3 ,  2 , -1 , -8 } ,
 { -1 , -3 , -3 , -3 , -5 , -2 , -3 , -4 , -2 ,  3 ,  5 , -3 ,  3 ,  2 , -2 ,
   -2 , -1 , -3 , -1 ,  2 , -3 ,  5 , -2 , -1 , -8 } ,
 {  0 ,  0 ,  1 ,  3 , -5 ,  3 ,  3 ,  0 ,  2 , -2 , -3 ,  0 , -2 , -5 ,  0 ,
    0 , -1 , -6 , -4 , -2 ,  2 , -2 ,  3 , -1 , -8 } ,
 { -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 ,
   -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -8 } ,
 { -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ,
   -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 , -8 ,  0 } } };
// clang-format on

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SequenceAlignment -------------------------*/
/*--------------------------------------------------------------------------*/
/// the optimal alignment of at most K sequences of amino acids
/** With \f$ n \leq K \f$ sequences \f$ s_0 , \ldots , s_{n-1} \f$ of
 * lengths \f$ L_0 , \ldots , L_{n-1} \f$, a state is the vector
 * \f$ ( p_0 , \ldots , p_{n-1} ) \f$ of the residues already aligned, the
 * lattice point of the formulation of the alignment as a shortest path
 * [see Ikeda and Imai, TCS 1999, and Yoshizumi, Miura and Ishida, AAAI
 * 2000]: the initial state is the origin, the goal is
 * \f$ ( L_0 , \ldots , L_{n-1} ) \f$, and a move advances by one the
 * positions of a nonempty set \f$ M \f$ of the sequences not yet ended,
 * i.e., it is a column of the alignment in which the sequence \f$ i \f$
 * has its residue \f$ s_i[ p_i ] \f$ if \f$ i \in M \f$ and a gap
 * otherwise. A column costs the sum over the pairs \f$ i < j \f$ of the
 * cost \f$ c( a_i , a_j ) \f$ of its two symbols (the sum of pairs),
 * \f[
 *  c( a , b ) = S_{\max} - S( a , b ) , \qquad c( - , - ) = 0 ,
 * \f]
 * \f$ S \f$ being the substitution matrix (pam250 by default) and
 * \f$ S_{\max} \f$ its largest score between two residues, so that every
 * cost is nonnegative and maximizing the score of the alignment minimizes
 * its cost; with PAM250, \f$ S_{\max} = 17 \f$ and a residue against a
 * gap costs \f$ 17 + 8 = 25 \f$. The shift by \f$ S_{\max} \f$ changes
 * the gap penalty: if the alignment induces on the pair \f$ i , j \f$
 * \f$ m \f$ columns of two residues and \f$ g \f$ of a residue and a
 * gap (the columns of two gaps dropped), then \f$ 2 m + g = L_i + L_j \f$
 * and the cost of the pair is
 * \f[
 *  \sum ( S_{\max} - S ) + ( S_{\max} - G ) g =
 *  S_{\max} \frac{ L_i + L_j }{ 2 } -
 *  \Bigl( \sum S - \bigl( \tfrac{ S_{\max} }{ 2 } - G \bigr) g \Bigr) ,
 * \f]
 * \f$ G \f$ being the score of a gap (-8 in PAM250) and the first sum
 * running over the \f$ m \f$ columns of two residues: a minimum cost is
 * a maximum score with the linear gap penalty
 * \f$ S_{\max} / 2 - G \f$, 16.5 for PAM250. A positive penalty
 * \f$ P \f$ given to the constructor replaces these costs by
 * \f$ 2 ( S_{\max} - S( a , b ) ) \f$ for two residues and
 * \f$ S_{\max} + 2 P \f$ for a residue and a gap, i.e., twice the costs
 * above with \f$ S_{\max} / 2 + P \f$ against a gap, so that by the
 * same identity a minimum cost is a maximum of \f$ \sum S - P g \f$ over
 * the pairs, the score with the linear gap penalty \f$ P \f$ (8 in the
 * experiments of Jinnai and Fukunaga, JAIR 2017, section 4.1.3), with
 * integer costs.
 *
 * The heuristic is the sum over the pairs of the cost
 * \f$ D_{ij}( p_i , p_j ) \f$ of an optimal alignment of the suffixes of
 * \f$ s_i \f$ and \f$ s_j \f$ from \f$ p_i \f$ and \f$ p_j \f$, computed
 * once for all by dynamic programming in \f$ O( L_i L_j ) \f$: it is
 * admissible because the alignment of all the sequences induces on each
 * pair an alignment no cheaper than \f$ D_{ij} \f$ (two gaps facing each
 * other cost 0 and are dropped), and consistent because each \f$ D_{ij}
 * \f$ is.
 *
 * The features of Zobrist hashing are the pairs \f$ ( i , p_i ) \f$, a
 * move changing those of the sequences that it advances; the abstract
 * features of Abstract Zobrist hashing are the pairs
 * \f$ ( i , \lfloor p_i / b \rfloor ) \f$, \f$ b \f$ being block() (4 by
 * default, as in Jinnai and Fukunaga, JAIR 2017, section 4.1.3), so that
 * the moves within a box of \f$ b^n \f$ lattice points keep the owner of
 * the state, and the abstract state of AHDA* is that box (the abstraction
 * of AHDA* of that paper keeps the position in the longest sequence only,
 * i.e., slices of the lattice instead of boxes). The
 * level of a state is \f$ \sum_i p_i \f$, which a move raises by the
 * number of the sequences it advances: the lattice is a levelled graph,
 * and HyperplaneDistribution gives it the hyperplane work distribution of
 * Kobayashi, Kishimoto and Watanabe (IJCAI 2011), whose planes
 * \f$ \sum_i p_i = c \f$ (or slabs of them) go to the workers in turn,
 * modulo their number. The boxes bound by \f$ 2^n - 1 \f$, the planes by
 * about \f$ n / d \f$, the workers that the successors of the states of a
 * worker go to, \f$ d \f$ being the thickness of the planes; the boxes
 * keep on one worker every move inside a box, the planes every move that
 * stays inside a slab. The nblocks of PBNF are
 * the boxes of side nblock_side(), each next to the boxes that differ by
 * one in some coordinates, all upwards or all downwards: the graph of the
 * boxes is directed, and SafePBNF wants it undirected. */

template< int K = 10 >
class SequenceAlignment {

 static_assert( ( K >= 2 ) && ( K <= 16 ) ,
                "SequenceAlignment: 2 to 16 sequences" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 /// the positions in the sequences, 0 beyond the n-th
 using state_type = std::array< std::uint16_t , K >;
 using cost_type = int;                          ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the sequences and the substitution matrix
 /** The sequences are strings of the symbols of amino_symbols but the
  * gap; \p block is the side of the boxes of AZH and AHDA*, \p
  * nblock_side that of the nblocks of PBNF, and a positive \p penalty the
  * linear gap penalty that replaces the gap scores of the matrix (see
  * above). Throws std::invalid_argument
  * if there are less than 2 or more than K sequences, if one is empty or
  * has 65535 residues or more, or if it has an unknown symbol. */

 explicit SequenceAlignment( const std::vector< std::string > & sequences ,
                             const ScoreMatrix & scores = pam250 ,
                             int block = 4 , int nblock_side = 16 ,
                             int penalty = 0 )
  : f_n( int( sequences.size() ) ) , f_block( block ) ,
    f_nside( nblock_side ) {
  if( ( f_n < 2 ) || ( f_n > K ) )
   throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: 2 "
                                 "to K sequences" ) );
  if( ( block <= 0 ) || ( nblock_side <= 0 ) )
   throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: "
                                 "boxes of side <= 0" ) );
  v_seq.resize( f_n );
  for( int i = 0 ; i < f_n ; ++i ) {
   const auto & s = sequences[ i ];
   if( s.empty() || ( s.size() >= 65535 ) )
    throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: "
                                  "empty or too long sequence" ) );
   for( char a : s ) {
    const int k = symbol( a );
    if( ( k < 0 ) || ( k == amino_gap ) )
     throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: "
                                   "unknown symbol" ) );
    v_seq[ i ].push_back( std::uint8_t( k ) );
    }
   f_goal[ i ] = std::uint16_t( s.size() );
   }

  int smax = scores[ 0 ][ 0 ];
  for( int a = 0 ; a < amino_gap ; ++a )
   for( int b = 0 ; b < amino_gap ; ++b )
    smax = std::max( smax , scores[ a ][ b ] );
  if( penalty < 0 )
   throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: "
                                 "negative gap penalty" ) );
  for( int a = 0 ; a < amino_count ; ++a )
   for( int b = 0 ; b < amino_count ; ++b )
    if( penalty == 0 )
     f_cost[ a ][ b ] = smax - scores[ a ][ b ];
    else
     f_cost[ a ][ b ] = ( a == amino_gap ) || ( b == amino_gap )
                        ? smax + 2 * penalty
                        : 2 * ( smax - scores[ a ][ b ] );
  f_cost[ amino_gap ][ amino_gap ] = 0;
  for( int a = 0 ; a < amino_count ; ++a )
   for( int b = 0 ; b < amino_count ; ++b )
    if( f_cost[ a ][ b ] < 0 )
     throw( std::invalid_argument( "SequenceAlignment::SequenceAlignment: "
                                   "a gap scores more than a residue" ) );

  std::size_t off = 0 , aoff = 0;
  std::uint64_t nb = 1;
  for( int i = 0 ; i < f_n ; ++i ) {
   v_off[ i ] = off;
   v_aoff[ i ] = aoff;
   off += std::size_t( f_goal[ i ] ) + 1;
   aoff += std::size_t( f_goal[ i ] / f_block ) + 1;
   v_nbs[ i ] = f_goal[ i ] / f_nside + 1;
   nb *= v_nbs[ i ];
   }
  f_features = off;
  f_afeatures = aoff;
  f_nblocks = std::size_t( nb );
  for( std::size_t k = 0 ; k < off ; ++k ) {
   int i = f_n - 1;
   while( v_off[ i ] > k )
    --i;
   v_aof.push_back( v_aoff[ i ] + ( k - v_off[ i ] ) / f_block );
   }
  pairwise_tables();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads the sequences of a file in the FASTA format
 /** Each sequence follows a line starting with '>', its name, and may span
  * several lines; the symbols are read as upper case, and the spaces and
  * the characters '-' and '.' of an aligned file are skipped. This is the
  * format of the files .tfa of BAliBASE 3.0 of Thompson et al. (Proteins
  * 2005). Throws std::invalid_argument if no sequence is read. */

 [[nodiscard]] static std::vector< std::string > read_fasta(
                                                       std::istream & in ) {
  std::vector< std::string > v;
  std::string line;
  while( std::getline( in , line ) ) {
   if( ! line.empty() && ( line[ 0 ] == '>' ) ) {
    v.emplace_back();
    continue;
    }
   if( v.empty() )
    continue;
   for( char c : line )
    if( ! std::isspace( static_cast< unsigned char >( c ) ) &&
        ( c != '-' ) && ( c != '.' ) )
     v.back().push_back( char( std::toupper(
                                static_cast< unsigned char >( c ) ) ) );
   }
  if( v.empty() )
   throw( std::invalid_argument( "SequenceAlignment::read_fasta: no "
                                 "sequence" ) );
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads a substitution matrix: the 625 scores, row by row, in the order
 /// of amino_symbols

 [[nodiscard]] static ScoreMatrix read_scores( std::istream & in ) {
  ScoreMatrix m{};
  for( auto & row : m )
   for( auto & x : row )
    if( ! ( in >> x ) )
     throw( std::invalid_argument( "SequenceAlignment::read_scores: bad "
                                   "matrix" ) );
  return( m );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// random sequences of the 20 amino acids
 /** The random instances of Jinnai and Fukunaga (JAIR 2017, section
  * 4.1.3): \p n sequences, n drawn uniformly in \f$ \{ 4 , \ldots , 9
  * \} \f$ if it is 0, each of a length drawn uniformly among the integers
  * \f$ l \f$ with \f$ 0.9 T / n < l < 1.1 T / n \f$, \f$ T \f$ being \p
  * total (5000 there), and each residue drawn uniformly among the first 20
  * symbols of amino_symbols, all from the seed. Throws
  * std::invalid_argument if no length is in the range or n exceeds K. */

 [[nodiscard]] static std::vector< std::string > random_sequences(
                      std::uint64_t seed , int n = 0 , int total = 5000 ) {
  std::mt19937_64 gen( seed );
  if( n == 0 )
   n = 4 + int( random_below( gen , 6 ) );
  const double mean = double( total ) / std::max( n , 1 );
  const int lo = int( std::floor( 0.9 * mean ) ) + 1;
  const int hi = int( std::ceil( 1.1 * mean ) ) - 1;
  if( ( n < 2 ) || ( n > K ) || ( lo < 1 ) || ( hi < lo ) )
   throw( std::invalid_argument( "SequenceAlignment::random_sequences: no "
                                 "such sequences" ) );
  std::vector< std::string > v( n );
  for( auto & s : v ) {
   s.resize( std::size_t( lo ) + std::size_t( random_below(
                                  gen , std::uint64_t( hi - lo + 1 ) ) ) );
   for( char & c : s )
    c = amino_symbols[ random_below( gen , 20 ) ];
   }
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the index of the symbol a in amino_symbols, -1 if it is not there

 [[nodiscard]] static int symbol( char a ) noexcept {
  const char * p = std::strchr( amino_symbols , a );
  return( ( a == '\0' ) || ! p ? -1 : int( p - amino_symbols ) );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/
 /// the origin of the lattice

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( state_type{} );
  }

 /// true if every sequence is ended

 [[nodiscard]] bool is_goal( const state_type & s ) const noexcept {
  return( s == f_goal );
  }

 /// the hash of a state: its positions, mixed

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const noexcept {
  std::uint64_t h = 0;
  for( int i = 0 ; i < K ; i += 4 ) {
   std::uint64_t w = 0;
   for( int j = i ; j < std::min( i + 4 , K ) ; ++j )
    w |= std::uint64_t( s[ j ] ) << ( 16 * ( j - i ) );
   h = ADS::splitmix64( h ^ w );
   }
  return( h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the sum of the costs of the optimal pairwise alignments of the
 /// suffixes

 [[nodiscard]] cost_type heuristic( const state_type & s ) const noexcept {
  cost_type h = 0;
  std::size_t k = 0;
  for( int i = 0 ; i < f_n ; ++i )
   for( int j = i + 1 ; j < f_n ; ++j , ++k )
    h += v_pair[ k ][ std::size_t( s[ i ] ) * ( f_goal[ j ] + 1u ) +
                      s[ j ] ];
  return( h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) ) for each column that can follow s, in
 /// decreasing order of the set of the sequences it advances

 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  std::uint32_t open = 0;
  int a[ K ];
  for( int i = 0 ; i < f_n ; ++i )
   if( s[ i ] < f_goal[ i ] ) {
    open |= std::uint32_t( 1 ) << i;
    a[ i ] = v_seq[ i ][ s[ i ] ];
    }
  for( std::uint32_t m = open ; m ; m = ( m - 1 ) & open ) {
   state_type t = s;
   int x[ K ];
   for( int i = 0 ; i < f_n ; ++i )
    if( m & ( std::uint32_t( 1 ) << i ) ) {
     ++t[ i ];
     x[ i ] = a[ i ];
     }
    else
     x[ i ] = amino_gap;
   cost_type c = 0;
   for( int i = 0 ; i < f_n ; ++i )
    for( int j = i + 1 ; j < f_n ; ++j )
     c += f_cost[ x[ i ] ][ x[ j ] ];
   f( t , c , heuristic( t ) );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the positions 0 to \f$ L_i \f$ of each
 /// sequence

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( f_features );
  }

 /// calls f for the feature \f$ ( i , p_i ) \f$ of each sequence

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const {
  for( int i = 0 ; i < f_n ; ++i )
   f( v_off[ i ] + s[ i ] );
  }

 /// calls f for the features of s and of t on the sequences that the move
 /// from s to t advances

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const {
  for( int i = 0 ; i < f_n ; ++i )
   if( s[ i ] != t[ i ] ) {
    f( v_off[ i ] + s[ i ] );
    f( v_off[ i ] + t[ i ] );
    }
  }

 /// the number of abstract features: the boxes of b positions of each
 /// sequence

 [[nodiscard]] std::size_t num_abstract_features( void ) const noexcept {
  return( f_afeatures );
  }

 /// the abstract feature \f$ ( i , \lfloor p / b \rfloor ) \f$ of the
 /// feature \f$ ( i , p ) \f$

 [[nodiscard]] std::size_t abstract_feature( std::size_t k ) const noexcept
 {
  return( v_aof[ k ] );
  }

/*--------------------------- STATE ABSTRACTIONS ---------------------------*/
 /// the level of s, the number of the residues aligned

 [[nodiscard]] std::uint64_t level( const state_type & s ) const noexcept {
  std::uint64_t l = 0;
  for( int i = 0 ; i < f_n ; ++i )
   l += s[ i ];
  return( l );
  }

 /// the level of the goal, the total length of the sequences

 [[nodiscard]] std::uint64_t max_level( void ) const noexcept {
  return( level( f_goal ) );
  }

 /// the abstract state of AHDA*: the box of side b of s

 [[nodiscard]] std::uint64_t abstract_state( const state_type & s ) const
  noexcept {
  std::uint64_t a = 0;
  for( int i = 0 ; i < f_n ; ++i )
   a = a * ( std::uint64_t( f_goal[ i ] / f_block ) + 1 ) +
       s[ i ] / f_block;
  return( a );
  }

 /// the number of nblocks

 [[nodiscard]] std::size_t num_nblocks( void ) const noexcept {
  return( f_nblocks );
  }

 /// the nblock of PBNF: the box of side nblock_side() of s

 [[nodiscard]] std::size_t nblock( const state_type & s ) const noexcept {
  std::size_t b = 0;
  for( int i = 0 ; i < f_n ; ++i )
   b = b * v_nbs[ i ] + s[ i ] / f_nside;
  return( b );
  }

 /// calls f( b' ) for each box b' that differs from b by one in a
 /// nonempty set of coordinates, all upwards or all downwards

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const {
  std::size_t c[ K ] , w[ K ];
  std::size_t r = b , weight = 1;
  for( int i = f_n - 1 ; i >= 0 ; --i ) {
   c[ i ] = r % v_nbs[ i ];
   r /= v_nbs[ i ];
   w[ i ] = weight;
   weight *= v_nbs[ i ];
   }
  const std::uint32_t all = ( std::uint32_t( 1 ) << f_n ) - 1;
  for( int dir : { 1 , -1 } ) {
   std::uint32_t ok = 0;
   for( int i = 0 ; i < f_n ; ++i )
    if( dir > 0 ? c[ i ] + 1 < v_nbs[ i ] : c[ i ] > 0 )
     ok |= std::uint32_t( 1 ) << i;
   for( std::uint32_t m = ok & all ; m ; m = ( m - 1 ) & ok ) {
    std::size_t t = b;
    for( int i = 0 ; i < f_n ; ++i )
     if( m & ( std::uint32_t( 1 ) << i ) )
      t = dir > 0 ? t + w[ i ] : t - w[ i ];
    f( t );
    }
   }
  }

 /// the side of the nblocks

 [[nodiscard]] int nblock_side( void ) const noexcept { return( f_nside ); }

 /// the side of the boxes of AZH and AHDA*

 [[nodiscard]] int block( void ) const noexcept { return( f_block ); }

/*----------------------------- OTHER METHODS ------------------------------*/
 /// the number of sequences

 [[nodiscard]] int size( void ) const noexcept { return( f_n ); }

 /// the length of the sequence i

 [[nodiscard]] int length( int i ) const noexcept { return( f_goal[ i ] ); }

 /// the cost of the pair of symbols a and b, by their indices

 [[nodiscard]] cost_type pair_cost( int a , int b ) const noexcept {
  return( f_cost[ a ][ b ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the rows of the alignment of a path of states, '-' for a gap

 [[nodiscard]] std::vector< std::string > alignment(
                            const std::vector< state_type > & path ) const {
  std::vector< std::string > rows( f_n );
  for( std::size_t k = 1 ; k < path.size() ; ++k )
   for( int i = 0 ; i < f_n ; ++i )
    rows[ i ].push_back( path[ k ][ i ] != path[ k - 1 ][ i ]
                         ? amino_symbols[ v_seq[ i ][ path[ k - 1 ][ i ] ] ]
                         : '-' );
  return( rows );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the tables D_ij of the optimal alignments of the suffixes of the pairs,
 // filled backwards from ( L_i , L_j )

 void pairwise_tables( void ) {
  for( int i = 0 ; i < f_n ; ++i )
   for( int j = i + 1 ; j < f_n ; ++j ) {
    const std::size_t li = f_goal[ i ] , lj = f_goal[ j ] , w = lj + 1;
    std::vector< cost_type > d( ( li + 1 ) * w );
    for( std::size_t p = li + 1 ; p-- > 0 ; )
     for( std::size_t q = lj + 1 ; q-- > 0 ; ) {
      cost_type best = p == li && q == lj
                      ? 0 : std::numeric_limits< cost_type >::max();
      if( p < li )
       best = std::min( best , d[ ( p + 1 ) * w + q ] +
                               f_cost[ v_seq[ i ][ p ] ][ amino_gap ] );
      if( q < lj )
       best = std::min( best , d[ p * w + q + 1 ] +
                               f_cost[ amino_gap ][ v_seq[ j ][ q ] ] );
      if( ( p < li ) && ( q < lj ) )
       best = std::min( best , d[ ( p + 1 ) * w + q + 1 ] +
                         f_cost[ v_seq[ i ][ p ] ][ v_seq[ j ][ q ] ] );
      d[ p * w + q ] = best;
      }
    v_pair.push_back( std::move( d ) );
    }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 int f_n;                                       ///< the sequences
 int f_block;                          ///< the side of the AZH boxes
 int f_nside;                          ///< the side of the nblocks
 std::vector< std::vector< std::uint8_t > > v_seq;  ///< the symbols
 state_type f_goal{};                           ///< the lengths
 std::array< std::array< cost_type , amino_count > , amino_count > f_cost{};
                                                ///< the costs of the pairs
 std::vector< std::vector< cost_type > > v_pair;  ///< the tables D_ij
 std::array< std::size_t , K > v_off{};   ///< the first feature of each
 std::array< std::size_t , K > v_aoff{};  ///< the first abstract one
 std::array< std::size_t , K > v_nbs{};   ///< the nblocks along each
 std::vector< std::size_t > v_aof;        ///< the abstract feature of each
 std::size_t f_features = 0;              ///< the features
 std::size_t f_afeatures = 0;             ///< the abstract features
 std::size_t f_nblocks = 0;               ///< the nblocks

 };  // end( class( SequenceAlignment ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SequenceAlignment */

/*--------------------------------------------------------------------------*/
/*---------------------- End File SequenceAlignment.h ----------------------*/
/*--------------------------------------------------------------------------*/
