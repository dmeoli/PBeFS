/*--------------------------------------------------------------------------*/
/*-------------------------- File PatternEntries.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class SparseEntries of the set of data structures of
 * ADS [see Structures.h], an encoding of the entries of a pattern table
 * [see PatternDatabase.h] on the rank and select of ADS; PatternDatabase.h
 * includes it beside its own encodings.
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

#ifndef __PBeFS_ads_PatternEntries
 #define __PBeFS_ads_PatternEntries
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include <utility>

#include <vector>

#include "ADS/succinct/RankSelect.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SparseEntries --------------------------*/
/*--------------------------------------------------------------------------*/
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
  std::vector< std::uint8_t > x;
  std::uint64_t z = 0;
  for( std::uint64_t r = 0 ; r < h.size() ; ++r ) {
   const int mr = m( r );
   if( ( h[ r ] < mr ) || ( ( h[ r ] - mr ) % 2 ) )
    throw( std::logic_error( "SparseEntries::encode: h - m is not even "
                             "and nonnegative" ) );
   const int e = std::min( ( h[ r ] - mr ) / 2 , 15 );
   if( e == 0 )
    continue;
   b.set( r );
   if( z % 2 == 0 )
    x.push_back( std::uint8_t( e ) );
   else
    x.back() = std::uint8_t( x.back() | ( e << 4 ) );
   ++z;
   }
  SparseEntries s;
  s.f_bits = ADS::RankSelect( std::move( b ) );
  s.v_nibbles = std::move( x );
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
  SparseEntries s;
  s.f_bits = ADS::RankSelect( std::move( b ) );
  s.v_nibbles.assign( data.begin() + std::ptrdiff_t( 8 * words ) ,
                      data.end() );
  return( s );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the rank r, whose placement has Manhattan distance m

 [[nodiscard]] int value( std::uint64_t r , int m ) const {
  if( ! f_bits[ r ] )
   return( m );
  const std::size_t i = f_bits.rank1( r );
  return( m + 2 * ( ( v_nibbles[ i >> 1 ] >> ( 4 * ( i & 1 ) ) ) & 0xF ) );
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
          v_nibbles.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes written in a file

 [[nodiscard]] std::vector< std::uint8_t > payload( void ) const {
  std::vector< std::uint8_t > d;
  const auto & w = f_bits.bits().words();
  const std::size_t words = ( f_bits.size() + 63 ) / 64;
  d.reserve( 8 * words + v_nibbles.size() );
  for( std::size_t i = 0 ; i < words ; ++i )
   for( int k = 0 ; k < 8 ; ++k )
    d.push_back( std::uint8_t( w[ i ] >> ( 8 * k ) ) );
  d.insert( d.end() , v_nibbles.begin() , v_nibbles.end() );
  return( d );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 ADS::RankSelect f_bits;            ///< B, with its rank and select
 std::vector< std::uint8_t > v_nibbles;   ///< the nonzero excesses

 };  // end( class( SparseEntries ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_PatternEntries */

/*--------------------------------------------------------------------------*/
/*------------------------ End File PatternEntries.h -----------------------*/
/*--------------------------------------------------------------------------*/
