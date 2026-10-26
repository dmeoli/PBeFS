/*--------------------------------------------------------------------------*/
/*---------------------- File MissionariesCannibals.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class MissionariesCannibals, which implements the
 * Domain concept [see Domain.h] for the missionaries and cannibals problem
 * of exercise 3.9 of Russell and Norvig, "Artificial Intelligence: A
 * Modern Approach" (3rd edition, 2010), with m missionaries, c cannibals
 * and a boat of b seats.
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

#ifndef __PBeFS_MissionariesCannibals
 #define __PBeFS_MissionariesCannibals
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS MissionariesCannibals -----------------------*/
/*--------------------------------------------------------------------------*/
/// the missionaries and cannibals problem
/** Everybody is on the left bank at the start, with the boat; a crossing
 * carries from 1 to b people, and on neither bank may the missionaries be
 * outnumbered by the cannibals, if any missionary is there. Every crossing
 * costs 1, and the heuristic is the least number of crossings that would
 * carry the people left over the river with the boat always full, which is
 * admissible. A state packs the missionaries on the left bank in the bits
 * 0 , ... , 15, the cannibals there in the bits 16 , ... , 31, and the
 * bank of the boat in the bit 32 (0 = left); the word is also its perfect
 * hash.
 *
 * The features of Zobrist hashing are the number of the missionaries on
 * the left bank (m + 1 features), that of the cannibals there (c + 1) and
 * the bank of the boat (2). */

class MissionariesCannibals {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint64_t;               ///< the packed state
 using cost_type = int;                          ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of the problem with m missionaries, c cannibals and a boat
 /// of b seats, m and c at most 65535

 MissionariesCannibals( int m , int c , int b )
  : f_m( m ) , f_c( c ) , f_b( b ) {
  if( ( m < 0 ) || ( c < 0 ) || ( m > 65535 ) || ( c > 65535 ) ||
      ( b < 1 ) )
   throw( std::invalid_argument( "MissionariesCannibals::"
                                 "MissionariesCannibals: bad sizes" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the missionaries on the left bank in the state s

 [[nodiscard]] static int missionaries( std::uint64_t s ) noexcept {
  return( int( s & 0xFFFF ) );
  }

 /// the cannibals on the left bank in the state s

 [[nodiscard]] static int cannibals( std::uint64_t s ) noexcept {
  return( int( ( s >> 16 ) & 0xFFFF ) );
  }

 /// true if the boat is on the left bank in the state s

 [[nodiscard]] static bool boat_left( std::uint64_t s ) noexcept {
  return( ( ( s >> 32 ) & 1 ) == 0 );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] std::uint64_t initial_state( void ) const noexcept {
  return( pack( f_m , f_c , true ) );
  }

 [[nodiscard]] bool is_goal( std::uint64_t s ) const noexcept {
  return( ( missionaries( s ) == 0 ) && ( cannibals( s ) == 0 ) );
  }

 [[nodiscard]] std::uint64_t hash( std::uint64_t s ) const noexcept {
  return( mix64( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// with n people on the left bank, a full boat carries b of them over and
 /// brings 1 back: at least 2 ceil( ( n - b ) / ( b - 1 ) ) + 1 crossings
 /// if n > b, 1 if 0 < n <= b, one more if the boat is on the right bank

 [[nodiscard]] int heuristic( std::uint64_t s ) const noexcept {
  const int n = missionaries( s ) + cannibals( s );
  if( n == 0 )
   return( 0 );
  int h = 1;
  if( n > f_b )
   h = f_b > 1 ? 2 * ( ( n - f_b + ( f_b - 1 ) - 1 ) / ( f_b - 1 ) ) + 1
               : 2 * n - 1;
  return( boat_left( s ) ? h : h + 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h( t ) ) for each safe crossing of dm missionaries
 /// and dc cannibals, 1 <= dm + dc <= b, in increasing dm and then dc

 template< class F >
 void for_each_successor( std::uint64_t s , int , F && f ) const {
  const bool left = boat_left( s );
  const int m = missionaries( s ) , c = cannibals( s );
  const int am = left ? m : f_m - m , ac = left ? c : f_c - c;  // aboard
  for( int dm = 0 ; dm <= std::min( am , f_b ) ; ++dm )
   for( int dc = dm ? 0 : 1 ; dc <= std::min( ac , f_b - dm ) ; ++dc ) {
    const int nm = left ? m - dm : m + dm , nc = left ? c - dc : c + dc;
    if( safe( nm , nc ) && safe( f_m - nm , f_c - nc ) ) {
     const std::uint64_t t = pack( nm , nc , ! left );
     f( t , 1 , heuristic( t ) );
     }
    }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the counts on the left bank and the boat

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( std::size_t( f_m ) + f_c + 4 );
  }

 /// calls f( m ), f( M + 1 + c ) and f( M + C + 2 + boat ) for the m
 /// missionaries and c cannibals on the left bank of s

 template< class F >
 void for_each_feature( std::uint64_t s , F && f ) const {
  f( std::size_t( missionaries( s ) ) );
  f( std::size_t( f_m ) + 1 + cannibals( s ) );
  f( std::size_t( f_m ) + f_c + ( boat_left( s ) ? 2 : 3 ) );
  }

/*------------------------------ OTHER HASHES ------------------------------*/
 /// the perfect hash: the packed state

 [[nodiscard]] std::uint64_t perfect_hash( std::uint64_t s ) const noexcept {
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the state with m missionaries and c cannibals on the left bank
 static std::uint64_t pack( int m , int c , bool left ) noexcept {
  return( std::uint64_t( m ) | ( std::uint64_t( c ) << 16 ) |
          ( std::uint64_t( left ? 0 : 1 ) << 32 ) );
  }

 // true if m missionaries and c cannibals on a bank are safe
 static bool safe( int m , int c ) noexcept {
  return( ( m == 0 ) || ( m >= c ) );
  }

 int f_m;                                        ///< the missionaries
 int f_c;                                        ///< the cannibals
 int f_b;                                        ///< the seats of the boat

 };  // end( class( MissionariesCannibals ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MissionariesCannibals */

/*--------------------------------------------------------------------------*/
/*-------------------- End File MissionariesCannibals.h --------------------*/
/*--------------------------------------------------------------------------*/
