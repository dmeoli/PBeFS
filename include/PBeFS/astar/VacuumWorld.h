/*--------------------------------------------------------------------------*/
/*--------------------------- File VacuumWorld.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class VacuumWorld, which implements the Domain
 * concept [see Domain.h] for the vacuum world of section 3.2.1 of Russell
 * and Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition,
 * 2010), on a row of n cells instead of 2.
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

#ifndef __PBeFS_VacuumWorld
 #define __PBeFS_VacuumWorld
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <bit>

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
/*--------------------------- CLASS VacuumWorld ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the vacuum world on a row of at most 32 cells
/** The robot moves Left or Right, or Sucks the dirt of its cell, each
 * action costing 1, and the goal is a clean row; the heuristic, the number
 * of the dirty cells, is consistent. A state packs the dirty cells in the
 * bits \f$ 0 , \ldots , 31 \f$ and the cell of the robot in the bits
 * \f$ 32 , \ldots , 37 \f$,
 * and the word is also its perfect hash.
 *
 * The features of Zobrist hashing are the dirt of the cell \f$ i \f$, for
 * \f$ i < n \f$, and the robot in the cell \f$ r \f$, numbered
 * \f$ n + r \f$. The abstraction keeps the dirt and forgets the robot: the
 * abstract features of Abstract Zobrist hashing are the dirty cells (the
 * \f$ n \f$ features of the robot being one abstract feature, which every
 * state has), and the abstract state of AHDA* is the set of the dirty cells,
 * so that a move keeps the owner and only a Suck changes it. The nblocks of
 * PBNF are the cells of the robot, the neighbours of a cell being its
 * abstract successors, and the heuristic of PA*SE between two states is the
 * number of the cells dirty in the first and clean in the second. */

class VacuumWorld {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint64_t;               ///< the packed state
 using cost_type = int;                          ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of the row of n cells, with the robot in the cell start
 /// and the dirty cells given by the bits of dirty

 VacuumWorld( int n , int start , std::uint32_t dirty )
  : f_n( n ) , f_start( start ) , f_dirty( dirty ) {
  if( ( n < 1 ) || ( n > 32 ) || ( start < 0 ) || ( start >= n ) )
   throw( std::invalid_argument( "VacuumWorld::VacuumWorld: bad sizes" ) );
  if( ( n < 32 ) && ( dirty >> n ) )
   throw( std::invalid_argument( "VacuumWorld::VacuumWorld: dirt out of "
                                 "the row" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cell of the robot in the state s

 [[nodiscard]] static int robot( std::uint64_t s ) noexcept {
  return( int( s >> 32 ) );
  }

 /// the dirty cells in the state s

 [[nodiscard]] static std::uint32_t dirt( std::uint64_t s ) noexcept {
  return( std::uint32_t( s ) );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] std::uint64_t initial_state( void ) const noexcept {
  return( pack( f_start , f_dirty ) );
  }

 [[nodiscard]] bool is_goal( std::uint64_t s ) const noexcept {
  return( dirt( s ) == 0 );
  }

 [[nodiscard]] std::uint64_t hash( std::uint64_t s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 [[nodiscard]] int heuristic( std::uint64_t s ) const noexcept {
  return( std::popcount( dirt( s ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h( t ) ) for Suck, Left and Right, if possible

 template< class F >
 void for_each_successor( std::uint64_t s , int h , F && f ) const {
  const int r = robot( s );
  const std::uint32_t d = dirt( s );
  if( d & ( std::uint32_t( 1 ) << r ) )                       // Suck
   f( pack( r , d & ~( std::uint32_t( 1 ) << r ) ) , 1 , h - 1 );
  if( r > 0 )                                                 // Left
   f( pack( r - 1 , d ) , 1 , h );
  if( r + 1 < f_n )                                           // Right
   f( pack( r + 1 , d ) , 1 , h );
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the dirty cells and the cells of the robot

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( 2 * std::size_t( f_n ) );
  }

 /// calls f( i ) for each dirty cell i and f( n + r ) for the robot in r

 template< class F >
 void for_each_feature( std::uint64_t s , F && f ) const {
  for( std::uint32_t d = dirt( s ) ; d ; d &= d - 1 )
   f( std::size_t( std::countr_zero( d ) ) );
  f( std::size_t( f_n + robot( s ) ) );
  }

 /// the number of abstract features: the dirty cells and the robot

 [[nodiscard]] std::size_t num_abstract_features( void ) const noexcept {
  return( std::size_t( f_n ) + 1 );
  }

 /// the abstract feature of a feature: the dirty cell, or n for the robot

 [[nodiscard]] std::size_t abstract_feature( std::size_t i ) const noexcept
 {
  return( i < std::size_t( f_n ) ? i : std::size_t( f_n ) );
  }

/*--------------------------- STATE ABSTRACTIONS ---------------------------*/
 /// the abstract state of AHDA*: the dirty cells

 [[nodiscard]] std::uint64_t abstract_state( std::uint64_t s ) const noexcept
 {
  return( dirt( s ) );
  }

 /// the number of nblocks: the cells

 [[nodiscard]] std::size_t num_nblocks( void ) const noexcept {
  return( std::size_t( f_n ) );
  }

 /// the nblock of PBNF: the cell of the robot

 [[nodiscard]] std::size_t nblock( std::uint64_t s ) const noexcept {
  return( std::size_t( robot( s ) ) );
  }

 /// calls f( b' ) for each cell b' next to the cell b

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const {
  if( b > 0 )
   f( b - 1 );
  if( b + 1 < std::size_t( f_n ) )
   f( b + 1 );
  }

/*---------------------- OTHER HEURISTICS AND HASHES -----------------------*/
 /// the heuristic of PA*SE: the cells dirty in s and clean in t

 [[nodiscard]] int pairwise_heuristic( std::uint64_t s , std::uint64_t t )
  const noexcept {
  return( std::popcount( dirt( s ) & ~dirt( t ) ) );
  }

 /// the perfect hash: the packed state

 [[nodiscard]] std::uint64_t perfect_hash( std::uint64_t s ) const noexcept {
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of cells

 [[nodiscard]] int size( void ) const noexcept { return( f_n ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the state with the robot in the cell r and the dirty cells d
 static std::uint64_t pack( int r , std::uint32_t d ) noexcept {
  return( ( std::uint64_t( r ) << 32 ) | d );
  }

 int f_n;                                        ///< the number of cells
 int f_start;                                    ///< the cell of the robot
 std::uint32_t f_dirty;                          ///< the dirty cells

 };  // end( class( VacuumWorld ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_VacuumWorld */

/*--------------------------------------------------------------------------*/
/*------------------------- End File VacuumWorld.h -------------------------*/
/*--------------------------------------------------------------------------*/
