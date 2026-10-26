/*--------------------------------------------------------------------------*/
/*---------------------------- File WaterJugs.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class WaterJugs, which implements the Domain concept
 * [see Domain.h] for the water jugs problem of exercise 3.6(d) of Russell
 * and Norvig, "Artificial Intelligence: A Modern Approach" (3rd edition,
 * 2010), which is the PourProblem of aima-python
 * (github.com/aimacode/aima-python).
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

#ifndef __PBeFS_WaterJugs
 #define __PBeFS_WaterJugs
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <stdexcept>

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
/*---------------------------- CLASS WaterJugs -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the water jugs problem, with up to 4 jugs of at most 65535 units
/** The jugs are empty at the start, and each move fills a jug from the
 * faucet, empties it on the ground, or pours it into another jug until
 * either the first is empty or the second is full; the goal is a jug
 * holding exactly the given level. Every move costs 1, and the heuristic is
 * 0 at a goal and 1 elsewhere, which is consistent. A state packs the
 * level of the jug \f$ j \f$ in the bits
 * \f$ 16 j , \ldots , 16 j + 15 \f$, and the word is
 * also its perfect hash.
 *
 * The features of Zobrist hashing are "the jug j holds l units", one for
 * each jug and each level from 0 to its capacity. */

class WaterJugs {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint64_t;               ///< the packed levels
 using cost_type = int;                          ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of the jugs of the given capacities, the goal being a jug
 /// holding goal units

 WaterJugs( std::vector< int > capacity , int goal )
  : v_cap( std::move( capacity ) ) , f_goal( goal ) {
  if( v_cap.empty() || ( v_cap.size() > 4 ) )
   throw( std::invalid_argument( "WaterJugs::WaterJugs: 1 to 4 jugs" ) );
  for( int c : v_cap ) {
   if( ( c <= 0 ) || ( c > 65535 ) )
    throw( std::invalid_argument( "WaterJugs::WaterJugs: bad capacity" ) );
   v_first.push_back( f_features );
   f_features += std::size_t( c ) + 1;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the level of the jug j in the state s

 [[nodiscard]] static int level( std::uint64_t s , std::size_t j ) noexcept {
  return( int( ( s >> ( 16 * j ) ) & 0xFFFF ) );
  }

 /// the number of jugs

 [[nodiscard]] std::size_t jugs( void ) const noexcept {
  return( v_cap.size() );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] std::uint64_t initial_state( void ) const noexcept {
  return( 0 );
  }

 [[nodiscard]] bool is_goal( std::uint64_t s ) const noexcept {
  for( std::size_t j = 0 ; j < v_cap.size() ; ++j )
   if( level( s , j ) == f_goal )
    return( true );
  return( false );
  }

 [[nodiscard]] std::uint64_t hash( std::uint64_t s ) const noexcept {
  return( mix64( s ) );
  }

 [[nodiscard]] int heuristic( std::uint64_t s ) const noexcept {
  return( is_goal( s ) ? 0 : 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h( t ) ) for each fill, empty and pour that changes s

 template< class F >
 void for_each_successor( std::uint64_t s , int , F && f ) const {
  const auto with = [ & ]( std::uint64_t t , std::size_t j , int l ) {
   return( ( t & ~( std::uint64_t( 0xFFFF ) << ( 16 * j ) ) ) |
           ( std::uint64_t( l ) << ( 16 * j ) ) );
   };
  const auto move = [ & ]( std::uint64_t t ) {
   if( t != s )
    f( t , 1 , heuristic( t ) );
   };
  for( std::size_t i = 0 ; i < v_cap.size() ; ++i ) {
   move( with( s , i , v_cap[ i ] ) );                       // fill i
   move( with( s , i , 0 ) );                                // empty i
   for( std::size_t j = 0 ; j < v_cap.size() ; ++j )         // i into j
    if( i != j ) {
     const int q = std::min( level( s , i ) , v_cap[ j ] - level( s , j ) );
     move( with( with( s , i , level( s , i ) - q ) , j ,
                 level( s , j ) + q ) );
     }
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: one per jug and level

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( f_features );
  }

 /// calls f( i ) for the feature i of the level of each jug

 template< class F >
 void for_each_feature( std::uint64_t s , F && f ) const {
  for( std::size_t j = 0 ; j < v_cap.size() ; ++j )
   f( v_first[ j ] + std::size_t( level( s , j ) ) );
  }

/*------------------------------ OTHER HASHES ------------------------------*/
 /// the perfect hash: the packed state

 [[nodiscard]] std::uint64_t perfect_hash( std::uint64_t s ) const noexcept {
  return( s );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 std::vector< int > v_cap;                       ///< the capacities
 int f_goal;                                     ///< the level to reach
 std::vector< std::size_t > v_first;  ///< the first feature of each jug
 std::size_t f_features = 0;                     ///< the number of features

 };  // end( class( WaterJugs ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_WaterJugs */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File WaterJugs.h --------------------------*/
/*--------------------------------------------------------------------------*/
