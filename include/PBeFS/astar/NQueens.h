/*--------------------------------------------------------------------------*/
/*----------------------------- File NQueens.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class NQueens, which implements the Domain concept
 * [see Domain.h] for the n-queens problem in the incremental formulation
 * of section 3.2.1 of Russell and Norvig, "Artificial Intelligence: A
 * Modern Approach" (3rd edition, 2010).
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

#ifndef __PBeFS_NQueens
 #define __PBeFS_NQueens
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <cstring>

#include <stdexcept>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS NQueens ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the n-queens problem, with at most 32 queens
/** A queen is added in the leftmost empty column, in a row where no queen
 * already on the board attacks it; a state is the number \f$ k \f$ of the
 * queens on the board and the row of the queen of each column \f$ c < k \f$,
 * the other rows being 0, with \f$ n \leq 32 \f$. Each queen costs 1, and
 * the heuristic, the number of the empty columns, is exact on the states
 * that lead to a solution: every state has \f$ f = n \f$, and A* with the
 * ties broken in favour of the larger g behaves as a depth-first search.
 *
 * The feature \f$ c n + r \f$ of Zobrist hashing is the queen of the
 * column \f$ c \f$ in the row \f$ r \f$, for \f$ c < k \f$; the states
 * do not fit a 64-bit word for \f$ n > 12 \f$, and the Domain gives no
 * perfect hash. */

class NQueens {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 /// a board: the queens placed and the row of each of them
 struct Board {
  std::array< std::uint8_t , 32 > row;   ///< the row of each column
  std::uint8_t k;                        ///< the queens on the board

  bool operator==( const Board & ) const = default;
  };

 using state_type = Board;                       ///< a board
 using cost_type = int;                          ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of the problem with \f$ n \f$ queens,
 /// \f$ 1 \leq n \leq 32 \f$

 explicit NQueens( int n ) : f_n( n ) {
  if( ( n < 1 ) || ( n > 32 ) )
   throw( std::invalid_argument( "NQueens::NQueens: n is 1 to 32" ) );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] Board initial_state( void ) const noexcept {
  return( Board{ {} , 0 } );
  }

 [[nodiscard]] bool is_goal( const Board & s ) const noexcept {
  return( s.k == f_n );
  }

 [[nodiscard]] std::uint64_t hash( const Board & s ) const noexcept {
  std::uint64_t w[ 4 ];
  std::memcpy( w , s.row.data() , sizeof( w ) );
  std::uint64_t z = ADS::splitmix64( s.k );
  for( std::uint64_t x : w )
   z = ADS::splitmix64( z ^ x );
  return( z );
  }

 [[nodiscard]] int heuristic( const Board & s ) const noexcept {
  return( f_n - s.k );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , 1 , h - 1 ) for each row of the column \f$ k \f$ free of
 /// attacks

 template< class F >
 void for_each_successor( const Board & s , int h , F && f ) const {
  const int k = s.k;
  for( int r = 0 ; r < f_n ; ++r ) {
   bool free = true;
   for( int c = 0 ; free && ( c < k ) ; ++c ) {
    const int rc = s.row[ c ];
    free = ( rc != r ) && ( std::abs( rc - r ) != k - c );
    }
   if( free ) {
    Board t = s;
    t.row[ k ] = std::uint8_t( r );
    t.k = std::uint8_t( k + 1 );
    f( t , 1 , h - 1 );
    }
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: one per pair ( column , row )

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( std::size_t( f_n ) * f_n );
  }

 /// calls f( c n + r ) for the queen of each column \f$ c < k \f$, in
 /// the row \f$ r \f$

 template< class F >
 void for_each_feature( const Board & s , F && f ) const {
  for( int c = 0 ; c < s.k ; ++c )
   f( std::size_t( c ) * f_n + s.row[ c ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the size of the board

 [[nodiscard]] int size( void ) const noexcept { return( f_n ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 int f_n;                                        ///< the size of the board

 };  // end( class( NQueens ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_NQueens */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File NQueens.h ---------------------------*/
/*--------------------------------------------------------------------------*/
