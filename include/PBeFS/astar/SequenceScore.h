/*--------------------------------------------------------------------------*/
/*-------------------------- File SequenceScore.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SequenceScore, which implements the
 * Domain concept [see Domain.h] for the sliding-tile puzzle with the
 * sequence score heuristic of Nilsson, an inadmissible heuristic.
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

#ifndef __PBeFS_SequenceScore
 #define __PBeFS_SequenceScore
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdlib>

#include "PBeFS/astar/SlidingTile.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SequenceScore ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the sliding-tile puzzle with the sequence score heuristic of Nilsson
/** The heuristic is h = M + W S, M being the Manhattan distance and S the
 * sequence score: on the 3 x 3 board S counts 1 if the center does not hold
 * its goal tile and 2 for each tile of the border not followed, clockwise,
 * by the tile that follows it in the goal; on the other boards S counts 1 if
 * the first cell does not hold its goal tile and 2 for each tile not
 * followed, in row-major order, by the tile that follows it in the goal
 * (Nilsson, "Problem-solving methods in artificial intelligence", 1971).
 * The heuristic is not admissible: A* finds a solution expanding fewer
 * nodes, but not necessarily an optimal one. */

template< int W , int H >
class SequenceScore {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using Base = SlidingTile< W , H >;                  ///< the puzzle
 using state_type = typename Base::state_type;       ///< a packed board
 using cost_type = typename Base::cost_type;         ///< the cost type
 using board_type = typename Base::board_type;       ///< a board

 static constexpr int cells = Base::cells;           ///< the cells

 /// constructor, taking the start and the goal boards

 explicit SequenceScore( const board_type & start ,
                         const board_type & goal = Base::korf_goal() )
  : f_base( start , goal ) , v_goal( goal ) {
  for( int c = 0 ; c < cells ; ++c )
   v_goal_cell[ goal[ c ] ] = c;
  if constexpr( ( W == 3 ) && ( H == 3 ) ) {
   const int ring[ 8 ] = { 0 , 1 , 2 , 5 , 8 , 7 , 6 , 3 };
   for( int k = 0 ; k < 8 ; ++k )
    v_next[ ring[ k ] ] = ring[ ( k + 1 ) % 8 ];
   v_next[ 4 ] = -1;
   }
  else {
   for( int c = 0 ; c < cells - 1 ; ++c )
    v_next[ c ] = c + 1;
   v_next[ cells - 1 ] = -1;
   }
  // the tile that follows t: the next of its goal cell, wrapping around to
  // the first tile of the sequence after the last one
  for( int t = 1 ; t < cells ; ++t ) {
   const int c = v_goal_cell[ t ];
   int n = v_next[ c ];
   if( n < 0 )
    n = ( W == 3 ) && ( H == 3 ) ? 0 : 1;
   v_follower[ t ] = goal[ n ];
   }
  }

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_base.initial_state() );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( f_base.is_goal( s ) );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( f_base.hash( s ) );
  }

 /// M + W S

 [[nodiscard]] cost_type heuristic( state_type s ) const noexcept {
  const int center = ( W == 3 ) && ( H == 3 ) ? 4 : 0;
  int m = 0 , sc = Base::tile( s , center ) != v_goal[ center ] ? 1 : 0;
  for( int c = 0 ; c < cells ; ++c ) {
   const int t = Base::tile( s , c );
   if( t == 0 )
    continue;
   const int g = v_goal_cell[ t ];
   m += std::abs( g % W - c % W ) + std::abs( g / W - c / W );
   if( ( v_next[ c ] >= 0 ) &&
       ( Base::tile( s , v_next[ c ] ) != v_follower[ t ] ) )
    sc += 2;
   }
  return( m + W * sc );
  }

 template< class F >
 void for_each_successor( state_type s , cost_type , F && f ) const {
  f_base.for_each_successor( s , 0 , [ & ]( state_type t , cost_type c ,
                                            cost_type ) {
   f( t , c , heuristic( t ) );
   } );
  }

 /// the underlying puzzle
 [[nodiscard]] const Base & base( void ) const noexcept { return( f_base ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 Base f_base;                                    ///< the puzzle
 board_type v_goal;                              ///< the goal board
 int v_goal_cell[ cells ] = {};                  ///< the goal cell of a tile
 int v_next[ cells ] = {};       ///< the next cell in the sequence, or -1
 int v_follower[ cells ] = {};          ///< the goal follower of a tile

 };  // end( class( SequenceScore ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SequenceScore */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SequenceScore.h ------------------------*/
/*--------------------------------------------------------------------------*/
