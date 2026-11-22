/*--------------------------------------------------------------------------*/
/*------------------------ File MultiObjectiveGrid.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template MultiObjectiveGrid, pathfinding on a
 * grid whose moves cost a vector of K components, which implements the
 * MultiObjectiveDomain concept [see MultiObjective.h].
 *
 * The cells are those of a GridMap [see Grid.h], the blocked ones (of cost
 * 0) excluded, and the moves go to the 4 or the 8 neighboring free cells,
 * a diagonal one only if both the cells beside it are free. The i-th cost
 * of the move between the cells a and b is an integer drawn uniformly in
 * \f$ [ \mathit{lo} , \mathit{hi} ] \f$ by a hash of the unordered pair
 * \f$ \{ a , b \} \f$, of i and of a seed, so that a move costs as much as
 * its inverse and the grid needs no memory for the costs: these are the
 * grids of Ren et al. (SoCS 2022), four-connected with costs in
 * \f$ [ 1 , 10 ] \f$, there on the empty grid of 10 x 10 cells and on maps
 * of the benchmarks of Sturtevant (2012). The heuristic is
 * \f$ h_i( u ) = \mathit{lo} \, \delta( u , t ) \f$, \f$ \delta \f$ being the
 * Manhattan distance with 4 moves and the Chebyshev distance with 8, which
 * a move lowers by 1 at most, so that h is consistent; set_ideal_heuristic()
 * replaces it with the ideal point of each cell [see MultiObjectiveGraph.h].
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

#ifndef __PBeFS_MultiObjectiveGrid
 #define __PBeFS_MultiObjectiveGrid
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <stdexcept>

#include <vector>

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/MultiObjectiveGraph.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MultiObjectiveGrid ------------------------*/
/*--------------------------------------------------------------------------*/
/// a grid whose moves cost K random integers, searched from a cell to a
/// cell
/** The states are the cells \f$ y W + x \f$ of a map of width W [see the
 * file header]; the constructor throws std::invalid_argument if the start
 * or the goal is blocked or outside the map, if the moves are neither 4
 * nor 8, or if \f$ 0 \leq \mathit{lo} \leq \mathit{hi} \f$ does not
 * hold. */

template< class Cost , std::size_t K >
class MultiObjectiveGrid {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint32_t;                  ///< a cell
 using cost_type = Cost;                            ///< a component
 static constexpr std::size_t objectives = K;       ///< the objectives
 using Vec = std::array< Cost , K >;                ///< a vector of costs

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor: the map, the moves (4 or 8), the start ( sx , sy ), the
 /// goal ( gx , gy ), the range of the costs and their seed

 MultiObjectiveGrid( const GridMap & map , int moves , int sx , int sy ,
                     int gx , int gy , Cost lo = 1 , Cost hi = 10 ,
                     std::uint64_t seed = 1 )
  : f_map( map ) , f_moves( moves ) , f_lo( lo ) , f_hi( hi ) ,
    f_seed( seed ) {
  if( ( map.cost( sx , sy ) == 0 ) || ( map.cost( gx , gy ) == 0 ) ||
      ( ( moves != 4 ) && ( moves != 8 ) ) || ( lo < Cost( 0 ) ) ||
      ( hi < lo ) )
   throw( std::invalid_argument( "MultiObjectiveGrid: bad start, goal, "
                                 "moves or costs" ) );
  f_start = cell( sx , sy );
  f_goal = cell( gx , gy );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// replaces the heuristic with the ideal point of each cell, by K
 /// Dijkstra's algorithms from the goal

 void set_ideal_heuristic( void ) {
  const std::size_t n = std::size_t( f_map.width() ) * f_map.height();
  v_h.assign( n , Vec{} );
  for( std::size_t i = 0 ; i < K ; ++i ) {
   // the moves are their own inverses, at the same cost
   const auto d = distances_to< Cost , K >( n , f_goal , i ,
                                            [ & ]( state_type v , auto && f ) {
    for_each_move( v , f );
    } );
   for( std::size_t u = 0 ; u < n ; ++u )
    v_h[ u ][ i ] = d[ u ];
   }
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_goal );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 [[nodiscard]] Vec heuristic( state_type s ) const noexcept {
  if( ! v_h.empty() )
   return( v_h[ s ] );
  const int w = f_map.width();
  const int dx = std::abs( int( s % w ) - int( f_goal % w ) );
  const int dy = std::abs( int( s / w ) - int( f_goal / w ) );
  Vec h;
  h.fill( f_lo * Cost( f_moves == 4 ? dx + dy : std::max( dx , dy ) ) );
  return( h );
  }

 template< class F >
 void for_each_successor( state_type s , const Vec & , F && f ) const {
  for_each_move( s , [ & ]( state_type t , const Vec & c ) {
   f( t , c , heuristic( t ) );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cell ( x , y )

 [[nodiscard]] state_type cell( int x , int y ) const noexcept {
  return( state_type( y ) * state_type( f_map.width() ) + state_type( x ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the vector of the costs of the move between the cells a and b

 [[nodiscard]] Vec move_cost( state_type a , state_type b ) const noexcept {
  const std::uint64_t pair = ( std::uint64_t( std::min( a , b ) ) << 32 ) |
                             std::max( a , b );
  Vec c;
  for( std::size_t i = 0 ; i < K ; ++i ) {
   const std::uint64_t r = ADS::splitmix64( ADS::splitmix64( pair ^ f_seed )
                                            + i );
   c[ i ] = f_lo + Cost( r % std::uint64_t( f_hi - f_lo + 1 ) );
   }
  return( c );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // calls f( t , c ) for each move from s to t of the costs c
 template< class F >
 void for_each_move( state_type s , F && f ) const {
  const int w = f_map.width() , x = int( s % w ) , y = int( s / w );
  static constexpr int dx[ 8 ] = { 1 , -1 , 0 , 0 , 1 , 1 , -1 , -1 };
  static constexpr int dy[ 8 ] = { 0 , 0 , 1 , -1 , 1 , -1 , 1 , -1 };
  for( int k = 0 ; k < f_moves ; ++k ) {
   const int u = x + dx[ k ] , v = y + dy[ k ];
   if( f_map.cost( u , v ) == 0 )
    continue;
   if( ( k >= 4 ) && ( ( f_map.cost( u , y ) == 0 ) ||
                       ( f_map.cost( x , v ) == 0 ) ) )
    continue;
   const state_type t = cell( u , v );
   f( t , move_cost( s , t ) );
   }
  }

 const GridMap & f_map;                             ///< the map
 int f_moves;                                       ///< 4 or 8
 Cost f_lo;                                         ///< the least cost
 Cost f_hi;                                         ///< the largest cost
 std::uint64_t f_seed;                              ///< the seed
 state_type f_start = 0;                            ///< the start
 state_type f_goal = 0;                             ///< the goal
 std::vector< Vec > v_h;                            ///< the ideal points

 };  // end( class( MultiObjectiveGrid ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MultiObjectiveGrid */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MultiObjectiveGrid.h ---------------------*/
/*--------------------------------------------------------------------------*/
