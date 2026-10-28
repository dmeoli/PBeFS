/*--------------------------------------------------------------------------*/
/*-------------------------- File GridAdjacency.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template GridAdjacency of the set of data
 * structures of ADS [see Structures.h], the moves from a cell of a GridMap
 * and their costs, on the GridGraph of ADS, whose adjacency is computed
 * from a copy of the costs of the cells and never stored, a move paying
 * the cost of the cell it leaves; Grid.h includes it in place of its own.
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

#ifndef __PBeFS_ads_GridAdjacency
 #define __PBeFS_ads_GridAdjacency
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>

#include <optional>

#include <type_traits>

#include "ADS/graph/GridGraph.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS GridAdjacency ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the moves from a cell of a GridMap and their costs, on the GridGraph of
/// ADS
/** The interface is that of the GridAdjacency of Grid.h, and so are the
 * moves and their costs (MoveCost::LEFT, with \f$ \sqrt{2} \f$ to a
 * corner); the cells are visited in increasing order of index. The row of
 * the cell s = ( x , y ) starts at s - x, so a neighbour t is in the row
 * above, in the same one or in the row below by comparing t with that
 * start, and its column follows with no division. */

template< int Moves >
class GridAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// the cost of a move: an integer with 4 moves, a real with 8
 using cost_type = std::conditional_t< Moves == 4 , int , double >;

 /// constructor, copying the costs of the cells of the map

 explicit GridAdjacency( const GridMap & map )
  : f_w( map.width() ) ,
    f_g( map.height() , map.width() , map.costs() ,
         Moves == 4 ? ADS::GridMoves::FOUR : ADS::GridMoves::EIGHT ,
         ADS::Corners::NO_CUT , cost_type( 1 ) ,
         Moves == 4 ? std::optional< cost_type >() :
                      std::optional< cost_type >( cost_type( grid_sqrt2 ) ) ,
         ADS::MoveCost::LEFT ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , x' , y' , c ) for each free cell t = ( x' , y' ) next to
 /// the free cell s = ( x , y ), c being the cost of the move

 template< class F >
 void for_each( std::uint32_t s , int x , int y , F && f ) const {
  const int u = int( s ) , row = u - x;
  f_g.forEachEdge( u , [ & ]( int v , cost_type c ) {
   const int dy = v < row ? -1 : ( v >= row + f_w ? 1 : 0 );
   f( std::uint32_t( v ) , v - u - dy * f_w + x , y + dy , c );
   } );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 int f_w;                                         ///< the width of the map
 ADS::GridGraph< cost_type , std::uint8_t > f_g;  ///< the graph of the map

 };  // end( class( GridAdjacency ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_GridAdjacency */

/*--------------------------------------------------------------------------*/
/*------------------------ End File GridAdjacency.h ------------------------*/
/*--------------------------------------------------------------------------*/
