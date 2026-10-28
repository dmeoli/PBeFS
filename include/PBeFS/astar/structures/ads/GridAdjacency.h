/*--------------------------------------------------------------------------*/
/*-------------------------- File GridAdjacency.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template GridAdjacency of the set of data
 * structures of ADS [see Structures.h], the free cells next to a cell of a
 * GridMap, on the GridGraph of ADS, whose adjacency is computed from a copy
 * of the costs of the cells and never stored; Grid.h includes it in place
 * of its own.
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
/// the free cells next to a cell of a GridMap, on the GridGraph of ADS
/** The interface is that of the GridAdjacency of Grid.h, and so are the
 * moves; the cells are visited in increasing order of index. The row of
 * the cell s = ( x , y ) starts at s - x, so a neighbour t is in the row
 * above, in the same one or in the row below by comparing t with that
 * start, and its column follows with no division. */

template< int Moves >
class GridAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor, copying the costs of the cells of the map

 explicit GridAdjacency( const GridMap & map )
  : f_w( map.width() ) ,
    f_g( map.height() , map.width() , map.costs() ,
         Moves == 4 ? ADS::GridMoves::FOUR : ADS::GridMoves::EIGHT ,
         ADS::Corners::NO_CUT ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , x' , y' , d ) for each free cell t = ( x' , y' ) next to
 /// the free cell s = ( x , y ), d being true if the move is to a corner

 template< class F >
 void for_each( std::uint32_t s , int x , int y , F && f ) const {
  const int u = int( s ) , row = u - x;
  f_g.forEachAdj( u , [ & ]( int v ) {
   const int dy = v < row ? -1 : ( v >= row + f_w ? 1 : 0 );
   const int dx = v - u - dy * f_w;
   f( std::uint32_t( v ) , x + dx , y + dy , ( dx != 0 ) && ( dy != 0 ) );
   } );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 int f_w;                                         ///< the width of the map
 ADS::GridGraph< double , std::uint8_t > f_g;     ///< the graph of the map

 };  // end( class( GridAdjacency ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ads_GridAdjacency */

/*--------------------------------------------------------------------------*/
/*------------------------ End File GridAdjacency.h ------------------------*/
/*--------------------------------------------------------------------------*/
