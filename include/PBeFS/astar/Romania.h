/*--------------------------------------------------------------------------*/
/*----------------------------- File Romania.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class Romania, the road map of Romania of figure 3.2
 * of Russell and Norvig, "Artificial Intelligence: A Modern Approach" (3rd
 * edition, 2010), as an ExplicitGraph [see ExplicitGraph.h].
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

#ifndef __PBeFS_Romania
 #define __PBeFS_Romania
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cmath>

#include <cstdint>

#include <stdexcept>

#include <string>

#include "PBeFS/astar/ExplicitGraph.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Romania ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a trip on the road map of Romania
/** The cities are Arad, Bucharest, Craiova, Drobeta, Eforie, Fagaras,
 * Giurgiu, Hirsova, Iasi, Lugoj, Mehadia, Neamt, Oradea, Pitesti, Rimnicu
 * (Vilcea), Sibiu, Timisoara, Urziceni, Vaslui and Zerind, and the roads
 * cost their length in km. Towards Bucharest the heuristic is the
 * straight-line distance of figure 3.22 of the book; towards any other
 * city it is the distance between the coordinates of the cities in
 * aima-python (github.com/aimacode/aima-python), which on this map never
 * exceeds the length of a road and is therefore consistent. The costs are
 * real numbers, since the latter heuristic is not an integer. */

class Romania : public ExplicitGraph< double > {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 /// the cities, in alphabetical order
 enum City : std::uint32_t { Arad , Bucharest , Craiova , Drobeta , Eforie ,
                             Fagaras , Giurgiu , Hirsova , Iasi , Lugoj ,
                             Mehadia , Neamt , Oradea , Pitesti , Rimnicu ,
                             Sibiu , Timisoara , Urziceni , Vaslui , Zerind ,
                             NCities };

 /// a road and its length in km
 struct Road { City a , b; int km; };

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of the trip from the city from to the city to

 Romania( City from = Arad , City to = Bucharest )
  : ExplicitGraph< double >( NCities , from , to ) {
  for( const Road & r : roads() )
   add_edge( r.a , r.b , r.km );
  for( std::uint32_t c = 0 ; c < NCities ; ++c )
   set_heuristic( c , to == Bucharest
                       ? double( c_to_bucharest[ c ] )
                       : std::hypot( c_xy[ c ][ 0 ] - c_xy[ to ][ 0 ] ,
                                     c_xy[ c ][ 1 ] - c_xy[ to ][ 1 ] ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the name of the city c

 [[nodiscard]] static const char * name( std::uint32_t c ) {
  static const char * const n[ NCities ] = {
   "Arad" , "Bucharest" , "Craiova" , "Drobeta" , "Eforie" , "Fagaras" ,
   "Giurgiu" , "Hirsova" , "Iasi" , "Lugoj" , "Mehadia" , "Neamt" ,
   "Oradea" , "Pitesti" , "Rimnicu" , "Sibiu" , "Timisoara" , "Urziceni" ,
   "Vaslui" , "Zerind" };
  if( c >= NCities )
   throw( std::invalid_argument( "Romania::name: no such city" ) );
  return( n[ c ] );
  }

 /// the city called name, NCities if there is none

 [[nodiscard]] static City find( const std::string & city_name ) {
  for( std::uint32_t c = 0 ; c < NCities ; ++c )
   if( city_name == name( c ) )
    return( City( c ) );
  return( NCities );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the roads of figure 3.2 of the book

 [[nodiscard]] static const std::array< Road , 23 > & roads( void ) {
  static const std::array< Road , 23 > r = { {
   { Oradea , Zerind , 71 } , { Oradea , Sibiu , 151 } ,
   { Zerind , Arad , 75 } , { Arad , Sibiu , 140 } ,
   { Arad , Timisoara , 118 } ,
   { Timisoara , Lugoj , 111 } , { Lugoj , Mehadia , 70 } ,
   { Mehadia , Drobeta , 75 } , { Drobeta , Craiova , 120 } ,
   { Craiova , Rimnicu , 146 } , { Craiova , Pitesti , 138 } ,
   { Sibiu , Rimnicu , 80 } , { Sibiu , Fagaras , 99 } ,
   { Rimnicu , Pitesti , 97 } , { Fagaras , Bucharest , 211 } ,
   { Pitesti , Bucharest , 101 } , { Bucharest , Giurgiu , 90 } ,
   { Bucharest , Urziceni , 85 } , { Urziceni , Hirsova , 98 } ,
   { Hirsova , Eforie , 86 } , { Urziceni , Vaslui , 142 } ,
   { Vaslui , Iasi , 92 } , { Iasi , Neamt , 87 } } };
  return( r );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the straight-line distance to Bucharest (figure 3.22)
 static constexpr int c_to_bucharest[ NCities ] = {
  366 , 0 , 160 , 242 , 161 , 176 , 77 , 151 , 226 , 244 , 241 , 234 , 380 ,
  100 , 193 , 253 , 329 , 80 , 199 , 374 };

 /// the coordinates of the cities in aima-python
 static constexpr int c_xy[ NCities ][ 2 ] = {
  { 91 , 492 } , { 400 , 327 } , { 253 , 288 } , { 165 , 299 } ,
  { 562 , 293 } , { 305 , 449 } , { 375 , 270 } , { 534 , 350 } ,
  { 473 , 506 } , { 165 , 379 } , { 168 , 339 } , { 406 , 537 } ,
  { 131 , 571 } , { 320 , 368 } , { 233 , 410 } , { 207 , 457 } ,
  { 94 , 410 } , { 456 , 350 } , { 509 , 444 } , { 108 , 531 } };

 };  // end( class( Romania ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Romania */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File Romania.h ---------------------------*/
/*--------------------------------------------------------------------------*/
