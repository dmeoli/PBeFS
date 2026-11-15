/*--------------------------------------------------------------------------*/
/*---------------------------- File symbolic.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The program that packaging/test builds against the installed library of
 * the symbolic searches, PBeFS::Symbolic, when it is there: the 8-puzzle of
 * figure 3.28 of Russell and Norvig, whose optimal solution has 26 moves,
 * solved by BDDA* bidirectionally, and the 181440 states that its start
 * reaches.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/symbolic.h"

#include <iostream>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const BddRuntime rt( 2 );
 const SymbolicSlidingTile< 3 , 3 > d(
  EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" ) );
 SymbolicOptions o;
 o.direction = SymbolicDirection::bidirectional;
 const auto s = symbolic_search( d , o );
 const auto r = symbolic_reachable( d );
 std::cout << "BDDA* " << s.cost << ", states " << r.states << std::endl;
 return( ( s.cost == 26 ) && ( r.states == 181440 ) ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File symbolic.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
