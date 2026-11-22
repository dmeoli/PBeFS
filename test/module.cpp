/*--------------------------------------------------------------------------*/
/*---------------------------- File module.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test of the C++20 module pbefs, built when PBEFS_BUILD_MODULE is ON [see
 * CMakeLists.txt].
 *
 * The program includes no header of PBeFS: it imports the module, and with
 * what that exports it solves an 8-puzzle three moves away from the goal of
 * Korf by A* and by IDA*, and a random knapsack by the sequential
 * Branch-and-Bound, and checks the costs against the optimal ones (3, and
 * the value of the dynamic program of the knapsack). A name that the module
 * fails to export stops the build, and a wrong result the run.
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

#include <iostream>

import pbefs;

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 int failed = 0;
 auto check = [ & ]( bool ok , const char * what ) {
  if( ! ok ) {
   std::cerr << "FAILED " << what << std::endl;
   ++failed;
   }
  };

 // the blank moved right, right and down from the goal of Korf
 const EightPuzzle d( EightPuzzle::parse( "1 2 5 3 4 0 6 7 8" ) );
 const auto a = astar( d );
 check( a.solved() && ( a.cost == 3 ) , "astar" );
 const auto i = idastar( d );
 check( i.solved() && ( i.cost == 3 ) , "idastar" );

 // a weakly correlated knapsack of 20 items
 const Knapsack k = Knapsack::random( 1 , 20 , 100 , 1 );
 const auto b = branch_and_bound( k );
 check( ( b.status == BBStatus::optimal ) &&
        ( long( b.value ) == k.optimum() ) , "branch_and_bound" );

 if( failed == 0 )
  std::cout << "the module pbefs: all the checks passed" << std::endl;
 return( failed == 0 ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File module.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
