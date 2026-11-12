/*--------------------------------------------------------------------------*/
/*--------------------------- File Solution.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the measures of a search computed from its counters:
 * the load balance and the communication overhead of SearchStats.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <numeric>

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS OF SearchStats -------------------------*/
/*--------------------------------------------------------------------------*/

double SearchStats::load_balance( void ) const
{
 if( expanded_by.empty() )
  return( 1 );
 const double sum = std::accumulate( expanded_by.begin() ,
                                     expanded_by.end() , 0.0 );
 if( sum == 0 )
  return( 1 );
 const auto mx = *std::max_element( expanded_by.begin() ,
                                    expanded_by.end() );
 return( double( mx ) * double( expanded_by.size() ) / sum );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double SearchStats::communication_overhead( void ) const
{
 return( generated == 0 ? 0 : double( sent ) / double( generated ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File Solution.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
