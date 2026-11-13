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

#include <cstdint>

#include <numeric>

#include <vector>

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS OF SearchStats -------------------------*/
/*--------------------------------------------------------------------------*/

// the largest of v over its mean, 1 if v is empty or all 0

static double max_over_mean( const std::vector< std::uint64_t > & v )
{
 if( v.empty() )
  return( 1 );
 const double sum = std::accumulate( v.begin() , v.end() , 0.0 );
 if( sum == 0 )
  return( 1 );
 const auto mx = *std::max_element( v.begin() , v.end() );
 return( double( mx ) * double( v.size() ) / sum );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double SearchStats::load_balance( void ) const
{
 return( max_over_mean( expanded_by ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double SearchStats::storage_balance( void ) const
{
 return( max_over_mean( stored_by ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double SearchStats::communication_overhead( void ) const
{
 return( generated == 0 ? 0 : double( sent ) / double( generated ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double SearchStats::search_overhead( std::uint64_t sequential ) const
{
 return( sequential == 0 ? 0 : double( expanded ) / double( sequential ) - 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File Solution.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
