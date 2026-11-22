/*--------------------------------------------------------------------------*/
/*------------------------- File ads-backends.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The programs by which makefile-c tells whether a backend of ADS is
 * there: compiled with PBEFS_TRY_HIGHWAY, PBEFS_TRY_HWLOC or
 * PBEFS_TRY_BOOST, each compiles and links only if Highway (its prefetch
 * and its sort, in libhwy and libhwy_contrib), hwloc 2 or the headers of
 * Boost 1.83 or later are found.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/

#if defined( PBEFS_TRY_HIGHWAY )

#include <cstdint>

#include <hwy/cache_control.h>

#include <hwy/contrib/sort/vqsort.h>

int main( void )
{
 std::int64_t k[ 2 ] = { 2 , 1 };
 hwy::Prefetch( k );
 hwy::VQSort( k , 2 , hwy::SortAscending() );
 return( int( k[ 0 ] ) - 1 );
 }

#elif defined( PBEFS_TRY_HWLOC )

#include <hwloc.h>

static_assert( HWLOC_API_VERSION >= 0x00020000 );

int main( void )
{
 hwloc_topology_t t;
 if( hwloc_topology_init( &t ) )
  return( 1 );
 hwloc_topology_destroy( t );
 return( 0 );
 }

#elif defined( PBEFS_TRY_BOOST )

#include <boost/unordered/unordered_flat_map.hpp>

#include <boost/version.hpp>

static_assert( BOOST_VERSION >= 108300 );

int main( void )
{
 boost::unordered_flat_map< int , int > m;
 return( int( m.size() ) );
 }

#endif

/*--------------------------------------------------------------------------*/
/*------------------------ End File ads-backends.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
