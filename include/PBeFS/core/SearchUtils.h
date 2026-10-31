/*--------------------------------------------------------------------------*/
/*--------------------------- File SearchUtils.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the small helpers shared by the searches: a stopwatch,
 * the placement of the workers on the CPUs and their run on an executor,
 * a prefetch, and the random numbers of the workers; the mixing of the
 * bits of a word, the backoff of the waiting workers, the placement on the
 * CPUs and the allocator of huge pages are those of ADS.
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

#ifndef __PBeFS_SearchUtils
 #define __PBeFS_SearchUtils
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <chrono>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <vector>

#include "ADS/concurrent/Backoff.h"

#include "ADS/dictionary/RandomHashing.h"

#include "ADS/memory/HugePageAllocator.h"

#include "ADS/parallel/Affinity.h"

#include "ADS/parallel/Executor.h"

#include "ADS/parallel/ThreadPool.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CONSTANTS AND FUNCTIONS -------------------------*/
/*--------------------------------------------------------------------------*/

/// the size of a cache line, used to keep the data of the workers apart
inline constexpr std::size_t cache_line = 64;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// pins the calling thread to the i-th CPU it may run on, if PBEFS_PIN is set
/** If the environment variable PBEFS_PIN is set to a non-empty value other
 * than 0, binds the calling thread to the i-th CPU (modulo their number) of
 * the affinity mask of the process, where the system allows it [see
 * ADS::pinThread()]; otherwise does nothing. The workers of every search
 * call it at their start, so that they are placed the same way on every
 * executor. */

inline void pin_worker( unsigned i ) noexcept
{
 static const bool pin = [] {
  const char * v = std::getenv( "PBEFS_PIN" );
  return( v && *v && ( v[ 0 ] != '0' ) );
  }();
 if( pin )
  static_cast< void >( ADS::pinThread( i ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the parameters shared by the parallel searches
/** The number of workers, and, for the searches that send their nodes to
 * their owners, the messages in a batch, the steps (expansions or
 * evaluations) between two flushes of the outboxes and the seed of the
 * owners; each search sets its own defaults of batch and flush_period. */

struct ParallelParams {

 unsigned threads = 1;          ///< the number of workers
 std::size_t batch = 64;        ///< nodes packed in a message
 unsigned flush_period = 64;    ///< steps between flushes of outboxes
 std::uint64_t seed = 1;        ///< seed of the owners

 /// sets to 1 any parameter that is 0
 void normalize( void ) noexcept {
  if( threads == 0 )
   threads = 1;
  if( batch == 0 )
   batch = 1;
  if( flush_period == 0 )
   flush_period = 1;
  }

 };  // end( struct( ParallelParams ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// runs the workers f( 0 ) , ... , f( p - 1 ) at once on an Executor
/** Calls f( i ) for i = 0 , ... , p - 1 by the run() of *ex, or, if ex is
 * nullptr, of an Executor of p threads built for the purpose, each call on
 * a thread of its own after pin_worker( i ) but on the thread of the caller,
 * which stays where it is, and returns when all the calls have returned. */

template< ADS::WorkerExecutor Executor , class F >
void run_workers( Executor * ex , unsigned p , F && f )
{
 const auto body = [ & f ]( std::size_t i ) {
  if( ! ADS::isRunCaller() )
   pin_worker( unsigned( i ) );
  f( unsigned( i ) );
  };
 if( ex )
  ex->run( p , body );
 else {
  Executor own( p );
  own.run( p , body );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// hints the processor to bring the cache line of p in, for reading

inline void prefetch( const void * p ) noexcept
{
 #if defined( __GNUC__ ) || defined( __clang__ )
  __builtin_prefetch( p , 0 , 1 );
 #else
  ( void ) p;
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a real number in [ 0 , 1 ) from the next word of a 64-bit generator
/** The distributions of \<random\> are not the same on every standard
 * library, while the words of std::mt19937_64 are; this and random_below()
 * draw from those words directly, so that a seed gives the same instance on
 * every platform. */

template< class Gen >
[[nodiscard]] double random_unit( Gen & gen )
{
 return( double( std::uint64_t( gen() ) >> 11 ) * 0x1.0p-53 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an integer in [ 0 , n ) from the next word of a 64-bit generator, n > 0

template< class Gen >
[[nodiscard]] std::uint64_t random_below( Gen & gen , std::uint64_t n )
{
 return( std::uint64_t( gen() ) % n );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Stopwatch -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a wall-clock stopwatch, started at construction

class Stopwatch {

 public:

 Stopwatch( void ) : f_start( std::chrono::steady_clock::now() ) {}

 /// the seconds elapsed since the construction
 [[nodiscard]] double seconds( void ) const {
  return( std::chrono::duration< double >(
                       std::chrono::steady_clock::now() - f_start ).count() );
  }

 private:

 std::chrono::steady_clock::time_point f_start;     ///< the start time

 };  // end( class( Stopwatch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SearchUtils */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SearchUtils.h -------------------------*/
/*--------------------------------------------------------------------------*/
