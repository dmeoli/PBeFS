/*--------------------------------------------------------------------------*/
/*--------------------------- File MultiAStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template MultiAStar, a portfolio of A*
 * searches: each worker runs A* on its own, with its own order of the
 * ties, and the first one to end gives the Solution.
 *
 * The workers share nothing but a flag, set by the first worker that ends:
 * the worker \f$ i > 0 \f$ generates the successors of each state in an
 * order drawn by a generator of its own, seeded by the seed of the options
 * and by \f$ i \f$, so that the nodes of equal f and g enter its OPEN list
 * in another order and the ties are broken otherwise, while the worker 0
 * runs the plain A*; once the flag is set, the next state that a worker
 * selects for expansion ends its search, whose result is dropped. Every
 * worker is optimal, and the portfolio pays when the expansions of A* in
 * the last f layer depend much on the ties, at the price of \f$ p \f$ times
 * its memory; with one worker it is A*.
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

#ifndef __PBeFS_MultiAStar
 #define __PBeFS_MultiAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstdint>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/HDAStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS MultiAStar ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a portfolio of A* searches with their own ties, the first one winning
/** The Solution is that of the worker that ended first, with the counters
 * of all the workers summed, expanded_by and stored_by holding those of
 * each; the node_limit holds for each worker. */

template< Domain D , ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class MultiAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the MultiAStar)

 explicit MultiAStar( const D & domain , const ParallelOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  std::atomic< int > winner{ -1 };
  std::vector< SolutionT > sol( p );
  run_workers( f_ex , p , [ & ]( unsigned i ) {
   const View v( f_domain , winner , i == 0 ? 0 :
                 ADS::splitmix64( f_options.seed ^
                                  ( std::uint64_t( i ) << 32 ) ) | 1 );
   sol[ i ] = AStar< View , Open >( v , f_options ).solve();
   int none = -1;
   winner.compare_exchange_strong( none , int( i ) );
   } );

  SearchStats st;
  for( const SolutionT & x : sol ) {
   const SearchStats & s = x.stats;
   st.expanded += s.expanded;
   st.reexpanded += s.reexpanded;
   st.generated += s.generated;
   st.duplicates += s.duplicates;
   st.stored += s.stored;
   st.expanded_by.push_back( s.expanded );
   st.stored_by.push_back( s.stored );
   }
  SolutionT res = std::move( sol[ unsigned( winner.load() ) ] );
  res.stats = std::move( st );
  res.stats.seconds = clock.seconds();
  return( res );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the Domain seen by a worker: the successors in the order of its
 /// generator (in that of the Domain if its seed is 0), and every state a
 /// goal once another worker has ended
 struct View {

  using state_type = State;
  using cost_type = Cost;

  View( const D & dom , const std::atomic< int > & w , std::uint64_t s )
   : d( dom ) , winner( w ) , seed( s ) , rng( s ) {}

  const D & d;
  const std::atomic< int > & winner;
  std::uint64_t seed;

  [[nodiscard]] State initial_state( void ) const {
   return( d.initial_state() );
   }
  [[nodiscard]] bool is_goal( const State & s ) const {
   return( ( winner.load( std::memory_order_relaxed ) >= 0 ) ||
           d.is_goal( s ) );
   }
  [[nodiscard]] Cost heuristic( const State & s ) const {
   return( Cost( d.heuristic( s ) ) );
   }
  [[nodiscard]] std::uint64_t hash( const State & s ) const {
   return( d.hash( s ) );
   }

  template< class F >
  void for_each_successor( const State & s , Cost h , F && f ) const {
   if( seed == 0 ) {
    d.for_each_successor( s , h , std::forward< F >( f ) );
    return;
    }
   kids.clear();
   d.for_each_successor( s , h , [ & ]( const State & t , Cost c ,
                                        Cost ht ) {
    kids.push_back( Kid{ t , c , ht , no_op } );
    } );
   shuffle();
   for( const Kid & k : kids )
    f( k.s , k.c , k.h );
   }

  template< class F > requires OperatorDomain< D >
  void for_each_successor( const State & s , Cost h , std::uint8_t op ,
                           F && f ) const {
   if( seed == 0 ) {
    d.for_each_successor( s , h , op , std::forward< F >( f ) );
    return;
    }
   kids.clear();
   d.for_each_successor( s , h , op , [ & ]( const State & t , Cost c ,
                                             Cost ht , std::uint8_t ot ) {
    kids.push_back( Kid{ t , c , ht , ot } );
    } );
   shuffle();
   for( const Kid & k : kids )
    f( k.s , k.c , k.h , k.op );
   }

  // a successor
  struct Kid {
   State s;
   Cost c;
   Cost h;
   std::uint8_t op;
   };

  // the successors in a random order, by Fisher and Yates on a xorshift
  void shuffle( void ) const {
   for( std::size_t i = kids.size() ; i > 1 ; --i ) {
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    std::swap( kids[ i - 1 ] , kids[ rng % i ] );
    }
   }

  mutable std::vector< Kid > kids;              // the successors of a state
  mutable std::uint64_t rng;                    // the generator

  };  // end( struct( View ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 ParallelOptions f_options;                            ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( MultiAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MultiAStar */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MultiAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
