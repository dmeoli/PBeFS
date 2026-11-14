/*--------------------------------------------------------------------------*/
/*----------------------------- File PWAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template PWAStar, the Parallel Weighted A*
 * that is the baseline of the edge-based parallel searches of Mukherjee,
 * Aine and Likhachev (SoCS 2022, section "Evaluation"; RA-L 2022): weighted
 * A* whose expansion evaluates the edges of the state in parallel.
 *
 * With the notation of EPASE.h, \f$ f' = g + w h \f$:
 *
 * \code{.unparsed}
 *  1  procedure PwA*(s₀, p)
 *  2     OPEN ← {s₀}, with g(s₀) ← 0;  CLOSED ← ∅
 *  3     while OPEN ≠ ∅ do
 *  4        s ← the state of OPEN of least f' (larger g on ties)
 *  5        if s is a goal then return the path to s
 *  6        add s to CLOSED
 *  7        the p workers evaluate the edges of s, each taking the next
 *  8           action not yet taken                  ▷ in parallel
 *  9        for each (v, c(s, v)) found, in the order of the actions do
 * 10           if v ∉ CLOSED and g(s) + c(s, v) < g(v) then
 * 11              g(v) ← g(s) + c(s, v);  parent(v) ← s;  put v in OPEN
 * \endcode
 *
 * The order of the expansions, and hence the path, the expansions and the
 * edges evaluated, are those of the sequential weighted A* without
 * reopening, whatever the number of workers: the cost is within a factor
 * \f$ w \f$ of the optimum for a consistent heuristic, optimal for
 * \f$ w = 1 \f$. The speedup is bounded by the branching factor, since
 * the workers evaluate the edges of one state at a time.
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

#ifndef __PBeFS_PWAStar
 #define __PBeFS_PWAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <compare>

#include <condition_variable>

#include <cstddef>

#include <cstdint>

#include <memory>

#include <mutex>

#include <vector>

#include "PBeFS/astar/HDAStar.h"

#include "ADS/priority/DaryHeap.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS PWAStarOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of PWAStar: those of ParallelOptions, whose weight
/// \f$ w \geq 1 \f$ is that of the heuristic, \f$ f' = g + w h \f$

struct PWAStarOptions : public ParallelOptions {

 };  // end( struct( PWAStarOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS PWAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Parallel Weighted A*: the edges of each expansion evaluated in parallel
/** The worker 0 runs the weighted A* and expands; at each expansion the
 * workers, the worker 0 among them, take the actions of the state from a
 * counter and evaluate their edges, and the worker 0 applies the results
 * once all are in. The statistics count in evaluated the edges evaluated,
 * and in expanded_by the edges evaluated by each worker. */

template< EdgeDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class PWAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the PWAStar)

 explicit PWAStar( const D & domain , const PWAStarOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( ! ( f_options.weight >= 1 ) )
   f_options.weight = 1;
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
  Shared sh;
  std::vector< SearchStats > st( p );
  SolutionT sol;
  run_workers( f_ex , p , [ & ]( unsigned i ) {
   if( i == 0 )
    master( sh , st , sol );
   else
    helper( sh , st[ i ] );
   } );
  for( const auto & x : st ) {
   sol.stats.evaluated += x.evaluated;
   sol.stats.expanded_by.push_back( x.evaluated );
   }
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // a mutex whatever the Executor: the workers that wait sleep on a condition
 // variable, while those of a spin lock (the lock of FastFlow) would take
 // the cores of the workers that evaluate the edges
 using Lock = std::mutex;

 // an item of OPEN: f', then larger g, then the order of the pushes, so
 // that the order does not depend on the addresses; stale if the node has
 // been closed or improved since
 struct Item {
  double f = 0;
  double ng = 0;
  std::uint64_t seq = 0;
  NodeT * node = nullptr;

  friend auto operator<=>( const Item & , const Item & ) = default;
  };

 // what the edge of an action gave
 struct Result {
  bool feasible = false;
  State t{};
  Cost c{};
  Cost ht{};
  };

 // the round of evaluations shared by the workers
 struct Shared {
  Lock lock;                                    // protects round and done
  std::condition_variable start;                // a new round, or done
  std::condition_variable end;                  // the round is over
  std::uint64_t round = 0;                      // the current round
  unsigned active = 0;                          // helpers in the round
  bool done = false;                            // the search is over
  State state{};                                // the state expanded
  Cost h{};                                     // its heuristic
  std::vector< std::uint32_t > actions;         // its actions
  std::vector< Result > results;                // one per action
  std::atomic< std::size_t > next{ 0 };         // the next action to take
  std::atomic< std::size_t > left{ 0 };         // the actions not done
  };

 // takes and evaluates actions until none is left
 void claim( Shared & sh , SearchStats & st ) const {
  const std::size_t n = sh.actions.size();
  for( std::size_t i ; ( i = sh.next.fetch_add( 1 ) ) < n ; ) {
   Result & r = sh.results[ i ];
   f_domain.evaluate_edge( sh.state , sh.h , sh.actions[ i ] ,
                           [ & ]( const State & t , Cost c , Cost ht ) {
    r = Result{ true , t , c , ht };
    } );
   ++st.evaluated;
   sh.left.fetch_sub( 1 , std::memory_order_acq_rel );
   }
  }

 // a worker other than 0: evaluates in each round; the worker 0 changes
 // the round only when no helper is in it
 void helper( Shared & sh , SearchStats & st ) const {
  std::unique_lock< Lock > guard( sh.lock );
  std::uint64_t seen = 0;
  while( true ) {
   while( ( sh.round == seen ) && ! sh.done )
    sh.start.wait( guard );
   if( sh.done )
    break;
   seen = sh.round;
   ++sh.active;
   guard.unlock();
   claim( sh , st );
   guard.lock();
   --sh.active;
   sh.end.notify_one();
   }
  }

 // the worker 0: the weighted A*, with a round of evaluations per
 // expansion
 void master( Shared & sh , std::vector< SearchStats > & st ,
              SolutionT & sol ) const {
  const double w = f_options.weight;
  const unsigned p = f_options.threads;
  NodePool< NodeT > pool;
  NodeTable< NodeT > table;
  ADS::DaryHeap< Item > open;
  SearchStats & s0 = st[ 0 ];
  std::uint64_t pushes = 0;
  const auto push = [ & ]( NodeT * x ) {
   open.push( Item{ double( x->g ) + w * double( x->h ) , - double( x->g ) ,
                    pushes++ , x } );
   };
  {
   const State r = f_domain.initial_state();
   NodeT * x = pool.allocate();
   *x = NodeT{ r , Cost( 0 ) , Cost( f_domain.heuristic( r ) ) , nullptr ,
               false , false };
   table.find_or_insert( r , f_domain.hash( r ) , [ & ] { return( x ); } );
   push( x );
   }
  NodeT * goal = nullptr;
  bool limit = false;
  while( ! open.empty() ) {
   const Item it = open.pop();
   NodeT * s = it.node;
   if( s->closed || ( - double( s->g ) != it.ng ) )
    continue;
   if( f_domain.is_goal( s->state ) ) {
    goal = s;
    break;
    }
   s->closed = s->expanded = true;
   ++s0.expanded;
   {
    std::unique_lock< Lock > guard( sh.lock , std::defer_lock );
    if( p > 1 ) {
     guard.lock();
     while( sh.active )
      sh.end.wait( guard );
     }
    sh.actions.clear();
    f_domain.for_each_action( s->state , [ & ]( std::size_t a ) {
     sh.actions.push_back( std::uint32_t( a ) );
     } );
    sh.results.assign( sh.actions.size() , Result{} );
    sh.state = s->state;
    sh.h = s->h;
    sh.next.store( 0 );
    sh.left.store( sh.actions.size() );
    if( ( p > 1 ) && ( sh.actions.size() > 1 ) ) {
     ++sh.round;
     sh.start.notify_all();
     }
    if( p > 1 )
     guard.unlock();
    claim( sh , s0 );
    if( p > 1 ) {
     guard.lock();
     while( sh.left.load() || sh.active )
      sh.end.wait( guard );
     }
    }
   for( const Result & r : sh.results ) {
    if( ! r.feasible )
     continue;
    ++s0.generated;
    const Cost g = s->g + r.c;
    auto [ x , fresh ] = table.find_or_insert( r.t , f_domain.hash( r.t ) ,
                                               [ & ] {
     NodeT * y = pool.allocate();
     *y = NodeT{ r.t , g , r.ht , s , false , false };
     return( y );
     } );
    if( fresh )
     push( x );
    else
     if( ! x->closed && ( g < x->g ) ) {
      x->g = g;
      x->parent = s;
      push( x );
      }
     else
      ++s0.duplicates;
    }
   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    limit = true;
    break;
    }
   }
  {
   std::unique_lock< Lock > guard( sh.lock );
   sh.done = true;
   sh.start.notify_all();
   }
  SearchStats & t = sol.stats;
  t.expanded = s0.expanded;
  t.generated = s0.generated;
  t.duplicates = s0.duplicates;
  t.stored = pool.size();
  t.stored_by.push_back( pool.size() );
  if( limit )
   sol.status = Status::node_limit;
  else
   if( goal ) {
    sol.status = Status::solved;
    sol.cost = goal->g;
    sol.path = trace_path( goal );
    }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 PWAStarOptions f_options;                             ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( PWAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PWAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File PWAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
