/*--------------------------------------------------------------------------*/
/*----------------------------- File MQAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template MQAStar, MultiQueue A*: the workers
 * share one CLOSED list, as in SPAStar, and an OPEN list that is a relaxed
 * priority queue, the MultiQueue of Rihani, Sanders and Dementiev, in place
 * of a priority queue under one lock.
 *
 * The search is that of SPAStar.h, with the notation of the page
 * "Best-first search and A*" of the documentation, but for OPEN and the
 * test of termination: OPEN is the MultiQueue of ADS, whose pop takes the
 * better of the least entries of two queues drawn at random, and therefore
 * not always the least of all, and the counter \f$ m \f$ of the entries
 * that are in OPEN or being expanded, which is atomic, takes the place of
 * busy and of the lock \f$ L \f$:
 *
 * \code{.unparsed}
 *  1  procedure MQA*(s₀, p)
 *  2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {(f(s₀), g(s₀), s₀)}
 *  3     z ← ∞;  m ← 1;  done ← false
 *  4     in parallel, on each worker i = 0, ..., p − 1 do
 *  5        while not done do
 *  6           e ← a pop of OPEN, possibly nothing
 *  7           if no entry was popped then
 *  8              if m = 0 then done ← true
 *  9              else wait a little
 * 10              continue
 * 11           (f, g, u) ← e;  F ← ∅
 * 12           if f < z then expand u as at lines 13 to 25 of SPA*
 * 13           m ← m + |F|;  push the entries of F in OPEN
 * 14           m ← m − 1
 * 15     return the path to the goal of cost z, if z < ∞
 * \endcode
 *
 * An entry popped with \f$ f \ge z \f$ is dropped, since the entries of
 * OPEN are not ordered across its queues and a better one may still be
 * there. The test at line 8 is exact: \f$ m \f$ grows before the pushes
 * of line 13 and falls only at line 14, after them, so that \f$ m = 0 \f$
 * means that OPEN is empty and that no worker is expanding, whatever the
 * order of the pops; then every node generated with \f$ f < z \f$ was
 * expanded, which is the condition of the page "Parallel A*" for an
 * optimal \f$ z \f$. The order of the pops only changes how many nodes
 * are expanded before \f$ z \f$ is optimal, and how many are reopened.
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

#ifndef __PBeFS_MQAStar
 #define __PBeFS_MQAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <memory>

#include <mutex>

#include <optional>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/HDAStar.h"

#include "ADS/concurrent/MultiQueue.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS MQAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// MultiQueue A*: one relaxed OPEN list and one CLOSED list shared by the
/// workers
/** The OPEN list of MQAStar, a priority queue under one lock, is the
 * bottleneck of that search when the expansions are fast; here it is the
 * MultiQueue of ADS (Rihani, Sanders and Dementiev, 2015), with \f$ 2 p \f$
 * sequential queues each under its own lock, so that the workers rarely
 * wait for one another. A pop takes the better of the least entries of two
 * queues drawn at random: the expected number of entries of OPEN with a
 * smaller f is less than \f$ \frac{5}{6} \cdot 2 p \f$, and therefore the
 * order of the expansions departs from that of A* by little more than it
 * does in MQAStar, where \f$ p \f$ nodes are expanded at once anyway.
 * Postnikova et al. (2022) measure the MultiQueue as the scheduler of
 * parallel searches on graphs.
 *
 * CLOSED is split into 256 shards as in MQAStar, the incumbent is atomic,
 * and the search ends when the atomic counter of the entries in OPEN or
 * being expanded falls to 0 [see the pseudocode above]. The stickiness of
 * the MultiQueue is 1, i.e., each operation draws its two queues, which
 * keeps the bound on the rank error above. A node may be expanded before
 * its g is optimal, and then reexpanded, as in MQAStar.
 *
 * \see A. Postnikova, N. Koval, G. Nadiradze, D. Alistarh, Multi-queues
 *      can be state-of-the-art priority schedulers, ACM SIGPLAN Symposium
 *      on Principles and Practice of Parallel Programming (PPoPP),
 *      353-367, 2022, doi:10.1145/3503221.3508432 */

template< Domain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class MQAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the MQAStar)

 explicit MQAStar( const D & domain , const ParallelOptions & options = {} )
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
  Shared sh( p , f_options.seed );
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh );

  {  // the root, inserted by the first worker before the start
   const State s0 = f_domain.initial_state();
   NodeT * root = w[ 0 ]->update( s0 , Cost( 0 ) ,
                                  Cost( f_domain.heuristic( s0 ) ) ,
                                  nullptr );
   auto hd = sh.open.handle();
   hd.push( root->g + root->h , OpenEntry< NodeT >{ root , root->g } );
   }

  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  if( sh.limit_hit.load() )
   sol.status = Status::node_limit;
  else
   if( sh.goal ) {
    sol.status = Status::solved;
    sol.cost = sh.goal->g;
    sol.path = trace_path( sh.goal );
    }
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 static constexpr unsigned c_shards = 256;         ///< shards of CLOSED

 /// the OPEN list: f as the key, the node and its g as the item
 using Mq = ADS::MultiQueue< Cost , OpenEntry< NodeT > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a shard of CLOSED and its lock

 struct alignas( cache_line ) Shard {
  Lock lock;
  NodeTable< NodeT > table;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {
  explicit Shared( unsigned p , std::uint64_t seed )
   : open( 2 * std::size_t( p ) , 1 , seed ) {}

  Mq open;                                      // the OPEN list
  // the entries in OPEN or being expanded
  alignas( cache_line ) std::atomic< std::size_t > pending{ 1 };
  std::unique_ptr< Shard[] > shards{ new Shard[ c_shards ] };
  alignas( cache_line ) std::atomic< Cost > incumbent{
                                                   infinite_cost< Cost >() };
  std::atomic< bool > done{ false };
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  NodeT * goal = nullptr;                       // the best goal found
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const MQAStar & a , Shared & s ) : f_a( a ) , f_sh( s ) {}

  // the shard of a hash value
  Shard & shard( std::uint64_t h ) {
   return( f_sh.shards[ h >> 56 ] );
   }

  // the main loop
  void run( void ) {
   ADS::Backoff backoff;
   auto hd = f_sh.open.handle();
   std::vector< OpenEntry< NodeT > > fresh;
   while( ! f_sh.done.load( std::memory_order_acquire ) ) {
    const auto e = hd.tryPop();
    if( ! e ) {
     if( f_sh.pending.load( std::memory_order_acquire ) == 0 ) {
      f_sh.done.store( true , std::memory_order_release );
      break;
      }
     backoff.pause();
     continue;
     }
    backoff.reset();
    fresh.clear();
    if( e->key < f_sh.incumbent.load() )
     expand( e->item , fresh );
    f_sh.pending.fetch_add( fresh.size() , std::memory_order_relaxed );
    for( const auto & x : fresh )
     hd.push( x.g + x.node->h , x );
    f_sh.pending.fetch_sub( 1 , std::memory_order_acq_rel );
    }
   }

  // expands the node of e unless the entry is stale; the nodes to put in
  // OPEN, with the g they had under the lock of their shard, are appended
  // to fresh
  void expand( const OpenEntry< NodeT > & e ,
               std::vector< OpenEntry< NodeT > > & fresh ) {
   NodeT * n = e.node;
   State s;
   Cost g , h;
   {
    Shard & sd = shard( f_a.f_domain.hash( n->state ) );
    std::lock_guard< Lock > guard( sd.lock );
    if( e.stale() )
     return;
    n->closed = true;
    if( n->expanded )
     ++f_st.reexpanded;
    n->expanded = true;
    s = n->state;
    g = n->g;
    h = n->h;
    }
   ++f_st.expanded;

   if( f_a.f_domain.is_goal( s ) ) {
    std::lock_guard< std::mutex > guard( f_sh.goal_mutex );
    if( g < f_sh.incumbent.load() ) {
     f_sh.incumbent.store( g );
     f_sh.goal = n;
     }
    return;
    }

   f_a.f_domain.for_each_successor( s , h , [ & ]( const State & t , Cost c ,
                                                   Cost ht ) {
    ++f_st.generated;
    const Cost gt = g + c;
    if( gt + ht >= f_sh.incumbent.load( std::memory_order_relaxed ) )
     return;
    if( NodeT * x = update( t , gt , ht , n ) )
     fresh.push_back( OpenEntry< NodeT >{ x , gt } );
    } );

   if( f_a.f_options.node_limit &&
       ( f_pool.size() >= f_a.f_options.node_limit ) ) {
    f_sh.limit_hit.store( true );
    f_sh.done.store( true , std::memory_order_release );
    }
   }

  // records the path to t of cost g through parent; returns the node of t
  // if it has to go into OPEN, nullptr if it is a duplicate
  NodeT * update( const State & t , Cost g , Cost h , NodeT * parent ) {
   const std::uint64_t hash = f_a.f_domain.hash( t );
   Shard & sd = shard( hash );
   std::lock_guard< Lock > guard( sd.lock );
   auto [ x , fresh ] = sd.table.find_or_insert( t , hash , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ t , g , h , parent , false , false };
    return( y );
    } );
   if( fresh )
    return( x );
   if( g < x->g ) {
    x->g = g;
    x->parent = parent;
    x->closed = false;
    return( x );
    }
   ++f_st.duplicates;
   return( nullptr );
   }

  const MQAStar & f_a;                          // the search
  Shared & f_sh;                                // the shared data
  NodePool< NodeT > f_pool;                     // the nodes it allocated
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 ParallelOptions f_options;                            ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( MQAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MQAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File MQAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
