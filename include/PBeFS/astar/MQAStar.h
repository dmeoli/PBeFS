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
 * 13           m ← m + |F| − 1;  push the entries of F in OPEN
 * 14     return the path to the goal of cost z, if z < ∞
 * \endcode
 *
 * An entry popped with \f$ f \ge z \f$ is dropped, since the entries of
 * OPEN are not ordered across its queues and a better one may still be
 * there. The test at line 8 is exact: at line 13, \f$ m \f$ counts the
 * entries of \f$ F \f$ and drops the one popped in one atomic step,
 * before the pushes, so that \f$ m = 0 \f$ means that OPEN is empty and
 * that no worker is expanding, whatever the order of the pops; then every
 * node generated with \f$ f < z \f$ was expanded, which is the condition
 * of the page "Parallel A*" for an optimal \f$ z \f$. One atomic step
 * per expansion, rather than one before and one after the pushes, halves
 * the contention on the line of \f$ m \f$, which all the workers share.
 * The order of the pops only changes how many nodes are expanded before
 * \f$ z \f$ is optimal, and how many are reopened. The weight, the
 * duplicate dropping and the anytime search change the test of line 12
 * and those of the expansion as in SPAStar, an entry popped being dropped
 * if \f$ f' \geq z \f$ or if \f$ z \f$ prunes it [see Weight.h].
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

#ifndef __PBeFS_MQAStar
 #define __PBeFS_MQAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstdint>

#include <memory>

#include <mutex>

#include <optional>

#include <ranges>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/ConcurrentClosed.h"

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
/** The OPEN list of SPAStar, a priority queue under one lock, is the
 * bottleneck of that search when the expansions are fast; here it is the
 * MultiQueue of ADS (Rihani, Sanders and Dementiev, 2015), with \f$ 2 p \f$
 * sequential queues each under its own lock, so that the workers rarely
 * wait for one another. A pop takes the better of the least entries of two
 * queues drawn at random: in the sequential model with single pushes the
 * expected number of entries of OPEN with a smaller f is less than
 * \f$ \frac{5}{6} \cdot 2 p \f$ (the pushes of all the successors at once
 * relax it as a stickiness would), and therefore the
 * order of the expansions departs from that of A* by little more than it
 * does in SPAStar, where \f$ p \f$ nodes are expanded at once anyway.
 * Postnikova et al. (2022) measure the MultiQueue as the scheduler of
 * parallel searches on graphs.
 *
 * CLOSED is the ConcurrentClosed of SPAStar, the incumbent is atomic,
 * and the search ends when the atomic counter of the entries in OPEN or
 * being expanded falls to 0 [see the pseudocode above]. The stickiness of
 * the MultiQueue is 1, i.e., each operation draws its two queues, which
 * keeps the bound on the rank error above. A node may be expanded before
 * its g is optimal, and then reexpanded, as in SPAStar.
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
  : f_domain( domain ) , f_options( options ) , f_w( options ) {
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
  Shared sh( f_domain , p , f_options.seed , clock );
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh );

  {  // the root, inserted by the first worker before the start
   const State s0 = f_domain.initial_state();
   sh.h0 = Cost( f_domain.heuristic( s0 ) );
   NodeT * root = w[ 0 ]->update( s0 , Cost( 0 ) , sh.h0 , nullptr , no_op );
   auto hd = sh.open.handle();
   hd.push( key( f_w , root->g , root->h ) ,
            OpenEntry< NodeT >{ root , root->g } );
   }

  if( f_w.unit() )               // the weight 1 known at compile time
   run_workers( f_ex , p , [ & ]( unsigned i ) {
    w[ i ]->run( UnitWeight< Cost >() );
    } );
  else
   run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run( f_w ); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.stored += x->f_local.stored();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_local.stored() );
   }
  sol.improvements = std::move( sh.improvements );
  finish_solution( sol , sh.goal , sh.limit_hit.load() , sh.h0 ,
                   f_w , f_domain );
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 using Closed = ConcurrentClosed< D , Lock >;      ///< the CLOSED list

 /// the key of OPEN: the key of w [see Weight], and with integer costs
 /// that key in the high half of a word and the complement of g in the low
 /// one, so that the ties go to the deeper node, as in AStar
 static constexpr bool kPacked = std::is_integral_v< Cost >;
 using Key = std::conditional_t< kPacked , std::uint64_t , Cost >;

 template< class W >
 static Key key( const W & w , Cost g , Cost h ) {
  if constexpr( kPacked )
   return( ( std::uint64_t( w.key( g , h ) ) << 32 ) |
           ( ~std::uint64_t( g ) & 0xFFFFFFFFu ) );
  else
   return( w.key( g , h ) );
  }

 static Cost f_of( Key k ) {
  if constexpr( kPacked )
   return( Cost( k >> 32 ) );
  else
   return( k );
  }

 /// the OPEN list: the key, the node and its g as the item
 using Mq = ADS::MultiQueue< Key , OpenEntry< NodeT > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {
  Shared( const D & d , unsigned p , std::uint64_t seed ,
          const Stopwatch & c )
   : open( 2 * std::size_t( p ) , 1 , seed ) , closed( d ) , clock( c ) {}

  Mq open;                                      // the OPEN list
  // the entries in OPEN or being expanded
  alignas( cache_line ) std::atomic< std::size_t > pending{ 1 };
  Closed closed;                                // the CLOSED list
  alignas( cache_line ) std::atomic< Cost > incumbent{
                                                   infinite_cost< Cost >() };
  std::atomic< bool > done{ false };
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  NodeT * goal = nullptr;                       // the best goal found
  // the incumbents of an anytime search, under goal_mutex
  std::vector< typename SolutionT::Improvement > improvements;
  const Stopwatch & clock;                      // started at the start
  Cost h0{};                                    // the h of the root
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const MQAStar & a , Shared & s ) : f_a( a ) , f_sh( s ) {}

  // the main loop
  template< class W >
  void run( const W & w ) {
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
    if( w.live( f_of( e->key ) , f_sh.incumbent.load() ) )
     expand( w , e->item , fresh );
    // the entries of F counted and the one popped dropped in one step,
    // before the pushes, so that m never misses an entry [see line 13]
    if( fresh.size() != 1 )
     f_sh.pending.fetch_add( fresh.size() - 1 , std::memory_order_acq_rel );
    auto entries = fresh | std::views::transform(
                            [ & ]( const OpenEntry< NodeT > & x ) {
                             return( typename Mq::Entry{
                                      key( w , x.g , x.node->h ) ,
                                      x } ); } );
    hd.pushAll( entries.begin() , entries.end() );
    }
   }

  // expands the node of e unless the entry is stale; the nodes to put in
  // OPEN, with the g they had under their lock, are appended to fresh
  template< class W >
  void expand( const W & w , const OpenEntry< NodeT > & e ,
               std::vector< OpenEntry< NodeT > > & fresh ) {
   NodeT * n = e.node;
   State s;
   Cost g , h;
   std::uint8_t op;
   {
    std::lock_guard< Lock > guard( f_sh.closed.lock( n ) );
    if( e.stale() || w.prunes( e.g , n->h , f_sh.incumbent.load(
                                               std::memory_order_relaxed ) ) )
     return;
    n->closed = true;
    if( n->expanded )
     ++f_st.reexpanded;
    n->expanded = true;
    s = n->state;
    g = n->g;
    h = n->h;
    op = n->op;
    }
   ++f_st.expanded;

   if( f_a.f_domain.is_goal( s ) ) {
    std::lock_guard< std::mutex > guard( f_sh.goal_mutex );
    if( g < f_sh.incumbent.load() ) {
     f_sh.incumbent.store( g );
     f_sh.goal = n;
     if( f_a.f_w.anytime() )
      f_sh.improvements.push_back( { g , f_sh.h0 , f_st.expanded ,
                                     f_sh.clock.seconds() } );
     }
    return;
    }

   // the successors first, but the parent [see OperatorDomain], with the
   // places of CLOSED that their updates read prefetched, so that the
   // misses of all of them overlap
   f_succ.clear();
   for_each_child( f_a.f_domain , s , h , op , [ & ]( const State & t ,
                                                      Cost c , Cost ht ,
                                                      std::uint8_t ot ) {
    ++f_st.generated;
    const Cost gt = g + c;
    if( w.prunes( gt , ht , f_sh.incumbent.load(
                                           std::memory_order_relaxed ) ) )
     return;
    const std::uint64_t hash = f_a.f_domain.hash( t );
    f_sh.closed.prefetch( hash );
    f_succ.push_back( { t , gt , ht , hash , ot } );
    } );
   for( const auto & v : f_succ )
    if( NodeT * x = f_sh.closed.update( v.state , v.hash , v.g , v.h , n ,
                                        v.op , f_local , f_st.duplicates ,
                                        w ) )
     fresh.push_back( OpenEntry< NodeT >{ x , v.g } );

   if( f_a.f_options.node_limit &&
       ( f_local.stored() >= f_a.f_options.node_limit ) ) {
    f_sh.limit_hit.store( true );
    f_sh.done.store( true , std::memory_order_release );
    }
   }

  // records the path to t of cost g through parent and op; returns the
  // node of t if it has to go into OPEN, nullptr if it is a duplicate
  NodeT * update( const State & t , Cost g , Cost h , NodeT * parent ,
                  std::uint8_t op ) {
   return( f_sh.closed.update( t , f_a.f_domain.hash( t ) , g , h ,
                               parent , op , f_local , f_st.duplicates ,
                               f_a.f_w ) );
   }

  const MQAStar & f_a;                          // the search
  Shared & f_sh;                                // the shared data
  typename Closed::Local f_local;               // the nodes it allocated
  std::vector< typename Closed::Successor > f_succ;  // of the expansion
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 ParallelOptions f_options;                            ///< the options
 Weight< Cost > f_w;                                   ///< the weight
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( MQAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MQAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File MQAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
