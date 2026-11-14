/*--------------------------------------------------------------------------*/
/*---------------------------- File HDLPAStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template HDLPAStar, Hash Distributed Lifelong
 * Planning A*: the incremental search of LPAStar [see LPAStar.h] run by the
 * workers of HDAStar [see HDAStar.h], each of which owns the states that a
 * distribution assigns to it, so that a search is repaired in parallel
 * after a change of the costs; and for the class template HDDStarLite, the
 * same on the reversed Domain with a target that moves, i.e., a parallel
 * D* Lite.
 *
 * The worker \f$ i \f$ owns the states \f$ v \f$ with \f$ \mathrm{owner}( v
 * ) = i \f$, and is the only one that writes their \f$ g \f$, \f$
 * \mathrm{rhs} \f$, parent and place in its queue \f$ U_i \f$; the others
 * read \f$ g( v ) \f$ only, by an atomic load, when they compute an rhs
 * over all the predecessors (a rescan). The goals are joined to a virtual
 * target \f$ t^* \f$ by edges of cost 0, whose owner keeps the last g of
 * each goal in an indexed heap and publishes the least, \f$ B \f$, in an
 * atomic word, so that the test of the end costs a load. Three messages
 * cross the workers, all with absolute values: \f$ \mathrm{Relax}( v , x ,
 * u ) \f$, after an overconsistent expansion of \f$ u \f$, lowers \f$
 * \mathrm{rhs}( v ) \f$ to \f$ x = g( u ) + c( u , v ) \f$ with the parent
 * \f$ u \f$ if \f$ x \f$ is smaller; \f$ \mathrm{Invalidate}( v , u ) \f$,
 * after an underconsistent expansion of \f$ u \f$, makes the owner of
 * \f$ v \f$ rescan \f$ v \f$ if its parent is \f$ u \f$; and
 * \f$ \mathrm{Wake} \f$, sent by the owner of \f$ t^* \f$ to all the others
 * when \f$ B \f$ rises, since the mailbox of ADS makes a passive worker
 * active only by a message. A state is eligible when it is inconsistent
 * and its key is below \f$ [ B + k_m ; B ] \f$, or it is a goal and its
 * key is not above it:
 *
 * \code{.unparsed}
 *  1  procedure WORKER(i, φ)                       ▷ one phase of one episode
 *  2     while the mailbox is not quiescent do
 *  3        for each message m received do RECEIVE(m)
 *  4        if active then
 *  5           u ← the least eligible state of U_i, the stale entries
 *                 dropped and, if φ = 1, the overconsistent ones expanded
 *                 underconsistent in this episode put aside (frozen)
 *  6           if u exists then EXPAND(u, φ) else become passive
 *  7  procedure EXPAND(u, φ)
 *  8     if g(u) > rhs(u) then                         ▷ overconsistent
 *  9        g(u) ← rhs(u), by an atomic store, then the sends
 * 10        SEND Relax(v, g(u) + c(u, v), u) to each successor v,
 *              and Relax(t*, g(u), u) if u is a goal
 * 11     else                                          ▷ underconsistent
 * 12        g(u) ← ∞;  raised(u) ← the episode
 * 13        SEND Invalidate(v, u) to each successor v, and
 *              Invalidate(t*, u) if u is a goal;  RESCAN(u)
 * 14  procedure RECEIVE(m)                             ▷ on the owner of m.v
 * 15     if m.v = t* then update the mirror and B;  if B rose and φ = 1,
 *           SEND Wake to the other workers
 * 16     else if m = Relax(v, x, u), v ≠ s₀ and x < rhs(v) then
 *           rhs(v) ← x;  parent(v) ← u;  fix v in U
 * 17     else if m = Invalidate(v, u) and parent(v) = u then RESCAN(v)
 * 18  procedure RESCAN(v)
 * 19     rhs(v) ← 0 if v = s₀, otherwise the least atomic load of g(w) plus
 *           c(w, v) over the predecessors w, parent(v) ← that w;  fix v in U
 * 20  procedure EPISODE(changes)                       ▷ between barriers
 * 21     each worker applies the changes of the edges into its states, as
 *           LPAStar::edge_changed(), and the owner of t* those of the goals
 * 22     restart, or repair: if an increase met a parent edge or the key of
 *           t* rose, WORKER(·, 1) to quiescence and a barrier; then the
 *           frozen states back in U and WORKER(·, 2) to quiescence
 * \endcode
 *
 * A SEND to the worker itself is applied at once, to another one goes in
 * its outbox; a Relax creates the state it names if it is new, an
 * Invalidate of a state unknown is dropped. The first episode, and every
 * restart, is phase 2 alone from an empty table, i.e., HDA* on the pair
 * \f$ ( g , \mathrm{rhs} ) \f$ with no pruning at insertion. The rule of
 * phase 1, that a state expanded underconsistent waits for phase 2 to be
 * expanded overconsistent ("freeze after raise"), is what makes an episode
 * end: without it two workers may count to infinity, as asynchronous
 * distance-vector routing does, on a cycle whose states have become
 * unreachable. The page of theory "Bounded-suboptimal, anytime and
 * incremental search" states the invariants and proves that at the end of
 * each episode \f$ B \f$ is the cost of a shortest path and the parents
 * from the goal of least g give one, whatever the interleaving of the
 * workers, with each state expanded underconsistent at most once per
 * episode.
 *
 * The repair or restart of an episode is chosen from the changes, the
 * prediction of the cost of the repair and the cost of the last restart
 * [see IncrementalOptions], and a repair that exceeds a budget is stopped
 * and restarted, so that an episode costs at most the repair or about
 * twice the restart.
 *
 * One mailbox of ADS serves all the phases, re-armed at each barrier; the
 * table of each worker is a ConcurrentFlatTable of ADS, which the others
 * read while its owner inserts, emptied by a restart in \f$ O( 1 ) \f$ by
 * its epochs, the nodes of the old table being kept in a list of free
 * nodes for the new one (the arena of the nodes has no epochs). The
 * messages are not coalesced in the outboxes, which the proofs allow, and
 * the audit is that of the search.
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

#ifndef __PBeFS_HDLPAStar
 #define __PBeFS_HDLPAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <barrier>

#include <cmath>

#include <cstdint>

#include <memory>

#include <optional>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <unordered_map>

#include <vector>

#include "ADS/concurrent/ConcurrentFlatTable.h"

#include "ADS/parallel/Channels.h"

#include "ADS/parallel/Mailbox.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/LPAStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS IncrementalOptions ------------------------*/
/*--------------------------------------------------------------------------*/
/// the choice between repairing a search and restarting it
/** With the policy adaptive an episode restarts when the changes that
 * matter (the decreases that lower an rhs and the increases of a parent
 * edge) are more than rho times \f$ \hat S \f$, the mean of the
 * expansions of the restarts, or when the expansions that the repair is
 * predicted to cost exceed gamma times \f$ \hat S \f$, the
 * prediction summing over the heads of the changes an estimate learned
 * from the past repairs, by the place of the head between the start and
 * the goal; a repair that expands more than beta times \f$ \hat S \f$ is
 * stopped and restarted, so that with beta = 1 and \f$ \hat S \f$ exact an
 * episode costs at most \f$ \min( R , 2 S ) \f$ expansions, \f$ R \f$
 * the repair and \f$ S \f$ the restart (the deterministic bound of ski
 * rental). The means weigh the last value by mean_weight. The policies
 * repair and restart always do that. freeze is the rule of phase 1:
 * Freeze::all holds back every state expanded underconsistent in the
 * episode until phase 2, so that none is expanded underconsistent twice;
 * Freeze::raised_parent only those whose parent has been expanded
 * underconsistent too, which ends as well, since a state is lowered in
 * phase 1 only through a parent not yet raised, and expands fewer states.
 * audit checks the invariants of the search at the end of every phase, and
 * makes solve() throw std::logic_error if one fails. */

struct IncrementalOptions {

 /// the choice of an episode
 enum class Policy {
  adaptive ,                                     ///< by the prediction
  repair ,                                       ///< always repair
  restart                                        ///< always restart
  };

 /// the states that phase 1 does not lower
 enum class Freeze {
  all ,                                          ///< every one raised
  raised_parent                                  ///< those whose parent
  };                                             ///< was raised too

 Policy policy = Policy::adaptive;               ///< the choice
 Freeze freeze = Freeze::raised_parent;          ///< the rule of phase 1
 double rho = 0.5;                               ///< the bound on changes
 double gamma = 1;                               ///< on the predicted work
 double beta = 1;                                ///< on the repair
 double mean_weight = 0.5;                       ///< of the last value
 bool audit = false;                             ///< check the invariants

 };  // end( struct( IncrementalOptions ) )

/*--------------------------------------------------------------------------*/
/// the options of HDLPAStar: those of the parallel searches and the choice
/// between repair and restart

struct HDLPAOptions : public ParallelOptions , public IncrementalOptions {
 };  // end( struct( HDLPAOptions ) )

/*--------------------------------------------------------------------------*/

namespace detail {

 /// the queue of the keys of LPA*: buckets for integer costs, a heap
 /// otherwise
 template< class T , class Cost >
 struct lpa_default_queue {
  using type = LPAHeap< T , Cost >;
  };

 template< class T , class Cost > requires std::is_integral_v< Cost >
 struct lpa_default_queue< T , Cost > {
  using type = LPABuckets< T , Cost >;
  };

 }  // end( namespace detail )

/// LPABuckets for integer costs, LPAHeap otherwise

template< class T , class Cost >
using LPADefaultQueue = typename detail::lpa_default_queue< T , Cost >::type;

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS HDLPAStar -----------------------------*/
/*--------------------------------------------------------------------------*/
/// Hash Distributed Lifelong Planning A*: an incremental search repaired in
/// parallel by the owners of the states
/** The use is that of LPAStar: the first solve() searches from scratch, the
 * caller changes the Domain and tells each change by edge_changed() and
 * goal_changed() (or take_events() for an EventDomain), and the next
 * solve() repairs the search, or restarts it [see IncrementalOptions]; the
 * changes of one edge are told once, from the cost before the first to the
 * cost after the last. shift_keys() adds to \f$ k_m \f$ and tells that the
 * heuristic changed [see HDDStarLite], restart() makes the next solve()
 * start from scratch. The distribution must give each state always the
 * same owner (not RandomDistribution); Queue is that of the keys of each
 * worker, LPABuckets for integer costs and LPAHeap otherwise by default.
 * The counters of a Solution are those of its episode: under counts the
 * underconsistent expansions, reexpanded the expansions of a state already
 * expanded in the episode, sent the messages to another worker, messages
 * their batches, rescans and remote_reads the rescans and the loads of a g
 * of another worker that they make, phases the phases run, restarts 1 if
 * the episode restarted, delta the states whose g changed. The constructor
 * throws std::invalid_argument if the options give a weight other than 1. */

template< PredecessorDomain D , class Dist = StateHashDistribution< D > ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          template< class , class > class Queue = LPADefaultQueue >
class HDLPAStar {

 static_assert( ! requires { Dist::per_worker; } ,
                "HDLPAStar: the owner of a state must never change" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain and the distribution (which must
 /// outlive the HDLPAStar)

 HDLPAStar( const D & domain , const Dist & dist ,
            const HDLPAOptions & options = {} )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) {
  if( options.weight != 1 )
   throw( std::invalid_argument( "HDLPAStar: no weight" ) );
  f_options.normalize();
  f_p = f_options.threads;
  f_s0 = domain.initial_state();
  f_s0_owner = f_dist.owner( f_s0 , f_p );
  f_t_owner = unsigned( ADS::splitmix64( f_options.seed ) % f_p );
  for( unsigned i = 0 ; i < f_p ; ++i )
   v_w.push_back( std::make_unique< Worker >( *this , i ) );
  f_mail = std::make_unique< Mailbox >( f_p , f_options.batch );
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// finds the shortest path with the current costs, repairing or
 /// restarting the last search after the changes told

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  ++f_episode;
  f_target_rose = f_target_rose || f_km_changed;
  f_km_changed = false;
  f_audit_error.clear();
  for( auto & w : v_w )
   w->begin_episode();
  for( const Event & e : v_events ) {
   if( e.goal )
    v_w[ f_dist.owner( e.u , f_p ) ]->v_ev.push_back( e );
   else
    v_w[ f_dist.owner( e.v , f_p ) ]->v_ev.push_back( e );
   }
  v_events.clear();

  std::barrier sync( std::ptrdiff_t( f_p ) , Step{ this } );
  f_sync = & sync;
  f_stage = Stage::events;
  run_workers( f_ex , f_p , [ & ]( unsigned i ) { v_w[ i ]->episode(); } );
  f_sync = nullptr;
  if( ! f_audit_error.empty() )
   throw( std::logic_error( "HDLPAStar::solve: " + f_audit_error ) );

  SolutionT sol;
  SearchStats & st = sol.stats;
  Worker & t = *v_w[ f_t_owner ];
  for( auto & x : v_w ) {
   x->end_episode();
   const SearchStats & y = x->f_st;
   st.expanded += y.expanded;
   st.reexpanded += y.reexpanded;
   st.under += y.under;
   st.generated += y.generated;
   st.sent += y.sent;
   st.relaxes += y.relaxes;
   st.invalidates += y.invalidates;
   st.wakes += y.wakes;
   st.rescans += y.rescans;
   st.remote_reads += y.remote_reads;
   st.delta += y.delta;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( y.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  st.messages = batches() - f_batches;
  st.phases = f_phases;
  st.restarts = f_restarted ? 1 : 0;
  std::uint64_t rexp = 0;
  for( auto & x : v_w )
   rexp += x->f_rexp;
  learn( st.expanded , rexp );
  if( ! t.f_mu.empty() ) {
   const LNode * goal = t.v_mnode[ std::size_t( t.f_mu.top().id ) ];
   sol.status = Status::solved;
   sol.cost = goal->g.load( std::memory_order_relaxed );
   sol.bound = sol.cost;
   std::size_t n = 0;
   for( const LNode * x = goal ; x ; x = x->par ) {
    sol.path.push_back( x->state );
    if( x->state == f_s0 )
     break;
    if( ++n > st.stored )
     throw( std::logic_error( "HDLPAStar::solve: a cycle of parents" ) );
    }
   std::reverse( sol.path.begin() , sol.path.end() );
   }
  st.seconds = clock.seconds();
  return( sol );
  }

/*------------------------------ THE CHANGES -------------------------------*/
 /// tells that the cost of the edge ( u , v ) has changed from from to to,
 /// \f$ \infty \f$ standing for an absent edge

 void edge_changed( const State & u , const State & v , Cost from , Cost to )
 {
  if( from == to )
   return;
  if constexpr( cost_changes_of< D >() == CostChanges::decrease_only )
   if( to > from )
    throw( std::invalid_argument( "HDLPAStar: an increase of a Domain"
                                  " that only decreases" ) );
  if constexpr( cost_changes_of< D >() == CostChanges::increase_only )
   if( to < from )
    throw( std::invalid_argument( "HDLPAStar: a decrease of a Domain"
                                  " that only increases" ) );
  v_events.push_back( Event{ u , v , from , to , false , false } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells that t has become a goal of the Domain (added) or no longer is

 void goal_changed( const State & t , bool added ) {
  v_events.push_back( Event{ t , t , Cost() , Cost() , true , added } );
  if( ! added )
   f_target_rose = true;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// edge_changed() and goal_changed() for each event of the EventDomain d

 template< class E >
 void take_events( E & d ) requires EventDomain< E > {
  d.drain_events( [ & ]( const auto & e ) {
   if constexpr( requires { e.from; } )
    edge_changed( e.u , e.v , e.from , e.to );
   else
    goal_changed( e.t , e.added );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds delta to \f$ k_m \f$ and computes h again at the next use of each
 /// state: the target of D* Lite has moved [see HDDStarLite]

 void shift_keys( Cost delta ) {
  f_km = add_cost( f_km , delta );
  ++f_hv;
  f_km_changed = true;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// makes the next solve() start from scratch

 void restart( void ) noexcept { f_fresh = true; }

/*-------------------------------- QUERIES ---------------------------------*/
 /// the most underconsistent expansions of one state in the last episode,
 /// at most 1 with Freeze::all [see theorem 3 of the page of theory]

 [[nodiscard]] std::uint32_t max_raises( void ) const {
  std::uint32_t m = 0;
  for( const auto & w : v_w )
   m = std::max( m , w->f_raises );
  return( m );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Key = LPAKey< Cost >;

 /// a state of a worker; g is read by the others, the rest only by the
 /// owner during a phase, the parents by all after the end of the episode
 struct LNode {
  State state;
  std::atomic< Cost > g;                          // read by the others
  Cost rhs;
  Cost h;
  Key q;                                          // the key pushed last
  Cost g0;                                        // g at the episode start
  LNode * par;
  Cost tm;                                        // min( g , rhs ) told
  std::uint32_t key;                              // of the distribution
  std::uint32_t hv;                               // the version of h
  std::atomic< std::uint32_t > raised;            // the episode of a raise
  std::uint32_t xep;                              // of the last expansion
  std::uint32_t gep;                              // g0 is of this episode
  std::uint32_t rep;                              // raises is of it
  std::uint32_t raises;                           // underconsistent ones
  int gid;                                        // the index of a goal
  bool inq;
  bool goal;
  };

 /// a change, of an edge ( u , v ) or of the goal u
 struct Event {
  State u;
  State v;
  Cost from;
  Cost to;
  bool goal;
  bool added;
  };

 /// the kinds of the messages
 enum Kind : std::uint8_t { RELAX , INVAL , WAKE , TARGET };

 /// a message: Relax( v , x , u ) with h( v ), Invalidate( v , u ), Wake,
 /// or the new \f$ \min( g , \mathrm{rhs} ) \f$ of the goal u to t*
 struct Msg {
  State v;
  Cost x;
  Cost h;
  LNode * u;
  std::uint32_t key;
  Kind kind;
  };

 using Mailbox = ADS::Mailbox< Msg , ADS::DefaultChannels< Msg , Executor > >;

 /// the hash of a state pointed to, by the Domain
 struct ByState {
  using is_avalanching = void;
  const D * d;
  std::uint64_t operator()( const State * s ) const {
   return( d->hash( *s ) );
   }
  };

 /// the equality of the states pointed to
 struct SameState {
  bool operator()( const State * a , const State * b ) const {
   return( *a == *b );
   }
  };

 /// the nodes of a worker, which the others look up while it inserts,
 /// emptied in O( 1 ) at a restart
 using Table = ADS::ConcurrentFlatTable< const State * , LNode * , ByState ,
                                         SameState , ADS::EpochClear<> >;

 /// what the barrier does next
 enum class Stage { events , phase1 , phase2 , restart , done };

 /// the step at a barrier, run by one worker while the others wait
 struct Step {
  HDLPAStar * a;
  void operator()( void ) noexcept { a->at_barrier(); }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its nodes, its queue, its goals and, on the owner of t*, the
 // mirror of the goals

 struct alignas( cache_line ) Worker {

  Worker( HDLPAStar & a , unsigned i )
   : f_a( a ) , f_id( i ) , f_dist( a.f_dist ) ,
    f_table( 0 , ByState{ & a.f_domain } ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  // the counters of a new episode
  void begin_episode( void ) {
   f_st = SearchStats();
   v_ev.clear();
   f_up = f_down = f_rexp = 0;
   f_raises = 0;
   std::fill( std::begin( f_hist ) , std::end( f_hist ) , 0 );
   }

  // the states whose g changed in the episode
  void end_episode( void ) {
   for( LNode * n : v_touched )
    if( n->g.load( std::memory_order_relaxed ) != n->g0 )
     ++f_st.delta;
   v_touched.clear();
   }

  // the episode, in the stages that the barriers choose
  void episode( void ) {
   apply_events();
   for( ; ; ) {
    f_a.f_sync->arrive_and_wait();
    switch( f_a.f_stage ) {
     case( Stage::restart ):
      f_phase = 2;
      if( f_id == f_a.f_s0_owner ) {
       LNode * s = node( f_a.f_s0 , Cost( f_a.f_domain.heuristic( f_a.f_s0 ) )
                         , distribution_key< D >( f_dist , f_a.f_s0 ) );
       s->rhs = Cost( 0 );
       requeue( s );
       }
      run_phase( 2 , false );
      break;
     case( Stage::phase1 ):
      run_phase( 1 , true );
      break;
     case( Stage::phase2 ):
      f_phase = 2;
      for( LNode * n : v_frozen )
       requeue( n );
      v_frozen.clear();
      run_phase( 2 , true );
      break;
     default:
      return;
     }
    }
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  // the changes of the edges into the states of the worker and, on the
  // owner of t*, of the goals [see line 21]

  void apply_events( void ) {
   if( f_a.f_fresh )
    return;
   for( const Event & e : v_ev ) {
    if( e.goal ) {
     goal_of_owner( e.u , e.added );
     continue;
     }
    if( e.to < e.from ) {
     if constexpr( cost_changes_of< D >() != CostChanges::increase_only ) {
      LNode * un = nullptr;
      const Cost gu = read_g( e.u , un );
      if( gu == infinite_cost< Cost >() )
       continue;
      LNode * vn = node( e.v , Cost( f_a.f_domain.heuristic( e.v ) ) ,
                         distribution_key< D >( f_dist , e.v ) );
      const Cost x = add_cost( gu , e.to );
      if( ( e.v != f_a.f_s0 ) && ( x < vn->rhs ) ) {
       vn->rhs = x;
       vn->par = un;
       requeue( vn );
       ++f_down;
       count_head( vn );
       }
      }
     }
    else {
     LNode * vn = find( e.v );
     if( ! vn || ! vn->par )
      continue;
     LNode * un = nullptr;
     static_cast< void >( read_g( e.u , un ) );
     if( un && ( vn->par == un ) ) {
      rescan( vn );
      ++f_up;
      count_head( vn );
      }
     }
    }
   }

  // the head n of a change that matters, in the histogram by its place
  // between the start and the goal [see learn()]
  void count_head( const LNode * n ) {
   ++f_hist[ HDLPAStar::bucket( std::min( n->g.load(
                                std::memory_order_relaxed ) , n->rhs ) ,
                                f_a.f_B.load( std::memory_order_relaxed ) ) ];
   }

  // the goal t added or removed, on the owner of t
  void goal_of_owner( const State & t , bool added ) {
   LNode * n = added ? node( t , Cost( f_a.f_domain.heuristic( t ) ) ,
                             distribution_key< D >( f_dist , t ) )
                     : find( t );
   if( ! n || ( n->goal == added ) )
    return;
   if( added && ( n->gid < 0 ) ) {
    n->gid = int( v_goals.size() );
    v_goals.push_back( n );
    }
   n->goal = added;
   if( ! added )
    erase_heap_key( f_gu , n->gid );
   requeue( n );
   }


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  // forgets the states, whose nodes become free, the queue, the goals
  // and the mirror; at a barrier only
  void clear( void ) {
   f_table.clear();
   v_free.insert( v_free.end() , v_nodes.begin() , v_nodes.end() );
   v_nodes.clear();
   v_goals.clear();
   f_u.clear();
   f_gu.reset( f_gu.ids() );
   v_frozen.clear();
   v_touched.clear();
   if( f_id == f_a.f_t_owner ) {
    f_mu.reset( f_mu.ids() );
    f_mid.clear();
    v_mnode.clear();
    }
   }

  // the node of s in the table of the worker, nullptr if absent
  LNode * find( const State & s ) {
   const auto r = f_table.find( & s , f_a.f_domain.hash( s ) );
   return( r ? *r : nullptr );
   }

  // the node of s of the worker, created with g = rhs = ∞, the heuristic h
  // and the key k of the distribution, from the free ones if any
  LNode * node( const State & s , Cost h , std::uint32_t k ) {
   const std::uint64_t hs = f_a.f_domain.hash( s );
   if( const auto r = f_table.find( & s , hs ) )
    return( *r );
   LNode * n;
   if( v_free.empty() )
    n = f_pool.allocate();
   else {
    n = v_free.back();
    v_free.pop_back();
    }
   v_nodes.push_back( n );
   n->state = s;
   n->gid = -1;
   n->g.store( infinite_cost< Cost >() , std::memory_order_relaxed );
   n->rhs = infinite_cost< Cost >();
   n->h = h;
   n->hv = f_a.f_hv;
   n->par = nullptr;
   n->key = k;
   n->raised.store( 0 , std::memory_order_relaxed );
   n->tm = infinite_cost< Cost >();
   n->xep = n->gep = n->rep = n->raises = 0;
   n->inq = false;
   n->goal = bool( f_a.f_domain.is_goal( s ) );
   if( n->goal && ( n->gid < 0 ) ) {
    n->gid = int( v_goals.size() );
    v_goals.push_back( n );
    }
   static_cast< void >( f_table.findOrInsert( & n->state , hs , n ) );
   return( n );
   }

  // g of s, read in the table of its owner, ∞ if absent; n its node
  Cost read_g( const State & s , LNode * & n ) {
   const unsigned o = f_dist.owner( s , f_a.f_p );
   if( o == f_id ) {
    n = find( s );
    return( n ? n->g.load( std::memory_order_relaxed )
              : infinite_cost< Cost >() );
    }
   ++f_st.remote_reads;
   const auto r = f_a.v_w[ o ]->f_table.find( & s , f_a.f_domain.hash( s ) );
   n = nullptr;
   if( ! r )
    return( infinite_cost< Cost >() );
   n = *r;
   return( n->g.load( std::memory_order_acquire ) );
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  // h of n, computed again if the heuristic changed since
  Cost heur( LNode * n ) {
   if( n->hv != f_a.f_hv ) {
    n->h = Cost( f_a.f_domain.heuristic( n->state ) );
    n->hv = f_a.f_hv;
    }
   return( n->h );
   }

  // the key of n
  Key key( LNode * n ) {
   const Cost m = std::min( n->g.load( std::memory_order_relaxed ) ,
                            n->rhs );
   if( m == infinite_cost< Cost >() )
    return( Key{ m , m } );
   return( Key{ add_cost( Cost( m + heur( n ) ) , f_a.f_km ) , m } );
   }

  // g( n ) ← x, by an atomic store
  void set_g( LNode * n , Cost x ) {
   if( n->gep != f_a.f_episode ) {
    n->gep = f_a.f_episode;
    n->g0 = n->g.load( std::memory_order_relaxed );
    v_touched.push_back( n );
    }
   n->g.store( x , std::memory_order_release );
   }

  // n in U if inconsistent, with its key, and a goal in the heap of the
  // inconsistent goals
  void requeue( LNode * n ) {
   const Cost g = n->g.load( std::memory_order_relaxed );
   const bool inc = g != n->rhs;
   if( inc ) {
    const Key k = key( n );
    if( ! n->inq || ! ( k == n->q ) ) {
     f_u.push( n , k );
     n->q = k;
     }
    }
   n->inq = inc;
   if( n->goal ) {
    if( inc )
     set_heap_key( f_gu , n->gid , std::min( g , n->rhs ) );
    else
     erase_heap_key( f_gu , n->gid );
    }
   const Cost m = n->goal ? std::min( g , n->rhs ) : infinite_cost< Cost >();
   if( m != n->tm ) {
    n->tm = m;
    tell_target( n , m );
    }
   }

  // the new min( g , rhs ) m of the goal n to the owner of t*, at once if
  // it is the worker, at the next barrier if no phase runs
  void tell_target( LNode * n , Cost m ) {
   if( f_phase == 0 )
    v_tpend.emplace_back( n , m );
   else
    if( f_id == f_a.f_t_owner )
     target_set( n , m , true );
    else {
     ++f_st.relaxes;
     send( f_a.f_t_owner , Msg{ n->state , m , Cost() , n , 0 , TARGET } );
     }
   }

  // whether n is frozen in phase 1: overconsistent, raised in this
  // episode and, with Freeze::raised_parent, its parent raised too; never
  // with one worker, which expands in the order of LPA*, whose episodes
  // end with no freeze
  bool frozen( const LNode * n ) const {
   if( f_a.f_p == 1 )
    return( false );
   const std::uint32_t e = f_a.f_episode;
   if( ( n->raised.load( std::memory_order_relaxed ) != e ) ||
       ! ( n->g.load( std::memory_order_relaxed ) > n->rhs ) )
    return( false );
   if( f_a.f_options.freeze == IncrementalOptions::Freeze::all )
    return( true );
   return( n->par &&
           ( n->par->raised.load( std::memory_order_relaxed ) == e ) );
   }

  // the least eligible state of U, removed; nullptr if none [see line 5]
  LNode * next( Cost b , int phase ) {
   const Key kt{ add_cost( b , f_a.f_km ) , b };
   while( ! f_u.empty() ) {
    const auto & e = f_u.top();
    LNode * n = e.x;
    if( ! n->inq || ! ( e.k == n->q ) ) {
     f_u.pop();
     continue;
     }
    const Key k = key( n );
    if( n->q < k ) {
     f_u.pop();
     f_u.push( n , k );
     n->q = k;
     continue;
     }
    if( ( phase == 1 ) && frozen( n ) ) {
     f_u.pop();
     put_aside( n );
     continue;
     }
    if( ! ( k < kt ) )
     break;
    f_u.pop();
    take( n );
    return( n );
    }
   while( ! f_gu.empty() && ( f_gu.top().key <= b ) ) {
    LNode * n = v_goals[ std::size_t( f_gu.top().id ) ];
    if( ( phase == 1 ) && frozen( n ) ) {
     put_aside( n );
     continue;
     }
    take( n );
    return( n );
    }
   return( nullptr );
   }

  // n out of U and of the heap of the goals
  void take( LNode * n ) {
   n->inq = false;
   if( n->goal )
    erase_heap_key( f_gu , n->gid );
   }

  // n frozen until phase 2
  void put_aside( LNode * n ) {
   take( n );
   v_frozen.push_back( n );
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  // RESCAN [see line 19]
  void rescan( LNode * v ) {
   ++f_st.rescans;
   v->rhs = infinite_cost< Cost >();
   v->par = nullptr;
   if( v->state == f_a.f_s0 )
    v->rhs = Cost( 0 );
   else
    f_a.f_domain.for_each_predecessor( v->state , [ & ]( const State & w ,
                                                         Cost c ) {
     LNode * wn = nullptr;
     const Cost x = add_cost( read_g( w , wn ) , c );
     if( x < v->rhs ) {
      v->rhs = x;
      v->par = wn;
      }
     } );
   requeue( v );
   }

  // the message m, on its owner [see lines 14 to 17]
  void receive( const Msg & m ) {
   switch( m.kind ) {
    case( RELAX ): {
     LNode * n = node( m.v , m.h , m.key );
     if( ( m.v != f_a.f_s0 ) && ( m.x < n->rhs ) ) {
      n->rhs = m.x;
      n->par = m.u;
      requeue( n );
      }
     break;
     }
    case( INVAL ): {
     LNode * n = find( m.v );
     if( n && ( n->par == m.u ) )
      rescan( n );
     break;
     }
    case( TARGET ):
     target_set( m.u , m.x , true );
     break;
    default:
     break;
    }
   }

  // the mirror of the goal t set to x, on the owner of t*, and B with it;
  // a rise of B woken up in phase 1 if wake [see line 15]
  void target_set( LNode * t , Cost x , bool wake ) {
   auto it = f_mid.find( t );
   int id;
   if( it == f_mid.end() ) {
    if( x == infinite_cost< Cost >() )
     return;
    id = int( v_mnode.size() );
    f_mid.emplace( t , id );
    v_mnode.push_back( t );
    }
   else
    id = it->second;
   if( x == infinite_cost< Cost >() )
    erase_heap_key( f_mu , id );
   else
    set_heap_key( f_mu , id , x );
   const Cost b = f_mu.empty() ? infinite_cost< Cost >() : f_mu.top().key;
   const Cost old = f_a.f_B.load( std::memory_order_relaxed );
   if( b == old )
    return;
   f_a.f_B.store( b , std::memory_order_relaxed );
   if( wake && ( b > old ) && ( f_phase == 1 ) )
    for( unsigned j = 0 ; j < f_a.f_p ; ++j )
     if( j != f_id ) {
      f_a.f_mail->send( f_id , j , Msg{ t->state , Cost() , Cost() ,
                                        nullptr , 0 , WAKE } );
      ++f_st.wakes;
      ++f_st.sent;
      }
   }

  // the message m sent, or applied at once if it is for the worker
  void send( unsigned to , const Msg & m ) {
   if( to == f_id ) {
    receive( m );
    return;
    }
   f_a.f_mail->send( f_id , to , m );
   ++f_st.sent;
   }

  // EXPAND [see lines 7 to 13]
  void expand( LNode * u ) {
   ++f_st.expanded;
   if( u->xep == f_a.f_episode )
    ++f_st.reexpanded;
   u->xep = f_a.f_episode;
   const unsigned p = f_a.f_p;
   const auto owner = successor_owners< D >( f_dist , u->state , u->key ,
                                             p );
   const Cost g = u->g.load( std::memory_order_relaxed );
   if( g > u->rhs ) {                            // overconsistent
    const Cost gu = u->rhs;
    set_g( u , gu );
    f_a.f_domain.for_each_successor( u->state , heur( u ) ,
                                     [ & ]( const State & t , Cost c ,
                                            Cost ht ) {
     ++f_st.generated;
     if( t == f_a.f_s0 )
      return;
     const auto [ dst , k ] = owner( t );
     ++f_st.relaxes;
     send( dst , Msg{ t , add_cost( gu , c ) , ht , u , k , RELAX } );
     } );
    }
   else {                                         // underconsistent
    ++f_st.under;
    if( u->rep != f_a.f_episode ) {
     u->rep = f_a.f_episode;
     u->raises = 0;
     }
    f_raises = std::max( f_raises , ++u->raises );
    set_g( u , infinite_cost< Cost >() );
    u->raised.store( f_a.f_episode , std::memory_order_relaxed );
    f_a.f_domain.for_each_successor( u->state , heur( u ) ,
                                     [ & ]( const State & t , Cost ,
                                            Cost ht ) {
     ++f_st.generated;
     if( t == f_a.f_s0 )
      return;
     const auto [ dst , k ] = owner( t );
     ++f_st.invalidates;
     send( dst , Msg{ t , Cost() , ht , u , k , INVAL } );
     } );
    rescan( u );
    }
   requeue( u );
   }

  // WORKER [see lines 1 to 6], within the budget of a repair if budget
  void run_phase( int phase , bool budget ) {
   f_phase = phase;
   Mailbox & mail = *f_a.f_mail;
   bool active = true;
   unsigned since = 0;
   std::uint64_t spent = 0;
   ADS::Backoff backoff;
   for( ; ; ) {
    v_in.clear();
    if( mail.poll( f_id , [ & ]( Msg && m ) { v_in.push_back( m ); } ) ) {
     active = true;
     for( const Msg & m : v_in )
      receive( m );
     }
    if( active ) {
     const Cost b = f_a.f_B.load( std::memory_order_relaxed );
     if( LNode * u = next( b , phase ) ) {
      backoff.reset();
      expand( u );
      ++spent;
      if( ! budget )
       ++f_rexp;
      if( ++since >= f_a.f_options.flush_period ) {
       since = 0;
       mail.flush( f_id );
       if( budget ) {
        const std::uint64_t all = f_a.f_spent.fetch_add(
         spent , std::memory_order_relaxed ) + spent;
        spent = 0;
        if( double( all ) > f_a.f_budget ) {
         f_a.f_aborted.store( true , std::memory_order_relaxed );
         mail.stop();
         }
        }
       }
      continue;
      }
     }
    if( mail.quiescent( f_id ) )
     break;
    active = false;
    backoff.pause();
    }
   if( budget )
    f_a.f_spent.fetch_add( spent , std::memory_order_relaxed );
   f_phase = 0;
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  HDLPAStar & f_a;                              // the search
  unsigned f_id;                                // the index of the worker
  Dist f_dist;                                  // its distribution
  NodePool< LNode > f_pool;                     // its nodes
  Table f_table;                                // their table
  Queue< LNode * , Cost > f_u;                  // the queue U
  ADS::IndexedDaryHeap< Cost > f_gu;            // its inconsistent goals
  std::vector< LNode * > v_goals;               // its goals, by gid
  std::vector< LNode * > v_frozen;              // put aside in phase 1
  std::vector< LNode * > v_touched;             // g0 of the episode set
  std::vector< LNode * > v_nodes;               // in the table
  std::vector< LNode * > v_free;                // not in the table
  std::vector< Event > v_ev;                    // the changes it applies
  std::vector< Msg > v_in;                      // the last inbox
  std::vector< std::pair< LNode * , Cost > > v_tpend;  // for t*, at the
                                                // barrier
  ADS::IndexedDaryHeap< Cost > f_mu;            // on the owner of t*: the
  std::unordered_map< const LNode * , int > f_mid;   // mirror of the goals
  std::vector< LNode * > v_mnode;               // and their nodes
  SearchStats f_st;                             // its counters
  std::uint64_t f_up = 0;                       // increases of a parent
  std::uint64_t f_down = 0;                     // decreases that lower
  std::uint64_t f_rexp = 0;                     // expansions of a restart
  std::uint32_t f_raises = 0;                   // most of one state
  std::uint64_t f_hist[ 8 ] = {};               // heads by place
  int f_phase = 0;                              // the phase running
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the bucket of a head of value m (its least of g and rhs) between the
 // start and B, in 8 parts
 static std::size_t bucket( Cost m , Cost b ) {
  if( ( b == infinite_cost< Cost >() ) || ( m == infinite_cost< Cost >() ) ||
      ! ( b > Cost( 0 ) ) )
   return( 7 );
  return( std::min( std::size_t( 7 ) ,
                    std::size_t( 8.0 * double( m ) / double( b ) ) ) );
  }

 // what the barrier does next, decided by one worker while the others wait
 void at_barrier( void ) noexcept {
  const bool aborted = f_aborted.exchange( false ,
                                           std::memory_order_relaxed );
  if( ( f_stage != Stage::events ) && ! aborted && f_options.audit &&
      f_audit_error.empty() )
   audit( f_stage == Stage::phase1 );
  Stage next = Stage::done;
  switch( f_stage ) {
   case( Stage::events ): {
    for( auto & w : v_w ) {
     for( const auto & [ n , m ] : w->v_tpend )
      v_w[ f_t_owner ]->target_set( n , m , false );
     w->v_tpend.clear();
     }
    f_restarted = false;
    f_phases = 0;
    f_batches = batches();
    f_spent.store( 0 , std::memory_order_relaxed );
    std::uint64_t up = 0 , down = 0;
    for( auto & w : v_w ) {
     up += w->f_up;
     down += w->f_down;
     }
    const Cost b = f_B.load( std::memory_order_relaxed );
    const bool rose = f_target_rose || ( f_B_end < b );
    f_target_rose = false;
    if( f_fresh || decide_restart( up + down ) )
     next = Stage::restart;
    else
     if( ( up > 0 ) || rose )
      next = Stage::phase1;
     else
      next = Stage::phase2;
    f_budget = ( f_S > 0 ) && ( f_options.policy ==
                                IncrementalOptions::Policy::adaptive )
               ? f_options.beta * f_S : HUGE_VAL;
    break;
    }
   case( Stage::phase1 ):
    next = aborted ? Stage::restart : Stage::phase2;
    break;
   case( Stage::phase2 ):
    next = aborted ? Stage::restart : Stage::done;
    break;
   default:
    next = Stage::done;
   }
  if( ! f_audit_error.empty() )
   next = Stage::done;
  if( f_mail->stopped() )
   static_cast< void >( f_mail->leftovers( []( Msg && ) {} ) );
  if( next == Stage::restart ) {
   for( auto & w : v_w )
    w->clear();
   f_fresh = false;
   f_restarted = true;
   f_B.store( infinite_cost< Cost >() , std::memory_order_relaxed );
   }
  if( next != Stage::done ) {
   if( f_new_mail )
    f_new_mail = false;                      // all active, as built
   else
    f_mail->rearm();
   ++f_phases;
   }
  else
   f_B_end = f_B.load( std::memory_order_relaxed );
  f_stage = next;
  }

 // the batches pushed so far
 std::uint64_t batches( void ) const {
  std::uint64_t n = 0;
  for( unsigned i = 0 ; i < f_p ; ++i )
   n += f_mail->batches( i );
  return( n );
  }

 // whether the changes that matter, n of them, and the prediction of the
 // repair ask for a restart [see IncrementalOptions]
 bool decide_restart( std::uint64_t n ) {
  using P = IncrementalOptions::Policy;
  if( f_options.policy == P::restart )
   return( true );
  if( ( f_options.policy == P::repair ) || ( f_S <= 0 ) )
   return( false );
  if( double( n ) > f_options.rho * f_S )
   return( true );
  double r = 0;
  for( std::size_t b = 0 ; b < 8 ; ++b ) {
   std::uint64_t c = 0;
   for( auto & w : v_w )
    c += w->f_hist[ b ];
   v_heads[ b ] = c;
   r += double( c ) * f_r[ b ];
   }
  return( r > f_options.gamma * f_S );
  }

 // the means of the restarts and of the repairs, after an episode of e
 // expansions, r of them by a restart
 void learn( std::uint64_t e , std::uint64_t r ) {
  const double w = f_options.mean_weight;
  if( f_restarted ) {
   f_S = f_S > 0 ? ( 1 - w ) * f_S + w * double( r ) : double( r );
   return;
   }
  std::uint64_t h = 0;
  for( std::size_t b = 0 ; b < 8 ; ++b ) {
   std::uint64_t c = 0;
   for( auto & x : v_w )
    c += x->f_hist[ b ];
   v_heads[ b ] = c;
   h += c;
   }
  if( h == 0 )
   return;
  const double per = double( e ) / double( h );
  for( std::size_t b = 0 ; b < 8 ; ++b )
   if( v_heads[ b ] > 0 )
    f_r[ b ] = f_r[ b ] > 0 ? ( 1 - w ) * f_r[ b ] + w * per : per;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the invariants at a quiescent point: rhs exact with its parent, every
 // state reached from a state of finite g in a table, the mirror equal to
 // the least of g and rhs of the goals and its least B to the least g, U
 // the inconsistent states, and none of them below the key of t* but the
 // frozen ones after phase 1; the first failure in f_audit_error
 void audit( bool after_phase1 ) {
  auto fail = [ & ]( const std::string & s ) {
   if( f_audit_error.empty() )
    f_audit_error = "audit: " + s;
   };
  Worker & tw = *v_w[ f_t_owner ];
  const Cost b = f_B.load( std::memory_order_relaxed );
  const Key kt{ add_cost( b , f_km ) , b };
  Cost least = infinite_cost< Cost >() , least_m = least;
  std::vector< std::uint64_t > reads;
  for( auto & w : v_w )
   reads.push_back( w->f_st.remote_reads );
  for( auto & w : v_w )
   for( LNode * v : w->v_nodes ) {
    const Cost g = v->g.load( std::memory_order_relaxed );
    Cost best = infinite_cost< Cost >();
    bool par_ok = false;
    if( v->state == f_s0 ) {
     best = Cost( 0 );
     par_ok = true;
     }
    else
     f_domain.for_each_predecessor( v->state , [ & ]( const State & s ,
                                                      Cost c ) {
      LNode * n = nullptr;
      const Cost x = add_cost( w->read_g( s , n ) , c );
      best = std::min( best , x );
      if( n && ( n == v->par ) && ( x == v->rhs ) )
       par_ok = true;
      } );
    if( best != v->rhs )
     fail( "rhs not exact" );
    if( ( best != infinite_cost< Cost >() ) && ! par_ok )
     fail( "parent not exact" );
    if( g != infinite_cost< Cost >() )
     f_domain.for_each_successor( v->state , w->heur( v ) ,
                                  [ & ]( const State & t , Cost , Cost ) {
      LNode * n = nullptr;
      static_cast< void >( w->read_g( t , n ) );
      if( ! n && ! ( t == f_s0 ) )
       fail( "a successor of finite rhs not in a table" );
      } );
    const Cost m = std::min( g , v->rhs );
    if( v->goal && ( m != infinite_cost< Cost >() ) ) {
     least = std::min( least , g );
     least_m = std::min( least_m , m );
     const auto it = tw.f_mid.find( v );
     if( ( it == tw.f_mid.end() ) || ! tw.f_mu.contains( it->second ) ||
         ( tw.f_mu.key( it->second ) != m ) )
      fail( "the mirror of a goal not its least of g and rhs" );
     }
    if( g == v->rhs )
     continue;
    const bool fz = after_phase1 && w->frozen( v );
    if( ! v->inq && ! fz )
     fail( "an inconsistent state not in U" );
    if( fz )
     continue;
    const Key k = w->key( v );
    if( ( k < kt ) || ( v->goal && ! ( kt < k ) ) )
     fail( "an inconsistent state below the key of t*" );
    }
  if( least_m != b )
   fail( "B not the least of g and rhs of the goals" );
  if( ! after_phase1 && ( least != b ) )
   fail( "B not the least g of the goals at the end" );
  for( unsigned i = 0 ; i < f_p ; ++i )
   v_w[ i ]->f_st.remote_reads = reads[ i ];
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                    ///< the work distribution
 HDLPAOptions f_options;                               ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr
 unsigned f_p = 1;                                     ///< the workers
 State f_s0{};                                         ///< the start
 unsigned f_s0_owner = 0;                              ///< its owner
 unsigned f_t_owner = 0;                               ///< the owner of t*
 std::vector< std::unique_ptr< Worker > > v_w;         ///< the workers
 std::vector< Event > v_events;                        ///< changes told
 std::unique_ptr< Mailbox > f_mail;                    ///< the messages
 bool f_new_mail = true;                               ///< not used yet
 std::barrier< Step > * f_sync = nullptr;              ///< of the episode
 Stage f_stage = Stage::events;                        ///< at the barrier
 alignas( cache_line ) std::atomic< Cost > f_B{
                                    infinite_cost< Cost >() }; ///< g( t* )
 alignas( cache_line ) std::atomic< std::uint64_t > f_spent{ 0 };
                                              ///< expansions of the repair
 std::atomic< bool > f_aborted{ false };              ///< repair stopped
 double f_budget = HUGE_VAL;                          ///< of the repair
 Cost f_B_end = infinite_cost< Cost >();              ///< B at the last end
 Cost f_km = Cost( 0 );                               ///< k_m
 std::uint32_t f_hv = 0;                              ///< version of h
 std::uint32_t f_episode = 0;                         ///< of the solve()
 bool f_fresh = true;                                 ///< restart next
 bool f_restarted = false;                            ///< in this episode
 bool f_target_rose = false;                          ///< a goal removed
 bool f_km_changed = false;                           ///< k_m shifted
 std::uint64_t f_batches = 0;                         ///< of the episode
 std::uint64_t f_phases = 0;                          ///< of the episode
 double f_S = 0;                              ///< mean of the restarts
 double f_r[ 8 ] = {};                        ///< expansions per head
 std::uint64_t v_heads[ 8 ] = {};             ///< heads by bucket
 std::string f_audit_error;                   ///< the first failure

 };  // end( class( HDLPAStar ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS HDDStarLite ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a parallel D* Lite: HDLPAStar on the reversed Domain with a robot that
/// moves
/** HDLPAStar on the ReverseDomain of D [see LPAStar.h], from the goal of the
 * robot to its state, as DStarLite is LPAStar on it: solve() returns the
 * shortest path from the robot to its goal, the caller tells the changes of
 * the edges of D, and move_to( r ) adds the estimate of the move to
 * \f$ k_m \f$ at the next barrier and makes r the target. The distribution
 * is one of D, which gives the owner of a state of the reversed Domain as
 * well. */

template< class D , class Dist = StateHashDistribution< D > ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          template< class , class > class Queue = LPADefaultQueue >
 requires PredecessorDomain< D > && PairwiseDomain< D >
class HDDStarLite {

 public:

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// constructor, taking D and the distribution (which must outlive the
 /// HDDStarLite), the goal and the state of the robot

 HDDStarLite( D & d , const Dist & dist , const State & goal ,
              const State & robot , const HDLPAOptions & options = {} )
  : f_r( d , goal , robot ) , f_s( f_r , dist , options ) {}

 /// sets the Executor that runs the workers [see HDLPAStar]

 void set_executor( Executor & ex ) noexcept { f_s.set_executor( ex ); }

 /// the shortest path from the robot to its goal, path.front() being the
 /// robot

 [[nodiscard]] SolutionT solve( void ) {
  SolutionT s = f_s.solve();
  std::reverse( s.path.begin() , s.path.end() );
  return( s );
  }

 /// tells that the cost of the edge ( u , v ) of D has changed

 void edge_changed( const State & u , const State & v , Cost from , Cost to )
 {
  f_s.edge_changed( v , u , from , to );
  }

 /// edge_changed() for each event of D

 void take_events( void ) requires EventDomain< D > {
  f_s.take_events( f_r );
  }

 /// moves the robot to r

 void move_to( const State & r ) {
  const State a = f_r.robot();
  if( a == r )
   return;
  const Cost delta = Cost( f_r.domain().pairwise_heuristic( a , r ) );
  f_r.set_robot( r );
  f_s.goal_changed( a , false );
  f_s.goal_changed( r , true );
  f_s.shift_keys( delta );
  }

 /// makes the next solve() start from scratch

 void restart( void ) noexcept { f_s.restart(); }

 /// the state of the robot

 [[nodiscard]] const State & robot( void ) const { return( f_r.robot() ); }

 private:

 ReverseDomain< D > f_r;                                ///< the space
 HDLPAStar< ReverseDomain< D > , Dist , Executor , Queue > f_s;
                                                        ///< the search

 };  // end( class( HDDStarLite ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_HDLPAStar */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File HDLPAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
