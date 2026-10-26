/*--------------------------------------------------------------------------*/
/*----------------------------- File HDAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template HDAStar, Hash Distributed A* on a
 * shared-memory machine, whose workers run on a runtime [see
 * ThreadRuntime.h] and exchange nodes through its channels.
 *
 * The search follows section 3 of Kishimoto, Fukunaga and Botea (AIJ 2013)
 * and section 4.1 of Fukunaga et al. (2017), with the notation of the page
 * "Best-first search and A*" of the documentation; the worker \f$ i \f$ owns
 * the states \f$ v \f$ with \f$ \mathrm{owner}( v ) = i \f$, which only it
 * stores, in its lists \f$ \mathrm{OPEN}_i \f$ and
 * \f$ \mathrm{CLOSED}_i \f$, the incumbent \f$ z \f$ is atomic, and a
 * message \f$ ( v , g' , u ) \f$ carries a state with the cost and the
 * predecessor of the path that reached it:
 *
 * \code{.unparsed}
 *  1  procedure HDA*(s₀, p, owner)
 *  2     z ← ∞;  send (s₀, 0, nil) to owner(s₀)
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        while the search is not over do
 *  5           for each message (v, g', u) received do INSERT(v, g', u)
 *  6           if OPEN_i has a node of f < z then
 *  7              u ← the first node of OPEN_i, removed (stale entries
 *                     dropped);  close u
 *  8              if u is a goal then
 *  9                 z ← min(z, g(u)), recording u if g(u) < z
 * 10              else for each successor v of u do
 * 11                 g' ← g(u) + c(u, v)
 * 12                 if g' + h(v) ≥ z then continue
 * 13                 if owner(v) = i then INSERT(v, g', u)
 * 14                 else send (v, g', u) to owner(v)
 * 15              every k expansions, send the outboxes not yet full
 * 16           else
 * 17              send all the outboxes;  WAIT(i)      ▷ see ThreadMailbox.h
 * 18     return the path to the goal of cost z, if z < ∞
 * 19  procedure INSERT(v, g', u)                       ▷ on the owner i of v
 * 20     if g' + h(v) ≥ z then return
 * 21     if v ∉ CLOSED_i then
 * 22        add v to CLOSED_i, with g(v) ← g' and parent(v) ← u
 * 23     else if g' < g(v) then
 * 24        g(v) ← g';  parent(v) ← u;  v is not closed
 * 25     else return                                         ▷ a duplicate
 * 26     push (g' + h(v), g', v) in OPEN_i
 * \endcode
 *
 * The messages travel in batches of batch nodes, and an outbox is sent when
 * it is full (line 14), every flush_period expansions (line 15) and before the
 * worker waits (line 17); WAIT ends when mail arrives or when the end is
 * detected by the three counters of ThreadMailbox, i.e., when no worker has
 * a node of \f$ f < z \f$ and no message is in transit, which is the
 * condition of the page "Parallel A*" for an optimal \f$ z \f$. The parents
 * cross the workers, and the path is followed only after all the workers
 * have finished.
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

#ifndef __PBeFS_HDAStar
 #define __PBeFS_HDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <memory>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/core/ThreadMailbox.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ParallelOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the parallel searches: those of the sequential search,
/// and the parameters of the parallel ones [see ParallelParams], with
/// batches of 64 nodes flushed every 64 expansions

struct ParallelOptions : public SearchOptions , public ParallelParams {
 };  // end( struct( ParallelOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS HDAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Hash Distributed A* on a shared-memory machine
/** Every worker runs A* on its own OPEN and CLOSED lists, and a global
 * function Dist assigns every state to the worker that owns it: a node
 * generated by a worker is sent to its owner, which checks it against its
 * closed list and possibly adds it to its OPEN list [see section 4.1 of
 * Fukunaga et al., 2017, and Kishimoto, Fukunaga and Botea, AIJ 2013].
 * Nodes are sent asynchronously, packed in batches of batch nodes, and an
 * outbox not yet full is flushed every flush_period expansions and before
 * the worker waits. The goal test is done at the expansion, the cost of
 * the best goal found (the incumbent) prunes every node with f not
 * smaller, and the search ends when no worker has a node with f smaller
 * than the incumbent and no node is travelling, which is detected with
 * the three counters of ThreadMailbox [see ThreadMailbox.h].
 *
 * The node of a state belongs to its owner, which is the only worker that
 * writes it; parent pointers cross the workers, and they are followed only
 * after all the workers have finished. */

template< Domain D , class Dist = StateHashDistribution< D > ,
          class Runtime = ThreadRuntime ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class HDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// a node travelling to its owner
 struct Msg {
  State s;
  Cost g;
  Cost h;
  NodeT * parent;
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain and the distribution (which must
 /// outlive the HDAStar)

 HDAStar( const D & domain , const Dist & dist ,
          const ParallelOptions & options = {} )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) {
  f_options.normalize();
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  Shared sh( p );

  // the root goes to its owner as the first message
  {
   const State s0 = f_domain.initial_state();
   std::vector< Msg > b{ Msg{ s0 , Cost( 0 ) ,
                              Cost( f_domain.heuristic( s0 ) ) , nullptr } };
   sh.mail.seed( f_dist.owner( s0 , p ) , b );
   }

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );

  Runtime::run( p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.sent += x->f_st.sent;
   st.messages += x->f_out.messages();
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

 using Mailbox = ThreadMailbox< Msg , Runtime >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  explicit Shared( unsigned np ) : p( np ) , mail( np ) {}

  unsigned p;                                   // the number of workers
  Mailbox mail;                                 // the messages
  alignas( cache_line ) std::atomic< Cost > incumbent{
                                                   infinite_cost< Cost >() };
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  NodeT * goal = nullptr;                       // the best goal found
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its OPEN and CLOSED lists, its outboxes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const HDAStar & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_dist( a.f_dist ) ,
     f_out( s.mail , i , a.f_options.batch ) {}

  // the main loop
  void run( void ) {
   unsigned since_flush = 0;
   while( ! f_sh.mail.done() ) {
    if( std::size_t n = f_sh.mail.receive( f_id , v_in ,
                                           f_open.empty() ) ) {
     // the hashes and the prefetches of the whole batch first, so that
     // the cache misses of the lookups overlap
     v_hash.resize( n );
     for( std::size_t k = 0 ; k < n ; ++k ) {
      v_hash[ k ] = f_a.f_domain.hash( v_in[ k ].s );
      f_table.prefetch( v_hash[ k ] );
      }
     for( std::size_t k = 0 ; k < n ; ++k )
      insert( v_in[ k ] , v_hash[ k ] );
     f_sh.mail.processed( n );
     }

    if( NodeT * n = next() ) {
     expand( n );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      f_out.flush( false );
      }
     if( f_a.f_options.node_limit &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      f_sh.limit_hit.store( true );
      f_sh.mail.stop();
      }
     continue;
     }

    f_out.flush( true );     // nothing to expand: send everything, then
    f_sh.mail.wait( f_id );  // wait for mail or for the end
    }
   }

  // the best node of OPEN with f < incumbent, removed; nullptr if none
  NodeT * next( void ) {
   while( ! f_open.empty() ) {
    const OpenEntry< NodeT > e = f_open.top();
    if( e.stale() ) {
     f_open.pop();
     continue;
     }
    if( e.node->g + e.node->h >= f_sh.incumbent.load(
                                                std::memory_order_relaxed ) )
     return( nullptr );
    f_open.pop();
    return( e.node );
    }
   return( nullptr );
   }

  // expands n, or records it as the incumbent if it is a goal
  void expand( NodeT * n ) {
   n->closed = true;
   ++f_st.expanded;
   if( n->expanded )
    ++f_st.reexpanded;
   n->expanded = true;

   if( f_a.f_domain.is_goal( n->state ) ) {
    std::lock_guard< std::mutex > lock( f_sh.goal_mutex );
    if( n->g < f_sh.incumbent.load() ) {
     f_sh.incumbent.store( n->g );
     f_sh.goal = n;
     }
    return;
    }

   const unsigned p = f_sh.p;
   f_a.f_domain.for_each_successor( n->state , n->h ,
                                    [ & ]( const State & t , Cost c ,
                                           Cost ht ) {
    ++f_st.generated;
    const Cost g = n->g + c;
    if( g + ht >= f_sh.incumbent.load( std::memory_order_relaxed ) )
     return;
    const unsigned dst = f_dist.owner( t , p );
    if( dst == f_id ) {
     insert( Msg{ t , g , ht , n } , f_a.f_domain.hash( t ) );
     return;
     }
    f_out.push( dst , Msg{ t , g , ht , n } );
    ++f_st.sent;
    } );
   }

  // adds the node of m, whose state has hash h, to OPEN, unless it is a
  // duplicate
  void insert( const Msg & m , std::uint64_t h ) {
   if( m.g + m.h >= f_sh.incumbent.load( std::memory_order_relaxed ) )
    return;
   auto [ x , fresh ] = f_table.find_or_insert( m.s , h , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ m.s , m.g , m.h , m.parent , false , false };
    return( y );
    } );
   if( fresh )
    f_open.push( x );
   else
    if( m.g < x->g ) {
     x->g = m.g;
     x->parent = m.parent;
     x->closed = false;
     f_open.push( x );
     }
    else
     ++f_st.duplicates;
   }

  const HDAStar & f_a;                          // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  Dist f_dist;                                  // its distribution
  NodePool< NodeT > f_pool;                     // its nodes
  NodeTable< NodeT > f_table;                   // its CLOSED list
  Open f_open;                                  // its OPEN list
  typename Mailbox::Outboxes f_out;             // its outboxes
  std::vector< Msg > v_in;                      // the last inbox received
  std::vector< std::uint64_t > v_hash;          // the hashes of v_in
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                    ///< the work distribution
 ParallelOptions f_options;                            ///< the options

 };  // end( class( HDAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_HDAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File HDAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
