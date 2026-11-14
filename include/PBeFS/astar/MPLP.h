/*--------------------------------------------------------------------------*/
/*------------------------------ File MPLP.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template MPLP, the Massively Parallelized Lazy
 * Planning of Mukherjee, Aine and Likhachev (RA-L 2022).
 *
 * The search is lazy: the successors of a state are generated without
 * evaluating their edges [see LazyEdgeDomain], each edge taking a lower
 * bound \f$ \underline{c} \f$ on its cost, and the edges are evaluated by
 * the other workers while the search goes on. With
 * \f$ c( e ) \f$ the cost of the edge \f$ e \f$ if evaluated (\f$ \infty \f$
 * if it does not exist) and \f$ \underline{c}( e ) \f$ otherwise:
 *
 * \code{.unparsed}
 *  1  procedure MPLP(s₀, p)                        ▷ on the worker 0
 *  2     z ← −∞;  Π ← ∅
 *  3     repeat
 *  4        π ← weighted A* from s₀ on the costs c, each edge read at the
 *  5            expansion of its state, the new edges queued for evaluation
 *  6        if there is no π then return no path
 *  7        z ← max(z, c(π));  if π ∉ Π then add π to Π, and its edges
 *  8            not evaluated to the front of the queue
 *  9        wait until an edge evaluated has changed its cost
 * 10     until a path is accepted
 * 11  procedure EVALUATE                           ▷ on the workers 1, …
 * 12     loop: take the first edge e of the queue, evaluate it
 * 13        for each π ∈ Π through e do
 * 14           if e does not exist then drop π
 * 15           else if all the edges of π are evaluated and c(π) ≤ z then
 * 16              accept π
 * \endcode
 *
 * The edges discovered go to the back of the queue, those of a path found to
 * its front, which is the priority 1 and 2 of the paper. Since
 * \f$ \underline{c} \leq c \f$, every path of the graph is a path of the
 * graph searched at line 4 and is not longer there; hence, if no path is
 * found, none exists (line 6), and each \f$ \pi \f$ found costs at most
 * \f$ w \, C^* \f$, the weighted A* without reopening being within
 * \f$ w \f$ of the optimum of a graph whose edges are not longer, and the
 * heuristic being consistent with \f$ \underline{c} \f$ (Lemma 3 of the
 * paper). So \f$ z \leq w \, C^* \f$, and a path accepted at line 16 costs
 * at most \f$ w \, C^* \f$ (Theorem 5); the search is complete on a finite
 * graph, since every edge is eventually evaluated and the paths without
 * repetition are finitely many (Theorem 2).
 *
 * The paper runs the search, the check of the paths and the hand-out of
 * the edges on three threads; here the check of a path is done by the
 * worker that evaluates its last edge, and the workers take the edges from
 * the queue themselves, so that \f$ p \f$ workers evaluate on \f$ p - 1 \f$
 * of them instead of \f$ p - 3 \f$. The search does not repeat while no
 * evaluation has changed a cost (line 9), since it would find the same path;
 * the paper repeats it at once. With one worker the search evaluates the
 * edges of each path it finds before searching again, which is the lazy
 * shortest path of Dellin and Srinivasa (ICAPS 2016).
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

#ifndef __PBeFS_MPLP
 #define __PBeFS_MPLP
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cmath>

#include <compare>

#include <condition_variable>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <memory>

#include <mutex>

#include <type_traits>

#include <unordered_map>

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
/*---------------------------- CLASS MPLPOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of MPLP: those of ParallelOptions, whose weight
/// \f$ w \geq 1 \f$ is that of the heuristic of the lazy searches, the
/// cost being within a factor \f$ w \f$ of the optimum

struct MPLPOptions : public ParallelOptions {

 };  // end( struct( MPLPOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS MPLP -------------------------------*/
/*--------------------------------------------------------------------------*/
/// Massively Parallelized Lazy Planning
/** The worker 0 repeats the lazy weighted A* and the other workers
 * evaluate the edges [see the file]. The statistics count in expanded the
 * expansions of all the lazy searches, in iterations the searches, in
 * generated the edges discovered, in evaluated the edges evaluated, and in
 * expanded_by the edges evaluated by each worker. */

template< LazyEdgeDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class MPLP {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the MPLP)

 explicit MPLP( const D & domain , const MPLPOptions & options = {} )
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
  run_workers( f_ex , p , [ & ]( unsigned i ) {
   if( i == 0 )
    search( sh , st[ 0 ] );
   else
    evaluator( sh , st[ i ] );
   } );
  SolutionT sol;
  SearchStats & t = sol.stats;
  t = st[ 0 ];
  t.expanded_by.clear();
  t.evaluated = 0;
  for( const auto & x : st ) {
   t.evaluated += x.evaluated;
   t.expanded_by.push_back( x.evaluated );
   }
  t.stored = sh.n_nodes;
  t.stored_by.assign( 1 , sh.n_nodes );
  if( sh.solution >= 0 ) {
   const Path & pi = sh.paths[ std::size_t( sh.solution ) ];
   sol.status = Status::solved;
   sol.cost = pi.cost;
   sol.path.push_back( pi.edges.empty() ? sh.root->state
                                        : pi.edges.back()->src->state );
   for( auto it = pi.edges.rbegin() ; it != pi.edges.rend() ; ++it )
    sol.path.push_back( ( *it )->succ->state );
   }
  else
   if( sh.limit_hit )
    sol.status = Status::node_limit;
  t.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 // the state of an edge
 static constexpr std::uint8_t c_queued = 0;     // waiting for evaluation
 static constexpr std::uint8_t c_running = 1;    // being evaluated
 static constexpr std::uint8_t c_done = 2;       // evaluated, it exists
 static constexpr std::uint8_t c_missing = 3;    // evaluated, it does not

 struct MNode;

 // an edge, discovered at the first expansion of its state
 struct MEdge {
  MNode * src = nullptr;                        // its state
  MNode * succ = nullptr;                       // its successor
  Cost low{};                                   // the lower bound
  Cost cost{};                                  // the cost, once evaluated
  std::uint32_t action = 0;                     // the action
  bool front = false;                           // at the front of the queue
  std::atomic< std::uint8_t > status{ c_queued };   // read by the search
  std::vector< std::uint32_t > waiting;         // the paths through it
  };

 // a node: one per state, made by the search only
 struct MNode {
  State state{};                                // the state
  Cost h{};                                     // the heuristic
  Cost g{};                                     // g in the search sid
  MNode * parent = nullptr;                     // the parent, and
  MEdge * pedge = nullptr;                      // the edge from it
  std::uint64_t sid = 0;                        // the search of g
  std::uint64_t closed = 0;                     // the search that closed it
  MEdge * edges = nullptr;                      // the edges, once expanded
  std::uint32_t n_edges = 0;                    // how many
  bool expanded = false;                        // ever expanded
  };

 // a path found: its edges from the goal back, those not yet evaluated,
 // whether one is missing, and its cost once all are evaluated
 struct Path {
  std::vector< MEdge * > edges;
  std::uint32_t left = 0;
  bool dropped = false;
  Cost cost{};
  };

 // an item of OPEN of the lazy search: f', larger g, the order of the
 // pushes
 struct Item {
  double f = 0;
  double ng = 0;
  std::uint64_t seq = 0;
  MNode * node = nullptr;

  friend auto operator<=>( const Item & , const Item & ) = default;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {
  Lock lock;                                    // protects all but nodes
  std::condition_variable_any work;             // an edge queued, or done
  std::condition_variable_any changed;          // a cost changed, or done
  std::deque< MEdge * > front;                  // the edges of the paths
  std::deque< MEdge * > back;                   // the edges discovered
  std::vector< Path > paths;                    // Π
  std::unordered_map< std::uint64_t , std::vector< std::uint32_t > > seen;
  Cost bound{};                                 // z
  bool has_bound = false;
  std::uint64_t epoch = 0;                      // costs changed so far
  long solution = -1;                           // the path accepted
  bool done = false;                            // the search is over
  bool limit_hit = false;                       // the node limit is hit
  MNode * root = nullptr;                       // the initial state
  std::size_t n_nodes = 0;                      // the nodes made
  // the nodes and the edges, made and read by the search alone, kept to
  // the end for the path
  NodePool< MNode > pool{ std::size_t( 1 ) << 12 };
  NodeTable< MNode > table;
  std::vector< std::unique_ptr< MEdge[] > > arena;
  };

 // ends the search (lock held)
 static void stop( Shared & sh ) {
  sh.done = true;
  sh.work.notify_all();
  sh.changed.notify_all();
  }

 // accepts the path i if all its edges are evaluated and it is not longer
 // than z (lock held)
 void check( Shared & sh , std::uint32_t i ) const {
  Path & pi = sh.paths[ i ];
  if( pi.dropped || pi.left )
   return;
  Cost c = Cost( 0 );
  for( const MEdge * e : pi.edges )
   c += e->cost;
  pi.cost = c;
  bool ok = c <= sh.bound;
  if constexpr( std::is_floating_point_v< Cost > )
   ok = c <= sh.bound + 1e-9 * std::max( Cost( 1 ) , std::abs( sh.bound ) );
  if( ok ) {
   sh.solution = long( i );
   stop( sh );
   }
  }

 // the evaluation of e, without the lock, and its outcome, under it
 void evaluate( Shared & sh , MEdge * e , SearchStats & st ,
                std::unique_lock< Lock > & guard ) const {
  e->status.store( c_running , std::memory_order_relaxed );
  guard.unlock();
  bool found = false;
  Cost c{};
  f_domain.evaluate_edge( e->src->state , e->src->h , e->action ,
                          [ & ]( const State & , Cost ce , Cost ) {
   found = true;
   c = ce;
   } );
  ++st.evaluated;
  guard.lock();
  if( found ) {
   e->cost = c;
   e->status.store( c_done , std::memory_order_release );
   }
  else
   e->status.store( c_missing , std::memory_order_release );
  if( ! found || ( c != e->low ) ) {
   ++sh.epoch;
   sh.changed.notify_all();
   }
  for( std::uint32_t i : e->waiting ) {
   Path & pi = sh.paths[ i ];
   if( pi.dropped )
    continue;
   if( ! found )
    pi.dropped = true;
   else
    if( --pi.left == 0 )
     check( sh , i );
   }
  e->waiting.clear();
  }

 // the next edge to evaluate, or nullptr (lock held)
 static MEdge * next_edge( Shared & sh ) {
  for( auto * q : { & sh.front , & sh.back } )
   while( ! q->empty() ) {
    MEdge * e = q->front();
    q->pop_front();
    if( e->status.load( std::memory_order_relaxed ) == c_queued )
     return( e );
    }
  return( nullptr );
  }

 // a worker other than 0: evaluates the edges of the queue
 void evaluator( Shared & sh , SearchStats & st ) const {
  std::unique_lock< Lock > guard( sh.lock );
  while( ! sh.done ) {
   MEdge * e = next_edge( sh );
   if( ! e ) {
    sh.work.wait( guard );
    continue;
    }
   evaluate( sh , e , st , guard );
   }
  }

 // the path to goal registered in Π (lock held), with its edges not
 // evaluated in todo; false if it was already there
 bool record( Shared & sh , MNode * goal ,
              std::vector< MEdge * > & todo ) const {
  todo.clear();
  if( ! sh.has_bound || ( sh.bound < goal->g ) )
   sh.bound = goal->g;
  sh.has_bound = true;
  Path pi;
  std::uint64_t key = 0x9E3779B97F4A7C15u;
  for( MNode * x = goal ; x->pedge ; x = x->parent ) {
   pi.edges.push_back( x->pedge );
   key = ADS::splitmix64( key ^ reinterpret_cast< std::uintptr_t >(
                                                              x->pedge ) );
   }
  auto & same = sh.seen[ key ];
  for( std::uint32_t i : same )
   if( sh.paths[ i ].edges == pi.edges ) {
    check( sh , i );                  // z may have grown
    return( false );
    }
  const auto id = std::uint32_t( sh.paths.size() );
  same.push_back( id );
  for( MEdge * e : pi.edges ) {
   const auto s = e->status.load( std::memory_order_relaxed );
   if( s == c_missing )
    pi.dropped = true;
   else
    if( s != c_done ) {
     ++pi.left;
     e->waiting.push_back( id );
     todo.push_back( e );
     if( ( s == c_queued ) && ! e->front ) {
      e->front = true;
      sh.front.push_back( e );
      }
     }
   }
  sh.paths.push_back( std::move( pi ) );
  if( ! todo.empty() )
   sh.work.notify_all();
  check( sh , id );
  return( true );
  }

 // the edges of x, discovered by the lazy successors (the search alone
 // makes the nodes; the queue under the lock)
 void discover( Shared & sh , MNode * x , SearchStats & st ) const {
  auto & table = sh.table;
  std::vector< std::uint32_t > acts;
  f_domain.for_each_action( x->state , [ & ]( std::size_t a ) {
   acts.push_back( std::uint32_t( a ) );
   } );
  sh.arena.push_back( std::make_unique< MEdge[] >( acts.size() ) );
  MEdge * es = sh.arena.back().get();
  std::uint32_t n = 0;
  for( std::uint32_t a : acts )
   f_domain.lazy_edge( x->state , x->h , a ,
                       [ & ]( const State & t , Cost c , Cost ht ) {
    MEdge & e = es[ n++ ];
    e.src = x;
    e.low = c;
    e.action = a;
    e.succ = table.find_or_insert( t , f_domain.hash( t ) , [ & ] {
     MNode * y = sh.pool.allocate();
     y->state = t;
     y->h = ht;
     return( y );
     } ).first;
    } );
  x->edges = es;
  x->n_edges = n;
  st.generated += n;
  if( n ) {
   std::unique_lock< Lock > guard( sh.lock );
   for( std::uint32_t i = 0 ; i < n ; ++i )
    sh.back.push_back( es + i );
   sh.n_nodes = table.size();
   sh.work.notify_all();
   }
  }

 // the worker 0: the lazy weighted A*, repeated
 void search( Shared & sh , SearchStats & st ) const {
  const unsigned p = f_options.threads;
  const double w = f_options.weight;
  auto & table = sh.table;
  ADS::DaryHeap< Item > open;
  std::vector< MEdge * > todo;
  {
   const State r = f_domain.initial_state();
   MNode * x = sh.pool.allocate();
   x->state = r;
   x->h = Cost( f_domain.heuristic( r ) );
   table.find_or_insert( r , f_domain.hash( r ) , [ & ] { return( x ); } );
   sh.root = x;
   sh.n_nodes = 1;
   }
  std::unique_lock< Lock > guard( sh.lock , std::defer_lock );
  for( std::uint64_t sid = 1 ; ; ++sid ) {
   guard.lock();
   if( sh.done )
    break;
   const std::uint64_t epoch = sh.epoch;
   guard.unlock();
   ++st.iterations;
   // the lazy weighted A*, each edge read once, at its expansion
   open.clear();
   MNode * root = sh.root;
   root->g = Cost( 0 );
   root->sid = sid;
   root->parent = nullptr;
   root->pedge = nullptr;
   std::uint64_t pushes = 0;
   open.push( Item{ w * double( root->h ) , 0 , pushes++ , root } );
   MNode * goal = nullptr;
   while( ! open.empty() ) {
    const Item it = open.pop();
    MNode * x = it.node;
    if( ( x->closed == sid ) || ( - double( x->g ) != it.ng ) )
     continue;
    if( f_domain.is_goal( x->state ) ) {
     goal = x;
     break;
     }
    x->closed = sid;
    ++st.expanded;
    if( x->expanded )
     ++st.reexpanded;
    x->expanded = true;
    if( ! x->edges ) {
     discover( sh , x , st );
     if( f_options.node_limit && ( table.size() >= f_options.node_limit ) ) {
      guard.lock();
      sh.limit_hit = true;
      stop( sh );
      guard.unlock();
      return;
      }
     }
    for( std::uint32_t i = 0 ; i < x->n_edges ; ++i ) {
     MEdge & e = x->edges[ i ];
     const auto s = e.status.load( std::memory_order_acquire );
     if( s == c_missing )
      continue;
     MNode * y = e.succ;
     if( y->closed == sid )
      continue;
     const Cost g = x->g + ( s == c_done ? e.cost : e.low );
     if( ( y->sid != sid ) || ( g < y->g ) ) {
      y->sid = sid;
      y->g = g;
      y->parent = x;
      y->pedge = & e;
      open.push( Item{ double( g ) + w * double( y->h ) , - double( g ) ,
                       pushes++ , y } );
      }
     }
    }
   guard.lock();
   if( sh.done )
    break;
   if( ! goal ) {                     // no path, then none exists
    stop( sh );
    break;
    }
   record( sh , goal , todo );
   if( p == 1 ) {                     // the edges of the path evaluated
    for( MEdge * e : todo )
     if( ! sh.done &&
         ( e->status.load( std::memory_order_relaxed ) == c_queued ) )
      evaluate( sh , e , st , guard );
    guard.unlock();
    continue;
    }
   while( ! sh.done && ( sh.epoch == epoch ) )
    sh.changed.wait( guard );
   guard.unlock();
   }
  if( guard.owns_lock() )
   guard.unlock();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 MPLPOptions f_options;                                ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( MPLP ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MPLP */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File MPLP.h ----------------------------*/
/*--------------------------------------------------------------------------*/
