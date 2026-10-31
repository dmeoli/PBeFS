/*--------------------------------------------------------------------------*/
/*---------------------------- File FarmAStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template FarmAStar, a centralized parallel A*
 * written as a farm with feedback on an executor of ADS: a master holds OPEN
 * and CLOSED, and the workers expand the nodes it sends them.
 *
 * The search is the centralized scheme of section 3.2 of Fukunaga et al.
 * (2017), written as a farm with feedback, with the notation of the page
 * "Best-first search and A*" of the documentation: the master alone holds
 * OPEN, CLOSED, the nodes and the incumbent \f$ z \f$, out counts the tasks
 * at the workers, and INSERT is the one of HDA* [see HDAStar.h] on the lists
 * of the master:
 *
 * \code{.unparsed}
 *  1  procedure MASTER(s₀, p)
 *  2     z ← ∞;  out ← 0;  INSERT(s₀, 0, nil)
 *  3     loop
 *  4        DISPATCH()
 *  5        if out = 0 and OPEN has no node of f < z then return
 *  6        receive a task T back from a worker;  out ← out − 1
 *  7        for each successor v of a node u of T, with c(u, v) and h(v) do
 *  8           if g(u) + c(u, v) + h(v) < z then INSERT(v, g(u) + c(u, v), u)
 *  9  procedure DISPATCH()
 * 10     while out < p and OPEN has a node of f < z do
 * 11        T ← ∅
 * 12        while |T| < batch and OPEN has a node of f < z do
 * 13           u ← the first node of OPEN, removed;  close u
 * 14           if u is a goal then z ← min(z, g(u)), recording u if g(u) < z
 * 15           else T ← T ∪ {u}
 * 16        if T = ∅ then return
 * 17        send T to a free worker;  out ← out + 1
 * 18  procedure WORKER(T)                       ▷ on each of the p workers
 * 19     for each u ∈ T, add to T the successors v of u, with c(u, v), h(v)
 * 20     send T back to the master
 * \endcode
 *
 * The test of line 5 is exact, since the master alone changes OPEN and
 * \f$ z \f$ and a node is either in OPEN or in a task that is out; the
 * workers only generate the successors, and therefore the scheme pays when
 * the generation costs much more than the bookkeeping of line 8, which the
 * master does alone.
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

#ifndef __PBeFS_FarmAStar
 #define __PBeFS_FarmAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <memory>

#include <stdexcept>

#include <tuple>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/HDAStar.h"

#include "ADS/parallel/Farm.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS FarmAStar -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a master-worker parallel A* on a farm with feedback
/** The master of the farm alone owns OPEN, CLOSED and
 * the nodes, so it needs no lock. It sends the best nodes of OPEN, batch at
 * a time, to the workers that are free (the farm schedules on demand); a
 * worker generates the successors of its nodes with their costs and
 * heuristic values and sends them back on the feedback channel, and the
 * master checks them against CLOSED and puts the new and improved ones in
 * OPEN. It is the centralized approach of Kumar, Ramesh and Rao (1988), the
 * same as SPAStar but with the master in place of the locks, close to the
 * K-best-first search of Felner, Kraus and Korf (2003) without its
 * synchronization at each round. A node may be expanded before its g is
 * optimal, and then re-expanded; the incumbent prunes, and the search ends
 * when no node of OPEN has f below the incumbent and no node is out at a
 * worker. Since the master does all the bookkeeping, FarmAStar pays only
 * when the expansions are slow.
 *
 * The farm runs threads + 1 threads, those of an Executor of ADS [see
 * ADS::feedbackFarm()]: with PBEFS_PIN the master takes the first CPU of the
 * affinity mask and the worker \f$ i \f$ the CPU \f$ i + 1 \f$ of it,
 * counting from 0 [see pin_worker()], but for the thread of the caller. */

template< Domain D , ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class FarmAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the FarmAStar);
 /// threads is the number of workers, batch the nodes per task

 explicit FarmAStar( const D & domain , const ParallelOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the master and the workers, which must
 /// outlive the searches and have at least threads + 1 threads; without
 /// it, each search builds one of threads + 1 threads

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  const unsigned w = f_options.threads;
  Master master( *this );
  std::vector< char > pinned( w , 0 );
  const auto start = [ & ]( auto && emit , std::size_t ) {
   if( ! ADS::isRunCaller() )
    pin_worker( 0 );
   master.start( emit );
   };
  const auto back = [ & ]( Task * t , auto && emit ) {
   return( master.back( t , emit ) );
   };
  const auto work = [ & ]( Task * t , std::size_t i ) {
   if( ! pinned[ i ] ) {                   // each index on one thread
    pinned[ i ] = 1;
    if( ! ADS::isRunCaller() )
     pin_worker( unsigned( i + 1 ) );
    }
   expand( t );
   return( t );
   };
  if( f_ex )
   ADS::feedbackFarm< Task * , Task * >( *f_ex , w , start , back , work );
  else {
   Executor ex( w + 1 );
   ADS::feedbackFarm< Task * , Task * >( ex , w , start , back , work );
   }

  SolutionT sol;
  sol.stats = master.f_st;
  sol.stats.stored = master.f_pool.size();
  sol.stats.expanded_by = { sol.stats.expanded };
  sol.stats.stored_by = { sol.stats.stored };
  if( master.f_limit_hit )
   sol.status = Status::node_limit;
  else
   if( master.f_goal ) {
    sol.status = Status::solved;
    sol.cost = master.f_goal->g;
    sol.path = trace_path( master.f_goal );
    }
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // a task: nodes to expand, and then their successors
 struct Task {
  std::vector< NodeT * > nodes;                           // the nodes
  std::vector< std::tuple< std::size_t , State , Cost , Cost > > succ;
                         // per successor: its node, state, edge cost, h
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the master: OPEN, CLOSED, the nodes and the counters

 struct Master {

  explicit Master( const FarmAStar & a ) : f_a( a ) {}

  Master( const Master & ) = delete;
  Master & operator=( const Master & ) = delete;

  ~Master() {
   for( Task * t : v_free )
    delete t;
   }

  // the start: the root, and the first tasks
  template< class E >
  void start( E && emit ) {
   const State s0 = f_a.f_domain.initial_state();
   insert( s0 , Cost( 0 ) , Cost( f_a.f_domain.heuristic( s0 ) ) ,
           nullptr , f_a.f_domain.hash( s0 ) );
   dispatch( emit );
   }

  // a task back, with the successors of its nodes; false once nothing is
  // out at the workers and nothing is left to send
  template< class E >
  bool back( Task * t , E && emit ) {
   --f_out;
   v_hash.resize( t->succ.size() );       // prefetch the batch first
   for( std::size_t i = 0 ; i < t->succ.size() ; ++i ) {
    v_hash[ i ] = f_a.f_domain.hash( std::get< 1 >( t->succ[ i ] ) );
    f_table.prefetch( v_hash[ i ] );
    }
   for( std::size_t i = 0 ; i < t->succ.size() ; ++i ) {
    const auto & [ k , s , c , h ] = t->succ[ i ];
    ++f_st.generated;
    NodeT * parent = t->nodes[ k ];
    const Cost g = parent->g + c;
    if( g + h < f_incumbent )
     insert( s , g , h , parent , v_hash[ i ] );
    }
   if( f_a.f_options.node_limit &&
       ( f_pool.size() >= f_a.f_options.node_limit ) )
    f_limit_hit = true;
   v_free.push_back( t );
   dispatch( emit );
   return( ( f_out > 0 ) || ( ! f_limit_hit && has_work() ) );
   }

  // sends tasks while there are nodes and free workers
  template< class E >
  void dispatch( E && emit ) {
   while( ! f_limit_hit && ( f_out < f_a.f_options.threads ) ) {
    Task * t = nullptr;
    while( NodeT * n = next() ) {
     n->closed = true;
     ++f_st.expanded;
     if( n->expanded )
      ++f_st.reexpanded;
     n->expanded = true;
     if( f_a.f_domain.is_goal( n->state ) ) {
      if( n->g < f_incumbent ) {
       f_incumbent = n->g;
       f_goal = n;
       }
      continue;
      }
     if( ! t ) {
      if( v_free.empty() )
       t = new Task;
      else {
       t = v_free.back();
       v_free.pop_back();
       }
      t->nodes.clear();
      t->succ.clear();
      }
     t->nodes.push_back( n );
     if( t->nodes.size() >= f_a.f_options.batch )
      break;
     }
    if( ! t )
     return;
    ++f_out;
    emit( t );
    }
   }

  // true if OPEN has a node with f below the incumbent
  bool has_work( void ) {
   while( ! f_open.empty() && f_open.top().stale() )
    f_open.pop();
   return( ! f_open.empty() && ( f_open.top_f() < f_incumbent ) );
   }

  // the best node of OPEN with f below the incumbent, removed, or nullptr
  NodeT * next( void ) {
   return( has_work() ? f_open.pop().node : nullptr );
   }

  // the path to s, whose hash is hash, of cost g through parent
  void insert( const State & s , Cost g , Cost h , NodeT * parent ,
               std::uint64_t hash ) {
   auto [ x , fresh ] = f_table.find_or_insert( s , hash , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ s , g , h , parent , false , false };
    return( y );
    } );
   if( fresh )
    f_open.push( x );
   else
    if( g < x->g ) {
     x->g = g;
     x->parent = parent;
     x->closed = false;
     f_open.push( x );
     }
    else
     ++f_st.duplicates;
   }

  const FarmAStar & f_a;                        // the search
  NodePool< NodeT > f_pool;                     // the nodes
  NodeTable< NodeT > f_table;                   // CLOSED
  Open f_open;                                  // OPEN
  std::vector< Task * > v_free;                 // tasks to reuse
  std::vector< std::uint64_t > v_hash;          // hashes of a task back
  unsigned f_out = 0;                           // tasks at the workers
  Cost f_incumbent = infinite_cost< Cost >();   // the best cost found
  NodeT * f_goal = nullptr;                     // the best goal found
  bool f_limit_hit = false;                     // the node limit is hit
  SearchStats f_st;                             // the counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the work of a worker: the successors of the nodes of a task

 void expand( Task * t ) const {
  for( std::size_t k = 0 ; k < t->nodes.size() ; ++k ) {
   const NodeT * n = t->nodes[ k ];
   f_domain.for_each_successor( n->state , n->h ,
                                [ & ]( const State & s , Cost c , Cost h ) {
    t->succ.emplace_back( k , s , c , h );
    } );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 ParallelOptions f_options;                            ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( FarmAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_FarmAStar */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File FarmAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
