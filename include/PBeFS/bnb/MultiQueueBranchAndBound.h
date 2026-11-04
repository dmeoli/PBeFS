/*--------------------------------------------------------------------------*/
/*-------------------- File MultiQueueBranchAndBound.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MultiQueueBranchAndBound, the parallel Branch-and-Bound
 * (B&B) whose workers share the open nodes in the MultiQueue of ADS, a
 * relaxed priority queue of \f$ 2 p \f$ sequential queues each under its
 * own lock, as the workers of MQAStar do. The problem, the options and the
 * Solution are those of BranchAndBound.h and ParallelBranchAndBound.h, and
 * the design is discussed in docs/branch-and-bound.md.
 *
 * The search is that of SharedBranchAndBound [see ParallelBranchAndBound.h]
 * but for the queue and the test of the end: a pop takes the better of the
 * first nodes of two queues drawn at random, and therefore not always the
 * first of all, and the counter \f$ m \f$ of the nodes that are in the
 * queue or held by a worker, which is atomic, takes the place of busy and
 * of the lock \f$ L \f$:
 *
 * \code{.unparsed}
 *  1  procedure MULTIQUEUE(R)
 *  2     Q ← R;  m ← |R|;  done ← false
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        d ← nil
 *  5        while not done do
 *  6           if d ≠ nil then n ← d;  d ← nil
 *  7           else
 *  8              n ← a pop of Q, possibly nothing
 *  9              if no node was popped then
 * 10                 if m = 0 then done ← true else wait a little
 * 11                 continue
 * 12              if PRUNABLE(the bound of n) then
 * 13                 n is a leaf (pruned);  m ← m − 1;  continue
 * 14           if a limit is reached then n is a leaf (open);  done ← true
 * 15           K ← PROCESS(n);  d ← DIVE(K)
 * 16           m ← m + |K| − 1;  push the nodes of K but d in one queue of Q
 * \endcode
 *
 * The test at line 10 is exact: at line 16, \f$ m \f$ counts the children
 * and drops the node processed in one atomic step, before the pushes (the
 * child that dives is held by the worker, and counted as such), so that
 * \f$ m = 0 \f$ means that \f$ Q \f$ is empty and that no worker holds a
 * node, whatever the order of the pops. Since the first node of \f$ Q \f$
 * is not known, a node is pruned when popped (line 12) rather than with
 * all the others, and DIVE compares the children with the bound of the
 * node processed instead of the least bound of \f$ Q \f$. The children of
 * a node go into one queue under one lock, as by a stickiness of their
 * number.
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

#ifndef __PBeFS_MultiQueueBranchAndBound
 #define __PBeFS_MultiQueueBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstdint>

#include <iterator>

#include <limits>

#include <memory>

#include <optional>

#include <utility>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "ADS/concurrent/MultiQueue.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CONCEPT MultiQueueValue -------------------------*/
/*--------------------------------------------------------------------------*/
/// the values that the MultiQueue of ADS takes as keys: arithmetic, with a
/// greatest value, and atomic with no lock

template< class V >
concept MultiQueueValue = std::numeric_limits< V >::is_specialized &&
                          std::atomic< V >::is_always_lock_free;

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS MultiQueueBranchAndBound ---------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B on the MultiQueue of ADS
/** The workers pop the nodes from the MultiQueue of ADS, evaluate and
 * branch them, and push the children into one of its queues under one lock
 * (with depth-first or plunging, the child that dives is kept by the
 * worker); no lock is shared by all the workers, and two of them wait for
 * each other only when they draw the same queue. The order of extraction
 * is relaxed: a pop takes a node whose rank among the open ones is less
 * than \f$ \frac{5}{6} \cdot 2 p \f$ in expectation, and the nodes of
 * equal key come out in no given order. The search ends when the atomic
 * counter of the nodes in the queue or held by a worker falls to 0 [see
 * the pseudocode above]. The values have to be arithmetic and atomic with
 * no lock [see MultiQueueValue], since the least key of each queue is read
 * with no lock, and the nodes default constructible, since the buffers of
 * the queues are arrays. */

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool >
 requires MultiQueueValue< bb_value_t< P > >
class MultiQueueBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search)

 explicit MultiQueueBranchAndBound( const P & problem ,
                                    const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution [see
 /// BranchAndBoundSearch::solve()]

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // a node in the queue: the node and its depth, its bound beside them
 struct Held {
  Node node;          // the node
  Value bound;        // its bound, in the minimization
  unsigned depth;     // its depth in the tree
  };

 using Mq = ADS::MultiQueue< Value , Held >;
 using Item = typename BBQueue< Node , Value >::Item;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const unsigned p = f_options.threads;
  Stopwatch clock;
  const Value start = f_options.incumbent ?
                      detail::bb_signed< P >( *f_options.incumbent ) : inf;
  Shared sh( p , start , f_options.seed );
  {
   auto hd = sh.q.handle();
   std::size_t n = 0;
   for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) ) {
    const Value key = detail::bb_key( f_p , l.node , l.bound , l.depth ,
                                       f_options.selection );
    hd.push( key , Held{ std::move( l.node ) , l.bound , l.depth } );
    ++n;
    }
   sh.held.store( n );
   }

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , clock );
  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run( i ); } );

  // the nodes left open by a limit bound the optimum too
  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( f_options.keep_frontier );
  for( auto & x : w ) {
   sol.stats.add( x->f_st );
   fr.merge( std::move( x->f_fr ) );
   }
  auto hd = sh.q.handle();
  for( auto e = hd.tryPop() ; e ; e = hd.tryPop() )
   fr.add( std::move( e->item.node ) , e->item.bound , e->item.depth ,
           BBLeafKind::open );
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( unsigned np , Value z , std::uint64_t seed )
   : q( 2 * std::size_t( np ) , 1 , seed ) , inc( z ) {}

  Mq q;                                         // the open nodes
  BBIncumbent< Node , Value > inc;              // the incumbent
  // the nodes in q or held by a worker
  alignas( cache_line ) std::atomic< std::size_t > held{ 0 };
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  alignas( cache_line ) std::atomic< bool > done{ false };
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its counters and its leaves

 struct alignas( cache_line ) Worker {

  Worker( const MultiQueueBranchAndBound & a , Shared & s ,
          const Stopwatch & c )
   : f_a( a ) , f_sh( s ) , f_clock( c ) ,
     f_fr( a.f_options.keep_frontier ) {}

  // the main loop
  void run( unsigned i ) {
   auto w = detail::bb_worker( f_a.f_p , i );
   const auto & o = f_a.f_options;
   const Selection sel = o.selection;
   auto hd = f_sh.q.handle();
   ADS::Backoff backoff;
   std::optional< Item > dive;               // the node of the dive
   std::vector< Item > kids;                 // the children of a node
   std::vector< typename Mq::Entry > out;    // the children pushed
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };

   while( ! f_sh.done.load( std::memory_order_acquire ) ) {
    std::optional< Item > it;
    if( dive ) {
     it.emplace( std::move( *dive ) );
     dive.reset();
     }
    else {
     auto e = hd.tryPop();
     if( ! e ) {
      if( f_sh.held.load( std::memory_order_acquire ) == 0 ) {
       f_sh.done.store( true , std::memory_order_release );
       break;
       }
      backoff.pause();
      continue;
      }
     backoff.reset();
     if( detail::bb_prunable( e->item.bound , inc() , o ) ) {
      f_fr.add( std::move( e->item.node ) , e->item.bound , e->item.depth ,
                BBLeafKind::pruned );
      ++f_st.pruned_extracted;
      f_sh.held.fetch_sub( 1 , std::memory_order_acq_rel );
      continue;
      }
     it.emplace( Item{ std::move( e->item.node ) , e->item.bound ,
                       e->item.depth } );
     }

    if( stop() ) {                     // the node stays open
     f_fr.add( std::move( it->node ) , it->bound , it->depth ,
               BBLeafKind::open );
     return;
     }

    kids.clear();
    detail::bb_process< P >( it->node , it->bound , it->depth , w , inc ,
                             improve , [ & ]( Node && c , Value b ,
                                              unsigned d ) {
                              kids.push_back( Item{ std::move( c ) , b , d } );
                              } , o , f_st , f_fr );
    if( o.node_limit )
     f_sh.evaluated.fetch_add( 1 );

    const std::size_t k = detail::bb_dive( kids , it->bound , inc() , o );
    out.clear();
    for( std::size_t j = 0 ; j < kids.size() ; ++j )
     if( j == k )
      dive.emplace( std::move( kids[ j ] ) );
     else
      out.push_back( typename Mq::Entry{
       detail::bb_key( f_a.f_p , kids[ j ].node , kids[ j ].bound ,
                       kids[ j ].depth , sel ) ,
       Held{ std::move( kids[ j ].node ) , kids[ j ].bound ,
             kids[ j ].depth } } );
    // the children counted and the node processed dropped in one step,
    // before the pushes, so that held never misses a node [see line 16]
    if( kids.size() != 1 ) {
     const std::size_t m = f_sh.held.fetch_add( kids.size() - 1 ,
                                                std::memory_order_acq_rel );
     f_st.max_open = std::max( f_st.max_open , m + kids.size() - 1 );
     }
    hd.pushAll( std::make_move_iterator( out.begin() ) ,
                std::make_move_iterator( out.end() ) );
    }
   if( dive )                          // stopped by a limit of another
    f_fr.add( std::move( dive->node ) , dive->bound , dive->depth ,
              BBLeafKind::open );
   }

  // true if a limit is reached, which ends the search
  bool stop( void ) {
   if( ! detail::bb_limit( f_a.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.done.store( true , std::memory_order_release );
   return( true );
   }

  const MultiQueueBranchAndBound & f_a;         // the search
  Shared & f_sh;                                // the shared data
  const Stopwatch & f_clock;                    // the clock of the search
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options
 Executor * f_ex = nullptr;                   ///< given, or nullptr

 };  // end( class( MultiQueueBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MultiQueueBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------ End File MultiQueueBranchAndBound.h -------------------*/
/*--------------------------------------------------------------------------*/
