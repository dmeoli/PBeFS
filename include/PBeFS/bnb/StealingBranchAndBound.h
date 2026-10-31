/*--------------------------------------------------------------------------*/
/*--------------------- File StealingBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for StealingBranchAndBound, the parallel Branch-and-Bound
 * (B&B) whose workers keep the children they generate in a queue of their
 * own and take about half of the queue of another worker only when theirs is
 * empty, i.e., by work stealing. The problem, the options and the Solution
 * are those of BranchAndBound.h and ParallelBranchAndBound.h, and the
 * design is discussed in docs/branch-and-bound.md.
 *
 * The search is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation, the procedures of the sequential
 * search [see BranchAndBound.h] and DIVE of the parallel ones [see
 * ParallelBranchAndBound.h]; each worker \f$ i \f$ has a queue
 * \f$ Q_i \f$ in the frontier of work stealing of ADS:
 *
 * \code{.unparsed}
 *  1  procedure STEALING(R)
 *  2     the k-th node of R in Q_(k mod p)
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        d ← nil
 *  5        loop
 *  6           if d ≠ nil then n ← d;  d ← nil
 *  7           else n ← NEXT(i);  if n = nil then return
 *  8           if a limit is reached then
 *  9              n is a leaf (open);  stop all the workers;  return
 * 10           K ← PROCESS(n)
 * 11           d ← DIVE(K);  push the nodes of K but d in Q_i
 * 12  procedure NEXT(i)
 * 13     if Q_i ≠ ∅ and PRUNABLE(the least bound of Q_i) then
 * 14        the nodes of Q_i are leaves (pruned);  Q_i ← ∅
 * 15     if Q_i ≠ ∅ then return the first node of Q_i, removed
 * 16     loop
 * 17        if stopped, or no worker holds a node or works on one, return nil
 * 18        v ← a worker ≠ i at random
 * 19        if Q_v holds two nodes at least then
 * 20           move about half of them, but the first, to Q_i
 * 21           if some were moved, return the first node of Q_i, removed
 * 22        wait a little
 * \endcode
 *
 * Lines 15 to 22 are the extraction of the frontier, whose test of the end
 * at line 17 counts the workers that hold a node or work on one and the
 * thieves that are trying a steal, and is therefore exact; the owner of a
 * queue extracts from it with no lock, and a thief waits for the next
 * operation of the victim, i.e., at most the evaluation of a node. A worker
 * steals only when its own queue is empty (line 16), and therefore it stays
 * in its own part of the tree; the selection is then best-first within
 * each queue only.
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

#ifndef __PBeFS_StealingBranchAndBound
 #define __PBeFS_StealingBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstdint>

#include <memory>

#include <optional>

#include <utility>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "ADS/concurrent/StealablePriorityQueue.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS StealingBranchAndBound ----------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B by work stealing
/** Each worker extracts the nodes from its own queue and pushes there the
 * children it generates (with depth-first or plunging, the child that dives
 * is kept), so that a worker goes on in the part of the tree where it is,
 * which is what matters when moving a worker from a node to a distant one
 * costs (e.g., when the evaluator of a worker keeps the changes along the
 * path from the root to its last node). The queues are those of the
 * frontier of work stealing of ADS: the owner of a queue extracts its nodes
 * in the exact order of the selection, with no lock, and a worker whose
 * queue is empty steals about half of the queue of another one, drawn at
 * random, the half that keeps out the best node of the victim. The
 * selection is then best-first within each queue only, and the order of
 * extraction may be farther from the sequential one than that of
 * SharedBranchAndBound. The incumbent is shared, and the search ends when
 * no worker holds a node nor works on one, which the frontier tells; the
 * counter sent of the statistics holds the steals. The nodes have to be
 * copyable with moves that do not throw, since a steal moves them to the
 * thief. */

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class StealingBranchAndBound {

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

 explicit StealingBranchAndBound( const P & problem ,
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
 /// BranchAndBoundSearch::solve()]; the node \f$ i \f$ goes to the worker
 /// \f$ i \bmod p \f$

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // a node in a queue: the node and its depth, its bound being the bound of
 // the pair
 struct Held {
  Node node;          // the node
  unsigned depth;     // its depth in the tree
  };

 // the order of the selection: the key, then the push, reversed (the last
 // pushed first) unless first in, first out
 using Order = std::pair< Value , std::uint64_t >;

 using Frontier = ADS::StealingFrontier< Order , Held , std::less<> , Value >;
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

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i , clock );
  std::size_t n = 0;
  for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) )
   w[ n++ % p ]->push( Item{ std::move( l.node ) , l.bound , l.depth } );

  const auto body = [ & ]( std::size_t i ) {
   if( ! ADS::isRunCaller() )
    pin_worker( unsigned( i ) );
   w[ i ]->run();
   };
  if( f_ex )
   sh.frontier.run( *f_ex , body );
  else {
   Executor ex( p );
   sh.frontier.run( ex , body );
   }

  // the nodes left open by a limit bound the optimum too
  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( f_options.keep_frontier );
  for( auto & x : w ) {
   sol.stats.add( x->f_st );
   fr.merge( std::move( x->f_fr ) );
   if( x->f_dive )
    fr.add( std::move( x->f_dive->node ) , x->f_dive->bound ,
            x->f_dive->depth , BBLeafKind::open );
   }
  sh.frontier.leftovers( [ & ]( auto && e ) {
   fr.add( std::move( e.item.node ) , e.bound , e.item.depth ,
           BBLeafKind::open );
   } );
  sol.stats.sent = sh.frontier.steals();
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( unsigned np , Value z , std::uint64_t seed )
   : frontier( np , seed ) , inc( z ) {}

  Frontier frontier;                            // the queues
  BBIncumbent< Node , Value > inc;              // the incumbent
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its counters, its dive and its leaves; its queue is in the
 // frontier

 struct alignas( cache_line ) Worker {

  Worker( const StealingBranchAndBound & a , Shared & s , unsigned i ,
          const Stopwatch & c )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_clock( c ) ,
     f_fifo( a.f_options.selection == Selection::breadth_first ) ,
     f_fr( a.f_options.keep_frontier ) {}

  // the main loop
  void run( void ) {
   auto w = detail::bb_worker( f_a.f_p , f_id );
   const auto & o = f_a.f_options;
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };

   for( ;; ) {
    std::optional< Item > it = next();
    if( ! it )
     return;
    if( stop() ) {                     // the node stays open
     f_fr.add( std::move( it->node ) , it->bound , it->depth ,
               BBLeafKind::open );
     return;
     }

    v_kids.clear();
    detail::bb_process< P >( it->node , it->bound , it->depth , w , inc ,
                             improve , [ & ]( Node && c , Value b ,
                                              unsigned d ) {
                              v_kids.push_back( Item{ std::move( c ) , b ,
                                                      d } );
                              } , o , f_st , f_fr );
    f_sh.evaluated.fetch_add( 1 );

    const std::size_t k = detail::bb_dive( v_kids , best_bound() , inc() ,
                                           o );
    for( std::size_t j = 0 ; j < v_kids.size() ; ++j )
     if( j == k )
      f_dive.emplace( std::move( v_kids[ j ] ) );
     else
      push( std::move( v_kids[ j ] ) );
    }
   }

  // pushes a node in the queue of the worker
  void push( Item && it ) {
   const Value key = detail::bb_key( f_a.f_p , it.node , it.bound ,
                                      it.depth , f_a.f_options.selection );
   const std::uint64_t seq = ++f_seq;
   f_sh.frontier.push( f_id , Order( key , f_fifo ? seq : ~seq ) , it.bound ,
                       Held{ std::move( it.node ) , it.depth } );
   f_st.max_open = std::max( f_st.max_open ,
                             f_sh.frontier.queue( f_id ).size() );
   }

  // the least bound of the nodes of the queue of the worker, +infinity if
  // there is none
  [[nodiscard]] Value best_bound( void ) const {
   const auto b = f_sh.frontier.minBound( f_id );
   return( b ? *b : BBValueTraits< Value >::infinity() );
   }

  // the next node: the dive, the best of the queue, or the best of what is
  // stolen; none when the search is over, i.e., when no worker holds a
  // node or a limit is reached; the queue is dropped when none of its
  // nodes can improve the incumbent
  std::optional< Item > next( void ) {
   if( f_dive ) {
    std::optional< Item > it( std::move( f_dive ) );
    f_dive.reset();
    return( it );
    }
   if( const auto b = f_sh.frontier.minBound( f_id ) ;
       b && detail::bb_prunable( *b , f_sh.inc.value() , f_a.f_options ) )
    f_st.pruned_extracted += f_sh.frontier.drain( f_id , [ & ]( auto && e ) {
     f_fr.add( std::move( e.item.node ) , e.bound , e.item.depth ,
               BBLeafKind::pruned );
     } );
   auto e = f_sh.frontier.pop( f_id );
   if( ! e )
    return( std::nullopt );
   return( Item{ std::move( e->item.node ) , e->bound , e->item.depth } );
   }

  // true if a limit is reached, which ends the search
  bool stop( void ) {
   if( ! detail::bb_limit( f_a.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.frontier.stop();
   return( true );
   }

  const StealingBranchAndBound & f_a;           // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  const Stopwatch & f_clock;                    // the clock of the search
  bool f_fifo;                                  // first in, first out
  std::uint64_t f_seq = 0;                      // the pushes so far
  std::optional< Item > f_dive;                 // the node of its dive
  std::vector< Item > v_kids;                   // the children of a node
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options
 Executor * f_ex = nullptr;                   ///< given, or nullptr

 };  // end( class( StealingBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_StealingBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------- End File StealingBranchAndBound.h --------------------*/
/*--------------------------------------------------------------------------*/
