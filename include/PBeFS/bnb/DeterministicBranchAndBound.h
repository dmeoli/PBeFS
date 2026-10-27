/*--------------------------------------------------------------------------*/
/*------------------- File DeterministicBranchAndBound.h -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for DeterministicBranchAndBound, the parallel Branch-and-Bound
 * (B&B) by synchronous rounds, whose outcome does not depend on the timing
 * of the workers. The problem, the options and the Solution are those of
 * BranchAndBound.h and ParallelBranchAndBound.h, and the design is discussed
 * in docs/branch-and-bound.md.
 *
 * The search is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation and the procedures PROCESS and
 * PRUNABLE of the sequential search [see BranchAndBound.h], which a worker
 * runs here with the incumbent \f$ z_0 \f$ of the start of the round and
 * with an incumbent of its node alone; \f$ r \f$ is the number of nodes of a
 * round (round, \f$ p \f$ if 0), and the procedure NEXTROUND is run by one
 * worker while the others wait at a barrier:
 *
 * \code{.unparsed}
 *  1  procedure DETERMINISTIC(R)
 *  2     Q ← R;  z ← the given incumbent, or ∞;  N ← ()
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        loop
 *  5           barrier, at which one worker runs NEXTROUND()
 *  6           if done then return
 *  7           for j = i, i + p, i + 2p, ... while j < |N| do
 *  8              (K_j, v_j) ← PROCESS(N_j) against z₀, v_j being the best
 *                    feasible value of N_j below z₀, if any
 *  9  procedure NEXTROUND()
 * 10     for j = 0, ..., |N| − 1 do                  ▷ in the order of N
 * 11        if v_j < z then z ← v_j;  N_j is the node of z
 * 12        for each m ∈ K_j do
 * 13           if PRUNABLE(β(m)) then m is a leaf (pruned) else push m in Q
 * 14     if a limit is reached, Q = ∅ or PRUNABLE(the least bound of Q) then
 * 15        done ← true;  return
 * 16     N ← the first min(r, evaluations left) nodes of Q, removed
 * 17     z₀ ← z
 * \endcode
 *
 * The nodes of a round, their order and the incumbent \f$ z_0 \f$ against
 * which they are evaluated are fixed before the round starts (lines 16 and
 * 17), and the outcomes are merged in the order of the nodes (line 10);
 * hence, the nodes evaluated, the incumbent and the bound depend only on
 * \f$ r \f$, and neither on the timing nor on the number of the workers,
 * provided that the evaluation of a node by the worker \f$ i \f$ depends
 * only on the node and on the nodes evaluated by \f$ i \f$ before. The node
 * limit is exact (line 16), while the time limit, checked at line 14, makes
 * the end depend on the timing; no dive is taken.
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

#ifndef __PBeFS_DeterministicBranchAndBound
 #define __PBeFS_DeterministicBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <barrier>

#include <exception>

#include <mutex>

#include <optional>

#include <utility>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------- CLASS DeterministicBranchAndBound --------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B by synchronous rounds
/** At each round the round best open nodes (by the selection of the
 * options, round being the number of workers if 0) are extracted from one
 * queue, and the worker \f$ i \f$ evaluates and branches the nodes
 * \f$ i , i + p , i + 2 p , \ldots \f$ of the round, in this order, against
 * the incumbent of the start of the round. When all of them are done, the
 * outcomes are merged in the order of the nodes of the round: the feasible
 * values update the incumbent, the children that cannot improve it are
 * pruned, and the others are pushed in the queue. Hence, the nodes
 * evaluated, the incumbent and the bound depend neither on the timing of the
 * workers nor on their number, but only on round, provided that the
 * evaluation of a node by the worker i depends only on the node and on the
 * nodes evaluated by that worker before (which make_worker() allows), and
 * the search can be repeated exactly, as the experiments of a paper need.
 *
 * The price is the wait of each round for its slowest node, which is
 * small when the evaluation dominates and its time varies little among the
 * nodes; also, a node is evaluated against an incumbent that may be one
 * round old. The node limit is exact, the last round being shortened to
 * it, while the time limit is checked between two rounds and therefore
 * makes the end depend on the timing. A dive is not taken, the plunging
 * selection extracting the nodes by their bound as best-first does. */

template< BranchAndBound P , class Runtime = ThreadRuntime >
class DeterministicBranchAndBound {

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

 explicit DeterministicBranchAndBound( const P & problem ,
                                       const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  }

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

 using Item = typename BBQueue< Node , Value >::Item;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a node of a round, and what its evaluation gives

 struct Slot {
  Item it;                                      // the node
  std::vector< Item > kids;                     // its children
  std::optional< std::pair< Value , Node > > best;  // its best value
  BBStats st;                                   // its counters
  detail::BBFrontier< Node , Value > fr;        // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const Options & o = f_options;
  const unsigned p = o.threads;
  const std::size_t round = o.round ? o.round : p;
  Stopwatch clock;

  BBQueue< Node , Value > q( o.selection );
  SolutionT sol;
  BBStats & st = sol.stats;
  st.evaluated_by.assign( p , 0 );
  const auto push = [ & ]( Item && it ) {
   const Value key = detail::bb_key( f_p , it.node , it.bound , it.depth ,
                                      o.selection );
   q.push( std::move( it.node ) , it.bound , it.depth , key );
   st.max_open = std::max( st.max_open , q.size() );
   };
  for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) )
   push( Item{ std::move( l.node ) , l.bound , l.depth } );

  BBIncumbent< Node , Value > inc( o.incumbent ?
                                   detail::bb_signed< P >( *o.incumbent ) :
                                   inf );
  detail::BBFrontier< Node , Value > fr( o.keep_frontier );
  std::vector< Slot > slots( round );
  std::size_t nslots = 0;               // the nodes of the round
  Value snapshot = inc.value();         // the incumbent of the round
  bool done = false , limit = false , late = false;
  std::exception_ptr error;
  std::mutex error_m;

  // merges the outcomes of the round just ended, in the order of its
  // nodes, and extracts the nodes of the next one; it is run by one worker
  // while the others wait
  const auto next_round = [ & ]( void ) noexcept {
   try {
    for( std::size_t j = 0 ; j < nslots ; ++j ) {
     Slot & s = slots[ j ];
     s.st.incumbents = 0;              // counted here, by the merge
     add_counters( st , s.st );
     st.evaluated_by[ j % p ] += s.st.evaluated;
     if( s.best && inc.template improve< P >( s.best->first ,
                                              s.best->second , o ) )
      ++st.incumbents;
     fr.merge( std::move( s.fr ) );
     for( Item & k : s.kids )
      if( detail::bb_prunable( k.bound , inc.value() , o ) ) {
       ++st.pruned_generated;
       fr.add( std::move( k.node ) , k.bound , k.depth ,
               BBLeafKind::pruned );
       }
      else
       push( std::move( k ) );
     }
    nslots = 0;

    if( error ) {
     done = true;
     return;
     }
    if( o.node_limit && ( st.evaluated >= o.node_limit ) ) {
     limit = done = true;
     return;
     }
    if( ( o.time_limit > 0 ) && ( clock.seconds() >= o.time_limit ) ) {
     limit = late = done = true;
     return;
     }
    if( q.empty() ) {
     done = true;
     return;
     }
    if( detail::bb_prunable( q.best_bound() , inc.value() , o ) ) {
     st.pruned_extracted += fr.add( q , BBLeafKind::pruned );
     done = true;                      // nothing can improve
     return;
     }

    std::size_t r = round;
    if( o.node_limit )
     r = std::min( r , std::size_t( o.node_limit - st.evaluated ) );
    while( ( nslots < r ) && ! q.empty() ) {
     Slot & s = slots[ nslots++ ];
     s.it = q.pop();
     s.kids.clear();
     s.best.reset();
     s.st = BBStats();
     s.fr = detail::BBFrontier< Node , Value >( o.keep_frontier );
     }
    snapshot = inc.value();
    }
   catch( ... ) {
    std::lock_guard< std::mutex > guard( error_m );
    if( ! error )
     error = std::current_exception();
    done = true;
    }
   };

  std::barrier sync( std::ptrdiff_t( p ) , next_round );
  Runtime::run( p , [ & ]( unsigned i ) {
   auto w = detail::bb_worker( f_p , i );
   for( ;; ) {
    sync.arrive_and_wait();
    if( done )
     return;
    for( std::size_t j = i ; j < nslots ; j += p ) {
     Slot & s = slots[ j ];
     try {
      detail::bb_process< P >( s.it.node , s.it.bound , s.it.depth , w ,
                               [ & ] { return( snapshot ); } ,
                               [ & ]( Value v , const Node & n ) {
                                if( ( v >= snapshot ) ||
                                    ( s.best && ( v >= s.best->first ) ) )
                                 return( false );
                                s.best.emplace( v , n );
                                return( true );
                                } ,
                               [ & ]( Node && c , Value b , unsigned d ) {
                                s.kids.push_back( Item{ std::move( c ) , b ,
                                                        d } );
                                } , o , s.st , s.fr );
      }
     catch( ... ) {
      std::lock_guard< std::mutex > guard( error_m );
      if( ! error )
       error = std::current_exception();
      }
     }
    }
   } );
  if( error )
   std::rethrow_exception( error );

  fr.add( q , BBLeafKind::open );       // left open by a limit
  detail::bb_finish< P >( sol , inc , fr , limit , late );
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // adds the counters of b to those of a, but for evaluated_by

 static void add_counters( BBStats & a , const BBStats & b ) {
  a.evaluated += b.evaluated;
  a.branched += b.branched;
  a.generated += b.generated;
  a.infeasible += b.infeasible;
  a.pruned_generated += b.pruned_generated;
  a.pruned_extracted += b.pruned_extracted;
  a.pruned_evaluated += b.pruned_evaluated;
  a.incumbents += b.incumbents;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options

 };  // end( class( DeterministicBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_DeterministicBranchAndBound */

/*--------------------------------------------------------------------------*/
/*----------------- End File DeterministicBranchAndBound.h -----------------*/
/*--------------------------------------------------------------------------*/
