/*--------------------------------------------------------------------------*/
/*--------------------- File SelfSplitBranchAndBound.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for SelfSplitBranchAndBound, the parallel Branch-and-Bound
 * (B&B) whose workers share nothing but the nodes of a first sampling of
 * the tree, as in the SelfSplit of Fischetti, Monaci and Salvagnin (2014).
 * The problem, the options and the Solution are those of BranchAndBound.h
 * and ParallelBranchAndBound.h, and the design is discussed in
 * docs/branch-and-bound.md.
 *
 * The search is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation and the procedures of the
 * sequential search [see BranchAndBound.h]; \f$ \sigma \f$ is the number of
 * open nodes per worker at which the sampling stops (split, 4 if 0):
 *
 * \code{.unparsed}
 *  1  procedure SELFSPLIT(R)
 *  2     Q ← R;  z ← the given incumbent, or ∞     ▷ the sampling, by worker 0
 *  3     while Q ≠ ∅ and |Q| < σ p do
 *  4        if a limit is reached then break
 *  5        if PRUNABLE(the least bound of Q) then
 *  6           the nodes of Q are leaves (pruned);  Q ← ∅;  break
 *  7        n ← the first node of Q, removed;  Q ← Q ∪ PROCESS(n)
 *  8     for k = 0, 1, ..., while Q ≠ ∅ do           ▷ the split
 *  9        add the first node of Q, removed, to R_(k mod p)
 * 10     in parallel, on each worker i = 0, ..., p − 1 do
 * 11        S_i ← B&B(R_i) on the evaluator of i, from the incumbent z, with
 *              its share of the node limit and the time left
 * 12     z ← min(z, the values of the S_i), the node of the least i on ties
 * 13     return z, its node, the leaves of the sampling and of the S_i, and
 *           the least of their bounds
 * \endcode
 *
 * The workers neither share the incumbent nor exchange nodes between lines
 * 10 and 11, and the nodes of line 9 are those of the sampling in the order
 * of the selection; hence, the outcome depends only on \f$ \sigma \f$ and
 * on \f$ p \f$, provided that the evaluation of a node by a worker depends
 * only on the node and on the nodes evaluated by that worker before. The
 * sampling uses the selection of the options, with best-first in place of
 * plunging, and the search is correct since the subsets of the nodes of the
 * split cover what the sampling left open.
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

#ifndef __PBeFS_SelfSplitBranchAndBound
 #define __PBeFS_SelfSplitBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

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
/*--------------------- CLASS SelfSplitBranchAndBound ----------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B by a split of the tree after a sampling
/** The tree is first explored sequentially (the *sampling*), by the
 * selection of the options and without dives, until the open nodes are split
 * per worker (4 if 0) or the tree is solved. The open nodes are then
 * extracted in the order of the selection, the \f$ k \f$-th going to the
 * worker \f$ k \bmod p \f$, and each worker searches the subtrees of its
 * nodes sequentially [see BranchAndBoundSearch], with its own evaluator, its
 * own incumbent (starting from that of the sampling) and its share of the
 * node limit. The Solution has the best incumbent of the workers (that of
 * the lowest worker among equal values), the least of their bounds, and
 * their leaves and counters together with those of the sampling.
 *
 * The workers neither share the incumbent nor exchange nodes, and
 * therefore the outcome depends only on split and on the number p of the
 * workers, provided that the evaluation of a node by a worker depends only
 * on the node and on the nodes that worker evaluated before; the sampling
 * is done once, by the evaluator of the worker 0. The price is the balance
 * of the load, which is that of the subtrees of the split, and an
 * incumbent found by a worker that does not prune the nodes of the others.
 * As in the deterministic search, the time limit makes the end depend on
 * the timing. */

template< BranchAndBound P , class Runtime = ThreadRuntime >
class SelfSplitBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using SolutionT = BBSolution< Node >;          ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search)

 explicit SelfSplitBranchAndBound( const P & problem ,
                                   const ParallelBBOptions & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  if( f_options.split == 0 )
   f_options.split = 4;
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution [see
 /// BranchAndBoundSearch::solve()]

 [[nodiscard]] SolutionT solve( std::vector< BBLeaf< Node > > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Item = typename BBQueue< Node >::Item;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 SolutionT run( std::optional< std::vector< BBLeaf< Node > > > seeds ) const
 {
  constexpr double sgn = detail::bb_sign< P >();
  constexpr double inf = std::numeric_limits< double >::infinity();
  const ParallelBBOptions & o = f_options;
  const unsigned p = o.threads;
  Stopwatch clock;

  // the sampling, by the evaluator of the worker 0
  SolutionT sol;
  BBStats & st = sol.stats;
  const Selection sel = o.selection == Selection::plunging ?
                        Selection::best_first : o.selection;
  BBQueue< Node > q( sel );
  const auto push = [ & ]( Item && it ) {
   const double key = detail::bb_key( f_p , it.node , it.bound , it.depth ,
                                      sel );
   q.push( std::move( it.node ) , it.bound , it.depth , key );
   st.max_open = std::max( st.max_open , q.size() );
   };
  for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) )
   push( Item{ std::move( l.node ) , l.bound , l.depth } );

  double inc = o.incumbent ? sgn * *o.incumbent : inf;
  std::optional< Node > inc_node;
  detail::BBFrontier< Node > fr( o.keep_frontier );
  bool limit = false , late = false;
  {
   auto w = detail::bb_worker( f_p , 0u );
   const auto improve = [ & ]( double v , const Node & n ) {
    if( v >= inc )
     return( false );
    inc = v;
    inc_node = n;
    if( o.on_incumbent )
     o.on_incumbent( sgn * v );
    return( true );
    };
   std::vector< Item > kids;
   while( ! q.empty() && ( q.size() < o.split * p ) ) {
    if( o.node_limit && ( st.evaluated >= o.node_limit ) ) {
     limit = true;
     break;
     }
    if( ( o.time_limit > 0 ) && ( clock.seconds() >= o.time_limit ) ) {
     limit = late = true;
     break;
     }
    if( detail::bb_prunable( q.best_bound() , inc , o ) ) {
     st.pruned_extracted += fr.add( q , BBLeafKind::pruned );
     break;                            // nothing can improve
     }
    Item it = q.pop();
    kids.clear();
    detail::bb_process< P >( it.node , it.bound , it.depth , w ,
                             [ & ] { return( inc ); } , improve ,
                             [ & ]( Node && c , double b , unsigned d ) {
                              kids.push_back( Item{ std::move( c ) , b ,
                                                    d } );
                              } , o , st , fr );
    for( Item & k : kids )
     push( std::move( k ) );
    }
  }

  // the split: the k-th open node, in the order of the selection, to the
  // worker k mod p, unless a limit or the incumbent ends the search here
  std::vector< std::vector< BBLeaf< Node > > > part( p );
  if( ! limit && ! q.empty() &&
      detail::bb_prunable( q.best_bound() , inc , o ) )
   st.pruned_extracted += fr.add( q , BBLeafKind::pruned );
  std::uint64_t rest = 0;                       // the evaluations left
  double time = 0;                              // the time left
  if( o.node_limit && ! limit )
   rest = o.node_limit - st.evaluated;
  if( ( o.time_limit > 0 ) && ! limit ) {
   time = o.time_limit - clock.seconds();
   if( time <= 0 )
    limit = late = true;
   }
  if( ! limit )
   for( std::size_t k = 0 ; ! q.empty() ; ++k ) {
    Item it = q.pop();
    const unsigned i = unsigned( k % p );
    if( o.node_limit && ( share( rest , p , i ) == 0 ) ) {
     limit = true;                      // a node left open by the limit
     fr.add( std::move( it.node ) , it.bound , it.depth , BBLeafKind::open );
     }
    else
     part[ i ].push_back( BBLeaf< Node >{ std::move( it.node ) ,
                                          sgn * it.bound , it.depth ,
                                          BBLeafKind::open } );
    }
  fr.add( q , BBLeafKind::open );        // left open by a limit
  st.evaluated_by.assign( p , 0 );
  st.evaluated_by[ 0 ] = st.evaluated;

  // the workers, each with its share of the limits and an incumbent that
  // calls on_incumbent only when it improves on all of them
  std::vector< std::optional< SolutionT > > out( p );
  std::mutex m;
  double best = inc;                     // the best value reported
  std::exception_ptr error;
  Runtime::run( p , [ & ]( unsigned i ) {
   if( part[ i ].empty() )
    return;
   try {
    BBOptions wo = o;
    if( std::isfinite( inc ) )
     wo.incumbent = sgn * inc;
    wo.node_limit = o.node_limit ? share( rest , p , i ) : 0;
    wo.time_limit = time;
    if( o.on_incumbent )
     wo.on_incumbent = [ & ]( double v ) {
      std::lock_guard< std::mutex > guard( m );
      if( sgn * v < best ) {
       best = sgn * v;
       o.on_incumbent( v );
       }
      };
    out[ i ] = BranchAndBoundSearch< P >( f_p , wo , i ).solve(
                                                  std::move( part[ i ] ) );
    }
   catch( ... ) {
    std::lock_guard< std::mutex > guard( m );
    if( ! error )
     error = std::current_exception();
    }
   } );
  if( error )
   std::rethrow_exception( error );

  // the merge, in the order of the workers
  double low = std::min( inc , fr.low );
  fr.take( sol.frontier , sgn );
  for( unsigned i = 0 ; i < p ; ++i ) {
   if( ! out[ i ] )
    continue;
   SolutionT & s = *out[ i ];
   if( s.feasible() && ( sgn * s.value < inc ) ) {
    inc = sgn * s.value;
    inc_node = std::move( s.node );
    }
   low = std::min( low , sgn * s.bound );
   if( s.status == BBStatus::time_limit )
    limit = late = true;
   else
    if( s.status == BBStatus::node_limit )
     limit = true;
   st.add( s.stats );                    // which appends its evaluations
   st.evaluated_by.pop_back();
   st.evaluated_by[ i ] += s.stats.evaluated;
   sol.frontier.insert( sol.frontier.end() ,
                        std::make_move_iterator( s.frontier.begin() ) ,
                        std::make_move_iterator( s.frontier.end() ) );
   }

  if( limit )
   sol.status = late ? BBStatus::time_limit : BBStatus::node_limit;
  else
   sol.status = std::isfinite( inc ) ? BBStatus::optimal
                                     : BBStatus::infeasible;
  sol.node = std::move( inc_node );
  if( std::isfinite( inc ) )
   sol.value = sgn * inc;
  low = std::min( low , inc );
  sol.bound = sgn * low;
  sol.gap = std::isfinite( inc ) ? inc - low : inf;
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the evaluations of the worker i of p out of the rest of the node limit

 static std::uint64_t share( std::uint64_t rest , unsigned p , unsigned i )
  noexcept {
  return( rest / p + ( i < rest % p ? 1 : 0 ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 ParallelBBOptions f_options;                         ///< the options

 };  // end( class( SelfSplitBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SelfSplitBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------- End File SelfSplitBranchAndBound.h -------------------*/
/*--------------------------------------------------------------------------*/
