/*--------------------------------------------------------------------------*/
/*------------------------- File BranchAndBound.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the best-first search with an incumbent, i.e., for
 * Branch-and-Bound (B&B): the concept BranchAndBound of a problem, its
 * options and its Solution, the queue of the open nodes, and the
 * sequential search BranchAndBoundSearch. The parallel searches are in
 * ParallelBranchAndBound.h, and the design is discussed in
 * docs/branch-and-bound.md.
 *
 * The values (bounds, feasible values, incumbent and absolute tolerance)
 * are of the value_type of the problem, double if it gives none: any signed
 * type with a total order for which BBValueTraits says what infinity is and
 * how large the tolerance of a gap is, e.g., a floating point type, an
 * integer one (whose bounds are then integer, which prunes the nodes whose
 * bound equals the incumbent) or a type of exact arithmetic for which the
 * caller specializes BBValueTraits; a maximization problem is solved by
 * negating them internally, so that the code below only minimizes.
 *
 * The sequential search is the following, with the notation of the page
 * "Branch-and-bound as best-first search" of the documentation and the
 * conventions of the pseudocode of the page "Best-first search and A*":
 * \f$ \beta( n ) \f$ is the bound of the node \f$ n \f$, inherited from its
 * parent until \f$ n \f$ is evaluated, \f$ z \f$ the incumbent value, Q the
 * queue of the open nodes, ordered by the key of the selection, \f$ \gamma \f$
 * the fraction plunge_gap of the dive, and the leaves are the nodes that the
 * search does not branch, which give the bound on the optimum:
 *
 * \code{.unparsed}
 *  1  procedure B&B(R)          ▷ the nodes to start from, {root} by default
 *  2     z ← the given incumbent, or ∞;  Q ← R;  d ← nil  ▷ d: the dive
 *  3     while Q ≠ ∅ or d ≠ nil do
 *  4        if the node or the time limit is reached then break
 *  5        βmin ← the least bound of Q and of d
 *  6        if PRUNABLE(βmin) then the nodes of Q and d are leaves;  break
 *  7        n ← d if d ≠ nil, else the first node of Q, removed;  d ← nil
 *  8        K ← PROCESS(n)
 *  9        if the selection is plunging and K ≠ ∅ then
 * 10           k ← the node of K of least bound
 * 11           βmin ← min(the least bound of Q, β(k))
 * 12           if z = ∞ or β(k) ≤ βmin + γ (z − βmin) then d ← k, taken from K
 * 13        push the nodes of K in Q
 * 14     the nodes left in Q and d are leaves             ▷ open by a limit
 * 15     return z, its node, and the bound min(z, least bound of the leaves)
 * 16  procedure PROCESS(n)
 * 17     if PRUNABLE(β(n)) then n is a leaf (pruned);  return ∅
 * 18     cutoff(z);  (β', infeasible, v) ← evaluate(n)
 * 19     if infeasible then n is a leaf (infeasible);  return ∅
 * 20     β(n) ← max(β(n), β')
 * 21     if v is given and v < z then z ← v;  n is the node of z
 * 22     if PRUNABLE(β(n)) then n is a leaf (pruned or feasible);  return ∅
 * 23     K ← ∅;  cutoff(z)
 * 24     for each child m of n, with its estimate e, given by branch(n) do
 * 25        β(m) ← max(β(n), e)
 * 26        if PRUNABLE(β(m)) then m is a leaf (pruned)
 * 27        else K ← K ∪ {m}
 * 28     if n has no child then n is a leaf
 * 29     return K
 * 30  procedure PRUNABLE(β)
 * 31     if z = ∞ then return β = ∞
 * 32     return β ≥ z − max(τ_a, τ_r |z|)
 * \endcode
 *
 * A node is pruned when it is generated (line 26), when it is extracted
 * (line 17, by an incumbent found after its generation) and after its
 * evaluation (line 22, by its own bound), with the absolute and relative
 * tolerances \f$ \tau_a \f$ and \f$ \tau_r \f$ (abs_gap and rel_gap); every
 * element of the set below \f$ z - \max \{ \tau_a , \tau_r |z| \} \f$ stays
 * in the subset of an open node, and hence the search may stop at line 6 as
 * soon as no open node can improve \f$ z \f$ by more than the tolerance. The
 * parallel searches use the same PROCESS, each with the incumbent it sees.
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

#ifndef __PBeFS_BranchAndBound
 #define __PBeFS_BranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <functional>

#include <limits>

#include <optional>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/core/SearchUtils.h"

#include "ADS/priority/IndexedHeapPair.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the sense of the optimization
enum class Sense {
 minimize ,     ///< the least value is sought
 maximize       ///< the largest value is sought
 };

/// the order in which the open nodes are extracted
enum class Selection {
 best_first ,     ///< the least bound first
 depth_first ,    ///< the deepest node first
 breadth_first ,  ///< the shallowest node first, the oldest among those
 best_estimate ,  ///< the least estimate( n ) of the problem first
 plunging         ///< a dive into the children, then best-first
 };

/// how a B&B ended
enum class BBStatus {
 optimal ,      ///< the incumbent is optimal within the tolerance
 infeasible ,   ///< no feasible solution exists (nor was given)
 node_limit ,   ///< the limit on the evaluated nodes was reached first
 time_limit     ///< the time limit was reached first
 };

/// the outcome of the evaluation of a node, with values of type V
template< class V = double >
struct BasicEvaluation {
 V bound = V( 0 );                ///< a bound on the values of the node
 bool infeasible = false;         ///< true if the node has no solution
 std::optional< V > feasible;     ///< a feasible value, n being its
                                  ///< solution, if there is one
 };

/// the outcome of the evaluation of a node, with values of type double
using Evaluation = BasicEvaluation<>;

/// why a node is a leaf of the tree explored
enum class BBLeafKind {
 open ,         ///< left open by a limit
 pruned ,       ///< its bound cannot improve the incumbent
 feasible ,     ///< evaluated with a feasible value, and not branched
 infeasible     ///< found infeasible
 };

/// what the search looks for
/** In the mode optimize the search seeks an optimal value, pruning the
 * nodes that cannot improve the incumbent. In the mode decide it seeks a
 * feasible value at least as good as the target of the options (any
 * feasible value if none is given): it prunes the nodes whose bound is
 * worse than the target, and it stops as soon as one is found. In the mode
 * enumerate it visits the whole tree, pruning only the infeasible nodes,
 * and counts the feasible values [see BBStats::solutions]; the incumbent is
 * then the best of them, and it prunes nothing. */
//! [branch-and-bound-sig-mode]
enum class BBMode {
 optimize ,     ///< an optimal value
 decide ,       ///< a value within the target, stopping at the first one
 enumerate      ///< every node, counting the feasible ones
 };
//! [branch-and-bound-sig-mode]

/// a leaf of the tree explored, or a node to start a search from
/** The bound is in the sense of the problem: \f$ -\infty \f$
 * (\f$ +\infty \f$ for a maximization) if none is known, \f$ +\infty \f$
 * (\f$ -\infty \f$) for an infeasible leaf. A pruned or open leaf may not
 * have been evaluated, its bound being then that of its parent or its
 * estimate. */

template< class Node , class V = double >
struct BBLeaf {
 Node node;                                ///< the node
 V bound;                                  ///< its bound
 unsigned depth = 0;                       ///< its depth in the tree
 BBLeafKind kind = BBLeafKind::open;       ///< why it is a leaf
 };

/// a child of a node, with its estimate [see BBEvaluator]
//! [branch-and-bound-sig-child]
template< class Node , class V = double >
struct BBChild {
 Node node;          ///< the child
 V estimate;         ///< a bound of it that is cheap to compute
 };
//! [branch-and-bound-sig-child]

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BBValueTraits ---------------------------*/
/*--------------------------------------------------------------------------*/
/// what the B&B needs to know of the type V of its values
/** The primary template serves the signed arithmetic types: infinity() is
 * the infinity of a floating point type and the largest value of an integer
 * one (whose opposite is then the least infinite value), none() the value
 * of a Solution with no feasible value (a quiet NaN, or infinity() if V has
 * none), and tolerance( z , a , r ) is \f$ \max \{ a , r |z| \} \f$, rounded
 * down for an integer type. within( b , l , z , g ) tells whether a bound
 * \f$ b \f$ is at most \f$ l + \gamma ( z - l ) \f$, false if \f$ l \f$ is
 * not finite, computed in V for a floating point type and in long double
 * for an integer one, which cannot overflow. A caller with another type of
 * values (e.g., the rationals of an exact arithmetic) specializes the
 * class with the same members. */

template< class V >
struct BBValueTraits {

 static_assert( std::is_arithmetic_v< V > && std::is_signed_v< V > ,
                "BBValueTraits: specialize it for a type that is not a "
                "signed arithmetic one" );

 /// the infinite value
 [[nodiscard]] static constexpr V infinity( void ) noexcept {
  if constexpr( std::numeric_limits< V >::has_infinity )
   return( std::numeric_limits< V >::infinity() );
  else
   return( std::numeric_limits< V >::max() );
  }

 /// the value of a Solution with no feasible value
 [[nodiscard]] static constexpr V none( void ) noexcept {
  if constexpr( std::numeric_limits< V >::has_quiet_NaN )
   return( std::numeric_limits< V >::quiet_NaN() );
  else
   return( infinity() );
  }

 /// true if v is finite
 [[nodiscard]] static constexpr bool finite( V v ) noexcept {
  return( ( v > - infinity() ) && ( v < infinity() ) );
  }

 /// the tolerance on the gap of the finite incumbent z, by the absolute
 /// tolerance a and the relative one r
 [[nodiscard]] static V tolerance( V z , V a , double r ) noexcept {
  if constexpr( std::is_floating_point_v< V > )
   return( std::max( a , V( r * std::abs( z ) ) ) );
  else {
   using L = long double;
   return( std::max( a , V( std::floor( L( r ) * std::abs( L( z ) ) ) ) ) );
   }
  }

 /// true if b is at most l + g ( z - l ), z being finite
 [[nodiscard]] static bool within( V b , V l , V z , double g ) noexcept {
  if( ! finite( l ) )
   return( false );
  if constexpr( std::is_floating_point_v< V > )
   return( b <= l + g * ( z - l ) );
  else {
   using L = long double;
   return( L( b ) <= L( l ) + L( g ) * ( L( z ) - L( l ) ) );
   }
  }

 };  // end( struct( BBValueTraits ) )

/*--------------------------------------------------------------------------*/
/// a type of values of a B&B: signed, totally ordered, made from an int,
/// and with its BBValueTraits

template< class V >
concept BBValue = std::totally_ordered< V > && std::copyable< V > &&
 std::constructible_from< V , int > &&
 requires( const V & a , double r ) {
  { - a } -> std::convertible_to< V >;
  { a - a } -> std::convertible_to< V >;
  { BBValueTraits< V >::infinity() } -> std::convertible_to< V >;
  { BBValueTraits< V >::none() } -> std::convertible_to< V >;
  { BBValueTraits< V >::finite( a ) } -> std::convertible_to< bool >;
  { BBValueTraits< V >::tolerance( a , a , r ) } -> std::convertible_to< V >;
  { BBValueTraits< V >::within( a , a , a , r ) }
   -> std::convertible_to< bool >;
  };

namespace detail {
 /// a callable accepting a child and its estimate, used only to state the
 /// concepts
 template< class N , class V >
 struct bb_sink {
  void operator()( const N & , V ) const noexcept {}
  };

 /// the value_type of the problem P, double if it gives none
 template< class P >
 struct bb_value { using type = double; };

 template< class P > requires requires { typename P::value_type; }
 struct bb_value< P > { using type = typename P::value_type; };
 }  // end( namespace detail )

/// the type of the values of the B&B problem P: its value_type, or double
template< class P >
using bb_value_t = typename detail::bb_value< P >::type;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the children of a node, given one at a time
/** next() gives the next child of the node with its estimate, or nothing
 * when there are no more; a search that stops early (e.g., a depth-first
 * one that prunes the rest of the children) never asks for the others, so
 * that they are never built. */

//! [branch-and-bound-sig-cursor]
template< class C , class N , class V = double >
concept BBCursor = requires( C & c ) {
 { c.next() } -> std::convertible_to< std::optional< BBChild< N , V > > >;
 };
//! [branch-and-bound-sig-cursor]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// what evaluates and branches the nodes of a B&B
/** evaluate( n ) returns the Evaluation of the node n, possibly storing in
 * n what the branching needs, and either branch( n , sink ) calls
 * sink( c , e ) for every child c of the evaluated node n, e being a bound
 * of c that is cheap to compute (typically the bound of n), or
 * children( n ) returns a BBCursor that gives the children one at a time
 * (the lazy children), so that a depth-first search builds only those it
 * visits; an evaluator may give both, and the searches that take all the
 * children at once use branch( n , sink ) when it is there.
 *
 * Optionally, the evaluator may give child( n , k ), which returns the
 * child k (counted from 0) of n as an std::optional of BBChild, empty if n
 * has fewer children, without n being evaluated: the searches that rebuild
 * a node from the path of the choices that lead to it [see
 * SkeletonBranchAndBound] then follow the path by it alone, and otherwise
 * they evaluate each node of the path and take its child k; either way, the
 * children of a node must come in an order that depends on the node alone,
 * and child( n , k ) must give the k-th of them.
 *
 * Optionally, the evaluator may give cutoff( z ), which the searches call
 * before each evaluate( n ) and branch( n , sink ) with the incumbent value
 * z (in the sense of the problem, infinite when there is none, and always
 * in the modes decide and enumerate [see BBMode]), e.g., to stop the
 * relaxation of a node as soon as it is proven not to improve z, the node
 * being then found infeasible. A worker should take z from here rather
 * than from a value it shares with the others: in the deterministic search
 * z is the incumbent of the start of the round [see
 * DeterministicBranchAndBound], so that the evaluation does not depend on
 * what the other workers find meanwhile. */

template< class W , class N , class V = double >
concept BBEvaluator = requires( W & w , N & n ) {
 { w.evaluate( n ) } -> std::convertible_to< BasicEvaluation< V > >;
 } && ( requires( W & w , const N & cn ) {
         w.branch( cn , detail::bb_sink< N , V >{} );
         } ||
        requires( W & w , const N & cn ) {
         { w.children( cn ) } -> BBCursor< N , V >;
         } );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an evaluator that gives the lazy children [see BBEvaluator]

//! [branch-and-bound-sig-lazy]
template< class W , class N , class V = double >
concept LazyBBEvaluator = BBEvaluator< W , N , V > &&
 requires( W & w , const N & cn ) {
  { w.children( cn ) } -> BBCursor< N , V >;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an evaluator that gives the child k of a node [see BBEvaluator]

template< class W , class N , class V = double >
concept IndexedBBEvaluator = BBEvaluator< W , N , V > &&
 requires( W & w , const N & cn , std::size_t k ) {
  { w.child( cn , k ) }
   -> std::convertible_to< std::optional< BBChild< N , V > > >;
  };
//! [branch-and-bound-sig-lazy]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a B&B problem
/** A type P is a BranchAndBound problem if it provides the copyable type
 * node_type, the constant sense, the method root() giving the node of the
 * whole set, optionally the type value_type of its values [see BBValue],
 * double if it gives none, and either the methods evaluate( n ) and
 * branch( n , sink ) of a BBEvaluator, which have to be safe to call
 * concurrently in a parallel search, or make_worker( i ), which returns the
 * BBEvaluator of the worker i (e.g., one with its own copy of a model and
 * of a solver).
 *
 * Optionally, P may give hash( n ), which the distributed search uses to
 * send n to its owner, and estimate( n ), the estimate of the best
 * feasible value under n that orders the nodes with
 * Selection::best_estimate. */

//! [branch-and-bound-sig-concept]
template< class P >
concept BranchAndBound = requires( const P & p , unsigned i ) {
 typename P::node_type;
 requires std::copyable< typename P::node_type >;
 requires BBValue< bb_value_t< P > >;
 { P::sense } -> std::convertible_to< Sense >;
 { p.root() } -> std::convertible_to< typename P::node_type >;
 } && ( BBEvaluator< const P , typename P::node_type , bb_value_t< P > > ||
        requires( const P & p , unsigned i ) {
         { p.make_worker( i ) }
          -> BBEvaluator< typename P::node_type , bb_value_t< P > >;
         } );
//! [branch-and-bound-sig-concept]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a B&B problem that gives a hash of its nodes

template< class P >
concept HashedBranchAndBound = BranchAndBound< P > &&
 requires( const P & p , const typename P::node_type & n ) {
  { p.hash( n ) } -> std::convertible_to< std::uint64_t >;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a B&B problem that estimates the best feasible value under a node

template< class P >
concept EstimatingBranchAndBound = BranchAndBound< P > &&
 requires( const P & p , const typename P::node_type & n ) {
  { p.estimate( n ) } -> std::convertible_to< bb_value_t< P > >;
  };

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BBOptions -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of a B&B with values of type V
/** A node is pruned when its bound cannot improve the incumbent \f$ z \f$
 * by more than the tolerance \f$ \max \{ \tau_a , \tau_r |z| \} \f$,
 * \f$ \tau_a \f$ being abs_gap and \f$ \tau_r \f$ rel_gap, and the search
 * ends when no open node can. With Selection::plunging, the best child of
 * the node just branched is extracted next (a dive) if its bound is at most
 * \f$ \underline{\beta} + \gamma ( z - \underline{\beta} ) \f$,
 * \f$ \gamma \f$ being plunge_gap and \f$ \underline{\beta} \f$ the least
 * open bound, and always while there is no incumbent. With keep_frontier,
 * the Solution holds the leaves of the tree explored [see BBSolution], from
 * which a later search can start again [see BranchAndBoundSearch::solve()],
 * and with keep_open_only too it holds only those left open by a limit,
 * which is what a checkpoint needs [see BBCheckpoint.h]. With prune_level,
 * the children of a node are taken to come in the order of their estimates,
 * the best first, so that once a child is pruned when it is generated so
 * are the children after it, which the depth-first searches with the lazy
 * children then never build [see BBEvaluator]. The mode says what the
 * search looks for [see BBMode], and target is the value of the mode
 * decide. */

//! [branch-and-bound-sig-options]
template< class V = double >
struct BasicBBOptions {

 Selection selection = Selection::best_first;   ///< the order of the nodes
 BBMode mode = BBMode::optimize;   ///< what the search looks for
 std::optional< V > target;     ///< the value sought by the mode decide
 V abs_gap = V( 0 );           ///< the absolute tolerance on the gap
 double rel_gap = 0;           ///< the relative tolerance on the gap
 double plunge_gap = 0.5;      ///< how far from the best bound a dive goes
 std::optional< V > incumbent;  ///< a known feasible value, if any
 /// a value that the search reads at each extraction and prunes by when it
 /// is better than its incumbent, e.g., the incumbent of other searches
 /// running at the same time; the Solution has only the feasible values
 /// that the search itself finds (or incumbent), so that it is infeasible
 /// when none of them is better than the values read, and its bound is
 /// still that of the subtrees searched
 std::function< V( void ) > shared_incumbent;
 std::size_t node_limit = 0;   ///< evaluated nodes (0 = no limit)
 double time_limit = 0;        ///< seconds (0 = no limit)
 bool keep_frontier = false;   ///< keep the leaves of the tree explored
 bool keep_open_only = false;  ///< ... only those left open by a limit
 bool prune_level = false;     ///< the children come by their estimates
 /// called with every new incumbent value (in a parallel search, by the
 /// worker that finds it, under a lock)
 std::function< void( V ) > on_incumbent;

 };  // end( struct( BasicBBOptions ) )
//! [branch-and-bound-sig-options]

/// the options of a B&B with values of type double
using BBOptions = BasicBBOptions<>;

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS BBStats ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the counters of a B&B
/** A node is pruned when it is generated, when it is extracted (by the
 * incumbent found after its generation) or after its evaluation (by its
 * own bound); the vector has one entry per worker. A feasible value
 * counts among solutions whenever an evaluation gives one, which in the
 * mode enumerate is the number of the feasible nodes of the tree, and
 * recomputed counts the evaluations made only to rebuild a node from the
 * path that leads to it [see SkeletonBranchAndBound], which evaluated
 * does not count; sent and messages count the nodes and the batches that
 * travel between the workers, or the steals of the searches by stealing. */

struct BBStats {

 std::uint64_t evaluated = 0;         ///< nodes evaluated
 std::uint64_t branched = 0;          ///< nodes branched
 std::uint64_t generated = 0;         ///< children generated
 std::uint64_t infeasible = 0;        ///< nodes found infeasible
 std::uint64_t pruned_generated = 0;  ///< children pruned when generated
 std::uint64_t pruned_extracted = 0;  ///< nodes pruned when extracted
 std::uint64_t pruned_evaluated = 0;  ///< nodes pruned by their bound
 std::uint64_t incumbents = 0;        ///< improvements of the incumbent
 std::uint64_t solutions = 0;         ///< feasible values found
 std::uint64_t recomputed = 0;        ///< evaluations that rebuild nodes
 std::uint64_t sent = 0;              ///< nodes sent to another worker
 std::uint64_t messages = 0;          ///< batches of nodes sent
 std::size_t max_open = 0;            ///< the most open nodes of a queue
 double seconds = 0;                  ///< wall-clock time of the search
 std::vector< std::uint64_t > evaluated_by;  ///< evaluations per worker

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the counters of a worker

 void add( const BBStats & o ) {
  evaluated += o.evaluated;
  branched += o.branched;
  generated += o.generated;
  infeasible += o.infeasible;
  pruned_generated += o.pruned_generated;
  pruned_extracted += o.pruned_extracted;
  pruned_evaluated += o.pruned_evaluated;
  incumbents += o.incumbents;
  solutions += o.solutions;
  recomputed += o.recomputed;
  sent += o.sent;
  messages += o.messages;
  max_open = std::max( max_open , o.max_open );
  evaluated_by.push_back( o.evaluated );
  }

 };  // end( struct( BBStats ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BBSolution ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the outcome of a B&B
/** value is the incumbent value and bound a bound on the optimum, both in
 * the sense of the problem, and the gap, value minus bound in absolute
 * value, is valid whatever
 * the status; value is none() of BBValueTraits if no feasible value is
 * known, and node is the incumbent node, empty if no feasible node was
 * found (even when an initial incumbent value was given).
 *
 * With BBOptions::keep_frontier, frontier holds the leaves of the tree
 * explored, in no particular order: the nodes pruned, feasible or
 * infeasible, and those left open by a limit. Their subtrees cover the
 * subtrees of the nodes the search started from, so that a search started
 * from them [see BranchAndBoundSearch::solve()] explores the same region,
 * e.g., after a change of the data that leaves the branching valid. With
 * BBOptions::keep_open_only, frontier holds only the leaves left open, and
 * closed_bound is the least bound of the others (pruned by the tolerance),
 * so that the bound of the region they cover is known without them.
 *
 * In the mode decide [see BBMode] the status is optimal if a value within
 * the target was found, value being then that value, and infeasible
 * otherwise. */

template< class Node , class V = double >
struct BBSolution {

 BBStatus status = BBStatus::infeasible;   ///< how the search ended
 std::optional< Node > node;               ///< the incumbent node
 V value = BBValueTraits< V >::none();     ///< its value
 V bound = BBValueTraits< V >::none();     ///< the bound
 V gap = BBValueTraits< V >::infinity();   ///< the gap
 /// the least bound of the leaves that are not open
 V closed_bound = BBValueTraits< V >::none();
 BBStats stats;                            ///< the counters
 std::vector< BBLeaf< Node , V > > frontier;  ///< the leaves, if kept

 /// true if a feasible value is known
 [[nodiscard]] bool feasible( void ) const noexcept {
  return( BBValueTraits< V >::finite( value ) );
  }

 };  // end( struct( BBSolution ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS BBQueue ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the open nodes of a B&B, ordered by a key
/** The nodes are kept in slots with a free list, since a node is destroyed
 * when it is branched, and the indexed heap pair of ADS on the slots gives
 * the next one, the least key first and, among equal keys, the last pushed
 * (the first pushed with Selection::breadth_first), and the least bound of
 * the open nodes, which the termination test needs when the key of the
 * selection is not the bound. The heaps grow by doubling with the slots. */

template< class Node , class V = double >
class BBQueue {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// an open node
 struct Item {
  Node node;          ///< the node
  V bound;            ///< its bound, in the minimization
  unsigned depth;     ///< its depth in the tree
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, taking the selection the keys are made by

 explicit BBQueue( Selection s = Selection::best_first )
  : f_fifo( s == Selection::breadth_first ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds a node with its bound, its depth and its key

 void push( Node && n , V bound , unsigned depth , V key ) {
  std::size_t i;
  if( v_free.empty() ) {
   i = v_slot.size();
   v_slot.emplace_back();
   if( i >= std::size_t( f_q.ids() ) )
    grow( 2 * i + 16 );
   }
  else {
   i = v_free.back();
   v_free.pop_back();
   }
  const std::uint64_t seq = ++f_seq;
  v_slot[ i ].emplace( Item{ std::move( n ) , bound , depth } );
  f_q.push( int( i ) , Order( key , f_fifo ? seq : ~seq ) , bound );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the next node, which stays in the queue (there must be one)

 [[nodiscard]] const Item & top( void ) const {
  return( *v_slot[ std::size_t( f_q.top() ) ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes and returns the next node

 Item pop( void ) {
  const auto i = std::size_t( f_q.pop() );
  Item it = std::move( *v_slot[ i ] );
  v_slot[ i ].reset();
  v_free.push_back( i );
  return( it );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the least bound of the open nodes (\f$ +\infty \f$ if there is none)

 [[nodiscard]] V best_bound( void ) const {
  return( f_q.empty() ? BBValueTraits< V >::infinity() : f_q.minBound() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes all the nodes, returning how many there were

 std::size_t clear( void ) {
  const std::size_t n = f_q.size();
  v_slot.clear();
  v_free.clear();
  f_q.clear();
  return( n );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes all the nodes, giving each to f( item ), and returns how many
 /// there were

 template< class F >
 std::size_t drain( F && f ) {
  for( auto & x : v_slot )
   if( x )
    f( std::move( *x ) );
  return( clear() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( item ) for every open node, in no particular order

 template< class F >
 void for_each( F && f ) const {
  for( const auto & x : v_slot )
   if( x )
    f( *x );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if there is no open node

 [[nodiscard]] bool empty( void ) const noexcept { return( f_q.empty() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of open nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_q.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the order of the selection: the key, then the push, reversed (the
 /// last pushed first) unless first in, first out
 using Order = std::pair< V , std::uint64_t >;

 /// the heaps of the keys and of the bounds, on the slots
 using Queue = ADS::IndexedHeapPair< Order , V >;

 // the heaps on n slots, with the nodes they hold
 void grow( std::size_t n ) {
  Queue q{ int( n ) };
  for( std::size_t j = 0 ; j < v_slot.size() ; ++j )
   if( v_slot[ j ] && f_q.contains( int( j ) ) )
    q.push( int( j ) , f_q.key( int( j ) ) , f_q.bound( int( j ) ) );
  f_q = std::move( q );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool f_fifo;                          ///< first in, first out on ties
 std::uint64_t f_seq = 0;              ///< the pushes so far
 std::vector< std::optional< Item > > v_slot;   ///< the slots
 std::vector< std::size_t > v_free;    ///< the free slots
 Queue f_q;                            ///< the heaps of the open nodes

 };  // end( class( BBQueue ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// the evaluator of a problem without make_worker(): the problem itself

template< class P >
struct BBSelf {
 using Node = typename P::node_type;
 const P & p;
 BasicEvaluation< bb_value_t< P > > evaluate( Node & n ) {
  return( p.evaluate( n ) );
  }
 template< class F >
 void branch( const Node & n , F && f )
  requires requires( const P & q ) { q.branch( n , f ); } {
  p.branch( n , std::forward< F >( f ) );
  }
 auto children( const Node & n )
  requires requires( const P & q ) { q.children( n ); } {
  return( p.children( n ) );
  }
 auto child( const Node & n , std::size_t k )
  requires requires( const P & q ) { q.child( n , k ); } {
  return( p.child( n , k ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// adds the counters of b to those of a, but for max_open, seconds and
/// evaluated_by

inline void bb_add_counters( BBStats & a , const BBStats & b )
{
 a.evaluated += b.evaluated;
 a.branched += b.branched;
 a.generated += b.generated;
 a.infeasible += b.infeasible;
 a.pruned_generated += b.pruned_generated;
 a.pruned_extracted += b.pruned_extracted;
 a.pruned_evaluated += b.pruned_evaluated;
 a.incumbents += b.incumbents;
 a.solutions += b.solutions;
 a.recomputed += b.recomputed;
 a.sent += b.sent;
 a.messages += b.messages;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the evaluator of the worker i: make_worker( i ), or the problem itself

template< BranchAndBound P >
auto bb_worker( const P & p , unsigned i ) {
 if constexpr( requires { p.make_worker( i ); } )
  return( p.make_worker( i ) );
 else {
  ( void ) i;
  return( BBSelf< P >{ p } );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the value v turned from the sense of P to the minimization, or back: v
/// itself for a minimization, - v for a maximization

template< BranchAndBound P , class V >
constexpr V bb_signed( const V & v ) {
 if( P::sense == Sense::minimize )
  return( v );
 return( - v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// true if the incumbent inc (in the minimization) is what the mode decide
/// seeks: a feasible value, within the target if one is given

template< BranchAndBound P , class V >
bool bb_met( const V & inc , const BasicBBOptions< V > & o )
{
 return( BBValueTraits< V >::finite( inc ) &&
         ( ( ! o.target ) || ( inc <= bb_signed< P >( *o.target ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// true if a node of the given bound cannot improve the incumbent inc (both
/// in the minimization) by more than the tolerance; in the mode enumerate,
/// if the bound is infinite, and in the mode decide, if it is worse than
/// the target or the target is met [see BBMode]

template< BranchAndBound P , class V >
bool bb_prunable( const V & bound , const V & inc ,
                  const BasicBBOptions< V > & o )
{
 using Tr = BBValueTraits< V >;
 if( o.mode == BBMode::enumerate )
  return( bound == Tr::infinity() );
 if( o.mode == BBMode::decide )
  return( bb_met< P >( inc , o ) || ( bound == Tr::infinity() ) ||
          ( o.target && ( bound > bb_signed< P >( *o.target ) ) ) );
 if( ! Tr::finite( inc ) )
  return( bound == Tr::infinity() );
 return( bound >= inc - Tr::tolerance( inc , o.abs_gap , o.rel_gap ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the status of a search that ended without a limit with the incumbent z
/// (in the minimization)

template< BranchAndBound P , class V >
BBStatus bb_status( const V & z , const BasicBBOptions< V > & o )
{
 if( o.mode == BBMode::decide )
  return( bb_met< P >( z , o ) ? BBStatus::optimal : BBStatus::infeasible );
 return( BBValueTraits< V >::finite( z ) ? BBStatus::optimal
                                         : BBStatus::infeasible );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the leaves of the part of the tree explored by a worker
/** low is the least bound (in the minimization) of the nodes pruned by the
 * tolerance or left open, so that \f$ \min \{ \mathrm{low} , z \} \f$,
 * \f$ z \f$ being the incumbent, is a bound on the optimum, and closed the
 * least of those that are not open; the leaves themselves are kept only if
 * keep, and only the open ones if open_only. */

template< class Node , class V = double >
struct BBFrontier {

 explicit BBFrontier( bool k = false , bool o = false )
  : keep( k ) , open_only( o ) {}

 explicit BBFrontier( const BasicBBOptions< V > & o )
  : keep( o.keep_frontier ) , open_only( o.keep_open_only ) {}

 /// true if a leaf of the given kind is kept
 [[nodiscard]] bool kept( BBLeafKind kind ) const noexcept {
  return( keep && ( ( ! open_only ) || ( kind == BBLeafKind::open ) ) );
  }

 /// lowers low, and closed unless kind is open, by the bound b of a leaf
 void note( V b , BBLeafKind kind ) {
  low = std::min( low , b );
  if( kind != BBLeafKind::open )
   closed = std::min( closed , b );
  }

 /// adds the leaf n of bound b (in the minimization), lowering low unless
 /// the subtree of n holds no better solution than its own
 void add( Node && n , V b , unsigned depth , BBLeafKind kind ,
           bool lower = true ) {
  if( lower )
   note( b , kind );
  if( kept( kind ) )
   leaves.push_back( BBLeaf< Node , V >{ std::move( n ) , b , depth ,
                                         kind } );
  }

 /// adds all the nodes of q as leaves of the given kind, returning their
 /// number
 std::size_t add( BBQueue< Node , V > & q , BBLeafKind kind ) {
  if( ! q.empty() )
   note( q.best_bound() , kind );
  if( ! kept( kind ) )
   return( q.clear() );
  return( q.drain( [ & ]( auto && it ) {
   leaves.push_back( BBLeaf< Node , V >{ std::move( it.node ) , it.bound ,
                                     it.depth , kind } );
   } ) );
  }

 /// appends the leaves of another worker
 void merge( BBFrontier && o ) {
  low = std::min( low , o.low );
  closed = std::min( closed , o.closed );
  leaves.insert( leaves.end() , std::make_move_iterator( o.leaves.begin() ) ,
                 std::make_move_iterator( o.leaves.end() ) );
  o.leaves.clear();
  }

 /// moves the leaves, with their bounds turned to the sense of P, into v
 template< class P >
 void take( std::vector< BBLeaf< Node , V > > & v ) {
  for( auto & l : leaves )
   l.bound = bb_signed< P >( l.bound );
  v = std::move( leaves );
  leaves.clear();
  }

 bool keep;                                     ///< true if leaves is kept
 bool open_only;                                ///< ... only the open ones
 V low = BBValueTraits< V >::infinity();        ///< the least bound
 V closed = BBValueTraits< V >::infinity();     ///< ... of those not open
 std::vector< BBLeaf< Node , V > > leaves;      ///< the leaves

 };  // end( struct( BBFrontier ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the nodes a search starts from: the given leaves, with their bounds
/// turned to the minimization, or the root if none given

template< BranchAndBound P , class Seeds >
std::vector< BBLeaf< typename P::node_type , bb_value_t< P > > > bb_starts(
                                                   const P & p ,
                                                          Seeds && seeds )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 std::vector< BBLeaf< Node , V > > v;
 if( ! seeds )
  v.push_back( BBLeaf< Node , V >{ p.root() ,
                                   - BBValueTraits< V >::infinity() ,
                               0 , BBLeafKind::open } );
 else {
  v = std::move( *seeds );
  for( auto & l : v )
   l.bound = bb_signed< P >( l.bound );
  }
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the key of a child in the queue, by the selection

template< BranchAndBound P >
bb_value_t< P > bb_key( const P & p , const typename P::node_type & n ,
                       const bb_value_t< P > & bound ,
               unsigned depth , Selection s )
{
 switch( s ) {
  case( Selection::depth_first ):
   return( - bb_value_t< P >( int( depth ) ) );
  case( Selection::breadth_first ):
   return( bb_value_t< P >( int( depth ) ) );
  case( Selection::best_estimate ):
   if constexpr( EstimatingBranchAndBound< P > )
    return( bb_signed< P >( bb_value_t< P >( p.estimate( n ) ) ) );
   else
    return( bound );
  default:
   return( bound );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// calls cutoff( z ) of the evaluator w, if it has one, with the incumbent
/// inc (in the minimization) turned to the sense of P, or with infinity in
/// the modes decide and enumerate

template< BranchAndBound P , class W , class V >
void bb_cutoff( W & w , const V & inc , const BasicBBOptions< V > & o )
{
 if constexpr( requires { w.cutoff( inc ); } )
  w.cutoff( bb_signed< P >( o.mode == BBMode::optimize ?
                            inc : BBValueTraits< V >::infinity() ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a cursor on the children that an evaluator gives all at once by
/// branch( n , sink ), for the searches that take them one at a time

template< class Node , class V >
struct BBEagerCursor {

 std::optional< BBChild< Node , V > > next( void ) {
  if( at == kids.size() )
   return( std::nullopt );
  return( std::move( kids[ at++ ] ) );
  }

 std::vector< BBChild< Node , V > > kids;       ///< the children
 std::size_t at = 0;                            ///< the next one
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// gives each child of the evaluated node n, with its estimate (in the
/// sense of P), to f( child , estimate ): by branch( n , f ) if the
/// evaluator has it, by the cursor of children( n ) otherwise

template< BranchAndBound P , class W , class F >
void bb_branch( W & w , const typename P::node_type & n , F && f )
{
 if constexpr( requires { w.branch( n , f ); } )
  w.branch( n , f );
 else {
  auto c = w.children( n );
  while( auto x = c.next() )
   f( std::move( x->node ) , bb_value_t< P >( x->estimate ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the cursor of the children of the evaluated node n: children( n ) if
/// the evaluator has it, else one on all those given by branch( n , sink )

template< BranchAndBound P , class W >
auto bb_children( W & w , const typename P::node_type & n )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 if constexpr( LazyBBEvaluator< W , Node , V > )
  return( w.children( n ) );
 else {
  BBEagerCursor< Node , V > c;
  w.branch( n , [ & ]( auto && x , V e ) {
   c.kids.push_back( BBChild< Node , V >{
                      Node( std::forward< decltype( x ) >( x ) ) , e } );
   } );
  return( c );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// prunes an extracted node by the incumbent, evaluates it, updates the
/// incumbent and prunes it by its bound; returns the bound of n (in the
/// minimization) if it is to be branched, nothing if it is a leaf
/** inc is the incumbent value (in the minimization), improve( v , n ) offers
 * the feasible value v of the node n, and n is moved to fr if it is a leaf
 * [see BBFrontier]; kind is set to the kind of leaf that n is if it has no
 * child, feasible if it gave a feasible value and pruned otherwise. */

template< BranchAndBound P , class W , class Inc , class Improve >
std::optional< bb_value_t< P > > bb_evaluate( typename P::node_type & n ,
                                              bb_value_t< P > bound ,
                                              unsigned depth , W & w ,
                                              Inc && inc ,
                                              Improve && improve ,
                    const BasicBBOptions< bb_value_t< P > > & o ,
                                              BBStats & st ,
                BBFrontier< typename P::node_type , bb_value_t< P > > & fr ,
                                              BBLeafKind & kind )
{
 using V = bb_value_t< P >;
 const V inf = BBValueTraits< V >::infinity();
 if( bb_prunable< P >( bound , inc() , o ) ) {
  ++st.pruned_extracted;
  fr.add( std::move( n ) , bound , depth , BBLeafKind::pruned );
  return( std::nullopt );
  }
 bb_cutoff< P >( w , inc() , o );
 const BasicEvaluation< V > e = w.evaluate( n );
 ++st.evaluated;
 if( e.infeasible ) {
  ++st.infeasible;
  fr.add( std::move( n ) , inf , depth , BBLeafKind::infeasible , false );
  return( std::nullopt );
  }
 const V b = std::max( bound , bb_signed< P >( e.bound ) );
 if( e.feasible ) {
  ++st.solutions;
  if( improve( bb_signed< P >( *e.feasible ) , n ) )
   ++st.incumbents;
  }
 kind = e.feasible ? BBLeafKind::feasible : BBLeafKind::pruned;
 if( bb_prunable< P >( b , inc() , o ) ) {
  ++st.pruned_evaluated;
  fr.add( std::move( n ) , b , depth , kind );
  return( std::nullopt );
  }
 return( b );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the bound (in the minimization) of a child of estimate est (in the
/// sense of P) of a node of bound b, if the child survives, nothing if it
/// is pruned, its bound lowering then the bounds of fr unless fr keeps it,
/// which is left to the caller; cut is set when the child is pruned with
/// prune_level, and then all the children after it are pruned too

template< BranchAndBound P >
std::optional< bb_value_t< P > > bb_child_bound( bb_value_t< P > b ,
                                                 bb_value_t< P > est ,
                                                 const bb_value_t< P > & inc ,
                    const BasicBBOptions< bb_value_t< P > > & o ,
                                                 bool & cut , BBStats & st ,
                BBFrontier< typename P::node_type , bb_value_t< P > > & fr )
{
 ++st.generated;
 const bb_value_t< P > cb = std::max( b , bb_signed< P >( est ) );
 if( cut || bb_prunable< P >( cb , inc , o ) ) {
  ++st.pruned_generated;
  cut = cut || o.prune_level;
  if( ! fr.kept( BBLeafKind::pruned ) )
   fr.note( cb , BBLeafKind::pruned );
  return( std::nullopt );
  }
 return( cb );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// processes an extracted node: prunes it by the incumbent, evaluates it,
/// updates the incumbent, prunes it by its bound, and branches it, giving
/// each child that survives to push( child , bound , depth )
/** inc is the incumbent value (in the minimization), improve( v , n ) offers
 * the feasible value v of the node n, and the leaves (n itself, or its
 * children pruned when generated) go to fr [see BBFrontier], n being moved
 * there if it is one. */

template< BranchAndBound P , class W , class Inc , class Improve ,
          class Push >
void bb_process( typename P::node_type & n , bb_value_t< P > bound ,
                 unsigned depth ,
                 W & w , Inc && inc , Improve && improve , Push && push ,
                 const BasicBBOptions< bb_value_t< P > > & o , BBStats & st ,
                 BBFrontier< typename P::node_type , bb_value_t< P > > & fr )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 BBLeafKind kind = BBLeafKind::pruned;
 const std::optional< V > ob = bb_evaluate< P >( n , bound , depth , w , inc ,
                                                 improve , o , st , fr ,
                                                 kind );
 if( ! ob )
  return;
 const V b = *ob;
 ++st.branched;
 const std::uint64_t g = st.generated;
 bool cut = false;
 bb_cutoff< P >( w , inc() , o );
 bb_branch< P >( w , std::as_const( n ) , [ & ]( auto && c , V est ) {
  const std::optional< V > cb = bb_child_bound< P >( b , est , inc() , o ,
                                                     cut , st , fr );
  if( ! cb ) {
   if( fr.kept( BBLeafKind::pruned ) )
    fr.add( Node( std::forward< decltype( c ) >( c ) ) ,
            std::max( b , bb_signed< P >( est ) ) , depth + 1 ,
            BBLeafKind::pruned );
   return;
   }
  push( Node( std::forward< decltype( c ) >( c ) ) , *cb , depth + 1 );
  } );
 if( st.generated == g )                // no children: its subtree is n
  fr.add( std::move( n ) , b , depth , kind , false );
 }

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BranchAndBoundSearch -----------------------*/
/*--------------------------------------------------------------------------*/
/// the sequential best-first search with an incumbent
/** The nodes are extracted by the Selection of the options and evaluated
 * only when they are extracted, a child inheriting until then the bound of
 * its parent or its own estimate, whichever is larger; a node is pruned
 * when it is generated, when it is extracted and after its evaluation, and
 * the search ends when no open node can improve the incumbent by more than
 * the tolerance [see docs/branch-and-bound.md]. */

template< BranchAndBound P >
class BranchAndBoundSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicBBOptions< Value >;       ///< the options type
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search) and
 /// the worker whose evaluator the search uses [see make_worker()], so that
 /// sequential searches can run at the same time on different workers

 explicit BranchAndBoundSearch( const P & problem ,
                                const Options & options = {} ,
                                unsigned worker = 0 )
  : f_p( problem ) , f_options( options ) , f_worker( worker ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution
 /** The nodes (the kind of a BBLeaf is not looked at) are open nodes with
  * their bounds and depths, evaluated only when extracted, e.g., the
  * leaves of a previous search [see BBSolution::frontier] whose bounds are
  * still valid, or set to the unknown bound where they are not; with none,
  * the search ends at once. The bound of the Solution is that of their
  * subtrees, so that the caller combines it with the bounds of whatever
  * nodes it leaves out. */

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  using Tr = BBValueTraits< Value >;
  const Value inf = Tr::infinity();
  const Options & o = f_options;
  const Selection sel = o.selection;
  Stopwatch clock;
  SolutionT sol;
  BBStats & st = sol.stats;
  auto w = detail::bb_worker( f_p , f_worker );

  Value inc = o.incumbent ? detail::bb_signed< P >( *o.incumbent ) : inf;
  Value cut = inc;          // what prunes: inc, or shared_incumbent if less
  std::optional< Node > inc_node;
  const auto improve = [ & ]( Value v , const Node & n ) {
   if( v >= cut )
    return( false );
   inc = cut = v;
   inc_node = n;
   if( o.on_incumbent )
    o.on_incumbent( detail::bb_signed< P >( v ) );
   return( true );
   };

  using Item = typename BBQueue< Node , Value >::Item;
  BBQueue< Node , Value > q( sel );
  std::optional< Item > dive;              // the node of the dive, if any
  std::vector< Item > kids;                // the children of a node
  detail::BBFrontier< Node , Value > fr( o );  // the leaves

  const auto push = [ & ]( Item && it ) {
   const Value key = detail::bb_key( f_p , it.node , it.bound , it.depth ,
                                      sel );
   q.push( std::move( it.node ) , it.bound , it.depth , key );
   st.max_open = std::max( st.max_open , q.size() );
   };

  for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) )
   push( Item{ std::move( l.node ) , l.bound , l.depth } );
  sol.status = BBStatus::optimal;

  for( ;; ) {
   if( o.shared_incumbent )
    cut = std::min( cut , detail::bb_signed< P >( o.shared_incumbent() ) );
   if( q.empty() && ! dive )          // the search is over, at a limit too
    break;
   if( o.node_limit && ( st.evaluated >= o.node_limit ) ) {
    sol.status = BBStatus::node_limit;
    break;
    }
   if( ( o.time_limit > 0 ) && ( clock.seconds() >= o.time_limit ) ) {
    sol.status = BBStatus::time_limit;
    break;
    }
   const Value best = std::min( q.best_bound() ,
                                 dive ? dive->bound : inf );
   if( detail::bb_prunable< P >( best , cut , o ) ) {   // nothing can improve
    st.pruned_extracted += fr.add( q , BBLeafKind::pruned );
    if( dive ) {
     ++st.pruned_extracted;
     fr.add( std::move( dive->node ) , dive->bound , dive->depth ,
             BBLeafKind::pruned );
     dive.reset();
     }
    break;
    }

   Item it = dive ? std::move( *dive ) : q.pop();
   dive.reset();
   kids.clear();
   detail::bb_process< P >( it.node , it.bound , it.depth , w ,
                            [ & ] { return( cut ); } , improve ,
                            [ & ]( Node && c , Value b , unsigned d ) {
                             kids.push_back( Item{ std::move( c ) , b , d } );
                             } , o , st , fr );
   if( kids.empty() )
    continue;

   // with plunging, the best child dives if it is close to the best bound
   std::size_t k = kids.size();
   if( sel == Selection::plunging ) {
    std::size_t c = 0;
    for( std::size_t j = 1 ; j < kids.size() ; ++j )
     if( kids[ j ].bound <= kids[ c ].bound )
      c = j;
    const Value b0 = std::min( q.best_bound() , kids[ c ].bound );
    if( ! Tr::finite( cut ) ||
        Tr::within( kids[ c ].bound , b0 , cut , o.plunge_gap ) )
     k = c;
    }
   for( std::size_t j = 0 ; j < kids.size() ; ++j )
    if( j == k )
     dive.emplace( std::move( kids[ j ] ) );
    else
     push( std::move( kids[ j ] ) );
   }

  // the open nodes left by a limit bound the optimum too
  fr.add( q , BBLeafKind::open );
  if( dive )
   fr.add( std::move( dive->node ) , dive->bound , dive->depth ,
           BBLeafKind::open );
  if( sol.status == BBStatus::optimal )
   sol.status = detail::bb_status< P >( inc , o );
  sol.node = std::move( inc_node );
  if( Tr::finite( inc ) )
   sol.value = detail::bb_signed< P >( inc );
  const Value lb = std::min( inc , fr.low );
  sol.bound = detail::bb_signed< P >( lb );
  sol.gap = Tr::finite( inc ) && Tr::finite( lb ) ? Value( inc - lb ) : inf;
  sol.closed_bound = detail::bb_signed< P >( fr.closed );
  fr.template take< P >( sol.frontier );
  st.evaluated_by = { st.evaluated };
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                      ///< the problem
 Options f_options;                                   ///< the options
 unsigned f_worker;                                   ///< the worker

 };  // end( class( BranchAndBoundSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the sequential B&B of a problem, with the given options

template< BranchAndBound P >
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
 branch_and_bound( const P & problem ,
                   const BasicBBOptions< bb_value_t< P > > & o = {} )
{
 return( BranchAndBoundSearch< P >( problem , o ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the sequential B&B of the subtrees of the given nodes [see
/// BranchAndBoundSearch::solve()], with the given options

template< BranchAndBound P >
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
 branch_and_bound( const P & problem ,
                   std::vector< BBLeaf< typename P::node_type ,
                                        bb_value_t< P > > > seeds ,
                   const BasicBBOptions< bb_value_t< P > > & o = {} )
{
 return( BranchAndBoundSearch< P >( problem , o ).solve(
                                                       std::move( seeds ) ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BranchAndBound */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BranchAndBound.h ------------------------*/
/*--------------------------------------------------------------------------*/
