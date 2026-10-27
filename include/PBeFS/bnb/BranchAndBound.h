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
 *         Universita' di Pisa \n
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
/// what evaluates and branches the nodes of a B&B
/** evaluate( n ) returns the Evaluation of the node n, possibly storing in
 * n what the branching needs, and branch( n , sink ) calls sink( c , e )
 * for every child c of the evaluated node n, e being a bound of c that is
 * cheap to compute (typically the bound of n).
 *
 * Optionally, the evaluator may give cutoff( z ), which the searches call
 * before each evaluate( n ) and branch( n , sink ) with the incumbent value
 * z (in the sense of the problem, infinite when there is none), e.g., to
 * stop the relaxation of a node as soon as it is proven not to improve z,
 * the node being then found infeasible. A worker should take z from here
 * rather than from a value it shares with the others: in the deterministic
 * search z is the incumbent of the start of the round [see
 * DeterministicBranchAndBound], so that the evaluation does not depend on
 * what the other workers find meanwhile. */

template< class W , class N , class V = double >
concept BBEvaluator = requires( W & w , N & n , const N & cn ) {
 { w.evaluate( n ) } -> std::convertible_to< BasicEvaluation< V > >;
 w.branch( cn , detail::bb_sink< N , V >{} );
 };

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
 * which a later search can start again [see BranchAndBoundSearch::solve()]. */

template< class V = double >
struct BasicBBOptions {

 Selection selection = Selection::best_first;   ///< the order of the nodes
 V abs_gap = V( 0 );           ///< the absolute tolerance on the gap
 double rel_gap = 0;           ///< the relative tolerance on the gap
 double plunge_gap = 0.5;      ///< how far from the best bound a dive goes
 std::optional< V > incumbent;  ///< a known feasible value, if any
 std::size_t node_limit = 0;   ///< evaluated nodes (0 = no limit)
 double time_limit = 0;        ///< seconds (0 = no limit)
 bool keep_frontier = false;   ///< keep the leaves of the tree explored
 /// called with every new incumbent value (in a parallel search, by the
 /// worker that finds it, under a lock)
 std::function< void( V ) > on_incumbent;

 };  // end( struct( BasicBBOptions ) )

/// the options of a B&B with values of type double
using BBOptions = BasicBBOptions<>;

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS BBStats ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the counters of a B&B
/** A node is pruned when it is generated, when it is extracted (by the
 * incumbent found after its generation) or after its evaluation (by its
 * own bound); the vector has one entry per worker. */

struct BBStats {

 std::uint64_t evaluated = 0;         ///< nodes evaluated
 std::uint64_t branched = 0;          ///< nodes branched
 std::uint64_t generated = 0;         ///< children generated
 std::uint64_t infeasible = 0;        ///< nodes found infeasible
 std::uint64_t pruned_generated = 0;  ///< children pruned when generated
 std::uint64_t pruned_extracted = 0;  ///< nodes pruned when extracted
 std::uint64_t pruned_evaluated = 0;  ///< nodes pruned by their bound
 std::uint64_t incumbents = 0;        ///< improvements of the incumbent
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
 * e.g., after a change of the data that leaves the branching valid. */

template< class Node , class V = double >
struct BBSolution {

 BBStatus status = BBStatus::infeasible;   ///< how the search ended
 std::optional< Node > node;               ///< the incumbent node
 V value = BBValueTraits< V >::none();     ///< its value
 V bound = BBValueTraits< V >::none();     ///< the bound
 V gap = BBValueTraits< V >::infinity();   ///< the gap
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
 * when it is branched, and a binary heap of the pairs ( key , index ) gives
 * the next one: the least key first and, among equal keys, the last pushed
 * (the first pushed with Selection::breadth_first). When the key of the
 * selection is not the bound, a second heap of the bounds, whose entries of
 * the removed nodes are deleted lazily, gives the least bound of the open
 * nodes, which the termination test needs. */

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
  : f_track( ( s == Selection::depth_first ) ||
             ( s == Selection::breadth_first ) ||
             ( s == Selection::best_estimate ) ) ,
    f_fifo( s == Selection::breadth_first ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds a node with its bound, its depth and its key

 void push( Node && n , V bound , unsigned depth , V key ) {
  std::size_t i;
  if( v_free.empty() ) {
   i = v_slot.size();
   v_slot.emplace_back();
   }
  else {
   i = v_free.back();
   v_free.pop_back();
   }
  const std::uint64_t seq = ++f_seq;
  v_slot[ i ].item.emplace( Item{ std::move( n ) , bound , depth } );
  v_slot[ i ].seq = seq;
  v_heap.push_back( Key{ key , seq , i } );
  std::push_heap( v_heap.begin() , v_heap.end() , order() );
  if( f_track ) {
   v_bounds.push_back( Key{ bound , seq , i } );
   std::push_heap( v_bounds.begin() , v_bounds.end() , after );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes and returns the next node

 Item pop( void ) {
  std::pop_heap( v_heap.begin() , v_heap.end() , order() );
  const std::size_t i = v_heap.back().idx;
  v_heap.pop_back();
  Item it = std::move( *v_slot[ i ].item );
  v_slot[ i ].item.reset();
  v_slot[ i ].seq = 0;
  v_free.push_back( i );
  return( it );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the least bound of the open nodes (\f$ +\infty \f$ if there is none)

 [[nodiscard]] V best_bound( void ) {
  if( ! f_track )
   return( v_heap.empty() ? BBValueTraits< V >::infinity()
                          : v_slot[ v_heap.front().idx ].item->bound );
  while( ! v_bounds.empty() &&
         ( v_slot[ v_bounds.front().idx ].seq != v_bounds.front().seq ) ) {
   std::pop_heap( v_bounds.begin() , v_bounds.end() , after );
   v_bounds.pop_back();
   }
  return( v_bounds.empty() ? BBValueTraits< V >::infinity()
                           : v_bounds.front().key );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes all the nodes, returning how many there were

 std::size_t clear( void ) {
  const std::size_t n = v_heap.size();
  v_slot.clear();
  v_free.clear();
  v_heap.clear();
  v_bounds.clear();
  return( n );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes all the nodes, giving each to f( item ), and returns how many
 /// there were

 template< class F >
 std::size_t drain( F && f ) {
  for( const Key & k : v_heap )
   f( std::move( *v_slot[ k.idx ].item ) );
  return( clear() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if there is no open node

 [[nodiscard]] bool empty( void ) const noexcept { return( v_heap.empty() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of open nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_heap.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a slot: the node (empty if free) and the number of its push (0 if free)
 struct Slot {
  std::optional< Item > item;
  std::uint64_t seq = 0;
  };

 /// an entry of a heap
 struct Key {
  V key;
  std::uint64_t seq;
  std::size_t idx;
  };

 // the order of the heaps: a comes after b if its key is larger or, for
 // equal keys, if it was pushed earlier
 static bool after( const Key & a , const Key & b ) noexcept {
  return( ( a.key > b.key ) || ( ( a.key == b.key ) && ( a.seq < b.seq ) ) );
  }

 // the same, but with the earlier pushed first for equal keys
 static bool after_fifo( const Key & a , const Key & b ) noexcept {
  return( ( a.key > b.key ) || ( ( a.key == b.key ) && ( a.seq > b.seq ) ) );
  }

 // the order of the heap of the keys
 [[nodiscard]] auto order( void ) const noexcept {
  return( f_fifo ? after_fifo : after );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool f_track;                         ///< true if v_bounds is kept
 bool f_fifo;                          ///< first in, first out on ties
 std::uint64_t f_seq = 0;              ///< the pushes so far
 std::vector< Slot > v_slot;           ///< the slots
 std::vector< std::size_t > v_free;    ///< the free slots
 std::vector< Key > v_heap;            ///< the heap of the keys
 std::vector< Key > v_bounds;          ///< the heap of the bounds

 };  // end( class( BBQueue ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// the evaluator of a problem without make_worker(): the problem itself

template< class P >
struct BBSelf {
 const P & p;
 BasicEvaluation< bb_value_t< P > > evaluate( typename P::node_type & n ) {
  return( p.evaluate( n ) );
  }
 template< class F >
 void branch( const typename P::node_type & n , F && f ) {
  p.branch( n , std::forward< F >( f ) );
  }
 };

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
/// true if a node of the given bound cannot improve the incumbent inc (both
/// in the minimization) by more than the tolerance

template< class V >
bool bb_prunable( const V & bound , const V & inc ,
                  const BasicBBOptions< V > & o )
{
 using Tr = BBValueTraits< V >;
 if( ! Tr::finite( inc ) )
  return( bound == Tr::infinity() );
 return( bound >= inc - Tr::tolerance( inc , o.abs_gap , o.rel_gap ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the leaves of the part of the tree explored by a worker
/** low is the least bound (in the minimization) of the nodes pruned by the
 * tolerance or left open, so that \f$ \min \{ \mathrm{low} , z \} \f$,
 * \f$ z \f$ being the incumbent, is a bound on the optimum; the leaves
 * themselves are kept only if keep. */

template< class Node , class V = double >
struct BBFrontier {

 explicit BBFrontier( bool k = false ) : keep( k ) {}

 /// adds the leaf n of bound b (in the minimization), lowering low unless
 /// the subtree of n holds no better solution than its own
 void add( Node && n , V b , unsigned depth , BBLeafKind kind ,
           bool lower = true ) {
  if( lower )
   low = std::min( low , b );
  if( keep )
   leaves.push_back( BBLeaf< Node , V >{ std::move( n ) , b , depth ,
                                         kind } );
  }

 /// adds all the nodes of q as leaves of the given kind, returning their
 /// number
 std::size_t add( BBQueue< Node , V > & q , BBLeafKind kind ) {
  low = std::min( low , q.best_bound() );
  if( ! keep )
   return( q.clear() );
  return( q.drain( [ & ]( auto && it ) {
   leaves.push_back( BBLeaf< Node , V >{ std::move( it.node ) , it.bound ,
                                     it.depth , kind } );
   } ) );
  }

 /// appends the leaves of another worker
 void merge( BBFrontier && o ) {
  low = std::min( low , o.low );
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
 V low = BBValueTraits< V >::infinity();        ///< the least bound
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
 const V inf = BBValueTraits< V >::infinity();
 if( bb_prunable( bound , inc() , o ) ) {
  ++st.pruned_extracted;
  fr.add( std::move( n ) , bound , depth , BBLeafKind::pruned );
  return;
  }
 const auto cutoff = [ & ]( void ) {
  if constexpr( requires { w.cutoff( inc() ); } )
   w.cutoff( bb_signed< P >( inc() ) );
  };
 cutoff();
 const BasicEvaluation< V > e = w.evaluate( n );
 ++st.evaluated;
 if( e.infeasible ) {
  ++st.infeasible;
  fr.add( std::move( n ) , inf , depth , BBLeafKind::infeasible , false );
  return;
  }
 const V b = std::max( bound , bb_signed< P >( e.bound ) );
 if( e.feasible && improve( bb_signed< P >( *e.feasible ) , n ) )
  ++st.incumbents;
 const BBLeafKind kind = e.feasible ? BBLeafKind::feasible
                                    : BBLeafKind::pruned;
 if( bb_prunable( b , inc() , o ) ) {
  ++st.pruned_evaluated;
  fr.add( std::move( n ) , b , depth , kind );
  return;
  }
 ++st.branched;
 const std::uint64_t g = st.generated;
 cutoff();
 w.branch( std::as_const( n ) , [ & ]( auto && c , V est ) {
  ++st.generated;
  const V cb = std::max( b , bb_signed< P >( est ) );
  if( bb_prunable( cb , inc() , o ) ) {
   ++st.pruned_generated;
   if( fr.keep )
    fr.add( Node( std::forward< decltype( c ) >( c ) ) , cb , depth + 1 ,
            BBLeafKind::pruned );
   else
    fr.low = std::min( fr.low , cb );
   return;
   }
  push( Node( std::forward< decltype( c ) >( c ) ) , cb , depth + 1 );
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
  std::optional< Node > inc_node;
  const auto improve = [ & ]( Value v , const Node & n ) {
   if( v >= inc )
    return( false );
   inc = v;
   inc_node = n;
   if( o.on_incumbent )
    o.on_incumbent( detail::bb_signed< P >( v ) );
   return( true );
   };

  using Item = typename BBQueue< Node , Value >::Item;
  BBQueue< Node , Value > q( sel );
  std::optional< Item > dive;              // the node of the dive, if any
  std::vector< Item > kids;                // the children of a node
  detail::BBFrontier< Node , Value > fr( o.keep_frontier );  // the leaves

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
   if( q.empty() && ! dive )
    break;
   if( detail::bb_prunable( best , inc , o ) ) {   // nothing can improve
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
                            [ & ] { return( inc ); } , improve ,
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
    if( ! Tr::finite( inc ) ||
        Tr::within( kids[ c ].bound , b0 , inc , o.plunge_gap ) )
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
  if( ( sol.status == BBStatus::optimal ) && ! Tr::finite( inc ) )
   sol.status = BBStatus::infeasible;
  sol.node = std::move( inc_node );
  if( Tr::finite( inc ) )
   sol.value = detail::bb_signed< P >( inc );
  const Value lb = std::min( inc , fr.low );
  sol.bound = detail::bb_signed< P >( lb );
  sol.gap = Tr::finite( inc ) && Tr::finite( lb ) ? Value( inc - lb ) : inf;
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
