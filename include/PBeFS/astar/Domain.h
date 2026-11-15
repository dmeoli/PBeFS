/*--------------------------------------------------------------------------*/
/*----------------------------- File Domain.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the concepts a search problem has to satisfy in order to
 * be searched by the algorithms of the library, and for the few type
 * helpers built upon them.
 *
 * A domain describes an implicit state-space graph in the sense of
 * Edelkamp and Schroedl: an initial state, a goal predicate, and an
 * expansion procedure that enumerates the successors of a state together
 * with the cost of the edge and the heuristic value of the successor. All
 * the searches are templates over the domain, so that no virtual call is
 * paid in the inner loop [see Burns et al., SoCS 2012].
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

#ifndef __PBeFS_Domain
 #define __PBeFS_Domain
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <span>

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
/*--------------------------------- EVENTS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// the cost of the edge ( u , v ) changed from from to to
/** An absent edge costs \f$ \infty \f$ [see infinite_cost()], so that an
 * edge that appears is a change from \f$ \infty \f$ and one that
 * disappears a change to \f$ \infty \f$, and a state blocked or freed is
 * the change of each edge into it and out of it [see Ramalingam and Reps,
 * TR 1087, 1992, footnote 1]. The old cost is in the event because the
 * Domain no longer knows it, and because the repair of an incremental
 * search treats a decrease in \f$ O( 1 ) \f$ and an increase only when the
 * edge was the one of the parent [see figure 7 of Koenig, Likhachev and
 * Furcy, AIJ 2004]. */

template< class State , class Cost >
struct EdgeCostChanged {
 State u;                                         ///< the tail
 State v;                                         ///< the head
 Cost from;                                       ///< the cost before
 Cost to;                                         ///< the cost after
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the state t became a goal (added) or stopped being one

template< class State >
struct GoalChanged {
 State t;                                         ///< the state
 bool added;                                      ///< whether now a goal
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the direction of the changes of costs that a Domain makes: an
/// incremental search compiles out the repair of the direction that never
/// happens

enum class CostChanges {
 increase_only ,                                  ///< costs only go up
 decrease_only ,                                  ///< costs only go down
 both                                             ///< either way
 };

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// a callable accepting a successor, used only to state the concepts
 template< class State , class Cost >
 struct successor_sink {
  void operator()( const State & , Cost , Cost ) const noexcept {}
  };

 /// a callable accepting a successor and its operator, used only to state
 /// the concepts
 template< class State , class Cost >
 struct op_sink {
  void operator()( const State & , Cost , Cost , std::uint8_t ) const
   noexcept {}
  };

 /// a callable accepting a predecessor, used only to state the concepts
 template< class State , class Cost >
 struct pred_sink {
  void operator()( const State & , Cost ) const noexcept {}
  };

 /// a callable accepting a successor without its heuristic, used only to
 /// state the concepts
 template< class State , class Cost >
 struct state_sink {
  void operator()( const State & , Cost ) const noexcept {}
  };

 /// a callable accepting an index, used only to state the concepts
 struct index_sink {
  void operator()( std::size_t ) const noexcept {}
  };

 /// a callable accepting the events, used only to state the concepts
 template< class State , class Cost >
 struct event_sink {
  void operator()( const EdgeCostChanged< State , Cost > & ) const noexcept {}
  void operator()( const GoalChanged< State > & ) const noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// a search domain
/** A type D is a Domain if it provides:
 *
 * - the type state_type, a regular value type compared with ==;
 *
 * - the type cost_type, an arithmetic type for g, h and the edge costs;
 *
 * - initial_state(), the root of the search;
 *
 * - is_goal( s ), the goal predicate;
 *
 * - heuristic( s ), the estimate h( s ) of the cost from s to a goal;
 *
 * - hash( s ), a well spread 64-bit hash of s, used by the closed list;
 *
 * - for_each_successor( s , h_s , f ), which calls f( t , c , h_t ) for
 *   every successor t of s, c being the cost of the edge ( s , t ) and h_t
 *   = h( t ); the value h_s = h( s ) is passed in so that h_t can be
 *   computed incrementally. */

template< class D >
concept Domain = requires( const D & d , const typename D::state_type & s ,
                           typename D::cost_type h ) {
 typename D::state_type;
 typename D::cost_type;
 requires std::is_arithmetic_v< typename D::cost_type >;
 requires std::equality_comparable< typename D::state_type >;
 requires std::copyable< typename D::state_type >;
 { d.initial_state() } -> std::convertible_to< typename D::state_type >;
 { d.is_goal( s ) } -> std::convertible_to< bool >;
 { d.heuristic( s ) } -> std::convertible_to< typename D::cost_type >;
 { d.hash( s ) } -> std::convertible_to< std::uint64_t >;
 d.for_each_successor( s , h ,
                       detail::successor_sink< typename D::state_type ,
                                               typename D::cost_type >{} );
 };

/*--------------------------------------------------------------------------*/
/// a Domain whose states are described by a set of features
/** A state is a set of propositions x_0 , ... , x_n taken from a universe
 * of num_features() features, and for_each_feature( s , f ) calls f( i )
 * for every feature i that holds in s. This is what Zobrist hashing needs
 * [see Kishimoto, Fukunaga and Botea, AIJ 2013]. */

template< class D >
concept FeatureDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.num_features() } -> std::convertible_to< std::size_t >;
  d.for_each_feature( s , detail::index_sink{} );
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain that also gives the features changed by a move
/** for_each_feature_change( s , t , f ) calls f( i ) for every feature i
 * that holds in exactly one of s and of t, a successor of s: since a
 * Zobrist hash is the exclusive or of the entries of the features, the
 * hash of t is that of s with the entries of these features added, and it
 * costs as many lookups as the features that the move changes rather than
 * as the features of t [see Zobrist, 1970]. */

template< class D >
concept IncrementalFeatureDomain = FeatureDomain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  d.for_each_feature_change( s , s , detail::index_sink{} );
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain that also provides a feature projection
/** abstract_feature( i ) maps the feature i to an abstract feature in
 * [ 0 , num_abstract_features() ): this is the projection A of Abstract
 * Zobrist hashing [see Jinnai and Fukunaga, JAIR 2017]. */

template< class D >
concept AbstractFeatureDomain = FeatureDomain< D > &&
 requires( const D & d , std::size_t i ) {
  { d.num_abstract_features() } -> std::convertible_to< std::size_t >;
  { d.abstract_feature( i ) } -> std::convertible_to< std::size_t >;
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain whose features are split into groups
/** feature_group( i ) maps the feature i to a group in
 * [ 0 , num_feature_groups() ), the features of a group being mutually
 * exclusive, as the values of a variable of SAS+ (an atom group): a move
 * that deletes a feature of a group adds another one of the same group.
 * The automatic distributions [see AutoDistribution.h] partition the
 * transition graph of each group [see Jinnai and Fukunaga, JAIR 2017,
 * section 6]; without groups, they infer them from the moves. */

template< class D >
concept FeatureGroupDomain = FeatureDomain< D > &&
 requires( const D & d , std::size_t i ) {
  { d.num_feature_groups() } -> std::convertible_to< std::size_t >;
  { d.feature_group( i ) } -> std::convertible_to< std::size_t >;
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain that lists its ground transitions
/** for_each_feature_transition( f ) calls f( i , j ) for every ground
 * operator and every feature i that it deletes and j that it adds in the
 * same group: the arcs of the domain transition graphs of SAS+, each
 * counted as many times as the operators that realize it, which are the
 * weights of the cut of GRAZHDA* [see Jinnai and Fukunaga, JAIR 2017,
 * section 6.2.1]. Without it, a FeatureModel counts the distinct operators
 * seen in its sample [see AutoDistribution.h]. */

template< class D >
concept FeatureTransitionDomain = FeatureDomain< D > &&
 requires( const D & d ) {
  d.for_each_feature_transition( []( std::size_t , std::size_t ) {} );
  };

/*--------------------------------------------------------------------------*/
/// a Domain with a state abstraction
/** abstract_state( s ) returns the image of s in an abstract space, coded
 * as an integer: it is the projection by which AHDA* keeps the successors
 * of a state on the same processor [see Burns et al., JAIR 2010]. */

template< class D >
concept AbstractStateDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.abstract_state( s ) } -> std::convertible_to< std::uint64_t >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain whose state space is partitioned into nblocks
/** Besides an abstraction with num_nblocks() abstract states, given by
 * nblock( s ), the Domain enumerates with for_each_nblock_successor( b , f )
 * the abstract successors of the nblock b, which define the duplicate
 * detection scope of Parallel Best-NBlock-First search [see Burns et al.,
 * JAIR 2010]. */

template< class D >
concept NBlockDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s , std::size_t b ) {
  { d.num_nblocks() } -> std::convertible_to< std::size_t >;
  { d.nblock( s ) } -> std::convertible_to< std::size_t >;
  d.for_each_nblock_successor( b , detail::index_sink{} );
  };

/*--------------------------------------------------------------------------*/
/// a Domain with a heuristic between any two states
/** pairwise_heuristic( s , t ) estimates the cost of a path from s to t; it
 * has to be forward-backward consistent, as required by PA*SE [see
 * Phillips, Likhachev and Koenig, ICAPS 2014]. */

template< class D >
concept PairwiseDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.pairwise_heuristic( s , s ) }
   -> std::convertible_to< typename D::cost_type >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain that estimates the number of edges from a state to a goal
/** distance( s ) estimates the number of edges of a cheapest path from s
 * to a goal, the \f$ d \f$ of the searches that prefer the nodes closer to
 * a goal in edges rather than in cost [see Thayer and Ruml, IJCAI 2011];
 * it need not be admissible, and it is \f$ h \f$ when every edge costs 1. */

template< class D >
concept DistanceDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.distance( s ) } -> std::convertible_to< typename D::cost_type >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain whose heuristic is computed for many states at once
/** heuristic_batch( s , h ) writes in h[ i ] the heuristic of the state
 * s[ i ], for the spans s and h of equal size: a learned heuristic on an
 * accelerator, or any heuristic whose cost per state falls with the
 * number of states, pays once per batch [see Agostinelli et al., Nature
 * Machine Intelligence 2019]. The searches that use it ignore the h that
 * for_each_successor passes with each successor, which may then be 0. */

template< class D >
concept BatchHeuristicDomain = Domain< D > &&
 requires( const D & d , std::span< const typename D::state_type > s ,
           std::span< typename D::cost_type > h ) {
  d.heuristic_batch( s , h );
  };

/*--------------------------------------------------------------------------*/
/// the abstract feature of a feature that a projection ignores, i.e., whose
/// entry in the table of Zobrist hashing is 0

inline constexpr std::uint32_t no_abstract_feature = 0xFFFFFFFF;

/*--------------------------------------------------------------------------*/
/// the code of no operator: the op of the root, and of a state generated
/// by an operator that has no code; given as op, it skips no successor

inline constexpr std::uint8_t no_op = 0xFF;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a Domain that knows the operator by which a state was generated
/** for_each_successor( s , h_s , op , f ) calls f( t , c , h_t , op_t ) for
 * every successor t of s except the one that the inverses of op lead to,
 * i.e., the parent, op being the code of the operator that generated s and
 * op_t that of the operator that generates t. With op = no_op (at the
 * root, or when the operator has no code) no successor is skipped, and a
 * Domain may give op_t = no_op to an operator it gives no code, e.g., one
 * without an inverse or beyond the codes below no_op [see Planning.h].
 * Since the costs are not negative, going back to the parent never gives
 * it a cheaper path, and the searches skip it without looking it up [see
 * Burns et al., SoCS 2012]. */

template< class D >
concept OperatorDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::uint8_t op ) {
  d.for_each_successor( s , h , op , detail::op_sink<
                         typename D::state_type , typename D::cost_type >{} );
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an OperatorDomain that prunes the moves a tree search need not try
/** for_each_tree_successor( s , h , op , f ) calls f( t , c , h_t , op_t )
 * as for_each_successor( s , h , op , f ), but it may skip more of the
 * successors: those of every sequence of two moves, op and op_t, that the
 * Domain declares redundant because another sequence, which it keeps,
 * reaches the same state at a cost not larger, such as two turns of the
 * same face of Rubik's cube or two moves of the same disk of the Towers of
 * Hanoi [see Korf, AAAI 1997]. Every path is then matched by a path of the
 * same end and of a cost not larger whose moves are all kept, so that a
 * search without duplicate detection (IDA*) stays complete and optimal;
 * the searches with a closed list do not use it, since the path kept may
 * go through a state closed with another operator. */

template< class D >
concept TreeOperatorDomain = OperatorDomain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::uint8_t op ) {
  d.for_each_tree_successor( s , h , op , detail::op_sink<
                              typename D::state_type ,
                              typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// a Domain whose moves all raise a level
/** level( s ) is a nonnegative integer that every move raises, and
 * max_level() the largest one, that of the goals: the state space is a
 * levelled graph, as the lattice of the alignment of several sequences,
 * whose level is the number of residues aligned. This is what the
 * hyperplane work distribution of Kobayashi, Kishimoto and Watanabe
 * (IJCAI 2011) needs [see HyperplaneDistribution]. */

template< class D >
concept LevelDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.level( s ) } -> std::convertible_to< std::uint64_t >;
  { d.max_level() } -> std::convertible_to< std::uint64_t >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain with a perfect hash of its states
/** perfect_hash( s ) is an injective map of the states to integers, used
 * by the HDA*[P] work distribution of Burns et al. (JAIR 2010). */

template< class D >
concept PerfectHashDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.perfect_hash( s ) } -> std::convertible_to< std::uint64_t >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain that enumerates the predecessors of a state
/** for_each_predecessor( s , f ) calls f( t , c ) for every state t with
 * an edge ( t , s ), c being its cost; the edges are those that
 * for_each_successor() gives at the same time. */

template< class D >
concept PredecessorDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  d.for_each_predecessor( s , detail::pred_sink< typename D::state_type ,
                                                 typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// a PredecessorDomain that tells the changes of its edges and goals
/** drain_events( f ) calls f( e ) for every event e since the last call, an
 * EdgeCostChanged or a GoalChanged, and forgets them; the events of one
 * edge are coalesced into one, from the cost before the first change to
 * the cost after the last, and those that change nothing are dropped [see
 * Ramalingam and Reps, TR 1087, 1992, p. 2]. A Domain may declare the
 * static member cost_changes, of type CostChanges [see cost_changes_of()].
 * The heuristic must stay consistent for every cost the edges take. */

template< class D >
concept EventDomain = PredecessorDomain< D > &&
 requires( D & d ) {
  d.drain_events( detail::event_sink< typename D::state_type ,
                                      typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// the direction of the changes of costs of D: its member cost_changes if
/// it declares one, CostChanges::both otherwise

template< class D >
[[nodiscard]] constexpr CostChanges cost_changes_of( void ) noexcept
{
 if constexpr( requires { { D::cost_changes } ->
                           std::convertible_to< CostChanges >; } )
  return( D::cost_changes );
 else
  return( CostChanges::both );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a Domain that also enumerates the predecessors of a state, with the
/// heuristic of the backward search
/** Besides the Domain, D provides:
 *
 * - goal_state(), the one state for which is_goal() holds;
 *
 * - backward_heuristic( s ), an estimate of the cost of a path from the
 *   initial state to s;
 *
 * - for_each_predecessor( s , hb_s , f ), which calls f( t , c , hb_t ) for
 *   every state t with an edge ( t , s ) of cost c, hb_t being
 *   backward_heuristic( t ), computed incrementally from hb_s.
 *
 * This is what a bidirectional search needs: the backward search from the
 * goal follows the edges against their direction, guided by the backward
 * heuristic [see Chen et al., IJCAI 2017, section 3]. */

template< class D >
concept BidirectionalDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h ) {
  { d.goal_state() } -> std::convertible_to< typename D::state_type >;
  { d.backward_heuristic( s ) }
   -> std::convertible_to< typename D::cost_type >;
  d.for_each_predecessor( s , h , detail::successor_sink<
                           typename D::state_type , typename D::cost_type >{}
                          );
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a BidirectionalDomain that knows the operator by which a predecessor
/// was generated
/** for_each_predecessor( s , hb , op , f ) calls f( t , c , hb_t , op_t )
 * for every predecessor t of s but the one that the inverse of op leads
 * to, op being the operator by which the backward search generated s
 * (no_op at the goal), as OperatorDomain does for the successors. */

template< class D >
concept BackwardOperatorDomain = BidirectionalDomain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::uint8_t op ) {
  d.for_each_predecessor( s , h , op , detail::op_sink<
                           typename D::state_type , typename D::cost_type >{}
                          );
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a Domain with a positive lower bound on the costs of its edges
/** min_edge_cost() is the \f$ \epsilon > 0 \f$ of the bound
 * \f$ g_F( u ) + g_B( v ) + \epsilon \f$ of the bidirectional searches
 * [see Holte et al., AIJ 2017], 0 being allowed when no such bound is
 * known. */

template< class D >
concept MinCostDomain = Domain< D > && requires( const D & d ) {
 { d.min_edge_cost() } -> std::convertible_to< typename D::cost_type >;
 };

/*--------------------------------------------------------------------------*/
/// a Domain with more than one heuristic
/** num_heuristics() is the number \f$ m \geq 1 \f$ of heuristics, and
 * heuristic( s , i ) the value \f$ h_i( s ) \f$ of the i-th one, for
 * \f$ 0 \leq i < m \f$, \f$ h_0 \f$ being heuristic( s ). The greedy
 * searches alternate their queues over them, break the ties of one by
 * another, or keep the Pareto front of their vectors [see Röger and
 * Helmert, ICAPS 2010]. */

template< class D >
concept MultiHeuristicDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s , std::size_t i ) {
  { d.num_heuristics() } -> std::convertible_to< std::size_t >;
  { d.heuristic( s , i ) } -> std::convertible_to< typename D::cost_type >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain that prefers some of the successors of a state
/** preferences( s ) returns a predicate p, computed once per expansion of
 * s, such that p( t ) is true if the successor t of s is preferred: the
 * helpful actions of the relaxed plan of FF [see Hoffmann and Nebel, JAIR
 * 2001], or the preferred operators of a heuristic [see Richter and
 * Helmert, ICAPS 2009], which the greedy searches put into queues of
 * their own and enforced hill-climbing may expand alone. */

template< class D >
concept PreferredDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  { d.preferences( s )( s ) } -> std::convertible_to< bool >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain that generates the successors of a state without evaluating
/// them
/** for_each_successor_state( s , f ) calls f( t , c ) for every successor
 * t of s, c being the cost of the edge, without computing h( t ): the
 * search with deferred evaluation then computes h( t ) only for the
 * states that it selects [see Richter and Helmert, ICAPS 2009]. */

template< class D >
concept DeferredDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  d.for_each_successor_state( s , detail::state_sink<
                               typename D::state_type ,
                               typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain whose goal is a conjunction of features
/** for_each_goal_feature( f ) calls f( i ) for every feature i of the goal,
 * a state being a goal if and only if it has all of them; the serialized
 * width search achieves them one at a time, and the width searches count
 * those still to achieve, \f$ \#g \f$ [see Lipovetzky and Geffner, ECAI
 * 2012 and AAAI 2017]. */

template< class D >
concept GoalFeatureDomain = FeatureDomain< D > &&
 requires( const D & d ) {
  d.for_each_goal_feature( detail::index_sink{} );
  };

/*--------------------------------------------------------------------------*/
/// a FeatureDomain that gives the relevant features of a state
/** for_each_relevant_feature( s , f ) calls f( i ) for every feature i of a
 * relaxed plan from s, the set R whose features reached on the path to a
 * state are counted by \f$ \#r \f$ in best-first width search [see
 * Lipovetzky and Geffner, AAAI 2017]. */

template< class D >
concept RelevantFeatureDomain = FeatureDomain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  d.for_each_relevant_feature( s , detail::index_sink{} );
  };

/*--------------------------------------------------------------------------*/
/// an OperatorDomain that generates the successors of a given
/// \f$ \Delta f \f$
/** for_each_successor_delta( s , h_s , op , delta , f ) calls f( t , c ,
 * h_t , op_t ) for those successors t of s, but the one that the inverse of
 * op leads to, whose \f$ \Delta f = c + h_t - h_s \f$ is delta (with real
 * costs, up to a slack such as that of cost_slack()), and returns the least
 * \f$ \Delta f \f$ above it of the others, or infinite_cost() if there is
 * none: the operator selection function of Enhanced Partial
 * Expansion A*, which never generates a successor that it does not put in
 * OPEN [see Goldenberg et al., JAIR 2014]. */

template< class D >
concept PartialExpansionDomain = OperatorDomain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::uint8_t op ) {
  { d.for_each_successor_delta( s , h , op , h ,
                                detail::op_sink< typename D::state_type ,
                                               typename D::cost_type >{} ) }
   -> std::convertible_to< typename D::cost_type >;
  };

/*--------------------------------------------------------------------------*/
/// a Domain whose edges are evaluated one at a time
/** for_each_action( s , f ) calls f( a ) for every action a of s, an index
 * smaller than \f$ 2^{32} - 1 \f$, without evaluating the edge
 * \f$ ( s , a ) \f$, which costs little; evaluate_edge( s , h_s , a , f )
 * evaluates it, which may cost much (a collision check, a simulation), and
 * calls f( t , c , h_t ) if the edge exists, t being the state that a leads
 * to and c the cost of the edge, nothing otherwise. The successors given by
 * for_each_successor( s , h_s , f ) are those given by evaluate_edge on the
 * actions of s, so that a search by edges and one by states search the
 * same graph [see Mukherjee, Aine and Likhachev, SoCS 2022]. */

template< class D >
concept EdgeDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::size_t a ) {
  d.for_each_action( s , detail::index_sink{} );
  d.evaluate_edge( s , h , a , detail::successor_sink<
                    typename D::state_type , typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// an EdgeDomain whose edges are cheap or expensive to evaluate
/** expensive_edge( s , a ) tells, without evaluating it, whether the edge
 * \f$ ( s , a ) \f$ is expensive to evaluate: GePA*SE evaluates the cheap
 * edges of a state together, and each expensive one apart [see Mukherjee
 * and Likhachev, SoCS 2023]. */

template< class D >
concept HeterogeneousEdgeDomain = EdgeDomain< D > &&
 requires( const D & d , const typename D::state_type & s , std::size_t a ) {
  { d.expensive_edge( s , a ) } -> std::convertible_to< bool >;
  };

/*--------------------------------------------------------------------------*/
/// an EdgeDomain whose successors are known before the edges are evaluated
/** lazy_edge( s , h_s , a , f ) calls f( t , c , h_t ), at little cost, if
 * the edge \f$ ( s , a ) \f$ may exist, t being the state that a leads to
 * and c a lower bound on the cost of the edge, and nothing if the edge
 * surely does not exist; evaluate_edge then gives the same t and a cost not
 * smaller than c, or no edge. The heuristic has to be consistent also with
 * the costs c, as the lazy searches use them [see Mukherjee, Aine and
 * Likhachev, RA-L 2022]. */

template< class D >
concept LazyEdgeDomain = EdgeDomain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::size_t a ) {
  d.lazy_edge( s , h , a , detail::successor_sink<
                typename D::state_type , typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the state type of a Domain
template< Domain D >
using state_t = typename D::state_type;

/// the cost type of a Domain
template< Domain D >
using cost_t = typename D::cost_type;

/// the operator that generated a state, kept in a message [see
/// OperatorDomain]; an empty type, with op always no_op, when D is not an
/// OperatorDomain, so that a message with [[no_unique_address]] on it does
/// not grow

template< class D >
struct OperatorOf {
 std::uint8_t op = no_op;                         ///< the operator
 constexpr OperatorOf( void ) = default;
 constexpr OperatorOf( std::uint8_t o ) : op( o ) {}
 };

template< class D > requires( ! OperatorDomain< D > )
struct OperatorOf< D > {
 static constexpr std::uint8_t op = no_op;        ///< no operator
 constexpr OperatorOf( void ) = default;
 constexpr OperatorOf( std::uint8_t ) {}
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the successors of s but its parent, when D knows the operators
/** Calls f( t , c , h_t , op_t ) for every successor t of s, s having been
 * generated by the operator op: without the one that the inverse of op
 * leads to and with op_t the operator of t if D is an OperatorDomain, with
 * all of them and op_t = no_op otherwise. */

template< Domain D , class F >
void for_each_child( const D & d , const state_t< D > & s , cost_t< D > h ,
                     std::uint8_t op , F && f ) {
 if constexpr( OperatorDomain< D > )
  d.for_each_successor( s , h , op , std::forward< F >( f ) );
 else
  d.for_each_successor( s , h , [ & ]( const state_t< D > & t , cost_t< D > c ,
                                       cost_t< D > ht ) {
   f( t , c , ht , no_op );
   } );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the successors of s that a tree search tries
/** As for_each_child(), with the moves pruned by for_each_tree_successor()
 * if D is a TreeOperatorDomain [see TreeOperatorDomain]. */

template< Domain D , class F >
void for_each_tree_child( const D & d , const state_t< D > & s ,
                          cost_t< D > h , std::uint8_t op , F && f ) {
 if constexpr( TreeOperatorDomain< D > )
  d.for_each_tree_successor( s , h , op , std::forward< F >( f ) );
 else
  for_each_child( d , s , h , op , std::forward< F >( f ) );
 }

/*--------------------------------------------------------------------------*/
/// the cost of the path from path.front() to path.back(), each step by the
/// cheapest edge between its two states (0 if the path is empty)
/** A search returns the cost of its path with it, and the g of a goal may
 * exceed the cost of the path that its parents give, when a node of the
 * path has been reached through a cheaper path after the goal was
 * generated: this is the cost of the path itself. */

template< Domain D >
[[nodiscard]] cost_t< D > path_cost( const D & d ,
                                     const std::vector< state_t< D > > & path )
{
 cost_t< D > sum( 0 );
 for( std::size_t i = 1 ; i < path.size() ; ++i ) {
  bool found = false;
  cost_t< D > best( 0 );
  d.for_each_successor( path[ i - 1 ] , cost_t< D >( d.heuristic(
                         path[ i - 1 ] ) ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > c ,
                               cost_t< D > ) {
   if( ( t == path[ i ] ) && ( ! found || ( c < best ) ) ) {
    found = true;
    best = c;
    }
   } );
  sum += best;
  }
 return( sum );
 }

/*--------------------------------------------------------------------------*/
/// the predecessors of s but its child, when D knows the operators
/** Calls f( t , c , hb_t , op_t ) for every predecessor t of s, s having
 * been generated by the backward search by the operator op, as
 * for_each_child() does for the successors. */

template< BidirectionalDomain D , class F >
void for_each_parent( const D & d , const state_t< D > & s , cost_t< D > h ,
                      std::uint8_t op , F && f ) {
 if constexpr( BackwardOperatorDomain< D > )
  d.for_each_predecessor( s , h , op , std::forward< F >( f ) );
 else
  d.for_each_predecessor( s , h , [ & ]( const state_t< D > & t ,
                                         cost_t< D > c , cost_t< D > ht ) {
   f( t , c , ht , no_op );
   } );
 }

/*--------------------------------------------------------------------------*/
/// the number of heuristics of a Domain: num_heuristics() if it is a
/// MultiHeuristicDomain, 1 otherwise

template< Domain D >
[[nodiscard]] std::size_t heuristics_of( const D & d )
{
 if constexpr( MultiHeuristicDomain< D > )
  return( std::size_t( d.num_heuristics() ) );
 else {
  static_cast< void >( d );
  return( 1 );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the i-th heuristic of a Domain at s [see MultiHeuristicDomain], which
/// is heuristic( s ) for i = 0

template< Domain D >
[[nodiscard]] cost_t< D > heuristic_of( const D & d , const state_t< D > & s ,
                                        std::size_t i )
{
 if constexpr( MultiHeuristicDomain< D > )
  if( i != 0 )
   return( cost_t< D >( d.heuristic( s , i ) ) );
 static_cast< void >( i );
 return( cost_t< D >( d.heuristic( s ) ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Domain */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File Domain.h ----------------------------*/
/*--------------------------------------------------------------------------*/
