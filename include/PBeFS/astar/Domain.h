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
 *         Universita' di Pisa \n
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

#include <type_traits>

#include <utility>

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

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

 /// a callable accepting an index, used only to state the concepts
 struct index_sink {
  void operator()( std::size_t ) const noexcept {}
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
/// the operator of the root, which no state was generated by

inline constexpr std::uint8_t no_op = 0xFF;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a Domain that knows the operator by which a state was generated
/** for_each_successor( s , h_s , op , f ) calls f( t , c , h_t , op_t ) for
 * every successor t of s except the one that the inverse of op leads to,
 * op being the operator that generated s (no_op at the root) and op_t <
 * no_op the one that generates t. Since the inverse move costs as much as
 * the move, going back to the parent never gives a cheaper path, and the
 * searches skip it without looking it up [see Burns et al., SoCS 2012]. */

template< class D >
concept OperatorDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ,
           typename D::cost_type h , std::uint8_t op ) {
  d.for_each_successor( s , h , op , detail::op_sink<
                         typename D::state_type , typename D::cost_type >{} );
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

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Domain */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File Domain.h ----------------------------*/
/*--------------------------------------------------------------------------*/
