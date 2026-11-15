/*--------------------------------------------------------------------------*/
/*-------------------------- File SymbolicDomain.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the concepts a search problem has to satisfy in order to
 * be searched symbolically, i.e., on sets of states represented by binary
 * decision diagrams [see SymbolicSearch.h], and for the descriptions of the
 * operators and of the heuristics that such a problem gives.
 *
 * A SymbolicDomain describes its states as the values of finitely many
 * variables, each with a finite range, and its transitions as operators
 * with a cost, a precondition (values that some variables must have) and
 * effects, each of which sets a variable to a constant, as in the SAS+
 * formalism of Bäckström and Nebel, or to the value that another variable
 * had before the move, which describes in a few effects the moves of the
 * permutation puzzles; the variables that no effect writes keep their
 * value. The domain knows nothing of the diagrams: the search compiles the
 * description into them [see SymbolicModel.h].
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

#ifndef __PBeFS_SymbolicDomain
 #define __PBeFS_SymbolicDomain
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <limits>

#include <type_traits>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the value of a variable, a fact \f$ v = d \f$

struct SymbolicFact {
 std::uint32_t var;                     ///< the variable
 std::uint32_t value;                   ///< its value
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an effect: var takes value, or the value that the variable source had

struct SymbolicEffect {
 std::uint32_t var;                     ///< the variable written
 std::uint32_t value = 0;               ///< the constant, if not a copy
 bool copy = false;                     ///< true if var takes source
 std::uint32_t source = 0;              ///< the variable copied

 /// var takes the constant value
 [[nodiscard]] static SymbolicEffect assign( std::uint32_t var ,
                                             std::uint32_t value ) {
  return( SymbolicEffect{ var , value , false , 0 } );
  }

 /// var takes the value that source had
 [[nodiscard]] static SymbolicEffect copy_of( std::uint32_t var ,
                                              std::uint32_t source ) {
  return( SymbolicEffect{ var , 0 , true , source } );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an operator: applicable where every fact of pre holds, it applies all
/// its effects at once (each reading the values before the move) and
/// costs cost; a variable is written by one effect at most

template< class Cost >
struct SymbolicOperator {
 Cost cost = 1;                          ///< the cost, \f$ \geq 0 \f$
 std::vector< SymbolicFact > pre;        ///< the precondition
 std::vector< SymbolicEffect > effects;  ///< the effects
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a term of an additive heuristic
/** The term on the variables vars (at most a few) has the value
 * table[ i ] on the values \f$ d_0 , \ldots , d_{k-1} \f$ of vars[ 0 ] ,
 * ... , vars[ k - 1 ], i being \f$ \sum_j d_j \prod_{l < j} |D_l| \f$
 * (the first variable the fastest), and the heuristic of a state is the sum
 * of all its terms; a value equal to infinite_term< Cost >() marks the
 * states from which no goal is reachable. */

template< class Cost >
struct SymbolicTerm {
 std::vector< std::uint32_t > vars;     ///< the variables of the term
 std::vector< Cost > table;             ///< its values
 };

/// the value of a term that marks a dead end

template< class Cost >
[[nodiscard]] constexpr Cost infinite_term( void ) noexcept {
 return( std::numeric_limits< Cost >::max() );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// a callable accepting an operator, used only to state the concepts
 template< class Cost >
 struct operator_sink {
  void operator()( const SymbolicOperator< Cost > & ) const noexcept {}
  };

 /// a callable accepting a partial assignment, used only to state the
 /// concepts
 struct facts_sink {
  void operator()( const std::vector< SymbolicFact > & ) const noexcept {}
  };

 /// a callable accepting a term, used only to state the concepts
 template< class Cost >
 struct term_sink {
  void operator()( const SymbolicTerm< Cost > & ) const noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// a domain searched on sets of states
/** A type D is a SymbolicDomain if it provides:
 *
 * - the type cost_type, an integral type for the costs, which are \f$
 *   \geq 0 \f$;
 *
 * - num_variables(), the number n of the variables, and range( v ), the
 *   number \f$ | D_v | \geq 1 \f$ of the values \f$ 0 , \ldots , | D_v | - 1
 *   \f$ of the variable v;
 *
 * - initial_values(), the value of each variable in the initial state;
 *
 * - for_each_goal( f ), which calls f( p ) with partial assignments p, a
 *   state being a goal if it satisfies at least one of them;
 *
 * - for_each_operator( f ), which calls f( o ) with every operator. */

template< class D >
concept SymbolicDomain = requires( const D & d , std::uint32_t v ) {
 typename D::cost_type;
 requires std::is_integral_v< typename D::cost_type >;
 { d.num_variables() } -> std::convertible_to< std::size_t >;
 { d.range( v ) } -> std::convertible_to< std::size_t >;
 { d.initial_values() } -> std::convertible_to< std::vector< std::uint32_t > >;
 d.for_each_goal( detail::facts_sink{} );
 d.for_each_operator( detail::operator_sink< typename D::cost_type >{} );
 };

/*--------------------------------------------------------------------------*/
/// a SymbolicDomain with an additive heuristic towards the goal
/** for_each_heuristic_term( f ) calls f( t ) with each term of a
 * heuristic h, the sum of the terms, which must be admissible, and
 * consistent for BDDA* to be optimal [see SymbolicSearch.h]. */

template< class D >
concept SymbolicHeuristicDomain = SymbolicDomain< D > &&
 requires( const D & d ) {
  d.for_each_heuristic_term(
   detail::term_sink< typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// a SymbolicDomain with an additive heuristic towards the initial state
/** for_each_backward_heuristic_term( f ) gives the terms of an estimate of
 * the cost from the initial state to a state, which the backward searches
 * use as SymbolicHeuristicDomain the forward ones. */

template< class D >
concept SymbolicBackwardHeuristicDomain = SymbolicDomain< D > &&
 requires( const D & d ) {
  d.for_each_backward_heuristic_term(
   detail::term_sink< typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/// a SymbolicDomain whose states are also those of a Domain
/** decode( x ) gives the state of the Domain whose variables have the
 * values x, so that a path found symbolically can be checked, or used, on
 * the explicit Domain [see Domain.h]. */

template< class D >
concept DecodableSymbolicDomain = SymbolicDomain< D > &&
 requires( const D & d , const std::vector< std::uint32_t > & x ) {
  typename D::state_type;
  { d.decode( x ) } -> std::convertible_to< typename D::state_type >;
  };

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SymbolicDomain */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SymbolicDomain.h ------------------------*/
/*--------------------------------------------------------------------------*/
