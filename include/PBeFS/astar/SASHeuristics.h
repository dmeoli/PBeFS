/*--------------------------------------------------------------------------*/
/*-------------------------- File SASHeuristics.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the heuristics of classical planning on a SASTask [see
 * SASTask.h]: the class RelaxedHeuristic, which gives \f$ h^{\max} \f$
 * (admissible), \f$ h^{\rm add} \f$ and \f$ h^{\rm FF} \f$ (not
 * admissible) from the relaxation that ignores the delete lists, the class
 * LMCut, the landmark-cut heuristic (admissible), the class
 * SASPatternDatabase, the pattern database of a projection of the task
 * onto a subset of its variables, and the class SASMergeAndShrink, the
 * merge-and-shrink abstraction (both admissible and consistent).
 *
 * Each heuristic is evaluated on the values of the variables in a state,
 * by const methods that may be called by several threads at once; the
 * value dead_end marks a state from which no goal is reachable.
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

#ifndef __PBeFS_SASHeuristics
 #define __PBeFS_SASHeuristics
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <limits>

#include <vector>

#include "PBeFS/astar/SASTask.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the value of a heuristic on a dead end

inline constexpr int dead_end = std::numeric_limits< int >::max();

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS RelaxedHeuristic -------------------------*/
/*--------------------------------------------------------------------------*/
/// the heuristics of the delete relaxation
/** With the unary operators \f$ pre( o ) \to q \f$, one for each operator
 * \f$ o \f$ and each \f$ q \f$ of its add list, the cost \f$ h( p ) \f$
 * of a fact is 0 if it holds in the state \f$ s \f$ and otherwise
 * \f[
 *  h( q ) = \min_{ o : q \in add( o ) }
 *           \{ c( o ) + \textstyle{\bigoplus}_{ p \in pre( o ) } h( p ) \} ,
 * \f]
 * computed by the generalized Dijkstra algorithm, \f$ \bigoplus \f$ being
 * the maximum for \f$ h^{\max} \f$ and the sum for \f$ h^{\rm add} \f$,
 * which give the same aggregation over the goal [see Bonet and Geffner,
 * AIJ 2001]. \f$ h^{\rm FF} \f$ is the cost of the relaxed plan traced
 * back from the goal through the operators that attain the minima of
 * \f$ h^{\rm add} \f$, each counted once [see Hoffmann and Nebel, JAIR
 * 2001, and Keyder and Geffner, ECAI 2008]. Only \f$ h^{\max} \f$ is
 * admissible. */

class RelaxedHeuristic {

 public:

 //! [planning-sig-relaxed]
 /// the heuristic of the relaxation
 enum class Kind {
  hmax ,     ///< the maximum over the preconditions
  hadd ,     ///< the sum over the preconditions
  hff        ///< the cost of the relaxed plan of h^add
  };

 /// constructor, on the STRIPS view of a task, which must outlive it

 RelaxedHeuristic( const StripsTask & task , Kind kind );

 /// the value on the state of the values s, one per variable

 [[nodiscard]] int operator()( const std::vector< int > & s ) const;
 //! [planning-sig-relaxed]

 private:

 struct Unary {
  int eff;                                   ///< the fact added
  int op;                                    ///< the operator
  int cost;                                  ///< its cost
  int npre;                                  ///< its preconditions
  };

 const StripsTask & f_task;                  ///< the task
 Kind f_kind;                                ///< the heuristic
 std::vector< Unary > v_unary;               ///< the unary operators
 std::vector< int > v_first;                 ///< CSR of the unary ops
 std::vector< int > v_of;                    ///< by precondition
 std::vector< int > v_free;                  ///< with no preconditions

 };  // end( class( RelaxedHeuristic ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LMCut -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the landmark-cut heuristic of Helmert and Domshlak (ICAPS 2009)
/** Each round computes \f$ h^{\max} \f$ with the current costs, chooses
 * for each operator a precondition of largest \f$ h^{\max} \f$, and finds
 * in the justification graph so obtained the goal zone, i.e., the facts
 * that reach the goal through operators of cost zero, and the cut, i.e.,
 * the operators that lead into the goal zone from the facts reachable from
 * the state without entering it. The cut is a disjunctive action landmark:
 * the least cost \f$ m \f$ of its operators is added to the heuristic and
 * subtracted from the cost of each of them, until \f$ h^{\max} \f$ of the
 * goal is 0. The sum of the \f$ m \f$ is admissible, and it dominates
 * \f$ h^{\max} \f$. */

class LMCut {

 public:

 //! [planning-sig-lmcut]
 /// constructor, on the STRIPS view of a task, which must outlive it

 explicit LMCut( const StripsTask & task );

 /// the value on the state of the values s, one per variable

 [[nodiscard]] int operator()( const std::vector< int > & s ) const;
 //! [planning-sig-lmcut]

 private:

 struct Op {
  std::vector< int > pre;                    ///< the preconditions
  std::vector< int > add;                    ///< the add list
  int cost;                                  ///< the cost
  };

 const StripsTask & f_task;                  ///< the task
 int f_true;                                 ///< the fact that always holds
 int f_goal;                                 ///< the fact of the goal
 std::vector< Op > v_op;                     ///< with the goal operator
 std::vector< std::vector< int > > v_pre_of; ///< the ops of a precondition
 std::vector< std::vector< int > > v_add_of; ///< the ops adding a fact

 };  // end( class( LMCut ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SASPatternDatabase ------------------------*/
/*--------------------------------------------------------------------------*/
/// the pattern database of a projection of a SAS+ task
/** The projection onto the pattern \f$ P \f$, a set of variables, keeps
 * the facts of the variables of \f$ P \f$ in the states, in the goal and
 * in the operators; the distances of all its abstract states to the
 * abstract goal, computed once by a Dijkstra search backward from the goal
 * states (each operator regressed, the variables of its effects with no
 * value before taking any value), are admissible and consistent estimates
 * of the distances of the states [see Edelkamp, ECP 2001]. The abstract
 * state of the values \f$ s \f$ is the index \f$ \sum_{ v \in P } s_v
 * m_v \f$, \f$ m_v \f$ being the product of the ranges of the variables
 * of \f$ P \f$ before \f$ v \f$. */

class SASPatternDatabase {

 public:

 //! [planning-sig-pdb]
 /// constructor, building the distances of the projection of task onto
 /// the given variables

 SASPatternDatabase( const SASTask & task , std::vector< int > pattern );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a pattern of at most max_size abstract states
 /** The goal variables are added while the size allows, then, in the
  * order of their indices, the variables of the preconditions of the
  * operators that change a variable of the pattern (its parents in the
  * causal graph), as long as one of them fits. */

 [[nodiscard]] static std::vector< int > greedy_pattern(
                               const SASTask & task , std::size_t max_size );
 //! [planning-sig-pdb]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the distance of the abstract state of the values s, one per variable

 [[nodiscard]] int operator()( const std::vector< int > & s ) const {
  std::size_t i = 0;
  for( std::size_t k = 0 ; k < v_pattern.size() ; ++k )
   i += std::size_t( s[ std::size_t( v_pattern[ k ] ) ] ) * v_mult[ k ];
  return( v_dist[ i ] );
  }

 /// the pattern

 [[nodiscard]] const std::vector< int > & pattern( void ) const noexcept {
  return( v_pattern );
  }

 /// the number of abstract states

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_dist.size() );
  }

 private:

 std::vector< int > v_pattern;              ///< the variables
 std::vector< std::size_t > v_mult;         ///< their multipliers
 std::vector< int > v_dist;                 ///< the distances

 };  // end( class( SASPatternDatabase ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SASMergeAndShrink -------------------------*/
/*--------------------------------------------------------------------------*/
/// the merge-and-shrink heuristic of a SAS+ task
/** An abstraction is a transition system whose states stand for sets of
 * states of the task, with a transition labelled by an operator \f$ o \f$
 * from \f$ a \f$ to \f$ b \f$ if \f$ o \f$ leads from a state of
 * \f$ a \f$ to one of \f$ b \f$; the distance of the abstract state of a
 * state to an abstract goal is an admissible and consistent estimate of its
 * distance to a goal [see Helmert, Haslum, Hoffmann and Nissim, JACM 2014].
 * The atomic abstraction of a variable has its values as states. The
 * abstraction is built by a linear merge: the variables are taken in the
 * order of the causal graph and of the goal (the first goal variable, then
 * repeatedly the variable of least index that is a parent, in the causal
 * graph, of one already taken, or else the next goal variable), and each
 * one is merged into the current abstraction by the synchronized product
 * with its atomic abstraction. Before each product the current abstraction
 * is shrunk to at most \f$ \lfloor N / | D_v | \rfloor \f$ states, so that
 * the product has at most \f$ N \f$, by the greedy bisimulation of
 * Nissim, Hoffmann and Helmert (IJCAI 2011): starting from the partition
 * of the states by their distance \f$ h \f$ to the goal, a block is split
 * by the labels and the blocks of the transitions \f$ ( a , o , b ) \f$
 * with \f$ h( a ) = c( o ) + h( b ) \f$, as long as the number of blocks
 * stays within the bound; a greedy bisimulation keeps \f$ h \f$ exact, and
 * any partition keeps it admissible. After each product the states not
 * reachable from the initial one or that reach no goal are pruned (the
 * latter are dead ends). The variables that are not ancestors of the goal
 * in the causal graph are left out, which is a further abstraction.
 *
 * A state of the task is mapped to its abstract state by a table per merge
 * step, from the pair ( abstract state so far , value of the variable
 * merged ) to the abstract state after the step. */

class SASMergeAndShrink {

 public:

 //! [planning-sig-mas]
 /// constructor, building the abstraction of the task, of at most
 /// max_states states

 SASMergeAndShrink( const SASTask & task , std::size_t max_states );
 //! [planning-sig-mas]

 /// the distance of the abstract state of the values s, one per variable;
 /// dead_end if it reaches no abstract goal

 [[nodiscard]] int operator()( const std::vector< int > & s ) const {
  int a = 0;
  for( std::size_t k = 0 ; k < v_vars.size() ; ++k ) {
   a = v_table[ k ][ std::size_t( a ) * std::size_t( v_range[ k ] ) +
                     std::size_t( s[ std::size_t( v_vars[ k ] ) ] ) ];
   if( a < 0 )
    return( dead_end );
   }
  return( v_dist[ std::size_t( a ) ] );
  }

 /// the variables in the order of the merge

 [[nodiscard]] const std::vector< int > & variables( void ) const noexcept {
  return( v_vars );
  }

 /// the number of abstract states

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_dist.size() );
  }

 private:

 std::vector< int > v_vars;                 ///< the variables merged
 std::vector< int > v_range;                ///< their ranges
 std::vector< std::vector< int > > v_table; ///< the tables of the steps
 std::vector< int > v_dist;                 ///< the abstract distances

 };  // end( class( SASMergeAndShrink ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SASHeuristics */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SASHeuristics.h ------------------------*/
/*--------------------------------------------------------------------------*/
