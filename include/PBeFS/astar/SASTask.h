/*--------------------------------------------------------------------------*/
/*----------------------------- File SASTask.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class SASTask, a classical planning task in the
 * SAS+ formalism of Bäckström and Nebel (Computational Intelligence,
 * 1995), read from the file output.sas written by the translator of Fast
 * Downward [see Helmert, JAIR 2006], and for the class StripsTask, its
 * view as a STRIPS task over the facts \f$ v = d \f$.
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

#ifndef __PBeFS_SASTask
 #define __PBeFS_SASTask
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <iosfwd>

#include <string>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS SASTask -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a SAS+ planning task
/** The task has the variables \f$ v = 0 , \ldots , n - 1 \f$, the variable
 * \f$ v \f$ with the values \f$ \{ 0 , \ldots , | D_v | - 1 \} \f$, an
 * initial state assigning a value to each variable, a goal given by a
 * partial assignment, and operators. An operator \f$ o \f$ has the prevail
 * conditions, the values that some variables must have and keep, and the
 * effects \f$ ( v , p , q ) \f$, which set \f$ v \f$ to \f$ q \f$ and ask
 * \f$ v = p \f$ before, unless \f$ p = -1 \f$; it is applicable in the
 * states that satisfy all of them, and costs \f$ c( o ) \geq 0 \f$, or 1
 * if the task has no metric.
 *
 * read() parses the version 3 of the format output.sas of the translator
 * of Fast Downward; it throws std::invalid_argument on a malformed file and
 * on what the searches of the library do not support: the axioms (derived
 * variables), the conditional effects and the negative costs. */

class SASTask {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 //! [planning-sig-sastask]
 /// the fact \f$ v = d \f$
 struct Fact {
  int var;                                    ///< the variable
  int val;                                    ///< the value
  };

 /// the effect that sets var to post, asking var = pre if pre >= 0
 struct Effect {
  int var;                                    ///< the variable
  int pre;                                    ///< the value before, or -1
  int post;                                   ///< the value after
  };

 /// an operator
 struct Operator {
  std::string name;                           ///< the name
  std::vector< Fact > prevail;                ///< the prevail conditions
  std::vector< Effect > effects;              ///< the effects
  int cost;                                   ///< the cost
  };

/*----------------------------- READING A TASK -----------------------------*/
 /// the task in the stream, in the format output.sas

 [[nodiscard]] static SASTask read( std::istream & in );

 /// the task in the file of the given name

 [[nodiscard]] static SASTask read_file( const std::string & name );

/*---------------------------- THE PRECONDITIONS ---------------------------*/
 /// the preconditions of the operator o: its prevail conditions and the
 /// values asked by its effects, sorted by variable

 [[nodiscard]] std::vector< Fact > preconditions( std::size_t o ) const;
 //! [planning-sig-sastask]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the values s, one per variable, satisfy the goal

 [[nodiscard]] bool satisfies_goal( const std::vector< int > & s ) const {
  for( const Fact & f : goal )
   if( s[ f.var ] != f.val )
    return( false );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of variables

 [[nodiscard]] std::size_t num_vars( void ) const noexcept {
  return( range.size() );
  }

/*--------------------------- PUBLIC DATA MEMBERS --------------------------*/

 //! [planning-sig-sastask-data]
 std::vector< std::string > var_names;   ///< the names of the variables
 std::vector< int > range;               ///< | D_v | of each variable
 std::vector< std::vector< std::string > > values;
                                         ///< the names of the values
 std::vector< std::vector< Fact > > mutexes;   ///< the mutex groups
 std::vector< int > init;                ///< the initial value of each v
 std::vector< Fact > goal;               ///< the goal
 std::vector< Operator > ops;            ///< the operators
 bool use_costs = false;                 ///< true if the costs are read
 //! [planning-sig-sastask-data]

 };  // end( class( SASTask ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS StripsTask ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the STRIPS view of a SASTask
/** The fact \f$ v = d \f$ is the proposition numbered offset[ v ] + d; an
 * operator has the preconditions, the add list (its effects \f$ v = q \f$)
 * and the delete list (\f$ v = p \f$ for an effect with \f$ p \geq 0 \f$,
 * all the other values of \f$ v \f$ otherwise). */

class StripsTask {

 public:

 //! [planning-sig-strips]
 /// a STRIPS operator, the facts sorted
 struct Operator {
  std::vector< int > pre;                      ///< the preconditions
  std::vector< int > add;                      ///< the add list
  std::vector< int > del;                      ///< the delete list
  int cost;                                    ///< the cost
  };

 /// the view of the task t

 explicit StripsTask( const SASTask & t );
 //! [planning-sig-strips]

 /// the proposition of the fact v = d

 [[nodiscard]] int fact( int v , int d ) const noexcept {
  return( offset[ std::size_t( v ) ] + d );
  }

 //! [planning-sig-strips-data]
 int num_facts = 0;                            ///< the propositions
 std::vector< int > offset;                    ///< the first one of each v
 std::vector< Operator > ops;                  ///< the operators
 std::vector< int > init;                      ///< the initial facts
 std::vector< int > goal;                      ///< the goal facts
 //! [planning-sig-strips-data]

 };  // end( class( StripsTask ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SASTask */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File SASTask.h ---------------------------*/
/*--------------------------------------------------------------------------*/
