/*--------------------------------------------------------------------------*/
/*--------------------------- File SymbolicSAS.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SymbolicSAS, which gives a SAS+
 * planning task as a SymbolicDomain [see SymbolicDomain.h].
 *
 * The task is any type with the public members of a SAS+ task in the
 * format of the translator of Fast Downward [see Helmert, JAIR 2006]: range,
 * the number of values of each variable; init, the initial value of each;
 * goal, the facts of the goal, each with the members var and val; ops, the
 * operators, each with prevail (facts), effects (each with var, pre, -1
 * if none, and post) and cost; and use_costs, false if every operator
 * costs 1. The operators become operators of assignments, their
 * preconditions the prevail conditions and the values that the effects
 * ask; conditional effects and axioms are not part of the format taken.
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

#ifndef __PBeFS_SymbolicSAS
 #define __PBeFS_SymbolicSAS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include <vector>

#include "PBeFS/symbolic/SymbolicDomain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS SymbolicSAS ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a SAS+ task as a SymbolicDomain
/** The task must outlive the domain; the values of a state are those of
 * the variables of the task, so that a path is a sequence of states of the
 * task. */

template< class Task >
class SymbolicSAS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using cost_type = int;                         ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the domain of the task t

 explicit SymbolicSAS( const Task & t ) : f_t( t ) {
  if( t.init.size() != t.range.size() )
   throw( std::invalid_argument( "SymbolicSAS: an initial value for each "
                                 "variable" ) );
  }

/*-------------------------- METHODS OF THE DOMAIN -------------------------*/

 [[nodiscard]] std::size_t num_variables( void ) const noexcept {
  return( f_t.range.size() );
  }

 [[nodiscard]] std::size_t range( std::uint32_t v ) const {
  return( std::size_t( f_t.range[ v ] ) );
  }

 [[nodiscard]] std::vector< std::uint32_t > initial_values( void ) const {
  std::vector< std::uint32_t > x;
  for( auto d : f_t.init )
   x.push_back( std::uint32_t( d ) );
  return( x );
  }

 template< class F >
 void for_each_goal( F && f ) const {
  std::vector< SymbolicFact > p;
  for( const auto & g : f_t.goal )
   p.push_back( { std::uint32_t( g.var ) , std::uint32_t( g.val ) } );
  f( p );
  }

 template< class F >
 void for_each_operator( F && f ) const {
  for( const auto & op : f_t.ops ) {
   SymbolicOperator< cost_type > o;
   o.cost = f_t.use_costs ? cost_type( op.cost ) : cost_type( 1 );
   for( const auto & p : op.prevail )
    o.pre.push_back( { std::uint32_t( p.var ) , std::uint32_t( p.val ) } );
   for( const auto & e : op.effects ) {
    if( e.pre >= 0 )
     o.pre.push_back( { std::uint32_t( e.var ) , std::uint32_t( e.pre ) } );
    o.effects.push_back( SymbolicEffect::assign( std::uint32_t( e.var ) ,
                                                 std::uint32_t( e.post ) ) );
    }
   f( o );
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const Task & f_t;                              ///< the task

 };  // end( class( SymbolicSAS ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SymbolicSAS */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SymbolicSAS.h -------------------------*/
/*--------------------------------------------------------------------------*/
