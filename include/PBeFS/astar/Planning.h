/*--------------------------------------------------------------------------*/
/*---------------------------- File Planning.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template Planning, which implements the Domain
 * concept [see Domain.h] for classical planning on a SAS+ task [see
 * SASTask.h], with the heuristics of SASHeuristics.h, the features of
 * Zobrist hashing, an abstraction of the features for Abstract Zobrist
 * hashing, the nblocks of PBNF and the operators of OperatorDomain, and
 * for its base PlanningBase, which holds what does not depend on the size
 * of the state.
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

#ifndef __PBeFS_Planning
 #define __PBeFS_Planning
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cstddef>

#include <cstdint>

#include <memory>

#include <mutex>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/SASHeuristics.h"

#include "PBeFS/astar/SASTask.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the heuristic of a Planning domain

enum class PlanningHeuristic {
 blind ,        ///< 0, admissible
 goal_count ,   ///< the goal facts that do not hold, not admissible
 hmax ,         ///< h^max, admissible
 hadd ,         ///< h^add, not admissible
 hff ,          ///< h^FF, not admissible
 pdb ,          ///< a pattern database, admissible
 lmcut          ///< the landmark-cut heuristic, admissible
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the heuristic of the given name: blind, goalcount, hmax, hadd, hff, pdb
/// or lmcut; throws std::invalid_argument on another name

[[nodiscard]] PlanningHeuristic planning_heuristic( const std::string & n );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the options of a Planning domain

struct PlanningOptions {

 PlanningHeuristic heuristic = PlanningHeuristic::lmcut;  ///< the heuristic

 /// the pattern of the pattern database, the greedy one of
 /// SASPatternDatabase if empty
 std::vector< int > pattern;

 /// the most abstract states of the greedy pattern
 std::size_t pdb_size = std::size_t( 1 ) << 16;

 /// the fraction of the variables, those of largest fluency, that
 /// Abstract Zobrist hashing leaves out [see PlanningBase]; 0 leaves out
 /// none
 double ignored = 0.3;

 /// the most nblocks of the abstraction of PBNF and of AHDA*
 std::size_t max_nblocks = 10000;

 /// if positive, the largest fraction of the facts that the variables of
 /// the abstraction of the nblocks may hold
 double nblock_facts = 0;

 };  // end( struct( PlanningOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS PlanningBase --------------------------*/
/*--------------------------------------------------------------------------*/
/// the part of a Planning domain that does not depend on the state type
/** It packs the variables into words, each in a field of the bits its
 * range needs and none across two words, builds the successor generator,
 * the heuristic and the abstractions, and evaluates the heuristic on the
 * values of the variables. The successor generator is a decision tree over
 * the variables in the order of their indices [see Helmert, JAIR 2006]: a
 * node tests a variable, has a child for each of its values and one for
 * the operators that ask nothing of it, and lists the operators whose
 * preconditions have all been tested on the way to it.
 *
 * The features of Zobrist hashing are the facts \f$ v = d \f$. The
 * abstraction of Abstract Zobrist hashing is that of FAZHDA* [see Jinnai
 * and Fukunaga, ICAPS 2016, section 3, and JAIR 2017, section 4.2]: the
 * fraction PlanningOptions::ignored of the variables of largest fluency
 * (the fraction of the operators that change them) is left out, all their
 * facts mapped to one abstract feature, so that their moves never change
 * the owner of a state, and the values of each other variable are split in
 * two abstract features by GreedyAFG, i.e., its undirected domain
 * transition graph is bisected by starting from a value of least degree
 * and adding the value with the most edges into the first half until it
 * holds half of the values. With ignored = 0 the abstraction is that of
 * GAZHDA*.
 *
 * The nblocks of PBNF, which are also the abstract states of AHDA*, are the
 * states of the projection onto a set of variables, chosen by the greedy
 * abstraction of structured duplicate detection [see Zhou and Hansen,
 * AAAI 2006, and Burns et al., JAIR 2010]: the variable added at each step
 * is the one whose abstract graph has the least largest number of
 * successors of a node (the larger graph on a tie), as long as the graph
 * has at most PlanningOptions::max_nblocks nodes and, if
 * PlanningOptions::nblock_facts is positive, the variables hold at most
 * that fraction of the facts, which is the dynamic threshold of DAHDA*
 * [see Jinnai and Fukunaga, ICAPS 2016, section 5]. The abstract graph is
 * then taken undirected, its edges in both directions, since the
 * interference scopes of SafePBNF count the successors of an nblock among
 * its predecessors; the abstraction is built at its first use.
 *
 * The operator that generated a state is given as a code of 8 bits:
 * \f$ o' \f$ is an inverse of \f$ o \f$ if it changes the same variables,
 * sets each of them back to the value that \f$ o \f$ asks before, and its
 * preconditions hold after \f$ o \f$ whatever the state, so that
 * \f$ o' ( o ( s ) ) = s \f$. The first 255 operators with an inverse get
 * the codes \f$ 0 , \ldots , 254 \f$, the others no_op, and a state
 * generated by a coded operator skips the inverses of that operator
 * [see OperatorDomain]. */

class PlanningBase {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, packing the variables into at most max_words words of 64
 /// bits; throws std::invalid_argument if they do not fit

 PlanningBase( SASTask task , const PlanningOptions & o ,
               std::size_t max_words );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the words of 64 bits that the state of the task t needs

 [[nodiscard]] static std::size_t words_for( const SASTask & t );

/*------------------------------ THE TASK ----------------------------------*/

 /// the task

 [[nodiscard]] const SASTask & task( void ) const noexcept {
  return( *p_task );
  }

 /// its STRIPS view

 [[nodiscard]] const StripsTask & strips( void ) const noexcept {
  return( *p_strips );
  }

 /// the words of the packed state actually used

 [[nodiscard]] std::size_t words( void ) const noexcept { return( f_words ); }

 /// the heuristic of the values s, one per variable; dead_end if no goal
 /// can be reached

 [[nodiscard]] int evaluate( const std::vector< int > & s ) const;

 /// the pattern database, if the heuristic is one

 [[nodiscard]] const SASPatternDatabase * pdb( void ) const noexcept {
  return( p_pdb.get() );
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: the facts

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( std::size_t( p_strips->num_facts ) );
  }

 /// the number of abstract features

 [[nodiscard]] std::size_t num_abstract_features( void ) const noexcept {
  return( f_num_abstract );
  }

 /// the abstract feature of the fact i

 [[nodiscard]] std::size_t abstract_feature( std::size_t i ) const noexcept
 {
  return( v_abstract[ i ] );
  }

/*-------------------------------- NBLOCKS ---------------------------------*/
 /// the number of nblocks

 [[nodiscard]] std::size_t num_nblocks( void ) const {
  return( nblocks().first.size() - 1 );
  }

 /// the variables of the abstraction of the nblocks

 [[nodiscard]] const std::vector< int > & nblock_vars( void ) const {
  return( nblocks().vars );
  }

 /// calls f( b' ) for each nblock b' reached from b by one operator, or
 /// reaching b by one operator

 template< class F >
 void for_each_nblock_successor( std::size_t b , F && f ) const {
  const NBlocks & x = nblocks();
  for( auto k = x.first[ b ] ; k < x.first[ b + 1 ] ; ++k )
   f( std::size_t( x.succ[ k ] ) );
  }

/*------------------------------- OPERATORS --------------------------------*/
 /// the code of the operator o [see OperatorDomain], no_op if none

 [[nodiscard]] std::uint8_t code( std::size_t o ) const noexcept {
  return( v_code[ o ] );
  }

 /// true if the operator o is an inverse of the operator of code c

 [[nodiscard]] bool undoes( std::uint8_t c , std::uint32_t o ) const
  noexcept {
  if( c == no_op )
   return( false );
  for( auto k = v_inv_first[ c ] ; k < v_inv_first[ c + 1u ] ; ++k )
   if( v_inv[ k ] == o )
    return( true );
  return( false );
  }

/*----------------------- PROTECTED PART OF THE CLASS ----------------------*/

 protected:

 /// a node of the successor generator
 struct GenNode {
  int var;                      ///< the variable tested, -1 at a leaf
  std::uint32_t ops_begin;      ///< the operators listed here, in v_gen_ops
  std::uint32_t ops_end;
  std::uint32_t child_begin;    ///< the children by value, in v_gen_child
  int dont_care;                ///< the child asking nothing, or -1
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( o ) for each operator o applicable where the variable v has
 /// the value get( v ), from the node n of the generator

 template< class Get , class F >
 void visit( int n , const Get & get , F & f ) const {
  const GenNode & x = v_gen[ std::size_t( n ) ];
  for( auto k = x.ops_begin ; k < x.ops_end ; ++k )
   f( v_gen_ops[ k ] );
  if( x.var < 0 )
   return;
  const int c = v_gen_child[ x.child_begin + unsigned( get( x.var ) ) ];
  if( c >= 0 )
   visit( c , get , f );
  if( x.dont_care >= 0 )
   visit( x.dont_care , get , f );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 std::shared_ptr< const SASTask > p_task;          ///< the task
 std::shared_ptr< const StripsTask > p_strips;     ///< its STRIPS view
 PlanningHeuristic f_heuristic;                    ///< the heuristic
 std::shared_ptr< const RelaxedHeuristic > p_relaxed;
                                                   ///< h^max, h^add, h^FF
 std::shared_ptr< const LMCut > p_lmcut;           ///< LM-cut
 std::shared_ptr< const SASPatternDatabase > p_pdb;
                                                   ///< the PDB

 std::size_t f_words;                              ///< the words used
 std::vector< std::uint32_t > v_word;              ///< the word of each v
 std::vector< std::uint32_t > v_shift;             ///< its first bit
 std::vector< std::uint64_t > v_mask;              ///< its mask, unshifted
 std::vector< std::vector< int > > v_word_vars;    ///< the vars of a word
 std::vector< std::uint64_t > v_init;              ///< the packed init

 std::vector< GenNode > v_gen;                     ///< the generator
 std::vector< std::uint32_t > v_gen_ops;           ///< its operators
 std::vector< int > v_gen_child;                   ///< its children

 std::vector< std::uint32_t > v_eff_first;         ///< CSR of the effects
 std::vector< std::uint32_t > v_eff_var;           ///< the variable
 std::vector< std::uint64_t > v_eff_post;          ///< the value after
 std::vector< int > v_cost;                        ///< the cost of each o

 std::vector< std::uint8_t > v_code;               ///< the code of each o
 std::vector< std::uint32_t > v_inv_first;         ///< CSR of the inverses
 std::vector< std::uint32_t > v_inv;               ///< of each code

 std::vector< std::size_t > v_abstract;            ///< by fact
 std::size_t f_num_abstract;                       ///< abstract features

 /// the abstraction of the nblocks, built at its first use
 struct NBlocks {
  std::once_flag once;                             ///< for the building
  std::size_t max_nblocks;                         ///< the options
  double max_facts;
  std::vector< int > vars;                         ///< the variables
  std::vector< std::size_t > mult;                 ///< their multipliers
  std::vector< std::size_t > first;                ///< CSR of the graph
  std::vector< std::uint32_t > succ;
  };

 /// the abstraction of the nblocks, built if not yet

 [[nodiscard]] const NBlocks & nblocks( void ) const {
  std::call_once( p_nb->once , [ this ] { build_nblocks(); } );
  return( *p_nb );
  }

 /// builds the abstraction of the nblocks

 void build_nblocks( void ) const;

 std::shared_ptr< NBlocks > p_nb;                  ///< the nblocks

 };  // end( class( PlanningBase ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Planning -----------------------------*/
/*--------------------------------------------------------------------------*/
/// classical planning on a SAS+ task, the state in W words of 64 bits
/** A state is the array of W words into which the values of the variables
 * are packed [see PlanningBase], compared and hashed as a whole; the cost
 * of a move is that of its operator, and a successor on which the
 * heuristic finds a dead end is not generated. The class is a Domain, an
 * IncrementalFeatureDomain, an AbstractFeatureDomain, an
 * AbstractStateDomain, an NBlockDomain and an OperatorDomain, so that every
 * search and every distribution of the library runs on it. */

template< std::size_t W = 4 >
class Planning : public PlanningBase {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::array< std::uint64_t , W >;  ///< a packed state
 using cost_type = int;                               ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the task and the options; throws
 /// std::invalid_argument if the state does not fit in W words

 explicit Planning( SASTask task , const PlanningOptions & o = {} )
  : PlanningBase( std::move( task ) , o , W ) {
  for( std::size_t i = 0 ; i < W ; ++i )
   f_init[ i ] = i < f_words ? v_init[ i ] : 0;
  }

/*----------------------------- PACKED STATES ------------------------------*/
 /// the value of the variable v in the state s

 [[nodiscard]] int value( const state_type & s , int v ) const noexcept {
  const auto x = std::size_t( v );
  return( int( ( s[ v_word[ x ] ] >> v_shift[ x ] ) & v_mask[ x ] ) );
  }

 /// the values of the variables in s

 [[nodiscard]] std::vector< int > unpack( const state_type & s ) const {
  std::vector< int > x( p_task->num_vars() );
  for( std::size_t v = 0 ; v < x.size() ; ++v )
   x[ v ] = value( s , int( v ) );
  return( x );
  }

 /// the state of the values x, one per variable

 [[nodiscard]] state_type pack( const std::vector< int > & x ) const {
  state_type s{};
  for( std::size_t v = 0 ; v < x.size() ; ++v )
   s[ v_word[ v ] ] |= std::uint64_t( x[ v ] ) << v_shift[ v ];
  return( s );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_init );
  }

 [[nodiscard]] bool is_goal( const state_type & s ) const noexcept {
  for( const auto & g : p_task->goal )
   if( value( s , g.var ) != g.val )
    return( false );
  return( true );
  }

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const noexcept {
  std::uint64_t h = ADS::splitmix64( s[ 0 ] );
  for( std::size_t i = 1 ; i < W ; ++i )
   h = ADS::splitmix64( h ^ s[ i ] );
  return( h );
  }

 /// the heuristic of the options; dead_end if no goal is reachable

 [[nodiscard]] int heuristic( const state_type & s ) const {
  thread_local std::vector< int > x;
  x.resize( p_task->num_vars() );
  for( std::size_t v = 0 ; v < x.size() ; ++v )
   x[ v ] = value( s , int( v ) );
  return( evaluate( x ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( t , c , h( t ) ) for each successor t of s that is not a
 /// dead end

 template< class F >
 void for_each_successor( const state_type & s , int h , F && f ) const {
  for_each_successor( s , h , no_op , [ & ]( const state_type & t , int c ,
                                             int ht , std::uint8_t ) {
   f( t , c , ht );
   } );
  }

 /// calls f( t , c , h( t ) , op_t ) for each successor t of s that is not
 /// a dead end, except those of the inverses of the operator of code op

 template< class F >
 void for_each_successor( const state_type & s , int , std::uint8_t op ,
                          F && f ) const {
  const auto get = [ & ]( int v ) { return( value( s , v ) ); };
  auto each = [ & ]( std::uint32_t o ) {
   if( undoes( op , o ) )
    return;
   const state_type t = apply( s , o );
   const int ht = heuristic( t );
   if( ht != dead_end )
    f( t , v_cost[ o ] , ht , v_code[ o ] );
   };
  visit( 0 , get , each );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the state reached from s by the operator o, applicable in s

 [[nodiscard]] state_type apply( const state_type & s , std::uint32_t o )
  const noexcept {
  state_type t = s;
  for( auto k = v_eff_first[ o ] ; k < v_eff_first[ o + 1 ] ; ++k ) {
   const auto v = v_eff_var[ k ];
   auto & w = t[ v_word[ v ] ];
   w = ( w & ~( v_mask[ v ] << v_shift[ v ] ) ) |
       ( v_eff_post[ k ] << v_shift[ v ] );
   }
  return( t );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the names of the operators of a cheapest plan along the states of
 /// path, which must be a path of the Domain

 [[nodiscard]] std::vector< std::string > plan(
                                 const std::vector< state_type > & path ) const
 {
  std::vector< std::string > p;
  for( std::size_t i = 1 ; i < path.size() ; ++i ) {
   const state_type & s = path[ i - 1 ];
   long best = -1;
   const auto get = [ & ]( int v ) { return( value( s , v ) ); };
   auto each = [ & ]( std::uint32_t o ) {
    if( ( apply( s , o ) == path[ i ] ) &&
        ( ( best < 0 ) || ( v_cost[ o ] < v_cost[ std::size_t( best ) ] ) ) )
     best = long( o );
    };
   visit( 0 , get , each );
   if( best < 0 )
    throw( std::invalid_argument( "Planning::plan: not a path" ) );
   p.push_back( p_task->ops[ std::size_t( best ) ].name );
   }
  return( p );
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// calls f( i ) for the fact i of each variable in s

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const {
  for( std::size_t v = 0 ; v < v_word.size() ; ++v )
   f( std::size_t( p_strips->offset[ v ] + value( s , int( v ) ) ) );
  }

 /// calls f( i ) for each fact i that holds in exactly one of s and t

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const {
  for( std::size_t w = 0 ; w < f_words ; ++w ) {
   const std::uint64_t d = s[ w ] ^ t[ w ];
   if( ! d )
    continue;
   for( int v : v_word_vars[ w ] )
    if( ( d >> v_shift[ std::size_t( v ) ] ) & v_mask[ std::size_t( v ) ] ) {
     const int o = p_strips->offset[ std::size_t( v ) ];
     f( std::size_t( o + value( s , v ) ) );
     f( std::size_t( o + value( t , v ) ) );
     }
   }
  }

/*------------------------------ ABSTRACTIONS ------------------------------*/
 /// the nblock of s: its index in the projection onto nblock_vars()

 [[nodiscard]] std::size_t nblock( const state_type & s ) const {
  const NBlocks & x = nblocks();
  std::size_t b = 0;
  for( std::size_t k = 0 ; k < x.vars.size() ; ++k )
   b += std::size_t( value( s , x.vars[ k ] ) ) * x.mult[ k ];
  return( b );
  }

 /// the abstract state of AHDA*: the nblock

 [[nodiscard]] std::uint64_t abstract_state( const state_type & s ) const
 {
  return( nblock( s ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 state_type f_init;                               ///< the initial state

 };  // end( class( Planning ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Planning */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Planning.h -------------------------*/
/*--------------------------------------------------------------------------*/
