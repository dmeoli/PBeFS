/*--------------------------------------------------------------------------*/
/*------------------------ File PartialExpansion.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class templates PEAStar and EPEAStar, Partial
 * Expansion A* [see Yoshizumi, Miura and Ishida, AAAI 2000] and Enhanced
 * Partial Expansion A* [see Felner et al., AAAI 2012, and Goldenberg et
 * al., JAIR 2014].
 *
 * Both keep with a node \f$ n \f$ a stored value \f$ F( n ) \geq f( n ) \f$,
 * \f$ f( n ) \f$ at first, by which OPEN is ordered. When \f$ n \f$ is
 * selected, only its successors \f$ c \f$ with \f$ f( c ) \leq F( n ) +
 * C \f$ (PEA*, C being a cutoff) or with \f$ f( c ) = F( n ) \f$ (EPEA*)
 * are put in OPEN, and \f$ n \f$ goes back to OPEN with \f$ F( n ) \f$ the
 * least f of the others, if any, or is closed otherwise:
 *
 * \code{.unparsed}
 *  1  procedure PEA*(s₀)                     ▷ EPEA* at the lines marked E
 *  2     F(s₀) ← f(s₀);  push s₀ in OPEN, by F and then by larger g
 *  3     while OPEN ≠ ∅ do
 *  4        n ← the first node of OPEN, removed;  if n is a goal then
 *              return the path to n
 *  5        generate the successors c of n         ▷ E: those with
 *              f(c) − f(n) = F(n) − f(n) only, by the operator selection
 *              function, which also gives the least larger Δf
 *  6        for each c with f(c) ≤ F(n) + C do
 *  7           if c is new or reached through a cheaper path then
 *                 F(c) ← f(c);  push c in OPEN
 *  8        F' ← the least f(c) > F(n) + C of the others (E: f(n) plus
 *              the least larger Δf), ∞ if none
 *  9        if F' < ∞ then F(n) ← F';  push n in OPEN again
 * 10        else close n
 * \endcode
 *
 * With an admissible h the first goal selected is optimal, since a node is
 * in OPEN with \f$ F( n ) \f$ at most the f of every successor not yet put
 * there, so that the least F of OPEN is a lower bound of the cost of the
 * paths not yet explored, as the least f is in A*; the closed node is
 * reopened if reached through a cheaper path. PEA* saves memory, the
 * surplus nodes (those with \f$ f > C^* \f$) never entering OPEN, at the
 * price of generating the successors of a node at each of its
 * selections; EPEA* saves also that time, generating each successor at
 * most once per selection of \f$ n \f$ and only when it enters OPEN, on a
 * PartialExpansionDomain. The counter expanded counts the selections, and
 * reexpanded those of a node selected before.
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

#ifndef __PBeFS_PartialExpansion
 #define __PBeFS_PartialExpansion
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

#include "ADS/priority/DaryHeap.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS PartialOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of PEAStar

struct PartialOptions : public SearchOptions {

 /// the cutoff C of PEA*: the successors with \f$ f \leq F( n ) + C \f$
 /// enter OPEN at a selection of n
 double cutoff = 0;

 };  // end( struct( PartialOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS PartialExpansionBase -----------------------*/
/*--------------------------------------------------------------------------*/
/// what PEAStar and EPEAStar share: the nodes, OPEN and CLOSED
/** Expand( n , run ) of the derived class D puts the chosen successors of
 * n in OPEN by run.consider() and returns the new F of n, infinite if n
 * is to be closed. */

template< Domain D , class Derived >
class PartialExpansionBase {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

 /// a node, with its stored value F
 struct PNode {
  State state;
  Cost g;
  Cost h;
  PNode * parent;
  bool closed;
  bool expanded;
  std::uint8_t op;
  Cost F;
  };

 /// an entry of OPEN
 struct Item {
  Cost F;
  Cost g;
  PNode * node;
  };

 /// smaller F, then larger g
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.F < b.F ) || ( ( a.F == b.F ) && ( a.g > b.g ) ) );
   }
  };

 /// the work space of a search
 struct Run {
  NodePool< PNode > pool;
  NodeTable< PNode > table;
  ADS::DaryHeap< Item , Better > open;
  SearchStats * st;

  explicit Run( std::size_t slots ) : table( slots ) {}

  /// the successor t of n, of edge cost c and heuristic ht, enters OPEN
  /// if new or reached through a cheaper path
  void consider( PNode * n , const State & t , Cost c , Cost ht ,
                 std::uint8_t op , std::uint64_t hash ) {
   const Cost g = n->g + c;
   auto [ m , fresh ] = table.find_or_insert( t , hash , [ & ] {
    PNode * x = pool.allocate();
    *x = PNode{ t , g , ht , n , false , false , op , g + ht };
    return( x );
    } );
   if( fresh ) {
    open.push( Item{ m->F , g , m } );
    return;
    }
   if( g < m->g ) {
    m->g = g;
    m->parent = n;
    m->op = op;
    m->closed = false;
    m->F = g + m->h;
    open.push( Item{ m->F , g , m } );
    }
   else
    ++st->duplicates;
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 PartialExpansionBase( const D & domain , const SearchOptions & options )
  : f_domain( domain ) , f_base( options ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the search of the file comment

 SolutionT search( void ) const {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  Run r( f_base.expected_nodes
         ? NodeTable< PNode >::slots_for( f_base.expected_nodes ) : 1024 );
  r.st = &st;
  const State s0 = f_domain.initial_state();
  const Cost h0 = Cost( f_domain.heuristic( s0 ) );
  PNode * root = r.table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   PNode * x = r.pool.allocate();
   *x = PNode{ s0 , Cost( 0 ) , h0 , nullptr , false , false , no_op , h0 };
   return( x );
   } ).first;
  r.open.push( Item{ root->F , root->g , root } );
  const PNode * goal = nullptr;
  while( ! r.open.empty() ) {
   const Item e = r.open.pop();
   PNode * n = e.node;
   if( n->closed || ( n->g != e.g ) || ( n->F != e.F ) )
    continue;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( f_domain.is_goal( n->state ) ) {
    goal = n;
    break;
    }
   const Cost F = static_cast< const Derived & >( *this ).expand( n , r );
   if( F == infinite_cost< Cost >() )
    n->closed = true;
   else {
    n->F = F;
    r.open.push( Item{ F , n->g , n } );
    }
   if( f_base.node_limit && ( r.pool.size() >= f_base.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }
  if( goal ) {
   sol.status = Status::solved;
   sol.cost = goal->g;
   sol.path = trace_path( goal );
   }
  st.stored = r.pool.size();
  st.evaluated = st.generated + 1;
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_base;                               ///< the options

 };  // end( class( PartialExpansionBase ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS PEAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Partial Expansion A*, with the cutoff C of PartialOptions
/** The search of the file comment on any Domain: every selection of a node
 * generates all its successors (counted in generated), and those above the
 * cutoff are discarded. */

template< Domain D >
class PEAStar : public PartialExpansionBase< D , PEAStar< D > > {

 using Base = PartialExpansionBase< D , PEAStar< D > >;
 friend Base;

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using typename Base::State;
 using typename Base::Cost;
 using typename Base::SolutionT;

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit PEAStar( const D & domain , const PartialOptions & options = {} )
  : Base( domain , options ) , f_cutoff( options.cutoff ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const { return( this->search() ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// puts in OPEN the successors of n with f at most F( n ) + C, returning
 /// the least f of the others
 Cost expand( typename Base::PNode * n , typename Base::Run & r ) const {
  Cost next = infinite_cost< Cost >();
  const double top = double( n->F ) + f_cutoff;
  for_each_child( this->f_domain , n->state , n->h , n->op ,
                  [ & ]( const State & t , Cost c , Cost ht ,
                         std::uint8_t op ) {
   ++r.st->generated;
   const Cost f = n->g + c + ht;
   if( double( f ) <= top )
    r.consider( n , t , c , ht , op , this->f_domain.hash( t ) );
   else
    if( f < next )
     next = f;
   } );
  return( next );
  }

 double f_cutoff;                                    ///< the cutoff C

 };  // end( class( PEAStar ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS EPEAStar -----------------------------*/
/*--------------------------------------------------------------------------*/
/// Enhanced Partial Expansion A*, on the operator selection function of a
/// PartialExpansionDomain
/** The search of the file comment: a selection of a node generates only
 * its successors of \f$ \Delta f = F( n ) - f( n ) \f$, which are counted
 * in generated, so that generated is the number of nodes put in OPEN or
 * found duplicates. */

template< PartialExpansionDomain D >
class EPEAStar : public PartialExpansionBase< D , EPEAStar< D > > {

 using Base = PartialExpansionBase< D , EPEAStar< D > >;
 friend Base;

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using typename Base::State;
 using typename Base::Cost;
 using typename Base::SolutionT;

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit EPEAStar( const D & domain , const SearchOptions & options = {} )
  : Base( domain , options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const { return( this->search() ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// puts in OPEN the successors of n with \f$ \Delta f = F( n ) - f( n )
 /// \f$, returning f( n ) plus the least larger \f$ \Delta f \f$
 Cost expand( typename Base::PNode * n , typename Base::Run & r ) const {
  const Cost f = n->g + n->h;
  const Cost next = Cost( this->f_domain.for_each_successor_delta(
   n->state , n->h , n->op , Cost( n->F - f ) ,
   [ & ]( const State & t , Cost c , Cost ht , std::uint8_t op ) {
    ++r.st->generated;
    r.consider( n , t , c , ht , op , this->f_domain.hash( t ) );
    } ) );
  if( next == infinite_cost< Cost >() )
   return( next );
  return( f + next );
  }

 };  // end( class( EPEAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with PEAStar

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > pea_star(
 const D & domain , const PartialOptions & options = {} )
{
 return( PEAStar< D >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a PartialExpansionDomain with EPEAStar

template< PartialExpansionDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > epea_star(
 const D & domain , const SearchOptions & options = {} )
{
 return( EPEAStar< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PartialExpansion */

/*--------------------------------------------------------------------------*/
/*---------------------- End File PartialExpansion.h -----------------------*/
/*--------------------------------------------------------------------------*/
