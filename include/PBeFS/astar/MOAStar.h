/*--------------------------------------------------------------------------*/
/*----------------------------- File MOAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the sequential multi-objective searches: the class
 * template MOAStar, the framework of BOA*, EMOA* and LTMOA* with their
 * lazy variants, and the class template NAMOAStarDR, NAMOA* with
 * dimensionality reduction.
 *
 * MOAStar is the search framework of Ren et al. (SoCS 2022, Algorithm 1),
 * the generalization of BOA* of Hernández et al. (AIJ 2023): OPEN holds
 * labels, i.e., paths \f$ x \f$ to a state \f$ s( x ) \f$ of cost
 * \f$ g( x ) \f$, ordered by the lexicographic order of
 * \f$ f( x ) = g( x ) + h( s( x ) ) \f$, each state \f$ s \f$ has the
 * frontier \f$ G^{Tr}( s ) \f$ of the truncated g of the paths expanded
 * there, and \f$ S \f$ is the frontier of the truncated f of the solutions
 * found, with the notation of the page "Multi-objective search" of the
 * documentation:
 *
 * \code{.unparsed}
 *  1  procedure MOA*(s₀)
 *  2     OPEN ← {x₀}, with s(x₀) = s₀, g(x₀) = 0;  sols ← ∅
 *  3     while OPEN ≠ ∅ do
 *  4        x ← the first label of OPEN by the lexicographic order of f;
 *              remove it from OPEN
 *  5        if s(x) is not a goal and Tr(g(x)) ∈ dom G^Tr(s(x)) then
 *              continue                          ▷ the frontier check
 *  6        if Tr(f(x)) ∈ dom S then continue    ▷ the solution check
 *  7        if s(x) is a goal then
 *  8           add Tr(f(x)) to S, removing what it dominates;
 *              add x to sols;  continue
 *  9        add Tr(g(x)) to G^Tr(s(x)), removing what it dominates
 * 10        for each successor t of s(x) with no infinite h(t) do
 * 11           y ← the label of t with g(y) = g(x) + c(s(x), t), parent x
 * 12           if not lazy and (the check of line 5 or 6 prunes y) then
 *                 continue
 * 13           push y in OPEN
 * 14     return sols
 * \endcode
 *
 * where \f$ v \in \mathrm{dom}\, F \f$ means that a vector of F weakly
 * dominates v. The frontier, a template parameter, decides how lines 5, 6,
 * 8 and 9 are done: ScalarFrontier, the number \f$ g_2^{\min} \f$ of BOA*
 * for two objectives, ArrayFrontier, the array of LTMOA* (Hernández et
 * al., IJCAI 2023), TreeFrontier, the balanced tree of EMOA* (with the
 * check in logarithmic time of TOA* for three objectives) [see
 * MultiObjective.h]. With the option lazy the checks of line 12 are not
 * done, which is LazyLTMOA* with the array (and lazy EMOA* with the tree).
 * The ties of f go to the label pushed first, and the paths to a goal are
 * not expanded, since every path through a goal is weakly dominated by its
 * prefix to that goal.
 *
 * NAMOAStarDR is NAMOA* of Mandow and Pérez de la Cruz (JACM 2010) with the
 * dimensionality reduction of Pulido, Mandow and Pérez de la Cruz (C&OR
 * 2015, NAMOA*dr), for any number of objectives: besides \f$ G^{Tr}( s )
 * \f$, which it consults only when it generates a path, each state keeps
 * \f$ G_{op}( s ) \f$, the g of its labels in OPEN, against which a new
 * label is checked in full and which loses the labels that the new one
 * dominates, removed from OPEN as well; the label extracted is then never
 * dominated at its state, and only the solution check of line 6 is done
 * there.
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

#ifndef __PBeFS_MOAStar
 #define __PBeFS_MOAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cstddef>

#include <cstdint>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/MultiObjective.h"

#include "PBeFS/astar/NodeTable.h"

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

/// the frontier of BOA* for two objectives, that of LTMOA* otherwise
template< class Cost , std::size_t K >
using AutoFrontier = std::conditional_t< K == 2 , ScalarFrontier< Cost , 2 > ,
                                         ArrayFrontier< Cost , K > >;

/// the frontier of EMOA*, which for two objectives is that of BOA* [see
/// Ren et al., SoCS 2022]
template< class Cost , std::size_t K >
using EMOAFrontier = std::conditional_t< K == 2 ,
                                         ScalarFrontier< Cost , 2 > ,
                                         TreeFrontier< Cost , K > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the multi-objective searches, by name

enum class MOAlgorithm {
 boa ,          ///< BOA*, two objectives only
 emoa ,         ///< EMOA* (TOA* for three objectives)
 ltmoa ,        ///< LTMOA*
 namoa_dr       ///< NAMOA*dr
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// u + v, component by component

template< class T , std::size_t K >
[[nodiscard]] constexpr std::array< T , K > add_costs(
 const std::array< T , K > & u , const std::array< T , K > & v ) noexcept
{
 std::array< T , K > w;
 for( std::size_t i = 0 ; i < K ; ++i )
  w[ i ] = u[ i ] + v[ i ];
 return( w );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the states of the path that ends with the label x, following the
/// parents of the labels, each with the member rec, whose state is that of
/// the label

template< class Label >
[[nodiscard]] auto trace_labels( const Label * x )
{
 std::vector< std::remove_cvref_t< decltype( x->rec->state ) > > p;
 for( ; x ; x = x->parent )
  p.push_back( x->rec->state );
 std::reverse( p.begin() , p.end() );
 return( p );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the front of sol from the labels of the solutions found, in the order p,
/// with their paths if paths, and its status: limit if the search stopped
/// at its limit

template< class S , class Label , std::size_t K >
void finish_front( S & sol , const std::vector< Label * > & goals ,
                   bool limit , const std::array< unsigned , K > & p ,
                   bool paths )
{
 sol.status = limit ? Status::node_limit
                    : ( goals.empty() ? Status::unsolvable
                                      : Status::solved );
 for( const Label * x : goals ) {
  typename S::Path q;
  q.cost = unpermute( x->g , p );
  if( paths )
   q.states = trace_labels( x );
  sol.front.push_back( std::move( q ) );
  }
 std::sort( sol.front.begin() , sol.front.end() ,
            []( const auto & a , const auto & b ) {
             return( a.cost < b.cost );
             } );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS MOAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the sequential multi-objective A* of the framework of BOA*
/** The search of the file header, with the frontier F (ScalarFrontier for
 * BOA*, ArrayFrontier for LTMOA*, TreeFrontier for EMOA*) and the options
 * of MOSearchOptions; it returns a cost-unique Pareto-optimal set if the
 * heuristic is consistent [see the page "Multi-objective search"]. The
 * labels pruned when extracted are recycled. */

template< MultiObjectiveDomain D ,
          template< class , std::size_t > class F = AutoFrontier >
class MOAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr std::size_t K = objectives_of< D >;  ///< the objectives
 using State = mo_state_t< D >;                      ///< the state type
 using Cost = mo_scalar_t< D >;                      ///< a component
 using Vec = mo_cost_t< D >;                         ///< a vector of costs
 using Frontier = F< Cost , K >;                     ///< the frontier
 using SolutionT = pareto_solution_t< D >;           ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the MOAStar); throws
 /// std::invalid_argument if the order is not a permutation

 explicit MOAStar( const D & domain , const MOSearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_order( objective_order< K >( options.order ) ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Pareto front

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  SolutionT sol;
  MOSearchStats & st = sol.stats;
  NodePool< Rec > recs;
  NodeTable< Rec > table( f_options.expected_states
                          ? NodeTable< Rec >::slots_for(
                                            f_options.expected_states )
                          : 1024 );
  NodePool< Label > labels;
  std::vector< Label * > spare;                  // the labels to recycle
  MOHeap< Item , Better > open;
  Frontier sols;                                 // S
  std::vector< Label * > goals;                  // the solutions
  std::uint64_t seq = 0;

  const auto record = [ & ]( const State & s , const Vec & h ) {
   return( table.find_or_insert( s , f_domain.hash( s ) , [ & ] {
    Rec * r = recs.allocate();
    r->state = s;
    r->h = permute( h , f_order );
    r->goal = f_domain.is_goal( s );
    return( r );
    } ).first );
   };
  const auto label = [ & ]( Rec * r , const Vec & g , Label * parent ) {
   Label * x;
   if( spare.empty() )
    x = labels.allocate();
   else {
    x = spare.back();
    spare.pop_back();
    }
   *x = Label{ r , g , parent };
   return( x );
   };

  const State s0 = f_domain.initial_state();
  Rec * r0 = record( s0 , f_domain.heuristic( s0 ) );
  if( ! dead_end( r0->h ) ) {
   open.push( Item{ r0->h , label( r0 , Vec{} , nullptr ) , seq++ } );
   ++st.inserted;
   }

  bool limit = false;
  while( ! open.empty() ) {
   const Item it = open.pop();
   ++st.extracted;
   Label * x = it.label;
   Rec * r = x->rec;
   if( ( ! r->goal && r->front.dominated( x->g , st.checks ) ) ||
       sols.dominated( it.f , st.checks ) ) {
    spare.push_back( x );
    continue;
    }
   if( r->goal ) {
    sols.insert( it.f , st.checks );
    goals.push_back( x );
    ++st.solutions;
    continue;
    }
   r->front.insert( x->g , st.checks );
   ++st.expanded;

   f_domain.for_each_successor( r->state , unpermute( r->h , f_order ) ,
                                [ & ]( const State & t , const Vec & c ,
                                       const Vec & ht ) {
    ++st.generated;
    Rec * q = record( t , ht );
    if( dead_end( q->h ) )
     return;
    const Vec g = add_costs( x->g , permute( c , f_order ) );
    const Vec f = add_costs( g , q->h );
    if( ! f_options.lazy &&
        ( ( ! q->goal && q->front.dominated( g , st.checks ) ) ||
          sols.dominated( f , st.checks ) ) )
     return;
    open.push( Item{ f , label( q , g , x ) , seq++ } );
    ++st.inserted;
    } );

   if( f_options.node_limit && ( labels.size() >= f_options.node_limit ) ) {
    limit = true;
    break;
    }
   }

  finish_front( sol , goals , limit , f_order , f_options.paths );
  st.stored = labels.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 struct Rec;

 /// a path: the record of its last state, its g and its predecessor
 struct Label {
  Rec * rec;
  Vec g;
  Label * parent;
  };

 /// the record of a state: the state, its h in the order of the search,
 /// whether it is a goal, its frontier
 struct Rec {
  State state;
  Vec h;
  bool goal;
  Frontier front;
  };

 /// an entry of OPEN: f, the label and the order of the pushes
 struct Item {
  Vec f;
  Label * label;
  std::uint64_t seq;
  };

 /// the order of OPEN: the lexicographic order of f, then the first pushed
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.seq < b.seq ) ) );
   }
  };

 const D & f_domain;                                 ///< the Domain
 MOSearchOptions f_options;                          ///< the options
 std::array< unsigned , K > f_order;                 ///< the order

 };  // end( class( MOAStar ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS NAMOAStarDR ----------------------------*/
/*--------------------------------------------------------------------------*/
/// NAMOA* with dimensionality reduction
/** The search of the file header, with the frontier F for the truncated
 * \f$ G_{cl}( s ) \f$ (by default ScalarFrontier for two objectives and
 * ArrayFrontier for more) and the array \f$ G_{op}( s ) \f$ of the labels
 * of s in OPEN; a label removed from \f$ G_{op} \f$ because a new one
 * dominates it is marked, and dropped when it reaches the top of OPEN. The
 * option lazy is ignored, since the checks at the generation are those
 * that make the search what it is. */

template< MultiObjectiveDomain D ,
          template< class , std::size_t > class F = AutoFrontier >
class NAMOAStarDR {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr std::size_t K = objectives_of< D >;  ///< the objectives
 using State = mo_state_t< D >;                      ///< the state type
 using Cost = mo_scalar_t< D >;                      ///< a component
 using Vec = mo_cost_t< D >;                         ///< a vector of costs
 using Frontier = F< Cost , K >;                     ///< the frontier
 using SolutionT = pareto_solution_t< D >;           ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the NAMOAStarDR);
 /// throws std::invalid_argument if the order is not a permutation

 explicit NAMOAStarDR( const D & domain ,
                       const MOSearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_order( objective_order< K >( options.order ) ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Pareto front

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  SolutionT sol;
  MOSearchStats & st = sol.stats;
  NodePool< Rec > recs;
  NodeTable< Rec > table( f_options.expected_states
                          ? NodeTable< Rec >::slots_for(
                                            f_options.expected_states )
                          : 1024 );
  NodePool< Label > labels;
  MOHeap< Item , Better > open;
  Frontier sols;
  std::vector< Label * > goals;
  std::uint64_t seq = 0;

  const auto record = [ & ]( const State & s , const Vec & h ) {
   return( table.find_or_insert( s , f_domain.hash( s ) , [ & ] {
    Rec * r = recs.allocate();
    r->state = s;
    r->h = permute( h , f_order );
    r->goal = f_domain.is_goal( s );
    return( r );
    } ).first );
   };

  const State s0 = f_domain.initial_state();
  Rec * r0 = record( s0 , f_domain.heuristic( s0 ) );
  if( ! dead_end( r0->h ) ) {
   Label * x0 = labels.allocate();
   *x0 = Label{ r0 , Vec{} , nullptr , false };
   r0->open.push_back( x0 );
   open.push( Item{ r0->h , x0 , seq++ } );
   ++st.inserted;
   }

  bool limit = false;
  while( ! open.empty() ) {
   const Item it = open.pop();
   Label * x = it.label;
   if( x->dead ) {
    ++st.duplicates;
    continue;
    }
   ++st.extracted;
   Rec * r = x->rec;
   auto & op = r->open;
   const auto at = std::find( op.begin() , op.end() , x );
   *at = op.back();
   op.pop_back();
   if( sols.dominated( it.f , st.checks ) )
    continue;
   if( r->goal ) {
    sols.insert( it.f , st.checks );
    goals.push_back( x );
    ++st.solutions;
    continue;
    }
   r->front.insert( x->g , st.checks );
   ++st.expanded;

   f_domain.for_each_successor( r->state , unpermute( r->h , f_order ) ,
                                [ & ]( const State & t , const Vec & c ,
                                       const Vec & ht ) {
    ++st.generated;
    Rec * q = record( t , ht );
    if( dead_end( q->h ) )
     return;
    const Vec g = add_costs( x->g , permute( c , f_order ) );
    if( ! q->goal && q->front.dominated( g , st.checks ) )
     return;
    for( const Label * z : q->open ) {
     ++st.checks;
     if( weakly_dominates( z->g , g ) )
      return;
     }
    const Vec f = add_costs( g , q->h );
    if( sols.dominated( f , st.checks ) )
     return;
    auto & qo = q->open;
    for( std::size_t i = 0 ; i < qo.size() ; ) {
     ++st.checks;
     if( weakly_dominates( g , qo[ i ]->g ) ) {
      qo[ i ]->dead = true;
      qo[ i ] = qo.back();
      qo.pop_back();
      }
     else
      ++i;
     }
    Label * y = labels.allocate();
    *y = Label{ q , g , x , false };
    qo.push_back( y );
    open.push( Item{ f , y , seq++ } );
    ++st.inserted;
    } );

   if( f_options.node_limit && ( labels.size() >= f_options.node_limit ) ) {
    limit = true;
    break;
    }
   }

  finish_front( sol , goals , limit , f_order , f_options.paths );
  st.stored = labels.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 struct Rec;

 /// a path: the record of its last state, its g, its predecessor, and
 /// whether a label of its state that dominates it has removed it
 struct Label {
  Rec * rec;
  Vec g;
  Label * parent;
  bool dead;
  };

 /// the record of a state: the state, its h in the order of the search,
 /// whether it is a goal, its truncated G_cl and its G_op
 struct Rec {
  State state;
  Vec h;
  bool goal;
  Frontier front;
  std::vector< Label * > open;
  };

 /// an entry of OPEN: f, the label and the order of the pushes
 struct Item {
  Vec f;
  Label * label;
  std::uint64_t seq;
  };

 /// the order of OPEN: the lexicographic order of f, then the first pushed
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.seq < b.seq ) ) );
   }
  };

 const D & f_domain;                                 ///< the Domain
 MOSearchOptions f_options;                          ///< the options
 std::array< unsigned , K > f_order;                 ///< the order

 };  // end( class( NAMOAStarDR ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a MultiObjectiveDomain with the algorithm a; throws
/// std::invalid_argument for BOA* on more than two objectives

template< MultiObjectiveDomain D >
[[nodiscard]] pareto_solution_t< D > mo_search(
 const D & domain , MOAlgorithm a , const MOSearchOptions & options = {} )
{
 constexpr std::size_t K = objectives_of< D >;
 switch( a ) {
  case( MOAlgorithm::boa ):
   if constexpr( K == 2 )
    return( MOAStar< D , ScalarFrontier >( domain , options ).solve() );
   else
    throw( std::invalid_argument( "mo_search: BOA* is for two "
                                  "objectives" ) );
  case( MOAlgorithm::emoa ):
   return( MOAStar< D , EMOAFrontier >( domain , options ).solve() );
  case( MOAlgorithm::ltmoa ):
   return( MOAStar< D , ArrayFrontier >( domain , options ).solve() );
  default:
   return( NAMOAStarDR< D >( domain , options ).solve() );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the MOAlgorithm of a name (boa, emoa, ltmoa, namoa-dr); throws
/// std::invalid_argument for another one

[[nodiscard]] inline MOAlgorithm mo_algorithm( const std::string & name )
{
 if( name == "boa" )
  return( MOAlgorithm::boa );
 if( name == "emoa" )
  return( MOAlgorithm::emoa );
 if( name == "ltmoa" )
  return( MOAlgorithm::ltmoa );
 if( name == "namoa-dr" )
  return( MOAlgorithm::namoa_dr );
 throw( std::invalid_argument( "mo_algorithm: unknown " + name ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MOAStar */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File MOAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
