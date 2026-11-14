/*--------------------------------------------------------------------------*/
/*----------------------------- File IDAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template IDAStar, the iterative-deepening A*
 * of Korf (1985), which finds an optimal path with a memory linear in its
 * length: a depth-first search bounded by \f$ f = g + h \f$, repeated with
 * the bound raised each time to the least \f$ f \f$ that exceeded it.
 *
 * With the notation of the page "Best-first search and A*" of the
 * documentation, and \f$ h \f$ admissible:
 *
 * \code{.unparsed}
 *  1  procedure IDA*(s₀)
 *  2     β ← h(s₀)
 *  3     loop
 *  4        β' ← ∞
 *  5        if DFS(s₀, 0, β) then return the path
 *  6        if β' = ∞ then return unsolvable
 *  7        β ← β'
 *
 *  1  function DFS(u, g, β)
 *  2     if g + h(u) > β then β' ← min(β', g + h(u)); return false
 *  3     if u is a goal then return true
 *  4     for each successor v of u but the parent of u do
 *  5        push v on the path
 *  6        if DFS(v, g + c(u, v), β) then return true
 *  7        pop v from the path
 *  8     return false
 * \endcode
 *
 * The parent is skipped by its operator in an OperatorDomain, with the
 * redundant pairs of moves of a TreeOperatorDomain [see Domain.h], and
 * otherwise by comparing the successor with the state two places back on
 * the path; no other duplicate is detected, so that the search fits the
 * domains whose graph has few short cycles, such as the sliding-tile
 * puzzles, where IDA* with a pattern database is the method of Korf and
 * Felner (2002) for the 24-puzzle. The costs must be positive
 * on every cycle, or an iteration may not end; and when no goal is
 * reachable the search proves it only on an acyclic graph, since on one
 * with cycles the bound grows without end, unless the node limit stops it.
 * With a bound given, the search is the cost-bounded depth-first search,
 * i.e., the first call of line 5 with that \f$ \beta \f$ and no other: its
 * path, if any, has no node of \f$ f > \beta \f$, but need not be optimal.
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

#ifndef __PBeFS_IDAStar
 #define __PBeFS_IDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Deepening.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/astar/Weight.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS IDAOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of IDAStar: those of the searches, and the bounds of
/// IDA*-CR [see Deepening.h]

struct IDAOptions : public SearchOptions , public DeepeningOptions {

 IDAOptions( void ) = default;

 /// the options of the searches, with the classical bounds
 IDAOptions( const SearchOptions & o ) : SearchOptions( o ) {}

 };  // end( struct( IDAOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS IDAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the iterative-deepening A* search of a Domain
/** The counters of the Solution are summed over the iterations: expanded
 * counts every expansion, those of the states already expanded in an
 * earlier iteration included, and stored is the length of the longest path
 * held, i.e., the memory of the search in states. The node_limit of the
 * SearchOptions bounds the expansions, and the search stops with
 * Status::node_limit when it reaches it; expected_nodes is not used. With
 * a weight \f$ w > 1 \f$ the thresholds are on \f$ f' = g + w h \f$ [see
 * Weight], and the path found costs at most \f$ w C^* \f$, since every node
 * of an optimal path has \f$ f' \leq w C^* \f$; drop_duplicates does
 * nothing, there being no duplicate detection, and the constructor throws
 * std::invalid_argument if the options ask for an anytime search. With
 * IDAOptions::growth above 1 the bounds are those of IDA*-CR, and the
 * iteration that finds the first goal with a bound above the least f
 * pruned before is completed as a branch-and-bound if IDAOptions::optimal
 * [see Deepening.h]; if it is not, the path need not be within the bound
 * \f$ w C^* \f$, and Solution::bound is the one that the least \f$ f' \f$
 * pruned before gives [see Weight::key_lower_bound]. */

template< Domain D >
class IDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the IDAStar)

 explicit IDAStar( const D & domain , const IDAOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_w( options.weight ) {
  if( options.anytime && ! f_w.unit() )
   throw( std::invalid_argument( "IDAStar: no anytime search" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( search( f_w.key( Cost( 0 ) ,
                           Cost( f_domain.heuristic(
                            f_domain.initial_state() ) ) ) , true ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost-bounded depth-first search: the first path found whose
 /// nodes have f (f' with a weight) not above bound, which need not be one
 /// of least cost, i.e., one iteration of IDA*; Status::unsolvable if there
 /// is none

 [[nodiscard]] SolutionT solve( Cost bound ) const {
  return( search( f_w.scale( bound ) , false ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the iterations of IDA* from bound, or only the first one
 SolutionT search( Cost bound , bool iterate ) const
 {
  Stopwatch clock;
  Run r( f_domain , f_options.node_limit , f_w );
  const State s0 = f_domain.initial_state();
  const Cost h0 = f_domain.heuristic( s0 );
  r.path.push_back( s0 );
  SolutionT sol;
  const bool cr = iterate && ( f_options.growth > 1 );
  Cost low = f_w.key( Cost( 0 ) , h0 );   // the least f' pruned before
  for( ; ; ) {
   r.cut.reset( bound , cr );
   r.bnb = cr && f_options.optimal && ( bound > low );
   r.best = infinite_cost< Cost >();
   const std::uint64_t before = r.st.expanded;
   if( ( f_w.key( Cost( 0 ) , h0 ) <= bound ) &&
       r.dfs( s0 , Cost( 0 ) , h0 , no_op , bound ) ) {
    sol.status = Status::solved;
    sol.cost = r.cost;
    sol.bound = bound > low ? std::min( r.cost , f_w.key_lower_bound( low ) )
                            : f_w.lower_bound( r.cost );
    sol.path = r.path;
    break;
    }
   if( r.best < infinite_cost< Cost >() ) {    // the branch-and-bound
    sol.status = Status::solved;
    sol.cost = r.best;
    sol.bound = f_w.lower_bound( r.best );
    sol.path = std::move( r.best_path );
    break;
    }
   if( r.stop ) {
    sol.status = Status::node_limit;
    break;
    }
   low = r.cut.low();
   if( ! iterate || ( low == infinite_cost< Cost >() ) ) {
    sol.status = Status::unsolvable;
    break;
    }
   bound = r.cut.next( f_options.growth , r.st.expanded - before );
   }
  r.st.expanded_by = { r.st.expanded };
  r.st.stored_by = { r.st.stored };
  r.st.seconds = clock.seconds();
  sol.stats = std::move( r.st );
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the state of a search: the path, the next bound and the counters
 struct Run {

  const D & d;                                    ///< the Domain
  std::size_t limit;                              ///< the expansion limit
  std::vector< State > path;                      ///< the path from s0
  CutoffHistogram< Cost > cut;                    ///< the f pruned
  Cost cost{};                                    ///< the cost of a goal
  bool bnb = false;                               ///< goes on at a goal
  Cost best{};                                    ///< the incumbent if bnb
  std::vector< State > best_path;                 ///< its path
  bool stop = false;                              ///< the limit is reached
  SearchStats st;                                 ///< the counters
  const Weight< Cost > & w;                       ///< the weight

  Run( const D & domain , std::size_t node_limit , const Weight< Cost > & x )
   : d( domain ) , limit( node_limit ) , w( x ) {}

  // the depth-first search below u, of cost g and heuristic h with
  // g + h within the bound, generated by the operator op; true if a goal is
  // found, the path then ending at it
  bool dfs( const State & u , Cost g , Cost h , std::uint8_t op ,
            Cost bound ) {
   if( d.is_goal( u ) ) {
    if( ! bnb ) {
     cost = g;
     return( true );
     }
    if( g < best ) {                    // a better incumbent, and on
     best = g;
     best_path = path;
     }
    return( false );
    }
   if( limit && ( st.expanded >= limit ) ) {
    stop = true;
    return( false );
    }
   ++st.expanded;
   bool found = false;
   const auto child = [ & ]( const State & v , Cost c , Cost hv ,
                             std::uint8_t opv ) {
    ++st.generated;
    const Cost k = w.key( g + c , hv );
    if( k > bound ) {                            // cut: a candidate bound
     cut.add( k );
     return;
     }
    if( bnb && w.prunes( g + c , hv , best ) )   // not below the incumbent
     return;
    path.push_back( v );
    st.stored = std::max< std::uint64_t >( st.stored , path.size() );
    if( dfs( v , g + c , hv , opv , bound ) )
     found = true;
    else
     path.pop_back();
    };
   if constexpr( OperatorDomain< D > )
    for_each_tree_child( d , u , h , op , [ & ]( const State & v , Cost c ,
                                                 Cost hv , std::uint8_t opv ) {
     if( ! found && ! stop )
      child( v , c , hv , opv );
     } );
   else
    d.for_each_successor( u , h , [ & ]( const State & v , Cost c ,
                                         Cost hv ) {
     if( found || stop )
      return;
     if( ( path.size() >= 2 ) && ( v == path[ path.size() - 2 ] ) ) {
      ++st.generated;                            // back to the parent
      return;
      }
     child( v , c , hv , no_op );
     } );
   return( found );
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 IDAOptions f_options;                               ///< the options
 Weight< Cost > f_w;                                 ///< the weight

 };  // end( class( IDAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with IDA*

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > idastar(
 const D & domain , const IDAOptions & options = {} )
{
 return( IDAStar< D >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain by the depth-first search bounded by f not above
/// bound [see IDAStar::solve( Cost )]

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > cost_bounded_dfs(
 const D & domain , cost_t< D > bound , const SearchOptions & options = {} )
{
 return( IDAStar< D >( domain , options ).solve( bound ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_IDAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File IDAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
