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
 * The parent is skipped by its operator in an OperatorDomain [see
 * Domain.h], and otherwise by comparing the successor with the state two
 * places back on the path; no other duplicate is detected, so that the
 * search fits the domains whose graph has few short cycles, such as the
 * sliding-tile puzzles, where IDA* with a pattern database is the method
 * of Korf and Felner (2002) for the 24-puzzle. The costs must be positive
 * on every cycle, or an iteration may not end; and when no goal is
 * reachable the search proves it only on an acyclic graph, since on one
 * with cycles the bound grows without end, unless the node limit stops it.
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

#ifndef __PBeFS_IDAStar
 #define __PBeFS_IDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS IDAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the iterative-deepening A* search of a Domain
/** The counters of the Solution are summed over the iterations: expanded
 * counts every expansion, those of the states already expanded in an
 * earlier iteration included, and stored is the length of the longest path
 * held, i.e., the memory of the search in states. The node_limit of the
 * SearchOptions bounds the expansions, and the search stops with
 * Status::node_limit when it reaches it; expected_nodes is not used. */

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

 explicit IDAStar( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  Run r( f_domain , f_options.node_limit );
  const State s0 = f_domain.initial_state();
  const Cost h0 = f_domain.heuristic( s0 );
  r.path.push_back( s0 );
  Cost bound = h0;
  SolutionT sol;
  for( ; ; ) {
   r.next = infinite_cost< Cost >();
   if( r.dfs( s0 , Cost( 0 ) , h0 , no_op , bound ) ) {
    sol.status = Status::solved;
    sol.cost = r.cost;
    sol.path = r.path;
    break;
    }
   if( r.stop ) {
    sol.status = Status::node_limit;
    break;
    }
   if( r.next == infinite_cost< Cost >() ) {
    sol.status = Status::unsolvable;
    break;
    }
   bound = r.next;
   }
  r.st.expanded_by = { r.st.expanded };
  r.st.stored_by = { r.st.stored };
  r.st.seconds = clock.seconds();
  sol.stats = std::move( r.st );
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the state of a search: the path, the next bound and the counters
 struct Run {

  const D & d;                                    ///< the Domain
  std::size_t limit;                              ///< the expansion limit
  std::vector< State > path;                      ///< the path from s0
  Cost next{};                                    ///< the next bound
  Cost cost{};                                    ///< the cost of a goal
  bool stop = false;                              ///< the limit is reached
  SearchStats st;                                 ///< the counters

  Run( const D & domain , std::size_t node_limit )
   : d( domain ) , limit( node_limit ) {}

  // the depth-first search below u, of cost g and heuristic h with
  // g + h within the bound, generated by the operator op; true if a goal is
  // found, the path then ending at it
  bool dfs( const State & u , Cost g , Cost h , std::uint8_t op ,
            Cost bound ) {
   if( d.is_goal( u ) ) {
    cost = g;
    return( true );
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
    if( g + c + hv > bound ) {                   // cut: a candidate bound
     next = std::min( next , g + c + hv );
     return;
     }
    path.push_back( v );
    st.stored = std::max< std::uint64_t >( st.stored , path.size() );
    if( dfs( v , g + c , hv , opv , bound ) )
     found = true;
    else
     path.pop_back();
    };
   if constexpr( OperatorDomain< D > )
    d.for_each_successor( u , h , op , [ & ]( const State & v , Cost c ,
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
 SearchOptions f_options;                            ///< the options

 };  // end( class( IDAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with IDA*

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > idastar(
 const D & domain , const SearchOptions & options = {} )
{
 return( IDAStar< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_IDAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File IDAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
