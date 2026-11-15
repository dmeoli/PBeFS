/*--------------------------------------------------------------------------*/
/*-------------------------- File BreadthFirst.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template BreadthFirst, the breadth-first
 * search of a Domain with duplicate detection, the baseline that ignores
 * both the costs and the heuristic.
 *
 * With the notation of the page "Best-first search and A*" of the
 * documentation, and the goal test at the generation:
 *
 * \code{.unparsed}
 *  1  procedure BFS(s₀)
 *  2     if s₀ is a goal then return the path (s₀)
 *  3     CLOSED ← {s₀};  L ← (s₀)
 *  4     while L ≠ ∅ do
 *  5        L' ← ()
 *  6        for each u of L, in order, do
 *  7           for each successor v ∉ CLOSED of u do
 *  8              add v to CLOSED, with parent(v) ← u
 *  9              if v is a goal then return the path to v
 * 10              append v to L'
 * 11        L ← L'
 * 12     return unsolvable
 * \endcode
 *
 * The path has the least number of edges, and hence the least cost when
 * the edges cost the same; otherwise its cost is that of its edges, which
 * need not be the least. The search stops at the first goal generated (line
 * 9), but the Domain gives all the successors of its parent, and those
 * that follow the goal are counted in generated without being looked up or
 * stored.
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

#ifndef __PBeFS_BreadthFirst
 #define __PBeFS_BreadthFirst
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BreadthFirst ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the breadth-first search of a Domain
/** The layers L and L' of the pseudocode are two vectors of nodes, and
 * CLOSED is a NodeTable; the heuristic is computed, by the Domain, only
 * because the successors come with it. */

template< Domain D >
class BreadthFirst {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the BreadthFirst)

 explicit BreadthFirst( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  NodePool< NodeT > pool;
  NodeTable< NodeT > table( f_options.expected_nodes
                            ? NodeTable< NodeT >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  std::vector< NodeT * > layer , next;
  NodeT * goal = nullptr;

  const State s0 = f_domain.initial_state();
  NodeT * root = pool.allocate();
  *root = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) ,
                 nullptr , false , false };
  table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   return( root );
   } );
  if( f_domain.is_goal( s0 ) )
   goal = root;
  else
   layer.push_back( root );

  bool limit = false;
  while( ! layer.empty() && ! goal && ! limit ) {
   next.clear();
   for( NodeT * n : layer ) {
    n->closed = n->expanded = true;
    ++st.expanded;
    for_each_child( f_domain , n->state , n->h , n->op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
     ++st.generated;
     if( goal )                        // a successor after the goal: only
      return;                          // counted
     auto [ m , fresh ] = table.find_or_insert( t , f_domain.hash( t ) ,
                                                [ & ] {
      NodeT * x = pool.allocate();
      *x = NodeT{ t , n->g + c , ht , n , false , false , op };
      return( x );
      } );
     if( ! fresh ) {
      ++st.duplicates;
      return;
      }
     if( f_domain.is_goal( t ) )
      goal = m;
     else
      next.push_back( m );
     } );
    if( goal )
     break;
    if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
     limit = true;
     break;
     }
    }
   std::swap( layer , next );
   }

  if( goal ) {
   sol.status = Status::solved;
   sol.cost = goal->g;
   sol.path = trace_path( goal );
   }
  else
   if( limit )
    sol.status = Status::node_limit;
  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options

 };  // end( class( BreadthFirst ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain breadth-first

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > breadth_first(
 const D & domain , const SearchOptions & options = {} )
{
 return( BreadthFirst< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BreadthFirst */

/*--------------------------------------------------------------------------*/
/*------------------------ End File BreadthFirst.h -------------------------*/
/*--------------------------------------------------------------------------*/
