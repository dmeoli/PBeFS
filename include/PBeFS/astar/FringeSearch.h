/*--------------------------------------------------------------------------*/
/*-------------------------- File FringeSearch.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template FringeSearch, the search of
 * Björnsson, Enzenberger, Holte and Schaeffer (CIG 2005), between IDA*
 * and A*: the frontier of IDA* is kept in a list, so that an iteration
 * starts from the nodes where the one before stopped instead of from the
 * root, and the g of every state seen is kept in a cache, so that a state
 * is not expanded again through a path that is not cheaper.
 *
 * \code{.unparsed}
 *  1  procedure FRINGE(s₀)
 *  2     L ← ⟨s₀⟩;  C[s₀] ← (0, nil);  θ ← h(s₀)
 *  3     while L ≠ ⟨⟩ do
 *  4        θ' ← ∞
 *  5        for each node n of L, from the first to the last do
 *  6           if g(n) + h(n) > θ then θ' ← min(θ', g(n) + h(n));  continue
 *  7           if n is a goal then return the path to n
 *  8           for each successor c of n do
 *  9              g' ← g(n) + c(n, c)
 * 10              if c ∈ C and g' ≥ g(c) then continue
 * 11              remove c from L if it is there;  C[c] ← (g', n)
 * 12              insert c in L just after n
 * 13           remove n from L
 * 14        θ ← θ'
 * \endcode
 *
 * With an admissible h the first goal found within the threshold θ is
 * optimal: as in IDA*, θ is at each iteration a lower bound of
 * \f$ C^* \f$, and every state of a path within θ that the iteration
 * meets is expanded in it (the children are visited right after their
 * parent, which is the depth-first order of IDA*). The memory is that of
 * A*, the time per node smaller, since the list costs \f$ O( 1 ) \f$ per
 * operation instead of the \f$ O( \log n ) \f$ of a heap; the counter
 * expanded counts the expansions (line 8), reexpanded those of a state
 * already expanded, and the number of iterations is iterations().
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

#ifndef __PBeFS_FringeSearch
 #define __PBeFS_FringeSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <stdexcept>

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
/*--------------------------- CLASS FringeSearch ---------------------------*/
/*--------------------------------------------------------------------------*/
/// Fringe Search, the frontier of IDA* kept between the iterations

template< Domain D >
class FringeSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit FringeSearch( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  if( options.weight != 1 )
   throw( std::invalid_argument( "FringeSearch: no weight" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void )
 {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  NodePool< FNode > pool;
  NodeTable< FNode > table( f_options.expected_nodes
                            ? NodeTable< FNode >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  const State s0 = f_domain.initial_state();
  FNode * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   FNode * x = pool.allocate();
   *x = FNode{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               no_op , false , false , nullptr , nullptr };
   return( x );
   } ).first;
  // the list L, doubly linked through the nodes, with a sentinel
  FNode head{};
  head.next = head.prev = &head;
  const auto link_after = [ & ]( FNode * a , FNode * x ) {
   x->prev = a;
   x->next = a->next;
   a->next->prev = x;
   a->next = x;
   x->listed = true;
   };
  const auto unlink = [ & ]( FNode * x ) {
   x->prev->next = x->next;
   x->next->prev = x->prev;
   x->listed = false;
   };
  link_after( &head , root );
  std::vector< FNode * > kids;
  Cost theta = root->h;
  f_iterations = 0;
  const FNode * goal = nullptr;
  bool limit = false;
  while( ( head.next != &head ) && ! goal && ! limit ) {
   ++f_iterations;
   Cost next = infinite_cost< Cost >();
   for( FNode * n = head.next ; n != &head ; ) {
    const Cost f = n->g + n->h;
    if( f > theta ) {
     if( f < next )
      next = f;
     n = n->next;
     continue;
     }
    if( f_domain.is_goal( n->state ) ) {
     goal = n;
     break;
     }
    ++st.expanded;
    if( n->expanded )
     ++st.reexpanded;
    n->expanded = true;
    kids.clear();
    for_each_child( f_domain , n->state , n->h , n->op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
     ++st.generated;
     const Cost g = n->g + c;
     bool fresh = false;
     FNode * m = table.find_or_insert( t , f_domain.hash( t ) , [ & ] {
      fresh = true;
      FNode * x = pool.allocate();
      *x = FNode{ t , g , ht , n , op , false , false , nullptr , nullptr };
      return( x );
      } ).first;
     if( ! fresh ) {
      if( g >= m->g ) {
       ++st.duplicates;
       return;
       }
      m->g = g;
      m->parent = n;
      m->op = op;
      if( m->listed )
       unlink( m );
      }
     kids.push_back( m );
     } );
    // the children right after n, in the order of their generation
    for( auto it = kids.rbegin() ; it != kids.rend() ; ++it )
     link_after( n , *it );
    FNode * after = n->next;
    unlink( n );
    n = after;
    if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
     limit = true;
     break;
     }
    }
   theta = next;
   }
  if( goal ) {
   sol.status = Status::solved;
   sol.cost = goal->g;
   sol.bound = sol.cost;
   sol.path = trace_path( goal );
   }
  else
   if( limit )
    sol.status = Status::node_limit;
  st.stored = pool.size();
  st.evaluated = st.stored;
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the iterations of the last search

 [[nodiscard]] std::size_t iterations( void ) const noexcept {
  return( f_iterations );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a node of the cache C, linked in the list L if listed
 struct FNode {
  State state;
  Cost g;
  Cost h;
  FNode * parent;
  std::uint8_t op;
  bool listed;
  bool expanded;
  FNode * prev;
  FNode * next;
  };

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 std::size_t f_iterations = 0;                       ///< iterations

 };  // end( class( FringeSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with FringeSearch

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > fringe_search(
 const D & domain , const SearchOptions & options = {} )
{
 return( FringeSearch< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_FringeSearch */

/*--------------------------------------------------------------------------*/
/*------------------------ End File FringeSearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
