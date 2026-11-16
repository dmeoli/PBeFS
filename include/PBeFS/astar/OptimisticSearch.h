/*--------------------------------------------------------------------------*/
/*------------------------- File OptimisticSearch.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template OptimisticSearch, the bounded
 * suboptimal search of Thayer and Ruml (ICAPS 2008): it searches with a
 * weight larger than the bound asked, and then expands by \f$ f \f$ only as
 * much as is needed to prove that the path found is within the bound.
 *
 * With the bound \f$ w \f$, the search orders OPEN by
 * \f$ f'_o = g + w_o h \f$, \f$ w_o \f$ being the optimistic weight
 * (\f$ 2 w - 1 \f$ by default, as in the paper), and keeps the same nodes
 * in FQ, ordered by \f$ f \f$; \f$ z \f$ is the incumbent, and the notation is
 * that of AStar.h:
 *
 * \code{.unparsed}
 *  1  procedure OPTIMISTIC(s₀, w, w_o)
 *  2     z ← ∞;  push s₀ in OPEN (by f'_o) and in FQ (by f)
 *  3     while FQ ≠ ∅ and w · (the least f of FQ) < z do
 *  4        if z = ∞ or OPEN has a node of f'_o < z then
 *  5           u ← the first node of OPEN, removed
 *  6        else u ← the first node of FQ, removed          ▷ the cleanup
 *  7        if f(u) ≥ z then continue                        ▷ pruned
 *  8        close u
 *  9        if u is a goal then z ← g(u), recording u;  continue
 * 10        for each successor v of u with g(u) + c(u, v) + h(v) < z do
 * 11           insert v as A* does, in OPEN and in FQ
 * 12     return the path of cost z
 * \endcode
 *
 * The stale entries are dropped as they surface in either list, a node
 * removed from one of them being closed, and therefore stale, in the
 * other. Every node of an optimal path with \f$ g = g^* \f$ that is not
 * closed is in FQ (or pruned, when \f$ z \leq C^* \f$), so that the least
 * \f$ f \f$ of FQ is at most \f$ C^* \f$, and the test of line 3 stops only
 * when \f$ z \leq w C^* \f$, or when FQ is empty and \f$ z \f$ is optimal;
 * the least \f$ f \f$ of FQ is the lower bound of the Solution.
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

#ifndef __PBeFS_OptimisticSearch
 #define __PBeFS_OptimisticSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <stdexcept>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/OpenList.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/astar/Weight.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS OptimisticOptions -------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of OptimisticSearch: the weight of the SearchOptions is the
/// bound \f$ w \f$, and optimism the weight \f$ w_o \f$ of the order

//! [bounded-suboptimal-sig-optimistic]
struct OptimisticOptions : public SearchOptions {

 /// the optimistic weight \f$ w_o \geq 1 \f$; 0 stands for \f$ 2 w - 1 \f$
 double optimism = 0;

 };  // end( struct( OptimisticOptions ) )
//! [bounded-suboptimal-sig-optimistic]

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS OptimisticSearch -------------------------*/
/*--------------------------------------------------------------------------*/
/// the optimistic search of Thayer and Ruml
/** The path found costs at most \f$ w C^* \f$, \f$ w \f$ being the weight of
 * the options; each path that lowers the incumbent is recorded in
 * Solution::improvements with the least \f$ f \f$ of FQ at that moment.
 * The closed nodes reached through a cheaper path are reopened, and the
 * goal is tested when it is selected. With the node limit, the Solution
 * keeps the best path found, with Status::node_limit. anytime and
 * drop_duplicates of the options are not read; the constructor throws
 * std::invalid_argument if a weight is below 1. */

template< Domain D ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class OptimisticSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the search)

 explicit OptimisticSearch( const D & domain ,
                            const OptimisticOptions & options = {} )
  : f_domain( domain ) , f_options( options ) , f_wb( options.weight ) ,
    f_wo( options.optimism == 0 ? 2 * options.weight - 1
                                : options.optimism ) {}

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
  Open open , fq;
  std::vector< Kid > kids;
  NodeT * goal = nullptr;
  Cost z = infinite_cost< Cost >();
  bool limit_hit = false;

  const State s0 = f_domain.initial_state();
  NodeT * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * n = pool.allocate();
   *n = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               false , false };
   return( n );
   } ).first;
  open.push( root , Cost( 0 ) , f_wo.key( Cost( 0 ) , root->h ) );
  fq.push( root );

  // drops the stale entries at the top of q
  const auto clean = []( Open & q ) {
   while( ! q.empty() && q.top().stale() )
    q.pop();
   };

  for( ; ; ) {
   clean( fq );
   if( fq.empty() || f_wb.prunes( fq.top_f() , Cost( 0 ) , z ) )
    break;
   clean( open );
   const bool optimistic = ! open.empty() && f_wo.live( open.top_f() , z );
   const OpenEntry< NodeT > e = optimistic ? open.pop() : fq.pop();
   NodeT * n = e.node;
   if( e.g + n->h >= z )
    continue;
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;

   if( f_domain.is_goal( n->state ) ) {
    z = n->g;
    goal = n;
    clean( fq );
    sol.improvements.push_back( { z , fq.empty() ? z
                                                 : std::min( z , fq.top_f() ) ,
                                  st.expanded , clock.seconds() } );
    continue;
    }

   kids.clear();
   for_each_child( f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    const std::uint64_t hash = f_domain.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    } );
   st.generated += kids.size();
   for( const Kid & k : kids ) {
    const Cost g = n->g + k.c;
    if( g + k.h >= z )
     continue;
    auto [ m , fresh ] = table.find_or_insert( k.s , k.hash , [ & ] {
     NodeT * x = pool.allocate();
     *x = NodeT{ k.s , g , k.h , n , false , false , k.op };
     return( x );
     } );
    if( ! fresh ) {
     if( ! cheaper( g , m->g ) ) {
      ++st.duplicates;
      continue;
      }
     m->g = g;                     // a cheaper path: reopen
     m->parent = n;
     m->op = k.op;
     m->closed = false;
     }
    open.push( m , g , f_wo.key( g , k.h ) );
    fq.push( m , g );
    }

   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    limit_hit = true;
    break;
    }
   }

  clean( fq );
  if( limit_hit )
   sol.status = Status::node_limit;
  else
   if( goal )
    sol.status = Status::solved;
  if( goal ) {
   sol.bound = fq.empty() ? z : std::min( z , fq.top_f() );
   sol.path = trace_path( goal );
   sol.cost = std::min( z , Cost( path_cost( f_domain , sol.path ) ) );
   settle_improvements( sol );
   }
  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a successor: its state, the cost of its edge, its h, its hash
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint64_t hash;
  std::uint8_t op;
  };

 const D & f_domain;                                 ///< the Domain
 OptimisticOptions f_options;                        ///< the options
 Weight< Cost > f_wb;                                ///< the bound w
 Weight< Cost > f_wo;                                ///< the weight w_o

 };  // end( class( OptimisticSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with the optimistic search, choosing the OPEN list at
/// run time as astar() does

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > optimistic_search(
 const D & domain , const OptimisticOptions & options = {} ,
 OpenKind open = OpenKind::automatic )
{
 using NodeT = Node< state_t< D > , cost_t< D > >;
 if( open == OpenKind::heap_low_g )
  return( OptimisticSearch< D , LowGOpenList< HeapOpenList< NodeT > > >(
           domain , options ).solve() );
 if constexpr( std::is_integral_v< cost_t< D > > )
  if( open != OpenKind::heap )
   return( OptimisticSearch< D , BucketOpenList< NodeT > >( domain ,
                                                            options )
           .solve() );
 return( OptimisticSearch< D , HeapOpenList< NodeT > >( domain , options )
         .solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_OptimisticSearch */

/*--------------------------------------------------------------------------*/
/*----------------------- End File OptimisticSearch.h ----------------------*/
/*--------------------------------------------------------------------------*/
