/*--------------------------------------------------------------------------*/
/*------------------------------ File AStar.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template AStar, the sequential A* search of
 * Hart, Nilsson and Raphael (1968).
 *
 * The search follows Algorithm 1 of Fukunaga et al. (2017), with the
 * notation and the conventions of the page "Best-first search and A*" of
 * the documentation: CLOSED is the table of the generated states, each one
 * with its node and a flag that tells whether the node is closed, i.e.,
 * expanded with its current \f$ g \f$, and OPEN is a priority queue of
 * entries \f$ ( f , g , u ) \f$ ordered by increasing \f$ f \f$ and, for
 * equal \f$ f \f$, by decreasing \f$ g \f$:
 *
 * \code{.unparsed}
 *  1  procedure A*(s₀)
 *  2     CLOSED ← {s₀}, with g(s₀) ← 0 and parent(s₀) ← nil
 *  3     OPEN ← {(f(s₀), g(s₀), s₀)}
 *  4     while OPEN ≠ ∅ do
 *  5        (f, g, u) ← the first entry of OPEN; remove it from OPEN
 *  6        if u is closed or g ≠ g(u) then continue      ▷ a stale entry
 *  7        close u
 *  8        if u is a goal then return the path to u, of cost g(u)
 *  9        for each successor v of u do
 * 10           g' ← g(u) + c(u, v)
 * 11           if v ∉ CLOSED then
 * 12              add v to CLOSED, with g(v) ← g' and parent(v) ← u
 * 13              push (g' + h(v), g', v) in OPEN
 * 14           else if g' < g(v) then           ▷ a cheaper path: reopen v
 * 15              g(v) ← g';  parent(v) ← u;  v is not closed
 * 16              push (g' + h(v), g', v) in OPEN
 * 17     return unsolvable
 * \endcode
 *
 * With SearchOptions::goal_at_generation (off by default) the search also
 * keeps the incumbent \f$ \hat{u} \f$, the goal of least \f$ g \f$ among
 * those generated with \f$ h = 0 \f$ (every goal, if \f$ h \f$ is
 * admissible), and line 5 is preceded by
 *
 * \code{.unparsed}
 *  4' if û exists and no entry of OPEN has f < g(û) then
 *        return the path to û, of cost g(û)
 * \endcode
 *
 * where the stale entries at the top of OPEN are dropped first. If
 * \f$ h \f$ is admissible the returned path is optimal: as long as the
 * search runs, OPEN holds an entry \f$ ( f , g , v ) \f$ of a node of an
 * optimal path with \f$ g = g^*( v ) \f$, whence \f$ f \leq C^* \f$
 * [see Hart, Nilsson and Raphael, 1968, lemma 1 and its corollary, which
 * hold with the reopening of lines 14 to 16], and therefore
 * \f$ g( \hat{u} ) \leq f_{\min} \leq C^* \f$ at line 4', while
 * \f$ g( \hat{u} ) \geq C^* \f$ since it is the cost of a path. If the
 * parent of a goal on an optimal path has \f$ f < C^* \f$ (e.g., with the
 * blind heuristic and positive costs), that goal is the incumbent before
 * any node of \f$ f = C^* \f$ is selected, and the search expands only
 * nodes of \f$ f < C^* \f$, which every A* expands; the goal test of
 * line 8 alone also expands the nodes of \f$ f = C^* \f$ that come before
 * a goal in the order of the ties. A goal selected at line 8 still ends the
 * search. The goals that improve \f$ \hat{u} \f$ are counted in
 * SearchStats::goals_generated, and a successor is tested only if its
 * \f$ h \f$ is 0, so that a domain whose heuristic is informative pays
 * nothing for it. With a weight \f$ w > 1 \f$ line 4' compares the keys
 * \f$ f' = g + w h \f$ with \f$ g( \hat{u} ) \f$ on their scale [see
 * Weight::live], and the path keeps the bound \f$ w C^* \f$: the node of
 * the optimal path in OPEN has \f$ f' \leq w C^* \f$, whence
 * \f$ g( \hat{u} ) \leq f'_{\min} \leq w C^* \f$.
 *
 * Neither OPEN list decreases a key: line 16 pushes a second entry of
 * \f$ v \f$, and the first one, which records the old \f$ g \f$, is
 * dropped at line 6 when it reaches the top. The goal test at line 8, when
 * a node is selected and not when it is generated, makes the first goal
 * selected an optimal one if \f$ h \f$ is admissible, and the reopening
 * at lines 14 to 16 keeps the search optimal when \f$ h \f$ is admissible
 * but not consistent; with a consistent \f$ h \f$ no closed node is ever
 * reopened. The successors of \f$ u \f$ are generated first, with the slots
 * of their lookups in CLOSED prefetched, and looked up afterwards, so that
 * the cache misses of the lookups overlap.
 *
 * With a weight \f$ w > 1 \f$ [see SearchOptions] this is weighted A*
 * [see Pohl, AIJ 1970]: OPEN is ordered by \f$ f' = g + w h \f$ at lines 13
 * and 16, the first goal selected costs at most \f$ w C^* \f$, and with
 * duplicate dropping the test of line 14 is that of Weight::improves(),
 * which reopens a closed node only when the saving exceeds what the weight
 * allows [see Weight.h].
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

#ifndef __PBeFS_AStar
 #define __PBeFS_AStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <stdexcept>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/OpenList.h"

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
/*------------------------------ CLASS AStar -------------------------------*/
/*--------------------------------------------------------------------------*/
/// sequential A* over a Domain, with the OPEN list Open
/** This is Algorithm 1 of Fukunaga et al. (2017): the goal test is done
 * when a node is selected for expansion, so that the first goal expanded is
 * optimal if h is admissible. A node reached through a cheaper path is
 * reopened even if it is closed, so that the search stays admissible with
 * an admissible but inconsistent heuristic; with a consistent one no node
 * is ever reopened. The weight and the duplicate dropping of the
 * SearchOptions make it weighted A*, whose path costs at most
 * \f$ w C^* \f$; the anytime search is AnytimeAStar, and the constructor
 * throws std::invalid_argument if the options ask for it.
 *
 * Open defaults to BucketOpenList for integer costs and to HeapOpenList
 * otherwise [see OpenList.h]. */

template< Domain D ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class AStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the AStar)

 explicit AStar( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) , f_w( options ) {
  if( f_w.anytime() )
   throw( std::invalid_argument( "AStar: the anytime search is "
                                 "AnytimeAStar" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  // the weight 1 known at compile time, so that A* pays nothing for it
  return( f_w.unit() ? run( UnitWeight< Cost >() ) : run( f_w ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the search with the rules of w [see Weight and UnitWeight]

 template< class W >
 [[nodiscard]] SolutionT run( const W & w ) const
 {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  NodePool< NodeT > pool;
  NodeTable< NodeT > table( f_options.expected_nodes
                            ? NodeTable< NodeT >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  Open open;
  std::vector< Kid > kids;                 // the successors of a node

  const State s0 = f_domain.initial_state();
  NodeT * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * n = pool.allocate();
   *n = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               false , false };
   return( n );
   } ).first;
  open.push( root , root->g , w.key( root->g , root->h ) );
  NodeT * best = nullptr;                  // the incumbent goal, if any

  for( ; ; ) {
   if( best ) {                            // line 4'
    while( ! open.empty() && open.top().stale() )
     ( void ) open.pop();
    if( open.empty() || ! w.live( open.top_f() , best->g ) ) {
     sol.status = Status::solved;
     sol.bound = f_w.lower_bound( best->g );
     sol.path = trace_path( best );
     sol.cost = std::min( best->g ,
                          Cost( path_cost( f_domain , sol.path ) ) );
     break;
     }
    }
   if( open.empty() )
    break;
   const OpenEntry< NodeT > e = open.pop();
   if( e.stale() )
    continue;

   NodeT * n = e.node;
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;

   if( f_domain.is_goal( n->state ) ) {
    sol.status = Status::solved;
    sol.bound = f_w.lower_bound( n->g );
    sol.path = trace_path( n );
    sol.cost = std::min( n->g , Cost( path_cost( f_domain , sol.path ) ) );
    break;
    }

   // first the successors, with the slots of their lookups prefetched,
   // then the lookups, whose cache misses thus overlap
   kids.clear();
   const auto kid = [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
    const std::uint64_t hash = f_domain.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    };
   for_each_child( f_domain , n->state , n->h , n->op , kid );
   st.generated += kids.size();
   for( const Kid & k : kids ) {
    const Cost g = n->g + k.c;
    auto [ m , fresh ] = table.find_or_insert( k.s , k.hash , [ & ] {
     NodeT * x = pool.allocate();
     *x = NodeT{ k.s , g , k.h , n , false , false , k.op };
     return( x );
     } );
    if( fresh )
     open.push( m , g , w.key( g , k.h ) );
    else
     if( w.improves( m->g , g , m->closed ,
                     [ & ] { return( n->h ); } , k.h ) ) {
      m->g = g;           // a cheaper path: reopen
      m->parent = n;
      m->op = k.op;
      m->closed = false;
      open.push( m , g , w.key( g , k.h ) );
      }
     else {
      ++st.duplicates;
      continue;
      }
    if( f_options.goal_at_generation && ( k.h == Cost( 0 ) ) &&
        ( ( ! best ) || ( m->g < best->g ) ) && f_domain.is_goal( k.s ) ) {
     best = m;
     ++st.goals_generated;
     }
    }

   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }

  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

 /// a successor: its state, the cost of its edge, its h, its hash
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint64_t hash;
  std::uint8_t op;
  };

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 Weight< Cost > f_w;                                 ///< the weight

 };  // end( class( AStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with A*, choosing the OPEN list at run time
/** Searches \p domain with A*; the OPEN list is the bucket queue for
 * integer costs and the heap otherwise, unless \p open asks for the heap,
 * or for the heap whose ties go to the smaller g [see LowGOpenList]. */

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > astar(
 const D & domain , OpenKind open = OpenKind::automatic ,
 const SearchOptions & options = {} )
{
 using NodeT = Node< state_t< D > , cost_t< D > >;
 if( open == OpenKind::heap_low_g )
  return( AStar< D , LowGOpenList< HeapOpenList< NodeT > > >( domain ,
                                                              options )
          .solve() );
 if constexpr( std::is_integral_v< cost_t< D > > )
  if( open != OpenKind::heap )
   return( AStar< D , BucketOpenList< NodeT > >( domain , options ).solve()
           );
 return( AStar< D , HeapOpenList< NodeT > >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AStar */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File AStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
