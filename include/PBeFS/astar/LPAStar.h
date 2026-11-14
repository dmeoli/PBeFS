/*--------------------------------------------------------------------------*/
/*----------------------------- File LPAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template LPAStar, the Lifelong Planning A* of
 * Koenig, Likhachev and Furcy (AIJ 2004): an incremental search, which
 * after a change of the costs of some edges finds the new shortest path by
 * repairing the g values of the previous search instead of searching again.
 *
 * Every state \f$ s \f$ has the value \f$ g( s ) \f$ and the one-step
 * lookahead \f$ \mathrm{rhs}( s ) = \min_{ s' \in \mathrm{pred}( s ) }
 * g( s' ) + c( s' , s ) \f$ (0 at the start \f$ s_0 \f$); \f$ s \f$ is
 * locally consistent if \f$ g( s ) = \mathrm{rhs}( s ) \f$, and the queue U
 * holds the inconsistent states with the key \f$ k( s ) = [ \min( g ,
 * \mathrm{rhs} ) + h ; \min( g , \mathrm{rhs} ) ] \f$, compared
 * lexicographically. The goals are those of the Domain, joined to a virtual
 * goal by edges of cost 0, so that \f$ m \f$, the least g of the goals
 * generated, stands for the g of the only goal of the paper, and the test of
 * its consistency for that of the goals that may give \f$ m \f$:
 *
 * \code{.unparsed}
 *  1  procedure UPDATEVERTEX(u)
 *  2     if u ≠ s₀ then rhs(u) ← min over the predecessors s' of u of
 *           g(s') + c(s', u), parent(u) ← the s' of the minimum
 *  3     U ← U ∖ {u};  if g(u) ≠ rhs(u) then insert u in U with k(u)
 *  4  procedure COMPUTESHORTESTPATH()
 *  5     while U ≠ ∅ and (the least key of U < [m; m] or a goal t with
 *           min(g(t), rhs(t)) ≤ m has g(t) ≠ rhs(t)) do
 *  6        u ← the state of the least key, removed from U
 *  7        if g(u) > rhs(u) then                       ▷ overconsistent
 *  8           g(u) ← rhs(u)
 *  9           for each successor v ≠ s₀ of u do
 * 10              if g(u) + c(u, v) < rhs(v) then
 * 11                 rhs(v) ← g(u) + c(u, v);  parent(v) ← u
 * 12              fix the place of v in U, as at line 3
 * 13        else                                         ▷ underconsistent
 * 14           g(u) ← ∞
 * 15           for each v ∈ {u} ∪ the successors of u with parent(v) = u
 *                 or v = u do UPDATEVERTEX(v)
 * 16  procedure LPA*()
 * 17     rhs(s₀) ← 0;  U ← {s₀};  COMPUTESHORTESTPATH()
 * 18     forever: wait for edge costs to change;  for each edge (u, v)
 *           changed, UPDATEVERTEX(v);  COMPUTESHORTESTPATH()
 * \endcode
 *
 * This is the optimized version of figure 7 of the paper (section 6), with
 * lines 10 and 11 in place of a full UPDATEVERTEX for an overconsistent
 * state, and line 15 recomputing rhs only where the parent was u. With a
 * consistent heuristic a search expands each state at most twice, once
 * overconsistent and once underconsistent [see Theorem 4 of the paper]. The
 * states are created when they are first reached, with \f$ g =
 * \mathrm{rhs} = \infty \f$, or when update() names them: a state never
 * reached has no predecessor of finite g, so that its rhs is \f$ \infty \f$
 * as if it had been created, until an edge into it appears.
 * U has no decrease-key: a state moved or removed leaves its entry, which
 * is stale when the state is no longer in U or its key has changed.
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

#ifndef __PBeFS_LPAStar
 #define __PBeFS_LPAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstdint>

#include <stdexcept>

#include <vector>

#include "ADS/priority/DaryHeap.h"

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
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// a callable accepting a predecessor, used only to state the concept
 template< class State , class Cost >
 struct pred_sink {
  void operator()( const State & , Cost ) const noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// a Domain that enumerates the predecessors of a state
/** for_each_predecessor( s , f ) calls f( t , c ) for every state t with
 * an edge ( t , s ), c being its cost; the edges are those that
 * for_each_successor() gives at the same time. */

template< class D >
concept PredecessorDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s ) {
  d.for_each_predecessor( s , detail::pred_sink< typename D::state_type ,
                                                 typename D::cost_type >{} );
  };

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS LPAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Lifelong Planning A*: a search that is repaired after a change of costs
/** The first solve() is a search from scratch; then the caller changes the
 * costs of some edges in the Domain (which LPAStar reads, and does not
 * copy), calls update( v ) for the head v of each edge changed, and the
 * next solve() repairs the search and returns the shortest path with the
 * new costs. The start is the initial state of the Domain when the LPAStar
 * is built, and the heuristic of a state must not change; it must stay
 * consistent with the new costs for the repair to expand few states, and
 * admissible for the path to be optimal. The counters of a Solution are
 * those of its solve(), but stored, which counts all the states created;
 * reexpanded counts the underconsistent expansions (line 13).
 * The constructor throws std::invalid_argument if the options give a weight
 * other than 1. */

template< PredecessorDomain D >
class LPAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the LPAStar)

 explicit LPAStar( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_table( options.expected_nodes
             ? NodeTable< LNode >::slots_for( options.expected_nodes )
             : 1024 ) {
  if( options.weight != 1 )
   throw( std::invalid_argument( "LPAStar: no weight" ) );
  const State s0 = domain.initial_state();
  f_start = node( s0 , Cost( domain.heuristic( s0 ) ) );
  f_start->rhs = Cost( 0 );
  enqueue( f_start );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// finds the shortest path with the current costs, repairing the last
 /// search after the calls of update()

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  SolutionT sol;
  f_st = SearchStats();
  compute();
  if( f_limit_hit )
   sol.status = Status::node_limit;
  else {
   const LNode * goal = best_goal();
   if( goal && ( goal->g < infinite_cost< Cost >() ) ) {
    sol.status = Status::solved;
    sol.cost = goal->g;
    sol.bound = goal->g;
    for( const LNode * n = goal ; n ; n = n->parent ) {
     sol.path.push_back( n->state );
     if( n == f_start )
      break;
     if( sol.path.size() > f_pool.size() )
      throw( std::logic_error( "LPAStar::solve: a cycle of parents" ) );
     }
    std::reverse( sol.path.begin() , sol.path.end() );
    }
   }
  f_st.stored = f_pool.size();
  f_st.expanded_by = { f_st.expanded };
  f_st.stored_by = { f_st.stored };
  f_st.seconds = clock.seconds();
  sol.stats = f_st;
  return( sol );
  }

/*---------------------------- THE CHANGES ---------------------------------*/
 /// tells that the costs of some edges entering v have changed (or that
 /// such edges have appeared or disappeared): UPDATEVERTEX( v ), v being
 /// created if it is new

 void update( const State & v ) {
  update_vertex( node( v , Cost( f_domain.heuristic( v ) ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a state with its g, rhs, h and parent, and whether it is in U
 struct LNode {
  State state;
  Cost g;
  Cost rhs;
  Cost h;
  LNode * parent;
  bool inq;
  bool goal;
  };

 /// a key of U: \f$ [ k_1 ; k_2 ] \f$
 struct Key {
  Cost k1;
  Cost k2;
  bool operator==( const Key & ) const = default;
  [[nodiscard]] bool operator<( const Key & o ) const noexcept {
   return( ( k1 < o.k1 ) || ( ( k1 == o.k1 ) && ( k2 < o.k2 ) ) );
   }
  };

 /// an entry of U
 struct Item {
  Key k;
  LNode * n;
  };

 /// the order of U: the smaller key first
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( a.k < b.k );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 // the key of n
 static Key key( const LNode * n ) noexcept {
  const Cost m = std::min( n->g , n->rhs );
  if( m == infinite_cost< Cost >() )
   return( Key{ m , m } );
  return( Key{ Cost( m + n->h ) , m } );
  }

 // the node of s, created with g = rhs = ∞ and the heuristic h
 LNode * node( const State & s , Cost h ) {
  auto [ n , fresh ] = f_table.find_or_insert( s , f_domain.hash( s ) ,
                                               [ & ] {
   LNode * x = f_pool.allocate();
   *x = LNode{ s , infinite_cost< Cost >() , infinite_cost< Cost >() , h ,
               nullptr , false , bool( f_domain.is_goal( s ) ) };
   return( x );
   } );
  if( fresh && n->goal )
   v_goals.push_back( n );
  return( n );
  }

 // the goal of least g, or nullptr
 const LNode * best_goal( void ) const {
  const LNode * b = nullptr;
  for( const LNode * x : v_goals )
   if( ! b || ( x->g < b->g ) )
    b = x;
  return( b );
  }

 // puts n in U if it is inconsistent, with its key [see line 3]
 void enqueue( LNode * n ) {
  n->inq = n->g != n->rhs;
  if( n->inq )
   f_u.push( Item{ key( n ) , n } );
  }

 // UPDATEVERTEX [see lines 1 to 3]
 void update_vertex( LNode * n ) {
  if( n != f_start ) {
   n->rhs = infinite_cost< Cost >();
   n->parent = nullptr;
   f_domain.for_each_predecessor( n->state , [ & ]( const State & t ,
                                                    Cost c ) {
    LNode * p = f_table.find( t , f_domain.hash( t ) );
    if( ! p || ( p->g == infinite_cost< Cost >() ) )
     return;
    if( p->g + c < n->rhs ) {
     n->rhs = p->g + c;
     n->parent = p;
     }
    } );
   }
  enqueue( n );
  }

 // the least valid entry of U, the stale ones dropped; false if none
 bool top( Item & it ) {
  while( ! f_u.empty() ) {
   it = f_u.top();
   if( it.n->inq && ( it.k == key( it.n ) ) )
    return( true );
   f_u.pop();
   }
  return( false );
  }

 // COMPUTESHORTESTPATH [see lines 4 to 15]
 void compute( void ) {
  Item it;
  while( top( it ) ) {
   const LNode * goal = best_goal();
   const Cost m = goal ? goal->g : infinite_cost< Cost >();
   bool repair = false;               // a goal of g or rhs within m is
   for( const LNode * x : v_goals )   // inconsistent [see line 5]
    repair = repair || ( ( x->g != x->rhs ) &&
                         ( std::min( x->g , x->rhs ) <= m ) );
   if( ! ( it.k < Key{ m , m } ) && ! repair )
    return;
   f_u.pop();
   LNode * u = it.n;
   u->inq = false;
   ++f_st.expanded;
   v_succ.clear();
   f_domain.for_each_successor( u->state , u->h , [ & ]( const State & t ,
                                                         Cost c , Cost ht ) {
    v_succ.push_back( { node( t , ht ) , c } );
    } );
   f_st.generated += v_succ.size();
   if( u->g > u->rhs ) {                       // overconsistent
    u->g = u->rhs;
    for( auto [ v , c ] : v_succ ) {
     if( v == f_start )
      continue;
     if( u->g + c < v->rhs ) {
      v->rhs = u->g + c;
      v->parent = u;
      enqueue( v );
      }
     }
    }
   else {                                       // underconsistent
    ++f_st.reexpanded;
    u->g = infinite_cost< Cost >();
    update_vertex( u );
    for( auto [ v , c ] : v_succ )
     if( v->parent == u )
      update_vertex( v );
    }
   if( f_options.node_limit && ( f_pool.size() >= f_options.node_limit ) ) {
    f_limit_hit = true;
    return;
    }
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 NodePool< LNode > f_pool;                           ///< the states
 NodeTable< LNode > f_table;                         ///< their table
 ADS::DaryHeap< Item , Better , 4 > f_u;             ///< the queue U
 std::vector< LNode * > v_goals;                     ///< the goals reached
 std::vector< std::pair< LNode * , Cost > > v_succ;  ///< of an expansion
 LNode * f_start = nullptr;                          ///< the start
 SearchStats f_st;                                   ///< the counters
 bool f_limit_hit = false;                           ///< the limit is hit

 };  // end( class( LPAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_LPAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File LPAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
