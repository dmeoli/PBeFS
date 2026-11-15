/*--------------------------------------------------------------------------*/
/*----------------------------- File LPAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template LPAStar, the Lifelong Planning A* of
 * Koenig, Likhachev and Furcy (AIJ 2004): an incremental search, which
 * after a change of the costs of some edges finds the new shortest path by
 * repairing the g values of the previous search instead of searching again;
 * for the class template DStarLite, the D* Lite of Koenig and Likhachev
 * (AAAI 2002), which is LPAStar on the reversed Domain with a target that
 * moves; for the class template RobotDomain, that reversed Domain; and
 * for the queues of their keys, LPAHeap and LPABuckets.
 *
 * Every state \f$ s \f$ has the value \f$ g( s ) \f$ and the one-step
 * lookahead \f$ \mathrm{rhs}( s ) = \min_{ s' \in \mathrm{pred}( s ) }
 * g( s' ) + c( s' , s ) \f$ (0 at the start \f$ s_0 \f$); \f$ s \f$ is
 * locally consistent if \f$ g( s ) = \mathrm{rhs}( s ) \f$, and the queue U
 * holds the inconsistent states with the key \f$ k( s ) = [ \min( g ,
 * \mathrm{rhs} ) + h + k_m ; \min( g , \mathrm{rhs} ) ] \f$, compared
 * lexicographically, \f$ k_m \f$ being 0 until the target moves. The goals
 * are those of the Domain, joined to a virtual target \f$ t^* \f$ by edges
 * of cost 0, so that \f$ B = \min_t g( t ) \f$ over the goals is the g of
 * \f$ t^* \f$, kept consistent at once, whose key is \f$ [ B + k_m ; B ]
 * \f$; a goal has \f$ h = 0 \f$, and an inconsistent goal of key equal to
 * that of \f$ t^* \f$ has to be repaired too, since the edge to \f$ t^* \f$
 * costs 0:
 *
 * \code{.unparsed}
 *  1  procedure UPDATEVERTEX(u)
 *  2     if u ≠ s₀ then rhs(u) ← min over the predecessors s' of u of
 *           g(s') + c(s', u), parent(u) ← the s' of the minimum
 *  3     U ← U ∖ {u};  if g(u) ≠ rhs(u) then insert u in U with k(u)
 *  4  procedure COMPUTESHORTESTPATH()
 *  5     while the least key of U < [B + k_m; B] or an inconsistent goal t
 *           has min(g(t), rhs(t)) ≤ B do
 *  6        u ← that state, removed from U;  if its key in U is below
 *           k(u) then insert u with k(u) and go back to line 5
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
 * 16  procedure EDGECHANGED(u, v, from, to)
 * 17     if to < from and v ≠ s₀ and g(u) + to < rhs(v) then
 *           rhs(v) ← g(u) + to;  parent(v) ← u;  fix the place of v in U
 * 18     if to > from and parent(v) = u then UPDATEVERTEX(v)
 * 19  procedure LPA*()
 * 20     rhs(s₀) ← 0;  U ← {s₀};  COMPUTESHORTESTPATH()
 * 21     forever: wait for edge costs to change;  for each edge (u, v)
 *           changed, EDGECHANGED(u, v, from, to);  COMPUTESHORTESTPATH()
 * \endcode
 *
 * This is the optimized version of figure 7 of the paper (section 6), with
 * lines 10 and 11 in place of a full UPDATEVERTEX for an overconsistent
 * state, line 15 recomputing rhs only where the parent was u, and lines 17
 * and 18 those of {33'}-{44'}, which need the old cost of the edge: a
 * decrease costs \f$ O( 1 ) \f$, an increase a scan of the predecessors
 * only where the edge was the one of the parent. With a consistent
 * heuristic a search expands each state at most twice, once overconsistent
 * and once underconsistent [see Theorem 4 of the paper]. The states are
 * created when they are first reached, with \f$ g = \mathrm{rhs} = \infty
 * \f$, or when a change names them: a state never reached has no
 * predecessor of finite g, so that its rhs is \f$ \infty \f$ as if it had
 * been created, until an edge into it appears. The goals of finite g are
 * in an indexed heap by g, whose least key is B, and the inconsistent ones
 * in another by \f$ \min( g , \mathrm{rhs} ) \f$, so that the test of line
 * 5 costs \f$ O( 1 ) \f$ whatever the number of goals.
 *
 * U has no decrease-key: a state moved or removed leaves its entry, which
 * is stale when the state is no longer in U or when its key in U is no
 * longer the last one it was pushed with; the latter is kept in the state,
 * so that line 6, that of D* Lite [see figure 5 of Koenig and Likhachev,
 * IEEE T-RO 2005, lines {13''}-{14''}], finds an entry whose key the growth of
 * \f$ k_m \f$ and the change of the heuristic have raised, and pushes it
 * again.
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

#include <type_traits>

#include <utility>

#include <vector>

#include "ADS/priority/BucketQueue.h"

#include "ADS/priority/DaryHeap.h"

#include "ADS/priority/IndexedDaryHeap.h"

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
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// \f$ a + b \f$, \f$ \infty \f$ if either is [see infinite_cost()]

template< class Cost >
[[nodiscard]] constexpr Cost add_cost( Cost a , Cost b ) noexcept
{
 if( ( a == infinite_cost< Cost >() ) || ( b == infinite_cost< Cost >() ) )
  return( infinite_cost< Cost >() );
 return( Cost( a + b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// sets the key of id in the ADS::IndexedDaryHeap h to k, inserting it if
/// absent, the ids grown to id + 1 if fewer

template< class Heap , class K >
void set_heap_key( Heap & h , int id , const K & k )
{
 if( id >= h.ids() )
  h.grow( id + 1 );
 if( h.contains( id ) )
  h.update( id , k );
 else
  h.push( id , k );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// removes id from h if it is there

template< class Heap >
void erase_heap_key( Heap & h , int id )
{
 if( ( id >= 0 ) && ( id < h.ids() ) && h.contains( id ) )
  h.erase( id );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LPAKey ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a key \f$ [ k_1 ; k_2 ] \f$ of LPA*, compared lexicographically

template< class Cost >
struct LPAKey {

 Cost k1;                                          ///< the first key
 Cost k2;                                          ///< the second key

 /// the equality of both keys
 bool operator==( const LPAKey & ) const = default;

 /// the lexicographic order
 [[nodiscard]] bool operator<( const LPAKey & o ) const noexcept {
  return( ( k1 < o.k1 ) || ( ( k1 == o.k1 ) && ( k2 < o.k2 ) ) );
  }

 /// the lexicographic order or the equality
 [[nodiscard]] bool operator<=( const LPAKey & o ) const noexcept {
  return( ! ( o < *this ) );
  }

 };  // end( struct( LPAKey ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS LPAHeap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the queue of the keys of LPA*, a 4-ary heap with no decrease-key
/** push( x , k ), top() and pop() of the entries \f$ ( k , x ) \f$, least
 * key first, the ties of the key in no particular order. */

template< class T , class Cost >
class LPAHeap {

 public:

 /// an entry
 struct Entry {
  LPAKey< Cost > k;                                ///< its key
  T x;                                             ///< its item
  };

 /// pushes x with the key k
 void push( T x , LPAKey< Cost > k ) { f_h.push( Entry{ k , x } ); }

 /// the entry of least key
 [[nodiscard]] const Entry & top( void ) const { return( f_h.top() ); }

 /// removes the entry of least key
 void pop( void ) { static_cast< void >( f_h.pop() ); }

 /// whether there is no entry
 [[nodiscard]] bool empty( void ) const { return( f_h.empty() ); }

 /// the number of the entries
 [[nodiscard]] std::size_t size( void ) const { return( f_h.size() ); }

 /// removes every entry
 void clear( void ) { f_h.clear(); }

 private:

 // the smaller key first
 struct Better {
  bool operator()( const Entry & a , const Entry & b ) const noexcept {
   return( a.k < b.k );
   }
  };

 ADS::DaryHeap< Entry , Better , 4 > f_h;          ///< the heap

 };  // end( class( LPAHeap ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS LPABuckets ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the queue of the keys of LPA* with integer costs, by buckets
/** The ADS::BucketQueue of the least first key and, among them, of the
 * least second key (ADS::LeastSecondKey), the order of theorem 6 of
 * Koenig, Likhachev and Furcy (AIJ 2004), the last pushed first among
 * equal keys; the least first key is found by a scan that only moves
 * forward between two pushes of a smaller one, so that the pushes below
 * the last key extracted, which the decreases of an incremental search
 * make, cost nothing more. */

template< class T , class Cost >
 requires std::is_integral_v< Cost >
class LPABuckets {

 public:

 /// an entry
 struct Entry {
  LPAKey< Cost > k;                                ///< its key
  T x;                                             ///< its item
  };

 /// pushes x with the key k, \f$ k_1 , k_2 \geq 0 \f$
 void push( T x , LPAKey< Cost > k ) { f_q.push( x , k.k1 , k.k2 ); }

 /// the entry of least key
 [[nodiscard]] Entry top( void ) {
  const Cost k1 = f_q.topF();
  const Cost k2 = f_q.topG();
  return( Entry{ LPAKey< Cost >{ k1 , k2 } , f_q.top() } );
  }

 /// removes the entry of least key
 void pop( void ) { static_cast< void >( f_q.pop() ); }

 /// whether there is no entry
 [[nodiscard]] bool empty( void ) const { return( f_q.empty() ); }

 /// the number of the entries
 [[nodiscard]] std::size_t size( void ) const { return( f_q.size() ); }

 /// removes every entry
 void clear( void ) { f_q.clear(); }

 private:

 ADS::BucketQueue< T , Cost , ADS::LeastSecondKey > f_q;   ///< the queue

 };  // end( class( LPABuckets ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS RobotDomain ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the Domain D with its edges reversed, from a goal to a target that moves
/** The successors of a state are its predecessors in D and vice versa, with
 * the same costs; the initial state is the goal \f$ s_g \f$ of the robot,
 * the only goal is the state \f$ s_r \f$ of the robot, and the heuristic of
 * \f$ s \f$ is pairwise_heuristic( \f$ s_r , s \f$ ) of D, an estimate of
 * the cost from \f$ s_r \f$ to \f$ s \f$ in D, which has to be consistent
 * and to satisfy the triangle inequality, so that the keys stored before a
 * move of the robot stay lower bounds once \f$ k_m \f$ grows by the
 * estimate of the move [see Koenig and Likhachev, IEEE T-RO 2005, figure
 * 5, line {30''}].
 * The events of D are handed with their edges reversed; the hashes, the
 * features, the abstract states and the perfect hash are those of D. This
 * is the search space of D* Lite [see DStarLite]. */

template< class D >
 requires PredecessorDomain< D > && PairwiseDomain< D >
class RobotDomain {

 public:

 using state_type = state_t< D >;                 ///< the state type
 using cost_type = cost_t< D >;                   ///< the cost type

 /// constructor, taking D (which must outlive the RobotDomain), the goal
 /// of the robot and its state

 RobotDomain( D & d , const state_type & goal , const state_type & robot )
  : f_d( & d ) , f_goal( goal ) , f_robot( robot ) {}

 /// the goal of the robot

 [[nodiscard]] state_type initial_state( void ) const { return( f_goal ); }

 /// true if s is the state of the robot

 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( s == f_robot );
  }

 /// the hash of s in D

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d->hash( s ) );
  }

 /// the estimate of the cost from the robot to s in D

 [[nodiscard]] cost_type heuristic( const state_type & s ) const {
  return( cost_type( f_d->pairwise_heuristic( f_robot , s ) ) );
  }

 /// the estimate of the cost from t to s in D

 [[nodiscard]] cost_type pairwise_heuristic( const state_type & s ,
                                             const state_type & t ) const {
  return( cost_type( f_d->pairwise_heuristic( t , s ) ) );
  }

 /// calls f( t , c , h( t ) ) for each predecessor t of s in D, c being
 /// the cost of ( t , s )

 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const
 {
  f_d->for_each_predecessor( s , [ & ]( const state_type & t , cost_type c ) {
   f( t , c , heuristic( t ) );
   } );
  }

 /// calls f( t , c ) for each successor t of s in D

 template< class F >
 void for_each_predecessor( const state_type & s , F && f ) const {
  f_d->for_each_successor( s , cost_type( f_d->heuristic( s ) ) ,
                           [ & ]( const state_type & t , cost_type c ,
                                  cost_type ) {
   f( t , c );
   } );
  }

 /// the edge events of D since the last call, reversed

 template< class F >
 void drain_events( F && f ) requires EventDomain< D > {
  f_d->drain_events( [ & ]( const auto & e ) {
   if constexpr( requires { e.from; } )
    f( EdgeCostChanged< state_type , cost_type >{ e.v , e.u , e.from ,
                                                  e.to } );
   } );
  }

 /// the number of features of D

 [[nodiscard]] std::size_t num_features( void ) const
  requires FeatureDomain< D > {
  return( f_d->num_features() );
  }

 /// the features of s in D

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const
  requires FeatureDomain< D > {
  f_d->for_each_feature( s , std::forward< F >( f ) );
  }

 /// the perfect hash of s in D

 [[nodiscard]] std::uint64_t perfect_hash( const state_type & s ) const
  requires PerfectHashDomain< D > {
  return( f_d->perfect_hash( s ) );
  }

 /// the abstract state of s in D

 [[nodiscard]] std::uint64_t abstract_state( const state_type & s ) const
  requires AbstractStateDomain< D > {
  return( f_d->abstract_state( s ) );
  }

 /// the state of the robot

 [[nodiscard]] const state_type & robot( void ) const { return( f_robot ); }

 /// moves the robot to r

 void set_robot( const state_type & r ) { f_robot = r; }

 /// the Domain reversed

 [[nodiscard]] D & domain( void ) const { return( *f_d ); }

 private:

 D * f_d;                                          ///< the Domain
 state_type f_goal;                                ///< the goal
 state_type f_robot;                               ///< the robot

 };  // end( class( RobotDomain ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS LPAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Lifelong Planning A*: a search that is repaired after a change of costs
/** The first solve() is a search from scratch; then the caller changes the
 * costs of some edges in the Domain (which LPAStar reads, and does not
 * copy), tells LPAStar of each change by edge_changed( u , v , from , to )
 * (or, for an EventDomain, by take_events() on the Domain), and of each
 * goal added or removed by goal_changed(), and the next solve() repairs the
 * search and returns the shortest path with the new costs. The changes of
 * one edge are told once, from the cost before the first to the cost after
 * the last, after all of them are made in the Domain; update( v ), which
 * recomputes rhs( v ) over all the predecessors, serves when the old costs
 * are not known. The start is the initial state of the Domain when the
 * LPAStar is built; the heuristic must stay consistent with the new costs
 * for the repair to expand few states, and admissible for the path to be
 * optimal. shift_keys( delta ) adds delta to \f$ k_m \f$ and tells that
 * the heuristic has changed, so that each state computes its h again: the
 * move of the target of D* Lite [see DStarLite]. The counters of a Solution
 * are those of its solve() and of the changes told before it, but stored,
 * which counts all the states created; reexpanded and under count the
 * underconsistent expansions (line 13), rescans the calls of
 * UPDATEVERTEX, delta the states whose g differs from the one before. Queue
 * is LPAHeap (the default) or, for integer costs, LPABuckets. The
 * constructor throws std::invalid_argument if the options give a weight
 * other than 1. */

template< PredecessorDomain D ,
          template< class , class > class Queue = LPAHeap >
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
 /// search after the changes told

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  SolutionT sol;
  compute();
  if( f_limit_hit )
   sol.status = Status::node_limit;
  else
   if( ! f_mu.empty() ) {
    const LNode * goal = v_goals[ std::size_t( f_mu.top().id ) ];
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
  for( LNode * n : v_touched ) {
   if( n->g != n->g0 )
    ++f_st.delta;
   n->touched = false;
   }
  v_touched.clear();
  f_st.stored = f_pool.size();
  f_st.expanded_by = { f_st.expanded };
  f_st.stored_by = { f_st.stored };
  f_st.seconds = clock.seconds();
  sol.stats = f_st;
  f_st = SearchStats();
  return( sol );
  }

/*------------------------------ THE CHANGES -------------------------------*/
 /// tells that the cost of the edge ( u , v ) has changed from from to to,
 /// \f$ \infty \f$ standing for an absent edge [see lines 16 to 18]

 void edge_changed( const State & u , const State & v , Cost from , Cost to )
 {
  if( from == to )
   return;
  LNode * un = f_table.find( u , f_domain.hash( u ) );
  if( to < from ) {
   if( ! un || ( un->g == infinite_cost< Cost >() ) )
    return;
   LNode * vn = node( v , Cost( f_domain.heuristic( v ) ) );
   const Cost x = add_cost( un->g , to );
   if( ( vn != f_start ) && ( x < vn->rhs ) ) {
    vn->rhs = x;
    vn->parent = un;
    enqueue( vn );
    }
   return;
   }
  LNode * vn = f_table.find( v , f_domain.hash( v ) );
  if( vn && un && ( vn->parent == un ) )
   update_vertex( vn );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// edge_changed() for each event of the EventDomain d, and goal_changed()
 /// for each goal added or removed

 template< class E >
 void take_events( E & d ) requires EventDomain< E > {
  d.drain_events( [ & ]( const auto & e ) {
   if constexpr( requires { e.from; } )
    edge_changed( e.u , e.v , e.from , e.to );
   else
    goal_changed( e.t , e.added );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells that the costs of some edges entering v have changed, with no
 /// old cost: UPDATEVERTEX( v ), v being created if it is new

 void update( const State & v ) {
  update_vertex( node( v , Cost( f_domain.heuristic( v ) ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// tells that t has become a goal of the Domain (added) or no longer is

 void goal_changed( const State & t , bool added ) {
  LNode * n = added ? node( t , Cost( f_domain.heuristic( t ) ) )
                    : f_table.find( t , f_domain.hash( t ) );
  if( ! n || ( n->goal == added ) )
   return;
  if( added && ( n->gid < 0 ) ) {
   n->gid = int( v_goals.size() );
   v_goals.push_back( n );
   }
  n->goal = added;
  if( added )
   set_g( n , n->g );
  else {
   erase_heap_key( f_mu , n->gid );
   erase_heap_key( f_gu , n->gid );
   }
  enqueue( n );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds delta to \f$ k_m \f$ and computes h again at the next use of each
 /// state: the target of D* Lite has moved [see DStarLite]

 void shift_keys( Cost delta ) {
  f_km = add_cost( f_km , delta );
  ++f_hv;
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Key = LPAKey< Cost >;

 /// a state with its g, rhs, h and parent, the key it was last pushed in
 /// U with, its g before the solve() and whether it is in U and a goal
 struct LNode {
  State state;
  Cost g;
  Cost rhs;
  Cost h;
  Key q;                                          // the key pushed last
  Cost g0;                                        // g before the solve()
  LNode * parent;
  std::uint32_t hv;                               // the version of h
  int gid;                                        // the index of a goal
  bool inq;
  bool goal;
  bool touched;                                   // g0 is set
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 // h of n, computed again if the heuristic changed since
 Cost heur( LNode * n ) {
  if( n->hv != f_hv ) {
   n->h = Cost( f_domain.heuristic( n->state ) );
   n->hv = f_hv;
   }
  return( n->h );
  }

 // the key of n
 Key key( LNode * n ) {
  const Cost m = std::min( n->g , n->rhs );
  if( m == infinite_cost< Cost >() )
   return( Key{ m , m } );
  return( Key{ add_cost( Cost( m + heur( n ) ) , f_km ) , m } );
  }

 // the node of s, created with g = rhs = ∞ and the heuristic h
 LNode * node( const State & s , Cost h ) {
  auto [ n , fresh ] = f_table.find_or_insert( s , f_domain.hash( s ) ,
                                               [ & ] {
   LNode * x = f_pool.allocate();
   *x = LNode{ s , infinite_cost< Cost >() , infinite_cost< Cost >() , h ,
               Key{} , Cost() , nullptr , f_hv , -1 , false ,
               bool( f_domain.is_goal( s ) ) , false };
   return( x );
   } );
  if( fresh && n->goal ) {
   n->gid = int( v_goals.size() );
   v_goals.push_back( n );
   }
  return( n );
  }

 // g( n ) ← x, kept in the heap of the goals
 void set_g( LNode * n , Cost x ) {
  if( ! n->touched ) {
   n->touched = true;
   n->g0 = n->g;
   v_touched.push_back( n );
   }
  n->g = x;
  if( n->goal ) {
   if( x == infinite_cost< Cost >() )
    erase_heap_key( f_mu , n->gid );
   else
    set_heap_key( f_mu , n->gid , x );
   }
  }

 // puts n in U if it is inconsistent, with its key [see line 3], and a
 // goal in the heap of the inconsistent ones
 void enqueue( LNode * n ) {
  const bool inc = n->g != n->rhs;
  if( inc ) {
   const Key k = key( n );
   if( ! n->inq || ! ( k == n->q ) ) {
    f_u.push( n , k );
    n->q = k;
    }
   }
  n->inq = inc;
  if( n->goal ) {
   if( inc )
    set_heap_key( f_gu , n->gid , std::min( n->g , n->rhs ) );
   else
    erase_heap_key( f_gu , n->gid );
   }
  }

 // UPDATEVERTEX [see lines 1 to 3]
 void update_vertex( LNode * n ) {
  ++f_st.rescans;
  if( n != f_start ) {
   n->rhs = infinite_cost< Cost >();
   n->parent = nullptr;
   f_domain.for_each_predecessor( n->state , [ & ]( const State & t ,
                                                    Cost c ) {
    LNode * p = f_table.find( t , f_domain.hash( t ) );
    if( ! p )
     return;
    const Cost x = add_cost( p->g , c );
    if( x < n->rhs ) {
     n->rhs = x;
     n->parent = p;
     }
    } );
   }
  enqueue( n );
  }

 // the least valid entry of U, the stale ones dropped and those whose key
 // has risen pushed again [see line 6]; nullptr if none, k its key
 LNode * top( Key & k ) {
  while( ! f_u.empty() ) {
   const auto & e = f_u.top();
   LNode * n = e.x;
   if( ! n->inq || ! ( e.k == n->q ) ) {
    f_u.pop();
    continue;
    }
   k = key( n );
   if( n->q < k ) {
    f_u.pop();
    f_u.push( n , k );
    n->q = k;
    continue;
    }
   return( n );
   }
  return( nullptr );
  }

 // COMPUTESHORTESTPATH [see lines 4 to 15]
 void compute( void ) {
  for( ; ; ) {
   const Cost b = f_mu.empty() ? infinite_cost< Cost >() : f_mu.top().key;
   const Key kt{ add_cost( b , f_km ) , b };
   Key k;
   LNode * u = top( k );
   if( u && ( k < kt ) )
    f_u.pop();
   else
    if( ! f_gu.empty() && ( f_gu.top().key <= b ) )
     u = v_goals[ std::size_t( f_gu.top().id ) ];
    else
     return;
   u->inq = false;
   if( u->goal )
    erase_heap_key( f_gu , u->gid );
   ++f_st.expanded;
   v_succ.clear();
   f_domain.for_each_successor( u->state , heur( u ) ,
                                [ & ]( const State & t , Cost c , Cost ht ) {
    v_succ.push_back( { node( t , ht ) , c } );
    } );
   f_st.generated += v_succ.size();
   if( u->g > u->rhs ) {                       // overconsistent
    set_g( u , u->rhs );
    for( auto [ v , c ] : v_succ ) {
     if( v == f_start )
      continue;
     const Cost x = add_cost( u->g , c );
     if( x < v->rhs ) {
      v->rhs = x;
      v->parent = u;
      enqueue( v );
      }
     }
    }
   else {                                       // underconsistent
    ++f_st.reexpanded;
    ++f_st.under;
    set_g( u , infinite_cost< Cost >() );
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
 Queue< LNode * , Cost > f_u;                        ///< the queue U
 ADS::IndexedDaryHeap< Cost > f_mu;                  ///< g of the goals
 ADS::IndexedDaryHeap< Cost > f_gu;                  ///< inconsistent goals
 std::vector< LNode * > v_goals;                     ///< the goals reached
 std::vector< std::pair< LNode * , Cost > > v_succ;  ///< of an expansion
 std::vector< LNode * > v_touched;                   ///< g0 set
 LNode * f_start = nullptr;                          ///< the start
 SearchStats f_st;                                   ///< the counters
 Cost f_km = Cost( 0 );                              ///< k_m
 std::uint32_t f_hv = 0;                             ///< version of h
 bool f_limit_hit = false;                           ///< the limit is hit

 };  // end( class( LPAStar ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS DStarLite -----------------------------*/
/*--------------------------------------------------------------------------*/
/// D* Lite: the shortest path of a robot that moves, repaired after each
/// move and change of costs
/** LPAStar on the RobotDomain of D, from the goal of the robot to its
 * state [see Koenig and Likhachev, AAAI 2002, figure 4, the optimized
 * version]: solve() returns the shortest path from the robot to its goal
 * with the current costs, which the caller changes in D and tells as to
 * LPAStar, with the edges of D; move_to( r ) moves the robot to r, which
 * adds the estimate of the move to \f$ k_m \f$ and makes r the target of
 * the search. The heuristic of D must be pairwise [see RobotDomain]. */

template< class D , template< class , class > class Queue = LPAHeap >
 requires PredecessorDomain< D > && PairwiseDomain< D >
class DStarLite {

 public:

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// constructor, taking D (which must outlive the DStarLite), the goal and
 /// the state of the robot

 DStarLite( D & d , const State & goal , const State & robot ,
            const SearchOptions & options = {} )
  : f_r( d , goal , robot ) , f_lpa( f_r , options ) {}

 /// the shortest path from the robot to its goal, path.front() being the
 /// robot

 [[nodiscard]] SolutionT solve( void ) {
  SolutionT s = f_lpa.solve();
  std::reverse( s.path.begin() , s.path.end() );
  return( s );
  }

 /// tells that the cost of the edge ( u , v ) of D has changed

 void edge_changed( const State & u , const State & v , Cost from , Cost to )
 {
  f_lpa.edge_changed( v , u , from , to );
  }

 /// edge_changed() for each event of D

 void take_events( void ) requires EventDomain< D > {
  f_lpa.take_events( f_r );
  }

 /// moves the robot to r

 void move_to( const State & r ) {
  const State a = f_r.robot();
  if( a == r )
   return;
  const Cost delta = Cost( f_r.domain().pairwise_heuristic( a , r ) );
  f_r.set_robot( r );
  f_lpa.goal_changed( a , false );
  f_lpa.goal_changed( r , true );
  f_lpa.shift_keys( delta );
  }

 /// the state of the robot

 [[nodiscard]] const State & robot( void ) const { return( f_r.robot() ); }

 private:

 RobotDomain< D > f_r;                             ///< the search space
 LPAStar< RobotDomain< D > , Queue > f_lpa;        ///< the search

 };  // end( class( DStarLite ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_LPAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File LPAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
