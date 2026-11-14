/*--------------------------------------------------------------------------*/
/*--------------------------- File AnytimeAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template AnytimeAStar, the sequential anytime
 * searches: Anytime Weighted A* (AWA*) of Hansen and Zhou (JAIR 2007) and
 * Anytime Repairing A* (ARA*) of Likhachev, Gordon and Thrun (NIPS 2003).
 * Both find a first path quickly with a weight \f$ w > 1 \f$ on the
 * heuristic and improve it until they prove it optimal, recording each path
 * with the lower bound on \f$ C^* \f$ known when it is found.
 *
 * AWA* is weighted A* that does not stop at the first goal, with the
 * notation of AStar.h, \f$ z \f$ the incumbent and FQ a second OPEN list
 * ordered by \f$ f \f$ over the same nodes, whose least \f$ f \f$ is the
 * lower bound:
 *
 * \code{.unparsed}
 *  1  procedure AWA*(s₀, w)
 *  2     z ← ∞;  push s₀ in OPEN (by f' = g + w h) and in FQ (by f)
 *  3     while OPEN ≠ ∅ do
 *  4        u ← the first node of OPEN, removed (stale entries dropped)
 *  5        if f(u) ≥ z then continue                       ▷ pruned
 *  6        close u
 *  7        if u is a goal then
 *  8           z ← g(u);  record (z, min(z, min f of FQ));  continue
 *  9        for each successor v of u with g(u) + c(u, v) + h(v) < z do
 * 10           insert v as A* does (lines 10 to 16 of AStar.h), in OPEN
 *                 and in FQ
 * 11     return the path of cost z, optimal
 * \endcode
 *
 * The goal is tested when it is selected, as in AStar, whereas Algorithm 1
 * of Hansen and Zhou tests it when it is generated and does not put it in
 * OPEN; the first goal selected costs at most \f$ w C^* \f$, as in weighted
 * A*;
 * every node of \f$ f < z \f$ is eventually expanded, with reopening, so
 * that the search ends with \f$ z = C^* \f$ [see Theorem 1 of Hansen and
 * Zhou]. A node of FQ is stale when it is closed or its g changed, as in
 * OPEN, and since a node on an optimal path with \f$ g = g^* \f$ is always
 * in OPEN (or pruned, or \f$ z = C^* \f$), \f$ \min( z , \min_{\rm FQ} f )
 * \leq C^* \f$ at line 8.
 *
 * ARA* runs weighted A* with a decreasing schedule of weights
 * \f$ w_0 > w_1 > \dots \f$, ending with 1, and reuses the work of each
 * iteration in the next. Within an iteration a node is expanded at most
 * once: a closed node whose g decreases goes to the list INCONS instead of
 * OPEN, and at the end of the iteration OPEN \f$ \cup \f$ INCONS, ordered by
 * the next weight, is the OPEN of the next one, CLOSED being emptied; the
 * goal is tested when it is generated, and \f$ z \f$ is its g:
 *
 * \code{.unparsed}
 *  1  procedure IMPROVEPATH(w)
 *  2     while OPEN has a node of f' < z do
 *  3        u ← the first node of OPEN, removed (stale entries dropped)
 *  4        if f(u) ≥ z then continue                       ▷ pruned
 *  5        close u
 *  6        for each successor v of u with g' = g(u) + c(u, v), g' + h(v) < z
 *  7           if v is new or g' < g(v) then
 *  8              g(v) ← g';  parent(v) ← u
 *  9              if v is a goal then z ← g', recording v
 * 10              else if v is closed then add v to INCONS
 * 11              else push (g' + w h(v), g', v) in OPEN
 * 12  procedure ARA*(s₀, w₀, w₁, ...)
 * 13     z ← ∞;  OPEN ← {s₀};  i ← 0
 * 14     loop
 * 15        IMPROVEPATH(w_i)
 * 16        L ← min(z, min f over OPEN ∪ INCONS)
 * 17        if z decreased in IMPROVEPATH then record (z, L)
 * 18        if L ≥ z or OPEN ∪ INCONS = ∅ then return
 * 19        i ← i + 1 (w_i = 1 once the schedule ends)
 * 20        OPEN ← OPEN ∪ INCONS, by g + w_i h;  INCONS ← ∅;  CLOSED ← ∅
 * \endcode
 *
 * By Theorem 1 of Likhachev et al., with a consistent heuristic, at the
 * end of IMPROVEPATH(w) every state s with \f$ g( s ) + w h( s ) \f$ not
 * larger than the least key of OPEN has \f$ g( s ) \leq w g^*( s ) \f$, and
 * in particular \f$ z \leq w C^* \f$; the bound of line 16 holds with any
 * admissible heuristic, since a node of an optimal path with
 * \f$ g = g^* \f$ is in OPEN or in INCONS, or has been pruned with
 * \f$ f \geq z \f$, and \f$ z / L \f$ is the suboptimality \f$ \epsilon' \f$
 * of the paper. OPEN \f$ \cup \f$ INCONS is found at line 16 by emptying
 * OPEN, whose live entries are then pushed again with the keys of the next
 * weight, so that no second queue is needed.
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

#ifndef __PBeFS_AnytimeAStar
 #define __PBeFS_AnytimeAStar
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
/*-------------------------- CLASS AnytimeOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of AnytimeAStar

struct AnytimeOptions : public SearchOptions {

 /// the weights of the iterations of ARA*, decreasing, the iterations
 /// after the last one having the weight 1; empty for AWA*, which keeps the
 /// weight of the SearchOptions for the whole search
 std::vector< double > schedule;

 };  // end( struct( AnytimeOptions ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS AnytimeAStar ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the anytime weighted searches AWA* and ARA*
/** The search is AWA* with the weight of the options if the schedule is
 * empty, ARA* with the schedule otherwise. Each path that lowers the
 * incumbent is recorded in Solution::improvements, with the lower bound on
 * \f$ C^* \f$ known at that moment; the search ends with an optimal path,
 * bound being its cost, unless it reaches the node limit, in which case the
 * Solution keeps the best path found, with Status::node_limit and the last
 * bound known. drop_duplicates and anytime of the options are not read:
 * ARA* has its own rule for the closed nodes (INCONS), and AWA* reopens
 * them to reach the optimum. The constructor throws std::invalid_argument
 * if a weight is below 1 or the schedule is not decreasing. */

template< Domain D ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class AnytimeAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the AnytimeAStar)

 explicit AnytimeAStar( const D & domain , const AnytimeOptions & options =
                        {} )
  : f_domain( domain ) , f_options( options ) {
  for( std::size_t i = 0 ; i < f_options.schedule.size() ; ++i )
   if( ! ( f_options.schedule[ i ] >= 1 ) ||
       ( i && ! ( f_options.schedule[ i ] < f_options.schedule[ i - 1 ] ) ) )
    throw( std::invalid_argument( "AnytimeAStar: the schedule must decrease"
                                  " and stay at least 1" ) );
  static_cast< void >( Weight< Cost >( f_options.weight ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Run r( *this );
  if( f_options.schedule.empty() )
   r.awa();
  else
   r.ara();
  return( r.finish() );
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the state of a search

 struct Run {

  explicit Run( const AnytimeAStar & a )
   : f_a( a ) , d( a.f_domain ) ,
     table( a.f_options.expected_nodes
            ? NodeTable< NodeT >::slots_for( a.f_options.expected_nodes )
            : 1024 ) {
   const State s0 = d.initial_state();
   root = table.find_or_insert( s0 , d.hash( s0 ) , [ & ] {
    NodeT * n = pool.allocate();
    *n = NodeT{ s0 , Cost( 0 ) , Cost( d.heuristic( s0 ) ) , nullptr ,
                false , false };
    return( n );
    } ).first;
   bound = root->h;
   }

  // true if the node limit is reached
  bool limit( void ) {
   if( f_a.f_options.node_limit &&
       ( pool.size() >= f_a.f_options.node_limit ) )
    limit_hit = true;
   return( limit_hit );
   }

  // closes n and counts its expansion
  void close( NodeT * n ) {
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   }

  // the successors of n, with the slots of their lookups prefetched
  void children( NodeT * n ) {
   kids.clear();
   for_each_child( d , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    const std::uint64_t hash = d.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    } );
   st.generated += kids.size();
   }

  // records a path to the goal n of cost z, with the bound lb
  void record( NodeT * n , Cost lb ) {
   z = n->g;
   goal = n;
   bound = std::min( lb , z );
   imp.push_back( { z , bound , st.expanded , clock.seconds() } );
   }

  // the least f of FQ, its stale entries dropped (infinite if empty)
  Cost least_f( void ) {
   while( ! fq.empty() && fq.top().stale() )
    fq.pop();
   return( fq.empty() ? infinite_cost< Cost >() : fq.top_f() );
   }

  // AWA* [see lines 1 to 11 of the pseudocode]
  void awa( void ) {
   const Weight< Cost > w( f_a.f_options.weight , true );
   open.push( root , Cost( 0 ) , w.key( Cost( 0 ) , root->h ) );
   fq.push( root );
   while( ! open.empty() ) {
    const OpenEntry< NodeT > e = open.pop();
    if( e.stale() || ( e.g + e.node->h >= z ) )
     continue;
    NodeT * n = e.node;
    close( n );
    if( d.is_goal( n->state ) ) {
     record( n , least_f() );
     continue;
     }
    children( n );
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
      if( ! ( g < m->g ) ) {
       ++st.duplicates;
       continue;
       }
      m->g = g;
      m->parent = n;
      m->op = k.op;
      m->closed = false;
      }
     open.push( m , g , w.key( g , k.h ) );
     fq.push( m , g );
     }
    if( limit() )
     return;
    }
   optimal = true;
   }

  // IMPROVEPATH of ARA* with the weight w [see lines 1 to 11]
  void improve( const Weight< Cost > & w ) {
   while( ! open.empty() ) {
    const OpenEntry< NodeT > e = open.top();
    if( e.stale() ) {
     open.pop();
     continue;
     }
    if( ! w.live( open.top_f() , z ) )
     return;
    open.pop();
    NodeT * n = e.node;
    if( e.g + n->h >= z )
     continue;
    close( n );
    v_closed.push_back( n );
    children( n );
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
      if( ! ( g < m->g ) ) {
       ++st.duplicates;
       continue;
       }
      m->g = g;
      m->parent = n;
      m->op = k.op;
      }
     if( d.is_goal( k.s ) ) {
      z = g;
      goal = m;
      }
     else
      if( m->closed )
       v_incons.push_back( m );
      else
       open.push( m , g , w.key( g , k.h ) );
     }
    if( limit() )
     return;
    }
   }

  // ARA* [see lines 12 to 20 of the pseudocode]
  void ara( void ) {
   const auto & sch = f_a.f_options.schedule;
   std::size_t i = 0;
   Weight< Cost > w( sch[ 0 ] );
   if( d.is_goal( root->state ) ) {
    z = Cost( 0 );
    goal = root;
    }
   else
    open.push( root , Cost( 0 ) , w.key( Cost( 0 ) , root->h ) );
   Cost reported = infinite_cost< Cost >();
   for( ; ; ) {
    improve( w );
    if( limit_hit )
     break;
    // OPEN ∪ INCONS, OPEN emptied, and L
    std::vector< NodeT * > next;
    Cost lb = z;
    while( ! open.empty() ) {
     const OpenEntry< NodeT > e = open.pop();
     if( e.stale() )
      continue;
     next.push_back( e.node );
     lb = std::min( lb , Cost( e.g + e.node->h ) );
     }
    for( NodeT * n : v_incons )
     lb = std::min( lb , n->f() );
    bound = lb;
    if( goal && ( z < reported ) ) {
     imp.push_back( { z , lb , st.expanded , clock.seconds() } );
     reported = z;
     }
    if( ( lb >= z ) || ( next.empty() && v_incons.empty() ) ) {
     optimal = true;
     return;
     }
    if( ++i < sch.size() )
     w = Weight< Cost >( sch[ i ] );
    else
     w = Weight< Cost >( 1 );
    for( NodeT * n : v_closed )
     n->closed = false;
    v_closed.clear();
    for( NodeT * n : next )
     open.push( n , n->g , w.key( n->g , n->h ) );
    for( NodeT * n : v_incons )
     open.push( n , n->g , w.key( n->g , n->h ) );
    v_incons.clear();
    }
   }

  // the Solution
  SolutionT finish( void ) {
   SolutionT sol;
   if( limit_hit )
    sol.status = Status::node_limit;
   else
    if( goal )
     sol.status = Status::solved;
   if( goal ) {
    sol.bound = ( optimal && ! limit_hit ) ? z : bound;
    sol.path = trace_path( goal );
    sol.cost = std::min( z , Cost( path_cost( d , sol.path ) ) );
    }
   sol.improvements = std::move( imp );
   settle_improvements( sol );
   st.stored = pool.size();
   st.expanded_by = { st.expanded };
   st.stored_by = { st.stored };
   st.seconds = clock.seconds();
   sol.stats = std::move( st );
   return( sol );
   }

  const AnytimeAStar & f_a;                     // the search
  const D & d;                                  // the Domain
  Stopwatch clock;                              // started at the start
  NodePool< NodeT > pool;                       // the nodes
  NodeTable< NodeT > table;                     // CLOSED
  Open open;                                    // OPEN, by f'
  Open fq;                                      // the same nodes by f
  std::vector< Kid > kids;                      // the successors of a node
  std::vector< NodeT * > v_closed;              // closed in the iteration
  std::vector< NodeT * > v_incons;              // INCONS
  std::vector< typename SolutionT::Improvement > imp;  // the paths found
  SearchStats st;                               // the counters
  NodeT * root = nullptr;                       // the node of s0
  NodeT * goal = nullptr;                       // the best goal
  Cost z = infinite_cost< Cost >();             // its cost
  Cost bound{};                                 // the last lower bound
  bool optimal = false;                         // z is proven optimal
  bool limit_hit = false;                       // the node limit is hit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 AnytimeOptions f_options;                           ///< the options

 };  // end( class( AnytimeAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with AWA* (empty schedule) or ARA*, choosing the OPEN
/// list at run time as astar() does

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > anytime_astar(
 const D & domain , const AnytimeOptions & options = {} ,
 OpenKind open = OpenKind::automatic )
{
 using NodeT = Node< state_t< D > , cost_t< D > >;
 if( open == OpenKind::heap_low_g )
  return( AnytimeAStar< D , LowGOpenList< HeapOpenList< NodeT > > >(
           domain , options ).solve() );
 if constexpr( std::is_integral_v< cost_t< D > > )
  if( open != OpenKind::heap )
   return( AnytimeAStar< D , BucketOpenList< NodeT > >( domain , options )
           .solve() );
 return( AnytimeAStar< D , HeapOpenList< NodeT > >( domain , options )
         .solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AnytimeAStar */

/*--------------------------------------------------------------------------*/
/*------------------------ End File AnytimeAStar.h -------------------------*/
/*--------------------------------------------------------------------------*/
