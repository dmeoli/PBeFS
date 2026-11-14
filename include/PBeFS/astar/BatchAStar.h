/*--------------------------------------------------------------------------*/
/*--------------------------- File BatchAStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template BatchAStar, the batch weighted A*:
 * A* that takes the \f$ K \f$ best nodes of OPEN at each cycle, expands
 * them all, and computes the heuristic of their new successors in one
 * batch, with the f of a node \f$ g + w h \f$.
 *
 * With the notation of the page "Best-first search and A*" of the
 * documentation, \f$ K \ge 1 \f$ and \f$ w \ge 1 \f$:
 *
 * \code{.unparsed}
 *  1  procedure BATCHA*(s₀, K, w)
 *  2     CLOSED ← {s₀};  OPEN ← {(g(s₀) + w h(s₀), s₀)};  z ← ∞
 *  3     while OPEN has a node of f < z do
 *  4        B ← the (at most) K first nodes of OPEN of f < z, removed
 *  5        N ← ∅
 *  6        for each u ∈ B, closed, do
 *  7           if u is a goal then z ← min(z, g(u)), recording u;  continue
 *  8           for each successor v of u, g' ← g(u) + c(u, v) do
 *  9              if v ∉ CLOSED then add v to CLOSED and to N
 * 10              else if g' < g(v) then v is not closed
 * 11              else continue
 * 12              g(v) ← g';  parent(v) ← u
 * 13        compute h(v) for all the v ∈ N at once
 * 14        push in OPEN the nodes of lines 9 and 10 of f < z
 * 15     return the path to the goal of cost z, if z < ∞
 * \endcode
 *
 * With \f$ K = 1 \f$ the search is A* with the heuristic \f$ w h \f$: the
 * first goal selected sets \f$ z \f$ to its \f$ f \f$, which no node of
 * OPEN is below, and the loop ends. With \f$ K > 1 \f$ a goal may be
 * selected together with nodes of smaller f, whose expansions can still
 * lead to a cheaper goal, which is why the goal does not stop the search
 * at once; at the end every node of \f$ f < z \f$ has been expanded, so
 * that with \f$ w = 1 \f$ and \f$ h \f$ admissible \f$ z \f$ is optimal
 * whatever \f$ K \f$, and with \f$ w > 1 \f$ it is at most \f$ w \f$
 * times the optimum, since a node of the optimal path either is expanded or
 * has \f$ g + w h \le w C^* \f$ not below \f$ z \f$. A larger \f$ K \f$
 * expands more nodes, those of the last batch beyond the optimum, but calls
 * the heuristic \f$ K \f$ times less often.
 * This is the batch weighted A* by which Agostinelli et al. (2019) search
 * with a learned heuristic, written there with \f$ f = \lambda g + h \f$,
 * i.e., \f$ w = 1 / \lambda \f$.
 *
 * Line 13 calls heuristic_batch() of a BatchHeuristicDomain once per
 * cycle, on the new states only (a duplicate keeps the h it has); on
 * another Domain, the h given by for_each_successor is used, and the search
 * differs from A* only by the order of the expansions.
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

#ifndef __PBeFS_BatchAStar
 #define __PBeFS_BatchAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include <cstdint>

#include <span>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/AStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BatchOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of BatchAStar, whose weight \f$ w \f$ of the heuristic is
/// that of SearchOptions (1 if smaller)

struct BatchOptions : public SearchOptions {

 /// the nodes expanded at each cycle, 1 if 0
 std::size_t batch = 1;

 };  // end( struct( BatchOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BatchAStar ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the batch weighted A* over a Domain, with the OPEN list Open
/** A node keeps \f$ w h \f$ in its field h, rounded down with integer
 * costs, so that the OPEN list orders by \f$ g + w h \f$; the Solution
 * counts in generated the successors, and the heuristic is computed for
 * stored of them, one per node. */

template< Domain D ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class BatchAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the BatchAStar)

 explicit BatchAStar( const D & domain , const BatchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  if( f_options.batch == 0 )
   f_options.batch = 1;
  if( ! ( f_options.weight >= 1 ) )
   f_options.weight = 1;
  }

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
  Open open;
  std::vector< NodeT * > batch , fresh , push;
  std::vector< State > states;
  std::vector< Cost > hs;

  const State s0 = f_domain.initial_state();
  NodeT * root = pool.allocate();
  *root = NodeT{ s0 , Cost( 0 ) , inflate( heuristic( s0 , states , hs ) ) ,
                 nullptr , false , false };
  table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   return( root );
   } );
  open.push( root );

  NodeT * goal = nullptr;
  Cost z = infinite_cost< Cost >();
  for( ;; ) {
   batch.clear();                                       // line 4
   while( batch.size() < f_options.batch ) {
    while( ! open.empty() && open.top().stale() )
     open.pop();
    if( open.empty() || ( open.top_f() >= z ) )
     break;
    NodeT * n = open.pop().node;
    n->closed = true;
    batch.push_back( n );
    }
   if( batch.empty() )                                  // line 3
    break;

   fresh.clear();                                       // lines 5 to 12
   push.clear();
   for( NodeT * n : batch ) {
    ++st.expanded;
    if( n->expanded )
     ++st.reexpanded;
    n->expanded = true;
    if( f_domain.is_goal( n->state ) ) {
     if( n->g < z ) {
      z = n->g;
      goal = n;
      }
     continue;
     }
    // the Domain computes the h of the successors from the true one
    const Cost hn = ( f_options.weight == 1 ) || BatchHeuristicDomain< D >
                    ? n->h : Cost( f_domain.heuristic( n->state ) );
    for_each_child( f_domain , n->state , hn , n->op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
     ++st.generated;
     const Cost g = n->g + c;
     auto [ m , isnew ] = table.find_or_insert( t , f_domain.hash( t ) ,
                                                [ & ] {
      NodeT * x = pool.allocate();
      *x = NodeT{ t , g , ht , n , false , false , op };
      return( x );
      } );
     if( isnew ) {
      fresh.push_back( m );
      return;
      }
     if( g >= m->g ) {
      ++st.duplicates;
      return;
      }
     m->g = g;
     m->parent = n;
     m->op = op;
     m->closed = false;
     push.push_back( m );
     } );
    }

   heuristics( fresh , states , hs );                   // line 13
   for( NodeT * m : fresh )                             // line 14
    if( m->g + m->h < z )
     open.push( m );
   for( NodeT * m : push )
    if( m->g + m->h < z )
     open.push( m );

   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }

  if( ( sol.status != Status::node_limit ) && goal ) {
   sol.status = Status::solved;
   sol.cost = goal->g;
   sol.path = trace_path( goal );
   }

  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// \f$ w h \f$, rounded down with integer costs
 [[nodiscard]] Cost inflate( Cost h ) const {
  if( f_options.weight == 1 )
   return( h );
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( std::floor( f_options.weight * double( h ) ) ) );
  else
   return( Cost( f_options.weight * h ) );
  }

 /// the heuristic of s, through heuristic_batch() if there is one
 Cost heuristic( const State & s , std::vector< State > & states ,
                 std::vector< Cost > & hs ) const {
  if constexpr( BatchHeuristicDomain< D > ) {
   states.assign( 1 , s );
   hs.assign( 1 , Cost( 0 ) );
   f_domain.heuristic_batch( std::span< const State >( states ) ,
                             std::span< Cost >( hs ) );
   return( hs[ 0 ] );
   }
  else
   return( Cost( f_domain.heuristic( s ) ) );
  }

 /// the field h of the nodes of v, \f$ w h \f$, computed in one batch for
 /// a BatchHeuristicDomain, from the h given with the successors otherwise
 void heuristics( const std::vector< NodeT * > & v ,
                  std::vector< State > & states ,
                  std::vector< Cost > & hs ) const {
  if constexpr( BatchHeuristicDomain< D > ) {
   if( v.empty() )
    return;
   states.clear();
   for( const NodeT * m : v )
    states.push_back( m->state );
   hs.assign( v.size() , Cost( 0 ) );
   f_domain.heuristic_batch( std::span< const State >( states ) ,
                             std::span< Cost >( hs ) );
   for( std::size_t i = 0 ; i < v.size() ; ++i )
    v[ i ]->h = inflate( hs[ i ] );
   }
  else
   for( NodeT * m : v )
    m->h = inflate( m->h );
  }

 const D & f_domain;                                 ///< the Domain
 BatchOptions f_options;                             ///< the options

 };  // end( class( BatchAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with the batch weighted A*

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > batch_astar(
 const D & domain , const BatchOptions & options = {} )
{
 return( BatchAStar< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BatchAStar */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BatchAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
