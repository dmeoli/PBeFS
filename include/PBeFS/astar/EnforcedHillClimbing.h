/*--------------------------------------------------------------------------*/
/*---------------------- File EnforcedHillClimbing.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template EnforcedHillClimbing, the local search
 * of the planner FF [see Hoffmann and Nebel, JAIR 2001].
 *
 * From the current state \f$ s \f$, a breadth-first search looks for a
 * state of heuristic value smaller than \f$ h( s ) \f$, which becomes the
 * current state, the path to it being appended to the plan; the search
 * ends at a goal, and fails if a breadth-first search exhausts the states
 * reachable from \f$ s \f$ without finding a better one:
 *
 * \code{.unparsed}
 *  1  procedure EHC(s₀)
 *  2     s ← s₀;  plan ← ⟨s₀⟩
 *  3     while s is not a goal do
 *  4        Q ← ⟨s⟩;  visited ← {s}
 *  5        loop
 *  6           if Q = ∅ then return failure      ▷ or a complete search
 *  7           u ← the first state of Q, removed
 *  8           for each successor v ∉ visited of u (the preferred ones
 *                only, or first, if so asked) do
 *  9              visited ← visited ∪ {v};  parent(v) ← u
 * 10              if v is a goal or h(v) < h(s) then
 * 11                 append the path from s to v to plan;  s ← v
 * 12                 continue at line 3
 * 13              append v to Q
 * 14     return plan
 * \endcode
 *
 * The values of h are kept over the whole search, so that a state is
 * evaluated once even if several breadth-first searches meet it. On a
 * PreferredDomain the successors may be restricted to the preferred ones
 * (the helpful actions of FF), or ordered with them first; FF falls back on
 * a complete best-first search when the search with the helpful actions
 * fails, which the option fallback reproduces with GreedySearch.
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

#ifndef __PBeFS_EnforcedHillClimbing
 #define __PBeFS_EnforcedHillClimbing
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <vector>

#include "PBeFS/astar/GreedySearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// how enforced hill-climbing uses the preferred successors

enum class PreferredUsage {
 ignore ,     ///< as the others
 prune ,      ///< only the preferred successors are generated
 first        ///< the preferred successors are visited first
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the options of EnforcedHillClimbing

struct EHCOptions : public SearchOptions {

 /// the use of the preferred successors [see PreferredDomain]
 PreferredUsage preferred = PreferredUsage::ignore;

 /// if the search fails, a GreedySearch from the initial state with the
 /// options greedy, as FF does
 bool fallback = false;

 /// the options of the search of the fallback
 GreedyOptions greedy;

 };  // end( struct( EHCOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS EnforcedHillClimbing -----------------------*/
/*--------------------------------------------------------------------------*/
/// enforced hill-climbing, the local search of FF
/** The search of the file comment. The counters are those of SearchStats,
 * summed over the breadth-first searches (and the fallback, if it runs);
 * evaluated counts each state once, since the values are kept, and
 * expanded counts a state once per breadth-first search that expands it.
 * The node limit bounds the states evaluated. */

template< Domain D >
class EnforcedHillClimbing {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit EnforcedHillClimbing( const D & domain ,
                                const EHCOptions & options = {} )
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
  std::vector< std::uint32_t > stamp;     // the phase that visited a node
  std::deque< NodeT * > q;
  std::vector< Kid > kids;

  // the node of t, evaluated with ht if new
  const auto node = [ & ]( const State & t , Cost ht , bool & fresh ) {
   fresh = false;
   NodeT * x = table.find_or_insert( t , f_domain.hash( t ) , [ & ] {
    fresh = true;
    NodeT * y = pool.allocate();
    *y = NodeT{ t , Cost( 0 ) , ht , nullptr , false , false , no_op ,
                std::uint32_t( pool.size() - 1 ) };
    stamp.push_back( 0 );
    return( y );
    } ).first;
   if( fresh )
    ++st.evaluated;
   return( x );
   };

  const State s0 = f_domain.initial_state();
  bool fresh;
  NodeT * s = node( s0 , Cost( f_domain.heuristic( s0 ) ) , fresh );
  std::vector< State > plan{ s0 };
  Cost cost( 0 );
  std::uint32_t phase = 0;
  bool failed = s->h == infinite_cost< Cost >();
  while( ( ! failed ) && ! f_domain.is_goal( s->state ) ) {
   ++phase;
   q.clear();
   q.push_back( s );
   stamp[ s->index ] = phase;
   s->g = Cost( 0 );
   s->parent = nullptr;
   NodeT * better = nullptr;
   while( ( ! better ) && ! q.empty() ) {
    NodeT * u = q.front();
    q.pop_front();
    ++st.expanded;
    kids.clear();
    for_each_child( f_domain , u->state , u->h , u->op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
     kids.push_back( Kid{ t , c , ht , op , false } );
     } );
    st.generated += kids.size();
    if constexpr( PreferredDomain< D > )
     if( f_options.preferred != PreferredUsage::ignore ) {
      const auto p = f_domain.preferences( u->state );
      for( Kid & k : kids )
       k.preferred = bool( p( k.s ) );
      if( f_options.preferred == PreferredUsage::prune )
       std::erase_if( kids , []( const Kid & k ) { return( ! k.preferred ); }
                      );
      else
       std::stable_partition( kids.begin() , kids.end() ,
                              []( const Kid & k ) { return( k.preferred ); }
                              );
      }
    for( const Kid & k : kids ) {
     NodeT * v = node( k.s , k.h , fresh );
     if( stamp[ v->index ] == phase ) {
      ++st.duplicates;
      continue;
      }
     stamp[ v->index ] = phase;
     v->g = u->g + k.c;
     v->parent = u;
     v->op = k.op;
     if( v->h == infinite_cost< Cost >() )
      continue;                                 // a dead end
     if( f_domain.is_goal( v->state ) || ( v->h < s->h ) ) {
      better = v;
      break;
      }
     q.push_back( v );
     }
    if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) )
     break;
    }
   if( ! better ) {
    if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) )
     sol.status = Status::node_limit;
    failed = true;
    break;
    }
   // the path from s to better, appended to the plan
   std::vector< State > piece;
   for( const NodeT * x = better ; x != s ; x = x->parent )
    piece.push_back( x->state );
   plan.insert( plan.end() , piece.rbegin() , piece.rend() );
   cost += better->g;
   s = better;
   s->op = better->op;
   }

  if( ! failed ) {
   sol.status = Status::solved;
   sol.cost = cost;
   sol.path = std::move( plan );
   }
  else
   if( f_options.fallback && ( sol.status != Status::node_limit ) ) {
    GreedyOptions g = f_options.greedy;
    g.node_limit = f_options.node_limit;
    g.expected_nodes = f_options.expected_nodes;
    SolutionT b = GreedySearch< D >( f_domain , g ).solve();
    b.stats.expanded += st.expanded;
    b.stats.generated += st.generated;
    b.stats.duplicates += st.duplicates;
    b.stats.evaluated += st.evaluated;
    b.stats.stored += pool.size();
    sol = std::move( b );
    sol.stats.expanded_by = { sol.stats.expanded };
    sol.stats.stored_by = { sol.stats.stored };
    sol.stats.seconds = clock.seconds();
    return( sol );
    }
  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using NodeT = GreedyNode< State , Cost >;         ///< the node type

 /// a successor
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint8_t op;
  bool preferred;
  };

 const D & f_domain;                                 ///< the Domain
 EHCOptions f_options;                               ///< the options

 };  // end( class( EnforcedHillClimbing ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with EnforcedHillClimbing and the given options

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > ehc(
 const D & domain , const EHCOptions & options = {} )
{
 return( EnforcedHillClimbing< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_EnforcedHillClimbing */

/*--------------------------------------------------------------------------*/
/*-------------------- End File EnforcedHillClimbing.h ---------------------*/
/*--------------------------------------------------------------------------*/
