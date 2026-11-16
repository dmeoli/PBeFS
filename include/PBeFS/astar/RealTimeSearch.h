/*--------------------------------------------------------------------------*/
/*------------------------- File RealTimeSearch.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template RealTimeSearch, the agent-centered
 * real-time searches LRTA* [see Korf, AIJ 1990] and LSS-LRTA* [see Koenig
 * and Sun, JAAMAS 2009], which interleave a search of bounded size with
 * the moves of an agent, and learn a heuristic H on the way.
 *
 * From the state s of the agent, an A* ordered by \f$ g + H \f$ expands at
 * most L states (the lookahead), stopping earlier if a goal is first in
 * OPEN; then the values H of the states it expanded are raised by a
 * Dijkstra search backward from OPEN, i.e., to the least
 * \f$ c( u , v ) + H( v ) \f$ over the paths to OPEN within the local
 * space, and the agent moves along the path to the state of least
 * \f$ g + H \f$ in OPEN:
 *
 * \code{.unparsed}
 *  1  procedure LSS-LRTA*(s₀, L)
 *  2     s ← s₀;  H(u) ← h(u) for each state u not yet met
 *  3     while s is not a goal do
 *  4        A* from s on g + H, at most L expansions, stopped when a goal
 *             is the first of OPEN;  CLOSED_s and OPEN_s its lists
 *  5        H(u) ← ∞ for each u ∈ CLOSED_s
 *  6        while some u ∈ CLOSED_s is not yet taken do     ▷ Dijkstra
 *  7           v ← the state of least H in OPEN_s ∪ CLOSED_s not yet taken,
 *                 whose H is then final
 *  8           for each u ∈ CLOSED_s with v a successor of u do
 *  9              H(u) ← min(H(u), c(u, v) + H(v))
 * 10        move the agent along the path to the first state of OPEN_s
 * \endcode
 *
 * With L = 1 this is LRTA*, whose update is
 * \f$ H( s ) = \min_{t} c( s , t ) + H( t ) \f$ and whose move is to the
 * best successor. With an admissible h, H stays admissible (Koenig and
 * Sun, theorem 1); on a finite graph where a goal is reachable from every
 * state, each trial reaches a goal, and trials repeated with the H
 * learned converge to an optimal path (theorem 3 there, and Korf for
 * LRTA*). solve() runs the trials of RealTimeOptions, returning the path
 * of the last one, which may visit a state more than once; trials()
 * gives the cost of each.
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

#ifndef __PBeFS_RealTimeSearch
 #define __PBeFS_RealTimeSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <functional>

#include <queue>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS RealTimeOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of RealTimeSearch

//! [memory-bounded-sig-realtimeoptions]
struct RealTimeOptions : public SearchOptions {

 /// the expansions L of the local search (1: LRTA*)
 std::size_t lookahead = 1;

 /// the trials, each from the initial state with the H learned so far
 std::size_t trials = 1;

 /// the trials go on, up to trials, until one leaves H unchanged
 bool until_convergence = false;

 /// the most moves of a trial (0: no limit)
 std::size_t max_moves = 0;

 };  // end( struct( RealTimeOptions ) )
//! [memory-bounded-sig-realtimeoptions]

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS RealTimeSearch --------------------------*/
/*--------------------------------------------------------------------------*/
/// LSS-LRTA*, and LRTA* for a lookahead of 1

template< Domain D >
class RealTimeSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit RealTimeSearch( const D & domain ,
                          const RealTimeOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_H( 1024 , Hasher{ &domain } ) {
  if( f_options.lookahead == 0 )
   f_options.lookahead = 1;
  if( f_options.trials == 0 )
   f_options.trials = 1;
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// runs the trials, returning the path of the last one

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  v_trials.clear();
  for( std::size_t k = 0 ; k < f_options.trials ; ++k ) {
   f_changed = false;
   SolutionT x = trial( st );
   v_trials.push_back( x.solved() ? x.cost : infinite_cost< Cost >() );
   sol.status = x.status;
   sol.cost = x.cost;
   sol.path = std::move( x.path );
   if( ! sol.solved() )
    break;
   if( f_options.until_convergence && ! f_changed )
    break;
   }
  st.stored = f_H.size();
  st.evaluated = f_H.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost of the path of each trial of the last solve()

 [[nodiscard]] const std::vector< Cost > & trials( void ) const noexcept {
  return( v_trials );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the heuristic learned for s, h( s ) if never updated

 [[nodiscard]] Cost learned( const State & s ) const {
  const auto it = f_H.find( s );
  return( it == f_H.end() ? Cost( f_domain.heuristic( s ) ) : it->second );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the hash of a state by the Domain
 struct Hasher {
  const D * d;
  std::size_t operator()( const State & s ) const {
   return( std::size_t( d->hash( s ) ) );
   }
  };

 /// a state of the local search
 struct Local {
  Cost g;
  State parent;
  bool has_parent;
  bool closed;
  std::vector< std::pair< State , Cost > > succ;   // if closed
  };

 using LocalMap = std::unordered_map< State , Local , Hasher >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// H of s, with h( s ) recorded if s is new

 Cost H( const State & s ) {
  auto it = f_H.find( s );
  if( it == f_H.end() )
   it = f_H.emplace( s , Cost( f_domain.heuristic( s ) ) ).first;
  return( it->second );
  }

 /// sets H of s to v
 void set_H( const State & s , Cost v ) { f_H[ s ] = v; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// one trial from the initial state

 SolutionT trial( SearchStats & st ) {
  SolutionT sol;
  State s = f_domain.initial_state();
  std::vector< State > path{ s };
  Cost cost( 0 );
  std::size_t moves = 0;
  using QItem = std::pair< Cost , std::uint64_t >;
  std::vector< State > by_id;
  while( ! f_domain.is_goal( s ) ) {
   if( f_options.max_moves && ( moves >= f_options.max_moves ) ) {
    sol.status = Status::node_limit;
    return( sol );
    }
   // the local A*, ordered by g + H and then by larger g, then FIFO
   LocalMap loc( 64 , Hasher{ &f_domain } );
   struct OItem {
    Cost f , g;
    std::uint64_t seq;
    State s;
    };
   const auto worse = []( const OItem & a , const OItem & b ) {
    if( a.f != b.f )
     return( a.f > b.f );
    if( a.g != b.g )
     return( a.g < b.g );
    return( a.seq > b.seq );
    };
   std::priority_queue< OItem , std::vector< OItem > , decltype( worse ) >
    open( worse );
   std::uint64_t seq = 0;
   loc.emplace( s , Local{ Cost( 0 ) , s , false , false , {} } );
   open.push( OItem{ H( s ) , Cost( 0 ) , seq++ , s } );
   std::size_t expanded = 0;
   while( ! open.empty() ) {
    const OItem e = open.top();
    Local & l = loc.at( e.s );
    if( l.closed || ( l.g != e.g ) ) {
     open.pop();
     continue;
     }
    if( f_domain.is_goal( e.s ) || ( expanded >= f_options.lookahead ) )
     break;
    open.pop();
    l.closed = true;
    ++expanded;
    ++st.expanded;
    std::vector< std::pair< State , Cost > > succ;
    f_domain.for_each_successor( e.s , Cost( 0 ) ,
                                 [ & ]( const State & t , Cost c , Cost ) {
     succ.emplace_back( t , c );
     } );
    st.generated += succ.size();
    for( const auto & [ t , c ] : succ ) {
     const Cost g = e.g + c;
     auto it = loc.find( t );
     if( it == loc.end() ) {
      loc.emplace( t , Local{ g , e.s , true , false , {} } );
      open.push( OItem{ g + H( t ) , g , seq++ , t } );
      }
     else
      if( ( ! it->second.closed ) && cheaper( g , it->second.g ) ) {
       it->second.g = g;
       it->second.parent = e.s;
       it->second.has_parent = true;
       open.push( OItem{ g + H( t ) , g , seq++ , t } );
       }
      else
       ++st.duplicates;
     }
    loc.at( e.s ).succ = std::move( succ );
    }
   if( open.empty() ) {                        // no way out: a dead end
    sol.status = Status::unsolvable;
    return( sol );
    }
   const State target = open.top().s;
   // the Dijkstra update of H on the closed states, from the open ones
   std::unordered_map< State , std::vector< std::pair< State , Cost > > ,
                       Hasher > pred( 64 , Hasher{ &f_domain } );
   std::size_t left = 0;
   std::vector< std::pair< State , Cost > > before;
   for( auto & [ u , l ] : loc )
    if( l.closed ) {
     ++left;
     before.emplace_back( u , H( u ) );
     set_H( u , infinite_cost< Cost >() );
     for( const auto & [ v , c ] : l.succ )
      pred[ v ].emplace_back( u , c );
     }
   by_id.clear();
   std::priority_queue< QItem , std::vector< QItem > ,
                        std::greater< QItem > > q;
   for( auto & [ u , l ] : loc )
    if( ! l.closed ) {
     by_id.push_back( u );
     q.push( QItem{ H( u ) , by_id.size() - 1 } );
     }
   while( left && ! q.empty() ) {
    const auto [ hv , id ] = q.top();
    q.pop();
    const State v = by_id[ id ];
    if( hv != H( v ) )
     continue;
    if( loc.at( v ).closed )
     --left;                                 // v is final
    const auto pit = pred.find( v );
    if( pit == pred.end() )
     continue;
    for( const auto & [ u , c ] : pit->second ) {
     const Cost cand = ( hv == infinite_cost< Cost >() )
                       ? hv : Cost( c + hv );
     const Cost hu = H( u );
     if( cand < hu ) {
      set_H( u , cand );
      by_id.push_back( u );
      q.push( QItem{ cand , by_id.size() - 1 } );
      }
     }
    }
   for( const auto & [ u , v ] : before )
    if( H( u ) != v )
     f_changed = true;
   // the agent moves along the path to target
   std::vector< State > piece;
   for( State x = target ; ; ) {
    const Local & l = loc.at( x );
    if( ! l.has_parent )
     break;
    piece.push_back( x );
    x = l.parent;
    }
   std::reverse( piece.begin() , piece.end() );
   State prev = s;
   for( const State & x : piece ) {
    Cost best = infinite_cost< Cost >();
    f_domain.for_each_successor( prev , Cost( 0 ) ,
                                 [ & ]( const State & t , Cost c , Cost ) {
     if( ( t == x ) && ( c < best ) )
      best = c;
     } );
    cost = cost + best;
    path.push_back( x );
    prev = x;
    ++moves;
    }
   s = target;
   }
  sol.status = Status::solved;
  sol.cost = cost;
  sol.path = std::move( path );
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 RealTimeOptions f_options;                          ///< the options
 std::unordered_map< State , Cost , Hasher > f_H;    ///< the H learned
 bool f_changed = false;                             ///< H changed in the
                                                     ///< trial
 std::vector< Cost > v_trials;                       ///< the trials

 };  // end( class( RealTimeSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with LRTA* (lookahead 1) or LSS-LRTA*

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > lss_lrta_star(
 const D & domain , const RealTimeOptions & options = {} )
{
 return( RealTimeSearch< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_RealTimeSearch */

/*--------------------------------------------------------------------------*/
/*----------------------- End File RealTimeSearch.h ------------------------*/
/*--------------------------------------------------------------------------*/
