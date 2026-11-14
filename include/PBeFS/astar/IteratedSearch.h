/*--------------------------------------------------------------------------*/
/*------------------------- File IteratedSearch.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template IteratedSearch, a sequence of greedy
 * or weighted searches, each one bounded by the cost of the best path that
 * those before it have found, as the planners LAMA [see Richter and
 * Westphal, JAIR 2010] and Fast Downward run them, and as Restarting
 * Weighted A* restarts its search with smaller and smaller weights [see
 * Richter, Thayer and Ruml, ICAPS 2010].
 *
 * Each phase is a GreedySearch with its own GreedyOptions, whose bound is
 * the cost \f$ c \f$ of the incumbent: the successors with
 * \f$ g \geq c \f$ (or \f$ g + h \geq c \f$ for an admissible h) are
 * pruned, so that a phase that finds a path finds a cheaper one. The
 * phases may share the values of the heuristics [see HeuristicCache], so
 * that a state is evaluated once over the whole sequence, which is what
 * makes the restarts of Richter, Thayer and Ruml cheap.
 *
 * \code{.unparsed}
 *  1  procedure ITERATED(phases)
 *  2     c ← ∞;  best ← none
 *  3     for each phase P (the last one again while it improves, if so
 *          asked) do
 *  4        π ← the path found by P with the bound c, if any
 *  5        if π exists then c ← cost(π);  best ← π
 *  6        else if a path is known and failing ends the search then stop
 *  7     return best
 * \endcode
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

#ifndef __PBeFS_IteratedSearch
 #define __PBeFS_IteratedSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <limits>

#include <stdexcept>

#include <vector>

#include "PBeFS/astar/GreedySearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS IteratedOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of IteratedSearch

struct IteratedOptions {

 /// the phases, in order; the node limit of each is its own
 std::vector< GreedyOptions > phases;

 /// each phase is bounded by the cost of the incumbent
 bool pass_bound = true;

 /// the bound prunes by \f$ g + h \f$ (for an admissible h) instead of g
 bool prune_by_f = false;

 /// the last phase is repeated as long as it finds a cheaper path
 bool repeat_last = false;

 /// a phase that finds no path does not end the search, once a path is
 /// known (before, the search always goes on)
 bool continue_on_fail = false;

 /// the phases share the values of the heuristics [see HeuristicCache]
 bool share_cache = true;

 /// the most phases run, the repetitions included (0: no limit)
 std::size_t max_phases = 0;

 };  // end( struct( IteratedOptions ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS IteratedSearch --------------------------*/
/*--------------------------------------------------------------------------*/
/// a sequence of greedy or weighted searches, each bounded by the best path
/// found before it
/** solve() returns the best path, with the counters summed over the
 * phases, and the status of a phase that hits its node limit if no path
 * was found; phases() gives the cost and the counters of each phase run
 * (a cost of infinite_cost() for a phase that found no path). */

template< Domain D >
class IteratedSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// what a phase did
 struct Phase {
  Cost cost;                    ///< the cost of its path, or infinite
  SearchStats stats;            ///< its counters
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search); throws
 /// std::invalid_argument if there is no phase

 IteratedSearch( const D & domain , const IteratedOptions & options )
  : f_domain( domain ) , f_options( options ) {
  if( f_options.phases.empty() )
   throw( std::invalid_argument( "IteratedSearch: no phase" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// runs the phases, returning the best Solution

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  v_phases.clear();
  SolutionT best;
  SearchStats total;
  HeuristicCache< D > cache( f_domain );
  double bound = std::numeric_limits< double >::infinity();
  const std::size_t np = f_options.phases.size();
  bool limit = false;
  for( std::size_t i = 0 ; ; ) {
   if( f_options.max_phases && ( v_phases.size() >= f_options.max_phases ) )
    break;
   GreedyOptions o = f_options.phases[ i ];
   if( f_options.pass_bound ) {
    o.bound = bound;
    o.prune_by_f = f_options.prune_by_f;
    }
   GreedySearch< D > s( f_domain , o );
   if( f_options.share_cache )
    s.set_cache( &cache );
   SolutionT x = s.solve();
   add( total , x.stats );
   v_phases.push_back( Phase{ x.solved() ? x.cost : infinite_cost< Cost >() ,
                              x.stats } );
   limit = limit || ( x.status == Status::node_limit );
   const bool better = x.solved() &&
                       ( ( ! best.solved() ) || ( x.cost < best.cost ) );
   if( better ) {
    bound = double( x.cost );
    best = std::move( x );
    }
   else
    if( best.solved() && ! f_options.continue_on_fail )
     break;
   if( i + 1 < np )
    ++i;
   else
    if( ! ( f_options.repeat_last && better ) )
     break;
   }
  if( ( ! best.solved() ) && limit )
   best.status = Status::node_limit;
  best.stats = total;
  if( f_options.share_cache )
   best.stats.evaluated = cache.misses();
  best.stats.expanded_by = { best.stats.expanded };
  best.stats.stored_by = { best.stats.stored };
  best.stats.seconds = clock.seconds();
  return( best );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the phases run by the last solve(), in order

 [[nodiscard]] const std::vector< Phase > & phases( void ) const noexcept {
  return( v_phases );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// adds the counters of s to t
 static void add( SearchStats & t , const SearchStats & s ) {
  t.expanded += s.expanded;
  t.reexpanded += s.reexpanded;
  t.generated += s.generated;
  t.duplicates += s.duplicates;
  t.stored += s.stored;
  t.evaluated += s.evaluated;
  }

 const D & f_domain;                                 ///< the Domain
 IteratedOptions f_options;                          ///< the options
 std::vector< Phase > v_phases;                      ///< the phases run

 };  // end( class( IteratedSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the phases of the planner LAMA on a Domain with m heuristics
/** A lazy GBFS, then lazy searches ordered by \f$ g + w h \f$ with
 * \f$ w = 5 , 3 , 2 , 1 \f$, each one reopening the closed nodes, the last
 * one repeated while it improves; every phase alternates a queue on each
 * heuristic and, if preferred, one of the preferred successors on each
 * heuristic, boosted by 1000 [see Richter and Westphal, JAIR 2010]. */

[[nodiscard]] inline IteratedOptions lama_phases( std::size_t m ,
                                                  bool preferred )
{
 GreedyOptions g;
 g.evaluation = GreedyEvaluation::lazy;
 g.queues.clear();
 for( std::size_t i = 0 ; i < m ; ++i ) {
  g.queues.push_back( QueueSpec{ QueueKind::best , { i } , false } );
  if( preferred )
   g.queues.push_back( QueueSpec{ QueueKind::best , { i } , true } );
  }
 g.boost = preferred ? 1000 : 0;
 g.preferred_first = preferred;
 IteratedOptions it;
 it.phases.push_back( g );
 for( double w : { 5.0 , 3.0 , 2.0 , 1.0 } ) {
  GreedyOptions x = g;
  x.weight = w;
  x.reopen = true;
  it.phases.push_back( x );
  }
 it.repeat_last = true;
 it.continue_on_fail = true;
 return( it );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_IteratedSearch */

/*--------------------------------------------------------------------------*/
/*----------------------- End File IteratedSearch.h ------------------------*/
/*--------------------------------------------------------------------------*/
