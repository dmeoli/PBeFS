/*--------------------------------------------------------------------------*/
/*---------------------- File RacingBranchAndBound.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the racing ramp-up of the Branch-and-Bound (B&B): several
 * settings search the tree from the root at the same time, each with a
 * budget of nodes, sharing their incumbents, and the search goes on, by any
 * of the parallel schemes, from the open nodes of the one whose bound is the
 * best, i.e., the tree of the winner is kept, as in the racing ramp-up of
 * the Ubiquity Generator of Shinano et al. (2016). The problem, the options
 * and the Solution are those of BranchAndBound.h and
 * ParallelBranchAndBound.h, and the design is discussed in
 * docs/branch-and-bound.md.
 *
 * The ramp-up is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation and the procedures of the
 * sequential search [see BranchAndBound.h]; \f$ s_k \f$ is the setting
 * \f$ k \f$ of BasicParallelBBOptions::race and \f$ N \f$ the budget of
 * nodes of a racer:
 *
 * \code{.unparsed}
 *  1  procedure RACE(R)
 *  2     in parallel, for each setting k do
 *  3        S_k ← B&B(R) with the selection of s_k, the evaluator of the
 *              worker k, at most N evaluations, pruning by the best
 *              incumbent of all the racers, and keeping its open leaves
 *  4     if some S_k ended without a limit then return the first of them
 *  5     w ← the k of the best bound of S_k, the fewest open leaves on ties
 *  6     z ← the best value of all the S_k
 *  7     S ← the scheme from the open leaves of S_w, with the incumbent z
 *  8     return S, its bound lowered by the closed leaves of S_w
 * \endcode
 *
 * The search is correct since the open and the closed leaves of S_w cover
 * R, and each closed one was pruned by an incumbent not better than z; the
 * racers that lose give only their incumbents. Racing pays when the
 * settings differ much in how early they find a good incumbent or how fast
 * they raise the bound, which is typical of the first nodes of a B&B, when
 * there are fewer open nodes than workers; the counters of the Solution
 * include the evaluations of all the racers. In the mode enumerate there is
 * nothing to race, and the scheme starts from R at once.
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

#ifndef __PBeFS_RacingBranchAndBound
 #define __PBeFS_RacingBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <exception>

#include <mutex>

#include <optional>

#include <utility>

#include <variant>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the B&B of the problem started by a race of the settings of o.race,
/// which goes on by search( problem , options , starts ), e.g., a call of
/// parallel_branch_and_bound() with a scheme, from the open leaves of the
/// winner, with the incumbent of all the racers in the options
/** The racer k runs BranchAndBoundSearch with the evaluator of the worker k
 * and the options o, but for the selection, the dive and prune_level of
 * o.race[ k ] and at most o.race_nodes evaluations (100 times the workers
 * if 0); search gets o without the race, with the incumbent of the race and
 * what is left of the node limit. The racers run at once on an executor of
 * as many threads as there are settings. */

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Search >
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
 racing_branch_and_bound( const P & problem ,
                          const BasicParallelBBOptions< bb_value_t< P > > & o ,
                          Search && search ,
                          std::optional< std::vector< BBLeaf<
                           typename P::node_type , bb_value_t< P > > > >
                           seeds = std::nullopt )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 using Tr = BBValueTraits< V >;
 using Leaf = BBLeaf< Node , V >;
 using SolutionT = BBSolution< Node , V >;

 BasicParallelBBOptions< V > rest = o;
 rest.race.clear();
 if( o.race.empty() || ( o.mode == BBMode::enumerate ) ) {
  std::vector< Leaf > from;
  if( seeds )
   from = std::move( *seeds );
  else
   from.push_back( Leaf{ problem.root() ,
                         detail::bb_signed< P >( - Tr::infinity() ) , 0 ,
                         BBLeafKind::open } );
  return( search( problem , rest , std::move( from ) ) );
  }
 Stopwatch clock;
 const unsigned k = unsigned( o.race.size() );
 const std::size_t budget = o.race_nodes ? o.race_nodes :
                            100 * std::size_t( o.threads ? o.threads : 1 );

 // the racers, pruning by the best incumbent of all of them
 BBIncumbent< std::monostate , V > best( o.incumbent ?
                                         detail::bb_signed< P >(
                                          *o.incumbent ) :
                                         Tr::infinity() );
 std::vector< std::optional< SolutionT > > out( k );
 std::exception_ptr error;
 std::mutex m;
 run_workers< Executor >( nullptr , k , [ & ]( unsigned i ) {
  try {
   BasicBBOptions< V > ro = o;
   ro.selection = o.race[ i ].selection;
   ro.plunge_gap = o.race[ i ].plunge_gap;
   ro.prune_level = o.race[ i ].prune_level;
   ro.node_limit = budget;
   if( o.node_limit )
    ro.node_limit = std::min( budget , std::size_t( o.node_limit ) );
   ro.keep_frontier = true;
   ro.keep_open_only = ! ( o.keep_frontier && ! o.keep_open_only );
   ro.on_incumbent = [ & ]( V v ) {
    ( void ) best.template improve< P >( detail::bb_signed< P >( v ) ,
                                         std::monostate() , o );
    };
   ro.shared_incumbent = [ & ] {
    return( detail::bb_signed< P >( best.value() ) );
    };
   std::vector< Leaf > from;
   if( seeds )
    from = *seeds;
   else
    from.push_back( Leaf{ problem.root() ,
                          detail::bb_signed< P >( - Tr::infinity() ) , 0 ,
                          BBLeafKind::open } );
   out[ i ] = BranchAndBoundSearch< P >( problem , ro , i ).solve(
                                                           std::move( from ) );
   }
  catch( ... ) {
   std::lock_guard< std::mutex > guard( m );
   if( ! error )
    error = std::current_exception();
   }
  } );
 if( error )
  std::rethrow_exception( error );

 // the winner: a racer that ended, or the best bound, the fewest open
 // leaves on ties; the incumbent: the best of all
 BBStats st;
 std::size_t w = k;
 for( std::size_t i = 0 ; i < k ; ++i ) {
  const SolutionT & s = *out[ i ];
  detail::bb_add_counters( st , s.stats );
  if( ( s.status != BBStatus::node_limit ) &&
      ( s.status != BBStatus::time_limit ) ) {
   if( ( w == k ) || ( out[ w ]->status == BBStatus::node_limit ) )
    w = i;
   continue;
   }
  if( w == k )
   w = i;
  else
   if( out[ w ]->status == BBStatus::node_limit ) {
    const V bi = detail::bb_signed< P >( s.bound );
    const V bw = detail::bb_signed< P >( out[ w ]->bound );
    if( ( bi > bw ) || ( ( bi == bw ) && ( s.frontier.size() <
                                           out[ w ]->frontier.size() ) ) )
     w = i;
    }
  }
 V z = o.incumbent ? detail::bb_signed< P >( *o.incumbent ) : Tr::infinity();
 std::optional< Node > zn;
 for( std::size_t i = 0 ; i < k ; ++i )
  if( out[ i ]->node && ( detail::bb_signed< P >( out[ i ]->value ) < z ) ) {
   z = detail::bb_signed< P >( out[ i ]->value );
   zn = out[ i ]->node;
   }

 // the scheme, from the open leaves of the winner, unless it ended; its
 // closed leaves are leaves of the search
 SolutionT & win = *out[ w ];
 std::vector< Leaf > open , closed;
 for( Leaf & l : win.frontier )
  ( l.kind == BBLeafKind::open ? open : closed ).push_back( std::move( l ) );
 const V wclosed = detail::bb_signed< P >( win.closed_bound );
 if( Tr::finite( z ) )
  rest.incumbent = detail::bb_signed< P >( z );
 if( o.node_limit )
  rest.node_limit = st.evaluated >= o.node_limit ? 1 :
                    o.node_limit - std::size_t( st.evaluated );
 if( o.time_limit > 0 )
  rest.time_limit = std::max( 1e-9 , o.time_limit - clock.seconds() );
 SolutionT sol;
 if( open.empty() ) {
  sol.status = detail::bb_status< P >( z , o );
  if( Tr::finite( z ) )
   sol.value = detail::bb_signed< P >( z );
  sol.bound = detail::bb_signed< P >( z );
  sol.closed_bound = detail::bb_signed< P >( Tr::infinity() );
  }
 else
  sol = search( problem , rest , std::move( open ) );

 // the node of the incumbent of the race if the scheme found none better
 // (its value is then that of the race), and the bound lowered by the
 // closed leaves of the winner
 if( ! sol.node && zn )
  sol.node = std::move( zn );
 const V sz = sol.feasible() ? detail::bb_signed< P >( sol.value )
                             : Tr::infinity();
 const V lb = std::min( std::min( wclosed , sz ) ,
                        detail::bb_signed< P >( sol.bound ) );
 sol.bound = detail::bb_signed< P >( lb );
 sol.closed_bound = detail::bb_signed< P >( std::min(
                     detail::bb_signed< P >( sol.closed_bound ) , wclosed ) );
 sol.gap = Tr::finite( sz ) && Tr::finite( lb ) ? V( sz - lb )
                                               : Tr::infinity();
 if( o.keep_frontier )
  sol.frontier.insert( sol.frontier.end() ,
                       std::make_move_iterator( closed.begin() ) ,
                       std::make_move_iterator( closed.end() ) );
 detail::bb_add_counters( st , sol.stats );
 st.evaluated_by = std::move( sol.stats.evaluated_by );
 st.max_open = std::max( st.max_open , sol.stats.max_open );
 sol.stats = std::move( st );
 sol.stats.seconds = clock.seconds();
 return( sol );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_RacingBranchAndBound */

/*--------------------------------------------------------------------------*/
/*-------------------- End File RacingBranchAndBound.h ---------------------*/
/*--------------------------------------------------------------------------*/
