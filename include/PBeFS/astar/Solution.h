/*--------------------------------------------------------------------------*/
/*---------------------------- File Solution.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the outcome of a search: the class template Solution, the
 * statistics SearchStats collected on the way, and the SearchOptions shared
 * by all the searches.
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

#ifndef __PBeFS_Solution
 #define __PBeFS_Solution
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// how a search ended
enum class Status {
 solved ,       ///< a path has been found, optimal or within the bound
                ///< of the weight [see SearchOptions]
 unsolvable ,   ///< OPEN emptied without reaching a goal
 node_limit ,   ///< the limit on the stored nodes was reached first;
                ///< an anytime search keeps the best path found
 disk_limit     ///< the limit on the bytes on the disk was reached first
 };

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SearchStats ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the counters of a search
/** The vectors have one entry per worker, a single one for a sequential
 * search. The field sent counts the nodes handed to another worker, so
 * that the communication overhead of Jinnai and Fukunaga (JAIR 2017) is
 * \f$ \mathrm{CO} = \mathrm{sent} / \mathrm{generated} \f$, and the load
 * balance is \f$ \mathrm{LB} = \max_i E_i / ( \frac{1}{p} \sum_i E_i ) \f$,
 * \f$ E_i \f$ being the expansions of the worker \f$ i \f$ in
 * expanded_by. */

struct SearchStats {

 std::uint64_t expanded = 0;      ///< expansions, re-expansions included
 std::uint64_t reexpanded = 0;    ///< expansions of a node already expanded
 std::uint64_t generated = 0;     ///< successors generated
 std::uint64_t duplicates = 0;    ///< generated states already known with
                                  ///< a \f$ g \f$ not larger
 std::uint64_t sent = 0;          ///< nodes sent to another worker
 std::uint64_t messages = 0;      ///< batches of nodes sent
 std::uint64_t stored = 0;        ///< distinct nodes stored at the end
 std::uint64_t evaluated = 0;     ///< heuristic evaluations, counted by
                                  ///< the searches that defer or cache
                                  ///< them (as many as the nodes stored
                                  ///< in the others)
 double seconds = 0;              ///< wall-clock time of the search
 std::vector< std::uint64_t > expanded_by;  ///< expansions of each worker
 std::vector< std::uint64_t > stored_by;    ///< nodes stored by each worker

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the load balance \f$ \mathrm{LB} \f$, the largest over the mean of
 /// the expansions per worker (1 is a perfect balance)

 [[nodiscard]] double load_balance( void ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the communication overhead \f$ \mathrm{CO} \f$, the nodes sent to
 /// another worker over the nodes generated

 [[nodiscard]] double communication_overhead( void ) const;

 };  // end( struct( SearchStats ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Solution -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a path from the initial state to a goal, its cost and a bound
/** path.front() is the initial state and path.back() the goal; the path
 * is empty unless status == Status::solved, or Status::node_limit after
 * an anytime search has found a goal. bound is a lower bound on the cost
 * \f$ C^* \f$ of an optimal path, which the search has proven: the cost
 * itself when the path is optimal, \f$ \lceil c / w \rceil \f$ (or
 * better) for a path of cost \f$ c \f$ found with the weight \f$ w \f$,
 * so that \f$ c / \mathrm{bound} \f$ bounds the suboptimality. The
 * anytime searches record in improvements each path that lowers the
 * incumbent, as it is found. */

template< class State , class Cost >
struct Solution {

 /// a path found by an anytime search: its cost, the lower bound on
 /// \f$ C^* \f$ known at that moment, the expansions done so far (by the
 /// worker that found it, in a parallel search) and the time
 struct Improvement {
  Cost cost;                             ///< the cost of the path
  Cost bound;                            ///< the lower bound on C*
  std::uint64_t expanded;                ///< the expansions so far
  double seconds;                        ///< the time since the start
  };

 Status status = Status::unsolvable;     ///< how the search ended
 Cost cost{};                            ///< the cost of the path
 Cost bound{};                           ///< a lower bound on C*
 std::vector< State > path;              ///< the states along the path
 std::vector< Improvement > improvements;  ///< the paths of an anytime
                                           ///< search, by time
 SearchStats stats;                      ///< the counters of the search

 /// true if a path has been found
 [[nodiscard]] bool solved( void ) const noexcept {
  return( status == Status::solved );
  }

 };  // end( struct( Solution ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SearchOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options shared by all the searches

struct SearchOptions {

 /// stop with Status::node_limit after storing this many nodes (0 = no
 /// limit); in a parallel search the limit holds for each worker
 std::size_t node_limit = 0;

 /// the number of nodes the search is expected to store (0 = unknown): the
 /// closed list is sized for them at the start, instead of being doubled
 /// (and rehashed) as it fills; in a parallel search it is per worker
 std::size_t expected_nodes = 0;

 /// the weight \f$ w \geq 1 \f$ of the heuristic: OPEN is ordered by
 /// \f$ f' = g + w h \f$, and the path found costs at most \f$ w C^* \f$
 /// [see Weight]; 1 is A*; the greedy searches take any \f$ w \geq 0 \f$,
 /// infinite (h alone) by default [see GreedyOptions]
 double weight = 1;

 /// with a weight above 1, a closed node reached through a cheaper path is
 /// reopened only if the saving exceeds what the weight allows, instead of
 /// always, which keeps the bound \f$ w C^* \f$ [see Weight::improves]
 bool drop_duplicates = false;

 /// with a weight above 1, the search goes on after the first goal, the
 /// incumbent pruning the nodes by \f$ f \f$ instead of \f$ f' \f$, until
 /// it proves the incumbent optimal, recording each better path found [see
 /// Solution::improvements]
 bool anytime = false;

 };  // end( struct( SearchOptions ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Solution */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Solution.h ---------------------------*/
/*--------------------------------------------------------------------------*/
