/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndBoundScheme.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the choice of the Branch-and-Bound (B&B) by its scheme:
 * BBScheme names the searches of the core of PBeFS, bb_scheme() reads it
 * from its name (e.g., from a file of parameters), and
 * parallel_branch_and_bound() runs the search of the scheme of the options
 * on a runtime, so that a caller chooses the search by a value rather than
 * by a type. The searches that need an external library, i.e., the farm of
 * FastFlow and the B&B over MPI, are not among them, and they are used by
 * their own classes [see FarmBranchAndBound.h and MPIBranchAndBound.h].
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_BranchAndBoundScheme
 #define __PBeFS_BranchAndBoundScheme
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <optional>

#include <stdexcept>

#include <string>

#include <string_view>

#include <utility>

#include <vector>

#include "PBeFS/bnb/DeterministicBranchAndBound.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "PBeFS/bnb/SelfSplitBranchAndBound.h"

#include "PBeFS/bnb/StealingBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- ENUM BBScheme ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the schemes of the B&B of the core
/** Each has a name, given by bb_scheme_name() and read by bb_scheme(). */

enum class BBScheme {
 sequential ,      ///< one worker [see BranchAndBoundSearch]
 shared ,          ///< one shared queue [see SharedBranchAndBound]
 distributed ,     ///< a queue per worker [see DistributedBranchAndBound]
 deterministic ,   ///< synchronous rounds [see DeterministicBranchAndBound]
 selfsplit ,       ///< a split of a sampling [see SelfSplitBranchAndBound]
 stealing          ///< work stealing [see StealingBranchAndBound]
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {
 /// the names of the schemes, in the order of BBScheme
 inline constexpr std::array< std::string_view , 6 > bb_scheme_names = {
  "sequential" , "shared" , "distributed" , "deterministic" , "selfsplit" ,
  "stealing" };
 }  // end( namespace detail )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the name of the scheme s

[[nodiscard]] inline std::string_view bb_scheme_name( BBScheme s ) {
 return( detail::bb_scheme_names.at( std::size_t( s ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the scheme of the given name; throws std::invalid_argument if there is
/// none

[[nodiscard]] inline BBScheme bb_scheme( std::string_view name ) {
 for( std::size_t i = 0 ; i < detail::bb_scheme_names.size() ; ++i )
  if( detail::bb_scheme_names[ i ] == name )
   return( BBScheme( i ) );
 throw( std::invalid_argument( "bb_scheme: unknown scheme " +
                               std::string( name ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the B&B of the problem by the scheme s on the given runtime, with the
/// given options, from the given nodes [see BranchAndBoundSearch::solve()]
/// or from the root if none

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool >
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
 parallel_branch_and_bound( const P & problem , BBScheme s ,
                            const BasicParallelBBOptions< bb_value_t< P > > &
                             o = {} ,
                            std::optional< std::vector< BBLeaf<
                             typename P::node_type , bb_value_t< P > > > >
                             seeds = std::nullopt )
{
 const auto solve = [ & ]( const auto & b ) {
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 switch( s ) {
  case( BBScheme::sequential ):
   return( solve( BranchAndBoundSearch< P >( problem , o ) ) );
  case( BBScheme::shared ):
   return( solve( SharedBranchAndBound< P , Executor >( problem , o ) ) );
  case( BBScheme::distributed ):
   return( solve( DistributedBranchAndBound< P , Executor >( problem ,
                                                             o ) ) );
  case( BBScheme::deterministic ):
   return( solve( DeterministicBranchAndBound< P , Executor >( problem ,
                                                               o ) ) );
  case( BBScheme::selfsplit ):
   return( solve( SelfSplitBranchAndBound< P , Executor >( problem , o ) ) );
  default:
   return( solve( StealingBranchAndBound< P , Executor >( problem , o ) ) );
  }
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BranchAndBoundScheme */

/*--------------------------------------------------------------------------*/
/*-------------------- End File BranchAndBoundScheme.h ---------------------*/
/*--------------------------------------------------------------------------*/
