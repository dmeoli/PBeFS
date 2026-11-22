/*--------------------------------------------------------------------------*/
/*-------------------------- File BBCheckpoint.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the checkpoint and the restart of the Branch-and-Bound
 * (B&B): the bytes of the nodes (which the B&B over MPI sends too), the
 * file of a checkpoint, holding the open nodes of a search with their
 * bounds and depths, the incumbent and its node, and the bound of the
 * region closed so far, and checkpointed_branch_and_bound(), which runs any
 * search in segments of a given number of evaluations, writing a
 * checkpoint after each, and which, given the file of an earlier run,
 * starts from it. The problem, the options and the Solution are those of
 * BranchAndBound.h and ParallelBranchAndBound.h, and the design is
 * discussed in docs/branch-and-bound.md.
 *
 * A segment is a search from the open nodes of the previous one (from the
 * root at first) with the incumbent found so far, the node limit of the
 * segment, and BBOptions::keep_open_only, so that it keeps only the leaves
 * left open by the limit, and the least bound of the others in
 * BBSolution::closed_bound. Its open leaves, with the closed region, cover
 * the region it started from, and therefore, by induction, the open nodes
 * of the last checkpoint and the closed regions of all the segments cover
 * the tree; the bound of the search is the least of the bounds of the
 * closed regions, of the open nodes and of the incumbent. A run stopped at
 * any time (e.g., by the end of the time granted by a batch system) is
 * resumed from the last checkpoint, losing at most the work of a segment,
 * by any scheme and with any number of workers; each segment costs the
 * synchronization of a search that ends and one that starts, which is
 * small when a segment is long w.r.t. the time an evaluation takes.
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

#ifndef __PBeFS_BBCheckpoint
 #define __PBeFS_BBCheckpoint
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>

#include <cstdio>

#include <cstring>

#include <fstream>

#include <iterator>

#include <optional>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a B&B problem whose nodes and values can be written as bytes
/** The values are trivially copyable, and are written as their bytes;
 * either the node is trivially copyable too, or the problem gives
 * serialize( n , out ), which appends the bytes of the node n to the
 * std::vector< char > out, and deserialize( in ), which rebuilds a node
 * from the bytes at the const char * in and moves in past them. A problem
 * may also give deserialize( in , end ), which throws if the bytes from in
 * to end hold no whole node, and which bb_read_checkpoint() then calls, so
 * that it reads nothing past the end of a truncated or corrupted file; with
 * deserialize( in ) alone, it checks only after the node that it did not
 * go past the end. */

//! [branch-and-bound-sig-serializable]
template< class P >
concept SerializableBranchAndBound = BranchAndBound< P > &&
 std::is_trivially_copyable_v< bb_value_t< P > > &&
 ( std::is_trivially_copyable_v< typename P::node_type > ||
   requires( const P & p , const typename P::node_type & n ,
             std::vector< char > & out , const char * & in ) {
    p.serialize( n , out );
    { p.deserialize( in ) } -> std::convertible_to< typename P::node_type >;
    } );
//! [branch-and-bound-sig-serializable]

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// appends the bytes of x to out

template< class T >
void bb_put_bytes( const T & x , std::vector< char > & out )
{
 const std::size_t n = out.size();
 out.resize( n + sizeof( T ) );
 std::memcpy( out.data() + n , & x , sizeof( T ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads an object from the bytes at in, and moves in past them

template< class T >
T bb_get_bytes( const char * & in )
{
 T x;
 std::memcpy( static_cast< void * >( & x ) , in , sizeof( T ) );
 in += sizeof( T );
 return( x );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads an object from the bytes at in, which end at end, and moves in
/// past them; throws std::runtime_error( who ) if fewer than sizeof( T )
/// bytes are left

template< class T >
T bb_get_bytes( const char * & in , const char * end ,
                const std::string & who )
{
 if( std::size_t( end - in ) < sizeof( T ) )
  throw( std::runtime_error( who ) );
 return( bb_get_bytes< T >( in ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// appends the node n to out

template< SerializableBranchAndBound P >
void bb_put_node( const P & p , const typename P::node_type & n ,
                  std::vector< char > & out )
{
 if constexpr( requires { p.serialize( n , out ); } )
  p.serialize( n , out );
 else
  bb_put_bytes( n , out );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads a node from the bytes at in, and moves in past them

template< SerializableBranchAndBound P >
typename P::node_type bb_get_node( const P & p , const char * & in )
{
 if constexpr( requires { p.deserialize( in ); } )
  return( p.deserialize( in ) );
 else
  return( bb_get_bytes< typename P::node_type >( in ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads a node from the bytes at in, which end at end, and moves in past
/// them; throws std::runtime_error( who ) if they hold no whole node [see
/// SerializableBranchAndBound]

template< SerializableBranchAndBound P >
typename P::node_type bb_get_node( const P & p , const char * & in ,
                                   const char * end ,
                                   const std::string & who )
{
 if constexpr( requires { p.deserialize( in , end ); } )
  return( p.deserialize( in , end ) );
 else
  if constexpr( requires { p.deserialize( in ); } ) {
   if( in >= end )
    throw( std::runtime_error( who ) );
   auto n = p.deserialize( in );
   if( in > end )
    throw( std::runtime_error( who ) );
   return( n );
   }
  else
   return( bb_get_bytes< typename P::node_type >( in , end , who ) );
 }

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BBCheckpoint --------------------------*/
/*--------------------------------------------------------------------------*/
/// the state of a B&B saved in a checkpoint, in the sense of the problem
/** open are the open nodes, with their bounds and depths, incumbent the
 * incumbent value and node its node (either may be missing), closed the
 * least bound of the region already closed (\f$ +\infty \f$, or
 * \f$ -\infty \f$ for a maximization, if none), and stats the counters of
 * the segments so far. */

//! [branch-and-bound-sig-checkpoint]
template< class Node , class V = double >
struct BBCheckpoint {
 std::vector< BBLeaf< Node , V > > open;       ///< the open nodes
 std::optional< V > incumbent;                 ///< the incumbent value
 std::optional< Node > node;                   ///< ... and its node
 V closed = BBValueTraits< V >::infinity();    ///< the bound of the closed
 BBStats stats;                                ///< the counters so far
 };
//! [branch-and-bound-sig-checkpoint]

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {
 /// the first bytes of a checkpoint
 inline constexpr char bb_magic[ 8 ] = { 'P' , 'B' , 'e' , 'F' , 'S' , 'c' ,
                                         'k' , '1' };
 }  // end( namespace detail )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// writes the checkpoint c of the problem p to the file, first to a
/// temporary one that then takes its name, so that a run stopped while
/// writing leaves the previous checkpoint whole; throws std::runtime_error
/// if it cannot be written

template< SerializableBranchAndBound P >
void bb_write_checkpoint( const P & p ,
                          const BBCheckpoint< typename P::node_type ,
                                              bb_value_t< P > > & c ,
                          const std::string & file )
{
 using V = bb_value_t< P >;
 std::vector< char > b( std::begin( detail::bb_magic ) ,
                        std::end( detail::bb_magic ) );
 detail::bb_put_bytes( std::uint32_t( sizeof( V ) ) , b );
 detail::bb_put_bytes( std::uint8_t( c.incumbent ? 1 : 0 ) , b );
 detail::bb_put_bytes( c.incumbent ? *c.incumbent : V( 0 ) , b );
 detail::bb_put_bytes( std::uint8_t( c.node ? 1 : 0 ) , b );
 if( c.node )
  detail::bb_put_node( p , *c.node , b );
 detail::bb_put_bytes( c.closed , b );
 const std::uint64_t k[ 12 ] = {
  c.stats.evaluated , c.stats.branched , c.stats.generated ,
  c.stats.infeasible , c.stats.pruned_generated , c.stats.pruned_extracted ,
  c.stats.pruned_evaluated , c.stats.incumbents , c.stats.solutions ,
  c.stats.recomputed , c.stats.sent , c.stats.messages };
 for( const auto x : k )
  detail::bb_put_bytes( x , b );
 detail::bb_put_bytes( std::uint64_t( c.open.size() ) , b );
 for( const auto & l : c.open ) {
  detail::bb_put_bytes( l.bound , b );
  detail::bb_put_bytes( std::uint32_t( l.depth ) , b );
  detail::bb_put_node( p , l.node , b );
  }
 const std::string tmp = file + ".tmp";
 {
  std::ofstream out( tmp , std::ios::binary | std::ios::trunc );
  out.write( b.data() , std::streamsize( b.size() ) );
  if( ! out )
   throw( std::runtime_error( "bb_write_checkpoint: cannot write " + tmp ) );
  }
 std::remove( file.c_str() );
 if( std::rename( tmp.c_str() , file.c_str() ) != 0 )
  throw( std::runtime_error( "bb_write_checkpoint: cannot rename " + tmp ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the checkpoint of the problem p in the file, or nothing if there is no
/// file; throws std::runtime_error if the file is not a checkpoint of the
/// values of P, or if it is truncated or corrupted, before reading past
/// its end [see SerializableBranchAndBound]

template< SerializableBranchAndBound P >
std::optional< BBCheckpoint< typename P::node_type , bb_value_t< P > > >
 bb_read_checkpoint( const P & p , const std::string & file )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 std::ifstream in( file , std::ios::binary );
 if( ! in )
  return( std::nullopt );
 const std::vector< char > b( ( std::istreambuf_iterator< char >( in ) ) ,
                              std::istreambuf_iterator< char >() );
 const std::size_t h = sizeof( detail::bb_magic );
 if( ( b.size() < h + 4 ) ||
     ( std::memcmp( b.data() , detail::bb_magic , h ) != 0 ) )
  throw( std::runtime_error( "bb_read_checkpoint: not a checkpoint " +
                             file ) );
 const char * at = b.data() + h;
 const char * const end = b.data() + b.size();
 const std::string bad = "bb_read_checkpoint: truncated or corrupted " +
                         file;
 if( detail::bb_get_bytes< std::uint32_t >( at , end , bad ) != sizeof( V ) )
  throw( std::runtime_error( "bb_read_checkpoint: other values in " +
                             file ) );
 BBCheckpoint< Node , V > c;
 const auto has_inc = detail::bb_get_bytes< std::uint8_t >( at , end , bad );
 const V inc = detail::bb_get_bytes< V >( at , end , bad );
 if( has_inc > 1 )
  throw( std::runtime_error( bad ) );
 if( has_inc )
  c.incumbent = inc;
 const auto has_node = detail::bb_get_bytes< std::uint8_t >( at , end , bad );
 if( has_node > 1 )
  throw( std::runtime_error( bad ) );
 if( has_node )
  c.node.emplace( detail::bb_get_node( p , at , end , bad ) );
 c.closed = detail::bb_get_bytes< V >( at , end , bad );
 std::uint64_t k[ 12 ];
 for( auto & x : k )
  x = detail::bb_get_bytes< std::uint64_t >( at , end , bad );
 c.stats.evaluated = k[ 0 ];
 c.stats.branched = k[ 1 ];
 c.stats.generated = k[ 2 ];
 c.stats.infeasible = k[ 3 ];
 c.stats.pruned_generated = k[ 4 ];
 c.stats.pruned_extracted = k[ 5 ];
 c.stats.pruned_evaluated = k[ 6 ];
 c.stats.incumbents = k[ 7 ];
 c.stats.solutions = k[ 8 ];
 c.stats.recomputed = k[ 9 ];
 c.stats.sent = k[ 10 ];
 c.stats.messages = k[ 11 ];
 const std::uint64_t n = detail::bb_get_bytes< std::uint64_t >( at , end ,
                                                                bad );
 // an open node takes its bound and its depth at least
 if( n > std::uint64_t( end - at ) / ( sizeof( V ) + 4 ) )
  throw( std::runtime_error( bad ) );
 c.open.reserve( std::size_t( n ) );
 for( std::uint64_t i = 0 ; i < n ; ++i ) {
  const V bound = detail::bb_get_bytes< V >( at , end , bad );
  const auto depth = detail::bb_get_bytes< std::uint32_t >( at , end , bad );
  c.open.push_back( BBLeaf< Node , V >{ detail::bb_get_node( p , at , end ,
                                                             bad ) ,
                                        bound , depth , BBLeafKind::open } );
  }
 if( at != end )
  throw( std::runtime_error( "bb_read_checkpoint: bad size of " + file ) );
 return( c );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the B&B of the problem by search( problem , options , starts ), e.g., a
/// call of parallel_branch_and_bound() with a scheme, run in segments of
/// every evaluations, a checkpoint being written to file after each; if
/// file holds a checkpoint, the search starts from it
/** The options of each segment are o, with the incumbent found so far, the
 * node limit of the segment (and what is left of that of o, which counts
 * the evaluations of the earlier runs too) and keep_open_only; the file is
 * left with the last checkpoint, which has no open node if the search
 * ended. The counters of the Solution are those of all the segments, of
 * the earlier runs included, and its frontier is that of the last one. */

template< SerializableBranchAndBound P , class Search >
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
 checkpointed_branch_and_bound( const P & problem ,
                                const BasicParallelBBOptions<
                                 bb_value_t< P > > & o ,
                                Search && search , const std::string & file ,
                                std::size_t every )
{
 using Node = typename P::node_type;
 using V = bb_value_t< P >;
 using Tr = BBValueTraits< V >;
 using Leaf = BBLeaf< Node , V >;
 Stopwatch clock;
 const V inf = detail::bb_signed< P >( Tr::infinity() );   // in the sense
 const auto better = [ & ]( V a , V b ) {                  // of P
  return( detail::bb_signed< P >( a ) < detail::bb_signed< P >( b ) );
  };

 BBCheckpoint< Node , V > c;
 if( auto r = bb_read_checkpoint( problem , file ) )
  c = std::move( *r );
 else {
  c.closed = inf;
  c.open.push_back( Leaf{ problem.root() , - inf , 0 , BBLeafKind::open } );
  if( o.incumbent )
   c.incumbent = o.incumbent;
  }

 BBSolution< Node , V > sol;
 sol.status = BBStatus::optimal;
 for( ;; ) {
  if( c.open.empty() ||
      ( o.node_limit && ( c.stats.evaluated >= o.node_limit ) ) )
   break;
  BasicParallelBBOptions< V > so = o;
  so.keep_frontier = true;
  so.keep_open_only = true;
  so.incumbent = c.incumbent;
  so.node_limit = every;
  if( o.node_limit ) {
   const std::size_t left = o.node_limit - std::size_t( c.stats.evaluated );
   so.node_limit = every ? std::min( every , left ) : left;
   }
  if( o.time_limit > 0 )
   so.time_limit = std::max( 1e-9 , o.time_limit - clock.seconds() );
  sol = search( problem , so , std::move( c.open ) );
  c.open.clear();
  detail::bb_add_counters( c.stats , sol.stats );
  c.stats.max_open = std::max( c.stats.max_open , sol.stats.max_open );
  if( sol.node && ( ! c.incumbent || better( sol.value , *c.incumbent ) ) ) {
   c.incumbent = sol.value;
   c.node = sol.node;
   }
  if( better( sol.closed_bound , c.closed ) )
   c.closed = sol.closed_bound;
  const bool limited = ( sol.status == BBStatus::node_limit ) ||
                       ( sol.status == BBStatus::time_limit );
  if( limited )
   for( auto & l : sol.frontier )
    c.open.push_back( std::move( l ) );
  sol.frontier.clear();
  bb_write_checkpoint( problem , c , file );
  if( ( ! limited ) || ( sol.status == BBStatus::time_limit ) )
   break;
  }
 if( c.open.empty() )
  bb_write_checkpoint( problem , c , file );

 // the Solution of all the segments
 BBSolution< Node , V > out;
 const bool ended = c.open.empty();
 out.status = ended ? ( c.incumbent ? BBStatus::optimal
                                    : BBStatus::infeasible )
                    : ( sol.status == BBStatus::time_limit ?
                        BBStatus::time_limit : BBStatus::node_limit );
 if( ended && ( o.mode == BBMode::decide ) && c.incumbent )
  out.status = detail::bb_status< P >( detail::bb_signed< P >(
                                        *c.incumbent ) , o );
 V lb = detail::bb_signed< P >( c.closed );
 for( const auto & l : c.open )
  lb = std::min( lb , detail::bb_signed< P >( l.bound ) );
 if( c.incumbent ) {
  out.value = *c.incumbent;
  lb = std::min( lb , detail::bb_signed< P >( *c.incumbent ) );
  }
 out.node = c.node;
 out.bound = detail::bb_signed< P >( lb );
 out.closed_bound = c.closed;
 out.gap = c.incumbent && Tr::finite( lb ) ?
           V( detail::bb_signed< P >( *c.incumbent ) - lb ) : Tr::infinity();
 out.stats = c.stats;
 out.stats.seconds = clock.seconds();
 out.frontier = std::move( c.open );
 return( out );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BBCheckpoint */

/*--------------------------------------------------------------------------*/
/*------------------------ End File BBCheckpoint.h -------------------------*/
/*--------------------------------------------------------------------------*/
