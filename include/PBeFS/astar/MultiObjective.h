/*--------------------------------------------------------------------------*/
/*-------------------------- File MultiObjective.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the parts shared by the multi-objective searches: the
 * concept MultiObjectiveDomain, the vectors of costs and their dominance,
 * the outcome of a search (ParetoSolution, the Pareto front of the paths
 * from the initial state to a goal), its options, and the frontiers, the
 * sets of truncated vectors by which the searches decide whether a path is
 * dominated.
 *
 * A path \f$ \pi \f$ costs a vector \f$ c( \pi ) \in \mathbb{R}^K_{\geq 0}
 * \f$, the sum of the vectors of its edges; \f$ u \preceq v \f$ (u weakly
 * dominates v) if \f$ u_i \leq v_i \f$ for every i, and \f$ u \prec v \f$
 * (u dominates v) if moreover \f$ u \neq v \f$. A path to a goal is
 * Pareto-optimal if no path to a goal dominates it, and the searches
 * return a cost-unique Pareto-optimal set, i.e., one path for each vector of
 * the Pareto front \f$ \mathcal{F} \f$, the set of the costs of the
 * Pareto-optimal paths [see Mandow and Pérez de la Cruz, JACM 2010]. The
 * truncation \f$ \mathrm{Tr}( v ) = ( v_2 , \ldots , v_K ) \f$ drops the
 * first component: a search that extracts its paths in the lexicographic
 * order of \f$ f = g + h \f$ with a consistent h extracts them by
 * non-decreasing \f$ f_1 \f$, so that a path extracted later is weakly
 * dominated by one extracted earlier as soon as its truncation is [see
 * Pulido, Mandow and Pérez de la Cruz, C&OR 2015].
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

#ifndef __PBeFS_MultiObjective
 #define __PBeFS_MultiObjective
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <numeric>

#include <stdexcept>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/astar/MOStructures.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// a callable accepting a successor with its vectors, used only to state
 /// the concepts
 template< class State , class Vec >
 struct mo_successor_sink {
  void operator()( const State & , const Vec & , const Vec & ) const
   noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// a multi-objective search domain
/** A type D is a MultiObjectiveDomain if it provides:
 *
 * - the type state_type, a regular value type compared with ==;
 *
 * - the type cost_type, an arithmetic type for the components of the
 *   vectors of g, h and the edge costs, and the constant objectives, the
 *   number \f$ K \geq 2 \f$ of components;
 *
 * - initial_state(), is_goal( s ) and hash( s ), as a Domain [see
 *   Domain.h];
 *
 * - heuristic( s ), the vector \f$ h( s ) \f$ of
 *   std::array< cost_type , objectives >;
 *
 * - for_each_successor( s , h_s , f ), which calls f( t , c , h_t ) for
 *   every successor t of s, c being the vector of the edge ( s , t ) and
 *   h_t = h( t ).
 *
 * The costs are nonnegative, and the searches assume h consistent, i.e.,
 * \f$ h( s ) \preceq c( s , t ) + h( t ) \f$ for every edge and
 * \f$ h( s ) = 0 \f$ at every goal, each component a consistent heuristic
 * of its own objective; a component equal to infinite_cost() marks a state
 * from which no goal is reachable, which the searches do not generate. */

template< class D >
concept MultiObjectiveDomain =
 requires( const D & d , const typename D::state_type & s ,
           const std::array< typename D::cost_type , D::objectives > & v ) {
  typename D::state_type;
  typename D::cost_type;
  requires std::is_arithmetic_v< typename D::cost_type >;
  requires ( std::size_t( D::objectives ) >= 2 );
  requires std::equality_comparable< typename D::state_type >;
  requires std::copyable< typename D::state_type >;
  { d.initial_state() } -> std::convertible_to< typename D::state_type >;
  { d.is_goal( s ) } -> std::convertible_to< bool >;
  { d.hash( s ) } -> std::convertible_to< std::uint64_t >;
  { d.heuristic( s ) } -> std::convertible_to<
                          std::array< typename D::cost_type ,
                                      D::objectives > >;
  d.for_each_successor( s , v , detail::mo_successor_sink<
                         typename D::state_type , std::array<
                          typename D::cost_type , D::objectives > >{} );
  };

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the state type of a MultiObjectiveDomain
template< MultiObjectiveDomain D >
using mo_state_t = typename D::state_type;

/// the type of a component of the costs of a MultiObjectiveDomain
template< MultiObjectiveDomain D >
using mo_scalar_t = typename D::cost_type;

/// the number of objectives of a MultiObjectiveDomain
template< MultiObjectiveDomain D >
inline constexpr std::size_t objectives_of = std::size_t( D::objectives );

/// the vector of the costs of a MultiObjectiveDomain
template< MultiObjectiveDomain D >
using mo_cost_t = std::array< mo_scalar_t< D > , objectives_of< D > >;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// true if \f$ u_i \leq v_i \f$ for every \f$ i \geq \mathit{from} \f$
/// (\f$ u \preceq v \f$ for from = 0)

template< class T , std::size_t K >
[[nodiscard]] constexpr bool weakly_dominates( const std::array< T , K > & u ,
                                               const std::array< T , K > & v ,
                                               std::size_t from = 0 )
 noexcept
{
 for( std::size_t i = from ; i < K ; ++i )
  if( v[ i ] < u[ i ] )
   return( false );
 return( true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// true if \f$ u \prec v \f$

template< class T , std::size_t K >
[[nodiscard]] constexpr bool dominates( const std::array< T , K > & u ,
                                        const std::array< T , K > & v )
 noexcept
{
 return( weakly_dominates( u , v ) && ( u != v ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// true if some component of v is infinite_cost(), the mark of a state
/// from which no goal is reachable

template< class T , std::size_t K >
[[nodiscard]] constexpr bool dead_end( const std::array< T , K > & v )
 noexcept
{
 for( const T & x : v )
  if( x == infinite_cost< T >() )
   return( true );
 return( false );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the vectors of v that no other one weakly dominates, each once, in
/// lexicographic order: the Pareto front of v

template< class T , std::size_t K >
[[nodiscard]] std::vector< std::array< T , K > > pareto_front(
 std::vector< std::array< T , K > > v )
{
 std::sort( v.begin() , v.end() );
 v.erase( std::unique( v.begin() , v.end() ) , v.end() );
 // in lexicographic order no vector dominates one before it
 std::vector< std::array< T , K > > front;
 for( const auto & x : v )
  if( std::none_of( front.begin() , front.end() ,
                    [ & ]( const auto & y ) {
                     return( weakly_dominates( y , x ) );
                     } ) )
   front.push_back( x );
 return( front );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the permutation of the objectives order (the identity if empty), which
/// has to be one of \f$ 0 , \ldots , K - 1 \f$; throws
/// std::invalid_argument otherwise

template< std::size_t K >
[[nodiscard]] std::array< unsigned , K > objective_order(
 const std::vector< unsigned > & order )
{
 std::array< unsigned , K > p;
 std::iota( p.begin() , p.end() , 0u );
 if( order.empty() )
  return( p );
 if( order.size() != K )
  throw( std::invalid_argument( "objective_order: not K objectives" ) );
 std::copy( order.begin() , order.end() , p.begin() );
 auto q = p;
 std::sort( q.begin() , q.end() );
 for( std::size_t i = 0 ; i < K ; ++i )
  if( q[ i ] != i )
   throw( std::invalid_argument( "objective_order: not a permutation" ) );
 return( p );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the vector v in the order p, i.e., \f$ ( v_{p_0} , \ldots ,
/// v_{p_{K-1}} ) \f$

template< class T , std::size_t K >
[[nodiscard]] constexpr std::array< T , K > permute(
 const std::array< T , K > & v , const std::array< unsigned , K > & p )
 noexcept
{
 std::array< T , K > w;
 for( std::size_t i = 0 ; i < K ; ++i )
  w[ i ] = v[ p[ i ] ];
 return( w );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the vector w back from the order p, the inverse of permute()

template< class T , std::size_t K >
[[nodiscard]] constexpr std::array< T , K > unpermute(
 const std::array< T , K > & w , const std::array< unsigned , K > & p )
 noexcept
{
 std::array< T , K > v;
 for( std::size_t i = 0 ; i < K ; ++i )
  v[ p[ i ] ] = w[ i ];
 return( v );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MOSearchStats ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the counters of a multi-objective search
/** Those of SearchStats, where expanded counts the paths (the labels)
 * whose successors are generated, the paths to a goal excluded, generated
 * the successors generated, stored the labels allocated, and besides:
 * extracted, the labels taken from OPEN, pruned or not; inserted, the
 * labels put into OPEN, the root included; checks, the dominance checks,
 * i.e., the comparisons of two vectors (or of two components, for the
 * frontiers of one component); solutions, the paths to a goal found, which
 * a parallel search may find more than once. */

struct MOSearchStats : public SearchStats {

 std::uint64_t extracted = 0;     ///< labels taken from OPEN
 std::uint64_t inserted = 0;      ///< labels put into OPEN
 std::uint64_t checks = 0;        ///< dominance checks
 std::uint64_t solutions = 0;     ///< paths to a goal found

 };  // end( struct( MOSearchStats ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ParetoSolution --------------------------*/
/*--------------------------------------------------------------------------*/
/// the Pareto front of the paths from the initial state to a goal
/** front has one path for each vector of the Pareto front, in the
 * lexicographic order of the costs, which are given in the order of the
 * objectives of the Domain whatever the order of the search; a path starts
 * at the initial state and ends at a goal (it is empty if the search was
 * asked not to keep the paths). status is Status::solved if the front is
 * complete, also when it is empty because no goal is reachable
 * (Status::unsolvable then), and Status::node_limit if the search stopped
 * at its limit, with the paths found so far, each Pareto-optimal but for
 * HDMOAStar, whose paths are then only mutually non-dominated. */

template< class State , class Cost , std::size_t K >
struct ParetoSolution {

 /// a Pareto-optimal path and its cost
 struct Path {
  std::array< Cost , K > cost;           ///< the vector of its costs
  std::vector< State > states;           ///< the states along it
  };

 Status status = Status::unsolvable;     ///< how the search ended
 std::vector< Path > front;              ///< the Pareto-optimal paths
 MOSearchStats stats;                    ///< the counters of the search

 /// true if the front is complete and not empty
 [[nodiscard]] bool solved( void ) const noexcept {
  return( status == Status::solved );
  }

 /// the vectors of the front, in lexicographic order
 [[nodiscard]] std::vector< std::array< Cost , K > > costs( void ) const {
  std::vector< std::array< Cost , K > > c;
  c.reserve( front.size() );
  for( const Path & p : front )
   c.push_back( p.cost );
  return( c );
  }

 };  // end( struct( ParetoSolution ) )

/// the ParetoSolution of a MultiObjectiveDomain
template< MultiObjectiveDomain D >
using pareto_solution_t = ParetoSolution< mo_state_t< D > , mo_scalar_t< D > ,
                                          objectives_of< D > >;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MOSearchOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the multi-objective searches

struct MOSearchOptions {

 /// stop with Status::node_limit after storing this many labels (0 = no
 /// limit); in a parallel search the limit holds for each worker
 std::size_t node_limit = 0;

 /// the number of states the search is expected to reach (0 = unknown),
 /// for which the table of the states is sized at the start
 std::size_t expected_states = 0;

 /// no dominance check when a path is generated, only when it is
 /// extracted: the lazy searches of Hernández et al. (IJCAI 2023), which
 /// put into OPEN more paths and check fewer
 bool lazy = false;

 /// the order of the objectives in the lexicographic order of OPEN, a
 /// permutation of \f$ 0 , \ldots , K - 1 \f$, the first one the primary
 /// objective; empty for \f$ 0 , \ldots , K - 1 \f$
 std::vector< unsigned > order;

 /// keep the paths of the front, not only their costs
 bool paths = true;

 };  // end( struct( MOSearchOptions ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the truncation \f$ \mathrm{Tr}( v ) = ( v_2 , \ldots , v_K ) \f$

template< class T , std::size_t K >
[[nodiscard]] constexpr std::array< T , K - 1 > truncate(
 const std::array< T , K > & v ) noexcept
{
 std::array< T , K - 1 > k;
 for( std::size_t i = 1 ; i < K ; ++i )
  k[ i - 1 ] = v[ i ];
 return( k );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ScalarFrontier --------------------------*/
/*--------------------------------------------------------------------------*/
/// the frontier of a state for two objectives: the least second component
/** Since the paths are extracted by non-decreasing \f$ g_1 \f$ at a state,
 * the truncated vectors of the paths expanded there are numbers, and one of
 * them, their minimum \f$ g_2^{\min} \f$, weakly dominates them all: the
 * dominance checks of BOA* in constant time [see Hernández et al., AIJ
 * 2023]. dominated( v , n ) tells whether \f$ \mathrm{Tr}( v ) \f$ is
 * weakly dominated by the frontier, insert( v , n ) adds it, which removes
 * what it dominates; n counts the checks. */

template< class Cost , std::size_t K >
class ScalarFrontier {

 static_assert( K == 2 , "ScalarFrontier is for two objectives" );

 public:

 using Vec = std::array< Cost , K >;              ///< a vector of costs

 /// the name of the frontier in the names of the searches
 static constexpr const char * name = "scalar";

 [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const
  noexcept {
  ++n;
  return( f_min <= v[ 1 ] );
  }

 void insert( const Vec & v , std::uint64_t & ) noexcept { f_min = v[ 1 ]; }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_min == infinite_cost< Cost >() ? 0 : 1 );
  }

 private:

 Cost f_min = infinite_cost< Cost >();            ///< g_2^min

 };  // end( class( ScalarFrontier ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ArrayFrontier ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the frontier of a state as an array of truncated vectors
/** The non-dominated truncated vectors, in no order: a check scans them
 * all, an insertion removes the dominated ones by moving the last vector
 * into their place, the frontier of LTMOA* [see Hernández et al., IJCAI
 * 2023]. */

template< class Cost , std::size_t K >
class ArrayFrontier {

 public:

 using Vec = std::array< Cost , K >;              ///< a vector of costs

 /// the name of the frontier in the names of the searches
 static constexpr const char * name = "array";

 [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const
  noexcept {
  const Key b = truncate( v );
  for( const Key & u : v_keys ) {
   ++n;
   if( weakly_dominates( u , b ) )
    return( true );
   }
  return( false );
  }

 void insert( const Vec & v , std::uint64_t & n ) {
  const Key b = truncate( v );
  for( std::size_t i = 0 ; i < v_keys.size() ; ) {
   ++n;
   if( weakly_dominates( b , v_keys[ i ] ) ) {
    v_keys[ i ] = v_keys.back();
    v_keys.pop_back();
    }
   else
    ++i;
   }
  v_keys.push_back( b );
  }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_keys.size() );
  }

 private:

 using Key = std::array< Cost , K - 1 >;          ///< a truncated vector

 std::vector< Key > v_keys;                       ///< the vectors

 };  // end( class( ArrayFrontier ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS TreeFrontier ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the frontier of a state in a balanced search tree
/** The non-dominated truncated vectors in a LexSet, the frontier of EMOA*
 * [see Ren et al., SoCS 2022]: since \f$ u \preceq v \f$ implies
 * \f$ u \leq_{lex} v \f$, a check looks only at the vectors
 * \f$ \leq_{lex} \mathrm{Tr}( v ) \f$ and an insertion removes only some of
 * the vectors \f$ \geq_{lex} \mathrm{Tr}( v ) \f$. For three objectives the
 * truncated vectors have two components, the second one decreasing as the
 * first one increases, so that \f$ \mathrm{Tr}( v ) \f$ is dominated if
 * and only if the largest vector \f$ \leq_{lex} \mathrm{Tr}( v ) \f$ is,
 * one check in \f$ O( \log n ) \f$, and the vectors it dominates follow it
 * in a run (TOA* of the same paper). */

template< class Cost , std::size_t K >
class TreeFrontier {

 static_assert( K >= 3 , "TreeFrontier is for three objectives or more" );

 public:

 using Vec = std::array< Cost , K >;              ///< a vector of costs

 /// the name of the frontier in the names of the searches
 static constexpr const char * name = "tree";

 [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const {
  const Key b = truncate( v );
  if constexpr( K == 3 ) {
   const Key * k = f_set.floor( b );
   ++n;
   return( k && ( ( *k )[ 1 ] <= b[ 1 ] ) );
   }
  else
   return( f_set.any_upto( b , [ & ]( const Key & k ) {
    ++n;
    return( weakly_dominates( k , b ) );
    } ) );
  }

 void insert( const Vec & v , std::uint64_t & n ) {
  const Key b = truncate( v );
  f_set.erase_from( b , [ & ]( const Key & k ) {
   ++n;
   return( weakly_dominates( b , k ) );
   } , K == 3 );
  f_set.insert( b );
  }

 [[nodiscard]] std::size_t size( void ) const { return( f_set.size() ); }

 private:

 using Key = std::array< Cost , K - 1 >;          ///< a truncated vector

 LexSet< Key > f_set;                             ///< the vectors

 };  // end( class( TreeFrontier ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS MOBucketQueue ----------------------------*/
/*--------------------------------------------------------------------------*/
/// an OPEN list of buckets by the first component of f, for integer costs
/** The interface of MOHeap [see MOStructures.h], for the items whose
 * member f is a vector of integers: a bucket for each value of
 * \f$ f_1 \f$, each an MOHeap by Better, on a ring that doubles when a
 * value too far from the least one arrives, as the bucket queues of Ahmadi
 * et al. (ICAPS 2025) but with the ties of \f$ f_1 \f$ kept in the order
 * Better, which keeps the dimensionality reduction exact. The pushes must
 * not go below the first one nor below the last \f$ f_1 \f$ extracted, as
 * in the searches with a consistent heuristic; the items come out in the
 * order of Better if Better orders first by \f$ f_1 \f$. */

template< class Item , class Better >
class MOBucketQueue {

 using Cost = std::remove_cvref_t< decltype( std::declval< Item >().f[ 0 ] ) >;

 static_assert( std::is_integral_v< Cost > ,
                "MOBucketQueue requires integer costs" );

 public:

 void push( const Item & x ) {
  const Cost k = x.f[ 0 ];
  if( v_ring.empty() ) {                 // the first push: the least f_1
   v_ring.resize( 64 );
   f_low = k;
   }
  while( std::size_t( k - f_low ) >= v_ring.size() )
   grow();
  v_ring[ ( f_at + std::size_t( k - f_low ) ) & ( v_ring.size() - 1 ) ]
   .push( x );
  ++f_size;
  }

 [[nodiscard]] const Item & top( void ) {
  advance();
  return( v_ring[ f_at ].top() );
  }

 Item pop( void ) {
  advance();
  --f_size;
  return( v_ring[ f_at ].pop() );
  }

 [[nodiscard]] bool empty( void ) const noexcept { return( f_size == 0 ); }

 [[nodiscard]] std::size_t size( void ) const noexcept { return( f_size ); }

 private:

 // the bucket of f_low the first one not empty
 void advance( void ) {
  while( v_ring[ f_at ].empty() ) {
   f_at = ( f_at + 1 ) & ( v_ring.size() - 1 );
   ++f_low;
   }
  }

 // the ring doubled, its buckets in order from f_low
 void grow( void ) {
  std::vector< MOHeap< Item , Better > > r( 2 * v_ring.size() );
  for( std::size_t i = 0 ; i < v_ring.size() ; ++i )
   std::swap( r[ i ] , v_ring[ ( f_at + i ) & ( v_ring.size() - 1 ) ] );
  v_ring.swap( r );
  f_at = 0;
  }

 std::vector< MOHeap< Item , Better > > v_ring;   ///< the buckets
 std::size_t f_at = 0;                             ///< the bucket of f_low
 Cost f_low{};                                     ///< the least f_1
 std::size_t f_size = 0;                           ///< the items

 };  // end( class( MOBucketQueue ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS LayerIndex -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a set of full vectors in layers of their first component, for the
/// dominance checks of vectors in any order
/** The vectors are kept by increasing first component in blocks of B to
 * 4 B of them, the layers, and a Fenwick tree over the layers keeps, for
 * each of its nodes, the truncations of the vectors of the layers it
 * covers in a TreeFrontier. A vector u of a layer before the one of v,
 * whose first component is at most \f$ v_1 \f$, weakly dominates v as soon
 * as \f$ \mathrm{Tr}( u ) \f$ weakly dominates \f$ \mathrm{Tr}( v ) \f$:
 * a check asks the \f$ O( \log L ) \f$ nodes that cover the layers before
 * that of v, each a check of TreeFrontier (one comparison for three
 * objectives), and compares v with the vectors of its own layer whose
 * first component is not larger. This is the indexed pruning of Xu et al.
 * (2026) on the prefixes of a Fenwick tree [see Fenwick, Software:
 * Practice and Experience 1994]. The vectors that a later one dominates are
 * kept in the layers: a check answers whether a vector added weakly
 * dominates v, which a dominated vector does only if the one that dominates
 * it does too. dominated( v , n ) and insert( v , n ) count the
 * comparisons in n. */

template< class Cost , std::size_t K >
class LayerIndex {

 static_assert( K >= 3 , "LayerIndex is for three objectives or more" );

 public:

 using Vec = std::array< Cost , K >;              ///< a vector of costs

 [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const {
  if( v_layer.empty() )
   return( false );
  const std::size_t b = layer_of( v );
  for( std::size_t i = b ; i > 0 ; i -= i & ( ~i + 1 ) )
   if( v_node[ i ].dominated( v , n ) )
    return( true );
  for( const Vec & z : v_layer[ b ] ) {
   if( v[ 0 ] < z[ 0 ] )
    break;
   ++n;
   if( weakly_dominates( z , v ) )
    return( true );
   }
  return( false );
  }

 void insert( const Vec & v , std::uint64_t & n ) {
  if( v_layer.empty() ) {
   v_layer.emplace_back();
   build( n );
   }
  std::size_t b = layer_of( v );
  if( ( b + 1 == v_layer.size() ) && ( v_layer[ b ].size() >= 2 * B ) &&
      ( v_layer[ b ].back()[ 0 ] <= v[ 0 ] ) ) {
   // a new layer at the end, its node made of the layers it covers
   v_layer.emplace_back();
   ++b;
   if( b + 1 >= v_node.size() ) {
    v_layer.back().push_back( v );
    build( n );
    return;
    }
   const std::size_t i = b + 1;
   for( std::size_t j = i - ( i & ( ~i + 1 ) ) ; j < b ; ++j )
    for( const Vec & z : v_layer[ j ] )
     add( v_node[ i ] , z , n );
   }
  auto & L = v_layer[ b ];
  L.insert( std::upper_bound( L.begin() , L.end() , v ,
                              []( const Vec & a , const Vec & c ) {
                               return( a[ 0 ] < c[ 0 ] );
                               } ) , v );
  if( L.size() > 4 * B ) {
   build( n );
   return;
   }
  for( std::size_t i = b + 1 ; i < v_node.size() ; i += i & ( ~i + 1 ) )
   add( v_node[ i ] , v , n );
  }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  std::size_t s = 0;
  for( const auto & L : v_layer )
   s += L.size();
  return( s );
  }

 private:

 static constexpr std::size_t B = 8;              ///< the least layer

 // the layer of v: the last one whose first vector is not above v_1
 [[nodiscard]] std::size_t layer_of( const Vec & v ) const {
  const auto it = std::upper_bound( v_layer.begin() + 1 , v_layer.end() ,
                                    v[ 0 ] , []( Cost c , const auto & L ) {
                                     return( c < L.front()[ 0 ] );
                                     } );
  return( std::size_t( it - v_layer.begin() ) - 1 );
  }

 // v added to the frontier F, unless it dominates v
 static void add( TreeFrontier< Cost , K > & F , const Vec & v ,
                  std::uint64_t & n ) {
  if( ! F.dominated( v , n ) )
   F.insert( v , n );
  }

 // the layers cut again in B vectors each, the nodes made anew
 void build( std::uint64_t & n ) {
  std::vector< Vec > all;
  for( const auto & L : v_layer )
   all.insert( all.end() , L.begin() , L.end() );
  v_layer.clear();
  for( std::size_t k = 0 ; k < all.size() ; k += B )
   v_layer.emplace_back( all.begin() + std::ptrdiff_t( k ) ,
                         all.begin() + std::ptrdiff_t(
                          std::min( all.size() , k + B ) ) );
  if( v_layer.empty() )
   v_layer.emplace_back();
  std::size_t cap = 2;
  while( cap < 2 * v_layer.size() + 1 )
   cap *= 2;
  v_node.clear();
  v_node.resize( cap );
  // the node i covers the layers i - lowbit( i ) , ... , i - 1, so that
  // the layer b is in the nodes b + 1 , b + 1 + lowbit( b + 1 ) , ...
  for( std::size_t b = 0 ; b < v_layer.size() ; ++b )
   for( const Vec & z : v_layer[ b ] )
    for( std::size_t i = b + 1 ; i < cap ; i += i & ( ~i + 1 ) )
     add( v_node[ i ] , z , n );
  }

 std::vector< std::vector< Vec > > v_layer;       ///< the layers
 std::vector< TreeFrontier< Cost , K > > v_node;  ///< the Fenwick tree

 };  // end( class( LayerIndex ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BlockIndex -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a set of full vectors with the minima of its blocks, for the dominance
/// checks of vectors in any order
/** The vectors are kept by increasing first component in blocks of at most
 * 2 B of them, each block with the minimum of its vectors component by
 * component, and a segment tree over the blocks keeps the minima of its
 * ranges: a check descends only into the nodes whose minimum weakly
 * dominates the vector, since no vector below the others can, and scans
 * the blocks it reaches, the block-minimum index of Xu et al. (2026) on a
 * segment tree. An insertion compares no vector. The vectors that a later
 * one dominates are kept, as in LayerIndex. dominated( v , n ) counts in n
 * one check per node and per vector compared. */

template< class Cost , std::size_t K >
class BlockIndex {

 public:

 using Vec = std::array< Cost , K >;              ///< a vector of costs

 [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const {
  return( ! v_block.empty() && down( 1 , 0 , f_cap , v , n ) );
  }

 void insert( const Vec & v , std::uint64_t & ) {
  // the last block that starts with a first component not above v_1
  std::size_t b = std::size_t( std::upper_bound(
                   v_block.begin() , v_block.end() , v[ 0 ] ,
                   []( Cost c , const Block & k ) {
                    return( c < k.v.front()[ 0 ] );
                    } ) - v_block.begin() );
  b = b ? b - 1 : 0;
  if( v_block.empty() )
   v_block.emplace_back();
  Block & k = v_block[ b ];
  k.v.insert( std::upper_bound( k.v.begin() , k.v.end() , v ,
                                []( const Vec & a , const Vec & c ) {
                                 return( a[ 0 ] < c[ 0 ] );
                                 } ) , v );
  if( k.v.size() == 1 )
   k.min = v;
  else
   for( std::size_t i = 0 ; i < K ; ++i )
    k.min[ i ] = std::min( k.min[ i ] , v[ i ] );
  if( k.v.size() > 2 * B ) {
   Block r;
   r.v.assign( k.v.begin() + B , k.v.end() );
   k.v.resize( B );
   k.min = minimum( k.v );
   r.min = minimum( r.v );
   v_block.insert( v_block.begin() + std::ptrdiff_t( b ) + 1 ,
                   std::move( r ) );
   build();
   }
  else
   if( v_block.size() > f_cap )
    build();
   else
    for( std::size_t x = f_cap + b ; x >= 1 ; x /= 2 )
     v_tree[ x ] = x >= f_cap ? v_block[ x - f_cap ].min
                              : meet( 2 * x , 2 * x + 1 );
  }

 [[nodiscard]] std::size_t size( void ) const noexcept {
  std::size_t s = 0;
  for( const Block & k : v_block )
   s += k.v.size();
  return( s );
  }

 private:

 static constexpr std::size_t B = 16;             ///< half a block

 /// a block: its vectors by increasing first component, their minimum
 struct Block {
  std::vector< Vec > v;
  Vec min;
  };

 [[nodiscard]] static Vec minimum( const std::vector< Vec > & w ) {
  Vec m = w.front();
  for( const Vec & x : w )
   for( std::size_t i = 0 ; i < K ; ++i )
    m[ i ] = std::min( m[ i ] , x[ i ] );
  return( m );
  }

 // true if the node x covers no block
 [[nodiscard]] bool empty_node( std::size_t x ) const {
  while( x < f_cap )
   x *= 2;
  return( x - f_cap >= v_block.size() );
  }

 // the minimum of the nodes x and y, a node beyond the blocks being empty
 [[nodiscard]] Vec meet( std::size_t x , std::size_t y ) const {
  if( empty_node( x ) )
   return( v_tree[ y ] );
  if( empty_node( y ) )
   return( v_tree[ x ] );
  Vec m;
  for( std::size_t i = 0 ; i < K ; ++i )
   m[ i ] = std::min( v_tree[ x ][ i ] , v_tree[ y ][ i ] );
  return( m );
  }

 // the tree over the blocks, rebuilt
 void build( void ) {
  f_cap = 1;
  while( f_cap < v_block.size() )
   f_cap *= 2;
  v_tree.assign( 2 * f_cap , Vec{} );
  for( std::size_t b = 0 ; b < v_block.size() ; ++b )
   v_tree[ f_cap + b ] = v_block[ b ].min;
  for( std::size_t x = f_cap - 1 ; x >= 1 ; --x )
   v_tree[ x ] = meet( 2 * x , 2 * x + 1 );
  }

 // true if a vector of the blocks under the node x, which covers
 // [ lo , lo + w ), weakly dominates v
 bool down( std::size_t x , std::size_t lo , std::size_t w , const Vec & v ,
            std::uint64_t & n ) const {
  if( lo >= v_block.size() )
   return( false );
  ++n;
  if( ! weakly_dominates( v_tree[ x ] , v ) )
   return( false );
  if( w == 1 ) {
   for( const Vec & z : v_block[ lo ].v ) {
    if( v[ 0 ] < z[ 0 ] )
     return( false );
    ++n;
    if( weakly_dominates( z , v ) )
     return( true );
    }
   return( false );
   }
  return( down( 2 * x , lo , w / 2 , v , n ) ||
          down( 2 * x + 1 , lo + w / 2 , w / 2 , v , n ) );
  }

 std::vector< Block > v_block;                    ///< the blocks
 std::vector< Vec > v_tree = std::vector< Vec >( 2 );  ///< the minima
 std::size_t f_cap = 1;                           ///< the leaves of the tree

 };  // end( class( BlockIndex ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MultiObjective */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MultiObjective.h -----------------------*/
/*--------------------------------------------------------------------------*/
