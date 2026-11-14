/*--------------------------------------------------------------------------*/
/*-------------------------- File ExplicitGraph.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ExplicitAdjacency, the arcs of a
 * weighted graph given one by one, and for the class template
 * ExplicitGraph, which implements the Domain concept [see Domain.h] for
 * such a graph, with a heuristic given node by node, such as the map of
 * Romania of Russell and Norvig. The adjacency is that of the set of data
 * structures [see Structures.h]: lists of the library for the sets own and
 * std, the graph in compressed rows of Boost for the set boost, and the
 * GraphBuilder of ADS, frozen into its ArrayGraph, for the set ads.
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

#ifndef __PBeFS_ExplicitGraph
 #define __PBeFS_ExplicitGraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <limits>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

#if defined( PBEFS_STRUCTURES_STD )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ExplicitAdjacency -------------------------*/
/*--------------------------------------------------------------------------*/
/// the arcs of a weighted graph, in a list per node
/** The nodes are \f$ 0 , \ldots , n - 1 \f$; the arcs leaving a node are
 * visited in the order in which they were added. */

template< class Cost >
class ExplicitAdjacency {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// constructor of n nodes and no arc

 explicit ExplicitAdjacency( std::size_t n ) : v_adj( n ) {}

 /// adds the arc ( u , v ) of cost c

 void add( std::uint32_t u , std::uint32_t v , Cost c ) {
  v_adj[ u ].emplace_back( v , c );
  }

 /// calls f( v , c ) for each arc ( u , v ) of cost c

 template< class F >
 void for_each( std::uint32_t u , F && f ) const {
  for( const auto & [ v , c ] : v_adj[ u ] )
   f( v , c );
  }

 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_adj.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 std::vector< std::vector< std::pair< std::uint32_t , Cost > > > v_adj;
                                                   ///< the adjacency lists

 };  // end( class( ExplicitAdjacency ) )

#endif  /* PBEFS_STRUCTURES_STD */

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*------------------ THE ADJACENCY OF THE OTHER STRUCTURES -----------------*/
/*--------------------------------------------------------------------------*/

#if defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/ExplicitAdjacency.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/ExplicitAdjacency.h"
#endif

/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ExplicitGraph ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a weighted graph given explicitly, searched from a source to a target
/** The nodes are \f$ 0 , \ldots , n - 1 \f$; add_edge() adds an arc (or
 * an edge) with its cost, and set_heuristic() gives \f$ h \f$ node by node
 * (0 by default). The states are the nodes, and the features of Zobrist
 * hashing and the perfect hash too. The arcs are kept in an
 * ExplicitAdjacency, whose order of the arcs leaving a node is that of
 * their insertion in every set but boost; in the set ads, a loop ( u , u )
 * is not kept, as it never shortens a path. The graph is a
 * BidirectionalDomain: the arcs are kept reversed too, for the
 * predecessors [see PredecessorDomain], and set_backward_heuristic() gives
 * the heuristic of the backward search, an estimate of the cost from the
 * source (0 by default). */

template< class Cost = int >
class ExplicitGraph {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint32_t;          ///< a node
 using cost_type = Cost;                    ///< the cost type

/*--------------------- CONSTRUCTOR AND INITIALIZATION ---------------------*/

 /// constructor of a graph with n nodes and no arc, searched from s to t

 ExplicitGraph( std::size_t n , state_type s , state_type t )
  : f_source( s ) , f_target( t ) , f_adj( n ) , f_radj( n ) ,
    v_h( n , Cost( 0 ) ) , v_hb( n , Cost( 0 ) ) {
  if( ( s >= n ) || ( t >= n ) )
   throw( std::invalid_argument( "ExplicitGraph::ExplicitGraph: source or "
                                 "target out of range" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the arc ( u , v ) of cost c, and ( v , u ) too unless directed

 void add_edge( state_type u , state_type v , Cost c ,
                bool directed = false ) {
  if( ( u >= f_adj.size() ) || ( v >= f_adj.size() ) || ( c < 0 ) )
   throw( std::invalid_argument( "ExplicitGraph::add_edge: bad arc" ) );
  f_adj.add( u , v , c );
  f_radj.add( v , u , c );
  if( ! directed ) {
   f_adj.add( v , u , c );
   f_radj.add( u , v , c );
   }
  if( ( u != v ) && ( c < f_min ) )
   f_min = c;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the heuristic value of the node u

 void set_heuristic( state_type u , Cost h ) { v_h.at( u ) = h; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the backward heuristic value of the node u, an estimate of the
 /// cost of a path from the source to u (0 by default)

 void set_backward_heuristic( state_type u , Cost h ) { v_hb.at( u ) = h; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the source of the search

 void set_source( state_type s ) {
  if( s >= f_adj.size() )
   throw( std::invalid_argument( "ExplicitGraph::set_source: bad node" ) );
  f_source = s;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the target of the search, whose heuristic, if any, has to be set
 /// again

 void set_target( state_type t ) {
  if( t >= f_adj.size() )
   throw( std::invalid_argument( "ExplicitGraph::set_target: bad node" ) );
  f_target = t;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the target of the search

 [[nodiscard]] state_type target( void ) const noexcept {
  return( f_target );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_source );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_target );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 [[nodiscard]] Cost heuristic( state_type s ) const noexcept {
  return( v_h[ s ] );
  }

 template< class F >
 void for_each_successor( state_type s , Cost , F && f ) const {
  f_adj.for_each( s , [ & ]( state_type v , Cost c ) {
   f( v , c , v_h[ v ] );
   } );
  }

 /// calls f( u , c ) for each arc ( u , s ) of cost c

 template< class F >
 void for_each_predecessor( state_type s , F && f ) const {
  f_radj.for_each( s , std::forward< F >( f ) );
  }

/*--------------------- METHODS OF THE BACKWARD SEARCH ---------------------*/
 [[nodiscard]] state_type goal_state( void ) const noexcept {
  return( f_target );
  }

 [[nodiscard]] Cost backward_heuristic( state_type s ) const noexcept {
  return( v_hb[ s ] );
  }

 template< class F >
 void for_each_predecessor( state_type s , Cost , F && f ) const {
  for_each_predecessor( s , [ & ]( state_type u , Cost c ) {
   f( u , c , v_hb[ u ] );
   } );
  }

 /// the least cost of an arc that is not a loop (the largest value of Cost
 /// while there is none)
 [[nodiscard]] Cost min_edge_cost( void ) const noexcept { return( f_min ); }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( f_adj.size() );
  }

 template< class F >
 void for_each_feature( state_type s , F && f ) const {
  f( std::size_t( s ) );
  }

/*------------------------------ OTHER HASHES ------------------------------*/
 /// the perfect hash: the node

 [[nodiscard]] std::uint64_t perfect_hash( state_type s ) const noexcept {
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_adj.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 state_type f_source;                                  ///< the source
 state_type f_target;                                  ///< the target
 ExplicitAdjacency< Cost > f_adj;                     ///< the arcs
 ExplicitAdjacency< Cost > f_radj;                    ///< the arcs reversed
 std::vector< Cost > v_h;                     ///< the heuristic, by node
 std::vector< Cost > v_hb;           ///< the backward heuristic, by node
 Cost f_min = std::numeric_limits< Cost >::max();  ///< the cheapest arc

 };  // end( class( ExplicitGraph ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ExplicitGraph */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ExplicitGraph.h ------------------------*/
/*--------------------------------------------------------------------------*/
