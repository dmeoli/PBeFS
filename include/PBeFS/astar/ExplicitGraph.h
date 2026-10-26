/*--------------------------------------------------------------------------*/
/*-------------------------- File ExplicitGraph.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ExplicitGraph, which implements the
 * Domain concept [see Domain.h] for a weighted graph given by adjacency
 * lists, with a heuristic given node by node, such as the map of Romania
 * of Russell and Norvig.
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

#ifndef __PBeFS_ExplicitGraph
 #define __PBeFS_ExplicitGraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
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
 * hashing and the perfect hash too. */

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
  : f_source( s ) , f_target( t ) , v_adj( n ) , v_h( n , Cost( 0 ) ) {
  if( ( s >= n ) || ( t >= n ) )
   throw( std::invalid_argument( "ExplicitGraph::ExplicitGraph: source or "
                                 "target out of range" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the arc ( u , v ) of cost c, and ( v , u ) too unless directed

 void add_edge( state_type u , state_type v , Cost c ,
                bool directed = false ) {
  if( ( u >= v_adj.size() ) || ( v >= v_adj.size() ) || ( c < 0 ) )
   throw( std::invalid_argument( "ExplicitGraph::add_edge: bad arc" ) );
  v_adj[ u ].emplace_back( v , c );
  if( ! directed )
   v_adj[ v ].emplace_back( u , c );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the heuristic value of the node u

 void set_heuristic( state_type u , Cost h ) { v_h.at( u ) = h; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the source of the search

 void set_source( state_type s ) {
  if( s >= v_adj.size() )
   throw( std::invalid_argument( "ExplicitGraph::set_source: bad node" ) );
  f_source = s;
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_source );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_target );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( mix64( s ) );
  }

 [[nodiscard]] Cost heuristic( state_type s ) const noexcept {
  return( v_h[ s ] );
  }

 template< class F >
 void for_each_successor( state_type s , Cost , F && f ) const {
  for( const auto & [ v , c ] : v_adj[ s ] )
   f( v , c , v_h[ v ] );
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( v_adj.size() );
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
  return( v_adj.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 state_type f_source;                                  ///< the source
 state_type f_target;                                  ///< the target
 std::vector< std::vector< std::pair< state_type , Cost > > > v_adj;
                                                   ///< the adjacency lists
 std::vector< Cost > v_h;                     ///< the heuristic, by node

 };  // end( class( ExplicitGraph ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ExplicitGraph */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ExplicitGraph.h ------------------------*/
/*--------------------------------------------------------------------------*/
