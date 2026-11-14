/*--------------------------------------------------------------------------*/
/*----------------------- File SubgraphIsomorphism.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for SubgraphIsomorphism, the (non-induced) subgraph
 * isomorphism problem as a problem of the Branch-and-Bound (B&B) in the
 * modes decide and enumerate [see BBMode]: whether a pattern graph is a
 * subgraph of a target graph, and in how many ways. It is an example of the
 * B&B layer, outside the engine.
 *
 * A node gives each vertex \f$ p \f$ of the pattern a domain
 * \f$ D_p \f$ of target vertices, a bitset, the assigned vertices having a
 * single value. At the root \f$ D_p \f$ holds the target vertices of degree
 * at least that of \f$ p \f$. The evaluation finds the node infeasible if
 * a domain is empty or if the union of the domains of the unassigned
 * vertices is smaller than their number (they must take distinct values),
 * feasible (of value 0) if all the domains are single values that make an
 * embedding (and infeasible if they do not), and otherwise
 * picks the unassigned vertex of the smallest domain (of the largest
 * degree on ties); its children give it each value of its domain, by
 * increasing number, removing the value from the other domains and
 * restricting the domains of its neighbours to the neighbours of the value
 * (forward checking). The children come one at a time [see BBCursor]. An
 * embedding maps the edges of the pattern to edges of the target; with
 * induced, also its non-edges to non-edges, the non-neighbours of a vertex
 * being then restricted to the non-neighbours of its value.
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

#ifndef __PBeFS_SubgraphIsomorphism
 #define __PBeFS_SubgraphIsomorphism
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstdint>

#include <fstream>

#include <optional>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/bnb/MaxClique.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SubgraphIsomorphism ------------------------*/
/*--------------------------------------------------------------------------*/
/// the subgraph isomorphism of a pattern graph into a target graph

struct SubgraphIsomorphism {

 using value_type = int;          ///< the type of the values

 /// a node: the domains of the pattern vertices, and the vertex to branch
 /// on once it is evaluated
 struct node_type {
  std::vector< BitSet > dom;        ///< the domains
  int var = -1;                     ///< the vertex to branch on
  };

 static constexpr Sense sense = Sense::minimize;

 BitGraph pattern;                  ///< the pattern
 BitGraph target;                   ///< the target
 bool induced = false;              ///< non-edges to non-edges too

 /// the graph of a file of the format of LAD: the number of vertices, then
 /// for each vertex its degree and its neighbours (from 0); an edge given
 /// in one direction is taken in both; throws std::runtime_error if the
 /// file cannot be read
 static BitGraph lad( const std::string & file ) {
  std::ifstream in( file );
  std::size_t n;
  if( ! ( in >> n ) )
   throw( std::runtime_error( "SubgraphIsomorphism::lad: cannot read " +
                              file ) );
  BitGraph g( n );
  for( std::size_t v = 0 ; v < n ; ++v ) {
   std::size_t d , u;
   in >> d;
   for( std::size_t k = 0 ; k < d ; ++k ) {
    in >> u;
    if( ! in || ( u >= n ) )
     throw( std::runtime_error( "SubgraphIsomorphism::lad: bad file " +
                                file ) );
    g.add_edge( v , u );
    }
   }
  return( g );
  }

 node_type root( void ) const {
  node_type r;
  r.dom.assign( pattern.n , BitSet( target.n ) );
  for( std::size_t p = 0 ; p < pattern.n ; ++p ) {
   const std::size_t dp = pattern.degree( p );
   for( std::size_t t = 0 ; t < target.n ; ++t )
    if( target.degree( t ) >= dp )
     r.dom[ p ].set( t );
   }
  return( r );
  }

 BasicEvaluation< int > evaluate( node_type & x ) const {
  BasicEvaluation< int > e;
  BitSet all( target.n );
  std::size_t open = 0 , best = 0;
  x.var = -1;
  for( std::size_t p = 0 ; p < pattern.n ; ++p ) {
   const std::size_t c = x.dom[ p ].count();
   if( c == 0 ) {
    e.infeasible = true;
    return( e );
    }
   if( c == 1 )
    continue;
   ++open;
   for( std::size_t j = 0 ; j < all.w.size() ; ++j )
    all.w[ j ] |= x.dom[ p ].w[ j ];
   if( ( x.var < 0 ) || ( c < best ) ||
       ( ( c == best ) && ( pattern.degree( p ) >
                            pattern.degree( std::size_t( x.var ) ) ) ) ) {
    x.var = int( p );
    best = c;
    }
   }
  if( all.count() < open ) {
   e.infeasible = true;
   return( e );
   }
  if( x.var < 0 ) {                    // all single: an embedding?
   if( ! check( x ) ) {
    e.infeasible = true;
    return( e );
    }
   e.feasible = 0;
   }
  return( e );
  }

 /// the child of x that gives the value t to its vertex
 [[nodiscard]] node_type assign( const node_type & x ,
                                 std::size_t t ) const {
  node_type c;
  c.dom = x.dom;
  const std::size_t v = std::size_t( x.var );
  c.dom[ v ] = BitSet( target.n );
  c.dom[ v ].set( t );
  for( std::size_t q = 0 ; q < pattern.n ; ++q ) {
   if( q == v )
    continue;
   c.dom[ q ].reset( t );
   if( pattern.edge( v , q ) )
    c.dom[ q ] &= target.adj[ t ];
   else
    if( induced )
     c.dom[ q ].remove( target.adj[ t ] );
   }
  return( c );
  }

 /// the children of an evaluated node, one at a time
 struct Cursor {
  const SubgraphIsomorphism * s;
  const node_type * x;
  BitSet left;
  std::optional< BBChild< node_type , int > > next( void ) {
   if( ( x->var < 0 ) || left.empty() )
    return( std::nullopt );
   const std::size_t t = left.first();
   left.reset( t );
   return( BBChild< node_type , int >{ s->assign( *x , t ) , 0 } );
   }
  };

 Cursor children( const node_type & x ) const {
  return( Cursor{ this , & x , x.var < 0 ? BitSet( target.n )
                                         : x.dom[ std::size_t( x.var ) ] } );
  }

 /// true if the domains of a feasible node are an embedding
 bool check( const node_type & x ) const {
  std::vector< std::size_t > f( pattern.n );
  for( std::size_t p = 0 ; p < pattern.n ; ++p ) {
   if( x.dom[ p ].count() != 1 )
    return( false );
   f[ p ] = x.dom[ p ].first();
   }
  for( std::size_t p = 0 ; p < pattern.n ; ++p )
   for( std::size_t q = p + 1 ; q < pattern.n ; ++q ) {
    if( f[ p ] == f[ q ] )
     return( false );
    if( pattern.edge( p , q ) && ! target.edge( f[ p ] , f[ q ] ) )
     return( false );
    if( induced && ! pattern.edge( p , q ) && target.edge( f[ p ] , f[ q ] ) )
     return( false );
    }
  return( true );
  }

 };  // end( struct( SubgraphIsomorphism ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SubgraphIsomorphism */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SubgraphIsomorphism.h ---------------------*/
/*--------------------------------------------------------------------------*/
