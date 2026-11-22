/*--------------------------------------------------------------------------*/
/*------------------------ File MaxCommonSubgraph.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MaxCommonSubgraph, the maximum common induced subgraph of
 * two undirected graphs as a problem of the Branch-and-Bound (B&B), with
 * the partition of the vertices into label classes and the bound of the
 * McSplit algorithm of McCreesh, Prosser and Trimble (2017). It is an
 * example of the B&B layer, outside the engine.
 *
 * A node is a matching \f$ M \f$ of vertices of \f$ G \f$ to vertices of
 * \f$ H \f$, which induce isomorphic subgraphs, and a partition of the
 * vertices left into classes \f$ ( L_c , R_c ) \f$, those of \f$ L_c \f$
 * in \f$ G \f$ and of \f$ R_c \f$ in \f$ H \f$ having the same adjacency
 * to the matched vertices, so that only a vertex of \f$ L_c \f$ can be
 * matched to one of \f$ R_c \f$. Its bound is
 * \f$ |M| + \sum_c \min \{ |L_c| , |R_c| \} \f$, and \f$ |M| \f$ is a
 * feasible value. The evaluation picks the class of the least
 * \f$ \max \{ |L_c| , |R_c| \} \f$ and in it the vertex \f$ v \f$ of
 * \f$ G \f$ of the largest degree; the children match \f$ v \f$ to each
 * vertex \f$ w \f$ of \f$ R_c \f$, by decreasing degree, splitting each
 * class into the neighbors and the non-neighbors of \f$ v \f$ and of
 * \f$ w \f$, and the last child leaves \f$ v \f$ unmatched; the estimate of
 * a child is its own bound. The children come one at a time [see
 * BBCursor].
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

#ifndef __PBeFS_MaxCommonSubgraph
 #define __PBeFS_MaxCommonSubgraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstdint>

#include <optional>

#include <utility>

#include <vector>

#include "PBeFS/bnb/MaxClique.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MaxCommonSubgraph -------------------------*/
/*--------------------------------------------------------------------------*/
/// the maximum common induced subgraph of two graphs, by McSplit

struct MaxCommonSubgraph {

 using value_type = int;          ///< the type of the values

 /// a label class: the vertices of G and of H left in it
 struct Class {
  std::vector< std::uint16_t > left;   ///< of G
  std::vector< std::uint16_t > right;  ///< of H
  };

 /// a node: the matching, the classes, and the class and the vertex to
 /// branch on once it is evaluated
 struct node_type {
  std::vector< std::pair< std::uint16_t , std::uint16_t > > match;
  std::vector< Class > classes;        ///< the label classes
  int cls = -1;                        ///< the class to branch on
  int v = -1;                          ///< the vertex of G in it
  };

 static constexpr Sense sense = Sense::maximize;

 BitGraph g;                       ///< the first graph
 BitGraph h;                       ///< the second graph

 node_type root( void ) const {
  node_type r;
  Class c;
  for( std::size_t v = 0 ; v < g.n ; ++v )
   c.left.push_back( std::uint16_t( v ) );
  for( std::size_t w = 0 ; w < h.n ; ++w )
   c.right.push_back( std::uint16_t( w ) );
  if( ! c.left.empty() && ! c.right.empty() )
   r.classes.push_back( std::move( c ) );
  return( r );
  }

 /// the bound of a node
 [[nodiscard]] static int bound_of( const node_type & x ) {
  std::size_t b = x.match.size();
  for( const Class & c : x.classes )
   b += std::min( c.left.size() , c.right.size() );
  return( int( b ) );
  }

 BasicEvaluation< int > evaluate( node_type & x ) const {
  BasicEvaluation< int > e;
  e.bound = bound_of( x );
  e.feasible = int( x.match.size() );
  x.cls = x.v = -1;
  std::size_t best = 0;
  for( std::size_t k = 0 ; k < x.classes.size() ; ++k ) {
   const std::size_t s = std::max( x.classes[ k ].left.size() ,
                                   x.classes[ k ].right.size() );
   if( ( x.cls < 0 ) || ( s < best ) ) {
    x.cls = int( k );
    best = s;
    }
   }
  if( x.cls >= 0 ) {
   const Class & c = x.classes[ std::size_t( x.cls ) ];
   for( const auto v : c.left )
    if( ( x.v < 0 ) ||
        ( g.degree( v ) > g.degree( std::size_t( x.v ) ) ) )
     x.v = v;
   }
  return( e );
  }

 /// the child of x that matches its vertex v to w, or leaves v unmatched
 /// if w is none
 [[nodiscard]] node_type split( const node_type & x ,
                                std::optional< std::size_t > w ) const {
  node_type c;
  c.match = x.match;
  const std::size_t v = std::size_t( x.v );
  if( ! w ) {
   c.classes = x.classes;
   auto & l = c.classes[ std::size_t( x.cls ) ].left;
   l.erase( std::find( l.begin() , l.end() , std::uint16_t( v ) ) );
   if( l.empty() )
    c.classes.erase( c.classes.begin() + x.cls );
   return( c );
   }
  c.match.emplace_back( std::uint16_t( v ) , std::uint16_t( *w ) );
  for( const Class & k : x.classes ) {
   Class a , b;                        // adjacent, non-adjacent
   for( const auto u : k.left )
    if( u != v )
     ( g.edge( v , u ) ? a : b ).left.push_back( u );
   for( const auto u : k.right )
    if( u != *w )
     ( h.edge( *w , u ) ? a : b ).right.push_back( u );
   if( ! a.left.empty() && ! a.right.empty() )
    c.classes.push_back( std::move( a ) );
   if( ! b.left.empty() && ! b.right.empty() )
    c.classes.push_back( std::move( b ) );
   }
  return( c );
  }

 /// the children of an evaluated node, one at a time
 struct Cursor {
  const MaxCommonSubgraph * s;
  const node_type * x;
  std::vector< std::uint16_t > ws;     // the vertices of H, by degree
  std::size_t k = 0;
  bool skip = false;                   // the child without v is done
  std::optional< BBChild< node_type , int > > next( void ) {
   if( x->cls < 0 )
    return( std::nullopt );
   std::optional< std::size_t > w;
   if( k < ws.size() )
    w = ws[ k++ ];
   else
    if( skip )
     return( std::nullopt );
    else
     skip = true;
   node_type c = s->split( *x , w );
   const int b = bound_of( c );
   return( BBChild< node_type , int >{ std::move( c ) , b } );
   }
  };

 Cursor children( const node_type & x ) const {
  Cursor c{ this , & x , {} };
  if( x.cls >= 0 ) {
   c.ws = x.classes[ std::size_t( x.cls ) ].right;
   std::stable_sort( c.ws.begin() , c.ws.end() ,
                     [ & ]( std::uint16_t a , std::uint16_t b ) {
                      return( h.degree( a ) > h.degree( b ) );
                      } );
   }
  return( c );
  }

 /// true if the matching of x induces isomorphic subgraphs
 bool check( const node_type & x ) const {
  for( std::size_t i = 0 ; i < x.match.size() ; ++i )
   for( std::size_t j = i + 1 ; j < x.match.size() ; ++j ) {
    const auto [ a , b ] = x.match[ i ];
    const auto [ c , d ] = x.match[ j ];
    if( ( a == c ) || ( b == d ) ||
        ( g.edge( a , c ) != h.edge( b , d ) ) )
     return( false );
    }
  return( true );
  }

 /// the modular product of g and h, whose largest clique has the size of
 /// the largest common induced subgraph: the vertex v h.n + w is the pair
 /// (v, w), adjacent to (v', w') if v differs from v', w from w', and the
 /// edge (v, v') of g is there exactly when (w, w') of h is
 [[nodiscard]] BitGraph modular_product( void ) const {
  BitGraph p( g.n * h.n );
  for( std::size_t v = 0 ; v < g.n ; ++v )
   for( std::size_t w = 0 ; w < h.n ; ++w )
    for( std::size_t a = v + 1 ; a < g.n ; ++a )
     for( std::size_t b = 0 ; b < h.n ; ++b )
      if( ( b != w ) && ( g.edge( v , a ) == h.edge( w , b ) ) )
       p.add_edge( v * h.n + w , a * h.n + b );
  return( p );
  }

 };  // end( struct( MaxCommonSubgraph ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MaxCommonSubgraph */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MaxCommonSubgraph.h ----------------------*/
/*--------------------------------------------------------------------------*/
