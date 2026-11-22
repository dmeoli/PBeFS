/*--------------------------------------------------------------------------*/
/*---------------------------- File MaxClique.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MaxClique, the maximum clique of an undirected graph as a
 * problem of the Branch-and-Bound (B&B), with the bound of a greedy
 * coloring computed on bitsets as in the algorithm of Tomita and Kameda
 * (2007) and in BBMC of San Segundo et al. (2011), and a reader of the
 * graphs of the second DIMACS challenge. It is an example of the B&B
 * layer, outside the engine.
 *
 * The vertices are numbered by decreasing degree (the lower number first on
 * ties), and a node is a clique \f$ C \f$ with the set \f$ P \f$ of the
 * vertices adjacent to all those of \f$ C \f$. Its evaluation colors
 * \f$ P \f$ greedily, a color class at a time, each taking the vertices of
 * \f$ P \f$ in their order that are not adjacent to those already in the
 * class; with \f$ k \f$ colors, \f$ |C| + k \f$ bounds the cliques of the
 * node, and \f$ |C| \f$ is a feasible value. The children are made in the
 * reverse order of the coloring: the child of the vertex \f$ v \f$ of
 * color \f$ c \f$ adds \f$ v \f$ to \f$ C \f$, restricts \f$ P \f$ to the
 * neighbors of \f$ v \f$ that come before it, and has the estimate
 * \f$ |C| + c \f$, so that the estimates do not increase and a search with
 * BBOptions::prune_level stops at the first child pruned. The children come
 * one at a time [see BBCursor], built only when the search asks for them.
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

#ifndef __PBeFS_MaxClique
 #define __PBeFS_MaxClique
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cstdint>

#include <fstream>

#include <numeric>

#include <optional>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/bnb/BranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS BitSet ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a set of the integers below a given size, as words of 64 bits

struct BitSet {

 BitSet( void ) = default;

 explicit BitSet( std::size_t n ) : w( ( n + 63 ) / 64 , 0 ) {}

 void set( std::size_t i ) { w[ i / 64 ] |= std::uint64_t( 1 ) << ( i % 64 ); }

 void reset( std::size_t i ) {
  w[ i / 64 ] &= ~( std::uint64_t( 1 ) << ( i % 64 ) );
  }

 [[nodiscard]] bool test( std::size_t i ) const {
  return( ( w[ i / 64 ] >> ( i % 64 ) ) & 1 );
  }

 [[nodiscard]] bool empty( void ) const {
  for( const auto x : w )
   if( x )
    return( false );
  return( true );
  }

 [[nodiscard]] std::size_t count( void ) const {
  std::size_t c = 0;
  for( const auto x : w )
   c += std::size_t( std::popcount( x ) );
  return( c );
  }

 /// the least element, the size of the universe if none
 [[nodiscard]] std::size_t first( void ) const {
  for( std::size_t j = 0 ; j < w.size() ; ++j )
   if( w[ j ] )
    return( 64 * j + std::size_t( std::countr_zero( w[ j ] ) ) );
  return( 64 * w.size() );
  }

 BitSet & operator&=( const BitSet & o ) {
  for( std::size_t j = 0 ; j < w.size() ; ++j )
   w[ j ] &= o.w[ j ];
  return( *this );
  }

 /// removes the elements of o
 void remove( const BitSet & o ) {
  for( std::size_t j = 0 ; j < w.size() ; ++j )
   w[ j ] &= ~o.w[ j ];
  }

 bool operator==( const BitSet & ) const = default;

 std::vector< std::uint64_t > w;      ///< the words

 };  // end( struct( BitSet ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BitGraph ------------------------------*/
/*--------------------------------------------------------------------------*/
/// an undirected graph without loops, as the bitsets of the neighbors

struct BitGraph {

 std::size_t n = 0;                    ///< the vertices
 std::vector< BitSet > adj;            ///< the neighbors of each one

 explicit BitGraph( std::size_t nv = 0 )
  : n( nv ) , adj( nv , BitSet( nv ) ) {}

 void add_edge( std::size_t u , std::size_t v ) {
  if( u == v )
   return;
  adj[ u ].set( v );
  adj[ v ].set( u );
  }

 [[nodiscard]] bool edge( std::size_t u , std::size_t v ) const {
  return( adj[ u ].test( v ) );
  }

 [[nodiscard]] std::size_t degree( std::size_t v ) const {
  return( adj[ v ].count() );
  }

 /// the graph of a file of the DIMACS challenges: lines "p edge n m" (or
 /// "p col n m") and "e u v", the vertices counted from 1; throws
 /// std::runtime_error if it cannot be read
 static BitGraph dimacs( const std::string & file ) {
  std::ifstream in( file );
  if( ! in )
   throw( std::runtime_error( "BitGraph::dimacs: cannot read " + file ) );
  BitGraph g;
  std::string line;
  while( std::getline( in , line ) ) {
   std::istringstream s( line );
   std::string t;
   s >> t;
   if( t == "p" ) {
    std::string kind;
    std::size_t nv , m;
    s >> kind >> nv >> m;
    g = BitGraph( nv );
    }
   else
    if( t == "e" ) {
     std::size_t u , v;
     s >> u >> v;
     if( ( u < 1 ) || ( v < 1 ) || ( u > g.n ) || ( v > g.n ) )
      throw( std::runtime_error( "BitGraph::dimacs: bad edge in " + file ) );
     g.add_edge( u - 1 , v - 1 );
     }
   }
  return( g );
  }

 /// the graph with the vertices renumbered by decreasing degree, the lower
 /// number first on ties; order[ i ] is the old number of the vertex i
 BitGraph by_degree( std::vector< std::size_t > & order ) const {
  order.resize( n );
  std::iota( order.begin() , order.end() , std::size_t( 0 ) );
  std::vector< std::size_t > d( n );
  for( std::size_t v = 0 ; v < n ; ++v )
   d[ v ] = degree( v );
  std::stable_sort( order.begin() , order.end() ,
                    [ & ]( std::size_t a , std::size_t b ) {
                     return( d[ a ] > d[ b ] );
                     } );
  BitGraph h( n );
  for( std::size_t i = 0 ; i < n ; ++i )
   for( std::size_t j = i + 1 ; j < n ; ++j )
    if( edge( order[ i ] , order[ j ] ) )
     h.add_edge( i , j );
  return( h );
  }

 };  // end( struct( BitGraph ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS MaxClique -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the maximum clique of a graph, by the bound of a greedy coloring

struct MaxClique {

 using value_type = int;          ///< the type of the values

 /// a clique, the vertices that extend it, and their coloring once the
 /// node is evaluated
 struct node_type {
  std::vector< std::uint32_t > clique;   ///< the vertices of C
  BitSet p;                              ///< the candidates P
  std::vector< std::uint32_t > order;    ///< P in the order of coloring
  std::vector< std::uint32_t > color;    ///< ... and their colors
  };

 static constexpr Sense sense = Sense::maximize;

 std::vector< std::size_t > label; ///< the number of each vertex in the file
 BitGraph g;                       ///< the graph, by decreasing degree

 MaxClique( void ) = default;

 /// the problem of the graph h, renumbered by decreasing degree
 explicit MaxClique( const BitGraph & h ) : g( h.by_degree( label ) ) {}

 node_type root( void ) const {
  node_type n;
  n.p = BitSet( g.n );
  for( std::size_t v = 0 ; v < g.n ; ++v )
   n.p.set( v );
  return( n );
  }

 /// colors P greedily, the bound being the size of the clique plus the
 /// colors, its size a feasible value
 BasicEvaluation< int > evaluate( node_type & n ) const {
  n.order.clear();
  n.color.clear();
  BitSet left = n.p;
  std::uint32_t c = 0;
  while( ! left.empty() ) {
   ++c;
   BitSet q = left;
   while( ! q.empty() ) {
    const std::size_t v = q.first();
    left.reset( v );
    q.reset( v );
    q.remove( g.adj[ v ] );
    n.order.push_back( std::uint32_t( v ) );
    n.color.push_back( c );
    }
   }
  BasicEvaluation< int > e;
  e.bound = int( n.clique.size() + c );
  e.feasible = int( n.clique.size() );
  return( e );
  }

 /// the children of an evaluated node, the last colored first; the cursor
 /// refers to the node, which must outlive it
 struct Cursor {
  const MaxClique * mc;
  const node_type * n;
  BitSet p;
  std::size_t k;
  std::optional< BBChild< node_type , int > > next( void ) {
   if( k == 0 )
    return( std::nullopt );
   --k;
   const std::uint32_t v = n->order[ k ];
   node_type c;
   c.clique = n->clique;
   c.clique.push_back( v );
   c.p = p;
   c.p &= mc->g.adj[ v ];
   p.reset( v );
   return( BBChild< node_type , int >{ std::move( c ) ,
                                       int( n->clique.size() +
                                            n->color[ k ] ) } );
   }
  };

 Cursor children( const node_type & n ) const {
  return( Cursor{ this , & n , n.p , n.order.size() } );
  }

 /// the clique of a node, numbered as in the file (from 0)
 std::vector< std::size_t > vertices( const node_type & n ) const {
  std::vector< std::size_t > v;
  for( const auto x : n.clique )
   v.push_back( label.empty() ? x : label[ x ] );
  return( v );
  }

 /// true if the vertices of n make a clique
 bool check( const node_type & n ) const {
  for( std::size_t i = 0 ; i < n.clique.size() ; ++i )
   for( std::size_t j = i + 1 ; j < n.clique.size() ; ++j )
    if( ! g.edge( n.clique[ i ] , n.clique[ j ] ) )
     return( false );
  return( true );
  }

 };  // end( struct( MaxClique ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MaxClique */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File MaxClique.h --------------------------*/
/*--------------------------------------------------------------------------*/
