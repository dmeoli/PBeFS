/*--------------------------------------------------------------------------*/
/*---------------------------- File OpenList.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the priority queues of the OPEN list: the class templates
 * HeapOpenList, a 4-ary heap for any cost type, and BucketOpenList, a
 * two-level bucket queue for integer costs.
 *
 * Both queues order the entries by increasing f and break the ties in
 * favour of the larger g, which reaches the goal earlier inside the last f
 * layer [see Burns et al., SoCS 2012]. Neither supports decrease-key: a node
 * reached through a cheaper path is pushed again, and an entry is stale when
 * the node is closed or its g no longer matches the g of the entry. The
 * caller discards the stale entries as they surface.
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

#ifndef __PBeFS_OpenList
 #define __PBeFS_OpenList
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/Structures.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the kind of OPEN list a search uses
enum class OpenKind {
 automatic ,  ///< buckets for integer costs, heap otherwise
 heap ,       ///< HeapOpenList
 buckets      ///< BucketOpenList (integer costs only)
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// an entry of the OPEN list: the node and the g it had when pushed

template< class NodeT >
struct OpenEntry {

 NodeT * node;               ///< the node
 decltype( NodeT::g ) g;     ///< the g of the node when it was pushed

 /// true if the node has been closed or improved since it was pushed
 [[nodiscard]] bool stale( void ) const noexcept {
  return( node->closed || ( node->g != g ) );
  }

 };  // end( struct( OpenEntry ) )

#ifdef PBEFS_OWN_STRUCTURES

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS HeapOpenList ---------------------------*/
/*--------------------------------------------------------------------------*/
/// an Arity-ary min-heap of nodes
/** A d-ary heap with d = 4 has half the depth of a binary heap, and the
 * children of a node lie in one cache line. The items keep f and g next to
 * the pointer, so that comparisons never touch the nodes. */

template< class NodeT , unsigned Arity = 4 >
class HeapOpenList {

 using Cost = decltype( NodeT::g );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// pushes a node with its current f and g

 void push( NodeT * node ) { push( node , node->g ); }

 /// pushes a node with the given g (and \f$ f = g + h \f$), not reading
 /// node->g

 void push( NodeT * node , Cost g ) {
  v_heap.push_back( { g + node->h , g , node } );
  sift_up( v_heap.size() - 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the entry with the smallest f (and the largest g among those)

 [[nodiscard]] OpenEntry< NodeT > top( void ) const noexcept {
  return( OpenEntry< NodeT >{ v_heap.front().node , v_heap.front().g } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the f of top()

 [[nodiscard]] Cost top_f( void ) const noexcept {
  return( v_heap.front().f );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes and returns top()

 OpenEntry< NodeT > pop( void ) {
  OpenEntry< NodeT > e{ v_heap.front().node , v_heap.front().g };
  v_heap.front() = v_heap.back();
  v_heap.pop_back();
  if( ! v_heap.empty() )
   sift_down( 0 );
  return( e );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept {
  return( v_heap.empty() );
  }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_heap.size() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// an item of the heap
 struct Item {
  Cost f;
  Cost g;
  NodeT * node;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the order of the heap: smaller f first, then larger g

 static bool better( const Item & a , const Item & b ) noexcept {
  return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.g > b.g ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void sift_up( std::size_t i ) {
  Item x = v_heap[ i ];
  while( i > 0 ) {
   std::size_t p = ( i - 1 ) / Arity;
   if( ! better( x , v_heap[ p ] ) )
    break;
   v_heap[ i ] = v_heap[ p ];
   i = p;
   }
  v_heap[ i ] = x;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void sift_down( std::size_t i ) {
  Item x = v_heap[ i ];
  const std::size_t n = v_heap.size();
  for( ;; ) {
   std::size_t c = i * Arity + 1;
   if( c >= n )
    break;
   std::size_t best = c;
   const std::size_t last = c + Arity < n ? c + Arity : n;
   for( std::size_t k = c + 1 ; k < last ; ++k )
    if( better( v_heap[ k ] , v_heap[ best ] ) )
     best = k;
   if( ! better( v_heap[ best ] , x ) )
    break;
   v_heap[ i ] = v_heap[ best ];
   i = best;
   }
  v_heap[ i ] = x;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 std::vector< Item > v_heap;                              ///< the heap

 };  // end( class( HeapOpenList ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BucketOpenList --------------------------*/
/*--------------------------------------------------------------------------*/
/// a two-level bucket queue for integer costs
/** The outer level is an array indexed by f, the inner one an array of LIFO
 * lists indexed by g, scanned from the largest g: all operations take
 * \f$ O( 1 ) \f$ amortized time, besides the empty buckets skipped by a
 * pop [see Dial, CACM 1969, and the page "Best-first search and A*" of the
 * documentation]. This is the nested bucket queue of
 * Burns et al. (SoCS 2012), with which A* on the 15-puzzle runs 2.4 times
 * faster than with a binary heap. The lists of an f layer are released when
 * the layer empties, so the memory follows the band of the f values that
 * are in OPEN. */

template< class NodeT >
class BucketOpenList {

 using Cost = decltype( NodeT::g );

 static_assert( std::is_integral_v< Cost > ,
                "BucketOpenList requires an integer cost type" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 BucketOpenList( void ) : f_min_f( 0 ) , f_size( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// pushes a node with its current f and g

 void push( NodeT * node ) { push( node , node->g ); }

 /// pushes a node with the given g (and \f$ f = g + h \f$), not reading
 /// node->g

 void push( NodeT * node , Cost gn ) {
  const auto f = static_cast< std::size_t >( gn + node->h );
  const auto g = static_cast< std::size_t >( gn );
  if( f >= v_layers.size() )
   v_layers.resize( f + 1 );
  Layer & l = v_layers[ f ];
  if( g >= l.by_g.size() )
   l.by_g.resize( g + 1 );
  l.by_g[ g ].push_back( node );
  if( ( l.count++ == 0 ) || ( g > l.top_g ) )
   l.top_g = g;
  if( ( f_size++ == 0 ) || ( f < f_min_f ) )
   f_min_f = f;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the entry with the smallest f (and the largest g among those)

 [[nodiscard]] OpenEntry< NodeT > top( void ) {
  Layer & l = front();
  return( OpenEntry< NodeT >{ l.by_g[ l.top_g ].back() ,
                              static_cast< Cost >( l.top_g ) } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the f of top()

 [[nodiscard]] Cost top_f( void ) {
  front();
  return( static_cast< Cost >( f_min_f ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes and returns top()

 OpenEntry< NodeT > pop( void ) {
  Layer & l = front();
  OpenEntry< NodeT > e{ l.by_g[ l.top_g ].back() ,
                        static_cast< Cost >( l.top_g ) };
  l.by_g[ l.top_g ].pop_back();
  --f_size;
  if( --l.count == 0 ) {
   std::vector< std::vector< NodeT * > >().swap( l.by_g );
   l.top_g = 0;
   }
  return( e );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if there is no entry

 [[nodiscard]] bool empty( void ) const noexcept { return( f_size == 0 ); }

 /// the number of entries, the stale ones included

 [[nodiscard]] std::size_t size( void ) const noexcept { return( f_size ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// the nodes with the same f, by g
 struct Layer {
  std::vector< std::vector< NodeT * > > by_g;
  std::size_t top_g = 0;
  std::size_t count = 0;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the first non-empty list of the first non-empty layer

 Layer & front( void ) {
  while( v_layers[ f_min_f ].count == 0 )
   ++f_min_f;
  Layer & l = v_layers[ f_min_f ];
  while( l.by_g[ l.top_g ].empty() )
   --l.top_g;
  return( l );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 std::vector< Layer > v_layers;                    ///< the layers, by f
 std::size_t f_min_f;              ///< no non-empty layer comes before
 std::size_t f_size;                 ///< the number of entries in OPEN

 };  // end( class( BucketOpenList ) )

#endif  /* PBEFS_OWN_STRUCTURES */

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

// the structures of another set, in place of those of the library
#if defined( PBEFS_STRUCTURES_STD )
 #include "PBeFS/astar/structures/std/OpenList.h"
#elif defined( PBEFS_STRUCTURES_BOOST )
 #include "PBeFS/astar/structures/boost/OpenList.h"
#elif defined( PBEFS_STRUCTURES_ADS )
 #include "PBeFS/astar/structures/ads/OpenList.h"
#endif

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_OpenList */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File OpenList.h ---------------------------*/
/*--------------------------------------------------------------------------*/
