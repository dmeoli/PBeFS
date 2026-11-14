/*--------------------------------------------------------------------------*/
/*---------------------------- File NBlockMap.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NBlockMap, the nblocks of an
 * NBlockDomain built at their first use, with their duplicate detection
 * scope and their interference scope, shared by the searches that
 * partition the state space into nblocks: SafePBNF, PSDD and BFPSDD.
 *
 * The duplicate detection scope of the nblock \f$ b \f$ is \f$ b \f$ and
 * its abstract successors \f$ S( b ) \f$, the only nblocks in which the
 * children of a node of \f$ b \f$ can fall; its interference scope
 * \f$ I( b ) \f$ is made of the nblocks whose duplicate detection scope
 * meets that of \f$ b \f$, i.e., the predecessors of \f$ b \f$ and of its
 * successors, \f$ b \f$ excluded. With an undirected abstract graph, the
 * predecessors are the successors, and
 * \f[
 *   I( b ) = \bigl( S( b ) \cup \textstyle\bigcup_{s \in S( b )} S( s )
 *            \bigr) \setminus \{ b \} ,
 * \f]
 * so that two nblocks can be in use at once if and only if neither is in
 * the interference scope of the other [see Zhou and Hansen, AAAI 2007, and
 * Burns et al., JAIR 2010].
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

#ifndef __PBeFS_NBlockMap
 #define __PBeFS_NBlockMap
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstddef>

#include <memory>

#include <mutex>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS NBlockMap -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the nblocks of an NBlockDomain, built at their first use
/** NB is the nblock of a search, default constructible, with the members
 * std::vector< std::size_t > scope and interferes, which block() fills
 * with the duplicate detection scope (the successors and the nblock
 * itself, last) and with the interference scope of the nblock. Any thread
 * may call find() and block(), the latter serializing the builds; an nblock
 * never moves once built, and lives as long as the map. */

template< NBlockDomain D , class NB >
class NBlockMap {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the map)

 explicit NBlockMap( const D & d )
  : f_d( d ) , v_b( new std::atomic< NB * >[ d.num_nblocks() ] ) {
  for( std::size_t i = 0 ; i < d.num_nblocks() ; ++i )
   v_b[ i ].store( nullptr , std::memory_order_relaxed );
  }

/*-------------------------------- READING ---------------------------------*/
 /// the nblock b if it has been built, nullptr otherwise

 [[nodiscard]] NB * find( std::size_t b ) const {
  return( v_b[ b ].load( std::memory_order_acquire ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the nblock b, built at its first use with its scopes

 NB & block( std::size_t b ) {
  if( NB * x = find( b ) )
   return( *x );
  std::lock_guard< std::mutex > guard( f_build );
  if( NB * x = v_b[ b ].load( std::memory_order_relaxed ) )
   return( *x );
  auto n = std::make_unique< NB >();
  std::vector< std::size_t > succ , iset;
  f_d.for_each_nblock_successor( b , [ & ]( std::size_t s ) {
   succ.push_back( s );
   } );
  for( std::size_t s : succ ) {
   iset.push_back( s );
   f_d.for_each_nblock_successor( s , [ & ]( std::size_t t ) {
    iset.push_back( t );
    } );
   }
  std::sort( iset.begin() , iset.end() );
  iset.erase( std::unique( iset.begin() , iset.end() ) , iset.end() );
  iset.erase( std::remove( iset.begin() , iset.end() , b ) , iset.end() );
  if constexpr( requires { n->id = b; } )
   n->id = b;
  n->interferes = std::move( iset );
  n->scope = std::move( succ );
  n->scope.push_back( b );
  v_b[ b ].store( n.get() , std::memory_order_release );
  NB & ref = *n;
  v_own.push_back( std::move( n ) );
  return( ref );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( n ) on every nblock n built so far, by increasing index

 template< class F >
 void for_each( F && f ) const {
  for( std::size_t i = 0 ; i < f_d.num_nblocks() ; ++i )
   if( NB * x = find( i ) )
    f( *x );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const D & f_d;                                   ///< the Domain
 std::unique_ptr< std::atomic< NB * >[] > v_b;    ///< the nblocks
 std::vector< std::unique_ptr< NB > > v_own;      ///< their ownership
 std::mutex f_build;                              ///< serializes block()

 };  // end( class( NBlockMap ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_NBlockMap */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NBlockMap.h --------------------------*/
/*--------------------------------------------------------------------------*/
