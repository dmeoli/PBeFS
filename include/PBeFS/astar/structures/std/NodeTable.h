/*--------------------------------------------------------------------------*/
/*---------------------------- File NodeTable.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NodeTable of the set of data structures
 * of the standard library [see Structures.h], the closed list of A*, on
 * std::unordered_map, a table of chains; NodeTable.h includes it in place of
 * its own.
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

#ifndef __PBeFS_std_NodeTable
 #define __PBeFS_std_NodeTable
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <unordered_map>

#include <utility>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS NodeTable -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the closed list, on std::unordered_map, a table of chains
/** The interface is the one described in NodeTable.h. */
template< class NodeT >
class NodeTable {

 using State = decltype( NodeT::state );

 public:

 /// the number of slots that holds n nodes without growing

 [[nodiscard]] static std::size_t slots_for( std::size_t n ) noexcept {
  return( n / 7 * 10 + 16 );
  }

 /// constructor, taking the initial number of slots (reserved for 7 / 10
 /// of them, as many nodes as the table of the library would hold)

 explicit NodeTable( std::size_t capacity = 1024 ) {
  f_map.reserve( capacity / 10 * 7 );
  }

 /// looks up a state, inserting the node made by make() if it is absent

 template< class Make >
 std::pair< NodeT * , bool > find_or_insert( const State & state ,
                                             std::uint64_t hash ,
                                             Make && make ) {
  auto [ it , fresh ] = f_map.try_emplace( Key{ state , hash } , nullptr );
  if( fresh )
   it->second = make();
  return( std::pair( it->second , fresh ) );
  }

 /// looks up a state, returning nullptr if it has never been generated

 [[nodiscard]] NodeT * find( const State & state , std::uint64_t hash ) const
 {
  const auto it = f_map.find( Key{ state , hash } );
  return( it == f_map.end() ? nullptr : it->second );
  }

 /// nothing: std::unordered_map has no way to prefetch a bucket

 void prefetch( std::uint64_t ) const noexcept {}

 /// the number of stored nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_map.size() );
  }

 /// an estimate of the memory of the table: the buckets, and a list node
 /// (key, pointer, next pointer, cached hash) per entry

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_map.bucket_count() * sizeof( void * ) +
          f_map.size() * ( sizeof( Key ) + 3 * sizeof( void * ) ) );
  }

 private:

 /// a key: the state and its hash, computed once by the search
 struct Key {
  State s;
  std::uint64_t h;
  bool operator==( const Key & o ) const { return( s == o.s ); }
  };

 /// the hash of a key: the one computed by the search
 struct Hash {
  std::size_t operator()( const Key & k ) const noexcept {
   return( std::size_t( k.h ) );
   }
  };

 std::unordered_map< Key , NodeT * , Hash > f_map;           ///< the table

 };  // end( class( NodeTable ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_std_NodeTable */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File NodeTable.h --------------------------*/
/*--------------------------------------------------------------------------*/
