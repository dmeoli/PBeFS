/*--------------------------------------------------------------------------*/
/*------------------------ File AbstractionPDB.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the pattern databases of any abstraction whose states are
 * ranked, i.e., indexed by the integers \f$ [ 0 , N ) \f$: the concept
 * AbstractSpace, the construction of the table of the distances to the
 * abstract goals by a backward search on the workers of an executor, and
 * the class template AbstractionTable, which stores the distances in 4 or 8
 * bits each and saves them to a file.
 *
 * The construction is the algorithm of Dial (1969) on the table itself: the
 * entries of value d are final when the layer d is scanned, since every
 * edge costs at least 1, and each of them lowers its predecessors to
 * d + c. The table is cut in blocks of \f$ 2^{14} \f$ entries, each with the
 * least value not yet scanned, so that a layer skips the blocks that have
 * none of its entries, which are most of them in the first and last
 * layers:
 *
 * \code{.unparsed}
 *  1  procedure DISTANCES(S)
 *  2     h[r] ← ∞ for each r ∈ [0, N);  h[g] ← 0 for each goal g
 *  3     lo[B] ← least h of the block B;  next[B] ← ∞
 *  4     d ← 0
 *  5     while d < ∞ do
 *  6        for each block B with lo[B] ≤ d do          ▷ in parallel
 *  7           m ← ∞
 *  8           for each r ∈ B do
 *  9              if h[r] = d then
 * 10                 for each predecessor r' of r, by an edge of cost c, do
 * 11                    if d + c < h[r'] then          ▷ by compare-and-swap
 * 12                       h[r'] ← d + c
 *                          next[B(r')] ← min(next[B(r')], d + c)
 * 13              else if d < h[r] < ∞ then m ← min(m, h[r])
 * 14           lo[B] ← m
 * 15        lo[B] ← min(lo[B], next[B]);  next[B] ← ∞ for each block B
 * 16        d ← min_B lo[B]
 * \endcode
 *
 * The writes of a layer go to next (line 12) and not to lo, which the
 * workers read during the layer: a block scanned at line 14 forgets the
 * entries lowered after it read them, and line 15 brings them back. With
 * unit costs this is the breadth-first search by layers of the pattern
 * databases of the sliding-tile puzzles [see PatternDatabase.h]; with
 * integer costs it is a uniform-cost search. The entries that no goal
 * reaches keep the value 255, and an abstract state that a state of the
 * Domain maps to is never one of them.
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

#ifndef __PBeFS_AbstractionPDB
 #define __PBeFS_AbstractionPDB
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <filesystem>

#include <fstream>

#include <functional>

#include <iterator>

#include <memory>

#include <stdexcept>

#include <string>

#include <thread>

#include <vector>

#include "PBeFS/core/SearchUtils.h"

#include "ADS/succinct/PackedArray.h"

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

 /// a callable accepting a rank, used only to state the concepts
 struct rank_sink {
  void operator()( std::uint64_t ) const noexcept {}
  };

 /// a callable accepting a rank and the cost of an edge, used only to
 /// state the concepts
 struct rank_cost_sink {
  void operator()( std::uint64_t , unsigned ) const noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// an abstract space whose states are the ranks [ 0 , size() )
/** for_each_goal( f ) calls f( r ) for each abstract goal r, and
 * for_each_predecessor( r , f ) calls f( r' , c ) for each edge ( r' , r ) of
 * cost c, an integer from 1 to 254; the edges are those of the abstraction
 * of the moves of a Domain, so that the distance of the image of a state is
 * at most the cost of its path to a goal. Both are const and may be called
 * by many threads at once. */

//! [heuristics-sig-abstractspace]
template< class S >
concept AbstractSpace = requires( const S & s , std::uint64_t r ) {
 { s.size() } -> std::convertible_to< std::uint64_t >;
 s.for_each_goal( detail::rank_sink{} );
 s.for_each_predecessor( r , detail::rank_cost_sink{} );
 };
//! [heuristics-sig-abstractspace]

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the distances of the states of space to its goals, one byte each
/** Runs DISTANCES [see the file header] on \p threads workers of Executor
 * (all the cores if 0); an entry that no goal reaches is 255. Throws
 * std::logic_error if a distance would exceed 254. */

template< ADS::WorkerExecutor Executor = ADS::ThreadPool , AbstractSpace S >
[[nodiscard]] std::vector< std::uint8_t > abstract_distances(
 const S & space , unsigned threads = 0 )
{
 if( threads == 0 )
  threads = std::max( 1u , std::thread::hardware_concurrency() );
 constexpr std::uint8_t unseen = 0xFF;
 constexpr unsigned shift = 14;
 const std::uint64_t n = space.size();
 const std::uint64_t nb = ( n + ( std::uint64_t( 1 ) << shift ) - 1 ) >> shift;
 std::vector< std::uint8_t > h( n , unseen );
 std::vector< std::uint8_t > lo( nb , unseen );
 std::unique_ptr< std::atomic< std::uint8_t >[] > next(
  new std::atomic< std::uint8_t >[ nb ] );
 for( std::uint64_t b = 0 ; b < nb ; ++b )
  next[ b ].store( unseen , std::memory_order_relaxed );
 space.for_each_goal( [ & ]( std::uint64_t g ) {
  h[ g ] = 0;
  lo[ g >> shift ] = 0;
  } );

 for( unsigned d = 0 ; d < unseen ; ) {
  std::atomic< std::uint64_t > cursor{ 0 };
  std::atomic< bool > deep{ false };
  run_workers< Executor >( nullptr , threads , [ & ]( unsigned ) {
   for( ; ; ) {
    const std::uint64_t b = cursor.fetch_add( 1 , std::memory_order_relaxed );
    if( b >= nb )
     break;
    if( lo[ b ] > d )
     continue;
    unsigned m = unseen;
    const std::uint64_t e = std::min( n , ( b + 1 ) << shift );
    for( std::uint64_t r = b << shift ; r < e ; ++r ) {
     const unsigned v = std::atomic_ref< std::uint8_t >( h[ r ] ).load(
                         std::memory_order_relaxed );
     if( v == d )
      space.for_each_predecessor( r , [ & ]( std::uint64_t q , unsigned c ) {
       const unsigned nv = d + c;
       if( nv >= unseen ) {
        deep.store( true , std::memory_order_relaxed );
        return;
        }
       std::atomic_ref< std::uint8_t > x( h[ q ] );
       for( std::uint8_t o = x.load( std::memory_order_relaxed ) ; nv < o ; )
        if( x.compare_exchange_weak( o , std::uint8_t( nv ) ,
                                     std::memory_order_relaxed ) ) {
         std::atomic< std::uint8_t > & w = next[ q >> shift ];
         for( std::uint8_t t = w.load( std::memory_order_relaxed ) ;
              nv < t ; )
          if( w.compare_exchange_weak( t , std::uint8_t( nv ) ,
                                       std::memory_order_relaxed ) )
           break;
         break;
         }
       } );
     else
      if( ( v > d ) && ( v < m ) )
       m = v;
     }
    lo[ b ] = std::uint8_t( m );
    }
   } );
  if( deep.load() )
   throw( std::logic_error( "abstract_distances: a distance beyond 254" ) );
  unsigned least = unseen;
  for( std::uint64_t b = 0 ; b < nb ; ++b ) {
   lo[ b ] = std::min( lo[ b ] , next[ b ].load( std::memory_order_relaxed ) );
   next[ b ].store( unseen , std::memory_order_relaxed );
   least = std::min< unsigned >( least , lo[ b ] );
   }
  if( least <= d )
   throw( std::logic_error( "abstract_distances: an edge of cost 0" ) );
  d = least;
  }
 return( h );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS AbstractionTable -------------------------*/
/*--------------------------------------------------------------------------*/
/// the distances of an abstraction, B bits per entry (4 or 8)
/** The values above \f$ 2^B - 1 \f$ are stored as \f$ 2^B - 1 \f$, which
 * keeps the table admissible; with B = 4 the tables of Rubik's cube of
 * Korf (1997) take half a byte per entry, as in his paper. A table is
 * saved with a key that its owner computes from the abstraction (the
 * goal, the pattern, the ranking), and is loaded only with the same key. */

template< unsigned B = 8 >
class AbstractionTable {

 static_assert( ( B == 4 ) || ( B == 8 ) , "AbstractionTable: 4 or 8 bits" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 static constexpr unsigned bits = B;            ///< bits per entry

 /// the version of the format of the files
 static constexpr std::uint32_t file_version = 1;

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// the empty table
 AbstractionTable( void ) = default;

 /// the table of the distances h, as abstract_distances() gives them
 explicit AbstractionTable( const std::vector< std::uint8_t > & h )
  : f_entries( h.size() , 0 ) {
  for( std::size_t i = 0 ; i < h.size() ; ++i )
   f_entries.set( i , h[ i ] );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the rank r

 [[nodiscard]] unsigned operator[]( std::uint64_t r ) const noexcept {
  return( unsigned( f_entries.get( std::size_t( r ) ) ) );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the entries

 [[nodiscard]] std::uint64_t size( void ) const noexcept {
  return( f_entries.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the bytes of the entries

 [[nodiscard]] std::size_t bytes( void ) const noexcept {
  return( f_entries.bytes() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of the entries of each value, 0 to \f$ 2^B - 1 \f$

 [[nodiscard]] std::vector< std::uint64_t > histogram( void ) const {
  std::vector< std::uint64_t > c( std::size_t( 1 ) << B , 0 );
  for( std::size_t i = 0 ; i < f_entries.size() ; ++i )
   ++c[ std::size_t( f_entries.get( i ) ) ];
  return( c );
  }

/*------------------------------ PERSISTENCE -------------------------------*/
 /// writes the table to the file path with the given key
 /** The file holds, in little-endian order, the magic "PBeFSABS", the
  * version of the format, B, the key, the number of the entries, the
  * words of the entries and a checksum of all that precedes it; it is
  * written to a temporary file in the same directory, then renamed. */

 void save( const std::string & path , std::uint64_t key ) const {
  std::vector< std::uint8_t > b = header( key , f_entries.size() );
  for( std::uint64_t w : f_entries.words() )
   put( b , w , 8 );
  put( b , checksum( b.data() , b.size() ) , 8 );
  const std::string tmp = path + ".tmp" +
                          std::to_string( std::hash< std::thread::id >()(
                                           std::this_thread::get_id() ) );
  {
   std::ofstream out( tmp , std::ios::binary | std::ios::trunc );
   out.write( reinterpret_cast< const char * >( b.data() ) ,
              std::streamsize( b.size() ) );
   if( ! out )
    throw( std::runtime_error( "AbstractionTable::save: cannot write " +
                               tmp ) );
   }
  std::filesystem::rename( tmp , path );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// reads the table of n entries with the key from the file path
 /** Throws std::runtime_error if the file cannot be read, if its magic,
  * version, checksum or length are wrong, or if it holds another table. */

 [[nodiscard]] static AbstractionTable load( const std::string & path ,
                                             std::uint64_t key ,
                                             std::uint64_t n ) {
  std::ifstream in( path , std::ios::binary );
  if( ! in )
   throw( std::runtime_error( "AbstractionTable::load: cannot read " +
                              path ) );
  std::vector< std::uint8_t > b( ( std::istreambuf_iterator< char >( in ) ) ,
                                 std::istreambuf_iterator< char >() );
  const std::vector< std::uint8_t > head = header( key , n );
  const std::uint64_t words = ( n * B + 63 ) / 64;
  if( b.size() != head.size() + 8 * words + 8 )
   throw( std::runtime_error( "AbstractionTable::load: " + path +
                              " has a wrong length" ) );
  if( get( b.data() + b.size() - 8 , 8 ) !=
      checksum( b.data() , b.size() - 8 ) )
   throw( std::runtime_error( "AbstractionTable::load: " + path +
                              " is corrupted (checksum)" ) );
  if( ! std::equal( head.begin() , head.end() , b.begin() ) )
   throw( std::runtime_error( "AbstractionTable::load: " + path +
                              " holds another table" ) );
  AbstractionTable t;
  t.f_entries = ADS::PackedArray< B , true >( std::size_t( n ) , 0 );
  const std::uint8_t * p = b.data() + head.size();
  constexpr std::uint64_t per = 64 / B;
  for( std::uint64_t i = 0 ; i < n ; ++i ) {
   const std::uint64_t w = get( p + 8 * ( i / per ) , 8 );
   t.f_entries.set( std::size_t( i ) ,
                    ( w >> ( ( i % per ) * B ) ) & ( ( 1u << B ) - 1 ) );
   }
  return( t );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the table from the file dir / name with the key, or built by make()
 /// and saved there (dir created if needed); "" as dir saves nothing

 template< class Make >
 [[nodiscard]] static AbstractionTable cached( const std::string & dir ,
                                               const std::string & name ,
                                               std::uint64_t key ,
                                               std::uint64_t n ,
                                               Make && make ) {
  namespace fs = std::filesystem;
  if( dir.empty() )
   return( AbstractionTable( make() ) );
  const fs::path f = fs::path( dir ) / name;
  if( fs::exists( f ) )
   try {
    return( load( f.string() , key , n ) );
    }
   catch( std::runtime_error & ) {}            // built again below
  fs::create_directories( dir );
  AbstractionTable t( make() );
  t.save( f.string() , key );
  return( t );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the checksum of the n bytes at b, by words of 8 bytes

 [[nodiscard]] static std::uint64_t checksum( const std::uint8_t * b ,
                                              std::size_t n ) noexcept {
  std::uint64_t h = ADS::splitmix64( 0x5042654653414253ULL ^ n );
  std::size_t i = 0;
  for( ; i + 8 <= n ; i += 8 )
   h = ADS::splitmix64( h ^ get( b + i , 8 ) ) + i;
  if( i < n )
   h = ADS::splitmix64( h ^ get( b + i , int( n - i ) ) );
  return( h );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the first bytes of a file: magic, version, bits, key and entries
 [[nodiscard]] static std::vector< std::uint8_t > header( std::uint64_t key ,
                                                          std::uint64_t n ) {
  std::vector< std::uint8_t > b = { 'P' , 'B' , 'e' , 'F' , 'S' , 'A' , 'B' ,
                                    'S' };
  put( b , file_version , 4 );
  put( b , B , 4 );
  put( b , key , 8 );
  put( b , n , 8 );
  return( b );
  }

 // appends the k low bytes of x, the lowest first
 static void put( std::vector< std::uint8_t > & b , std::uint64_t x , int k ) {
  for( int i = 0 ; i < k ; ++i )
   b.push_back( std::uint8_t( x >> ( 8 * i ) ) );
  }

 // the integer of the k bytes at b, the lowest first
 [[nodiscard]] static std::uint64_t get( const std::uint8_t * b , int k )
  noexcept {
  std::uint64_t x = 0;
  for( int i = 0 ; i < k ; ++i )
   x |= std::uint64_t( b[ i ] ) << ( 8 * i );
  return( x );
  }

 ADS::PackedArray< B , true > f_entries;        ///< the entries

 };  // end( class( AbstractionTable ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AbstractionPDB */

/*--------------------------------------------------------------------------*/
/*----------------------- End File AbstractionPDB.h ------------------------*/
/*--------------------------------------------------------------------------*/
