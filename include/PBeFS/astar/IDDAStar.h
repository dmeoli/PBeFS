/*--------------------------------------------------------------------------*/
/*----------------------------- File IDDAStar.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template IDDAStar, the A* search with
 * immediate duplicate detection in external memory, A*-IDD, of Lin and
 * Fukunaga (AAAI 2018): the closed list is a hash table on the disk, read
 * at random, as fits a solid state drive, and the memory holds only a
 * compact table of the positions of its records.
 *
 * The OPEN list is a set of buckets \f$ ( f , g ) \f$, each a stack of
 * records on the disk with its top in memory, popped by least f and then
 * largest g. The CLOSED list has three parts: the records, in a file
 * written only at its end; a table in memory of \f$ S \f$ places of
 * \f$ \lceil \log_2 ( S + 1 ) \rceil \f$ bits each, holding the position
 * plus one of a record in the file (0 for an empty place), probed linearly
 * from the hash of the state; and \f$ P \f$ buffers in memory, one per
 * partition of the states by a second hash, which collect the states
 * closed until a buffer has \f$ B \f$ of them and is written to the file as
 * one block. Since each block holds the states of one partition, a table of
 * one partition number per block tells whether the record that a place
 * points to may be the state looked for, and the file is read only when it
 * may [Lin and Fukunaga, section 4]:
 *
 * \code{.unparsed}
 *  1  procedure A*-IDD(s₀)
 *  2     push (s₀, g = 0) in OPEN
 *  3     while OPEN ≠ ∅ do
 *  4        n ← pop the top of the bucket of least f and largest g
 *  5        if n is in the buffer of its partition or, probing the table
 *              and reading only the records of blocks of its partition,
 *              in the file, with a g not larger, then continue
 *  6        put n in the buffer of its partition (or update its record)
 *  7        if the buffer has B states, append them to the file as a
 *              block, and put their positions in the table
 *  8        if n is a goal then return the path, traced by the parents
 *  9        push each successor of n, with its g and parent, in OPEN
 * 10     return unsolvable
 * \endcode
 *
 * The duplicates are detected when a node is popped, before its
 * expansion, so that a state is expanded once with a consistent heuristic
 * and reopened, its record rewritten in place, when it is popped with a
 * smaller g otherwise. The memory is that of the table, about \f$ S \log_2
 * S / 8 \f$ bytes, of the buffers, \f$ P B \f$ records, and of the tops of
 * the buckets; the disk holds the records of the closed and of the open
 * states. A record holds the state of the parent, by which the path is
 * traced back from the goal. The states must be trivially copyable and the
 * costs integers (the concept DiskDomain); the search stops with
 * Status::node_limit when the table is 7/8 full or the node_limit of the
 * SearchOptions is reached, and with Status::disk_limit when the files
 * exceed IDDOptions::disk_limit.
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

#ifndef __PBeFS_IDDAStar
 #define __PBeFS_IDDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <filesystem>

#include <fstream>

#include <map>

#include <memory>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "ADS/dictionary/FlatHashTable.h"

#include "ADS/external/RecordFile.h"

#include "PBeFS/astar/DDDAStar.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS IDDOptions -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of IDDAStar, besides those of all the searches

//! [bidirectional-sig-iddoptions]
struct IDDOptions {

 /// where the files go, in a directory of their own removed at the end;
 /// the temporary directory of the system if empty
 std::filesystem::path directory;

 /// the places S of the table of the positions in memory, which bound
 /// the closed states to 7 S / 8; at least 8 blocks
 std::size_t slots = std::size_t( 1 ) << 20;

 /// the partitions P of the states, one buffer each
 std::size_t partitions = 100;

 /// the states B of a full buffer, i.e., of a block of the file
 std::size_t block = std::size_t( 1 ) << 10;

 /// the records of the top of an open bucket kept in memory
 std::size_t open_block = std::size_t( 1 ) << 10;

 /// the most bytes the files may hold on the disk (0 = no cap)
 std::uint64_t disk_limit = 0;

 };  // end( struct( IDDOptions ) )
//! [bidirectional-sig-iddoptions]

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS IDDAStar -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the A* search of a Domain with immediate duplicate detection on the disk
/** The counters of the Solution: expanded, reexpanded, generated,
 * duplicates (the states popped already closed with a g not larger),
 * stored (the closed states), sent (the bytes written to the disk) and
 * messages (the records of the closed file read at random, i.e., the
 * probes that reached the disk). */

template< DiskDomain D >
class IDDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the IDDAStar)

 explicit IDDAStar( const D & domain , const SearchOptions & options = {} ,
                    IDDOptions idd = {} )
  : f_domain( domain ) , f_options( options ) , f_idd( std::move( idd ) ) {
  if( ( f_idd.partitions == 0 ) || ( f_idd.block == 0 ) ||
      ( f_idd.open_block == 0 ) || ( f_idd.slots < 8 * f_idd.block ) )
   throw( std::invalid_argument( "IDDAStar: an option is 0, or the places "
                                 "are fewer than 8 blocks" ) );
  if( options.weight != 1 )
   throw( std::invalid_argument( "IDDAStar: no weight" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Stopwatch clock;
  Run r( f_domain , f_options , f_idd );
  SolutionT sol = r.search();
  sol.stats = r.st;
  sol.stats.expanded_by = { sol.stats.expanded };
  sol.stats.stored_by = { sol.stats.stored };
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a record: a state, its g, the state of its parent, its operator
 struct Rec {
  State s;
  State parent;
  Cost g;
  [[no_unique_address]] OperatorOf< D > o;
  };

 /// what a buffer keeps of a state
 struct Val {
  State parent;
  Cost g;
  [[no_unique_address]] OperatorOf< D > o;
  };

 struct Hash {
  const D * d;
  std::uint64_t operator()( const State & s ) const {
   return( d->hash( s ) );
   }
  };

 /// a stack of records on the disk, its top in memory
 struct Stack {
  std::filesystem::path path;           ///< the file of the bottom
  std::size_t on_disk = 0;              ///< the records in the file
  std::vector< Rec > top;               ///< the records in memory
  };

 /// a bucket, ordered by least f and then largest g
 struct Key {
  Cost f , g;
  bool operator<( const Key & o ) const {
   return( ( f < o.f ) || ( ( f == o.f ) && ( g > o.g ) ) );
   }
  };

 /// the state of a search
 struct Run {

  const D & d;
  const SearchOptions & opt;
  const IDDOptions & idd;
  ADS::TempDirectory dir;
  std::map< Key , Stack > open;           ///< the buckets
  unsigned bits;                          ///< of a place of the table
  std::vector< std::uint64_t > table;     ///< the places, packed
  std::size_t used = 0;                   ///< the places taken
  std::vector< std::uint32_t > blocks;    ///< the partition of each block
  std::vector< ADS::FlatHashTable< State , Val , Hash > > buf;  ///< buffers
  std::fstream closed;                    ///< the records of the closed
  std::uint64_t n_closed = 0;             ///< the records in that file
  std::uint64_t disk = 0;                 ///< the bytes on the disk
  SearchStats st;
  std::size_t counter = 0;                ///< files made

  Run( const D & domain , const SearchOptions & o , const IDDOptions & x )
   : d( domain ) , opt( o ) , idd( x ) , dir( x.directory , "pbefs-idd-" ) ,
     bits( unsigned( std::bit_width( std::uint64_t( x.slots ) ) ) ) ,
     table( ( x.slots * bits + 63 ) / 64 + 1 , 0 ) {
   buf.reserve( x.partitions );
   for( std::size_t p = 0 ; p < x.partitions ; ++p )
    buf.emplace_back( 0 , Hash{ &domain } );
   const auto path = dir / "closed";
   { std::ofstream make( path , std::ios::binary ); }
   closed.open( path , std::ios::in | std::ios::out | std::ios::binary );
   if( ! closed )
    throw( std::runtime_error( "IDDAStar: cannot open the closed file" ) );
   }

  // the place i of the table
  std::uint64_t place( std::size_t i ) const {
   const std::size_t b = i * bits , w = b / 64 , k = b % 64;
   std::uint64_t v = table[ w ] >> k;
   if( k + bits > 64 )
    v |= table[ w + 1 ] << ( 64 - k );
   return( v & ( ( std::uint64_t( 1 ) << bits ) - 1 ) );
   }

  void set_place( std::size_t i , std::uint64_t v ) {
   const std::size_t b = i * bits , w = b / 64 , k = b % 64;
   const std::uint64_t m = ( std::uint64_t( 1 ) << bits ) - 1;
   table[ w ] = ( table[ w ] & ~( m << k ) ) | ( v << k );
   if( k + bits > 64 ) {
    const unsigned r = unsigned( 64 - k );
    table[ w + 1 ] = ( table[ w + 1 ] & ~( m >> r ) ) | ( v >> r );
    }
   }

  std::size_t partition( const State & s ) const {
   return( std::size_t( ADS::splitmix64( d.hash( s ) ^ 0x9E3779B97F4A7C15ULL )
                        % idd.partitions ) );
   }

  Rec read( std::uint64_t i ) {
   Rec r;
   closed.seekg( std::streamoff( i * sizeof( Rec ) ) );
   closed.read( reinterpret_cast< char * >( &r ) , sizeof( Rec ) );
   ++st.messages;
   return( r );
   }

  void write( std::uint64_t i , const Rec & r ) {
   closed.seekp( std::streamoff( i * sizeof( Rec ) ) );
   closed.write( reinterpret_cast< const char * >( &r ) , sizeof( Rec ) );
   st.sent += sizeof( Rec );
   }

  // the position of s in the closed file, or ~0
  std::uint64_t find_closed( const State & s , std::size_t p ) {
   const std::size_t S = idd.slots;
   for( std::size_t i = std::size_t( d.hash( s ) % S ) ; ;
        i = ( i + 1 ) % S ) {
    const std::uint64_t v = place( i );
    if( v == 0 )
     return( ~std::uint64_t( 0 ) );
    const std::uint64_t at = v - 1;
    if( blocks[ at / idd.block ] != p )      // another partition: skipped
     continue;
    if( read( at ).s == s )
     return( at );
    }
   }

  // writes the full buffer p to the closed file as a block
  void flush( std::size_t p ) {
   std::vector< Rec > v;
   v.reserve( buf[ p ].size() );
   buf[ p ].forEach( [ & ]( const State & s , const Val & x ) {
    v.push_back( Rec{ s , x.parent , x.g , x.o } );
    } );
   closed.seekp( std::streamoff( n_closed * sizeof( Rec ) ) );
   closed.write( reinterpret_cast< const char * >( v.data() ) ,
                 std::streamsize( v.size() * sizeof( Rec ) ) );
   st.sent += v.size() * sizeof( Rec );
   disk += v.size() * sizeof( Rec );
   for( const Rec & r : v ) {
    const std::size_t S = idd.slots;
    std::size_t i = std::size_t( d.hash( r.s ) % S );
    while( place( i ) != 0 )
     i = ( i + 1 ) % S;
    set_place( i , ++n_closed );
    ++used;
    }
   blocks.push_back( std::uint32_t( p ) );
   buf[ p ].clear();
   }

  void push( const Rec & r , Cost h ) {
   auto it = open.find( Key{ Cost( r.g + h ) , r.g } );
   if( it == open.end() ) {
    Stack x;
    x.path = dir / ( "o" + std::to_string( counter++ ) );
    it = open.emplace( Key{ Cost( r.g + h ) , r.g } , std::move( x ) ).first;
    }
   Stack & k = it->second;
   k.top.push_back( r );
   if( k.top.size() >= 2 * idd.open_block ) {    // the bottom half out
    ADS::RecordFile< Rec > f( k.path );
    f.append( k.top.data() , idd.open_block );
    k.top.erase( k.top.begin() ,
                 k.top.begin() + std::ptrdiff_t( idd.open_block ) );
    k.on_disk += idd.open_block;
    st.sent += idd.open_block * sizeof( Rec );
    disk += idd.open_block * sizeof( Rec );
    }
   }

  bool pop( Rec & r ) {
   while( ! open.empty() ) {
    auto it = open.begin();
    Stack & k = it->second;
    if( k.top.empty() && ( k.on_disk > 0 ) ) {  // the last block back in
     const std::size_t n = std::min( k.on_disk , idd.open_block );
     k.top.resize( n );
     {
      std::ifstream in( k.path , std::ios::binary );
      in.seekg( std::streamoff( ( k.on_disk - n ) * sizeof( Rec ) ) );
      in.read( reinterpret_cast< char * >( k.top.data() ) ,
               std::streamsize( n * sizeof( Rec ) ) );
      }
     k.on_disk -= n;
     std::filesystem::resize_file( k.path , k.on_disk * sizeof( Rec ) );
     disk -= n * sizeof( Rec );
     }
    if( k.top.empty() ) {
     std::error_code e;
     std::filesystem::remove( k.path , e );
     open.erase( it );
     continue;
     }
    r = k.top.back();
    k.top.pop_back();
    return( true );
    }
   return( false );
   }

  // the record of the closed state s, from its buffer or the file
  bool lookup( const State & s , Rec & r ) {
   const std::size_t p = partition( s );
   if( const Val * v = buf[ p ].find( s ) ) {
    r = Rec{ s , v->parent , v->g , v->o };
    return( true );
    }
   const std::uint64_t at = find_closed( s , p );
   if( at == ~std::uint64_t( 0 ) )
    return( false );
   r = read( at );
   return( true );
   }

  SolutionT search( void ) {
   SolutionT sol;
   const State s0 = d.initial_state();
   push( Rec{ s0 , s0 , Cost( 0 ) , {} } , Cost( d.heuristic( s0 ) ) );
   Rec n;
   while( pop( n ) ) {
    const std::size_t p = partition( n.s );
    bool reopen = false;
    if( Val * v = buf[ p ].find( n.s ) ) {      // closed, in the buffer
     if( v->g <= n.g ) {
      ++st.duplicates;
      continue;
      }
     *v = Val{ n.parent , n.g , n.o };
     reopen = true;
     }
    else {
     const std::uint64_t at = find_closed( n.s , p );
     if( at != ~std::uint64_t( 0 ) ) {          // closed, in the file
      if( read( at ).g <= n.g ) {
       ++st.duplicates;
       continue;
       }
      write( at , n );
      reopen = true;
      }
     else {
      buf[ p ].insert( n.s , Val{ n.parent , n.g , n.o } );
      ++st.stored;
      if( buf[ p ].size() >= idd.block )
       flush( p );
      }
     }
    ++st.expanded;
    if( reopen )
     ++st.reexpanded;
    if( d.is_goal( n.s ) ) {
     sol.status = Status::solved;
     sol.cost = n.g;
     sol.bound = sol.cost;
     sol.path = trace( n );
     return( sol );
     }
    const Cost h = Cost( d.heuristic( n.s ) );
    for_each_child( d , n.s , h , n.o.op , [ & ]( const State & t , Cost c ,
                                                  Cost ht ,
                                                  std::uint8_t ot ) {
     ++st.generated;
     push( Rec{ t , n.s , Cost( n.g + c ) , ot } , ht );
     } );
    if( ( used * 8 >= idd.slots * 7 ) ||
        ( opt.node_limit && ( st.stored >= opt.node_limit ) ) ) {
     sol.status = Status::node_limit;
     return( sol );
     }
    if( idd.disk_limit && ( disk > idd.disk_limit ) ) {
     sol.status = Status::disk_limit;
     return( sol );
     }
    }
   sol.status = Status::unsolvable;
   return( sol );
   }

  // the path from the start to the goal record x
  std::vector< State > trace( Rec x ) {
   std::vector< State > path{ x.s };
   while( ! ( ( x.g == Cost( 0 ) ) && ( x.parent == x.s ) ) ) {
    Rec y;
    if( ! lookup( x.parent , y ) )
     throw( std::logic_error( "IDDAStar: the parent of a state is not "
                              "closed" ) );
    x = y;
    path.push_back( x.s );
    }
   std::reverse( path.begin() , path.end() );
   return( path );
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 IDDOptions f_idd;                                   ///< the options of IDD

 };  // end( class( IDDAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with A*-IDD

template< DiskDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > idd_astar(
 const D & domain , const SearchOptions & options = {} ,
 const IDDOptions & idd = {} )
{
 return( IDDAStar< D >( domain , options , idd ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_IDDAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File IDDAStar.h --------------------------*/
/*--------------------------------------------------------------------------*/
