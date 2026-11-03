/*--------------------------------------------------------------------------*/
/*---------------------------- File DDDAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template DDDAStar, the A* search with delayed
 * duplicate detection, whose OPEN and CLOSED lists are files on the disk:
 * the External A* of Edelkamp, Jabbar and Schrödl (2004), with the
 * hash-based delayed duplicate detection of Korf (2008) in place of the
 * external sorts.
 *
 * The nodes are kept in buckets \f$ ( g , h ) \f$, one per pair of values,
 * each bucket split into \f$ P \f$ files by the hash of the state, so that
 * two copies of a state are always in the same file of their buckets. With
 * a consistent \f$ h \f$ and positive integer costs, the buckets are
 * expanded by increasing \f$ f = g + h \f$ and, for equal \f$ f \f$, by
 * increasing \f$ g \f$ (the successors of a bucket of the same \f$ f \f$
 * have a larger \f$ g \f$), and a state reached again has the same
 * \f$ h \f$ and a \f$ g \f$ not smaller than that of its first copy; hence
 * its duplicates are within its bucket or in the closed buckets of the same
 * \f$ h \f$, and nowhere else:
 *
 * \code{.unparsed}
 *  1  procedure DDD-A*(s₀)
 *  2     OPEN(0, h(s₀)) ← {s₀}
 *  3     while some OPEN(g, h) is not empty do
 *  4        (g, h) ← the bucket of least g + h, and then of least g
 *  5        for p ← 0 to P − 1 do
 *  6           T ← the states of the file p of OPEN(g, h), without copies
 *  7           for each closed bucket (g', h) with g' < g do
 *  8              T ← T ∖ the file p of CLOSED(g', h)
 *  9           the file p of CLOSED(g, h) ← T
 * 10           if T has a goal then return the path to it, of cost g
 * 11           for each u in T and each successor v of u do
 * 12              append (v, hash(u)) to the file hash(v) mod P of
 *                 OPEN(g + c(u, v), h(v))
 * 13     return unsolvable
 * \endcode
 *
 * Only one file of one bucket is in memory at a time (line 6), in a hash
 * table, while the closed files are scanned (line 8) and the successors are
 * appended through a buffer per file of the buckets they go to (line 12),
 * so that the memory of the search is that of the largest file, about
 * \f$ 1 / P \f$ of the largest bucket, and not that of the states
 * generated. The files of the closed buckets are kept to the end, since a
 * duplicate may be in any of them with the same \f$ h \f$; a domain whose
 * duplicates are only a few layers back (the undirected graphs with unit
 * costs, as in the frontier search of Korf et al.) could drop the older
 * ones, which is not done here. Each record carries the hash of the parent
 * of its state, by which the path is traced back from the goal at the end:
 * the parent is looked for in the closed files of its partition, with a
 * smaller \f$ g \f$, and checked to have the state as a successor at the
 * cost of the difference of the two \f$ g \f$.
 *
 * The states must be trivially copyable, since they are written to the
 * files as their bytes, and the costs integers, since a bucket is a pair of
 * values (the concept DiskDomain); a successor of a smaller \f$ f \f$ than its parent (a heuristic
 * that is not consistent) or an edge of a cost that is not positive throws
 * std::logic_error, since the order of the buckets would no longer make the
 * first goal optimal.
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

#ifndef __PBeFS_DDDAStar
 #define __PBeFS_DDDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <filesystem>

#include <map>

#include <memory>

#include <set>

#include <stdexcept>

#include <string>

#include <tuple>

#include <type_traits>

#include <utility>

#include <vector>

#include "ADS/dictionary/FlatHashTable.h"

#include "ADS/external/RecordFile.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Domain that DDDAStar can search: integer costs, since a bucket is a
/// pair of values, and states that are trivially copyable, since they are
/// written to the disk as their bytes

template< class D >
concept DiskDomain = Domain< D > && std::is_integral_v< cost_t< D > > &&
 std::is_trivially_copyable_v< state_t< D > >;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of DDDAStar, besides those of all the searches

struct DDDOptions {
 /// the files of a bucket, i.e., the parts by the hash of the state: the
 /// largest bucket over this number is what the search holds in memory
 std::size_t partitions = 64;
 /// where the files go, in a directory of their own removed at the end;
 /// the temporary directory of the system if empty
 std::filesystem::path directory;
 /// the bytes of the buffer of each file being written, i.e., of each
 /// file of the buckets that the successors of a bucket go to
 std::size_t buffer_bytes = std::size_t( 1 ) << 16;
 };  // end( struct( DDDOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS DDDAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the A* search of a Domain with delayed duplicate detection on the disk
/** The counters of the Solution: expanded counts the states expanded, each
 * once, generated the successors, duplicates the copies dropped (within a
 * file and against the closed files), and stored the states written to the
 * closed files, i.e., the distinct states expanded. The node_limit of the
 * SearchOptions bounds stored, and the search stops with
 * Status::node_limit when it reaches it; expected_nodes is not used. */

template< DiskDomain D >
class DDDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the DDDAStar)

 explicit DDDAStar( const D & domain , const SearchOptions & options = {} ,
                    DDDOptions ddd = {} )
  : f_domain( domain ) , f_options( options ) , f_ddd( std::move( ddd ) )
 {
  if( f_ddd.partitions == 0 )
   throw( std::invalid_argument( "DDDAStar: no partition" ) );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  Run r( f_domain , f_options , f_ddd );
  SolutionT sol = r.search();
  sol.stats = std::move( r.st );
  sol.stats.expanded_by = { sol.stats.expanded };
  sol.stats.stored_by = { sol.stats.stored };
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a record of a file: a state and the hash of its parent
 struct Rec {
  State s;
  std::uint64_t parent;
  };

 /// the hash of the Domain, for the table of a file
 struct Hash {
  const D * d;
  std::uint64_t operator()( const State & s ) const {
   return( d->hash( s ) );
   }
  };

 using Table = ADS::FlatHashTable< State , std::uint64_t , Hash >;

 /// a bucket, ordered by f and then by g
 struct Key {
  Cost f;
  Cost g;
  bool operator<( const Key & o ) const {
   return( ( f < o.f ) || ( ( f == o.f ) && ( g < o.g ) ) );
   }
  };

 /// the state of a search: the directory, the buckets and the counters
 struct Run {

  const D & d;                                    ///< the Domain
  const SearchOptions & opt;                      ///< the options
  const DDDOptions & ddd;                         ///< the options of DDD
  ADS::TempDirectory dir;                         ///< where the files are
  std::size_t buffer;                             ///< records per buffer
  std::set< Key > open;                           ///< the open buckets
  std::map< Cost , std::vector< Cost > > closed;  ///< the g of the closed
                                                  ///< buckets, by h
  SearchStats st;                                 ///< the counters

  Run( const D & domain , const SearchOptions & o , const DDDOptions & x )
   : d( domain ) , opt( o ) , ddd( x ) , dir( x.directory , "pbefs-ddd-" ) ,
     buffer( std::max< std::size_t >( 1 , x.buffer_bytes / sizeof( Rec ) ) )
   {}

  // the file p of the open (o) or closed (c) bucket ( g , h )
  [[nodiscard]] ADS::RecordFile< Rec > file( char kind , Cost g , Cost h ,
                                             std::size_t p ) const {
   return( ADS::RecordFile< Rec >( dir / ( std::string( 1 , kind ) + "_" +
                                           std::to_string( g ) + "_" +
                                           std::to_string( h ) + "_" +
                                           std::to_string( p ) ) ) );
   }

  [[nodiscard]] std::size_t part( const State & s ) const {
   return( std::size_t( d.hash( s ) % ddd.partitions ) );
   }

  SolutionT search( void ) {
   SolutionT sol;
   const State s0 = d.initial_state();
   const Cost h0 = d.heuristic( s0 );
   {
    const Rec root{ s0 , d.hash( s0 ) };
    file( 'o' , Cost( 0 ) , h0 , part( s0 ) ).append( &root , 1 );
    open.insert( Key{ h0 , Cost( 0 ) } );
    }

   while( ! open.empty() ) {
    const Key k = *open.begin();
    open.erase( open.begin() );
    const Cost g = k.g , h = k.f - k.g;
    // the files of the buckets ( g' , h' ) that the successors go to, by
    // ( g' , h' , p ), with their writers, kept for this bucket only
    std::map< std::tuple< Cost , Cost , std::size_t > ,
              std::pair< std::unique_ptr< ADS::RecordFile< Rec > > ,
                         std::unique_ptr< ADS::RecordWriter< Rec > > > > out;
    const auto writer = [ & ]( Cost gv , Cost hv , std::size_t p )
     -> ADS::RecordWriter< Rec > & {
     auto & slot = out[ { gv , hv , p } ];
     if( ! slot.second ) {
      slot.first = std::make_unique< ADS::RecordFile< Rec > >(
                                               file( 'o' , gv , hv , p ) );
      slot.second = std::make_unique< ADS::RecordWriter< Rec > >(
                                               *slot.first , buffer );
      }
     return( *slot.second );
     };
    for( std::size_t p = 0 ; p < ddd.partitions ; ++p ) {
     auto in = file( 'o' , g , h , p );
     if( in.empty() ) {
      in.remove();
      continue;
      }
     Table t( in.size() , Hash{ &d } );
     in.forEach( [ & ]( const Rec & x ) {
      if( ! t.insert( x.s , x.parent ) )
       ++st.duplicates;                           // a copy in the file
      } );
     in.remove();
     const auto cl = closed.find( h );
     if( cl != closed.end() )
      for( Cost gc : cl->second )
       file( 'c' , gc , h , p ).forEach( [ & ]( const Rec & x ) {
        if( t.remove( x.s ) )
         ++st.duplicates;                         // closed already
        } );
     if( t.empty() )
      continue;

     std::vector< Rec > fresh;                    // the closed file p
     fresh.reserve( t.size() );
     t.forEach( [ & ]( const State & s , std::uint64_t par ) {
      fresh.push_back( Rec{ s , par } );
      } );
     file( 'c' , g , h , p ).append( fresh.data() , fresh.size() );
     st.stored += fresh.size();
     closed[ h ];                                 // the list of this h

     for( const Rec & x : fresh )
      if( d.is_goal( x.s ) ) {
       mark_closed( g , h );
       sol.status = Status::solved;
       sol.cost = g;
       sol.path = trace( x , g );
       return( sol );
       }

     for( const Rec & x : fresh ) {
      ++st.expanded;
      const std::uint64_t hu = d.hash( x.s );
      d.for_each_successor( x.s , h , [ & ]( const State & v , Cost c ,
                                             Cost hv ) {
       ++st.generated;
       if( c <= Cost( 0 ) )
        throw( std::logic_error( "DDDAStar: an edge of a cost that is not"
                                 " positive" ) );
       if( g + c + hv < k.f )
        throw( std::logic_error( "DDDAStar: the heuristic is not"
                                 " consistent" ) );
       writer( g + c , hv , part( v ) ).push( Rec{ v , hu } );
       } );
      }

     if( opt.node_limit && ( st.stored >= opt.node_limit ) ) {
      mark_closed( g , h );
      sol.status = Status::node_limit;
      return( sol );
      }
     }

    for( auto & [ key , fw ] : out ) {
     fw.second->flush();
     const auto [ gv , hv , p ] = key;
     open.insert( Key{ gv + hv , gv } );
     }
    mark_closed( g , h );
    }
   sol.status = Status::unsolvable;
   return( sol );
   }

  // the bucket ( g , h ) among the closed ones of h, if it has a file
  void mark_closed( Cost g , Cost h ) {
   auto it = closed.find( h );
   if( ( it != closed.end() ) &&
       ( std::find( it->second.begin() , it->second.end() , g ) ==
         it->second.end() ) )
    it->second.push_back( g );
   }

  // the path from s0 to the goal x, of cost g: each parent
  // is looked for in the closed files of its partition, among the states
  // of smaller g whose hash is the one recorded and that have the child as
  // a successor at the cost of the difference
  std::vector< State > trace( Rec x , Cost g ) {
   std::vector< State > path{ x.s };
   const State s0 = d.initial_state();
   while( ! ( ( g == Cost( 0 ) ) && ( x.s == s0 ) ) ) {
    const std::size_t p = std::size_t( x.parent % ddd.partitions );
    bool found = false;
    for( const auto & [ hc , gs ] : closed ) {
     for( Cost gc : gs ) {
      if( gc >= g )
       continue;
      file( 'c' , gc , hc , p ).forEach( [ & ]( const Rec & y ) {
       if( found || ( d.hash( y.s ) != x.parent ) )
        return;
       d.for_each_successor( y.s , hc , [ & ]( const State & v , Cost c ,
                                               Cost ) {
        if( ! found && ( v == x.s ) && ( gc + c == g ) ) {
         found = true;
         x = y;
         g = gc;
         }
        } );
       } );
      if( found )
       break;
      }
     if( found )
      break;
     }
    if( ! found )
     throw( std::logic_error( "DDDAStar: the parent of a state is not in"
                              " the closed files" ) );
    path.push_back( x.s );
    }
   std::reverse( path.begin() , path.end() );
   return( path );
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 DDDOptions f_ddd;                                   ///< the options of DDD

 };  // end( class( DDDAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with A* and delayed duplicate detection on the disk

template< DiskDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > ddd_astar(
 const D & domain , const SearchOptions & options = {} ,
 const DDDOptions & ddd = {} )
{
 return( DDDAStar< D >( domain , options , ddd ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_DDDAStar */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File DDDAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
