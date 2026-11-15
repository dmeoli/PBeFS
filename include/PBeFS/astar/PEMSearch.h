/*--------------------------------------------------------------------------*/
/*---------------------------- File PEMSearch.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template PEMSearch, the parallel
 * external-memory searches of the framework PEM-BiHS of Siag, Shperberg,
 * Felner and Sturtevant (ECAI 2024): PEM-A*, the reverse PEM-rA* (PEM-A*
 * on a ReverseDomain), PEMM of Sturtevant and Chen (IJCAI 2016) and
 * PEM-BAE*, whose open and closed lists are files on the disk, the nodes of
 * one bucket being loaded, checked for duplicates and expanded by the
 * workers of an Executor.
 *
 * A node of the direction D is kept in the bucket \f$ ( D , g_D , h_F ,
 * h_B ) \f$ (the heuristic of the other direction being 0 for PEM-A*),
 * which carries all that the priorities and the bounds need, and the
 * bucket is split into \f$ P \f$ partitions by the hash of the state, each
 * partition being a file per worker that writes to it. Two copies of a
 * state have the same \f$ h_F \f$ and \f$ h_B \f$ and the same partition,
 * so that the duplicates of a partition are in that partition of the
 * buckets with the same heuristics, and so are the copies of the state in
 * the other direction. The search is algorithm 1 of Siag et al.:
 *
 * \code{.unparsed}
 *  1  procedure PEM-BiHS(s, t)
 *  2     U ← ∞; OPEN_F ← {s}; OPEN_B ← {t}
 *  3     while OPEN_F ≠ ∅ ∧ OPEN_B ≠ ∅ ∧ U > LB do
 *  4        D ← the direction chosen; b ← the best bucket of OPEN_D
 *  5        for each partition p of b, in parallel, do
 *  6           T ← the states of the files of p, without copies
 *  7           T ← T ∖ the partition p of the closed buckets of D with the
 *                 same h_F and h_B and a smaller g
 *  8           write T to the partition p of CLOSED_D(b)
 *  9           U ← min(U, g + g') for each state of T found, with g', in
 *                 the partition p of the buckets of the other direction
 *                 with the same h_F and h_B
 * 10           for each u in T and each child v of u in D do
 *                 append v to the partition of v of its bucket of OPEN_D
 * 11     return U
 * \endcode
 *
 * Line 9 is the delayed solution detection of Sturtevant and Chen: a
 * meeting is found when one of the two copies is expanded, against the
 * open and the closed buckets of the other direction; PEM-A* instead
 * tests the goal when a node is generated, as Siag et al. do. Since the
 * meeting is found late, the bound of MM cannot use \f$ \epsilon \f$
 * (Sturtevant and Chen, section 3.3), and LB is \f$ \min f \f$ of the open
 * buckets for PEM-A*, \f$ \max\{ C , f_{\min F} , f_{\min B} , g_{\min F} +
 * g_{\min B} \} \f$ for PEMM, and for PEM-BAE* the largest of
 * \f$ \lceil ( b_{\min F} + b_{\min B} ) / 2 \rceil \f$, \f$ g_{\min F} +
 * g_{\min B} \f$, \f$ f_{\min F} + d_{\min B} \f$ and \f$ d_{\min F} +
 * f_{\min B} \f$, with \f$ d_D = g_D - h_{\bar D} \f$ (the bounds of
 * Alcázar et al., AAAI 2020, valid with consistent heuristics, which the
 * PEM-BAE* of Siag et al. computes), all from the keys of the open buckets.
 * The bucket chosen is the one of least priority (f, the priority of MM or
 * the b of BAE*) and of least g among those, the lower-g-first order of
 * Siag et al. (section 3.3), under which, with a consistent heuristic and
 * positive costs, a bucket gets no node once expanded; PEMM takes the
 * direction of least priority and PEM-BAE* alternates them.
 *
 * The states must be trivially copyable, the costs integers, the edge
 * costs positive and the heuristics consistent, as for DDDAStar, which
 * throws std::logic_error otherwise. PEMOptions::disk_limit caps the bytes
 * of the files on the disk: the search stops with Status::disk_limit when
 * they exceed it.
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

#ifndef __PBeFS_PEMSearch
 #define __PBeFS_PEMSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <filesystem>

#include <map>

#include <memory>

#include <mutex>

#include <set>

#include <stdexcept>

#include <string>

#include <tuple>

#include <utility>

#include <vector>

#include "ADS/dictionary/FlatHashTable.h"

#include "ADS/external/RecordFile.h"

#include "PBeFS/astar/Bidirectional.h"

#include "PBeFS/astar/DDDAStar.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the algorithm of a PEMSearch

enum class PEMKind {
 astar ,      ///< PEM-A*, forward only
 mm ,         ///< PEMM
 bae          ///< PEM-BAE*
 };

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS PEMOptions -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of PEMSearch, besides those of all the searches and of the
/// parallel ones

struct PEMOptions {

 /// the partitions of a bucket by the hash of the state
 std::size_t partitions = 16;

 /// where the files go, in a directory of their own removed at the end;
 /// the temporary directory of the system if empty
 std::filesystem::path directory;

 /// the bytes of the buffer of each file being written
 std::size_t buffer_bytes = std::size_t( 1 ) << 16;

 /// the most bytes the files may hold on the disk (0 = no cap)
 std::uint64_t disk_limit = 0;

 };  // end( struct( PEMOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS PEMSearch -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a parallel external-memory search of PEM-BiHS
/** The counters of the Solution: expanded counts the states expanded,
 * expanded_by by direction, generated the children, duplicates the copies
 * dropped, stored the states written to the closed files, sent the bytes
 * written to the disk, and messages the most bytes held on the disk at
 * once. The node_limit of the SearchOptions bounds stored; threads of the
 * ParallelParams is the number of workers. */

template< DiskDomain D , PEMKind K ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
 requires( ( K == PEMKind::astar ) || BidirectionalDomain< D > )
class PEMSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 PEMSearch( const D & domain , unsigned threads = 1 ,
            const SearchOptions & options = {} , PEMOptions pem = {} )
  : f_domain( domain ) , f_threads( std::max( 1u , threads ) ) ,
    f_options( options ) , f_pem( std::move( pem ) ) {
  if( f_pem.partitions == 0 )
   throw( std::invalid_argument( "PEMSearch: no partition" ) );
  if( options.weight != 1 )
   throw( std::invalid_argument( "PEMSearch: no weight" ) );
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers, used for all its buckets

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Stopwatch clock;
  Run r( *this );
  SolutionT sol = r.search();
  sol.stats = r.stats();
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 static constexpr bool c_bi = K != PEMKind::astar;   ///< two directions

 static constexpr Cost inf( void ) noexcept {
  return( infinite_cost< Cost >() );
  }

 /// a record of a file: a state, the hash of its parent and the operator
 /// that generated it
 struct Rec {
  State s;
  std::uint64_t parent;
  [[no_unique_address]] OperatorOf< D > o;
  };

 /// what the table of a partition keeps of a state
 struct Val {
  std::uint64_t parent;
  [[no_unique_address]] OperatorOf< D > o;
  };

 struct Hash {
  const D * d;
  std::uint64_t operator()( const State & s ) const {
   return( d->hash( s ) );
   }
  };

 using Table = ADS::FlatHashTable< State , Val , Hash >;

 /// a bucket: its direction, g and heuristics
 struct Bucket {
  unsigned d;
  Cost g , hf , hb;
  auto tie( void ) const { return( std::tie( d , g , hf , hb ) ); }
  bool operator<( const Bucket & o ) const { return( tie() < o.tie() ); }
  bool operator==( const Bucket & o ) const { return( tie() == o.tie() ); }
  };

 /// a record found on both sides, the best path so far
 struct Meet {
  Rec r[ 2 ];
  Cost g[ 2 ];
  };

 /// the state of a search
 struct Run {

  const PEMSearch & s;
  const D & d;
  ADS::TempDirectory dir;
  std::size_t buffer;
  std::set< Bucket > open;                      ///< the open buckets
  std::map< std::tuple< unsigned , Cost , Cost > , std::vector< Cost > >
   closed;                                      ///< closed g by d , hf , hb
  std::atomic< Cost > U{ inf() };
  std::mutex lock;                              ///< of U and the meeting
  Meet meet{};
  bool have_goal = false;                       ///< PEM-A*: goal record
  Rec goal{};                                   ///< its record
  Cost goal_g = 0;                              ///< and its g
  std::atomic< std::uint64_t > disk{ 0 };       ///< bytes on the disk
  std::uint64_t peak = 0;                       ///< the most of them
  std::atomic< std::uint64_t > written{ 0 };
  std::atomic< std::uint64_t > n_expanded{ 0 } , n_generated{ 0 } ,
                                n_duplicates{ 0 } , n_stored{ 0 };
  std::uint64_t by[ 2 ] = { 0 , 0 };
  unsigned turn = 1;                            ///< to alternate
  std::atomic< bool > bad{ false };             ///< an edge or h refused
  std::unique_ptr< Executor > own;              ///< if none was given
  Executor * ex;                                ///< runs the workers

  explicit Run( const PEMSearch & search )
   : s( search ) , d( search.f_domain ) ,
     dir( search.f_pem.directory , "pbefs-pem-" ) ,
     buffer( std::max< std::size_t >( 1 , search.f_pem.buffer_bytes /
                                           sizeof( Rec ) ) ) ,
     own( search.f_ex ? nullptr
                      : std::make_unique< Executor >( search.f_threads ) ) ,
     ex( search.f_ex ? search.f_ex : own.get() ) {}

  SearchStats stats( void ) const {
   SearchStats st;
   st.expanded = n_expanded.load();
   st.generated = n_generated.load();
   st.duplicates = n_duplicates.load();
   st.stored = n_stored.load();
   st.sent = written.load();
   st.messages = peak;
   st.expanded_by = c_bi ? std::vector< std::uint64_t >{ by[ 0 ] , by[ 1 ] }
                         : std::vector< std::uint64_t >{ by[ 0 ] };
   st.stored_by = { st.stored };
   return( st );
   }

  // the heuristics of s
  void heur( const State & x , Cost & hf , Cost & hb ) const {
   hf = Cost( d.heuristic( x ) );
   if constexpr( c_bi )
    hb = Cost( d.backward_heuristic( x ) );
   else
    hb = Cost( 0 );
   }

  // the file of the partition p of the open (o) or closed (c) bucket b,
  // of the worker w for an open one
  [[nodiscard]] ADS::RecordFile< Rec > file( char kind , const Bucket & b ,
                                             std::size_t p ,
                                             unsigned w = 0 ) const {
   return( ADS::RecordFile< Rec >( dir / ( std::string( 1 , kind ) +
                                           std::to_string( b.d ) + "_" +
                                           std::to_string( b.g ) + "_" +
                                           std::to_string( b.hf ) + "_" +
                                           std::to_string( b.hb ) + "_" +
                                           std::to_string( p ) + "_" +
                                           std::to_string( w ) ) ) );
   }

  [[nodiscard]] std::size_t part( const State & x ) const {
   return( std::size_t( d.hash( x ) % s.f_pem.partitions ) );
   }

  // the priority of a bucket
  [[nodiscard]] Cost key( const Bucket & b ) const {
   const Cost h = b.d == forward_side ? b.hf : b.hb;
   const Cost ho = b.d == forward_side ? b.hb : b.hf;
   if constexpr( K == PEMKind::astar )
    return( b.g + h );
   else if constexpr( K == PEMKind::mm )
    return( std::max( Cost( b.g + h ) , Cost( 2 * b.g ) ) );
   else
    return( Cost( 2 * b.g + h - ho ) );
   }

  // the best bucket of the direction d, of least priority, then g, then
  // h of d, nullptr if none
  const Bucket * best( unsigned side ) const {
   const Bucket * x = nullptr;
   const auto rank = [ & ]( const Bucket & b ) {
    return( std::make_tuple( key( b ) , b.g , side == forward_side ? b.hf
                                                                  : b.hb ) );
    };
   for( const Bucket & b : open )
    if( ( b.d == side ) && ( ! x || ( rank( b ) < rank( *x ) ) ) )
     x = &b;
   return( x );
   }

  // the lower bound on C*
  Cost bound( void ) const {
   Cost kmin[ 2 ] = { inf() , inf() } , fmin[ 2 ] = { inf() , inf() } ,
        gmin[ 2 ] = { inf() , inf() } , dmin[ 2 ] = { inf() , inf() };
   for( const Bucket & b : open ) {
    kmin[ b.d ] = std::min( kmin[ b.d ] , key( b ) );
    fmin[ b.d ] = std::min( fmin[ b.d ] ,
                            Cost( b.g + ( b.d == 0 ? b.hf : b.hb ) ) );
    gmin[ b.d ] = std::min( gmin[ b.d ] , b.g );
    dmin[ b.d ] = std::min( dmin[ b.d ] ,
                            Cost( b.g - ( b.d == 0 ? b.hb : b.hf ) ) );
    }
   [[maybe_unused]] const auto add = []( Cost a , Cost b ) {
    return( ( a == inf() ) || ( b == inf() ) ? inf() : Cost( a + b ) );
    };
   if constexpr( K == PEMKind::astar )
    return( kmin[ 0 ] );
   else if constexpr( K == PEMKind::mm )
    return( std::max( { std::min( kmin[ 0 ] , kmin[ 1 ] ) , fmin[ 0 ] ,
                        fmin[ 1 ] , add( gmin[ 0 ] , gmin[ 1 ] ) } ) );
   else {
    const Cost b = add( kmin[ 0 ] , kmin[ 1 ] );
    return( std::max( { b == inf() ? inf() : Cost( ( b + 1 ) / 2 ) ,
                        add( gmin[ 0 ] , gmin[ 1 ] ) ,
                        add( fmin[ 0 ] , dmin[ 1 ] ) ,
                        add( dmin[ 0 ] , fmin[ 1 ] ) } ) );
    }
   }

  // the root of the direction side
  void root( unsigned side , const State & x ) {
   Cost hf , hb;
   heur( x , hf , hb );
   const Bucket b{ side , Cost( 0 ) , hf , hb };
   const Rec r{ x , d.hash( x ) , {} };
   file( 'o' , b , part( x ) ).append( &r , 1 );
   disk += sizeof( Rec );
   written += sizeof( Rec );
   open.insert( b );
   }

  // offers a path of cost c through the records a (forward) and b
  void offer( Cost c , const Rec & a , Cost ga , const Rec & b , Cost gb ) {
   if( c >= U.load() )
    return;
   std::lock_guard< std::mutex > guard( lock );
   if( c < U.load() ) {
    U.store( c );
    meet = Meet{ { a , b } , { ga , gb } };
    }
   }

  SolutionT search( void ) {
   SolutionT sol;
   const State s0 = d.initial_state();
   if constexpr( c_bi ) {
    const State t = d.goal_state();
    if( s0 == t ) {
     sol.status = Status::solved;
     sol.cost = Cost( 0 );
     sol.bound = sol.cost;
     sol.path = { s0 };
     return( sol );
     }
    root( forward_side , s0 );
    root( backward_side , t );
    }
   else {
    if( d.is_goal( s0 ) ) {
     sol.status = Status::solved;
     sol.cost = Cost( 0 );
     sol.bound = sol.cost;
     sol.path = { s0 };
     return( sol );
     }
    root( forward_side , s0 );
    }
   for( ;; ) {
    const Bucket * nf = best( forward_side );
    const Bucket * nb = c_bi ? best( backward_side ) : nf;
    if( ! nf || ! nb )
     break;
    if( ( U.load() < inf() ) && ( U.load() <= bound() ) )
     break;
    unsigned side = forward_side;
    if constexpr( K == PEMKind::mm )
     side = key( *nb ) < key( *nf ) ? backward_side : forward_side;
    else if constexpr( K == PEMKind::bae )
     side = turn = 1 - turn;
    const Bucket b = side == forward_side ? *nf : *nb;
    open.erase( b );
    expand( b );
    if( bad.load() )
     throw( std::logic_error( "PEMSearch: an edge of a cost that is not "
                              "positive, or a heuristic not consistent" ) );
    peak = std::max( peak , disk.load() );
    if( s.f_pem.disk_limit && ( disk.load() > s.f_pem.disk_limit ) ) {
     sol.status = Status::disk_limit;
     return( sol );
     }
    if( s.f_options.node_limit && ( n_stored.load() >=
                                    s.f_options.node_limit ) ) {
     sol.status = Status::node_limit;
     return( sol );
     }
    }
   if( U.load() == inf() ) {
    sol.status = Status::unsolvable;
    return( sol );
    }
   sol.status = Status::solved;
   sol.cost = U.load();
   sol.bound = sol.cost;
   if constexpr( c_bi ) {
    sol.path = trace( meet.r[ 0 ] , meet.g[ 0 ] , forward_side );
    auto back = trace( meet.r[ 1 ] , meet.g[ 1 ] , backward_side );
    sol.path.insert( sol.path.end() , back.rbegin() + 1 , back.rend() );
    }
   else
    sol.path = trace( goal , goal_g , forward_side );
   return( sol );
   }

  // the closed buckets of the direction side with the heuristics of b
  const std::vector< Cost > * closed_of( unsigned side ,
                                         const Bucket & b ) const {
   const auto it = closed.find( { side , b.hf , b.hb } );
   return( it == closed.end() ? nullptr : &it->second );
   }

  // the open buckets of the direction side with the heuristics of b
  std::vector< Bucket > open_of( unsigned side , const Bucket & b ) const {
   std::vector< Bucket > v;
   for( const Bucket & x : open )
    if( ( x.d == side ) && ( x.hf == b.hf ) && ( x.hb == b.hb ) )
     v.push_back( x );
   return( v );
   }

  // expands the bucket b, its partitions in parallel
  void expand( const Bucket & b ) {
   const unsigned w = s.f_threads;
   const std::size_t P = s.f_pem.partitions;
   const unsigned o = 1 - b.d;
   const std::vector< Cost > * cl = closed_of( b.d , b );
   std::vector< Bucket > oo;                     // the other side's open
   const std::vector< Cost > * oc = nullptr;     // and closed buckets
   if constexpr( c_bi ) {
    oo = open_of( o , b );
    oc = closed_of( o , b );
    }
   const Cost key_b = key( b );
   std::atomic< std::size_t > next{ 0 };
   std::vector< std::set< Bucket > > made( w );
   std::vector< std::uint64_t > exp( w , 0 );
   std::atomic< bool > any{ false };
   const auto work = [ & ]( unsigned i ) {
    using Writer = ADS::RecordWriter< Rec >;
    std::map< std::pair< Bucket , std::size_t > ,
              std::pair< std::unique_ptr< ADS::RecordFile< Rec > > ,
                         std::unique_ptr< Writer > > > out;
    const auto writer = [ & ]( const Bucket & nb , std::size_t p )
     -> Writer & {
     auto & slot = out[ { nb , p } ];
     if( ! slot.second ) {
      slot.first = std::make_unique< ADS::RecordFile< Rec > >(
                                                file( 'o' , nb , p , i ) );
      slot.second = std::make_unique< Writer >( *slot.first , buffer );
      }
     return( *slot.second );
     };
    std::uint64_t gen = 0 , dup = 0 , wrote = 0;
    for( std::size_t p ; ( p = next++ ) < P ; ) {
     std::size_t n = 0;
     std::vector< ADS::RecordFile< Rec > > in;
     for( unsigned k = 0 ; k < w ; ++k ) {
      in.emplace_back( file( 'o' , b , p , k ) );
      n += in.back().size();
      }
     if( n == 0 )
      continue;
     Table t( n , Hash{ &d } );
     for( auto & f : in ) {
      f.forEach( [ & ]( const Rec & x ) {
       if( ! t.insert( x.s , Val{ x.parent , x.o } ) )
        ++dup;
       } );
      disk -= f.size() * sizeof( Rec );
      f.remove();
      }
     if( cl )
      for( Cost gc : *cl )
       if( gc < b.g )
        file( 'c' , Bucket{ b.d , gc , b.hf , b.hb } , p ).forEach(
         [ & ]( const Rec & x ) {
          if( t.remove( x.s ) )
           ++dup;
          } );
     if( t.empty() )
      continue;
     std::vector< Rec > fresh;
     fresh.reserve( t.size() );
     t.forEach( [ & ]( const State & x , const Val & v ) {
      fresh.push_back( Rec{ x , v.parent , v.o } );
      } );
     file( 'c' , b , p ).append( fresh.data() , fresh.size() );
     wrote += fresh.size() * sizeof( Rec );
     n_stored += fresh.size();
     any.store( true );
     if constexpr( c_bi ) {                      // delayed detection
      const auto look = [ & ]( const Bucket & ob ) {
       const auto scan = [ & ]( ADS::RecordFile< Rec > f ) {
        f.forEach( [ & ]( const Rec & x ) {
         if( const Val * v = t.find( x.s ) ) {
          const Rec mine{ x.s , v->parent , v->o };
          if( b.d == forward_side )
           offer( b.g + ob.g , mine , b.g , x , ob.g );
          else
           offer( b.g + ob.g , x , ob.g , mine , b.g );
          }
         } );
        };
       return( scan );
       };
      for( const Bucket & ob : oo )
       for( unsigned k = 0 ; k < w ; ++k )
        look( ob )( file( 'o' , ob , p , k ) );
      if( oc )
       for( Cost gc : *oc ) {
        const Bucket ob{ o , gc , b.hf , b.hb };
        look( ob )( file( 'c' , ob , p ) );
        }
      }
     for( const Rec & x : fresh ) {
      ++exp[ i ];
      const std::uint64_t hu = d.hash( x.s );
      const Cost hd = b.d == forward_side ? b.hf : b.hb;
      const auto child = [ & ]( const State & v , Cost c , Cost hv ,
                                std::uint8_t ov ) {
       ++gen;
       if( c <= Cost( 0 ) ) {
        bad.store( true );
        return;
        }
       const Cost g = b.g + c;
       if( g + hv >= U.load() )
        return;
       Cost hf , hb;
       if( b.d == forward_side ) {
        hf = hv;
        hb = Cost( 0 );
        if constexpr( c_bi )
         hb = Cost( d.backward_heuristic( v ) );
        }
       else {
        hb = hv;
        hf = Cost( d.heuristic( v ) );
        }
       const Bucket nb{ b.d , g , hf , hb };
       if( key( nb ) < key_b ) {                 // the order would break
        bad.store( true );
        return;
        }
       const Rec r{ v , hu , ov };
       if constexpr( ! c_bi )
        if( d.is_goal( v ) ) {                   // detected at generation
         std::lock_guard< std::mutex > guard( lock );
         if( g < U.load() ) {
          U.store( g );
          goal = r;
          goal_g = g;
          have_goal = true;
          }
         return;
         }
       writer( nb , part( v ) ).push( r );
       wrote += sizeof( Rec );
       made[ i ].insert( nb );
       };
      if( b.d == forward_side )
       for_each_child( d , x.s , hd , x.o.op , child );
      else
       if constexpr( c_bi )
        for_each_parent( d , x.s , hd , x.o.op , child );
      }
     }
    for( auto & [ k , fw ] : out )
     fw.second->flush();
    n_generated += gen;
    n_duplicates += dup;
    written += wrote;
    disk += wrote;
    };
   run_workers( ex , w , work );
   for( unsigned i = 0 ; i < w ; ++i ) {
    n_expanded += exp[ i ];
    by[ b.d ] += exp[ i ];
    for( const Bucket & nb : made[ i ] )
     open.insert( nb );
    }
   if( any.load() ) {
    auto & v = closed[ { b.d , b.hf , b.hb } ];
    if( std::find( v.begin() , v.end() , b.g ) == v.end() )
     v.push_back( b.g );
    }
   }

  // the path from the root of the direction side to the record x of g,
  // each parent looked for in the closed files of its partition, of
  // smaller g, with the child as a successor (or predecessor) at the cost
  // of the difference
  std::vector< State > trace( Rec x , Cost g , unsigned side ) {
   std::vector< State > path{ x.s };
   const State r = side == forward_side ? d.initial_state()
                                        : goal_state();
   while( ! ( ( g == Cost( 0 ) ) && ( x.s == r ) ) ) {
    const std::size_t p = std::size_t( x.parent % s.f_pem.partitions );
    bool found = false;
    for( const auto & [ k , gs ] : closed ) {
     if( std::get< 0 >( k ) != side )
      continue;
     const Cost hd = side == forward_side ? std::get< 1 >( k )
                                          : std::get< 2 >( k );
     for( Cost gc : gs ) {
      if( gc >= g )
       continue;
      file( 'c' , Bucket{ side , gc , std::get< 1 >( k ) ,
                          std::get< 2 >( k ) } , p ).forEach(
       [ & ]( const Rec & y ) {
        if( found || ( d.hash( y.s ) != x.parent ) )
         return;
        const auto test = [ & ]( const State & v , Cost c , Cost ) {
         if( ! found && ( v == x.s ) && ( gc + c == g ) ) {
          found = true;
          x = y;
          g = gc;
          }
         };
        if( side == forward_side )
         d.for_each_successor( y.s , hd , test );
        else
         if constexpr( c_bi )
          d.for_each_predecessor( y.s , hd , test );
        } );
      if( found )
       break;
      }
     if( found )
      break;
     }
    if( ! found )
     throw( std::logic_error( "PEMSearch: the parent of a state is not in "
                              "the closed files" ) );
    path.push_back( x.s );
    }
   std::reverse( path.begin() , path.end() );
   return( path );
   }

  State goal_state( void ) const {
   if constexpr( c_bi )
    return( d.goal_state() );
   else
    return( d.initial_state() );
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 unsigned f_threads;                                 ///< the workers
 SearchOptions f_options;                            ///< the options
 PEMOptions f_pem;                                   ///< the options of PEM
 Executor * f_ex = nullptr;                          ///< given, or nullptr

 };  // end( class( PEMSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PEMSearch */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File PEMSearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
