/*--------------------------------------------------------------------------*/
/*-------------------------- File TraceRecorder.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template TraceRecorder, an observer of the
 * searches [see SearchObserver.h] that records their events in a buffer
 * per worker and writes them as a trace, one JSON object per line.
 *
 * Each worker appends its events to its own buffer, which no other worker
 * touches, so that recording takes no lock and no atomic operation; the
 * buffers are made by on_start(), before the workers start, and merged
 * after the search. The time of an event is its step, i.e., its index
 * among the events of its worker, and not the wall-clock time, so that
 * the trace of a sequential search is the same at every run; the merged
 * trace orders the events by step and then by worker, which interleaves
 * the workers one event at a time and is not the order in which they
 * ran. A cap bounds the events recorded, divided evenly among the
 * workers: the events of a worker beyond its share are counted as
 * dropped, and their steps are skipped. A mask chooses the kinds of
 * events recorded, and every records only the events whose step is a
 * multiple of it.
 *
 * The states and the nodes are written by an id, a word given by the
 * function Id of the recorder (e.g., the hash of the Domain), 0 for the
 * objects it does not take; the trace thus does not depend on the type of
 * the states, and a reader that wants them decodes the ids.
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

#ifndef __PBeFS_TraceRecorder
 #define __PBeFS_TraceRecorder
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <iomanip>

#include <limits>

#include <ostream>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/core/SearchObserver.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the kinds of the events of a trace, one per hook [see SearchObserver.h]

enum class TraceEvent : std::uint8_t {
 expand ,       ///< a state expanded
 generate ,     ///< a state generated
 duplicate ,    ///< a state generated and already known
 open_push ,    ///< an entry pushed in OPEN
 open_pop ,     ///< an entry popped from OPEN
 message ,      ///< nodes handed to another worker
 threshold ,    ///< an iteration of an iterative-deepening search
 solution ,     ///< a goal, or a feasible value, that improves the
                ///< incumbent
 evaluate ,     ///< a node of B&B evaluated
 branch ,       ///< a node of B&B branched
 prune          ///< a node of B&B pruned
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the bit of the event kind e in a mask of event kinds

[[nodiscard]] constexpr std::uint32_t trace_bit( TraceEvent e ) noexcept
{
 return( std::uint32_t( 1 ) << unsigned( e ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the mask of all the event kinds

inline constexpr std::uint32_t trace_all = ( std::uint32_t( 1 ) << 11 ) - 1;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the name of the event kind e in a trace

[[nodiscard]] constexpr const char * trace_name( TraceEvent e ) noexcept
{
 constexpr const char * names[] = { "expand" , "generate" , "duplicate" ,
                                    "push" , "pop" , "message" ,
                                    "threshold" , "solution" , "evaluate" ,
                                    "branch" , "prune" };
 return( names[ unsigned( e ) ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS TraceRecord ----------------------------*/
/*--------------------------------------------------------------------------*/
/// an event of a trace
/** The fields that an event does not use are 0. By kind:
 *
 * - expand: id, g, h and f of the state;
 * - generate: id, g and h of the state, parent the id of its parent;
 * - duplicate: id and g of the state;
 * - open_push, open_pop: size the entries of OPEN after the event;
 * - message: other the worker that receives, size the nodes;
 * - threshold: g the threshold;
 * - solution: g the cost or the value;
 * - evaluate, branch, prune: g the bound, other the depth, id the id of
 *   the node (not for prune), size the children (for branch). */

struct TraceRecord {

 std::uint64_t step = 0;     ///< the index of the event in its worker
 std::uint64_t id = 0;       ///< the id of the state or of the node
 std::uint64_t parent = 0;   ///< the id of the parent
 std::uint64_t size = 0;     ///< entries, nodes or children
 double g = 0;               ///< g, or the bound, threshold or value
 double h = 0;               ///< h
 double f = 0;               ///< the key of OPEN
 std::uint32_t worker = 0;   ///< the worker
 std::uint32_t other = 0;    ///< the receiver or the depth
 TraceEvent type = TraceEvent::expand;  ///< the kind of the event

 };  // end( struct( TraceRecord ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the id function of a TraceRecorder that gives no id

struct TraceNoId {
 };  // end( struct( TraceNoId ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS TraceRecorder ---------------------------*/
/*--------------------------------------------------------------------------*/
/// an observer that records the events of a search, a buffer per worker
/** Id is called as id( x ) on the states or the nodes x of the events,
 * where it can take them, and must return their id as a word; the
 * recorder must outlive the searches that it observes, which may be more
 * than one, one after the other, each on_start() clearing what the
 * previous one recorded. */

template< class Id = TraceNoId >
class TraceRecorder {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the id function, the cap on the events recorded
 /// (0 = none), the mask of their kinds and the sampling period every
 /// (an event of each every, by step; 0 is 1)

 explicit TraceRecorder( Id id = Id() , std::size_t cap = 1 << 20 ,
                         std::uint32_t mask = trace_all ,
                         std::uint64_t every = 1 )
  : f_id( std::move( id ) ) , f_cap( cap ) , f_mask( mask ) ,
    f_every( every ? every : 1 ) {}

/*------------------------------- THE HOOKS --------------------------------*/

 /// the buffers of p workers, empty, before the search
 void on_start( unsigned p ) {
  if( p == 0 )
   p = 1;
  v_buf.clear();
  v_buf.resize( p );
  f_share = f_cap ? ( f_cap + p - 1 ) / p
                  : std::numeric_limits< std::size_t >::max();
  for( auto & b : v_buf )
   b.v.reserve( std::min< std::size_t >( f_share , 4096 ) );
  }

 /// nothing, the buffers being merged when the trace is read
 void on_finish( void ) noexcept {}

 /// [see SearchObserver.h]
 template< class S , class C >
 void on_expand( unsigned i , const S & s , C g , C h , C f ) {
  if( TraceRecord * r = slot( i , TraceEvent::expand ) ) {
   r->id = id_of( s );
   r->g = value( g );
   r->h = value( h );
   r->f = value( f );
   }
  }

 /// [see SearchObserver.h]
 template< class S , class C >
 void on_generate( unsigned i , const S & s , const S * parent , C g ,
                   C h ) {
  if( TraceRecord * r = slot( i , TraceEvent::generate ) ) {
   r->id = id_of( s );
   r->parent = parent ? id_of( *parent ) : 0;
   r->g = value( g );
   r->h = value( h );
   }
  }

 /// [see SearchObserver.h]
 template< class S , class C >
 void on_duplicate( unsigned i , const S & s , C g ) {
  if( TraceRecord * r = slot( i , TraceEvent::duplicate ) ) {
   r->id = id_of( s );
   r->g = value( g );
   }
  }

 /// [see SearchObserver.h]
 void on_open_push( unsigned i , std::size_t n ) {
  if( TraceRecord * r = slot( i , TraceEvent::open_push ) )
   r->size = n;
  }

 /// [see SearchObserver.h]
 void on_open_pop( unsigned i , std::size_t n ) {
  if( TraceRecord * r = slot( i , TraceEvent::open_pop ) )
   r->size = n;
  }

 /// [see SearchObserver.h]
 void on_message( unsigned i , unsigned j , std::size_t n ) {
  if( TraceRecord * r = slot( i , TraceEvent::message ) ) {
   r->other = j;
   r->size = n;
   }
  }

 /// [see SearchObserver.h]
 template< class C >
 void on_threshold( unsigned i , C b ) {
  if( TraceRecord * r = slot( i , TraceEvent::threshold ) )
   r->g = value( b );
  }

 /// [see SearchObserver.h]
 template< class C >
 void on_solution( unsigned i , C c ) {
  if( TraceRecord * r = slot( i , TraceEvent::solution ) )
   r->g = value( c );
  }

 /// [see SearchObserver.h]
 template< class N , class V >
 void on_evaluate( unsigned i , const N & n , V b , unsigned d ) {
  if( TraceRecord * r = slot( i , TraceEvent::evaluate ) ) {
   r->id = id_of( n );
   r->g = value( b );
   r->other = d;
   }
  }

 /// [see SearchObserver.h]
 template< class N , class V >
 void on_branch( unsigned i , const N & n , V b , unsigned d ,
                 std::size_t k ) {
  if( TraceRecord * r = slot( i , TraceEvent::branch ) ) {
   r->id = id_of( n );
   r->g = value( b );
   r->other = d;
   r->size = k;
   }
  }

 /// [see SearchObserver.h]
 template< class V >
 void on_prune( unsigned i , V b , unsigned d ) {
  if( TraceRecord * r = slot( i , TraceEvent::prune ) ) {
   r->g = value( b );
   r->other = d;
   }
  }

/*-------------- METHODS FOR READING THE DATA OF THE TRACE -----------------*/

 /// the number of the workers of the last search
 [[nodiscard]] unsigned workers( void ) const noexcept {
  return( unsigned( v_buf.size() ) );
  }

 /// the events recorded by the worker i, by step
 [[nodiscard]] const std::vector< TraceRecord > & events( unsigned i )
  const {
  return( v_buf[ i ].v );
  }

 /// the events recorded, of all the workers
 [[nodiscard]] std::size_t size( void ) const noexcept {
  std::size_t n = 0;
  for( const auto & b : v_buf )
   n += b.v.size();
  return( n );
  }

 /// the events dropped by the cap, of all the workers
 [[nodiscard]] std::uint64_t dropped( void ) const noexcept {
  std::uint64_t n = 0;
  for( const auto & b : v_buf )
   n += b.dropped;
  return( n );
  }

 /// the events of kind e recorded by the worker i, or by all of them if i
 /// is the number of the workers
 [[nodiscard]] std::uint64_t count( TraceEvent e , unsigned i ) const {
  std::uint64_t n = 0;
  for( unsigned j = 0 ; j < v_buf.size() ; ++j )
   if( ( i == v_buf.size() ) || ( i == j ) )
    n += std::uint64_t( std::count_if( v_buf[ j ].v.begin() ,
                                       v_buf[ j ].v.end() ,
                                       [ e ]( const TraceRecord & r ) {
                                        return( r.type == e );
                                        } ) );
  return( n );
  }

 /// the events of all the workers, by step and then by worker
 [[nodiscard]] std::vector< TraceRecord > merged( void ) const {
  std::vector< TraceRecord > all;
  all.reserve( size() );
  for( const auto & b : v_buf )
   all.insert( all.end() , b.v.begin() , b.v.end() );
  std::stable_sort( all.begin() , all.end() ,
                    []( const TraceRecord & a , const TraceRecord & b ) {
                     return( a.step < b.step );
                     } );
  return( all );
  }

 /// writes the trace to out: a header line, with the number of the
 /// workers, of the events recorded and dropped and the cap, then the
 /// merged events, one per line, each with the fields of its kind [see
 /// TraceRecord], the ids in hexadecimal (a word does not fit the numbers
 /// of JSON) and the values not finite as null

 void write( std::ostream & out ) const {
  out << "{\"trace\":\"PBeFS\",\"version\":1,\"workers\":" << workers()
      << ",\"events\":" << size() << ",\"dropped\":" << dropped()
      << ",\"cap\":" << f_cap << ",\"every\":" << f_every << "}\n";
  for( const TraceRecord & r : merged() )
   write( out , r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes the event r to out as a line of JSON

 static void write( std::ostream & out , const TraceRecord & r ) {
  out << "{\"e\":\"" << trace_name( r.type ) << "\",\"w\":" << r.worker
      << ",\"t\":" << r.step;
  const auto id = [ & ]( const char * k , std::uint64_t x ) {
   out << ",\"" << k << "\":\"" << std::hex << x << std::dec << "\"";
   };
  const auto num = [ & ]( const char * k , double x ) {
   out << ",\"" << k << "\":";
   if( ! std::isfinite( x ) )
    out << "null";
   else
    if( ( x == std::floor( x ) ) && ( std::abs( x ) < 9.0e15 ) )
     out << std::int64_t( x );
    else
     out << std::setprecision( 17 ) << x << std::setprecision( 6 );
   };
  switch( r.type ) {
   case( TraceEvent::expand ):
    id( "id" , r.id ); num( "g" , r.g ); num( "h" , r.h ); num( "f" , r.f );
    break;
   case( TraceEvent::generate ):
    id( "id" , r.id ); id( "parent" , r.parent ); num( "g" , r.g );
    num( "h" , r.h );
    break;
   case( TraceEvent::duplicate ):
    id( "id" , r.id ); num( "g" , r.g );
    break;
   case( TraceEvent::open_push ): case( TraceEvent::open_pop ):
    out << ",\"size\":" << r.size;
    break;
   case( TraceEvent::message ):
    out << ",\"to\":" << r.other << ",\"size\":" << r.size;
    break;
   case( TraceEvent::threshold ): case( TraceEvent::solution ):
    num( "g" , r.g );
    break;
   case( TraceEvent::evaluate ):
    id( "id" , r.id ); num( "bound" , r.g ); out << ",\"depth\":" << r.other;
    break;
   case( TraceEvent::branch ):
    id( "id" , r.id ); num( "bound" , r.g ); out << ",\"depth\":" << r.other
                                                 << ",\"size\":" << r.size;
    break;
   case( TraceEvent::prune ):
    num( "bound" , r.g ); out << ",\"depth\":" << r.other;
    break;
   }
  out << "}\n";
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the buffer of a worker, on cache lines of its own
 struct alignas( cache_line ) Buffer {
  std::vector< TraceRecord > v;               // the events recorded
  std::uint64_t step = 0;                     // the steps so far
  std::uint64_t dropped = 0;                  // the events dropped
  };

 // the record of the next event of kind e of the worker i, its step and
 // kind set, or nullptr if the event is not recorded
 TraceRecord * slot( unsigned i , TraceEvent e ) {
  if( ! ( f_mask & trace_bit( e ) ) )
   return( nullptr );
  Buffer & b = v_buf[ i ];
  const std::uint64_t t = b.step++;
  if( t % f_every )
   return( nullptr );
  if( b.v.size() >= f_share ) {
   ++b.dropped;
   return( nullptr );
   }
  TraceRecord & r = b.v.emplace_back();
  r.step = t;
  r.worker = i;
  r.type = e;
  return( & r );
  }

 // the id of x, 0 if Id does not take it
 template< class X >
 std::uint64_t id_of( const X & x ) const {
  if constexpr( std::is_invocable_v< const Id & , const X & > )
   return( std::uint64_t( f_id( x ) ) );
  else
   return( 0 );
  }

 // a cost or a value as a double, NaN if it is not convertible
 template< class C >
 static double value( const C & c ) {
  if constexpr( std::is_convertible_v< C , double > )
   return( double( c ) );
  else
   return( std::numeric_limits< double >::quiet_NaN() );
  }

 Id f_id;                                     // the id function
 std::size_t f_cap;                           // the cap on the events
 std::uint32_t f_mask;                        // the kinds recorded
 std::uint64_t f_every;                       // the sampling period
 std::size_t f_share = 0;                     // the cap of a worker
 std::vector< Buffer > v_buf;                 // the buffers

 };  // end( class( TraceRecorder ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_TraceRecorder */

/*--------------------------------------------------------------------------*/
/*------------------------ End File TraceRecorder.h ------------------------*/
/*--------------------------------------------------------------------------*/
