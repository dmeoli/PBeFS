/*--------------------------------------------------------------------------*/
/*-------------------------- File ThreadMailbox.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for ThreadMailbox, the messages of the searches whose workers
 * run on a runtime of shared memory [see ThreadRuntime.h] and send their
 * nodes to the worker that owns them: the batches of the outboxes, the
 * inboxes of the channels of the runtime, and the detection of the end.
 *
 * The end of a search whose workers exchange messages is detected with three
 * shared counters, as follows (cf. the page "Parallel A*" of the
 * documentation): pending counts the messages sent and not yet processed,
 * idle the workers waiting, and epoch the times a worker stopped waiting;
 * \f$ p \f$ is the number of workers:
 *
 * \code{.unparsed}
 *  1  procedure SEND(B, j)                 ▷ the batch B, to the worker j
 *  2     pending ← pending + |B|          ▷ before B becomes visible
 *  3     append B to the inbox of j
 *  4  procedure PROCESSED(n)               ▷ after n messages are processed
 *  5     pending ← pending − n
 *  6  procedure WAIT(i)                    ▷ the worker i has nothing to do
 *  7     idle ← idle + 1
 *  8     loop
 *  9        if done then return
 * 10        if the inbox of i is not empty then
 * 11           idle ← idle − 1;  epoch ← epoch + 1;  return
 * 12        e ← epoch
 * 13        if idle = p and pending = 0 and epoch = e then
 * 14           done ← true;  return
 * 15        wait a little
 * \endcode
 *
 * The test at line 13 is sound: at the instant at which idle is read every
 * worker is waiting, and pending \f$ = 0 \f$ read after it means that every
 * message sent before that instant has been processed; a message sent later
 * would have to come from a worker that stopped waiting, which raises epoch
 * (line 11) after leaving idle, and then epoch differs from \f$ e \f$. A
 * send that does not wait for the lock of the inbox and finds it taken
 * lowers pending again and keeps the batch in its outbox.
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

#ifndef __PBeFS_ThreadMailbox
 #define __PBeFS_ThreadMailbox
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <vector>

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/core/ThreadRuntime.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ThreadMailbox ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the messages of the workers of a search in shared memory
/** The workers send messages of type Msg (say, nodes) to each other in
 * batches, through the channels of the Runtime [see ThreadRuntime], and
 * the end of the search is detected with three shared counters: pending,
 * the messages sent and not yet processed (raised before a batch becomes
 * visible, lowered after it has been processed [see processed()]), idle,
 * the workers waiting, and epoch, raised by a worker that stops waiting
 * after it has left idle. A waiting worker declares the end when epoch
 * does not change while it reads idle == p and pending == 0: at the moment
 * idle was read nobody had work, and no message could be in transit, since
 * consuming it would have required somebody to stop waiting, i.e., to
 * change epoch [see wait()]. The search may also be ended at once, say by
 * a limit [see stop()].
 *
 * Each worker has its Outboxes, one per destination, a batch being sent,
 * if the lock of the inbox of the destination is free, once it holds
 * batch messages, and all of them before the worker waits [see flush()]. */

template< class Msg , class Runtime = ThreadRuntime >
class ThreadMailbox {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using Channels = typename Runtime::template channels_type< Msg >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constructor, for p workers

 explicit ThreadMailbox( unsigned p ) : f_p( p ) , f_channels( p ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of workers

 [[nodiscard]] unsigned workers( void ) const noexcept { return( f_p ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves the messages b into the inbox of dst before the workers start
 /// (say, the root), counting them as pending; b is left empty

 void seed( unsigned dst , std::vector< Msg > & b ) {
  if( b.empty() )
   return;
  f_pending.fetch_add( std::int64_t( b.size() ) );
  f_channels.send( 0 , dst , b );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves the inbox of the worker i into in, waiting for the lock if
 /// block, and returns the number of messages; the worker calls
 /// processed() once it has processed them

 std::size_t receive( unsigned i , std::vector< Msg > & in , bool block ) {
  return( f_channels.receive( i , in , block ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// records that n messages received have been processed

 void processed( std::size_t n ) noexcept {
  if( n )
   f_pending.fetch_sub( std::int64_t( n ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the search is over

 [[nodiscard]] bool done( void ) const noexcept {
  return( f_done.load( std::memory_order_acquire ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// ends the search at once, whatever the messages left (say, by a limit)

 void stop( void ) noexcept {
  f_done.store( true , std::memory_order_release );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// waits until mail arrives for the worker i or the search ends, which
 /// it may detect itself [see the class comment]

 void wait( unsigned i ) {
  f_idle.fetch_add( 1 );
  Backoff backoff;
  for( ;; ) {
   if( done() )
    return;
   if( f_channels.has_mail( i ) ) {
    f_idle.fetch_sub( 1 );     // first leave idle, then change epoch
    f_epoch.fetch_add( 1 );
    return;
    }
   const auto e = f_epoch.load();
   if( ( f_idle.load() == f_p ) && ( f_pending.load() == 0 ) &&
       ( f_epoch.load() == e ) ) {
    stop();
    return;
    }
   backoff.pause();
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the outboxes of a worker

 class Outboxes {

  public:

  /// constructor, for the worker id of the mailbox m, with batches of
  /// batch messages
  Outboxes( ThreadMailbox & m , unsigned id , std::size_t batch )
   : f_m( m ) , f_id( id ) , f_batch( batch ) , v_out( m.workers() ) {
   for( auto & o : v_out )
    o.reserve( f_batch );
   }

  /// adds the message x to the outbox of dst, which is sent if it is full
  /// and the lock of the inbox of dst is free
  void push( unsigned dst , Msg && x ) {
   auto & o = v_out[ dst ];
   o.push_back( std::move( x ) );
   if( o.size() >= f_batch )
    post( dst , false );
   }

  /// sends the non-empty outboxes, waiting for the locks if block
  void flush( bool block ) {
   for( unsigned d = 0 ; d < f_m.workers() ; ++d )
    if( ! v_out[ d ].empty() )
     post( d , block );
   }

  /// the batches sent so far
  [[nodiscard]] std::uint64_t messages( void ) const noexcept {
   return( f_messages );
   }

  /// the outboxes, one per destination, with what is still unsent
  [[nodiscard]] std::vector< std::vector< Msg > > & boxes( void ) noexcept {
   return( v_out );
   }

  private:

  // sends the outbox of dst, waiting for the lock if block; true if sent
  bool post( unsigned dst , bool block ) {
   auto & o = v_out[ dst ];
   const auto k = std::int64_t( o.size() );
   f_m.f_pending.fetch_add( k );
   if( block )
    f_m.f_channels.send( f_id , dst , o );
   else
    if( ! f_m.f_channels.try_send( f_id , dst , o ) ) {
     f_m.f_pending.fetch_sub( k );
     return( false );
     }
   ++f_messages;
   o.reserve( f_batch );
   return( true );
   }

  ThreadMailbox & f_m;                          ///< the mailbox
  unsigned f_id;                                ///< the worker
  std::size_t f_batch;                          ///< messages in a batch
  std::vector< std::vector< Msg > > v_out;      ///< one per destination
  std::uint64_t f_messages = 0;                 ///< batches sent

  };  // end( class( Outboxes ) )

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 unsigned f_p;                                  ///< the number of workers
 Channels f_channels;                           ///< the inboxes
 alignas( cache_line ) std::atomic< std::int64_t > f_pending{ 0 };
 alignas( cache_line ) std::atomic< unsigned > f_idle{ 0 };
 alignas( cache_line ) std::atomic< std::uint64_t > f_epoch{ 0 };
 alignas( cache_line ) std::atomic< bool > f_done{ false };

 };  // end( class( ThreadMailbox ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ParallelParams ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the parameters shared by the parallel searches
/** The number of workers, and, for the searches that send their nodes to
 * their owners, the messages in a batch, the steps (expansions or
 * evaluations) between two flushes of the outboxes and the seed of the
 * owners; each search sets its own defaults of batch and flush_period. */

struct ParallelParams {

 unsigned threads = 1;          ///< the number of workers
 std::size_t batch = 64;        ///< nodes packed in a message
 unsigned flush_period = 64;    ///< steps between flushes of outboxes
 std::uint64_t seed = 1;        ///< seed of the owners

 /// sets to 1 any parameter that is 0
 void normalize( void ) noexcept {
  if( threads == 0 )
   threads = 1;
  if( batch == 0 )
   batch = 1;
  if( flush_period == 0 )
   flush_period = 1;
  }

 };  // end( struct( ParallelParams ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ThreadMailbox */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ThreadMailbox.h ------------------------*/
/*--------------------------------------------------------------------------*/
