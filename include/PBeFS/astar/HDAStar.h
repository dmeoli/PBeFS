/*--------------------------------------------------------------------------*/
/*----------------------------- File HDAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template HDAStar, Hash Distributed A* on a
 * shared-memory machine, whose workers run on an executor of ADS and
 * exchange nodes through the mailbox of ADS on its channels.
 *
 * The search follows section 3 of Kishimoto, Fukunaga and Botea (AIJ 2013)
 * and section 4.1 of Fukunaga et al. (2017), with the notation of the page
 * "Best-first search and A*" of the documentation; the worker \f$ i \f$ owns
 * the states \f$ v \f$ with \f$ \mathrm{owner}( v ) = i \f$, which only it
 * stores, in its lists \f$ \mathrm{OPEN}_i \f$ and
 * \f$ \mathrm{CLOSED}_i \f$, the incumbent \f$ z \f$ is atomic, and a
 * message \f$ ( v , g' , u ) \f$ carries a state with the cost and the
 * predecessor of the path that reached it:
 *
 * \code{.unparsed}
 *  1  procedure HDA*(s₀, p, owner)
 *  2     z ← ∞;  INSERT(s₀, 0, nil) on owner(s₀)
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        while the search is not over do
 *  5           for each message (v, g', u) received do INSERT(v, g', u)
 *  6           if OPEN_i has a node of f < z then
 *  7              u ← the first node of OPEN_i, removed (stale entries
 *                     dropped);  close u
 *  8              if u is a goal then
 *  9                 z ← min(z, g(u)), recording u if g(u) < z
 * 10              else for each successor v of u do
 * 11                 g' ← g(u) + c(u, v)
 * 12                 if g' + h(v) ≥ z then continue
 * 13                 if owner(v) = i then INSERT(v, g', u)
 * 14                 else send (v, g', u) to owner(v)
 * 15              every k expansions, send the outboxes not yet full
 * 16           else
 * 17              send all the outboxes;  wait for mail or for the end
 * 18     return the path to the goal of cost z, if z < ∞
 * 19  procedure INSERT(v, g', u)                       ▷ on the owner i of v
 * 20     if g' + h(v) ≥ z then return
 * 21     if v ∉ CLOSED_i then
 * 22        add v to CLOSED_i, with g(v) ← g' and parent(v) ← u
 * 23     else if g' < g(v) then
 * 24        g(v) ← g';  parent(v) ← u;  v is not closed
 * 25     else return                                         ▷ a duplicate
 * 26     push (g' + h(v), g', v) in OPEN_i
 * \endcode
 *
 * The messages travel in batches of batch nodes, and an outbox is sent when
 * it is full (line 14), every flush_period expansions (line 15) and before the
 * worker waits (line 17); the wait ends when mail arrives or when the
 * mailbox detects the end, i.e., when no worker has a node of
 * \f$ f < z \f$ and no message is in transit, which is the condition of
 * the page "Parallel A*" for an optimal \f$ z \f$. The parents
 * cross the workers, and the path is followed only after all the workers
 * have finished.
 *
 * With a weight \f$ w > 1 \f$ [see SearchOptions and Weight.h] OPEN is
 * ordered by \f$ f' = g + w h \f$, the test of line 6 asks for a node of
 * \f$ f' < z \f$, and lines 12 and 20 prune the nodes of \f$ w f \geq z \f$
 * (\f$ f' \geq z \f$ with duplicate dropping, whose test replaces that of
 * line 23): the search ends with \f$ z \leq w C^* \f$, whatever the order
 * in which the workers expand. The anytime search (anytime weighted HDA*)
 * tests at line 6 that OPEN is not empty, drops at line 7 the nodes of
 * \f$ f \geq z \f$, and prunes those of \f$ f \geq z \f$ at lines 12 and
 * 20, so that it ends with an optimal \f$ z \f$ after recording the
 * incumbents found on the way.
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

#ifndef __PBeFS_HDAStar
 #define __PBeFS_HDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <memory>

#include <mutex>

#include <optional>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/core/SearchUtils.h"

#include "ADS/parallel/Channels.h"

#include "ADS/parallel/Mailbox.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ParallelOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the parallel searches: those of the sequential search,
/// and the parameters of the parallel ones [see ParallelParams], with
/// batches of 64 nodes flushed every 64 expansions

struct ParallelOptions : public SearchOptions , public ParallelParams {
 };  // end( struct( ParallelOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS HDAOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of HDAStar: those of the parallel searches, the synchronous
/// sends of PRA* and the lending of the nodes to the workers left behind

struct HDAOptions : public ParallelOptions {

 /// true to push every node to its owner as soon as it is generated, as
 /// in PRA*, instead of packing the nodes in batches
 bool sync_send = false;

 /// if positive, a worker lends the node it is about to expand to the
 /// worker whose best f is the least among those not smaller than its f
 /// plus outsource, which expands it at once; 0 for none
 double outsource = 0;

 HDAOptions( void ) = default;

 /// the options of the parallel searches, with no synchronous send and no
 /// lending
 HDAOptions( const ParallelOptions & o ) : ParallelOptions( o ) {}

 };  // end( struct( HDAOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS LockedChannelsOf -------------------------*/
/*--------------------------------------------------------------------------*/
/// the channels of HDAStar whose inboxes are vectors under the lock of the
/// Executor, which a worker takes, waiting for it, to receive [see
/// ADS::LockedChannels]

template< ADS::WorkerExecutor Executor >
struct LockedChannelsOf {
 /// the family for the messages Msg
 template< class Msg >
 using type = ADS::LockedChannels< Msg , typename Executor::lock_type >;
 };  // end( struct( LockedChannelsOf ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS HDAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Hash Distributed A* on a shared-memory machine
/** Every worker runs A* on its own OPEN and CLOSED lists, and a global
 * function Dist assigns every state to the worker that owns it: a node
 * generated by a worker is sent to its owner, which checks it against its
 * closed list and possibly adds it to its OPEN list [see section 4.1 of
 * Fukunaga et al., 2017, and Kishimoto, Fukunaga and Botea, AIJ 2013].
 * Nodes are sent asynchronously, packed in batches of batch nodes, and an
 * outbox not yet full is flushed every flush_period expansions and before
 * the worker waits. The goal test is done at the expansion, the cost of
 * the best goal found (the incumbent) prunes every node with f not
 * smaller, and the search ends when no worker has a node with f smaller
 * than the incumbent and no node is travelling, which the mailbox of ADS
 * detects.
 *
 * The node of a state belongs to its owner, which is the only worker that
 * writes it; parent pointers cross the workers, and they are followed only
 * after all the workers have finished.
 *
 * The channels between the workers are the family ChannelsOf, by default
 * the one of ADS for the Executor, which never waits; with
 * LockedChannelsOf< Executor > an inbox is a vector under a lock, which the
 * receiver takes waiting for it, and with the option sync_send every node
 * is pushed alone as soon as it is generated. The two together are PRA*,
 * the Parallel Retracting A* of Evett et al. (1995) without the retraction,
 * whose communication is synchronous at both ends, and with the
 * AbstractStateDistribution they are APRA* [see Burns et al., JAIR 2010,
 * who compare the four combinations of synchronous and asynchronous sends
 * and receives]. With the option outsource positive, a worker publishes
 * the f of the node it is about to expand, and lends the node to the
 * worker whose published f is the least among those at least outsource
 * above it, if there is one: the borrower expands it at once and sends its
 * successors to their owners, while the node stays in the CLOSED list of
 * its owner; this moves work from a worker whose nodes are much better
 * than those of the others, which the distribution alone cannot do. A lent
 * node is never lent again. */

template< Domain D , class Dist = StateHashDistribution< D > ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > ,
          class ChannelsOf = ADS::ChannelsOf< Executor > >
class HDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// a node travelling to its owner
 struct Msg {
  State s;
  Cost g;
  Cost h;
  NodeT * parent;
  [[no_unique_address]] OperatorOf< D > o;        // the operator of s
  std::uint32_t key;                               // the key of s
  bool lent = false;            // s is lent, parent being its own node
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain and the distribution (which must
 /// outlive the HDAStar)

 HDAStar( const D & domain , const Dist & dist ,
          const HDAOptions & options = {} )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) ,
    f_w( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  Shared sh( p , f_options.sync_send ? 1 : f_options.batch , clock );

  // the root, which its owner inserts before the first expansion
  const State s0 = f_domain.initial_state();
  sh.root.emplace( Msg{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) ,
                        nullptr , no_op ,
                        distribution_key< D >( f_dist , s0 ) } );
  sh.h0 = sh.root->h;
  sh.root_owner = f_dist.owner( s0 , p );

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );

  // the weight 1 known at compile time, so that A* pays nothing for it
  if( f_w.unit() )
   run_workers( f_ex , p , [ & ]( unsigned i ) {
    w[ i ]->run( UnitWeight< Cost >() );
    } );
  else
   run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run( f_w ); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.sent += x->f_st.sent;
   st.messages += sh.mail.batches( x->f_id );
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  sol.improvements = std::move( sh.improvements );
  finish_solution( sol , sh.goal , sh.limit_hit.load() , sh.h0 ,
                   f_w , f_domain );
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Mailbox = ADS::Mailbox< Msg ,
                               typename ChannelsOf::template type< Msg > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( unsigned np , std::size_t batch , const Stopwatch & c )
   : p( np ) , mail( np , batch ) , clock( c ) , best( new Best[ np ] ) {}

  // the best f of a worker, which it publishes when it lends nodes
  struct alignas( cache_line ) Best {
   std::atomic< Cost > f{ infinite_cost< Cost >() };
   };

  unsigned p;                                   // the number of workers
  Mailbox mail;                                 // the messages
  std::optional< Msg > root;                    // the root
  unsigned root_owner = 0;                      // its owner
  alignas( cache_line ) std::atomic< Cost > incumbent{
                                                   infinite_cost< Cost >() };
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  NodeT * goal = nullptr;                       // the best goal found
  // the incumbents of an anytime search, under goal_mutex
  std::vector< typename SolutionT::Improvement > improvements;
  const Stopwatch & clock;                      // started at the start
  Cost h0{};                                    // the h of the root
  std::unique_ptr< Best[] > best;               // the best f of each one
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its OPEN and CLOSED lists, its outboxes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const HDAStar & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_dist( a.f_dist ) {}

  // the main loop, with the rules of w [see Weight and UnitWeight]
  template< class W >
  void run( const W & w ) {
   if( f_id == f_sh.root_owner )
    insert( w , *f_sh.root , f_a.f_domain.hash( f_sh.root->s ) );
   unsigned since_flush = 0;
   ADS::Backoff backoff;
   while( ! f_sh.mail.stopped() ) {
    //! [hash-distributed-batch]
    v_in.clear();
    f_sh.mail.poll( f_id , [ & ]( Msg && m ) {
     v_in.push_back( std::move( m ) );
     } );
    if( const std::size_t n = v_in.size() ) {
     // the hashes and the prefetches of the whole batch first, so that
     // the cache misses of the lookups overlap
     v_hash.resize( n );
     for( std::size_t k = 0 ; k < n ; ++k ) {
      v_hash[ k ] = f_a.f_domain.hash( v_in[ k ].s );
      f_table.prefetch( v_hash[ k ] );
      }
     for( std::size_t k = 0 ; k < n ; ++k )
      if( v_in[ k ].lent )
       expand_lent( w , v_in[ k ] );
      else
       insert( w , v_in[ k ] , v_hash[ k ] );
     }
    //! [hash-distributed-batch]

    if( NodeT * n = next( w ) ) {
     backoff.reset();
     if( ( f_a.f_options.outsource <= 0 ) || ! lend( n ) )
      expand( w , n );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      f_sh.mail.flush( f_id );
      }
     if( f_a.f_options.node_limit &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      f_sh.limit_hit.store( true );
      f_sh.mail.stop();
      }
     continue;
     }

    // nothing to expand: send everything, and end if no worker has work
    // and no message is in transit
    if( f_sh.mail.quiescent( f_id ) )
     return;
    backoff.pause();
    }
   }

  // the best node of OPEN with f' < incumbent (any node, if anytime) not
  // pruned, removed; nullptr if none; its f, or the infinity, published if
  // the nodes are lent
  template< class W >
  NodeT * next( const W & w ) {
   NodeT * n = nullptr;
   while( ! f_open.empty() ) {
    const OpenEntry< NodeT > e = f_open.top();
    if( e.stale() ) {
     f_open.pop();
     continue;
     }
    const Cost z = f_sh.incumbent.load( std::memory_order_relaxed );
    if( ! w.live( w.key( e.g , e.node->h ) , z ) )   // the key of e
     break;
    f_open.pop();
    if( ! w.unit() && w.prunes( e.g , e.node->h , z ) )
     continue;                  // with the weight 1, live is f < z already
    n = e.node;
    break;
    }
   if( f_a.f_options.outsource > 0 )
    f_sh.best[ f_id ].f.store( n ? n->g + n->h : infinite_cost< Cost >() ,
                               std::memory_order_relaxed );
   return( n );
   }

  // lends n, which is not a goal, to the worker whose published f is the
  // least among those not smaller than f( n ) plus the gap; false if no
  // worker is that far behind
  bool lend( NodeT * n ) {
   const double f = double( n->g + n->h ) + f_a.f_options.outsource;
   unsigned to = f_id;
   Cost best = infinite_cost< Cost >();
   for( unsigned j = 0 ; j < f_sh.p ; ++j ) {
    if( j == f_id )
     continue;
    const Cost fj = f_sh.best[ j ].f.load( std::memory_order_relaxed );
    if( ( double( fj ) >= f ) && ( ( to == f_id ) || ( fj < best ) ) ) {
     best = fj;
     to = j;
     }
    }
   if( ( to == f_id ) || f_a.f_domain.is_goal( n->state ) )
    return( false );
   n->closed = true;
   n->expanded = true;
   f_sh.mail.send( f_id , to , Msg{ n->state , n->g , n->h , n , n->op ,
                                    n->key , true } );
   f_sh.mail.flush( f_id );
   ++f_st.sent;
   return( true );
   }

  // expands the node lent by its owner, which keeps it in its CLOSED list
  template< class W >
  void expand_lent( const W & w , const Msg & m ) {
   if( w.prunes( m.g , m.h , f_sh.incumbent.load(
                                         std::memory_order_relaxed ) ) )
    return;
   ++f_st.expanded;
   generate( w , m.s , m.g , m.h , m.o.op , m.key , m.parent );
   }

  // expands n, or records it as the incumbent if it is a goal
  template< class W >
  void expand( const W & w , NodeT * n ) {
   n->closed = true;
   ++f_st.expanded;
   if( n->expanded )
    ++f_st.reexpanded;
   n->expanded = true;

   if( f_a.f_domain.is_goal( n->state ) ) {
    std::lock_guard< std::mutex > lock( f_sh.goal_mutex );
    if( n->g < f_sh.incumbent.load() ) {
     f_sh.incumbent.store( n->g );
     f_sh.goal = n;
     if( f_a.f_w.anytime() )
      f_sh.improvements.push_back( { n->g , f_sh.h0 , f_st.expanded ,
                                     f_sh.clock.seconds() } );
     }
    return;
    }

   //! [hash-distributed-expand]
   generate( w , n->state , n->g , n->h , n->op , n->key , n );
   }

  // sends the successors of the state s of the node n, of cost g and
  // heuristic h, generated by op, of key key, to their owners
  template< class W >
  void generate( const W & w , const State & s , Cost gs , Cost h ,
                 std::uint8_t op , std::uint32_t key , NodeT * n ) {
   const auto owner = successor_owners< D >( f_dist , s , key , f_sh.p );
   for_each_child( f_a.f_domain , s , h , op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t ot ) {
    ++f_st.generated;
    const Cost g = gs + c;
    if( w.prunes( g , ht , f_sh.incumbent.load(
                                    std::memory_order_relaxed ) ) )
     return;
    const auto [ dst , key ] = owner( t );
    if( dst == f_id ) {
     insert( w , Msg{ t , g , ht , n , ot , key } ,
             f_a.f_domain.hash( t ) );
     return;
     }
    f_sh.mail.send( f_id , dst , Msg{ t , g , ht , n , ot , key } );
    ++f_st.sent;
    } );
   //! [hash-distributed-expand]
   }

  // adds the node of m, whose state has hash h, to OPEN, unless it is a
  // duplicate
  //! [hash-distributed-insert]
  template< class W >
  void insert( const W & w , const Msg & m , std::uint64_t h ) {
   if( w.prunes( m.g , m.h , f_sh.incumbent.load(
                                         std::memory_order_relaxed ) ) )
    return;
   auto [ x , fresh ] = f_table.find_or_insert( m.s , h , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ m.s , m.g , m.h , m.parent , false , false , m.o.op ,
                m.key };
    return( y );
    } );
   if( fresh ) {
    f_open.push( x , m.g , w.key( m.g , m.h ) );
    return;
    }
   if( w.improves( x->g , m.g , x->closed ,
                   [ & ] { return( m.parent->h ); } , m.h ) ) {
    x->g = m.g;
    x->parent = m.parent;
    x->op = m.o.op;
    x->closed = false;
    f_open.push( x , m.g , w.key( m.g , m.h ) );
    }
   else
    ++f_st.duplicates;
   }
  //! [hash-distributed-insert]

  const HDAStar & f_a;                          // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  Dist f_dist;                                  // its distribution
  NodePool< NodeT > f_pool;                     // its nodes
  NodeTable< NodeT > f_table;                   // its CLOSED list
  Open f_open;                                  // its OPEN list
  std::vector< Msg > v_in;                      // the last inbox received
  std::vector< std::uint64_t > v_hash;          // the hashes of v_in
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                    ///< the work distribution
 HDAOptions f_options;                                 ///< the options
 Weight< Cost > f_w;                                   ///< the weight
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( HDAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_HDAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File HDAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
