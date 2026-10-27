/*--------------------------------------------------------------------------*/
/*--------------------- File ParallelBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the parallel best-first searches with an incumbent, i.e.,
 * the parallel Branch-and-Bound (B&B) on a runtime (the C++ threads,
 * OpenMP or FastFlow): SharedBranchAndBound, whose workers extract the
 * nodes from one queue as the workers of SPA* do, and
 * DistributedBranchAndBound, whose workers have a queue each and send the
 * children to their owners as the workers of HDA* do. The problem, the
 * options and the Solution are those of BranchAndBound.h, and the design
 * is discussed in docs/branch-and-bound.md.
 *
 * The two searches are the following, with the notation of the page
 * "Parallel branch-and-bound" of the documentation and the procedures
 * PROCESS and PRUNABLE of the sequential search [see BranchAndBound.h]; the
 * incumbent \f$ z \f$ is shared (BBIncumbent), DIVE(K) is the child of
 * \f$ K \f$ that the worker goes on with (the last one with depth-first, the
 * one of least bound with plunging if it is close enough to the least open
 * bound, none otherwise), and \f$ p \f$ is the number of workers:
 *
 * \code{.unparsed}
 *  1  procedure SHARED(R)       ▷ one queue Q, under the lock L, for all
 *  2     Q ← R;  busy ← 0;  done ← false
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        d ← nil
 *  5        loop
 *  6           if d ≠ nil then n ← d;  d ← nil
 *  7           else n ← NEXT();  if n = nil then return
 *  8           if a limit is reached then
 *  9              n is a leaf (open);  busy ← busy − 1;  done ← true;  return
 * 10           K ← PROCESS(n)
 * 11           lock L;  d ← DIVE(K);  push the nodes of K but d in Q
 * 12           if d = nil then busy ← busy − 1
 * 13           unlock L
 * 14  procedure NEXT()
 * 15     loop
 * 16        if done then return nil
 * 17        if Q = ∅ then                             ▷ read without the lock
 * 18           if busy = 0 then
 * 19              lock L;  if Q = ∅ and busy = 0 then done ← true;  unlock L
 * 20           else wait a little
 * 21           continue
 * 22        lock L
 * 23        if PRUNABLE(the least bound of Q) then
 * 24           the nodes of Q are leaves (pruned);  Q ← ∅;  unlock L
 * 25           continue
 * 26        busy ← busy + 1;  n ← the first node of Q, removed;  unlock L
 * 27        return n
 * 28  procedure DISTRIBUTED(R)  ▷ a queue Q_i per worker, nodes to the owners
 * 29     send each node of R to its owner
 * 30     in parallel, on each worker i = 0, ..., p − 1 do
 * 31        d ← nil
 * 32        while the search is not over do
 * 33           push the nodes received in Q_i
 * 34           if d ≠ nil then n ← d;  d ← nil
 * 35           else if Q_i ≠ ∅ and not PRUNABLE(the least bound of Q_i) then
 * 36              n ← the first node of Q_i, removed
 * 37           else
 * 38              the nodes of Q_i are leaves (pruned);  Q_i ← ∅
 * 39              send all the outboxes;  WAIT(i);  continue
 * 40           if a limit is reached then n is a leaf (open); stop all; return
 * 41           K ← PROCESS(n);  d ← DIVE(K)
 * 42           for each m ∈ K but d do
 * 43              j ← hash(m) mod p, or a random worker without hash(m)
 * 44              if j = i then push m in Q_i else send m to j
 * 45           every k evaluations, send the outboxes not yet full
 * \endcode
 *
 * The counter busy of the shared search is raised under the lock when a node
 * is extracted (line 26) and lowered only after the children of the node
 * have been pushed (line 12), so that \f$ Q = \emptyset \f$ and busy
 * \f$ = 0 \f$, read under the lock at line 19, mean that no node is open
 * anywhere. The distributed search detects its end by WAIT, i.e., by the
 * three counters of ThreadMailbox [see ThreadMailbox.h]; since a tree has no
 * duplicates, the owner of line 43 only balances the load.
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

#ifndef __PBeFS_ParallelBranchAndBound
 #define __PBeFS_ParallelBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <memory>

#include <mutex>

#include <random>

#include <type_traits>

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/core/ThreadMailbox.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ParallelBBOptions -------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the parallel B&B with values of type V
/** Those of the sequential search, and the parameters of the parallel ones
 * [see ParallelParams], with batches of 16 nodes flushed at every
 * evaluation, since an evaluation is usually "expensive"; node_limit
 * bounds the evaluations of all the workers together, batch is used by the
 * searches that send nodes (the distributed, the farm and the MPI ones),
 * flush_period by the distributed and the MPI ones, and seed gives the
 * owners of the nodes when the problem gives no hash( n ), while round is
 * used by the deterministic search only and split by the SelfSplit one. */

template< class V = double >
struct BasicParallelBBOptions : public BasicBBOptions< V > ,
                                public ParallelParams {

 BasicParallelBBOptions( void ) {
  batch = 16;
  flush_period = 1;
  }

 /// the nodes evaluated at each round of the deterministic search, 0 for
 /// one per worker [see DeterministicBranchAndBound]
 std::size_t round = 0;

 /// the open nodes per worker at which the sampling of the SelfSplit search
 /// stops, 0 for 4 [see SelfSplitBranchAndBound]
 std::size_t split = 0;

 };  // end( struct( BasicParallelBBOptions ) )

/// the options of the parallel B&B with values of type double
using ParallelBBOptions = BasicParallelBBOptions<>;

namespace detail {

/// true if std::atomic< V > exists and is always lock-free
template< class V , bool = std::is_trivially_copyable_v< V > >
struct bb_lock_free : std::false_type {};

template< class V >
struct bb_lock_free< V , true >
 : std::bool_constant< std::atomic< V >::is_always_lock_free > {};

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BBIncumbent ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the incumbent shared by the workers of a parallel B&B
/** The value (in the minimization) is atomic when V has a lock-free
 * std::atomic, so that the workers read it without a lock when they prune;
 * an improvement takes the lock, checks the value again, and stores the
 * node and calls the callback under it. For another V (e.g., the rationals
 * of an exact arithmetic) the value is read under the lock too. */

template< class Node , class V = double >
class BBIncumbent {

 public:

 /// constructor, taking the initial value (infinity() of BBValueTraits if
 /// none)

 explicit BBIncumbent( V v ) : f_value( v ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the current value

 [[nodiscard]] V value( void ) const {
  if constexpr( c_atomic )
   return( f_value.load( std::memory_order_relaxed ) );
  else {
   std::lock_guard< std::mutex > guard( f_m );
   return( f_value );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// offers the value v of the node n of the problem P; true if it is a new
 /// incumbent

 template< class P >
 bool improve( V v , const Node & n , const BasicBBOptions< V > & o ) {
  if( c_atomic && ( v >= value() ) )
   return( false );
  std::lock_guard< std::mutex > guard( f_m );
  if( v >= stored() )
   return( false );
  if constexpr( c_atomic )
   f_value.store( v , std::memory_order_relaxed );
  else
   f_value = v;
  f_node = n;
  if( o.on_incumbent )
   o.on_incumbent( detail::bb_signed< P >( v ) );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// moves out the node of the incumbent, if any

 std::optional< Node > take_node( void ) { return( std::move( f_node ) ); }

 private:

 /// true if the value is a lock-free std::atomic
 static constexpr bool c_atomic = detail::bb_lock_free< V >::value;

 // the value, read under the lock
 [[nodiscard]] V stored( void ) const {
  if constexpr( c_atomic )
   return( f_value.load( std::memory_order_relaxed ) );
  else
   return( f_value );
  }

 alignas( cache_line ) std::conditional_t< c_atomic , std::atomic< V > , V >
  f_value;                                      ///< the value
 mutable std::mutex f_m;                        ///< protects the rest
 std::optional< Node > f_node;                  ///< the node

 };  // end( class( BBIncumbent ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// fills the Solution of a parallel B&B from the incumbent and the leaves
/// of all the workers, merged in fr

template< BranchAndBound P >
void bb_finish( BBSolution< typename P::node_type , bb_value_t< P > > & sol ,
                BBIncumbent< typename P::node_type , bb_value_t< P > > & inc ,
                BBFrontier< typename P::node_type , bb_value_t< P > > & fr ,
                bool limit , bool time )
{
 using V = bb_value_t< P >;
 using Tr = BBValueTraits< V >;
 const V z = inc.value();
 if( limit )
  sol.status = time ? BBStatus::time_limit : BBStatus::node_limit;
 else
  sol.status = Tr::finite( z ) ? BBStatus::optimal : BBStatus::infeasible;
 sol.node = inc.take_node();
 if( Tr::finite( z ) )
  sol.value = bb_signed< P >( z );
 const V lb = std::min( z , fr.low );
 sol.bound = bb_signed< P >( lb );
 sol.gap = Tr::finite( z ) && Tr::finite( lb ) ? V( z - lb )
                                               : Tr::infinity();
 fr.template take< P >( sol.frontier );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the index of the child that the worker goes on with, or kids.size() if
/// none: with depth-first the last child, which the queue would extract
/// next anyway, and with plunging the child of least bound, if it is close
/// enough to the best bound

template< class Item , class V >
std::size_t bb_dive( const std::vector< Item > & kids , V best , V inc ,
                     const BasicBBOptions< V > & o )
{
 if( ( o.selection == Selection::depth_first ) && ! kids.empty() )
  return( kids.size() - 1 );
 if( ( o.selection != Selection::plunging ) || kids.empty() )
  return( kids.size() );
 std::size_t c = 0;
 for( std::size_t j = 1 ; j < kids.size() ; ++j )
  if( kids[ j ].bound <= kids[ c ].bound )
   c = j;
 using Tr = BBValueTraits< V >;
 const V b0 = std::min( best , kids[ c ].bound );
 if( ! Tr::finite( inc ) ||
     Tr::within( kids[ c ].bound , b0 , inc , o.plunge_gap ) )
  return( c );
 return( kids.size() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// true if the node or the time limit of o is reached, given the
/// evaluations of all the workers and the clock of the search; then limit
/// is set, and time too if the time limit is reached

template< class V >
bool bb_limit( const BasicBBOptions< V > & o , std::uint64_t evaluated ,
               const Stopwatch & clock , std::atomic< bool > & limit ,
               std::atomic< bool > & time )
{
 const bool nodes = o.node_limit && ( evaluated >= o.node_limit );
 const bool late = ( o.time_limit > 0 ) && ( clock.seconds() >= o.time_limit );
 if( ! ( nodes || late ) )
  return( false );
 if( late )
  time.store( true );
 limit.store( true );
 return( true );
 }

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SharedBranchAndBound -----------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B with one queue shared by the workers
/** The workers extract the nodes from one queue under a lock, evaluate and
 * branch them outside of it, and push the children back under it (with
 * depth-first or plunging, the child that dives is kept by the worker). The
 * search ends when the queue is empty and no worker holds a node; the open
 * nodes are dropped at once when none of them can improve the incumbent. The
 * order of extraction is close to the sequential one, and the lock is cheap
 * when the evaluation dominates, which is the usual case of a B&B. */

template< BranchAndBound P , class Runtime = ThreadRuntime >
class SharedBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search)

 explicit SharedBranchAndBound( const P & problem ,
                                const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution [see
 /// BranchAndBoundSearch::solve()]

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const unsigned p = f_options.threads;
  const Selection sel = f_options.selection;
  Stopwatch clock;
  const Value start = f_options.incumbent ?
                      detail::bb_signed< P >( *f_options.incumbent ) : inf;
  Shared sh( start , sel );
  for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) ) {
   const Value key = detail::bb_key( f_p , l.node , l.bound , l.depth ,
                                      sel );
   sh.q.push( std::move( l.node ) , l.bound , l.depth , key );
   }
  sh.open.store( sh.q.size() );
  sh.max_open = sh.q.size();

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , clock );
  Runtime::run( p , [ & ]( unsigned i ) { w[ i ]->run( i ); } );

  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( f_options.keep_frontier );
  fr.add( sh.q , BBLeafKind::open );          // left open by a limit
  for( auto & x : w ) {
   sol.stats.add( x->f_st );
   fr.merge( std::move( x->f_fr ) );
   }
  sol.stats.max_open = sh.max_open;
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 using Lock = typename Runtime::lock_type;
 using Item = typename BBQueue< Node , Value >::Item;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( Value z , Selection s ) : inc( z ) , q( s ) {}

  BBIncumbent< Node , Value > inc;             // the incumbent
  alignas( cache_line ) Lock lock;              // protects q and max_open
  BBQueue< Node , Value > q;                   // the open nodes
  std::size_t max_open = 0;                    // the most open nodes
  alignas( cache_line ) std::atomic< std::size_t > open{ 0 };  // q.size()
  alignas( cache_line ) std::atomic< unsigned > busy{ 0 };  // holding a node
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  alignas( cache_line ) std::atomic< bool > done{ false };
  std::atomic< bool > limit{ false };          // stopped by a limit
  std::atomic< bool > time{ false };           // ... the time limit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker

 struct alignas( cache_line ) Worker {

  Worker( const SharedBranchAndBound & a , Shared & s , const Stopwatch & c )
   : f_a( a ) , f_sh( s ) , f_clock( c ) ,
     f_fr( a.f_options.keep_frontier ) {}

  // the main loop
  void run( unsigned i ) {
   auto w = detail::bb_worker( f_a.f_p , i );
   const auto & o = f_a.f_options;
   const Selection sel = o.selection;
   std::optional< Item > dive;               // the node of the dive
   std::vector< Item > kids;                 // the children of a node
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };

   for( ;; ) {
    std::optional< Item > it;
    if( dive ) {
     it.emplace( std::move( *dive ) );
     dive.reset();
     }
    else
     if( ! ( it = next() ) )
      return;

    if( stop() ) {                     // the node stays open
     f_fr.add( std::move( it->node ) , it->bound , it->depth ,
               BBLeafKind::open );
     f_sh.busy.fetch_sub( 1 );
     return;
     }

    kids.clear();
    detail::bb_process< P >( it->node , it->bound , it->depth , w , inc ,
                             improve , [ & ]( Node && c , Value b ,
                                              unsigned d ) {
                              kids.push_back( Item{ std::move( c ) , b , d } );
                              } , o , f_st , f_fr );
    f_sh.evaluated.fetch_add( 1 );

    std::lock_guard< Lock > guard( f_sh.lock );
    const std::size_t k = detail::bb_dive( kids , f_sh.q.best_bound() ,
                                           inc() , o );
    for( std::size_t j = 0 ; j < kids.size() ; ++j )
     if( j == k )
      dive.emplace( std::move( kids[ j ] ) );
     else {
      const Value key = detail::bb_key( f_a.f_p , kids[ j ].node ,
                                         kids[ j ].bound , kids[ j ].depth ,
                                         sel );
      f_sh.q.push( std::move( kids[ j ].node ) , kids[ j ].bound ,
                   kids[ j ].depth , key );
      }
    f_sh.open.store( f_sh.q.size() );
    f_sh.max_open = std::max( f_sh.max_open , f_sh.q.size() );
    if( ! dive )
     f_sh.busy.fetch_sub( 1 );
    }
   }

  // the next node of the queue, counting the worker as busy; none when the
  // search is over
  std::optional< Item > next( void ) {
   Backoff backoff;
   for( ;; ) {
    if( f_sh.done.load( std::memory_order_acquire ) )
     return( std::nullopt );
    if( f_sh.open.load() == 0 ) {      // wait without the lock
     if( f_sh.busy.load() == 0 ) {
      std::lock_guard< Lock > guard( f_sh.lock );
      if( f_sh.q.empty() && ( f_sh.busy.load() == 0 ) )
       f_sh.done.store( true , std::memory_order_release );
      continue;
      }
     backoff.pause();
     continue;
     }
    std::lock_guard< Lock > guard( f_sh.lock );
    if( f_sh.q.empty() )
     continue;
    const Value best = f_sh.q.best_bound();
    if( detail::bb_prunable( best , f_sh.inc.value() , f_a.f_options ) ) {
     f_st.pruned_extracted += f_fr.add( f_sh.q , BBLeafKind::pruned );
     f_sh.open.store( 0 );                       // nothing can improve
     continue;
     }
    f_sh.busy.fetch_add( 1 );
    Item it = f_sh.q.pop();
    f_sh.open.store( f_sh.q.size() );
    return( it );
    }
   }

  // true if a limit is reached, which ends the search
  bool stop( void ) {
   if( ! detail::bb_limit( f_a.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.done.store( true , std::memory_order_release );
   return( true );
   }

  const SharedBranchAndBound & f_a;             // the search
  Shared & f_sh;                                // the shared data
  const Stopwatch & f_clock;                    // the clock of the search
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options

 };  // end( class( SharedBranchAndBound ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS DistributedBranchAndBound ---------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B with a queue per worker and the children sent to
/// their owners
/** Each worker extracts the nodes from its own queue, and sends each child,
 * in batches, to the worker that owns it:
 * \f$ \mathrm{hash}( n ) \bmod p \f$ if the problem gives hash( n ), a
 * random worker otherwise (with depth-first or plunging,
 * the child that dives is kept). Since a tree has no duplicates, the owners
 * only balance the load. The incumbent is shared, and the termination is
 * detected as in HDAStar, by the three counters of ThreadMailbox [see
 * ThreadMailbox.h]. */

template< BranchAndBound P , class Runtime = ThreadRuntime >
class DistributedBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

 /// a node travelling to its owner
 struct Msg {
  Node node;
  Value bound;
  unsigned depth;
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search)

 explicit DistributedBranchAndBound( const P & problem ,
                                     const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution [see
 /// BranchAndBoundSearch::solve()]; each node goes to its owner, which is
 /// the worker \f$ i \bmod p \f$ for the node \f$ i \f$ if the problem
 /// gives no hash( n )

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const unsigned p = f_options.threads;
  Stopwatch clock;
  const Value start = f_options.incumbent ?
                      detail::bb_signed< P >( *f_options.incumbent ) : inf;
  Shared sh( p , start );
  {
   std::vector< std::vector< Msg > > b( p );
   std::size_t i = 0;
   for( auto & l : detail::bb_starts( f_p , std::move( seeds ) ) ) {
    unsigned dst;
    if constexpr( HashedBranchAndBound< P > )
     dst = unsigned( std::uint64_t( f_p.hash( l.node ) ) % p );
    else
     dst = unsigned( i % p );
    ++i;
    b[ dst ].push_back( Msg{ std::move( l.node ) , l.bound , l.depth } );
    }
   for( unsigned d = 0 ; d < p ; ++d )
    sh.mail.seed( d , b[ d ] );
   }

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i , clock );
  Runtime::run( p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  // the nodes left open by a limit bound the optimum too
  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( f_options.keep_frontier );
  std::vector< Msg > in;
  const auto open = [ & ]( Msg & m ) {
   fr.add( std::move( m.node ) , m.bound , m.depth , BBLeafKind::open );
   };
  for( unsigned i = 0 ; i < p ; ++i ) {
   Worker & x = *w[ i ];
   x.f_st.messages = x.f_out.messages();
   sol.stats.add( x.f_st );
   fr.merge( std::move( x.f_fr ) );
   fr.add( x.f_q , BBLeafKind::open );
   if( x.f_dive )
    fr.add( std::move( x.f_dive->node ) , x.f_dive->bound ,
            x.f_dive->depth , BBLeafKind::open );
   for( auto & o : x.f_out.boxes() )
    for( Msg & m : o )
     open( m );
   sh.mail.receive( i , in , true );
   for( Msg & m : in )
    open( m );
   }
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 using Mailbox = ThreadMailbox< Msg , Runtime >;
 using Item = typename BBQueue< Node , Value >::Item;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( unsigned np , Value z ) : p( np ) , mail( np ) , inc( z ) {}

  unsigned p;                                   // the number of workers
  Mailbox mail;                                 // the messages
  BBIncumbent< Node , Value > inc;              // the incumbent
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its queue, its outboxes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const DistributedBranchAndBound & a , Shared & s , unsigned i ,
          const Stopwatch & c )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_clock( c ) ,
     f_q( a.f_options.selection ) ,
     f_gen( a.f_options.seed + i ) ,
     f_out( s.mail , i , a.f_options.batch ) ,
     f_fr( a.f_options.keep_frontier ) {}

  // the main loop
  void run( void ) {
   auto w = detail::bb_worker( f_a.f_p , f_id );
   const auto & o = f_a.f_options;
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };
   unsigned since_flush = 0;

   while( ! f_sh.mail.done() ) {
    if( std::size_t n = f_sh.mail.receive( f_id , v_in ,
                                           f_q.empty() && ! f_dive ) ) {
     for( Msg & m : v_in )
      push( Item{ std::move( m.node ) , m.bound , m.depth } );
     f_sh.mail.processed( n );
     }

    if( std::optional< Item > it = next() ) {
     if( stop() ) {                   // the node stays open
      f_fr.add( std::move( it->node ) , it->bound , it->depth ,
                BBLeafKind::open );
      return;
      }
     v_kids.clear();
     detail::bb_process< P >( it->node , it->bound , it->depth , w , inc ,
                              improve , [ & ]( Node && c , Value b ,
                                               unsigned d ) {
                               v_kids.push_back( Item{ std::move( c ) , b ,
                                                       d } );
                               } , o , f_st , f_fr );
     f_sh.evaluated.fetch_add( 1 );
     const std::size_t k = detail::bb_dive( v_kids , f_q.best_bound() ,
                                            inc() , o );
     for( std::size_t j = 0 ; j < v_kids.size() ; ++j )
      if( j == k )
       f_dive.emplace( std::move( v_kids[ j ] ) );
      else
       send( std::move( v_kids[ j ] ) );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      f_out.flush( false );
      }
     continue;
     }

    f_out.flush( true );     // nothing to evaluate: send everything, then
    f_sh.mail.wait( f_id );  // wait for mail or for the end
    }
   }

  // pushes a node in the queue
  void push( Item && it ) {
   const Value key = detail::bb_key( f_a.f_p , it.node , it.bound ,
                                      it.depth , f_a.f_options.selection );
   f_q.push( std::move( it.node ) , it.bound , it.depth , key );
   f_st.max_open = std::max( f_st.max_open , f_q.size() );
   }

  // the next node: the dive, or the best of the queue; none if there is
  // none, or if none can improve the incumbent (which empties the queue)
  std::optional< Item > next( void ) {
   if( f_dive ) {
    std::optional< Item > it( std::move( f_dive ) );
    f_dive.reset();
    return( it );
    }
   if( f_q.empty() )
    return( std::nullopt );
   const Value best = f_q.best_bound();
   if( detail::bb_prunable( best , f_sh.inc.value() , f_a.f_options ) ) {
    f_st.pruned_extracted += f_fr.add( f_q , BBLeafKind::pruned );
    return( std::nullopt );
    }
   return( f_q.pop() );
   }

  // sends a child to its owner, or keeps it if the worker owns it
  void send( Item && it ) {
   const unsigned p = f_sh.p;
   unsigned dst;
   if constexpr( HashedBranchAndBound< P > )
    dst = unsigned( std::uint64_t( f_a.f_p.hash( it.node ) ) % p );
   else
    dst = unsigned( random_below( f_gen , p ) );
   if( dst == f_id ) {
    push( std::move( it ) );
    return;
    }
   f_out.push( dst , Msg{ std::move( it.node ) , it.bound , it.depth } );
   ++f_st.sent;
   }

  // true if a limit is reached, which ends the search
  bool stop( void ) {
   if( ! detail::bb_limit( f_a.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.mail.stop();
   return( true );
   }

  const DistributedBranchAndBound & f_a;        // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  const Stopwatch & f_clock;                    // the clock of the search
  BBQueue< Node , Value > f_q;                  // its open nodes
  std::optional< Item > f_dive;                 // the node of its dive
  std::mt19937_64 f_gen;                        // the random owners
  typename Mailbox::Outboxes f_out;             // its outboxes
  std::vector< Msg > v_in;                      // the last inbox received
  std::vector< Item > v_kids;                   // the children of a node
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options

 };  // end( class( DistributedBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ParallelBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------- End File ParallelBranchAndBound.h --------------------*/
/*--------------------------------------------------------------------------*/
