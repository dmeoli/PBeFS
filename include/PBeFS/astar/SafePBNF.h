/*--------------------------------------------------------------------------*/
/*---------------------------- File SafePBNF.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SafePBNF, the livelock-free Parallel
 * Best-NBlock-First search of Burns, Lemons, Ruml and Zhou (JAIR 2010).
 *
 * The search follows the pseudocode of appendix A of Burns et al. (JAIR
 * 2010), with the notation of the page "Best-first search and A*" of the
 * documentation. Each nblock \f$ b \f$ has its lists
 * \f$ \mathrm{OPEN}( b ) \f$ and \f$ \mathrm{CLOSED}( b ) \f$, its
 * interference scope
 * \f$ I( b ) \f$ (the nblocks that, if in use, keep \f$ b \f$ from being
 * free), the counters \f$ \sigma( b ) \f$ of the nblocks in use and
 * \f$ \sigma_h( b ) \f$ of the hot nblocks whose interference scopes hold
 * \f$ b \f$, and the value \f$ f( b ) \f$, the least \f$ f \f$ of its OPEN
 * list, which its worker publishes after each expansion; one lock protects
 * the nblock graph, i.e., the counters, the hot flags and the free list,
 * which is ordered by \f$ f( b ) \f$, and \f$ m \f$ is the minimum number of
 * expansions in an nblock before a switch:
 *
 * \code{.unparsed}
 *  1  procedure THREADSEARCH()                        ▷ on each worker
 *  2     b ← nil
 *  3     loop
 *  4        b ← NEXTNBLOCK(b);  if b = nil then return
 *  5        k ← 0
 *  6        while not SHOULDSWITCH(b, k) do
 *  7           u ← the first node of OPEN(b), removed (stale entries dropped)
 *  8           if f(u) ≥ z then OPEN(b) ← ∅;  break
 *  9           close u
 * 10           if u is a goal then z ← min(z, g(u)), recording u if g(u) < z
 * 11           else for each successor v of u with g(u) + c(u,v) + h(v) < z do
 * 12              insert v in the lists of its nblock, as A* does
 * 13           k ← k + 1;  publish f(b)
 * 14  procedure SHOULDSWITCH(b, k)
 * 15     if OPEN(b) = ∅ then return true
 * 16     if k < m then return false
 * 17     k ← 0;  lock the graph
 * 18     b' ← the nblock of I(b) of least published f
 * 19     if f(free list) < f(b) or f(b') < f(b) then
 * 20        if f(b') < f(free list) then SETHOT(b')
 * 21        return true
 * 22     SETCOLD each hot nblock of I(b);  return false
 * 23  procedure NEXTNBLOCK(b)
 * 24     if b ≠ nil, OPEN(b) ≠ ∅ and b was not just set hot then
 * 25        if the graph is locked by another worker then return b
 * 26     lock the graph
 * 27     if b ≠ nil then
 * 28        if f(b) < z, f(b) ≤ f(b'), f(b) ≤ f(free list) and b was not
 *              just set hot then return b
 * 29        RELEASE(b)
 * 30     drop from the free list the nblocks with no node of f < z
 * 31     while the free list is empty do
 * 32        if no nblock is in use then done ← true;  wake all;  return nil
 * 33        sleep until woken
 * 34     b ← the first nblock of the free list;  ACQUIRE(b);  return b
 * 35  procedure ACQUIRE(b)
 * 36     b leaves the free list and is in use
 * 37     for each n ∈ I(b) do σ(n) ← σ(n) + 1;  n leaves the free list
 * 38  procedure RELEASE(b)
 * 39     b is not in use
 * 40     for each n ∈ I(b) do
 * 41        σ(n) ← σ(n) − 1
 * 42        if σ(n) = 0 then
 * 43           if σ_h(n) = 0 and n is hot then SETCOLD(n)
 * 44           FREE(n)
 * 45     FREE(b)
 * 46  procedure SETHOT(b)
 * 47     if b is hot, σ(b) = 0, or a hot nblock of I(b) has a smaller f
 *           then return
 * 48     b is hot
 * 49     for each n ∈ I(b) do
 * 50        if n is hot then SETCOLD(n)
 * 51        n leaves the free list;  σ_h(n) ← σ_h(n) + 1
 * 52  procedure SETCOLD(b)
 * 53     b is not hot
 * 54     for each n ∈ I(b) do
 * 55        σ_h(n) ← σ_h(n) − 1
 * 56        if σ(n) = 0 and σ_h(n) = 0 then
 * 57           if n is hot then SETCOLD(n)
 * 58           FREE(n)
 * 59  procedure FREE(n)
 * 60     if n is not in use, σ(n) = σ_h(n) = 0 and f(n) < z then
 * 61        put n in the free list, with the key f(n);  wake all
 * \endcode
 *
 * The worker of an nblock expands its nodes and inserts their children
 * (line 12) without any lock, since the nblocks of the children lie in the
 * duplicate detection scope of \f$ b \f$, which no other worker can hold
 * while \f$ b \f$ is in use. The hot nblocks (lines 20 and 46 to 51) prevent
 * the livelock of plain PBNF: an nblock that interferes with a hot one leaves
 * the free list and is released at the next check of its worker, so that
 * every hot nblock eventually becomes free. The search ends at line 32, when
 * no nblock is in use and none is free, the nblocks with no node below the
 * incumbent having left the free list at line 30. The value
 * \f$ f( b' ) \f$ read at lines 18 and 28 is the one published by the other
 * workers, which may be slightly out of date; this changes the order of
 * expansion, and not the cost of the path found.
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

#ifndef __PBeFS_SafePBNF
 #define __PBeFS_SafePBNF
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <condition_variable>

#include <memory>

#include <mutex>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/astar/HDAStar.h"

#include "ADS/tree/BTreeMap.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS PBNFOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of SafePBNF

struct PBNFOptions : public ParallelOptions {

 /// expansions in an nblock before a worker may switch to a better one
 unsigned min_expansions = 32;

 };  // end( struct( PBNFOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS SafePBNF -----------------------------*/
/*--------------------------------------------------------------------------*/
/// Safe Parallel Best-NBlock-First search
/** The state space is partitioned into nblocks by the abstraction of the
 * Domain, each nblock with its own OPEN and CLOSED lists. The duplicate
 * detection scope of an nblock is the nblock and its abstract successors,
 * the only nblocks its children can fall in; a worker acquires an nblock
 * whose scope does not overlap the scope of any nblock in use, and then
 * expands its nodes and checks their children for duplicates with no lock
 * at all [see Zhou and Hansen, AAAI 2007, for parallel structured duplicate
 * detection]. The nblocks that can be acquired form the free list, ordered
 * by the best f of their OPEN lists: a worker keeps its nblock while it is
 * not worse than the best free one and than the best nblock its use
 * prevents, after min_expansions expansions at least, and otherwise
 * releases it and acquires the best free one.
 *
 * The interference scope of an nblock b is made of the predecessors, the
 * successors and the predecessors of the successors of b: the nblocks that
 * cannot be acquired while b is in use. \f$ \sigma( b ) \f$ counts the
 * nblocks in use whose interference scope contains b. A worker whose nblock
 * interferes with a better one marks it hot, and \f$ \sigma_h( b ) \f$
 * counts the hot nblocks with b in their scope: a nblock interfering with a
 * hot one is not free, and its worker releases it, so that every hot nblock
 * eventually becomes free, which is what prevents the livelock of plain PBNF
 * in an infinite state space. This is the pseudocode of appendix A of Burns
 * et al. (JAIR 2010); one lock protects the nblock graph (\f$ \sigma \f$,
 * \f$ \sigma_h \f$, the hot flags and the free list), and a worker with no
 * free nblock waits on a condition variable.
 *
 * The incumbent prunes the nodes with f not smaller, a goal is tested when
 * it is selected for expansion, and the search ends when no nblock is in
 * use and the free list is empty. The best f of the nblocks in use by
 * other workers, read to decide whether to switch, is a value they publish
 * after each expansion, which may be slightly out of date. */

template< NBlockDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Open = std::conditional_t<
           std::is_integral_v< cost_t< D > > ,
           BucketOpenList< Node< state_t< D > , cost_t< D > > > ,
           HeapOpenList< Node< state_t< D > , cost_t< D > > > > >
class SafePBNF {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the SafePBNF)

 explicit SafePBNF( const D & domain , const PBNFOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( f_options.min_expansions == 0 )
   f_options.min_expansions = 1;
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
  Graph g( f_domain );
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , g );

  {  // the root, in its nblock, which is the first free one
   const State s0 = f_domain.initial_state();
   const std::size_t b = f_domain.nblock( s0 );
   w[ 0 ]->insert( g.block( b ) , s0 , Cost( 0 ) ,
                   Cost( f_domain.heuristic( s0 ) ) , nullptr ,
                   f_domain.hash( s0 ) );
   g.make_free( b );
   }

  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  if( g.limit_hit.load() )
   sol.status = Status::node_limit;
  else
   if( g.goal ) {
    sol.status = Status::solved;
    sol.cost = g.goal->g;
    sol.path = trace_path( g.goal );
    }
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // an nblock: its OPEN and CLOSED lists and its place in the graph

 struct NBlock {
  Open open;                                 // its OPEN list
  NodeTable< NodeT > table{ 16 };            // its CLOSED list
  std::vector< std::size_t > interferes;     // its interference scope
  std::vector< std::size_t > scope;          // its successors (and itself)
  unsigned sigma = 0;                        // nblocks in use interfering
  unsigned sigma_h = 0;                      // hot nblocks interfering
  bool hot = false;                          // marked hot
  bool inuse = false;                        // acquired by a worker
  bool in_free = false;                      // in the free list
  Cost free_key{};                           // its key in the free list
  std::atomic< Cost > best{ infinite_cost< Cost >() };  // published best f
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the nblock graph, shared by the workers

 struct Graph {

  explicit Graph( const D & d )
   : f_d( d ) , v_b( new std::atomic< NBlock * >[ d.num_nblocks() ] ) {
   for( std::size_t i = 0 ; i < d.num_nblocks() ; ++i )
    v_b[ i ].store( nullptr , std::memory_order_relaxed );
   }

  // the nblock b if it has been built, nullptr otherwise
  NBlock * find( std::size_t b ) const {
   return( v_b[ b ].load( std::memory_order_acquire ) );
   }

  // the nblock b, built at its first use with its scopes
  NBlock & block( std::size_t b ) {
   if( NBlock * x = find( b ) )
    return( *x );
   std::lock_guard< std::mutex > guard( f_build );
   if( NBlock * x = v_b[ b ].load( std::memory_order_relaxed ) )
    return( *x );
   auto n = std::make_unique< NBlock >();
   std::vector< std::size_t > succ , iset;
   f_d.for_each_nblock_successor( b , [ & ]( std::size_t s ) {
    succ.push_back( s );
    } );
   // undirected abstract graphs: the predecessors are the successors
   for( std::size_t s : succ ) {
    iset.push_back( s );
    f_d.for_each_nblock_successor( s , [ & ]( std::size_t t ) {
     iset.push_back( t );
     } );
    }
   std::sort( iset.begin() , iset.end() );
   iset.erase( std::unique( iset.begin() , iset.end() ) , iset.end() );
   iset.erase( std::remove( iset.begin() , iset.end() , b ) , iset.end() );
   n->interferes = std::move( iset );
   n->scope = std::move( succ );
   n->scope.push_back( b );
   v_b[ b ].store( n.get() , std::memory_order_release );
   NBlock & ref = *n;
   v_own.push_back( std::move( n ) );
   return( ref );
   }

  // the best f of the OPEN list of b, the stale entries dropped; only for
  // an nblock that no worker can touch, or that the caller holds
  Cost best_f( NBlock & n ) {
   while( ! n.open.empty() && n.open.top().stale() )
    n.open.pop();
   return( n.open.empty() ? infinite_cost< Cost >() : n.open.top_f() );
   }

  // free: not in use, no interference, nodes below the incumbent
  bool is_free( NBlock & n ) {
   return( ! n.inuse && ( n.sigma == 0 ) && ( n.sigma_h == 0 ) &&
           ( best_f( n ) < incumbent.load() ) );
   }

  // puts b in the free list if it is free (lock held)
  void make_free( std::size_t b ) {
   NBlock & n = block( b );
   if( n.in_free || ! is_free( n ) )
    return;
   n.free_key = best_f( n );
   n.best.store( n.free_key , std::memory_order_relaxed );
   free.insert( { n.free_key , b } , true );
   n.in_free = true;
   wake.notify_all();
   }

  // takes b out of the free list (lock held)
  void make_busy( std::size_t b ) {
   NBlock & n = block( b );
   if( n.in_free ) {
    free.remove( { n.free_key , b } );
    n.in_free = false;
    }
   }

  // the f of the best free nblock (lock held)
  Cost best_free( void ) const {
   return( free.empty() ? infinite_cost< Cost >()
                        : ( *free.begin() ).first.first );
   }

  // the nblock of best published f among the interferers of n, or none
  std::pair< Cost , std::size_t > best_scope( const NBlock & n ) {
   std::pair< Cost , std::size_t > best{ infinite_cost< Cost >() ,
                                         std::size_t( -1 ) };
   for( std::size_t i : n.interferes )
    if( const NBlock * x = find( i ) ) {
     const Cost f = x->best.load( std::memory_order_relaxed );
     if( f < best.first )
      best = { f , i };
     }
   return( best );
   }

  // SetHot of appendix A (lock held)
  void set_hot( std::size_t b ) {
   NBlock & n = block( b );
   if( n.hot || ( n.sigma == 0 ) )
    return;
   for( std::size_t i : n.interferes )
    if( const NBlock * x = find( i ) ; x && x->hot &&
        ( x->best.load( std::memory_order_relaxed ) <
          n.best.load( std::memory_order_relaxed ) ) )
     return;                       // a better hot nblock interferes
   n.hot = true;
   for( std::size_t i : n.interferes ) {
    NBlock & m = block( i );
    if( m.hot )
     set_cold( i );
    if( m.in_free )
     make_busy( i );
    ++m.sigma_h;
    }
   }

  // SetCold of appendix A (lock held)
  void set_cold( std::size_t b ) {
   NBlock & n = block( b );
   n.hot = false;
   for( std::size_t i : n.interferes ) {
    NBlock & m = block( i );
    --m.sigma_h;
    if( ( m.sigma == 0 ) && ( m.sigma_h == 0 ) ) {
     if( m.hot )
      set_cold( i );
     make_free( i );
     }
    }
   }

  // Release of appendix A (lock held)
  void release( std::size_t b ) {
   NBlock & n = block( b );
   n.inuse = false;
   --inuse;
   for( std::size_t i : n.interferes ) {
    NBlock & m = block( i );
    if( --m.sigma == 0 ) {
     if( ( m.sigma_h == 0 ) && m.hot )
      set_cold( i );
     make_free( i );
     }
    }
   make_free( b );
   }

  // acquires b, taken from the free list (lock held)
  void acquire( std::size_t b ) {
   make_busy( b );
   NBlock & n = block( b );
   n.inuse = true;
   ++inuse;
   for( std::size_t i : n.interferes ) {
    NBlock & m = block( i );
    if( m.in_free )
     make_busy( i );
    ++m.sigma;
    }
   }

  // drops from the free list the nblocks with no node below the incumbent
  void prune_free( void ) {
   const Cost c = incumbent.load();
   for( auto it = free.lowerBound( { c , 0 } ) ; it != free.end() ; ) {
    block( ( *it ).first.second ).in_free = false;
    it = free.erase( it );
    }
   }

  const D & f_d;                                // the Domain
  std::unique_ptr< std::atomic< NBlock * >[] > v_b;  // the nblocks
  std::vector< std::unique_ptr< NBlock > > v_own;    // their ownership
  std::mutex f_build;                           // serializes block()
  Lock lock;                                    // protects the graph
  std::condition_variable_any wake;             // a free nblock, or the end
  ADS::BTreeMap< std::pair< Cost , std::size_t > , bool > free;  // free list
  unsigned inuse = 0;                           // the nblocks in use
  bool done = false;                            // the search is over
  std::atomic< Cost > incumbent{ infinite_cost< Cost >() };
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  NodeT * goal = nullptr;                       // the best goal found
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const SafePBNF & a , Graph & g ) : f_a( a ) , f_g( g ) {}

  // ThreadSearch of appendix A
  void run( void ) {
   std::size_t b = none;
   bool just_hot = false;
   for( ;; ) {
    b = next_nblock( b , just_hot );
    just_hot = false;
    if( b == none )
     return;
    NBlock & n = f_g.block( b );
    unsigned exp = 0;
    while( ! should_switch( b , n , exp , just_hot ) ) {
     const OpenEntry< NodeT > e = n.open.pop();
     if( e.stale() )
      continue;
     NodeT * m = e.node;
     if( m->g + m->h >= f_g.incumbent.load() ) {
      while( ! n.open.empty() )          // all the others are not better
       n.open.pop();
      break;
      }
     expand( n , m );
     ++exp;
     n.best.store( n.open.empty() ? infinite_cost< Cost >()
                                  : n.open.top_f() ,
                   std::memory_order_relaxed );
     if( f_a.f_options.node_limit &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      std::lock_guard< Lock > guard( f_g.lock );
      f_g.limit_hit.store( true );
      f_g.done = true;
      f_g.wake.notify_all();
      return;
      }
     }
    }
   }

  // ShouldSwitch of appendix A
  bool should_switch( std::size_t b , NBlock & n , unsigned & exp ,
                      bool & just_hot ) {
   while( ! n.open.empty() && n.open.top().stale() )
    n.open.pop();
   if( n.open.empty() )
    return( true );
   if( exp < f_a.f_options.min_expansions )
    return( false );
   exp = 0;
   const Cost mine = n.open.top_f();
   std::lock_guard< Lock > guard( f_g.lock );
   if( f_g.done )
    return( true );
   const auto [ fs , bs ] = f_g.best_scope( n );
   const Cost ff = f_g.best_free();
   if( ( ff < mine ) || ( fs < mine ) ) {
    if( ( bs != none ) && ( fs < ff ) ) {
     f_g.set_hot( bs );
     just_hot = true;
     }
    return( true );
    }
   for( std::size_t i : n.interferes )
    if( const NBlock * x = f_g.find( i ) ; x && x->hot )
     f_g.set_cold( i );
   (void) b;
   return( false );
   }

  // NextNBlock of appendix A
  std::size_t next_nblock( std::size_t b , bool just_hot ) {
   std::unique_lock< Lock > guard( f_g.lock , std::defer_lock );
   NBlock * n = b == none ? nullptr : & f_g.block( b );
   if( ! n || n->open.empty() || just_hot )
    guard.lock();
   else
    if( ! guard.try_lock() )
     return( b );
   if( n ) {
    const Cost mine = f_g.best_f( *n );
    if( ! f_g.done && ( mine < f_g.incumbent.load() ) &&
        ( mine <= f_g.best_scope( *n ).first ) &&
        ( mine <= f_g.best_free() ) && ! just_hot )
     return( b );
    n->best.store( mine , std::memory_order_relaxed );
    f_g.release( b );
    }
   f_g.prune_free();
   for( ;; ) {
    if( f_g.done )
     return( none );
    if( ! f_g.free.empty() )
     break;
    if( f_g.inuse == 0 ) {
     f_g.done = true;
     f_g.wake.notify_all();
     return( none );
     }
    f_g.wake.wait( guard );
    f_g.prune_free();
    }
   const std::size_t m = ( *f_g.free.begin() ).first.second;
   f_g.acquire( m );
   return( m );
   }

  // expands m, which is in n, or records it as the incumbent
  void expand( NBlock & , NodeT * m ) {
   m->closed = true;
   ++f_st.expanded;
   if( m->expanded )
    ++f_st.reexpanded;
   m->expanded = true;
   if( f_a.f_domain.is_goal( m->state ) ) {
    std::lock_guard< std::mutex > guard( f_g.goal_mutex );
    if( m->g < f_g.incumbent.load() ) {
     f_g.incumbent.store( m->g );
     f_g.goal = m;
     }
    return;
    }
   // first the successors, with their slots prefetched, then the lookups
   v_kids.clear();
   f_a.f_domain.for_each_successor( m->state , m->h ,
                                    [ & ]( const State & t , Cost c ,
                                           Cost ht ) {
    ++f_st.generated;
    const Cost g = m->g + c;
    if( g + ht >= f_g.incumbent.load( std::memory_order_relaxed ) )
     return;
    NBlock & cb = f_g.block( f_a.f_domain.nblock( t ) );
    const std::uint64_t hash = f_a.f_domain.hash( t );
    cb.table.prefetch( hash );
    v_kids.push_back( Kid{ t , g , ht , hash , & cb } );
    } );
   for( const Kid & k : v_kids )
    if( insert( *k.cb , k.s , k.g , k.h , m , k.hash ) )
     k.cb->best.store( std::min( k.cb->best.load( std::memory_order_relaxed ) ,
                                 k.g + k.h ) , std::memory_order_relaxed );
   }

  // records the path to t, whose hash is hash, of cost g in the nblock cb
  // (which the worker holds in its scope); true if t goes into OPEN
  bool insert( NBlock & cb , const State & t , Cost g , Cost h ,
               NodeT * parent , std::uint64_t hash ) {
   auto [ x , fresh ] = cb.table.find_or_insert( t , hash , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ t , g , h , parent , false , false };
    return( y );
    } );
   if( ! fresh ) {
    if( g >= x->g ) {
     ++f_st.duplicates;
     return( false );
     }
    x->g = g;
    x->parent = parent;
    x->closed = false;
    }
   cb.open.push( x );
   return( true );
   }

  static constexpr std::size_t none = std::size_t( -1 );

  // a successor: its state, g, h and hash, and its nblock
  struct Kid {
   State s;
   Cost g;
   Cost h;
   std::uint64_t hash;
   NBlock * cb;
   };

  std::vector< Kid > v_kids;                    // the successors of a node

  const SafePBNF & f_a;                         // the search
  Graph & f_g;                                  // the nblock graph
  NodePool< NodeT > f_pool;                     // the nodes it allocated
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 PBNFOptions f_options;                                ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( SafePBNF ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SafePBNF */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File SafePBNF.h ---------------------------*/
/*--------------------------------------------------------------------------*/
