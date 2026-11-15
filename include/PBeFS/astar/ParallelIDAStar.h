/*--------------------------------------------------------------------------*/
/*------------------------- File ParallelIDAStar.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ParallelIDAStar, the parallel
 * iterative-deepening A* that searches a frontier of the state space: a
 * best-first search (or a tree of fixed depth) stops at a set of open nodes,
 * the units, and the iterations of IDA* [see IDAStar.h] then run below the
 * units on the workers of an executor of ADS.
 *
 * With the notation of the page "Parallel A*" of the documentation, the
 * frontier is built by one of three rules (cf. FrontierRule): A* until
 * OPEN holds N live nodes (astar, as in the AIDA* of Reinefeld and
 * Schnecke, 1994, and in the A*+IDA* of Bu and Korf, 2019), uniform-cost
 * search until the same (breadth, a breadth-first frontier with unit
 * costs), or the tree of the smallest depth with N leaves, without
 * duplicate detection (depth, the work units of the parallel IDA* of
 * Sturtevant and of Siag et al., 2024). Each unit u has a threshold
 * \f$ t_u \f$, initially \f$ f( u ) \f$, and the cost-bounded depth-first
 * search of u with the bound \f$ \beta \geq t_u \f$ returns, besides a
 * goal of cost at most \f$ \beta \f$ if one exists below u, the least
 * \f$ f \f$ of the nodes it cuts, which becomes \f$ t_u \f$:
 *
 * \code{.unparsed}
 *  1  procedure PIDA*(s₀, p, N)
 *  2     build the units U by the frontier rule (return if it ends the search)
 *  3     t_u ← f(u) for each u ∈ U
 *  4     loop                                       ▷ the synchronous scheme
 *  5        β ← min_u t_u;  if β = ∞ then return unsolvable
 *  6        A ← {u ∈ U : t_u ≤ β}       ▷ all the units without unit bounds
 *  7        in parallel, the workers take the units of A
 *              (by work stealing, or from a shared cursor):
 *  8           (goal, t'_u) ← DFS(u, g(u), β);  t_u ← t'_u
 *  9           if goal then stop all the workers;  return its path
 *
 *  1  procedure AIDA*(s₀, p, N)                    ▷ the asynchronous scheme
 *  2     build U;  z ← ∞;  Q ← U by increasing t_u
 *  3     in parallel, each worker repeats until Q = ∅ and no worker searches:
 *  4        u ← the first unit of Q, removed;  if t_u ≥ z then drop it
 *  5        (goal, t'_u) ← DFS(u, g(u), t_u) cutting also f ≥ z
 *  6        if goal then z ← min(z, its cost), recording its path
 *  7        t_u ← t'_u;  if t_u < z then put u back in Q
 *  8     return the path of cost z, or unsolvable
 * \endcode
 *
 * With unit bounds (cf. PIDAOptions::unit_bounds), a unit whose threshold
 * exceeds \f$ \beta \f$ sits the iteration out: its search with
 * \f$ \beta \f$ would expand the same nodes as its last one, since no node
 * below it has \f$ f \f$ in \f$ ( \beta_{\mathrm{last}} , t_u ) \f$, so that
 * the sum over the iterations never exceeds that of the classical scheme,
 * which searches every unit at every bound. The units are scheduled by
 * work stealing, each worker taking first those it searched before and
 * stealing half of the range of another worker when it has none, or from a
 * cursor shared by all, in both cases in the decreasing order of the
 * expansions of their last search (the longest processing time first).
 *
 * Within a unit the depth-first search keeps an explicit stack of the
 * successors not yet tried, and with splitting (cf. PIDAOptions::split) a
 * worker that sees an idle one gives it half of the successors not yet
 * tried at the shallowest level of its stack, as in the parallel IDA* of
 * Rao and Kumar (1987); the pieces of a unit report their least
 * cut \f$ f \f$ to the unit, whose threshold is the least of them. With a
 * table of the frontier (rules astar and breadth), the first table_depth
 * levels below a unit are looked up in it, and a node found there with a
 * \f$ g \f$ not larger is pruned.
 *
 * The page "Parallel A*" proves the scheme optimal: with an admissible
 * \f$ h \f$, the open node of largest \f$ g \f$ on an optimal path, whose
 * \f$ g \f$ is optimal, keeps a threshold not above \f$ C^* \f$ until the
 * goal is found below it, and no pruning cuts the optimal suffix below it,
 * so that \f$ \beta \leq C^* \f$ at every synchronous iteration and every
 * goal found is optimal; in the asynchronous scheme, z is optimal once
 * every threshold is at least z.
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

#ifndef __PBeFS_ParallelIDAStar
 #define __PBeFS_ParallelIDAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <functional>

#include <memory>

#include <mutex>

#include <queue>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/Solution.h"

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
/// how the frontier of a ParallelIDAStar is built

enum class FrontierRule : std::uint8_t {
 astar ,     ///< A* until OPEN holds the units wanted
 breadth ,   ///< uniform-cost search (by g) until the same
 depth       ///< the tree of the least depth with that many leaves
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// how the units of an iteration are given to the workers

enum class UnitSchedule : std::uint8_t {
 steal ,     ///< a range of units per worker, half stolen when it is empty
 shared      ///< one cursor over the units, shared by all the workers
 };

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS PIDAOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of ParallelIDAStar
/** Those of the parallel searches, of which batch and flush_period are not
 * used, and node_limit bounds the expansions of each worker in the
 * iterations; deterministic needs the synchronous scheme and no splitting,
 * table_depth a frontier rule with a table, and the constructor of
 * ParallelIDAStar throws std::invalid_argument otherwise. */

struct PIDAOptions : public ParallelOptions {

 /// the units wanted (0: 128 per worker)
 std::size_t frontier = 0;

 /// the rule that builds the frontier
 FrontierRule rule = FrontierRule::astar;

 /// the scheduling of the units of an iteration
 UnitSchedule schedule = UnitSchedule::steal;

 /// a unit is searched only when the bound reaches its threshold; false
 /// searches every unit at every bound
 bool unit_bounds = true;

 /// no barrier between the iterations: the units are searched by
 /// increasing threshold, each with its own (AIDA*)
 bool asynchronous = false;

 /// with a value s > 0, a worker that sees an idle one gives it half of
 /// the successors not yet tried at the shallowest level of its stack,
 /// those whose f is at least s below the bound; 0 disables the splitting
 double split = 0;

 /// the levels below a unit whose nodes are looked up in the table of the
 /// frontier, and pruned if found there with a g not larger
 unsigned table_depth = 0;

 /// the path found is that of the first goal in the order of the units,
 /// the same at every run
 bool deterministic = false;

 };  // end( struct( PIDAOptions ) )

template< Domain D >
class ParallelIDAStarMPI;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ParallelIDAStar --------------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel iterative-deepening A* below a frontier of units
/** The counters of the Solution add those of the frontier (expanded,
 * generated, duplicates of its table) to those of the iterations, whose
 * expansions by worker are expanded_by; duplicates also counts the nodes
 * pruned by the table in the iterations, sent the units stolen and the
 * nodes given away by splitting, messages the steals and the splits, and
 * stored the nodes of the frontier plus the longest path of each worker.
 * iterations() and units() tell the iterations (the searches of units in
 * the asynchronous scheme) and the units of the last solve(). */

template< Domain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class ParallelIDAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the search)

 explicit ParallelIDAStar( const D & domain ,
                           const PIDAOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( f_options.weight != 1 )
   throw( std::invalid_argument( "ParallelIDAStar: no weight" ) );
  if( f_options.frontier == 0 )
   f_options.frontier = std::size_t( 128 ) * f_options.threads;
  if( f_options.deterministic &&
      ( f_options.asynchronous || ( f_options.split > 0 ) ) )
   throw( std::invalid_argument( "ParallelIDAStar: deterministic needs the "
                                 "synchronous scheme and no splitting" ) );
  if( ( f_options.table_depth > 0 ) &&
      ( f_options.rule == FrontierRule::depth ) )
   throw( std::invalid_argument( "ParallelIDAStar: the rule depth has no "
                                 "table" ) );
  if( f_options.split < 0 )
   throw( std::invalid_argument( "ParallelIDAStar: negative split" ) );
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
  Stopwatch clock;
  Frontier fr( f_options.frontier );
  SolutionT sol;
  f_iterations = 0;
  f_units = 0;
  if( f_options.rule == FrontierRule::depth )
   build_depth( fr , sol );
  else
   build_best_first( fr , sol );
  f_units = fr.units.size();
  if( ! fr.done )
   search( fr , sol );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the iterations of the last solve(), or its searches of units in the
 /// asynchronous scheme

 [[nodiscard]] std::uint64_t iterations( void ) const noexcept {
  return( f_iterations );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the units of the last solve()

 [[nodiscard]] std::size_t units( void ) const noexcept { return( f_units ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 template< Domain > friend class ParallelIDAStarMPI;     // the frontier

/*---------------------------- PRIVATE TYPES -------------------------------*/

 using Lock = typename Executor::lock_type;

 /// a unit: a node of the frontier and its threshold
 struct Unit {
  const NodeT * node;                    // its node, with the path to it
  Cost threshold;                        // the bound of its next search
  std::uint64_t work = 0;                // expansions of its last search
  };

 /// the frontier: its nodes, its table and its units
 struct Frontier {
  NodePool< NodeT > pool;
  NodeTable< NodeT > table;
  std::vector< Unit > units;
  bool done = false;                     // the frontier ended the search
  explicit Frontier( std::size_t n )
   : table( NodeTable< NodeT >::slots_for( 4 * n + 16 ) ) {}
  };

 /// a successor in the stack of the depth-first search
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint8_t op;
  };

 /// a level of the stack: the successors [ begin , end ) of the node of
 /// cost g, the next one to try at next
 struct Frame {
  std::size_t begin , end , next;
  Cost g;
  };

 /// a part of a unit given away by splitting: a successor not tried, with
 /// the states from the unit to it
 struct Piece {
  std::uint32_t unit;
  Cost bound;
  Cost g;
  Cost h;
  std::uint8_t op;
  std::vector< State > path;
  };

 /// the best goal found
 struct Best {
  Cost cost = infinite_cost< Cost >();
  std::uint64_t index = ~std::uint64_t( 0 );    // its unit
  std::vector< State > path;
  };

 /// the range [ b , e ) of positions of a worker, packed in a word
 struct alignas( cache_line ) Range {
  std::atomic< std::uint64_t > r{ 0 };
  };

/*-------------------------------- FRONTIER --------------------------------*/
 // the frontier by A* (by f) or uniform-cost search (by g), with
 // reopening: it stops when the live nodes of OPEN are N, ending the search
 // if a goal is expanded first or OPEN empties

 void build_best_first( Frontier & fr , SolutionT & sol ) const {
  SearchStats & st = sol.stats;
  const bool by_f = f_options.rule == FrontierRule::astar;
  struct Entry {
   Cost key;
   Cost g;
   std::uint64_t seq;
   NodeT * n;
   };
  // the entry first popped: least key, then (by f) largest g, then oldest
  const auto later = [ by_f ]( const Entry & a , const Entry & b ) {
   if( a.key != b.key )
    return( a.key > b.key );
   if( by_f && ( a.g != b.g ) )
    return( a.g < b.g );
   return( a.seq > b.seq );
   };
  std::priority_queue< Entry , std::vector< Entry > , decltype( later ) >
   open( later );
  std::uint64_t seq = 0;
  const auto push = [ & ]( NodeT * n ) {
   open.push( Entry{ by_f ? n->g + n->h : n->g , n->g , seq++ , n } );
   };
  const State s0 = f_domain.initial_state();
  NodeT * root = fr.table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * n = fr.pool.allocate();
   *n = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               false , false };
   return( n );
   } ).first;
  push( root );
  std::size_t live = 1;
  while( ! open.empty() && ( live < f_options.frontier ) ) {
   const Entry e = open.top();
   open.pop();
   NodeT * n = e.n;
   if( n->closed || ( n->g != e.g ) )
    continue;
   n->closed = true;
   --live;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( f_domain.is_goal( n->state ) ) {
    sol.status = Status::solved;
    sol.cost = n->g;
    sol.bound = sol.cost;
    sol.path = trace_path( n );
    fr.done = true;
    break;
    }
   for_each_child( f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    ++st.generated;
    const Cost g = n->g + c;
    auto [ m , fresh ] = fr.table.find_or_insert( t , f_domain.hash( t ) ,
                                                  [ & ] {
     NodeT * x = fr.pool.allocate();
     *x = NodeT{ t , g , ht , n , false , false , op };
     return( x );
     } );
    if( fresh ) {
     push( m );
     ++live;
     }
    else
     if( g < m->g ) {
      if( m->closed )
       ++live;
      m->g = g;
      m->parent = n;
      m->op = op;
      m->closed = false;
      push( m );
      }
     else
      ++st.duplicates;
    } );
   }
  if( ! fr.done ) {
   while( ! open.empty() ) {              // the live entries, in order
    const Entry e = open.top();
    open.pop();
    if( ! e.n->closed && ( e.n->g == e.g ) )   // one per live node
     fr.units.push_back( Unit{ e.n , e.n->g + e.n->h } );
    }
   if( fr.units.empty() ) {
    sol.status = Status::unsolvable;
    fr.done = true;
    }
   }
  st.stored = fr.pool.size();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the frontier as the tree of the least depth with N leaves, without
 // duplicate detection: a goal is a leaf, never expanded, so that the
 // iterations find the cheapest goal

 void build_depth( Frontier & fr , SolutionT & sol ) const {
  SearchStats & st = sol.stats;
  const State s0 = f_domain.initial_state();
  NodeT * root = fr.pool.allocate();
  *root = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
                 false , false };
  std::vector< NodeT * > level{ root } , next;
  while( level.size() < f_options.frontier ) {
   next.clear();
   bool grew = false;
   for( NodeT * n : level ) {
    if( f_domain.is_goal( n->state ) ) {
     next.push_back( n );
     continue;
     }
    ++st.expanded;
    n->expanded = true;
    grew = true;
    for_each_tree_child( f_domain , n->state , n->h , n->op ,
                         [ & ]( const State & t , Cost c , Cost ht ,
                                std::uint8_t op ) {
     ++st.generated;
     if constexpr( ! OperatorDomain< D > )
      if( n->parent && ( t == n->parent->state ) )
       return;
     NodeT * x = fr.pool.allocate();
     *x = NodeT{ t , n->g + c , ht , n , false , false , op };
     next.push_back( x );
     } );
    }
   if( ! grew )
    break;                                // only goals left
   level.swap( next );
   if( level.empty() )
    break;
   }
  for( NodeT * n : level )
   fr.units.push_back( Unit{ n , n->g + n->h } );
  if( fr.units.empty() ) {
   sol.status = Status::unsolvable;
   fr.done = true;
   }
  st.stored = fr.pool.size();
  }

/*-------------------------------- THE RUN ---------------------------------*/
 /// the state shared by the workers

 struct Shared {

  const Frontier & fr;
  std::vector< Unit > units;
  const unsigned p;

  // the next bound of each unit in its current search, the searches of
  // it not ended, and its expansions in them
  std::unique_ptr< std::atomic< Cost >[] > next;
  std::unique_ptr< std::atomic< int >[] > pending;
  std::unique_ptr< std::atomic< std::uint64_t >[] > work;

  // the units of the iteration, by position, and their claims
  std::vector< std::uint32_t > order;
  std::atomic< std::size_t > cursor{ 0 };
  std::unique_ptr< Range[] > ranges;
  Cost beta{};

  // the searches running or waiting (units claimed and pieces), the idle
  // workers that want a piece, the pieces queued
  std::atomic< std::int64_t > outstanding{ 0 };
  std::atomic< int > want{ 0 };
  std::atomic< std::size_t > queued{ 0 };
  Lock lock;                             // the pieces and the queue Q
  std::deque< Piece > pieces;
  std::priority_queue< std::pair< Cost , std::uint32_t > ,
                       std::vector< std::pair< Cost , std::uint32_t > > ,
                       std::greater<> > q;

  // the end: a goal (or the limit) stops the workers; in the synchronous
  // scheme, the units after cut are abandoned
  std::atomic< bool > stop{ false };
  std::atomic< bool > limit{ false };
  std::atomic< std::uint64_t > cut{ ~std::uint64_t( 0 ) };
  std::atomic< Cost > z;
  Lock best_lock;
  Best best;

  // the barrier between the iterations
  std::atomic< unsigned > arrived{ 0 };
  std::atomic< unsigned > generation{ 0 };
  bool finished = false;
  std::uint64_t iterations = 0;

  Shared( const Frontier & f , unsigned workers )
   : fr( f ) , units( f.units ) , p( workers ) ,
     next( new std::atomic< Cost >[ f.units.size() ] ) ,
     pending( new std::atomic< int >[ f.units.size() ] ) ,
     work( new std::atomic< std::uint64_t >[ f.units.size() ] ) ,
     ranges( new Range[ workers ] ) , z( infinite_cost< Cost >() ) {
   for( std::size_t i = 0 ; i < units.size() ; ++i ) {
    next[ i ].store( infinite_cost< Cost >() , std::memory_order_relaxed );
    pending[ i ].store( 0 , std::memory_order_relaxed );
    work[ i ].store( 0 , std::memory_order_relaxed );
    }
   }

  };  // end( struct( Shared ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the state of a worker

 struct alignas( cache_line ) Worker {

  unsigned id = 0;
  std::vector< Kid > kids;               // the successors on the stack
  std::vector< Frame > frames;           // the stack
  std::vector< State > path;             // the states from the unit
  std::uint64_t expanded = 0 , generated = 0 , duplicates = 0 ,
                sent = 0 , messages = 0 , longest = 0;

  };  // end( struct( Worker ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the iterations below the units of fr, on p workers

 void search( Frontier & fr , SolutionT & sol ) const {
  const unsigned p = f_options.threads;
  Shared sh( fr , p );
  std::vector< Worker > workers( p );
  for( unsigned i = 0 ; i < p ; ++i )
   workers[ i ].id = i;
  if( f_options.asynchronous )
   for( std::uint32_t i = 0 ; i < sh.units.size() ; ++i )
    sh.q.emplace( sh.units[ i ].threshold , i );
  else
   next_iteration( sh );
  if( ! sh.finished )
   run_workers< Executor >( f_ex , p , [ & ]( unsigned i ) {
    if( f_options.asynchronous )
     run_asynchronous( sh , workers[ i ] );
    else
     run_synchronous( sh , workers[ i ] );
    } );

  SearchStats & st = sol.stats;
  st.expanded_by.assign( p , 0 );
  st.stored_by.assign( p , 0 );
  for( const Worker & w : workers ) {
   st.expanded += w.expanded;
   st.generated += w.generated;
   st.duplicates += w.duplicates;
   st.sent += w.sent;
   st.messages += w.messages;
   st.stored += w.longest;
   st.expanded_by[ w.id ] = w.expanded;
   st.stored_by[ w.id ] = w.longest;
   }
  f_iterations = sh.iterations;
  if( sh.best.cost < infinite_cost< Cost >() ) {
   sol.status = Status::solved;
   sol.cost = sh.best.cost;
   sol.bound = sol.cost;
   sol.path = std::move( sh.best.path );
   }
  else
   sol.status = sh.limit.load() ? Status::node_limit : Status::unsolvable;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the synchronous scheme: the iterations, separated by a barrier whose
 // last worker prepares the next one

 void run_synchronous( Shared & sh , Worker & w ) const {
  for( ; ; ) {
   ADS::Backoff backoff;
   bool idle = false;
   for( ; ; ) {
    if( take_piece( sh , w , idle ) || take_unit( sh , w , idle ) ) {
     backoff.reset();
     continue;
     }
    if( sh.outstanding.load() == 0 )
     break;
    if( ! idle ) {
     idle = true;
     sh.want.fetch_add( 1 );
     }
    backoff.pause();
    }
   if( idle )
    sh.want.fetch_sub( 1 );
   const unsigned gen = sh.generation.load( std::memory_order_acquire );
   if( sh.arrived.fetch_add( 1 , std::memory_order_acq_rel ) + 1 == sh.p ) {
    next_iteration( sh );
    sh.arrived.store( 0 , std::memory_order_relaxed );
    sh.generation.store( gen + 1 , std::memory_order_release );
    }
   else {
    ADS::Backoff wait;
    while( sh.generation.load( std::memory_order_acquire ) == gen )
     wait.pause();
    }
   if( sh.finished )
    break;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the asynchronous scheme: the units by increasing threshold, until Q is
 // empty and no search runs

 void run_asynchronous( Shared & sh , Worker & w ) const {
  ADS::Backoff backoff;
  bool idle = false;
  for( ; ; ) {
   if( sh.stop.load( std::memory_order_relaxed ) )
    break;
   if( take_piece( sh , w , idle ) ) {
    backoff.reset();
    continue;
    }
   sh.lock.lock();
   while( ! sh.q.empty() &&
          ( sh.q.top().first >= sh.z.load( std::memory_order_relaxed ) ) )
    sh.q.pop();                          // all the rest is at least z
   if( ! sh.q.empty() ) {
    const std::uint32_t u = sh.q.top().second;
    sh.q.pop();
    sh.outstanding.fetch_add( 1 );
    ++sh.iterations;
    sh.lock.unlock();
    if( idle ) {
     idle = false;
     sh.want.fetch_sub( 1 );
     }
    start_unit( sh , w , u , sh.units[ u ].threshold );
    backoff.reset();
    continue;
    }
   const bool over = sh.outstanding.load() == 0;
   sh.lock.unlock();
   if( over )
    break;
   if( ! idle ) {
    idle = true;
    sh.want.fetch_add( 1 );
    }
   backoff.pause();
   }
  if( idle )
   sh.want.fetch_sub( 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the step of the last worker at the barrier: the thresholds of the units
 // searched, the next bound, and the units of the next iteration

 void next_iteration( Shared & sh ) const {
  const std::size_t n = sh.units.size();
  for( std::uint32_t u : sh.order ) {
   sh.units[ u ].threshold = sh.next[ u ].load( std::memory_order_relaxed );
   sh.units[ u ].work = sh.work[ u ].load( std::memory_order_relaxed );
   }
  if( sh.stop.load() || ( sh.best.cost < infinite_cost< Cost >() ) ) {
   sh.finished = true;
   return;
   }
  Cost beta = infinite_cost< Cost >();
  for( const Unit & u : sh.units )
   beta = std::min( beta , u.threshold );
  if( beta == infinite_cost< Cost >() ) {
   sh.finished = true;
   return;
   }
  sh.beta = beta;
  ++sh.iterations;
  std::vector< std::uint32_t > a;
  for( std::uint32_t u = 0 ; u < n ; ++u )
   if( f_options.unit_bounds ? sh.units[ u ].threshold <= beta
                             : sh.units[ u ].threshold <
                               infinite_cost< Cost >() )
    a.push_back( u );
  std::stable_sort( a.begin() , a.end() , [ & ]( std::uint32_t x ,
                                                 std::uint32_t y ) {
   return( sh.units[ x ].work > sh.units[ y ].work );
   } );
  sh.order.resize( a.size() );
  if( f_options.schedule == UnitSchedule::shared ) {
   std::copy( a.begin() , a.end() , sh.order.begin() );
   sh.cursor.store( 0 );
   }
  else {                                  // dealt to the workers in turn
   std::size_t pos = 0;
   for( unsigned i = 0 ; i < sh.p ; ++i ) {
    const std::size_t b = pos;
    for( std::size_t j = i ; j < a.size() ; j += sh.p )
     sh.order[ pos++ ] = a[ j ];
    sh.ranges[ i ].r.store( ( std::uint64_t( b ) << 32 ) |
                            std::uint64_t( pos ) );
    }
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // claims a unit of the iteration and searches it; false if none is left

 bool take_unit( Shared & sh , Worker & w , bool & idle ) const {
  sh.outstanding.fetch_add( 1 );
  std::size_t pos = ~std::size_t( 0 );
  if( f_options.schedule == UnitSchedule::shared ) {
   const std::size_t c = sh.cursor.fetch_add( 1 );
   if( c < sh.order.size() )
    pos = c;
   }
  else
   pos = claim( sh , w );
  if( pos == ~std::size_t( 0 ) ) {
   sh.outstanding.fetch_sub( 1 );
   return( false );
   }
  if( idle ) {
   idle = false;
   sh.want.fetch_sub( 1 );
   }
  const std::uint32_t u = sh.order[ pos ];
  if( f_options.deterministic &&
      ( u > sh.cut.load( std::memory_order_relaxed ) ) ) {
   sh.next[ u ].store( sh.units[ u ].threshold , std::memory_order_relaxed );
   sh.outstanding.fetch_sub( 1 );
   return( true );                         // abandoned, its bound kept
   }
  start_unit( sh , w , u , sh.beta );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the position of a unit taken from the range of w, or stolen from the
 // range of another worker (its upper half), or ~0 if all are empty

 std::size_t claim( Shared & sh , Worker & w ) const {
  std::atomic< std::uint64_t > & mine = sh.ranges[ w.id ].r;
  for( std::uint64_t x = mine.load() ; ; ) {
   const std::uint64_t b = x >> 32 , e = x & 0xFFFFFFFFu;
   if( b >= e )
    break;
   if( mine.compare_exchange_weak( x , ( ( b + 1 ) << 32 ) | e ) )
    return( std::size_t( b ) );
   }
  for( unsigned k = 1 ; k < sh.p ; ++k ) {
   std::atomic< std::uint64_t > & other = sh.ranges[ ( w.id + k ) % sh.p ].r;
   for( std::uint64_t x = other.load() ; ; ) {
    const std::uint64_t b = x >> 32 , e = x & 0xFFFFFFFFu;
    if( b >= e )
     break;
    const std::uint64_t m = e - ( e - b + 1 ) / 2;
    if( other.compare_exchange_weak( x , ( b << 32 ) | m ) ) {
     mine.store( ( ( m + 1 ) << 32 ) | e );
     ++w.messages;
     w.sent += e - m;
     return( std::size_t( m ) );
     }
    }
   }
  return( ~std::size_t( 0 ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // takes a piece given away by another worker and searches it; false if
 // there is none

 bool take_piece( Shared & sh , Worker & w , bool & idle ) const {
  if( sh.queued.load( std::memory_order_relaxed ) == 0 )
   return( false );
  sh.lock.lock();
  if( sh.pieces.empty() ) {
   sh.lock.unlock();
   return( false );
   }
  Piece pc = std::move( sh.pieces.front() );
  sh.pieces.pop_front();
  sh.queued.fetch_sub( 1 );
  sh.lock.unlock();
  if( idle ) {
   idle = false;
   sh.want.fetch_sub( 1 );
   }
  w.path = std::move( pc.path );
  ++w.generated;                          // never tried by its donor
  const Cost nb = dfs( sh , w , pc.unit , pc.g , pc.h , pc.op , pc.bound ,
                       nullptr );
  finish( sh , pc.unit , nb );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the search of the unit u with the bound b, claimed by w

 void start_unit( Shared & sh , Worker & w , std::uint32_t u , Cost b ) const
 {
  const NodeT * n = sh.units[ u ].node;
  sh.next[ u ].store( infinite_cost< Cost >() , std::memory_order_relaxed );
  sh.work[ u ].store( 0 , std::memory_order_relaxed );
  sh.pending[ u ].store( 1 );
  w.path.assign( 1 , n->state );
  const Cost nb = dfs( sh , w , u , n->g , n->h , n->op , b ,
                       n->parent ? & n->parent->state : nullptr );
  finish( sh , u , nb );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the end of a search of the unit u (or of one of its pieces), whose
 // least cut f is nb: the last one puts the unit back in Q (asynchronous)

 void finish( Shared & sh , std::uint32_t u , Cost nb ) const {
  atomic_min( sh.next[ u ] , nb );
  if( sh.pending[ u ].fetch_sub( 1 ) == 1 )
   if( f_options.asynchronous ) {
    Unit & x = sh.units[ u ];
    x.threshold = sh.next[ u ].load();
    x.work = sh.work[ u ].load();
    if( x.threshold < sh.z.load() ) {
     sh.lock.lock();
     sh.q.emplace( x.threshold , u );
     sh.lock.unlock();
     }
    }
  sh.outstanding.fetch_sub( 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the cost-bounded depth-first search below the last state of w.path,
 // generated by op with cost g and heuristic h, in the unit u with the
 // bound b; above is the state before the unit (nullptr if none or if
 // w.path holds it); returns the least f of the nodes cut

 Cost dfs( Shared & sh , Worker & w , std::uint32_t u , Cost g , Cost h ,
           std::uint8_t op , Cost b , const State * above ) const {
  const std::size_t off = w.path.size() - 1;     // the root in w.path
  const std::size_t base = off;           // its depth below the unit
  Cost nb = infinite_cost< Cost >();
  std::uint64_t mine = 0;                 // expansions of this search
  const bool async = f_options.asynchronous;
  const auto over_z = [ & ]( Cost f ) {
   return( async && ( f >= sh.z.load( std::memory_order_relaxed ) ) );
   };

  // the root: the end, the cut, the table, the goal, then its expansion
  if( sh.stop.load( std::memory_order_relaxed ) ||
      ( f_options.deterministic &&
        ( u > sh.cut.load( std::memory_order_relaxed ) ) ) )
   return( nb );
  if( ( g + h > b ) || over_z( g + h ) )
   return( g + h );
  if( ( base > 0 ) && pruned( sh , w , w.path.back() , g , base ) )
   return( nb );
  if( f_domain.is_goal( w.path.back() ) ) {
   goal( sh , w , u , g );
   return( nb );
   }
  w.kids.clear();
  w.frames.clear();
  if( ! expand( sh , w , g , h , op , mine ) )
   return( nb );

  while( ! w.frames.empty() ) {
   Frame & F = w.frames.back();
   if( F.next == F.end ) {
    w.kids.resize( F.begin );
    w.frames.pop_back();
    w.path.pop_back();
    continue;
    }
   const Kid k = w.kids[ F.next++ ];
   ++w.generated;
   const std::size_t at = off + w.frames.size() - 1;    // k's parent
   if constexpr( ! OperatorDomain< D > ) {
    const State * gp = at > 0 ? & w.path[ at - 1 ] : above;
    if( gp && ( k.s == *gp ) )
     continue;                            // back to the parent
    }
   const Cost gk = F.g + k.c , fk = gk + k.h;
   if( fk > b ) {
    nb = std::min( nb , fk );
    continue;
    }
   if( over_z( fk ) )
    continue;
   const std::size_t depth = base + w.frames.size();
   if( pruned( sh , w , k.s , gk , depth ) )
    continue;
   if( f_domain.is_goal( k.s ) ) {
    w.path.push_back( k.s );
    goal( sh , w , u , gk );
    w.path.pop_back();
    if( ! async )
     break;
    continue;
    }
   if( sh.stop.load( std::memory_order_relaxed ) ||
       ( f_options.deterministic &&
         ( u > sh.cut.load( std::memory_order_relaxed ) ) ) )
    break;
   if( f_options.node_limit && ( w.expanded >= f_options.node_limit ) ) {
    sh.limit.store( true );
    sh.stop.store( true );
    break;
    }
   if( ( f_options.split > 0 ) &&
       ( sh.want.load( std::memory_order_relaxed ) >
         int( sh.queued.load( std::memory_order_relaxed ) ) ) )
    give( sh , w , u , b );
   w.path.push_back( k.s );
   expand( sh , w , gk , k.h , k.op , mine );
   }
  w.frames.clear();
  sh.work[ u ].fetch_add( mine , std::memory_order_relaxed );
  return( nb );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // expands the last state of w.path, of cost g, pushing a frame of its
 // successors; true unless the search has to stop

 bool expand( Shared & sh , Worker & w , Cost g , Cost h , std::uint8_t op ,
              std::uint64_t & mine ) const {
  if( f_options.node_limit && ( w.expanded >= f_options.node_limit ) ) {
   sh.limit.store( true );
   sh.stop.store( true );
   return( false );
   }
  ++w.expanded;
  ++mine;
  w.longest = std::max< std::uint64_t >( w.longest , w.path.size() );
  const std::size_t b = w.kids.size();
  for_each_tree_child( f_domain , w.path.back() , h , op ,
                       [ & ]( const State & t , Cost c , Cost ht ,
                              std::uint8_t o ) {
   w.kids.push_back( Kid{ t , c , ht , o } );
   } );
  w.frames.push_back( Frame{ b , w.kids.size() , b , g } );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // true if s, of cost g at the given depth below its unit, is in the
 // table of the frontier with a g not larger

 bool pruned( Shared & sh , Worker & w , const State & s , Cost g ,
              std::size_t depth ) const {
  if( depth > f_options.table_depth )
   return( false );
  const NodeT * m = sh.fr.table.find( s , f_domain.hash( s ) );
  if( m && ( m->g <= g ) ) {
   ++w.duplicates;
   return( true );
   }
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a goal of cost g at the end of w.path, in the unit u

 void goal( Shared & sh , Worker & w , std::uint32_t u , Cost g ) const {
  std::lock_guard< Lock > guard( sh.best_lock );
  const bool better = f_options.deterministic ? u < sh.best.index
                                              : g < sh.best.cost;
  if( better ) {
   sh.best.cost = g;
   sh.best.index = u;
   sh.best.path = trace_path( sh.units[ u ].node );
   sh.best.path.insert( sh.best.path.end() , w.path.begin() + 1 ,
                        w.path.end() );
   }
  if( f_options.asynchronous )
   atomic_min( sh.z , g );
  else
   if( f_options.deterministic )
    atomic_min( sh.cut , std::uint64_t( u ) );
   else
    sh.stop.store( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // gives the idle workers half of the successors not tried at the
 // shallowest level that has some, those of f at least split below b

 void give( Shared & sh , Worker & w , std::uint32_t u , Cost b ) const {
  const std::size_t off = w.path.size() - w.frames.size();
  for( std::size_t l = 0 ; l < w.frames.size() ; ++l ) {
   Frame & F = w.frames[ l ];
   const std::size_t left = F.end - F.next;
   const bool top = l + 1 == w.frames.size();
   if( left < ( top ? 2u : 1u ) )
    continue;
   std::size_t count = top ? left / 2 : ( left + 1 ) / 2;
   while( ( count > 0 ) &&
          ( double( b ) - double( F.g + w.kids[ F.end - count ].c +
                                  w.kids[ F.end - count ].h )
            < f_options.split ) )
    --count;
   if( count == 0 )
    continue;
   std::vector< Piece > out;
   out.reserve( count );
   for( std::size_t i = F.end - count ; i < F.end ; ++i ) {
    const Kid & k = w.kids[ i ];
    Piece pc{ u , b , F.g + k.c , k.h , k.op ,
              std::vector< State >( w.path.begin() ,
                                    w.path.begin() + std::ptrdiff_t(
                                     off + l + 1 ) ) };
    pc.path.push_back( k.s );
    out.push_back( std::move( pc ) );
    }
   F.end -= count;
   sh.pending[ u ].fetch_add( int( count ) );
   sh.outstanding.fetch_add( std::int64_t( count ) );
   sh.lock.lock();
   for( Piece & pc : out )
    sh.pieces.push_back( std::move( pc ) );
   sh.queued.fetch_add( count );
   sh.lock.unlock();
   w.sent += count;
   ++w.messages;
   return;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // x ← min( x , v ), atomically

 template< class T >
 static void atomic_min( std::atomic< T > & x , T v ) noexcept {
  for( T c = x.load( std::memory_order_relaxed ) ; v < c ; )
   if( x.compare_exchange_weak( c , v ) )
    break;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 PIDAOptions f_options;                              ///< the options
 Executor * f_ex = nullptr;                          ///< the Executor
 mutable std::uint64_t f_iterations = 0;     ///< iterations of last solve
 mutable std::size_t f_units = 0;            ///< units of last solve

 };  // end( class( ParallelIDAStar ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with ParallelIDAStar on the threads of ADS

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > parallel_idastar(
 const D & domain , const PIDAOptions & options = {} )
{
 return( ParallelIDAStar< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ParallelIDAStar */

/*--------------------------------------------------------------------------*/
/*----------------------- End File ParallelIDAStar.h -----------------------*/
/*--------------------------------------------------------------------------*/
