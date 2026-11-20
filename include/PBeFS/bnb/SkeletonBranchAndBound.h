/*--------------------------------------------------------------------------*/
/*--------------------- File SkeletonBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for SkeletonBranchAndBound, the parallel Branch-and-Bound
 * (B&B) whose workers search depth-first and share the work by one of the
 * four rules of the search skeletons of Archibald et al. (2018, 2020): the
 * depth-bounded spawning, the budget spawning, the stack stealing and the
 * ordered tasks. The problem, the options and the Solution are those of
 * BranchAndBound.h and ParallelBranchAndBound.h, and the design is
 * discussed in docs/branch-and-bound.md.
 *
 * A worker searches a task (a node and its subtree) depth-first with a
 * stack of frames, each holding an evaluated node and the cursor of its
 * children [see BBCursor], which are built only when the search reaches
 * them; with the lazy children of the evaluator [see LazyBBEvaluator] a
 * child that is pruned with BBOptions::prune_level ends its siblings, which
 * are never built. With the notation of the page "Parallel branch-and-bound"
 * of the documentation and the procedures of the sequential search [see
 * BranchAndBound.h], EVALUATE(n) being lines 17 to 22 of PROCESS, the search
 * is the following; \f$ d \f$ is the spawn depth, \f$ b \f$ the budget and
 * \f$ P_i \f$ the pool of the tasks of the worker \f$ i \f$:
 *
 * \code{.unparsed}
 *  1  procedure SKELETON(R)
 *  2     with ordered: T ← the nodes of depth d under R, by discrepancy;
 *           the nodes above them evaluated and branched
 *     else the k-th node of R in P_(k mod p)
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        while t ← NEXT(i) ≠ nil do DFS(t)
 *  5  procedure DFS(t)
 *  6     S ← ();  VISIT(t);  c ← 0                 ▷ c: the backtracks
 *  7     while S ≠ () do
 *  8        if a limit is reached then the rest of S are leaves (open); stop
 *  9        with stack stealing, if a thief j asks, move to P_j the next
 *              child (all of them with steal_all) of the lowest frame of S,
 *              the top excepted, that has one
 * 10        with budget, if c = b then move the children left in the
 *              frames of S, the top excepted, to P_i;  c ← 0
 * 11        (m, e) ← the next child of the top frame of S, if any
 * 12        if there is none then pop S;  c ← c + 1;  continue
 * 13        β(m) ← max(β(top), e);  if PRUNABLE(β(m)) then m is a leaf
 *              (pruned), and so are the next children with prune_level
 * 14        else VISIT(m)
 * 15  procedure VISIT(n)
 * 16     if EVALUATE(n) makes n a leaf then return
 * 17     if depth-bounded and depth(n) < d then move all the children of n
 *           that survive to P_i
 * 18     else push n and the cursor of its children on S
 * 19  procedure NEXT(i)
 * 20     loop
 * 21        if ordered then return the next task of T, or nil
 * 22        if P_i ≠ ∅ then return its deepest task, removed
 * 23        if a pool P_v ≠ ∅ then return its shallowest task, removed
 * 24        with stack stealing, ask a busy worker v for its work
 * 25        if no task is left anywhere then return nil
 * \endcode
 *
 * The tasks are counted when they are made and when their subtree is done,
 * so that NEXT ends (line 25) when the count is 0, which is exact; the
 * steal of line 9 is answered by the victim at its next step, i.e., after
 * at most the evaluation of a node, and a worker that waits for an answer
 * answers the requests it gets with nothing. With BasicParallelBBOptions::
 * indexed, a task is not a node but the path of the choices (the index of
 * the child at each level) from the node of R it comes from, and a worker
 * rebuilds it before VISIT, by child( n , k ) of the evaluator if it has
 * one [see IndexedBBEvaluator], or else by evaluating each node of the
 * path, which BBStats::recomputed counts; a task then takes a few bytes
 * per level, whatever the node. The ordered tasks are taken in the order of
 * T by all the workers, the worker 0 included, so that, as in the Ordered
 * skeleton, the tasks done by the time a task is started are always those
 * before it; by discrepancy, the priority of a task is the sum of the
 * indices of its choices, ties broken left to right.
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

#ifndef __PBeFS_SkeletonBranchAndBound
 #define __PBeFS_SkeletonBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstdint>

#include <deque>

#include <exception>

#include <memory>

#include <mutex>

#include <numeric>

#include <optional>

#include <random>

#include <utility>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ ENUM BBRule -------------------------------*/
/*--------------------------------------------------------------------------*/
/// how the workers of SkeletonBranchAndBound share the work

//! [branch-and-bound-sig-rule]
enum class BBRule {
 depth_bounded ,   ///< the children of the nodes above a depth are tasks
 budget ,          ///< after a budget of backtracks, the rest are tasks
 stack_stealing ,  ///< an idle worker takes children from a busy one
 ordered           ///< the nodes of a depth, by a fixed priority
 };
//! [branch-and-bound-sig-rule]

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS SkeletonBranchAndBound ----------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel depth-first B&B of the search skeletons
/** Each worker searches its tasks depth-first, visiting the children of a
 * node in the order the evaluator gives them, and the work is shared by
 * the rule of the options [see BBRule and BasicParallelBBOptions]:
 * depth_bounded spawns as tasks the children of the nodes of depth less
 * than spawn_depth, budget spawns the children left on the stack of a
 * worker each time it has backtracked budget times, stack_stealing lets an
 * idle worker take from a busy one the next child of its lowest frame (or
 * all its children, with steal_all), and ordered makes tasks of the nodes
 * of depth spawn_depth, which all the workers take in the order of their
 * discrepancy (or left to right). A task is the node itself, or, with
 * indexed, the path of the choices that lead to it, which the worker that
 * takes it rebuilds. The incumbent is shared; the selection of the options
 * is not used, the order being depth-first, and with prune_level the lazy
 * children pruned when generated end their siblings [see BBOptions]. The
 * counter sent of the statistics holds the tasks that moved between the
 * workers, and messages the successful steals. */

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class SkeletonBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search), the
 /// rule by which the work is shared and the options

 explicit SkeletonBranchAndBound( const P & problem ,
                                  BBRule rule = BBRule::stack_stealing ,
                                  const Options & options = {} )
  : f_p( problem ) , f_rule( rule ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  return( run( std::nullopt ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes, returning the Solution [see
 /// BranchAndBoundSearch::solve()]; the node \f$ i \f$ goes to the worker
 /// \f$ i \bmod p \f$

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Item = typename BBQueue< Node , Value >::Item;
 using W = decltype( detail::bb_worker( std::declval< const P & >() , 0u ) );
 using Cursor = decltype( detail::bb_children< P >(
                                       std::declval< W & >() ,
                                       std::declval< const Node & >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a task: a node, or the path to it from the start it comes from

 struct Task {
  std::optional< Node > node;          // the node, unless indexed
  Value bound;                         // its bound, in the minimization
  unsigned depth;                      // its depth
  std::uint32_t start = 0;             // the start of the path
  std::vector< std::uint32_t > path;   // the choices from there
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the tasks of a worker, by their depth: the owner takes the deepest
 // first, a thief the shallowest, and each level is first in, first out

 struct Pool {

  void push( Task && t ) {
   std::lock_guard< std::mutex > guard( m );
   if( t.depth >= level.size() )
    level.resize( t.depth + 1 );
   level[ t.depth ].push_back( std::move( t ) );
   size.fetch_add( 1 , std::memory_order_release );
   }

  std::optional< Task > pop( bool deepest ) {
   if( size.load( std::memory_order_acquire ) == 0 )
    return( std::nullopt );
   std::lock_guard< std::mutex > guard( m );
   const std::size_t n = level.size();
   for( std::size_t j = 0 ; j < n ; ++j ) {
    auto & q = level[ deepest ? n - 1 - j : j ];
    if( ! q.empty() ) {
     std::optional< Task > t( std::in_place , std::move( q.front() ) );
     q.pop_front();
     size.fetch_sub( 1 , std::memory_order_release );
     return( t );
     }
    }
   return( std::nullopt );
   }

  // calls f on each task, the deepest first, and empties the pool
  template< class F >
  void drain( F && f ) {
   std::lock_guard< std::mutex > guard( m );
   for( std::size_t j = level.size() ; j-- > 0 ; ) {
    for( Task & t : level[ j ] )
     f( t );
    level[ j ].clear();
    }
   size.store( 0 , std::memory_order_release );
   }

  std::mutex m;                                  // protects level
  std::vector< std::deque< Task > > level;       // the tasks, by depth
  std::atomic< std::size_t > size{ 0 };          // their number
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the request of a steal: 0 none, j + 1 from the worker j, or claimed

 static constexpr unsigned c_claimed = ~0u;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {

  Shared( unsigned np , Value z ) : p( np ) , inc( z ) , pool( np ) ,
   ask( np ) , answered( np ) , answer( np ) {}

  unsigned p;                                   // the number of workers
  BBIncumbent< Node , Value > inc;              // the incumbent
  std::vector< Pool > pool;                     // the tasks
  std::vector< Leaf > starts;                   // the starts of the paths
  std::vector< Task > ordered;                  // the ordered tasks
  alignas( cache_line ) std::atomic< std::size_t > next{ 0 };  // ... next
  alignas( cache_line ) std::atomic< std::int64_t > pending{ 0 };  // tasks
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< bool > stop{ false };            // a limit was reached
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  std::vector< std::atomic< unsigned > > ask;   // the steal requests
  std::vector< std::atomic< bool > > answered;  // a steal was answered
  std::vector< std::vector< Task > > answer;    // ... with these tasks
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const Options & o = f_options;
  const unsigned p = o.threads;
  Stopwatch clock;
  Shared sh( p , o.incumbent ? detail::bb_signed< P >( *o.incumbent ) :
                               inf );
  sh.starts = detail::bb_starts( f_p , std::move( seeds ) );

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i , clock );

  if( f_rule == BBRule::ordered )
   w[ 0 ]->order();                    // the tasks, by the worker 0
  else
   for( std::size_t k = 0 ; k < sh.starts.size() ; ++k ) {
    const Leaf & l = sh.starts[ k ];
    Task t{ std::nullopt , l.bound , l.depth , std::uint32_t( k ) , {} };
    if( ! o.indexed )
     t.node.emplace( l.node );
    sh.pending.fetch_add( 1 );
    sh.pool[ k % p ].push( std::move( t ) );
    }

  std::exception_ptr error;
  std::mutex error_m;
  run_workers( f_ex , p , [ & ]( unsigned i ) {
   try {
    w[ i ]->run();
    }
   catch( ... ) {
    std::lock_guard< std::mutex > guard( error_m );
    if( ! error )
     error = std::current_exception();
    sh.stop.store( true );
    }
   } );
  if( error )
   std::rethrow_exception( error );

  // the tasks left by a limit are open leaves
  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( o );
  for( unsigned i = 0 ; i < p ; ++i ) {
   w[ i ]->leave();
   sol.stats.add( w[ i ]->f_st );
   fr.merge( std::move( w[ i ]->f_fr ) );
   }
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() , o );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its evaluator, its stack, its counters and its leaves

 struct alignas( cache_line ) Worker {

  // a frame of the stack: an evaluated node and the cursor of its
  // children, which may refer to the node, since a frame does not move
  struct Frame {
   Frame( Node && n , Value b , unsigned d , BBLeafKind k )
    : node( std::move( n ) ) , bound( b ) , depth( d ) , kind( k ) {}
   Node node;                         // the node
   Value bound;                       // its bound
   unsigned depth;                    // its depth
   BBLeafKind kind;                   // its kind, if it has no child
   std::optional< Cursor > cur;       // its children
   std::uint32_t k = 0;               // the index of the next one
   bool any = false;                  // true if one was generated
   bool cut = false;                  // true if the rest are pruned
   };

  Worker( const SkeletonBranchAndBound & a , Shared & s , unsigned i ,
          const Stopwatch & c )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_clock( c ) ,
     f_w( detail::bb_worker( a.f_p , i ) ) ,
     f_gen( a.f_options.seed + i ) , f_fr( a.f_options ) {}

  // the incumbent and its improvement
  Value inc( void ) const { return( f_sh.inc.value() ); }

  bool improve( Value v , const Node & n ) {
   return( f_sh.inc.template improve< P >( v , n , f_a.f_options ) );
   }

  // the main loop
  void run( void ) {
   while( std::optional< Task > t = next() ) {
    dfs( std::move( *t ) );
    f_sh.pending.fetch_sub( 1 , std::memory_order_acq_rel );
    }
   answer_nothing();
   }

  // true if a limit is reached, which stops all the workers
  bool stop( void ) {
   if( f_sh.stop.load( std::memory_order_relaxed ) )
    return( true );
   if( ! detail::bb_limit( f_a.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.stop.store( true );
   return( true );
   }

  // the node of a task: its own, or the one rebuilt from its path
  Node rebuild( Task & t ) {
   if( t.node )
    return( std::move( *t.node ) );
   Node n = f_sh.starts[ t.start ].node;
   for( const std::uint32_t k : t.path ) {
    if constexpr( IndexedBBEvaluator< W , Node , Value > ) {
     auto c = f_w.child( std::as_const( n ) , std::size_t( k ) );
     n = std::move( c->node );
     }
    else {
     ( void ) f_w.evaluate( n );
     ++f_st.recomputed;
     auto cur = detail::bb_children< P >( f_w , std::as_const( n ) );
     for( std::uint32_t j = 0 ; j < k ; ++j )
      ( void ) cur.next();
     auto c = cur.next();
     n = std::move( c->node );
     }
    }
   return( n );
   }

  // a task of the child of index k of the frame of the stack at level j
  // (the last, at the top, if j is the size of the path), of bound b
  Task task_of( Node && c , Value b , unsigned depth , std::size_t j ,
                std::uint32_t k ) {
   Task t{ std::nullopt , b , depth , f_start , {} };
   if( f_a.f_options.indexed ) {
    t.path = v_base;
    t.path.insert( t.path.end() , v_path.begin() ,
                   v_path.begin() + std::ptrdiff_t( j ) );
    t.path.push_back( k );
    }
   else
    t.node.emplace( std::move( c ) );
   return( t );
   }

  // spawns a task in the pool of the worker
  void spawn( Task && t ) {
   f_sh.pending.fetch_add( 1 , std::memory_order_acq_rel );
   f_sh.pool[ f_id ].push( std::move( t ) );
   }

  // the depth-first search of a task
  void dfs( Task && t ) {
   if( stop() ) {
    Node n = rebuild( t );
    f_fr.add( std::move( n ) , t.bound , t.depth , BBLeafKind::open );
    return;
    }
   Node n = rebuild( t );
   f_start = t.start;
   v_base = std::move( t.path );
   v_path.clear();
   v_stack.clear();
   std::uint64_t backtracks = 0;
   visit( std::move( n ) , t.bound , t.depth );
   const auto & o = f_a.f_options;
   const BBRule rule = f_a.f_rule;
   while( ! v_stack.empty() ) {
    if( stop() ) {
     leave_stack();
     return;
     }
    if( rule == BBRule::stack_stealing )
     answer();
    if( ( rule == BBRule::budget ) && ( backtracks >= o.budget ) ) {
     backtracks = 0;
     for( std::size_t j = 0 ; j + 1 < v_stack.size() ; ++j )
      while( std::optional< Task > s = take( j ) )
       spawn( std::move( *s ) );
     }
    Frame & f = v_stack.back();
    std::optional< BBChild< Node , Value > > c;
    if( ! f.cut )
     c = f.cur->next();
    if( ! c ) {
     if( ! f.any )                     // no children: its subtree is n
      f_fr.add( std::move( f.node ) , f.bound , f.depth , f.kind , false );
     else
      if( f.cut )
       drain_pruned( f );
     v_stack.pop_back();
     if( ! v_path.empty() )
      v_path.pop_back();
     ++backtracks;
     continue;
     }
    f.any = true;
    const std::uint32_t k = f.k++;
    const Value fb = f.bound;
    const unsigned d = f.depth + 1;
    const std::optional< Value > cb = detail::bb_child_bound< P >(
     fb , Value( c->estimate ) , inc() , o , f.cut , f_st , f_fr );
    if( ! cb ) {
     if( f_fr.kept( BBLeafKind::pruned ) )
      f_fr.add( std::move( c->node ) ,
                std::max( fb , detail::bb_signed< P >(
                                Value( c->estimate ) ) ) ,
                d , BBLeafKind::pruned );
     continue;
     }
    const std::size_t before = v_stack.size();
    v_path.push_back( k );
    visit( std::move( c->node ) , *cb , d );
    if( v_stack.size() == before )     // a leaf, or its children spawned
     v_path.pop_back();
    }
   }

  // evaluates a node, and pushes it on the stack if it is to be branched,
  // or spawns its children with depth_bounded above the spawn depth
  void visit( Node && n , Value bound , unsigned depth ) {
   const auto & o = f_a.f_options;
   BBLeafKind kind = BBLeafKind::pruned;
   const std::uint64_t e = f_st.evaluated;
   const std::optional< Value > b = detail::bb_evaluate< P >( n , bound ,
                                    depth , f_w , [ & ] { return( inc() ); } ,
                                    [ & ]( Value v , const Node & x ) {
                                     return( improve( v , x ) );
                                     } , o , f_st , f_fr , kind );
   if( f_st.evaluated != e )
    f_sh.evaluated.fetch_add( 1 , std::memory_order_relaxed );
   if( ! b )
    return;
   ++f_st.branched;
   detail::bb_cutoff< P >( f_w , inc() , o );
   Frame & f = v_stack.emplace_back( std::move( n ) , *b , depth , kind );
   f.cur.emplace( detail::bb_children< P >( f_w , std::as_const( f.node ) ) );
   if( ( f_a.f_rule == BBRule::depth_bounded ) &&
       ( depth < o.spawn_depth ) ) {
    const std::size_t j = v_stack.size() - 1;
    while( std::optional< Task > s = take( j ) )
     spawn( std::move( *s ) );
    if( ! f.any )
     f_fr.add( std::move( f.node ) , f.bound , f.depth , f.kind , false );
    else
     if( f.cut )
      drain_pruned( f );
    v_stack.pop_back();
    }
   }

  // the next child that survives of the frame at level j, as a task, or
  // none if there is no more; the pruned ones are leaves
  std::optional< Task > take( std::size_t j ) {
   Frame & f = v_stack[ j ];
   const auto & o = f_a.f_options;
   while( ! f.cut ) {
    std::optional< BBChild< Node , Value > > c = f.cur->next();
    if( ! c )
     return( std::nullopt );
    f.any = true;
    const std::uint32_t k = f.k++;
    const std::optional< Value > cb = detail::bb_child_bound< P >(
     f.bound , Value( c->estimate ) , inc() , o , f.cut , f_st , f_fr );
    if( cb )
     return( task_of( std::move( c->node ) , *cb , f.depth + 1 , j , k ) );
    if( f_fr.kept( BBLeafKind::pruned ) )
     f_fr.add( std::move( c->node ) ,
               std::max( f.bound , detail::bb_signed< P >(
                                    Value( c->estimate ) ) ) ,
               f.depth + 1 , BBLeafKind::pruned );
    }
   return( std::nullopt );
   }

  // the children left in a frame cut by prune_level: with the leaves kept,
  // each is a leaf (pruned); without, they are never built
  void drain_pruned( Frame & f ) {
   if( ! f_fr.kept( BBLeafKind::pruned ) )
    return;
   while( std::optional< BBChild< Node , Value > > c = f.cur->next() ) {
    ++f_st.generated;
    ++f_st.pruned_generated;
    f_fr.add( std::move( c->node ) ,
              std::max( f.bound , detail::bb_signed< P >(
                                   Value( c->estimate ) ) ) ,
              f.depth + 1 , BBLeafKind::pruned );
    }
   }

  // after a limit, the children left on the stack are open leaves
  void leave_stack( void ) {
   while( ! v_stack.empty() ) {
    Frame & f = v_stack.back();
    while( ! f.cut ) {
     std::optional< BBChild< Node , Value > > c = f.cur->next();
     if( ! c )
      break;
     f.any = true;
     const Value cb = std::max( f.bound , detail::bb_signed< P >(
                                           Value( c->estimate ) ) );
     f_fr.add( std::move( c->node ) , cb , f.depth + 1 , BBLeafKind::open );
     }
    if( ! f.any )
     f_fr.add( std::move( f.node ) , f.bound , f.depth , f.kind , false );
    v_stack.pop_back();
    }
   v_path.clear();
   }

  // after the end, the tasks left in the pool and in the answer are open
  // leaves
  void leave( void ) {
   f_sh.pool[ f_id ].drain( [ & ]( Task & t ) {
    Node n = rebuild( t );
    f_fr.add( std::move( n ) , t.bound , t.depth , BBLeafKind::open );
    } );
   for( Task & t : f_sh.answer[ f_id ] ) {
    Node n = rebuild( t );
    f_fr.add( std::move( n ) , t.bound , t.depth , BBLeafKind::open );
    }
   f_sh.answer[ f_id ].clear();
   if( f_a.f_rule == BBRule::ordered && ( f_id == 0 ) )
    for( std::size_t k = std::min( f_sh.next.load() ,
                                   f_sh.ordered.size() ) ;
         k < f_sh.ordered.size() ; ++k ) {
     Task & t = f_sh.ordered[ k ];
     Node n = rebuild( t );
     f_fr.add( std::move( n ) , t.bound , t.depth , BBLeafKind::open );
     }
   }

  // answers a steal request with the next child (or all the children,
  // with steal_all) of the lowest frame that has one, the top excepted
  void answer( void ) {
   std::atomic< unsigned > & a = f_sh.ask[ f_id ];
   unsigned r = a.load( std::memory_order_acquire );
   if( ( r == 0 ) || ( r == c_claimed ) ||
       ! a.compare_exchange_strong( r , c_claimed ,
                                    std::memory_order_acq_rel ) )
    return;
   const unsigned thief = r - 1;
   std::vector< Task > & out = f_sh.answer[ thief ];
   for( std::size_t j = 0 ; ( j + 1 < v_stack.size() ) && out.empty() ;
        ++j )
    while( std::optional< Task > t = take( j ) ) {
     out.push_back( std::move( *t ) );
     if( ! f_a.f_options.steal_all )
      break;
     }
   if( ! out.empty() ) {
    f_sh.pending.fetch_add( std::int64_t( out.size() ) ,
                            std::memory_order_acq_rel );
    f_st.sent += out.size();
    ++f_st.messages;
    }
   f_sh.answered[ thief ].store( true , std::memory_order_release );
   a.store( 0 , std::memory_order_release );
   }

  // answers with nothing a request that comes while the worker has no
  // stack
  void answer_nothing( void ) {
   std::atomic< unsigned > & a = f_sh.ask[ f_id ];
   unsigned r = a.load( std::memory_order_acquire );
   if( ( r == 0 ) || ( r == c_claimed ) ||
       ! a.compare_exchange_strong( r , c_claimed ,
                                    std::memory_order_acq_rel ) )
    return;
   f_sh.answered[ r - 1 ].store( true , std::memory_order_release );
   a.store( 0 , std::memory_order_release );
   }

  // asks the worker v for some of its work, waiting for the answer; true
  // if some came, which is then in the pool
  bool steal_stack( unsigned v ) {
   unsigned zero = 0;
   if( ! f_sh.ask[ v ].compare_exchange_strong( zero , f_id + 1 ,
                                                std::memory_order_acq_rel ) )
    return( false );
   ADS::Backoff backoff;
   while( ! f_sh.answered[ f_id ].load( std::memory_order_acquire ) ) {
    answer_nothing();
    if( f_sh.stop.load() ||
        ( f_sh.pending.load( std::memory_order_acquire ) == 0 ) ) {
     unsigned mine = f_id + 1;
     if( f_sh.ask[ v ].compare_exchange_strong( mine , 0 ,
                                                std::memory_order_acq_rel ) )
      return( false );                 // withdrawn before an answer
     }
    backoff.pause();
    }
   f_sh.answered[ f_id ].store( false , std::memory_order_relaxed );
   std::vector< Task > & in = f_sh.answer[ f_id ];
   const bool got = ! in.empty();
   for( Task & t : in )
    f_sh.pool[ f_id ].push( std::move( t ) );
   in.clear();
   return( got );
   }

  // the next task: an ordered one, one of the worker, one of the others,
  // or one stolen from a stack; none when no task is left or a limit is
  // reached
  std::optional< Task > next( void ) {
   const unsigned p = f_sh.p;
   if( f_a.f_rule == BBRule::ordered ) {
    if( stop() )
     return( std::nullopt );
    const std::size_t k = f_sh.next.fetch_add( 1 );
    if( k >= f_sh.ordered.size() )
     return( std::nullopt );
    return( std::move( f_sh.ordered[ k ] ) );
    }
   ADS::Backoff backoff;
   for( ;; ) {
    answer_nothing();
    if( std::optional< Task > t = f_sh.pool[ f_id ].pop( true ) )
     return( t );
    if( f_sh.stop.load() ||
        ( f_sh.pending.load( std::memory_order_acquire ) == 0 ) )
     return( std::nullopt );
    if( p > 1 ) {
     const unsigned v = unsigned( ( f_id + 1 +
                                    random_below( f_gen , p - 1 ) ) % p );
     if( std::optional< Task > t = f_sh.pool[ v ].pop( false ) ) {
      ++f_st.sent;
      ++f_st.messages;
      return( t );
      }
     if( ( f_a.f_rule == BBRule::stack_stealing ) && steal_stack( v ) )
      continue;
     }
    backoff.pause();
    }
   }

  // the ordered tasks: the nodes of depth spawn_depth under the starts
  // (a start of depth spawn_depth or more being a task itself, of priority
  // 0), the nodes above them evaluated and branched, sorted by their
  // discrepancy (or left to right), ties in the order of generation; it
  // runs on the worker 0 alone, before the workers start, and does not
  // check the node and the time limits, so that all the nodes above
  // spawn_depth are evaluated even past a limit (they are counted in
  // evaluated, and the limit stops the workers before their first task).
  // The tasks are not counted in pending, which the ordered rule does not
  // read: NEXT ends when the tasks are used up, and those not taken when a
  // limit stops the workers are made open leaves by leave() of the worker 0
  void order( void ) {
   const auto & o = f_a.f_options;
   std::vector< Task > all;               // the tasks
   std::vector< std::uint64_t > prio;     // ... and their priorities
   for( std::size_t s = 0 ; s < f_sh.starts.size() ; ++s ) {
    f_start = std::uint32_t( s );
    v_base.clear();
    v_path.clear();
    v_stack.clear();
    const Leaf & l = f_sh.starts[ s ];
    if( l.depth >= o.spawn_depth ) {
     Task t{ std::nullopt , l.bound , l.depth , std::uint32_t( s ) , {} };
     if( ! o.indexed )
      t.node.emplace( l.node );
     all.push_back( std::move( t ) );
     prio.push_back( 0 );
     continue;
     }
    std::vector< std::uint64_t > disc;   // the discrepancy of each frame
    visit( Node( l.node ) , l.bound , l.depth );
    if( ! v_stack.empty() )
     disc.push_back( 0 );
    while( ! v_stack.empty() ) {
     Frame & f = v_stack.back();
     std::optional< BBChild< Node , Value > > c;
     if( ! f.cut )
      c = f.cur->next();
     if( ! c ) {
      if( ! f.any )
       f_fr.add( std::move( f.node ) , f.bound , f.depth , f.kind , false );
      else
       if( f.cut )
        drain_pruned( f );
      v_stack.pop_back();
      disc.pop_back();
      if( ! v_path.empty() )
       v_path.pop_back();
      continue;
      }
     f.any = true;
     const std::uint32_t k = f.k++;
     const unsigned d = f.depth + 1;
     const std::optional< Value > cb = detail::bb_child_bound< P >(
      f.bound , Value( c->estimate ) , inc() , o , f.cut , f_st , f_fr );
     if( ! cb ) {
      if( f_fr.kept( BBLeafKind::pruned ) )
       f_fr.add( std::move( c->node ) ,
                 std::max( f.bound , detail::bb_signed< P >(
                                      Value( c->estimate ) ) ) ,
                 d , BBLeafKind::pruned );
      continue;
      }
     const std::uint64_t dk = disc.back() + ( o.discrepancy ? k : 0 );
     if( d >= o.spawn_depth ) {
      all.push_back( task_of( std::move( c->node ) , *cb , d ,
                              v_stack.size() - 1 , k ) );
      prio.push_back( dk );
      continue;
      }
     const std::size_t before = v_stack.size();
     v_path.push_back( k );
     visit( std::move( c->node ) , *cb , d );
     if( v_stack.size() == before )
      v_path.pop_back();
     else
      disc.push_back( dk );
     }
    }
   std::vector< std::size_t > by( all.size() );
   std::iota( by.begin() , by.end() , std::size_t( 0 ) );
   std::stable_sort( by.begin() , by.end() ,
                     [ & ]( std::size_t a , std::size_t b ) {
                      return( prio[ a ] < prio[ b ] );
                      } );
   f_sh.ordered.reserve( all.size() );
   for( const std::size_t j : by )
    f_sh.ordered.push_back( std::move( all[ j ] ) );
   }

  const SkeletonBranchAndBound & f_a;           // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  const Stopwatch & f_clock;                    // the clock of the search
  W f_w;                                        // its evaluator
  std::mt19937_64 f_gen;                        // the random victims
  std::deque< Frame > v_stack;                  // the depth-first stack
  std::uint32_t f_start = 0;                    // the start of the task
  std::vector< std::uint32_t > v_base;          // the path to the task
  std::vector< std::uint32_t > v_path;          // ... and from it
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 BBRule f_rule;                                       ///< the rule
 Options f_options;                                   ///< the options
 Executor * f_ex = nullptr;                           ///< given, or nullptr

 };  // end( class( SkeletonBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SkeletonBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------- End File SkeletonBranchAndBound.h --------------------*/
/*--------------------------------------------------------------------------*/
