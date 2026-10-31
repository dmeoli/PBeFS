/*--------------------------------------------------------------------------*/
/*----------------------- File FarmBranchAndBound.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for FarmBranchAndBound, the master-worker Branch-and-Bound
 * (B&B) on a farm with feedback on an executor of ADS: the master holds the
 * queue of the open nodes, and the workers evaluate and branch the nodes it
 * sends. The problem, the options and the Solution are those of
 * BranchAndBound.h and ParallelBranchAndBound.h, and the design is discussed
 * in docs/branch-and-bound.md.
 *
 * The search is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation, the procedures of the sequential
 * search [see BranchAndBound.h] and DIVE of the parallel ones [see
 * ParallelBranchAndBound.h]; the master alone holds the queue Q, out counts
 * the tasks at the workers, and the incumbent \f$ z \f$ is shared:
 *
 * \code{.unparsed}
 *  1  procedure MASTER(R, p)
 *  2     Q ← R;  out ← 0
 *  3     loop
 *  4        DISPATCH()
 *  5        if out = 0 and (a limit is reached or Q = ∅) then return
 *  6        receive a task T back from a worker;  out ← out − 1
 *  7        push the children of T in Q
 *  8  procedure DISPATCH()
 *  9     while no limit is reached, out < p and Q ≠ ∅ do
 * 10        if PRUNABLE(the least bound of Q) then
 * 11           the nodes of Q are leaves (pruned);  Q ← ∅;  return
 * 12        T ← the first batch nodes of Q, removed
 * 13        send T to a free worker;  out ← out + 1
 * 14  procedure WORKER(T)                       ▷ on each of the p workers
 * 15     for each n ∈ T do
 * 16        while n ≠ nil do                  ▷ the node, then its dive
 * 17           K ← PROCESS(n);  n ← DIVE(K), unless the node limit is reached
 * 18           add the nodes of K but n to the children of T
 * 19     send T back to the master
 * \endcode
 *
 * The test of line 5 is exact, since the master alone changes Q and every
 * open node is either in Q or in a task that is out; the order of extraction
 * is the one of the selection, and the master only moves nodes, which suits
 * an evaluation that is "expensive" w.r.t. the transfer of a node.
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

#ifndef __PBeFS_FarmBranchAndBound
 #define __PBeFS_FarmBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>

#include <memory>

#include <stdexcept>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

#include "ADS/parallel/Farm.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS FarmBranchAndBound ------------------------*/
/*--------------------------------------------------------------------------*/
/// the master-worker B&B on a farm with feedback
/** The master extracts the nodes from its queue by the Selection and sends
 * them, batch at a time, to the workers that are free; a worker evaluates
 * and branches each node, offering its feasible value to the shared
 * incumbent, and sends back the children that survive the pruning (with
 * depth-first or plunging, it goes on with the child that dives in the same
 * task). The search ends when the queue is empty, or no node of it can
 * improve the incumbent, and no task is out at a worker. Since the master
 * only moves nodes, the order of extraction is the one of the selection, and
 * the scheme suits an evaluation that is expensive w.r.t. the transfer of a
 * node.
 *
 * The farm runs threads + 1 threads, those of an Executor of ADS [see
 * ADS::feedbackFarm()]: with PBEFS_PIN the master takes the first CPU of the
 * affinity mask and the worker \f$ i \f$ the CPU \f$ i + 1 \f$ of it,
 * counting from 0 [see pin_worker()], but for the thread of the caller. */

template< BranchAndBound P , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class FarmBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using Leaf = BBLeaf< Node , Value >;           ///< the leaf type
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search);
 /// threads is the number of workers, batch the nodes per task

 explicit FarmBranchAndBound( const P & problem ,
                              const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the master and the workers, which must
 /// outlive the searches and have at least threads + 1 threads; without
 /// it, each search builds one of threads + 1 threads

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

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

 using Item = typename BBQueue< Node , Value >::Item;

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  Stopwatch clock;
  const Value start = f_options.incumbent ?
                      detail::bb_signed< P >( *f_options.incumbent ) : inf;
  Shared sh( start , clock );
  Master master( *this , sh , detail::bb_starts( f_p , std::move( seeds ) ) );
  const unsigned n = f_options.threads;
  std::vector< std::unique_ptr< Worker > > w( n );
  for( unsigned i = 0 ; i < n ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );
  const auto begin = [ & ]( auto && emit , std::size_t ) {
   if( ! ADS::isRunCaller() )
    pin_worker( 0 );
   master.start( emit );
   };
  const auto back = [ & ]( Task * t , auto && emit ) {
   return( master.back( t , emit ) );
   };
  const auto work = [ & ]( Task * t , std::size_t i ) {
   return( w[ i ]->process( t ) );
   };
  if( f_ex )
   ADS::feedbackFarm< Task * , Task * >( *f_ex , n , begin , back , work );
  else {
   Executor ex( n + 1 );
   ADS::feedbackFarm< Task * , Task * >( ex , n , begin , back , work );
   }

  SolutionT sol;
  // left open by a limit
  detail::BBFrontier< Node , Value > & fr = master.f_fr;
  fr.add( master.f_q , BBLeafKind::open );
  for( auto & x : w ) {
   sol.stats.add( x->f_st );
   fr.merge( std::move( x->f_fr ) );
   }
  sol.stats.pruned_extracted += master.f_st.pruned_extracted;
  sol.stats.max_open = master.f_st.max_open;
  sol.stats.messages = master.f_st.messages;
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

 // a task: the nodes to evaluate, and then the children that survive
 struct Task {
  std::vector< Item > nodes;                    // the nodes
  std::vector< Item > kids;                     // their children
  Value best = Value( 0 );                      // the best open bound
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the master and the workers

 struct Shared {

  Shared( Value z , const Stopwatch & c ) : inc( z ) , clock( c ) {}

  BBIncumbent< Node , Value > inc;              // the incumbent
  const Stopwatch & clock;                      // the clock of the search
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the master: the queue of the open nodes

 struct Master {

  Master( const FarmBranchAndBound & a , Shared & s ,
          std::vector< Leaf > && starts )
   : f_a( a ) , f_sh( s ) , f_q( a.f_options.selection ) ,
     f_fr( a.f_options.keep_frontier ) , v_starts( std::move( starts ) ) {}

  Master( const Master & ) = delete;
  Master & operator=( const Master & ) = delete;

  ~Master() {
   for( Task * t : v_free )
    delete t;
   }

  // the start: the first nodes, and the first tasks
  template< class E >
  void start( E && emit ) {
   for( auto & l : v_starts )
    enqueue( Item{ std::move( l.node ) , l.bound , l.depth } );
   v_starts.clear();
   check_limits();
   dispatch( emit );
   }

  // a task back, with the children of its nodes; false once nothing is
  // out at the workers and nothing is left to send
  template< class E >
  bool back( Task * t , E && emit ) {
   --f_out;
   for( Item & it : t->kids )
    enqueue( std::move( it ) );
   t->kids.clear();
   v_free.push_back( t );
   check_limits();
   dispatch( emit );
   return( ( f_out > 0 ) || ( ! f_sh.limit.load() && ! f_q.empty() ) );
   }

  // pushes a node in the queue
  void enqueue( Item && it ) {
   const Value key = detail::bb_key( f_a.f_p , it.node , it.bound ,
                                      it.depth , f_a.f_options.selection );
   f_q.push( std::move( it.node ) , it.bound , it.depth , key );
   f_st.max_open = std::max( f_st.max_open , f_q.size() );
   }

  // stops the dispatch if a limit is reached
  void check_limits( void ) {
   ( void ) detail::bb_limit( f_a.f_options , f_sh.evaluated.load() ,
                              f_sh.clock , f_sh.limit , f_sh.time );
   }

  // sends tasks while there are nodes and free workers
  template< class E >
  void dispatch( E && emit ) {
   while( ! f_sh.limit.load() && ( f_out < f_a.f_options.threads ) &&
          ! f_q.empty() ) {
    const Value best = f_q.best_bound();
    if( detail::bb_prunable( best , f_sh.inc.value() , f_a.f_options ) ) {
     f_st.pruned_extracted += f_fr.add( f_q , BBLeafKind::pruned );
     return;                                    // nothing can improve
     }
    Task * t;
    if( v_free.empty() )
     t = new Task;
    else {
     t = v_free.back();
     v_free.pop_back();
     }
    t->nodes.clear();
    while( ! f_q.empty() && ( t->nodes.size() < f_a.f_options.batch ) )
     t->nodes.push_back( f_q.pop() );
    t->best = std::min( best , f_q.best_bound() );
    ++f_out;
    ++f_st.messages;
    emit( t );
    }
   }

  const FarmBranchAndBound & f_a;               // the search
  Shared & f_sh;                                // the shared data
  BBQueue< Node , Value > f_q;                  // the open nodes
  std::vector< Task * > v_free;                 // tasks to reuse
  unsigned f_out = 0;                           // tasks at the workers
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  std::vector< Leaf > v_starts;       // the nodes to start from
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: evaluates and branches the nodes of a task

 struct Worker {

  Worker( const FarmBranchAndBound & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) ,
     f_fr( a.f_options.keep_frontier ) {}

  // the nodes of t evaluated and branched, the children in t
  Task * process( Task * t ) {
   if( ! f_w ) {                     // the first task, on its own thread
    if( ! ADS::isRunCaller() )
     pin_worker( f_id + 1 );
    f_w.emplace( detail::bb_worker( f_a.f_p , f_id ) );
    }
   const auto & o = f_a.f_options;
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };
   for( Item & first : t->nodes ) {
    std::optional< Item > it( std::move( first ) );
    while( it ) {                  // the node, then the child that dives
     v_kids.clear();
     detail::bb_process< P >( it->node , it->bound , it->depth , *f_w , inc ,
                              improve , [ & ]( Node && c , Value b ,
                                               unsigned d ) {
                               v_kids.push_back( Item{ std::move( c ) , b ,
                                                       d } );
                               } , o , f_st , f_fr );
     f_sh.evaluated.fetch_add( 1 );
     it.reset();
     const std::size_t k = detail::bb_dive( v_kids , t->best , inc() , o );
     for( std::size_t j = 0 ; j < v_kids.size() ; ++j )
      if( j == k )
       it.emplace( std::move( v_kids[ j ] ) );
      else
       t->kids.push_back( std::move( v_kids[ j ] ) );
     if( it && o.node_limit &&
         ( f_sh.evaluated.load() >= o.node_limit ) ) {
      t->kids.push_back( std::move( *it ) );     // the dive stops here
      it.reset();
      }
     }
    }
   t->nodes.clear();
   return( t );
   }

  using Evaluator = decltype( detail::bb_worker( std::declval< const P & >() ,
                                                 0u ) );

  const FarmBranchAndBound & f_a;               // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  std::optional< Evaluator > f_w;               // its evaluator
  std::vector< Item > v_kids;                   // the children of a node
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options
 Executor * f_ex = nullptr;                   ///< given, or nullptr

 };  // end( class( FarmBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_FarmBranchAndBound */

/*--------------------------------------------------------------------------*/
/*--------------------- End File FarmBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
