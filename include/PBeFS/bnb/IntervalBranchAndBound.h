/*--------------------------------------------------------------------------*/
/*--------------------- File IntervalBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for IntervalBranchAndBound, the parallel Branch-and-Bound
 * (B&B) whose work is an interval of the leaves of the tree, numbered in a
 * mixed radix, and which an idle worker shares by taking the second half of
 * the interval of another one, as the Integer-Vector-Matrix (IVM) of Mezmaz,
 * Leroy, Melab and Tuyttens (2014) and of Gmys et al. (2016) for the
 * permutation problems. The problem, the options and the Solution are those
 * of BranchAndBound.h and ParallelBranchAndBound.h, and the design is
 * discussed in docs/branch-and-bound.md.
 *
 * The tree must have the same number \f$ r_d \f$ of children at all the
 * nodes of depth \f$ d \f$ (\f$ n - d \f$ for the permutations of \f$ n \f$
 * elements), given by arity( d ) of the problem, 0 at the depth \f$ D \f$
 * of the leaves, and the evaluator must give the child \f$ k \f$ of a node
 * without its evaluation [see IndexedBBEvaluator]. A leaf is then the
 * vector \f$ ( x_0 , \ldots , x_{D-1} ) \f$ of the indices of the children
 * on its path, \f$ 0 \le x_d < r_d \f$, and the leaves are the integers of
 * \f$ [ 0 , \prod_d r_d ) \f$ in that mixed radix, the leaves of a node of
 * depth \f$ d \f$ making the interval of the numbers that begin with its
 * path and go on with \f$ D - d \f$ digits. A worker holds the interval
 * \f$ [ a , b ) \f$ of the leaves it has still to explore and the nodes of
 * the path to \f$ a \f$ (the vector and the matrix of the IVM), and the
 * search is the following, with the procedures of the sequential search
 * [see BranchAndBound.h]:
 *
 * \code{.unparsed}
 *  1  procedure INTERVAL()
 *  2     the worker 0 holds [0, ∏ r_d), the others nothing
 *  3     in parallel, on each worker i = 0, ..., p − 1 do
 *  4        loop
 *  5           while a_i < b_i do
 *  6              n ← the node of the path to a_i at the current depth d
 *  7              if d = D or PRUNABLE(β(n)) after EVALUATE(n) then
 *  8                 a_i ← the first leaf after the subtree of n
 *  9              else go down to the child a_i[d] of n
 * 10           if no worker holds leaves then return
 * 11           v ← a random worker;  with the lock of v, if b_v − a_v ≥ 2
 *                 then m ← a_v + ⌊(b_v − a_v)/2⌋;  e ← b_v;  b_v ← m;
 *                 then, with the lock of i alone, [a_i, b_i) ← [m, e)
 * \endcode
 *
 * The evaluation of a node of line 7 counts as such only if its interval
 * begins at \f$ a_i \f$ (the first leaf of its subtree is the next of the
 * worker); otherwise the node is an ancestor of the leaf at which the
 * worker took its interval, which the worker it came from evaluates too,
 * and it counts among the recomputed nodes. Hence, each node is evaluated
 * once, by the worker whose interval holds its first leaf, and with the
 * incumbent fixed (e.g., given as the optimum) the nodes evaluated are those
 * of the sequential depth-first search, whatever the steals. A steal moves
 * two vectors of \f$ D \f$ digits, whatever the node, and the thief builds
 * the path by child( n , k ). After a limit, the leaves left in the
 * intervals are covered by the fewest nodes, which are open leaves of the
 * Solution, with the estimate of their parent as their bound. The search
 * starts from the root only.
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

#ifndef __PBeFS_IntervalBranchAndBound
 #define __PBeFS_IntervalBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstdint>

#include <exception>

#include <memory>

#include <mutex>

#include <optional>

#include <random>

#include <vector>

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a B&B problem whose nodes of a depth have the same number of children,
/// given by arity( d ), and whose evaluators give the child k of a node

template< class P >
concept RadixBranchAndBound = BranchAndBound< P > &&
 requires( const P & p , unsigned d ) {
  { p.arity( d ) } -> std::convertible_to< std::size_t >;
  } &&
 IndexedBBEvaluator< decltype( detail::bb_worker( std::declval< const P & >() ,
                                                  0u ) ) ,
                     typename P::node_type , bb_value_t< P > >;

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS IntervalBranchAndBound ----------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel B&B by the work stealing of intervals of leaves (IVM)
/** Each worker explores depth-first the leaves of its interval, and an
 * idle one takes the second half of the interval of another worker, drawn
 * at random [see IntervalBranchAndBound.h]. The incumbent is shared; the
 * counter sent of the statistics holds the steals, and recomputed the
 * evaluations of the ancestors of the leaves at which the intervals were
 * taken. The selection of the options is not used. */

template< RadixBranchAndBound P ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
class IntervalBranchAndBound {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using Node = typename P::node_type;            ///< the node type
 using Value = bb_value_t< P >;                 ///< the value type
 using Options = BasicParallelBBOptions< Value >; ///< the options
 using SolutionT = BBSolution< Node , Value >;  ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the problem (which must outlive the search)

 explicit IntervalBranchAndBound( const P & problem ,
                                  const Options & options = {} )
  : f_p( problem ) , f_options( options ) {
  f_options.normalize();
  for( unsigned d = 0 ; problem.arity( d ) > 0 ; ++d )
   v_radix.push_back( std::uint64_t( problem.arity( d ) ) );
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const Value inf = BBValueTraits< Value >::infinity();
  const Options & o = f_options;
  const unsigned p = o.threads;
  Stopwatch clock;
  Shared sh( p , o.incumbent ? detail::bb_signed< P >( *o.incumbent ) :
                               inf );
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i , clock );
  w[ 0 ]->f_a.assign( v_radix.size() , 0 );   // the whole tree
  w[ 0 ]->f_b.assign( v_radix.size() , 0 );
  if( v_radix.empty() )
   w[ 0 ]->f_b.push_back( 1 );                 // the root is the leaf
  else
   w[ 0 ]->f_b[ 0 ] = v_radix[ 0 ];
  w[ 0 ]->f_held = true;
  sh.active.store( 1 );

  std::exception_ptr error;
  std::mutex m;
  run_workers( f_ex , p , [ & ]( unsigned i ) {
   try {
    w[ i ]->run();
    }
   catch( ... ) {
    std::lock_guard< std::mutex > guard( m );
    if( ! error )
     error = std::current_exception();
    sh.stop.store( true );
    }
   } );
  if( error )
   std::rethrow_exception( error );

  SolutionT sol;
  detail::BBFrontier< Node , Value > fr( o );
  for( auto & x : w ) {
   x->leave();
   sol.stats.add( x->f_st );
   fr.merge( std::move( x->f_fr ) );
   }
  detail::bb_finish< P >( sol , sh.inc , fr , sh.limit.load() ,
                          sh.time.load() , o );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using W = decltype( detail::bb_worker( std::declval< const P & >() , 0u ) );
 using Digits = std::vector< std::uint64_t >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 struct Worker;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {
  Shared( unsigned np , Value z ) : p( np ) , inc( z ) {}
  unsigned p;                                   // the number of workers
  BBIncumbent< Node , Value > inc;              // the incumbent
  alignas( cache_line ) std::atomic< unsigned > active{ 0 };  // holding
  alignas( cache_line ) std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< bool > stop{ false };            // a limit was reached
  std::atomic< bool > limit{ false };           // stopped by a limit
  std::atomic< bool > time{ false };            // ... the time limit
  std::vector< Worker * > workers;              // all of them
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its interval, its path, its counters and its leaves

 struct alignas( cache_line ) Worker {

  Worker( const IntervalBranchAndBound & s , Shared & sh , unsigned i ,
          const Stopwatch & c )
   : f_s( s ) , f_sh( sh ) , f_id( i ) , f_clock( c ) ,
     f_w( detail::bb_worker( s.f_p , i ) ) , f_gen( s.f_options.seed + i ) ,
     f_fr( s.f_options ) {
   if( sh.workers.size() < sh.p )
    sh.workers.resize( sh.p , nullptr );
   sh.workers[ i ] = this;
   }

  const Digits & radix( void ) const { return( f_s.v_radix ); }

  // a < b, lexicographically
  static bool less( const Digits & a , const Digits & b ) {
   return( std::lexicographical_compare( a.begin() , a.end() , b.begin() ,
                                         b.end() ) );
   }

  // a + b, and a - b (a >= b), in the mixed radix, digit 0 unbounded
  Digits add( const Digits & a , const Digits & b ) const {
   Digits c( a.size() );
   std::uint64_t carry = 0;
   for( std::size_t d = a.size() ; d-- > 0 ; ) {
    std::uint64_t x = a[ d ] + b[ d ] + carry;
    carry = 0;
    if( ( d > 0 ) && ( x >= radix()[ d ] ) ) {
     x -= radix()[ d ];
     carry = 1;
     }
    c[ d ] = x;
    }
   return( c );
   }

  Digits sub( const Digits & a , const Digits & b ) const {
   Digits c( a.size() );
   std::uint64_t borrow = 0;
   for( std::size_t d = a.size() ; d-- > 0 ; ) {
    const std::uint64_t y = b[ d ] + borrow;
    if( a[ d ] >= y ) {
     c[ d ] = a[ d ] - y;
     borrow = 0;
     }
    else {
     c[ d ] = a[ d ] + radix()[ d ] - y;
     borrow = 1;
     }
    }
   return( c );
   }

  // a / 2, and whether it is 0
  Digits half( const Digits & a , bool & zero ) const {
   Digits c( a.size() );
   std::uint64_t rem = 0;
   zero = true;
   for( std::size_t d = 0 ; d < a.size() ; ++d ) {
    const std::uint64_t x = a[ d ] + rem * ( d > 0 ? radix()[ d ] : 0 );
    c[ d ] = x / 2;
    rem = x % 2;
    zero = zero && ( c[ d ] == 0 );
    }
   return( c );
   }

  // true if a limit is reached, which stops all the workers
  bool stop( void ) {
   if( f_sh.stop.load( std::memory_order_relaxed ) )
    return( true );
   if( ! detail::bb_limit( f_s.f_options , f_sh.evaluated.load() , f_clock ,
                           f_sh.limit , f_sh.time ) )
    return( false );
   f_sh.stop.store( true );
   return( true );
   }

  // true if the interval holds no leaf, under the lock
  bool empty( void ) {
   std::lock_guard< std::mutex > guard( f_m );
   return( ! less( f_a , f_b ) );
   }

  // the main loop
  void run( void ) {
   ADS::Backoff backoff;
   for( ;; ) {
    if( f_held ) {
     explore();
     f_held = false;
     f_sh.active.fetch_sub( 1 , std::memory_order_acq_rel );
     }
    if( f_sh.stop.load() ||
        ( f_sh.active.load( std::memory_order_acquire ) == 0 ) )
     return;
    if( steal() )
     backoff.reset();
    else
     backoff.pause();
    }
   }

  // takes the second half of the interval of a random worker
  bool steal( void ) {
   const unsigned p = f_sh.p;
   if( p < 2 )
    return( false );
   const unsigned v = unsigned( ( f_id + 1 +
                                  random_below( f_gen , p - 1 ) ) % p );
   Worker & x = *f_sh.workers[ v ];
   Digits mid , end;
   {
    std::lock_guard< std::mutex > guard( x.f_m );
    if( ! less( x.f_a , x.f_b ) )
     return( false );
    bool zero;
    const Digits h = half( sub( x.f_b , x.f_a ) , zero );
    if( zero )
     return( false );
    mid = add( x.f_a , h );
    end = x.f_b;
    x.f_b = mid;
    // counted while the victim, which holds leaves, is counted too
    f_sh.active.fetch_add( 1 , std::memory_order_acq_rel );
    }
   {                                   // the thief's interval is empty here
    std::lock_guard< std::mutex > mine( f_m );
    f_a = std::move( mid );
    f_b = std::move( end );
    }
   f_held = true;
   ++f_st.sent;
   ++f_st.messages;
   return( true );
   }

  // the depth-first search of the interval
  void explore( void ) {
   const auto & o = f_s.f_options;
   const std::size_t D = radix().size();
   v_node.clear();
   v_bound.clear();
   v_node.push_back( f_s.f_p.root() );
   v_bound.push_back( - BBValueTraits< Value >::infinity() );
   std::size_t d = 0;                  // the depth of the current node
   const auto inc = [ & ] { return( f_sh.inc.value() ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    return( f_sh.inc.template improve< P >( v , n , o ) );
    };
   for( ;; ) {
    if( stop() )
     return;
    Digits a;
    {
     std::lock_guard< std::mutex > guard( f_m );
     if( ! less( f_a , f_b ) )
      return;
     a = f_a;
     }
    // the node at depth d is on the path to a: its own if a has zeros
    // from d on, an ancestor of the start of the interval otherwise
    bool own = true;
    for( std::size_t j = d ; j < D ; ++j )
     own = own && ( a[ j ] == 0 );
    Node & n = v_node[ d ];
    const std::uint64_t e0 = f_st.evaluated;
    BBLeafKind kind = BBLeafKind::pruned;
    BBStats tmp;
    BBStats & st = own ? f_st : tmp;
    const std::optional< Value > b = detail::bb_evaluate< P >( n ,
                                     v_bound[ d ] , unsigned( d ) , f_w ,
                                     inc , improve , o , st , f_fr , kind );
    if( ! own )
     f_st.recomputed += tmp.evaluated;
    else
     if( f_st.evaluated != e0 )
      f_sh.evaluated.fetch_add( 1 , std::memory_order_relaxed );
    if( ! b || ( d == D ) ) {
     if( ! skip( d ) )
      return;
     d = sibling();
     continue;
     }
    v_bound[ d ] = *b;
    if( down( d , a[ d ] ) )            // down to the child a[ d ]
     ++d;
    else {
     if( ! skip( d + 1 ) )
      return;
     d = sibling();
     }
    }
   }

  // moves a past the subtree of the node at depth d of its path; false if
  // the whole tree is done
  bool skip( std::size_t d ) {
   std::lock_guard< std::mutex > guard( f_m );
   const std::size_t D = radix().size();
   for( std::size_t j = d ; j < D ; ++j )
    f_a[ j ] = 0;
   if( d == 0 ) {                       // the root: nothing is left
    f_a = f_b;
    return( false );
    }
   std::size_t j = d - 1;
   ++f_a[ j ];
   while( ( j > 0 ) && ( f_a[ j ] == radix()[ j ] ) ) {
    f_a[ j ] = 0;
    ++f_a[ --j ];
    }
   f_skipped = j;                       // the depth of the last change
   return( true );
   }

  // after skip(), the depth of the node to visit next, whose path is the
  // new a: the nodes above the last digit changed are kept, and the
  // children down from there are built, skipping those pruned by their
  // estimate; 0 if the interval is done
  std::size_t sibling( void ) {
   for( ;; ) {
    const std::size_t j = f_skipped;
    std::uint64_t k;
    {
     std::lock_guard< std::mutex > guard( f_m );
     if( ! less( f_a , f_b ) )
      return( 0 );
     k = f_a[ j ];
     }
    v_node.resize( j + 1 );
    v_bound.resize( j + 1 );
    if( down( j , k ) )
     return( j + 1 );
    if( ! skip( j + 1 ) )
     return( 0 );
    }
   }

  // builds the child k of the node at depth d as the node at depth d + 1;
  // false if it is pruned by its estimate (and then it is a leaf)
  bool down( std::size_t d , std::uint64_t k ) {
   const auto & o = f_s.f_options;
   auto c = f_w.child( std::as_const( v_node[ d ] ) , std::size_t( k ) );
   bool cut = false;
   const std::optional< Value > cb = detail::bb_child_bound< P >(
    v_bound[ d ] , Value( c->estimate ) , f_sh.inc.value() , o , cut ,
    f_st , f_fr );
   if( ! cb ) {
    if( f_fr.kept( BBLeafKind::pruned ) )
     f_fr.add( std::move( c->node ) ,
               std::max( v_bound[ d ] , detail::bb_signed< P >(
                                         Value( c->estimate ) ) ) ,
               unsigned( d + 1 ) , BBLeafKind::pruned );
    return( false );
    }
   v_node.resize( d + 2 );
   v_bound.resize( d + 2 );
   v_node[ d + 1 ] = std::move( c->node );
   v_bound[ d + 1 ] = *cb;
   return( true );
   }

  // after a limit, the leaves left in the interval as the fewest nodes
  void leave( void ) {
   const std::size_t D = radix().size();
   while( less( f_a , f_b ) ) {
    // the shallowest node whose leaves begin at a and end before b
    std::size_t d = D;
    while( ( d > 0 ) && ( f_a[ d - 1 ] == 0 ) ) {
     Digits last = f_a;
     for( std::size_t j = d - 1 ; j < D ; ++j )
      last[ j ] = radix()[ j ] - 1;
     if( ! less( last , f_b ) )
      break;
     --d;
     }
    Node n = f_s.f_p.root();
    Value bound = - BBValueTraits< Value >::infinity();
    for( std::size_t j = 0 ; j < d ; ++j ) {
     auto c = f_w.child( std::as_const( n ) , std::size_t( f_a[ j ] ) );
     bound = std::max( bound , detail::bb_signed< P >(
                                Value( c->estimate ) ) );
     n = std::move( c->node );
     }
    f_fr.add( std::move( n ) , bound , unsigned( d ) , BBLeafKind::open );
    if( ! skip( d ) )
     break;
    }
   }

  const IntervalBranchAndBound & f_s;           // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  const Stopwatch & f_clock;                    // the clock of the search
  W f_w;                                        // its evaluator
  std::mt19937_64 f_gen;                        // the random victims
  std::mutex f_m;                               // protects the interval
  Digits f_a;                                   // the next leaf
  Digits f_b;                                   // the end of the interval
  bool f_held = false;                          // true if it has leaves
  std::size_t f_skipped = 0;                    // the digit last changed
  std::vector< Node > v_node;                   // the nodes of the path
  std::vector< Value > v_bound;                 // their bounds
  BBStats f_st;                                 // its counters
  detail::BBFrontier< Node , Value > f_fr;      // its leaves
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                                   ///< the options
 Digits v_radix;                                      ///< the arities
 Executor * f_ex = nullptr;                           ///< given, or nullptr

 };  // end( class( IntervalBranchAndBound ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_IntervalBranchAndBound */

/*--------------------------------------------------------------------------*/
/*------------------- End File IntervalBranchAndBound.h --------------------*/
/*--------------------------------------------------------------------------*/
