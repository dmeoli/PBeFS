/*--------------------------------------------------------------------------*/
/*------------------------------ File EPASE.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template EPASE, the Edge-based Parallel A* for
 * Slow Evaluations of Mukherjee, Aine and Likhachev (SoCS 2022) with its
 * weighted variant w-ePA*SE, the generalized GePA*SE of Mukherjee and
 * Likhachev (SoCS 2023), and their anytime variants A-ePA*SE and A-GePA*SE
 * (Yang, Mukherjee and Likhachev, SoCS 2023).
 *
 * The search works on the edges of an EdgeDomain [see Domain.h] instead of
 * its states. OPEN holds the states whose g has improved, each standing for
 * all its edges (the dummy edge of the paper); READY holds the edges of the
 * states being expanded (the set BE) still to be evaluated. A state, an
 * edge and the states of BE all have the key \f$ f' = g + w h \f$ of their
 * state, \f$ g \f$ being, for BE, the one at the expansion; \f$ h( t , s ) \f$
 * is the pairwise heuristic of the Domain, \f$ \varepsilon \geq w \f$ the
 * inflation of the independence rule and \f$ K \f$ the states of OPEN
 * examined (all, if 0). One lock \f$ L \f$ protects all the sets:
 *
 * \code{.unparsed}
 *  1  procedure SELECT                              ▷ under L
 *  2     r ← the first edge of READY, if any
 *  3     for each s among the first K states of OPEN with key(s) < key(r) do
 *  4        if SAFE(s) and (s is not a goal or ENDS(s)) then return s
 *  5     return r, if any
 *  6  procedure SAFE(s)
 *  7     for each t ∈ BE ∪ OPEN with f'(t) < f'(s) do
 *  8        if g(s) − g(t) > ε h(t, s) then return false
 *  9     return true
 * 10  procedure EXPAND(s)                           ▷ under L
 * 11     move s from OPEN to BE
 * 12     for each action a of s do
 * 13        if a is cheap (GePA*SE) then add a to the batch B of s
 * 14        else add the edge (s, a) to READY
 * 15     if B ≠ ∅ then add B to READY
 * 16  procedure EVALUATE(s, A)          ▷ the edges A of s, outside of L
 * 17     S ← the successors (v, c(s, v)) of the edges of A
 * 18     lock L
 * 19     for each (v, c(s, v)) ∈ S do
 * 20        g' ← g(s) + c(s, v)
 * 21        if v is new then add v to OPEN, with g(v) ← g', parent(v) ← s
 * 22        else if v ∈ OPEN and g' < g(v) then
 * 23           g(v) ← g';  parent(v) ← s;  reorder v in OPEN
 * 24     if all the edges of s are evaluated then move s from BE to CLOSED
 * 25     unlock L
 * 26  procedure ePA*SE(s₀, p)
 * 27     OPEN ← {s₀}, with g(s₀) ← 0;  READY ← ∅;  BE ← ∅
 * 28     on the worker 0 (the dispatcher), lock L
 * 29     while not done do
 * 30        x ← SELECT, an edge only if a worker is idle
 * 31        if there is none then
 * 32           if OPEN = READY = BE = ∅ then done ← true
 * 33           else wait on L until a set changes or a worker is idle
 * 34        else if x is a state and a goal then record x;  done ← true
 * 35        else if x is a state then EXPAND(x)
 * 36        else remove x from READY;  give it to the idle worker of
 * 37             least index, which runs EVALUATE on it and waits again
 * \endcode
 *
 * where ENDS(s) holds if no state of OPEN or BE has a key smaller than
 * \f$ g( s ) \f$. A state \f$ s \f$ is independent of \f$ t \f$ when
 * \f$ g( s ) - g( t ) \leq \varepsilon \, h( t , s ) \f$. A state of OPEN is
 * expanded only when it is safe (line 4), and an edge only when its state is
 * in BE (line 5), i.e., when the state was safe: this is the safety of the
 * edges of Equations 2 and 3 of the first paper, where the edges of a state
 * in BE are always safe, since their g no longer changes. The states with a
 * key not smaller than that of \f$ s \f$ are independent of \f$ s \f$ when
 * \f$ w \leq \varepsilon \f$ and \f$ h( t ) \leq h( t , s ) + h( s ) \f$
 * (Theorem 2 of the first paper for OPEN, the same argument for BE, as in
 * section "Algorithm" of GePA*SE), so that line 7 checks only those of
 * smaller key. Hence, when a state is expanded its g is within a factor
 * \f$ \varepsilon \f$ of the optimum (Theorem 1 of the first paper), each
 * state is expanded at most once and each edge evaluated at most once, and
 * the cost of the path found is within \f$ \varepsilon \f$ of the optimum
 * (Theorem 3, with ENDS for the domains with several goals, where the
 * safety of a goal bounds its g only by the cost of a cheapest path to that
 * goal), optimal for \f$ \varepsilon = 1 \f$; a smaller
 * \f$ \varepsilon \f$ is raised to \f$ w \f$. The proofs, for this search
 * and for the anytime one, are on the page of the parallel searches of the
 * documentation.
 *
 * Only the dispatcher scans OPEN, and it wakes one worker per edge, the
 * idle one of least index, so that the workers beyond those that the
 * independent edges keep busy never run, which is the thread management of
 * the first paper (its threads spawned as needed); the states are expanded
 * by the dispatcher itself, which costs little, also while all the workers
 * are busy. Without the dispatcher, every worker selects, expands and
 * evaluates under the lock as in PA*SE, which is the variant "without
 * thread management" of the paper; with one worker the search is the
 * sequential eA* of the paper.
 *
 * GePA*SE evaluates the cheap edges of a state [see HeterogeneousEdgeDomain]
 * together, in one batch of READY given to one worker, and each expensive
 * edge apart: with no expensive edge it is PA*SE on the edges, with no cheap
 * one it is ePA*SE.
 *
 * The anytime search (A-ePA*SE, A-GePA*SE) repairs the solution with the
 * weights \f$ w_0 > w_1 > \dots > w_k = 1 \f$, \f$ w_{i+1} = \max( 1 , w_i -
 * \delta ) \f$, and \f$ \varepsilon = w_i \f$, as ARA* does [see Likhachev,
 * Gordon and Thrun, NIPS 2003]: a state of BE or CLOSED expanded in the
 * current iteration whose g improves goes to INCON instead of OPEN, and a
 * goal selected at line 34 ends the iteration, which waits for the edges
 * being evaluated, publishes the path with the bound \f$ w_i \f$, moves
 * INCON and BE back to OPEN with the keys of \f$ w_{i+1} \f$, and goes on;
 * the edges already evaluated are kept with their states, so that an edge
 * is evaluated at most once over all the iterations.
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

#ifndef __PBeFS_EPASE
 #define __PBeFS_EPASE
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <compare>

#include <condition_variable>

#include <cstddef>

#include <cstdint>

#include <initializer_list>

#include <memory>

#include <mutex>

#include <optional>

#include <tuple>

#include <type_traits>

#include <utility>

#include <vector>

#include "PBeFS/astar/HDAStar.h"

#include "ADS/priority/DaryHeap.h"

#include "ADS/tree/BTreeMap.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS EPASEOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of EPASE

struct EPASEOptions : public ParallelOptions {

 /// the inflation \f$ \varepsilon \geq 1 \f$ of the independence rule,
 /// raised to the weight \f$ w \f$ of SearchOptions if smaller: the cost
 /// is within a factor \f$ \varepsilon \f$ of the optimum; w is that of
 /// the keys, \f$ f' = g + w h \f$ (w-ePA*SE if above 1), and the first
 /// weight of the anytime search
 double eps = 1;

 /// GePA*SE: the cheap edges of a state evaluated together [see
 /// HeterogeneousEdgeDomain]; all the edges are expensive if the Domain
 /// does not tell them apart
 bool generalized = false;

 /// one worker selects and expands, the others evaluate, woken one at a
 /// time (the thread management of ePA*SE); otherwise every worker
 /// selects, expands and evaluates
 bool dispatcher = true;

 /// the states of OPEN examined for a safe one, from the first (0: all)
 std::size_t candidates = 64;

 /// the decrease \f$ \delta > 0 \f$ of the weight of the anytime search
 /// (SearchOptions::anytime), which lowers the weight by it down to 1, with
 /// \f$ \varepsilon \f$ equal to the weight, after each path found, and
 /// records in Solution::improvements the best path of each iteration with
 /// the bound \f$ c / w_i \f$ on the optimal cost
 double weight_step = 0.5;

 /// the seconds after which the anytime search stops with the best path
 /// found (0: no limit)
 double time_limit = 0;

 };  // end( struct( EPASEOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS EPASE -------------------------------*/
/*--------------------------------------------------------------------------*/
/// Edge-based Parallel A* for Slow Evaluations, and its variants
/** Each worker evaluates one edge (or one batch of cheap edges) at a time,
 * outside the lock, while the dispatcher expands the safe states and hands
 * out the edges of the states it has expanded [see the file]. Compared with
 * PASE, an edge waits only for its own evaluation rather than for those of
 * the other edges of its state, the edges of a state are evaluated in
 * parallel, and an edge whose state stops being the best is not evaluated
 * at all: hence, EPASE evaluates fewer edges and keeps more workers busy
 * when the evaluations are slow and the branching factor is smaller than the
 * workers, as in motion planning with collision checks.
 *
 * The statistics count in expanded the states expanded (more than once
 * only in the anytime search, which counts the repeated ones in
 * reexpanded), in evaluated the edges evaluated, in generated the edges
 * that exist, in iterations the iterations of the anytime search, and in
 * expanded_by the expansions and evaluations of each worker, which shows
 * how many workers have run. The constructor raises \f$ \varepsilon \f$ to
 * \f$ w \f$; the anytime search sets it to the weight of each iteration. */

template< PairwiseDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
 requires EdgeDomain< D >
class EPASE {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the EPASE)

 explicit EPASE( const D & domain , const EPASEOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( ! ( f_options.weight >= 1 ) )
   f_options.weight = 1;
  if( ! ( f_options.eps >= f_options.weight ) )
   f_options.eps = f_options.weight;
  if( ! ( f_options.weight_step > 0 ) )
   f_options.weight_step = 0.5;
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution (the best one found by the
 /// anytime search)

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Shared sh( p );
  sh.w = f_options.weight;
  sh.eps = f_options.anytime ? sh.w : f_options.eps;
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );

  {
   const State s0 = f_domain.initial_state();
   ENode * root = w[ 0 ]->f_pool.allocate();
   root->state = s0;
   root->g = Cost( 0 );
   root->h = Cost( f_domain.heuristic( s0 ) );
   root->where = c_open;
   sh.table.find_or_insert( s0 , f_domain.hash( s0 ) ,
                            [ & ] { return( root ); } );
   sh.open.insert( w[ 0 ]->key( root ) , true );
   }

  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.evaluated += x->f_st.evaluated;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded + x->f_st.evaluated );
   st.stored_by.push_back( x->f_pool.size() );
   }
  st.iterations = sh.iterations;
  if( ! sh.best_path.empty() ) {
   sol.status = Status::solved;
   sol.cost = sh.best_cost;
   sol.bound = lower_bound_of( sol.cost , sh.eps );
   sol.path = std::move( sh.best_path );
   }
  else
   if( sh.limit_hit )
    sol.status = Status::node_limit;
  sol.improvements = std::move( sh.improvements );
  st.seconds = sh.clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 // the lower bound c / w on the optimal cost of a path of cost c found
 // within the factor w, rounded up for integer costs
 static Cost lower_bound_of( Cost c , double w ) noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( std::ceil( double( c ) / w - 1e-9 ) ) );
  else
   return( Cost( c / w ) );
  }

 // where a node is
 static constexpr std::uint8_t c_open = 1;       // in OPEN
 static constexpr std::uint8_t c_be = 2;         // in BE
 static constexpr std::uint8_t c_closed = 3;     // expanded, all evaluated

 // the state of an edge
 static constexpr std::uint8_t c_todo = 0;       // not evaluated
 static constexpr std::uint8_t c_feasible = 1;   // evaluated, it exists
 static constexpr std::uint8_t c_infeasible = 2; // evaluated, it does not

 // the edge of READY that stands for the cheap edges of a state, before
 // its expensive edges, which are 1 + their index, as GePA*SE evaluates
 // the cheap edges at the expansion
 static constexpr std::uint32_t c_batch = 0;

 struct ENode;

 // an edge of a state, made at its first expansion
 struct Edge {
  ENode * succ = nullptr;                       // the successor, if any
  Cost c{};                                     // the cost, if evaluated
  std::uint32_t action = 0;                     // the action
  std::uint8_t status = c_todo;                 // evaluated or not
  bool cheap = false;                           // in the batch of GePA*SE
  };

 // a node: one per state
 struct ENode {
  State state{};                                // the state
  Cost g{};                                     // the best g found
  Cost h{};                                     // the heuristic
  ENode * parent = nullptr;                     // the parent on the path
  Cost v{};                                     // g at the last expansion
  double bkey = 0;                              // f' at the last expansion
  Edge * edges = nullptr;                       // the edges, once expanded
  std::uint32_t n_edges = 0;                    // how many
  std::uint32_t pending = 0;                    // edges still to evaluate
  std::uint32_t iter = 0;                       // iteration of the last
                                                // expansion (0: never)
  std::uint8_t where = 0;                       // OPEN, BE or CLOSED
  bool incon = false;                           // in INCON
  };

 // a key of OPEN and of BE: f', then larger g, then the address
 struct Key {
  double f = 0;
  double ng = 0;
  std::uintptr_t node = 0;

  ENode * get( void ) const { return( reinterpret_cast< ENode * >( node ) ); }

  friend auto operator<=>( const Key & , const Key & ) = default;
  };

 // an item of READY: the key of its state, the state and the edge (1 +
 // its index, or c_batch)
 struct Item {
  double f = 0;
  double ng = 0;
  std::uintptr_t node = 0;
  std::uint32_t edge = 0;

  ENode * get( void ) const { return( reinterpret_cast< ENode * >( node ) ); }

  friend auto operator<=>( const Item & , const Item & ) = default;
  };

 // the edges of the expanded states, in blocks that never move
 struct EdgeArena {
  std::vector< std::unique_ptr< Edge[] > > blocks;
  std::size_t used = 0;
  std::size_t cap = 0;

  Edge * get( std::size_t k ) {
   if( used + k > cap ) {
    cap = std::max< std::size_t >( 4096 , k );
    blocks.push_back( std::make_unique< Edge[] >( cap ) );
    used = 0;
    }
   Edge * e = blocks.back().get() + used;
   used += k;
   return( e );
   }
  };

 // what a worker is given by the dispatcher, and the wait for it
 struct alignas( cache_line ) Slot {
  Item job;                                     // the edge or batch
  bool has = false;                             // job is given
  bool idle = false;                            // waiting for a job
  std::condition_variable_any cv;               // job given, or done
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers, all protected by lock

 struct Shared {

  explicit Shared( unsigned p ) : slot( new Slot[ p ] ) , n_slots( p ) {}

  Lock lock;                                    // protects everything
  std::condition_variable_any changed;          // a set or a worker changed
  ADS::BTreeMap< Key , bool > open;            // OPEN
  ADS::BTreeMap< Key , bool > be;              // BE, by the key of expansion
  ADS::DaryHeap< Item > ready;                  // READY
  std::vector< ENode * > incon;                 // INCON
  NodeTable< ENode > table;                     // all the states
  EdgeArena edges;                              // the edges
  std::unique_ptr< Slot[] > slot;               // a slot per worker
  unsigned n_slots;                             // how many
  unsigned busy = 0;                            // edges being evaluated
  double w = 1;                                 // the weight
  double eps = 1;                               // the inflation
  std::uint32_t iter = 1;                       // the current iteration
  std::uint64_t iterations = 0;                 // iterations ended
  bool done = false;                            // the search is over
  bool draining = false;                        // iteration ending
  bool limit_hit = false;                       // the node limit is hit
  ENode * goal = nullptr;                       // the goal of the iteration
  Cost best_cost{};                             // the best path found
  std::vector< State > best_path;
  std::vector< typename SolutionT::Improvement > improvements;  // anytime
  Stopwatch clock;                              // started at the start
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes, its counters and its buffers

 struct alignas( cache_line ) Worker {

  Worker( const EPASE & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) {}

  // a successor found by an evaluation: the edge, whether it exists, the
  // successor, the cost of the edge, its heuristic and its hash
  struct Result {
   std::uint32_t edge;
   bool feasible;
   State t;
   Cost c;
   Cost ht;
   std::uint64_t hash;
   };

  // the key of x in OPEN, with its current g
  Key key( const ENode * x ) const {
   return( Key{ double( x->g ) + f_sh.w * double( x->h ) , - double( x->g ) ,
                reinterpret_cast< std::uintptr_t >( x ) } );
   }

  // the key of x in BE, with its g at the expansion
  Key bkey( const ENode * x ) const {
   return( Key{ x->bkey , - double( x->v ) ,
                reinterpret_cast< std::uintptr_t >( x ) } );
   }

  // true if s is independent of t, whose g is gt: g( s ) - gt <= eps
  // h( t , s ), up to a relative 1e-9 for real costs
  bool independent( const ENode * s , const ENode * t , Cost gt ) const {
   const double lhs = double( s->g - gt );
   const double rhs = f_sh.eps *
    double( f_a.f_domain.pairwise_heuristic( t->state , s->state ) );
   if constexpr( std::is_floating_point_v< Cost > )
    return( lhs <= rhs + 1e-9 * std::max( 1.0 ,
                                          std::abs( double( s->g ) ) ) );
   else
    return( lhs <= rhs );
   }

  // true if s, of key fs, is safe (lock held): independent of the states
  // of BE (with their g at the expansion, which bounds their edges) and of
  // OPEN with a smaller f', the others being independent of s
  bool safe( const ENode * s , double fs ) const {
   for( const auto & k : f_sh.be ) {
    if( ! ( k.first.f < fs ) )
     break;
    const ENode * t = k.first.get();
    if( ! independent( s , t , t->v ) )
     return( false );
    }
   for( const auto & k : f_sh.open ) {
    if( ! ( k.first.f < fs ) )
     break;
    const ENode * t = k.first.get();
    if( ! independent( s , t , t->g ) )
     return( false );
    }
   return( true );
   }

  // SELECT (lock held): the safe state of OPEN before the first edge of
  // READY, if any; otherwise that edge, removed, if take; otherwise nothing
  using Pick = std::pair< ENode * , std::optional< Item > >;

  Pick select( bool take ) {
   const bool has_edge = ! f_sh.ready.empty();
   Item r;
   if( has_edge )
    r = f_sh.ready.top();
   std::size_t n = 0;
   const std::size_t k_max = f_a.f_options.candidates;
   for( const auto & k : f_sh.open ) {
    if( has_edge && ! ( std::tie( k.first.f , k.first.ng ) <
                        std::tie( r.f , r.ng ) ) )
     break;
    if( k_max && ( n++ == k_max ) )
     break;
    ENode * s = k.first.get();
    if( safe( s , k.first.f ) &&
        ( ! f_a.f_domain.is_goal( s->state ) || may_end( s ) ) )
     return( Pick( s , std::nullopt ) );
    }
   if( has_edge && take )
    return( Pick( nullptr , f_sh.ready.pop() ) );
   return( Pick( nullptr , std::nullopt ) );
   }

  // true if the goal s may end the search (lock held): no state of OPEN or
  // BE has a key smaller than g( s ), so that every cheapest path to a goal
  // has a state of key at least g( s ), and g( s ) <= eps C*; the safety
  // of s alone bounds g( s ) by eps g*( s ), which is not eps C* when there
  // are several goals
  bool may_end( const ENode * s ) const {
   const double gs = double( s->g );
   return( ! ( ( *f_sh.open.begin() ).first.f < gs ) &&
           ( f_sh.be.empty() || ! ( ( *f_sh.be.begin() ).first.f < gs ) ) );
   }

  // a better g' for the known node x, reached from s (lock held)
  void improve( ENode * s , ENode * x , Cost gx ) {
   if( ! ( gx < x->g ) ) {
    ++f_st.duplicates;
    return;
    }
   if( x->where == c_open ) {
    f_sh.open.remove( key( x ) );
    x->g = gx;
    x->parent = s;
    f_sh.open.insert( key( x ) , true );
    return;
    }
   if( ! f_a.f_options.anytime ) {     // in BE or CLOSED, never improved
    ++f_st.duplicates;
    return;
    }
   x->g = gx;
   x->parent = s;
   if( x->iter == f_sh.iter ) {        // expanded in this iteration
    if( ! x->incon ) {
     x->incon = true;
     f_sh.incon.push_back( x );
     }
    }
   else {                              // closed in an earlier one
    x->where = c_open;
    f_sh.open.insert( key( x ) , true );
    }
   }

  // moves s from BE to CLOSED (lock held)
  void close( ENode * s ) {
   f_sh.be.remove( bkey( s ) );
   s->where = c_closed;
   }

  // EXPAND (lock held): s from OPEN to BE, its edges to READY
  void expand( ENode * s ) {
   f_sh.open.remove( key( s ) );
   s->where = c_be;
   s->v = s->g;
   s->bkey = key( s ).f;
   if( s->iter )
    ++f_st.reexpanded;
   s->iter = f_sh.iter;
   ++f_st.expanded;
   f_sh.be.insert( bkey( s ) , true );
   if( ! s->edges ) {
    f_actions.clear();
    f_a.f_domain.for_each_action( s->state , [ & ]( std::size_t a ) {
     f_actions.push_back( std::uint32_t( a ) );
     } );
    s->n_edges = std::uint32_t( f_actions.size() );
    s->edges = f_sh.edges.get( f_actions.size() );
    for( std::uint32_t i = 0 ; i < s->n_edges ; ++i ) {
     Edge & e = s->edges[ i ];
     e.action = f_actions[ i ];
     if constexpr( HeterogeneousEdgeDomain< D > )
      e.cheap = f_a.f_options.generalized &&
                ! f_a.f_domain.expensive_edge( s->state , e.action );
     }
    }
   s->pending = 0;
   std::uint32_t cheap = 0;
   const auto ptr = reinterpret_cast< std::uintptr_t >( s );
   for( std::uint32_t i = 0 ; i < s->n_edges ; ++i ) {
    Edge & e = s->edges[ i ];
    if( e.status == c_feasible )       // evaluated in an earlier iteration
     improve( s , e.succ , s->g + e.c );
    else
     if( e.status == c_todo ) {
      ++s->pending;
      if( e.cheap )
       ++cheap;
      else
       f_sh.ready.push( Item{ s->bkey , - double( s->v ) , ptr , i + 1 } );
      }
    }
   if( cheap )
    f_sh.ready.push( Item{ s->bkey , - double( s->v ) , ptr , c_batch } );
   if( s->pending == 0 )
    close( s );
   }

  // the edges of the job and the state to evaluate them from (lock held)
  void take( const Item & job ) {
   ENode * s = job.get();
   f_s = s;
   f_state = s->state;
   f_h = s->h;
   f_todo.clear();
   if( job.edge != c_batch )
    f_todo.emplace_back( job.edge - 1 , s->edges[ job.edge - 1 ].action );
   else
    for( std::uint32_t i = 0 ; i < s->n_edges ; ++i )
     if( s->edges[ i ].cheap && ( s->edges[ i ].status == c_todo ) )
      f_todo.emplace_back( i , s->edges[ i ].action );
   }

  // the evaluation of the edges taken, without the lock
  void evaluate( void ) {
   f_res.clear();
   for( const auto & [ i , a ] : f_todo ) {
    bool found = false;
    f_a.f_domain.evaluate_edge( f_state , f_h , a ,
                                [ & ]( const State & t , Cost c , Cost ht ) {
     found = true;
     f_res.push_back( Result{ i , true , t , c , ht ,
                              f_a.f_domain.hash( t ) } );
     } );
    if( ! found )
     f_res.push_back( Result{ i , false , State{} , Cost( 0 ) , Cost( 0 ) ,
                              0 } );
    ++f_st.evaluated;
    }
   }

  // the results of the evaluation applied (lock held)
  void apply( void ) {
   ENode * s = f_s;
   for( const auto & r : f_res )
    if( r.feasible )
     f_sh.table.prefetch( r.hash );
   for( const auto & r : f_res ) {
    Edge & e = s->edges[ r.edge ];
    if( ! r.feasible ) {
     e.status = c_infeasible;
     continue;
     }
    ++f_st.generated;
    const Cost gt = s->g + r.c;
    auto [ x , fresh ] = f_sh.table.find_or_insert( r.t , r.hash , [ & ] {
     ENode * y = f_pool.allocate();
     y->state = r.t;
     y->g = gt;
     y->h = r.ht;
     y->parent = s;
     y->where = c_open;
     return( y );
     } );
    e.succ = x;
    e.c = r.c;
    e.status = c_feasible;
    if( fresh )
     f_sh.open.insert( key( x ) , true );
    else
     improve( s , x , gt );
    }
   s->pending -= std::uint32_t( f_res.size() );
   if( ( s->pending == 0 ) && ( s->where == c_be ) )
    close( s );
   if( f_a.f_options.node_limit &&
       ( f_pool.size() >= f_a.f_options.node_limit ) ) {
    f_sh.limit_hit = true;
    stop();
    }
   }

  // the end of the search (lock held)
  void stop( void ) {
   f_sh.done = true;
   for( unsigned i = 0 ; i < f_sh.n_slots ; ++i )
    f_sh.slot[ i ].cv.notify_all();
   f_sh.changed.notify_all();
   }

  // a goal s selected safe (lock held): the end of the search, or of the
  // iteration of the anytime one
  void reach( ENode * s ) {
   if( f_sh.best_path.empty() || ( s->g < f_sh.best_cost ) ) {
    f_sh.goal = s;
    }
   if( ! f_a.f_options.anytime ) {
    f_sh.best_cost = s->g;
    f_sh.best_path = trace_path( s );
    stop();
    return;
    }
   f_sh.draining = true;
   }

  // the end of an iteration of the anytime search, no edge being evaluated
  // (lock held): the path published, the weight lowered, INCON and BE back
  // into OPEN with the new keys
  void next_iteration( void ) {
   ENode * goal = f_sh.goal;
   if( f_sh.best_path.empty() || ( goal->g < f_sh.best_cost ) ) {
    f_sh.best_cost = goal->g;
    f_sh.best_path = trace_path( goal );
    }
   ++f_sh.iterations;
   f_sh.improvements.push_back( { f_sh.best_cost ,
                                  lower_bound_of( f_sh.best_cost , f_sh.w ) ,
                                  f_st.expanded , f_sh.clock.seconds() } );
   if( ( f_sh.w <= 1 ) || out_of_time() ) {
    stop();
    return;
    }
   f_sh.w = std::max( 1.0 , f_sh.w - f_a.f_options.weight_step );
   f_sh.eps = f_sh.w;
   f_nodes.clear();
   for( const auto & k : f_sh.open )
    f_nodes.push_back( k.first.get() );
   for( const auto & k : f_sh.be )
    f_nodes.push_back( k.first.get() );
   for( ENode * x : f_sh.incon )
    if( x->where == c_closed )         // those in BE are taken above
     f_nodes.push_back( x );
   f_sh.open.clear();
   f_sh.be.clear();
   f_sh.ready.clear();
   f_sh.incon.clear();
   ++f_sh.iter;
   for( ENode * x : f_nodes ) {
    x->where = c_open;
    x->incon = false;
    f_sh.open.insert( key( x ) , true );
    }
   f_sh.draining = false;
   f_sh.changed.notify_all();
   }

  // true if the anytime search has run out of time
  bool out_of_time( void ) const {
   return( f_a.f_options.anytime && ( f_a.f_options.time_limit > 0 ) &&
           ( f_sh.clock.seconds() >= f_a.f_options.time_limit ) );
   }

  // true if nothing is left to do (lock held)
  bool exhausted( void ) const {
   return( f_sh.open.empty() && f_sh.ready.empty() && f_sh.be.empty() &&
           ( f_sh.busy == 0 ) );
   }

  // the idle worker of least index, or 0 if none (lock held)
  unsigned idle_worker( void ) const {
   for( unsigned i = 1 ; i < f_sh.n_slots ; ++i )
    if( f_sh.slot[ i ].idle )
     return( i );
   return( 0 );
   }

  // the work of the worker
  void run( void ) {
   std::unique_lock< Lock > guard( f_sh.lock );
   if( ( f_sh.n_slots == 1 ) || ! f_a.f_options.dispatcher )
    run_all( guard );
   else
    if( f_id == 0 )
     run_dispatcher( guard );
    else
     run_evaluator( guard );
   }

  // the dispatcher: selects, expands, gives the edges to the idle workers
  void run_dispatcher( std::unique_lock< Lock > & guard ) {
   while( ! f_sh.done ) {
    if( out_of_time() ) {
     stop();
     break;
     }
    if( f_sh.draining ) {
     if( f_sh.busy == 0 )
      next_iteration();
     else
      f_sh.changed.wait( guard );
     continue;
     }
    const unsigned i = idle_worker();
    auto [ s , e ] = select( i != 0 );
    if( s ) {
     if( f_a.f_domain.is_goal( s->state ) )
      reach( s );
     else
      expand( s );
     continue;
     }
    if( ! e ) {
     if( exhausted() ) {
      stop();
      break;
      }
     f_sh.changed.wait( guard );
     continue;
     }
    Slot & x = f_sh.slot[ i ];
    x.job = *e;
    x.has = true;
    x.idle = false;
    ++f_sh.busy;
    x.cv.notify_one();
    }
   }

  // a worker that evaluates what the dispatcher gives it
  void run_evaluator( std::unique_lock< Lock > & guard ) {
   Slot & x = f_sh.slot[ f_id ];
   x.idle = true;
   f_sh.changed.notify_all();
   while( true ) {
    while( ! x.has && ! f_sh.done )
     x.cv.wait( guard );
    if( f_sh.done )
     break;
    x.has = false;
    take( x.job );
    guard.unlock();
    evaluate();
    guard.lock();
    if( ! f_sh.done )
     apply();
    --f_sh.busy;
    x.idle = true;
    f_sh.changed.notify_all();
    }
   }

  // a worker that selects, expands and evaluates, as in PA*SE
  void run_all( std::unique_lock< Lock > & guard ) {
   while( ! f_sh.done ) {
    if( out_of_time() ) {
     stop();
     break;
     }
    if( f_sh.draining ) {
     if( f_sh.busy == 0 )
      next_iteration();
     else
      f_sh.changed.wait( guard );
     continue;
     }
    auto [ s , e ] = select( true );
    if( s ) {
     if( f_a.f_domain.is_goal( s->state ) )
      reach( s );
     else {
      expand( s );
      f_sh.changed.notify_all();
      }
     continue;
     }
    if( ! e ) {
     if( exhausted() ) {
      stop();
      break;
      }
     f_sh.changed.wait( guard );
     continue;
     }
    ++f_sh.busy;
    take( *e );
    guard.unlock();
    evaluate();
    guard.lock();
    if( ! f_sh.done )
     apply();
    --f_sh.busy;
    f_sh.changed.notify_all();
    }
   }

  const EPASE & f_a;                            // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  NodePool< ENode > f_pool;                     // the nodes it allocated
  SearchStats f_st;                             // its counters
  ENode * f_s = nullptr;                        // the state of the job
  State f_state{};                              // its state, copied
  Cost f_h{};                                   // and its heuristic
  std::vector< std::pair< std::uint32_t , std::uint32_t > > f_todo;
  std::vector< Result > f_res;                  // the results of the job
  std::vector< std::uint32_t > f_actions;       // the actions of a state
  std::vector< ENode * > f_nodes;               // the nodes moved to OPEN
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 EPASEOptions f_options;                               ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( EPASE ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_EPASE */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File EPASE.h ----------------------------*/
/*--------------------------------------------------------------------------*/
