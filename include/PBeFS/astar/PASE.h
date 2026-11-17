/*--------------------------------------------------------------------------*/
/*------------------------------ File PASE.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template PASE, Parallel A* for Slow Expansions
 * of Phillips, Likhachev and Koenig (ICAPS 2014).
 *
 * The search follows Algorithm 1 of Phillips, Likhachev and Koenig (ICAPS
 * 2014), with the notation of the page "Best-first search and A*" of the
 * documentation: BE is the set of the states being expanded, OPEN is an
 * ordered set by increasing \f$ f' = g + w h \f$ and, for equal \f$ f' \f$,
 * by decreasing \f$ g \f$, \f$ h( t , s ) \f$ is the pairwise heuristic of
 * the Domain, \f$ w \geq 1 \f$ the weight of the heuristic (1 for PA*SE,
 * above 1 for wPA*SE), \f$ \varepsilon_p \geq w \f$ the inflation of the
 * independence rule, \f$ K \f$ the number of states of OPEN examined (all,
 * if 0), and one lock \f$ L \f$ protects OPEN, BE and CLOSED:
 *
 * \code{.unparsed}
 *  1  procedure PA*SE(s₀, p)
 *  2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {s₀};  BE ← ∅
 *  3     in parallel, on each worker do
 *  4        lock L
 *  5        while not done do
 *  6           s ← the first state of OPEN, among its first K, that is SAFE
 *                  and, if a goal, ENDS
 *  7           if there is none then
 *  8              if OPEN = ∅ and BE = ∅ then done ← true;  break
 *  9              wait on L until OPEN or BE change;  continue
 * 10           remove s from OPEN
 * 11           if s is a goal then record s;  done ← true;  break
 * 12           add s to BE;  close s;  unlock L
 * 13           S ← the successors of s, with their costs   ▷ the slow part
 * 14           lock L
 * 15           for each (v, c(s, v)) ∈ S do
 * 16              g' ← g(s) + c(s, v)
 * 17              if v ∉ CLOSED then
 * 18                 add v to CLOSED and OPEN, with g(v) ← g', parent(v) ← s
 * 19              else if v is not closed and g' < g(v) then
 * 20                 g(v) ← g';  parent(v) ← s;  reorder v in OPEN
 * 21           remove s from BE;  wake the waiting workers
 * 22        unlock L
 * 23  procedure SAFE(s)
 * 24     for each t ∈ BE with f'(t) < f'(s) do
 * 25        if g(s) − g(t) > ε_p h(t, s) then return false
 * 26     for each t ∈ OPEN with f'(t) < f'(s) do
 * 27        if g(s) − g(t) > ε_p h(t, s) then return false
 * 28     return true
 * 29  procedure ENDS(s)
 * 30     return no state t of OPEN ∪ BE has g(t) + ε_p h(t) < g(s)
 * \endcode
 *
 * A state \f$ s \f$ is independent of \f$ t \f$ when
 * \f$ g( s ) - g( t ) \leq \varepsilon_p \, h( t , s ) \f$, since then no
 * path through \f$ t \f$ can lower \f$ g( s ) \f$ below
 * \f$ g( s ) / \varepsilon_p \f$. If \f$ f'( s ) \leq f'( t ) \f$, then
 * \f$ g( s ) - g( t ) \leq w ( h( t ) - h( s ) ) \leq w \, h( t , s ) \leq
 * \varepsilon_p \, h( t , s ) \f$, provided that
 * \f$ h( t ) \leq h( t , s ) + h( s ) \f$ (which a heuristic computed as
 * the pairwise one to the nearest goal satisfies): by Theorem 2 of the paper
 * only the states of OPEN with a smaller \f$ f' \f$ have to be checked (line
 * 26), and by the same argument only those of BE (line 24), whose \f$ g \f$
 * no longer changes. With \f$ \varepsilon_p = 1 \f$ a state is expanded
 * only with its optimal \f$ g \f$, hence at most once, and a closed state
 * is never improved (line 19); the goal is returned as soon as it is safe
 * (line 11), with an optimal cost, or within a factor \f$ \varepsilon_p \f$
 * of the optimum otherwise (Theorems 1 and 3 of the paper, which hold for
 * \f$ w \leq \varepsilon_p \f$; a smaller \f$ \varepsilon_p \f$ is
 * raised to \f$ w \f$). Theorem 3 bounds \f$ g \f$ of the goal by
 * \f$ \varepsilon_p \f$ times the cost of a cheapest path to that goal,
 * which with several goals may exceed \f$ C^* \f$; hence a goal is taken
 * only if also no state \f$ t \f$ of OPEN or BE has
 * \f$ g( t ) + \varepsilon_p h( t ) \f$ smaller than its \f$ g \f$ (line 6):
 * a cheapest path to any goal has a first state in OPEN or BE, with
 * \f$ g + \varepsilon_p h \leq \varepsilon_p C^* \f$ by Theorem 1 and the
 * admissibility of \f$ h \f$. Meanwhile, a state with
 * \f$ g + \varepsilon_p h \geq z \f$, \f$ z \f$ the least \f$ g \f$ of a goal
 * in OPEN, is not expanded: if it is that first state, then
 * \f$ z \leq \varepsilon_p C^* \f$ already. The scan at line 6 may be
 * limited to the first \f$ K \f$ states of OPEN (64 by default), which
 * bounds the time spent under the lock: a state further down would have
 * to be checked against a longer prefix of OPEN, and with BE empty the
 * first state is always safe, so that the search goes on; the limit
 * changes the order of the expansions, not the bound on the cost. For real
 * costs the test at lines 25 and 27 allows a relative error of
 * \f$ 10^{-9} \f$, so that the rounding does not make the goal dependent
 * on the states of equal \f$ f' \f$. With integer
 * costs the weight is the fraction \f$ w = a / b \f$, \f$ b = 1 \f$ if
 * \f$ w \f$ is an integer and 256 otherwise (\f$ a \f$ rounded up, so that
 * the weight is not smaller than the one given), and the key of a state is
 * \f$ b g + a h \f$, exact, packed with \f$ g \f$ in one word [see
 * packed_fg]; then \f$ b g + a h < 2^{32} \f$ is required.
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

#ifndef __PBeFS_PASE
 #define __PBeFS_PASE
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

#include <stdexcept>

#include <tuple>

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
/*--------------------------- CLASS PASEOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of PASE

struct PASEOptions : public ParallelOptions {

 /// the inflation \f$ \varepsilon_p \geq 1 \f$ of the independence rule:
 /// with \f$ \varepsilon_p = 1 \f$ the search is optimal, otherwise the
 /// cost is within a factor \f$ \varepsilon_p \f$ of the optimum; it is
 /// raised to the weight \f$ w \f$ of SearchOptions if smaller, which
 /// orders OPEN by \f$ f' = g + w h \f$ (wPA*SE if above 1)
 double eps = 1;

 /// the states of OPEN examined for a safe one, from the first (0: all)
 std::size_t candidates = 64;

 /// once a goal of cost \f$ z \f$ is in OPEN, no state with
 /// \f$ g + \varepsilon_p h \geq z \f$ is expanded, which keeps the bound
 bool prune = true;

 };  // end( struct( PASEOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS PASE -------------------------------*/
/*--------------------------------------------------------------------------*/
/// Parallel A* for Slow Expansions
/** A state \f$ s \f$ is independent of a state \f$ s' \f$ if
 * \f$ g( s ) - g( s' ) \leq \varepsilon_p \, h( s' , s ) \f$,
 * \f$ h( s' , s ) \f$ being the forward-backward consistent pairwise
 * heuristic of the Domain: then no path through \f$ s' \f$ can lower
 * \f$ g( s ) \f$ below \f$ g( s ) / \varepsilon_p \f$. A state is safe to
 * expand if it is independent of all the states being expanded (the set BE)
 * with a smaller \f$ f' = g + w h \f$ and of all the states of OPEN with a
 * smaller \f$ f' \f$ [see Theorem 2 of the paper]. Each worker takes from
 * OPEN the safe state of smallest \f$ f' \f$, expands
 * it outside the lock and puts back the successors whose g has improved; a
 * state is expanded at most once, and when it is expanded its g is optimal
 * (with \f$ \varepsilon_p = 1 \f$). The search ends when the goal is safe to
 * expand, or when OPEN and BE are both empty. This is Algorithm 1 of the
 * paper, and with \f$ w > 1 \f$ the wPA*SE of section 4 of the paper;
 * OPEN is an ordered set scanned in increasing \f$ f' \f$, and a worker
 * that finds no safe state waits until OPEN or BE change.
 *
 * Finding a safe state costs a scan of the states of OPEN with a smaller f
 * and of BE, all under one lock: PASE pays only when the expansions are
 * slow, as in motion planning with collision checks. The suboptimality is
 * that of \f$ \varepsilon_p \f$, which the constructor raises to the
 * weight. Observer is called at the events of the search [see
 * SearchObserver.h], the key of OPEN being \f$ g + w h \f$ (rounded down
 * for integer costs), the improvement of a state of OPEN being a pop and
 * a push, and the goal that ends the search its solution, which is not an
 * expansion; the root is inserted by the worker 0. */

template< PairwiseDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          class Observer = NullObserver >
class PASE {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the PASE)

 explicit PASE( const D & domain , const PASEOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( ! ( f_options.weight >= 1 ) )
   f_options.weight = 1;
  if constexpr( std::is_integral_v< Cost > ) {
   const double w = f_options.weight;
   if( w == std::floor( w ) )
    f_wa = std::uint64_t( w );
   else {
    f_wb = 256;
    f_wa = std::uint64_t( std::ceil( w * 256 - 1e-9 ) );
    }
   f_options.weight = double( f_wa ) / double( f_wb );
   }
  if( ! ( f_options.eps >= f_options.weight ) )
   f_options.eps = f_options.weight;
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the observer of the searches, which must outlive them

 void set_observer( Observer & o ) noexcept { f_obs.set( o ); }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  Observer & obs = f_obs.get( "PASE" );
  obs.on_start( p );
  Shared sh;
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i , obs );

  {
   const State s0 = f_domain.initial_state();
   NodeT * root = w[ 0 ]->f_pool.allocate();
   *root = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) ,
                  nullptr , false , false };
   sh.table.find_or_insert( s0 , f_domain.hash( s0 ) ,
                            [ & ] { return( root ); } );
   sh.open.insert( Key( root , *this ) , true );
   obs.on_open_push( 0 , sh.open.size() );
   if( f_domain.is_goal( s0 ) )
    sh.z = Cost( 0 );
   }

  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  if( sh.limit_hit )
   sol.status = Status::node_limit;
  else
   if( sh.goal ) {
    sol.status = Status::solved;
    sol.cost = sh.goal->g;
    if constexpr( std::is_integral_v< Cost > )
     sol.bound = Cost( std::ceil( double( sol.cost ) / f_options.eps -
                                  1e-9 ) );
    else
     sol.bound = sol.cost / f_options.eps;
    sol.path = trace_path( sh.goal );
    }
  st.seconds = clock.seconds();
  obs.on_finish();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 // the f' of x, b g + a h for integer costs, g + w h otherwise
 auto fw( const NodeT * x ) const noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( f_wb * std::uint64_t( x->g ) + f_wa * std::uint64_t( x->h ) );
  else
   return( x->g + Cost( f_options.weight ) * x->h );
  }

 // the f' of x on the scale of the costs, rounded down for integer costs
 Cost key( const NodeT * x ) const noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( fw( x ) / f_wb ) );
  else
   return( fw( x ) );
  }

 // g( x ) on the scale of fw()
 auto fw0( const NodeT * x ) const noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( f_wb * std::uint64_t( x->g ) );
  else
   return( x->g );
  }

 // OPEN and BE: by f', then by larger g, then by address; with integer
 // costs f' and g are one word [see packed_fg], compared at once
 struct Key {
  std::conditional_t< std::is_integral_v< Cost > , std::uint64_t ,
                      std::pair< Cost , Cost > > fg{};
  std::uintptr_t node = 0;

  Key( void ) = default;

  Key( const NodeT * x , const PASE & a )
   : node( reinterpret_cast< std::uintptr_t >( x ) ) {
   if constexpr( std::is_integral_v< Cost > )
    fg = packed_fg( x->g , Cost( a.fw( x ) - std::uint64_t( x->g ) ) );
   else
    fg = { a.fw( x ) , - x->g };
   }

  NodeT * get( void ) const { return( reinterpret_cast< NodeT * >( node ) ); }

  friend auto operator<=>( const Key & , const Key & ) = default;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers, all protected by lock

 struct Shared {
  Lock lock;                                    // protects everything
  std::condition_variable_any changed;          // OPEN or BE changed
  ADS::BTreeMap< Key , bool > open;            // OPEN
  ADS::BTreeMap< Key , bool > be;              // being expanded
  NodeTable< NodeT > table;                     // all the states
  bool done = false;                            // the search is over
  bool limit_hit = false;                       // the node limit is hit
  NodeT * goal = nullptr;                       // the goal, once safe
  Cost z = infinite_cost< Cost >();             // the least g of a goal
                                                // in OPEN
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const PASE & a , Shared & s , unsigned i , Observer & o )
   : f_a( a ) , f_sh( s ) , f_obs( o ) , f_id( i ) {}

  // true if s is independent of t: g( s ) - g( t ) <= eps h( t , s ), up
  // to a relative 1e-9 for real costs, whose rounding would otherwise make
  // the goal dependent on the states with f' = f'( goal )
  bool independent( const NodeT * s , const NodeT * t ) const {
   const double lhs = double( s->g - t->g );
   const double rhs = f_a.f_options.eps *
    double( f_a.f_domain.pairwise_heuristic( t->state , s->state ) );
   if constexpr( std::is_floating_point_v< Cost > )
    return( lhs <= rhs + 1e-9 * std::max( 1.0 ,
                                          std::abs( double( s->g ) ) ) );
   else
    return( lhs <= rhs );
   }

  // true if x is pruned (lock held): not a goal, and g + eps h not
  // smaller than the least g of a goal in OPEN, so that no goal cheaper
  // than that one is reached through x
  bool pruned( const NodeT * x ) const {
   return( f_a.f_options.prune && ! f_a.f_domain.is_goal( x->state ) &&
           ! ( bound( x ) < double( f_sh.z ) ) );
   }

  // true if s is safe to expand (lock held): independent of the states
  // of BE and of OPEN with a smaller f', which come first in both, the
  // others being independent of s by Theorem 2 of the paper; the states of
  // equal f' are not checked, so that a layer of equal f' costs no scan,
  // nor the states of OPEN pruned. A state s expanded with
  // g( s ) > eps g*( s ) depends on the first state q in OPEN of a cheapest
  // path to it; if q is pruned and s is on a cheapest path to a goal, then
  // so is q, and z <= g( q ) + eps h( q ) <= eps C*, which the goal of g z
  // in OPEN keeps for the end [see may_end()]. With BE empty the first
  // state of OPEN not pruned is then safe, and if a goal it has g z and
  // may end, so that the workers never all wait with OPEN not empty
  bool safe( const NodeT * s ) const {
   const auto fs = f_a.fw( s );
   for( const auto & k : f_sh.be ) {
    const NodeT * t = k.first.get();
    if( ! ( f_a.fw( t ) < fs ) )
     break;
    if( ! independent( s , t ) )
     return( false );
    }
   for( const auto & k : f_sh.open ) {
    const NodeT * t = k.first.get();
    if( ! ( f_a.fw( t ) < fs ) )
     break;
    if( ! pruned( t ) && ! independent( s , t ) )
     return( false );
    }
   return( true );
   }

  // true if the goal s may end the search (lock held): no state t of OPEN
  // or BE has g( t ) + eps h( t ) < g( s ), so that the first state t of a
  // cheapest path to a goal, which has g( t ) + eps h( t ) <= eps C*, gives
  // g( s ) <= eps C*; the states of f' not smaller than g( s ) satisfy it,
  // since w <= eps. The safety of s alone bounds g( s ) by eps g*( s ),
  // which is not eps C* when there are several goals
  bool may_end( const NodeT * s ) const {
   const auto gs = f_a.fw0( s );
   for( const auto * set : { & f_sh.open , & f_sh.be } )
    for( const auto & k : *set ) {
     const NodeT * t = k.first.get();
     if( ! ( f_a.fw( t ) < gs ) )
      break;
     if( bound( t ) < double( s->g ) )
      return( false );
     }
   return( true );
   }

  // g + eps h of x: a cheapest path to a goal whose first state in OPEN or
  // BE is x costs at least this over eps
  double bound( const NodeT * x ) const {
   return( double( x->g ) + f_a.f_options.eps * double( x->h ) );
   }

  // Algorithm 1 of the paper
  void run( void ) {
   std::vector< std::tuple< State , Cost , Cost , std::uint64_t ,
                            std::uint8_t > > succ;
   std::unique_lock< Lock > guard( f_sh.lock );
   while( ! f_sh.done ) {
    // the safe state of smallest f' among the first candidates of OPEN
    // not pruned (all, if 0), if any: further down, each candidate would
    // have to be checked against the larger and larger part of OPEN before
    // it, under the lock; with BE empty the first one is safe, and the
    // search goes on
    NodeT * s = nullptr;
    std::size_t n = 0;
    const std::size_t k_max = f_a.f_options.candidates;
    for( const auto & k : f_sh.open ) {
     NodeT * x = k.first.get();
     if( pruned( x ) )
      continue;                         // no cheaper goal through x
     if( k_max && ( n++ == k_max ) )
      break;
     const bool goal = f_a.f_domain.is_goal( x->state );
     if( safe( x ) && ( ! goal || may_end( x ) ) ) {
      s = x;
      break;
      }
     }
    if( ! s ) {
     if( f_sh.open.empty() && f_sh.be.empty() ) {
      f_sh.done = true;
      f_sh.changed.notify_all();
      break;
      }
     f_sh.changed.wait( guard );
     continue;
     }
    f_sh.open.remove( Key( s , f_a ) );
    f_obs.on_open_pop( f_id , f_sh.open.size() );
    if( f_a.f_domain.is_goal( s->state ) ) {
     f_obs.on_solution( f_id , s->g );
     f_sh.goal = s;               // safe, hence of optimal g
     f_sh.done = true;
     f_sh.changed.notify_all();
     break;
     }
    f_sh.be.insert( Key( s , f_a ) , true );
    s->closed = s->expanded = true;
    const State st = s->state;
    const Cost g = s->g , h = s->h;
    const std::uint8_t op = s->op;
    guard.unlock();

    ++f_st.expanded;
    if constexpr( observes< Observer > )
     f_obs.on_expand( f_id , st , g , h , f_a.key( s ) );
    succ.clear();
    for_each_child( f_a.f_domain , st , h , op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t ot ) {
     succ.emplace_back( t , g + c , ht , f_a.f_domain.hash( t ) , ot );
     } );

    // the hashes computed without the lock, and the places of all the
    // successors prefetched under it before any lookup, so that their
    // misses overlap and the lock is held for about one of them
    guard.lock();
    for( const auto & sc : succ )
     f_sh.table.prefetch( std::get< 3 >( sc ) );
    for( const auto & sc : succ ) {    // no structured binding: OpenMP
     const State & t = std::get< 0 >( sc );     // cannot capture one
     const Cost gt = std::get< 1 >( sc ) , ht = std::get< 2 >( sc );
     ++f_st.generated;
     f_obs.on_generate( f_id , t , & st , gt , ht );
     auto [ x , fresh ] = f_sh.table.find_or_insert( t ,
                                                     std::get< 3 >( sc ) ,
                                                     [ & ] {
      NodeT * y = f_pool.allocate();
      *y = NodeT{ t , gt , ht , s , false , false , std::get< 4 >( sc ) };
      return( y );
      } );
     if( fresh ) {
      f_sh.open.insert( Key( x , f_a ) , true );
      f_obs.on_open_push( f_id , f_sh.open.size() );
      }
     else
      if( ! x->closed && ( gt < x->g ) ) {
       f_sh.open.remove( Key( x , f_a ) );
       f_obs.on_open_pop( f_id , f_sh.open.size() );
       x->g = gt;
       x->parent = s;
       x->op = std::get< 4 >( sc );
       f_sh.open.insert( Key( x , f_a ) , true );
       f_obs.on_open_push( f_id , f_sh.open.size() );
       }
      else {
       ++f_st.duplicates;
       f_obs.on_duplicate( f_id , t , gt );
       continue;
       }
     if( ( x->g < f_sh.z ) && f_a.f_domain.is_goal( t ) )
      f_sh.z = x->g;
     }
    f_sh.be.remove( Key( s , f_a ) );
    if( f_a.f_options.node_limit &&
        ( f_pool.size() >= f_a.f_options.node_limit ) ) {
     f_sh.limit_hit = true;
     f_sh.done = true;
     }
    f_sh.changed.notify_all();
    }
   }

  const PASE & f_a;                             // the search
  Shared & f_sh;                                // the shared data
  NodePool< NodeT > f_pool;                     // the nodes it allocated
  SearchStats f_st;                             // its counters
  [[no_unique_address]] observer_ref< Observer > f_obs;  // the observer
  unsigned f_id;                                // the index of the worker
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 PASEOptions f_options;                                ///< the options
 std::uint64_t f_wa = 1;                               ///< w = a / b, the a
 std::uint64_t f_wb = 1;                               ///< w = a / b, the b
 Executor * f_ex = nullptr;                            ///< given, or nullptr
 [[no_unique_address]] ObserverSlot< Observer > f_obs;  ///< the observer

 };  // end( class( PASE ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PASE */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File PASE.h -----------------------------*/
/*--------------------------------------------------------------------------*/
