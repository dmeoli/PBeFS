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
 * ordered set by increasing \f$ f \f$ and, for equal \f$ f \f$, by
 * decreasing \f$ g \f$, \f$ h( t , s ) \f$ is the pairwise heuristic of the
 * Domain, \f$ \varepsilon_p \geq 1 \f$ the inflation of the independence
 * rule, and one lock \f$ L \f$ protects OPEN, BE and CLOSED:
 *
 * \code{.unparsed}
 *  1  procedure PA*SE(s₀, p)
 *  2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {s₀};  BE ← ∅
 *  3     in parallel, on each worker do
 *  4        lock L
 *  5        while not done do
 *  6           s ← the first state of OPEN, among its first 64, that is SAFE
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
 * 24     for each t ∈ BE do
 * 25        if g(s) − g(t) > ε_p h(t, s) then return false
 * 26     for each t ∈ OPEN with f(t) < f(s) do
 * 27        if g(s) − g(t) > ε_p h(t, s) then return false
 * 28     return true
 * \endcode
 *
 * A state \f$ s \f$ is independent of \f$ t \f$ when
 * \f$ g( s ) - g( t ) \leq \varepsilon_p \, h( t , s ) \f$, since then no
 * path through \f$ t \f$ can lower \f$ g( s ) \f$ below
 * \f$ g( s ) / \varepsilon_p \f$; by Theorem 2 of the paper only the states
 * of OPEN with a smaller \f$ f \f$ have to be checked (line 26). With
 * \f$ \varepsilon_p = 1 \f$ a state is expanded only with its optimal
 * \f$ g \f$, hence at most once, and a closed state is never improved
 * (line 19); the goal is returned as soon as it is safe (line 11), with an
 * optimal cost, or within a factor \f$ \varepsilon_p \f$ of the optimum
 * otherwise (Theorem 3 of the paper). The scan at line 6 is limited to the
 * first 64 states of OPEN, which bounds the time spent under the lock: a
 * state further down would have to be checked against a longer prefix of
 * OPEN, and with BE empty the first state is always safe, so that the search
 * goes on. For real costs the test at lines 25 and 27 allows a relative
 * error of \f$ 10^{-9} \f$, so that the rounding does not make the goal
 * dependent on the states of equal \f$ f \f$.
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

#ifndef __PBeFS_PASE
 #define __PBeFS_PASE
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <condition_variable>

#include <memory>

#include <mutex>

#include <tuple>

#include <type_traits>

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
 /// cost is within a factor \f$ \varepsilon_p \f$ of the optimum
 double eps = 1;

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
 * and of all the states of OPEN with a smaller f [see Theorem 2 of the
 * paper]. Each worker takes from OPEN the safe state of smallest f, expands
 * it outside the lock and puts back the successors whose g has improved; a
 * state is expanded at most once, and when it is expanded its g is optimal
 * (with \f$ \varepsilon_p = 1 \f$). The search ends when the goal is safe to
 * expand, or when OPEN and BE are both empty. This is Algorithm 1 of the
 * paper; OPEN is an ordered set scanned in increasing f, and a worker that
 * finds no safe state waits until OPEN or BE change.
 *
 * Finding a safe state costs a scan of the states of OPEN with a smaller f
 * and of BE, all under one lock: PASE pays only when the expansions are
 * slow, as in motion planning with collision checks. */

template< PairwiseDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
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
  if( f_options.eps < 1 )
   f_options.eps = 1;
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
  Shared sh;
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh );

  {
   const State s0 = f_domain.initial_state();
   NodeT * root = w[ 0 ]->f_pool.allocate();
   *root = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) ,
                  nullptr , false , false };
   sh.table.find_or_insert( s0 , f_domain.hash( s0 ) ,
                            [ & ] { return( root ); } );
   sh.open.insert( Key( root->f() , - root->g , root ) , true );
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
    sol.path = trace_path( sh.goal );
    }
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Lock = typename Executor::lock_type;

 // OPEN: by f, then by larger g, then by address
 using Key = std::tuple< Cost , Cost , NodeT * >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers, all protected by lock

 struct Shared {
  Lock lock;                                    // protects everything
  std::condition_variable_any changed;          // OPEN or BE changed
  ADS::BTreeMap< Key , bool > open;            // OPEN
  std::vector< NodeT * > be;                    // being expanded
  NodeTable< NodeT > table;                     // all the states
  bool done = false;                            // the search is over
  bool limit_hit = false;                       // the node limit is hit
  NodeT * goal = nullptr;                       // the goal, once safe
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its pool of nodes and its counters

 struct alignas( cache_line ) Worker {

  Worker( const PASE & a , Shared & s ) : f_a( a ) , f_sh( s ) {}

  // the states of OPEN that a worker examines for a safe one
  static constexpr std::size_t c_candidates = 64;

  // true if s is independent of t: g( s ) - g( t ) <= eps h( t , s ), up
  // to a relative 1e-9 for real costs, whose rounding would otherwise make
  // the goal dependent on the states with f = f( goal )
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

  // true if s is safe to expand (lock held): independent of BE and of the
  // states of OPEN with a smaller f, which come first in OPEN, as in the
  // definition of safety of the paper; the states of equal f are not
  // checked, so that a layer of equal f costs a scan of BE per state
  // instead of a scan of the layer
  bool safe( const NodeT * s ) const {
   for( const NodeT * t : f_sh.be )
    if( ! independent( s , t ) )
     return( false );
   const Cost fs = s->f();
   for( const auto & k : f_sh.open ) {
    const NodeT * t = std::get< 2 >( k.first );
    if( ! ( t->f() < fs ) )
     break;
    if( ! independent( s , t ) )
     return( false );
    }
   return( true );
   }

  // Algorithm 1 of the paper
  void run( void ) {
   std::vector< std::tuple< State , Cost , Cost > > succ;
   std::unique_lock< Lock > guard( f_sh.lock );
   while( ! f_sh.done ) {
    // the safe state of smallest f among the first c_candidates of OPEN,
    // if any: further down, each candidate would have to be checked
    // against the larger and larger part of OPEN before it, under the lock;
    // with BE empty the first state is safe, and the search goes on
    NodeT * s = nullptr;
    std::size_t n = 0;
    for( const auto & k : f_sh.open ) {
     if( n++ == c_candidates )
      break;
     if( safe( std::get< 2 >( k.first ) ) ) {
      s = std::get< 2 >( k.first );
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
    f_sh.open.remove( Key( s->f() , - s->g , s ) );
    if( f_a.f_domain.is_goal( s->state ) ) {
     f_sh.goal = s;               // safe, hence of optimal g
     f_sh.done = true;
     f_sh.changed.notify_all();
     break;
     }
    f_sh.be.push_back( s );
    s->closed = s->expanded = true;
    const State st = s->state;
    const Cost g = s->g , h = s->h;
    guard.unlock();

    ++f_st.expanded;
    succ.clear();
    f_a.f_domain.for_each_successor( st , h , [ & ]( const State & t ,
                                                      Cost c , Cost ht ) {
     succ.emplace_back( t , g + c , ht );
     } );

    guard.lock();
    for( const auto & sc : succ ) {    // no structured binding: OpenMP
     const State & t = std::get< 0 >( sc );     // cannot capture one
     const Cost gt = std::get< 1 >( sc ) , ht = std::get< 2 >( sc );
     ++f_st.generated;
     auto [ x , fresh ] = f_sh.table.find_or_insert( t ,
                                                     f_a.f_domain.hash( t ) ,
                                                     [ & ] {
      NodeT * y = f_pool.allocate();
      *y = NodeT{ t , gt , ht , s , false , false };
      return( y );
      } );
     if( fresh )
      f_sh.open.insert( Key( x->f() , - x->g , x ) , true );
     else
      if( ! x->closed && ( gt < x->g ) ) {
       f_sh.open.remove( Key( x->f() , - x->g , x ) );
       x->g = gt;
       x->parent = s;
       f_sh.open.insert( Key( x->f() , - x->g , x ) , true );
       }
      else
       ++f_st.duplicates;
     }
    f_sh.be.erase( std::find( f_sh.be.begin() , f_sh.be.end() , s ) );
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
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 PASEOptions f_options;                                ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( PASE ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PASE */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File PASE.h -----------------------------*/
/*--------------------------------------------------------------------------*/
