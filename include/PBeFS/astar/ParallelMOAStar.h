/*--------------------------------------------------------------------------*/
/*------------------------- File ParallelMOAStar.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the parallel multi-objective searches on the workers of
 * an executor of ADS: the class template MOAStarPar, the searches on
 * different orders of the objectives that share their upper bounds, and the
 * class template HDMOAStar, the multi-objective search distributed by the
 * hash of the states.
 *
 * MOAStarPar is the parallel MOA* of Ahmadi et al. (ICAPS 2025): each
 * worker \f$ j \f$ runs MOA* [see MOAStar.h] with its own order of the
 * objectives, whose first one is its primary objective \f$ i_j \f$, and its
 * own OPEN, frontiers and solutions; the workers publish, for each
 * objective \f$ i \f$, the largest \f$ F_i \f$ of the \f$ f_i \f$ that a
 * worker of primary objective i has extracted, and the upper bounds
 * \f$ \bar f_i \f$. A worker keeps the costs of its solutions, in its
 * order, in the sets \f$ V_1 , \ldots , V_{K-1} \f$, \f$ V_\lambda \f$
 * holding vectors of which only the components from the \f$ \lambda \f$-th
 * on count, \f$ V_1 \f$ being the frontier S of MOA*:
 *
 * \code{.unparsed}
 *  1  procedure MOA*-par(s₀)
 *  2     F ← 0;  f̄ ← ∞
 *  3     in parallel, on each worker j with its order (i_j first) do
 *  4        MOA*(s₀) in the order of j, where
 *  5           extracting x:  F_{i_j} ← max(F_{i_j}, f_{i_j}(x));
 *                 if f_{i_j}(x) ≥ f̄_{i_j} then stop all the workers
 *  6           the solution check of x (and of a successor, if not lazy):
 *                 f_i(x) ≥ f̄_i for some i, or the components from λ on
 *                 of f(x) weakly dominated by a vector of V_λ, some λ
 *  7           a solution z:  add z to V₁;  UPDATE
 *  8     return the union of the solutions of the workers, each cost once
 *  9  procedure UPDATE                              ▷ on the worker j
 * 10     for λ = 1, ..., K − 2 do
 * 11        i ← the objective of the component λ in the order of j
 * 12        for each v ∈ V_λ with v_i ≤ F_i do
 * 13           move v from V_λ to V_{λ+1}, unless a vector there
 *                 dominates it, removing those that it dominates
 * 14     if V_{K−1} ≠ ∅ then f̄_i ← min(f̄_i, min V_{K−1}), i the last
 *              objective in the order of j
 * \endcode
 *
 * The workers are as many as the threads, in the cyclic orders of the
 * objectives first, as in the paper, and then in the other permutations, in
 * lexicographic order, up to \f$ K! \f$; the search ends as soon as one
 * worker ends, and its theorem, proven on the page "Multi-objective search"
 * of the documentation, is that the union is the Pareto front, each vector
 * found by at least one worker. UPDATE is run at every solution, as in the
 * paper, and also every flush_period extractions when F has grown since,
 * which makes the truncations of line 12 earlier.
 *
 * HDMOAStar is the hash-distributed MOA*, in which a function of the hash
 * of a state assigns it to the worker that owns it, as in HDA* [see
 * HDAStar.h]: the owner keeps the labels of its states in its OPEN, by the
 * lexicographic order of f, and their frontiers, and the paths generated
 * by a worker travel in batches to the owner of their last state through
 * the mailbox of ADS. Since the workers do not extract the labels in one
 * global order, neither the dimensionality reduction nor the permanence of
 * a label expanded holds: the frontier of a state keeps the full g of the
 * labels expanded there, a label is checked in full at its state and against
 * the solutions, and the owner of a goal keeps its solutions mutually non-
 * dominated, removing those that a later one dominates, and sends each new
 * one to every worker, which prunes with its copy of the solutions. This is
 * a label-correcting search, which ends, when no worker has a label and no
 * message is travelling, with the Pareto front [see the page
 * "Multi-objective search"]. Since a label expanded early may be dominated
 * by one that arrives later, a worker expands a label only if its
 * \f$ f_1 \f$ is within a window of the least \f$ f_1 \f$ of the tops of
 * all the OPEN lists, which every worker publishes [see
 * ParallelMOOptions::window]: this keeps the workers close to the
 * lexicographic order of the sequential search, and decides only the order
 * of the expansions, not the front.
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

#ifndef __PBeFS_ParallelMOAStar
 #define __PBeFS_ParallelMOAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <memory>

#include <numeric>

#include <vector>

#include "PBeFS/astar/MOAStar.h"

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
/*------------------------ CLASS ParallelMOOptions -------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the parallel multi-objective searches: those of the
/// sequential ones and the parameters of the parallel ones [see
/// ParallelParams]; order is ignored by MOAStarPar, whose workers have their
/// own orders

struct ParallelMOOptions : public MOSearchOptions , public ParallelParams {

 /// with false, MOAStarPar shares no bound among its workers (the variant
 /// without upper bounds of Ahmadi et al., ICAPS 2025)
 bool share_bounds = true;

 /// HDMOAStar expands a label only if its \f$ f_1 \f$ (in the order of the
 /// search) is at most \f$ L + \mathit{window} \f$, L being the least
 /// \f$ f_1 \f$ of the tops of the OPEN lists of the workers; negative for
 /// no such limit
 double window = 0;

 };  // end( struct( ParallelMOOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the first n orders of K objectives of MOAStarPar: the K cyclic ones, then
/// the others in lexicographic order, n at most \f$ K! \f$

template< std::size_t K >
[[nodiscard]] std::vector< std::array< unsigned , K > > objective_orders(
 std::size_t n )
{
 std::vector< std::array< unsigned , K > > o;
 for( std::size_t j = 0 ; ( j < K ) && ( o.size() < n ) ; ++j ) {
  std::array< unsigned , K > p;
  for( std::size_t i = 0 ; i < K ; ++i )
   p[ i ] = unsigned( ( j + i ) % K );
  o.push_back( p );
  }
 std::array< unsigned , K > p;
 std::iota( p.begin() , p.end() , 0u );
 do
  if( ( o.size() < n ) && ( std::find( o.begin() , o.end() , p ) ==
                            o.end() ) )
   o.push_back( p );
 while( ( o.size() < n ) && std::next_permutation( p.begin() , p.end() ) );
 return( o );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS MOAStarPar ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel MOA* on different orders of the objectives
/** The search of the file header, each worker with the frontier F of its
 * states (ArrayFrontier by default, LTMOA*-par; TreeFrontier for EMOA*,
 * ScalarFrontier for BOA* with two objectives, whose two workers are then
 * BOA* in the orders ( 1 , 2 ) and ( 2 , 1 ) sharing their bounds), the
 * options of ParallelMOOptions, threads being the number of the workers (at
 * most \f$ K! \f$), and the lazy checks if asked; OPEN is a Q, as in
 * MOAStar. The labels pruned when extracted are recycled. */

template< MultiObjectiveDomain D ,
          template< class , std::size_t > class F = AutoFrontier ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool ,
          template< class , class > class Q = MOHeap >
class MOAStarPar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr std::size_t K = objectives_of< D >;  ///< the objectives
 using State = mo_state_t< D >;                      ///< the state type
 using Cost = mo_scalar_t< D >;                      ///< a component
 using Vec = mo_cost_t< D >;                         ///< a vector of costs
 using Frontier = F< Cost , K >;                     ///< the frontier
 using SolutionT = pareto_solution_t< D >;           ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the MOAStarPar)

 explicit MOAStarPar( const D & domain ,
                      const ParallelMOOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Pareto front

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  std::size_t fact = 1;
  for( std::size_t i = 2 ; i <= K ; ++i )
   fact *= i;
  const auto orders = objective_orders< K >( std::min< std::size_t >(
                                         f_options.threads , fact ) );
  const unsigned p = unsigned( orders.size() );
  Shared sh;
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned j = 0 ; j < p ; ++j )
   w[ j ] = std::make_unique< Worker >( *this , sh , orders[ j ] );
  run_workers( f_ex , p , [ & ]( unsigned j ) { w[ j ]->run(); } );

  SolutionT sol;
  MOSearchStats & st = sol.stats;
  std::vector< const Label * > goals;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.generated += x->f_st.generated;
   st.extracted += x->f_st.extracted;
   st.inserted += x->f_st.inserted;
   st.checks += x->f_st.checks;
   st.solutions += x->f_st.solutions;
   st.stored += x->f_labels.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_labels.size() );
   goals.insert( goals.end() , x->v_goals.begin() , x->v_goals.end() );
   }

  // the union of the solutions, each cost once
  sol.status = sh.limit.load() ? Status::node_limit
                               : ( goals.empty() ? Status::unsolvable
                                                 : Status::solved );
  for( const Label * x : goals ) {
   typename SolutionT::Path q;
   q.cost = unpermute( x->g , *x->rec->order );
   if( f_options.paths )
    q.states = trace_labels( x );
   sol.front.push_back( std::move( q ) );
   }
  std::sort( sol.front.begin() , sol.front.end() ,
             []( const auto & a , const auto & b ) {
              return( a.cost < b.cost );
              } );
  sol.front.erase( std::unique( sol.front.begin() , sol.front.end() ,
                                []( const auto & a , const auto & b ) {
                                 return( a.cost == b.cost );
                                 } ) , sol.front.end() );
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 struct Rec;

 /// a path: the record of its last state, its g in the order of its
 /// worker and its predecessor
 struct Label {
  Rec * rec;
  Vec g;
  Label * parent;
  };

 /// the record of a state for a worker: the state, its h in the order of
 /// the worker, the order, whether it is a goal, its frontier
 struct Rec {
  State state;
  Vec h;
  const std::array< unsigned , K > * order;
  bool goal;
  Frontier front;
  };

 /// an entry of OPEN: f, the label and the order of the pushes
 struct Item {
  Vec f;
  Label * label;
  std::uint64_t seq;
  };

 /// the order of OPEN: the lexicographic order of f, then the first pushed
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.seq < b.seq ) ) );
   }
  };

 /// a value of the workers on a cache line of its own
 struct alignas( cache_line ) Bound {
  std::atomic< Cost > value;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers: F, the upper bounds, the end

 struct Shared {
  Shared( void ) {
   for( std::size_t i = 0 ; i < K ; ++i ) {
    explored[ i ].value.store( Cost( 0 ) );
    ub[ i ].value.store( infinite_cost< Cost >() );
    }
   }
  std::array< Bound , K > explored;        // F_i
  std::array< Bound , K > ub;              // the upper bounds
  alignas( cache_line ) std::atomic< bool > done{ false };
  std::atomic< bool > limit{ false };
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its order, OPEN, states, labels, solutions and counters

 struct alignas( cache_line ) Worker {

  Worker( const MOAStarPar & a , Shared & s ,
          const std::array< unsigned , K > & o )
   : f_a( a ) , f_sh( s ) , f_order( o ) , v_level( K ) {}

  // the search of the worker
  void run( void ) {
   const D & d = f_a.f_domain;
   const State s0 = d.initial_state();
   Rec * r0 = record( s0 , d.heuristic( s0 ) );
   if( ! dead_end( r0->h ) ) {
    f_open.push( Item{ r0->h , label( r0 , Vec{} , nullptr ) , f_seq++ } );
    ++f_st.inserted;
    }
   const unsigned i0 = f_order[ 0 ];
   unsigned since = 0;
   while( ! f_open.empty() &&
          ! f_sh.done.load( std::memory_order_relaxed ) ) {
    const Item it = f_open.pop();
    ++f_st.extracted;
    Label * x = it.label;
    Rec * r = x->rec;
    if( f_a.f_options.share_bounds ) {
     raise( f_sh.explored[ i0 ].value , it.f[ 0 ] );
     if( it.f[ 0 ] >= f_sh.ub[ i0 ].value.load(
                                         std::memory_order_acquire ) )
      break;                     // all the labels left are out of bounds
     if( ++since >= f_a.f_options.flush_period ) {
      since = 0;
      update();
      }
     }
    if( ( ! r->goal && r->front.dominated( x->g , f_st.checks ) ) ||
        solution_check( it.f ) ) {
     v_spare.push_back( x );
     continue;
     }
    if( r->goal ) {
     add_solution( it.f );
     v_goals.push_back( x );
     ++f_st.solutions;
     continue;
     }
    r->front.insert( x->g , f_st.checks );
    ++f_st.expanded;
    d.for_each_successor( r->state , unpermute( r->h , f_order ) ,
                          [ & ]( const State & t , const Vec & c ,
                                 const Vec & ht ) {
     ++f_st.generated;
     Rec * q = record( t , ht );
     if( dead_end( q->h ) )
      return;
     const Vec g = add_costs( x->g , permute( c , f_order ) );
     const Vec f = add_costs( g , q->h );
     if( ! f_a.f_options.lazy &&
         ( ( ! q->goal && q->front.dominated( g , f_st.checks ) ) ||
           solution_check( f ) ) )
      return;
     f_open.push( Item{ f , label( q , g , x ) , f_seq++ } );
     ++f_st.inserted;
     } );
    if( f_a.f_options.node_limit &&
        ( f_labels.size() >= f_a.f_options.node_limit ) ) {
     f_sh.limit.store( true );
     break;
     }
    }
   f_sh.done.store( true );
   }

  // F_i <- max( F_i , v ), released, so that a worker that reads it sees
  // what happened before [see the page "Multi-objective search"]
  static void raise( std::atomic< Cost > & a , Cost v ) {
   Cost old = a.load( std::memory_order_relaxed );
   while( ( old < v ) &&
          ! a.compare_exchange_weak( old , v , std::memory_order_release ,
                                     std::memory_order_relaxed ) )
    ;
   }

  // ub <- min( ub , v ), released as raise()
  static void lower( std::atomic< Cost > & a , Cost v ) {
   Cost old = a.load( std::memory_order_relaxed );
   while( ( v < old ) &&
          ! a.compare_exchange_weak( old , v , std::memory_order_release ,
                                     std::memory_order_relaxed ) )
    ;
   }

  // true if the solution check of line 6 prunes f
  bool solution_check( const Vec & f ) {
   if( f_a.f_options.share_bounds )
    for( std::size_t m = 0 ; m < K ; ++m ) {
     ++f_st.checks;
     if( f[ m ] >= f_sh.ub[ f_order[ m ] ].value.load(
                                           std::memory_order_acquire ) )
      return( true );
     }
   for( std::size_t l = K - 1 ; l >= 1 ; --l )
    for( const Vec & u : v_level[ l ] ) {
     ++f_st.checks;
     if( weakly_dominates( u , f , l ) )
      return( true );
     }
   return( false );
   }

  // adds v to V_l, unless it is dominated there, removing what it
  // dominates; false if dominated
  bool add_level( std::size_t l , const Vec & v ) {
   auto & V = v_level[ l ];
   for( const Vec & u : V ) {
    ++f_st.checks;
    if( weakly_dominates( u , v , l ) )
     return( false );
    }
   for( std::size_t k = 0 ; k < V.size() ; ) {
    ++f_st.checks;
    if( weakly_dominates( v , V[ k ] , l ) ) {
     V[ k ] = V.back();
     V.pop_back();
     }
    else
     ++k;
    }
   V.push_back( v );
   return( true );
   }

  // a solution of cost f (in the order of the worker): into V_1, then
  // UPDATE
  void add_solution( const Vec & f ) {
   add_level( 1 , f );
   if( f_a.f_options.share_bounds )
    update();
   }

  // UPDATE of the file header
  void update( void ) {
   for( std::size_t l = 1 ; l + 1 < K ; ++l ) {
    const Cost Fi = f_sh.explored[ f_order[ l ] ].value.load(
                                               std::memory_order_acquire );
    auto & V = v_level[ l ];
    for( std::size_t k = 0 ; k < V.size() ; )
     if( V[ k ][ l ] <= Fi ) {
      const Vec v = V[ k ];
      V[ k ] = V.back();
      V.pop_back();
      static_cast< void >( add_level( l + 1 , v ) );
      }
     else
      ++k;
    }
   const auto & last = v_level[ K - 1 ];
   if( ! last.empty() ) {
    Cost m = last[ 0 ][ K - 1 ];
    for( const Vec & v : last )
     m = std::min( m , v[ K - 1 ] );
    lower( f_sh.ub[ f_order[ K - 1 ] ].value , m );
    }
   }

  // the record of the state s, of heuristic h, made if new
  Rec * record( const State & s , const Vec & h ) {
   return( f_table.find_or_insert( s , f_a.f_domain.hash( s ) , [ & ] {
    Rec * r = f_recs.allocate();
    r->state = s;
    r->h = permute( h , f_order );
    r->order = & f_order;
    r->goal = f_a.f_domain.is_goal( s );
    return( r );
    } ).first );
   }

  // a label, recycled if possible
  Label * label( Rec * r , const Vec & g , Label * parent ) {
   Label * x;
   if( v_spare.empty() )
    x = f_labels.allocate();
   else {
    x = v_spare.back();
    v_spare.pop_back();
    }
   *x = Label{ r , g , parent };
   return( x );
   }

  const MOAStarPar & f_a;                       // the search
  Shared & f_sh;                                // the shared data
  std::array< unsigned , K > f_order;           // the order
  NodePool< Rec > f_recs;                       // the records
  NodeTable< Rec > f_table{ 1024 };             // the states
  NodePool< Label > f_labels;                   // the labels
  std::vector< Label * > v_spare;               // the labels to recycle
  Q< Item , Better > f_open;                    // OPEN
  std::uint64_t f_seq = 0;                      // the pushes
  std::vector< std::vector< Vec > > v_level;    // V_1 , ... , V_{K-1}
  std::vector< Label * > v_goals;               // the solutions
  MOSearchStats f_st;                           // the counters
  };

 const D & f_domain;                                   ///< the Domain
 ParallelMOOptions f_options;                          ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( MOAStarPar ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS HDMOAStar -----------------------------*/
/*--------------------------------------------------------------------------*/
/// the multi-objective search distributed by the hash of the states
/** The search of the file header: the owner of a state is
 * \f$ \mathrm{mix}( \mathrm{hash}( s ) \oplus \mathit{seed} ) \bmod p \f$,
 * the paths travel in batches of batch messages, an outbox not yet full is
 * flushed every flush_period extractions and before the worker waits, and
 * each worker orders its OPEN by the lexicographic order of f in the order
 * of the options. With lazy, the paths are checked only when extracted.
 * The labels of the parents cross the workers, and the paths are followed
 * only after all the workers have finished. */

template< MultiObjectiveDomain D ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
class HDMOAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 static constexpr std::size_t K = objectives_of< D >;  ///< the objectives
 using State = mo_state_t< D >;                      ///< the state type
 using Cost = mo_scalar_t< D >;                      ///< a component
 using Vec = mo_cost_t< D >;                         ///< a vector of costs
 using SolutionT = pareto_solution_t< D >;           ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the HDMOAStar);
 /// throws std::invalid_argument if the order is not a permutation

 explicit HDMOAStar( const D & domain ,
                     const ParallelMOOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_order( objective_order< K >( options.order ) ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Pareto front

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  const unsigned p = f_options.threads;
  Shared sh( p , f_options.batch );
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );
  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  MOSearchStats & st = sol.stats;
  std::vector< const Label * > goals;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.extracted += x->f_st.extracted;
   st.inserted += x->f_st.inserted;
   st.checks += x->f_st.checks;
   st.solutions += x->f_st.solutions;
   st.sent += x->f_st.sent;
   st.messages += sh.mail.batches( x->f_id );
   st.stored += x->f_labels.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_labels.size() );
   goals.insert( goals.end() , x->v_goals.begin() , x->v_goals.end() );
   }

  // the solutions of the owners of the goals, which may dominate one
  // another if there are several goals, filtered
  std::sort( goals.begin() , goals.end() ,
             []( const Label * a , const Label * b ) {
              return( a->g < b->g );
              } );
  std::vector< const Label * > front;
  for( const Label * x : goals )
   if( std::none_of( front.begin() , front.end() ,
                     [ & ]( const Label * y ) {
                      return( weakly_dominates( y->g , x->g ) );
                      } ) )
    front.push_back( x );
  sol.status = sh.limit.load() ? Status::node_limit
                               : ( front.empty() ? Status::unsolvable
                                                 : Status::solved );
  for( const Label * x : front ) {
   typename SolutionT::Path q;
   q.cost = unpermute( x->g , f_order );
   if( f_options.paths )
    q.states = trace_labels( x );
   sol.front.push_back( std::move( q ) );
   }
  std::sort( sol.front.begin() , sol.front.end() ,
             []( const auto & a , const auto & b ) {
              return( a.cost < b.cost );
              } );
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 struct Rec;

 /// ( 0 , v ), whose truncation is v
 [[nodiscard]] static std::array< Cost , K + 1 > lift( const Vec & v ) {
  std::array< Cost , K + 1 > w;
  w[ 0 ] = Cost( 0 );
  std::copy( v.begin() , v.end() , w.begin() + 1 );
  return( w );
  }

 /// a set of vectors, kept twice: in full, for the checks of any vector,
 /// in the tree of EMOA* of the vectors ( 0 , v ) for two objectives, where
 /// the vectors form a staircase and a check takes O( log n ), in a
 /// LayerIndex for three, in a BlockIndex for more, and truncated, in the
 /// frontier of EMOA* of the vectors v, which answers exactly for a vector
 /// whose first component is not smaller than that of any vector added [see
 /// the page "Multi-objective search"], the case of the labels that arrive
 /// in the order of f_1; add( v ) needs v not dominated
 struct VecSet {
  std::conditional_t< K == 2 , TreeFrontier< Cost , K + 1 > ,
   std::conditional_t< K == 3 , LayerIndex< Cost , K > ,
                       BlockIndex< Cost , K > > > full;  // the vectors
  EMOAFrontier< Cost , K > tr;               // their truncations
  Cost top{};                                 // the largest v_1 added
  bool any = false;                           // a vector was added
  [[nodiscard]] bool dominated( const Vec & v , std::uint64_t & n ) const {
   if( ! any )
    return( false );
   if( v[ 0 ] >= top )
    return( tr.dominated( v , n ) );
   if constexpr( K == 2 )
    return( full.dominated( lift( v ) , n ) );
   else
    return( full.dominated( v , n ) );
   }
  void add( const Vec & v , std::uint64_t & n ) {
   if constexpr( K == 2 )
    full.insert( lift( v ) , n );
   else
    full.insert( v , n );
   if( ! tr.dominated( v , n ) )
    tr.insert( v , n );
   if( ! any || ( top < v[ 0 ] ) )
    top = v[ 0 ];
   any = true;
   }
  };

 /// a path: the record of its last state, its g and its predecessor
 struct Label {
  Rec * rec;
  Vec g;
  Label * parent;
  };

 /// the record of a state at its owner: the state, its h in the order of
 /// the search, whether it is a goal, the g of the labels expanded there
 struct Rec {
  State state;
  Vec h;
  bool goal;
  VecSet closed;
  };

 /// a message: a path to its owner (s, g in the order of the search, h as
 /// the Domain gives it, parent), or the cost f of a solution to every
 /// worker
 struct Msg {
  State s;
  Vec g;
  Vec h;
  Label * parent;
  bool solution;
  };

 /// an entry of OPEN: f, the label and the order of the pushes
 struct Item {
  Vec f;
  Label * label;
  std::uint64_t seq;
  };

 /// the order of OPEN: the lexicographic order of f, then the first pushed
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.f < b.f ) || ( ( a.f == b.f ) && ( a.seq < b.seq ) ) );
   }
  };

 using Mailbox = ADS::Mailbox< Msg , ADS::DefaultChannels< Msg , Executor > >;

 /// the f_1 of the top of the OPEN of a worker, on a cache line of its own
 struct alignas( cache_line ) Low {
  std::atomic< Cost > value{ infinite_cost< Cost >() };
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data shared by the workers

 struct Shared {
  Shared( unsigned np , std::size_t batch )
   : p( np ) , mail( np , batch ) , low( np ) {}
  unsigned p;                                   // the number of workers
  Mailbox mail;                                 // the messages
  std::vector< Low > low;                       // the f_1 of their tops
  std::atomic< bool > limit{ false };
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its OPEN, its states, its labels, its copy of the solutions

 struct alignas( cache_line ) Worker {

  Worker( const HDMOAStar & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) {}

  // the main loop
  void run( void ) {
   const D & d = f_a.f_domain;
   const State s0 = d.initial_state();
   if( owner( s0 ) == f_id )
    insert( Msg{ s0 , Vec{} , d.heuristic( s0 ) , nullptr , false } );
   unsigned since = 0;
   ADS::Backoff backoff;
   while( ! f_sh.mail.stopped() ) {
    v_in.clear();
    f_sh.mail.poll( f_id , [ & ]( Msg && m ) {
     v_in.push_back( std::move( m ) );
     } );
    for( const Msg & m : v_in )
     if( m.solution )
      static_cast< void >( add_vector( v_known , m.g ) );
     else
      insert( m );
    const Cost top = f_open.empty() ? infinite_cost< Cost >()
                                    : f_open.top().f[ 0 ];
    // the least f_1 of the labels of the worker, in OPEN or in its outboxes
    f_sh.low[ f_id ].value.store( std::min( top , f_pending ) ,
                                  std::memory_order_relaxed );
    if( ! f_open.empty() ) {
     if( ! in_window( top ) ) {
      // a label of a smaller f_1 may still dominate this one: wait
      flush();
      backoff.pause();
      continue;
      }
     backoff.reset();
     extract();
     if( ++since >= f_a.f_options.flush_period ) {
      since = 0;
      flush();
      }
     if( f_a.f_options.node_limit &&
         ( f_labels.size() >= f_a.f_options.node_limit ) ) {
      f_sh.limit.store( true );
      f_sh.mail.stop();
      }
     continue;
     }
    flush();
    f_sh.low[ f_id ].value.store( f_pending , std::memory_order_relaxed );
    if( f_sh.mail.quiescent( f_id ) )
     return;
    backoff.pause();
    }
   }

  // the outboxes sent; f_pending forgotten if they all left
  void flush( void ) {
   if( f_sh.mail.flush( f_id ) )
    f_pending = infinite_cost< Cost >();
   }

  // true if a label of f_1 = top may be expanded [see window]
  [[nodiscard]] bool in_window( Cost top ) const {
   const double w = f_a.f_options.window;
   if( w < 0 )
    return( true );
   Cost L = top;
   for( const Low & l : f_sh.low )
    L = std::min( L , l.value.load( std::memory_order_relaxed ) );
   return( double( top ) <= double( L ) + w );
   }

  // the owner of s
  [[nodiscard]] unsigned owner( const State & s ) const {
   return( unsigned( ADS::splitmix64( f_a.f_domain.hash( s ) ^
                                      f_a.f_options.seed ) % f_sh.p ) );
   }

  // true if a vector of V weakly dominates v
  bool covered( const VecSet & V , const Vec & v ) {
   return( V.dominated( v , f_st.checks ) );
   }

  // adds v to V, removing the vectors it dominates, unless one dominates
  // it; false then
  bool add_vector( VecSet & V , const Vec & v ) {
   if( V.dominated( v , f_st.checks ) )
    return( false );
   V.add( v , f_st.checks );
   return( true );
   }

  // a path that has reached its owner: into OPEN, unless pruned
  void insert( const Msg & m ) {
   Rec * r = f_table.find_or_insert( m.s , f_a.f_domain.hash( m.s ) , [ & ] {
    Rec * q = f_recs.allocate();
    q->state = m.s;
    q->h = permute( m.h , f_a.f_order );
    q->goal = f_a.f_domain.is_goal( m.s );
    return( q );
    } ).first;
   if( dead_end( r->h ) )
    return;
   const Vec & g = m.g;
   const Vec f = add_costs( g , r->h );
   if( ! f_a.f_options.lazy && ( ( ! r->goal && covered( r->closed , g ) ) ||
                                 covered( v_known , f ) ) ) {
    ++f_st.duplicates;
    return;
    }
   Label * x;
   if( v_spare.empty() )
    x = f_labels.allocate();
   else {
    x = v_spare.back();
    v_spare.pop_back();
    }
   *x = Label{ r , g , m.parent };
   f_open.push( Item{ f , x , f_seq++ } );
   ++f_st.inserted;
   }

  // the first label of OPEN: pruned, a solution, or expanded
  void extract( void ) {
   const Item it = f_open.pop();
   ++f_st.extracted;
   Label * x = it.label;
   Rec * r = x->rec;
   if( ( ! r->goal && covered( r->closed , x->g ) ) ||
       covered( v_known , it.f ) ) {
    v_spare.push_back( x );
    return;
    }
   if( r->goal ) {
    v_known.add( it.f , f_st.checks );                 // checked above
    for( std::size_t k = 0 ; k < v_goals.size() ; ) {
     ++f_st.checks;
     if( weakly_dominates( x->g , v_goals[ k ]->g ) ) {
      v_goals[ k ] = v_goals.back();
      v_goals.pop_back();
      }
     else
      ++k;
     }
    v_goals.push_back( x );
    ++f_st.solutions;
    for( unsigned j = 0 ; j < f_sh.p ; ++j )
     if( j != f_id ) {
      f_sh.mail.send( f_id , j , Msg{ r->state , it.f , Vec{} , nullptr ,
                                      true } );
      ++f_st.sent;
      }
    return;
    }
   r->closed.add( x->g , f_st.checks );               // checked above
   ++f_st.expanded;
   const auto & order = f_a.f_order;
   f_a.f_domain.for_each_successor( r->state , unpermute( r->h , order ) ,
                                    [ & ]( const State & t , const Vec & c ,
                                           const Vec & ht ) {
    ++f_st.generated;
    if( dead_end( ht ) )
     return;
    const Vec g = add_costs( x->g , permute( c , order ) );
    if( ! f_a.f_options.lazy &&
        covered( v_known , add_costs( g , permute( ht , order ) ) ) )
     return;
    const Msg m{ t , g , ht , x , false };
    const unsigned dst = owner( t );
    if( dst == f_id )
     insert( m );
    else {
     f_sh.mail.send( f_id , dst , m );
     ++f_st.sent;
     f_pending = std::min( f_pending , Cost( g[ 0 ] + ht[ order[ 0 ] ] ) );
     }
    } );
   }

  const HDMOAStar & f_a;                        // the search
  Shared & f_sh;                                // the shared data
  unsigned f_id;                                // the index of the worker
  NodePool< Rec > f_recs;                       // the records
  NodeTable< Rec > f_table{ 1024 };             // its states
  NodePool< Label > f_labels;                   // its labels
  std::vector< Label * > v_spare;               // the labels to recycle
  MOHeap< Item , Better > f_open;               // its OPEN
  std::uint64_t f_seq = 0;                      // the pushes
  Cost f_pending = infinite_cost< Cost >();     // the least f_1 sent and
                                                // still in its outboxes
  VecSet v_known;                               // the solutions known
  std::vector< Label * > v_goals;               // its solutions
  std::vector< Msg > v_in;                      // the last inbox
  MOSearchStats f_st;                           // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 ParallelMOOptions f_options;                          ///< the options
 std::array< unsigned , K > f_order;                   ///< the order
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( HDMOAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ParallelMOAStar */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ParallelMOAStar.h ------------------------*/
/*--------------------------------------------------------------------------*/
