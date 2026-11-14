/*--------------------------------------------------------------------------*/
/*------------------------------- File TDS.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template TDS, the Transposition-table Driven
 * Scheduling of Romein, Plaat, Bal and Schaeffer (AAAI 1999): a parallel
 * IDA* in which every node is sent to the worker that owns the entry of
 * its state in the transposition table, which looks it up there and, if
 * the state is not already being searched with a bound at least as large,
 * pushes it on its own stack.
 *
 * With the notation of the page "Best-first search and A*" of the
 * documentation, \f$ \beta \f$ the bound of the iteration, \f$ T_i \f$ the
 * part of the table of the worker \f$ i \f$, each entry with the g with
 * which its state was last pushed in the iteration, and a message
 * \f$ ( v , g' , u ) \f$ a state with the cost and the predecessor of its
 * path:
 *
 * \code{.unparsed}
 *  1  procedure TDS(s₀, p, owner)
 *  2     β ← h(s₀)
 *  3     loop
 *  4        β' ← ∞;  INSERT(s₀, 0, nil) on owner(s₀)
 *  5        in parallel, on each worker i, until no worker has work and
 *              no message is in transit, or a goal is found, do
 *  6           for each message (v, g', u) received do INSERT(v, g', u)
 *  7           if the stack of i is not empty then
 *  8              u ← pop the stack (stale entries dropped)
 *  9              if u is a goal then record u;  stop all
 * 10              for each successor v of u, g' ← g(u) + c(u, v) do
 * 11                 if g' + h(v) > β then β' ← min(β', g' + h(v))
 * 12                 else if owner(v) = i then INSERT(v, g', u)
 * 13                 else send (v, g', u) to owner(v)
 * 14        if a goal was recorded then return the path to it
 * 15        if β' = ∞ then return unsolvable
 * 16        β ← β'
 * 17  procedure INSERT(v, g', u)                       ▷ on the owner i of v
 * 18     if T_i has v from this iteration with g(v) ≤ g' then return
 * 19     g(v) ← g';  parent(v) ← u;  push v on the stack of i
 * \endcode
 *
 * Line 18 is what TDS adds to IDA*: since the bound is the same for the
 * whole iteration, a state already pushed in it with a g not larger has
 * been, or is being, searched with a remaining budget
 * \f$ \beta - g \f$ at least as large, and its subtree is not searched again
 * by any worker (Romein et al., p. 3). The goal ends the search at once: in
 * the first iteration whose \f$ \beta \f$ is not below the optimum
 * \f$ C^* \f$, \f$ \beta = C^* \f$, the bounds being the least f pruned, so
 * that a goal of the iteration has \f$ g \le C^* \f$. Each iteration is a
 * computation of the mailbox of ADS, which detects its end by the deficit
 * counting of HDAStar [see HDAStar.h]; the least f pruned is the least
 * over the workers, read after the iteration.
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

#ifndef __PBeFS_TDS
 #define __PBeFS_TDS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <memory>

#include <mutex>

#include <optional>

#include <vector>

#include "PBeFS/astar/Deepening.h"

#include "PBeFS/astar/HDAStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS TDSOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of TDS: those of the parallel searches, and the bounds of
/// IDA*-CR [see Deepening.h]

struct TDSOptions : public ParallelOptions , public DeepeningOptions {

 TDSOptions( void ) = default;

 /// the options of the parallel searches, with the classical bounds
 TDSOptions( const ParallelOptions & o ) : ParallelOptions( o ) {}

 };  // end( struct( TDSOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS TDS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Transposition-table Driven Scheduling, a parallel IDA*
/** The transposition table is split among the workers by the distribution
 * Dist, as the CLOSED list of HDAStar, and each worker searches depth first
 * from its own stack; the nodes travel in batches, as in HDAStar, with the
 * outboxes flushed every flush_period expansions and before the worker
 * waits. An entry of the table is a node, kept from an iteration to the
 * next and marked with the iteration that last pushed it, so that the
 * table holds every state generated within the bounds, as the CLOSED list
 * of A*, and the node_limit (per worker) bounds it; with no replacement of
 * the entries, no subtree is searched twice in an iteration, which is the
 * case of a table large enough in Romein et al. (1999).
 *
 * The parents cross the workers, and the path is followed only after all
 * the workers have finished; its cost is that of its edges, recomputed
 * from the Domain, since the parent of a node on it may have been improved
 * in the same iteration by a cheaper path, which only lowers the cost. The
 * counters of the Solution are summed over the iterations. With
 * TDSOptions::growth above 1 the bounds are those of IDA*-CR, and an
 * iteration whose bound is above the least f pruned before does not stop
 * at the first goal if TDSOptions::optimal: the goals lower an incumbent,
 * which prunes the nodes of f not below it, and the iteration ends when no
 * worker has work, with an optimal path [see Deepening.h]. */

template< Domain D , class Dist = StateHashDistribution< D > ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
class TDS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain and the distribution (which must
 /// outlive the TDS)

 TDS( const D & domain , const Dist & dist ,
      const TDSOptions & options = {} )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers, for all its iterations

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , i );

  const State s0 = f_domain.initial_state();
  const Cost h0 = Cost( f_domain.heuristic( s0 ) );
  const unsigned root = f_dist.owner( s0 , p );
  std::optional< Executor > own;
  Executor * ex = f_ex ? f_ex : & own.emplace( p );
  Shared sh;
  SolutionT sol;
  const bool cr = f_options.growth > 1;
  Cost low = h0;                       // the least f pruned before
  std::uint64_t done = 0;              // the expansions before
  CutoffHistogram< Cost > cut;         // the f pruned, and the step
  for( sh.bound = h0 ; ; ) {
   ++sh.iteration;
   sh.bnb = cr && f_options.optimal && ( sh.bound > low );
   sh.z.store( infinite_cost< Cost >() );
   cut.reset( sh.bound , cr );
   for( auto & x : w )
    x->f_cut.reset( sh.bound , cr );
   Mailbox mail( p , f_options.batch );
   sh.mail = & mail;
   const Msg r{ s0 , Cost( 0 ) , h0 , nullptr , no_op ,
                distribution_key< D >( f_dist , s0 ) };
   run_workers( ex , p , [ & ]( unsigned i ) {
    w[ i ]->run( sh , i == root ? & r : nullptr );
    } );
   std::uint64_t all = 0;
   for( auto & x : w ) {
    cut.merge( x->f_cut );
    all += x->f_st.expanded;
    }
   if( sh.limit_hit.load() ) {
    sol.status = Status::node_limit;
    break;
    }
   if( sh.goal ) {
    sol.status = Status::solved;
    sol.path = trace_path( static_cast< const Node< State , Cost > * >(
                            sh.goal ) );
    sol.cost = path_cost( sol.path );
    break;
    }
   low = cut.low();
   if( low == infinite_cost< Cost >() )
    break;
   sh.bound = cut.next( f_options.growth , all - done );
   done = all;
   }

  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.reexpanded += x->f_st.reexpanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.sent += x->f_st.sent;
   st.messages += x->f_st.messages;
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // an entry of the table: a node, with the iteration that last pushed it

 struct TNode : public Node< State , Cost > {
  unsigned iteration = 0;
  };

 // a node travelling to its owner
 struct Msg {
  State s;
  Cost g;
  Cost h;
  TNode * parent;
  [[no_unique_address]] OperatorOf< D > o;        // the operator of s
  std::uint32_t key;                               // the key of s
  };

 using Mailbox = ADS::Mailbox< Msg , ADS::DefaultChannels< Msg , Executor > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the data of an iteration shared by the workers

 struct Shared {
  Mailbox * mail = nullptr;                     // the messages
  unsigned iteration = 0;                       // the iteration
  Cost bound{};                                 // its bound
  bool bnb = false;                             // goes on at a goal
  std::atomic< Cost > z{ infinite_cost< Cost >() };  // the incumbent if bnb
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;                        // protects goal
  TNode * goal = nullptr;                       // the goal found
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a worker: its part of the table, its stack and its counters

 struct alignas( cache_line ) Worker {

  Worker( const TDS & a , unsigned i )
   : f_a( a ) , f_id( i ) , f_dist( a.f_dist ) {}

  // one iteration, inserting the root r if it is given
  void run( Shared & sh , const Msg * r ) {
   f_sh = & sh;
   v_stack.clear();
   Mailbox & mail = *sh.mail;
   if( r )
    insert( *r , f_a.f_domain.hash( r->s ) );
   unsigned since_flush = 0;
   ADS::Backoff backoff;
   while( ! mail.stopped() ) {
    mail.poll( f_id , [ & ]( Msg && m ) {
     insert( m , f_a.f_domain.hash( m.s ) );
     } );
    if( TNode * n = pop() ) {
     backoff.reset();
     expand( n );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      mail.flush( f_id );
      }
     if( f_a.f_options.node_limit &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      sh.limit_hit.store( true );
      mail.stop();
      }
     continue;
     }
    if( mail.quiescent( f_id ) )
     break;
    backoff.pause();
    }
   f_st.messages += mail.batches( f_id );
   }

  // the top of the stack, the stale entries dropped; nullptr if none
  TNode * pop( void ) {
   while( ! v_stack.empty() ) {
    const auto [ n , g ] = v_stack.back();
    v_stack.pop_back();
    if( ( n->g == g ) &&
        ( n->g + n->h < f_sh->z.load( std::memory_order_relaxed ) ) )
     return( n );
    }
   return( nullptr );
   }

  // expands n, or records it and stops all if it is a goal
  void expand( TNode * n ) {
   ++f_st.expanded;
   if( n->expanded )
    ++f_st.reexpanded;
   n->expanded = true;
   if( f_a.f_domain.is_goal( n->state ) ) {
    std::lock_guard< std::mutex > lock( f_sh->goal_mutex );
    if( ! f_sh->bnb ) {
     if( ! f_sh->goal )
      f_sh->goal = n;
     f_sh->mail->stop();
     }
    else
     if( n->g < f_sh->z.load() ) {
      f_sh->z.store( n->g );
      f_sh->goal = n;
      }
    return;
    }
   const auto owner = successor_owners< D >(
    f_dist , n->state , n->key , unsigned( f_a.f_options.threads ) );
   const Cost bound = f_sh->bound;
   for_each_child( f_a.f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t ot ) {
    ++f_st.generated;
    const Cost g = n->g + c;
    if( g + ht > bound ) {
     f_cut.add( g + ht );
     return;
     }
    if( g + ht >= f_sh->z.load( std::memory_order_relaxed ) )
     return;
    const auto [ dst , key ] = owner( t );
    if( dst == f_id ) {
     insert( Msg{ t , g , ht , n , ot , key } , f_a.f_domain.hash( t ) );
     return;
     }
    f_sh->mail->send( f_id , dst , Msg{ t , g , ht , n , ot , key } );
    ++f_st.sent;
    } );
   }

  // INSERT: pushes the node of m, whose state has hash h, unless its state
  // was pushed in this iteration with a g not larger
  void insert( const Msg & m , std::uint64_t h ) {
   const unsigned it = f_sh->iteration;
   auto [ x , fresh ] = f_table.find_or_insert( m.s , h , [ & ] {
    TNode * y = f_pool.allocate();
    y->state = m.s;
    y->g = m.g;
    y->h = m.h;
    y->parent = m.parent;
    y->closed = y->expanded = false;
    y->op = m.o.op;
    y->key = m.key;
    y->iteration = it;
    return( y );
    } );
   if( ! fresh ) {
    if( ( x->iteration == it ) && ( x->g <= m.g ) ) {
     ++f_st.duplicates;
     return;
     }
    x->g = m.g;
    x->parent = m.parent;
    x->op = m.o.op;
    x->iteration = it;
    }
   v_stack.emplace_back( x , m.g );
   }

  const TDS & f_a;                              // the search
  unsigned f_id;                                // the index of the worker
  Dist f_dist;                                  // its distribution
  Shared * f_sh = nullptr;                      // the iteration
  NodePool< TNode > f_pool;                     // its entries
  NodeTable< TNode > f_table;                   // its part of the table
  std::vector< std::pair< TNode * , Cost > > v_stack;  // its stack
  CutoffHistogram< Cost > f_cut;                // the f above the bound
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the cost of the edges of a path, each the cheapest from a state to the
 // next one

 Cost path_cost( const std::vector< State > & path ) const {
  Cost sum = Cost( 0 );
  for( std::size_t i = 1 ; i < path.size() ; ++i ) {
   Cost best = infinite_cost< Cost >();
   f_domain.for_each_successor( path[ i - 1 ] ,
                                Cost( f_domain.heuristic( path[ i - 1 ] ) ) ,
                                [ & ]( const State & t , Cost c , Cost ) {
    if( ( t == path[ i ] ) && ( c < best ) )
     best = c;
    } );
   sum += best;
   }
  return( sum );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                                  ///< the distribution
 TDSOptions f_options;                                 ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( TDS ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_TDS */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File TDS.h -----------------------------*/
/*--------------------------------------------------------------------------*/
