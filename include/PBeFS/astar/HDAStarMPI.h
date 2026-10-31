/*--------------------------------------------------------------------------*/
/*--------------------------- File HDAStarMPI.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template HDAStarMPI, Hash Distributed A* on a
 * distributed-memory machine, with one MPI process per worker.
 *
 * The search is HDA* [see HDAStar.h] on the process \f$ i \f$ of
 * \f$ p \f$, with the notation of the page "Best-first search and A*" of the
 * documentation; a node records the state of its parent instead of a
 * pointer, since the parent may be stored by another process, and
 * \f$ z \f$ is the incumbent known to the process:
 *
 * \code{.unparsed}
 *  1  procedure HDA*-MPI(s₀) on the process i
 *  2     z ← ∞;  if owner(s₀) = i then INSERT(s₀, 0, nil)
 *  3     loop
 *  4        for each message m received do
 *  5           if m is a batch of nodes then INSERT each of them
 *  6           else if m is an incumbent z' then z ← min(z, z')
 *  7           else if m is a stop then stop expanding
 *  8        if OPEN_i has a node of f < z then
 *  9           expand it as at lines 7 to 15 of HDA*, a goal that lowers z
 *                 sending z to all the other processes
 * 10        else send all the outboxes
 * 11        if OPEN_i has no node of f < z and the mailbox detects the end
 * 12           then break
 * 13     the sends in progress complete
 * 14     the process of the least goal broadcasts it, and the path is rebuilt
 *           from it, each state broadcast with its parent by its owner
 * \endcode
 *
 * The nodes and the incumbents travel through the mailbox of the processes
 * of ADS, packed in batches of batch messages and sent with non-blocking
 * sends, and the end is detected by the four-counter method of Mattern,
 * which by line 11 a process joins only when it has no node of
 * \f$ f < z \f$ and nothing to send. A process that reaches the node limit
 * stops the mailbox of all the processes (line 7), after which they only
 * wait for the end. The path is
 * rebuilt by asking the owner of each state for its parent (line 14), which
 * is why the distribution has to be deterministic.
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

#ifndef __PBeFS_HDAStarMPI
 #define __PBeFS_HDAStarMPI
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <stdexcept>

#include <type_traits>

#include <variant>

#include <vector>

#include <mpi.h>

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/core/SearchUtils.h"

#include "ADS/distributed/MPIMailbox.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS HDAStarMPI ----------------------------*/
/*--------------------------------------------------------------------------*/
/// Hash Distributed A* with one MPI process per worker
/** The original HDA* of Kishimoto, Fukunaga and Botea (AIJ 2013): each
 * process owns the states that Dist assigns to it, keeps their nodes in its
 * own memory, and sends the nodes it generates for the other processes in
 * batches with non-blocking sends, checking for incoming batches between
 * two expansions. State and Cost have to be trivially copyable, since the
 * nodes travel as bytes, and Dist deterministic, since the path is rebuilt
 * by asking the owner of each state for its parent.
 *
 * A node records the state of its parent instead of a pointer. A process
 * that improves the incumbent sends its cost to all the others. The end is
 * detected with the four-counter method of Mattern by the mailbox of ADS, a
 * process joining a wave only when it has no node with f below the
 * incumbent and no batch to send. If a process reaches the node limit, it
 * tells the others, and all of them stop expanding and only drain the
 * messages until the end is detected.
 *
 * MPI must be initialized by the caller, with at least MPI_THREAD_FUNNELED
 * if the Domain uses threads; all the processes of the communicator call
 * solve(), and all of them get the same Solution. */

template< Domain D , class Dist = StateHashDistribution< D > >
class HDAStarMPI {

 static_assert( std::is_trivially_copyable_v< state_t< D > > ,
                "HDAStarMPI: the states have to be trivially copyable" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain, the distribution and the communicator

 HDAStarMPI( const D & domain , const Dist & dist ,
             const ParallelOptions & options = {} ,
             MPI_Comm comm = MPI_COMM_WORLD )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) ,
    f_comm( comm ) {
  f_options.normalize();
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution (on every process)

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  Process w( *this );
  w.run();
  SolutionT sol;
  w.gather( sol );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // a node: the parent is a state, owned by another process maybe
 struct MNode {
  State state;
  Cost g;
  Cost h;
  State parent;
  bool has_parent;
  bool closed;
  bool expanded;
  };

 // a node travelling to its owner
 struct Msg {
  State s;
  Cost g;
  Cost h;
  State parent;
  bool has_parent;
  };

 // a cost of a goal, sent to every process
 struct Incumbent {
  Cost c;
  };

 // the messages: the nodes and the incumbents
 using Mail = std::variant< Msg , Incumbent >;

 // the OPEN list: buckets for integer costs, the heap otherwise
 using Open_t = std::conditional_t< std::is_integral_v< Cost > ,
                                    BucketOpenList< MNode > ,
                                    HeapOpenList< MNode > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the search on one process

 struct Process {

  explicit Process( const HDAStarMPI & a )
   : f_a( a ) , f_comm( a.f_comm ) ,
     f_mail( a.f_comm , a.f_options.batch ) , f_dist( a.f_dist ) {
   MPI_Comm_rank( f_comm , & f_rank );
   int np;
   MPI_Comm_size( f_comm , & np );
   f_p = unsigned( np );
   }

  // the main loop
  void run( void ) {
   const State s0 = f_a.f_domain.initial_state();
   if( f_dist.owner( s0 , f_p ) == unsigned( f_rank ) )
    insert( Msg{ s0 , Cost( 0 ) , Cost( f_a.f_domain.heuristic( s0 ) ) , s0 ,
                 false } , f_a.f_domain.hash( s0 ) );

   unsigned since_flush = 0;
   for( ;; ) {
    receive();
    if( MNode * n = next() ) {
     expand( n );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      f_mail.flush();
      }
     if( f_a.f_options.node_limit && ! f_mail.stopped() &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      f_limit_hit = true;
      f_mail.stop();
      }
     continue;
     }
    if( f_mail.quiescent() )      // nothing to expand nor to send
     break;
    }
   }

  // the best node of OPEN with f < incumbent, not removed; nullptr if none
  MNode * next_peek( void ) {
   if( f_mail.stopped() )
    return( nullptr );
   while( ! f_open.empty() ) {
    const OpenEntry< MNode > e = f_open.top();
    if( e.node->closed || ( e.node->g != e.g ) ) {
     f_open.pop();
     continue;
     }
    if( e.node->g + e.node->h >= f_incumbent )
     return( nullptr );
    return( e.node );
    }
   return( nullptr );
   }

  // the best node of OPEN with f < incumbent, removed; nullptr if none
  MNode * next( void ) {
   MNode * n = next_peek();
   if( n )
    f_open.pop();
   return( n );
   }

  // expands n, or records it as the incumbent if it is a goal
  void expand( MNode * n ) {
   n->closed = true;
   ++f_st.expanded;
   if( n->expanded )
    ++f_st.reexpanded;
   n->expanded = true;

   if( f_a.f_domain.is_goal( n->state ) ) {
    if( n->g < f_incumbent ) {
     f_incumbent = n->g;
     f_goal = n;
     for( unsigned q = 0 ; q < f_p ; ++q )
      if( q != unsigned( f_rank ) )
       f_mail.send( q , Mail( Incumbent{ n->g } ) );
     f_mail.flush();
     }
    return;
    }

   f_a.f_domain.for_each_successor( n->state , n->h ,
                                    [ & ]( const State & t , Cost c ,
                                           Cost ht ) {
    ++f_st.generated;
    const Cost g = n->g + c;
    if( g + ht >= f_incumbent )
     return;
    const Msg m{ t , g , ht , n->state , true };
    const unsigned dst = f_dist.owner( t , f_p );
    if( dst == unsigned( f_rank ) ) {
     insert( m , f_a.f_domain.hash( t ) );
     return;
     }
    f_mail.send( dst , Mail( m ) );
    ++f_st.sent;
    } );
   }

  // adds the node of m, whose state has hash h, to OPEN, unless it is a
  // duplicate
  void insert( const Msg & m , std::uint64_t h ) {
   if( f_mail.stopped() || ( m.g + m.h >= f_incumbent ) )
    return;
   auto [ x , fresh ] = f_table.find_or_insert( m.s , h , [ & ] {
    MNode * y = f_pool.allocate();
    *y = MNode{ m.s , m.g , m.h , m.parent , m.has_parent , false , false };
    return( y );
    } );
   if( fresh )
    f_open.push( x );
   else
    if( m.g < x->g ) {
     x->g = m.g;
     x->parent = m.parent;
     x->has_parent = m.has_parent;
     x->closed = false;
     f_open.push( x );
     }
    else
     ++f_st.duplicates;
   }

  // receives and processes all the messages that have arrived, the nodes
  // of a batch hashed and prefetched first
  void receive( void ) {
   v_msg.clear();
   f_mail.poll( [ & ]( const Mail & x ) {
    if( const Msg * m = std::get_if< Msg >( & x ) )
     v_msg.push_back( *m );
    else
     f_incumbent = std::min( f_incumbent , std::get< Incumbent >( x ).c );
    } );
   const std::size_t k = v_msg.size();
   v_hash.resize( k );
   for( std::size_t i = 0 ; i < k ; ++i ) {
    v_hash[ i ] = f_a.f_domain.hash( v_msg[ i ].s );
    f_table.prefetch( v_hash[ i ] );
    }
   for( std::size_t i = 0 ; i < k ; ++i )
    insert( v_msg[ i ] , v_hash[ i ] );
   }

  // fills the Solution on every process: counters, status, path
  void gather( SolutionT & sol ) {
   SearchStats & st = sol.stats;
   unsigned long long mine[ 7 ] = { f_st.expanded , f_st.reexpanded ,
                                    f_st.generated , f_st.duplicates ,
                                    f_st.sent , f_mail.batches() ,
                                    f_pool.size() };
   std::vector< unsigned long long > all( 7 * std::size_t( f_p ) );
   MPI_Allgather( mine , 7 , MPI_UNSIGNED_LONG_LONG , all.data() , 7 ,
                  MPI_UNSIGNED_LONG_LONG , f_comm );
   for( unsigned q = 0 ; q < f_p ; ++q ) {
    const auto * x = all.data() + 7 * std::size_t( q );
    st.expanded += x[ 0 ];
    st.reexpanded += x[ 1 ];
    st.generated += x[ 2 ];
    st.duplicates += x[ 3 ];
    st.sent += x[ 4 ];
    st.messages += x[ 5 ];
    st.stored += x[ 6 ];
    st.expanded_by.push_back( x[ 0 ] );
    st.stored_by.push_back( x[ 6 ] );
    }

   int limit = f_limit_hit ? 1 : 0 , any_limit = 0;
   MPI_Allreduce( & limit , & any_limit , 1 , MPI_INT , MPI_MAX , f_comm );
   if( any_limit ) {
    sol.status = Status::node_limit;
    return;
    }

   struct { double c; int r; } best{ f_goal ? double( f_goal->g )
                                            : infinite_cost< double >() ,
                                     f_rank } , win{};
   MPI_Allreduce( & best , & win , 1 , MPI_DOUBLE_INT , MPI_MINLOC ,
                  f_comm );
   if( win.c == infinite_cost< double >() )
    return;                                    // no goal: unsolvable

   // the path, from the goal back to the root, one owner at a time
   struct Step {
    State s;
    Cost g;
    State parent;
    bool has_parent;
    } step{};
   if( f_rank == win.r )
    step = Step{ f_goal->state , f_goal->g , f_goal->parent ,
                 f_goal->has_parent };
   MPI_Bcast( & step , int( sizeof( Step ) ) , MPI_BYTE , win.r , f_comm );
   sol.status = Status::solved;
   sol.cost = step.g;
   std::vector< State > path{ step.s };
   while( step.has_parent ) {
    const State s = step.parent;
    const int owner = int( f_dist.owner( s , f_p ) );
    if( f_rank == owner ) {
     const MNode * x = f_table.find( s , f_a.f_domain.hash( s ) );
     if( ! x )
      throw( std::logic_error( "HDAStarMPI::solve: lost node" ) );
     step = Step{ x->state , x->g , x->parent , x->has_parent };
     }
    MPI_Bcast( & step , int( sizeof( Step ) ) , MPI_BYTE , owner , f_comm );
    path.push_back( step.s );
    }
   sol.path.assign( path.rbegin() , path.rend() );
   }

  const HDAStarMPI & f_a;                       // the search
  MPI_Comm f_comm;                              // the communicator
  ADS::MPIMailbox< Mail > f_mail;               // its messages
  int f_rank;                                   // this process
  unsigned f_p;                                 // the number of processes
  Dist f_dist;                                  // the distribution
  NodePool< MNode > f_pool;                     // its nodes
  NodeTable< MNode > f_table;                   // its CLOSED list
  Open_t f_open;                                // its OPEN list
  std::vector< Msg > v_msg;                     // the nodes received
  std::vector< std::uint64_t > v_hash;          // their hashes
  Cost f_incumbent = infinite_cost< Cost >();   // the best cost known
  MNode * f_goal = nullptr;                     // its best goal
  bool f_limit_hit = false;                     // it hit the node limit
  SearchStats f_st;                             // its counters
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                    ///< the work distribution
 ParallelOptions f_options;                            ///< the options
 MPI_Comm f_comm;                                      ///< the processes

 };  // end( class( HDAStarMPI ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_HDAStarMPI */

/*--------------------------------------------------------------------------*/
/*------------------------- End File HDAStarMPI.h --------------------------*/
/*--------------------------------------------------------------------------*/
