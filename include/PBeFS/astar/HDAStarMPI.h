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
 *  7           else if m is an abort then stop expanding
 *  8        if OPEN_i has a node of f < z then
 *  9           expand it as at lines 7 to 15 of HDA*, a goal that lowers z
 *                 sending z to all the other processes
 * 10        else send all the outboxes
 * 11        if ENDED(OPEN_i has no node of f < z and the outboxes are empty)
 * 12           then break                              ▷ see MPIMailbox.h
 * 13     wait for the sends in progress
 * 14     the process of the least goal broadcasts it, and the path is rebuilt
 *           from it, each state broadcast with its parent by its owner
 * \endcode
 *
 * The nodes travel as bytes, packed in batches of batch nodes and sent with
 * non-blocking sends, and the end is detected by the four-counter method of
 * Mattern [see MPIMailbox], which by line 11 a process joins only when it
 * has no node of \f$ f < z \f$ and nothing to send. A process that reaches
 * the node limit sends an abort to the others (line 7), after which all of
 * them only drain their messages until the end is detected. The path is
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

#include <cstring>

#include <stdexcept>

#include <type_traits>

#include <vector>

#include <mpi.h>

#include "PBeFS/astar/HDAStar.h"

#include "PBeFS/core/MPIMailbox.h"

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
 * detected with the four-counter method of Mattern [see MPIMailbox], a
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

 enum Tag { c_tag_work = 7001 , c_tag_incumbent , c_tag_abort };

 // the OPEN list: buckets for integer costs, the heap otherwise
 using Open_t = std::conditional_t< std::is_integral_v< Cost > ,
                                    BucketOpenList< MNode > ,
                                    HeapOpenList< MNode > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the search on one process

 struct Process {

  explicit Process( const HDAStarMPI & a )
   : f_a( a ) , f_comm( a.f_comm ) , f_mail( a.f_comm ) ,
     f_dist( a.f_dist ) {
   MPI_Comm_rank( f_comm , & f_rank );
   int np;
   MPI_Comm_size( f_comm , & np );
   f_p = unsigned( np );
   v_out.resize( f_p );
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
    f_mail.test();

    if( MNode * n = next() ) {
     expand( n );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      flush();
      }
     if( f_a.f_options.node_limit && ! f_abort &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      f_abort = f_limit_hit = true;
      for( unsigned q = 0 ; q < f_p ; ++q )
       if( q != unsigned( f_rank ) )
        f_mail.post_value( q , c_tag_abort , Cost( 0 ) );
      }
     }
    else
     flush();

    if( f_mail.ended( idle() ) )
     break;
    }
   f_mail.finish();                      // all received: they complete
   }

  // true if there is nothing to expand and nothing to send
  bool idle( void ) {
   if( next_peek() )
    return( false );
   for( auto & o : v_out )
    if( ! o.empty() )
     return( false );
   return( true );
   }

  // the best node of OPEN with f < incumbent, not removed; nullptr if none
  MNode * next_peek( void ) {
   if( f_abort )
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
       f_mail.post_value( q , c_tag_incumbent , n->g );
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
    v_out[ dst ].push_back( m );
    ++f_st.sent;
    if( v_out[ dst ].size() >= f_a.f_options.batch )
     post( dst );
    } );
   }

  // adds the node of m, whose state has hash h, to OPEN, unless it is a
  // duplicate
  void insert( const Msg & m , std::uint64_t h ) {
   if( f_abort || ( m.g + m.h >= f_incumbent ) )
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

  // sends the outbox of dst
  void post( unsigned dst ) {
   auto & o = v_out[ dst ];
   std::vector< char > b( o.size() * sizeof( Msg ) );
   std::memcpy( b.data() , o.data() , b.size() );
   o.clear();
   f_mail.post( dst , c_tag_work , b );
   ++f_st.messages;
   }

  // sends the non-empty outboxes
  void flush( void ) {
   for( unsigned d = 0 ; d < f_p ; ++d )
    if( ! v_out[ d ].empty() )
     post( d );
   }

  // receives and processes all the messages that have arrived
  void receive( void ) {
   int tag;
   while( f_mail.receive( v_in , tag ) ) {
    if( tag == c_tag_work ) {
     const std::size_t k = v_in.size() / sizeof( Msg );
     v_msg.resize( k );
     v_hash.resize( k );
     std::memcpy( v_msg.data() , v_in.data() , k * sizeof( Msg ) );
     for( std::size_t i = 0 ; i < k ; ++i ) {   // prefetch the batch first
      v_hash[ i ] = f_a.f_domain.hash( v_msg[ i ].s );
      f_table.prefetch( v_hash[ i ] );
      }
     for( std::size_t i = 0 ; i < k ; ++i )
      insert( v_msg[ i ] , v_hash[ i ] );
     }
    else
     if( tag == c_tag_incumbent ) {
      const Cost c = MPIMailbox::value< Cost >( v_in );
      if( c < f_incumbent )
       f_incumbent = c;
      }
     else
      if( tag == c_tag_abort )
       f_abort = true;
    }
   }

  // fills the Solution on every process: counters, status, path
  void gather( SolutionT & sol ) {
   SearchStats & st = sol.stats;
   unsigned long long mine[ 7 ] = { f_st.expanded , f_st.reexpanded ,
                                    f_st.generated , f_st.duplicates ,
                                    f_st.sent , f_st.messages ,
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
  MPIMailbox f_mail;                            // its messages
  int f_rank;                                   // this process
  unsigned f_p;                                 // the number of processes
  Dist f_dist;                                  // the distribution
  NodePool< MNode > f_pool;                     // its nodes
  NodeTable< MNode > f_table;                   // its CLOSED list
  Open_t f_open;                                // its OPEN list
  std::vector< std::vector< Msg > > v_out;      // its outboxes
  std::vector< char > v_in;                     // the last message
  std::vector< Msg > v_msg;                     // its nodes
  std::vector< std::uint64_t > v_hash;          // their hashes
  Cost f_incumbent = infinite_cost< Cost >();   // the best cost known
  MNode * f_goal = nullptr;                     // its best goal
  bool f_abort = false;                         // stop expanding
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
