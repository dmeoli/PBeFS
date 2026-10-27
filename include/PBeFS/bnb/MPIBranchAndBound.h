/*--------------------------------------------------------------------------*/
/*------------------------ File MPIBranchAndBound.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MPIBranchAndBound, the distributed Branch-and-Bound (B&B)
 * with one MPI process per worker: each process has a queue of open nodes
 * and sends the children to their owners, as DistributedBranchAndBound does
 * with threads [see ParallelBranchAndBound.h] and HDAStarMPI with the
 * states of A* [see HDAStarMPI.h]. The problem, the options and the
 * Solution are those of BranchAndBound.h and ParallelBranchAndBound.h, and
 * the design is discussed in docs/branch-and-bound.md.
 *
 * The search is the following, with the notation of the page "Parallel
 * branch-and-bound" of the documentation, the procedures of the sequential
 * search [see BranchAndBound.h] and DIVE of the parallel ones [see
 * ParallelBranchAndBound.h], on the process \f$ i \f$ of \f$ p \f$, with
 * its queue \f$ Q_i \f$, whose least bound is written βmin(Q_i), and the
 * incumbent \f$ z \f$ known to it:
 *
 * \code{.unparsed}
 *  1  procedure MPI-B&B(R_i) on the process i       ▷ R_i: its own nodes
 *  2     Q_i ← R_i;  d ← nil
 *  3     loop
 *  4        for each message m received do
 *  5           if m is a batch of nodes then push them in Q_i
 *  6           else if m is an incumbent z' then z ← min(z, z')
 *  7           else if m is an abort then stop evaluating
 *  8        if not stopped and (d ≠ nil or not PRUNABLE(βmin(Q_i))) then
 *  9           n ← d, or the first node of Q_i, removed;  d ← nil
 * 10           K ← PROCESS(n), a value v < z sending v to all the others
 * 11           d ← DIVE(K)
 * 12           for each m ∈ K but d do
 * 13              j ← hash(m) mod p, or a random process without hash(m)
 * 14              if j = i then push m in Q_i else send m to j
 * 15           every k evaluations, send the outboxes not yet full
 * 16           if the share of the node limit or the time limit is reached
 *                 then send an abort to all the others;  stop evaluating
 * 17        else send all the outboxes
 * 18        if ENDED(stopped, or d = nil and PRUNABLE(βmin(Q_i)),
 *                  and the outboxes are empty) then break
 * 19     the process of the least z sends its node to all
 * \endcode
 *
 * When no node of \f$ Q_i \f$ can improve \f$ z \f$ (lines 8 and 18) the
 * nodes of \f$ Q_i \f$ become leaves (pruned), and after an abort the nodes
 * left in \f$ Q_i \f$ are leaves (open), which bound the optimum. The end is
 * detected by the four-counter method of Mattern [see MPIMailbox.h], which
 * by line 18 a process joins only when it has nothing to evaluate and nothing
 * to send; every process gets the same incumbent, value and bound, and keeps
 * the leaves that it owns.
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

#ifndef __PBeFS_MPIBranchAndBound
 #define __PBeFS_MPIBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstring>

#include <random>

#include <type_traits>

#include <vector>

#include <mpi.h>

#include "PBeFS/core/MPIMailbox.h"

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
/// a B&B problem whose nodes and values can travel between processes
/** The values are trivially copyable, and travel as their bytes; either the
 * node is trivially copyable too, or the problem gives
 * serialize( n , out ), which appends the bytes of the node n to the
 * std::vector< char > out, and deserialize( in ), which rebuilds a node
 * from the bytes at the const char * in and moves in past them. */

template< class P >
concept MPIBranchAndBound = BranchAndBound< P > &&
 std::is_trivially_copyable_v< bb_value_t< P > > &&
 ( std::is_trivially_copyable_v< typename P::node_type > ||
   requires( const P & p , const typename P::node_type & n ,
             std::vector< char > & out , const char * & in ) {
    p.serialize( n , out );
    { p.deserialize( in ) } -> std::convertible_to< typename P::node_type >;
    } );

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

/// appends the bytes of x to out

template< class T >
void bb_put_bytes( const T & x , std::vector< char > & out )
{
 const auto * b = reinterpret_cast< const char * >( & x );
 out.insert( out.end() , b , b + sizeof( T ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads an object from the bytes at in, and moves in past them

template< class T >
T bb_get_bytes( const char * & in )
{
 T x;
 std::memcpy( static_cast< void * >( & x ) , in , sizeof( T ) );
 in += sizeof( T );
 return( x );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// appends the node n to out

template< MPIBranchAndBound P >
void bb_put_node( const P & p , const typename P::node_type & n ,
                  std::vector< char > & out )
{
 if constexpr( requires { p.serialize( n , out ); } )
  p.serialize( n , out );
 else
  bb_put_bytes( n , out );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// reads a node from the bytes at in, and moves in past them

template< MPIBranchAndBound P >
typename P::node_type bb_get_node( const P & p , const char * & in )
{
 if constexpr( requires { p.deserialize( in ); } )
  return( p.deserialize( in ) );
 else
  return( bb_get_bytes< typename P::node_type >( in ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the values v of all the processes of comm, in the order of their ranks

template< class V >
std::vector< V > bb_all_values( const V & v , MPI_Comm comm )
{
 int np;
 MPI_Comm_size( comm , & np );
 std::vector< V > all( std::size_t( np ) , v );
 MPI_Allgather( & v , int( sizeof( V ) ) , MPI_BYTE , all.data() ,
                int( sizeof( V ) ) , MPI_BYTE , comm );
 return( all );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the least of the values v of the processes of comm, and the least rank
/// holding it

template< class V >
std::pair< V , int > bb_all_min( const V & v , MPI_Comm comm )
{
 const std::vector< V > all = bb_all_values( v , comm );
 int r = 0;
 for( int q = 1 ; q < int( all.size() ) ; ++q )
  if( all[ q ] < all[ r ] )
   r = q;
 return( std::pair( all[ r ] , r ) );
 }

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MPIBranchAndBound -------------------------*/
/*--------------------------------------------------------------------------*/
/// the distributed B&B with one MPI process per worker
/** Each process extracts the nodes from its own queue, and sends each child,
 * in batches and with non-blocking sends, to the process that owns it:
 * \f$ \mathrm{hash}( n ) \bmod p \f$ if the problem gives hash( n ), a
 * random process otherwise (with
 * depth-first or plunging, the child that dives is kept). A process that
 * improves the incumbent sends its value to all the others, and calls the
 * callback of the options; the node of the incumbent stays on it until the
 * end, when it is sent to all. The end is detected with the four-counter
 * method of Mattern [see MPIMailbox], a process joining a wave when it has
 * no node that can improve the incumbent and no batch to send. The node
 * limit is split evenly among the processes; a process that reaches its
 * share, or the time limit, tells the others, and all of them stop
 * evaluating and drain the messages until the end is detected, the nodes
 * left open bounding the optimum.
 *
 * MPI must be initialized by the caller, with at least MPI_THREAD_FUNNELED
 * if the problem uses threads; all the processes of the communicator call
 * solve(), and all of them get the same Solution, but for the frontier
 * [see BBSolution], each process keeping the leaves that it owns. Started
 * from given nodes, each process starts from its own ones, and these
 * leaves can be given back to the same process for a later search. */

template< MPIBranchAndBound P >
class MPIBranchAndBoundSearch {

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
 /// options (threads is ignored) and the communicator

 explicit MPIBranchAndBoundSearch( const P & problem ,
                                   const Options & options = {} ,
                                   MPI_Comm comm = MPI_COMM_WORLD )
  : f_p( problem ) , f_options( options ) , f_comm( comm ) {
  f_options.normalize();
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the tree from the root, returning the Solution (on every
 /// process)

 [[nodiscard]] SolutionT solve( void ) const {
  int rank;
  MPI_Comm_rank( f_comm , & rank );
  if( rank == 0 )
   return( run( std::nullopt ) );
  return( run( std::vector< Leaf >() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// searches the subtrees of the given nodes of this process, returning
 /// the Solution (on every process) [see BranchAndBoundSearch::solve()]

 [[nodiscard]] SolutionT solve( std::vector< Leaf > seeds ) const {
  return( run( std::move( seeds ) ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 SolutionT run( std::optional< std::vector< Leaf > > seeds ) const
 {
  Stopwatch clock;
  Process w( *this , clock );
  w.run( detail::bb_starts( f_p , std::move( seeds ) ) );
  SolutionT sol;
  w.gather( sol );
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

 using Item = typename BBQueue< Node , Value >::Item;

 enum Tag { c_tag_work = 7101 , c_tag_incumbent , c_tag_abort };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the search on one process

 struct Process {

  Process( const MPIBranchAndBoundSearch & a , const Stopwatch & c )
   : f_a( a ) , f_comm( a.f_comm ) , f_mail( a.f_comm ) , f_clock( c ) ,
     f_q( a.f_options.selection ) , f_fr( a.f_options.keep_frontier ) {
   MPI_Comm_rank( f_comm , & f_rank );
   int np;
   MPI_Comm_size( f_comm , & np );
   f_np = unsigned( np );
   v_out.resize( f_np );
   v_count.assign( f_np , 0 );
   f_gen.seed( a.f_options.seed + unsigned( f_rank ) );
   f_inc = a.f_options.incumbent ?
           detail::bb_signed< P >( *a.f_options.incumbent ) :
           BBValueTraits< Value >::infinity();
   const std::size_t lim = a.f_options.node_limit;
   f_limit = lim ? ( lim + f_np - 1 ) / f_np : 0;
   }

  // the main loop, from the given nodes
  void run( std::vector< Leaf > && starts ) {
   auto w = detail::bb_worker( f_a.f_p , unsigned( f_rank ) );
   const auto & o = f_a.f_options;
   const auto inc = [ & ] { return( f_inc ); };
   const auto improve = [ & ]( Value v , const Node & n ) {
    if( v >= f_inc )
     return( false );
    f_inc = f_node_value = v;
    f_node = n;
    if( o.on_incumbent )
     o.on_incumbent( detail::bb_signed< P >( v ) );
    for( unsigned q = 0 ; q < f_np ; ++q )
     if( q != unsigned( f_rank ) )
      f_mail.post_value( q , c_tag_incumbent , v );
    return( true );
    };

   for( auto & l : starts )
    push( Item{ std::move( l.node ) , l.bound , l.depth } );
   starts.clear();
   unsigned since_flush = 0;
   for( ;; ) {
    receive();
    f_mail.test();
    if( std::optional< Item > it = next() ) {
     v_kids.clear();
     detail::bb_process< P >( it->node , it->bound , it->depth , w , inc ,
                              improve , [ & ]( Node && c , Value b ,
                                               unsigned d ) {
                               v_kids.push_back( Item{ std::move( c ) , b ,
                                                       d } );
                               } , o , f_st , f_fr );
     const std::size_t k = detail::bb_dive( v_kids , f_q.best_bound() ,
                                            f_inc , o );
     for( std::size_t j = 0 ; j < v_kids.size() ; ++j )
      if( j == k )
       f_dive.emplace( std::move( v_kids[ j ] ) );
      else
       send( std::move( v_kids[ j ] ) );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      flush();
      }
     check_limits();
     }
    else
     flush();
    if( f_mail.ended( idle() ) )
     break;
    }
   f_mail.finish();                      // all received: they complete
   }

  // pushes a node in the queue
  void push( Item && it ) {
   const Value key = detail::bb_key( f_a.f_p , it.node , it.bound ,
                                      it.depth , f_a.f_options.selection );
   f_q.push( std::move( it.node ) , it.bound , it.depth , key );
   f_st.max_open = std::max( f_st.max_open , f_q.size() );
   }

  // true if there is nothing to evaluate and nothing to send
  bool idle( void ) {
   if( ! f_abort && ( f_dive || ! prunable() ) )
    return( false );
   for( auto & o : v_out )
    if( ! o.empty() )
     return( false );
   return( true );
   }

  // true if no node of the queue can improve the incumbent; then, without
  // an abort, the queue is emptied
  bool prunable( void ) {
   if( f_q.empty() )
    return( true );
   const Value best = f_q.best_bound();
   if( ! detail::bb_prunable( best , f_inc , f_a.f_options ) )
    return( false );
   f_st.pruned_extracted += f_fr.add( f_q , BBLeafKind::pruned );
   return( true );
   }

  // the next node: the dive, or the best of the queue; none if there is
  // none that can improve the incumbent, or after an abort
  std::optional< Item > next( void ) {
   if( f_abort )
    return( std::nullopt );
   if( f_dive ) {
    std::optional< Item > it( std::move( f_dive ) );
    f_dive.reset();
    return( it );
    }
   if( prunable() )
    return( std::nullopt );
   return( f_q.pop() );
   }

  // stops the evaluations of all the processes if a limit is reached
  void check_limits( void ) {
   const auto & o = f_a.f_options;
   const bool nodes = f_limit && ( f_st.evaluated >= f_limit );
   const bool time = ( o.time_limit > 0 ) &&
                     ( f_clock.seconds() >= o.time_limit );
   if( f_abort || ! ( nodes || time ) )
    return;
   f_abort = f_limit_hit = true;
   f_time_hit = time;
   for( unsigned q = 0 ; q < f_np ; ++q )
    if( q != unsigned( f_rank ) )
     f_mail.post_value( q , c_tag_abort , 0.0 );
   }

  // sends a child to its owner, or keeps it if the process owns it
  void send( Item && it ) {
   unsigned dst;
   if constexpr( HashedBranchAndBound< P > )
    dst = unsigned( std::uint64_t( f_a.f_p.hash( it.node ) ) % f_np );
   else
    dst = unsigned( random_below( f_gen , f_np ) );
   if( dst == unsigned( f_rank ) ) {
    push( std::move( it ) );
    return;
    }
   auto & o = v_out[ dst ];
   detail::bb_put_bytes( it.bound , o );
   detail::bb_put_bytes( it.depth , o );
   detail::bb_put_node( f_a.f_p , it.node , o );
   ++f_st.sent;
   if( ++v_count[ dst ] >= f_a.f_options.batch )
    post( dst );
   }

  // sends the outbox of dst
  void post( unsigned dst ) {
   f_mail.post( dst , c_tag_work , v_out[ dst ] );
   v_count[ dst ] = 0;
   ++f_st.messages;
   }

  // sends the non-empty outboxes
  void flush( void ) {
   for( unsigned d = 0 ; d < f_np ; ++d )
    if( ! v_out[ d ].empty() )
     post( d );
   }

  // receives and processes all the messages that have arrived
  void receive( void ) {
   int tag;
   while( f_mail.receive( v_in , tag ) ) {
    const char * in = v_in.data();
    if( tag == c_tag_work )
     while( in < v_in.data() + v_in.size() ) {
      const Value b = detail::bb_get_bytes< Value >( in );
      const unsigned d = detail::bb_get_bytes< unsigned >( in );
      push( Item{ detail::bb_get_node( f_a.f_p , in ) , b , d } );
      }
    else
     if( tag == c_tag_incumbent )
      f_inc = std::min( f_inc , detail::bb_get_bytes< Value >( in ) );
     else
      if( tag == c_tag_abort )
       f_abort = true;
    }
   }

  // fills the Solution on every process: counters, status, bound, node
  void gather( SolutionT & sol ) {
   const Value inf = BBValueTraits< Value >::infinity();
   BBStats & st = sol.stats;
   unsigned long long mine[ 11 ] = {
    f_st.evaluated , f_st.branched , f_st.generated , f_st.infeasible ,
    f_st.pruned_generated , f_st.pruned_extracted , f_st.pruned_evaluated ,
    f_st.incumbents , f_st.sent , f_st.messages , f_st.max_open };
   std::vector< unsigned long long > all( 11 * std::size_t( f_np ) );
   MPI_Allgather( mine , 11 , MPI_UNSIGNED_LONG_LONG , all.data() , 11 ,
                  MPI_UNSIGNED_LONG_LONG , f_comm );
   for( unsigned q = 0 ; q < f_np ; ++q ) {
    const auto * x = all.data() + 11 * std::size_t( q );
    BBStats s;
    s.evaluated = x[ 0 ];
    s.branched = x[ 1 ];
    s.generated = x[ 2 ];
    s.infeasible = x[ 3 ];
    s.pruned_generated = x[ 4 ];
    s.pruned_extracted = x[ 5 ];
    s.pruned_evaluated = x[ 6 ];
    s.incumbents = x[ 7 ];
    s.sent = x[ 8 ];
    s.messages = x[ 9 ];
    s.max_open = std::size_t( x[ 10 ] );
    st.add( s );
    }

   int flags[ 2 ] = { f_limit_hit ? 1 : 0 , f_time_hit ? 1 : 0 } , any[ 2 ];
   MPI_Allreduce( flags , any , 2 , MPI_INT , MPI_MAX , f_comm );

   // the bound: the nodes pruned by the tolerance and those left open
   f_fr.add( f_q , BBLeafKind::open );
   if( f_dive ) {
    f_fr.add( std::move( f_dive->node ) , f_dive->bound , f_dive->depth ,
              BBLeafKind::open );
    f_dive.reset();
    }
   using Tr = BBValueTraits< Value >;
   const Value glow = detail::bb_all_min( f_fr.low , f_comm ).first;

   // the incumbent and the process holding its node
   const auto [ wv , wr ] = detail::bb_all_min( f_node ? f_node_value : inf ,
                                                f_comm );
   const Value z = detail::bb_all_min( f_inc , f_comm ).first;
   if( Tr::finite( wv ) && ( wv <= z ) ) {
    std::vector< char > bytes;
    if( f_rank == wr )
     detail::bb_put_node( f_a.f_p , *f_node , bytes );
    unsigned long long n = bytes.size();
    MPI_Bcast( & n , 1 , MPI_UNSIGNED_LONG_LONG , wr , f_comm );
    bytes.resize( std::size_t( n ) );
    MPI_Bcast( bytes.data() , int( n ) , MPI_BYTE , wr , f_comm );
    const char * in = bytes.data();
    sol.node.emplace( detail::bb_get_node( f_a.f_p , in ) );
    }

   if( any[ 0 ] )
    sol.status = any[ 1 ] ? BBStatus::time_limit : BBStatus::node_limit;
   else
    sol.status = Tr::finite( z ) ? BBStatus::optimal : BBStatus::infeasible;
   if( Tr::finite( z ) )
    sol.value = detail::bb_signed< P >( z );
   const Value lb = std::min( z , glow );
   sol.bound = detail::bb_signed< P >( lb );
   sol.gap = Tr::finite( z ) && Tr::finite( lb ) ? Value( z - lb ) : inf;
   f_fr.template take< P >( sol.frontier );
   }

  const MPIBranchAndBoundSearch & f_a;          // the search
  MPI_Comm f_comm;                              // the communicator
  MPIMailbox f_mail;                            // its messages
  const Stopwatch & f_clock;                    // the clock of the search
  int f_rank = 0;                               // the rank of the process
  unsigned f_np = 1;                            // the number of processes
  BBQueue< Node , Value > f_q;                  // the open nodes
  std::optional< Item > f_dive;                 // the node of the dive
  std::vector< Item > v_kids;                   // the children of a node
  std::vector< std::vector< char > > v_out;     // the outboxes
  std::vector< std::size_t > v_count;           // the nodes in each one
  std::vector< char > v_in;                     // the last message
  std::mt19937_64 f_gen;                        // the random owners
  Value f_inc;                                  // the incumbent value
  std::optional< Node > f_node;                 // its node, if found here
  Value f_node_value = Value( 0 );              // the value of f_node
  std::size_t f_limit;                          // its share of node_limit
  bool f_abort = false;                         // stop evaluating
  bool f_limit_hit = false;                     // a limit was reached here
  bool f_time_hit = false;                      // ... the time limit
  BBStats f_st;                                 // the counters
  detail::BBFrontier< Node , Value > f_fr;      // the leaves owned here
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                           ///< the options
 MPI_Comm f_comm;                                     ///< the communicator

 };  // end( class( MPIBranchAndBoundSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MPIBranchAndBound */

/*--------------------------------------------------------------------------*/
/*---------------------- End File MPIBranchAndBound.h ----------------------*/
/*--------------------------------------------------------------------------*/
