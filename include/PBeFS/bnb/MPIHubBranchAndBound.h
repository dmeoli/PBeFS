/*--------------------------------------------------------------------------*/
/*---------------------- File MPIHubBranchAndBound.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for MPIHubBranchAndBoundSearch, the distributed
 * Branch-and-Bound (B&B) with one MPI process per worker organized as a
 * master, hubs and workers, as in the Abstract Library for Parallel Search
 * (ALPS) of Xu, Ralphs, Ladányi and Saltzman (2005): each process searches
 * its own subtrees, i.e., the open nodes it holds, a unit of work at a
 * time, and the work moves only when a hub tells a donor to give some of
 * its nodes to a receiver. The problem, the options and the Solution are
 * those of BranchAndBound.h and ParallelBranchAndBound.h, and the design is
 * discussed in docs/branch-and-bound.md.
 *
 * The processes are split in clusters of consecutive ranks, the first of
 * each being its hub, and the rank 0 is also the master; every process
 * searches. With \f$ u \f$ the unit of work (unit_work evaluations), the
 * workload of a process is its quantity, the number of its open nodes, and
 * its quality, the mean bound of its \f$ e \f$ best nodes (elite), and the
 * search is the following, the procedures of the sequential search being
 * those of BranchAndBound.h:
 *
 * \code{.unparsed}
 *  1  procedure HUB-B&B() on the process i
 *  2     Q_i ← {root} on the rank 0, ∅ elsewhere
 *  3     loop
 *  4        for each message m received do
 *  5           nodes: push them in Q_i;  incumbent z': z ← min(z, z')
 *  6           donate k to j: send to j up to k nodes of Q_i, the 2nd,
 *                 4th, ... best, keeping the best;  report to the hub
 *  7           (hub) report of a worker: BALANCE()
 *  8           (master) report of a hub, or a need: send to the hub of the
 *                 largest quantity "donate to j"
 *  9        evaluate up to u nodes of Q_i by PROCESS
 * 10        report the quantity and the quality of Q_i to the hub, after
 *              each unit and when Q_i becomes empty
 * 11        if Q_i = ∅ and the mailbox detects the end then break
 * 12  procedure BALANCE() on a hub
 * 13     R ← the workers of quantity ≤ need_work
 * 14     if R ≠ ∅ then D ← those of quantity > need_work   ▷ quantity
 * 15     else (quality, once every |cluster| reports) D ← those whose quality
 *           is better than the mean by more than the fraction
 *           donor_threshold, R those worse than it by more than
 *           receiver_threshold
 * 16     pair D and R, the best donor with the neediest receiver, and send
 *           to each donor "donate half of your nodes to its receiver"
 * 17     the receivers left with no donor are sent to the master
 * \endcode
 *
 * The end is detected by the four-counter method of Mattern in the mailbox
 * of ADS, as in MPIBranchAndBoundSearch: a process joins a wave only when
 * it has no node, and every message (nodes, incumbents, reports, orders)
 * is counted; a process with no node sends no report but the one of line
 * 10 when it becomes empty, and the hubs and the master answer only
 * reports, so that the messages stop when every queue is empty. A
 * receiver paired with a donor is not paired again until that donor
 * reports, which it does after each order. The incumbent and the end are
 * those of MPIBranchAndBoundSearch, and so are the Solution and the
 * requirements on the problem [see MPIBranchAndBound]; the counter sent of
 * the statistics holds the nodes donated, and messages the donations.
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

#ifndef __PBeFS_MPIHubBranchAndBound
 #define __PBeFS_MPIHubBranchAndBound
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstdint>

#include <optional>

#include <span>

#include <utility>

#include <vector>

#include "PBeFS/bnb/MPIBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BBHubParams ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the parameters of the balance of MPIHubBranchAndBoundSearch

//! [branch-and-bound-sig-hub]
struct BBHubParams {
 unsigned cluster = 4;             ///< the processes of a cluster
 std::size_t unit_work = 32;       ///< evaluations between two reports
 std::size_t need_work = 2;        ///< open nodes at which one is a receiver
 std::size_t elite = 4;            ///< the best nodes of the quality
 double donor_threshold = 0.1;     ///< the quality better than the mean by
 double receiver_threshold = 0.1;  ///< ... and worse, as fractions
 };
//! [branch-and-bound-sig-hub]

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MPIHubBranchAndBound -----------------------*/
/*--------------------------------------------------------------------------*/
/// the distributed B&B of a master, hubs and workers, with the balance of
/// the quantity and of the quality of their work [see MPIHubBranchAndBound.h]
/** MPI must be initialized by the caller; all the processes of the
 * communicator call solve(), and all of them get the same Solution but for
 * the frontier, as for MPIBranchAndBoundSearch. The node limit is split
 * evenly among the processes. */

template< MPIBranchAndBound P >
class MPIHubBranchAndBoundSearch {

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
 /// options (threads is ignored), the parameters of the balance and the
 /// communicator

 explicit MPIHubBranchAndBoundSearch( const P & problem ,
                                      const Options & options = {} ,
                                      const BBHubParams & hub = {} ,
                                      MPI_Comm comm = MPI_COMM_WORLD )
  : f_p( problem ) , f_options( options ) , f_hub( hub ) , f_comm( comm ) {
  f_options.normalize();
  if( f_hub.cluster == 0 )
   f_hub.cluster = 1;
  if( f_hub.unit_work == 0 )
   f_hub.unit_work = 1;
  if( f_hub.elite == 0 )
   f_hub.elite = 1;
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

 // the tags of the records
 enum Tag : std::uint32_t {
  c_work = 1 ,           // nodes
  c_incumbent = 2 ,      // an incumbent value
  c_report = 3 ,         // worker to hub: rank, quantity, quality
  c_donate = 4 ,         // hub to donor: receiver, how many
  c_need = 5 ,           // hub to master: a receiver
  c_hub = 6 ,            // hub to master: hub, its quantity
  c_ask = 7              // master to hub: a receiver to serve
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the search on one process

 struct Process {

  Process( const MPIHubBranchAndBoundSearch & a , const Stopwatch & c )
   : f_a( a ) , f_comm( a.f_comm ) , f_mail( a.f_comm , 1 ) , f_clock( c ) ,
     f_q( a.f_options.selection ) , f_fr( a.f_options ) {
   MPI_Comm_rank( f_comm , & f_rank );
   int np;
   MPI_Comm_size( f_comm , & np );
   f_np = unsigned( np );
   const unsigned cl = std::min( a.f_hub.cluster , f_np );
   f_hubrank = ( unsigned( f_rank ) / cl ) * cl;
   f_size = std::min( cl , f_np - f_hubrank );
   if( unsigned( f_rank ) == f_hubrank ) {
    v_qty.assign( f_size , 0 );
    v_qual.assign( f_size , BBValueTraits< Value >::infinity() );
    v_served.assign( f_size , -1 );
    }
   if( f_rank == 0 )
    v_hubqty.assign( ( f_np + cl - 1 ) / cl , 0 );
   f_inc = a.f_options.incumbent ?
           detail::bb_signed< P >( *a.f_options.incumbent ) :
           BBValueTraits< Value >::infinity();
   const std::size_t lim = a.f_options.node_limit;
   f_limit = lim ? ( lim + f_np - 1 ) / f_np : 0;
   }

  // sends a record of tag t made of the values xs
  template< class ... T >
  void say( unsigned to , std::uint32_t t , const T & ... xs ) {
   std::vector< char > b;
   ( detail::bb_put_bytes( xs , b ) , ... );
   f_mail.send( to , t , std::as_bytes( std::span( b ) ) );
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
      say( q , c_incumbent , v );
    return( true );
    };
   for( auto & l : starts )
    push( Item{ std::move( l.node ) , l.bound , l.depth } );
   starts.clear();
   bool reported_empty = false;
   for( ;; ) {
    receive();
    std::size_t done = 0;
    while( ( done < f_a.f_hub.unit_work ) && ! f_mail.stopped() ) {
     if( prunable() )
      break;
     Item it = f_q.pop();
     v_kids.clear();
     detail::bb_process< P >( it.node , it.bound , it.depth , w , inc ,
                              improve , [ & ]( Node && c , Value b ,
                                               unsigned d ) {
                               v_kids.push_back( Item{ std::move( c ) , b ,
                                                       d } );
                               } , o , f_st , f_fr );
     for( Item & k : v_kids )
      push( std::move( k ) );
     ++done;
     check_limits();
     }
    if( done > 0 ) {
     report();
     reported_empty = f_q.empty();
     continue;
     }
    if( ! reported_empty ) {           // empty: say it once
     report();
     reported_empty = true;
     }
    if( f_mail.quiescent() )
     break;
    }
   f_mail.leftovers( [ & ]( std::uint32_t t ,
                            std::span< const std::byte > b ) {
    take( t , b , true );
    } );
   }

  // pushes a node in the queue
  void push( Item && it ) {
   const Value key = detail::bb_key( f_a.f_p , it.node , it.bound ,
                                      it.depth , f_a.f_options.selection );
   f_q.push( std::move( it.node ) , it.bound , it.depth , key );
   f_st.max_open = std::max( f_st.max_open , f_q.size() );
   }

  // true if no node of the queue can improve the incumbent, which then
  // empties it
  bool prunable( void ) {
   if( f_q.empty() )
    return( true );
   if( ! detail::bb_prunable< P >( f_q.best_bound() , f_inc ,
                                   f_a.f_options ) )
    return( false );
   f_st.pruned_extracted += f_fr.add( f_q , BBLeafKind::pruned );
   return( true );
   }

  // stops all the processes if a limit is reached
  void check_limits( void ) {
   const auto & o = f_a.f_options;
   const bool nodes = f_limit && ( f_st.evaluated >= f_limit );
   const bool time = ( o.time_limit > 0 ) &&
                     ( f_clock.seconds() >= o.time_limit );
   if( f_mail.stopped() || ! ( nodes || time ) )
    return;
   f_limit_hit = true;
   f_time_hit = time;
   f_mail.stop();
   }

  // the quality of the queue: the mean bound of its elite nodes
  Value quality( void ) const {
   const std::size_t e = f_a.f_hub.elite;
   if( f_q.empty() )
    return( BBValueTraits< Value >::infinity() );
   if( e == 1 )
    return( f_q.best_bound() );
   std::vector< Value > b;
   f_q.for_each( [ & ]( const Item & it ) { b.push_back( it.bound ); } );
   const std::size_t k = std::min( e , b.size() );
   std::partial_sort( b.begin() , b.begin() + std::ptrdiff_t( k ) ,
                      b.end() );
   long double s = 0;
   for( std::size_t j = 0 ; j < k ; ++j )
    s += static_cast< long double >( b[ j ] );
   return( Value( s / static_cast< long double >( k ) ) );
   }

  // reports the workload to the hub
  void report( void ) {
   say( f_hubrank , c_report , std::uint32_t( f_rank ) ,
        std::uint64_t( f_q.size() ) , quality() );
   }

  // gives up to k nodes to the process r: the 2nd, 4th, ... best, keeping
  // the best
  void donate( unsigned r , std::uint64_t k ) {
   if( ( f_q.size() >= 2 ) && ( k > 0 ) ) {
    std::vector< Item > all;
    f_q.drain( [ & ]( Item && it ) { all.push_back( std::move( it ) ); } );
    std::stable_sort( all.begin() , all.end() ,
                      []( const Item & x , const Item & y ) {
                       return( x.bound < y.bound );
                       } );
    std::vector< char > b;
    std::uint64_t given = 0;
    for( std::size_t j = 0 ; j < all.size() ; ++j )
     if( ( j % 2 == 1 ) && ( given < k ) ) {
      detail::bb_put_bytes( all[ j ].bound , b );
      detail::bb_put_bytes( all[ j ].depth , b );
      detail::bb_put_node( f_a.f_p , all[ j ].node , b );
      ++given;
      }
     else
      push( std::move( all[ j ] ) );
    f_mail.send( r , c_work , std::as_bytes( std::span( b ) ) );
    f_st.sent += given;
    ++f_st.messages;
    }
   report();                           // the donor reports after an order
   }

  // the hub: a report of the worker of rank r
  void hub_report( unsigned r , std::uint64_t qty , Value qual ) {
   const unsigned j = r - f_hubrank;
   v_qty[ j ] = qty;
   v_qual[ j ] = qual;
   for( auto & s : v_served )          // its receivers may ask again
    if( s == int( r ) )
     s = -1;
   balance();
   std::uint64_t sum = 0;
   for( const auto x : v_qty )
    sum += x;
   say( 0 , c_hub , std::uint32_t( f_hubrank ) , sum );
   }

  // the hub: pairs donors and receivers of its cluster
  void balance( void ) {
   const auto & h = f_a.f_hub;
   std::vector< unsigned > rec , don;
   for( unsigned j = 0 ; j < f_size ; ++j )
    if( ( v_qty[ j ] <= h.need_work ) && ( v_served[ j ] < 0 ) )
     rec.push_back( j );
   if( ! rec.empty() ) {               // quantity
    for( unsigned j = 0 ; j < f_size ; ++j )
     if( v_qty[ j ] > h.need_work )
      don.push_back( j );
    std::stable_sort( don.begin() , don.end() , [ & ]( unsigned x ,
                                                       unsigned y ) {
     return( v_qual[ x ] < v_qual[ y ] );
     } );
    }
   else {                              // quality, once a round
    if( ++f_reports % f_size != 0 )
     return;
    long double mean = 0;
    std::size_t k = 0;
    for( unsigned j = 0 ; j < f_size ; ++j )
     if( BBValueTraits< Value >::finite( v_qual[ j ] ) ) {
      mean += static_cast< long double >( v_qual[ j ] );
      ++k;
      }
    if( k < 2 )
     return;
    mean /= static_cast< long double >( k );
    const long double scale = std::max( std::fabs( mean ) ,
                                        static_cast< long double >( 1 ) );
    for( unsigned j = 0 ; j < f_size ; ++j ) {
     if( ! BBValueTraits< Value >::finite( v_qual[ j ] ) ||
         ( v_served[ j ] >= 0 ) )
      continue;
     const long double d = ( static_cast< long double >( v_qual[ j ] ) -
                             mean ) / scale;
     if( ( d < - h.donor_threshold ) && ( v_qty[ j ] > h.need_work ) )
      don.push_back( j );
     else
      if( d > h.receiver_threshold )
       rec.push_back( j );
     }
    }
   const std::size_t pairs = std::min( rec.size() , don.size() );
   for( std::size_t k = 0 ; k < pairs ; ++k ) {
    const unsigned d = f_hubrank + don[ k ] , r = f_hubrank + rec[ k ];
    say( d , c_donate , std::uint32_t( r ) ,
         std::uint64_t( std::max< std::uint64_t >( 1 ,
                                                   v_qty[ don[ k ] ] / 2 ) ) );
    v_served[ rec[ k ] ] = int( d );
    v_qty[ don[ k ] ] /= 2;
    }
   for( std::size_t k = pairs ; k < rec.size() ; ++k )
    if( v_qty[ rec[ k ] ] == 0 ) {     // no donor here: ask the master
     v_served[ rec[ k ] ] = int( f_np );
     say( 0 , c_need , std::uint32_t( f_hubrank + rec[ k ] ) );
     }
   }

  // the master: a receiver r with no donor in its cluster
  void master_need( unsigned r ) {
   const unsigned c = std::min( f_a.f_hub.cluster , f_np );
   std::size_t best = v_hubqty.size();
   for( std::size_t k = 0 ; k < v_hubqty.size() ; ++k )
    if( ( k != r / c ) && ( v_hubqty[ k ] > f_a.f_hub.need_work ) &&
        ( ( best == v_hubqty.size() ) ||
          ( v_hubqty[ k ] > v_hubqty[ best ] ) ) )
     best = k;
   if( best < v_hubqty.size() )
    say( unsigned( best * c ) , c_ask , std::uint32_t( r ) );
   else
    say( ( r / c ) * c , c_ask , std::uint32_t( f_np ) + r );   // none
   }

  // the hub: serve the receiver r of another cluster, or (r ≥ p) release
  // the receiver r - p of its own, for which the master found no donor
  void hub_ask( unsigned r ) {
   if( r >= f_np ) {
    const unsigned j = r - f_np - f_hubrank;
    if( ( j < f_size ) && ( v_served[ j ] == int( f_np ) ) )
     v_served[ j ] = -1;
    return;
    }
   unsigned best = f_size;
   for( unsigned j = 0 ; j < f_size ; ++j )
    if( ( v_qty[ j ] > f_a.f_hub.need_work ) &&
        ( ( best == f_size ) || ( v_qty[ j ] > v_qty[ best ] ) ) )
     best = j;
   if( best == f_size )
    return;
   say( f_hubrank + best , c_donate , std::uint32_t( r ) ,
        std::uint64_t( std::max< std::uint64_t >( 1 , v_qty[ best ] / 2 ) ) );
   v_qty[ best ] /= 2;
   }

  // receives and processes all the messages that have arrived
  void receive( void ) {
   f_mail.poll( [ & ]( std::uint32_t t , std::span< const std::byte > b ) {
    take( t , b , false );
    } );
   }

  // a record of tag t; its nodes go among the open leaves if open
  void take( std::uint32_t t , std::span< const std::byte > b , bool open ) {
   const char * in = reinterpret_cast< const char * >( b.data() );
   const char * end = in + b.size();
   switch( t ) {
    case( c_incumbent ):
     f_inc = std::min( f_inc , detail::bb_get_bytes< Value >( in ) );
     return;
    case( c_work ):
     while( in < end ) {
      const Value v = detail::bb_get_bytes< Value >( in );
      const unsigned d = detail::bb_get_bytes< unsigned >( in );
      Node n = detail::bb_get_node( f_a.f_p , in );
      if( open )
       f_fr.add( std::move( n ) , v , d , BBLeafKind::open );
      else
       push( Item{ std::move( n ) , v , d } );
      }
     return;
    default:
     break;
    }
   if( open )
    return;
   switch( t ) {
    case( c_report ): {
     const auto r = detail::bb_get_bytes< std::uint32_t >( in );
     const auto q = detail::bb_get_bytes< std::uint64_t >( in );
     const auto v = detail::bb_get_bytes< Value >( in );
     hub_report( r , q , v );
     return;
     }
    case( c_donate ): {
     const auto r = detail::bb_get_bytes< std::uint32_t >( in );
     const auto k = detail::bb_get_bytes< std::uint64_t >( in );
     donate( r , k );
     return;
     }
    case( c_need ):
     master_need( detail::bb_get_bytes< std::uint32_t >( in ) );
     return;
    case( c_hub ): {
     const auto h = detail::bb_get_bytes< std::uint32_t >( in );
     const auto q = detail::bb_get_bytes< std::uint64_t >( in );
     v_hubqty[ h / std::min( f_a.f_hub.cluster , f_np ) ] = q;
     return;
     }
    case( c_ask ):
     hub_ask( detail::bb_get_bytes< std::uint32_t >( in ) );
     return;
    default:
     return;
    }
   }

  // fills the Solution on every process, as MPIBranchAndBoundSearch does
  void gather( SolutionT & sol ) {
   const Value inf = BBValueTraits< Value >::infinity();
   using Tr = BBValueTraits< Value >;
   BBStats & st = sol.stats;
   constexpr int k = 13;
   unsigned long long mine[ k ] = {
    f_st.evaluated , f_st.branched , f_st.generated , f_st.infeasible ,
    f_st.pruned_generated , f_st.pruned_extracted , f_st.pruned_evaluated ,
    f_st.incumbents , f_st.sent , f_st.messages , f_st.max_open ,
    f_st.solutions , f_st.recomputed };
   std::vector< unsigned long long > all( k * std::size_t( f_np ) );
   MPI_Allgather( mine , k , MPI_UNSIGNED_LONG_LONG , all.data() , k ,
                  MPI_UNSIGNED_LONG_LONG , f_comm );
   for( unsigned q = 0 ; q < f_np ; ++q ) {
    const auto * x = all.data() + k * std::size_t( q );
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
    s.solutions = x[ 11 ];
    s.recomputed = x[ 12 ];
    st.add( s );
    }
   int flags[ 2 ] = { f_limit_hit ? 1 : 0 , f_time_hit ? 1 : 0 } , any[ 2 ];
   MPI_Allreduce( flags , any , 2 , MPI_INT , MPI_MAX , f_comm );
   f_fr.add( f_q , BBLeafKind::open );
   const Value glow = detail::bb_all_min( f_fr.low , f_comm ).first;
   const Value gclosed = detail::bb_all_min( f_fr.closed , f_comm ).first;
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
    sol.status = detail::bb_status< P >( z , f_a.f_options );
   if( Tr::finite( z ) )
    sol.value = detail::bb_signed< P >( z );
   const Value lb = std::min( z , glow );
   sol.bound = detail::bb_signed< P >( lb );
   sol.gap = Tr::finite( z ) && Tr::finite( lb ) ? Value( z - lb ) : inf;
   sol.closed_bound = detail::bb_signed< P >( gclosed );
   f_fr.template take< P >( sol.frontier );
   }

  const MPIHubBranchAndBoundSearch & f_a;       // the search
  MPI_Comm f_comm;                              // the communicator
  ADS::MPIByteMailbox f_mail;                   // its messages
  const Stopwatch & f_clock;                    // the clock of the search
  int f_rank = 0;                               // the rank of the process
  unsigned f_np = 1;                            // the number of processes
  unsigned f_hubrank = 0;                       // the rank of its hub
  unsigned f_size = 1;                          // the size of its cluster
  BBQueue< Node , Value > f_q;                  // the open nodes
  std::vector< Item > v_kids;                   // the children of a node
  Value f_inc;                                  // the incumbent value
  std::optional< Node > f_node;                 // its node, if found here
  Value f_node_value = Value( 0 );              // the value of f_node
  std::size_t f_limit;                          // its share of node_limit
  bool f_limit_hit = false;                     // a limit was reached here
  bool f_time_hit = false;                      // ... the time limit
  std::vector< std::uint64_t > v_qty;           // hub: the quantities
  std::vector< Value > v_qual;                  // hub: the qualities
  std::vector< int > v_served;                  // hub: the donor of each
  std::vector< std::uint64_t > v_hubqty;        // master: those of the hubs
  std::size_t f_reports = 0;                    // hub: the reports so far
  BBStats f_st;                                 // the counters
  detail::BBFrontier< Node , Value > f_fr;      // the leaves owned here
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const P & f_p;                                       ///< the problem
 Options f_options;                                   ///< the options
 BBHubParams f_hub;                                   ///< the balance
 MPI_Comm f_comm;                                     ///< the communicator

 };  // end( class( MPIHubBranchAndBoundSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MPIHubBranchAndBound */

/*--------------------------------------------------------------------------*/
/*-------------------- End File MPIHubBranchAndBound.h ---------------------*/
/*--------------------------------------------------------------------------*/
