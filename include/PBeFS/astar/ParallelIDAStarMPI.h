/*--------------------------------------------------------------------------*/
/*----------------------- File ParallelIDAStarMPI.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ParallelIDAStarMPI, the synchronous
 * scheme of ParallelIDAStar [see ParallelIDAStar.h] on a distributed-memory
 * machine, with one MPI process per worker.
 *
 * Every process builds the same frontier, since its rules are
 * deterministic, so that the units need not travel; an iteration gives
 * each process a set of units, searches them, and ends with two reductions
 * over all the units, the minimum of the thresholds and the maximum of the
 * expansions (a unit searched by one process, the others giving infinity
 * and 0), after which every process knows the next bound and the next
 * assignment:
 *
 * \code{.unparsed}
 *  1  procedure PIDA*-MPI(s₀) on the process i of p
 *  2     build the units U by the frontier rule;  t_u ← f(u), w_u ← 0
 *  3     loop
 *  4        β ← min_u t_u;  if β = ∞ then return unsolvable
 *  5        A ← {u : t_u ≤ β}, by decreasing w_u (then by u)
 *  6        give each unit of A in turn to the process of least load,
 *              the load growing by max(w_u, 1)
 *  7        for each unit u given to i: (goal, t_u, w_u) ← DFS(u, β),
 *              a goal telling all the other processes to stop
 *  8        t ← min-reduce(t);  w ← max-reduce(w)  ▷ others: ∞ and 0
 *  9        if some process has a goal then broadcast its path; return it
 * \endcode
 *
 * The assignment of line 6 is the longest processing time first of Graham
 * (1969) on the expansions of the last search of each unit, which predict
 * those of the next one within the factor of growth of an iteration; the
 * processes compute it alike, from the same reduced vectors. A process that
 * finds a goal sends a message to the others, which poll for it every 1024
 * expansions and stop, since every goal found is optimal. State and Cost
 * have to be trivially copyable, since the path travels as bytes.
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

#ifndef __PBeFS_ParallelIDAStarMPI
 #define __PBeFS_ParallelIDAStarMPI
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstdint>

#include <limits>

#include <stdexcept>

#include <type_traits>

#include <vector>

#include <mpi.h>

#include "PBeFS/astar/ParallelIDAStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ParallelIDAStarMPI ------------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel iterative-deepening A* below a frontier, one MPI process
/// per worker
/** The options are those of ParallelIDAStar, of which threads, schedule,
 * asynchronous, split, table_depth and deterministic are not used (the
 * processes are those of the communicator); node_limit bounds the
 * expansions of each process. MPI must be initialized by the caller; all
 * the processes of the communicator call solve(), and all of them get the
 * same Solution, whose counters are summed over them. */

template< Domain D >
class ParallelIDAStarMPI {

 static_assert( std::is_trivially_copyable_v< state_t< D > > ,
                "ParallelIDAStarMPI: the states have to be trivially "
                "copyable" );

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain and the communicator

 ParallelIDAStarMPI( const D & domain , const PIDAOptions & options = {} ,
                     MPI_Comm comm = MPI_COMM_WORLD )
  : f_domain( domain ) , f_options( options ) , f_comm( comm ) {
  int np;
  MPI_Comm_size( f_comm , & np );
  f_options.threads = unsigned( np );
  f_options.table_depth = 0;
  f_options.deterministic = false;
  f_options.split = 0;
  f_options.asynchronous = false;
  f_options.normalize();
  if( f_options.frontier == 0 )
   f_options.frontier = std::size_t( 128 ) * f_options.threads;
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution (on every process)

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  f_iterations = 0;
  int rank , np;
  MPI_Comm_rank( f_comm , & rank );
  MPI_Comm_size( f_comm , & np );
  const ParallelIDAStar< D , ADS::ThreadPool > local( f_domain , f_options );
  typename ParallelIDAStar< D , ADS::ThreadPool >::Frontier fr(
   f_options.frontier );
  SolutionT sol;
  if( f_options.rule == FrontierRule::depth )
   local.build_depth( fr , sol );
  else
   local.build_best_first( fr , sol );
  if( fr.done ) {
   sol.stats.seconds = clock.seconds();
   return( sol );
   }

  const std::size_t n = fr.units.size();
  std::vector< double > t( n ) , next( n );
  std::vector< unsigned long long > w( n , 0 ) , work( n );
  for( std::size_t u = 0 ; u < n ; ++u )
   t[ u ] = double( fr.units[ u ].node->g + fr.units[ u ].node->h );
  Run r( *this , rank , np );
  bool found = false , limit = false;
  std::vector< State > best;
  Cost cost{};
  for( ; ; ) {
   const double beta = *std::min_element( t.begin() , t.end() );
   if( beta == std::numeric_limits< double >::infinity() )
    break;
   ++f_iterations;
   std::vector< std::size_t > a;
   for( std::size_t u = 0 ; u < n ; ++u )
    if( f_options.unit_bounds ? t[ u ] <= beta
                              : t[ u ] < std::numeric_limits< double >::
                                         infinity() )
     a.push_back( u );
   std::stable_sort( a.begin() , a.end() , [ & ]( std::size_t x ,
                                                  std::size_t y ) {
    return( w[ x ] > w[ y ] );
    } );
   std::vector< unsigned long long > load( std::size_t( np ) , 0 );
   next = t;                          // the units not searched keep theirs
   work = w;
   for( std::size_t u : a ) {
    const std::size_t q = std::size_t( std::min_element( load.begin() ,
                                                          load.end() ) -
                                       load.begin() );
    load[ q ] += std::max< unsigned long long >( w[ u ] , 1 );
    next[ u ] = std::numeric_limits< double >::infinity();   // the owner's
    work[ u ] = 0;
    if( q != std::size_t( rank ) )
     continue;
    if( r.stopped() ) {                   // its bound kept
     next[ u ] = t[ u ];
     continue;
     }
    const NodeT * x = fr.units[ u ].node;
    const std::uint64_t before = r.expanded;
    const Cost nb = r.search( x , Cost( beta ) );
    next[ u ] = nb == infinite_cost< Cost >()
                ? std::numeric_limits< double >::infinity() : double( nb );
    work[ u ] = r.expanded - before;
    if( r.goal_found && best.empty() ) {
     best = trace_path( x );
     best.insert( best.end() , r.path.begin() + 1 , r.path.end() );
     cost = r.goal_cost;
     r.tell_others();
     }
    }
   MPI_Allreduce( next.data() , t.data() , int( n ) , MPI_DOUBLE , MPI_MIN ,
                  f_comm );
   MPI_Allreduce( work.data() , w.data() , int( n ) ,
                  MPI_UNSIGNED_LONG_LONG , MPI_MAX , f_comm );
   r.drain();
   int lim = r.limit ? 1 : 0 , any = 0;
   MPI_Allreduce( & lim , & any , 1 , MPI_INT , MPI_MAX , f_comm );
   struct { double c; int r; } mine{ best.empty()
                                     ? std::numeric_limits< double >::
                                       infinity() : double( cost ) , rank } ,
                                   win{};
   MPI_Allreduce( & mine , & win , 1 , MPI_DOUBLE_INT , MPI_MINLOC ,
                  f_comm );
   if( win.c < std::numeric_limits< double >::infinity() ) {
    unsigned long long len = best.size();
    MPI_Bcast( & len , 1 , MPI_UNSIGNED_LONG_LONG , win.r , f_comm );
    best.resize( std::size_t( len ) );
    MPI_Bcast( best.data() , int( len * sizeof( State ) ) , MPI_BYTE , win.r ,
               f_comm );
    MPI_Bcast( & cost , int( sizeof( Cost ) ) , MPI_BYTE , win.r , f_comm );
    found = true;
    break;
    }
   if( any ) {
    limit = true;
    break;
    }
   }

  // the counters, summed over the processes
  SearchStats & st = sol.stats;
  unsigned long long mine[ 3 ] = { r.expanded , r.generated , r.longest };
  std::vector< unsigned long long > all( 3 * std::size_t( np ) );
  MPI_Allgather( mine , 3 , MPI_UNSIGNED_LONG_LONG , all.data() , 3 ,
                 MPI_UNSIGNED_LONG_LONG , f_comm );
  for( int q = 0 ; q < np ; ++q ) {
   st.expanded += all[ 3 * std::size_t( q ) ];
   st.generated += all[ 3 * std::size_t( q ) + 1 ];
   st.stored += all[ 3 * std::size_t( q ) + 2 ];
   st.expanded_by.push_back( all[ 3 * std::size_t( q ) ] );
   st.stored_by.push_back( all[ 3 * std::size_t( q ) + 2 ] );
   }
  if( found ) {
   sol.status = Status::solved;
   sol.cost = cost;
   sol.path = std::move( best );
   }
  else
   sol.status = limit ? Status::node_limit : Status::unsolvable;
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the iterations of the last solve()

 [[nodiscard]] std::uint64_t iterations( void ) const noexcept {
  return( f_iterations );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using NodeT = Node< State , Cost >;

 static constexpr int c_tag = 0x5049;           // the tag of a stop

 /// the searches of a process
 struct Run {

  const ParallelIDAStarMPI & a;
  int rank , np;
  std::vector< State > path;
  std::uint64_t expanded = 0 , generated = 0 , longest = 0;
  bool goal_found = false , limit = false , told = false;
  Cost goal_cost{};
  int heard = 0;                                // stops received
  std::vector< MPI_Request > sends;
  int dummy = 0;

  Run( const ParallelIDAStarMPI & s , int r , int p )
   : a( s ) , rank( r ) , np( p ) {}

  [[nodiscard]] bool stopped( void ) const noexcept {
   return( goal_found || limit || ( heard > 0 ) );
   }

  // the cost-bounded search below the unit x with the bound b; returns
  // the least f cut
  Cost search( const NodeT * x , Cost b ) {
   path.assign( 1 , x->state );
   Cost nb = infinite_cost< Cost >();
   if( x->g + x->h > b )
    return( x->g + x->h );
   const State * above = x->parent ? & x->parent->state : nullptr;
   dfs( x->g , x->h , x->op , b , above , nb );
   return( nb );
   }

  // true if a goal is found below the last state of path
  bool dfs( Cost g , Cost h , std::uint8_t op , Cost b , const State * above ,
            Cost & nb ) {
   const State u = path.back();
   if( a.f_domain.is_goal( u ) ) {
    goal_found = true;
    goal_cost = g;
    return( true );
    }
   if( a.f_options.node_limit && ( expanded >= a.f_options.node_limit ) ) {
    limit = true;
    return( false );
    }
   if( ( ( ++expanded & 1023 ) == 0 ) && poll() )
    return( false );
   longest = std::max< std::uint64_t >( longest , path.size() );
   bool done = false;
   const State * gp = path.size() >= 2 ? & path[ path.size() - 2 ] : above;
   for_each_tree_child( a.f_domain , u , h , op ,
                        [ & ]( const State & v , Cost c , Cost hv ,
                               std::uint8_t o ) {
    if( done || stopped() )
     return;
    ++generated;
    if constexpr( ! OperatorDomain< D > )
     if( gp && ( v == *gp ) )
      return;
    if( g + c + hv > b ) {
     nb = std::min( nb , g + c + hv );
     return;
     }
    path.push_back( v );
    if( dfs( g + c , hv , o , b , nullptr , nb ) )
     done = true;
    else
     path.pop_back();
    } );
   return( done );
   }

  // true if another process has found a goal
  bool poll( void ) {
   int flag = 0;
   MPI_Status s;
   MPI_Iprobe( MPI_ANY_SOURCE , c_tag , a.f_comm , & flag , & s );
   if( flag ) {
    int x;
    MPI_Recv( & x , 1 , MPI_INT , s.MPI_SOURCE , c_tag , a.f_comm ,
              MPI_STATUS_IGNORE );
    ++heard;
    }
   return( heard > 0 );
   }

  // tells the other processes that a goal is found
  void tell_others( void ) {
   if( told )
    return;
   told = true;
   for( int q = 0 ; q < np ; ++q )
    if( q != rank ) {
     sends.emplace_back();
     MPI_Isend( & dummy , 1 , MPI_INT , q , c_tag , a.f_comm ,
                & sends.back() );
     }
   }

  // receives the stops not yet received, so that none is left for the
  // next search, and completes the sends
  void drain( void ) {
   int t = told ? 1 : 0 , tellers = 0;
   MPI_Allreduce( & t , & tellers , 1 , MPI_INT , MPI_SUM , a.f_comm );
   const int expect = tellers - ( told ? 1 : 0 );
   while( heard < expect ) {
    int x;
    MPI_Recv( & x , 1 , MPI_INT , MPI_ANY_SOURCE , c_tag , a.f_comm ,
              MPI_STATUS_IGNORE );
    ++heard;
    }
   if( ! sends.empty() )
    MPI_Waitall( int( sends.size() ) , sends.data() , MPI_STATUSES_IGNORE );
   sends.clear();
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 PIDAOptions f_options;                              ///< the options
 MPI_Comm f_comm;                                    ///< the processes
 mutable std::uint64_t f_iterations = 0;     ///< iterations of last solve

 };  // end( class( ParallelIDAStarMPI ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ParallelIDAStarMPI */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ParallelIDAStarMPI.h ---------------------*/
/*--------------------------------------------------------------------------*/
