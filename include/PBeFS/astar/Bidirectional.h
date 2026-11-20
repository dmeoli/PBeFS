/*--------------------------------------------------------------------------*/
/*-------------------------- File Bidirectional.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for what the front-to-end bidirectional searches share: the
 * options BidirectionalOptions, the node BiNode, which holds the state with
 * its g, h and parent in both directions, the class template BiFronts,
 * the table of the nodes of both searches with the incumbent and the
 * expansion of a node in either direction, the lazy priority queue
 * BiQueue, and the adapters ReverseDomain, which swaps the start and the
 * goal of a BidirectionalDomain, and BlindDomain, which drops its
 * heuristics.
 *
 * A forward search from the start s and a backward search from the goal t
 * share one table of the states, each state with the cost
 * \f$ g_F \f$ of the best path found from s and the cost \f$ g_B \f$ of the
 * best path found to t, \f$ \infty \f$ until the search of that direction
 * reaches it. A state generated in a direction with a \f$ g \f$ finite in
 * the other one is a path from s to t of cost \f$ g_F + g_B \f$, which
 * lowers the incumbent \f$ U \f$, the cost of the best path found: this is
 * the immediate solution detection of Sturtevant and Chen (IJCAI 2016),
 * done with the same lookup that detects the duplicates, since the two
 * directions have one table. The searches stop when \f$ U \f$ is not
 * larger than a lower bound on the cost \f$ C^* \f$ of an optimal path,
 * which each algorithm computes from its open lists; the bounds and their
 * proofs are in the page "Bidirectional search" of the documentation. The
 * open lists hold no copy of a node but an entry with its \f$ g \f$, as
 * in AStar: an entry is stale when the node has been closed in that
 * direction or reached again with a smaller \f$ g \f$, and it is dropped
 * when it reaches the top.
 *
 * With admissible heuristics, a node whose \f$ g \f$ plus the cost of an
 * edge, or whose \f$ f \f$, is not smaller than \f$ U \f$ cannot lie on a
 * path cheaper than the incumbent, and it is not stored. With consistent
 * heuristics (BidirectionalOptions::consistent), moreover, a search does
 * not open a node that the other one has closed, the "nipping" of BS* [see
 * Kwa, AIJ 1989]. Without the nipping, a node expanded in one direction
 * has there its optimal \f$ g \f$; with it, this holds, while \f$ U \f$
 * is larger than \f$ C^* \f$, for the nodes of the optimal paths, so
 * that a node of such a path closed by the other direction has both its
 * optimal values of \f$ g \f$ and every path through it is known by
 * \f$ g_F + g_B \f$. A node off the optimal paths may instead be expanded
 * with a \f$ g \f$ larger than its optimal one, when its cheaper path
 * goes through a node not opened by the nipping, which costs some search
 * but not the optimality of the path returned. Both claims hold for every
 * algorithm of the library whose priority does not decrease along an edge
 * and increases with \f$ g \f$ at a fixed state, ties broken by the larger
 * \f$ g \f$, which is the case of A*, MM, fMM, NBS, DVCBS, GBFHS and BAE*
 * under consistent heuristics.
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

#ifndef __PBeFS_Bidirectional
 #define __PBeFS_Bidirectional
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <limits>

#include <stdexcept>

#include <type_traits>

#include <utility>

#include <vector>

#include "ADS/priority/BucketQueue.h"

#include "ADS/priority/DaryHeap.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the two directions of a bidirectional search, used as indices

//! [bidirectional-sig-sides]
enum Side : unsigned {
 forward_side = 0 ,   ///< the search from the start
 backward_side = 1    ///< the search from the goal
 };
//! [bidirectional-sig-sides]

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// how a bidirectional search chooses the direction of the next expansion

//! [bidirectional-sig-sidechoice]
enum class SideChoice {
 automatic ,    ///< the choice of the algorithm as published
 priority ,     ///< the direction whose open list has the least priority
 alternate ,    ///< forward and backward in turn
 cardinality    ///< the direction with fewer open nodes [Pohl, 1971]
 };
//! [bidirectional-sig-sidechoice]

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BidirectionalOptions -----------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the bidirectional searches, besides those of all the
/// searches
/** The node_limit bounds the nodes stored by both directions together. */

//! [bidirectional-sig-options]
struct BidirectionalOptions : public SearchOptions {

 /// the least cost \f$ \epsilon \f$ of an edge, by which the bound
 /// \f$ g_F + g_B + \epsilon \f$ is raised; negative to take
 /// min_edge_cost() of a MinCostDomain, and 0 for any other Domain
 double epsilon = -1;

 /// true if both heuristics are consistent: a search then does not open a
 /// node that the other has closed, and BAE* may use its bound; false
 /// asks only for admissible heuristics, and BAE* refuses it
 bool consistent = true;

 /// true to stop at the largest of all the lower bounds valid under the
 /// heuristics (those of A*, MM and, with consistent heuristics, BAE*)
 /// instead of the one of the algorithm alone
 bool all_bounds = false;

 /// the direction of the next expansion, for the algorithms that expand
 /// one node at a time
 SideChoice side = SideChoice::automatic;

 /// the fraction \f$ p \in ( 0 , 1 ) \f$ of the optimal path at which
 /// fMM( p ) meets, and the share of the forward search in the split of
 /// GBFHS
 double fraction = 0.5;

 /// for NBS and DVCBS, true to make ready at once the nodes whose
 /// \f$ f \f$ equals the lower bound, which is the version that covers the
 /// pairs needed to find all the optimal paths [Shperberg et al., AAAI
 /// 2019], false for the version of Chen et al. (IJCAI 2017), which makes
 /// them ready one at a time
 bool eager = false;

 };  // end( struct( BidirectionalOptions ) )
//! [bidirectional-sig-options]

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS BiNode ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a node of a bidirectional search: a state with its data in both
/// directions
/** The index 0 is the forward search and 1 the backward one [see Side].
 * g[ d ] is infinite until the search d reaches the state, and h[ d ] is
 * valid once g[ d ] is finite (or when the search asks for both
 * heuristics of every node). flags holds, for each direction d, whether
 * the node is closed (bit d), has ever been expanded (bit 2 + d), is
 * open (bit 4 + d), and is ready in the sense of NBS and DVCBS (bit 6 +
 * d). */

template< class State , class Cost >
struct BiNode {
 State state;                      ///< the state
 Cost g[ 2 ];                      ///< g_F and g_B
 Cost h[ 2 ];                      ///< h_F and h_B
 BiNode * parent[ 2 ];             ///< the predecessor on each path
 std::uint8_t op[ 2 ];             ///< the operator that generated it
 std::uint8_t flags;               ///< closed, expanded and open bits

 /// true if closed in the direction d
 [[nodiscard]] bool closed( unsigned d ) const noexcept {
  return( ( flags >> d ) & 1u );
  }

 /// true if open in the direction d
 [[nodiscard]] bool open( unsigned d ) const noexcept {
  return( ( flags >> ( 4 + d ) ) & 1u );
  }

 /// true if ready in the direction d
 [[nodiscard]] bool ready( unsigned d ) const noexcept {
  return( ( flags >> ( 6 + d ) ) & 1u );
  }

 /// sets or clears the bit b of the flags
 void flag( unsigned b , bool on ) noexcept {
  flags = std::uint8_t( on ? ( flags | ( 1u << b ) )
                           : ( flags & ~( 1u << b ) ) );
  }

 /// \f$ f_d = g_d + h_d \f$
 [[nodiscard]] Cost f( unsigned d ) const noexcept {
  return( g[ d ] + h[ d ] );
  }

 };  // end( struct( BiNode ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS BiQueue -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a lazy priority queue of the nodes of one direction, by a key and then
/// by the larger g
/** An entry is a node with the key and the g it had when pushed; the
 * queue is the BucketQueue of ADS when both are integers (and the keys
 * nonnegative), and its DaryHeap otherwise. Neither decreases a key: the
 * caller pushes a node again when its key changes and drops the stale
 * entries when they reach the top. */

template< class NodeT , class Key , class Cost >
class BiQueue {

 public:

 /// an entry of the queue
 struct Entry {
  NodeT * node;                      ///< the node
  Key key;                           ///< its key when pushed
  Cost g;                            ///< its g when pushed
  };

 /// pushes the node n with the key k and the g g
 void push( NodeT * n , Key k , Cost g ) {
  if constexpr( c_buckets )
   f_q.push( n , k , g );
  else
   f_q.push( Entry{ n , k , g } );
  }

 /// the entry of least key, of largest g among those
 [[nodiscard]] Entry top( void ) {
  if constexpr( c_buckets ) {
   const auto e = f_q.topEntry();
   return( Entry{ e.item , e.f , e.g } );
   }
  else
   return( f_q.top() );
  }

 /// removes and returns top()
 Entry pop( void ) {
  if constexpr( c_buckets ) {
   const auto e = f_q.popEntry();
   return( Entry{ e.item , e.f , e.g } );
   }
  else
   return( f_q.pop() );
  }

 /// true if there is no entry
 [[nodiscard]] bool empty( void ) const noexcept { return( f_q.empty() ); }

 /// the number of entries, the stale ones included
 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_q.size() );
  }

 /// drops the entries at the top for which stale( e ) holds; true if an
 /// entry is left
 template< class Stale >
 bool clean( Stale && stale ) {
  while( ! empty() && stale( top() ) )
   static_cast< void >( pop() );
  return( ! empty() );
  }

 private:

 /// true if the keys are integers, for the bucket queue
 static constexpr bool c_buckets = std::is_integral_v< Key > &&
                                   std::is_integral_v< Cost > &&
                                   std::is_same_v< Key , Cost >;

 /// the order of the heap: the smaller key, then the larger g
 struct Better {
  bool operator()( const Entry & a , const Entry & b ) const noexcept {
   return( ( a.key < b.key ) || ( ( a.key == b.key ) && ( a.g > b.g ) ) );
   }
  };

 /// the bucket queue, named only when the keys are integers
 template< bool B , class = void >
 struct Pick { using type = ADS::DaryHeap< Entry , Better >; };
 template< class V >
 struct Pick< true , V > { using type = ADS::BucketQueue< NodeT * , Key >; };

 typename Pick< c_buckets >::type f_q;                     ///< the queue

 };  // end( class( BiQueue ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS BiFronts ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the two fronts of a bidirectional search: the table of the nodes, the
/// incumbent and the expansion of a node in a direction
/** A search builds a BiFronts, calls start() and then expand() on the
 * nodes it selects, and finish() at the end. The counters: expanded and
 * expanded_by (forward, backward), generated, duplicates (states reached
 * again with a g not smaller), reexpanded, and stored, the nodes in the
 * table. */

template< BidirectionalDomain D >
class BiFronts {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = BiNode< State , Cost >;             ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain, the options and whether every node
 /// needs both its heuristics from its creation

 BiFronts( const D & domain , const BidirectionalOptions & options ,
           bool both_h )
  : f_d( domain ) , f_o( options ) , f_both( both_h ) ,
    f_table( options.expected_nodes
             ? NodeTable< NodeT >::slots_for( options.expected_nodes )
             : 1024 ) ,
    f_eps( epsilon_of( domain , options ) ) {
  if( options.weight != 1 )
   throw( std::invalid_argument( "BiFronts: no weight" ) );
  }

/*-------------------------------- SEARCHING -------------------------------*/
 /// creates the start and the goal nodes, and calls push( d , n ) for the
 /// root n of each direction d; if they are the same state, U is 0

 template< class Push >
 void start( Push && push ) {
  const State s = f_d.initial_state() , t = f_d.goal_state();
  NodeT * a = node( s , forward_side , Cost( 0 ) ,
                    Cost( f_d.heuristic( s ) ) );
  NodeT * b = node( t , backward_side , Cost( 0 ) ,
                    Cost( f_d.backward_heuristic( t ) ) );
  f_root[ 0 ] = a;
  f_root[ 1 ] = b;
  if( ( a->g[ 1 ] < inf() ) ) {
   U = Cost( 0 );
   meet = a;
   }
  open_node( a , forward_side );
  push( forward_side , a );
  if( b != a ) {
   open_node( b , backward_side );
   push( backward_side , b );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// expands the node n in the direction d, calling push( m , g ) for each
 /// node m whose g in the direction d is set or lowered from g (infinite
 /// for a node new in that direction) and which is to be opened, and
 /// lowering U on each path found

 template< class Push >
 void expand( unsigned d , NodeT * n , Push && push ) {
  const unsigned o = 1 - d;
  n->flags = std::uint8_t( n->flags | ( 1u << d ) );
  ++st.expanded;
  ++f_exp[ d ];
  if( n->flags & ( 4u << d ) )
   ++st.reexpanded;
  n->flags = std::uint8_t( n->flags | ( 4u << d ) );
  f_kids.clear();
  const auto kid = [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
   const std::uint64_t hash = f_d.hash( t );
   f_table.prefetch( hash );
   f_kids.push_back( Kid{ t , c , ht , hash , op } );
   };
  if( d == forward_side )
   for_each_child( f_d , n->state , n->h[ 0 ] , n->op[ 0 ] , kid );
  else
   for_each_parent( f_d , n->state , n->h[ 1 ] , n->op[ 1 ] , kid );
  st.generated += f_kids.size();
  const Cost gn = n->g[ d ];
  for( const Kid & k : f_kids ) {
   const Cost g = gn + k.c;
   if( ( g >= U ) || ( g + k.h >= U ) )      // no cheaper path through it
    continue;
   auto [ m , fresh ] = f_table.find_or_insert( k.s , k.hash , [ & ] {
    NodeT * x = f_pool.allocate();
    init( x , k.s , d , g , k.h );
    return( x );
    } );
   Cost old = inf();
   if( fresh ) {
    m->parent[ d ] = n;
    m->op[ d ] = k.op;
    }
   else {
    old = m->g[ d ];
    if( ! cheaper( g , m->g[ d ] ) ) {
     ++st.duplicates;
     continue;
     }
    if( m->g[ d ] == inf() )
     m->h[ d ] = k.h;
    m->g[ d ] = g;
    m->parent[ d ] = n;
    m->op[ d ] = k.op;
    m->flags = std::uint8_t( m->flags & ~( 1u << d ) );   // reopened
    if( m->g[ o ] < inf() && ( g + m->g[ o ] < U ) ) {
     U = g + m->g[ o ];
     meet = m;
     }
    }
   if( f_o.consistent && m->closed( o ) )      // every path through m is
    continue;                                 // known already
   open_node( m , d );
   push( m , old );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// marks the node n open in the direction d, counting it once

 void open_node( NodeT * n , unsigned d ) noexcept {
  if( ! n->open( d ) ) {
   n->flags = std::uint8_t( n->flags | ( 16u << d ) );
   ++f_open[ d ];
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// marks the node n no longer open in the direction d

 void leave( NodeT * n , unsigned d ) noexcept {
  if( n->open( d ) ) {
   n->flags = std::uint8_t( n->flags & ~( 16u << d ) );
   --f_open[ d ];
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the entry of n with the g g is stale in the direction d: the
 /// node has been closed there, or improved, or (with consistent
 /// heuristics) closed in the other direction; a node closed in the other
 /// direction leaves the open list of d

 bool stale( unsigned d , NodeT * n , Cost g ) noexcept {
  if( n->closed( d ) || ( n->g[ d ] != g ) )
   return( true );
  if( f_o.consistent && n->closed( 1 - d ) ) {
   leave( n , d );
   return( true );
   }
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the node limit has been reached

 [[nodiscard]] bool full( void ) const noexcept {
  return( f_o.node_limit && ( f_pool.size() >= f_o.node_limit ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Solution: the path through the meeting node if U is finite and the
 /// search was not stopped, with the counters

 [[nodiscard]] SolutionT finish( bool stopped , double seconds ) {
  SolutionT sol;
  if( stopped )
   sol.status = Status::node_limit;
  else
   if( U < inf() ) {
    sol.status = Status::solved;
    sol.cost = U;
    sol.bound = sol.cost;
    for( const NodeT * x = meet ; x ; x = x->parent[ 0 ] )
     sol.path.push_back( x->state );
    std::reverse( sol.path.begin() , sol.path.end() );
    for( const NodeT * x = meet->parent[ 1 ] ; x ; x = x->parent[ 1 ] )
     sol.path.push_back( x->state );
    }
   else
    sol.status = Status::unsolvable;
  st.stored = f_pool.size();
  st.expanded_by = { f_exp[ 0 ] , f_exp[ 1 ] };
  st.stored_by = { st.stored };
  st.seconds = seconds;
  sol.stats = st;
  return( sol );
  }

/*------------------------------- ACCESSORS --------------------------------*/
 /// the infinite cost, the g of a direction that has not reached a node

 [[nodiscard]] static constexpr Cost inf( void ) noexcept {
  return( infinite_cost< Cost >() );
  }

 /// the least cost of an edge used by the bounds
 [[nodiscard]] Cost epsilon( void ) const noexcept { return( f_eps ); }

 /// the number of nodes open in the direction d
 [[nodiscard]] std::size_t open_count( unsigned d ) const noexcept {
  return( f_open[ d ] );
  }

 /// the expansions in the direction d
 [[nodiscard]] std::uint64_t expanded( unsigned d ) const noexcept {
  return( f_exp[ d ] );
  }

 /// the options
 [[nodiscard]] const BidirectionalOptions & options( void ) const noexcept {
  return( f_o );
  }

 /// the Domain
 [[nodiscard]] const D & domain( void ) const noexcept { return( f_d ); }

 /// \f$ a + b \f$, or the infinite cost if either is infinite
 [[nodiscard]] static Cost add( Cost a , Cost b ) noexcept {
  if( ( a == inf() ) || ( b == inf() ) )
   return( inf() );
  return( a + b );
  }

/*--------------------------------- FIELDS ---------------------------------*/

 Cost U = inf();                  ///< the cost of the best path found
 NodeT * meet = nullptr;          ///< the node of the best path found
 SearchStats st;                  ///< the counters

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a successor or predecessor: its state, the cost of the edge, its h in
 /// the direction of the search, its hash and its operator
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint64_t hash;
  std::uint8_t op;
  };

 // the epsilon of the options or of the Domain
 static Cost epsilon_of( const D & d , const BidirectionalOptions & o ) {
  if( o.epsilon >= 0 )
   return( Cost( o.epsilon ) );
  if constexpr( MinCostDomain< D > )
   return( Cost( d.min_edge_cost() ) );
  else
   return( Cost( 0 ) );
  }

 // a node of state s reached first in the direction d with g and h
 void init( NodeT * x , const State & s , unsigned d , Cost g , Cost h ) {
  x->state = s;
  x->g[ d ] = g;
  x->g[ 1 - d ] = inf();
  x->h[ d ] = h;
  x->h[ 1 - d ] = f_both ? ( d == forward_side
                             ? Cost( f_d.backward_heuristic( s ) )
                             : Cost( f_d.heuristic( s ) ) ) : Cost( 0 );
  x->parent[ 0 ] = x->parent[ 1 ] = nullptr;
  x->op[ 0 ] = x->op[ 1 ] = no_op;
  x->flags = 0;
  }

 // the node of s, reached in the direction d with g and h, made if absent
 NodeT * node( const State & s , unsigned d , Cost g , Cost h ) {
  auto [ m , fresh ] = f_table.find_or_insert( s , f_d.hash( s ) , [ & ] {
   NodeT * x = f_pool.allocate();
   init( x , s , d , g , h );
   return( x );
   } );
  if( ! fresh ) {
   m->g[ d ] = g;
   m->h[ d ] = h;
   }
  return( m );
  }

 const D & f_d;                                    ///< the Domain
 const BidirectionalOptions & f_o;                 ///< the options
 bool f_both;                       ///< both heuristics at the creation
 NodePool< NodeT > f_pool;                         ///< the nodes
 NodeTable< NodeT > f_table;                       ///< the table
 Cost f_eps;                                       ///< the epsilon
 NodeT * f_root[ 2 ] = { nullptr , nullptr };      ///< start and goal
 std::uint64_t f_exp[ 2 ] = { 0 , 0 };             ///< expansions by side
 std::size_t f_open[ 2 ] = { 0 , 0 };              ///< open nodes by side
 std::vector< Kid > f_kids;                        ///< the generated ones

 };  // end( class( BiFronts ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS ReverseDomain --------------------------*/
/*--------------------------------------------------------------------------*/
/// a BidirectionalDomain with the start and the goal swapped
/** The search of ReverseDomain< D > goes from the goal of D to its start
 * along the edges reversed: its successors are the predecessors of D, its
 * heuristic the backward heuristic of D, and the other way round, so that
 * A* on it is the reverse A* (rA*) of Siag et al. (ECAI 2024), and a path
 * found on it, reversed, is a path of D. The operators are kept when D is
 * an OperatorDomain and a BackwardOperatorDomain. */

template< BidirectionalDomain D >
class ReverseDomain {

 public:

 using state_type = state_t< D >;               ///< the state type
 using cost_type = cost_t< D >;                 ///< the cost type

 /// constructor, taking the Domain (which must outlive the ReverseDomain)
 explicit ReverseDomain( const D & d ) : f_d( d ) {}

 [[nodiscard]] state_type initial_state( void ) const {
  return( f_d.goal_state() );
  }
 [[nodiscard]] state_type goal_state( void ) const {
  return( f_d.initial_state() );
  }
 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( s == f_d.initial_state() );
  }
 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d.hash( s ) );
  }
 [[nodiscard]] cost_type heuristic( const state_type & s ) const {
  return( f_d.backward_heuristic( s ) );
  }
 [[nodiscard]] cost_type backward_heuristic( const state_type & s ) const {
  return( f_d.heuristic( s ) );
  }
 template< class F >
 void for_each_successor( const state_type & s , cost_type h , F && f )
  const {
  f_d.for_each_predecessor( s , h , std::forward< F >( f ) );
  }
 template< class F >
 void for_each_successor( const state_type & s , cost_type h ,
                          std::uint8_t op , F && f ) const
  requires( BackwardOperatorDomain< D > && OperatorDomain< D > ) {
  f_d.for_each_predecessor( s , h , op , std::forward< F >( f ) );
  }
 template< class F >
 void for_each_predecessor( const state_type & s , cost_type h , F && f )
  const {
  f_d.for_each_successor( s , h , std::forward< F >( f ) );
  }
 template< class F >
 void for_each_predecessor( const state_type & s , cost_type h ,
                            std::uint8_t op , F && f ) const
  requires( BackwardOperatorDomain< D > && OperatorDomain< D > ) {
  f_d.for_each_successor( s , h , op , std::forward< F >( f ) );
  }
 [[nodiscard]] cost_type min_edge_cost( void ) const
  requires( MinCostDomain< D > ) {
  return( f_d.min_edge_cost() );
  }

 private:

 const D & f_d;                                    ///< the Domain

 };  // end( class( ReverseDomain ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BlindDomain ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a BidirectionalDomain with both heuristics replaced by 0
/** A* on BlindDomain< D > is Dijkstra's algorithm, and MM on it the
 * bidirectional brute-force search MM0 of Holte et al. (AAAI 2016). */

template< BidirectionalDomain D >
class BlindDomain {

 public:

 using state_type = state_t< D >;               ///< the state type
 using cost_type = cost_t< D >;                 ///< the cost type

 /// constructor, taking the Domain (which must outlive the BlindDomain)
 explicit BlindDomain( const D & d ) : f_d( d ) {}

 [[nodiscard]] state_type initial_state( void ) const {
  return( f_d.initial_state() );
  }
 [[nodiscard]] state_type goal_state( void ) const {
  return( f_d.goal_state() );
  }
 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( f_d.is_goal( s ) );
  }
 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d.hash( s ) );
  }
 [[nodiscard]] cost_type heuristic( const state_type & ) const {
  return( 0 );
  }
 [[nodiscard]] cost_type backward_heuristic( const state_type & ) const {
  return( 0 );
  }
 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const {
  f_d.for_each_successor( s , cost_type( 0 ) ,
                          [ & ]( const state_type & t , cost_type c ,
                                 cost_type ) {
   f( t , c , cost_type( 0 ) );
   } );
  }
 template< class F >
 void for_each_predecessor( const state_type & s , cost_type , F && f )
  const {
  f_d.for_each_predecessor( s , cost_type( 0 ) ,
                            [ & ]( const state_type & t , cost_type c ,
                                   cost_type ) {
   f( t , c , cost_type( 0 ) );
   } );
  }
 [[nodiscard]] cost_type min_edge_cost( void ) const
  requires( MinCostDomain< D > ) {
  return( f_d.min_edge_cost() );
  }

 private:

 const D & f_d;                                    ///< the Domain

 };  // end( class( BlindDomain ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Bidirectional */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Bidirectional.h ------------------------*/
/*--------------------------------------------------------------------------*/
