/*--------------------------------------------------------------------------*/
/*--------------------------- File DynamicGraph.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template DynamicGraph, which implements the
 * EventDomain concept [see Domain.h] for a directed graph whose arcs change
 * their costs, whose nodes are blocked and freed and whose goals change,
 * as the gridworlds of the experiments of Koenig, Likhachev and Furcy (AIJ
 * 2004, section 7), and tells each change as an event to the incremental
 * searches.
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

#ifndef __PBeFS_DynamicGraph
 #define __PBeFS_DynamicGraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <cstdlib>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS DynamicGraph ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a directed graph whose arc costs, blocked nodes and goals change
/** The nodes are \f$ 0 , \ldots , n - 1 \f$, the arcs are numbered in the
 * order in which they are added, and an arc costs its own cost unless one
 * of its ends is blocked, when it is absent (of cost \f$ \infty \f$). The
 * heuristic is the Manhattan distance to the nearest goal times scale() if
 * the nodes have positions, given by set_position(), and a value per node
 * otherwise (0 unless set_heuristic() gives it); since an incremental
 * search computes the h of a state once, a goal set changing under the
 * Manhattan distances asks for the fixed set of set_heuristic_goals(), the
 * union of the goals of every search, so that h stays consistent and 0 on
 * the goals of each of them. pairwise_heuristic() is
 * the Manhattan distance between two nodes times scale(). Every change of
 * set_cost(), set_blocked() and set_goal() is recorded, coalesced per arc
 * and per goal, and drain_events() hands the net changes to an incremental
 * search [see EventDomain]. grid() builds the four-connected gridworlds
 * with arcs both ways between adjacent cells, gridworld() and perturb()
 * those of the experiments of Koenig, Likhachev and Furcy and their
 * changes. */

template< class Cost = int >
class DynamicGraph {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint32_t;                     ///< a node
 using cost_type = Cost;                               ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor of n nodes, no arc and no goal, the start being s

 DynamicGraph( std::size_t n , state_type s )
  : f_start( s ) , v_out( n ) , v_in( n ) , v_blocked( n , 0 ) ,
    v_goal( n , 0 ) , v_gpend( n , -1 ) , v_h( n , Cost( 0 ) ) {
  if( s >= n )
   throw( std::invalid_argument( "DynamicGraph: no such start" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the four-connected gridworld of w by h cells, the cell ( x , y ) being
 /// the node \f$ y w + x \f$ at the position ( x , y ), with an arc of cost
 /// c from each cell to each adjacent one, added up, left, right and down
 /// for each cell in increasing order, the start s and no goal

 [[nodiscard]] static DynamicGraph grid( int w , int h , Cost c ,
                                         state_type s = 0 ) {
  DynamicGraph g( std::size_t( w ) * std::size_t( h ) , s );
  static constexpr int dx[ 4 ] = { 0 , -1 , 1 , 0 };
  static constexpr int dy[ 4 ] = { -1 , 0 , 0 , 1 };
  for( int y = 0 ; y < h ; ++y )
   for( int x = 0 ; x < w ; ++x ) {
    const state_type u = state_type( y * w + x );
    g.set_position( u , x , y );
    for( int k = 0 ; k < 4 ; ++k ) {
     const int nx = x + dx[ k ] , ny = y + dy[ k ];
     if( ( nx >= 0 ) && ( ny >= 0 ) && ( nx < w ) && ( ny < h ) )
      g.add_edge( u , state_type( ny * w + nx ) , c );
     }
    }
  return( g );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a gridworld of the experiments of Koenig, Likhachev and Furcy (AIJ
 /// 2004, section 7): grid( w , w ) with the start and the goal drawn
 /// uniformly from the cells by gen, the arcs of cost 1 or 2 drawn
 /// uniformly if block is 0, and of cost 1 with each cell but those two
 /// blocked with probability block otherwise [see random_below() and
 /// random_unit()]; the changes of the construction are not recorded

 template< class Gen >
 [[nodiscard]] static DynamicGraph gridworld( int w , double block ,
                                              Gen & gen ) {
  const auto n = std::uint64_t( w ) * std::uint64_t( w );
  const auto s = state_type( random_below( gen , n ) );
  const auto t = state_type( random_below( gen , n ) );
  auto g = grid( w , w , Cost( 1 ) , s );
  for( std::size_t e = 0 ; e < g.edges() ; ++e )
   g.v_cost[ e ] = block > 0 ? Cost( 1 )
                             : Cost( 1 + random_below( gen , 2 ) );
  g.set_goal( t , true );
  for( state_type v = 0 ; v < n ; ++v )
   if( ( v != s ) && ( v != t ) && ( random_unit( gen ) < block ) )
    g.set_blocked( v , true );
  g.drain_events( []( const auto & ) {} );
  return( g );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the changes of an episode of those experiments: if blocks is 0, a
 /// fraction frac of the arcs, drawn with replacement, given a cost of 1
 /// or 2; otherwise blocks free cells blocked and as many blocked cells
 /// freed, never the start nor a goal (there must be that many)

 template< class Gen >
 void perturb( double frac , int blocks , Gen & gen ) {
  if( blocks == 0 ) {
   const auto k = std::size_t( frac * double( edges() ) + 0.5 );
   for( std::size_t j = 0 ; j < k ; ++j )
    set_cost( std::size_t( random_below( gen , edges() ) ) ,
              Cost( 1 + random_below( gen , 2 ) ) );
   return;
   }
  for( const bool want : { true , false } )
   for( int j = 0 ; j < blocks ; ) {
    const auto v = state_type( random_below( gen , size() ) );
    if( ( blocked( v ) == want ) || ( v == f_start ) || is_goal( v ) )
     continue;
    set_blocked( v , want );
    ++j;
    }
  }

/*------------------------------ THE CHANGES -------------------------------*/
 /// adds the arc ( u , v ) of cost c, returning its index; not recorded as
 /// a change

 std::size_t add_edge( state_type u , state_type v , Cost c ) {
  if( ( u >= size() ) || ( v >= size() ) || ( c < 0 ) )
   throw( std::invalid_argument( "DynamicGraph::add_edge: bad arc" ) );
  const std::size_t e = v_tail.size();
  v_tail.push_back( u );
  v_head.push_back( v );
  v_cost.push_back( c );
  v_epend.push_back( -1 );
  v_out[ u ].push_back( e );
  v_in[ v ].push_back( e );
  return( e );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the cost of the arc e to c, \f$ \infty \f$ removing it

 void set_cost( std::size_t e , Cost c ) {
  const Cost from = cost( e );
  v_cost[ e ] = c;
  changed( e , from );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// blocks the node s (b true) or frees it: the arcs into it and out of it
 /// are absent while it is blocked

 void set_blocked( state_type s , bool b ) {
  if( bool( v_blocked[ s ] ) == b )
   return;
  std::vector< std::pair< std::size_t , Cost > > before;
  for( const std::size_t e : v_out[ s ] )
   before.emplace_back( e , cost( e ) );
  for( const std::size_t e : v_in[ s ] )
   before.emplace_back( e , cost( e ) );
  v_blocked[ s ] = b;
  for( const auto & [ e , from ] : before )
   changed( e , from );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// makes t a goal (b true) or no longer one

 void set_goal( state_type t , bool b ) {
  if( bool( v_goal[ t ] ) == b )
   return;
  if( v_gpend[ t ] < 0 ) {
   v_gpend[ t ] = int( v_gev.size() );
   v_gev.push_back( { t , ! b } );
   }
  v_goal[ t ] = b;
  if( b )
   v_goals.push_back( t );
  else
   v_goals.erase( std::find( v_goals.begin() , v_goals.end() , t ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( EdgeCostChanged ) for each arc whose cost differs from the
 /// one at the last call, and f( GoalChanged ) for each node that became
 /// or stopped being a goal since then, and forgets the changes

 template< class F >
 void drain_events( F && f ) {
  for( const auto & [ e , from ] : v_eev ) {
   v_epend[ e ] = -1;
   const Cost to = cost( e );
   if( to != from )
    f( EdgeCostChanged< state_type , Cost >{ v_tail[ e ] , v_head[ e ] ,
                                             from , to } );
   }
  v_eev.clear();
  for( const auto & [ t , was ] : v_gev ) {
   v_gpend[ t ] = -1;
   if( bool( v_goal[ t ] ) != was )
    f( GoalChanged< state_type >{ t , bool( v_goal[ t ] ) } );
   }
  v_gev.clear();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gives the node s the position ( x , y ) of the Manhattan distances

 void set_position( state_type s , int x , int y ) {
  if( v_x.empty() ) {
   v_x.assign( size() , 0 );
   v_y.assign( size() , 0 );
   }
  v_x[ s ] = x;
  v_y[ s ] = y;
  }

 /// the factor of the Manhattan distances, 1 by default

 void set_scale( Cost c ) { f_scale = c; }

 /// makes s the start, which is not a change that an incremental search
 /// can repair: a search from scratch, as A* from the robot of D* Lite

 void set_start( state_type s ) {
  if( s >= size() )
   throw( std::invalid_argument( "DynamicGraph: no such start" ) );
  f_start = s;
  }

 /// the Manhattan distances to the nodes t instead of the goals

 void set_heuristic_goals( std::vector< state_type > t ) {
  v_href = std::move( t );
  f_href = true;
  }

 /// gives the node s the heuristic h, when the nodes have no position

 void set_heuristic( state_type s , Cost h ) { v_h.at( s ) = h; }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/
 /// the start

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 /// true if s is a goal

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( v_goal[ s ] );
  }

 /// the hash of a node: its index, mixed

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 /// the Manhattan distance to the nearest goal times the scale, or the
 /// value given

 [[nodiscard]] Cost heuristic( state_type s ) const noexcept {
  if( v_x.empty() )
   return( v_h[ s ] );
  const auto & ts = f_href ? v_href : v_goals;
  Cost best = ts.empty() ? Cost( 0 ) : infinite_cost< Cost >();
  for( const state_type t : ts )
   best = std::min( best , pairwise_heuristic( s , t ) );
  return( best );
  }

 /// the Manhattan distance between s and t times the scale, 0 if the
 /// nodes have no position

 [[nodiscard]] Cost pairwise_heuristic( state_type s , state_type t ) const
  noexcept {
  if( v_x.empty() )
   return( Cost( 0 ) );
  return( Cost( f_scale * Cost( std::abs( v_x[ s ] - v_x[ t ] ) +
                                std::abs( v_y[ s ] - v_y[ t ] ) ) ) );
  }

 /// calls f( v , c , h( v ) ) for each arc ( s , v ) present, of cost c

 template< class F >
 void for_each_successor( state_type s , Cost , F && f ) const {
  for( const std::size_t e : v_out[ s ] ) {
   const Cost c = cost( e );
   if( c != infinite_cost< Cost >() )
    f( v_head[ e ] , c , heuristic( v_head[ e ] ) );
   }
  }

 /// calls f( u , c ) for each arc ( u , s ) present, of cost c

 template< class F >
 void for_each_predecessor( state_type s , F && f ) const {
  for( const std::size_t e : v_in[ s ] ) {
   const Cost c = cost( e );
   if( c != infinite_cost< Cost >() )
    f( v_tail[ e ] , c );
   }
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/
 /// the number of features: one per node

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( size() );
  }

 /// calls f( s ), the only feature of s

 template< class F >
 void for_each_feature( state_type s , F && f ) const {
  f( std::size_t( s ) );
  }

 /// the perfect hash: the node

 [[nodiscard]] std::uint64_t perfect_hash( state_type s ) const noexcept {
  return( s );
  }

/*-------------------------------- QUERIES ---------------------------------*/
 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( v_out.size() );
  }

 /// the number of arcs

 [[nodiscard]] std::size_t edges( void ) const noexcept {
  return( v_tail.size() );
  }

 /// the tail of the arc e

 [[nodiscard]] state_type tail( std::size_t e ) const {
  return( v_tail[ e ] );
  }

 /// the head of the arc e

 [[nodiscard]] state_type head( std::size_t e ) const {
  return( v_head[ e ] );
  }

 /// the cost of the arc e, \f$ \infty \f$ if absent or an end is blocked

 [[nodiscard]] Cost cost( std::size_t e ) const {
  if( v_blocked[ v_tail[ e ] ] || v_blocked[ v_head[ e ] ] )
   return( infinite_cost< Cost >() );
  return( v_cost[ e ] );
  }

 /// whether s is blocked

 [[nodiscard]] bool blocked( state_type s ) const { return( v_blocked[ s ] ); }

 /// the goals

 [[nodiscard]] const std::vector< state_type > & goals( void ) const {
  return( v_goals );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the arc e changed, its cost having been from before
 void changed( std::size_t e , Cost from ) {
  if( v_epend[ e ] < 0 ) {
   v_epend[ e ] = int( v_eev.size() );
   v_eev.emplace_back( e , from );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 state_type f_start;                                   ///< the start
 Cost f_scale = Cost( 1 );                     ///< of the Manhattan distance
 std::vector< state_type > v_tail;                     ///< the tails
 std::vector< state_type > v_head;                     ///< the heads
 std::vector< Cost > v_cost;                           ///< the costs
 std::vector< std::vector< std::size_t > > v_out;      ///< arcs out
 std::vector< std::vector< std::size_t > > v_in;       ///< arcs in
 std::vector< char > v_blocked;                        ///< blocked nodes
 std::vector< char > v_goal;                           ///< goal nodes
 std::vector< state_type > v_goals;                    ///< the goals
 std::vector< int > v_epend;                   ///< place in v_eev, or -1
 std::vector< int > v_gpend;                   ///< place in v_gev, or -1
 std::vector< std::pair< std::size_t , Cost > > v_eev; ///< arcs changed
 std::vector< std::pair< state_type , bool > > v_gev;  ///< goals changed
 std::vector< int > v_x;                               ///< the columns
 std::vector< int > v_y;                               ///< the rows
 std::vector< Cost > v_h;                      ///< the heuristic, by node
 std::vector< state_type > v_href;             ///< the goals of h
 bool f_href = false;                          ///< v_href given

 };  // end( class( DynamicGraph ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_DynamicGraph */

/*--------------------------------------------------------------------------*/
/*------------------------ End File DynamicGraph.h -------------------------*/
/*--------------------------------------------------------------------------*/
