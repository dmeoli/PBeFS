/*--------------------------------------------------------------------------*/
/*------------------------------ File GBFHS.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template GBFHS, the Generalized Breadth-First
 * Heuristic Search of Barley, Riddle, Linares López, Dobson and Pohl (SoCS
 * 2018), and for the bidirectional breadth-first search, which is GBFHS
 * without the heuristics.
 *
 * GBFHS raises a bound \f$ f_{Lim} \f$ on \f$ C^* \f$ one unit at a time
 * and, at each value, splits \f$ f_{Lim} - \epsilon + 1 \f$ between two
 * limits \f$ g_{Lim F} + g_{Lim B} \f$ on the depths of the two searches;
 * it then expands in the direction D every open node with
 * \f$ f_D \leq f_{Lim} \f$ and \f$ g_D < g_{Lim D} \f$, until none is left
 * or the incumbent \f$ U \f$ is not larger than \f$ f_{Lim} \f$:
 *
 * \code{.unparsed}
 *  1  procedure GBFHS(s, t)
 *  2     U ← ∞; f_Lim ← max(h_F(s), h_B(t), ε); g_LimF ← g_LimB ← 0
 *  3     while OPEN_F ≠ ∅ and OPEN_B ≠ ∅ do
 *  4        if U ≤ f_Lim then return U
 *  5        raise g_LimF and g_LimB so that g_LimF + g_LimB = f_Lim − ε + 1
 *  6        while some n in OPEN_D has f_D(n) ≤ f_Lim and g_D(n) < g_LimD do
 *  7           expand the one of least f_D, of largest g_D among those,
 *              of the direction of the least such f
 *  8           if U ≤ f_Lim then return U
 *  9        f_Lim ← f_Lim + 1
 * \endcode
 *
 * The costs are integers, \f$ \epsilon \geq 0 \f$ the least of them. After
 * the level \f$ f_{Lim} \f$ every path P of cost \f$ C \leq f_{Lim} \f$ has
 * been found: if the last node of P expanded forward is followed by
 * \f$ u \f$ and the first expanded backward preceded by \f$ v \f$, with
 * \f$ u \f$ before \f$ v \f$, then \f$ g_F( u ) \geq g_{Lim F} \f$,
 * \f$ g_B( v ) \geq g_{Lim B} \f$ and \f$ C \geq g_F( u ) + \epsilon +
 * g_B( v ) \geq f_{Lim} + 1 \f$, a contradiction, so the two searches meet
 * on P; with consistent heuristics every node on it is expanded with its
 * optimal g, since the nodes are expanded by increasing f, as by A*, and a
 * cheaper path to a node has all its nodes of an f smaller than that of
 * the node by consistency, and of a g within the limit. Hence
 * \f$ f_{Lim} \f$ stays a lower bound on \f$ C^* \f$, and GBFHS never
 * expands a node with \f$ f_D > C^* \f$ (Barley et al., section 4).
 *
 * The order of line 7 is ours: the expandable nodes of a level may be
 * expanded in any order, and the order of A* (least f, then largest g)
 * reaches the meeting point within the last level much earlier than the
 * order of g; on the 15-puzzle it expands 2.5 times fewer nodes. The nodes
 * of f within \f$ f_{Lim} \f$ but too deep wait in a queue by g until the
 * limit of their direction passes them.
 *
 * The split of line 5 is BidirectionalOptions::side: alternate (the
 * default) gives each new unit to the smaller limit, the backward one on a
 * tie, after the first split in halves; cardinality gives it to the
 * direction with fewer open nodes, Pohl's criterion; priority gives the
 * forward search the share fraction of the sum, as fMM( p ) does.
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

#ifndef __PBeFS_GBFHS
 #define __PBeFS_GBFHS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <type_traits>

#include "PBeFS/astar/Bidirectional.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS GBFHS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Generalized Breadth-First Heuristic Search of a BidirectionalDomain
/// with integer costs
/** The counters of the Solution are those of BiFronts. */

template< BidirectionalDomain D >
 requires std::is_integral_v< cost_t< D > >
class GBFHS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = BiNode< State , Cost >;             ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit GBFHS( const D & domain ,
                 const BidirectionalOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Stopwatch clock;
  Run r( *this );
  const bool stopped = r.search();
  return( r.fr.finish( stopped , clock.seconds() ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Queue = BiQueue< NodeT , Cost , Cost >;

 /// the state of a search
 struct Run {

  const GBFHS & s;                      ///< the search
  BiFronts< D > fr;                     ///< the two fronts
  Queue waiting[ 2 ];                   ///< by f, those above f_Lim
  Queue ready[ 2 ];                     ///< by f, those within f_Lim
  Queue deep[ 2 ];                      ///< by g, those beyond g_Lim
  Cost flim = 0;                        ///< the bound f_Lim
  Cost glim[ 2 ] = { 0 , 0 };           ///< the depth limits
  Cost eps = 0;                         ///< the least cost of an edge
  bool first = true;                    ///< no split done yet

  explicit Run( const GBFHS & search )
   : s( search ) , fr( search.f_domain , search.f_options , false ) ,
     eps( std::max( Cost( 0 ) , fr.epsilon() ) ) {}

  void push( unsigned d , NodeT * n ) {
   waiting[ d ].push( n , Cost( n->g[ d ] + n->h[ d ] ) , n->g[ d ] );
   }

  bool clean( Queue & q , unsigned d ) {
   return( q.clean( [ & ]( const auto & e ) {
            return( fr.stale( d , e.node , e.g ) ); } ) );
   }

  // moves to ready_d the nodes with f within f_Lim; true if d has an
  // open node
  bool settle( unsigned d ) {
   while( clean( waiting[ d ] , d ) && ( waiting[ d ].top().key <= flim ) ) {
    const auto e = waiting[ d ].pop();
    ready[ d ].push( e.node , e.key , e.g );
    }
   while( clean( deep[ d ] , d ) && ( deep[ d ].top().key < glim[ d ] ) ) {
    const auto e = deep[ d ].pop();
    ready[ d ].push( e.node , Cost( e.g + e.node->h[ d ] ) , e.g );
    }
   return( clean( ready[ d ] , d ) || clean( deep[ d ] , d ) ||
           ! waiting[ d ].empty() );
   }

  // raises the limits to the sum t
  void split( Cost t ) {
   Cost & a = glim[ 0 ];
   Cost & b = glim[ 1 ];
   Cost x = t - a - b;
   if( x <= 0 )
    return;
   switch( s.f_options.side ) {
    case SideChoice::priority: {
     const Cost fa = Cost( std::floor( s.f_options.fraction * double( t ) ) );
     const Cost na = std::clamp( fa , a , Cost( t - b ) );
     x -= na - a;
     a = na;
     b += x;
     break;
     }
    case SideChoice::cardinality:
     while( x-- > 0 )
      ++( fr.open_count( 1 ) < fr.open_count( 0 ) ? b : a );
     break;
    default:
     if( first ) {
      a += x / 2;
      b += x - x / 2;
      }
     else
      while( x-- > 0 )
       ++( a < b ? a : b );
    }
   first = false;
   }

  // the search; true if stopped by the node limit
  bool search( void ) {
   fr.start( [ & ]( unsigned d , NodeT * n ) { push( d , n ); } );
   const D & d = s.f_domain;
   flim = std::max( { Cost( d.heuristic( d.initial_state() ) ) ,
                      Cost( d.backward_heuristic( d.goal_state() ) ) ,
                      eps } );
   for( ;; ) {
    if( fr.U <= flim )
     return( false );
    if( ! settle( 0 ) || ! settle( 1 ) )
     return( false );
    split( Cost( flim - eps + 1 ) );
    for( ;; ) {
     settle( 0 );
     settle( 1 );
     bool ok[ 2 ];
     for( unsigned k = 0 ; k < 2 ; ++k ) {
      while( clean( ready[ k ] , k ) &&
             ( ready[ k ].top().g >= glim[ k ] ) ) {
       const auto e = ready[ k ].pop();
       deep[ k ].push( e.node , e.g , e.g );
       }
      ok[ k ] = clean( ready[ k ] , k );
      }
     if( ! ok[ 0 ] && ! ok[ 1 ] )
      break;
     const unsigned k = ! ok[ 0 ] ? 1 : ! ok[ 1 ] ? 0 :
                        ( ready[ 1 ].top().key < ready[ 0 ].top().key
                          ? 1 : 0 );
     NodeT * n = ready[ k ].pop().node;
     fr.leave( n , k );
     fr.expand( k , n , [ & ]( NodeT * m , Cost ) { push( k , m ); } );
     if( fr.full() )
      return( true );
     if( fr.U <= flim )
      return( false );
     }
    ++flim;
    }
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 BidirectionalOptions f_options;                     ///< the options

 };  // end( class( GBFHS ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with GBFHS

template< BidirectionalDomain D >
 requires std::is_integral_v< cost_t< D > >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > gbfhs(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( GBFHS< D >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with integer costs by the bidirectional breadth-first
/// search: GBFHS with no heuristic, each unit of depth given to the
/// direction with fewer open nodes [Pohl, 1971]

template< BidirectionalDomain D >
 requires std::is_integral_v< cost_t< D > >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > bibfs(
 const D & domain , BidirectionalOptions options = {} )
{
 options.side = SideChoice::cardinality;
 const BlindDomain< D > blind( domain );
 return( GBFHS< BlindDomain< D > >( blind , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_GBFHS */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File GBFHS.h ---------------------------*/
/*--------------------------------------------------------------------------*/
