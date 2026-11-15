/*--------------------------------------------------------------------------*/
/*-------------------------- File BiBestFirst.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template BiBestFirst, the front-to-end
 * bidirectional searches that expand one node at a time, the one of least
 * priority in the direction chosen: the bidirectional Dijkstra's algorithm,
 * MM of Holte et al. (AAAI 2016), its variant MMe of Sharon et al. (SoCS
 * 2016), the fractional MM, fMM( p ), of Shaham et al. (SoCS 2017), and
 * BAE* of Sadhukhan (2013), as defined by Alcázar et al. (AAAI 2020).
 *
 * With the notation of the page "Bidirectional search" of the
 * documentation, for a node n open in the direction D (F or B, \f$ \bar D
 * \f$ being the other one), the priorities are
 *
 * \f[ \mathrm{pr}_D( n ) = g_D( n ) \f]
 *
 * for the bidirectional Dijkstra's algorithm,
 * \f$ \mathrm{pr}_D( n ) = \max\{ f_D( n ) , 2 g_D( n ) \} \f$ for MM,
 * \f$ \max\{ f_D( n ) , 2 g_D( n ) + \epsilon \} \f$ for MMe,
 * \f$ \max\{ f_F( n ) , g_F( n ) / p \} \f$ and
 * \f$ \max\{ f_B( n ) , g_B( n ) / ( 1 - p ) \} \f$ for fMM( p ), and
 *
 * \f[ b_D( n ) = f_D( n ) + g_D( n ) - h_{\bar D}( n ) \f]
 *
 * for BAE*, the f of the node plus the error \f$ g_D - h_{\bar D} \f$ of
 * the heuristic of the other direction at it. The search is
 *
 * \code{.unparsed}
 *  1  procedure BiBestFirst(s, t)
 *  2     U ← ∞; OPEN_F ← {s}; OPEN_B ← {t}
 *  3     while OPEN_F ≠ ∅ and OPEN_B ≠ ∅ do
 *  4        if U ≤ LB then return U
 *  5        D ← the direction chosen
 *  6        n ← the node of least pr_D in OPEN_D, of largest g_D among those
 *  7        if f_D(n) < U then expand n in the direction D, and
 *              U ← min(U, g_D(m) + g_D̄(m)) for each m generated
 *  8     return U
 * \endcode
 *
 * with the lower bound LB on \f$ C^* \f$ of each algorithm: \f$ g_{\min F}
 * + g_{\min B} + \epsilon \f$ for Dijkstra's algorithm, \f$ \max\{ C ,
 * f_{\min F} , f_{\min B} , g_{\min F} + g_{\min B} + \epsilon \} \f$, \f$ C
 * \f$ being the least priority of the two open lists, for MM, MMe and fMM
 * [see Holte et al., AIJ 2017, section 4], and \f$ ( b_{\min F} +
 * b_{\min B} ) / 2 \f$ for BAE*, rounded up with integer costs, which is a
 * bound only if both heuristics are consistent [see Alcázar et al., AAAI
 * 2020]. With BidirectionalOptions::all_bounds the bound is the largest of
 * all of them that hold under the heuristics, with also \f$ f_{\min F} +
 * d_{\min B} \f$ and \f$ d_{\min F} + f_{\min B} \f$, \f$ d_D = g_D -
 * h_{\bar D} \f$, under consistent heuristics (the individual bounds of
 * Alcázar et al.). The direction of MM, MMe, fMM
 * and Dijkstra's algorithm is the one of least priority, the forward one on
 * a tie unless the backward top has a smaller g; BAE* alternates the two,
 * as in its definition and in PEM-BAE* [see Siag et al., ECAI 2024]; the
 * other choices of SideChoice are available to every algorithm.
 *
 * With SideChoice::priority (their default), MM never expands a node with
 * \f$ g_D( n ) > C^* / 2 \f$, and fMM( p ) none with
 * \f$ g_F( n ) > p C^* \f$ or \f$ g_B( n ) > ( 1 - p ) C^* \f$: every node
 * expanded has \f$ \mathrm{pr} \leq C^* \f$, since the bound \f$ C \f$
 * is not larger than \f$ C^* \f$; with the other choices of the side the
 * node expanded need not be one of least priority, and the property does
 * not hold. The proofs are in chapter 7 of the manual.
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

#ifndef __PBeFS_BiBestFirst
 #define __PBeFS_BiBestFirst
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <stdexcept>

#include <type_traits>

#include "PBeFS/astar/Bidirectional.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the priority, and the bound, of a BiBestFirst

//! [bidirectional-sig-bipriority]
enum class BiPriority {
 dijkstra ,   ///< g, the bidirectional Dijkstra's algorithm
 mm ,         ///< max( f , 2 g ), MM
 mme ,        ///< max( f , 2 g + epsilon ), MMe
 fmm ,        ///< max( f , g / p ) and max( f , g / ( 1 - p ) ), fMM( p )
 bae          ///< f + g - h of the other direction, BAE*
 };
//! [bidirectional-sig-bipriority]

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BiBestFirst ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a front-to-end bidirectional search that expands one node at a time
/** The counters of the Solution are those of BiFronts. BAE* throws
 * std::invalid_argument if BidirectionalOptions::consistent is false, and
 * fMM if the fraction is not in ( 0 , 1 ). */

template< BidirectionalDomain D , BiPriority P >
class BiBestFirst {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = BiNode< State , Cost >;             ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// the type of the priority: real for fMM, the cost otherwise
 using Key = std::conditional_t< P == BiPriority::fmm , double , Cost >;

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit BiBestFirst( const D & domain ,
                       const BidirectionalOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  if( ( P == BiPriority::bae ) && ! options.consistent )
   throw( std::invalid_argument( "BiBestFirst: BAE* needs consistent "
                                 "heuristics" ) );
  if( ( P == BiPriority::fmm ) &&
      ! ( ( options.fraction > 0 ) && ( options.fraction < 1 ) ) )
   throw( std::invalid_argument( "BiBestFirst: the fraction of fMM is not "
                                 "in ( 0 , 1 )" ) );
  }

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

 /// the queues of the bounds
 enum { by_f = 0 , by_g = 1 , by_pr = 2 , by_b = 3 , by_d = 4 };

 using Queue = BiQueue< NodeT , Key , Cost >;
 using CostQueue = BiQueue< NodeT , Cost , Cost >;

 /// the state of a search
 struct Run {

  const BiBestFirst & s;                ///< the search
  bool all;                             ///< all the bounds
  bool need[ 5 ];                       ///< which bound queues are kept
  BiFronts< D > fr;                     ///< the two fronts
  Queue open[ 2 ];                      ///< by priority
  CostQueue aux[ 5 ][ 2 ];              ///< by f, g, pr of MMe, b, d
  SideChoice side;                      ///< the choice of the direction
  unsigned turn = backward_side;        ///< the last side, to alternate

  explicit Run( const BiBestFirst & search )
   : s( search ) , all( search.f_options.all_bounds ) ,
     need{ all || mm_like() , all || mm_like() ,
           all && ( P != BiPriority::mme ) ,
           all && search.f_options.consistent && ( P != BiPriority::bae ) ,
           all && search.f_options.consistent } ,
     fr( search.f_domain , search.f_options ,
         ( P == BiPriority::bae ) || need[ by_d ] ) ,
     side( search.f_options.side ) {
   if( side == SideChoice::automatic )
    side = P == BiPriority::bae ? SideChoice::alternate
                                : SideChoice::priority;
   }

  static constexpr bool mm_like( void ) {
   return( ( P == BiPriority::mm ) || ( P == BiPriority::mme ) ||
           ( P == BiPriority::fmm ) );
   }

  // the priority of n in the direction d
  Key priority( unsigned d , const NodeT * n ) const {
   const Cost g = n->g[ d ] , h = n->h[ d ];
   if constexpr( P == BiPriority::dijkstra )
    return( g );
   else if constexpr( P == BiPriority::mm )
    return( std::max( Cost( g + h ) , Cost( 2 * g ) ) );
   else if constexpr( P == BiPriority::mme )
    return( std::max( Cost( g + h ) , Cost( 2 * g + fr.epsilon() ) ) );
   else if constexpr( P == BiPriority::fmm ) {
    const double p = d == forward_side ? s.f_options.fraction
                                       : 1 - s.f_options.fraction;
    return( std::max( double( g + h ) , double( g ) / p ) );
    }
   else
    return( Cost( 2 * g + h - n->h[ 1 - d ] ) );
   }

  // pushes n in the queues of the direction d
  void push( unsigned d , NodeT * n ) {
   const Cost g = n->g[ d ];
   open[ d ].push( n , priority( d , n ) , g );
   if( need[ by_f ] )
    aux[ by_f ][ d ].push( n , Cost( g + n->h[ d ] ) , g );
   if( need[ by_g ] )
    aux[ by_g ][ d ].push( n , g , g );
   if( need[ by_pr ] )
    aux[ by_pr ][ d ].push( n , std::max( Cost( g + n->h[ d ] ) ,
                                          Cost( 2 * g + fr.epsilon() ) ) ,
                            g );
   if( need[ by_b ] )
    aux[ by_b ][ d ].push( n , Cost( 2 * g + n->h[ d ] - n->h[ 1 - d ] ) ,
                           g );
   if( need[ by_d ] )
    aux[ by_d ][ d ].push( n , Cost( g - n->h[ 1 - d ] ) , g );
   }

  // the least key of the bound queue q of the direction d (after it has
  // been cleaned)
  template< class Q >
  bool least( Q & q , unsigned d , Cost & v ) {
   if( ! q.clean( [ & ]( const auto & e ) {
          return( fr.stale( d , e.node , e.g ) ); } ) )
    return( false );
   v = Cost( q.top().key );
   return( true );
   }

  // the lower bound on C*, given the open lists not empty
  Cost bound( void ) {
   const Cost eps = fr.epsilon();
   Cost lb = 0;
   const Key kf = open[ 0 ].top().key , kb = open[ 1 ].top().key;
   if constexpr( P == BiPriority::dijkstra )
    lb = BiFronts< D >::add( BiFronts< D >::add( Cost( kf ) , Cost( kb ) ) ,
                             eps );
   else if constexpr( P == BiPriority::bae )
    lb = half_up( Cost( kf ) , Cost( kb ) );
   else {
    if constexpr( P == BiPriority::fmm ) {
     const double c = std::min( kf , kb );
     lb = std::is_integral_v< Cost > ? Cost( std::ceil( c - 1e-9 ) )
                                     : Cost( c );
     }
    else
     lb = std::min( kf , kb );
    }
   Cost v[ 2 ];
   if( need[ by_f ] && least( aux[ by_f ][ 0 ] , 0 , v[ 0 ] ) &&
       least( aux[ by_f ][ 1 ] , 1 , v[ 1 ] ) )
    lb = std::max( { lb , v[ 0 ] , v[ 1 ] } );
   if( need[ by_g ] && least( aux[ by_g ][ 0 ] , 0 , v[ 0 ] ) &&
       least( aux[ by_g ][ 1 ] , 1 , v[ 1 ] ) )
    lb = std::max( lb , Cost( v[ 0 ] + v[ 1 ] + eps ) );
   if( need[ by_pr ] && least( aux[ by_pr ][ 0 ] , 0 , v[ 0 ] ) &&
       least( aux[ by_pr ][ 1 ] , 1 , v[ 1 ] ) )
    lb = std::max( lb , std::min( v[ 0 ] , v[ 1 ] ) );
   if( need[ by_b ] && least( aux[ by_b ][ 0 ] , 0 , v[ 0 ] ) &&
       least( aux[ by_b ][ 1 ] , 1 , v[ 1 ] ) )
    lb = std::max( lb , half_up( v[ 0 ] , v[ 1 ] ) );
   Cost w[ 2 ];
   if( need[ by_d ] && least( aux[ by_d ][ 0 ] , 0 , v[ 0 ] ) &&
       least( aux[ by_d ][ 1 ] , 1 , v[ 1 ] ) &&
       least( aux[ by_f ][ 0 ] , 0 , w[ 0 ] ) &&
       least( aux[ by_f ][ 1 ] , 1 , w[ 1 ] ) )
    lb = std::max( { lb , Cost( w[ 0 ] + v[ 1 ] ) ,
                     Cost( v[ 0 ] + w[ 1 ] ) } );
   return( lb );
   }

  // ( a + b ) / 2, rounded up with integer costs
  static Cost half_up( Cost a , Cost b ) {
   if constexpr( std::is_integral_v< Cost > )
    return( Cost( ( a + b + 1 ) / 2 ) );
   else
    return( ( a + b ) / 2 );
   }

  // the direction of the next expansion
  unsigned choose( void ) {
   switch( side ) {
    case SideChoice::alternate:
     turn = 1 - turn;
     return( turn );
    case SideChoice::cardinality:
     return( fr.open_count( 1 ) < fr.open_count( 0 ) ? backward_side
                                                    : forward_side );
    default: {
     const auto a = open[ 0 ].top() , b = open[ 1 ].top();
     if( b.key < a.key )
      return( backward_side );
     if( ( b.key == a.key ) && ( b.g < a.g ) )
      return( backward_side );
     return( forward_side );
     }
    }
   }

  // the search; true if stopped by the node limit
  bool search( void ) {
   fr.start( [ & ]( unsigned d , NodeT * n ) { push( d , n ); } );
   for( ;; ) {
    bool alive = true;
    for( unsigned d = 0 ; d < 2 ; ++d )
     alive = open[ d ].clean( [ & ]( const auto & e ) {
              return( fr.stale( d , e.node , e.g ) ); } ) && alive;
    if( ! alive )
     return( false );
    if( ( fr.U < BiFronts< D >::inf() ) && ( fr.U <= bound() ) )
     return( false );
    const unsigned d = choose();
    NodeT * n = open[ d ].pop().node;
    fr.leave( n , d );
    if( n->g[ d ] + n->h[ d ] >= fr.U ) {      // no cheaper path through n
     n->flags = std::uint8_t( n->flags | ( 1u << d ) );
     continue;
     }
    fr.expand( d , n , [ & ]( NodeT * m , Cost ) { push( d , m ); } );
    if( fr.full() )
     return( true );
    }
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 BidirectionalOptions f_options;                     ///< the options

 };  // end( class( BiBestFirst ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with the bidirectional Dijkstra's algorithm, the
/// heuristics being ignored (the MM0 of Holte et al.)

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > bidijkstra(
 const D & domain , const BidirectionalOptions & options = {} )
{
 const BlindDomain< D > blind( domain );
 return( BiBestFirst< BlindDomain< D > , BiPriority::dijkstra >(
          blind , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with MM

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > mm(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( BiBestFirst< D , BiPriority::mm >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with MMe

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > mme(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( BiBestFirst< D , BiPriority::mme >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with fMM( p ), p being options.fraction

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > fmm(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( BiBestFirst< D , BiPriority::fmm >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with BAE*

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > bae(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( BiBestFirst< D , BiPriority::bae >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BiBestFirst */

/*--------------------------------------------------------------------------*/
/*------------------------ End File BiBestFirst.h --------------------------*/
/*--------------------------------------------------------------------------*/
