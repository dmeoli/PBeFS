/*--------------------------------------------------------------------------*/
/*----------------------------- File Weight.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template Weight, the weight \f$ w \geq 1 \f$ of
 * the heuristic of the bounded-suboptimal and anytime searches, and the
 * rules by which they order OPEN, prune with the incumbent and reopen the
 * closed nodes.
 *
 * A weighted search orders OPEN by \f$ f' = g + w h \f$ in place of
 * \f$ f = g + h \f$ [see Pohl, AIJ 1970]. With integer costs the weight is
 * the exact fraction \f$ w = a / b \f$, and the key of a node is
 * \f$ b g + a h \f$, i.e., \f$ b f' \f$, so that no rounding changes the
 * order; scale( z ) puts a cost on the same scale. Let \f$ z \f$ be the
 * incumbent and \f$ C^* \f$ the optimal cost. A weighted search ends when
 * no node of OPEN has \f$ f' < z \f$, and it prunes a node when
 * \f$ w f \geq z \f$: if a node \f$ n \f$ on an optimal path, reached with
 * \f$ g( n ) = g^*( n ) \f$, is in OPEN or is pruned, then
 * \f$ z \leq f'( n ) \leq w f( n ) \leq w C^* \f$, whatever the order of
 * the expansions, so that the bound holds for the parallel searches as well
 * as for the sequential one [see Burns et al., JAIR 2010, for the pruning of
 * wPBNF].
 *
 * With duplicate dropping, a closed node \f$ v \f$ reached from \f$ u \f$
 * with \f$ g' < g( v ) \f$ is reopened only if
 * \f$ g( v ) - g' > ( w - 1 ) ( h( u ) - h( v ) ) \f$, and it keeps its
 * path otherwise. With a consistent heuristic this keeps the bound
 * \f$ w C^* \f$ when the nodes are pruned by \f$ f' \geq z \f$ (not by
 * \f$ w f \f$), in any order of expansion [see Proposition 6.4 of the
 * manual], and it reopens fewer nodes than the weak duplicate dropping of
 * section 5.2 of Burns et al. (JAIR 2010), which keeps the old path when
 * \f$ g( v ) \leq g( u ) + w c( u , v ) \f$, i.e., when the saving is at
 * most \f$ ( w - 1 ) c( u , v ) \f$, while \f$ ( w - 1 ) ( h( u ) - h( v ) )
 * \leq ( w - 1 ) c( u , v ) \f$. With an inconsistent heuristic the bound
 * may fail [see Example 6.5 of the manual], and the duplicates are better
 * reopened.
 *
 * An anytime search goes on after the first goal: it prunes the nodes with
 * \f$ f \geq z \f$ only, and it ends when OPEN is empty, i.e., when the
 * incumbent is optimal [see Hansen and Zhou, JAIR 2007].
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

#ifndef __PBeFS_Weight
 #define __PBeFS_Weight
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstdint>

#include <stdexcept>

#include <type_traits>

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/Solution.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS Weight ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the weight of the heuristic, and the rules of a weighted search
/** With integer costs, the weight is the fraction \f$ a / b \f$ with the
 * least \f$ b \leq 10^4 \f$ equal to it, or \f$ \lfloor 10^4 w \rfloor /
 * 10^4 \f$ if there is none, and the keys \f$ b g + a h \f$ have to fit in
 * the cost type. The weight 1 gives the keys and the rules of A*. */

template< class Cost >
class Weight {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the weight and the modes [see SearchOptions];
 /// throws std::invalid_argument if \f$ w < 1 \f$, or if both modes are
 /// asked, since a dropped duplicate may hold the only optimal path

 explicit Weight( double w = 1 , bool anytime = false , bool drop = false )
 {
  if( ! ( w >= 1 ) || ! std::isfinite( w ) )
   throw( std::invalid_argument( "Weight: the weight must be at least 1" ) );
  if( anytime && drop )
   throw( std::invalid_argument( "Weight: an anytime search cannot drop "
                                 "duplicates" ) );
  if constexpr( std::is_integral_v< Cost > ) {
   f_den = 10000;
   f_num = std::int64_t( std::floor( w * 1e4 ) );
   for( std::int64_t b = 1 ; b <= 10000 ; ++b ) {
    const double x = w * double( b );
    const double r = std::round( x );
    if( std::abs( x - r ) <= 1e-9 * x ) {
     f_num = std::int64_t( r );
     f_den = b;
     break;
     }
    }
   f_w = double( f_num ) / double( f_den );
   }
  else
   f_w = w;
  f_unit = f_w == 1;
  f_anytime = anytime && ! f_unit;
  f_drop = drop && ! f_unit;
  }

 /// constructor, taking the weight and the modes of the options

 explicit Weight( const SearchOptions & o )
  : Weight( o.weight , o.anytime , o.drop_duplicates ) {}

/*-------------------------------- THE KEYS --------------------------------*/
 /// the key of a node in OPEN: \f$ b g + a h \f$ with integer costs,
 /// \f$ g + w h \f$ otherwise

 [[nodiscard]] Cost key( Cost g , Cost h ) const noexcept {
  return( f_unit ? Cost( g + h ) : key_weighted( g , h ) );
  }

 /// key() with a weight above 1

 [[nodiscard]] Cost key_weighted( Cost g , Cost h ) const noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( f_den * g + f_num * h ) );
  else
   return( g + f_w * h );
  }

 /// the cost z on the scale of the keys, the infinite cost staying such

 [[nodiscard]] Cost scale( Cost z ) const noexcept {
  if( f_unit )
   return( z );
  if constexpr( std::is_integral_v< Cost > ) {
   if( z >= infinite_cost< Cost >() / Cost( f_den ) )
    return( infinite_cost< Cost >() );
   return( Cost( f_den * z ) );
   }
  else
   return( z );
  }

/*------------------------------- THE RULES --------------------------------*/
 /// false if no node with the key k, or larger, can lead to a path the
 /// search still looks for, z being the incumbent: \f$ k \geq b z \f$,
 /// unless the search is anytime

 [[nodiscard]] bool live( Cost k , Cost z ) const noexcept {
  return( f_anytime ? ( k < infinite_cost< Cost >() ) : ( k < scale( z ) ) );
  }

 /// true if the incumbent z prunes the node of the given g and h:
 /// \f$ f \geq z \f$ in an anytime search, \f$ f' \geq z \f$ with duplicate
 /// dropping, \f$ w f \geq z \f$ otherwise

 [[nodiscard]] bool prunes( Cost g , Cost h , Cost z ) const noexcept {
  if( f_unit || f_anytime )
   return( g + h >= z );
  return( prunes_weighted( g , h , z ) );
  }

 /// prunes() with a weight above 1 and no anytime search, apart so that
 /// the test of the weight 1 is all that the searches inline

 [[nodiscard]] bool prunes_weighted( Cost g , Cost h , Cost z ) const noexcept
 {
  if( f_drop )
   return( key( g , h ) >= scale( z ) );
  if constexpr( std::is_integral_v< Cost > ) {
   if( z == infinite_cost< Cost >() )
    return( false );
   return( f_num * ( std::int64_t( g ) + h ) >= f_den * std::int64_t( z ) );
   }
  else
   return( f_w * ( g + h ) >= z );
  }

 /// true if a node v of cost g_old, closed or not, reached through a path
 /// of cost g from a node u, gets the new path: if g is cheaper than
 /// \f$ g_{\mathrm{old}} \f$ [see cheaper()] and, when v is closed and
 /// duplicates are dropped,
 /// \f$ g_{\mathrm{old}} - g > ( w - 1 ) ( h( u ) - h( v ) ) \f$; h_u() is
 /// called only in the latter case, the parent being far away maybe

 template< class HU >
 [[nodiscard]] bool improves( Cost g_old , Cost g , bool closed , HU && h_u ,
                              Cost h_v ) const {
  if( ! cheaper( g , g_old ) )
   return( false );
  if( ! ( closed && f_drop ) )
   return( true );
  const Cost hu = h_u();
  if constexpr( std::is_integral_v< Cost > )
   return( f_den * ( std::int64_t( g_old ) - g ) >
           ( f_num - f_den ) * ( std::int64_t( hu ) - h_v ) );
  else
   return( g_old - g > ( f_w - 1 ) * ( hu - h_v ) );
  }

/*------------------------------- THE BOUNDS -------------------------------*/
 /// \f$ w f \f$, rounded down with integer costs: a node of cost
 /// \f$ f' \f$ has \f$ f' \leq w f \f$ if and only if \f$ f' \leq \f$
 /// times( f ); the infinite cost stays such

 [[nodiscard]] Cost times( Cost f ) const noexcept {
  if( f == infinite_cost< Cost >() )
   return( f );
  if constexpr( std::is_integral_v< Cost > ) {
   const std::int64_t x = f_num * std::int64_t( f ) / f_den;
   return( x >= std::int64_t( infinite_cost< Cost >() )
           ? infinite_cost< Cost >() : Cost( x ) );
   }
  else
   return( f_w * f );
  }

 /// the lower bound on \f$ C^* \f$ proven by a search that ended with the
 /// incumbent c: c itself if the search is anytime, which ends only when c
 /// is optimal, lower_bound( c ) otherwise

 [[nodiscard]] Cost proven( Cost c ) const noexcept {
  return( f_anytime ? c : lower_bound( c ) );
  }

 /// the lower bound on \f$ C^* \f$ that a path of cost \f$ c \leq w C^* \f$
 /// gives: \f$ \lceil b c / a \rceil \f$ with integer costs, \f$ c / w \f$
 /// otherwise

 [[nodiscard]] Cost lower_bound( Cost c ) const noexcept {
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( ( f_den * std::int64_t( c ) + f_num - 1 ) / f_num ) );
  else
   return( c / f_w );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the lower bound on \f$ C^* \f$ that a search gives when every path to
 /// a goal has a node of key at least k: \f$ \lceil k / a \rceil \f$ with
 /// integer costs, \f$ k / w \f$ otherwise, since every node of an optimal
 /// path has a key of at most \f$ a C^* \f$ (\f$ w C^* \f$)

 [[nodiscard]] Cost key_lower_bound( Cost k ) const noexcept {
  if( f_unit || ( k >= infinite_cost< Cost >() ) )
   return( k );
  if constexpr( std::is_integral_v< Cost > )
   return( Cost( ( std::int64_t( k ) + f_num - 1 ) / f_num ) );
  else
   return( k / f_w );
  }

/*---------------------------- READING THE DATA ----------------------------*/
 /// the weight, \f$ a / b \f$ with integer costs

 [[nodiscard]] double value( void ) const noexcept { return( f_w ); }

 /// true if the weight is 1

 [[nodiscard]] bool unit( void ) const noexcept { return( f_unit ); }

 /// true if the search is anytime (never with the weight 1, which makes
 /// the first goal optimal)

 [[nodiscard]] bool anytime( void ) const noexcept { return( f_anytime ); }

 /// true if the closed nodes are reopened by improves() only

 [[nodiscard]] bool drops( void ) const noexcept { return( f_drop ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 std::int64_t f_num = 1;                     ///< a, with integer costs
 std::int64_t f_den = 1;                     ///< b, with integer costs
 double f_w = 1;                             ///< the weight
 bool f_anytime = false;                     ///< the search is anytime
 bool f_drop = false;                        ///< duplicates are dropped
 bool f_unit = true;                         ///< the weight is 1

 };  // end( class( Weight ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS UnitWeight ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the weight 1 known at compile time: the interface of Weight, whose
/// rules are those of A* and cost nothing, so that a search written over
/// either compiles, with this one, to the code it would have without
/// weights

template< class Cost >
struct UnitWeight {

 /// \f$ g + h \f$
 [[nodiscard]] static constexpr Cost key( Cost g , Cost h ) noexcept {
  return( g + h );
  }

 /// z
 [[nodiscard]] static constexpr Cost scale( Cost z ) noexcept {
  return( z );
  }

 /// \f$ k < z \f$
 [[nodiscard]] static constexpr bool live( Cost k , Cost z ) noexcept {
  return( k < z );
  }

 /// \f$ g + h \geq z \f$
 [[nodiscard]] static constexpr bool prunes( Cost g , Cost h , Cost z )
  noexcept { return( g + h >= z ); }

 /// g cheaper than \f$ g_{\mathrm{old}} \f$ [see cheaper()]
 template< class HU >
 [[nodiscard]] static constexpr bool improves( Cost g_old , Cost g , bool ,
                                               HU && , Cost ) noexcept {
  return( cheaper( g , g_old ) );
  }

 /// true
 [[nodiscard]] static constexpr bool unit( void ) noexcept { return( true ); }

 };  // end( struct( UnitWeight ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// records in sol.improvements the cost of its path if it is below the
/// last one recorded, which happens when the path to the last goal has
/// become cheaper than its g [see path_cost()]

template< class S >
void settle_improvements( S & sol )
{
 if( ! sol.improvements.empty() &&
     ( sol.cost < sol.improvements.back().cost ) ) {
  auto x = sol.improvements.back();
  x.cost = sol.cost;
  x.bound = std::min( x.bound , sol.cost );
  sol.improvements.push_back( x );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// fills sol with the best goal found by a search on d that may have hit
/// the node limit: solved, with the bound proven by w, if it did not; else
/// node_limit, with the path to the goal and the bound h0 (the h of the
/// root) if the search is anytime; the cost is that of the path [see
/// path_cost()], and if it is below the last incumbent recorded it is
/// recorded too

template< class S , class NodeT , class Cost , class D >
void finish_solution( S & sol , const NodeT * goal , bool limit_hit ,
                      Cost h0 , const Weight< Cost > & w , const D & d )
{
 if( limit_hit )
  sol.status = Status::node_limit;
 if( ! goal || ( limit_hit && ! w.anytime() ) )
  return;
 if( limit_hit )
  sol.bound = h0;
 else {
  sol.status = Status::solved;
  sol.bound = w.proven( goal->g );
  }
 sol.path = trace_path( goal );
 sol.cost = std::min( goal->g , Cost( path_cost( d , sol.path ) ) );
 settle_improvements( sol );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Weight */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File Weight.h ----------------------------*/
/*--------------------------------------------------------------------------*/
