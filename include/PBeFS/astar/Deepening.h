/*--------------------------------------------------------------------------*/
/*---------------------------- File Deepening.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the bounds of the iterative-deepening searches (IDAStar,
 * IDPSDD and TDS): the options DeepeningOptions and the class template
 * CutoffHistogram, which chooses the next bound from the f of the nodes
 * that an iteration pruned.
 *
 * An iteration of bound \f$ \beta_k \f$ prunes the nodes of
 * \f$ f > \beta_k \f$, and the classical next bound is the least of their
 * f, \f$ \ell_k \f$, so that each iteration adds at least one f value; when
 * the f values are many and close (costs that are not integers, or grids,
 * whose number of nodes grows slowly with the bound), the iterations are
 * as many as the f values and the search expands the nodes of the first
 * ones again at each of them. IDA*-CR (Sarkar, Chakrabarti, Ghose and De
 * Sarkar, AIJ 1991) chooses instead a bound beyond which enough pruned nodes
 * lie: with the growth \f$ r > 1 \f$, \f$ \beta_{k+1} \f$ is the largest
 * f in the first buckets of the histogram of the f of the pruned nodes that
 * hold at least \f$ ( r - 1 ) E_k \f$ of them,
 * \f$ E_k \f$ being the expansions of the iteration, and never below
 * \f$ \ell_k \f$. When even all the pruned nodes are fewer than that and
 * the iteration expanded less than \f$ r \f$ times the one before, as on
 * a grid, where the children of the nodes of the last f layer are all just
 * above it and each layer adds few nodes, the bound grows instead by twice
 * its last step, \f$ \beta_{k+1} = \beta_k + 2 ( \beta_k - \beta_{k-1} )
 * \f$, so that the steps grow geometrically while the iterations stay
 * small; where the iterations grow by \f$ r \f$ anyway, as on the
 * sliding-tile puzzles, the bounds stay the classical ones. The
 * buckets are those of the binary exponent of the
 * excess \f$ f - \beta_k \f$, the excesses in \f$ [ 2^{x-1} , 2^x ) \f$ in
 * the bucket \f$ x \f$, so that the histogram has a fixed size whatever
 * the costs, and its buckets are fine near the bound and coarse far from
 * it.
 *
 * The cost of the first path found then has a bound: an iteration that
 * finds no goal expands every node of \f$ f \le \beta_k \f$, hence the
 * optimum \f$ C^* \f$ is at least \f$ \ell_k \f$, while the path found in
 * the iteration \f$ k + 1 \f$ costs \f$ c \le \beta_{k+1} \f$, so that
 * \f[
 *   c - C^* \le \beta_{k+1} - \ell_k ,
 * \f]
 * which is 0 when \f$ \beta_{k+1} = \ell_k \f$, i.e., for the classical
 * bound. With DeepeningOptions::optimal, the iteration that finds the first
 * goal goes on as a depth-first branch-and-bound, the nodes of f not below
 * the incumbent being pruned: every node of f not above
 * \f$ \min( \beta_{k+1} , z ) \f$ is then expanded, the nodes of an optimal
 * path have \f$ f \le C^* \le c \le \beta_{k+1} \f$, and the incumbent at
 * the end is \f$ C^* \f$; this completes the iteration, and needs no
 * further one.
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

#ifndef __PBeFS_Deepening
 #define __PBeFS_Deepening
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include "PBeFS/astar/Node.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DeepeningOptions -------------------------*/
/*--------------------------------------------------------------------------*/
/// the choice of the bounds of an iterative-deepening search

struct DeepeningOptions {

 /// the growth \f$ r \f$ of IDA*-CR: if larger than 1, the next bound is
 /// the one of CutoffHistogram, otherwise the least f pruned
 double growth = 0;

 /// true to complete the iteration of the first goal as a branch-and-bound
 /// when its bound is above the least f pruned by the iteration before, so
 /// that the path is optimal; false to stop at the first goal, whose cost
 /// is then within \f$ \beta_{k+1} - \ell_k \f$ of the optimum
 bool optimal = true;

 };  // end( struct( DeepeningOptions ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS CutoffHistogram --------------------------*/
/*--------------------------------------------------------------------------*/
/// the f of the nodes pruned by an iteration, and the next bound

template< class Cost >
class CutoffHistogram {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*-------------------------------- CHANGING --------------------------------*/
 /// empties the histogram for the iteration of bound b, which counts the
 /// f in buckets only if buckets is true

 void reset( Cost b , bool buckets ) noexcept {
  f_bound = b;
  f_low = infinite_cost< Cost >();
  f_buckets = buckets;
  if( buckets ) {
   v_n.fill( 0 );
   v_max.fill( Cost( 0 ) );
   }
  }

 /// forgets the last step of the bound, for a new search

 void restart( void ) noexcept {
  f_step = Cost( 0 );
  f_last = 0;
  }

 /// a node of f above the bound

 void add( Cost f ) noexcept {
  f_low = std::min( f_low , f );
  if( f_buckets ) {
   int x;
   static_cast< void >( std::frexp( double( f - f_bound ) , & x ) );
   const auto i = std::size_t( std::clamp( x + Offset , 0 ,
                                           int( Size ) - 1 ) );
   ++v_n[ i ];
   v_max[ i ] = std::max( v_max[ i ] , f );
   }
  }

 /// the nodes of another histogram of the same iteration

 void merge( const CutoffHistogram & o ) noexcept {
  f_low = std::min( f_low , o.f_low );
  if( f_buckets )
   for( std::size_t i = 0 ; i < Size ; ++i ) {
    v_n[ i ] += o.v_n[ i ];
    v_max[ i ] = std::max( v_max[ i ] , o.v_max[ i ] );
    }
  }

/*-------------------------------- READING ---------------------------------*/
 /// the least f pruned, the infinity if none

 [[nodiscard]] Cost low( void ) const noexcept { return( f_low ); }

 /// the next bound: the least f pruned if growth is not above 1, or the
 /// largest f of the first buckets that hold at least ( growth - 1 ) times
 /// expanded nodes pruned, or the bound plus twice the last step if fewer
 /// were pruned and expanded is less than growth times the expansions of
 /// the iteration before, never below the least f pruned; the infinity if
 /// no node was pruned

 Cost next( double growth , std::uint64_t expanded ) {
  if( ! f_buckets || ( growth <= 1 ) ||
      ( f_low == infinite_cost< Cost >() ) )
   return( f_low );
  const double target = ( growth - 1 ) * double( expanded );
  double cum = 0;
  std::size_t last = 0;
  for( std::size_t i = 0 ; i < Size ; ++i ) {
   if( v_n[ i ] == 0 )
    continue;
   last = i;
   cum += double( v_n[ i ] );
   if( cum >= target )
    break;
   }
  Cost edge = v_max[ last ];       // the largest f of the buckets taken
  const bool grew = double( expanded ) >= growth * double( f_last );
  if( ( cum < target ) && ! grew && ( f_step > Cost( 0 ) ) )
   edge = std::max( edge , Cost( f_bound + 2 * f_step ) );
  edge = std::max( f_low , edge );
  f_step = edge - f_bound;
  f_last = expanded;
  return( edge );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 static constexpr int Offset = 1075;                ///< exponent of 0
 static constexpr std::size_t Size = 2100;          ///< all the exponents

 Cost f_bound{};                                    ///< the bound
 Cost f_low = infinite_cost< Cost >();              ///< the least f
 Cost f_step{};                                     ///< the last step
 std::uint64_t f_last = 0;                          ///< its expansions
 bool f_buckets = false;                            ///< counts the f
 std::array< std::uint64_t , Size > v_n{};          ///< the buckets
 std::array< Cost , Size > v_max{};                 ///< their largest f

 };  // end( class( CutoffHistogram ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Deepening */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Deepening.h ---------------------------*/
/*--------------------------------------------------------------------------*/
