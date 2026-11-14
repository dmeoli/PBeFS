/*--------------------------------------------------------------------------*/
/*------------------------------ File RBFS.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template RBFS, the Recursive Best-First Search
 * of Korf (AIJ 1993), a best-first search in linear space.
 *
 * Every node \f$ n \f$ on the current path carries a backed-up value
 * \f$ F( n ) \geq f( n ) \f$, the least f of the frontier below it when it
 * was last explored; the search of the subtree of \f$ n \f$ goes on while
 * the best child does not exceed the bound B, the F of the best
 * alternative path elsewhere:
 *
 * \code{.unparsed}
 *  1  function RBFS(n, F(n), B)
 *  2     if f(n) > B then return f(n)
 *  3     if n is a goal then stop with the path to n
 *  4     if n has no successor then return ∞
 *  5     for each successor c of n do
 *  6        F(c) ← max(F(n), f(c)) if f(n) < F(n), else f(c)
 *  7     sort the successors by F
 *  8     while F(c₁) ≤ B and F(c₁) < ∞ do      ▷ c₁, c₂: the two best
 *  9        F(c₁) ← RBFS(c₁, F(c₁), min(B, F(c₂)));  sort again
 * 10     return F(c₁)
 * 11  RBFS(s₀, f(s₀), ∞)
 * \endcode
 *
 * Line 6 is what makes it efficient: a node explored before carries its
 * \f$ F > f \f$ down to its children, which are thus not explored again
 * below that value. With an admissible h the first goal found is optimal,
 * and with a consistent one RBFS expands the nodes in best-first order
 * when it expands a node for the first time (Korf, theorems 4.2 and 4.3);
 * a state reached by several paths is expanded once per path, as in IDA*,
 * but the parent of a node is never generated again on an
 * OperatorDomain. The counter expanded counts the calls that reach line 5.
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

#ifndef __PBeFS_RBFS
 #define __PBeFS_RBFS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS RBFS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// Recursive Best-First Search
/** The search of the file comment; the node limit bounds the expansions,
 * since nothing is stored but the current path and the successors of its
 * nodes. */

template< Domain D >
class RBFS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit RBFS( const D & domain , const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void )
 {
  Stopwatch clock;
  SolutionT sol;
  f_st = SearchStats();
  f_path.clear();
  f_found = false;
  f_limit = false;
  const State s0 = f_domain.initial_state();
  const Cost h0 = Cost( f_domain.heuristic( s0 ) );
  f_path.push_back( s0 );
  search( s0 , Cost( 0 ) , h0 , no_op , h0 , infinite_cost< Cost >() , 0 );
  if( f_found ) {
   sol.status = Status::solved;
   sol.cost = f_cost;
   sol.path = f_path;
   }
  else
   if( f_limit )
    sol.status = Status::node_limit;
  sol.stats = f_st;
  sol.stats.stored = f_peak;
  sol.stats.evaluated = f_st.generated + 1;
  sol.stats.expanded_by = { f_st.expanded };
  sol.stats.stored_by = { f_peak };
  sol.stats.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a successor, with its backed-up value
 struct Kid {
  State s;
  Cost g;
  Cost h;
  std::uint8_t op;
  Cost F;
  };

 /// RBFS( n , F , B ) of the file comment, n being s with g, h and the
 /// operator op; the successors of the nodes at depth d are in v_kids[ d ]
 Cost search( const State & s , Cost g , Cost h , std::uint8_t op , Cost F ,
              Cost B , std::size_t d ) {
  const Cost f = g + h;
  if( f > B )
   return( f );
  if( f_domain.is_goal( s ) ) {
   f_found = true;
   f_cost = g;
   return( f );
   }
  if( f_options.node_limit && ( f_st.expanded >= f_options.node_limit ) ) {
   f_limit = true;
   return( infinite_cost< Cost >() );
   }
  ++f_st.expanded;
  if( v_kids.size() <= d )
   v_kids.resize( d + 1 );
  std::vector< Kid > & k = v_kids[ d ];
  k.clear();
  for_each_child( f_domain , s , h , op ,
                  [ & ]( const State & t , Cost c , Cost ht ,
                         std::uint8_t ot ) {
   if( ht == infinite_cost< Cost >() )
    return;
   const Cost fc = g + c + ht;
   k.push_back( Kid{ t , g + c , ht , ot , ( f < F ) ? std::max( F , fc )
                                                      : fc } );
   } );
  f_st.generated += k.size();
  f_peak = std::max( f_peak , count_stored( d ) );
  if( k.empty() )
   return( infinite_cost< Cost >() );
  const auto by_F = []( const Kid & a , const Kid & b ) {
   return( a.F < b.F );
   };
  std::stable_sort( k.begin() , k.end() , by_F );
  for( ; ; ) {
   Kid & best = v_kids[ d ][ 0 ];
   if( ( best.F > B ) || ( best.F == infinite_cost< Cost >() ) )
    return( best.F );
   const Cost alt = v_kids[ d ].size() > 1 ? v_kids[ d ][ 1 ].F
                                           : infinite_cost< Cost >();
   f_path.push_back( best.s );
   const Kid b = best;
   const Cost r = search( b.s , b.g , b.h , b.op , b.F , std::min( B , alt ) ,
                          d + 1 );
   if( f_found || f_limit )
    return( r );
   f_path.pop_back();
   auto & v = v_kids[ d ];
   v[ 0 ].F = r;
   // the first successor back in its place among the others, sorted
   for( std::size_t i = 1 ; ( i < v.size() ) && ( v[ i ].F < v[ i - 1 ].F ) ;
        ++i )
    std::swap( v[ i ] , v[ i - 1 ] );
   }
  }

 /// the successors held on the current path, down to the depth d
 [[nodiscard]] std::size_t count_stored( std::size_t d ) const {
  std::size_t n = 0;
  for( std::size_t i = 0 ; i <= d ; ++i )
   n += v_kids[ i ].size();
  return( n );
  }

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 SearchStats f_st;                                   ///< the counters
 std::vector< std::vector< Kid > > v_kids;           ///< by depth
 std::vector< State > f_path;                        ///< the current path
 bool f_found = false;                               ///< a goal is found
 bool f_limit = false;                               ///< the limit hit
 Cost f_cost{};                                      ///< its cost
 std::size_t f_peak = 0;                             ///< the most stored

 };  // end( class( RBFS ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with RBFS

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > rbfs(
 const D & domain , const SearchOptions & options = {} )
{
 return( RBFS< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_RBFS */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File RBFS.h -----------------------------*/
/*--------------------------------------------------------------------------*/
