/*--------------------------------------------------------------------------*/
/*------------------------------ File main.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The program that the recipes of packaging/ build against the installed
 * library: the 8-puzzle of figure 3.28 of Russell and Norvig, whose optimal
 * solution has 26 moves, solved by the sequential A* and by HDA* on the
 * C++ threads, and a 0-1 knapsack of four items, whose optimum is 9, solved
 * by the sequential branch-and-bound.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/PBeFS.h"

#include <cstddef>

#include <iostream>

#include <limits>

#include <utility>

#include <vector>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- KNAPSACK ---------------------------------*/
/*--------------------------------------------------------------------------*/
// the 0-1 knapsack with the items in the order of their ratio, the bound
// of the fractional knapsack and the feasible value of the items taken

struct Knapsack {

 struct node_type {
  std::size_t k = 0;              // the items fixed
  long weight = 0;                // the weight of those taken
  long profit = 0;                // their profit
  };

 static constexpr Sense sense = Sense::maximize;

 std::vector< long > p{ 6 , 5 , 3 , 1 };   // the profits
 std::vector< long > w{ 3 , 3 , 2 , 1 };   // the weights
 long cap = 5;                             // the capacity

 node_type root( void ) const { return( node_type() ); }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e;
  if( n.weight > cap ) {
   e.infeasible = true;
   return( e );
   }
  long r = cap - n.weight;
  double b = double( n.profit );
  std::size_t j = n.k;
  for( ; ( j < p.size() ) && ( w[ j ] <= r ) ; ++j ) {
   r -= w[ j ];
   b += double( p[ j ] );
   }
  if( j < p.size() )
   b += double( p[ j ] ) * double( r ) / double( w[ j ] );
  e.bound = b;
  e.feasible = double( n.profit );
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  if( n.k == p.size() )
   return;
  // no cheap bound of the children: +infinity, none, for a maximization
  constexpr double none = std::numeric_limits< double >::infinity();
  node_type out = n , in = n;
  ++out.k;
  ++in.k;
  in.weight += w[ n.k ];
  in.profit += p[ n.k ];
  f( std::move( out ) , none );
  f( std::move( in ) , none );
  }

 };

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const EightPuzzle d( EightPuzzle::parse( "7 2 4 5 0 6 8 3 1" ) );
 const auto s = astar( d );
 ParallelOptions o;
 o.threads = 2;
 const ZobristDistribution< EightPuzzle > z( d );
 const auto h = HDAStar< EightPuzzle , decltype( z ) >( d , z , o ).solve();
 const auto b = BranchAndBoundSearch< Knapsack >( Knapsack() ).solve();
 std::cout << "A* " << s.cost << ", HDA* " << h.cost << ", B&B "
           << b.value << std::endl;
 return( ( s.cost == 26 ) && ( h.cost == 26 ) && ( b.value == 9 ) ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File main.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
