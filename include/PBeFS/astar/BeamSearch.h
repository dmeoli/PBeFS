/*--------------------------------------------------------------------------*/
/*--------------------------- File BeamSearch.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template BeamSearch, the breadth-first beam
 * search: a breadth-first search that keeps, of each layer, only the B
 * best nodes by f (or by h), so that its memory grows linearly with the
 * depth [see Lowerre, 1976, and Zhou and Hansen, ICAPS 2005, section 2].
 *
 * \code{.unparsed}
 *  1  procedure BEAM(s₀, B)
 *  2     layer ← ⟨s₀⟩;  CLOSED ← {s₀}
 *  3     while layer ≠ ⟨⟩ do
 *  4        next ← ⟨⟩
 *  5        for each node n of layer, each successor c of n do
 *  6           if c is a goal then return the path to c
 *  7           if c ∉ CLOSED or it is reached through a cheaper path then
 *                 CLOSED ← CLOSED ∪ {c};  append c to next
 *  8        layer ← the B best nodes of next, by key and then by order
 * \endcode
 *
 * The search is neither complete nor optimal, since the best path may
 * leave the beam; the nodes dropped at line 8 stay in CLOSED, as in the
 * beam search of Zhou and Hansen, so that a later layer does not pay them
 * again. With an infinite width it is a breadth-first search that
 * reopens a node reached by a cheaper path.
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

#ifndef __PBeFS_BeamSearch
 #define __PBeFS_BeamSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <vector>

#include "PBeFS/astar/GreedySearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BeamOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of BeamSearch

//! [memory-bounded-sig-beamoptions]
struct BeamOptions : public SearchOptions {

 std::size_t width = 100;     ///< the nodes kept in a layer, B
 bool greedy = false;         ///< the key is h instead of f

 };  // end( struct( BeamOptions ) )
//! [memory-bounded-sig-beamoptions]

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BeamSearch ----------------------------*/
/*--------------------------------------------------------------------------*/
/// breadth-first beam search

template< Domain D >
class BeamSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit BeamSearch( const D & domain , const BeamOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  if( f_options.width == 0 )
   f_options.width = 1;
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  NodePool< NodeT > pool;
  NodeTable< NodeT > table( f_options.expected_nodes
                            ? NodeTable< NodeT >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  const State s0 = f_domain.initial_state();
  NodeT * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * x = pool.allocate();
   *x = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               false , false , no_op , 0 };
   return( x );
   } ).first;
  const NodeT * goal = f_domain.is_goal( s0 ) ? root : nullptr;
  std::vector< NodeT * > layer{ root } , next;
  std::vector< std::pair< NodeT * , Cost > > cand;  // node, g on arrival
  f_layers = 0;
  while( ( ! goal ) && ! layer.empty() ) {
   ++f_layers;
   cand.clear();
   for( NodeT * n : layer ) {
    if( n->closed )
     continue;
    n->closed = true;
    ++st.expanded;
    if( n->expanded )
     ++st.reexpanded;
    n->expanded = true;
    for_each_child( f_domain , n->state , n->h , n->op ,
                    [ & ]( const State & t , Cost c , Cost ht ,
                           std::uint8_t op ) {
     if( goal )
      return;
     ++st.generated;
     if( ht == infinite_cost< Cost >() )
      return;
     const Cost g = n->g + c;
     bool fresh = false;
     NodeT * m = table.find_or_insert( t , f_domain.hash( t ) , [ & ] {
      fresh = true;
      NodeT * x = pool.allocate();
      *x = NodeT{ t , g , ht , n , false , false , op , 0 };
      return( x );
      } ).first;
     if( ! fresh ) {
      if( ! cheaper( g , m->g ) ) {
       ++st.duplicates;
       return;
       }
      m->g = g;
      m->parent = n;
      m->op = op;
      m->closed = false;
      }
     if( f_domain.is_goal( t ) ) {
      goal = m;
      return;
      }
     cand.emplace_back( m , g );
     } );
    if( goal )
     break;
    }
   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   // the best B, by key and then by order of generation
   next.clear();
   for( auto & [ m , g ] : cand )
    if( m->g == g )
     next.push_back( m );
   next.erase( std::unique( next.begin() , next.end() ) , next.end() );
   const auto key = [ & ]( const NodeT * x ) {
    return( f_options.greedy ? x->h : x->g + x->h );
    };
   if( next.size() > f_options.width ) {
    std::stable_sort( next.begin() , next.end() ,
                      [ & ]( const NodeT * a , const NodeT * b ) {
                       return( key( a ) < key( b ) ); } );
    next.resize( f_options.width );
    }
   layer.swap( next );
   }
  if( goal ) {
   sol.status = Status::solved;
   sol.path = trace_path( goal );
   // a node of the path may have got a cheaper g after the goal
   sol.cost = std::min( goal->g , Cost( path_cost( f_domain , sol.path ) ) );
   }
  st.stored = pool.size();
  st.evaluated = st.stored;
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the layers of the last search

 [[nodiscard]] std::size_t layers( void ) const noexcept {
  return( f_layers );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using NodeT = GreedyNode< State , Cost >;         ///< the node type

 const D & f_domain;                                 ///< the Domain
 BeamOptions f_options;                              ///< the options
 mutable std::size_t f_layers = 0;                   ///< layers

 };  // end( class( BeamSearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with BeamSearch

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > beam_search(
 const D & domain , const BeamOptions & options = {} )
{
 return( BeamSearch< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_BeamSearch */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BeamSearch.h --------------------------*/
/*--------------------------------------------------------------------------*/
