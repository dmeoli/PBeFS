/*--------------------------------------------------------------------------*/
/*------------------------- File PotentialSearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template PotentialSearch, the Dynamic
 * Potential Search (DPS) of Gilon, Felner and Stern (SoCS 2016), a bounded
 * suboptimal search with no parameter but the bound \f$ w \f$.
 *
 * DPS expands the node of OPEN of largest potential
 * \f$ \mathrm{ud}( n ) = ( w f_{\min} - g( n ) ) / h( n ) \f$,
 * \f$ f_{\min} \f$ being the least \f$ f \f$ of OPEN, i.e., the node most
 * likely to lead to a path of cost at most \f$ w f_{\min} \f$ if the
 * relative error of \f$ h \f$ is uniform; \f$ h( n ) = 0 \f$ is taken
 * as \f$ 10^{-30} \f$, so that such a node comes before all the others if
 * \f$ g( n ) < w f_{\min} \f$, after all of them if \f$ g( n ) >
 * w f_{\min} \f$, and, among those of \f$ h = 0 \f$, by the larger
 * \f$ w f_{\min} - g( n ) \f$. With the notation of AStar.h:
 *
 * \code{.unparsed}
 *  1  procedure DPS(s₀, w)
 *  2     OPEN ← {s₀};  F ← f(s₀)
 *  3     while OPEN ≠ ∅ do
 *  4        if f_min ≠ F then F ← f_min;  order OPEN again by ud
 *  5        u ← the node of OPEN of largest ud, removed
 *  6        close u;  if u is a goal then return its path
 *  7        insert the successors of u as A* does, by ud
 * \endcode
 *
 * OPEN is kept twice: by \f$ f \f$, for \f$ f_{\min} \f$, and by ud, a heap
 * rebuilt at line 4 from its live entries, both dropping their stale
 * entries lazily. The reopening keeps in OPEN a node \f$ n \f$ of an
 * optimal path with \f$ g( n ) = g^* \f$, hence \f$ f_{\min} \leq C^* \f$.
 * A goal \f$ t \f$ selected at line 5 with \f$ g( t ) > w f_{\min} \f$
 * comes after every node of \f$ h > 0 \f$, so that \f$ h( n ) = 0 \f$ and
 * \f$ w f_{\min} - g( n ) \leq w f_{\min} - g( t ) \f$, i.e.,
 * \f$ g( t ) \leq g^*( n ) \leq C^* \f$; in all cases the path costs at
 * most \f$ \max( w f_{\min} , C^* ) \leq w C^* \f$ (cf. section 4.1 of the
 * paper, where the node of largest potential is shown to be in the focal
 * list of the nodes of \f$ f \leq w f_{\min} \f$). The paper gives every
 * node of \f$ h = 0 \f$ an infinite potential (footnotes 4 and 6), which
 * would expand first a goal of \f$ g > w f_{\min} \f$; the order above
 * keeps the bound in this case too. \f$ f_{\min} \f$ is the lower bound of the
 * Solution. The rebuild of line 4 costs \f$ O( |\mathrm{OPEN}| \log
 * |\mathrm{OPEN}| ) \f$ (the paper keeps OPEN in buckets of equal g and h
 * and orders only the buckets, section 4.2), once
 * for each value of \f$ f_{\min} \f$, which with integer costs and a
 * consistent heuristic are at most \f$ C^* - h( s_0 ) + 1 \f$.
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

#ifndef __PBeFS_PotentialSearch
 #define __PBeFS_PotentialSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <limits>

#include <vector>

#include "PBeFS/astar/FocalSearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS PotentialSearch --------------------------*/
/*--------------------------------------------------------------------------*/
/// Dynamic Potential Search
/** The path found costs at most \f$ w C^* \f$, \f$ w \f$ being the weight
 * of the options, and bound is \f$ f_{\min} \f$ when the goal is selected.
 * The closed nodes reached through a cheaper path are reopened. anytime and
 * drop_duplicates of the options are not read; the constructor throws
 * std::invalid_argument if the weight is below 1. */

template< Domain D >
class PotentialSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the search)

 explicit PotentialSearch( const D & domain ,
                           const SearchOptions & options = {} )
  : f_domain( domain ) , f_options( options ) , f_w( options.weight ) {}

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
  LazyQueue< NodeT , Cost > byf;              // OPEN by f
  LazyQueue< NodeT , double > byud;           // OPEN by - ud
  std::vector< NodeT * > live;
  const double w = f_w.value();
  Cost fcur{};
  // minus the potential of m, for the min-heap, h = 0 standing for the
  // tiny h = 10^{-30}
  const auto key = [ & ]( const NodeT * m ) {
   const double room = w * double( fcur ) - double( m->g );
   return( - room / ( m->h == Cost( 0 ) ? 1e-30 : double( m->h ) ) );
   };

  const State s0 = f_domain.initial_state();
  NodeT * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * n = pool.allocate();
   *n = NodeT{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) , nullptr ,
               false , false };
   return( n );
   } ).first;
  fcur = root->f();
  byf.push( root , root->f() );
  byud.push( root , key( root ) );

  std::vector< Kid > kids;
  bool limit_hit = false;
  while( byf.clean() ) {
   const Cost fmin = byf.top().k;
   if( fmin != fcur ) {                        // the potentials change
    fcur = fmin;
    live.clear();
    while( ! byud.empty() ) {
     const auto e = byud.pop();
     if( ! e.stale() )
      live.push_back( e.node );
     }
    for( NodeT * m : live )
     byud.push( m , key( m ) );
    }
   byud.clean();
   NodeT * n = byud.pop().node;
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( f_domain.is_goal( n->state ) ) {
    sol.status = Status::solved;
    sol.bound = std::min( fmin , n->g );
    sol.path = trace_path( n );
    sol.cost = std::min( n->g , Cost( path_cost( f_domain , sol.path ) ) );
    break;
    }
   kids.clear();
   for_each_child( f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    const std::uint64_t hash = f_domain.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    } );
   st.generated += kids.size();
   for( const Kid & k : kids ) {
    const Cost g = n->g + k.c;
    auto [ m , fresh ] = table.find_or_insert( k.s , k.hash , [ & ] {
     NodeT * x = pool.allocate();
     *x = NodeT{ k.s , g , k.h , n , false , false , k.op };
     return( x );
     } );
    if( ! fresh ) {
     if( ! cheaper( g , m->g ) ) {
      ++st.duplicates;
      continue;
      }
     m->g = g;
     m->parent = n;
     m->op = k.op;
     m->closed = false;
     }
    byf.push( m , m->f() );
    byud.push( m , key( m ) );
    }
   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    limit_hit = true;
    break;
    }
   }

  if( limit_hit )
   sol.status = Status::node_limit;
  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a successor: its state, the cost of its edge, its h, its hash
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint64_t hash;
  std::uint8_t op;
  };

 const D & f_domain;                                 ///< the Domain
 SearchOptions f_options;                            ///< the options
 Weight< Cost > f_w;                                 ///< the bound w

 };  // end( class( PotentialSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PotentialSearch */

/*--------------------------------------------------------------------------*/
/*----------------------- End File PotentialSearch.h -----------------------*/
/*--------------------------------------------------------------------------*/
