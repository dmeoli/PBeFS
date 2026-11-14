/*--------------------------------------------------------------------------*/
/*--------------------------- File FocalSearch.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template FocalSearch, the bounded suboptimal
 * searches with a focal list: \f$ A^*_\epsilon \f$ of Pearl and Kim (IEEE
 * TPAMI 1982), which expands, among the nodes of OPEN with
 * \f$ f \leq w f_{\min} \f$, the one of a second order, and Explicit
 * Estimation Search (EES) of Thayer and Ruml (IJCAI 2011), which orders the
 * focal list by an inadmissible estimate of the distance to a goal and
 * chooses between three lists.
 *
 * In the focal search, with the notation of AStar.h, \f$ f_{\min} \f$ the
 * least \f$ f \f$ of OPEN and \f$ k \f$ the second order (the distance
 * \f$ d \f$ to a goal, \f$ h \f$, or \f$ g + ( 2 w - 1 ) h \f$):
 *
 * \code{.unparsed}
 *  1  procedure FOCAL(s₀, w)
 *  2     OPEN ← {s₀};  B ← w f(s₀);  FOCAL ← {s₀}
 *  3     while OPEN ≠ ∅ do
 *  4        B' ← w f_min
 *  5        if B' > B then move to FOCAL the nodes of OPEN with B < f ≤ B'
 *  6        B ← B'
 *  7        u ← the first node of FOCAL by k, removed
 *  8        if f(u) > B then put u back among the nodes not in FOCAL;
 *              continue
 *  9        close u;  if u is a goal then return its path
 * 10        for each successor v of u, inserted as A* does, put v also in
 *              FOCAL if f(v) ≤ B
 * \endcode
 *
 * OPEN is kept by \f$ f \f$, with the nodes not in FOCAL also in a list
 * PENDING by \f$ f \f$, from which line 5 takes them in order; all the lists
 * drop their stale entries lazily, and line 8 handles a \f$ f_{\min} \f$
 * that decreases, which only an inconsistent heuristic causes. Every node of
 * OPEN with \f$ f \leq B \f$ is in FOCAL, and the node of \f$ f_{\min} \f$
 * among them, so that FOCAL is not empty at line 7. A goal selected at line
 * 9 has \f$ g = f \leq w f_{\min} \leq w C^* \f$, since the reopening keeps
 * a node of an optimal path with \f$ g = g^* \f$ in OPEN, whose \f$ f \f$ is
 * at most \f$ C^* \f$; \f$ f_{\min} \f$ is the lower bound of the Solution.
 *
 * EES keeps OPEN by \f$ \hat f = g + \hat h \f$, FOCAL with the nodes of
 * \f$ \hat f \leq w \hat f_{\min} \f$ by \f$ \hat d \f$, and CLEANUP by
 * \f$ f \f$; the estimates are corrected online by the mean one-step errors
 * \f$ \bar\epsilon_h \f$ and \f$ \bar\epsilon_d \f$, measured at each
 * expansion on the best child (that of least \f$ f \f$): \f$ \epsilon_h =
 * c( u , v ) + h( v ) - h( u ) \f$, \f$ \epsilon_d = 1 + d( v ) - d( u ) \f$,
 * \f$ \hat d = d / ( 1 - \bar\epsilon_d ) \f$ and \f$ \hat h = h +
 * \bar\epsilon_h \hat d \f$ [see section 3 of the paper], the values when
 * the node is inserted:
 *
 * \code{.unparsed}
 *  1  procedure SELECTNODE()
 *  2     b_d ← the first node of FOCAL;  b_f̂ ← that of OPEN;
 *        b_f ← that of CLEANUP
 *  3     if f̂(b_d) ≤ w f(b_f) then return b_d
 *  4     if f̂(b_f̂) ≤ w f(b_f) then return b_f̂
 *  5     return b_f
 * \endcode
 *
 * FOCAL is kept as in the focal search, on \f$ \hat f \f$ in place of
 * \f$ f \f$. A goal has \f$ \hat h = 0 \f$, since \f$ h = d = 0 \f$ there,
 * so that a goal selected at line 3 or 4 has \f$ g = \hat f \leq w f( b_f )
 * \leq w C^* \f$, and one selected at line 5 is optimal; this holds whatever
 * the estimates, which only change the order.
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

#ifndef __PBeFS_FocalSearch
 #define __PBeFS_FocalSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <limits>

#include <stdexcept>

#include <vector>

#include "ADS/priority/DaryHeap.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/OpenList.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/astar/Weight.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// the order of the focal list of FocalSearch
enum class FocalOrder {
 distance ,    ///< the distance d of a DistanceDomain, h otherwise
 heuristic ,   ///< h
 optimistic ,  ///< \f$ g + ( 2 w - 1 ) h \f$
 explicit_estimation   ///< EES, with \f$ \hat d \f$ and three lists
 };

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS FocalOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of FocalSearch: the weight of the SearchOptions is the
/// bound \f$ w \f$

struct FocalOptions : public SearchOptions {

 FocalOrder order = FocalOrder::distance;     ///< the order of FOCAL

 };  // end( struct( FocalOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS LazyQueue ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a heap of nodes by a key of any type, then by the larger g, that drops
/// its stale entries lazily, as the OPEN lists do [see OpenList.h]; each
/// entry carries a second value, read back with the entry

template< class NodeT , class Key >
class LazyQueue {

 using Cost = decltype( NodeT::g );

 public:

 /// an entry: the key, the second value, the g when pushed, the node
 struct Item {
  Key k;
  double aux;
  Cost g;
  NodeT * node;
  /// true if the node has been closed or improved since it was pushed
  [[nodiscard]] bool stale( void ) const noexcept {
   return( node->closed || ( node->g != g ) );
   }
  };

 /// pushes the node with its current g, the key k and the value aux

 void push( NodeT * n , Key k , double aux = 0 ) {
  f_heap.push( Item{ k , aux , n->g , n } );
  }

 /// drops the stale entries at the top; true if an entry is left

 bool clean( void ) {
  while( ! f_heap.empty() && f_heap.top().stale() )
   static_cast< void >( f_heap.pop() );
  return( ! f_heap.empty() );
  }

 /// the first entry
 [[nodiscard]] const Item & top( void ) const { return( f_heap.top() ); }

 /// removes and returns the first entry
 Item pop( void ) { return( f_heap.pop() ); }

 /// true if there is no entry
 [[nodiscard]] bool empty( void ) const { return( f_heap.empty() ); }

 private:

 /// the order: the smaller key, then the larger g
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.k < b.k ) || ( ! ( b.k < a.k ) && ( a.g > b.g ) ) );
   }
  };

 ADS::DaryHeap< Item , Better , 4 > f_heap;               ///< the heap

 };  // end( class( LazyQueue ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS FocalSearch ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the focal search \f$ A^*_\epsilon \f$ and Explicit Estimation Search
/** The path found costs at most \f$ w C^* \f$, \f$ w \f$ being the weight
 * of the options, and bound is the least \f$ f \f$ of OPEN when the goal is
 * selected. The order of the options chooses the focal search (with
 * \f$ d \f$, \f$ h \f$ or \f$ g + ( 2 w - 1 ) h \f$ as the second order) or
 * EES; \f$ d \f$ is the distance of a DistanceDomain, \f$ h \f$ otherwise.
 * The closed nodes reached through a cheaper path are reopened, and the
 * goal is tested when it is selected. anytime and drop_duplicates of the
 * options are not read; the constructor throws std::invalid_argument if the
 * weight is below 1. */

template< Domain D >
class FocalSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the search)

 explicit FocalSearch( const D & domain , const FocalOptions & options = {} )
  : f_domain( domain ) , f_options( options ) , f_w( options.weight ) ,
    f_wo( 2 * options.weight - 1 ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Run r( *this );
  if( f_options.order == FocalOrder::explicit_estimation )
   r.ees();
  else
   r.focal();
  return( r.finish() );
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the state of a search

 struct Run {

  explicit Run( const FocalSearch & a )
   : f_a( a ) , d( a.f_domain ) ,
     table( a.f_options.expected_nodes
            ? NodeTable< NodeT >::slots_for( a.f_options.expected_nodes )
            : 1024 ) {
   const State s0 = d.initial_state();
   root = table.find_or_insert( s0 , d.hash( s0 ) , [ & ] {
    NodeT * n = pool.allocate();
    *n = NodeT{ s0 , Cost( 0 ) , Cost( d.heuristic( s0 ) ) , nullptr ,
                false , false };
    return( n );
    } ).first;
   }

  // the distance d of the state s, whose heuristic is h
  Cost dist( const State & s , Cost h ) const {
   if constexpr( DistanceDomain< D > ) {
    static_cast< void >( h );
    return( Cost( d.distance( s ) ) );
    }
   else
    return( h );
   }

  // the distance d of the node n
  Cost dist( const NodeT * n ) const { return( dist( n->state , n->h ) ); }

  // the key of n in FOCAL, in the focal search
  Cost focal_key( const NodeT * n ) const {
   switch( f_a.f_options.order ) {
    case FocalOrder::heuristic: return( n->h );
    case FocalOrder::optimistic: return( f_a.f_wo.key( n->g , n->h ) );
    default: return( dist( n ) );
    }
   }

  // closes n, counts its expansion, and generates its successors; false if
  // n is a goal, which is then recorded with the bound lb
  bool expand( NodeT * n , Cost lb ) {
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( d.is_goal( n->state ) ) {
    goal = n;
    bound = std::min( lb , n->g );
    return( false );
    }
   kids.clear();
   for_each_child( d , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    const std::uint64_t hash = d.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    } );
   st.generated += kids.size();
   return( true );
   }

  // the node of the successor k of n if it is new or improved, after the
  // insertion as A* does; nullptr if it is a duplicate
  NodeT * insert( NodeT * n , const Kid & k ) {
   const Cost g = n->g + k.c;
   auto [ m , fresh ] = table.find_or_insert( k.s , k.hash , [ & ] {
    NodeT * x = pool.allocate();
    *x = NodeT{ k.s , g , k.h , n , false , false , k.op };
    return( x );
    } );
   if( fresh )
    return( m );
   if( ! ( g < m->g ) ) {
    ++st.duplicates;
    return( nullptr );
    }
   m->g = g;
   m->parent = n;
   m->op = k.op;
   m->closed = false;
   return( m );
   }

  // true if the node limit is reached
  bool limit( void ) {
   if( f_a.f_options.node_limit &&
       ( pool.size() >= f_a.f_options.node_limit ) )
    limit_hit = true;
   return( limit_hit );
   }

  // the focal search [see the first pseudocode]
  void focal( void ) {
   LazyQueue< NodeT , Cost > open , pend , fq;
   Cost b = f_a.f_w.times( root->f() );
   const auto add = [ & ]( NodeT * m ) {
    open.push( m , m->f() );
    if( m->f() <= b )
     fq.push( m , focal_key( m ) );
    else
     pend.push( m , m->f() );
    };
   add( root );
   while( open.clean() ) {
    const Cost fmin = open.top().k;
    const Cost nb = f_a.f_w.times( fmin );
    if( nb > b )
     while( pend.clean() && ( pend.top().k <= nb ) ) {
      NodeT * m = pend.pop().node;
      fq.push( m , focal_key( m ) );
      }
    b = nb;
    NodeT * n;
    if( fq.clean() ) {
     const auto e = fq.pop();
     n = e.node;
     if( n->f() > b ) {                         // f_min has decreased
      pend.push( n , n->f() );
      continue;
      }
     }
    else
     n = open.pop().node;
    if( ! expand( n , fmin ) )
     return;
    for( const Kid & k : kids )
     if( NodeT * m = insert( n , k ) )
      add( m );
    if( limit() )
     return;
    }
   }

  // EES [see the second pseudocode]
  void ees( void ) {
   LazyQueue< NodeT , double > open , pend , fq;
   LazyQueue< NodeT , Cost > cleanup;
   const double w = f_a.f_w.value();
   double sum_h = 0 , sum_d = 0;
   std::uint64_t samples = 0;
   // f̂ and d̂ of m with the current errors
   const auto hat = [ & ]( const NodeT * m , double & fh , double & dh ) {
    const double eh = samples ? sum_h / double( samples ) : 0;
    double ed = samples ? sum_d / double( samples ) : 0;
    ed = std::clamp( ed , 0.0 , 1.0 - 1e-6 );
    dh = double( dist( m ) ) / ( 1 - ed );
    fh = double( m->g ) + double( m->h ) + std::max( eh , 0.0 ) * dh;
    };
   double b = 0;
   bool first = true;
   const auto add = [ & ]( NodeT * m ) {
    double fh , dh;
    hat( m , fh , dh );
    open.push( m , fh );
    cleanup.push( m , m->f() );
    if( ! first && ( fh <= b ) )
     fq.push( m , dh , fh );
    else
     pend.push( m , fh , dh );
    };
   add( root );
   while( cleanup.clean() ) {
    open.clean();
    const double fhmin = open.top().k;
    const double nb = w * fhmin;
    if( first || ( nb > b ) ) {
     b = nb;
     while( pend.clean() && ( pend.top().k <= b ) ) {
      const auto e = pend.pop();
      fq.push( e.node , e.aux , e.k );
      }
     }
    else
     b = nb;
    first = false;
    while( fq.clean() && ( fq.top().aux > b ) ) {   // f̂_min has decreased
     const auto e = fq.pop();
     pend.push( e.node , e.aux , e.k );
     }
    const Cost fmin = cleanup.top().k;
    const double wf = w * double( fmin );
    NodeT * n;
    if( fq.clean() && ( fq.top().aux <= wf ) )
     n = fq.pop().node;
    else
     if( fhmin <= wf )
      n = open.pop().node;
     else
      n = cleanup.pop().node;
    if( ! expand( n , fmin ) )
     return;
    // the one-step errors on the best child, of least c + h
    const Kid * best = nullptr;
    for( const Kid & k : kids ) {
     if( ! best || ( k.c + k.h < best->c + best->h ) )
      best = & k;
     if( NodeT * m = insert( n , k ) )
      add( m );
     }
    if( best ) {
     sum_h += double( best->c ) + double( best->h ) - double( n->h );
     sum_d += 1 + double( dist( best->s , best->h ) ) - double( dist( n ) );
     ++samples;
     }
    if( limit() )
     return;
    }
   }

  // the Solution
  SolutionT finish( void ) {
   SolutionT sol;
   if( limit_hit )
    sol.status = Status::node_limit;
   else
    if( goal ) {
     sol.status = Status::solved;
     sol.bound = bound;
     sol.path = trace_path( goal );
     sol.cost = std::min( goal->g , Cost( path_cost( d , sol.path ) ) );
     }
   st.stored = pool.size();
   st.expanded_by = { st.expanded };
   st.stored_by = { st.stored };
   st.seconds = clock.seconds();
   sol.stats = std::move( st );
   return( sol );
   }

  const FocalSearch & f_a;                      // the search
  const D & d;                                  // the Domain
  Stopwatch clock;                              // started at the start
  NodePool< NodeT > pool;                       // the nodes
  NodeTable< NodeT > table;                     // CLOSED
  std::vector< Kid > kids;                      // the successors of a node
  SearchStats st;                               // the counters
  NodeT * root = nullptr;                       // the node of s0
  NodeT * goal = nullptr;                       // the goal found
  Cost bound{};                                 // the lower bound
  bool limit_hit = false;                       // the node limit is hit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 FocalOptions f_options;                             ///< the options
 Weight< Cost > f_w;                                 ///< the bound w
 Weight< Cost > f_wo;                                ///< the weight 2w - 1

 };  // end( class( FocalSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_FocalSearch */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FocalSearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
