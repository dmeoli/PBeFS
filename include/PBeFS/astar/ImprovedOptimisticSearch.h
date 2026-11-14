/*--------------------------------------------------------------------------*/
/*--------------------- File ImprovedOptimisticSearch.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ImprovedOptimisticSearch, the Improved
 * Optimistic Search (IOS) of Chen, Sturtevant, Doyle and Ruml (SoCS 2019):
 * a weighted A* with no reopening finds a first path, and an A* restarted
 * from the root, which reopens each node of the first search at most once,
 * proves it within the bound or improves it.
 *
 * With the bound \f$ w \f$ and the optimistic weight \f$ w_o \f$
 * (\f$ 2 w - 1 \f$ by default), \f$ z \f$ the cost of the first path
 * \f$ P \f$, \f$ r \f$ the largest saving found on a node of \f$ P \f$ and
 * \f$ L \f$ the lower bound, with the notation of AStar.h:
 *
 * \code{.unparsed}
 *  1  procedure IOS(s₀, w, w_o)
 *  2     L ← 0;  weighted A* by g + w_o h, no closed node reopened, each
 *           node selected raising L to (g + w_o h) / w_o, until a goal is
 *           selected:  z ← its g, P ← its path
 *  3     every node is closed;  r ← 0;  push s₀ in OPEN by f
 *  4     while OPEN ≠ ∅ and z − r > w L do
 *  5        u ← the first node of OPEN by f, removed;  L ← max(L, f(u))
 *  6        close u;  if u is a goal then z ← g(u);  r ← 0;  P ← its path;
 *              break
 *  7        for each successor v of u, with g' = g(u) + c(u, v) do
 *  8           if v ∈ P and g' < its g on P then r ← max(r, its g on P − g')
 *  9           if v was not reached since line 3 then
 * 10              g(v) ← g';  parent(v) ← u;  push v in OPEN     ▷ reopened
 * 11           else insert v as A* does
 * 12     return P improved by r, of cost z − r
 * \endcode
 *
 * With a consistent heuristic the weighted A* selects only nodes with
 * \f$ g + w_o h \leq w_o C^* \f$ [see Theorem 1 of Likhachev et al., NIPS
 * 2003, in AnytimeAStar.h], so that the \f$ L \f$ of line 2 is at most
 * \f$ C^* \f$ (the improved termination of the paper, with
 * \f$ f' = g / w_o + h \f$); the A* of line 5, which sees every node anew
 * at line 9, selects only nodes with \f$ f \leq C^* \f$ until it selects a
 * goal, and the test of line 4 stops with a path of cost at most
 * \f$ w C^* \f$. The path of line 12 is the one to the node of \f$ P \f$
 * with the saving \f$ r \f$, followed by the rest of \f$ P \f$; this is
 * the solution-update policy of the paper. Algorithm 1 of the paper goes on
 * instead from the OPEN of the first search, ordered by \f$ f \f$, with no
 * reopening, and proves the bound also by its least \f$ f \f$; we restart
 * from the root because, when the closed nodes are not reopened, a node of
 * an optimal path may stay with \f$ g > g^* \f$, and the least \f$ f \f$
 * of OPEN is then not a lower bound on \f$ C^* \f$ (cf. Corollary 3 of
 * Burns et al., JAIR 2010). The phases are not interleaved, as in the
 * configuration that the paper recommends.
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

#ifndef __PBeFS_ImprovedOptimisticSearch
 #define __PBeFS_ImprovedOptimisticSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <type_traits>

#include <vector>

#include "PBeFS/astar/OptimisticSearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ImprovedOptimisticSearch ---------------------*/
/*--------------------------------------------------------------------------*/
/// the Improved Optimistic Search
/** The options are those of OptimisticSearch: the weight is the bound
 * \f$ w \f$, and optimism the weight \f$ w_o \f$ of the first search
 * (\f$ 2 w - 1 \f$ if 0). The heuristic has to be consistent for the bound
 * to hold [see above]. The first path and each saving are recorded in
 * Solution::improvements; with the node limit, the Solution keeps the best
 * path found, with Status::node_limit. The constructor throws
 * std::invalid_argument if a weight is below 1. */

template< Domain D >
class ImprovedOptimisticSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/

 /// constructor, taking the Domain (which must outlive the search)

 explicit ImprovedOptimisticSearch( const D & domain ,
                                    const OptimisticOptions & options = {} )
  : f_domain( domain ) , f_options( options ) , f_wb( options.weight ) ,
    f_wo( options.optimism == 0 ? 2 * options.weight - 1
                                : options.optimism ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Run r( *this );
  r.run();
  return( r.finish() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a node: that of Node, with the marks of the second search
 struct INode {
  State state;
  Cost g;
  Cost h;
  INode * parent;
  bool closed;
  bool expanded;
  std::uint8_t op = no_op;
  bool seen = false;         ///< reached by the second search
  bool on_path = false;      ///< on the first path P
  Cost g_path{};             ///< its g on P
  [[nodiscard]] Cost f( void ) const noexcept { return( g + h ); }
  };

 using Open = std::conditional_t< std::is_integral_v< Cost > ,
                                  BucketOpenList< INode > ,
                                  HeapOpenList< INode > >;

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

  explicit Run( const ImprovedOptimisticSearch & a )
   : f_a( a ) , d( a.f_domain ) ,
     table( a.f_options.expected_nodes
            ? NodeTable< INode >::slots_for( a.f_options.expected_nodes )
            : 1024 ) {}

  // closes n, counts its expansion and generates its successors
  void expand( INode * n ) {
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   kids.clear();
   for_each_child( d , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    const std::uint64_t hash = d.hash( t );
    table.prefetch( hash );
    kids.push_back( Kid{ t , c , ht , hash , op } );
    } );
   st.generated += kids.size();
   }

  // the node of k, reached from n, and whether it is new
  std::pair< INode * , bool > find( INode * n , const Kid & k , bool seen ) {
   return( table.find_or_insert( k.s , k.hash , [ & ] {
    INode * x = pool.allocate();
    *x = INode{ k.s , n->g + k.c , k.h , n , false , false , k.op , seen };
    return( x );
    } ) );
   }

  // true if the node limit is reached
  bool limit( void ) {
   if( f_a.f_options.node_limit &&
       ( pool.size() >= f_a.f_options.node_limit ) )
    limit_hit = true;
   return( limit_hit );
   }

  // the lower bound that a key of the first search gives: key / w_o
  Cost first_bound( Cost key ) const {
   const Weight< Cost > & wo = f_a.f_wo;
   if constexpr( std::is_integral_v< Cost > ) {
    // key = b g + a h, and ( g + w_o h ) / w_o = key / a, rounded up
    const Cost a = wo.key( Cost( 0 ) , Cost( 1 ) );
    return( Cost( ( key + a - 1 ) / a ) );
    }
   else
    return( key / wo.value() );
   }

  // the cost of the path of the incumbent, improved by the saving
  Cost incumbent( void ) const { return( z - saving ); }

  // records the incumbent
  void record( void ) {
   imp.push_back( { incumbent() , std::min( lb , incumbent() ) ,
                    st.expanded , clock.seconds() } );
   }

  void run( void ) {
   const State s0 = d.initial_state();
   INode * root = table.find_or_insert( s0 , d.hash( s0 ) , [ & ] {
    INode * x = pool.allocate();
    *x = INode{ s0 , Cost( 0 ) , Cost( d.heuristic( s0 ) ) , nullptr ,
                false , false };
    return( x );
    } ).first;

   // the weighted A*, no closed node reopened [line 2]
   const Weight< Cost > & wo = f_a.f_wo;
   {
    Open open;
    open.push( root , Cost( 0 ) , wo.key( Cost( 0 ) , root->h ) );
    INode * goal = nullptr;
    while( ! open.empty() ) {
     const OpenEntry< INode > e = open.pop();
     if( e.stale() )
      continue;
     INode * n = e.node;
     lb = std::max( lb , first_bound( wo.key( e.g , n->h ) ) );
     if( d.is_goal( n->state ) ) {
      n->closed = true;
      ++st.expanded;
      n->expanded = true;
      goal = n;
      break;
      }
     expand( n );
     for( const Kid & k : kids ) {
      auto [ m , fresh ] = find( n , k , false );
      const Cost g = n->g + k.c;
      if( ! fresh ) {
       if( m->closed || ! ( g < m->g ) ) {
        ++st.duplicates;
        continue;
        }
       m->g = g;
       m->parent = n;
       m->op = k.op;
       }
      open.push( m , g , wo.key( g , k.h ) );
      }
     if( limit() )
      return;
     }
    if( ! goal ) {
     done = true;                                // unsolvable
     return;
     }
    z = goal->g;
    for( INode * x = goal ; x ; x = x->parent ) {
     x->on_path = true;
     x->g_path = x->g;
     v_path.push_back( x );
     }
    std::reverse( v_path.begin() , v_path.end() );
    best = goal;
    record();
    }

   // the A* from the root, every node closed [lines 3 to 11]
   Open open;
   root->seen = true;
   root->closed = false;
   open.push( root );
   while( ! open.empty() && ( incumbent() > f_a.f_wb.times( lb ) ) ) {
    const OpenEntry< INode > e = open.pop();
    if( e.stale() )
     continue;
    INode * n = e.node;
    lb = std::max( lb , n->f() );
    if( d.is_goal( n->state ) ) {
     n->closed = true;
     ++st.expanded;
     n->expanded = true;
     if( n->g < incumbent() ) {
      z = n->g;
      saving = Cost( 0 );
      best = n;
      v_path.clear();
      record();
      }
     break;
     }
    expand( n );
    for( const Kid & k : kids ) {
     const Cost g = n->g + k.c;
     auto [ m , fresh ] = find( n , k , true );
     if( fresh ) {
      open.push( m );
      continue;
      }
     if( ! m->seen )                            // reopened, once
      m->seen = true;
     else
      if( ! ( g < m->g ) ) {
       ++st.duplicates;
       continue;
       }
     m->g = g;
     m->parent = n;
     m->op = k.op;
     m->closed = false;
     open.push( m );
     if( m->on_path && ( m->g_path - g > saving ) ) {
      saving = m->g_path - g;
      best = m;
      record();
      }
     }
    if( limit() )
     return;
    }
   done = true;
   }

  // the Solution
  SolutionT finish( void ) {
   SolutionT sol;
   if( limit_hit )
    sol.status = Status::node_limit;
   else
    if( best )
     sol.status = Status::solved;
   if( best ) {
    sol.cost = incumbent();
    sol.bound = std::min( lb , sol.cost );
    sol.improvements = std::move( imp );
    // P itself if no saving was found, since the second search may have
    // given its nodes other parents; else the path to best, and the rest
    // of P after it
    if( ( saving == Cost( 0 ) ) && ( ! v_path.empty() ) &&
        ( best == v_path.back() ) )
     for( const INode * x : v_path )
      sol.path.push_back( x->state );
    else {
     sol.path = trace_path( best );
     const auto it = std::find( v_path.begin() , v_path.end() , best );
     if( it != v_path.end() )
      for( auto x = it + 1 ; x != v_path.end() ; ++x )
       sol.path.push_back( ( *x )->state );
     }
    sol.cost = std::min( sol.cost , Cost( path_cost( d , sol.path ) ) );
    settle_improvements( sol );
    }
   st.stored = pool.size();
   st.expanded_by = { st.expanded };
   st.stored_by = { st.stored };
   st.seconds = clock.seconds();
   sol.stats = std::move( st );
   return( sol );
   }

  const ImprovedOptimisticSearch & f_a;         // the search
  const D & d;                                  // the Domain
  Stopwatch clock;                              // started at the start
  NodePool< INode > pool;                       // the nodes
  NodeTable< INode > table;                     // CLOSED
  std::vector< Kid > kids;                      // the successors of a node
  std::vector< INode * > v_path;                // the first path P
  std::vector< typename SolutionT::Improvement > imp;  // the paths found
  SearchStats st;                               // the counters
  INode * best = nullptr;                       // the node of the saving
  Cost z = infinite_cost< Cost >();             // the cost of P
  Cost saving{};                                // the saving r
  Cost lb{};                                    // the lower bound L
  bool done = false;                            // the search is over
  bool limit_hit = false;                       // the node limit is hit
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 OptimisticOptions f_options;                        ///< the options
 Weight< Cost > f_wb;                                ///< the bound w
 Weight< Cost > f_wo;                                ///< the weight w_o

 };  // end( class( ImprovedOptimisticSearch ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ImprovedOptimisticSearch */

/*--------------------------------------------------------------------------*/
/*------------------- End File ImprovedOptimisticSearch.h ------------------*/
/*--------------------------------------------------------------------------*/
