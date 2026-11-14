/*--------------------------------------------------------------------------*/
/*------------------------------ File DVCBS.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template DVCBS, the Dynamic Vertex Cover
 * Bidirectional Search of Shperberg, Felner, Sturtevant, Shimony and
 * Hayoun (AAAI 2019).
 *
 * DVCBS keeps the high level of NBS, a lower bound LB on \f$ C^* \f$ raised
 * until the incumbent \f$ U \f$ is not larger, and the queues waiting_D, by
 * \f$ f_D \f$, and ready_D, of the open nodes with \f$ f_D < \mathrm{LB}
 * \f$ (or \f$ \leq \f$ with BidirectionalOptions::eager), but the nodes
 * move from waiting_D to ready_D a value of \f$ f_D \f$ at a time, and
 * those of ready_D are grouped in clusters, one per value of \f$ g_D \f$. The
 * clusters are the vertices of the dynamic must-expand graph DGMX, a
 * forward cluster of \f$ g_F = a \f$ and a backward one of \f$ g_B = b \f$
 * being joined if \f$ a + b + \epsilon \leq \mathrm{LB} \f$, every pair of
 * their nodes having then \f$ \mathrm{lb} \leq \mathrm{LB} \f$; instead of
 * both ends of one edge, as NBS, DVCBS expands at each step one cluster of
 * a minimum vertex cover of DGMX, weighted by the sizes of the clusters
 * [Shperberg et al., algorithm 3].
 *
 * A minimum vertex cover of DGMX is a prefix of the forward clusters by
 * increasing \f$ g \f$ and a prefix of the backward ones (Shaham et al.,
 * SoCS 2017): with the forward clusters \f$ a_0 < a_1 < \cdots < a_{m-1}
 * \f$ and the backward ones \f$ b_0 < \cdots < b_{n-1} \f$ that have an
 * edge, covering the first \f$ i + 1 \f$ forward clusters leaves the edges
 * of \f$ a_{i+1} \f$, the cluster of least g not covered, to the backward
 * clusters \f$ b_j \f$ with \f$ a_{i+1} + b_j + \epsilon \leq \mathrm{LB}
 * \f$, a prefix that covers the edges of every later forward cluster too;
 * so the minimum is found in a pass over \f$ i = -1 , \ldots , m - 1 \f$,
 * linear in the number of clusters. Among the minimum covers, the cluster
 * expanded is the forward \f$ a_0 \f$ or the backward \f$ b_0 \f$ that is
 * in one of them, the one with fewer nodes, the backward on a tie: this is
 * the variant reported by Shperberg et al. (section 5.1), its nodes
 * expanded in the order they became ready. DVCBS has no bound of the kind
 * of NBS on its expansions (section 5.2); it expanded fewer nodes than NBS
 * on every domain of Shperberg et al. (table 1).
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

#ifndef __PBeFS_DVCBS
 #define __PBeFS_DVCBS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <map>

#include <vector>

#include "PBeFS/astar/Bidirectional.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS DVCBS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Dynamic Vertex Cover Bidirectional Search of a BidirectionalDomain
/** The counters of the Solution are those of BiFronts. A node of a cluster
 * whose \f$ f \f$ is not smaller than \f$ U \f$ when its turn comes is
 * closed without being expanded. */

template< BidirectionalDomain D >
class DVCBS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = BiNode< State , Cost >;             ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit DVCBS( const D & domain ,
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

 /// the ready nodes of one g: those whose g and ready flag still match are
 /// live, the others are skipped
 struct Cluster {
  std::vector< NodeT * > nodes;         ///< in the order they became ready
  std::size_t live = 0;                 ///< the nodes still in it
  };

 /// the state of a search
 struct Run {

  const DVCBS & s;                      ///< the search
  BiFronts< D > fr;                     ///< the two fronts
  Queue waiting[ 2 ];                   ///< by f
  std::map< Cost , Cluster > ready[ 2 ];   ///< the clusters, by g
  Cost lb = 0;                          ///< the lower bound LB
  std::vector< NodeT * > batch;         ///< the cluster being expanded

  explicit Run( const DVCBS & search )
   : s( search ) , fr( search.f_domain , search.f_options , false ) {}

  // pushes n in waiting_d; with old its previous g there, it leaves its
  // cluster if it was ready
  void push( unsigned d , NodeT * n , Cost old ) {
   if( n->ready( d ) )
    unready( d , n , old );
   waiting[ d ].push( n , Cost( n->g[ d ] + n->h[ d ] ) , n->g[ d ] );
   }

  // n , ready in d with the g g, leaves its cluster
  void unready( unsigned d , NodeT * n , Cost g ) {
   n->flag( 6 + d , false );
   const auto it = ready[ d ].find( g );
   if( ( it != ready[ d ].end() ) && ( --it->second.live == 0 ) )
    ready[ d ].erase( it );
   }

  // true if n is still in the cluster of g in d
  bool alive( unsigned d , const NodeT * n , Cost g ) const {
   return( n->ready( d ) && ( n->g[ d ] == g ) && ! n->closed( d ) &&
           ! ( s.f_options.consistent && n->closed( 1 - d ) ) );
   }

  // drops the stale entries at the top of waiting_d; true if not empty
  bool clean( unsigned d ) {
   return( waiting[ d ].clean( [ & ]( const auto & e ) {
            return( e.node->ready( d ) ||
                    fr.stale( d , e.node , e.g ) ); } ) );
   }

  // moves the cluster of least f of waiting_d, which is not empty, to
  // the clusters of ready_d
  void move( unsigned d ) {
   const Cost f = waiting[ d ].top().key;
   do {
    NodeT * n = waiting[ d ].pop().node;
    n->flag( 6 + d , true );
    Cluster & c = ready[ d ][ n->g[ d ] ];
    c.nodes.push_back( n );
    ++c.live;
    } while( clean( d ) && ( waiting[ d ].top().key == f ) );
   }

  // the weights of the clusters of d with an edge, given the least g of
  // the other direction
  void clusters( unsigned d , Cost other , std::vector< Cost > & g ,
                 std::vector< std::size_t > & w ) const {
   g.clear();
   w.clear();
   for( const auto & [ k , c ] : ready[ d ] ) {
    if( k + other + fr.epsilon() > lb )
     break;
    g.push_back( k );
    w.push_back( c.live );
    }
   }

  // the direction of the cluster to expand, of a minimum vertex cover of
  // DGMX
  unsigned cover( void ) {
   std::vector< Cost > a , b;
   std::vector< std::size_t > wa , wb;
   clusters( 0 , ready[ 1 ].begin()->first , a , wa );
   clusters( 1 , ready[ 0 ].begin()->first , b , wb );
   const Cost eps = fr.epsilon();
   std::size_t best = ~std::size_t( 0 );
   bool fwd = false , bwd = false;   // a_0 , b_0 in some minimum cover
   std::size_t sa = 0;               // the weight of a_0 .. a_i
   std::vector< std::size_t > pb( b.size() + 1 , 0 );   // prefix weights
   for( std::size_t j = 0 ; j < b.size() ; ++j )
    pb[ j + 1 ] = pb[ j ] + wb[ j ];
   std::size_t nb = b.size();        // backward clusters joined to a_i+1
   for( std::size_t i = 0 ; i <= a.size() ; ++i ) {   // i forward covered
    if( i < a.size() )
     while( ( nb > 0 ) && ( a[ i ] + b[ nb - 1 ] + eps > lb ) )
      --nb;
    const std::size_t k = i < a.size() ? nb : 0;
    const std::size_t wt = sa + pb[ k ];
    if( wt < best ) {
     best = wt;
     fwd = bwd = false;
     }
    if( wt == best ) {
     fwd = fwd || ( i > 0 );
     bwd = bwd || ( k > 0 );
     }
    if( i < a.size() )
     sa += wa[ i ];
    }
   if( fwd && ( ! bwd || ( wa[ 0 ] < wb[ 0 ] ) ) )
    return( forward_side );
   return( backward_side );
   }

  // finds the next cluster; false if a direction has no open node
  bool prepare( unsigned & side ) {
   const bool eager = s.f_options.eager;
   for( unsigned d = 0 ; d < 2 ; ++d )
    while( clean( d ) && ( eager ? waiting[ d ].top().key <= lb
                                 : waiting[ d ].top().key < lb ) )
     move( d );
   const Cost eps = fr.epsilon();
   for( ;; ) {
    bool w[ 2 ];
    for( unsigned d = 0 ; d < 2 ; ++d ) {
     w[ d ] = clean( d );
     if( ready[ d ].empty() && ! w[ d ] )
      return( false );
     }
    const bool r = ! ready[ 0 ].empty() && ! ready[ 1 ].empty();
    if( r && ( ready[ 0 ].begin()->first + ready[ 1 ].begin()->first +
               eps <= lb ) ) {
     side = cover();
     return( true );
     }
    bool moved = false;
    for( unsigned d : { 1u , 0u } )
     if( w[ d ] && ( waiting[ d ].top().key <= lb ) ) {
      move( d );
      moved = true;
      }
    if( moved )
     continue;
    Cost x = BiFronts< D >::inf();
    for( unsigned d = 0 ; d < 2 ; ++d )
     if( w[ d ] )
      x = std::min( x , waiting[ d ].top().key );
    if( r )
     x = std::min( x , Cost( ready[ 0 ].begin()->first +
                             ready[ 1 ].begin()->first + eps ) );
    lb = x;
    }
   }

  // expands the cluster of least g of the direction d
  void expand( unsigned d ) {
   const auto it = ready[ d ].begin();
   const Cost g = it->first;
   batch = it->second.nodes;
   const unsigned o = 1 - d;
   for( NodeT * n : batch ) {
    if( ! alive( d , n , g ) )
     continue;
    unready( d , n , g );
    fr.leave( n , d );
    if( s.f_options.consistent && n->ready( o ) )
     unready( o , n , n->g[ o ] );     // closed here, it leaves the other
    if( n->g[ d ] + n->h[ d ] >= fr.U ) {
     n->flag( d , true );
     continue;
     }
    fr.expand( d , n , [ & ]( NodeT * m , Cost old ) {
     push( d , m , old );
     } );
    if( fr.full() || ( fr.U <= lb ) )        // done: the rest is left
     return;
    }
   // the cluster may keep nodes no longer alive whose leaving was not
   // counted (closed in the other direction): it goes as a whole
   const auto at = ready[ d ].find( g );
   if( at != ready[ d ].end() ) {
    for( NodeT * n : at->second.nodes )
     if( n->ready( d ) && ( n->g[ d ] == g ) )
      n->flag( 6 + d , false );
    ready[ d ].erase( at );
    }
   }

  // the search; true if stopped by the node limit
  bool search( void ) {
   fr.start( [ & ]( unsigned d , NodeT * n ) {
    push( d , n , BiFronts< D >::inf() );
    } );
   for( ;; ) {
    unsigned side = 0;
    if( ! prepare( side ) )
     return( false );
    if( fr.U <= lb )
     return( false );
    expand( side );
    if( fr.full() )
     return( true );
    }
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 BidirectionalOptions f_options;                     ///< the options

 };  // end( class( DVCBS ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with DVCBS

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > dvcbs(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( DVCBS< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_DVCBS */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File DVCBS.h ---------------------------*/
/*--------------------------------------------------------------------------*/
