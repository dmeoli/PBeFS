/*--------------------------------------------------------------------------*/
/*------------------------------- File NBS.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template NBS, the Near-Optimal Bidirectional
 * Search of Chen, Holte, Zilles and Sturtevant (IJCAI 2017), which expands
 * the two nodes of a pair at a time.
 *
 * For a node u open in the forward direction and a node v open in the
 * backward one,
 *
 * \f[ \mathrm{lb}( u , v ) = \max\{ f_F( u ) , f_B( v ) , g_F( u ) +
 *     g_B( v ) + \epsilon \} \f]
 *
 * is a lower bound on the cost of a path from the start through u and v to
 * the goal, and NBS expands both nodes of a pair of least lb, of least
 * \f$ g_F( u ) \f$ and then of least \f$ g_B( v ) \f$ among those, until
 * the least lb is not smaller than the incumbent \f$ U \f$ [Chen et al.,
 * algorithm 1]. With consistent heuristics, Theorem 1 of Eckerle et al.
 * (ICAPS 2017) says that one of the two states of every pair with
 * \f$ \mathrm{lb} < C^* \f$ must be expanded by every admissible
 * front-to-end algorithm: these pairs are the edges of the must-expand
 * graph \f$ G_{MX} \f$, whose minimum vertex cover is a lower bound on the
 * expansions, and NBS, expanding both ends of an edge, expands at most
 * twice its size (Chen et al., theorem 3), no admissible front-to-end
 * algorithm doing better in the worst case (theorem 5).
 *
 * The pair is found as in section 7 of Chen et al., with a lower bound
 * \f$ C_{lb} \f$ on the least lb raised until it is reached: in each
 * direction D, a queue waiting_D of the open nodes by \f$ f_D \f$, the
 * larger \f$ g_D \f$ first on a tie, and a queue ready_D of those with
 * \f$ f_D \leq C_{lb} \f$ by \f$ g_D \f$, the smaller \f$ f_D \f$ first on
 * a tie:
 *
 * \code{.unparsed}
 *  1  procedure PrepareBest
 *  2     move from waiting_D to ready_D every node with f_D < C_lb
 *  3     loop
 *  4        if ready_D ∪ waiting_D = ∅ for some D then return none
 *  5        if ready_F.g + ready_B.g + ε ≤ C_lb then return the two tops
 *  6        if waiting_D.f ≤ C_lb for some D then
 *              move the top of waiting_D to ready_D (for each such D)
 *  7        else C_lb ← min(waiting_F.f, waiting_B.f,
 *                           ready_F.g + ready_B.g + ε)
 * \endcode
 *
 * and the search stops when \f$ U \leq C_{lb} \f$. With
 * BidirectionalOptions::eager, line 2 moves the nodes with
 * \f$ f_D \leq C_{lb} \f$, which is NBS_A of Shperberg et al. (AAAI 2019),
 * whose pairs cover those with \f$ \mathrm{lb} \leq C^* \f$, needed to find
 * all the optimal paths, and line 6 is then never taken. The
 * \f$ \epsilon \f$ is that of BidirectionalOptions; with 0, the bound is the
 * one of Chen et al. without it (the base case of Shperberg et al.).
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

#ifndef __PBeFS_NBS
 #define __PBeFS_NBS
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include "PBeFS/astar/Bidirectional.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS NBS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Near-Optimal Bidirectional Search of a BidirectionalDomain
/** The counters of the Solution are those of BiFronts. A node of a pair
 * whose \f$ f \f$ is not smaller than \f$ U \f$ when its turn comes is
 * closed without being expanded. */

template< BidirectionalDomain D >
class NBS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = BiNode< State , Cost >;             ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit NBS( const D & domain , const BidirectionalOptions & options = {} )
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

 /// an entry of ready: a node, its g and its f when it became ready
 struct Ready {
  NodeT * node;
  Cost g;
  Cost f;
  };

 /// the order of ready: the smaller g, then the smaller f
 struct Earlier {
  bool operator()( const Ready & a , const Ready & b ) const noexcept {
   return( ( a.g < b.g ) || ( ( a.g == b.g ) && ( a.f < b.f ) ) );
   }
  };

 /// the ready queue, a heap
 using ReadyQueue = ADS::DaryHeap< Ready , Earlier >;

 /// the state of a search
 struct Run {

  const NBS & s;                        ///< the search
  BiFronts< D > fr;                     ///< the two fronts
  Queue waiting[ 2 ];                   ///< by f
  ReadyQueue ready[ 2 ];                ///< by g, then f
  Cost clb = 0;                         ///< the lower bound C_lb

  explicit Run( const NBS & search )
   : s( search ) , fr( search.f_domain , search.f_options , false ) {}

  // pushes n in waiting_d
  void push( unsigned d , NodeT * n ) {
   waiting[ d ].push( n , Cost( n->g[ d ] + n->h[ d ] ) , n->g[ d ] );
   }

  // drops the stale entries at the tops; true if q is not empty
  bool clean( Queue & q , unsigned d ) {
   return( q.clean( [ & ]( const auto & e ) {
            return( fr.stale( d , e.node , e.g ) ); } ) );
   }

  bool clean( ReadyQueue & q , unsigned d ) {
   while( ! q.empty() && fr.stale( d , q.top().node , q.top().g ) )
    static_cast< void >( q.pop() );
   return( ! q.empty() );
   }

  // moves the top of waiting_d to ready_d
  void move( unsigned d ) {
   const auto e = waiting[ d ].pop();
   ready[ d ].push( Ready{ e.node , e.g , e.key } );
   }

  // finds the next pair; false if a direction has no open node
  bool prepare( void ) {
   const bool eager = s.f_options.eager;
   for( unsigned d = 0 ; d < 2 ; ++d )
    while( clean( waiting[ d ] , d ) &&
           ( eager ? waiting[ d ].top().key <= clb
                   : waiting[ d ].top().key < clb ) )
     move( d );
   const Cost eps = fr.epsilon();
   for( ;; ) {
    bool r[ 2 ] , w[ 2 ];
    for( unsigned d = 0 ; d < 2 ; ++d ) {
     r[ d ] = clean( ready[ d ] , d );
     w[ d ] = clean( waiting[ d ] , d );
     if( ! r[ d ] && ! w[ d ] )
      return( false );
     }
    if( r[ 0 ] && r[ 1 ] &&
        ( ready[ 0 ].top().g + ready[ 1 ].top().g + eps <= clb ) )
     return( true );
    bool moved = false;
    for( unsigned d : { 1u , 0u } )
     if( w[ d ] && ( waiting[ d ].top().key <= clb ) ) {
      move( d );
      moved = true;
      }
    if( moved )
     continue;
    Cost lb = BiFronts< D >::inf();
    for( unsigned d = 0 ; d < 2 ; ++d )
     if( w[ d ] )
      lb = std::min( lb , waiting[ d ].top().key );
    if( r[ 0 ] && r[ 1 ] )
     lb = std::min( lb , Cost( ready[ 0 ].top().g +
                               ready[ 1 ].top().g + eps ) );
    clb = lb;
    }
   }

  // expands the top of ready_d, unless its f is not below U
  void expand( unsigned d ) {
   if( ! clean( ready[ d ] , d ) )
    return;
   NodeT * n = ready[ d ].pop().node;
   fr.leave( n , d );
   if( n->g[ d ] + n->h[ d ] >= fr.U ) {
    n->flag( d , true );
    return;
    }
   fr.expand( d , n , [ & ]( NodeT * m , Cost ) { push( d , m ); } );
   }

  // the search; true if stopped by the node limit
  bool search( void ) {
   fr.start( [ & ]( unsigned d , NodeT * n ) { push( d , n ); } );
   for( ;; ) {
    if( ! prepare() )
     return( false );
    if( fr.U <= clb )
     return( false );
    expand( forward_side );
    expand( backward_side );
    if( fr.full() )
     return( true );
    }
   }

  };  // end( struct( Run ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 BidirectionalOptions f_options;                     ///< the options

 };  // end( class( NBS ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with NBS

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > nbs(
 const D & domain , const BidirectionalOptions & options = {} )
{
 return( NBS< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_NBS */

/*--------------------------------------------------------------------------*/
/*------------------------------ End File NBS.h ----------------------------*/
/*--------------------------------------------------------------------------*/
