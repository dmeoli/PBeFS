/*--------------------------------------------------------------------------*/
/*--------------------------- File WidthSearch.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the width-based searches of Lipovetzky and Geffner on a
 * FeatureDomain: the class template NoveltyTable, the tuples of features
 * seen so far, and the class templates IteratedWidth (IW(k)),
 * SerializedWidth (SIW) and BestFirstWidth (BFWS).
 *
 * The novelty \f$ w( s ) \f$ of a state \f$ s \f$, generated after the
 * states \f$ s_1 , \ldots , s_n \f$, is the size of the smallest tuple of
 * features true in \f$ s \f$ and in none of \f$ s_1 , \ldots , s_n \f$
 * [see Lipovetzky and Geffner, ECAI 2012]. A NoveltyTable of width k keeps
 * the tuples of at most k features seen and gives
 * \f$ \min( w( s ) , k + 1 ) \f$. Since the tuples of the parent of
 * \f$ s \f$ are in the table, only those with a feature that the move
 * adds can be new, which an IncrementalFeatureDomain gives without
 * scanning the features of \f$ s \f$:
 *
 * - IW(k) is a breadth-first search that prunes the states of novelty
 *   larger than k, so that it generates \f$ O( n^k ) \f$ states for n
 *   features, and finds a path to any goal of width at most k [see
 *   Lipovetzky and Geffner, ECAI 2012, theorem 4];
 *
 * - SIW, on a GoalFeatureDomain, runs IW(1), IW(2), ... from the current
 *   state until one reaches a state with all the goal features of the
 *   current state and one more, which becomes the current state, i.e., it
 *   achieves the goal features one at a time [see the same paper];
 *
 * - BFWS is a best-first search ordered by \f$ ( w , \#g , g ) \f$, with
 *   \f$ \#g \f$ the goal features not achieved and the novelty \f$ w \f$
 *   computed separately for each value of \f$ ( \#g , \#r ) \f$,
 *   \f$ \#r \f$ being the features of a relaxed plan reached on the path,
 *   the plan computed again at each state where \f$ \#g \f$ decreases
 *   [see Lipovetzky and Geffner, AAAI 2017, BFWS( \f$ f_5 \f$ )]; on a
 *   Domain without goal features, \f$ \#g \f$ is replaced by the
 *   heuristic.
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

#ifndef __PBeFS_WidthSearch
 #define __PBeFS_WidthSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <iterator>

#include <limits>

#include <memory>

#include <stdexcept>

#include <type_traits>

#include <unordered_map>

#include <unordered_set>

#include <vector>

#include "PBeFS/astar/GreedySearch.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS NoveltyTable ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the tuples of at most k features seen so far
/** The tuples of one feature are bits of a vector of n bits, those of two a
 * triangle of \f$ n ( n - 1 ) / 2 \f$ bits, and those of \f$ m > 2 \f$
 * the codes \f$ \sum_j i_j n^j \f$ of their sorted features in a hash
 * set for each m, which needs \f$ n^k < 2^{64} \f$. */

class NoveltyTable {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the number n of features and the width k >= 1;
 /// throws std::invalid_argument if k = 0 or \f$ n^k \geq 2^{64} \f$

 NoveltyTable( std::size_t n , unsigned k ) : f_n( n ) , f_k( k ) ,
  v_one( ( n + 63 ) / 64 ) {
  if( k == 0 )
   throw( std::invalid_argument( "NoveltyTable: width 0" ) );
  if( k >= 2 )
   v_two.assign( ( n * ( n - ( n ? 1 : 0 ) ) / 2 + 63 ) / 64 , 0 );
  if( k >= 3 ) {
   v_more.resize( k + 1 );
   double p = 1;
   for( unsigned i = 0 ; i < k ; ++i )
    p *= double( n );
   if( p >= 18446744073709551616.0 )
    throw( std::invalid_argument( "NoveltyTable: n^k too large" ) );
   }
  }

/*-------------------------------- QUERIES ---------------------------------*/
 /// the width of the table

 [[nodiscard]] unsigned width( void ) const noexcept { return( f_k ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// records the tuples of the sorted features F, those with a feature of
 /// N only (a subset of F) if not all, returning the size of the
 /// smallest one that was new, or k + 1 if none was

 unsigned record( const std::vector< std::size_t > & F ,
                  const std::vector< std::size_t > & N ) {
  unsigned w = f_k + 1;
  for( std::size_t a : N )
   if( set_bit( v_one , a ) )
    w = 1;
  if( f_k >= 2 )
   for( std::size_t a : N )
    for( std::size_t b : F )
     if( ( a != b ) && set_bit( v_two , pair_index( a , b ) ) )
      w = std::min( w , 2u );
  for( unsigned m = 3 ; m <= f_k ; ++m )
   if( record_tuples( F , N , m ) )
    w = std::min( w , m );
  return( w );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// forgets all the tuples

 void clear( void ) {
  std::fill( v_one.begin() , v_one.end() , 0 );
  std::fill( v_two.begin() , v_two.end() , 0 );
  for( auto & x : v_more )
   x.clear();
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// sets the bit i of v, returning true if it was not set
 static bool set_bit( std::vector< std::uint64_t > & v , std::size_t i ) {
  const std::uint64_t b = std::uint64_t( 1 ) << ( i & 63 );
  std::uint64_t & x = v[ i >> 6 ];
  if( x & b )
   return( false );
  x |= b;
  return( true );
  }

 /// the index of the pair { a , b } in the triangle
 [[nodiscard]] std::size_t pair_index( std::size_t a , std::size_t b ) const
  noexcept {
  if( a > b )
   std::swap( a , b );
  return( b * ( b - 1 ) / 2 + a );
  }

 /// records the tuples of m features of F with one of N at least,
 /// returning true if one was new
 bool record_tuples( const std::vector< std::size_t > & F ,
                     const std::vector< std::size_t > & N , unsigned m ) {
  if( F.size() < m )
   return( false );
  std::vector< std::size_t > idx( m );
  for( unsigned i = 0 ; i < m ; ++i )
   idx[ i ] = i;
  const bool all = N.size() == F.size();
  bool fresh = false;
  for( ; ; ) {
   bool hit = all;
   std::uint64_t code = 0;
   for( unsigned i = 0 ; i < m ; ++i ) {
    const std::size_t f = F[ idx[ i ] ];
    if( ! hit )
     hit = std::binary_search( N.begin() , N.end() , f );
    code = code * f_n + f;
    }
   if( hit && v_more[ m ].insert( code ).second )
    fresh = true;
   // the next combination of m indices in [ 0 , |F| )
   int i = int( m ) - 1;
   while( ( i >= 0 ) && ( idx[ i ] == F.size() - m + std::size_t( i ) ) )
    --i;
   if( i < 0 )
    break;
   ++idx[ i ];
   for( unsigned j = unsigned( i ) + 1 ; j < m ; ++j )
    idx[ j ] = idx[ j - 1 ] + 1;
   }
  return( fresh );
  }

 std::size_t f_n;                                ///< the features
 unsigned f_k;                                   ///< the width
 std::vector< std::uint64_t > v_one;             ///< the single features
 std::vector< std::uint64_t > v_two;             ///< the pairs
 std::vector< std::unordered_set< std::uint64_t > > v_more;  ///< the
                                                 ///< tuples of m > 2
                                                 ///< features, by m

 };  // end( class( NoveltyTable ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the sorted features of s into F

template< FeatureDomain D >
void features_of( const D & d , const state_t< D > & s ,
                  std::vector< std::size_t > & F )
{
 F.clear();
 d.for_each_feature( s , [ & ]( std::size_t i ) { F.push_back( i ); } );
 std::sort( F.begin() , F.end() );
 F.erase( std::unique( F.begin() , F.end() ) , F.end() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the features of t, sorted in F, that its parent s does not have, into N:
/// from the changes of the move on an IncrementalFeatureDomain, from the
/// features Fs of s otherwise

template< FeatureDomain D >
void added_features( const D & d , const state_t< D > & s ,
                     const state_t< D > & t ,
                     const std::vector< std::size_t > & Fs ,
                     const std::vector< std::size_t > & F ,
                     std::vector< std::size_t > & N )
{
 N.clear();
 if constexpr( IncrementalFeatureDomain< D > ) {
  static_cast< void >( Fs );
  d.for_each_feature_change( s , t , [ & ]( std::size_t i ) {
   if( std::binary_search( F.begin() , F.end() , i ) )
    N.push_back( i );
   } );
  std::sort( N.begin() , N.end() );
  N.erase( std::unique( N.begin() , N.end() ) , N.end() );
  }
 else {
  static_cast< void >( d );
  static_cast< void >( s );
  static_cast< void >( t );
  std::set_difference( F.begin() , F.end() , Fs.begin() , Fs.end() ,
                       std::back_inserter( N ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS WidthOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of the width-based searches

//! [best-first-sig-width-options]
struct WidthOptions : public SearchOptions {

 /// the width k of IW, the largest width of SIW, and the width of the
 /// novelty of BFWS
 unsigned width = 2;

 /// BFWS: the states of novelty larger than width are pruned (k-BFWS)
 bool prune = false;

 /// BFWS: orders by the heuristic instead of the goal features to achieve,
 /// and computes the novelty for each value of it
 bool use_heuristic = false;

 /// BFWS: counts the features of the relaxed plans reached [see
 /// RelevantFeatureDomain]
 bool relevant = true;

 };  // end( struct( WidthOptions ) )
 //! [best-first-sig-width-options]

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS IteratedWidth ---------------------------*/
/*--------------------------------------------------------------------------*/
/// IW(k), the breadth-first search that prunes the states of novelty
/// larger than k
/** The goal test is done when a state is generated; a state already
 * generated is a duplicate and is not evaluated again. The Solution has
 * status unsolvable if the search ends without reaching a goal, whether or
 * not it pruned a state; pruned() tells how many states the last search
 * pruned by novelty, which SIW uses to tell a dead end. search() runs IW
 * from a given state to the first state that satisfies a predicate. */

template< FeatureDomain D >
class IteratedWidth {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit IteratedWidth( const D & domain , const WidthOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {}

/*-------------------------------- SOLVING ---------------------------------*/
 /// IW( options.width ) from the initial state to a goal

 [[nodiscard]] SolutionT solve( void ) {
  return( search( f_domain.initial_state() , f_options.width ,
                  [ & ]( const State & s , const std::vector< std::size_t > &
                         ) { return( f_domain.is_goal( s ) ); } ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// IW(k) from s0 to the first state s, of sorted features F, with
 /// target( s , F )

 template< class Target >
 [[nodiscard]] SolutionT search( const State & s0 , unsigned k ,
                                 Target && target ) {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  f_pruned = 0;
  NoveltyTable nt( f_domain.num_features() , k );
  NodePool< NodeT > pool;
  NodeTable< NodeT > table( f_options.expected_nodes
                            ? NodeTable< NodeT >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  std::deque< NodeT * > q;
  std::vector< std::size_t > Fs , F , N;

  const auto finish = [ & ]( const NodeT * goal ) {
   if( goal ) {
    sol.status = Status::solved;
    sol.cost = goal->g;
    sol.path = trace_path( goal );
    }
   st.stored = pool.size();
   st.expanded_by = { st.expanded };
   st.stored_by = { st.stored };
   st.seconds = clock.seconds();
   return( sol );
   };

  NodeT * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   NodeT * x = pool.allocate();
   *x = NodeT{ s0 , Cost( 0 ) , Cost( 0 ) , nullptr , false , false , no_op ,
               0 };
   return( x );
   } ).first;
  features_of( f_domain , s0 , F );
  nt.record( F , F );
  ++st.evaluated;
  if( target( s0 , F ) )
   return( finish( root ) );
  q.push_back( root );
  while( ! q.empty() ) {
   NodeT * n = q.front();
   q.pop_front();
   n->closed = true;
   ++st.expanded;
   n->expanded = true;
   features_of( f_domain , n->state , Fs );
   const NodeT * goal = nullptr;
   for_each_child( f_domain , n->state , Cost( 0 ) , n->op ,
                   [ & ]( const State & t , Cost c , Cost , std::uint8_t op ) {
    if( goal )
     return;
    ++st.generated;
    const std::uint64_t hash = f_domain.hash( t );
    if( table.find( t , hash ) ) {
     ++st.duplicates;
     return;
     }
    features_of( f_domain , t , F );
    added_features( f_domain , n->state , t , Fs , F , N );
    ++st.evaluated;
    if( nt.record( F , N ) > k ) {
     ++f_pruned;
     return;
     }
    NodeT * m = table.find_or_insert( t , hash , [ & ] {
     NodeT * x = pool.allocate();
     *x = NodeT{ t , n->g + c , Cost( 0 ) , n , false , false , op , 0 };
     return( x );
     } ).first;
    if( target( t , F ) )
     goal = m;
    else
     q.push_back( m );
    } );
   if( goal )
    return( finish( goal ) );
   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }
  return( finish( nullptr ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the states that the last search pruned by novelty

 [[nodiscard]] std::uint64_t pruned( void ) const noexcept {
  return( f_pruned );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using NodeT = GreedyNode< State , Cost >;         ///< the node type

 const D & f_domain;                                 ///< the Domain
 WidthOptions f_options;                             ///< the options
 std::uint64_t f_pruned = 0;                         ///< states pruned

 };  // end( class( IteratedWidth ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SerializedWidth --------------------------*/
/*--------------------------------------------------------------------------*/
/// SIW, IW run to achieve the goal features one at a time
/** From the current state s, with the goal features A of s, runs IW(1),
 * IW(2), ... , IW( options.width ) until one reaches a state with all the
 * features of A and one more goal feature, which becomes the current
 * state; it fails if IW( options.width ) does not, or if an IW that
 * pruned no state does not (then no such state is reachable). The
 * counters are summed over the calls of IW, and the path is the
 * concatenation of their paths. */

template< GoalFeatureDomain D >
class SerializedWidth {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit SerializedWidth( const D & domain ,
                           const WidthOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_domain.for_each_goal_feature( [ & ]( std::size_t i ) {
   v_goal.push_back( i );
   } );
  std::sort( v_goal.begin() , v_goal.end() );
  v_goal.erase( std::unique( v_goal.begin() , v_goal.end() ) ,
                v_goal.end() );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  IteratedWidth< D > iw( f_domain , f_options );
  State s = f_domain.initial_state();
  std::vector< State > plan{ s };
  Cost cost( 0 );
  std::vector< std::size_t > F , A;
  bool failed = false;
  f_calls = 0;
  f_max_width = 0;
  while( ! f_domain.is_goal( s ) ) {
   features_of( f_domain , s , F );
   A.clear();
   std::set_intersection( F.begin() , F.end() , v_goal.begin() ,
                          v_goal.end() , std::back_inserter( A ) );
   const std::size_t na = A.size();
   // a state with all of A, and one more goal feature, or a goal
   const auto target = [ & ]( const State & t ,
                              const std::vector< std::size_t > & Ft ) {
    if( f_domain.is_goal( t ) )
     return( true );
    std::size_t in = 0;
    for( std::size_t g : v_goal )
     in += std::binary_search( Ft.begin() , Ft.end() , g ) ? 1 : 0;
    if( in <= na )
     return( false );
    return( std::includes( Ft.begin() , Ft.end() , A.begin() , A.end() ) );
    };
   bool found = false;
   for( unsigned k = 1 ; k <= f_options.width ; ++k ) {
    SolutionT x = iw.search( s , k , target );
    ++f_calls;
    f_max_width = std::max( f_max_width , k );
    st.expanded += x.stats.expanded;
    st.generated += x.stats.generated;
    st.duplicates += x.stats.duplicates;
    st.evaluated += x.stats.evaluated;
    st.stored = std::max( st.stored , x.stats.stored );
    if( x.status == Status::node_limit ) {
     sol.status = Status::node_limit;
     failed = true;
     break;
     }
    if( x.solved() ) {
     plan.insert( plan.end() , x.path.begin() + 1 , x.path.end() );
     cost += x.cost;
     s = x.path.back();
     found = true;
     break;
     }
    if( iw.pruned() == 0 )
     break;                          // nothing more is reachable
    }
   if( failed || ! found ) {
    failed = true;
    break;
    }
   }
  if( ! failed ) {
   sol.status = Status::solved;
   sol.cost = cost;
   sol.path = std::move( plan );
   }
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the calls of IW, and the largest width used, of the last search

 [[nodiscard]] std::size_t calls( void ) const noexcept { return( f_calls ); }

 [[nodiscard]] unsigned max_width( void ) const noexcept {
  return( f_max_width );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const D & f_domain;                                 ///< the Domain
 WidthOptions f_options;                             ///< the options
 std::vector< std::size_t > v_goal;                  ///< the goal features
 std::size_t f_calls = 0;                            ///< the calls of IW
 unsigned f_max_width = 0;                           ///< the largest width

 };  // end( class( SerializedWidth ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BestFirstWidth --------------------------*/
/*--------------------------------------------------------------------------*/
/// BFWS, the best-first search ordered by novelty, then by the goal
/// features to achieve (or the heuristic), then by g
/** Each successor generated is evaluated, its novelty recorded in the
 * table of its value of \f$ ( \#g , \#r ) \f$ (or \f$ ( h , \#r ) \f$),
 * and then looked up: a state already generated is a duplicate, unless
 * its g decreases and it is not closed, in which case it is pushed again.
 * \f$ \#r \f$ is 0 unless the Domain is a RelevantFeatureDomain and
 * options.relevant: then the relevant features are those of a relaxed
 * plan of the nearest ancestor where \f$ \#g \f$ decreased (or of the
 * root), and \f$ \#r \f$ counts those added on the path from it. The goal
 * test is done at the expansion; with options.prune, the states of
 * novelty larger than options.width are pruned. */

template< FeatureDomain D >
class BestFirstWidth {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit BestFirstWidth( const D & domain ,
                          const WidthOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  if constexpr( GoalFeatureDomain< D > )
   f_domain.for_each_goal_feature( [ & ]( std::size_t i ) {
    v_goal.push_back( i );
    } );
  else
   f_options.use_heuristic = true;
  std::sort( v_goal.begin() , v_goal.end() );
  v_goal.erase( std::unique( v_goal.begin() , v_goal.end() ) ,
                v_goal.end() );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  NodePool< WNode > pool;
  NodeTable< WNode > table( f_options.expected_nodes
                            ? NodeTable< WNode >::slots_for(
                                               f_options.expected_nodes )
                            : 1024 );
  ADS::DaryHeap< Item , Better > open;
  std::unordered_map< std::uint64_t , std::unique_ptr< NoveltyTable > >
   tables;
  std::unordered_map< const WNode * , std::vector< std::size_t > > relevant;
  std::vector< std::size_t > Fs , F , N , R;
  std::uint64_t seq = 0;
  const bool use_r = relevant_used();

  // the table of the novelty of a value of ( #g or h , #r )
  const auto table_of = [ & ]( Cost key , std::uint32_t r ) -> NoveltyTable &
  {
   std::uint64_t k;
   if constexpr( std::is_floating_point_v< Cost > )
    k = std::bit_cast< std::uint64_t >( double( key ) );
   else
    k = std::uint64_t( std::int64_t( key ) );
   const std::uint64_t code = ADS::splitmix64( k * 0x9E3779B97F4A7C15ULL +
                                               r );
   auto & p = tables[ code ];
   if( ! p )
    p = std::make_unique< NoveltyTable >( f_domain.num_features() ,
                                          f_options.width );
   return( *p );
   };

  const State s0 = f_domain.initial_state();
  WNode * root = table.find_or_insert( s0 , f_domain.hash( s0 ) , [ & ] {
   WNode * x = pool.allocate();
   *x = WNode{ s0 , Cost( 0 ) , Cost( 0 ) , nullptr , false , false , no_op ,
               0 , 0 , 1 };
   return( x );
   } ).first;
  features_of( f_domain , s0 , F );
  root->key = order_key( s0 , F );
  ++st.evaluated;
  root->w = table_of( root->key , 0 ).record( F , F );
  if( use_r )
   relevant_of( root , relevant );
  open.push( Item{ root->w , root->key , root->g , seq++ , root } );

  const WNode * goal = nullptr;
  while( ! open.empty() ) {
   const Item e = open.pop();
   WNode * n = e.node;
   if( n->closed || ( n->g != e.g ) )
    continue;
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( f_domain.is_goal( n->state ) ) {
    goal = n;
    break;
    }
   if( use_r && ( n->parent ) && ( n->key < n->parent->key ) &&
       ! relevant.count( n ) )
    relevant_of( n , relevant );
   features_of( f_domain , n->state , Fs );
   for_each_child( f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
    ++st.generated;
    const Cost g = n->g + c;
    features_of( f_domain , t , F );
    const Cost key = order_key_h( F , ht );
    if( key == infinite_cost< Cost >() )
     return;                                   // a dead end
    const std::uint32_t r = use_r ? count_relevant( n , t , F , relevant )
                                  : 0;
    if( ( key == n->key ) && ( r == n->r ) )
     added_features( f_domain , n->state , t , Fs , F , N );
    else
     N = F;
    ++st.evaluated;
    const unsigned w = table_of( key , r ).record( F , N );
    if( f_options.prune && ( w > f_options.width ) )
     return;
    bool fresh = false;
    WNode * m = table.find_or_insert( t , f_domain.hash( t ) , [ & ] {
     fresh = true;
     WNode * x = pool.allocate();
     *x = WNode{ t , g , ht , n , false , false , op , 0 , r , w };
     x->key = key;
     return( x );
     } ).first;
    if( fresh ) {
     open.push( Item{ w , key , g , seq++ , m } );
     return;
     }
    if( cheaper( g , m->g ) && ! m->closed ) {
     m->g = g;
     m->parent = n;
     m->op = op;
     m->r = r;
     m->w = std::min( m->w , w );
     open.push( Item{ m->w , m->key , g , seq++ , m } );
     }
    else
     ++st.duplicates;
    } );
   if( f_options.node_limit && ( pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }

  if( goal ) {
   sol.status = Status::solved;
   sol.cost = goal->g;
   sol.path = trace_path( goal );
   }
  st.stored = pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a node: that of the greedy searches, with the key of its order (the
 /// goal features to achieve, or h), its count of relevant features and
 /// its novelty
 struct WNode {
  State state;
  Cost g;
  Cost h;
  WNode * parent;
  bool closed;
  bool expanded;
  std::uint8_t op;
  Cost key;
  std::uint32_t r;
  unsigned w;
  };

 /// an entry of OPEN
 struct Item {
  unsigned w;
  Cost key;
  Cost g;
  std::uint64_t seq;
  WNode * node;
  };

 /// the order ( w , key , g ) and then the first pushed
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   if( a.w != b.w )
    return( a.w < b.w );
   if( a.key != b.key )
    return( a.key < b.key );
   if( a.g != b.g )
    return( a.g < b.g );
   return( a.seq < b.seq );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if \#r is counted

 [[nodiscard]] bool relevant_used( void ) const noexcept {
  if constexpr( RelevantFeatureDomain< D > )
   return( f_options.relevant );
  else
   return( false );
  }

 /// the key of the order of a state of sorted features F and heuristic h
 [[nodiscard]] Cost order_key_h( const std::vector< std::size_t > & F ,
                                 Cost h ) const {
  if( f_options.use_heuristic )
   return( h );
  Cost missing( 0 );
  for( std::size_t g : v_goal )
   if( ! std::binary_search( F.begin() , F.end() , g ) )
    missing = missing + Cost( 1 );
  return( missing );
  }

 /// the key of the order of s, of sorted features F
 [[nodiscard]] Cost order_key( const State & s ,
                               const std::vector< std::size_t > & F ) const {
  return( order_key_h( F , f_options.use_heuristic
                           ? Cost( f_domain.heuristic( s ) ) : Cost( 0 ) ) );
  }

 /// computes the relevant features of the node n
 void relevant_of( const WNode * n ,
                   std::unordered_map< const WNode * ,
                                       std::vector< std::size_t > > & rel )
  const {
  if constexpr( RelevantFeatureDomain< D > ) {
   auto & v = rel[ n ];
   v.clear();
   f_domain.for_each_relevant_feature( n->state , [ & ]( std::size_t i ) {
    v.push_back( i );
    } );
   std::sort( v.begin() , v.end() );
   v.erase( std::unique( v.begin() , v.end() ) , v.end() );
   }
  else {
   static_cast< void >( n );
   static_cast< void >( rel );
   }
  }

 /// \#r of the successor t, of sorted features Ft, of n: the relevant
 /// features of the nearest ancestor a with them added on the path from a
 /// to t
 std::uint32_t count_relevant( const WNode * n , const State & t ,
                               const std::vector< std::size_t > & Ft ,
                               const std::unordered_map< const WNode * ,
                                std::vector< std::size_t > > & rel ) const {
  const WNode * a = n;
  while( ! rel.count( a ) )
   a = a->parent;
  const auto & R = rel.at( a );
  std::vector< std::size_t > got , Fp , Fc , N;
  // the edge from n to t, then those from the parents of n up to a
  added_features( f_domain , n->state , t , features( n->state , Fp ) , Ft ,
                  N );
  for( std::size_t i : N )
   if( std::binary_search( R.begin() , R.end() , i ) )
    got.push_back( i );
  for( const WNode * c = n ; c != a ; c = c->parent ) {
   features( c->state , Fc );
   features( c->parent->state , Fp );
   added_features( f_domain , c->parent->state , c->state , Fp , Fc , N );
   for( std::size_t i : N )
    if( std::binary_search( R.begin() , R.end() , i ) )
     got.push_back( i );
   }
  std::sort( got.begin() , got.end() );
  return( std::uint32_t( std::unique( got.begin() , got.end() ) -
                         got.begin() ) );
  }

 /// the sorted features of s into F, returned
 const std::vector< std::size_t > & features( const State & s ,
                                              std::vector< std::size_t > & F )
  const {
  features_of( f_domain , s , F );
  return( F );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 WidthOptions f_options;                             ///< the options
 std::vector< std::size_t > v_goal;                  ///< the goal features

 };  // end( class( BestFirstWidth ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with IW( options.width )

template< FeatureDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > iw(
 const D & domain , const WidthOptions & options = {} )
{
 return( IteratedWidth< D >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with SIW

template< GoalFeatureDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > siw(
 const D & domain , const WidthOptions & options = {} )
{
 return( SerializedWidth< D >( domain , options ).solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// searches a Domain with BFWS

template< FeatureDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > bfws(
 const D & domain , const WidthOptions & options = {} )
{
 return( BestFirstWidth< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_WidthSearch */

/*--------------------------------------------------------------------------*/
/*------------------------- End File WidthSearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
