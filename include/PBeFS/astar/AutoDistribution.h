/*--------------------------------------------------------------------------*/
/*------------------------ File AutoDistribution.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the work distributions of hash-based parallel A* that
 * are built automatically from the structure of a FeatureDomain, and for
 * the estimates of their overheads by which one of them is chosen [see
 * Jinnai and Fukunaga, AAAI 2016 and JAIR 2017].
 *
 * The structure is that of a FeatureModel: a sample of the state space
 * near the root, expanded in the order of A*, which stands for the
 * workload graph of section 5 of the JAIR paper (the states of
 * \f$ f < f^* \f$ and their edges), and, from the moves of the sample, the
 * groups of the features (the atom groups, i.e., the variables of SAS+),
 * the transition graph of each group (its domain transition graph, DTG),
 * the operators, i.e., the sets of features that a move changes, and how
 * often each of them is applied. In the notation HDA*[ hash , abstraction ]
 * of the paper, the distributions built on it are
 *
 * - greedy_abstraction(), HDA*[Z,Afeature/DTGgreedy] (GAZHDA*): each DTG
 *   cut in two halves grown greedily from its node of least degree;
 *
 * - sparsest_cut_abstraction(), HDA*[Z,Afeature/DTGsparsity]
 *   (GRAZHDA* with sparsity): each DTG cut in the two sides of largest
 *   sparsity, or of largest estimated efficiency;
 *
 * - fluency_abstraction(), HDA*[Z,Afeature/DTGfluency] (FAZHDA*): the
 *   groups that the moves change most often ignored, the others cut as in
 *   GAZHDA*;
 *
 * - dynamic_state_abstraction(), HDA*[Z,Astate/SDDdynamic] (DAHDA*): the
 *   AHDA* of Burns et al. with the groups chosen greedily until they hold a
 *   given fraction of the features;
 *
 * - kernel_selection_table() (KZHDA*): the Zobrist table whose key, a
 *   linear map over GF(2) of the features, vanishes on the change vectors
 *   of the most frequent operators, as many as the sample says, of which
 *   Z, AZ, OZ and the cuts of the DTGs are particular cases;
 *
 * - operator_table(), HDA*[Zoperator] (OZHDA*): Zobrist hashing whose
 *   table is chosen so that the most frequent operators change no bit of
 *   the key, which decides the owner;
 *
 * and auto_distribution() chooses among them (and the ones given by the
 * caller) the one of largest estimated efficiency on the sample. All are
 * TableZobristDistribution [see Distribution.h], so that the key of a
 * successor is computed from that of its parent.
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

#ifndef __PBeFS_AutoDistribution
 #define __PBeFS_AutoDistribution
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <functional>

#include <iterator>

#include <limits>

#include <map>

#include <queue>

#include <stdexcept>

#include <string>

#include <tuple>

#include <unordered_map>

#include <unordered_set>

#include <utility>

#include <vector>

#include "PBeFS/astar/Distribution.h"

#include "PBeFS/astar/Domain.h"

#include "ADS/dictionary/RandomHashing.h"

#include "ADS/mfset/UnionFind.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS DistributionEstimate -----------------------*/
/*--------------------------------------------------------------------------*/
/// the overheads of a distribution estimated on a workload graph
/** With \f$ S_i \f$ the expanded states of the sample that the worker
 * \f$ i \f$ owns, \f$ S \f$ all of them and \f$ c \f$ the cost of sending
 * a node relative to that of generating it [see Jinnai and Fukunaga, JAIR
 * 2017, equations (3) to (6)]:
 *
 * - co, the fraction of the edges ( s , t ), s expanded and t generated,
 *   whose ends have different owners, i.e., the CO measured by HDA* (sent
 *   over generated) on the expansions of the sample; equation (3) counts
 *   instead the undirected edges between states of \f$ f < f^* \f$ only;
 *
 * - lb, \f$ \max_i | S_i | / ( | S | / p ) \f$;
 *
 * - so, \f$ \sum_i ( \max_j | S_j | - | S_i | ) / | S | = \mathrm{lb} - 1
 *   \f$, the expansions that the other workers make while the most loaded
 *   one expands its states, relative to \f$ | S | \f$ (equation (4) of the
 *   paper writes the same sum as \f$ p ( \mathrm{LB} - 1 ) \f$, which is
 *   its value relative to the mean \f$ | S | / p \f$ rather than to
 *   \f$ | S | \f$, and which equation (5) then adds to 1 as a ratio to
 *   the expansions of A*);
 *
 * - efficiency, \f$ 1 / ( ( 1 + c \, \mathrm{co} ) ( 1 + \mathrm{so} ) )
 *   \f$;
 *
 * - sparsity, \f$ \prod_i ( p | S_i | / | S | ) / \mathrm{co} \f$, the
 *   sparsity of the partition with the sizes scaled by \f$ p \f$, so that a
 *   balanced partition has sparsity \f$ 1 / \mathrm{co} \f$ (infinite if no
 *   edge is cut). */

struct DistributionEstimate {

 double co = 0;              ///< the communication overhead
 double lb = 1;              ///< the load balance
 double so = 0;              ///< the search overhead, lb - 1
 double efficiency = 1;      ///< the estimated efficiency
 double sparsity = 0;        ///< the scaled sparsity

 };  // end( struct( DistributionEstimate ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS FeatureGraph ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the transition structure of the features of a sample, without the
/// Domain: the part of a FeatureModel on which the distributions are
/// built, so that their construction is compiled once for all the Domains

class FeatureGraph {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 /// the group of a feature that no move changes
 static constexpr std::uint32_t no_group = 0xFFFFFFFF;

 /// an edge of a DTG, to the feature to, of weight w
 struct Arc {
  std::uint32_t to;
  double w;
  };

 /// an operator: the features that its moves change, sorted, and the
 /// number of its moves in the sample
 struct Operator {
  std::vector< std::uint32_t > features;
  double weight;
  };

/*---------------------------- THE SAMPLE GRAPH ----------------------------*/
 /// the number of expanded states, the first ones of states()

 [[nodiscard]] std::size_t num_expanded( void ) const noexcept {
  return( f_expanded );
  }

 /// the edges ( s , t ) of the sample, as indices in states(), s expanded

 [[nodiscard]] const std::vector< std::pair< std::uint32_t ,
                                             std::uint32_t > > &
  edges( void ) const noexcept { return( v_edge ); }

/*------------------------------ THE FEATURES ------------------------------*/
 /// the number of features

 [[nodiscard]] std::size_t num_features( void ) const noexcept {
  return( f_n );
  }

 /// the number of groups

 [[nodiscard]] std::size_t num_groups( void ) const noexcept {
  return( v_members.size() );
  }

 /// the group of the feature i, no_group if no move changes it

 [[nodiscard]] std::uint32_t group( std::size_t i ) const noexcept {
  return( v_group[ i ] );
  }

 /// true if the feature i holds in a state of the sample

 [[nodiscard]] bool seen( std::size_t i ) const noexcept {
  return( v_seen[ i ] );
  }

 /// the number of the expanded states of the sample in which the feature
 /// i holds

 [[nodiscard]] double mass( std::size_t i ) const noexcept {
  return( v_mass[ i ] );
  }

 /// the nodes of the DTG of the group g, increasing: its features that
 /// hold in a state of the sample or that a ground transition of a
 /// FeatureTransitionDomain changes

 [[nodiscard]] const std::vector< std::uint32_t > & members( std::size_t g )
  const noexcept { return( v_members[ g ] ); }

 /// the edges of the DTG of the group of the feature i at i, by
 /// increasing neighbour, each weighing the moves of the sample that cross
 /// it in either direction

 [[nodiscard]] const std::vector< Arc > & arcs( std::size_t i ) const
  noexcept { return( v_adj[ i ] ); }

 /// the edges of the DTG at i, by increasing neighbour, each weighing the
 /// ground operators that realize it in either direction: those of the
 /// Domain if it is a FeatureTransitionDomain, the distinct operators of
 /// the sample otherwise

 [[nodiscard]] const std::vector< Arc > & static_arcs( std::size_t i ) const
  noexcept { return( v_sadj[ i ] ); }

 /// the fluency of the group g: the fraction of the moves of the sample
 /// that change a feature of g

 [[nodiscard]] double fluency( std::size_t g ) const noexcept {
  return( f_moves > 0 ? v_changes[ g ] / f_moves : 0 );
  }

 /// the operators, by decreasing weight (and increasing features on ties)

 [[nodiscard]] const std::vector< Operator > & operators( void ) const
  noexcept { return( v_op ); }

 /// the words of the states of the sample in the groups: words()[ k G + g ]
 /// is the exclusive or of random words of the features of the group g
 /// that hold in the k-th state of the sample, G being num_groups()

 [[nodiscard]] const std::vector< std::uint64_t > & words( void ) const
  noexcept { return( v_word ); }

 /// the features of the k-th state of the sample, increasing

 [[nodiscard]] std::pair< const std::uint32_t * , const std::uint32_t * >
  state_features( std::size_t k ) const noexcept {
  return( std::pair( v_sfeat.data() + v_soff[ k ] ,
                     v_sfeat.data() + v_soff[ k + 1 ] ) );
  }

/*------------------------------- ESTIMATES --------------------------------*/
 /// the overheads among p workers on the sample of Zobrist hashing on the
 /// table r, one entry per feature, with c the cost of sending a node
 /// relative to that of generating it [see DistributionEstimate]
 /** It is FeatureModel::estimate() of TableZobristDistribution( d , r ),
  * computed from the features of the states of the sample without the
  * Domain. */

 [[nodiscard]] DistributionEstimate estimate_table(
  const std::vector< std::uint64_t > & r , unsigned p , double c = 1 ) const
 {
  const std::size_t ns = v_soff.empty() ? 0 : v_soff.size() - 1;
  std::vector< unsigned > own( ns );
  for( std::size_t k = 0 ; k < ns ; ++k ) {
   std::uint64_t z = 0;
   for( auto i = v_soff[ k ] ; i < v_soff[ k + 1 ] ; ++i )
    z ^= r[ v_sfeat[ i ] ];
   own[ k ] = key_owner( std::uint32_t( z >> 32 ) , p );
   }
  return( estimate_owners( own , p , c ) );
  }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

 /// constructor, for n features

 explicit FeatureGraph( std::size_t n )
  : f_n( n ) , v_group( n , no_group ) , v_seen( n , false ) ,
    v_mass( n , 0 ) , v_adj( n ) , v_sadj( n ) {}

 /// the estimate of the owners own[ k ] of the states of the sample among
 /// p workers

 [[nodiscard]] DistributionEstimate estimate_owners(
  const std::vector< unsigned > & own , unsigned p , double c ) const {
  DistributionEstimate e;
  if( ( p == 0 ) || ( f_expanded == 0 ) )
   return( e );
  std::vector< double > load( p , 0 );
  for( std::size_t k = 0 ; k < f_expanded ; ++k )
   load[ own[ k ] ] += 1;
  std::size_t cut = 0;
  for( const auto & [ s , t ] : v_edge )
   if( own[ s ] != own[ t ] )
    ++cut;
  const double mean = double( f_expanded ) / p;
  e.co = v_edge.empty() ? 0 : double( cut ) / double( v_edge.size() );
  e.lb = *std::max_element( load.begin() , load.end() ) / mean;
  e.so = e.lb - 1;
  e.efficiency = 1 / ( ( 1 + c * e.co ) * ( 1 + e.so ) );
  double prod = 1;
  for( double l : load )
   prod *= l / mean;
  e.sparsity = e.co > 0 ? prod / e.co
                        : std::numeric_limits< double >::infinity();
  return( e );
  }

 std::size_t f_n;                                      ///< the features
 std::size_t f_expanded = 0;             ///< the expanded ones
 std::vector< std::pair< std::uint32_t , std::uint32_t > > v_edge;
                                         ///< the edges of the sample
 std::vector< std::uint32_t > v_group;   ///< the group of each feature
 std::vector< bool > v_seen;             ///< the features of the sample
 std::vector< double > v_mass;           ///< the expanded states of each
                                         ///< feature
 std::vector< std::vector< std::uint32_t > > v_members;
                                         ///< the nodes of each DTG
 std::vector< std::vector< Arc > > v_adj;  ///< the DTG edges at a feature
 std::vector< std::vector< Arc > > v_sadj; ///< the same, static weights
 std::vector< double > v_changes;        ///< the moves changing a group
 double f_moves = 0;                     ///< the moves of the sample
 std::vector< Operator > v_op;           ///< the operators
 std::vector< std::uint64_t > v_word;   ///< the word of each state in
                                         ///< each group
 std::vector< std::uint32_t > v_sfeat;  ///< the features of the states
 std::vector< std::size_t > v_soff;     ///< where those of each start

 };  // end( class( FeatureGraph ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS FeatureModel ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a sample of the state space and the transition structure of the features
/** The sample is built by A* from the initial state, with the Domain's
 * heuristic and the ties broken by larger g, up to a given number of
 * expansions or to the first goal: the expanded states are those of
 * \f$ f < f^* \f$ if the search reaches the goal, and the nearest to the
 * root otherwise. For each move ( s , t ) of an expanded s, the features
 * of s that t lacks are deleted and those of t that s lacks are added; a
 * deleted and an added feature of the same group make an edge of the DTG
 * of the group, of weight the number of the moves that make it, and the
 * set of the features changed is an operator, of weight the number of its
 * moves. The groups are those of the Domain if it is a FeatureGroupDomain,
 * and otherwise the connected components of the graph that joins every
 * deleted feature of a move to every added one; a feature that no move
 * changes has no group. */

template< FeatureDomain D >
class FeatureModel : public FeatureGraph {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, sampling the Domain with at most \p expansions expansions

 explicit FeatureModel( const D & domain , std::size_t expansions = 32768 )
  : FeatureGraph( domain.num_features() ) , f_domain( domain ) {
  sample( expansions );
  }

/*---------------------------- THE SAMPLE GRAPH ----------------------------*/
 /// the Domain

 [[nodiscard]] const D & domain( void ) const noexcept { return( f_domain ); }

 /// the states of the sample, the expanded ones first

 [[nodiscard]] const std::vector< State > & states( void ) const noexcept {
  return( v_state );
  }

/*------------------------------- ESTIMATES --------------------------------*/
 /// the overheads of the distribution x among p workers on the sample,
 /// with c the cost of sending a node relative to that of generating it
 /** \p x is any distribution with a method owner( s , p ). */

 template< class X >
 [[nodiscard]] DistributionEstimate estimate( const X & x , unsigned p ,
                                              double c = 1 ) const {
  std::vector< unsigned > own( v_state.size() );
  for( std::size_t k = 0 ; k < v_state.size() ; ++k )
   own[ k ] = p ? unsigned( x.owner( v_state[ k ] , p ) ) : 0;
  return( estimate_owners( own , p , c ) );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the features of s, increasing

 std::vector< std::uint32_t > features( const State & s ) const {
  std::vector< std::uint32_t > v;
  f_domain.for_each_feature( s , [ & ]( std::size_t i ) {
   v.push_back( std::uint32_t( i ) );
   } );
  std::sort( v.begin() , v.end() );
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // A* up to the given expansions or to a goal, recording the sample and
 // the moves

 void sample( std::size_t expansions ) {
  struct Hash {
   const D * d;
   std::size_t operator()( const State & s ) const {
    return( std::size_t( d->hash( s ) ) );
    }
   };
  std::unordered_map< State , std::uint32_t , Hash > index(
                                                     16 , Hash{ &f_domain } );
  std::vector< Cost > g;
  std::vector< Cost > h;
  std::vector< bool > closed;
  // ( f , -g , index ) of the entries of OPEN, the smallest first
  using Entry = std::tuple< Cost , Cost , std::uint32_t >;
  std::priority_queue< Entry , std::vector< Entry > ,
                       std::greater< Entry > > open;
  const auto add = [ & ]( const State & s , Cost gs , Cost hs ) {
   const auto [ it , fresh ] = index.try_emplace( s ,
                                                  std::uint32_t( g.size() ) );
   const std::uint32_t k = it->second;
   if( fresh ) {
    g.push_back( gs );
    h.push_back( hs );
    closed.push_back( false );
    v_order.push_back( s );
    }
   else
    if( gs < g[ k ] ) {
     g[ k ] = gs;
     closed[ k ] = false;
     }
    else
     return( k );
   open.emplace( gs + hs , -gs , k );
   return( k );
   };

  const State s0 = f_domain.initial_state();
  add( s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) );
  std::vector< std::uint32_t > expanded;                 // in order
  std::vector< std::pair< std::uint32_t , std::uint32_t > > edge;
  std::map< std::vector< std::uint32_t > , double > ops;
  std::map< std::pair< std::uint32_t , std::uint32_t > , double > dtg ,
                                                                  sdtg;
  ADS::DisjointSets<> uf{ int( f_n ) };
  std::vector< std::uint32_t > del , ins , chg;

  while( ! open.empty() && ( expanded.size() < expansions ) ) {
   const auto [ f , mg , k ] = open.top();
   open.pop();
   if( closed[ k ] || ( -mg != g[ k ] ) )
    continue;
   closed[ k ] = true;
   expanded.push_back( k );
   const State s = v_order[ k ];
   if( f_domain.is_goal( s ) )
    break;
   const auto fs = features( s );
   for( std::uint32_t i : fs ) {
    v_seen[ i ] = true;
    v_mass[ i ] += 1;
    }
   const Cost gs = g[ k ];
   std::vector< std::tuple< State , Cost , Cost > > succ;
   f_domain.for_each_successor( s , h[ k ] , [ & ]( const State & t , Cost c ,
                                                     Cost ht ) {
    succ.emplace_back( t , c , ht );
    } );
   for( const auto & [ t , c , ht ] : succ ) {
    const std::uint32_t kt = add( t , gs + c , ht );
    edge.emplace_back( k , kt );
    const auto ft = features( t );
    for( std::uint32_t i : ft )
     v_seen[ i ] = true;
    del.clear();
    ins.clear();
    std::set_difference( fs.begin() , fs.end() , ft.begin() , ft.end() ,
                         std::back_inserter( del ) );
    std::set_difference( ft.begin() , ft.end() , fs.begin() , fs.end() ,
                         std::back_inserter( ins ) );
    if( del.empty() && ins.empty() )
     continue;
    ++f_moves;
    chg.clear();
    std::merge( del.begin() , del.end() , ins.begin() , ins.end() ,
                std::back_inserter( chg ) );
    const auto [ it , fresh ] = ops.try_emplace( chg , 0.0 );
    it->second += 1;
    for( std::uint32_t a : del )
     for( std::uint32_t b : ins ) {
      if constexpr( FeatureGroupDomain< D > ) {
       if( f_domain.feature_group( a ) != f_domain.feature_group( b ) )
        continue;
       }
      else
       uf.merge( int( a ) , int( b ) );
      dtg[ { std::min( a , b ) , std::max( a , b ) } ] += 1;
      if( fresh && ! FeatureTransitionDomain< D > )
       sdtg[ { std::min( a , b ) , std::max( a , b ) } ] += 1;
      }
    }
   }

  // the states, the expanded ones first, and the edges between them
  std::vector< std::uint32_t > pos( v_order.size() ,
                                    std::uint32_t( v_order.size() ) );
  for( std::uint32_t k : expanded ) {
   pos[ k ] = std::uint32_t( v_state.size() );
   v_state.push_back( v_order[ k ] );
   }
  f_expanded = v_state.size();
  for( std::uint32_t k = 0 ; k < v_order.size() ; ++k )
   if( pos[ k ] == v_order.size() ) {
    pos[ k ] = std::uint32_t( v_state.size() );
    v_state.push_back( v_order[ k ] );
    }
  v_order.clear();
  v_order.shrink_to_fit();
  v_edge.reserve( edge.size() );
  for( const auto & [ s , t ] : edge )
   v_edge.emplace_back( pos[ s ] , pos[ t ] );

  // the ground transitions of the Domain
  if constexpr( FeatureTransitionDomain< D > )
   f_domain.for_each_feature_transition( [ & ]( std::size_t a ,
                                                std::size_t b ) {
    if constexpr( ! FeatureGroupDomain< D > )
     uf.merge( int( a ) , int( b ) );
    const auto x = std::uint32_t( std::min( a , b ) );
    const auto y = std::uint32_t( std::max( a , b ) );
    sdtg[ { x , y } ] += 1;
    } );

  // the groups of the features that the moves or the transitions change
  std::vector< bool > moved( f_n , false );
  for( const auto & [ fe , w ] : ops )
   for( std::uint32_t i : fe )
    moved[ i ] = true;
  for( const auto & [ e , w ] : sdtg )
   moved[ e.first ] = moved[ e.second ] = true;
  std::vector< std::uint32_t > id( f_n , no_group );
  std::uint32_t ng = 0;
  for( std::size_t i = 0 ; i < f_n ; ++i ) {
   if( ! moved[ i ] )
    continue;
   std::size_t r;
   if constexpr( FeatureGroupDomain< D > )
    r = f_domain.feature_group( i );
   else
    r = std::size_t( uf.find( int( i ) ) );
   if( r >= id.size() )
    id.resize( r + 1 , no_group );
   if( id[ r ] == no_group )
    id[ r ] = ng++;
   v_group[ i ] = id[ r ];
   }
  std::vector< bool > node( f_n , false );
  for( const auto & [ e , w ] : sdtg )
   node[ e.first ] = node[ e.second ] = true;
  v_members.assign( ng , {} );
  for( std::size_t i = 0 ; i < f_n ; ++i )
   if( ( v_group[ i ] != no_group ) && ( v_seen[ i ] || node[ i ] ) )
    v_members[ v_group[ i ] ].push_back( std::uint32_t( i ) );
  const auto arcs_of = []( const auto & m , auto & adj ) {
   for( const auto & [ e , w ] : m ) {
    adj[ e.first ].push_back( Arc{ e.second , w } );
    adj[ e.second ].push_back( Arc{ e.first , w } );
    }
   for( auto & a : adj )
    std::sort( a.begin() , a.end() , []( const Arc & x , const Arc & y ) {
     return( x.to < y.to );
     } );
   };
  arcs_of( dtg , v_adj );
  arcs_of( sdtg , v_sadj );

  // the operators by decreasing weight, and the fluency of the groups
  v_changes.assign( ng , 0 );
  std::vector< std::uint32_t > gs;
  for( auto & [ fe , w ] : ops ) {
   gs.clear();
   for( std::uint32_t i : fe )
    gs.push_back( v_group[ i ] );
   std::sort( gs.begin() , gs.end() );
   gs.erase( std::unique( gs.begin() , gs.end() ) , gs.end() );
   for( std::uint32_t x : gs )
    v_changes[ x ] += w;
   v_op.push_back( Operator{ fe , w } );
   }
  std::stable_sort( v_op.begin() , v_op.end() ,
                    []( const Operator & a , const Operator & b ) {
                     return( a.weight > b.weight );
                     } );

  // the features of the states, and their words in the groups
  const ADS::ZobristTable r( f_n ? f_n : 1 , 1 , 0x5eed );
  v_word.assign( v_state.size() * ng , 0 );
  v_soff.assign( 1 , 0 );
  for( std::size_t k = 0 ; k < v_state.size() ; ++k ) {
   f_domain.for_each_feature( v_state[ k ] , [ & ]( std::size_t i ) {
    v_sfeat.push_back( std::uint32_t( i ) );
    if( const auto x = v_group[ i ] ; x != no_group )
     v_word[ k * ng + x ] ^= r( i , 0 );
    } );
   v_soff.push_back( v_sfeat.size() );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                   ///< the Domain
 std::vector< State > v_order;           ///< the states, as discovered
 std::vector< State > v_state;           ///< the states of the sample
 };  // end( class( FeatureModel ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the objective by which sparsest_cut_abstraction() cuts a DTG

enum class CutObjective {
 sparsity ,          ///< the largest \f$ |S_1| |S_2| / w( S_1 , S_2 ) \f$,
                     ///< with the static weights
 sampled_sparsity ,  ///< the same, with the weights of the sample
 efficiency          ///< the largest \f$ 1 / ( ( 1 + c \, CO ) \, LB ) \f$,
                     ///< with the weights of the sample
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// the abstraction that gives the side side[ k ] of the node members[ k ]
 /// of each cut DTG, ignores the groups in ignored and keeps every other
 /// feature by itself: sides holds, by group, the sides of its members
 /// (empty for a group kept feature by feature)

 inline FeatureAbstraction cut_abstraction(
  const FeatureGraph & m ,
  const std::vector< std::vector< unsigned char > > & sides ,
  const std::vector< bool > & ignored ) {
  std::vector< std::uint32_t > a( m.num_features() , no_abstract_feature );
  std::uint32_t next = 0;
  std::vector< std::uint32_t > side_id;
  for( std::size_t g = 0 ; g < m.num_groups() ; ++g ) {
   if( ignored[ g ] || sides[ g ].empty() )
    continue;
   side_id.assign( 2 , no_abstract_feature );
   const auto & v = m.members( g );
   for( std::size_t k = 0 ; k < v.size() ; ++k ) {
    std::uint32_t & x = side_id[ sides[ g ][ k ] ? 1 : 0 ];
    if( x == no_abstract_feature )
     x = next++;
    a[ v[ k ] ] = x;
    }
   }
  for( std::size_t i = 0 ; i < m.num_features() ; ++i ) {
   const std::uint32_t g = m.group( i );
   if( ( g != FeatureGraph::no_group ) && ignored[ g ] )
    continue;
   if( a[ i ] == no_abstract_feature )
    a[ i ] = next++;
   }
  return( FeatureAbstraction::of( std::move( a ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the halves of Algorithm 4 of Jinnai and Fukunaga (JAIR 2017) on the
 /// DTG of the group g: side 0 grown from the node of least degree by the
 /// node with the most edges to it, until it holds half of the nodes

 inline std::vector< unsigned char > greedy_cut( const FeatureGraph & m ,
                                          std::size_t g ) {
  const auto & v = m.members( g );
  const std::size_t n = v.size();
  std::vector< unsigned char > side( n , 1 );
  if( n < 2 )
   return( side );
  std::unordered_map< std::uint32_t , std::size_t > at;
  for( std::size_t k = 0 ; k < n ; ++k )
   at[ v[ k ] ] = k;
  std::vector< std::size_t > deg( n , 0 ) , link( n , 0 );
  for( std::size_t k = 0 ; k < n ; ++k )
   for( const auto & e : m.static_arcs( v[ k ] ) )
    if( at.count( e.to ) )
     ++deg[ k ];
  const auto take = [ & ]( std::size_t k ) {
   side[ k ] = 0;
   for( const auto & e : m.static_arcs( v[ k ] ) )
    if( const auto it = at.find( e.to ) ; it != at.end() )
     ++link[ it->second ];
   };
  take( std::size_t( std::min_element( deg.begin() , deg.end() ) -
                     deg.begin() ) );
  for( std::size_t in = 1 ; 2 * in < n ; ++in ) {
   std::size_t best = n;
   for( std::size_t k = 0 ; k < n ; ++k )
    if( side[ k ] && ( ( best == n ) || ( link[ k ] > link[ best ] ) ) )
     best = k;
   take( best );
   }
  return( side );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the two sides of the DTG of the group g of best objective, with the
 /// static weights (a node weighing 1) for CutObjective::sparsity, and
 /// with those of the sample (a node weighing as the expanded states of its
 /// feature) otherwise: by enumeration up to 20 nodes, and by local search
 /// from several starts above

 inline std::vector< unsigned char > best_cut( const FeatureGraph & m ,
                                        std::size_t g , CutObjective o ,
                                        double c ) {
  const auto & v = m.members( g );
  const std::size_t n = v.size();
  if( n < 2 )
   return( std::vector< unsigned char >( n , 0 ) );
  std::unordered_map< std::uint32_t , std::size_t > at;
  for( std::size_t k = 0 ; k < n ; ++k )
   at[ v[ k ] ] = k;
  std::vector< std::vector< std::pair< std::size_t , double > > > adj( n );
  std::vector< double > q( n );
  double total = 0 , mass = 0;
  const bool unit = o == CutObjective::sparsity;
  for( std::size_t k = 0 ; k < n ; ++k ) {
   q[ k ] = unit ? 1 : m.mass( v[ k ] );
   mass += q[ k ];
   for( const auto & e : unit ? m.static_arcs( v[ k ] ) : m.arcs( v[ k ] ) )
    if( const auto it = at.find( e.to ) ; it != at.end() ) {
     adj[ k ].emplace_back( it->second , e.w );
     total += e.w / 2;
     }
   }
  if( mass <= 0 ) {
   std::fill( q.begin() , q.end() , 1 );
   mass = double( n );
   }

  // the value of a cut of weight w with mass m1 on side 1 (the larger the
  // better), a cut of weight 0 counting as one of an infinitesimal weight
  const double eps = total > 0 ? total * 1e-9 : 1;
  const auto value = [ & ]( double m1 , double w ) {
   const double a = std::max( mass - m1 , 0.0 ) , b = std::max( m1 , 0.0 );
   if( o != CutObjective::efficiency )
    return( a * b / std::max( w , eps ) );
   const double co = total > 0 ? w / total : 0;
   const double lb = std::max( a , b ) / ( mass / 2 );
   return( 1 / ( ( 1 + c * co ) * lb ) );
   };
  // the change of the cut weight if the node k changes side
  const auto delta = [ & ]( const std::vector< unsigned char > & s ,
                            std::size_t k ) {
   double d = 0;
   for( const auto & [ u , w ] : adj[ k ] )
    d += s[ u ] == s[ k ] ? w : -w;
   return( d );
   };

  std::vector< unsigned char > best( n , 0 );
  double best_v = value( 0 , 0 );
  if( n <= 20 ) {
   // the Gray code on the nodes 1 , ... , n - 1, node 0 on side 0
   std::vector< unsigned char > s( n , 0 );
   double m1 = 0 , w = 0;
   for( std::uint64_t i = 1 ; i < ( std::uint64_t( 1 ) << ( n - 1 ) ) ; ++i )
   {
    const std::size_t k = std::size_t( std::countr_zero( i ) ) + 1;
    w += delta( s , k );
    s[ k ] ^= 1;
    m1 += s[ k ] ? q[ k ] : -q[ k ];
    if( const double x = value( m1 , w ) ; x > best_v ) {
     best_v = x;
     best = s;
     }
    }
   return( best );
   }

  // passes of Fiduccia and Mattheyses from the greedy halves and from
  // the halves grown in breadth-first order from up to 16 nodes: in a
  // pass every node moves once, the one of best value first even if the
  // value drops, and the best of the cuts met is kept; the passes are
  // repeated while they improve
  std::vector< std::vector< unsigned char > > start;
  start.push_back( greedy_cut( m , g ) );
  const std::size_t seeds = std::min< std::size_t >( n , 16 );
  for( std::size_t j = 0 ; j < seeds ; ++j ) {
   const std::size_t r = j * n / seeds;
   std::vector< unsigned char > s( n , 1 );
   std::vector< std::size_t > fifo{ r };
   s[ r ] = 0;
   for( std::size_t h = 0 , in = 1 ; ( h < fifo.size() ) && ( 2 * in < n ) ;
        ++h )
    for( const auto & [ u , w ] : adj[ fifo[ h ] ] )
     if( s[ u ] && ( 2 * in < n ) ) {
      s[ u ] = 0;
      ++in;
      fifo.push_back( u );
      }
   start.push_back( std::move( s ) );
   }
  std::vector< unsigned char > locked( n );
  std::vector< std::size_t > moved;
  for( auto & s : start ) {
   double m1 = 0 , w = 0;
   for( std::size_t k = 0 ; k < n ; ++k ) {
    if( s[ k ] )
     m1 += q[ k ];
    for( const auto & [ u , x ] : adj[ k ] )
     if( ( u > k ) && ( s[ u ] != s[ k ] ) )
      w += x;
    }
   double cur = value( m1 , w );
   for( bool better = true ; better ; ) {
    better = false;
    std::fill( locked.begin() , locked.end() , 0 );
    moved.clear();
    double pm1 = m1 , pw = w , bx = cur;
    std::size_t bl = 0;                       // the moves of the best cut
    for( std::size_t step = 0 ; step < n ; ++step ) {
     std::size_t bk = n;
     double vx = 0;
     for( std::size_t k = 0 ; k < n ; ++k ) {
      if( locked[ k ] )
       continue;
      const double y = s[ k ] ? pm1 - q[ k ] : pm1 + q[ k ];
      const double x = value( y , pw + delta( s , k ) );
      if( ( bk == n ) || ( x > vx ) ) {
       bk = k;
       vx = x;
       }
      }
     pw += delta( s , bk );
     pm1 += s[ bk ] ? -q[ bk ] : q[ bk ];
     s[ bk ] ^= 1;
     locked[ bk ] = 1;
     moved.push_back( bk );
     if( vx > bx * ( 1 + 1e-12 ) ) {
      bx = vx;
      bl = moved.size();
      }
     }
    // back to the best cut of the pass
    for( std::size_t j = moved.size() ; j > bl ; --j )
     s[ moved[ j - 1 ] ] ^= 1;
    if( bl > 0 ) {
     better = true;
     cur = bx;
     m1 = 0;
     w = 0;
     for( std::size_t k = 0 ; k < n ; ++k ) {
      if( s[ k ] )
       m1 += q[ k ];
      for( const auto & [ u , x ] : adj[ k ] )
       if( ( u > k ) && ( s[ u ] != s[ k ] ) )
        w += x;
      }
     }
    }
   if( cur > best_v ) {
    best_v = cur;
    best = s;
    }
   }
  return( best );
  }

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// GAZHDA*: each DTG cut in two halves by Algorithm 4 of Jinnai and
/// Fukunaga (JAIR 2017), the halves being the abstract features
/** Side 0 starts from the node of least degree and grows by the node with
 * the most edges to it until it holds \f$ \lceil n / 2 \rceil \f$ of the
 * \f$ n \f$ nodes, the ties going to the node of least index: the cut is
 * balanced (\f$ |S_2| \leq |S_1| \leq |S_2| + 1 \f$), whatever the edges
 * it cuts. The degrees are those of the DTG of the ground transitions
 * (FeatureModel::static_arcs()), as in the paper. */

[[nodiscard]] inline FeatureAbstraction greedy_abstraction(
 const FeatureGraph & m ) {
 std::vector< std::vector< unsigned char > > sides( m.num_groups() );
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g )
  sides[ g ] = detail::greedy_cut( m , g );
 return( detail::cut_abstraction( m , sides ,
                                  std::vector< bool >( m.num_groups() ) ) );
 }

/*--------------------------------------------------------------------------*/
/// GRAZHDA*: each DTG cut in the two sides of best objective, which are
/// the abstract features [see Jinnai and Fukunaga, JAIR 2017, section 6]
/** With CutObjective::sparsity the cut maximizes the sparsity
 * \f$ |S_1| |S_2| / w( S_1 , S_2 ) \f$ of equation (6) of the paper, with
 * its weights: \f$ |S_i| \f$ the nodes of \f$ S_i \f$ and \f$ w \f$ the
 * ground operators that realize the edges cut, in either direction
 * (FeatureModel::static_arcs()). With CutObjective::sampled_sparsity the
 * weights are those of the sample: \f$ w \f$ the moves of the sample that
 * cross the cut and \f$ |S_i| \f$ the expanded states whose feature of the
 * group is in \f$ S_i \f$ (FeatureModel::mass()), so that a cut weighs the
 * part of the state space that the search visits rather than the domain
 * of the variable. With CutObjective::efficiency the cut maximizes the
 * estimated efficiency \f$ 1 / ( ( 1 + c \, CO ) \, LB ) \f$ of the DTG
 * alone, with CO the fraction of the weight of the sample that is cut and
 * LB the larger mass over half of the total; the paper finds the
 * objective min( CO + LB ) comparable to the sparsity and leaves the
 * others open. A cut of weight 0 counts as one of weight \f$ 10^{-9} \f$
 * of the total, so that among the cuts of weight 0 the most balanced wins.
 *
 * The DTGs of up to 20 nodes are cut optimally, by enumerating the
 * \f$ 2^{n-1} - 1 \f$ cuts with the first node on side 0 in the order of a
 * Gray code, in which one node changes side at each step and the weight of
 * the cut changes by the edges at that node only, so that the enumeration
 * costs \f$ O( 2^n \bar{d} ) \f$ for the mean degree \f$ \bar{d} \f$ (the
 * paper uses a depth-first branch-and-bound on DTGs of fewer than 10
 * nodes). The larger DTGs are cut by the passes of Fiduccia and Mattheyses
 * (1982), which the paper suggests for them: in a pass each node changes
 * side once, the one that gives the best value first even if the value
 * drops, so that a pass can cross a valley that single improving moves
 * cannot, and the best cut met in the pass is kept; the passes start from
 * the halves of GAZHDA* and from halves grown in breadth-first order from
 * up to 16 nodes, and are repeated while they improve. The result is at
 * least as good as each start, with \f$ O( n^2 \bar{d} ) \f$ time per
 * pass. */

[[nodiscard]] inline FeatureAbstraction sparsest_cut_abstraction(
 const FeatureGraph & m , CutObjective o = CutObjective::sparsity ,
 double c = 1 ) {
 std::vector< std::vector< unsigned char > > sides( m.num_groups() );
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g )
  sides[ g ] = detail::best_cut( m , g , o , c );
 return( detail::cut_abstraction( m , sides ,
                                  std::vector< bool >( m.num_groups() ) ) );
 }

/*--------------------------------------------------------------------------*/
/// FAZHDA*: the fraction \p ignored of the groups of largest fluency
/// ignored, the others cut as by greedy_abstraction()
/** The fluency of a group is the fraction of the moves that change it [see
 * Jinnai and Fukunaga, JAIR 2017, section 4.2.2]: a group that nearly
 * every move changes would make nearly every move change the owner,
 * whatever the cut of the others. The paper counts the ground actions,
 * whereas FeatureModel::fluency() counts the moves of the sample, i.e.,
 * the ground actions weighed by how often the search applies them near
 * the root; the paper finds 0.3 the best fraction on the planning domains,
 * and \f$ \lfloor \mathrm{ignored} \cdot G \rfloor \f$ of the \f$ G \f$
 * groups are ignored, the ties going to the group of least index. */

[[nodiscard]] inline FeatureAbstraction fluency_abstraction(
 const FeatureGraph & m , double ignored = 0.3 ) {
 const std::size_t ng = m.num_groups();
 std::vector< std::size_t > by( ng );
 for( std::size_t g = 0 ; g < ng ; ++g )
  by[ g ] = g;
 std::stable_sort( by.begin() , by.end() ,
                   [ & ]( std::size_t a , std::size_t b ) {
                    return( m.fluency( a ) > m.fluency( b ) );
                    } );
 std::vector< bool > ign( ng , false );
 const auto k = std::size_t( std::max( 0.0 , ignored ) * double( ng ) );
 for( std::size_t j = 0 ; ( j < k ) && ( j < ng ) ; ++j )
  ign[ by[ j ] ] = true;
 std::vector< std::vector< unsigned char > > sides( ng );
 for( std::size_t g = 0 ; g < ng ; ++g )
  if( ! ign[ g ] )
   sides[ g ] = detail::greedy_cut( m , g );
 return( detail::cut_abstraction( m , sides , ign ) );
 }

namespace detail {

 /// the groups chosen greedily by the criterion of Zhou and Hansen (AAAI
 /// 2006) on the sample of m, until stop( features , abstract states ) of
 /// the groups chosen so far is true
 /** The abstract state of a state for a set P of groups is the set of its
  * features in the groups of P, coded by the exclusive or of random words
  * of those features; the abstract graph has the abstract states of the
  * states of the sample as nodes and an arc between the abstract states of
  * the two ends of each move of the sample that changes it. The criterion
  * is the largest out-degree of the abstract graph over its number of
  * nodes, which is evaluated for each group added to the current P, and
  * the group of least value is added (ties to the group of least fluency,
  * then of least index). stop is called before each addition with the
  * features of the groups that it would add and the abstract states that
  * the union would have. */

 template< class Stop >
 std::vector< bool > zhou_hansen_groups( const FeatureGraph & m ,
                                         Stop && stop ) {
  const std::size_t ng = m.num_groups();
  const auto & word = m.words();
  const std::size_t ns = ng ? word.size() / ng : 0;
  std::vector< std::uint64_t > key( ns , 0 ) , cand( ns );
  std::vector< bool > in( ng , false );
  std::size_t features = 0 , total = 0;
  for( std::size_t g = 0 ; g < ng ; ++g )
   total += m.members( g ).size();
  std::unordered_map< std::uint64_t , std::uint32_t > nodes;
  std::unordered_map< std::uint64_t , std::size_t > outdeg;
  struct PairHash {
   std::size_t operator()( const std::pair< std::uint64_t ,
                                            std::uint64_t > & x ) const {
    return( std::size_t( ADS::splitmix64( x.first ^
                                          ADS::splitmix64( x.second ) ) ) );
    }
   };
  std::unordered_set< std::pair< std::uint64_t , std::uint64_t > ,
                      PairHash > arcs;
  for( ;; ) {
   std::size_t best = ng , best_nodes = 0;
   double best_v = 0;
   for( std::size_t g = 0 ; g < ng ; ++g ) {
    if( in[ g ] || m.members( g ).empty() )
     continue;
    nodes.clear();
    outdeg.clear();
    arcs.clear();
    for( std::size_t k = 0 ; k < ns ; ++k ) {
     cand[ k ] = key[ k ] ^ word[ k * ng + g ];
     nodes.emplace( cand[ k ] , 0 );
     }
    std::size_t mx = 0;
    for( const auto & [ s , t ] : m.edges() )
     if( ( cand[ s ] != cand[ t ] ) &&
         arcs.emplace( cand[ s ] , cand[ t ] ).second )
      mx = std::max( mx , ++outdeg[ cand[ s ] ] );
    const double v = double( mx ) / double( nodes.size() );
    if( ( best == ng ) || ( v < best_v ) ||
        ( ( v == best_v ) && ( m.fluency( g ) < m.fluency( best ) ) ) ) {
     best = g;
     best_v = v;
     best_nodes = nodes.size();
     }
    }
   if( ( best == ng ) ||
       stop( features + m.members( best ).size() , best_nodes ) )
    break;
   in[ best ] = true;
   features += m.members( best ).size();
   for( std::size_t k = 0 ; k < ns ; ++k )
    key[ k ] ^= word[ k * ng + best ];
   if( features >= total )
    break;
   }
  return( in );
  }

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// DAHDA*: Zobrist hashing of the abstract state made of the groups chosen
/// one at a time until they hold the fraction \p ratio of the features
/** The groups are chosen by the greedy abstraction of Zhou and Hansen (AAAI
 * 2006), which adds the group that minimizes the largest out-degree of the
 * abstract graph over its number of nodes, measured here on the abstract
 * graph of the sample [see detail::zhou_hansen_groups()]; Burns et al.
 * (JAIR 2010) stop before the abstract graph exceeds a given size [see
 * state_abstraction()], whereas DAHDA* adds groups as long as they hold
 * at most the fraction \p ratio of the features (at least one group),
 * which keeps the load balanced on instances of any size [see Jinnai and
 * Fukunaga, JAIR 2017, appendix A, where 0.3 is the best fraction]. The
 * features of the other groups are ignored. */

[[nodiscard]] inline FeatureAbstraction dynamic_state_abstraction(
 const FeatureGraph & m , double ratio = 0.3 ) {
 std::size_t total = 0;
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g )
  total += m.members( g ).size();
 const double limit = ratio * double( total );
 std::size_t taken = 0;
 const auto in = detail::zhou_hansen_groups( m , [ & ]( std::size_t f ,
                                                        std::size_t ) {
  const bool over = ( taken > 0 ) && ( double( f ) > limit );
  taken = f;
  return( over );
  } );
 std::vector< bool > ign( m.num_groups() );
 for( std::size_t g = 0 ; g < ign.size() ; ++g )
  ign[ g ] = ! in[ g ];
 return( detail::cut_abstraction(
          m , std::vector< std::vector< unsigned char > >( ign.size() ) ,
          ign ) );
 }

/*--------------------------------------------------------------------------*/
/// AHDA*: as dynamic_state_abstraction(), but the groups are added as long
/// as the abstract graph of the sample has at most \p max_states nodes
/** It is the abstraction of Burns et al. (JAIR 2010) for planning, i.e.,
 * that of Zhou and Hansen (AAAI 2006) with the bound \f$ N_{\max} \f$ on
 * the abstract graph, here measured on the sample: the group whose
 * addition would exceed it is not added. */

[[nodiscard]] inline FeatureAbstraction state_abstraction(
 const FeatureGraph & m , std::size_t max_states = 1000 ) {
 const auto in = detail::zhou_hansen_groups( m , [ & ]( std::size_t ,
                                                        std::size_t n ) {
  return( n > max_states );
  } );
 std::vector< bool > ign( m.num_groups() );
 for( std::size_t g = 0 ; g < ign.size() ; ++g )
  ign[ g ] = ! in[ g ];
 return( detail::cut_abstraction(
          m , std::vector< std::vector< unsigned char > >( ign.size() ) ,
          ign ) );
 }

/*--------------------------------------------------------------------------*/
/// OZHDA*: a Zobrist table in which the most frequent operators change no
/// bit of the key
/** The table is drawn from the seed; then, for the operators by decreasing
 * weight, up to the fraction \p fraction of them, the first feature x of
 * the operator that no previous operator has fixed gets as the high 32
 * bits of \f$ R[ x ] \f$ the exclusive or of those of the other features
 * of the operator, and all of them are fixed [see Jinnai and Fukunaga,
 * JAIR 2017, section 2.4, and AAAI 2016]. Since an entry changes only
 * while none of the operators already processed holds its feature, the
 * high bits of \f$ \bigoplus_{x \in a} R[ x ] \f$ stay 0 for every
 * operator a that fixed a feature, and a move of a keeps the key, and then
 * the owner, of its parent. The paper processes all the actions, by
 * increasing number of preconditions, and sets the whole entry; here the
 * order is by frequency in the sample, the fraction (1 for all of them) is
 * a parameter, and only the 32 bits of the key are set, so that the low
 * bits stay random and the hash still tells the states apart. With
 * \p as_paper, the whole entry is set and the operators go by increasing
 * number of features changed, which stands for the preconditions of the
 * paper (a FeatureDomain gives no preconditions); ties keep the order by
 * frequency. Since each
 * operator joins its unfixed feature to the others, the operators of a
 * DTG join whole paths of it, and a large fraction balances the load
 * badly. */

[[nodiscard]] inline std::vector< std::uint64_t > operator_table(
 const FeatureGraph & m , double fraction = 0.2 ,
 std::uint64_t seed = 1 , bool as_paper = false ) {
 auto r = FeatureAbstraction::identity( m.num_features() ).table( seed );
 std::vector< bool > fixed( m.num_features() , false );
 std::vector< const FeatureGraph::Operator * > ops;
 for( const auto & a : m.operators() )
  ops.push_back( &a );
 if( as_paper )
  std::stable_sort( ops.begin() , ops.end() ,
                    []( const auto * a , const auto * b ) {
                     return( a->features.size() < b->features.size() );
                     } );
 const auto k = std::size_t( std::ceil( std::clamp( fraction , 0.0 , 1.0 ) *
                                        double( ops.size() ) ) );
 const std::uint64_t low = as_paper ? 0 : 0xFFFFFFFF;
 for( std::size_t j = 0 ; j < k ; ++j ) {
  const auto & fe = ops[ j ]->features;
  if( fe.size() < 2 )
   continue;
  const auto x = std::find_if( fe.begin() , fe.end() ,
                               [ & ]( std::uint32_t i ) {
                                return( ! fixed[ i ] );
                                } );
  if( x == fe.end() )
   continue;
  std::uint64_t hi = 0;
  for( std::uint32_t i : fe )
   if( i != *x )
    hi ^= r[ i ] & ~low;
  r[ *x ] = hi | ( r[ *x ] & low );
  for( std::uint32_t i : fe )
   fixed[ i ] = true;
  }
 return( r );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Gf2Basis ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a basis of a subspace V of \f$ GF(2)^n \f$, kept in echelon form
/** insert() reduces a vector by the rows, each of which has a pivot that
 * no previous row has, and appends what is left if it is not 0;
 * reduced() returns the rows in reduced row echelon form, where no row has
 * a bit at the pivot of another row. */

class Gf2Basis {

 public:

 /// a vector of \f$ GF(2)^n \f$, by words of 64 bits
 using Vector = std::vector< std::uint64_t >;

 /// the basis of \f$ \{ 0 \} \f$ in \f$ GF(2)^n \f$

 explicit Gf2Basis( std::size_t n ) : f_w( ( n + 63 ) / 64 ) {}

 /// the dimension of V

 [[nodiscard]] std::size_t rank( void ) const noexcept {
  return( v_row.size() );
  }

 /// adds to V the vector of the set bits \p bits; true if V grows

 bool insert( const std::vector< std::uint32_t > & bits ) {
  Vector v( f_w , 0 );
  for( std::uint32_t i : bits )
   v[ i / 64 ] ^= std::uint64_t( 1 ) << ( i % 64 );
  for( std::size_t r = 0 ; r < v_row.size() ; ++r )
   if( ( v[ v_piv[ r ] / 64 ] >> ( v_piv[ r ] % 64 ) ) & 1 )
    for( std::size_t k = 0 ; k < f_w ; ++k )
     v[ k ] ^= v_row[ r ][ k ];
  for( std::size_t k = 0 ; k < f_w ; ++k )
   if( v[ k ] ) {
    v_piv.push_back( std::uint32_t( k * 64 + std::countr_zero( v[ k ] ) ) );
    v_row.push_back( std::move( v ) );
    return( true );
    }
  return( false );
  }

 /// the rows in reduced row echelon form, and their pivots

 [[nodiscard]] std::pair< std::vector< Vector > ,
                         std::vector< std::uint32_t > > reduced( void ) const
 {
  auto row = v_row;
  // row r lacks the pivots of the rows before it; clearing from the last
  // row backwards the pivots of the later rows leaves those absent too
  for( std::size_t r = row.size() ; r-- > 0 ; )
   for( std::size_t q = r + 1 ; q < row.size() ; ++q )
    if( ( row[ r ][ v_piv[ q ] / 64 ] >> ( v_piv[ q ] % 64 ) ) & 1 )
     for( std::size_t k = 0 ; k < f_w ; ++k )
      row[ r ][ k ] ^= row[ q ][ k ];
  return( std::pair( std::move( row ) , v_piv ) );
  }

 private:

 std::size_t f_w;                          ///< the words of a vector
 std::vector< Vector > v_row;              ///< the rows
 std::vector< std::uint32_t > v_piv;       ///< the pivot of each row

 };  // end( class( Gf2Basis ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Zobrist table whose key vanishes on the subspace V of the change
/// vectors in basis
/** Each feature that is not a pivot of the reduced basis gets a random
 * entry; a pivot x of the row \f$ v \f$ gets as the high 32 bits of
 * \f$ R[ x ] \f$ the exclusive or of those of the other features of
 * \f$ v \f$, which are not pivots, and random low bits. Then
 * \f$ \bigoplus_{y \in v} R[ y ] \f$ has its high 32 bits 0 for every row,
 * hence for every vector of V, and the 32 rows of key bits are uniform in
 * the annihilator of V. */

[[nodiscard]] inline std::vector< std::uint64_t > kernel_table(
 std::size_t n , const Gf2Basis & basis , std::uint64_t seed = 1 ) {
 auto r = FeatureAbstraction::identity( n ).table( seed );
 const auto [ row , piv ] = basis.reduced();
 constexpr std::uint64_t low = 0xFFFFFFFF;
 for( std::size_t q = 0 ; q < row.size() ; ++q ) {
  std::uint64_t hi = 0;
  for( std::size_t k = 0 ; k < row[ q ].size() ; ++k )
   for( std::uint64_t w = row[ q ][ k ] ; w ; w &= w - 1 ) {
    const std::size_t x = k * 64 + std::size_t( std::countr_zero( w ) );
    if( x != piv[ q ] )
     hi ^= r[ x ] & ~low;
    }
  r[ piv[ q ] ] = hi | ( r[ piv[ q ] ] & low );
  }
 return( r );
 }

/*--------------------------------------------------------------------------*/
/// KZHDA*: the Zobrist table whose key vanishes on the change vectors of
/// operators of the sample chosen to give the largest estimated efficiency
/// among p workers
/** A state is the vector \f$ x( s ) \in GF(2)^n \f$ of its features, and
 * the key of Zobrist hashing is \f$ H x( s ) \f$, with \f$ H \f$ the
 * \f$ 32 \times n \f$ matrix of the high bits of the table: a move from s
 * to t with change vector \f$ \delta = x( s ) \oplus x( t ) \f$ keeps the
 * key, and then the owner, if and only if \f$ H \delta = 0 \f$ (with
 * \f$ p = 2^b \f$ workers the owner is the first b bits of the key, and it
 * is kept if and only if the first b rows vanish on \f$ \delta \f$). The
 * subspace \f$ V \f$ on which \f$ H \f$ vanishes [see kernel_table()] is
 * grown greedily over the \p candidates operators of largest weight, by
 * decreasing weight: an operator already in \f$ V \f$ costs nothing, and
 * another one is kept if the table of the larger \f$ V \f$ has a larger
 * value than the current one, and dropped otherwise. The value is the
 * smaller of the efficiencies estimated on the sample among p and among
 * 2 p workers [see FeatureGraph::estimate_table()]: the sample holds the
 * states nearest to the root, and a subspace that balances them among p
 * workers only by a few owner classes balances the later states badly,
 * which the estimate among 2 p workers exposes. The table of every
 * \f$ V \f$ is drawn from the same seed, a feature that is not a pivot
 * keeping its entry, so that two steps differ only in the entries that
 * the new operator forces. The start is \f$ V = \{ 0 \} \f$, Zobrist
 * hashing, or, if \p start is not nullptr, the subspace of the feature
 * projection *start, i.e., the span of the \f$ e_x \oplus e_y \f$ with
 * x and y of the same abstract feature and of the \f$ e_x \f$ with x
 * ignored, whose tables are those of Abstract Zobrist hashing with that
 * projection; the result is never estimated worse than the start. Each
 * step costs the reduction of the basis and one estimate, i.e., about the
 * features of the states of the sample. */

[[nodiscard]] inline std::vector< std::uint64_t > kernel_selection_table(
 const FeatureGraph & m , unsigned p , double c = 1 ,
 std::uint64_t seed = 1 , std::size_t candidates = 1024 ,
 const FeatureAbstraction * start = nullptr ) {
 const auto & ops = m.operators();
 Gf2Basis basis( m.num_features() );
 if( start ) {
  std::vector< std::uint32_t > last( start->num_abstract ,
                                     no_abstract_feature );
  for( std::uint32_t x = 0 ; x < start->abstract.size() ; ++x ) {
   const std::uint32_t a = start->abstract[ x ];
   if( a == no_abstract_feature )
    basis.insert( { x } );
   else {
    if( last[ a ] != no_abstract_feature )
     basis.insert( { last[ a ] , x } );
    last[ a ] = x;
    }
   }
  }
 auto best = kernel_table( m.num_features() , basis , seed );
 // the efficiency estimated among p and among 2 p workers, the smaller
 const auto value = [ & ]( const std::vector< std::uint64_t > & t ) {
  return( std::min( m.estimate_table( t , p , c ).efficiency ,
                    m.estimate_table( t , 2 * p , c ).efficiency ) );
  };
 double best_e = value( best );
 for( std::size_t j = 0 ; j < std::min( candidates , ops.size() ) ; ++j ) {
  Gf2Basis next = basis;
  if( ! next.insert( ops[ j ].features ) )
   continue;                                   // already in V
  auto t = kernel_table( m.num_features() , next , seed );
  if( const double e = value( t ) ;
      e > best_e ) {
   best_e = e;
   best = std::move( t );
   basis = std::move( next );
   }
  }
 return( best );
 }

/*--------------------------------------------------------------------------*/
/// the names of the distributions among which auto_distribution() chooses
/** Each method of this header, the ones with a parameter at a few values
 * around the default [see make_distribution()]. */

inline const std::vector< std::string > & auto_distribution_names( void ) {
 static const std::vector< std::string > n = {
  "Z" , "GAZ" , "GRAZ" , "GRAZS" , "GRAZE" , "FAZ:0.1" , "FAZ:0.3" ,
  "FAZ:0.5" ,
  "OZ:0.1" , "OZ:0.2" , "OZ:0.3" , "DAHDA:0.2" , "DAHDA:0.3" , "DAHDA:0.5" ,
  "AHDA:100" , "AHDA:1000" , "AHDA:10000" , "KZ" , "KZ:GRAZS" ,
  "KZ:GRAZE" , "KZ:DAHDA" };
 return( n );
 }

/*--------------------------------------------------------------------------*/
/// the distribution of the given name built on m
/** The names are Z (Zobrist hashing), GAZ (greedy_abstraction()), GRAZ
 * (sparsest_cut_abstraction() by the sparsity of the paper), GRAZS (by the
 * sparsity with the weights of the sample), GRAZE (by efficiency, with the
 * cost c of a message), FAZ (fluency_abstraction(), 0.3 of the groups
 * ignored), OZ (operator_table(), 0.2 of the operators), OZP (the same as
 * in the paper, all the operators by default), DAHDA
 * (dynamic_state_abstraction(), 0.3 of the features) and AHDA
 * (state_abstraction(), at most 1000 abstract states), the last five
 * followed, if they are not at their default, by a colon and the
 * parameter, e.g., OZ:0.1 or AHDA:100, and KZ (kernel_selection_table()
 * for p workers, from Zobrist hashing), or KZ:GRAZ, KZ:GRAZS, KZ:GRAZE and
 * KZ:DAHDA (from the subspace of that projection); the table is drawn from
 * the seed. Throws
 * std::invalid_argument on an unknown name. */

template< FeatureDomain D >
[[nodiscard]] TableZobristDistribution< D > make_distribution(
 const FeatureModel< D > & m , const std::string & name , double c = 1 ,
 std::uint64_t seed = 1 , unsigned p = 8 ) {
 const D & d = m.domain();
 const auto colon = name.find( ':' );
 const std::string base = name.substr( 0 , colon );
 const auto par = [ & ]( double dflt ) {
  return( colon == std::string::npos ? dflt
                                     : std::stod( name.substr( colon + 1 ) ) );
  };
 if( base == "Z" )
  return( TableZobristDistribution< D >(
           d , FeatureAbstraction::identity( m.num_features() ) , seed ,
           name ) );
 if( base == "GAZ" )
  return( TableZobristDistribution< D >( d , greedy_abstraction( m ) , seed ,
                                         name ) );
 if( base == "GRAZ" )
  return( TableZobristDistribution< D >( d , sparsest_cut_abstraction( m ) ,
                                         seed , name ) );
 if( base == "GRAZS" )
  return( TableZobristDistribution< D >(
           d , sparsest_cut_abstraction( m , CutObjective::sampled_sparsity ) ,
           seed , name ) );
 if( base == "GRAZE" )
  return( TableZobristDistribution< D >(
           d , sparsest_cut_abstraction( m , CutObjective::efficiency , c ) ,
           seed , name ) );
 if( base == "FAZ" )
  return( TableZobristDistribution< D >(
           d , fluency_abstraction( m , par( 0.3 ) ) , seed , name ) );
 if( base == "OZ" )
  return( TableZobristDistribution< D >(
           d , operator_table( m , par( 0.2 ) , seed ) , name ) );
 if( base == "KZ" ) {
  if( colon == std::string::npos )
   return( TableZobristDistribution< D >(
            d , kernel_selection_table( m , p , c , seed ) , name ) );
  const std::string from = name.substr( colon + 1 );
  FeatureAbstraction a;
  if( from == "GRAZ" )
   a = sparsest_cut_abstraction( m );
  else
   if( from == "GRAZS" )
    a = sparsest_cut_abstraction( m , CutObjective::sampled_sparsity );
   else
    if( from == "GRAZE" )
     a = sparsest_cut_abstraction( m , CutObjective::efficiency , c );
    else
     if( from == "DAHDA" )
      a = dynamic_state_abstraction( m );
     else
      throw( std::invalid_argument( "make_distribution: no start " + from ) );
  return( TableZobristDistribution< D >(
           d , kernel_selection_table( m , p , c , seed , 1024 , &a ) ,
           name ) );
  }
 if( base == "OZP" )
  return( TableZobristDistribution< D >(
           d , operator_table( m , par( 1 ) , seed , true ) , name ) );
 if( base == "AHDA" )
  return( TableZobristDistribution< D >(
           d , state_abstraction( m , std::size_t( par( 1000 ) ) ) , seed ,
           name ) );
 if( base == "DAHDA" )
  return( TableZobristDistribution< D >(
           d , dynamic_state_abstraction( m , par( 0.3 ) ) , seed , name ) );
 throw( std::invalid_argument( "make_distribution: unknown distribution " +
                               name ) );
 }

/*--------------------------------------------------------------------------*/
/// the distribution of largest estimated efficiency among p workers
/** Builds every distribution of auto_distribution_names() on m, adds those
 * in \p extra, estimates each one on the sample of m [see
 * FeatureModel::estimate()] with c the cost of sending a node relative to
 * that of generating it, and returns the one of largest estimated
 * efficiency (the first on ties), named "auto:" and its name; if \p est
 * is not nullptr, it receives its estimate. The choice rests on the
 * correlation of the estimated with the actual efficiency that Jinnai and
 * Fukunaga (JAIR 2017, section 5.2.1) find; besides choosing the method,
 * it tunes the fraction of FAZHDA*, OZHDA* and DAHDA*, which the paper
 * fixes once for all the instances. */

template< FeatureDomain D >
[[nodiscard]] TableZobristDistribution< D > auto_distribution(
 const FeatureModel< D > & m , unsigned p , double c = 1 ,
 const std::vector< TableZobristDistribution< D > > & extra = {} ,
 DistributionEstimate * est = nullptr , std::uint64_t seed = 1 ) {
 std::vector< TableZobristDistribution< D > > cand;
 for( const auto & n : auto_distribution_names() )
  cand.push_back( make_distribution( m , n , c , seed , p ) );
 for( const auto & x : extra )
  cand.push_back( x );
 std::size_t best = 0;
 DistributionEstimate be;
 for( std::size_t k = 0 ; k < cand.size() ; ++k ) {
  const auto e = m.estimate( cand[ k ] , p , c );
  if( ( k == 0 ) || ( e.efficiency > be.efficiency ) ) {
   best = k;
   be = e;
   }
  }
 if( est )
  *est = be;
 return( TableZobristDistribution< D >( m.domain() , cand[ best ].table() ,
                                        "auto:" + cand[ best ].name() ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AutoDistribution */

/*--------------------------------------------------------------------------*/
/*---------------------- End File AutoDistribution.h -----------------------*/
/*--------------------------------------------------------------------------*/
