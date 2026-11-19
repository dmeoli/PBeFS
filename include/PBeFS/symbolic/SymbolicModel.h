/*--------------------------------------------------------------------------*/
/*-------------------------- File SymbolicModel.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SymbolicModel, which compiles a
 * SymbolicDomain [see SymbolicDomain.h] into binary decision diagrams
 * [see Bdd.h]: the encoding of the states into bits, the initial state, the
 * goal, the transition relations partitioned by cost and merged up to a
 * size, and the sets of the states of each value of an additive heuristic.
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

#ifndef __PBeFS_SymbolicModel
 #define __PBeFS_SymbolicModel
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cstddef>

#include <cstdint>

#include <iterator>

#include <limits>

#include <map>

#include <stdexcept>

#include <string>

#include <type_traits>

#include <vector>

#include "PBeFS/symbolic/Bdd.h"

#include "PBeFS/symbolic/SymbolicDomain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SymbolicLayout --------------------------*/
/*--------------------------------------------------------------------------*/
/// the order of the variables of a SymbolicModel
/** order is a permutation of the variables, the first one taking the first
 * bits; if it is empty the order is that of the domain, or, with optimize,
 * the one that a local search finds for the least sum over the operators
 * (and the terms of the heuristic) of the square of the span of their
 * variables, in bits, which keeps close the variables that a move reads and
 * writes together [see docs/theory/symbolic-search.md]. */

struct SymbolicLayout {
 std::vector< std::uint32_t > order;    ///< the variables, from the first
 bool optimize = false;                 ///< search an order if none given
 };

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SymbolicModel ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a SymbolicDomain compiled into binary decision diagrams
/** The variable v of the domain, with \f$ | D_v | \f$ values, takes \f$ b_v
 * = \lceil \log_2 | D_v | \rceil \f$ bits, the most significant first, and
 * the bits of all the variables, in their order, are the bits \f$ k = 0 ,
 * 1 , \ldots \f$ of a state, the bit k of the state before a move being the
 * variable \f$ 2 k \f$ of Sylvan and that after it the variable \f$ 2 k + 1
 * \f$ [see Bdd.h]; the two copies of a bit are thus next to each other,
 * which keeps small the relations of the moves that leave a variable as it
 * is.
 *
 * The relation of an operator o is \f$ T_o( x , x' ) = \mathrm{pre}_o( x )
 * \wedge \bigwedge_{e} E_e( x , x' ) \wedge \bigwedge_{v \in R_o \setminus
 * W_o} ( x'_v = x_v ) \f$ on the variables \f$ R_o \f$ that o reads or
 * writes, \f$ W_o \f$ those it writes, \f$ E_e \f$ being \f$ x'_v = d \f$
 * for an assignment and \f$ x'_v = x_u \f$ for a copy; the other variables
 * do not appear, and the image leaves them unchanged. The operators of
 * equal cost are merged, in their order, into partitions \f$ T = \bigvee_o
 * T_o \f$ on the union of their variables, a partition being closed when
 * the next operator would bring it above max_nodes nodes (0 keeps one
 * operator per partition), the disjunctive partitioning bounded by size of
 * the symbolic planners [see Torralba et al., AIJ 2017].
 *
 * The heuristic \f$ h = \sum_t h_t \f$ of a SymbolicHeuristicDomain is
 * compiled into the sets \f$ H_k = \{ x : h( x ) = k \} \f$ by a dynamic
 * program over the terms, \f$ P_t( a ) = \bigvee_{b} P_{t-1}( a - b )
 * \wedge \{ x : h_t( x ) = b \} \f$: the diagram of \f$ P_t( a ) \f$ has at
 * most as many nodes on a level as there are partial sums, which keeps
 * small the sets of, e.g., the Manhattan distance [see Edelkamp and Reffel,
 * KI 1998]. A state with a term infinite is in no \f$ H_k \f$, and the
 * searches drop it. Without terms, \f$ H_0 \f$ is every state. */

template< SymbolicDomain D >
class SymbolicModel {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using cost_type = typename D::cost_type;   ///< the cost type

 /// a partition of the transition relation
 struct Partition {
  cost_type cost;                ///< the cost of its operators
  Bdd relation;                  ///< \f$ T( x , x' ) \f$
  Bdd vars;                      ///< the variables it reads or writes
  std::size_t operators = 0;     ///< the operators merged in it
  };

 using Levels = std::map< cost_type , Bdd >;  ///< \f$ H_k \f$ by k

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the model of d, the partitions up to max_nodes nodes; with unit_costs
 /// every operator costs 1, as in a breadth-first search, and without
 /// terms the heuristic terms of d are left out (h = 0)
 /** Throws std::invalid_argument if the domain is malformed: a value out
  * of its range, a negative cost, a variable written twice by one
  * operator, a copy into a variable of fewer bits. Needs a BddRuntime. */

 explicit SymbolicModel( const D & d , std::size_t max_nodes = 100000 ,
                         bool unit_costs = false , bool terms = true ,
                         const SymbolicLayout & layout = {} ) : f_d( d ) {
  const std::size_t n = d.num_variables();
  v_offset.resize( n , 0 );
  v_bits.resize( n , 0 );
  v_range.resize( n );
  for( std::size_t v = 0 ; v < n ; ++v ) {
   v_range[ v ] = d.range( std::uint32_t( v ) );
   if( v_range[ v ] == 0 )
    throw( std::invalid_argument( "SymbolicModel: a variable without "
                                  "values" ) );
   v_bits[ v ] = v_range[ v ] <= 1 ? 0 :
    std::size_t( std::bit_width( v_range[ v ] - 1 ) );
   }

  // the derived variables, not in the state, and the variables of the
  // state that each one depends on
  v_derived.assign( n , false );
  v_default.assign( n , 0 );
  v_layer.assign( n , 0 );
  v_depends.resize( n );
  if constexpr( SymbolicAxiomDomain< D > ) {
   for( std::size_t v = 0 ; v < n ; ++v )
    if( d.is_derived( std::uint32_t( v ) ) ) {
     if( v_range[ v ] != 2 )
      throw( std::invalid_argument( "SymbolicModel: a derived variable "
                                    "without two values" ) );
     v_derived[ v ] = true;
     v_bits[ v ] = 0;
     v_default[ v ] = d.default_value( std::uint32_t( v ) );
     check( std::uint32_t( v ) , v_default[ v ] );
     v_layer[ v ] = d.layer( std::uint32_t( v ) );
     }
   d.for_each_axiom( [ & ]( const SymbolicAxiom & a ) {
    check( a.var , a.value );
    if( ! v_derived[ a.var ] || ( a.value == v_default[ a.var ] ) )
     throw( std::invalid_argument( "SymbolicModel: an axiom of a variable "
                                   "not derived, or to its default" ) );
    for( const auto & f : a.cond )
     check( f.var , f.value );
    v_axioms.push_back( a );
    } );
   // the dependencies, to their closure
   for( bool more = true ; more ; ) {
    more = false;
    for( const auto & a : v_axioms )
     for( const auto & f : a.cond ) {
      std::vector< std::uint32_t > add;
      if( v_derived[ f.var ] )
       add = v_depends[ f.var ];
      else
       add = { f.var };
      for( auto u : add )
       if( std::find( v_depends[ a.var ].begin() ,
                      v_depends[ a.var ].end() , u ) ==
           v_depends[ a.var ].end() ) {
        v_depends[ a.var ].push_back( u );
        more = true;
        }
      }
    }
   }

  // the operators by cost, in their order
  d.for_each_operator( [ & ]( const SymbolicOperator< cost_type > & o ) {
   if constexpr( std::is_signed_v< cost_type > )
    if( o.cost < 0 )
     throw( std::invalid_argument( "SymbolicModel: a negative cost" ) );
   support( o , nullptr );
   v_ops[ unit_costs ? cost_type( 1 ) : o.cost ].push_back( o );
   ++f_operators;
   } );

  // the order of the variables, and the first bit of each
  v_order = layout.order;
  if( v_order.empty() ) {
   for( std::size_t v = 0 ; v < n ; ++v )
    v_order.push_back( std::uint32_t( v ) );
   if( layout.optimize )
    optimize_order( d , terms );
   }
  {
   std::vector< bool > seen( n , false );
   for( auto v : v_order ) {
    if( ( v >= n ) || seen[ v ] )
     throw( std::invalid_argument( "SymbolicModel: the order is not a "
                                   "permutation of the variables" ) );
    seen[ v ] = true;
    }
   if( v_order.size() != n )
    throw( std::invalid_argument( "SymbolicModel: the order is not a "
                                  "permutation of the variables" ) );
   }
  for( auto v : v_order ) {
   v_offset[ v ] = f_bits;
   f_bits += v_bits[ v ];
   }
  std::vector< std::uint32_t > cur;
  for( std::size_t k = 0 ; k < f_bits ; ++k )
   cur.push_back( std::uint32_t( 2 * k ) );
  v_state_vars = cur;
  f_state_set = Bdd::variables( cur );

  f_valid = Bdd::one();
  for( std::size_t v = 0 ; v < n ; ++v )
   if( ! v_derived[ v ] )
    f_valid &= below( v , v_range[ v ] );

  // the states where each derived variable has its other value, layer by
  // layer, each to its least fixed point
  v_other.assign( n , Bdd() );
  if( ! v_axioms.empty() ) {
   int top = 0;
   for( std::size_t v = 0 ; v < n ; ++v )
    if( v_derived[ v ] )
     top = std::max( top , v_layer[ v ] );
   for( int L = 0 ; L <= top ; ++L )
    for( bool more = true ; more ; ) {
     more = false;
     for( const auto & a : v_axioms ) {
      if( v_layer[ a.var ] != L )
       continue;
      Bdd C = f_valid;
      for( const auto & f : a.cond )
       C &= value( f.var , f.value , false );
      const Bdd X = v_other[ a.var ] | C;
      if( ! ( X == v_other[ a.var ] ) ) {
       v_other[ a.var ] = X;
       more = true;
       }
      }
     }
   }

  const auto x0 = d.initial_values();
  if( x0.size() != n )
   throw( std::invalid_argument( "SymbolicModel: an initial value for each "
                                 "variable" ) );
  f_init = state( x0 );

  f_goal = Bdd();
  d.for_each_goal( [ & ]( const std::vector< SymbolicFact > & p ) {
   Bdd g = Bdd::one();
   for( const auto & f : p )
    g &= value( check( f.var , f.value ) , f.value , false );
   f_goal |= g;
   } );
  f_goal &= f_valid;

  f_max_nodes = max_nodes;
  for( const auto & [ c , list ] : v_ops )
   merge( c , list , max_nodes , nullptr , v_part );

  v_forward = { { cost_type( 0 ) , Bdd::one() } };
  v_backward = v_forward;
  if( ! terms )
   return;
  if constexpr( SymbolicHeuristicDomain< D > ) {
   std::vector< SymbolicTerm< cost_type > > t;
   d.for_each_heuristic_term( [ & ]( const SymbolicTerm< cost_type > & x ) {
    t.push_back( x );
    } );
   v_forward = levels( t );
   f_has_forward = true;
   }
  if constexpr( SymbolicBackwardHeuristicDomain< D > ) {
   std::vector< SymbolicTerm< cost_type > > t;
   d.for_each_backward_heuristic_term(
    [ & ]( const SymbolicTerm< cost_type > & x ) { t.push_back( x ); } );
   v_backward = levels( t );
   f_has_backward = true;
   }
  }

/*------------------------------- THE STATES -------------------------------*/
 /// the domain

 [[nodiscard]] const D & domain( void ) const noexcept { return( f_d ); }

 /// the variables of the domain

 [[nodiscard]] std::size_t num_variables( void ) const noexcept {
  return( v_range.size() );
  }

 /// the order of the variables, from the first
 [[nodiscard]] const std::vector< std::uint32_t > & order( void ) const
  noexcept {
  return( v_order );
  }

 /// the bits of a state

 [[nodiscard]] std::size_t num_bits( void ) const noexcept {
  return( f_bits );
  }

 /// the variables of Sylvan of a state, sorted

 [[nodiscard]] const std::vector< std::uint32_t > & state_vars( void ) const
  noexcept {
  return( v_state_vars );
  }

 /// the set of the variables of a state, for count()

 [[nodiscard]] const Bdd & state_set( void ) const noexcept {
  return( f_state_set );
  }

 /// the states whose variable v has the value d, after the move if next

 [[nodiscard]] Bdd value( std::size_t v , std::size_t d , bool next ) const {
  if( v_derived[ v ] ) {
   if( next )
    throw( std::invalid_argument( "SymbolicModel: a derived variable "
                                  "written" ) );
   return( d == v_default[ v ] ? ! v_other[ v ] : v_other[ v ] );
   }
  const std::size_t b = v_bits[ v ];
  if( ( b < 64 ) && ( d >> b ) )
   return( Bdd() );
  std::vector< std::uint32_t > vars( b );
  std::vector< bool > bits( b );
  for( std::size_t j = 0 ; j < b ; ++j ) {
   vars[ j ] = std::uint32_t( 2 * ( v_offset[ v ] + j ) + ( next ? 1 : 0 ) );
   bits[ j ] = ( d >> ( b - 1 - j ) ) & 1;
   }
  return( Bdd::cube( vars , bits ) );
  }

 /// the state of the values x

 [[nodiscard]] Bdd state( const std::vector< std::uint32_t > & x ) const {
  std::vector< bool > bits( num_bits() );
  for( std::size_t v = 0 ; v < x.size() ; ++v ) {
   check( std::uint32_t( v ) , x[ v ] );
   const std::size_t b = v_bits[ v ];
   for( std::size_t j = 0 ; j < b ; ++j )
    bits[ v_offset[ v ] + j ] = ( x[ v ] >> ( b - 1 - j ) ) & 1;
   }
  return( Bdd::cube( v_state_vars , bits ) );
  }

 /// the values of the variables in one state of S, which is not empty

 [[nodiscard]] std::vector< std::uint32_t > decode( const Bdd & S ) const {
  const auto bits = S.pick( v_state_vars );
  if( bits.empty() )
   throw( std::invalid_argument( "SymbolicModel::decode: an empty set" ) );
  std::vector< std::uint32_t > x( num_variables() , 0 );
  for( std::size_t v = 0 ; v < x.size() ; ++v )
   for( std::size_t k = v_offset[ v ] ; k < v_offset[ v ] + v_bits[ v ] ;
        ++k )
    x[ v ] = ( x[ v ] << 1 ) | ( bits[ k ] ? 1 : 0 );
  // the derived variables, which are not in the state
  for( std::size_t v = 0 ; v < x.size() ; ++v )
   if( v_derived[ v ] )
    x[ v ] = v_other[ v ].disjoint( state( x ) ) ? v_default[ v ]
                                                 : 1 - v_default[ v ];
  return( x );
  }

 /// the states in S, as a number

 [[nodiscard]] double count( const Bdd & S ) const {
  return( S.count( f_state_set ) );
  }

 /// the initial state

 [[nodiscard]] const Bdd & initial( void ) const noexcept { return( f_init ); }

 /// the goal states

 [[nodiscard]] const Bdd & goal( void ) const noexcept { return( f_goal ); }

 /// the states whose every variable is within its range

 [[nodiscard]] const Bdd & valid( void ) const noexcept { return( f_valid ); }

/*---------------------------- THE RELATIONS -------------------------------*/
 /// the partitions, by increasing cost

 [[nodiscard]] const std::vector< Partition > & partitions( void ) const
  noexcept {
  return( v_part );
  }

 /// the operators of the domain

 [[nodiscard]] std::size_t num_operators( void ) const noexcept {
  return( f_operators );
  }

 /// the least cost of an operator, 0 if there is none

 [[nodiscard]] cost_type min_cost( void ) const noexcept {
  return( v_part.empty() ? cost_type( 0 ) : v_part.front().cost );
  }

 /// the nodes of all the partitions

 [[nodiscard]] std::size_t relation_nodes( void ) const {
  std::vector< Bdd > r;
  for( const auto & p : v_part )
   r.push_back( p.relation );
  return( nodes( r ) );
  }

/*---------------------------- THE HEURISTICS ------------------------------*/
 /// the sets \f$ H_k \f$ of the heuristic towards the goal (towards the
 /// initial state if backward)

 [[nodiscard]] const Levels & levels( bool backward ) const noexcept {
  return( backward ? v_backward : v_forward );
  }

 /// true if the domain gives the heuristic of the direction

 [[nodiscard]] bool has_heuristic( bool backward ) const noexcept {
  return( backward ? f_has_backward : f_has_forward );
  }

 /// the nodes of the sets of the heuristic of the direction

 [[nodiscard]] std::size_t heuristic_nodes( bool backward ) const {
  std::vector< Bdd > r;
  for( const auto & [ k , H ] : levels( backward ) )
   r.push_back( H );
  return( nodes( r ) );
  }

 /// true if the heuristic of the direction is consistent and zero on the
 /// roots of the other direction (the goals, or the initial state if
 /// backward), checked on the states that the search of the direction
 /// reaches (from the initial state by images, or from the goals by
 /// preimages if backward)
 /** For every partition of cost c and every k, the image of \f$ H_k \f$
  * (its preimage if backward) must not meet \f$ H_j \f$ for any \f$ j < k
  * - c \f$, i.e., \f$ h( x ) \leq c + h( y ) \f$ along every move from x to
  * y, and the target must lie in \f$ H_0 \f$ (or have no level, a dead
  * end being no target): together these give the admissibility [see
  * docs/theory/symbolic-search.md]. The states reachable are found first,
  * by a breadth-first search of the sets; the encodings that are not
  * states of the domain (two pancakes in one position, say) are left out,
  * since no search reaches them. */

 [[nodiscard]] bool consistent( bool backward ) const {
  const Levels & H = levels( backward );
  Bdd reached = backward ? f_goal : f_init , frontier = reached;
  while( ! frontier.empty() ) {
   Bdd next;
   for( const auto & p : v_part )
    next |= backward ? preimage( frontier , p.relation , p.vars ) & f_valid
                     : image( frontier , p.relation , p.vars );
   frontier = next - reached;
   reached |= frontier;
   }
  const Bdd target = ( backward ? f_init : f_goal ) & reached;
  for( const auto & [ k , Hk ] : H )
   if( ( k > 0 ) && ! target.disjoint( Hk ) )
    return( false );
  for( const auto & p : v_part )
   for( const auto & [ k , Hk ] : H ) {
    const Bdd X = Hk & reached;
    const Bdd Y = backward ? preimage( X , p.relation , p.vars )
                           : image( X , p.relation , p.vars );
    for( const auto & [ j , Hj ] : H ) {
     if( j + p.cost >= k )
      break;
     if( ! Y.disjoint( Hj ) )
      return( false );
     }
    }
  return( true );
  }

/*------------------------- THE ABSTRACT HEURISTICS ------------------------*/
 /// the heuristic of the direction raised to the largest of itself, of the
 /// pattern database of each pattern and of the perimeter of the given sets
 /** A pattern P is a set of variables, and the projection of a state x on
  * it keeps the values \f$ x_P \f$; the operators are projected by
  * dropping their facts and their effects outside P (a copy from a variable
  * outside P leaves its variable free). The pattern database is the
  * distance in the projected space from \f$ x_P \f$ to the projections of
  * the targets (the goals, or the initial state if backward), computed by
  * Dijkstra's algorithm on sets backward (forward if backward), and its
  * sets \f$ H_k \f$ are on the bits of P alone [see Edelkamp and
  * Kissmann, IJCAI 2009]. With perimeter > 0 the search in the original
  * space first expands that many sets from the targets: the states it
  * closes with cost g get h = g, the others h = R, R being the least g of
  * its frontier (the states reached and not closed), and, for each
  * pattern, the least over the frontier of g plus the abstract distance,
  * the projected search starting from the projections of the sets of the
  * frontier with their g; this is the perimeter abstraction of Felner and
  * Ofek (SARA 2007), with the frontier of SymBA* as the perimeter [see
  * Torralba et al., IPC 2014]. Each of these heuristics is consistent and
  * admissible [see docs/theory/symbolic-search.md], and so is their
  * largest. With additive, the pattern databases are added instead, which
  * needs the patterns disjoint and no operator writing in two of them
  * (std::invalid_argument otherwise), and their sum raises the heuristic
  * as a whole; the perimeter is then used alone. */

 void strengthen( bool backward ,
                  const std::vector< std::vector< std::uint32_t > > &
                   patterns , std::size_t perimeter = 0 ,
                  bool additive = false ) {
  if( additive ) {
   check_additive( patterns );
   const std::vector< std::pair< cost_type , Bdd > > seeds = {
    { cost_type( 0 ) , backward ? f_init : f_goal } };
   Levels S = { { cost_type( 0 ) , Bdd::one() } };
   for( const auto & P : patterns )
    S = add( S , abstraction( P , backward , seeds ) );
   raise( backward , S );
   if( perimeter > 0 )
    strengthen( backward , {} , perimeter , false );
   return;
   }
  if( perimeter == 0 ) {
   const std::vector< std::pair< cost_type , Bdd > > seeds = {
    { cost_type( 0 ) , backward ? f_init : f_goal } };
   for( const auto & P : patterns )
    raise( backward , abstraction( P , backward , seeds ) );
   return;
   }
  // the states closed get their distance, and the others R, or the least
  // over the states of the frontier F of g( F ) + the abstract distance to
  // F, for each pattern
  std::vector< std::pair< cost_type , Bdd > > frontier;
  cost_type R{};
  const auto closed = around( backward , perimeter , frontier , R );
  Levels L;
  Bdd inside;
  for( const auto & [ g , C ] : closed ) {
   auto & x = L[ g ];
   x |= C;
   inside |= C;
   }
  Levels M = L;
  if( R != std::numeric_limits< cost_type >::max() ) {
   auto & x = L[ R ];
   x |= f_valid - inside;
   }
  raise( backward , L );
  for( const auto & P : patterns ) {
   Levels A = M;
   for( const auto & [ k , X ] : abstraction( P , backward , frontier ) ) {
    const Bdd Y = X - inside;
    if( ! Y.empty() ) {
     auto & a = A[ k ];
     a |= Y;
     }
    }
   raise( backward , A );
   }
  }

 /// the relational products and the sets expanded by the searches of the
 /// abstract heuristics and of the perimeter

 [[nodiscard]] std::uint64_t abstraction_images( void ) const noexcept {
  return( f_abs_images );
  }

 [[nodiscard]] std::uint64_t abstraction_steps( void ) const noexcept {
  return( f_abs_steps );
  }

 /// the heuristic of the direction raised with the perimeter of a blind
 /// search from its targets that is under way, with the pattern that a
 /// greedy choice finds for it; true if a pattern is found
 /** closed holds the sets that the blind search (Dijkstra's algorithm on
  * sets from the targets, by preimages, or by images if backward) has
  * expanded, each with its g, and frontier its open sets without the closed
  * states, each with its g, the least g being R: every state at a distance
  * below R is then closed, as the perimeter of strengthen() needs. The
  * pattern starts from all the variables of the state and loses, one at a
  * time, the variable whose removal leaves the smallest projection of the
  * frontier, until the projected search from the frontier keeps every set
  * within budget nodes; this is the relaxation of the frontier by
  * abstraction of SymBA* [see Torralba et al., IPC 2014], and the heuristic
  * is \f$ \max( h , h_C , h_{C,P} ) \f$ [see
  * docs/theory/symbolic-search.md]. */

 bool strengthen_from( bool backward ,
                       const std::vector< std::pair< cost_type , Bdd > > &
                        closed ,
                       const std::vector< std::pair< cost_type , Bdd > > &
                        frontier , std::size_t budget ) {
  if( frontier.empty() )
   return( false );
  const cost_type R = frontier.front().first;
  Levels C;
  Bdd inside;
  for( const auto & [ g , X ] : closed ) {
   auto & x = C[ g ];
   x |= X;
   inside |= X;
   }
  std::vector< Bdd > all;
  for( const auto & [ g , X ] : frontier )
   all.push_back( X );
  const Bdd F = disjoin( all );
  std::vector< std::uint32_t > P;
  for( std::size_t v = 0 ; v < num_variables() ; ++v )
   if( ! v_derived[ v ] && ( v_bits[ v ] > 0 ) )
    P.push_back( std::uint32_t( v ) );
  // the bits of a set of variables, for the projections
  const auto bits_of = [ & ]( const std::vector< std::uint32_t > & Q ) {
   std::vector< std::uint32_t > b;
   for( auto v : Q )
    for( std::size_t k = v_offset[ v ] ; k < v_offset[ v ] + v_bits[ v ] ;
         ++k )
     b.push_back( std::uint32_t( 2 * k ) );
   return( Bdd::variables( b ) );
   };
  while( ! P.empty() ) {
   if( F.exists( bits_of( complement( P ) ) ).nodes() <= budget ) {
    const Levels A = abstraction( P , backward , frontier , budget );
    if( ! A.empty() ) {
     Levels L = C , M = C;
     for( const auto & [ k , X ] : A ) {
      const Bdd Y = X - inside;
      if( ! Y.empty() ) {
       auto & m = M[ k ];
       m |= Y;
       }
      }
     auto & r = L[ R ];
     r |= f_valid - inside;
     raise( backward , L );
     raise( backward , M );
     return( true );
     }
    }
   // the variable whose removal leaves the smallest projected frontier
   std::size_t best = 0 , least = std::numeric_limits< std::size_t >::max();
   for( std::size_t i = 0 ; i < P.size() ; ++i ) {
    auto Q = P;
    Q.erase( Q.begin() + std::ptrdiff_t( i ) );
    const std::size_t n = F.exists( bits_of( complement( Q ) ) ).nodes();
    if( n < least ) {
     least = n;
     best = i;
     }
    }
   P.erase( P.begin() + std::ptrdiff_t( best ) );
   }
  return( false );
  }

/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/

 private:

 /// the closed sets, with their g, of Dijkstra's algorithm from the targets
 /// of the direction (preimages, or images if backward) stopped after steps
 /// sets, in frontier its open sets (without the closed states) with their
 /// g, and in R the least of those g (the largest cost if none)
 std::vector< std::pair< cost_type , Bdd > > around(
  bool backward , std::size_t steps ,
  std::vector< std::pair< cost_type , Bdd > > & frontier , cost_type & R ) {
  std::map< cost_type , Bdd > open;
  open[ cost_type( 0 ) ] = backward ? f_init : f_goal;
  Bdd closed;
  std::vector< std::pair< cost_type , Bdd > > out;
  while( ! open.empty() && ( steps > 0 ) ) {
   const auto it = open.begin();
   const cost_type g = it->first;
   const Bdd S = it->second - closed;
   open.erase( it );
   if( S.empty() )
    continue;
   --steps;
   ++f_abs_steps;
   closed |= S;
   out.emplace_back( g , S );
   for( const auto & p : v_part ) {
    Bdd Y = ( backward ? image( S , p.relation , p.vars )
                       : preimage( S , p.relation , p.vars ) ) & f_valid;
    ++f_abs_images;
    Y -= closed;
    if( ! Y.empty() ) {
     auto & b = open[ g + p.cost ];
     b |= Y;
     }
    }
   }
  // every state at a distance below R is closed, each with its distance
  frontier.clear();
  for( const auto & [ g , X ] : open ) {
   const Bdd Y = X - closed;
   if( ! Y.empty() )
    frontier.emplace_back( g , Y );
   }
  R = frontier.empty() ? std::numeric_limits< cost_type >::max()
                       : frontier.front().first;
  return( out );
  }

 /// the sets \f$ H_k \f$ of the pattern database of P, from the seeds
 /// (empty if budget > 0 and a set of the abstract search exceeds budget
 /// nodes)
 Levels abstraction( const std::vector< std::uint32_t > & P , bool backward ,
                     const std::vector< std::pair< cost_type , Bdd > > &
                      seeds , std::size_t budget = 0 ) {
  std::vector< bool > in( num_variables() , false );
  for( auto v : P )
   in[ check( v , 0 ) ] = true;
  std::vector< Partition > part;
  for( const auto & [ c , list ] : v_ops )
   merge( c , list , f_max_nodes , & in , part );
  // the bits outside P, quantified away from the seeds
  std::vector< std::uint32_t > out;
  Bdd VP = Bdd::one();
  for( std::size_t v = 0 ; v < num_variables() ; ++v )
   if( in[ v ] )
    VP &= below( v , v_range[ v ] );
   else
    for( std::size_t k = v_offset[ v ] ; k < v_offset[ v ] + v_bits[ v ] ;
        ++k )
     out.push_back( std::uint32_t( 2 * k ) );
  const Bdd outside = Bdd::variables( out );
  std::map< cost_type , Bdd > open;
  for( const auto & [ g , X ] : seeds ) {
   auto & b = open[ g ];
   b |= X.exists( outside );
   }
  Bdd closed;
  Levels L;
  while( ! open.empty() ) {
   const auto it = open.begin();
   const cost_type g = it->first;
   const Bdd S = it->second - closed;
   open.erase( it );
   if( S.empty() )
    continue;
   ++f_abs_steps;
   if( ( budget > 0 ) && ( S.nodes() > budget ) )
    return( Levels() );
   closed |= S;
   auto & l = L[ g ];
   l |= S;
   for( const auto & p : part ) {
    Bdd Y = ( backward ? image( S , p.relation , p.vars )
                       : preimage( S , p.relation , p.vars ) ) & VP;
    ++f_abs_images;
    Y -= closed;
    if( ! Y.empty() ) {
     auto & b = open[ g + p.cost ];
     b |= Y;
     }
    }
   }
  return( L );
  }

 /// the variables of the state not in P (sorted)
 std::vector< std::uint32_t > complement( const std::vector< std::uint32_t > &
                                          P ) const {
  std::vector< std::uint32_t > Q;
  for( std::size_t v = 0 ; v < num_variables() ; ++v )
   if( ! v_derived[ v ] &&
       ( std::find( P.begin() , P.end() , std::uint32_t( v ) ) == P.end() ) )
    Q.push_back( std::uint32_t( v ) );
  return( Q );
  }

 /// std::invalid_argument unless the patterns are disjoint and no operator
 /// writes in two of them
 void check_additive( const std::vector< std::vector< std::uint32_t > > &
                      patterns ) const {
  std::vector< int > of( num_variables() , -1 );
  for( std::size_t i = 0 ; i < patterns.size() ; ++i )
   for( auto v : patterns[ i ] ) {
    if( of[ check( v , 0 ) ] >= 0 )
     throw( std::invalid_argument( "SymbolicModel: additive patterns that "
                                   "share a variable" ) );
    of[ v ] = int( i );
    }
  for( const auto & [ c , list ] : v_ops )
   for( const auto & o : list ) {
    int seen = -1;
    for( const auto & e : o.effects )
     if( of[ e.var ] >= 0 ) {
      if( ( seen >= 0 ) && ( seen != of[ e.var ] ) )
       throw( std::invalid_argument( "SymbolicModel: an operator writes in "
                                     "two additive patterns" ) );
      seen = of[ e.var ];
      }
    }
  }

 /// the sets of the sum of the heuristics of A and of B
 static Levels add( const Levels & A , const Levels & B ) {
  Levels S;
  for( const auto & [ a , X ] : A )
   for( const auto & [ b , Y ] : B ) {
    const Bdd Z = X & Y;
    if( ! Z.empty() ) {
     auto & z = S[ a + b ];
     z |= Z;
     }
    }
  return( S );
  }

 /// the heuristic of the direction raised to its largest with L
 void raise( bool backward , const Levels & L ) {
  Levels & H = backward ? v_backward : v_forward;
  Levels M;
  // the states with h of H below k, and with that of L below k
  Bdd belowH , belowL;
  auto i = H.begin();
  auto j = L.begin();
  while( ( i != H.end() ) || ( j != L.end() ) ) {
   const cost_type k = ( j == L.end() ) || ( ( i != H.end() ) &&
                                             ( i->first < j->first ) )
                       ? i->first : j->first;
   const bool hi = ( i != H.end() ) && ( i->first == k );
   const bool lj = ( j != L.end() ) && ( j->first == k );
   Bdd m;
   if( hi )
    m |= i->second & belowL;
   if( lj )
    m |= j->second & belowH;
   if( hi && lj )
    m |= i->second & j->second;
   if( ! m.empty() )
    M[ k ] = m;
   if( hi ) {
    belowH |= i->second;
    ++i;
    }
   if( lj ) {
    belowL |= j->second;
    ++j;
    }
   }
  H = std::move( M );
  ( backward ? f_has_backward : f_has_forward ) = true;
  }

 /// the order of the least sum of the squared spans, in bits, of the
 /// variables of each operator and of each term, by moving one variable at a
 /// time to its best place while that lowers the sum
 void optimize_order( const D & d , bool terms ) {
  std::vector< std::vector< std::uint32_t > > groups;
  for( const auto & [ c , list ] : v_ops )
   for( const auto & o : list ) {
    auto U = support( o , nullptr );
    if( U.size() > 1 )
     groups.push_back( std::move( U ) );
    }
  if constexpr( SymbolicHeuristicDomain< D > )
   if( terms )
    d.for_each_heuristic_term( [ & ]( const SymbolicTerm< cost_type > & t ) {
     if( t.vars.size() > 1 )
      groups.push_back( t.vars );
     } );
  const std::size_t n = v_order.size();
  const auto cost = [ & ]( const std::vector< std::uint32_t > & ord ) {
   std::vector< std::size_t > at( n );
   std::size_t b = 0;
   for( auto v : ord ) {
    at[ v ] = b;
    b += v_bits[ v ];
    }
   double sum = 0;
   for( const auto & U : groups ) {
    std::size_t lo = b , hi = 0;
    for( auto v : U ) {
     lo = std::min( lo , at[ v ] );
     hi = std::max( hi , at[ v ] + v_bits[ v ] );
     }
    sum += double( hi - lo ) * double( hi - lo );
    }
   return( sum );
   };
  double best = cost( v_order );
  for( bool better = true ; better ; ) {
   better = false;
   for( std::size_t i = 0 ; i < n ; ++i )
    for( std::size_t j = 0 ; j < n ; ++j ) {
     if( i == j )
      continue;
     auto ord = v_order;
     const auto v = ord[ i ];
     ord.erase( ord.begin() + std::ptrdiff_t( i ) );
     ord.insert( ord.begin() + std::ptrdiff_t( j ) , v );
     const double c = cost( ord );
     if( c < best ) {
      best = c;
      v_order = std::move( ord );
      better = true;
      }
     }
   }
  }

 /// d, after checking that it is a value of v
 std::uint32_t check( std::uint32_t v , std::uint32_t d ) const {
  if( ( v >= v_range.size() ) || ( d >= v_range[ v ] ) )
   throw( std::invalid_argument( "SymbolicModel: the value " +
                                 std::to_string( d ) + " of the variable " +
                                 std::to_string( v ) + " out of range" ) );
  return( v );
  }

 /// the values of v below r
 Bdd below( std::size_t v , std::size_t r ) const {
  const std::size_t b = v_bits[ v ];
  if( ( b >= 64 ) || ( r >= ( std::size_t( 1 ) << b ) ) )
   return( Bdd::one() );
  Bdd s;
  for( std::size_t d = 0 ; d < r ; ++d )
   s |= value( v , d , false );
  return( s );
  }

 /// \f$ x'_t = x_s \f$, the missing high bits of x_s taken as 0
 Bdd copy( std::size_t t , std::size_t s ) const {
  const std::size_t bt = v_bits[ t ];
  const std::size_t bs = v_bits[ s ];
  if( bs > bt )
   throw( std::invalid_argument( "SymbolicModel: a copy into a variable of "
                                 "fewer bits" ) );
  Bdd r = Bdd::one();
  for( std::size_t j = 0 ; j < bt ; ++j ) {
   const Bdd next = Bdd::variable( std::uint32_t( 2 * ( v_offset[ t ] + j )
                                                  + 1 ) );
   if( j < bt - bs )
    r &= ! next;
   else
    r &= next.equiv( Bdd::variable(
     std::uint32_t( 2 * ( v_offset[ s ] + j - ( bt - bs ) ) ) ) );
   }
  return( r );
  }

 /// the variables of Sylvan (both copies) of the variables in U
 Bdd vars_of( const std::vector< std::uint32_t > & U ) const {
  std::vector< std::uint32_t > w;
  for( auto v : U )
   for( std::size_t k = v_offset[ v ] ; k < v_offset[ v ] + v_bits[ v ] ;
        ++k ) {
    w.push_back( std::uint32_t( 2 * k ) );
    w.push_back( std::uint32_t( 2 * k + 1 ) );
    }
  return( Bdd::variables( w ) );
  }

 /// the variables that o reads or writes, sorted, and checked; with a
 /// pattern in (true for its variables), those of its projection, empty if
 /// the projection writes nothing
 std::vector< std::uint32_t > support( const SymbolicOperator< cost_type > &
                                       o , const std::vector< bool > * in )
  const {
  const auto keep = [ & ]( std::uint32_t v ) { return( ! in || ( *in )[ v ] );
   };
  std::vector< std::uint32_t > U , W;
  // a fact on a derived variable reads what the variable depends on, and
  // the projections drop it
  const auto reads = [ & ]( const SymbolicFact & f ) {
   check( f.var , f.value );
   if( v_derived[ f.var ] ) {
    if( ! in )
     U.insert( U.end() , v_depends[ f.var ].begin() ,
               v_depends[ f.var ].end() );
    }
   else
    if( keep( f.var ) )
     U.push_back( f.var );
   };
  for( const auto & f : o.pre )
   reads( f );
  for( const auto & e : o.effects ) {
   if( e.copy ) {
    check( e.var , 0 );
    check( e.source , 0 );
    if( ! e.table.empty() ) {
     if( e.table.size() != v_range[ e.source ] )
      throw( std::invalid_argument( "SymbolicModel: a map whose table is "
                                    "not of the range of its source" ) );
     for( auto d : e.table )
      check( e.var , d );
     }
    }
   else
    check( e.var , e.value );
   for( const auto & f : e.cond )
    check( f.var , f.value );
   if( v_derived[ e.var ] )
    throw( std::invalid_argument( "SymbolicModel: an operator writes a "
                                  "derived variable" ) );
   if( ! keep( e.var ) )
    continue;
   if( e.copy && keep( e.source ) )
    U.push_back( e.source );
   for( const auto & f : e.cond )
    reads( f );
   U.push_back( e.var );
   W.push_back( e.var );
   }
  if( in && W.empty() )
   return( std::vector< std::uint32_t >() );
  std::sort( W.begin() , W.end() );
  if( std::adjacent_find( W.begin() , W.end() ) != W.end() )
   throw( std::invalid_argument( "SymbolicModel: a variable written twice "
                                 "by an operator" ) );
  std::sort( U.begin() , U.end() );
  U.erase( std::unique( U.begin() , U.end() ) , U.end() );
  return( U );
  }

 /// \f$ T_o \f$ on the variables U, which hold its support
 /// (of its projection on the pattern in, if any: the facts and the
 /// effects outside the pattern dropped, and a copy from outside it leaving
 /// its variable free)
 Bdd relation( const SymbolicOperator< cost_type > & o ,
               const std::vector< std::uint32_t > & U ,
               const std::vector< bool > * in ) const {
  const auto keep = [ & ]( std::uint32_t v ) {
   return( v_derived[ v ] ? ! in : ( ! in || ( *in )[ v ] ) );
   };
  Bdd r = Bdd::one();
  for( const auto & f : o.pre )
   if( keep( f.var ) )
    r &= value( f.var , f.value , false );
  std::vector< bool > written( num_variables() , false );
  for( const auto & e : o.effects ) {
   if( ! keep( e.var ) )
    continue;
   Bdd E = Bdd::one();
   if( ! e.copy )
    E = value( e.var , e.value , true );
   else
    if( keep( e.source ) ) {
     if( e.table.empty() )
      E = copy( e.var , e.source );
     else {
      E = Bdd();
      for( std::size_t d = 0 ; d < e.table.size() ; ++d )
       E |= value( e.source , d , false ) & value( e.var , e.table[ d ] ,
                                                   true );
      }
     }
   // a conditional effect, or, if a fact of its condition is dropped by
   // the projection, an effect that may or may not happen
   if( ! e.cond.empty() ) {
    Bdd C = Bdd::one();
    bool dropped = false;
    for( const auto & f : e.cond )
     if( keep( f.var ) )
      C &= value( f.var , f.value , false );
     else
      dropped = true;
    const Bdd same = copy( e.var , e.var );
    E = dropped ? ( same | ( C & E ) ) : ( ( C & E ) | ( ( ! C ) & same ) );
    }
   r &= E;
   written[ e.var ] = true;
   }
  for( auto v : U )
   if( ! written[ v ] )
    r &= copy( v , v );
  return( r );
  }

 /// the frame \f$ x'_v = x_v \f$ of the variables of U not in V
 Bdd frame( const std::vector< std::uint32_t > & U ,
            const std::vector< std::uint32_t > & V ) const {
  Bdd r = Bdd::one();
  for( auto v : U )
   if( ! std::binary_search( V.begin() , V.end() , v ) )
    r &= copy( v , v );
  return( r );
  }

 /// the partitions of the operators of cost c (of their projections on
 /// the pattern in, if any), appended to out
 void merge( cost_type c ,
             const std::vector< SymbolicOperator< cost_type > > & list ,
             std::size_t max_nodes , const std::vector< bool > * in ,
             std::vector< Partition > & out ) const {
  std::vector< std::uint32_t > U;    // the variables of the open partition
  Bdd T;                             // its relation
  std::size_t count = 0;             // its operators
  const auto close = [ & ]( void ) {
   if( count > 0 )
    out.push_back( { c , T , vars_of( U ) , count } );
   };
  for( const auto & o : list ) {
   const auto Uo = support( o , in );
   if( Uo.empty() )
    continue;
   if( count > 0 ) {
    std::vector< std::uint32_t > V;
    std::set_union( U.begin() , U.end() , Uo.begin() , Uo.end() ,
                    std::back_inserter( V ) );
    const Bdd merged = ( T & frame( V , U ) ) | relation( o , V , in );
    if( ( max_nodes > 0 ) && ( merged.nodes() <= max_nodes ) ) {
     T = merged;
     U = V;
     ++count;
     continue;
     }
    close();
    }
   U = Uo;
   T = relation( o , U , in );
   count = 1;
   }
  close();
  }

 /// the sets \f$ H_k \f$ of the sum of the terms t
 Levels levels( const std::vector< SymbolicTerm< cost_type > > & t ) const {
  Levels P = { { cost_type( 0 ) , Bdd::one() } };
  for( const auto & term : t ) {
   std::size_t size = 1;
   for( auto v : term.vars ) {
    size *= v_range[ check( v , 0 ) ];
    if( v_derived[ v ] )
     throw( std::invalid_argument( "SymbolicModel: a term of the heuristic "
                                   "on a derived variable" ) );
    }
   if( term.table.size() != size )
    throw( std::invalid_argument( "SymbolicModel: a term of the heuristic "
                                  "with a table of the wrong size" ) );
   // the assignments of the variables of the term, by value
   Levels B;
   std::vector< std::size_t > x( term.vars.size() , 0 );
   for( std::size_t i = 0 ; i < size ; ++i ) {
    if( term.table[ i ] != infinite_term< cost_type >() ) {
     if constexpr( std::is_signed_v< cost_type > )
      if( term.table[ i ] < 0 )
       throw( std::invalid_argument( "SymbolicModel: a negative term" ) );
     Bdd a = Bdd::one();
     for( std::size_t j = 0 ; j < x.size() ; ++j )
      a &= value( term.vars[ j ] , x[ j ] , false );
     auto & s = B[ term.table[ i ] ];
     s |= a;
     }
    for( std::size_t j = 0 ; j < x.size() ; ++j )
     if( ++x[ j ] < v_range[ term.vars[ j ] ] )
      break;
     else
      x[ j ] = 0;
    }
   Levels Q;
   for( const auto & [ a , Pa ] : P )
    for( const auto & [ b , Bb ] : B ) {
     const Bdd s = Pa & Bb;
     if( ! s.empty() ) {
      auto & q = Q[ a + b ];
      q |= s;
      }
     }
   P = std::move( Q );
   }
  return( P );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_d;                                   ///< the domain
 std::vector< bool > v_derived;                   ///< v is derived
 std::vector< std::uint32_t > v_default;          ///< its default value
 std::vector< int > v_layer;                      ///< its layer
 std::vector< std::vector< std::uint32_t > > v_depends;  ///< what it reads
 std::vector< SymbolicAxiom > v_axioms;           ///< the axioms
 std::vector< Bdd > v_other;                      ///< where it is not the
                                                  ///< default
 std::vector< std::size_t > v_offset;             ///< first bit of each v
 std::vector< std::size_t > v_bits;               ///< bits of each v
 std::vector< std::uint32_t > v_order;            ///< the variables in order
 std::size_t f_bits = 0;                          ///< bits of a state
 std::vector< std::size_t > v_range;              ///< | D_v |
 std::vector< std::uint32_t > v_state_vars;       ///< the even variables
 Bdd f_state_set;                                 ///< their set
 Bdd f_valid;                                     ///< the valid states
 Bdd f_init;                                      ///< the initial state
 Bdd f_goal;                                      ///< the goal states
 std::map< cost_type , std::vector< SymbolicOperator< cost_type > > > v_ops;
                                                  ///< the operators by cost
 std::size_t f_max_nodes = 0;                     ///< the bound of a merge
 std::uint64_t f_abs_images = 0;                  ///< products of the
                                                  ///< abstract searches
 std::uint64_t f_abs_steps = 0;                   ///< their sets expanded
 std::vector< Partition > v_part;                 ///< the partitions
 std::size_t f_operators = 0;                     ///< the operators
 Levels v_forward;                                ///< \f$ H_k \f$ forward
 Levels v_backward;                               ///< \f$ H_k \f$ backward
 bool f_has_forward = false;                      ///< forward terms given
 bool f_has_backward = false;                     ///< backward terms given

 };  // end( class( SymbolicModel ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SymbolicModel */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SymbolicModel.h ------------------------*/
/*--------------------------------------------------------------------------*/
