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
 /// every operator costs 1, as in a breadth-first search
 /** Throws std::invalid_argument if the domain is malformed: a value out
  * of its range, a negative cost, a variable written twice by one
  * operator, a copy into a variable of fewer bits. Needs a BddRuntime. */

 explicit SymbolicModel( const D & d , std::size_t max_nodes = 100000 ,
                         bool unit_costs = false ) : f_d( d ) {
  const std::size_t n = d.num_variables();
  v_offset.resize( n + 1 , 0 );
  v_range.resize( n );
  for( std::size_t v = 0 ; v < n ; ++v ) {
   v_range[ v ] = d.range( std::uint32_t( v ) );
   if( v_range[ v ] == 0 )
    throw( std::invalid_argument( "SymbolicModel: a variable without "
                                  "values" ) );
   const std::size_t b = v_range[ v ] <= 1 ? 0 :
    std::size_t( std::bit_width( v_range[ v ] - 1 ) );
   v_offset[ v + 1 ] = v_offset[ v ] + b;
   }
  std::vector< std::uint32_t > cur;
  for( std::size_t k = 0 ; k < v_offset[ n ] ; ++k )
   cur.push_back( std::uint32_t( 2 * k ) );
  v_state_vars = cur;
  f_state_set = Bdd::variables( cur );

  f_valid = Bdd::one();
  for( std::size_t v = 0 ; v < n ; ++v )
   f_valid &= below( v , v_range[ v ] );

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

  // the operators by cost, in their order
  std::map< cost_type , std::vector< SymbolicOperator< cost_type > > > ops;
  d.for_each_operator( [ & ]( const SymbolicOperator< cost_type > & o ) {
   if constexpr( std::is_signed_v< cost_type > )
    if( o.cost < 0 )
     throw( std::invalid_argument( "SymbolicModel: a negative cost" ) );
   ops[ unit_costs ? cost_type( 1 ) : o.cost ].push_back( o );
   ++f_operators;
   } );
  for( const auto & [ c , list ] : ops )
   merge( c , list , max_nodes );

  if constexpr( SymbolicHeuristicDomain< D > ) {
   std::vector< SymbolicTerm< cost_type > > t;
   d.for_each_heuristic_term( [ & ]( const SymbolicTerm< cost_type > & x ) {
    t.push_back( x );
    } );
   v_forward = levels( t );
   f_has_forward = true;
   }
  else
   v_forward = { { cost_type( 0 ) , Bdd::one() } };
  if constexpr( SymbolicBackwardHeuristicDomain< D > ) {
   std::vector< SymbolicTerm< cost_type > > t;
   d.for_each_backward_heuristic_term(
    [ & ]( const SymbolicTerm< cost_type > & x ) { t.push_back( x ); } );
   v_backward = levels( t );
   f_has_backward = true;
   }
  else
   v_backward = { { cost_type( 0 ) , Bdd::one() } };
  }

/*------------------------------- THE STATES -------------------------------*/
 /// the domain

 [[nodiscard]] const D & domain( void ) const noexcept { return( f_d ); }

 /// the variables of the domain

 [[nodiscard]] std::size_t num_variables( void ) const noexcept {
  return( v_range.size() );
  }

 /// the bits of a state

 [[nodiscard]] std::size_t num_bits( void ) const noexcept {
  return( v_offset.back() );
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
  const std::size_t b = v_offset[ v + 1 ] - v_offset[ v ];
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
   const std::size_t b = v_offset[ v + 1 ] - v_offset[ v ];
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
   for( std::size_t k = v_offset[ v ] ; k < v_offset[ v + 1 ] ; ++k )
    x[ v ] = ( x[ v ] << 1 ) | ( bits[ k ] ? 1 : 0 );
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
 /// backward), checked on the valid states
 /** For every partition of cost c and every k, the image of \f$ H_k \f$
  * (its preimage if backward) must not meet \f$ H_j \f$ for any \f$ j < k
  * - c \f$, i.e., \f$ h( x ) \leq c + h( y ) \f$ along every move from x to
  * y, and the target must lie in \f$ H_0 \f$ (or have no level, a dead
  * end being no target): together these give the admissibility [see
  * docs/theory/symbolic-search.md]. */

 [[nodiscard]] bool consistent( bool backward ) const {
  const Levels & H = levels( backward );
  const Bdd & target = backward ? f_init : f_goal;
  for( const auto & [ k , Hk ] : H )
   if( ( k > 0 ) && ! target.disjoint( Hk ) )
    return( false );
  for( const auto & p : v_part )
   for( const auto & [ k , Hk ] : H ) {
    const Bdd Y = backward ? preimage( Hk & f_valid , p.relation , p.vars )
                           : image( Hk & f_valid , p.relation , p.vars );
    for( const auto & [ j , Hj ] : H ) {
     if( j + p.cost >= k )
      break;
     if( ! Y.disjoint( Hj ) )
      return( false );
     }
    }
  return( true );
  }

/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/

 private:

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
  const std::size_t b = v_offset[ v + 1 ] - v_offset[ v ];
  if( ( b >= 64 ) || ( r >= ( std::size_t( 1 ) << b ) ) )
   return( Bdd::one() );
  Bdd s;
  for( std::size_t d = 0 ; d < r ; ++d )
   s |= value( v , d , false );
  return( s );
  }

 /// \f$ x'_t = x_s \f$, the missing high bits of x_s taken as 0
 Bdd copy( std::size_t t , std::size_t s ) const {
  const std::size_t bt = v_offset[ t + 1 ] - v_offset[ t ];
  const std::size_t bs = v_offset[ s + 1 ] - v_offset[ s ];
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
   for( std::size_t k = v_offset[ v ] ; k < v_offset[ v + 1 ] ; ++k ) {
    w.push_back( std::uint32_t( 2 * k ) );
    w.push_back( std::uint32_t( 2 * k + 1 ) );
    }
  return( Bdd::variables( w ) );
  }

 /// the variables that o reads or writes, sorted, and checked
 std::vector< std::uint32_t > support( const SymbolicOperator< cost_type > &
                                       o ) const {
  std::vector< std::uint32_t > U , W;
  for( const auto & f : o.pre )
   U.push_back( check( f.var , f.value ) );
  for( const auto & e : o.effects ) {
   if( e.copy ) {
    check( e.var , 0 );
    U.push_back( check( e.source , 0 ) );
    }
   else
    check( e.var , e.value );
   U.push_back( e.var );
   W.push_back( e.var );
   }
  std::sort( W.begin() , W.end() );
  if( std::adjacent_find( W.begin() , W.end() ) != W.end() )
   throw( std::invalid_argument( "SymbolicModel: a variable written twice "
                                 "by an operator" ) );
  std::sort( U.begin() , U.end() );
  U.erase( std::unique( U.begin() , U.end() ) , U.end() );
  return( U );
  }

 /// \f$ T_o \f$ on the variables U, which hold its support
 Bdd relation( const SymbolicOperator< cost_type > & o ,
               const std::vector< std::uint32_t > & U ) const {
  Bdd r = Bdd::one();
  for( const auto & f : o.pre )
   r &= value( f.var , f.value , false );
  std::vector< bool > written( num_variables() , false );
  for( const auto & e : o.effects ) {
   r &= e.copy ? copy( e.var , e.source ) : value( e.var , e.value , true );
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

 /// the partitions of the operators of cost c
 void merge( cost_type c ,
             const std::vector< SymbolicOperator< cost_type > > & list ,
             std::size_t max_nodes ) {
  std::vector< std::uint32_t > U;    // the variables of the open partition
  Bdd T;                             // its relation
  std::size_t count = 0;             // its operators
  const auto close = [ & ]( void ) {
   if( count > 0 )
    v_part.push_back( { c , T , vars_of( U ) , count } );
   };
  for( const auto & o : list ) {
   const auto Uo = support( o );
   if( count > 0 ) {
    std::vector< std::uint32_t > V;
    std::set_union( U.begin() , U.end() , Uo.begin() , Uo.end() ,
                    std::back_inserter( V ) );
    const Bdd merged = ( T & frame( V , U ) ) | relation( o , V );
    if( ( max_nodes > 0 ) && ( merged.nodes() <= max_nodes ) ) {
     T = merged;
     U = V;
     ++count;
     continue;
     }
    close();
    }
   U = Uo;
   T = relation( o , U );
   count = 1;
   }
  close();
  }

 /// the sets \f$ H_k \f$ of the sum of the terms t
 Levels levels( const std::vector< SymbolicTerm< cost_type > > & t ) const {
  Levels P = { { cost_type( 0 ) , Bdd::one() } };
  for( const auto & term : t ) {
   std::size_t size = 1;
   for( auto v : term.vars )
    size *= v_range[ check( v , 0 ) ];
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
 std::vector< std::size_t > v_offset;             ///< first bit of each v
 std::vector< std::size_t > v_range;              ///< | D_v |
 std::vector< std::uint32_t > v_state_vars;       ///< the even variables
 Bdd f_state_set;                                 ///< their set
 Bdd f_valid;                                     ///< the valid states
 Bdd f_init;                                      ///< the initial state
 Bdd f_goal;                                      ///< the goal states
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
