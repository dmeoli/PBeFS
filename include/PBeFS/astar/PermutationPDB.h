/*--------------------------------------------------------------------------*/
/*------------------------ File PermutationPDB.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the pattern databases of the permutation puzzles: the
 * concept PermutationDomain, the abstract space PermutationSpace of a
 * pattern, the class template PermutationPDB of its table, built by
 * abstract_distances() [see AbstractionPDB.h], PermutationHeuristic, the
 * maximum or the sum of several of them, and PDBDomain, any Domain with
 * such a heuristic in place of (or besides) its own.
 *
 * A permutation puzzle has groups of positions, each holding as many
 * objects as positions; an object of a group with modulus m also has an
 * orientation in \f$ [ 0 , m ) \f$. A move is, for each group, a
 * permutation of the positions and a turn per position: the object that
 * comes to the position i is the one that was at \f$ \mathrm{from}[ i ]
 * \f$, and its orientation grows by \f$ \mathrm{turn}[ i ] \f$ modulo m.
 * The pancakes, the burned pancakes (m = 2), TopSpin and the corners
 * (m = 3) and edges (m = 2) of Rubik's cube are such puzzles; the
 * sliding-tile puzzles are not, since a move depends on where the blank
 * is, and keep their own tables [see PatternDatabase.h].
 *
 * A pattern is a set of k objects of a group, all the others being "don't
 * care" as in Culberson and Schaeffer (1998): its abstract state is the
 * sequence of the positions of the k objects, a k-permutation of the n
 * positions, with their k orientations, and its rank is
 *
 * \f[
 *   \mathrm{rank}( p , o ) = \mathrm{rank}_R( p ) \, m^{k'} +
 *   \sum_{j < k'} o_j m^j ,
 * \f]
 *
 * \f$ \mathrm{rank}_R \f$ being a Ranking [see Permutation.h] and
 * \f$ k' = k \f$, or \f$ k - 1 \f$ when the pattern holds the whole group
 * and the orientations of a group sum to 0 modulo m in every state, the
 * last one being then implied (the \f$ 8! \cdot 3^7 \f$ states of the
 * corners of Korf, 1997). A move is an edge of the abstract space from the
 * image of a state to the image of its successor, of the cost of the move,
 * or a loop when it moves no object of the pattern; the abstract
 * predecessor of \f$ ( p , o ) \f$ by a move has the object j at
 * \f$ \mathrm{from}[ p_j ] \f$ with the orientation
 * \f$ o_j - \mathrm{turn}[ p_j ] \f$. The distance of the image of a state
 * is thus at most the cost of any of its paths, and the maximum of the
 * tables of several patterns is admissible and consistent; their sum is
 * admissible when no move moves objects of two patterns, which holds for
 * no puzzle of this header but holds for the Towers of Hanoi [see
 * TowersOfHanoi.h], whose moves move one disk.
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

#ifndef __PBeFS_PermutationPDB
 #define __PBeFS_PermutationPDB
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <concepts>

#include <cstddef>

#include <cstdint>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/AbstractionPDB.h"

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Permutation.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// a group of the positions of a permutation puzzle

struct PermutationGroup {
 unsigned size = 0;         ///< the positions, and the objects
 unsigned modulus = 1;      ///< the orientations of an object (1: none)
 bool zero_sum = false;     ///< the orientations sum to 0 modulo modulus
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a move of a permutation puzzle, by group

struct PermutationMove {
 /// from[ g ][ i ]: the position whose object comes to i
 std::vector< std::vector< std::uint8_t > > from;
 /// turn[ g ][ i ]: what its orientation grows by
 std::vector< std::vector< std::uint8_t > > turn;
 unsigned cost = 1;         ///< the cost, 1 to 254
 };

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONCEPTS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail {

 /// a callable accepting a permutation, used only to state the concepts
 struct permutation_sink {
  void operator()( const std::uint8_t * , const std::uint8_t * ) const
   noexcept {}
  };

 }  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// a Domain that is a permutation puzzle
/** permutation_groups() and permutation_moves() describe the puzzle [see
 * the file header], every move of the Domain being one of the moves;
 * permutation_of( s , g , obj , ori ) writes in obj[ i ] and ori[ i ] the
 * object at the position i of the group g in the state s and its
 * orientation, and for_each_goal_permutation( g , f ) calls f( obj , ori )
 * for the group g of each goal (of several, if the goals are more). */

template< class D >
concept PermutationDomain = Domain< D > &&
 requires( const D & d , const typename D::state_type & s , unsigned g ,
           std::uint8_t * o , std::uint8_t * r ) {
  { d.permutation_groups() }
   -> std::convertible_to< std::vector< PermutationGroup > >;
  { d.permutation_moves() }
   -> std::convertible_to< std::vector< PermutationMove > >;
  d.permutation_of( s , g , o , r );
  d.for_each_goal_permutation( g , detail::permutation_sink{} );
  };

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS PermutationSpace -------------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract space of a pattern of a permutation puzzle
/** An AbstractSpace [see AbstractionPDB.h]: the states are the ranks of
 * the positions and orientations of the objects of the pattern [see the
 * file header], the goals the images of the goals of the puzzle, and the
 * edges the moves that move an object of the pattern. */

template< Ranking R = LexRanking >
class PermutationSpace {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the space of the objects of the group g of a puzzle, with its goals
 /** goals holds, for each goal, the n objects and the n orientations of
  * the group g, by position. Throws std::invalid_argument if g, a move,
  * the objects or the goals do not fit the groups, or if the number of
  * the abstract states does not fit 64 bits. */

 PermutationSpace( const std::vector< PermutationGroup > & groups ,
                   const std::vector< PermutationMove > & moves , unsigned g ,
                   const std::vector< std::uint8_t > & objects ,
                   const std::vector< std::pair< std::vector< std::uint8_t > ,
                                                 std::vector< std::uint8_t >
                                                 > > & goals )
  : f_group( check_group( groups , g ) ) , f_objects( objects ) ,
    f_rank( f_group.size , unsigned( objects.size() ) ) {
  const unsigned n = f_group.size , k = unsigned( objects.size() );
  std::vector< bool > seen( n , false );
  for( std::uint8_t x : objects ) {
   if( ( x >= n ) || seen[ x ] )
    throw( std::invalid_argument( "PermutationSpace: bad objects" ) );
   seen[ x ] = true;
   }
  if( k == 0 )
   throw( std::invalid_argument( "PermutationSpace: no object" ) );
  f_digits = ( ( k == n ) && f_group.zero_sum ) ? k - 1 : k;
  f_size = f_rank.size();
  for( unsigned j = 0 ; j < f_digits ; ++j ) {
   if( f_size > ~std::uint64_t( 0 ) / f_group.modulus )
    throw( std::invalid_argument( "PermutationSpace: more than 2^64 "
                                  "states" ) );
   f_size *= f_group.modulus;
   }
  for( const PermutationMove & m : moves ) {
   if( ( m.from.size() <= g ) || ( m.from[ g ].size() != n ) ||
       ( ( m.turn.size() > g ) && ! m.turn[ g ].empty() &&
         ( m.turn[ g ].size() != n ) ) || ( m.cost == 0 ) ||
       ( m.cost > 254 ) )
    throw( std::invalid_argument( "PermutationSpace: bad move" ) );
   Step st;
   st.cost = m.cost;
   st.from = m.from[ g ];
   st.turn = ( m.turn.size() > g ) && ! m.turn[ g ].empty()
             ? m.turn[ g ] : std::vector< std::uint8_t >( n , 0 );
   v_steps.push_back( std::move( st ) );
   }
  for( const auto & [ obj , ori ] : goals ) {
   if( ( obj.size() != n ) || ( ori.size() != n ) )
    throw( std::invalid_argument( "PermutationSpace: bad goal" ) );
   std::vector< std::uint8_t > loc( n ) , o( n );
   for( unsigned i = 0 ; i < n ; ++i ) {
    loc[ obj[ i ] ] = std::uint8_t( i );
    o[ obj[ i ] ] = ori[ i ];
    }
   v_goals.push_back( rank_objects( loc.data() , o.data() ) );
   }
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the number of the abstract states

 [[nodiscard]] std::uint64_t size( void ) const noexcept { return( f_size ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the objects of the pattern

 [[nodiscard]] const std::vector< std::uint8_t > & objects( void ) const
  noexcept {
  return( f_objects );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the group of the pattern

 [[nodiscard]] const PermutationGroup & group( void ) const noexcept {
  return( f_group );
  }

/*-------------------------------- RANKING ---------------------------------*/
 /// the rank of the abstract state whose object j is at p[ j ] with the
 /// orientation o[ j ], j indexing the objects of the pattern

 [[nodiscard]] std::uint64_t rank( const std::uint8_t * p ,
                                   const std::uint8_t * o ) const noexcept {
  std::uint64_t r = f_rank.rank( p );
  for( unsigned j = f_digits ; j-- > 0 ; )
   r = r * f_group.modulus + o[ j ];
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the rank of the image of a state, loc[ x ] and ori[ x ] being the
 /// position and the orientation of each object x of the group

 [[nodiscard]] std::uint64_t rank_objects( const std::uint8_t * loc ,
                                           const std::uint8_t * ori ) const
  noexcept {
  std::uint8_t p[ 256 ] , o[ 256 ];
  for( std::size_t j = 0 ; j < f_objects.size() ; ++j ) {
   p[ j ] = loc[ f_objects[ j ] ];
   o[ j ] = ori[ f_objects[ j ] ];
   }
  return( rank( p , o ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes in p and o the abstract state of rank r

 void unrank( std::uint64_t r , std::uint8_t * p , std::uint8_t * o ) const
  noexcept {
  unsigned sum = 0;
  for( unsigned j = 0 ; j < f_digits ; ++j ) {
   o[ j ] = std::uint8_t( r % f_group.modulus );
   sum += o[ j ];
   r /= f_group.modulus;
   }
  if( f_digits < f_objects.size() )
   o[ f_digits ] = std::uint8_t( ( f_group.modulus - sum % f_group.modulus )
                                 % f_group.modulus );
  f_rank.unrank( r , p );
  }

/*---------------------------- THE ABSTRACT SPACE --------------------------*/
 /// calls f( r ) for each abstract goal

 template< class F >
 void for_each_goal( F && f ) const {
  for( std::uint64_t g : v_goals )
   f( g );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( r' , c ) for each predecessor r' of r by a move of cost c
 /// that moves an object of the pattern

 template< class F >
 void for_each_predecessor( std::uint64_t r , F && f ) const {
  const std::size_t k = f_objects.size();
  const unsigned m = f_group.modulus;
  std::uint8_t p[ 256 ] , o[ 256 ] , q[ 256 ] , u[ 256 ];
  unrank( r , p , o );
  for( const Step & st : v_steps ) {
   bool moved = false;
   for( std::size_t j = 0 ; j < k ; ++j ) {
    q[ j ] = st.from[ p[ j ] ];
    u[ j ] = std::uint8_t( ( o[ j ] + m - st.turn[ p[ j ] ] % m ) % m );
    moved = moved || ( q[ j ] != p[ j ] ) || ( u[ j ] != o[ j ] );
    }
   if( moved )
    f( rank( q , u ) , st.cost );
   }
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a move restricted to the group
 struct Step {
  std::vector< std::uint8_t > from;
  std::vector< std::uint8_t > turn;
  unsigned cost = 1;
  };

 // the group g, checked
 [[nodiscard]] static PermutationGroup check_group(
  const std::vector< PermutationGroup > & groups , unsigned g ) {
  if( ( g >= groups.size() ) || ( groups[ g ].size == 0 ) ||
      ( groups[ g ].size > 256 ) || ( groups[ g ].modulus == 0 ) ||
      ( groups[ g ].modulus > 256 ) )
   throw( std::invalid_argument( "PermutationSpace: bad group" ) );
  return( groups[ g ] );
  }

 PermutationGroup f_group;                      ///< the group
 std::vector< std::uint8_t > f_objects;         ///< the pattern
 R f_rank;                                      ///< the ranking
 unsigned f_digits = 0;                         ///< orientations ranked
 std::uint64_t f_size = 0;                      ///< the abstract states
 std::vector< Step > v_steps;                   ///< the moves
 std::vector< std::uint64_t > v_goals;          ///< the abstract goals

 };  // end( class( PermutationSpace ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS PermutationPDB --------------------------*/
/*--------------------------------------------------------------------------*/
/// the pattern database of a pattern of a permutation puzzle
/** The table of a PermutationSpace, B bits per entry, built on threads
 * workers of Executor and, with a cache directory, saved there under a
 * name and a key that depend on the puzzle, the goals, the pattern, the
 * ranking and B, and loaded back the next time. */

template< Ranking R = LexRanking , unsigned B = 8 >
class PermutationPDB {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the database of the objects of the group g of the Domain d

 template< PermutationDomain D ,
           ADS::WorkerExecutor Executor = ADS::ThreadPool >
 PermutationPDB( const D & d , unsigned g ,
                 const std::vector< std::uint8_t > & objects ,
                 unsigned threads = 0 , const std::string & cache = "" ,
                 Executor * = nullptr )
  : f_group( g ) , f_space( d.permutation_groups() , d.permutation_moves() ,
                            g , objects , goals_of( d , g ) ) {
  const std::uint64_t key = key_of( d , g , objects );
  char hex[ 17 ];
  std::snprintf( hex , sizeof( hex ) , "%016llx" ,
                 static_cast< unsigned long long >( key ) );
  f_table = AbstractionTable< B >::cached(
   cache , std::string( "perm-" ) + R::name + "-" + std::to_string( B ) +
   "-" + hex + ".pdb" , key , f_space.size() , [ & ] {
    return( abstract_distances< Executor >( f_space , threads ) );
    } );
  }

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the image of a state whose objects x of the group are at
 /// loc[ x ] with the orientations ori[ x ]

 [[nodiscard]] unsigned value( const std::uint8_t * loc ,
                               const std::uint8_t * ori ) const noexcept {
  return( f_table[ f_space.rank_objects( loc , ori ) ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the value of the image of the state s of the Domain d

 template< PermutationDomain D >
 [[nodiscard]] unsigned value( const D & d ,
                               const typename D::state_type & s ) const {
  std::uint8_t obj[ 256 ] , ori[ 256 ] , loc[ 256 ] , o[ 256 ];
  d.permutation_of( s , f_group , obj , ori );
  for( unsigned i = 0 ; i < f_space.group().size ; ++i ) {
   loc[ obj[ i ] ] = std::uint8_t( i );
   o[ obj[ i ] ] = ori[ i ];
   }
  return( value( loc , o ) );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the group of the pattern

 [[nodiscard]] unsigned group( void ) const noexcept { return( f_group ); }

 /// the abstract space

 [[nodiscard]] const PermutationSpace< R > & space( void ) const noexcept {
  return( f_space );
  }

 /// the table

 [[nodiscard]] const AbstractionTable< B > & table( void ) const noexcept {
  return( f_table );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // the goals of the group g of d
 template< PermutationDomain D >
 [[nodiscard]] static std::vector< std::pair< std::vector< std::uint8_t > ,
                                              std::vector< std::uint8_t > > >
  goals_of( const D & d , unsigned g ) {
  const auto groups = d.permutation_groups();
  if( g >= groups.size() )
   throw( std::invalid_argument( "PermutationPDB: bad group" ) );
  const unsigned n = groups[ g ].size;
  std::vector< std::pair< std::vector< std::uint8_t > ,
                          std::vector< std::uint8_t > > > v;
  d.for_each_goal_permutation( g , [ & ]( const std::uint8_t * o ,
                                          const std::uint8_t * r ) {
   v.emplace_back( std::vector< std::uint8_t >( o , o + n ) ,
                   std::vector< std::uint8_t >( r , r + n ) );
   } );
  return( v );
  }

 // the key of the table: a checksum of the puzzle, the goals, the pattern,
 // the ranking and B
 template< PermutationDomain D >
 [[nodiscard]] static std::uint64_t key_of( const D & d , unsigned g ,
                                            const std::vector< std::uint8_t >
                                            & objects ) {
  std::vector< std::uint8_t > b;
  const auto add = [ & ]( std::uint64_t x ) {
   for( int i = 0 ; i < 8 ; ++i )
    b.push_back( std::uint8_t( x >> ( 8 * i ) ) );
   };
  for( const PermutationGroup & gr : d.permutation_groups() ) {
   add( gr.size );
   add( gr.modulus );
   add( gr.zero_sum );
   }
  for( const PermutationMove & m : d.permutation_moves() ) {
   add( m.cost );
   if( m.from.size() > g )
    b.insert( b.end() , m.from[ g ].begin() , m.from[ g ].end() );
   if( m.turn.size() > g )
    b.insert( b.end() , m.turn[ g ].begin() , m.turn[ g ].end() );
   }
  for( const auto & [ o , r ] : goals_of( d , g ) ) {
   b.insert( b.end() , o.begin() , o.end() );
   b.insert( b.end() , r.begin() , r.end() );
   }
  add( g );
  add( objects.size() );
  b.insert( b.end() , objects.begin() , objects.end() );
  for( const char * c = R::name ; *c ; ++c )
   b.push_back( std::uint8_t( *c ) );
  add( B );
  return( AbstractionTable< B >::checksum( b.data() , b.size() ) );
  }

 unsigned f_group;                              ///< the group
 PermutationSpace< R > f_space;                 ///< the abstract space
 AbstractionTable< B > f_table;                 ///< the table

 };  // end( class( PermutationPDB ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS PermutationHeuristic -----------------------*/
/*--------------------------------------------------------------------------*/
/// the maximum, or the sum, of pattern databases of a permutation puzzle
/** The state is read once per group, and each table is looked up. The sum
 * is admissible only if no move moves objects of two patterns [see the
 * file header]; the constructor does not check it. */

template< PermutationDomain D , Ranking R = LexRanking , unsigned B = 8 >
class PermutationHeuristic {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using PDB = PermutationPDB< R , B >;              ///< a table

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the heuristic of the tables pdbs of the Domain d (which must outlive
 /// it), their maximum or, if additive, their sum

 PermutationHeuristic( const D & d , std::vector< PDB > pdbs ,
                       bool additive = false )
  : f_d( d ) , v_pdbs( std::move( pdbs ) ) , f_additive( additive ) ,
    f_groups( unsigned( d.permutation_groups().size() ) ) {}

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the value of the state s

 [[nodiscard]] unsigned operator()( const typename D::state_type & s ) const
 {
  std::uint8_t obj[ 256 ] , ori[ 256 ] , loc[ 256 ] , o[ 256 ];
  unsigned v = 0;
  for( unsigned g = 0 ; g < f_groups ; ++g ) {
   bool read = false;
   for( const PDB & p : v_pdbs ) {
    if( p.group() != g )
     continue;
    if( ! read ) {
     f_d.permutation_of( s , g , obj , ori );
     for( unsigned i = 0 ; i < p.space().group().size ; ++i ) {
      loc[ obj[ i ] ] = std::uint8_t( i );
      o[ obj[ i ] ] = ori[ i ];
      }
     read = true;
     }
    const unsigned x = p.value( loc , o );
    v = f_additive ? v + x : std::max( v , x );
    }
   }
  return( v );
  }

/*---------------------- METHODS FOR READING THE DATA ----------------------*/
 /// the tables

 [[nodiscard]] const std::vector< PDB > & pdbs( void ) const noexcept {
  return( v_pdbs );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const D & f_d;                                 ///< the Domain
 std::vector< PDB > v_pdbs;                     ///< the tables
 bool f_additive;                               ///< sum, not maximum
 unsigned f_groups;                             ///< the groups

 };  // end( class( PermutationHeuristic ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS PDBDomain -----------------------------*/
/*--------------------------------------------------------------------------*/
/// a Domain whose heuristic is the maximum of its own and of H
/** Everything is that of D, but the heuristic, \f$ \max( h_D( s ) , H( s )
 * ) \f$ for a callable H of the states (such as a PermutationHeuristic),
 * which is admissible and consistent if both are; the heuristic of each
 * successor is computed anew, that of D included. */

template< Domain D , class H >
class PDBDomain {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = state_t< D >;                 ///< the state type
 using cost_type = cost_t< D >;                   ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain and the heuristic (which must outlive
 /// the PDBDomain)

 PDBDomain( const D & d , const H & h ) : f_d( d ) , f_h( h ) {}

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const {
  return( f_d.initial_state() );
  }

 [[nodiscard]] bool is_goal( const state_type & s ) const {
  return( f_d.is_goal( s ) );
  }

 [[nodiscard]] std::uint64_t hash( const state_type & s ) const {
  return( f_d.hash( s ) );
  }

 [[nodiscard]] cost_type heuristic( const state_type & s ) const {
  return( std::max( cost_type( f_d.heuristic( s ) ) ,
                    cost_type( f_h( s ) ) ) );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type , F && f ) const
 {
  f_d.for_each_successor( s , f_d.heuristic( s ) ,
                          [ & ]( const state_type & t , cost_type c ,
                                 cost_type ht ) {
   f( t , c , std::max( ht , cost_type( f_h( t ) ) ) );
   } );
  }

 template< class F >
 void for_each_successor( const state_type & s , cost_type ,
                          std::uint8_t op , F && f ) const
  requires OperatorDomain< D > {
  f_d.for_each_successor( s , f_d.heuristic( s ) , op ,
                          [ & ]( const state_type & t , cost_type c ,
                                 cost_type ht , std::uint8_t o ) {
   f( t , c , std::max( ht , cost_type( f_h( t ) ) ) , o );
   } );
  }

 template< class F >
 void for_each_tree_successor( const state_type & s , cost_type ,
                               std::uint8_t op , F && f ) const
  requires TreeOperatorDomain< D > {
  f_d.for_each_tree_successor( s , f_d.heuristic( s ) , op ,
                               [ & ]( const state_type & t , cost_type c ,
                                      cost_type ht , std::uint8_t o ) {
   f( t , c , std::max( ht , cost_type( f_h( t ) ) ) , o );
   } );
  }

/*---------------------- FEATURES FOR ZOBRIST HASHING ----------------------*/

 [[nodiscard]] std::size_t num_features( void ) const
  requires FeatureDomain< D > {
  return( f_d.num_features() );
  }

 template< class F >
 void for_each_feature( const state_type & s , F && f ) const
  requires FeatureDomain< D > {
  f_d.for_each_feature( s , std::forward< F >( f ) );
  }

 template< class F >
 void for_each_feature_change( const state_type & s , const state_type & t ,
                               F && f ) const
  requires IncrementalFeatureDomain< D > {
  f_d.for_each_feature_change( s , t , std::forward< F >( f ) );
  }

/*------------------------------ OTHER HASHES ------------------------------*/

 [[nodiscard]] std::uint64_t perfect_hash( const state_type & s ) const
  requires PerfectHashDomain< D > {
  return( f_d.perfect_hash( s ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the Domain

 [[nodiscard]] const D & domain( void ) const noexcept { return( f_d ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 const D & f_d;                                 ///< the Domain
 const H & f_h;                                 ///< the heuristic

 };  // end( class( PDBDomain ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_PermutationPDB */

/*--------------------------------------------------------------------------*/
/*----------------------- End File PermutationPDB.h ------------------------*/
/*--------------------------------------------------------------------------*/
