/*--------------------------------------------------------------------------*/
/*----------------------- File MultiObjectiveGraph.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template MultiObjectiveGraph, a directed graph
 * whose arcs cost a vector of K components, given explicitly, which
 * implements the MultiObjectiveDomain concept [see MultiObjective.h], for
 * the readers of the road graphs of the 9th DIMACS Implementation
 * Challenge with several costs, and for the ideal-point heuristic.
 *
 * The heuristic of the multi-objective searches on the road graphs is the
 * ideal point of each state, \f$ h_i( u ) = d_i( u , t ) \f$, the cost of a
 * cheapest path from u to the target t by the i-th objective alone,
 * computed by a backward Dijkstra's algorithm for each objective [see
 * Hernández et al., AIJ 2023, section 8.1]: it is consistent, since
 * \f$ d_i( u , t ) \leq c_i( u , v ) + d_i( v , t ) \f$, and it is
 * infinite at the nodes from which t is not reachable, which the searches
 * never generate.
 *
 * The road graphs have two costs, the distances (USA-road-d) and the times
 * (USA-road-t), in two files with the same arcs in the same order; the
 * multi-objective searches add more objectives computed from the graph, as
 * the papers do: the number of arcs, i.e., 1 for each arc [see Maristany de
 * las Casas et al., Networks 2023], a cost of 2 for the arcs whose end
 * nodes have on average 4 neighbours or more and of 1 otherwise [see Ren et
 * al., SoCS 2022], the average number of neighbours of the two end nodes
 * [see Ahmadi et al., ICAPS 2025], and a random integer in [ 1 , 100 ] [see
 * Hernández et al., IJCAI 2023]; the neighbours of a node are the distinct
 * nodes joined to it by an arc in either direction.
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

#ifndef __PBeFS_MultiObjectiveGraph
 #define __PBeFS_MultiObjectiveGraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cstddef>

#include <cstdint>

#include <istream>

#include <random>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "ADS/priority/DaryHeap.h"

#include "PBeFS/astar/ExplicitGraph.h"

#include "PBeFS/astar/MultiObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// an objective computed from a road graph [see the file header]

enum class ExtraCost {
 arcs ,         ///< 1 for each arc
 degree ,       ///< 2 if the end nodes have 4 neighbours on average, else 1
 avg_degree ,   ///< the average neighbours of the end nodes, rounded down
 random         ///< a random integer in [ 1 , 100 ]
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the costs of the cheapest paths to t of n nodes by the objective i, the
/// arcs ( u , v ) of cost c into v given by pred( v , f ), which calls
/// f( u , c ) for each of them; infinite_cost() for the nodes from which t
/// is not reachable

template< class Cost , std::size_t K , class Pred >
[[nodiscard]] std::vector< Cost > distances_to( std::size_t n ,
                                                std::uint32_t t ,
                                                std::size_t i ,
                                                Pred && pred )
{
 struct Entry {
  Cost d;
  std::uint32_t v;
  };
 struct Better {
  bool operator()( const Entry & a , const Entry & b ) const noexcept {
   return( a.d < b.d );
   }
  };
 std::vector< Cost > d( n , infinite_cost< Cost >() );
 ADS::DaryHeap< Entry , Better , 4 > heap;
 d[ t ] = Cost( 0 );
 heap.push( Entry{ Cost( 0 ) , t } );
 while( ! heap.empty() ) {
  const Entry e = heap.pop();
  if( e.d > d[ e.v ] )
   continue;
  pred( e.v , [ & ]( std::uint32_t u , const std::array< Cost , K > & c ) {
   const Cost x = e.d + c[ i ];
   if( x < d[ u ] ) {
    d[ u ] = x;
    heap.push( Entry{ x , u } );
    }
   } );
  }
 return( d );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS MultiObjectiveGraph ------------------------*/
/*--------------------------------------------------------------------------*/
/// a directed graph with K costs on each arc, searched from a source to a
/// target
/** The nodes are \f$ 0 , \ldots , n - 1 \f$; add_arc() adds an arc with its
 * vector of costs, set_heuristic() gives \f$ h \f$ node by node (0 by
 * default) and set_ideal_heuristic() the ideal point of every node. The
 * arcs are kept in two ExplicitAdjacency, forward and reversed, of the set
 * of data structures chosen [see Structures.h]; the parallel arcs are all
 * kept, and so are the loops but in the set ads. */

template< class Cost , std::size_t K >
class MultiObjectiveGraph {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using state_type = std::uint32_t;                  ///< a node
 using cost_type = Cost;                            ///< a component
 static constexpr std::size_t objectives = K;       ///< the objectives
 using Vec = std::array< Cost , K >;                ///< a vector of costs

/*--------------------- CONSTRUCTOR AND INITIALIZATION ---------------------*/

 /// constructor of a graph with n nodes and no arc, searched from s to t

 MultiObjectiveGraph( std::size_t n , state_type s , state_type t )
  : f_source( s ) , f_target( t ) , f_adj( n ) , f_radj( n ) ,
    v_h( n , Vec{} ) {
  if( ( s >= n ) || ( t >= n ) )
   throw( std::invalid_argument( "MultiObjectiveGraph: source or target "
                                 "out of range" ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// adds the arc ( u , v ) of the costs c, all nonnegative

 void add_arc( state_type u , state_type v , const Vec & c ) {
  if( ( u >= size() ) || ( v >= size() ) ||
      std::any_of( c.begin() , c.end() , []( Cost x ) {
       return( x < Cost( 0 ) );
       } ) )
   throw( std::invalid_argument( "MultiObjectiveGraph::add_arc: bad arc" ) );
  f_adj.add( u , v , c );
  f_radj.add( v , u , c );
  ++f_arcs;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the heuristic of the node u

 void set_heuristic( state_type u , const Vec & h ) { v_h.at( u ) = h; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the heuristic of every node to its ideal point for the target [see
 /// the file header]

 void set_ideal_heuristic( void ) {
  for( std::size_t i = 0 ; i < K ; ++i ) {
   const auto d = distances_to< Cost , K >( size() , f_target , i ,
                                            [ & ]( state_type v , auto && f ) {
    f_radj.for_each( v , f );
    } );
   for( std::size_t u = 0 ; u < size() ; ++u )
    v_h[ u ][ i ] = d[ u ];
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the source of the search

 void set_source( state_type s ) {
  if( s >= size() )
   throw( std::invalid_argument( "MultiObjectiveGraph::set_source: bad "
                                 "node" ) );
  f_source = s;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the target of the search, whose heuristic has to be set again

 void set_target( state_type t ) {
  if( t >= size() )
   throw( std::invalid_argument( "MultiObjectiveGraph::set_target: bad "
                                 "node" ) );
  f_target = t;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the target of the search

 [[nodiscard]] state_type target( void ) const noexcept {
  return( f_target );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_source );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( s == f_target );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 [[nodiscard]] const Vec & heuristic( state_type s ) const noexcept {
  return( v_h[ s ] );
  }

 template< class F >
 void for_each_successor( state_type s , const Vec & , F && f ) const {
  f_adj.for_each( s , [ & ]( state_type v , const Vec & c ) {
   f( v , c , v_h[ v ] );
   } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// calls f( u , c ) for each arc ( u , s ) of the costs c

 template< class F >
 void for_each_predecessor( state_type s , F && f ) const {
  f_radj.for_each( s , std::forward< F >( f ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of nodes

 [[nodiscard]] std::size_t size( void ) const noexcept {
  return( f_adj.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of arcs added

 [[nodiscard]] std::size_t arcs( void ) const noexcept { return( f_arcs ); }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 state_type f_source;                               ///< the source
 state_type f_target;                               ///< the target
 ExplicitAdjacency< Vec > f_adj;                    ///< the arcs
 ExplicitAdjacency< Vec > f_radj;                   ///< the arcs reversed
 std::vector< Vec > v_h;                            ///< the heuristic
 std::size_t f_arcs = 0;                            ///< the arcs added

 };  // end( class( MultiObjectiveGraph ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the road graph of the files .gr in, one per objective, with the same arcs
/// in the same order, and with the objectives extra after theirs, searched
/// from source to target (numbered from 1, as in the files)
/** in.size() + extra.size() must be K; the random objective is drawn from
 * std::mt19937_64 seeded with seed, arc by arc in the order of the files.
 * The heuristic is left to 0 [see set_ideal_heuristic()]. Throws
 * std::invalid_argument if the files are not a graph of the format of the
 * challenge [see read_dimacs_graph()], have different nodes or arcs, or if
 * the objectives are not K. */

template< class Cost , std::size_t K >
[[nodiscard]] MultiObjectiveGraph< Cost , K > read_dimacs_objectives(
 const std::vector< std::istream * > & in ,
 const std::vector< ExtraCost > & extra = {} ,
 std::uint32_t source = 1 , std::uint32_t target = 1 ,
 std::uint64_t seed = 1 )
{
 const std::size_t m = in.size();
 if( ( m == 0 ) || ( m + extra.size() != K ) )
  throw( std::invalid_argument( "read_dimacs_objectives: not K "
                                "objectives" ) );
 // the line p of each file
 std::size_t n = 0 , arcs = 0;
 for( std::istream * f : in ) {
  std::string line;
  std::size_t nf = 0 , mf = 0;
  while( std::getline( *f , line ) )
   if( ! line.empty() && ( line[ 0 ] == 'p' ) ) {
    std::istringstream is( line );
    std::string p , sp;
    if( ! ( is >> p >> sp >> nf >> mf ) || ( sp != "sp" ) || ( nf == 0 ) )
     throw( std::invalid_argument( "read_dimacs_objectives: bad line p" ) );
    break;
    }
  if( ( nf == 0 ) || ( n && ( ( nf != n ) || ( mf != arcs ) ) ) )
   throw( std::invalid_argument( "read_dimacs_objectives: no line p, or "
                                 "files of different graphs" ) );
  n = nf;
  arcs = mf;
  }
 if( ( source == 0 ) || ( source > n ) || ( target == 0 ) || ( target > n ) )
  throw( std::invalid_argument( "read_dimacs_objectives: bad source or "
                                "target" ) );

 // the arcs, line by line in all the files at once
 std::vector< std::pair< std::uint32_t , std::uint32_t > > ends;
 std::vector< std::array< Cost , K > > costs;
 ends.reserve( arcs );
 costs.reserve( arcs );
 for( ; ; ) {
  std::array< Cost , K > c{};
  std::uint64_t u0 = 0 , v0 = 0;
  std::size_t got = 0;
  for( std::size_t i = 0 ; i < m ; ++i ) {
   std::string line;
   while( std::getline( *in[ i ] , line ) && ( line.empty() ||
                                               ( line[ 0 ] != 'a' ) ) )
    ;
   if( line.empty() || ( line[ 0 ] != 'a' ) )
    continue;
   std::istringstream is( line );
   char a;
   std::uint64_t u , v;
   long long w;
   if( ! ( is >> a >> u >> v >> w ) || ( u == 0 ) || ( v == 0 ) ||
       ( u > n ) || ( v > n ) || ( w < 0 ) || ( got && ( ( u != u0 ) ||
                                                         ( v != v0 ) ) ) )
    throw( std::invalid_argument( "read_dimacs_objectives: bad arc " +
                                  line ) );
   u0 = u;
   v0 = v;
   c[ i ] = Cost( w );
   ++got;
   }
  if( got == 0 )
   break;
  if( got != m )
   throw( std::invalid_argument( "read_dimacs_objectives: files of "
                                 "different arcs" ) );
  ends.emplace_back( std::uint32_t( u0 - 1 ) , std::uint32_t( v0 - 1 ) );
  costs.push_back( c );
  }
 if( ends.size() < arcs )
  throw( std::invalid_argument( "read_dimacs_objectives: missing arcs" ) );

 // the distinct neighbours of each node, for the objectives of the degrees
 std::vector< std::uint32_t > deg( n , 0 );
 if( std::any_of( extra.begin() , extra.end() , []( ExtraCost e ) {
      return( ( e == ExtraCost::degree ) || ( e == ExtraCost::avg_degree ) );
      } ) ) {
  std::vector< std::pair< std::uint32_t , std::uint32_t > > e;
  e.reserve( 2 * ends.size() );
  for( const auto & [ u , v ] : ends )
   if( u != v ) {
    e.emplace_back( u , v );
    e.emplace_back( v , u );
    }
  std::sort( e.begin() , e.end() );
  e.erase( std::unique( e.begin() , e.end() ) , e.end() );
  for( const auto & x : e )
   ++deg[ x.first ];
  }

 std::mt19937_64 gen( seed );
 MultiObjectiveGraph< Cost , K > g( n , source - 1 , target - 1 );
 for( std::size_t a = 0 ; a < ends.size() ; ++a ) {
  const auto [ u , v ] = ends[ a ];
  for( std::size_t j = 0 ; j < extra.size() ; ++j ) {
   Cost & x = costs[ a ][ m + j ];
   switch( extra[ j ] ) {
    case( ExtraCost::arcs ):
     x = Cost( 1 );
     break;
    case( ExtraCost::degree ):
     x = Cost( deg[ u ] + deg[ v ] >= 8 ? 2 : 1 );
     break;
    case( ExtraCost::avg_degree ):
     x = Cost( ( deg[ u ] + deg[ v ] ) / 2 );
     break;
    default:
     x = Cost( 1 + random_below( gen , 100 ) );
    }
   }
  g.add_arc( u , v , costs[ a ] );
  }
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the ExtraCost of a name (arcs, degree, avg-degree, random); throws
/// std::invalid_argument for another one

[[nodiscard]] inline ExtraCost extra_cost( const std::string & name )
{
 if( name == "arcs" )
  return( ExtraCost::arcs );
 if( name == "degree" )
  return( ExtraCost::degree );
 if( name == "avg-degree" )
  return( ExtraCost::avg_degree );
 if( name == "random" )
  return( ExtraCost::random );
 throw( std::invalid_argument( "extra_cost: unknown " + name ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_MultiObjectiveGraph */

/*--------------------------------------------------------------------------*/
/*--------------------- End File MultiObjectiveGraph.h ---------------------*/
/*--------------------------------------------------------------------------*/
