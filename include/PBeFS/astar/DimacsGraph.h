/*--------------------------------------------------------------------------*/
/*--------------------------- File DimacsGraph.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the road graphs of the 9th DIMACS Implementation
 * Challenge on shortest paths (Demetrescu, Goldberg and Johnson, 2009):
 * the readers of their files of arcs (.gr) and of coordinates (.co) into an
 * ExplicitGraph [see ExplicitGraph.h], and the heuristic of the
 * coordinates.
 *
 * A file .gr has the line "p sp n m" and m lines "a u v w", the arc
 * ( u , v ) of the integer cost w, the nodes numbered from 1; a file .co
 * has the line "p aux sp co n" and n lines "v v x y", the longitude and the
 * latitude of v in millionths of a degree; the lines "c" are comments. Both
 * directions of a road are arcs of their own. The heuristic of the node u
 * for the target t is \f$ h( u ) = \lfloor \kappa \, \delta( u , t )
 * \rfloor \f$, \f$ \delta \f$ being the great-circle distance in meters and
 * \f$ \kappa = \min_{( u , v )} w( u , v ) / \delta( u , v ) \f$ over the
 * arcs of positive length: since \f$ \delta \f$ is a metric, \f$ h( u )
 * \leq \lfloor w( u , v ) + \kappa \, \delta( v , t ) \rfloor = w( u , v ) +
 * h( v ) \f$, i.e., h is consistent (and 0 at the target) whatever the
 * units of the costs, meters for the graphs USA-road-d and tenths of a
 * second for USA-road-t, as long as the data are; \f$ \kappa \f$ is the
 * largest factor for which this holds.
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

#ifndef __PBeFS_DimacsGraph
 #define __PBeFS_DimacsGraph
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <istream>

#include <limits>

#include <sstream>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/astar/ExplicitGraph.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the coordinates of the nodes of a road graph, in millionths of a degree

//! [domains-sig-dimacs]
struct DimacsCoordinates {
 std::vector< std::int32_t > lon;     ///< the longitude of each node
 std::vector< std::int32_t > lat;     ///< the latitude of each node
 };
//! [domains-sig-dimacs]

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the graph of the file .gr in, searched from source to target (both
/// numbered from 1, as in the file)
/** Throws std::invalid_argument if the file has no line "p sp", an arc out
 * of range or of negative cost, or fewer arcs than it says. */

template< class Cost = std::int64_t >
[[nodiscard]] ExplicitGraph< Cost > read_dimacs_graph(
 std::istream & in , std::uint32_t source = 1 , std::uint32_t target = 1 )
{
 std::string line;
 std::size_t n = 0 , m = 0 , read = 0;
 while( std::getline( in , line ) )
  if( ! line.empty() && ( line[ 0 ] == 'p' ) ) {
   std::istringstream is( line );
   std::string p , sp;
   if( ! ( is >> p >> sp >> n >> m ) || ( sp != "sp" ) || ( n == 0 ) )
    throw( std::invalid_argument( "read_dimacs_graph: bad line p" ) );
   break;
   }
 if( ( n == 0 ) || ( source == 0 ) || ( source > n ) || ( target == 0 ) ||
     ( target > n ) )
  throw( std::invalid_argument( "read_dimacs_graph: no line p, or bad "
                                "source or target" ) );
 ExplicitGraph< Cost > g( n , source - 1 , target - 1 );
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] != 'a' ) )
   continue;
  std::istringstream is( line );
  char a;
  std::uint64_t u , v;
  long long w;
  if( ! ( is >> a >> u >> v >> w ) || ( u == 0 ) || ( v == 0 ) ||
      ( u > n ) || ( v > n ) || ( w < 0 ) )
   throw( std::invalid_argument( "read_dimacs_graph: bad arc " + line ) );
  g.add_edge( std::uint32_t( u - 1 ) , std::uint32_t( v - 1 ) , Cost( w ) ,
              true );
  ++read;
  }
 if( read < m )
  throw( std::invalid_argument( "read_dimacs_graph: missing arcs" ) );
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the coordinates of the file .co in, of n nodes
/** Throws std::invalid_argument if the file has no line "p aux sp co" for
 * n nodes, or a node out of range. */

[[nodiscard]] inline DimacsCoordinates read_dimacs_coordinates(
 std::istream & in , std::size_t n )
{
 std::string line;
 std::size_t k = 0;
 while( std::getline( in , line ) )
  if( ! line.empty() && ( line[ 0 ] == 'p' ) ) {
   std::istringstream is( line );
   std::string p , aux , sp , co;
   if( ! ( is >> p >> aux >> sp >> co >> k ) || ( co != "co" ) )
    throw( std::invalid_argument( "read_dimacs_coordinates: bad line p" ) );
   break;
   }
 if( k != n )
  throw( std::invalid_argument( "read_dimacs_coordinates: not the nodes "
                                "of the graph" ) );
 DimacsCoordinates c;
 c.lon.assign( n , 0 );
 c.lat.assign( n , 0 );
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] != 'v' ) )
   continue;
  std::istringstream is( line );
  char a;
  std::uint64_t v;
  long long x , y;
  if( ! ( is >> a >> v >> x >> y ) || ( v == 0 ) || ( v > n ) )
   throw( std::invalid_argument( "read_dimacs_coordinates: bad line " +
                                 line ) );
  c.lon[ v - 1 ] = std::int32_t( x );
  c.lat[ v - 1 ] = std::int32_t( y );
  }
 return( c );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the great-circle distance in meters between the nodes u and v, by the
/// formula of the haversine on the sphere of radius 6371008.8 meters

[[nodiscard]] inline double dimacs_distance( const DimacsCoordinates & c ,
                                             std::size_t u , std::size_t v )
 noexcept
{
 constexpr double rad = 3.14159265358979323846 / 180e6;
 const double p1 = c.lat[ u ] * rad , p2 = c.lat[ v ] * rad;
 const double dp = p2 - p1 , dl = double( c.lon[ v ] - c.lon[ u ] ) * rad;
 const double a = std::sin( dp / 2 ) * std::sin( dp / 2 ) +
                  std::cos( p1 ) * std::cos( p2 ) *
                  std::sin( dl / 2 ) * std::sin( dl / 2 );
 return( 2 * 6371008.8 * std::asin( std::min( 1.0 , std::sqrt( a ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the factor \f$ \kappa \f$ of the heuristic of the coordinates [see the
/// file header], slightly lowered against the rounding of the distances

template< class Cost >
[[nodiscard]] double dimacs_scale( const ExplicitGraph< Cost > & g ,
                                   const DimacsCoordinates & c )
{
 double k = std::numeric_limits< double >::infinity();
 for( std::uint32_t u = 0 ; u < g.size() ; ++u )
  g.for_each_successor( u , Cost( 0 ) , [ & ]( std::uint32_t v , Cost w ,
                                               Cost ) {
   const double d = dimacs_distance( c , u , v );
   if( d > 0 )
    k = std::min( k , double( w ) / d );
   } );
 if( ! std::isfinite( k ) )
  k = 0;
 return( k * ( 1 - 1e-9 ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// sets the heuristic of the coordinates for the target of g, with the
/// factor kappa (dimacs_scale( g , c ) if negative); returns kappa

template< class Cost >
double set_dimacs_heuristic( ExplicitGraph< Cost > & g ,
                             const DimacsCoordinates & c ,
                             double kappa = -1 )
{
 if( ( c.lon.size() != g.size() ) || ( c.lat.size() != g.size() ) )
  throw( std::invalid_argument( "set_dimacs_heuristic: coordinates and "
                                "graph differ" ) );
 if( kappa < 0 )
  kappa = dimacs_scale( g , c );
 const std::size_t t = g.target();
 for( std::uint32_t u = 0 ; u < g.size() ; ++u )
  g.set_heuristic( u , Cost( std::floor( kappa * dimacs_distance( c , u ,
                                                                  t ) ) ) );
 return( kappa );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_DimacsGraph */

/*--------------------------------------------------------------------------*/
/*------------------------ End File DimacsGraph.h --------------------------*/
/*--------------------------------------------------------------------------*/
