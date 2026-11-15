/*--------------------------------------------------------------------------*/
/*---------------------------- File RobotNav2D.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class NavMap, the occupancy of a Moving AI map scaled
 * up, and for the class RobotNav2D, which implements the EdgeDomain,
 * HeterogeneousEdgeDomain, LazyEdgeDomain and PairwiseDomain concepts [see
 * Domain.h] for the navigation of a square robot on a NavMap by long
 * straight motions, the domain on which the edge-based parallel searches
 * were compared (section "2D Grid World" of Mukherjee and Likhachev, SoCS
 * 2023; section "Evaluation" of Yang, Mukherjee and Likhachev, SoCS 2023).
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

#ifndef __PBeFS_RobotNav2D
 #define __PBeFS_RobotNav2D
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <array>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <istream>

#include <random>

#include <stdexcept>

#include <string>

#include <vector>

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS NavMap ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a rectangular map of free and blocked cells
/** The cell \f$ ( x , y ) \f$, \f$ 0 \leq x < W \f$, \f$ 0 \leq y < H \f$,
 * is free or blocked; a map read from a file of the benchmarks of
 * Sturtevant (2012) is scaled by an integer factor \f$ k \f$, each cell
 * becoming \f$ k \times k \f$ cells. */

class NavMap {

 public:

 /// constructor of a map of width by height free cells

 NavMap( int width = 0 , int height = 0 )
  : f_width( width ) , f_height( height ) ,
    v_blocked( std::size_t( width ) * std::size_t( height ) , 0 ) {}

 /// reads a map in the format of Sturtevant (2012), scaled by scale; the
 /// cells '.', 'G' and 'S' are free, and 'T' as well if trees_free, as in
 /// the experiments of the edge-based parallel searches

 [[nodiscard]] static NavMap read_movingai( std::istream & in ,
                                            int scale = 1 ,
                                            bool trees_free = true ) {
  std::string word;
  int width = -1 , height = -1;
  while( ( in >> word ) && ( word != "map" ) ) {
   if( word == "height" )
    in >> height;
   else
    if( word == "width" )
     in >> width;
   }
  if( ( width <= 0 ) || ( height <= 0 ) || ( scale <= 0 ) )
   throw( std::invalid_argument( "NavMap::read_movingai: bad header" ) );
  NavMap m( width * scale , height * scale );
  for( int y = 0 ; y < height ; ++y ) {
   std::string row;
   if( ! ( in >> row ) || ( int( row.size() ) != width ) )
    throw( std::invalid_argument( "NavMap::read_movingai: bad row" ) );
   for( int x = 0 ; x < width ; ++x ) {
    const char c = row[ x ];
    const bool free = ( c == '.' ) || ( c == 'G' ) || ( c == 'S' ) ||
                      ( trees_free && ( c == 'T' ) );
    if( ! free )
     for( int dy = 0 ; dy < scale ; ++dy )
      for( int dx = 0 ; dx < scale ; ++dx )
       m.set_blocked( x * scale + dx , y * scale + dy , true );
    }
   }
  return( m );
  }

 /// the width \f$ W \f$
 [[nodiscard]] int width( void ) const noexcept { return( f_width ); }

 /// the height \f$ H \f$
 [[nodiscard]] int height( void ) const noexcept { return( f_height ); }

 /// true if \f$ ( x , y ) \f$ is in the map
 [[nodiscard]] bool inside( int x , int y ) const noexcept {
  return( ( x >= 0 ) && ( y >= 0 ) && ( x < f_width ) && ( y < f_height ) );
  }

 /// true if \f$ ( x , y ) \f$, in the map, is blocked
 [[nodiscard]] bool blocked( int x , int y ) const noexcept {
  return( v_blocked[ std::size_t( y ) * std::size_t( f_width ) +
                     std::size_t( x ) ] != 0 );
  }

 /// blocks or frees \f$ ( x , y ) \f$
 void set_blocked( int x , int y , bool b ) {
  if( ! inside( x , y ) )
   throw( std::invalid_argument( "NavMap::set_blocked: bad cell" ) );
  v_blocked[ std::size_t( y ) * std::size_t( f_width ) +
             std::size_t( x ) ] = b ? 1 : 0;
  }

 private:

 int f_width;                                  ///< the width
 int f_height;                                 ///< the height
 std::vector< std::uint8_t > v_blocked;        ///< 1 if blocked, by rows

 };  // end( class( NavMap ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS NavOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the parameters of RobotNav2D, by default those of the experiments

//! [domains-sig-navoptions]
struct NavOptions {

 int length = 25;                ///< the poses checked along a motion, which
                                 ///< moves by length - 1 cells
 int half = 16;                  ///< the footprint covers the cells
                                 ///< [ - half , half - 1 ] around the pose
 double radius = 25;             ///< a goal is closer than this to the goal
                                 ///< point
 bool factors = false;           ///< the cost of a motion multiplied by a
                                 ///< factor of its first cell
 std::uint64_t seed = 1;         ///< the seed of the factors
 unsigned repeat = 1;            ///< the check of a diagonal motion is done
                                 ///< this many times, so that it is the
                                 ///< expensive one if repeat > 1

 };  // end( struct( NavOptions ) )
//! [domains-sig-navoptions]

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS RobotNav2D ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a square robot moving on a NavMap by long straight motions
/** A state is the cell \f$ ( x , y ) \f$ of the center of the robot, whose
 * footprint is the square of the cells \f$ [ x - r , x + r - 1 ] \times
 * [ y - r , y + r - 1 ] \f$, \f$ r \f$ = half (32 cells a side by default).
 * The 8 actions move by \f$ L - 1 \f$ cells, \f$ L \f$ = length, in the
 * directions up, up-right, right, right-down, down, down-left, left and
 * left-up (in this order, y growing downwards); the edge exists if the
 * footprint is in the map and on free cells at each of the \f$ L \f$ poses
 * \f$ ( x , y ) + i \, d \f$, \f$ 0 \leq i < L \f$, which the evaluation
 * checks cell by cell, as a collision checker does, so that it is the
 * expensive part. Its cost is the length of the motion, \f$ L - 1 \f$ or
 * \f$ ( L - 1 ) \sqrt{2} \f$, times, with factors, a factor of the first
 * cell drawn uniformly from \f$ \{ 1 , \dots , 99 \} \f$ by the seed. The
 * goals are the states closer than \f$ R \f$ = radius to the goal point
 * \f$ \gamma \f$, and the heuristic is
 * \f$ h( s ) = \max\{ 0 , \| s - \gamma \| - R \} \f$: every motion costs
 * at least its length, so that \f$ h \f$ is admissible and consistent, the
 * pairwise heuristic \f$ h( s , t ) = \| s - t \| \f$ is forward-backward
 * consistent, and \f$ h( t ) \leq h( t , s ) + h( s ) \f$, which EPASE
 * needs. The lazy edge checks the last pose only, with the cost of the
 * motion. With repeat \f$ > 1 \f$ the check of the diagonal motions is done
 * repeat times, which makes them the expensive edges of GePA*SE (4 out of
 * 8, with a ratio of their times of about repeat, 30 or 300 in the paper).
 *
 * The experiments used the maps hrt201n, den501d, den520d, ht_chantry and
 * brc203d scaled by 5, 5, 5, 10 and 5, with 50 pairs of start and goal
 * points each, a heuristic weight of 50, and either no factors or the
 * random factors; the copies in data/robot-nav-2d come from the code of
 * those experiments (MIT license). Their heuristic was the distance to the
 * goal point, 0 within the radius, which is admissible but not consistent;
 * the one above is used instead. */

class RobotNav2D {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using state_type = std::uint32_t;     ///< x in the low half, y in the high
 using cost_type = double;             ///< the cost type

 /// constructor, taking the map (which must outlive the RobotNav2D), the
 /// start cell and the goal point

 RobotNav2D( const NavMap & map , int sx , int sy , int gx , int gy ,
             const NavOptions & o = {} )
  : f_map( map ) , f_o( o ) , f_start( pack( sx , sy ) ) , f_gx( gx ) ,
    f_gy( gy ) {
  if( ! map.inside( sx , sy ) || ( map.width() > 65535 ) ||
      ( map.height() > 65535 ) || ( o.length < 1 ) || ( o.half < 0 ) )
   throw( std::invalid_argument( "RobotNav2D: bad map, cell or options" ) );
  if( o.factors ) {
   std::mt19937_64 gen( o.seed );
   v_factor.resize( std::size_t( map.width() ) * std::size_t( map.height() ) );
   for( auto & f : v_factor )
    f = std::uint8_t( 1 + random_below( gen , 99 ) );
   }
  }

/*-------------------------------- READING ---------------------------------*/

 /// the pairs of start and goal points of the experiments, from the files
 /// nav2d_starts.txt (x, y and two more numbers per line) and
 /// nav2d_goals.txt (x, y per line) of a map, scaled and rounded to cells,
 /// as { sx , sy , gx , gy }

 [[nodiscard]] static std::vector< std::array< int , 4 > > read_points(
                       std::istream & starts , std::istream & goals ,
                       int scale ) {
  std::vector< std::array< int , 4 > > v;
  double sx , sy , c , l , gx , gy;
  while( ( starts >> sx >> sy >> c >> l ) && ( goals >> gx >> gy ) )
   v.push_back( { int( std::lround( sx * scale ) ) ,
                  int( std::lround( sy * scale ) ) ,
                  int( std::lround( gx * scale ) ) ,
                  int( std::lround( gy * scale ) ) } );
  return( v );
  }

/*------------------------- METHODS OF THE DOMAIN --------------------------*/

 /// the state of the cell ( x , y )
 [[nodiscard]] static state_type pack( int x , int y ) noexcept {
  return( state_type( x ) | ( state_type( y ) << 16 ) );
  }

 /// the x of s
 [[nodiscard]] static int x_of( state_type s ) noexcept {
  return( int( s & 0xFFFFu ) );
  }

 /// the y of s
 [[nodiscard]] static int y_of( state_type s ) noexcept {
  return( int( s >> 16 ) );
  }

 [[nodiscard]] state_type initial_state( void ) const noexcept {
  return( f_start );
  }

 [[nodiscard]] bool is_goal( state_type s ) const noexcept {
  return( to_goal( s ) < f_o.radius );
  }

 [[nodiscard]] std::uint64_t hash( state_type s ) const noexcept {
  return( ADS::splitmix64( s ) );
  }

 [[nodiscard]] cost_type heuristic( state_type s ) const noexcept {
  return( std::max( 0.0 , to_goal( s ) - f_o.radius ) );
  }

 [[nodiscard]] cost_type pairwise_heuristic( state_type s , state_type t )
  const noexcept {
  return( std::hypot( double( x_of( s ) - x_of( t ) ) ,
                      double( y_of( s ) - y_of( t ) ) ) );
  }

 /// the 8 actions
 template< class F >
 void for_each_action( state_type , F && f ) const {
  for( std::size_t a = 0 ; a < 8 ; ++a )
   f( a );
  }

 /// true if a is diagonal and its check is repeated
 [[nodiscard]] bool expensive_edge( state_type , std::size_t a )
  const noexcept { return( ( f_o.repeat > 1 ) && ( a & 1 ) ); }

 /// the motion a from s, checked at each pose
 template< class F >
 void evaluate_edge( state_type s , cost_type , std::size_t a ,
                     F && f ) const {
  const int x = x_of( s ) , y = y_of( s );
  const unsigned times = expensive_edge( s , a ) ? f_o.repeat : 1;
  bool ok = true;
  for( unsigned k = 0 ; ok && ( k < times ) ; ++k )
   for( int i = 0 ; ok && ( i < f_o.length ) ; ++i )
    ok = fits( x + i * c_dx[ a ] , y + i * c_dy[ a ] );
  if( ok )
   emit( s , a , f );
  }

 /// the motion a from s, its last pose checked
 template< class F >
 void lazy_edge( state_type s , cost_type , std::size_t a , F && f ) const {
  const int n = f_o.length - 1;
  if( fits( x_of( s ) + n * c_dx[ a ] , y_of( s ) + n * c_dy[ a ] ) )
   emit( s , a , f );
  }

 /// the successors of s, all the motions evaluated
 template< class F >
 void for_each_successor( state_type s , cost_type h , F && f ) const {
  for( std::size_t a = 0 ; a < 8 ; ++a )
   evaluate_edge( s , h , a , f );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 static constexpr std::array< int , 8 > c_dx = { 0 , 1 , 1 , 1 , 0 , -1 ,
                                                 -1 , -1 };
 static constexpr std::array< int , 8 > c_dy = { 1 , 1 , 0 , -1 , -1 , -1 ,
                                                 0 , 1 };

 // the distance of s from the goal point
 double to_goal( state_type s ) const noexcept {
  return( std::hypot( double( x_of( s ) - f_gx ) ,
                      double( y_of( s ) - f_gy ) ) );
  }

 // true if the footprint at ( x , y ) is in the map and free, checked cell
 // by cell
 bool fits( int x , int y ) const noexcept {
  const int r = f_o.half;
  if( ! f_map.inside( x - r , y - r ) ||
      ( r && ! f_map.inside( x + r - 1 , y + r - 1 ) ) )
   return( false );
  for( int v = y - r ; v < y + std::max( r , 1 ) ; ++v )
   for( int u = x - r ; u < x + std::max( r , 1 ) ; ++u )
    if( f_map.blocked( u , v ) )
     return( false );
  return( true );
  }

 // calls f with the successor of the motion a from s and its cost
 template< class F >
 void emit( state_type s , std::size_t a , F & f ) const {
  const int n = f_o.length - 1;
  const int x = x_of( s ) , y = y_of( s );
  const state_type t = pack( x + n * c_dx[ a ] , y + n * c_dy[ a ] );
  double c = ( a & 1 ) ? double( n ) * std::sqrt( 2.0 ) : double( n );
  if( ! v_factor.empty() )
   c *= v_factor[ std::size_t( y ) * std::size_t( f_map.width() ) +
                  std::size_t( x ) ];
  f( t , c , heuristic( t ) );
  }

 const NavMap & f_map;                         ///< the map
 NavOptions f_o;                               ///< the parameters
 state_type f_start;                           ///< the start
 int f_gx;                                     ///< the goal point
 int f_gy;
 std::vector< std::uint8_t > v_factor;         ///< the factors, if any

 };  // end( class( RobotNav2D ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_RobotNav2D */

/*--------------------------------------------------------------------------*/
/*------------------------- End File RobotNav2D.h --------------------------*/
/*--------------------------------------------------------------------------*/
