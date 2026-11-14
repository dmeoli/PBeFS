/*--------------------------------------------------------------------------*/
/*----------------------------- File HPAStar.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template HPAStar, Hierarchical Path-Finding A*
 * on a GridMap [see Botea, Müller and Schaeffer, Journal of Game
 * Development 1, 2004], with one level of abstraction.
 *
 * The map is split in clusters of \f$ C \times C \f$ cells. Along the
 * border of two adjacent clusters, every maximal run of pairs of free cells
 * facing each other is an entrance, which gets a transition (a pair of
 * abstract nodes, one on each side, joined by the move between them) in
 * its middle if it is shorter than a split length, and one at each end
 * otherwise. Inside a cluster, its abstract nodes are joined by the
 * shortest paths that stay in the cluster, computed once by Dijkstra's
 * algorithm and kept. A query inserts the start and the goal into the
 * abstract graph, joined to the nodes of their clusters, searches it with
 * A* on the distance of the Grid, and refines the abstract path by
 * concatenating the paths kept; when the start and the goal are in the
 * same cluster, the shortest path inside it is a candidate too.
 *
 * The path is then smoothed as in the original method: from each of its
 * cells, a straight walk in each of the directions of the moves that
 * meets a later cell of the path replaces the piece between them if it
 * costs less, so that the detours through the transitions are cut.
 *
 * The path is not optimal in general: it goes through the transitions, and it
 * is optimal among such paths before the smoothing. The search is also not
 * complete when a cluster is disconnected inside, since a pair of an entrance
 * may be reachable within its cluster while its transition is not, as in the
 * original method. The query costs a search of the abstract graph and of
 * two clusters, instead of a search of the map, which is the point.
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

#ifndef __PBeFS_HPAStar
 #define __PBeFS_HPAStar
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <queue>

#include <stdexcept>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/Grid.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS HPAOptions ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of HPAStar

struct HPAOptions {

 int cluster = 10;            ///< the side C of a cluster
 int split = 6;               ///< the entrances at least this long get two
                              ///< transitions
 bool smooth = true;          ///< the path is smoothed [see HPAStar]

 };  // end( struct( HPAOptions ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS HPAStar ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Hierarchical Path-Finding A* on a GridMap with 4 or 8 moves
/** The constructor builds the abstract graph, whose size is given by
 * nodes() and edges(), and the expansions of its Dijkstra searches by
 * preprocessing(); solve() answers a query, with the counters of the
 * insertion of the start and of the goal and of the abstract search
 * (expanded, generated), and the costs of the moves of the Grid. */

template< int Moves = 4 >
class HPAStar {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using GridT = Grid< Moves >;                      ///< the grid
 using State = typename GridT::state_type;         ///< a cell
 using Cost = typename GridT::cost_type;           ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the map (which must outlive the search) and the
 /// options; throws std::invalid_argument if the map has no free cell or
 /// the side of a cluster is not positive

 explicit HPAStar( const GridMap & map , const HPAOptions & options = {} )
  : f_map( map ) , f_o( options ) , f_grid( make_grid( map ) ) {
  if( f_o.cluster <= 0 )
   throw( std::invalid_argument( "HPAStar: clusters of side <= 0" ) );
  f_cw = ( map.width() + f_o.cluster - 1 ) / f_o.cluster;
  f_ch = ( map.height() + f_o.cluster - 1 ) / f_o.cluster;
  v_members.resize( std::size_t( f_cw ) * f_ch );
  build_entrances();
  for( std::size_t c = 0 ; c < v_members.size() ; ++c )
   build_intra( c );
  }

/*-------------------------------- QUERIES ---------------------------------*/
 /// the abstract nodes, the abstract edges, and the expansions of the
 /// preprocessing

 [[nodiscard]] std::size_t nodes( void ) const noexcept {
  return( v_cell.size() );
  }

 [[nodiscard]] std::size_t edges( void ) const noexcept {
  std::size_t n = 0;
  for( const auto & a : v_adj )
   n += a.size();
  return( n );
  }

 [[nodiscard]] std::uint64_t preprocessing( void ) const noexcept {
  return( f_pre );
  }

/*-------------------------------- SOLVING ---------------------------------*/
 /// the path from ( sx , sy ) to ( gx , gy ); throws std::invalid_argument
 /// if either is blocked

 [[nodiscard]] SolutionT solve( int sx , int sy , int gx , int gy ) const {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  if( ! f_map.cost( sx , sy ) || ! f_map.cost( gx , gy ) )
   throw( std::invalid_argument( "HPAStar::solve: blocked start or goal" ) );
  const State s = f_grid.index( sx , sy ) , g = f_grid.index( gx , gy );
  if( s == g ) {
   sol.status = Status::solved;
   sol.cost = Cost( 0 );
   sol.path = { s };
   st.seconds = clock.seconds();
   return( sol );
   }

  // the abstract graph with the start and the goal: their edges apart
  const std::size_t n = v_cell.size();
  const std::size_t S = n , G = n + 1;
  std::vector< Edge > out_s;                    // from S
  std::vector< std::vector< Edge > > to_g( n );  // node -> G
  const std::size_t cs = cluster_of( s ) , cg = cluster_of( g );
  {
   const auto dist = dijkstra( s , cs , st );
   for( std::size_t a : v_members[ cs ] )
    if( dist.reached( v_cell[ a ] ) )
     out_s.push_back( Edge{ a , dist.cost( v_cell[ a ] ) ,
                            dist.path( v_cell[ a ] ) } );
   }
  for( std::size_t a : v_members[ cg ] ) {
   const auto dist = dijkstra( v_cell[ a ] , cg , st );
   if( dist.reached( g ) )
    to_g[ a ].push_back( Edge{ G , dist.cost( g ) , dist.path( g ) } );
   }
  // the candidate inside one cluster
  Cost best = infinite_cost< Cost >();
  std::vector< State > best_path;
  if( cs == cg ) {
   const auto dist = dijkstra( s , cs , st );
   if( dist.reached( g ) ) {
    best = dist.cost( g );
    best_path = dist.path( g );
    }
   }
  // A* on the abstract graph, h the distance of the Grid
  struct Item {
   Cost f , g;
   std::size_t v;
   bool operator>( const Item & o ) const {
    return( ( f > o.f ) || ( ( f == o.f ) && ( g < o.g ) ) );
    }
   };
  std::vector< Cost > dist( n + 2 , infinite_cost< Cost >() );
  std::vector< std::size_t > par( n + 2 , std::size_t( -1 ) );
  std::vector< const std::vector< State > * > via( n + 2 , nullptr );
  std::vector< char > closed( n + 2 , 0 );
  std::priority_queue< Item , std::vector< Item > , std::greater< Item > > q;
  const auto h = [ & ]( std::size_t v ) {
   return( v == G ? Cost( 0 ) : f_grid.pairwise_heuristic( cell( v , s , g ) ,
                                                           g ) );
   };
  dist[ S ] = Cost( 0 );
  q.push( Item{ h( S ) , Cost( 0 ) , S } );
  while( ! q.empty() ) {
   const Item e = q.top();
   q.pop();
   if( closed[ e.v ] || ( e.g != dist[ e.v ] ) )
    continue;
   closed[ e.v ] = 1;
   ++st.expanded;
   if( e.v == G )
    break;
   if( e.g >= best )
    break;                              // the candidate is not worse
   const auto relax = [ & ]( const Edge & x ) {
    ++st.generated;
    const Cost d = e.g + x.cost;
    if( d < dist[ x.to ] ) {
     dist[ x.to ] = d;
     par[ x.to ] = e.v;
     via[ x.to ] = &x.path;
     q.push( Item{ d + h( x.to ) , d , x.to } );
     }
    };
   if( e.v == S )
    for( const Edge & x : out_s )
     relax( x );
   else {
    for( const Edge & x : v_adj[ e.v ] )
     relax( x );
    for( const Edge & x : to_g[ e.v ] )
     relax( x );
    }
   }
  if( ( dist[ G ] < best ) && ( dist[ G ] != infinite_cost< Cost >() ) ) {
   best = dist[ G ];
   std::vector< const std::vector< State > * > pieces;
   for( std::size_t v = G ; v != S ; v = par[ v ] )
    pieces.push_back( via[ v ] );
   best_path = { s };
   for( auto it = pieces.rbegin() ; it != pieces.rend() ; ++it )
    best_path.insert( best_path.end() , ( *it )->begin() + 1 ,
                      ( *it )->end() );
   }
  if( best != infinite_cost< Cost >() ) {
   if( f_o.smooth )
    best = smooth( best_path );
   sol.status = Status::solved;
   sol.cost = best;
   sol.path = std::move( best_path );
   }
  st.expanded_by = { st.expanded };
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// an abstract edge: its head, its cost and its path of cells, from the
 /// cell of its tail to that of its head
 struct Edge {
  std::size_t to;
  Cost cost;
  std::vector< State > path;
  };

 /// the distances of a Dijkstra search inside a cluster
 struct Dist {
  std::unordered_map< State , std::pair< Cost , State > > d;  // cell ->
                                                // ( distance , parent )
  State root;
  [[nodiscard]] bool reached( State c ) const { return( d.count( c ) ); }
  [[nodiscard]] Cost cost( State c ) const { return( d.at( c ).first ); }
  [[nodiscard]] std::vector< State > path( State c ) const {
   std::vector< State > v;
   for( State x = c ; ; x = d.at( x ).second ) {
    v.push_back( x );
    if( x == root )
     break;
    }
   std::reverse( v.begin() , v.end() );
   return( v );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a Grid on map, to read its moves and its distance
 static GridT make_grid( const GridMap & map ) {
  for( int y = 0 ; y < map.height() ; ++y )
   for( int x = 0 ; x < map.width() ; ++x )
    if( map.cost( x , y ) )
     return( GridT( map , x , y , x , y ) );
  throw( std::invalid_argument( "HPAStar: no free cell" ) );
  }

 /// the cluster of a cell
 [[nodiscard]] std::size_t cluster_of( State c ) const noexcept {
  const int x = int( c % State( f_map.width() ) );
  const int y = int( c / State( f_map.width() ) );
  return( std::size_t( y / f_o.cluster ) * f_cw + x / f_o.cluster );
  }

 /// the cell of the abstract node v, S and G being s and g
 [[nodiscard]] State cell( std::size_t v , State s , State g ) const {
  if( v < v_cell.size() )
   return( v_cell[ v ] );
  return( v == v_cell.size() ? s : g );
  }

 /// the abstract node of the cell c, made if new
 std::size_t node_of( State c ) {
  auto [ it , fresh ] = f_node.emplace( c , v_cell.size() );
  if( fresh ) {
   v_cell.push_back( c );
   v_adj.emplace_back();
   v_members[ cluster_of( c ) ].push_back( it->second );
   }
  return( it->second );
  }

 /// the cost of the move from the cell a to the cell b next to it
 [[nodiscard]] Cost move( State a , State b ) const {
  Cost c = infinite_cost< Cost >();
  f_grid.for_each_successor( a , Cost( 0 ) , [ & ]( State t , Cost x ,
                                                    Cost ) {
   if( t == b )
    c = x;
   } );
  return( c );
  }

 /// the transition between the cells a and b, on two sides of a border
 void transition( State a , State b ) {
  const std::size_t u = node_of( a ) , v = node_of( b );
  v_adj[ u ].push_back( Edge{ v , move( a , b ) , { a , b } } );
  v_adj[ v ].push_back( Edge{ u , move( b , a ) , { b , a } } );
  }

 /// the transitions of a run of pairs ( a_i , b_i ), i in [ 0 , n )
 template< class Pair >
 void entrance( int n , Pair && pair ) {
  if( n <= 0 )
   return;
  if( n < f_o.split ) {
   const auto [ a , b ] = pair( n / 2 );
   transition( a , b );
   return;
   }
  const auto [ a0 , b0 ] = pair( 0 );
  transition( a0 , b0 );
  const auto [ a1 , b1 ] = pair( n - 1 );
  transition( a1 , b1 );
  }

 /// the entrances of all the borders
 void build_entrances( void ) {
  const int W = f_map.width() , H = f_map.height() , C = f_o.cluster;
  // the vertical borders, between the columns x and x + 1
  for( int x = C - 1 ; x + 1 < W ; x += C )
   for( int y0 = 0 ; y0 < H ; y0 += C ) {
    const int y1 = std::min( H , y0 + C );
    for( int y = y0 ; y < y1 ; ) {
     if( ! ( f_map.cost( x , y ) && f_map.cost( x + 1 , y ) ) ) {
      ++y;
      continue;
      }
     int e = y;
     while( ( e < y1 ) && f_map.cost( x , e ) && f_map.cost( x + 1 , e ) )
      ++e;
     const int start = y;
     entrance( e - start , [ & ]( int i ) {
      return( std::pair( f_grid.index( x , start + i ) ,
                         f_grid.index( x + 1 , start + i ) ) );
      } );
     y = e;
     }
    }
  // the horizontal borders, between the rows y and y + 1
  for( int y = C - 1 ; y + 1 < H ; y += C )
   for( int x0 = 0 ; x0 < W ; x0 += C ) {
    const int x1 = std::min( W , x0 + C );
    for( int x = x0 ; x < x1 ; ) {
     if( ! ( f_map.cost( x , y ) && f_map.cost( x , y + 1 ) ) ) {
      ++x;
      continue;
      }
     int e = x;
     while( ( e < x1 ) && f_map.cost( e , y ) && f_map.cost( e , y + 1 ) )
      ++e;
     const int start = x;
     entrance( e - start , [ & ]( int i ) {
      return( std::pair( f_grid.index( start + i , y ) ,
                         f_grid.index( start + i , y + 1 ) ) );
      } );
     x = e;
     }
    }
  }

 /// Dijkstra's algorithm from the cell r, inside the cluster c
 Dist dijkstra( State r , std::size_t c , SearchStats & st ) const {
  Dist out;
  out.root = r;
  using QI = std::pair< Cost , State >;
  std::priority_queue< QI , std::vector< QI > , std::greater< QI > > q;
  out.d[ r ] = { Cost( 0 ) , r };
  q.push( { Cost( 0 ) , r } );
  std::unordered_map< State , char > done;
  while( ! q.empty() ) {
   const Cost d = q.top().first;
   const State u = q.top().second;
   q.pop();
   if( done[ u ] || ( d != out.d[ u ].first ) )
    continue;
   done[ u ] = 1;
   ++st.expanded;
   f_grid.for_each_successor( u , Cost( 0 ) , [ & ]( State t , Cost x ,
                                                     Cost ) {
    if( cluster_of( t ) != c )
     return;
    ++st.generated;
    const Cost nd = d + x;
    auto it = out.d.find( t );
    if( ( it == out.d.end() ) || ( nd < it->second.first ) ) {
     out.d[ t ] = { nd , u };
     q.push( { nd , t } );
     }
    } );
   }
  return( out );
  }

 /// smooths the path p, returning its cost: from each cell, the straight
 /// walks in the directions of the moves that meet a later cell of p
 /// replace the piece between the two if they cost less
 Cost smooth( std::vector< State > & p ) const {
  static constexpr int dx[ 8 ] = { 0 , -1 , 1 , 0 , -1 , 1 , -1 , 1 };
  static constexpr int dy[ 8 ] = { -1 , 0 , 0 , 1 , -1 , -1 , 1 , 1 };
  const int W = f_map.width();
  for( std::size_t i = 0 ; i + 2 < p.size() ; ++i ) {
   std::unordered_map< State , std::size_t > at;
   std::vector< Cost > upto( p.size() , Cost( 0 ) );   // cost from p[ i ]
   for( std::size_t j = i + 1 ; j < p.size() ; ++j ) {
    at[ p[ j ] ] = j;
    upto[ j ] = upto[ j - 1 ] + move( p[ j - 1 ] , p[ j ] );
    }
   std::size_t best_j = 0;
   std::vector< State > best_walk;
   Cost best_save = Cost( 0 );
   for( int k = 0 ; k < Moves ; ++k ) {
    std::vector< State > walk;
    Cost c( 0 );
    State x = p[ i ];
    for( ; ; ) {
     const int nx = int( x % State( W ) ) + dx[ k ];
     const int ny = int( x / State( W ) ) + dy[ k ];
     if( ! f_map.cost( nx , ny ) )
      break;
     const State y = f_grid.index( nx , ny );
     const Cost m = move( x , y );
     if( m == infinite_cost< Cost >() )
      break;                                 // a corner cut
     c = c + m;
     walk.push_back( y );
     x = y;
     const auto it = at.find( y );
     if( ( it != at.end() ) && ( it->second > i + 1 ) &&
         ( upto[ it->second ] - c > best_save ) ) {
      best_save = upto[ it->second ] - c;
      best_j = it->second;
      best_walk = walk;
      }
     }
    }
   if( best_j ) {
    std::vector< State > q( p.begin() , p.begin() + i + 1 );
    q.insert( q.end() , best_walk.begin() , best_walk.end() );
    q.insert( q.end() , p.begin() + best_j + 1 , p.end() );
    p.swap( q );
    }
   }
  Cost c( 0 );
  for( std::size_t j = 1 ; j < p.size() ; ++j )
   c = c + move( p[ j - 1 ] , p[ j ] );
  return( c );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the intra-cluster edges of the cluster c
 void build_intra( std::size_t c ) {
  SearchStats st;
  for( std::size_t a : v_members[ c ] ) {
   const auto d = dijkstra( v_cell[ a ] , c , st );
   for( std::size_t b : v_members[ c ] )
    if( ( a != b ) && d.reached( v_cell[ b ] ) )
     v_adj[ a ].push_back( Edge{ b , d.cost( v_cell[ b ] ) ,
                                 d.path( v_cell[ b ] ) } );
   }
  f_pre += st.expanded;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const GridMap & f_map;                              ///< the map
 HPAOptions f_o;                                     ///< the options
 GridT f_grid;                                       ///< its moves
 int f_cw = 0;                                       ///< clusters per row
 int f_ch = 0;                                       ///< and per column
 std::vector< State > v_cell;                        ///< cell of a node
 std::vector< std::vector< Edge > > v_adj;           ///< the edges
 std::unordered_map< State , std::size_t > f_node;   ///< node of a cell
 std::vector< std::vector< std::size_t > > v_members;  ///< nodes by
                                                     ///< cluster
 std::uint64_t f_pre = 0;                            ///< expansions of
                                                     ///< the preprocessing

 };  // end( class( HPAStar ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_HPAStar */

/*--------------------------------------------------------------------------*/
/*--------------------------- End File HPAStar.h ---------------------------*/
/*--------------------------------------------------------------------------*/
