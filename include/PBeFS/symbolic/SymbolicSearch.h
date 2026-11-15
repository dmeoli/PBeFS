/*--------------------------------------------------------------------------*/
/*------------------------- File SymbolicSearch.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template SymbolicSearch, the searches of a
 * SymbolicDomain on sets of states represented by binary decision diagrams
 * [see SymbolicModel.h]: the breadth-first search, Dijkstra's algorithm and
 * BDDA*, forward, backward and bidirectional, with the images of the
 * partitions of the transition relation computed one after the other, as
 * tasks of Lace or on the workers of an executor of ADS.
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

#ifndef __PBeFS_SymbolicSearch
 #define __PBeFS_SymbolicSearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cstddef>

#include <cstdint>

#include <limits>

#include <map>

#include <optional>

#include <stdexcept>

#include <utility>

#include <vector>

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/symbolic/Bdd.h"

#include "PBeFS/symbolic/SymbolicModel.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the direction of a symbolic search

enum class SymbolicDirection {
 forward ,          ///< from the initial state, by images
 backward ,         ///< from the goals, by preimages
 bidirectional      ///< both, the cheaper image first
 };

/// how the images of the partitions of a set are computed

enum class SymbolicImages {
 sequential ,       ///< one partition after the other, each relational
                    ///< product parallel inside Sylvan
 lace ,             ///< all the partitions at once as tasks of Lace
 executor           ///< the partitions spread over the workers of an
                    ///< executor of ADS, each submitting its products to
                    ///< Lace
 };

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SymbolicOptions --------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of a symbolic search

struct SymbolicOptions {
 /// the direction
 SymbolicDirection direction = SymbolicDirection::forward;
 /// the heuristic of the domain, if it has one, in each direction (BDDA*);
 /// without it the search is Dijkstra's algorithm
 bool heuristic = true;
 /// every operator costs 1: a breadth-first search, whose cost is the
 /// number of moves
 bool unit_costs = false;
 /// the heuristic terms of the domain, if it has them (without them h = 0
 /// but for the patterns and the perimeter)
 bool terms = true;
 /// the patterns (sets of variables) whose pattern databases raise the
 /// heuristic, with SymbolicModel::strengthen()
 std::vector< std::vector< std::uint32_t > > patterns;
 /// the sets of the perimeter around the targets that raise the
 /// heuristic and seed the pattern databases (0: none)
 std::size_t perimeter = 0;
 /// the nodes beyond which a partition of the transition relation is not
 /// merged with the next operator of the same cost (0: one operator per
 /// partition)
 std::size_t max_relation_nodes = 100000;
 /// how the images of the partitions are computed
 SymbolicImages images = SymbolicImages::sequential;
 /// the workers of the executor, with SymbolicImages::executor
 unsigned threads = 2;
 /// the search stops with Status::node_limit if, after a collection, the
 /// table of the nodes of Sylvan is filled above this fraction of its
 /// largest size
 double max_filled = 0.5;
 /// the path found is rebuilt from the sets expanded
 bool path = true;
 };  // end( struct( SymbolicOptions ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SymbolicStats ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the counters of a symbolic search, all of them independent of the
/// number of workers, since the diagrams are canonical

struct SymbolicStats {
 std::uint64_t steps = 0;           ///< sets expanded, i.e., buckets of
                                    ///< equal g and h, and each layer of
                                    ///< their moves of cost 0
 std::uint64_t forward_steps = 0;   ///< those of the forward direction
 std::uint64_t backward_steps = 0;  ///< those of the backward direction
 std::uint64_t images = 0;          ///< relational products of the
                                    ///< expansions, one per partition and
                                    ///< set expanded
 std::uint64_t path_images = 0;     ///< relational products to rebuild the
                                    ///< path
 double expanded = 0;               ///< states expanded
 std::uint64_t expanded_nodes = 0;  ///< nodes of the sets expanded, summed
 std::uint64_t max_set_nodes = 0;   ///< nodes of the largest set expanded
 std::uint64_t max_image_nodes = 0; ///< nodes of the largest image of a set
 std::uint64_t partitions = 0;      ///< partitions of the relation
 std::uint64_t relation_nodes = 0;  ///< their nodes
 std::uint64_t heuristic_nodes = 0; ///< nodes of the sets \f$ H_k \f$ used
 std::uint64_t abstraction_images = 0;  ///< relational products of the
                                        ///< perimeter and of the pattern
                                        ///< databases
 std::uint64_t abstraction_steps = 0;   ///< their sets expanded
 std::vector< std::uint64_t > images_by;  ///< relational products made by
                                          ///< each worker of the executor
 };  // end( struct( SymbolicStats ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SymbolicSolution --------------------------*/
/*--------------------------------------------------------------------------*/
/// the outcome of a symbolic search
/** path holds the values of the variables of each state, from the initial
 * state to a goal, if the search is solved and SymbolicOptions::path is
 * set; bound is a lower bound on the optimal cost (the cost itself if
 * solved). */

template< class Cost >
struct SymbolicSolution {
 Status status = Status::unsolvable;              ///< how the search ended
 Cost cost{};                                     ///< the optimal cost
 Cost bound{};                                    ///< a lower bound on C*
 std::vector< std::vector< std::uint32_t > > path;  ///< the states
 SymbolicStats stats;                             ///< the counters

 /// true if a path has been found
 [[nodiscard]] bool solved( void ) const noexcept {
  return( status == Status::solved );
  }
 };  // end( struct( SymbolicSolution ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SymbolicReachable -------------------------*/
/*--------------------------------------------------------------------------*/
/// the states reachable from the roots of a direction, by layers

struct SymbolicReachable {
 std::vector< double > layers;        ///< the states first reached after
                                      ///< i moves, layers[ 0 ] the roots
 double states = 0;                   ///< all of them
 std::uint64_t images = 0;            ///< relational products
 std::uint64_t max_set_nodes = 0;     ///< nodes of the largest layer
 std::uint64_t reached_nodes = 0;     ///< nodes of the set of all of them
 };  // end( struct( SymbolicReachable ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS SymbolicSearch ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the symbolic searches of a SymbolicDomain
/** A direction keeps its open list as sets \f$ O_{g,f} \f$ of the states
 * reached with cost g and \f$ f = g + h \f$, h of the direction, the
 * closed set C of the states expanded, the sets \f$ R_g \f$ of the states
 * reached with cost g, and the sequence of the sets expanded with their g.
 * A step takes the set of least f, and of least g among them, removes C
 * from it, expands it, i.e., computes its image (its preimage backward)
 * under the partitions of each cost c, adds the image, split by the sets
 * \f$ H_k \f$, to the sets \f$ O_{g+c,g+c+k} \f$, and repeats on the
 * image of cost 0 that stays in its own set until none is left; this is
 * BDDA* with buckets of g and h [Edelkamp and Reffel, KI 1998; Jensen et
 * al., AIJ 2008], which is Dijkstra's algorithm on sets with h = 0 and a
 * breadth-first search with unit costs, and it returns g when the set
 * expanded meets the target (the goals forward, the initial state
 * backward), which is optimal if h is consistent [see
 * docs/theory/symbolic-search.md].
 *
 * The bidirectional search keeps both directions, detects every state
 * reached by both, the incumbent U being the least \f$ g_F + g_B \f$ of
 * such a state, and stops when \f$ U \leq \max( f_F , f_B , g_F + g_B +
 * \epsilon ) \f$, \f$ f_D \f$ and \f$ g_D \f$ being the least f and g in
 * the open list of D and \f$ \epsilon \f$ the least cost of an operator;
 * a step expands the direction whose next set is estimated cheaper to
 * expand, its nodes times the ratio of the nodes of the image to those of
 * the set at its last step, as the choice by the estimated time of the next
 * step in SymBA* [Torralba et al., AIJ 2017], here by counts.
 *
 * A BddRuntime must be alive while the search exists. */

template< SymbolicDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class SymbolicSearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using cost_type = typename D::cost_type;              ///< the cost type
 using solution_type = SymbolicSolution< cost_type >;  ///< the outcome

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// the search of d with the options o, the model compiled here

 explicit SymbolicSearch( const D & d , const SymbolicOptions & o = {} )
  : f_o( check( o ) ) , f_model( d , o.max_relation_nodes , o.unit_costs ,
                                     o.terms ) {
  for( const auto & p : f_model.partitions() ) {
   v_rel.push_back( p.relation );
   v_vars.push_back( p.vars );
   }
  if( o.heuristic && ! o.unit_costs &&
      ( ! o.patterns.empty() || ( o.perimeter > 0 ) ) ) {
   if( o.direction != SymbolicDirection::backward )
    f_model.strengthen( false , o.patterns , o.perimeter );
   if( o.direction != SymbolicDirection::forward )
    f_model.strengthen( true , o.patterns , o.perimeter );
   }
  }

/*-------------------------------- OPTIONS ---------------------------------*/
 /// sets the Executor of SymbolicImages::executor, which must outlive the
 /// search; without it the search builds its own

 void set_executor( Executor & ex ) noexcept { f_ex = & ex; }

 /// the model of the domain

 [[nodiscard]] const SymbolicModel< D > & model( void ) const noexcept {
  return( f_model );
  }

/*------------------------------- THE SEARCH -------------------------------*/
 /// the search, in the direction of the options

 [[nodiscard]] solution_type solve( void ) {
  f_sol = solution_type();
  f_sol.stats.partitions = f_model.partitions().size();
  f_sol.stats.relation_nodes = f_model.relation_nodes();
  f_sol.stats.abstraction_images = f_model.abstraction_images();
  f_sol.stats.abstraction_steps = f_model.abstraction_steps();
  start_executor();
  if( f_o.direction == SymbolicDirection::bidirectional )
   bidirectional();
  else
   unidirectional( f_o.direction == SymbolicDirection::backward );
  f_own.reset();
  return( std::move( f_sol ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// all the states reachable from the initial state (from the goals if
 /// backward), by a breadth-first search of the sets, every move counted
 /// as one whatever its cost

 [[nodiscard]] SymbolicReachable reachable( bool backward = false ) {
  SymbolicReachable r;
  f_sol = solution_type();
  start_executor();
  Bdd reached = backward ? f_model.goal() : f_model.initial();
  Bdd frontier = reached;
  while( ! frontier.empty() ) {
   r.layers.push_back( f_model.count( frontier ) );
   r.max_set_nodes = std::max< std::uint64_t >( r.max_set_nodes ,
                                                frontier.nodes() );
   std::vector< Bdd > next;
   for( const auto & [ c , Y ] : successors( frontier , backward ) )
    next.push_back( Y );
   frontier = disjoin( next ) - reached;
   reached |= frontier;
   }
  f_own.reset();
  r.states = f_model.count( reached );
  r.images = f_sol.stats.images;
  r.reached_nodes = reached.nodes();
  return( r );
  }

/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/

 private:

/*----------------------------- PRIVATE TYPES ------------------------------*/

 using Levels = typename SymbolicModel< D >::Levels;
 using Key = std::pair< cost_type , cost_type >;      // ( f , g )

 /// one direction of the search
 struct Side {
  bool backward = false;                      // by preimages
  const Levels * H = nullptr;                 // the sets of its h
  Bdd root;                                   // its roots
  Bdd target;                                 // what it looks for
  std::map< Key , Bdd > open;                 // O_{g,f}, by ( f , g )
  Bdd closed;                                 // C
  std::map< cost_type , Bdd > reached;        // R_g
  std::vector< std::pair< cost_type , Bdd > > layers;  // the sets expanded
  double ratio = 1;                           // nodes of the last image
                                              // over those of its set
  };

 /// a state reached by both directions
 struct Meeting {
  Bdd states;                                 // the states
  cost_type g_f{};                            // their g forward
  cost_type g_b{};                            // and backward
  };

 static constexpr cost_type c_inf = std::numeric_limits< cost_type >::max();

/*---------------------------- PRIVATE METHODS -----------------------------*/

 static const SymbolicOptions & check( const SymbolicOptions & o ) {
  if( ! BddRuntime::running() )
   throw( std::logic_error( "SymbolicSearch: no BddRuntime is alive" ) );
  if( ( o.images == SymbolicImages::executor ) && ( o.threads == 0 ) )
   throw( std::invalid_argument( "SymbolicSearch: no threads" ) );
  return( o );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void start_executor( void ) {
  if( f_o.images != SymbolicImages::executor )
   return;
  f_sol.stats.images_by.assign( f_o.threads , 0 );
  if( ! f_ex )
   f_own.emplace( f_o.threads );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the image (preimage if backward) of S under the partitions of each cost

 std::map< cost_type , Bdd > successors( const Bdd & S , bool backward ) {
  const auto & part = f_model.partitions();
  const std::size_t n = part.size();
  std::vector< Bdd > out( n );
  switch( f_o.images ) {
   case( SymbolicImages::sequential ):
    for( std::size_t i = 0 ; i < n ; ++i )
     out[ i ] = backward ? preimage( S , v_rel[ i ] , v_vars[ i ] )
                         : image( S , v_rel[ i ] , v_vars[ i ] );
    break;
   case( SymbolicImages::lace ):
    images( S , v_rel , v_vars , backward , out );
    break;
   case( SymbolicImages::executor ): {
    const unsigned p = f_o.threads;
    Executor * ex = f_ex ? f_ex : & *f_own;
    auto & by = f_sol.stats.images_by;
    run_workers( ex , p , [ & ]( unsigned w ) {
     for( std::size_t i = w ; i < n ; i += p ) {
      image_into( S , v_rel[ i ] , v_vars[ i ] , backward , out[ i ] );
      ++by[ w ];
      }
     } );
    break;
    }
   }
  f_sol.stats.images += n;
  // the partitions of each cost are consecutive
  std::map< cost_type , Bdd > r;
  for( std::size_t i = 0 ; i < n ; ) {
   std::size_t j = i;
   while( ( j < n ) && ( part[ j ].cost == part[ i ].cost ) )
    ++j;
   std::vector< Bdd > group( out.begin() + std::ptrdiff_t( i ) ,
                             out.begin() + std::ptrdiff_t( j ) );
   Bdd Y = group.size() == 1 ? group[ 0 ] : disjoin( group );
   if( backward )
    Y &= f_model.valid();
   f_sol.stats.max_image_nodes = std::max< std::uint64_t >(
    f_sol.stats.max_image_nodes , Y.nodes() );
   r.emplace( part[ i ].cost , Y );
   i = j;
   }
  return( r );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the side of a direction, its roots in the open list

 Side make_side( bool backward ) {
  Side s;
  s.backward = backward;
  s.H = ( f_o.heuristic && ! f_o.unit_costs )
        ? & f_model.levels( backward ) : & v_zero;
  if( s.H != & v_zero )
   f_sol.stats.heuristic_nodes += f_model.heuristic_nodes( backward );
  s.root = backward ? f_model.goal() : f_model.initial();
  s.target = backward ? f_model.initial() : f_model.goal();
  insert( s , s.root , cost_type( 0 ) , c_inf );
  return( s );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // Y, reached with cost g, into the open list of s, split by h, but the
 // dead ends and the states with f >= U; returns what is inserted

 Bdd insert( Side & s , const Bdd & Y , cost_type g , cost_type U ) {
  std::vector< Bdd > in;
  for( const auto & [ k , Hk ] : *s.H ) {
   if( ( U != c_inf ) && ( g + k >= U ) )
    break;
   const Bdd X = Y & Hk;
   if( X.empty() )
    continue;
   auto & b = s.open[ Key( g + k , g ) ];
   b |= X;
   in.push_back( X );
   }
  Bdd all = in.size() == 1 ? in[ 0 ] : disjoin( in );
  if( ! all.empty() ) {
   auto & R = s.reached[ g ];
   R |= all;
   }
  return( all );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // false if the table of the nodes stays too full after a collection

 bool memory_ok( void ) {
  if( BddRuntime::filled() <= f_o.max_filled )
   return( true );
  BddRuntime::collect();
  return( BddRuntime::filled() <= f_o.max_filled );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the step of s on its first set, of key ( f , g ): it is expanded, and
 // then the layers of its moves of cost 0; meet( X , g ) is called on each
 // set X inserted with cost g, and found( F , g ) on each set F expanded
 // with cost g, which returns true to stop; U prunes as in insert()

 template< class Meet , class Found >
 bool step( Side & s , cost_type U , Meet && meet , Found && found ) {
  const auto it = s.open.begin();
  const cost_type f = it->first.first , g = it->first.second;
  Bdd frontier = it->second - s.closed;
  s.open.erase( it );
  const cost_type h = f - g;
  const auto level = s.H->find( h );
  while( ! frontier.empty() ) {
   s.closed |= frontier;
   s.layers.emplace_back( g , frontier );
   auto & st = f_sol.stats;
   ++st.steps;
   ++( s.backward ? st.backward_steps : st.forward_steps );
   const std::size_t fn = frontier.nodes();
   st.expanded += f_model.count( frontier );
   st.expanded_nodes += fn;
   st.max_set_nodes = std::max< std::uint64_t >( st.max_set_nodes , fn );
   if( found( frontier , g ) )
    return( true );
   Bdd same;
   std::vector< Bdd > all;
   for( const auto & [ c , Y0 ] : successors( frontier , s.backward ) ) {
    Bdd Y = Y0 - s.closed;
    if( Y.empty() )
     continue;
    all.push_back( Y );
    if( ( c == 0 ) && ( level != s.H->end() ) ) {
     same = Y & level->second;
     Y -= same;
     if( ! same.empty() ) {
      auto & R = s.reached[ g ];
      R |= same;
      meet( same , g );
      }
     }
    const Bdd X = insert( s , Y , g + c , U );
    if( ! X.empty() )
     meet( X , g + c );
    }
   s.ratio = double( nodes( all ) ) / double( std::max< std::size_t >( fn ,
                                                                       1 ) );
   frontier = same;
   }
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void unidirectional( bool backward ) {
  Side s = make_side( backward );
  while( ! s.open.empty() ) {
   const cost_type f = s.open.begin()->first.first;
   if( ! memory_ok() ) {
    f_sol.status = Status::node_limit;
    f_sol.bound = f;
    return;
    }
   cost_type at = 0;
   const bool done = step( s , c_inf , []( const Bdd & , cost_type ) {} ,
                           [ & ]( const Bdd & F , cost_type g ) {
    if( F.disjoint( s.target ) )
     return( false );
    at = g;
    return( true );
    } );
   if( done ) {
    f_sol.status = Status::solved;
    f_sol.cost = f_sol.bound = at;
    if( f_o.path ) {
     const Bdd t = one_state( s.layers.back().second & s.target );
     auto walk = rebuild( s , t , at , s.layers.size() - 1 );
     if( ! backward )
      std::reverse( walk.begin() , walk.end() );
     for( const auto & x : walk )
      f_sol.path.push_back( f_model.decode( x ) );
     }
    return;
    }
   }
  f_sol.status = Status::unsolvable;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void bidirectional( void ) {
  Side F = make_side( false );
  Side B = make_side( true );
  cost_type U = c_inf;
  std::optional< Meeting > best;
  const cost_type eps = f_model.min_cost();
  // the states X reached by d with cost g against the other direction
  const auto meet = [ & ]( const Side & other , bool forward , const Bdd & X ,
                           cost_type g ) {
   for( const auto & [ g2 , R ] : other.reached ) {
    if( ( U != c_inf ) && ( g + g2 >= U ) )
     break;
    if( ! X.disjoint( R ) ) {
     U = g + g2;
     best = forward ? Meeting{ X & R , g , g2 } : Meeting{ X & R , g2 , g };
     break;
     }
    }
   };
  meet( B , true , F.reached.begin()->second , cost_type( 0 ) );
  const auto least_g = []( const Side & s ) {
   cost_type m = c_inf;
   for( const auto & [ key , X ] : s.open )
    m = std::min( m , key.second );
   return( m );
   };
  cost_type lb = 0;
  while( ! F.open.empty() && ! B.open.empty() ) {
   lb = std::max( { lb , F.open.begin()->first.first ,
                    B.open.begin()->first.first ,
                    cost_type( least_g( F ) + least_g( B ) + eps ) } );
   if( U <= lb )
    break;
   if( ! memory_ok() ) {
    f_sol.status = Status::node_limit;
    f_sol.bound = lb;
    return;
    }
   const double eF = double( F.open.begin()->second.nodes() ) * F.ratio;
   const double eB = double( B.open.begin()->second.nodes() ) * B.ratio;
   const bool forward = eF <= eB;
   Side & s = forward ? F : B;
   const Side & other = forward ? B : F;
   step( s , U , [ & ]( const Bdd & X , cost_type g ) {
    meet( other , forward , X , g );
    } , []( const Bdd & , cost_type ) { return( false ); } );
   }
  if( U == c_inf ) {
   f_sol.status = Status::unsolvable;
   return;
   }
  f_sol.status = Status::solved;
  f_sol.cost = f_sol.bound = U;
  if( ! f_o.path )
   return;
  const Bdd m = one_state( best->states );
  auto front = rebuild( F , m , best->g_f , F.layers.size() );
  std::reverse( front.begin() , front.end() );
  const auto back = rebuild( B , m , best->g_b , B.layers.size() );
  for( const auto & x : front )
   f_sol.path.push_back( f_model.decode( x ) );
  for( std::size_t i = 1 ; i < back.size() ; ++i )
   f_sol.path.push_back( f_model.decode( back[ i ] ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a set with one state of X, which is not empty

 Bdd one_state( const Bdd & X ) const {
  return( f_model.state( f_model.decode( X ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the states from t, reached by s with cost g by a set expanded before
 // the i-th, back to a root of s: each one is in a set expanded before
 // that of the one before it, with g lower by the cost of the move

 std::vector< Bdd > rebuild( const Side & s , Bdd t , cost_type g ,
                             std::size_t i ) {
  const auto & part = f_model.partitions();
  std::vector< Bdd > walk{ t };
  while( t.disjoint( s.root ) ) {
   std::map< cost_type , Bdd > back;   // the moves into t, by cost
   bool moved = false;
   while( ( i > 0 ) && ! moved ) {
    const auto & [ gj , L ] = s.layers[ --i ];
    if( gj > g )
     continue;
    const cost_type c = g - gj;
    auto it = back.find( c );
    if( it == back.end() ) {
     std::vector< Bdd > into;
     for( std::size_t p = 0 ; p < part.size() ; ++p )
      if( part[ p ].cost == c ) {
       into.push_back( s.backward ? image( t , v_rel[ p ] , v_vars[ p ] )
                                  : preimage( t , v_rel[ p ] , v_vars[ p ] ) );
       ++f_sol.stats.path_images;
       }
     it = back.emplace( c , disjoin( into ) ).first;
     }
    const Bdd X = it->second & L;
    if( ! X.empty() ) {
     t = one_state( X );
     g = gj;
     walk.push_back( t );
     moved = true;
     }
    }
   if( ! moved )
    throw( std::logic_error( "SymbolicSearch: the path cannot be rebuilt" ) );
   }
  return( walk );
  }

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SymbolicOptions f_o;                        ///< the options
 SymbolicModel< D > f_model;                 ///< the diagrams of the domain
 std::vector< Bdd > v_rel;                   ///< the relations
 std::vector< Bdd > v_vars;                  ///< and their variables
 Levels v_zero = { { cost_type( 0 ) , Bdd::one() } };  ///< h = 0
 Executor * f_ex = nullptr;                  ///< given, or nullptr
 std::optional< Executor > f_own;            ///< built by the search
 solution_type f_sol;                        ///< the outcome

 };  // end( class( SymbolicSearch ) )

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// the symbolic search of d with the options o (a BddRuntime alive)

template< SymbolicDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
[[nodiscard]] SymbolicSolution< typename D::cost_type > symbolic_search(
 const D & d , const SymbolicOptions & o = {} )
{
 SymbolicSearch< D , Executor > s( d , o );
 return( s.solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the symbolic breadth-first search of d, every move costing 1

template< SymbolicDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
[[nodiscard]] SymbolicSolution< typename D::cost_type > symbolic_bfs(
 const D & d , SymbolicOptions o = {} )
{
 o.unit_costs = true;
 return( symbolic_search< D , Executor >( d , o ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the states of d reachable from the initial state (from the goals if
/// backward), by layers

template< SymbolicDomain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
[[nodiscard]] SymbolicReachable symbolic_reachable(
 const D & d , bool backward = false , const SymbolicOptions & o = {} )
{
 SymbolicSearch< D , Executor > s( d , o );
 return( s.reachable( backward ) );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SymbolicSearch */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SymbolicSearch.h -----------------------*/
/*--------------------------------------------------------------------------*/
