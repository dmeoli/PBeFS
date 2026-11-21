/*--------------------------------------------------------------------------*/
/*------------------------- File SASHeuristics.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the heuristics of classical planning: the classes
 * RelaxedHeuristic, LMCut, SASPatternDatabase, SASMergeAndShrink,
 * SASCartesianAbstraction and SASAdditiveCartesian.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <chrono>

#include <functional>

#include <map>

#include <queue>

#include <stdexcept>

#include <utility>

#include "PBeFS/astar/SASHeuristics.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

namespace {

/// a min-queue of the pairs ( cost , fact )
using MinQueue = std::priority_queue< std::pair< int , int > ,
                                      std::vector< std::pair< int , int > > ,
                                      std::greater<> >;

/// empties q, keeping its storage
void clear( MinQueue & q )
{
 while( ! q.empty() )
  q.pop();
 }

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS RelaxedHeuristic -------------------------*/
/*--------------------------------------------------------------------------*/

RelaxedHeuristic::RelaxedHeuristic( const StripsTask & task , Kind kind )
 : f_task( task ) , f_kind( kind )
{
 std::vector< std::vector< int > > of( std::size_t( task.num_facts ) );
 for( std::size_t o = 0 ; o < task.ops.size() ; ++o ) {
  const auto & op = task.ops[ o ];
  for( int q : op.add ) {
   const int u = int( v_unary.size() );
   v_unary.push_back( { q , int( o ) , op.cost , int( op.pre.size() ) } );
   if( op.pre.empty() )
    v_free.push_back( u );
   for( int p : op.pre )
    of[ std::size_t( p ) ].push_back( u );
   }
  }
 v_first.push_back( 0 );
 for( const auto & l : of ) {
  v_of.insert( v_of.end() , l.begin() , l.end() );
  v_first.push_back( int( v_of.size() ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

int RelaxedHeuristic::operator()( const std::vector< int > & s ) const
{
 thread_local std::vector< int > cost , left , agg , best;
 thread_local std::vector< char > done;
 thread_local MinQueue q;
 const auto nf = std::size_t( f_task.num_facts );
 cost.assign( nf , dead_end );
 best.assign( nf , -1 );
 done.assign( nf , 0 );
 left.resize( v_unary.size() );
 agg.assign( v_unary.size() , 0 );
 for( std::size_t u = 0 ; u < v_unary.size() ; ++u )
  left[ u ] = v_unary[ u ].npre;
 clear( q );

 const auto relax = [ & ]( int f , int c , int u ) {
  if( c < cost[ std::size_t( f ) ] ) {
   cost[ std::size_t( f ) ] = c;
   best[ std::size_t( f ) ] = u;
   q.emplace( c , f );
   }
  };
 for( std::size_t v = 0 ; v < s.size() ; ++v )
  relax( f_task.fact( int( v ) , s[ v ] ) , 0 , -1 );
 for( int u : v_free )
  relax( v_unary[ std::size_t( u ) ].eff , v_unary[ std::size_t( u ) ].cost ,
         u );

 // the goal facts still to be settled
 std::size_t goals = 0;
 for( int g : f_task.goal )
  if( cost[ std::size_t( g ) ] != 0 )
   ++goals;
 const bool sum = f_kind != Kind::hmax;
 while( ( goals > 0 ) && ! q.empty() ) {
  const auto [ c , f ] = q.top();
  q.pop();
  if( done[ std::size_t( f ) ] )
   continue;
  done[ std::size_t( f ) ] = 1;
  if( std::binary_search( f_task.goal.begin() , f_task.goal.end() , f ) &&
      ( c > 0 ) )
   --goals;
  for( int k = v_first[ std::size_t( f ) ] ;
       k < v_first[ std::size_t( f ) + 1 ] ; ++k ) {
   const auto u = std::size_t( v_of[ std::size_t( k ) ] );
   agg[ u ] = sum ? agg[ u ] + c : std::max( agg[ u ] , c );
   if( --left[ u ] == 0 )
    relax( v_unary[ u ].eff , agg[ u ] + v_unary[ u ].cost , int( u ) );
   }
  }

 int h = 0;
 for( int g : f_task.goal ) {
  const int c = cost[ std::size_t( g ) ];
  if( c == dead_end )
   return( dead_end );
  h = f_kind == Kind::hmax ? std::max( h , c ) : h + c;
  }
 if( f_kind != Kind::hff )
  return( h );

 // the relaxed plan, from the goal back through the best supporters
 thread_local std::vector< char > used;
 thread_local std::vector< int > stack;
 used.assign( f_task.ops.size() , 0 );
 std::fill( done.begin() , done.end() , 0 );
 stack.assign( f_task.goal.begin() , f_task.goal.end() );
 h = 0;
 while( ! stack.empty() ) {
  const auto f = std::size_t( stack.back() );
  stack.pop_back();
  if( done[ f ] || ( best[ f ] < 0 ) )
   continue;
  done[ f ] = 1;
  const auto o = std::size_t( v_unary[ std::size_t( best[ f ] ) ].op );
  if( used[ o ] )
   continue;
  used[ o ] = 1;
  h += f_task.ops[ o ].cost;
  stack.insert( stack.end() , f_task.ops[ o ].pre.begin() ,
                f_task.ops[ o ].pre.end() );
  }
 return( h );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS LMCut -------------------------------*/
/*--------------------------------------------------------------------------*/

LMCut::LMCut( const StripsTask & task )
 : f_task( task ) , f_true( task.num_facts ) , f_goal( task.num_facts + 1 )
{
 for( const auto & op : task.ops ) {
  Op x{ op.pre , op.add , op.cost };
  if( x.pre.empty() )
   x.pre.push_back( f_true );
  v_op.push_back( std::move( x ) );
  }
 v_op.push_back( { task.goal , { f_goal } , 0 } );
 if( v_op.back().pre.empty() )
  v_op.back().pre.push_back( f_true );
 v_pre_of.resize( std::size_t( task.num_facts ) + 2 );
 v_add_of.resize( std::size_t( task.num_facts ) + 2 );
 for( std::size_t o = 0 ; o < v_op.size() ; ++o ) {
  for( int p : v_op[ o ].pre )
   v_pre_of[ std::size_t( p ) ].push_back( int( o ) );
  for( int a : v_op[ o ].add )
   v_add_of[ std::size_t( a ) ].push_back( int( o ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

int LMCut::operator()( const std::vector< int > & s ) const
{
 thread_local std::vector< int > cost , cur , left , opc , pcf , stack;
 thread_local std::vector< char > zone , seen , cut , done;
 thread_local MinQueue q;
 const std::size_t nf = std::size_t( f_task.num_facts ) + 2;
 const std::size_t no = v_op.size();
 cur.resize( no );
 for( std::size_t o = 0 ; o < no ; ++o )
  cur[ o ] = v_op[ o ].cost;

 int h = 0;
 for( ; ; ) {
  // h^max with the current costs
  cost.assign( nf , dead_end );
  done.assign( nf , 0 );
  opc.assign( no , 0 );
  left.resize( no );
  for( std::size_t o = 0 ; o < no ; ++o )
   left[ o ] = int( v_op[ o ].pre.size() );
  clear( q );
  const auto relax = [ & ]( int f , int c ) {
   if( c < cost[ std::size_t( f ) ] ) {
    cost[ std::size_t( f ) ] = c;
    q.emplace( c , f );
    }
   };
  for( std::size_t v = 0 ; v < s.size() ; ++v )
   relax( f_task.fact( int( v ) , s[ v ] ) , 0 );
  relax( f_true , 0 );
  while( ! q.empty() ) {
   const auto [ c , f ] = q.top();
   q.pop();
   if( done[ std::size_t( f ) ] )
    continue;
   done[ std::size_t( f ) ] = 1;
   for( int o : v_pre_of[ std::size_t( f ) ] ) {
    opc[ std::size_t( o ) ] = std::max( opc[ std::size_t( o ) ] , c );
    if( --left[ std::size_t( o ) ] == 0 )
     for( int a : v_op[ std::size_t( o ) ].add )
      relax( a , opc[ std::size_t( o ) ] + cur[ std::size_t( o ) ] );
    }
   }
  const int hg = cost[ std::size_t( f_goal ) ];
  if( hg == dead_end )
   return( dead_end );
  if( hg == 0 )
   return( h );

  // the precondition choice function
  pcf.assign( no , -1 );
  for( std::size_t o = 0 ; o < no ; ++o )
   if( left[ o ] == 0 ) {
    int b = v_op[ o ].pre.front();
    for( int p : v_op[ o ].pre )
     if( cost[ std::size_t( p ) ] > cost[ std::size_t( b ) ] )
      b = p;
    pcf[ o ] = b;
    }

  // the goal zone: the facts reaching the goal through cost zero
  zone.assign( nf , 0 );
  zone[ std::size_t( f_goal ) ] = 1;
  stack.assign( 1 , f_goal );
  while( ! stack.empty() ) {
   const int g = stack.back();
   stack.pop_back();
   for( int o : v_add_of[ std::size_t( g ) ] )
    if( ( pcf[ std::size_t( o ) ] >= 0 ) && ( cur[ std::size_t( o ) ] == 0 )
        && ! zone[ std::size_t( pcf[ std::size_t( o ) ] ) ] ) {
     zone[ std::size_t( pcf[ std::size_t( o ) ] ) ] = 1;
     stack.push_back( pcf[ std::size_t( o ) ] );
     }
   }

  // the cut: the operators into the zone from the facts reached before it
  seen.assign( nf , 0 );
  cut.assign( no , 0 );
  stack.clear();
  for( std::size_t v = 0 ; v < s.size() ; ++v )
   stack.push_back( f_task.fact( int( v ) , s[ v ] ) );
  stack.push_back( f_true );
  for( int f : stack )
   seen[ std::size_t( f ) ] = 1;
  int m = dead_end;
  while( ! stack.empty() ) {
   const int f = stack.back();
   stack.pop_back();
   for( int o : v_pre_of[ std::size_t( f ) ] ) {
    if( pcf[ std::size_t( o ) ] != f )
     continue;
    for( int a : v_op[ std::size_t( o ) ].add )
     if( zone[ std::size_t( a ) ] ) {
      if( ! cut[ std::size_t( o ) ] ) {
       cut[ std::size_t( o ) ] = 1;
       m = std::min( m , cur[ std::size_t( o ) ] );
       }
      }
     else
      if( ! seen[ std::size_t( a ) ] ) {
       seen[ std::size_t( a ) ] = 1;
       stack.push_back( a );
       }
    }
   }
  if( ( m == dead_end ) || ( m <= 0 ) )
   throw( std::logic_error( "LMCut::operator(): an empty cut" ) );
  h += m;
  for( std::size_t o = 0 ; o < no ; ++o )
   if( cut[ o ] )
    cur[ o ] -= m;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SASPatternDatabase ------------------------*/
/*--------------------------------------------------------------------------*/

SASPatternDatabase::SASPatternDatabase( const SASTask & task ,
                                        std::vector< int > pattern )
 : v_pattern( std::move( pattern ) )
{
 std::sort( v_pattern.begin() , v_pattern.end() );
 v_pattern.erase( std::unique( v_pattern.begin() , v_pattern.end() ) ,
                  v_pattern.end() );
 const std::size_t n = task.num_vars();
 std::vector< int > where( n , -1 );     // the place of v in the pattern
 std::size_t size = 1;
 for( std::size_t k = 0 ; k < v_pattern.size() ; ++k ) {
  const int v = v_pattern[ k ];
  if( ( v < 0 ) || ( std::size_t( v ) >= n ) )
   throw( std::invalid_argument(
             "SASPatternDatabase::SASPatternDatabase: no such variable" ) );
  where[ std::size_t( v ) ] = int( k );
  v_mult.push_back( size );
  const auto r = std::size_t( task.range[ std::size_t( v ) ] );
  if( size > ( std::size_t( 1 ) << 40 ) / r )
   throw( std::invalid_argument(
             "SASPatternDatabase::SASPatternDatabase: too many states" ) );
  size *= r;
  }

 // the abstract operators, which change some variable of the pattern
 struct AOp {
  std::vector< std::pair< int , int > > pre;  ///< ( place , value )
  std::vector< std::pair< int , int > > eff;  ///< ( place , post )
  std::vector< int > any;                     ///< places with no value
  int cost;
  };
 std::vector< AOp > ops;
 for( std::size_t o = 0 ; o < task.ops.size() ; ++o ) {
  AOp a;
  a.cost = task.ops[ o ].cost;
  for( const auto & e : task.ops[ o ].effects ) {
   const int k = where[ std::size_t( e.var ) ];
   if( k < 0 )
    continue;
   a.eff.emplace_back( k , e.post );
   if( e.pre < 0 )
    a.any.push_back( k );
   }
  if( a.eff.empty() )
   continue;
  for( const auto & p : task.preconditions( o ) )
   if( where[ std::size_t( p.var ) ] >= 0 )
    a.pre.emplace_back( where[ std::size_t( p.var ) ] , p.val );
  ops.push_back( std::move( a ) );
  }

 const std::size_t np = v_pattern.size();
 const auto decode = [ & ]( std::size_t i , std::vector< int > & x ) {
  for( std::size_t k = 0 ; k < np ; ++k ) {
   const auto r = std::size_t( task.range[ std::size_t( v_pattern[ k ] ) ] );
   x[ k ] = int( i % r );
   i /= r;
   }
  };
 const auto encode = [ & ]( const std::vector< int > & x ) {
  std::size_t i = 0;
  for( std::size_t k = 0 ; k < np ; ++k )
   i += std::size_t( x[ k ] ) * v_mult[ k ];
  return( i );
  };

 // the goal states, then the Dijkstra search backward
 v_dist.assign( size , dead_end );
 MinQueue q;
 std::vector< int > x( np ) , y( np );
 for( std::size_t i = 0 ; i < size ; ++i ) {
  decode( i , x );
  bool goal = true;
  for( const auto & g : task.goal )
   if( ( where[ std::size_t( g.var ) ] >= 0 ) &&
       ( x[ std::size_t( where[ std::size_t( g.var ) ] ) ] != g.val ) )
    goal = false;
  if( goal ) {
   v_dist[ i ] = 0;
   q.emplace( 0 , int( i ) );
   }
  }
 while( ! q.empty() ) {
  const auto [ d , j ] = q.top();
  q.pop();
  if( d > v_dist[ std::size_t( j ) ] )
   continue;
  decode( std::size_t( j ) , x );
  for( const AOp & a : ops ) {
   bool ok = true;
   for( const auto & [ k , post ] : a.eff )
    ok = ok && ( x[ std::size_t( k ) ] == post );
   if( ! ok )
    continue;
   // the predecessor: the values asked before, the others as in j
   y = x;
   for( const auto & [ k , val ] : a.pre )
    y[ std::size_t( k ) ] = val;
   for( const auto & [ k , val ] : a.pre ) {    // a prevail must hold in j
    bool changed = false;
    for( const auto & e : a.eff )
     changed = changed || ( e.first == k );
    ok = ok && ( changed || ( x[ std::size_t( k ) ] == val ) );
    }
   if( ! ok )
    continue;
   // each variable of the effects with no value before takes any value
   const std::size_t base = encode( y );
   std::vector< int > any = a.any;
   const auto nd = d + a.cost;
   std::size_t total = 1;
   for( int k : any )
    total *= std::size_t( task.range[ std::size_t( v_pattern[
                                       std::size_t( k ) ] ) ] );
   for( std::size_t c = 0 ; c < total ; ++c ) {
    std::size_t i = base , r = c;
    for( int k : any ) {
     const auto rk = std::size_t( task.range[ std::size_t( v_pattern[
                                               std::size_t( k ) ] ) ] );
     i -= std::size_t( y[ std::size_t( k ) ] ) * v_mult[ std::size_t( k ) ];
     i += ( r % rk ) * v_mult[ std::size_t( k ) ];
     r /= rk;
     }
    if( nd < v_dist[ i ] ) {
     v_dist[ i ] = nd;
     q.emplace( nd , int( i ) );
     }
    }
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::vector< int > SASPatternDatabase::greedy_pattern( const SASTask & task ,
                                                       std::size_t max_size )
{
 std::vector< char > in( task.num_vars() , 0 );
 std::vector< int > p;
 std::size_t size = 1;
 const auto fits = [ & ]( int v ) {
  return( ! in[ std::size_t( v ) ] &&
          ( size * std::size_t( task.range[ std::size_t( v ) ] ) <=
            max_size ) );
  };
 const auto add = [ & ]( int v ) {
  in[ std::size_t( v ) ] = 1;
  p.push_back( v );
  size *= std::size_t( task.range[ std::size_t( v ) ] );
  };
 for( const auto & g : task.goal )
  if( fits( g.var ) )
   add( g.var );
 for( bool grown = ! p.empty() ; grown ; ) {
  grown = false;
  std::vector< char > parent( task.num_vars() , 0 );
  for( std::size_t o = 0 ; o < task.ops.size() ; ++o ) {
   bool touches = false;
   for( const auto & e : task.ops[ o ].effects )
    touches = touches || in[ std::size_t( e.var ) ];
   if( ! touches )
    continue;
   for( const auto & f : task.preconditions( o ) )
    parent[ std::size_t( f.var ) ] = 1;
   }
  for( std::size_t v = 0 ; v < task.num_vars() ; ++v )
   if( parent[ v ] && fits( int( v ) ) ) {
    add( int( v ) );
    grown = true;
    break;
    }
  }
 return( p );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS SASMergeAndShrink -------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/// a transition system: the transitions of each label (operator), or all
/// the self-loops if the label leaves the system unchanged
struct TransitionSystem {
 int n = 0;                                          ///< the states
 int init = 0;                                       ///< the initial one
 std::vector< char > goal;                           ///< the goal states
 std::vector< std::vector< std::pair< int , int > > > tr;  ///< by label
 std::vector< char > irrelevant;                     ///< by label
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the distances of the states of x to its goal states, by a Dijkstra
/// search backward, dead_end for those that reach none

std::vector< int > goal_distances( const TransitionSystem & x ,
                                   const std::vector< int > & cost )
{
 std::vector< std::vector< std::pair< int , int > > > pred(
                                                   std::size_t( x.n ) );
 for( std::size_t o = 0 ; o < x.tr.size() ; ++o )
  if( ! x.irrelevant[ o ] )
   for( const auto & [ a , b ] : x.tr[ o ] )
    if( a != b )
     pred[ std::size_t( b ) ].emplace_back( a , cost[ o ] );
 std::vector< int > d( std::size_t( x.n ) , dead_end );
 MinQueue q;
 for( int s = 0 ; s < x.n ; ++s )
  if( x.goal[ std::size_t( s ) ] ) {
   d[ std::size_t( s ) ] = 0;
   q.emplace( 0 , s );
   }
 while( ! q.empty() ) {
  const auto [ c , s ] = q.top();
  q.pop();
  if( c > d[ std::size_t( s ) ] )
   continue;
  for( const auto & [ a , w ] : pred[ std::size_t( s ) ] )
   if( c + w < d[ std::size_t( a ) ] ) {
    d[ std::size_t( a ) ] = c + w;
    q.emplace( c + w , a );
    }
  }
 return( d );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// x with its states renamed by m (-1 dropping a state) into nb states

TransitionSystem rename( const TransitionSystem & x ,
                         const std::vector< int > & m , int nb )
{
 TransitionSystem y;
 y.n = nb;
 y.init = x.init >= 0 ? m[ std::size_t( x.init ) ] : -1;
 y.goal.assign( std::size_t( nb ) , 0 );
 for( int s = 0 ; s < x.n ; ++s )
  if( ( m[ std::size_t( s ) ] >= 0 ) && x.goal[ std::size_t( s ) ] )
   y.goal[ std::size_t( m[ std::size_t( s ) ] ) ] = 1;
 y.irrelevant = x.irrelevant;
 y.tr.resize( x.tr.size() );
 for( std::size_t o = 0 ; o < x.tr.size() ; ++o ) {
  for( const auto & [ a , b ] : x.tr[ o ] )
   if( ( m[ std::size_t( a ) ] >= 0 ) && ( m[ std::size_t( b ) ] >= 0 ) )
    y.tr[ o ].emplace_back( m[ std::size_t( a ) ] , m[ std::size_t( b ) ] );
  std::sort( y.tr[ o ].begin() , y.tr[ o ].end() );
  y.tr[ o ].erase( std::unique( y.tr[ o ].begin() , y.tr[ o ].end() ) ,
                   y.tr[ o ].end() );
  }
 return( y );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the renaming that keeps the states of x reachable from its initial
/// state and reaching a goal; returns the number of states kept

int prune_map( const TransitionSystem & x , const std::vector< int > & cost ,
               std::vector< int > & m )
{
 std::vector< char > seen( std::size_t( x.n ) , 0 );
 if( x.init >= 0 ) {
  std::vector< std::vector< int > > succ( std::size_t( x.n ) );
  for( std::size_t o = 0 ; o < x.tr.size() ; ++o )
   if( ! x.irrelevant[ o ] )
    for( const auto & [ a , b ] : x.tr[ o ] )
     succ[ std::size_t( a ) ].push_back( b );
  std::vector< int > stack{ x.init };
  seen[ std::size_t( x.init ) ] = 1;
  while( ! stack.empty() ) {
   const int a = stack.back();
   stack.pop_back();
   for( int b : succ[ std::size_t( a ) ] )
    if( ! seen[ std::size_t( b ) ] ) {
     seen[ std::size_t( b ) ] = 1;
     stack.push_back( b );
     }
   }
  }
 const auto d = goal_distances( x , cost );
 m.assign( std::size_t( x.n ) , -1 );
 int k = 0;
 for( int s = 0 ; s < x.n ; ++s )
  if( seen[ std::size_t( s ) ] && ( d[ std::size_t( s ) ] != dead_end ) )
   m[ std::size_t( s ) ] = k++;
 return( k );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the greedy bisimulation of x with at most limit blocks, as a renaming;
/// every state of x reaches a goal; returns the number of blocks

int shrink_map( const TransitionSystem & x , const std::vector< int > & cost ,
                int limit , std::vector< int > & m )
{
 const auto h = goal_distances( x , cost );
 // the blocks of equal h, the largest values merged beyond the limit
 std::vector< int > values( h );
 std::sort( values.begin() , values.end() );
 values.erase( std::unique( values.begin() , values.end() ) ,
               values.end() );
 m.resize( std::size_t( x.n ) );
 for( int s = 0 ; s < x.n ; ++s )
  m[ std::size_t( s ) ] = std::min( int( std::lower_bound(
                                     values.begin() , values.end() ,
                                     h[ std::size_t( s ) ] ) -
                                    values.begin() ) , limit - 1 );
 int nb = std::min( int( values.size() ) , limit );
 if( int( values.size() ) >= limit )
  return( nb );

 // the transitions that keep h, by source
 std::vector< std::vector< std::pair< int , int > > > out(
                                                    std::size_t( x.n ) );
 for( std::size_t o = 0 ; o < x.tr.size() ; ++o )
  if( ! x.irrelevant[ o ] )
   for( const auto & [ a , b ] : x.tr[ o ] )
    if( h[ std::size_t( a ) ] == cost[ o ] + h[ std::size_t( b ) ] )
     out[ std::size_t( a ) ].emplace_back( int( o ) , b );

 for( ; ; ) {
  std::map< std::vector< int > , int > sig;
  std::vector< int > next( std::size_t( x.n ) );
  for( int s = 0 ; s < x.n ; ++s ) {
   std::vector< std::pair< int , int > > l;
   for( const auto & [ o , b ] : out[ std::size_t( s ) ] )
    l.emplace_back( o , m[ std::size_t( b ) ] );
   std::sort( l.begin() , l.end() );
   l.erase( std::unique( l.begin() , l.end() ) , l.end() );
   std::vector< int > key{ m[ std::size_t( s ) ] };
   for( const auto & [ o , b ] : l ) {
    key.push_back( o );
    key.push_back( b );
    }
   next[ std::size_t( s ) ] = sig.emplace( std::move( key ) ,
                                           int( sig.size() ) ).first->second;
   }
  const int k = int( sig.size() );
  if( ( k == nb ) || ( k > limit ) )
   return( nb );
  m = std::move( next );
  nb = k;
  }
 }

}  // end( anonymous namespace )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

SASMergeAndShrink::SASMergeAndShrink( const SASTask & task ,
                                      std::size_t max_states )
{
 const std::size_t n = task.num_vars() , no = task.ops.size();
 const int bound = int( std::clamp< std::size_t >( max_states , 1 ,
                                                   std::size_t( 1 ) << 24 ) );
 std::vector< int > cost( no );
 std::vector< std::vector< SASTask::Fact > > pre( no );
 for( std::size_t o = 0 ; o < no ; ++o ) {
  cost[ o ] = task.ops[ o ].cost;
  pre[ o ] = task.preconditions( o );
  }

 // the causal graph: the parents of each variable
 std::vector< std::vector< char > > parent( n , std::vector< char >( n , 0 ) );
 for( std::size_t o = 0 ; o < no ; ++o )
  for( const auto & e : task.ops[ o ].effects ) {
   for( const auto & p : pre[ o ] )
    parent[ std::size_t( e.var ) ][ std::size_t( p.var ) ] = 1;
   for( const auto & f : task.ops[ o ].effects )
    parent[ std::size_t( e.var ) ][ std::size_t( f.var ) ] = 1;
   }
 // the ancestors of the goal, then the order of the linear merge
 std::vector< char > relevant( n , 0 ) , taken( n , 0 );
 std::vector< int > stack;
 for( const auto & g : task.goal )
  if( ! relevant[ std::size_t( g.var ) ] ) {
   relevant[ std::size_t( g.var ) ] = 1;
   stack.push_back( g.var );
   }
 while( ! stack.empty() ) {
  const auto v = std::size_t( stack.back() );
  stack.pop_back();
  for( std::size_t u = 0 ; u < n ; ++u )
   if( parent[ v ][ u ] && ! relevant[ u ] ) {
    relevant[ u ] = 1;
    stack.push_back( int( u ) );
    }
  }
 std::size_t next_goal = 0;
 for( ; ; ) {
  int pick = -1;
  for( std::size_t u = 0 ; ( u < n ) && ( pick < 0 ) ; ++u )
   if( relevant[ u ] && ! taken[ u ] )
    for( int v : v_vars )
     if( parent[ std::size_t( v ) ][ u ] ) {
      pick = int( u );
      break;
      }
  while( ( pick < 0 ) && ( next_goal < task.goal.size() ) ) {
   const int g = task.goal[ next_goal++ ].var;
   if( ! taken[ std::size_t( g ) ] )
    pick = g;
   }
  if( pick < 0 )
   break;
  taken[ std::size_t( pick ) ] = 1;
  v_vars.push_back( pick );
  v_range.push_back( task.range[ std::size_t( pick ) ] );
  }

 // the linear merge, the current system x
 TransitionSystem x;
 x.n = 1;
 x.goal.assign( 1 , 1 );
 x.tr.assign( no , {} );
 x.irrelevant.assign( no , 1 );
 for( std::size_t k = 0 ; k < v_vars.size() ; ++k ) {
  const int v = v_vars[ k ] , r = v_range[ k ];
  std::vector< int > m;
  // shrink x so that the product has at most bound states
  if( ( k > 0 ) && ( x.n > std::max( 1 , bound / r ) ) ) {
   const int nb = shrink_map( x , cost , std::max( 1 , bound / r ) , m );
   for( int & e : v_table.back() )
    if( e >= 0 )
     e = m[ std::size_t( e ) ];
   x = rename( x , m , nb );
   }
  // the atomic system of v
  TransitionSystem y;
  y.n = r;
  y.init = task.init[ std::size_t( v ) ];
  y.goal.assign( std::size_t( r ) , 1 );
  for( const auto & g : task.goal )
   if( g.var == v )
    for( int d = 0 ; d < r ; ++d )
     y.goal[ std::size_t( d ) ] = d == g.val;
  y.tr.assign( no , {} );
  y.irrelevant.assign( no , 1 );
  for( std::size_t o = 0 ; o < no ; ++o ) {
   for( const auto & e : task.ops[ o ].effects )
    if( e.var == v ) {
     y.irrelevant[ o ] = 0;
     for( int d = 0 ; d < r ; ++d )
      if( ( e.pre < 0 ) || ( e.pre == d ) )
       y.tr[ o ].emplace_back( d , e.post );
     }
   if( y.irrelevant[ o ] )
    for( const auto & p : task.ops[ o ].prevail )
     if( p.var == v ) {
      y.irrelevant[ o ] = 0;
      y.tr[ o ].emplace_back( p.val , p.val );
      }
   }
  // the synchronized product, the state ( a , d ) numbered a r + d
  TransitionSystem z;
  z.n = x.n * r;
  z.init = x.init * r + y.init;
  z.goal.resize( std::size_t( z.n ) );
  for( int a = 0 ; a < x.n ; ++a )
   for( int d = 0 ; d < r ; ++d )
    z.goal[ std::size_t( a * r + d ) ] = x.goal[ std::size_t( a ) ] &&
                                         y.goal[ std::size_t( d ) ];
  z.tr.assign( no , {} );
  z.irrelevant.assign( no , 0 );
  for( std::size_t o = 0 ; o < no ; ++o ) {
   if( x.irrelevant[ o ] && y.irrelevant[ o ] ) {
    z.irrelevant[ o ] = 1;
    continue;
    }
   std::vector< std::pair< int , int > > lx = x.tr[ o ] , ly = y.tr[ o ];
   if( x.irrelevant[ o ] )
    for( int a = 0 ; a < x.n ; ++a )
     lx.emplace_back( a , a );
   if( y.irrelevant[ o ] )
    for( int d = 0 ; d < r ; ++d )
     ly.emplace_back( d , d );
   for( const auto & [ a , b ] : lx )
    for( const auto & [ c , d ] : ly )
     z.tr[ o ].emplace_back( a * r + c , b * r + d );
   }
  // the table of the step, then the pruning
  std::vector< int > t( std::size_t( x.n ) * std::size_t( r ) );
  for( std::size_t i = 0 ; i < t.size() ; ++i )
   t[ i ] = int( i );
  const int kept = prune_map( z , cost , m );
  for( int & e : t )
   e = m[ std::size_t( e ) ];
  v_table.push_back( std::move( t ) );
  x = rename( z , m , kept );
  if( kept == 0 ) {                      // the initial state is a dead end
   v_vars.resize( v_table.size() );
   v_range.resize( v_table.size() );
   break;
   }
  }
 v_dist = x.n > 0 ? goal_distances( x , cost ) : std::vector< int >();
 if( v_vars.empty() )
  v_dist.assign( 1 , 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS SASCartesianAbstraction ----------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/// the transition system of a Cartesian abstraction under refinement, with
/// its refinement hierarchy

class Refiner {

 public:

 /// a transition, from or to the other abstract state
 struct Arc {
  int op;
  int other;
  };

 /// the condition of an operator on a variable: the value asked, or -1,
 /// and the value given, or -1 if the variable is not changed
 struct Cond {
  int var;
  int pre;
  int post;
  };

 Refiner( const SASTask & task , std::span< const SASTask::Fact > goal ,
          std::span< const int > cost )
  : t( task ) , g( goal ) , c( cost ) , n( task.num_vars() )
 {
  offset.resize( n );
  int f = 0;
  for( std::size_t v = 0 ; v < n ; ++v ) {
   offset[ v ] = f;
   f += task.range[ v ];
   }
  words = ( std::size_t( f ) + 63 ) / 64;
  cond.resize( task.ops.size() );
  pre.resize( task.ops.size() );
  for( std::size_t o = 0 ; o < task.ops.size() ; ++o ) {
   for( const auto & p : task.ops[ o ].prevail )
    cond[ o ].push_back( { p.var , p.val , -1 } );
   for( const auto & e : task.ops[ o ].effects )
    cond[ o ].push_back( { e.var , e.pre , e.post } );
   std::sort( cond[ o ].begin() , cond[ o ].end() ,
              []( const Cond & x , const Cond & y ) {
               return( x.var < y.var );
               } );
   pre[ o ] = task.preconditions( o );
   }
  // the single abstract state, with every operator as a self-loop
  bits.assign( words , 0 );
  for( int i = 0 ; i < f ; ++i )
   bits[ std::size_t( i ) / 64 ] |= std::uint64_t( 1 ) << ( i % 64 );
  out.resize( 1 );
  in.resize( 1 );
  loops.resize( 1 );
  for( std::size_t o = 0 ; o < task.ops.size() ; ++o )
   loops[ 0 ].push_back( int( o ) );
  var.push_back( -1 );
  first.push_back( 0 );
  left.push_back( 0 );
  right.push_back( 0 );
  leaf.push_back( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int states( void ) const { return( int( out.size() ) ); }

 [[nodiscard]] bool has( int a , int v , int d ) const {
  const auto i = std::size_t( offset[ std::size_t( v ) ] + d );
  return( ( bits[ std::size_t( a ) * words + i / 64 ] >> ( i % 64 ) ) & 1 );
  }

 void set( int a , int v , int d , bool x ) {
  const auto i = std::size_t( offset[ std::size_t( v ) ] + d );
  auto & w = bits[ std::size_t( a ) * words + i / 64 ];
  const auto m = std::uint64_t( 1 ) << ( i % 64 );
  w = x ? ( w | m ) : ( w & ~m );
  }

 [[nodiscard]] int count( int a , int v ) const {
  int k = 0;
  for( int d = 0 ; d < t.range[ std::size_t( v ) ] ; ++d )
   k += has( a , v , d );
  return( k );
  }

 [[nodiscard]] bool is_goal( int a ) const {
  for( const auto & f : g )
   if( ! has( a , f.var , f.val ) )
    return( false );
  return( true );
  }

 /// the abstract state of the values s
 [[nodiscard]] int state_of( const std::vector< int > & s ) const {
  int k = 0;
  while( var[ std::size_t( k ) ] >= 0 )
   k = split[ std::size_t( first[ std::size_t( k ) ] +
                           s[ std::size_t( var[ std::size_t( k ) ] ) ] ) ] ?
       right[ std::size_t( k ) ] : left[ std::size_t( k ) ];
  return( left[ std::size_t( k ) ] );
  }

 /// true if, as far as the variable v is concerned, o leads from a state
 /// of a to one of b
 [[nodiscard]] bool allows( int o , int v , int a , int b ) const {
  for( const auto & x : cond[ std::size_t( o ) ] )
   if( x.var == v ) {
    if( ( x.pre >= 0 ) && ! has( a , v , x.pre ) )
     return( false );
    return( has( b , v , x.post >= 0 ? x.post : x.pre ) );
    }
  for( int d = 0 ; d < t.range[ std::size_t( v ) ] ; ++d )
   if( has( a , v , d ) && has( b , v , d ) )
    return( true );
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// splits a on v: the values in wanted go to a new abstract state

 void refine( int a , int v , const std::vector< char > & wanted ) {
  const int b = states();
  bits.resize( bits.size() + words );
  std::copy_n( bits.begin() + std::ptrdiff_t( std::size_t( a ) * words ) ,
               words , bits.begin() + std::ptrdiff_t( std::size_t( b ) *
                                                      words ) );
  for( int d = 0 ; d < t.range[ std::size_t( v ) ] ; ++d ) {
   const bool x = has( a , v , d );
   set( b , v , d , x && wanted[ std::size_t( d ) ] );
   set( a , v , d , x && ! wanted[ std::size_t( d ) ] );
   }
  out.emplace_back();
  in.emplace_back();
  loops.emplace_back();
  mark.resize( out.size() , 0 );

  // the hierarchy: the leaf of a becomes the node of the split
  const auto k = std::size_t( leaf[ std::size_t( a ) ] );
  var[ k ] = v;
  first[ k ] = int( split.size() );
  for( int d = 0 ; d < t.range[ std::size_t( v ) ] ; ++d )
   split.push_back( wanted[ std::size_t( d ) ] );
  for( int x : { a , b } ) {
   ( x == a ? left[ k ] : right[ k ] ) = int( var.size() );
   leaf.resize( out.size() );
   leaf[ std::size_t( x ) ] = int( var.size() );
   var.push_back( -1 );
   first.push_back( 0 );
   left.push_back( x );
   right.push_back( 0 );
   }

  // the transitions
  auto old_in = std::move( in[ std::size_t( a ) ] );
  auto old_out = std::move( out[ std::size_t( a ) ] );
  auto old_loops = std::move( loops[ std::size_t( a ) ] );
  in[ std::size_t( a ) ].clear();
  out[ std::size_t( a ) ].clear();
  loops[ std::size_t( a ) ].clear();
  arcs -= old_in.size() + old_out.size();
  const auto drop = [ & ]( std::vector< Arc > & l ) {
   std::erase_if( l , [ & ]( const Arc & x ) { return( x.other == a ); } );
   };
  const auto add = [ & ]( int x , int o , int y ) {
   out[ std::size_t( x ) ].push_back( { o , y } );
   in[ std::size_t( y ) ].push_back( { o , x } );
   ++arcs;
   };
  ++stamp;
  for( const auto & x : old_in )
   if( mark[ std::size_t( x.other ) ] != stamp ) {
    mark[ std::size_t( x.other ) ] = stamp;
    drop( out[ std::size_t( x.other ) ] );
    }
  for( const auto & x : old_in )
   for( int y : { a , b } )
    if( allows( x.op , v , x.other , y ) )
     add( x.other , x.op , y );
  ++stamp;
  for( const auto & x : old_out )
   if( mark[ std::size_t( x.other ) ] != stamp ) {
    mark[ std::size_t( x.other ) ] = stamp;
    drop( in[ std::size_t( x.other ) ] );
    }
  for( const auto & x : old_out )
   for( int y : { a , b } )
    if( allows( x.op , v , y , x.other ) )
     add( y , x.op , x.other );
  for( int o : old_loops )
   for( int x : { a , b } )
    for( int y : { a , b } )
     if( allows( o , v , x , y ) ) {
      if( x == y )
       loops[ std::size_t( x ) ].push_back( o );
      else
       add( x , o , y );
      }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// an optimal abstract plan from a, as the pairs ( operator , abstract
 /// state reached ), by A* with the heuristic h, which is consistent;
 /// false if a reaches no abstract goal

 bool plan( int a , const std::vector< int > & h ,
            std::vector< Arc > & p ) const {
  if( h[ std::size_t( a ) ] == dead_end )
   return( false );
  std::vector< int > d( out.size() , dead_end ) , op( out.size() , -1 ) ,
   from( out.size() , -1 );
  MinQueue q;
  d[ std::size_t( a ) ] = 0;
  q.emplace( h[ std::size_t( a ) ] , a );
  while( ! q.empty() ) {
   const auto [ fx , x ] = q.top();
   q.pop();
   const int dx = d[ std::size_t( x ) ];
   if( fx > dx + h[ std::size_t( x ) ] )
    continue;
   if( is_goal( x ) ) {
    p.clear();
    for( int y = x ; y != a ; y = from[ std::size_t( y ) ] )
     p.push_back( { op[ std::size_t( y ) ] , y } );
    std::reverse( p.begin() , p.end() );
    return( true );
    }
   for( const auto & e : out[ std::size_t( x ) ] ) {
    const int dy = dx + c[ std::size_t( e.op ) ];
    const int hy = h[ std::size_t( e.other ) ];
    if( ( hy != dead_end ) && ( dy < d[ std::size_t( e.other ) ] ) ) {
     d[ std::size_t( e.other ) ] = dy;
     op[ std::size_t( e.other ) ] = e.op;
     from[ std::size_t( e.other ) ] = x;
     q.emplace( dy + hy , e.other );
     }
    }
   }
  return( false );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the distances of the abstract states to the abstract goals

 [[nodiscard]] std::vector< int > goal_distances( void ) const {
  std::vector< int > d( out.size() , dead_end );
  MinQueue q;
  for( int x = 0 ; x < states() ; ++x )
   if( is_goal( x ) ) {
    d[ std::size_t( x ) ] = 0;
    q.emplace( 0 , x );
    }
  while( ! q.empty() ) {
   const auto [ dx , x ] = q.top();
   q.pop();
   if( dx > d[ std::size_t( x ) ] )
    continue;
   for( const auto & e : in[ std::size_t( x ) ] ) {
    const int dy = dx + c[ std::size_t( e.op ) ];
    if( dy < d[ std::size_t( e.other ) ] ) {
     d[ std::size_t( e.other ) ] = dy;
     q.emplace( dy , e.other );
     }
    }
   }
  return( d );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// one step of the refinement: traces the abstract plan p on the task and
 /// splits on its first flaw; returns the abstract state split, or -1 if p is
 /// a plan of the task

 int split_flaw( const std::vector< Arc > & p ) {
  std::vector< int > s = t.init;
  int a = state_of( s );
  // the variable of the candidates ( v , wanted ) whose abstract domain is
  // the smallest fraction of its domain, the least index on a tie
  std::vector< std::pair< int , std::vector< char > > > cand;
  const auto split_best = [ & ]( void ) {
   std::size_t best = 0;
   for( std::size_t i = 1 ; i < cand.size() ; ++i ) {
    const int u = cand[ i ].first , w = cand[ best ].first;
    const long lu = long( count( a , u ) ) * t.range[ std::size_t( w ) ] ,
               lw = long( count( a , w ) ) * t.range[ std::size_t( u ) ];
    if( ( lu < lw ) || ( ( lu == lw ) && ( u < w ) ) )
     best = i;
    }
   refine( a , cand[ best ].first , cand[ best ].second );
   return( a );
   };
  const auto only = [ & ]( int v , int d ) {
   std::vector< char > w( std::size_t( t.range[ std::size_t( v ) ] ) , 0 );
   w[ std::size_t( d ) ] = 1;
   return( w );
   };
  for( const auto & [ o , b ] : p ) {
   for( const auto & f : pre[ std::size_t( o ) ] )      // precondition
    if( s[ std::size_t( f.var ) ] != f.val )
     cand.emplace_back( f.var , only( f.var , f.val ) );
   if( ! cand.empty() )
    return( split_best() );
   std::vector< int > next = s;
   for( const auto & e : t.ops[ std::size_t( o ) ].effects )
    next[ std::size_t( e.var ) ] = e.post;
   for( std::size_t v = 0 ; v < n ; ++v )                // deviation
    if( ! has( b , int( v ) , next[ v ] ) ) {
     std::vector< char > w( std::size_t( t.range[ v ] ) );
     for( int d = 0 ; d < t.range[ v ] ; ++d )
      w[ std::size_t( d ) ] = has( a , int( v ) , d ) &&
                              has( b , int( v ) , d );
     cand.emplace_back( int( v ) , std::move( w ) );
     }
   if( ! cand.empty() )
    return( split_best() );
   s = std::move( next );
   a = b;
   }
  for( const auto & f : g )                              // goal
   if( s[ std::size_t( f.var ) ] != f.val )
    cand.emplace_back( f.var , only( f.var , f.val ) );
  return( cand.empty() ? -1 : split_best() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const SASTask & t;                              ///< the task
 std::span< const SASTask::Fact > g;             ///< the goal
 std::span< const int > c;                       ///< the costs
 std::size_t n;                                  ///< the variables
 std::vector< int > offset;                      ///< the first bit of v
 std::size_t words;                              ///< the words of a state
 std::vector< std::vector< Cond > > cond;        ///< by operator and var
 std::vector< std::vector< SASTask::Fact > > pre; ///< the preconditions
 std::vector< std::uint64_t > bits;              ///< the values of each a
 std::vector< std::vector< Arc > > out;          ///< the transitions out
 std::vector< std::vector< Arc > > in;           ///< the transitions in
 std::vector< std::vector< int > > loops;        ///< the self-loops
 std::size_t arcs = 0;                           ///< the transitions
 std::vector< int > mark;                        ///< the stamps of a scan
 int stamp = 0;                                  ///< the current one
 std::vector< int > var , first , left , right;  ///< the hierarchy
 std::vector< char > split;                      ///< the values split
 std::vector< int > leaf;                        ///< the leaf of each a

 };  // end( class( Refiner ) )

}  // end( anonymous namespace )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

SASCartesianAbstraction::SASCartesianAbstraction( const SASTask & task ,
                                                  const Limits & limits )
{
 std::vector< int > cost( task.ops.size() );
 for( std::size_t o = 0 ; o < cost.size() ; ++o )
  cost[ o ] = task.ops[ o ].cost;
 build( task , task.goal , cost , limits );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

SASCartesianAbstraction::SASCartesianAbstraction( const SASTask & task ,
                                std::span< const SASTask::Fact > goal ,
                                std::span< const int > cost ,
                                const Limits & limits )
{
 build( task , goal , cost , limits );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void SASCartesianAbstraction::build( const SASTask & task ,
                                     std::span< const SASTask::Fact > goal ,
                                     std::span< const int > cost ,
                                     const Limits & limits )
{
 if( cost.size() != task.ops.size() )
  throw( std::invalid_argument(
          "SASCartesianAbstraction::build: one cost per operator" ) );
 for( const auto & f : goal )
  if( ( f.var < 0 ) || ( std::size_t( f.var ) >= task.num_vars() ) ||
      ( f.val < 0 ) || ( f.val >= task.range[ std::size_t( f.var ) ] ) )
   throw( std::invalid_argument(
           "SASCartesianAbstraction::build: a goal fact out of range" ) );
 using Clock = std::chrono::steady_clock;
 const auto start = Clock::now();
 Refiner r( task , goal , cost );
 std::vector< Refiner::Arc > p;
 // the goal distances of an earlier abstraction, the new state of a split
 // taking those of the old one: the distances only grow with the splits,
 // and a transition after a split comes from one before it, whence they
 // stay admissible and consistent; they are computed again when the
 // states have grown by a quarter
 std::vector< int > h = r.goal_distances();
 std::size_t fresh = 1;
 for( ; ; ) {
  if( ( std::size_t( r.states() ) >= limits.max_states ) ||
      ( r.arcs >= limits.max_transitions ) ||
      ( ( limits.max_time > 0 ) && ( std::chrono::duration< double >(
           Clock::now() - start ).count() >= limits.max_time ) ) )
   break;
  if( std::size_t( r.states() ) * 4 >= fresh * 5 ) {
   h = r.goal_distances();
   fresh = std::size_t( r.states() );
   }
  if( ! r.plan( r.state_of( task.init ) , h , p ) )
   break;                                // the initial state is a dead end
  const int a = r.split_flaw( p );
  if( a < 0 ) {
   f_exact = true;
   break;
   }
  h.push_back( h[ std::size_t( a ) ] );
  }

 v_dist = r.goal_distances();
 v_saturated.assign( task.ops.size() , 0 );
 for( int a = 0 ; a < r.states() ; ++a )
  if( v_dist[ std::size_t( a ) ] != dead_end )
   for( const auto & e : r.out[ std::size_t( a ) ] )
    if( v_dist[ std::size_t( e.other ) ] != dead_end )
     v_saturated[ std::size_t( e.op ) ] = std::max(
      v_saturated[ std::size_t( e.op ) ] ,
      v_dist[ std::size_t( a ) ] - v_dist[ std::size_t( e.other ) ] );
 f_transitions = r.arcs;
 v_split = std::move( r.split );
 v_node.resize( r.var.size() );
 for( std::size_t k = 0 ; k < v_node.size() ; ++k )
  v_node[ k ] = { r.var[ k ] , r.first[ k ] , r.left[ k ] , r.right[ k ] };
 }

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS SASAdditiveCartesian -----------------------*/
/*--------------------------------------------------------------------------*/

SASAdditiveCartesian::SASAdditiveCartesian( const SASTask & task ,
                            const SASCartesianAbstraction::Limits & limits )
{
 using Clock = std::chrono::steady_clock;
 const auto start = Clock::now();
 std::vector< int > cost( task.ops.size() );
 for( std::size_t o = 0 ; o < cost.size() ; ++o )
  cost[ o ] = task.ops[ o ].cost;
 const std::size_t k = task.goal.size();
 std::size_t states = limits.max_states;
 v_abs.reserve( k );
 for( std::size_t i = 0 ; i < k ; ++i ) {
  auto l = limits;
  l.max_states = std::max< std::size_t >( 1 , states / ( k - i ) );
  if( limits.max_time > 0 ) {
   const double left = limits.max_time - std::chrono::duration< double >(
                                          Clock::now() - start ).count();
   l.max_time = std::max( left , 1e-9 ) / double( k - i );
   }
  v_abs.emplace_back( task , std::span( task.goal ).subspan( i , 1 ) ,
                      cost , l );
  const auto & a = v_abs.back();
  states -= std::min( states , a.size() );
  for( std::size_t o = 0 ; o < cost.size() ; ++o )
   cost[ o ] -= a.saturated_costs()[ o ];
  if( a( task.init ) == dead_end )
   break;
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File SASHeuristics.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
