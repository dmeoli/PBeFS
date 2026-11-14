/*--------------------------------------------------------------------------*/
/*------------------------- File SASHeuristics.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the heuristics of classical planning: the classes
 * RelaxedHeuristic, LMCut and SASPatternDatabase.
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

#include <functional>

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
/*---------------------- End File SASHeuristics.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
