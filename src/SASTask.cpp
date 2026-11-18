/*--------------------------------------------------------------------------*/
/*---------------------------- File SASTask.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the classes SASTask and StripsTask: the reading of a
 * file output.sas and the STRIPS view of the task.
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

#include <fstream>

#include <istream>

#include <sstream>

#include <stdexcept>

#include "PBeFS/astar/SASTask.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- THE FILE READER -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

// the lines of a file output.sas, each read whole, the blank ones skipped

class Lines {

 public:

 explicit Lines( std::istream & in ) : f_in( in ) {}

 // the next line, without the spaces at its ends

 std::string next( void ) {
  std::string l;
  while( std::getline( f_in , l ) ) {
   ++f_line;
   const auto a = l.find_first_not_of( " \t\r" );
   if( a == std::string::npos )
    continue;
   const auto b = l.find_last_not_of( " \t\r" );
   return( l.substr( a , b - a + 1 ) );
   }
  fail( "unexpected end of file" );
  return( l );
  }

 // the next line, which must be w

 void expect( const std::string & w ) {
  if( next() != w )
   fail( "expected " + w );
  }

 // the next line, which must be a number in [ lo , hi ]

 int number( int lo , int hi ) {
  std::istringstream is( next() );
  int x;
  if( ! ( is >> x ) || ( x < lo ) || ( x > hi ) )
   fail( "a number out of range" );
  return( x );
  }

 // the integers of the next line

 std::vector< int > numbers( void ) {
  std::istringstream is( next() );
  std::vector< int > v;
  for( int x ; is >> x ; )
   v.push_back( x );
  return( v );
  }

 [[noreturn]] void fail( const std::string & what ) const {
  throw( std::invalid_argument( "SASTask::read: line " +
                                std::to_string( f_line ) + ": " + what ) );
  }

 private:

 std::istream & f_in;
 int f_line = 0;

 };  // end( class( Lines ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the fact var val of the next line, checked against the task

SASTask::Fact read_fact( Lines & in , const SASTask & t )
{
 const auto v = in.numbers();
 if( ( v.size() != 2 ) || ( v[ 0 ] < 0 ) ||
     ( std::size_t( v[ 0 ] ) >= t.num_vars() ) || ( v[ 1 ] < 0 ) ||
     ( v[ 1 ] >= t.range[ std::size_t( v[ 0 ] ) ] ) )
  in.fail( "a fact out of range" );
 return( SASTask::Fact{ v[ 0 ] , v[ 1 ] } );
 }

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*---------------------------- READING A TASK ------------------------------*/
/*--------------------------------------------------------------------------*/

SASTask SASTask::read( std::istream & in )
{
 Lines l( in );
 SASTask t;
 l.expect( "begin_version" );
 if( l.number( 0 , 1000 ) != 3 )
  l.fail( "only the version 3 of the format is read" );
 l.expect( "end_version" );
 l.expect( "begin_metric" );
 t.use_costs = l.number( 0 , 1 ) == 1;
 l.expect( "end_metric" );

 const int n = l.number( 1 , 1 << 24 );
 for( int v = 0 ; v < n ; ++v ) {
  l.expect( "begin_variable" );
  t.var_names.push_back( l.next() );
  if( l.number( -1 , 1 << 24 ) != -1 )
   l.fail( "the axioms are not supported" );
  const int r = l.number( 1 , 1 << 24 );
  t.range.push_back( r );
  t.values.emplace_back();
  for( int d = 0 ; d < r ; ++d )
   t.values.back().push_back( l.next() );
  l.expect( "end_variable" );
  }

 const int m = l.number( 0 , 1 << 30 );
 for( int g = 0 ; g < m ; ++g ) {
  l.expect( "begin_mutex_group" );
  const int k = l.number( 0 , 1 << 30 );
  t.mutexes.emplace_back();
  for( int i = 0 ; i < k ; ++i )
   t.mutexes.back().push_back( read_fact( l , t ) );
  l.expect( "end_mutex_group" );
  }

 l.expect( "begin_state" );
 for( int v = 0 ; v < n ; ++v )
  t.init.push_back( l.number( 0 , t.range[ std::size_t( v ) ] - 1 ) );
 l.expect( "end_state" );

 l.expect( "begin_goal" );
 const int ng = l.number( 0 , n );
 for( int i = 0 ; i < ng ; ++i )
  t.goal.push_back( read_fact( l , t ) );
 l.expect( "end_goal" );

 const int no = l.number( 0 , 1 << 30 );
 t.ops.reserve( std::size_t( no ) );
 for( int o = 0 ; o < no ; ++o ) {
  l.expect( "begin_operator" );
  Operator op;
  op.name = l.next();
  const int np = l.number( 0 , n );
  for( int i = 0 ; i < np ; ++i )
   op.prevail.push_back( read_fact( l , t ) );
  const int ne = l.number( 0 , 1 << 24 );
  for( int i = 0 ; i < ne ; ++i ) {
   const auto e = l.numbers();
   if( e.empty() || ( e[ 0 ] != 0 ) )
    l.fail( "the conditional effects are not supported" );
   if( ( e.size() != 4 ) || ( e[ 1 ] < 0 ) || ( e[ 1 ] >= n ) )
    l.fail( "a malformed effect" );
   const int r = t.range[ std::size_t( e[ 1 ] ) ];
   if( ( e[ 2 ] < -1 ) || ( e[ 2 ] >= r ) || ( e[ 3 ] < 0 ) ||
       ( e[ 3 ] >= r ) )
    l.fail( "a value out of range" );
   op.effects.push_back( { e[ 1 ] , e[ 2 ] , e[ 3 ] } );
   }
  const int c = l.number( -( 1 << 30 ) , 1 << 30 );
  if( c < 0 )
   l.fail( "a negative cost" );
  op.cost = t.use_costs ? c : 1;
  l.expect( "end_operator" );
  t.ops.push_back( std::move( op ) );
  }

 if( l.number( 0 , 1 << 30 ) != 0 )
  l.fail( "the axioms are not supported" );
 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

SASTask SASTask::read_file( const std::string & name )
{
 std::ifstream in( name );
 if( ! in )
  throw( std::invalid_argument( "SASTask::read_file: cannot read " +
                                name ) );
 return( read( in ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- THE PRECONDITIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

std::vector< SASTask::Fact > SASTask::preconditions( std::size_t o ) const
{
 std::vector< Fact > p = ops[ o ].prevail;
 for( const Effect & e : ops[ o ].effects )
  if( e.pre >= 0 )
   p.push_back( { e.var , e.pre } );
 std::sort( p.begin() , p.end() , []( const Fact & a , const Fact & b ) {
  return( ( a.var < b.var ) || ( ( a.var == b.var ) && ( a.val < b.val ) ) );
  } );
 p.erase( std::unique( p.begin() , p.end() ,
                       []( const Fact & a , const Fact & b ) {
                        return( ( a.var == b.var ) && ( a.val == b.val ) );
                        } ) , p.end() );
 return( p );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- THE STRIPS VIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

StripsTask::StripsTask( const SASTask & t )
{
 for( int r : t.range ) {
  offset.push_back( num_facts );
  num_facts += r;
  }
 for( std::size_t v = 0 ; v < t.num_vars() ; ++v )
  init.push_back( fact( int( v ) , t.init[ v ] ) );
 for( const auto & g : t.goal )
  goal.push_back( fact( g.var , g.val ) );
 std::sort( goal.begin() , goal.end() );
 ops.reserve( t.ops.size() );
 for( std::size_t o = 0 ; o < t.ops.size() ; ++o ) {
  Operator op;
  for( const auto & p : t.preconditions( o ) )
   op.pre.push_back( fact( p.var , p.val ) );
  for( const auto & e : t.ops[ o ].effects ) {
   op.add.push_back( fact( e.var , e.post ) );
   if( e.pre >= 0 ) {
    if( e.pre != e.post )
     op.del.push_back( fact( e.var , e.pre ) );
    }
   else
    for( int d = 0 ; d < t.range[ std::size_t( e.var ) ] ; ++d )
     if( d != e.post )
      op.del.push_back( fact( e.var , d ) );
   }
  std::sort( op.pre.begin() , op.pre.end() );
  std::sort( op.add.begin() , op.add.end() );
  std::sort( op.del.begin() , op.del.end() );
  op.cost = t.ops[ o ].cost;
  ops.push_back( std::move( op ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File SASTask.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
