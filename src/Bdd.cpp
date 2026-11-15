/*--------------------------------------------------------------------------*/
/*------------------------------ File Bdd.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the class Bdd and of BddRuntime on Sylvan and Lace.
 *
 * Sylvan drops at each garbage collection the nodes that are reachable
 * neither from a registered address nor from the stacks of its tasks: a Bdd
 * registers the address of its handle with mtbdd_protect() for as long as
 * it lives. A result returned by Sylvan to a thread outside Lace is not
 * registered until that thread stores it in a Bdd, and a collection asked
 * for meanwhile by another thread could drop it; hence image_into(), the
 * one operation that many threads may call at once, runs as a task of Lace
 * that writes the result in the registered handle of its output before it
 * ends, and the collection, which stops every worker at the start of a
 * task, never sees the result unregistered. The same holds for the tasks of
 * images() and disjoin(), which write each partial result in a registered
 * handle.
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

#include <atomic>

#include <stdexcept>

#include <type_traits>

#include <sylvan.h>

#include <sylvan_int.h>

#include "PBeFS/symbolic/Bdd.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

using namespace sylvan;

/*--------------------------------------------------------------------------*/
/*--------------------------------- TASKS ----------------------------------*/
/*--------------------------------------------------------------------------*/
// the tasks of Lace are defined at the global scope, as its macros need

namespace {

/// what the tasks of images() share
struct ImageJob {
 BDD set;                              // the set of states
 const std::uint64_t * rel;            // the relations
 const std::uint64_t * vars;           // their variables
 std::uint64_t * out;                  // the registered outputs
 int backward;                         // preimages if not 0
 };

std::atomic< bool > g_running{ false };  // a BddRuntime is alive

/// a call of a callable, run as a task of Lace
struct Thunk {
 void ( * call )( void * );            // calls obj
 void * obj;                           // the callable
 };

}  // end( namespace )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the call of the thunk

VOID_TASK_1( pbefs_call , const Thunk * , t )
{
 t->call( t->obj );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// f() run on a worker of Lace: the operations of Sylvan that are not tasks
// (the creation of a node among them) take the region of the table of the
// worker that runs them, and need one; f must not throw

namespace {

template< class F >
void in_lace( F && f )
{
 using G = std::remove_reference_t< F >;
 const Thunk t{ []( void * o ) { ( *static_cast< G * >( o ) )(); } ,
                static_cast< void * >( & f ) };
 RUN( pbefs_call , & t );
 }

}  // end( namespace )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the image (preimage) of set under rel, written in *out

VOID_TASK_5( pbefs_image_into , BDD , set , BDD , rel , BDD , vars , int ,
             backward , std::uint64_t * , out )
{
 *out = backward ? sylvan_relprev( rel , set , vars )
                 : sylvan_relnext( set , rel , vars );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the images of the relations lo .. hi - 1 of the job, the two halves of
// the range in parallel

VOID_TASK_3( pbefs_images , const ImageJob * , job , std::size_t , lo ,
             std::size_t , hi )
{
 if( hi - lo == 1 ) {
  job->out[ lo ] = job->backward
   ? sylvan_relprev( job->rel[ lo ] , job->set , job->vars[ lo ] )
   : sylvan_relnext( job->set , job->rel[ lo ] , job->vars[ lo ] );
  return;
  }
 const std::size_t mid = lo + ( hi - lo ) / 2;
 SPAWN( pbefs_images , job , lo , mid );
 CALL( pbefs_images , job , mid , hi );
 SYNC( pbefs_images );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the disjunction of v[ lo .. hi - 1 ], registered, written in v[ lo ]

VOID_TASK_3( pbefs_disjoin , std::uint64_t * , v , std::size_t , lo ,
             std::size_t , hi )
{
 if( hi - lo == 1 )
  return;
 const std::size_t mid = lo + ( hi - lo ) / 2;
 SPAWN( pbefs_disjoin , v , lo , mid );
 CALL( pbefs_disjoin , v , mid , hi );
 SYNC( pbefs_disjoin );
 v[ lo ] = sylvan_or( v[ lo ] , v[ mid ] );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- METHODS OF Bdd -----------------------------*/
/*--------------------------------------------------------------------------*/

Bdd::Bdd( void ) noexcept : f_bdd( sylvan_false )
{
 mtbdd_protect( & f_bdd );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd::Bdd( std::uint64_t h ) noexcept : f_bdd( h )
{
 mtbdd_protect( & f_bdd );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd::Bdd( const Bdd & b ) noexcept : f_bdd( b.f_bdd )
{
 mtbdd_protect( & f_bdd );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd & Bdd::operator=( const Bdd & b ) noexcept
{
 f_bdd = b.f_bdd;
 return( *this );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd::~Bdd()
{
 mtbdd_unprotect( & f_bdd );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::one( void )
{
 return( Bdd( sylvan_true ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::variable( std::uint32_t v )
{
 BDD r = sylvan_false;
 in_lace( [ & ]( void ) { r = sylvan_ithvar( v ); } );
 return( Bdd( r ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::cube( const std::vector< std::uint32_t > & vars ,
               const std::vector< bool > & values )
{
 if( vars.size() != values.size() )
  throw( std::invalid_argument( "Bdd::cube: as many values as variables" ) );
 std::vector< std::uint32_t > v( vars );
 if( ! std::is_sorted( v.begin() , v.end() ) )
  throw( std::invalid_argument( "Bdd::cube: the variables not sorted" ) );
 std::vector< std::uint8_t > c( values.size() );
 for( std::size_t i = 0 ; i < c.size() ; ++i )
  c[ i ] = values[ i ] ? 1 : 0;
 BDD r = sylvan_false;
 in_lace( [ & ]( void ) {
  const BDD set = mtbdd_set_from_array( v.data() , v.size() );
  mtbdd_refs_push( set );
  r = sylvan_cube( set , c.data() );
  mtbdd_refs_pop( 1 );
  } );
 return( Bdd( r ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::variables( std::vector< std::uint32_t > vars )
{
 std::sort( vars.begin() , vars.end() );
 vars.erase( std::unique( vars.begin() , vars.end() ) , vars.end() );
 BDD r = sylvan_false;
 in_lace( [ & ]( void ) {
  r = mtbdd_set_from_array( vars.data() , vars.size() );
  } );
 return( Bdd( r ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::operator&( const Bdd & b ) const
{
 return( Bdd( sylvan_and( f_bdd , b.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::operator|( const Bdd & b ) const
{
 return( Bdd( sylvan_or( f_bdd , b.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::operator-( const Bdd & b ) const
{
 return( Bdd( sylvan_diff( f_bdd , b.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::operator!( void ) const
{
 return( Bdd( sylvan_not( f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::equiv( const Bdd & b ) const
{
 return( Bdd( sylvan_equiv( f_bdd , b.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd Bdd::exists( const Bdd & vars ) const
{
 return( Bdd( sylvan_exists( f_bdd , vars.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool Bdd::empty( void ) const noexcept
{
 return( f_bdd == sylvan_false );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool Bdd::full( void ) const noexcept
{
 return( f_bdd == sylvan_true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool Bdd::disjoint( const Bdd & b ) const
{
 return( sylvan_disjoint( f_bdd , b.f_bdd ) != 0 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::size_t Bdd::nodes( void ) const
{
 return( mtbdd_nodecount( f_bdd ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double Bdd::count( const Bdd & vars ) const
{
 return( sylvan_satcount( f_bdd , vars.f_bdd ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::vector< bool > Bdd::pick( const std::vector< std::uint32_t > & vars )
 const
{
 if( f_bdd == sylvan_false )
  return( std::vector< bool >() );
 std::vector< std::uint32_t > v( vars );
 std::vector< std::uint8_t > c( vars.size() + 1 , 0 );
 in_lace( [ & ]( void ) {
  const BDD set = mtbdd_set_from_array( v.data() , v.size() );
  mtbdd_refs_push( set );
  sylvan_sat_one( f_bdd , set , c.data() );
  mtbdd_refs_pop( 1 );
  } );
 std::vector< bool > r( vars.size() );
 for( std::size_t i = 0 ; i < r.size() ; ++i )
  r[ i ] = c[ i ] == 1;
 return( r );
 }

/*--------------------------------------------------------------------------*/
/*------------------- RELATIONAL PRODUCTS AND REDUCTIONS -------------------*/
/*--------------------------------------------------------------------------*/

Bdd PBeFS_di_unipi_it::image( const Bdd & S , const Bdd & T ,
                              const Bdd & vars )
{
 return( Bdd( sylvan_relnext( S.f_bdd , T.f_bdd , vars.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd PBeFS_di_unipi_it::preimage( const Bdd & S , const Bdd & T ,
                                 const Bdd & vars )
{
 return( Bdd( sylvan_relprev( T.f_bdd , S.f_bdd , vars.f_bdd ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PBeFS_di_unipi_it::image_into( const Bdd & S , const Bdd & T ,
                                    const Bdd & vars , bool backward ,
                                    Bdd & out )
{
 RUN( pbefs_image_into , S.f_bdd , T.f_bdd , vars.f_bdd ,
      backward ? 1 : 0 , & out.f_bdd );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PBeFS_di_unipi_it::images( const Bdd & S , const std::vector< Bdd > & T ,
                                const std::vector< Bdd > & vars ,
                                bool backward , std::vector< Bdd > & out )
{
 if( T.size() != vars.size() )
  throw( std::invalid_argument( "images: as many sets of variables as "
                                "relations" ) );
 out.assign( T.size() , Bdd() );
 if( T.empty() )
  return;
 std::vector< std::uint64_t > rel( T.size() ) , var( T.size() );
 for( std::size_t i = 0 ; i < T.size() ; ++i ) {
  rel[ i ] = T[ i ].f_bdd;
  var[ i ] = vars[ i ].f_bdd;
  }
 // the tasks write the results in registered handles
 std::vector< std::uint64_t > res( T.size() , sylvan_false );
 for( auto & r : res )
  mtbdd_protect( & r );
 const ImageJob job{ S.f_bdd , rel.data() , var.data() , res.data() ,
                     backward ? 1 : 0 };
 RUN( pbefs_images , & job , std::size_t( 0 ) , T.size() );
 for( std::size_t i = 0 ; i < T.size() ; ++i )
  out[ i ] = Bdd( res[ i ] );
 for( auto & r : res )
  mtbdd_unprotect( & r );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

Bdd PBeFS_di_unipi_it::disjoin( const std::vector< Bdd > & v )
{
 if( v.empty() )
  return( Bdd() );
 std::vector< std::uint64_t > w( v.size() );
 for( std::size_t i = 0 ; i < v.size() ; ++i ) {
  w[ i ] = v[ i ].f_bdd;
  mtbdd_protect( & w[ i ] );
  }
 RUN( pbefs_disjoin , w.data() , std::size_t( 0 ) , w.size() );
 const Bdd r( w[ 0 ] );
 for( auto & x : w )
  mtbdd_unprotect( & x );
 return( r );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::size_t PBeFS_di_unipi_it::nodes( const std::vector< Bdd > & v )
{
 std::vector< std::uint64_t > w( v.size() );
 for( std::size_t i = 0 ; i < v.size() ; ++i )
  w[ i ] = v[ i ].f_bdd;
 return( mtbdd_nodecount_more( w.data() , w.size() ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- METHODS OF BddRuntime --------------------------*/
/*--------------------------------------------------------------------------*/

BddRuntime::BddRuntime( unsigned workers , std::size_t memory )
{
 if( g_running.exchange( true ) )
  throw( std::logic_error( "BddRuntime: Sylvan is running already" ) );
 lace_start( workers , 0 );
 // the table of the nodes twice the cache, both an eighth of their largest
 // size at the start
 sylvan_set_limits( memory , 1 , 3 );
 sylvan_init_package();
 sylvan_init_mtbdd();
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

BddRuntime::~BddRuntime()
{
 sylvan_quit();
 lace_stop();
 g_running = false;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

unsigned BddRuntime::workers( void )
{
 return( lace_workers() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool BddRuntime::running( void ) noexcept
{
 return( g_running );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double BddRuntime::filled( void )
{
 std::size_t used = 0 , total = 0;
 sylvan_table_usage( & used , & total );
 const std::size_t mx = llmsset_get_max_size( ::sylvan::nodes );
 return( mx ? double( used ) / double( mx ) : 0 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void BddRuntime::collect( void )
{
 sylvan_gc();
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Bdd.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
