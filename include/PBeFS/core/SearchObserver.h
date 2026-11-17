/*--------------------------------------------------------------------------*/
/*------------------------- File SearchObserver.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for NullObserver, the observer of a search that observes
 * nothing, and for the hooks that a search calls on its observer.
 *
 * An observer is a template parameter of the searches that call it (the
 * last one, NullObserver by default), given to a search by set_observer()
 * as the executor is by set_executor() and kept in an ObserverSlot; the
 * search calls its hooks at the events below, each with the index of the
 * worker that calls it, \f$ 0 \f$ in a sequential search. The hooks of
 * NullObserver are empty and inline, and their arguments have no side
 * effects, so that they compile to nothing; the hooks inside the lambdas
 * that a search passes to the Domain are called only if observes< O >,
 * since a lambda captures what such a call names even when the call is
 * discarded. A worker calls the hooks with its own index only, and two
 * workers may call them at the same time: an observer that records is
 * thread-safe if it keeps apart what each worker records [see
 * TraceRecorder].
 *
 * The hooks of the best-first searches of a Domain, State and Cost being
 * those of the Domain, are
 *
 * - on_start( p ), before the search, with the number p of the workers
 *   (that of the indices the hooks may receive), and on_finish(), after
 *   it, both from the thread of the caller;
 * - on_expand( i , s , g , h , f ), when the worker i expands s with
 *   \f$ g \f$ and \f$ h \f$, \f$ f \f$ being the key of OPEN (\f$ g + h
 *   \f$, or \f$ g + w h \f$ with a weight), the goal test included, so
 *   that the expansions are as many as SearchStats::expanded;
 * - on_generate( i , s , parent , g , h ), when it generates s from the
 *   state *parent, as many as SearchStats::generated;
 * - on_duplicate( i , s , g ), when the state s generated with \f$ g \f$
 *   is already known with a \f$ g \f$ not larger, as many as
 *   SearchStats::duplicates;
 * - on_open_push( i , n ) and on_open_pop( i , n ), when the worker i
 *   pushes an entry in OPEN or pops one from it, n being the entries of
 *   that OPEN after the push or the pop;
 * - on_message( i , j , n ), when the worker i hands n nodes to the
 *   worker j;
 * - on_threshold( i , b ), when an iterative-deepening search starts the
 *   iteration of threshold b;
 * - on_solution( i , c ), when the worker i finds a goal of cost c that
 *   lowers the incumbent (with a weight, or with PA*SE, the goal
 *   accepted).
 *
 * Those of the Branch-and-Bound searches, Node and Value being those of
 * the problem and the values in its sense (of the minimization or of the
 * maximization), are on_start(), on_finish(), on_open_push(),
 * on_open_pop() and on_solution( i , v ) for a feasible value v that
 * improves the incumbent, and
 *
 * - on_evaluate( i , n , b , d ), when the worker i evaluates the node n
 *   of depth d, b being its bound after the evaluation;
 * - on_branch( i , n , b , d , k ), after it has branched the node n of
 *   bound b and depth d into k children that survive;
 * - on_prune( i , b , d ), when it prunes a node of bound b and depth d,
 *   extracted, evaluated or generated.
 *
 * A search calls only the hooks of its events: a search with no OPEN
 * calls no on_open_push(), a search of one worker no on_message().
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

#ifndef __PBeFS_SearchObserver
 #define __PBeFS_SearchObserver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <stdexcept>

#include <string>

#include <type_traits>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS NullObserver ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the observer of a search that observes nothing, the default one: each
/// hook [see SearchObserver.h] does nothing

struct NullObserver {

 void on_start( unsigned ) noexcept {}

 void on_finish( void ) noexcept {}

 template< class S , class C >
 void on_expand( unsigned , const S & , C , C , C ) noexcept {}

 template< class S , class C >
 void on_generate( unsigned , const S & , const S * , C , C ) noexcept {}

 template< class S , class C >
 void on_duplicate( unsigned , const S & , C ) noexcept {}

 void on_open_push( unsigned , std::size_t ) noexcept {}

 void on_open_pop( unsigned , std::size_t ) noexcept {}

 void on_message( unsigned , unsigned , std::size_t ) noexcept {}

 template< class C >
 void on_threshold( unsigned , C ) noexcept {}

 template< class C >
 void on_solution( unsigned , C ) noexcept {}

 template< class N , class V >
 void on_evaluate( unsigned , const N & , V , unsigned ) noexcept {}

 template< class N , class V >
 void on_branch( unsigned , const N & , V , unsigned , std::size_t )
  noexcept {}

 template< class V >
 void on_prune( unsigned , V , unsigned ) noexcept {}

 };  // end( struct( NullObserver ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS ObserverSlot ---------------------------*/
/*--------------------------------------------------------------------------*/
/// where a search keeps the observer of type O that set_observer() gives
/// it: a pointer, nullptr until then, so that get() throws
/// std::logic_error if the search runs before it

template< class O >
class ObserverSlot {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// sets the observer
 void set( O & o ) noexcept { f_o = &o; }

 /// the observer, or std::logic_error naming the search if there is none
 [[nodiscard]] O & get( const char * search ) const {
  if( ! f_o )
   throw( std::logic_error( std::string( search ) +
                            ": set_observer() before solve()" ) );
  return( *f_o );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 O * f_o = nullptr;                                ///< the observer

 };  // end( class( ObserverSlot ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the ObserverSlot of NullObserver, which keeps nothing, so that as a
/// member with [[no_unique_address]] it takes no room, and the search reads
/// no pointer

template<>
class ObserverSlot< NullObserver > {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// does nothing
 void set( NullObserver & ) noexcept {}

 /// a NullObserver, which has no state
 [[nodiscard]] NullObserver & get( const char * ) const noexcept {
  static NullObserver none;
  return( none );
  }

 };  // end( class( ObserverSlot< NullObserver > ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- TYPES AND CONSTANTS ---------------------------*/
/*--------------------------------------------------------------------------*/
/// true if O (a reference to it, possibly const) is not NullObserver

template< class O >
inline constexpr bool observes =
 ! std::is_same_v< std::remove_cvref_t< O > , NullObserver >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// how a part of a search (a worker, a recursion) holds the observer of
/// type O: a NullObserver, which takes no room as a member with
/// [[no_unique_address]], or a reference to any other

template< class O >
using observer_ref = std::conditional_t< observes< O > , O & ,
                                         NullObserver >;

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_SearchObserver */

/*--------------------------------------------------------------------------*/
/*----------------------- End File SearchObserver.h ------------------------*/
/*--------------------------------------------------------------------------*/
