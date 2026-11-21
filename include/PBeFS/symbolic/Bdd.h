/*--------------------------------------------------------------------------*/
/*------------------------------- File Bdd.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class Bdd, a reduced ordered binary decision diagram
 * of the library Sylvan [see van Dijk and van de Pol, STTT 2017], for the
 * class BddRuntime, which starts and stops Sylvan and its work-stealing
 * scheduler Lace, and for the few operations of the symbolic searches on
 * sets of states and on transition relations: the image, the preimage, the
 * images under several relations at once and the disjunction of many
 * diagrams.
 *
 * The header does not include those of Sylvan: a Bdd holds the 64-bit
 * handle of a node of Sylvan, which the implementation in src/Bdd.cpp
 * registers with Sylvan for as long as the Bdd lives, so that the garbage
 * collector of Sylvan keeps it [see BddRuntime]. The library is compiled
 * only if Sylvan is found (PBEFS_HAS_SYLVAN), and linked as the target
 * PBeFS::Symbolic.
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

#ifndef __PBeFS_Bdd
 #define __PBeFS_Bdd
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASS Bdd --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a binary decision diagram of Sylvan
/** A Bdd is a Boolean function of the variables \f$ 0 , 1 , \ldots \f$ of
 * Sylvan, which are ordered by their index; equal functions have the same
 * diagram, so that == compares the functions in constant time. Every
 * operation but image_into() must be called by one thread at a time, the
 * one that drives the search, while a BddRuntime is alive; Sylvan runs
 * each of them on the workers of Lace. */

class Bdd {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTORS ------------------------------*/
 /// the constant false

 Bdd( void ) noexcept;

 /// a copy, registered with the collector

 Bdd( const Bdd & b ) noexcept;

 /// a copy

 Bdd & operator=( const Bdd & b ) noexcept;

 /// the diagram released to the collector

 ~Bdd();

/*------------------------------- CONSTANTS --------------------------------*/
 /// the constant true

 [[nodiscard]] static Bdd one( void );

 /// the constant false

 [[nodiscard]] static Bdd zero( void ) { return( Bdd() ); }

 /// the variable v, true where it is true

 [[nodiscard]] static Bdd variable( std::uint32_t v );

 /// the conjunction of the literals of the variables vars (sorted by
 /// index), the variable vars[ i ] positive if values[ i ] is true

 [[nodiscard]] static Bdd cube( const std::vector< std::uint32_t > & vars ,
                                const std::vector< bool > & values );

 /// the set of the variables vars (in any order), as the quantifications
 /// and the counts take it

 [[nodiscard]] static Bdd variables( std::vector< std::uint32_t > vars );

/*----------------------------- THE OPERATIONS -----------------------------*/
 /// the conjunction

 [[nodiscard]] Bdd operator&( const Bdd & b ) const;

 /// the disjunction

 [[nodiscard]] Bdd operator|( const Bdd & b ) const;

 /// the difference, *this and not b

 [[nodiscard]] Bdd operator-( const Bdd & b ) const;

 /// the negation

 [[nodiscard]] Bdd operator!( void ) const;

 /// the equivalence, true where both are true or both are false

 [[nodiscard]] Bdd equiv( const Bdd & b ) const;

 Bdd & operator&=( const Bdd & b ) { return( *this = *this & b ); }

 Bdd & operator|=( const Bdd & b ) { return( *this = *this | b ); }

 Bdd & operator-=( const Bdd & b ) { return( *this = *this - b ); }

 /// \f$ \exists V \f$ of the function, V given by variables()

 [[nodiscard]] Bdd exists( const Bdd & vars ) const;

/*------------------------------ THE QUERIES -------------------------------*/
 /// true if the two are the same function

 [[nodiscard]] bool operator==( const Bdd & b ) const noexcept {
  return( f_bdd == b.f_bdd );
  }

 /// true if the function is false

 [[nodiscard]] bool empty( void ) const noexcept;

 /// true if the function is true

 [[nodiscard]] bool full( void ) const noexcept;

 /// true if *this and b have no common point

 [[nodiscard]] bool disjoint( const Bdd & b ) const;

 /// the internal nodes and the terminal of the diagram

 [[nodiscard]] std::size_t nodes( void ) const;

 /// the assignments of the variables vars (given by variables()) that
 /// satisfy the function, whose support must be in vars

 [[nodiscard]] double count( const Bdd & vars ) const;

 /// the values of the variables vars (sorted by index) in one assignment
 /// that satisfies the function, false for those it does not depend on;
 /// empty if the function is false

 [[nodiscard]] std::vector< bool > pick( const std::vector< std::uint32_t > &
                                         vars ) const;

 /// the handle of the node of Sylvan

 [[nodiscard]] std::uint64_t handle( void ) const noexcept {
  return( f_bdd );
  }

/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/

 private:

 explicit Bdd( std::uint64_t h ) noexcept;

 friend Bdd image( const Bdd & , const Bdd & , const Bdd & );
 friend Bdd preimage( const Bdd & , const Bdd & , const Bdd & );
 friend void image_into( const Bdd & , const Bdd & , const Bdd & , bool ,
                         Bdd & );
 friend void images( const Bdd & , const std::vector< Bdd > & ,
                     const std::vector< Bdd > & , bool ,
                     std::vector< Bdd > & );
 friend Bdd disjoin( const std::vector< Bdd > & );
 friend std::size_t nodes( const std::vector< Bdd > & );

 std::uint64_t f_bdd;                    ///< the handle, registered with
                                         ///< the collector at its address

 };  // end( class( Bdd ) )

/*--------------------------------------------------------------------------*/
/*------------------- RELATIONAL PRODUCTS AND REDUCTIONS -------------------*/
/*--------------------------------------------------------------------------*/
/** The variables of a state are the even ones of Sylvan, \f$ x_k = 2 k
 * \f$, and those of a successor the odd ones, \f$ x'_k = 2 k + 1 \f$: a
 * relation \f$ T( x , x' ) \f$ is a Bdd on both, a set of states \f$ S( x )
 * \f$ one on the even variables, and vars, given by Bdd::variables(), the
 * set of the variables \f$ x_k \f$ and \f$ x'_k \f$ that T reads or writes;
 * the variables of a state not in vars are left unchanged by T. */

/// the image \f$ \exists x : S( x ) \wedge T( x , x' ) \f$, renamed to x

[[nodiscard]] Bdd image( const Bdd & S , const Bdd & T , const Bdd & vars );

/// the preimage \f$ \exists x' : T( x , x' ) \wedge S( x' ) \f$

[[nodiscard]] Bdd preimage( const Bdd & S , const Bdd & T ,
                            const Bdd & vars );

/// the image of S (the preimage if backward) written in out, which may be
/// called by any number of threads at once [see the implementation]

void image_into( const Bdd & S , const Bdd & T , const Bdd & vars ,
                 bool backward , Bdd & out );

/// the images of S under every T[ i ] (the preimages if backward) in out,
/// computed at once as tasks of Lace, which balances them over its workers

void images( const Bdd & S , const std::vector< Bdd > & T ,
             const std::vector< Bdd > & vars , bool backward ,
             std::vector< Bdd > & out );

/// the disjunction of all of v, as a balanced tree of tasks of Lace

[[nodiscard]] Bdd disjoin( const std::vector< Bdd > & v );

/// the nodes of the diagrams of v, those shared counted once

[[nodiscard]] std::size_t nodes( const std::vector< Bdd > & v );

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BddRuntime ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the scope in which Sylvan and Lace run
/** Sylvan and Lace are one per process: the BddRuntime starts Lace with the
 * given number of workers (0 for one per core) and Sylvan with its table
 * of nodes and its cache within the given bytes, and stops them when it is
 * destroyed; a second BddRuntime alive at the same time throws
 * std::logic_error, and every Bdd must be destroyed before the BddRuntime.
 * The table grows up to its maximum, and a search that would fill it
 * stops before [see filled()]. */

class BddRuntime {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /// Lace with workers workers, Sylvan within memory bytes

 explicit BddRuntime( unsigned workers = 0 ,
                      std::size_t memory = std::size_t( 1 ) << 30 );

 BddRuntime( const BddRuntime & ) = delete;

 BddRuntime & operator=( const BddRuntime & ) = delete;

 /// Sylvan and Lace stopped

 ~BddRuntime();

 /// the workers of Lace

 [[nodiscard]] static unsigned workers( void );

 /// true if a BddRuntime is alive

 [[nodiscard]] static bool running( void ) noexcept;

 /// the fraction of the largest table of the nodes that is in use

 [[nodiscard]] static double filled( void );

 /// a garbage collection, which drops the nodes of no Bdd alive

 static void collect( void );

 };  // end( class( BddRuntime ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Bdd */

/*--------------------------------------------------------------------------*/
/*----------------------------- End File Bdd.h -----------------------------*/
/*--------------------------------------------------------------------------*/
