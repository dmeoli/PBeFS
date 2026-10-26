/*--------------------------------------------------------------------------*/
/*-------------------------- File Distribution.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the work distribution functions of hash-based parallel
 * A*: each one maps a state to its owner among p workers [see section 6 of
 * Fukunaga et al., 2017, and Jinnai and Fukunaga, JAIR 2017].
 *
 * A distribution is a class with a method owner( s , p ) returning the
 * worker in \f$ \{ 0 , \ldots , p - 1 \} \f$ that owns the state
 * \f$ s \f$. In the notation of Jinnai and
 * Fukunaga, HDA*[ hash , abstraction ]:
 *
 * - ZobristDistribution is HDA*[Z] (Kishimoto, Fukunaga and Botea);
 *
 * - AbstractZobristDistribution is HDA*[Z,Afeature], i.e., AZHDA*;
 *
 * - AbstractStateDistribution is HDA*[Z,Astate], i.e., AHDA* with Zobrist
 *   hashing of the abstract state (Burns et al.);
 *
 * - PerfectHashDistribution is HDA*[P], the perfect hash modulo p;
 *
 * - StateHashDistribution uses the hash of the closed list;
 *
 * - RandomDistribution sends a state to a random worker, so that the
 *   duplicates are not detected across the workers (Kumar, Ramesh and
 *   Rao; Karp and Zhang): it is the baseline that HDA* is compared to.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_Distribution
 #define __PBeFS_Distribution
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cstddef>

#include <cstdint>

#include <random>

#include <string>

#include <vector>

#if defined( _MSC_VER )
 #include <intrin.h>
#endif

#include "PBeFS/astar/Domain.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the worker in \f$ \{ 0 , \ldots , p - 1 \} \f$ of a 64-bit hash value
/** Maps \p h to \f$ \lfloor h p / 2^{64} \rfloor \f$, the multiply-high
 * reduction of Lemire, which uses the high bits of \f$ h \f$ and is
 * unbiased for a uniform \f$ h \f$. */

[[nodiscard]] inline unsigned reduce( std::uint64_t h , unsigned p ) noexcept
{
 #if defined( __SIZEOF_INT128__ )
  __extension__ using u128 = unsigned __int128;
  return( unsigned( ( u128( h ) * p ) >> 64 ) );
 #elif defined( _MSC_VER )
  return( unsigned( __umulh( h , p ) ) );
 #else
  return( unsigned( h % p ) );
 #endif
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a table of n random 64-bit strings, reproducible from the seed

[[nodiscard]] std::vector< std::uint64_t > random_table( std::size_t n ,
                                                     std::uint64_t seed );

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ZobristDistribution ------------------------*/
/*--------------------------------------------------------------------------*/
/// HDA*[Z]: Zobrist hashing of the features of the state
/** \f$ Z( s ) = R[ x_0 ] \oplus \cdots \oplus R[ x_n ] \f$, where
 * \f$ x_0 , \ldots , x_n \f$ are the features of \f$ s \f$, \f$ R \f$ a
 * table of random 64-bit strings and \f$ \oplus \f$ the exclusive or of
 * the bits [see Zobrist, 1970,
 * and Kishimoto, Fukunaga and Botea, AIJ 2013]. Zobrist hashing is simple
 * tabulation hashing, which is 3-independent: it spreads the states
 * uniformly over the workers, but nearly every successor lands on another
 * worker than its parent. */

template< FeatureDomain D >
class ZobristDistribution {

 public:

 /// constructor, drawing the table of D from the seed

 explicit ZobristDistribution( const D & domain , std::uint64_t seed = 1 )
  : f_domain( domain ) , v_r( random_table( domain.num_features() , seed ) )
 {}

 /// the Zobrist hash of s

 [[nodiscard]] std::uint64_t hash( const state_t< D > & s ) const {
  std::uint64_t z = 0;
  f_domain.for_each_feature( s , [ & ]( std::size_t i ) { z ^= v_r[ i ]; } );
  return( z );
  }

 /// the owner of s among p workers

 [[nodiscard]] unsigned owner( const state_t< D > & s , unsigned p ) const {
  return( reduce( hash( s ) , p ) );
  }

 /// the name in the notation of Jinnai and Fukunaga

 [[nodiscard]] static std::string name( void ) { return( "Z" ); }

 private:

 const D & f_domain;                                   ///< the Domain
 std::vector< std::uint64_t > v_r;                     ///< R, by feature

 };  // end( class( ZobristDistribution ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS AbstractZobristDistribution --------------------*/
/*--------------------------------------------------------------------------*/
/// HDA*[Z,Afeature]: Zobrist hashing of the abstract features (AZHDA*)
/** \f$ AZ( s ) = R[ A( x_0 ) ] \oplus \cdots \oplus R[ A( x_n ) ] \f$,
 * where \f$ A \f$ is the feature projection of the Domain [see Jinnai and
 * Fukunaga, JAIR 2017,
 * equation (2)]. A successor that changes no abstract feature keeps the
 * owner of its parent, which cuts the communication, while the Zobrist
 * hashing of all the abstract features keeps the load balanced. The table
 * is laid out by raw feature (Algorithm 3 of the paper), so that hashing
 * costs exactly as much as with ZobristDistribution. */

template< AbstractFeatureDomain D >
class AbstractZobristDistribution {

 public:

 /// constructor, drawing the table of the abstract features from the seed

 explicit AbstractZobristDistribution( const D & domain ,
                                       std::uint64_t seed = 1 )
  : f_domain( domain ) , v_r( domain.num_features() ) {
  const auto r0 = random_table( domain.num_abstract_features() , seed );
  for( std::size_t i = 0 ; i < v_r.size() ; ++i )
   v_r[ i ] = r0[ domain.abstract_feature( i ) ];
  }

 /// the abstract Zobrist hash of s

 [[nodiscard]] std::uint64_t hash( const state_t< D > & s ) const {
  std::uint64_t z = 0;
  f_domain.for_each_feature( s , [ & ]( std::size_t i ) { z ^= v_r[ i ]; } );
  return( z );
  }

 /// the owner of s among p workers

 [[nodiscard]] unsigned owner( const state_t< D > & s , unsigned p ) const {
  return( reduce( hash( s ) , p ) );
  }

 /// the name in the notation of Jinnai and Fukunaga

 [[nodiscard]] static std::string name( void ) { return( "AZ" ); }

 private:

 const D & f_domain;                                   ///< the Domain
 std::vector< std::uint64_t > v_r;           ///< R[ A( x ) ], by feature

 };  // end( class( AbstractZobristDistribution ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS AbstractStateDistribution ---------------------*/
/*--------------------------------------------------------------------------*/
/// HDA*[Z,Astate]: the abstract state of the Domain, hashed (AHDA*)
/** All the states with the same abstract state belong to the same worker
 * [see Burns et al., JAIR 2010]: the communication is minimal, but nothing
 * balances the load, and the search overhead grows. */

template< AbstractStateDomain D >
class AbstractStateDistribution {

 public:

 explicit AbstractStateDistribution( const D & domain ,
                                     std::uint64_t seed = 1 )
  : f_domain( domain ) , f_seed( mix64( seed ) ) {}

 /// the owner of s among p workers

 [[nodiscard]] unsigned owner( const state_t< D > & s , unsigned p ) const {
  return( reduce( mix64( f_domain.abstract_state( s ) ^ f_seed ) , p ) );
  }

 /// the name in the notation of Jinnai and Fukunaga

 [[nodiscard]] static std::string name( void ) { return( "A" ); }

 private:

 const D & f_domain;                                   ///< the Domain
 std::uint64_t f_seed;                         ///< mixed into the hash

 };  // end( class( AbstractStateDistribution ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS PerfectHashDistribution ----------------------*/
/*--------------------------------------------------------------------------*/
/// HDA*[P]: the perfect hash of the Domain modulo p
/** The distribution "HDA*" of Burns et al. (JAIR 2010): with the rank of
 * the permutation of the tiles, modulo p only the last few tiles matter,
 * and Jinnai and Fukunaga (JAIR 2017) find both a high search and a high
 * communication overhead. It is kept to reproduce that comparison. */

template< PerfectHashDomain D >
class PerfectHashDistribution {

 public:

 explicit PerfectHashDistribution( const D & domain , std::uint64_t = 1 )
  : f_domain( domain ) {}

 [[nodiscard]] unsigned owner( const state_t< D > & s , unsigned p ) const {
  return( unsigned( f_domain.perfect_hash( s ) % p ) );
  }

 [[nodiscard]] static std::string name( void ) { return( "P" ); }

 private:

 const D & f_domain;                                   ///< the Domain

 };  // end( class( PerfectHashDistribution ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS StateHashDistribution -----------------------*/
/*--------------------------------------------------------------------------*/
/// the hash of the closed list, remixed with a seed
/** A uniform distribution for any Domain, with no feature to define; it
 * behaves as Zobrist hashing, but its value cannot be updated
 * incrementally. */

template< Domain D >
class StateHashDistribution {

 public:

 explicit StateHashDistribution( const D & domain , std::uint64_t seed = 1 )
  : f_domain( domain ) , f_seed( mix64( seed ) ) {}

 [[nodiscard]] unsigned owner( const state_t< D > & s , unsigned p ) const {
  return( reduce( mix64( f_domain.hash( s ) ^ f_seed ) , p ) );
  }

 [[nodiscard]] static std::string name( void ) { return( "H" ); }

 private:

 const D & f_domain;                                   ///< the Domain
 std::uint64_t f_seed;                         ///< mixed into the hash

 };  // end( class( StateHashDistribution ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS RandomDistribution ------------------------*/
/*--------------------------------------------------------------------------*/
/// a random worker for each generated state
/** The randomized strategy of Kumar, Ramesh and Rao (1988) and of Karp and
 * Zhang (1993): a state is detected as a duplicate only if its copies
 * happen to reach the same worker. owner() is const but changes the state
 * of the generator: each worker must have its own copy. */

template< Domain D >
class RandomDistribution {

 public:

 explicit RandomDistribution( const D & , std::uint64_t seed = 1 )
  : f_gen( seed ) {}

 [[nodiscard]] unsigned owner( const state_t< D > & , unsigned p ) const {
  return( reduce( f_gen() , p ) );
  }

 [[nodiscard]] static std::string name( void ) { return( "R" ); }

 /// true: each worker needs its own copy
 static constexpr bool per_worker = true;

 private:

 mutable std::mt19937_64 f_gen;                        ///< the generator

 };  // end( class( RandomDistribution ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_Distribution */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Distribution.h -------------------------*/
/*--------------------------------------------------------------------------*/
