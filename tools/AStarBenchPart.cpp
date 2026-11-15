/*--------------------------------------------------------------------------*/
/*------------------------- File AStarBenchPart.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The searches of AStarBench on one Domain of AStarBenchDomains.h, the one
 * numbered PBEFS_BENCH_DOMAIN_INDEX (from 1), in one of its variants,
 * chosen by PBEFS_BENCH_PART: 0 instantiates solve_maybe_slow() and the
 * searches on the Domain itself, 1 the searches on its SlowDomain and 2
 * those on its SlowEdgeDomain [see AStarBenchSolve.h]. The build compiles
 * this file once for each Domain and part, so that each translation unit
 * holds the searches of a single Domain type, which g++ compiles in at
 * most about 3 GB, while a unit with all the variants of a few Domains
 * takes more than 12 GB.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AStarBenchSolve.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ THE DOMAIN --------------------------------*/
/*--------------------------------------------------------------------------*/

#if ! defined( PBEFS_BENCH_DOMAIN_INDEX ) || ! defined( PBEFS_BENCH_PART )
 #error "PBEFS_BENCH_DOMAIN_INDEX and PBEFS_BENCH_PART are both needed"
#endif

namespace {

/// a permutation puzzle with its pattern databases

template< class P >
using PermutationPDB = PDBDomain< P , PermutationHeuristic< P > >;

/// the Domains of AStarBenchDomains.h, after a void that numbers them from 1

template< class ... T > struct BenchDomains {};

#define PBEFS_BENCH_DOMAIN( ... ) , __VA_ARGS__

using AllBenchDomains = BenchDomains< void
#include "AStarBenchDomains.h"
 >;

#undef PBEFS_BENCH_DOMAIN

/// the Domain number i of L

template< std::size_t i , class L > struct BenchDomain;

template< std::size_t i , class H , class ... T >
struct BenchDomain< i , BenchDomains< H , T ... > > :
 BenchDomain< i - 1 , BenchDomains< T ... > > {};

template< class H , class ... T >
struct BenchDomain< 0 , BenchDomains< H , T ... > > { using type = H; };

}  // end( unnamed namespace )

using BenchD = BenchDomain< PBEFS_BENCH_DOMAIN_INDEX ,
                            AllBenchDomains >::type;

static_assert( ! std::is_void_v< BenchD > );

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

#if PBEFS_BENCH_PART == 0

// the variants are left to the units of the other parts

extern template Solution< state_t< SlowDomain< BenchD > > ,
                          cost_t< SlowDomain< BenchD > > >
 solve< SlowDomain< BenchD > >( const SlowDomain< BenchD > & );

extern template Solution< state_t< SlowEdgeDomain< BenchD > > ,
                          cost_t< SlowEdgeDomain< BenchD > > >
 solve< SlowEdgeDomain< BenchD > >( const SlowEdgeDomain< BenchD > & );

template Solution< state_t< BenchD > , cost_t< BenchD > >
 solve_maybe_slow< BenchD >( const BenchD & );

#elif PBEFS_BENCH_PART == 1

template Solution< state_t< SlowDomain< BenchD > > ,
                   cost_t< SlowDomain< BenchD > > >
 solve< SlowDomain< BenchD > >( const SlowDomain< BenchD > & );

#elif PBEFS_BENCH_PART == 2

template Solution< state_t< SlowEdgeDomain< BenchD > > ,
                   cost_t< SlowEdgeDomain< BenchD > > >
 solve< SlowEdgeDomain< BenchD > >( const SlowEdgeDomain< BenchD > & );

#else
 #error "PBEFS_BENCH_PART is 0, 1 or 2"
#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File AStarBenchPart.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
