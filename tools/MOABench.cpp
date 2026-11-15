/*--------------------------------------------------------------------------*/
/*--------------------------- File MOABench.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the multi-objective searches of PBeFS: solves pairs of
 * nodes of a road graph of DIMACS with several costs, or a grid whose moves
 * cost several random integers, with a chosen search, runtime and number of
 * workers, and prints one line of comma-separated values per run.
 *
 *     MOABench [options]
 *
 *     --graph   <file.gr>,<file.gr>,...  the files of the first objectives,
 *                                        with the same arcs in the same
 *                                        order (USA-road-d, USA-road-t)
 *     --extra   arcs | degree | avg-degree | random , ...
 *                                        the objectives computed from the
 *                                        graph after those of the files
 *                                        [see MultiObjectiveGraph.h] (none)
 *     --queries <file>                   the pairs "source target" of the
 *                                        runs, numbered from 1, one per line
 *     --pair    <source>,<target>        a single pair instead
 *     --grid    <W>x<H>:<moves>:<lo>:<hi>:<seed>:<K>
 *                                        a MultiObjectiveGrid of K
 *                                        objectives instead of a graph,
 *                                        free unless --map is given
 *     --map     <Moving AI map>          the blocked cells of the grid
 *     --points  <sx>,<sy>,<gx>,<gy>      its start and goal (the corners)
 *     --ideal   0 | 1   the ideal point as the heuristic of the grid   (1)
 *     --algo    boa | emoa | ltmoa | namoa-dr | par | hd           (ltmoa)
 *     --frontier scalar | array | tree   the frontier of par       (array)
 *     --lazy    0 | 1   no dominance check at the generation           (0)
 *     --order   <i>,<j>,...  the order of the objectives (0, 1, ...)
 *     --runtime thread | openmp | fastflow                        (thread)
 *     --threads <workers>                                              (1)
 *     --batch   <nodes of a message of hd>                            (64)
 *     --flush   <extractions between flushes or updates>              (64)
 *     --share   0 | 1   the bounds shared by the workers of par        (1)
 *     --seed    <seed of the random objective and of the owners>       (1)
 *     --limit   <labels stored at most, per worker>          (0: none)
 *     --header  print the header of the columns
 *
 * The columns are those of --header: the instance (source and target, or
 * the grid), the search, its frontier, lazy, the runtime, the workers, the
 * objectives, the status, the size of the front and the counters of
 * MOSearchStats [see MultiObjective.h], then the seconds of the search,
 * which do not include the reading and the heuristic.
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

#include <fstream>

#include <iostream>

#include <map>

#include <memory>

#include <sstream>

#include <stdexcept>

#include <string>

#include <utility>

#include <vector>

#include "PBeFS/astar/MOAStar.h"

#include "PBeFS/astar/MultiObjectiveGraph.h"

#include "PBeFS/astar/MultiObjectiveGrid.h"

#include "PBeFS/astar/ParallelMOAStar.h"

#ifdef PBEFS_HAS_OPENMP
 #include "ADS/parallel/OpenMPExecutor.h"
#endif

#ifdef PBEFS_HAS_FASTFLOW
 #include "ADS/parallel/FastFlowExecutor.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

namespace {

/*--------------------------------------------------------------------------*/
/*----------------------- TYPES AND GLOBAL VARIABLES -----------------------*/
/*--------------------------------------------------------------------------*/

// the options of the command line, by name
std::map< std::string , std::string > opt = {
 { "graph" , "" } , { "extra" , "" } , { "queries" , "" } , { "pair" , "" } ,
 { "grid" , "" } , { "map" , "" } , { "points" , "" } , { "ideal" , "1" } ,
 { "algo" , "ltmoa" } , { "frontier" , "array" } , { "lazy" , "0" } ,
 { "order" , "" } , { "runtime" , "thread" } , { "threads" , "1" } ,
 { "batch" , "64" } , { "flush" , "64" } , { "share" , "1" } ,
 { "seed" , "1" } , { "limit" , "0" } };

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the fields of s separated by c

std::vector< std::string > split( const std::string & s , char c )
{
 std::vector< std::string > v;
 std::istringstream in( s );
 for( std::string x ; std::getline( in , x , c ) ; )
  if( ! x.empty() )
   v.push_back( x );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of the parallel searches

ParallelMOOptions options( void )
{
 ParallelMOOptions o;
 o.lazy = opt[ "lazy" ] == "1";
 for( const auto & x : split( opt[ "order" ] , ',' ) )
  o.order.push_back( unsigned( std::stoul( x ) ) );
 o.paths = false;
 o.node_limit = std::stoull( opt[ "limit" ] );
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoull( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.seed = std::stoull( opt[ "seed" ] );
 o.share_bounds = opt[ "share" ] == "1";
 return( o );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel search of --algo on d, on the executor Executor

template< class Executor , class D >
pareto_solution_t< D > solve_parallel( const D & d )
{
 constexpr std::size_t K = objectives_of< D >;
 const auto o = options();
 Executor ex( o.threads );
 const std::string & a = opt[ "algo" ] , & f = opt[ "frontier" ];
 if( a == "hd" ) {
  HDMOAStar< D , Executor > s( d , o );
  s.set_executor( ex );
  return( s.solve() );
  }
 if( a != "par" )
  throw( std::invalid_argument( "unknown --algo " + a ) );
 if( f == "tree" ) {
  MOAStarPar< D , EMOAFrontier , Executor > s( d , o );
  s.set_executor( ex );
  return( s.solve() );
  }
 if( f == "scalar" ) {
  if constexpr( K == 2 ) {
   MOAStarPar< D , ScalarFrontier , Executor > s( d , o );
   s.set_executor( ex );
   return( s.solve() );
   }
  else
   throw( std::invalid_argument( "--frontier scalar is for two "
                                 "objectives" ) );
  }
 MOAStarPar< D , ArrayFrontier , Executor > s( d , o );
 s.set_executor( ex );
 return( s.solve() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the search of --algo on d

template< class D >
pareto_solution_t< D > solve( const D & d )
{
 const std::string & a = opt[ "algo" ];
 if( ( a == "par" ) || ( a == "hd" ) ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( solve_parallel< ADS::ThreadPool >( d ) );
  #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( solve_parallel< ADS::OpenMPExecutor >( d ) );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
   if( r == "fastflow" )
    return( solve_parallel< ADS::FastFlowExecutor >( d ) );
  #endif
  throw( std::invalid_argument( "unknown --runtime " + r ) );
  }
 const MOSearchOptions o = options();
 return( mo_search( d , mo_algorithm( a ) , o ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// one line of the results of the run of d labelled by instance

template< class D >
void run( const D & d , const std::string & instance )
{
 const auto s = solve( d );
 const MOSearchStats & st = s.stats;
 const std::string & a = opt[ "algo" ];
 const bool par = ( a == "par" ) || ( a == "hd" );
 std::cout << instance << "," << a << ","
           << ( a == "par" ? opt[ "frontier" ] : "-" ) << ","
           << opt[ "lazy" ] << "," << ( par ? opt[ "runtime" ] : "-" ) << ","
           << ( par ? opt[ "threads" ] : "1" ) << ","
           << objectives_of< D > << "," << int( s.status ) << ","
           << s.front.size() << "," << st.extracted << "," << st.expanded
           << "," << st.generated << "," << st.inserted << "," << st.checks
           << "," << st.solutions << "," << st.stored << "," << st.sent
           << "," << st.messages << "," << st.seconds << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the runs on the road graph of --graph and --extra with K objectives

template< std::size_t K >
void run_graph( void )
{
 std::vector< std::unique_ptr< std::ifstream > > files;
 std::vector< std::istream * > in;
 for( const auto & f : split( opt[ "graph" ] , ',' ) ) {
  files.push_back( std::make_unique< std::ifstream >( f ) );
  if( ! *files.back() )
   throw( std::invalid_argument( "cannot read " + f ) );
  in.push_back( files.back().get() );
  }
 std::vector< ExtraCost > extra;
 for( const auto & e : split( opt[ "extra" ] , ',' ) )
  extra.push_back( extra_cost( e ) );
 auto g = read_dimacs_objectives< long , K >( in , extra , 1 , 1 ,
                                              std::stoull( opt[ "seed" ] ) );

 std::vector< std::pair< std::uint32_t , std::uint32_t > > pairs;
 if( ! opt[ "pair" ].empty() ) {
  const auto p = split( opt[ "pair" ] , ',' );
  pairs.emplace_back( std::stoul( p.at( 0 ) ) , std::stoul( p.at( 1 ) ) );
  }
 else {
  std::ifstream q( opt[ "queries" ] );
  if( ! q )
   throw( std::invalid_argument( "cannot read " + opt[ "queries" ] ) );
  for( std::uint32_t s , t ; q >> s >> t ; )
   pairs.emplace_back( s , t );
  }
 for( const auto & [ s , t ] : pairs ) {
  if( ( s == 0 ) || ( t == 0 ) || ( s > g.size() ) || ( t > g.size() ) )
   throw( std::invalid_argument( "a pair out of the graph" ) );
  g.set_source( s - 1 );
  g.set_target( t - 1 );
  g.set_ideal_heuristic();
  run( g , std::to_string( s ) + " " + std::to_string( t ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the run on the grid of --grid with K objectives

template< std::size_t K >
void run_grid( const std::vector< std::string > & f )
{
 const auto wh = split( f.at( 0 ) , 'x' );
 const int w = std::stoi( wh.at( 0 ) ) , h = std::stoi( wh.at( 1 ) );
 GridMap map( w , h );
 if( ! opt[ "map" ].empty() ) {
  std::ifstream in( opt[ "map" ] );
  if( ! in )
   throw( std::invalid_argument( "cannot read " + opt[ "map" ] ) );
  map = GridMap::read_movingai( in );
  }
 std::vector< int > p{ 0 , 0 , map.width() - 1 , map.height() - 1 };
 if( ! opt[ "points" ].empty() ) {
  const auto x = split( opt[ "points" ] , ',' );
  for( std::size_t i = 0 ; i < 4 ; ++i )
   p[ i ] = std::stoi( x.at( i ) );
  }
 MultiObjectiveGrid< long , K > g( map , std::stoi( f.at( 1 ) ) , p[ 0 ] ,
                                   p[ 1 ] , p[ 2 ] , p[ 3 ] ,
                                   std::stol( f.at( 2 ) ) ,
                                   std::stol( f.at( 3 ) ) ,
                                   std::stoull( f.at( 4 ) ) );
 if( opt[ "ideal" ] == "1" )
  g.set_ideal_heuristic();
 run( g , opt[ "grid" ] );
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 try {
  for( int i = 1 ; i < argc ; ++i ) {
   const std::string a = argv[ i ];
   if( a == "--header" ) {
    std::cout << "instance,algo,frontier,lazy,runtime,threads,objectives,"
                 "status,front,extracted,expanded,generated,inserted,"
                 "checks,solutions,stored,sent,messages,seconds"
              << std::endl;
    continue;
    }
   if( ( a.rfind( "--" , 0 ) != 0 ) || ! opt.count( a.substr( 2 ) ) ||
       ( i + 1 >= argc ) )
    throw( std::invalid_argument( "bad option " + a ) );
   opt[ a.substr( 2 ) ] = argv[ ++i ];
   }
  std::size_t k;
  std::vector< std::string > grid;
  if( ! opt[ "grid" ].empty() ) {
   grid = split( opt[ "grid" ] , ':' );
   k = std::stoul( grid.at( 5 ) );
   }
  else
   k = split( opt[ "graph" ] , ',' ).size() +
       split( opt[ "extra" ] , ',' ).size();
  switch( k ) {
   case( 2 ): grid.empty() ? run_graph< 2 >() : run_grid< 2 >( grid ); break;
   case( 3 ): grid.empty() ? run_graph< 3 >() : run_grid< 3 >( grid ); break;
   case( 4 ): grid.empty() ? run_graph< 4 >() : run_grid< 4 >( grid ); break;
   case( 5 ): grid.empty() ? run_graph< 5 >() : run_grid< 5 >( grid ); break;
   default:
    throw( std::invalid_argument( "the objectives must be 2 to 5" ) );
   }
  }
 catch( std::exception & e ) {
  std::cerr << "MOABench: " << e.what() << std::endl;
  return( 1 );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File MOABench.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
