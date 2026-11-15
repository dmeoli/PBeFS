/*--------------------------------------------------------------------------*/
/*------------------------- File SymbolicBench.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * SymbolicBench, the benchmark of the symbolic searches [see
 * SymbolicSearch.h]: it solves one instance with one search and prints a
 * line of CSV with the cost and the counters, which do not depend on the
 * workers, and, with --astar, the expansions of A* on the same instance.
 *
 *     SymbolicBench [options]
 *       --domain tiles8 | tiles15 | pancake<N>   (N from 5 to 12)
 *       --instance <k>   the instance k of Korf (tiles15), the seed of a
 *                        random stack (pancake), the seed of a walk of 4 k
 *                        moves of the blank from the goal (tiles8)
 *       --data <dir>     the directory of the data (../data)
 *       --dir f | b | bi the direction (f)
 *       --blind          Dijkstra's algorithm, h = 0
 *       --bfs            every move costing 1, without h
 *       --reach          all the states reachable from the start, by layers
 *       --images seq | lace | thread | openmp | fastflow   (seq)
 *       --threads <p>    the workers of the executor (2)
 *       --lace <w>       the workers of Lace (1)
 *       --merge <n>      the nodes of a partition, 0 one per operator
 *       --memory <MB>    the memory of Sylvan (2048)
 *       --astar          also the expansions of A*
 *       --header         the header of the CSV first
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

#include <cstdlib>

#include <fstream>

#include <iostream>

#include <sstream>

#include <stdexcept>

#include <string>

#include "PBeFS/astar/AStar.h"

#include "PBeFS/symbolic.h"

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

/*--------------------------------------------------------------------------*/
/*-------------------------------- OPTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

struct Args {
 std::string domain = "tiles8";
 int instance = 1;
 std::string data = "../data";
 std::string dir = "f";
 bool blind = false , bfs = false , reach = false , astar = false ,
      header = false;
 std::string images = "seq";
 unsigned threads = 2 , lace = 1;
 long merge = -1;
 std::size_t memory = 2048;
 };

/*--------------------------------------------------------------------------*/
// the line of CSV of a solution

template< class Sol >
void print( const Args & a , const Sol & s , std::uint64_t astar_expanded )
{
 const auto & t = s.stats;
 std::cout << a.domain << ',' << a.instance << ',' << a.dir << ','
           << ( a.bfs ? "bfs" : a.blind ? "dijkstra" : "bdda" ) << ','
           << a.images << ',' << a.threads << ',' << a.lace << ','
           << ( s.solved() ? std::to_string( s.cost ) : "-" ) << ','
           << t.steps << ',' << t.forward_steps << ',' << t.backward_steps
           << ',' << t.images << ',' << t.path_images << ',' << t.expanded
           << ',' << t.expanded_nodes << ',' << t.max_set_nodes << ','
           << t.max_image_nodes << ',' << t.partitions << ','
           << t.relation_nodes << ',' << t.heuristic_nodes << ','
           << astar_expanded << std::endl;
 }

/*--------------------------------------------------------------------------*/
// the search of the symbolic domain sd, on the executor of a.images

template< class SD , class E >
void run( const Args & a , const SD & sd , const E & e )
{
 SymbolicOptions o;
 o.direction = a.dir == "b" ? SymbolicDirection::backward :
               a.dir == "bi" ? SymbolicDirection::bidirectional :
                               SymbolicDirection::forward;
 o.heuristic = ! a.blind;
 o.unit_costs = a.bfs;
 o.threads = a.threads;
 if( a.merge >= 0 )
  o.max_relation_nodes = std::size_t( a.merge );
 o.images = a.images == "seq" ? SymbolicImages::sequential :
            a.images == "lace" ? SymbolicImages::lace :
                                 SymbolicImages::executor;
 const auto go = [ & ]< class Executor >( void ) {
  if( a.reach ) {
   const auto r = symbolic_reachable< SD , Executor >( sd , false , o );
   std::cout << a.domain << ',' << a.instance << ",reach," << r.states
             << ',' << r.layers.size() << ',' << r.images << ','
             << r.max_set_nodes << ',' << r.reached_nodes << std::endl;
   return;
   }
  const auto s = symbolic_search< SD , Executor >( sd , o );
  print( a , s , a.astar ? astar( e ).stats.expanded : 0 );
  };
 if( a.images == "openmp" ) {
  #ifdef PBEFS_HAS_OPENMP
   go.template operator()< ADS::OpenMPExecutor >();
   return;
  #else
   throw( std::invalid_argument( "no OpenMP" ) );
  #endif
  }
 if( a.images == "fastflow" ) {
  #ifdef PBEFS_HAS_FASTFLOW
   go.template operator()< ADS::FastFlowExecutor >();
   return;
  #else
   throw( std::invalid_argument( "no FastFlow" ) );
  #endif
  }
 go.template operator()< ADS::ThreadPool >();
 }

/*--------------------------------------------------------------------------*/
// the pancakes with N pancakes, or N - 1 down to 5

template< int N >
void pancakes( const Args & a , int n )
{
 if constexpr( N >= 5 ) {
  if( n != N ) {
   pancakes< N - 1 >( a , n );
   return;
   }
  const auto s = PancakePuzzle< N >::random( std::uint64_t( a.instance ) );
  run( a , SymbolicPancake< N >( s ) , PancakePuzzle< N >( s ) );
  }
 else
  throw( std::invalid_argument( "from 5 to 12 pancakes" ) );
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 Args a;
 try {
  for( int i = 1 ; i < argc ; ++i ) {
   const std::string k = argv[ i ];
   const auto next = [ & ]( void ) -> std::string {
    if( i + 1 >= argc )
     throw( std::invalid_argument( "no value of " + k ) );
    return( argv[ ++i ] );
    };
   if( k == "--domain" ) a.domain = next();
   else if( k == "--instance" ) a.instance = std::stoi( next() );
   else if( k == "--data" ) a.data = next();
   else if( k == "--dir" ) a.dir = next();
   else if( k == "--blind" ) a.blind = true;
   else if( k == "--bfs" ) a.bfs = true;
   else if( k == "--reach" ) a.reach = true;
   else if( k == "--images" ) a.images = next();
   else if( k == "--threads" ) a.threads = unsigned( std::stoul( next() ) );
   else if( k == "--lace" ) a.lace = unsigned( std::stoul( next() ) );
   else if( k == "--merge" ) a.merge = std::stol( next() );
   else if( k == "--memory" ) a.memory = std::stoul( next() );
   else if( k == "--astar" ) a.astar = true;
   else if( k == "--header" ) a.header = true;
   else
    throw( std::invalid_argument( "unknown option " + k ) );
   }
  if( a.header )
   std::cout << "domain,instance,dir,search,images,threads,lace,cost,steps,"
                "forward_steps,backward_steps,images,path_images,expanded,"
                "expanded_nodes,max_set_nodes,max_image_nodes,partitions,"
                "relation_nodes,heuristic_nodes,astar_expanded" << std::endl;
  const BddRuntime rt( a.lace , a.memory << 20 );
  if( a.domain == "tiles8" ) {
   auto b = EightPuzzle::korf_goal();
   std::uint64_t x = std::uint64_t( a.instance );
   int z = 0;
   for( int k = 0 ; k < 4 * a.instance ; ++k ) {
    x = ADS::splitmix64( x );
    int c = z;
    switch( x % 4 ) {
     case( 0 ): c = z >= 3 ? z - 3 : z; break;
     case( 1 ): c = z % 3 ? z - 1 : z; break;
     case( 2 ): c = z % 3 < 2 ? z + 1 : z; break;
     default: c = z < 6 ? z + 3 : z;
     }
    std::swap( b[ z ] , b[ c ] );
    z = c;
    }
   run( a , SymbolicSlidingTile< 3 , 3 >( b ) , EightPuzzle( b ) );
   }
  else
   if( a.domain == "tiles15" ) {
    std::ifstream in( a.data + "/15-puzzle/korf100.txt" );
    int id = 0;
    for( std::string line ; std::getline( in , line ) ; ) {
     if( line.empty() || ( line[ 0 ] == '#' ) || ( ++id != a.instance ) )
      continue;
     std::istringstream is( line );
     FifteenPuzzle::board_type b{};
     for( auto & t : b ) {
      int v;
      is >> v;
      t = std::uint8_t( v );
      }
     run( a , SymbolicSlidingTile< 4 , 4 >( b ) , FifteenPuzzle( b ) );
     return( 0 );
     }
    throw( std::invalid_argument( "no instance of Korf in " + a.data ) );
    }
   else
    if( a.domain.rfind( "pancake" , 0 ) == 0 )
     pancakes< 12 >( a , std::stoi( a.domain.substr( 7 ) ) );
    else
     throw( std::invalid_argument( "unknown domain " + a.domain ) );
  }
 catch( std::exception & e ) {
  std::cerr << "SymbolicBench: " << e.what() << std::endl;
  return( 2 );
  }
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File SymbolicBench.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
