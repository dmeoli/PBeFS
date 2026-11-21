/*--------------------------------------------------------------------------*/
/*----------------------------- File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tester of PBeFS: each case, selected by the first argument, solves
 * instances of known optimal cost with one family of searches and checks
 * both the cost and the path, which must start at the initial state, end at
 * a goal, be made of edges of the Domain and cost what the Solution says.
 *
 * The references do not come from the library: a breadth-first search of
 * the whole state space of the 8-puzzle, the optimal costs of the 100
 * instances of Korf (1985), Dijkstra's algorithm on the grids and on the
 * boards, the dynamic program on the lattice of the alignments, the one of
 * Held and Karp on the cities, and the distances on the map of Romania of
 * Russell and Norvig; the symbolic searches are checked against the
 * explicit ones [see symbolic.cpp].
 *
 *     test astar | grid | graph | limits | reopen | ddd | idd
 *     test msa | board | files
 *     test bidir | bidir-thread | bidir-openmp | bidir-fastflow
 *     test pem-thread | pem-openmp | pem-fastflow
 *     test pdb | pdb-search
 *     test hda-thread | hda-openmp | hda-fastflow
 *     test farm-thread | farm-openmp | farm-fastflow
 *     test spa-thread | spa-openmp
 *     test mqa-thread | mqa-openmp
 *     test pbnf-thread | pbnf-openmp | pbnf-fastflow | pbnf-options
 *     test pase-thread
 *     test pida-thread | pida-openmp | pida-fastflow
 *     test bfs | batch | psdd-thread | psdd-openmp | psdd-fastflow
 *     test tds-thread | tds-openmp | tds-fastflow
 *     test multi-thread | multi-openmp | multi-fastflow
 *     test planning | planning-thread | planning-openmp | planning-fastflow
 *     test distributions | hda-owners
 *     test lpa | hdlpa-thread | hdlpa-openmp | hdlpa-fastflow
 *     test shift | weighted-wrapper | domains | manual
 *     test permutation | abstraction | pancake | topspin | hanoi | dimacs
 *     test rubik | rubik-pdb
 *     test gbfs | gbfs-queues | ehc | iterated | width | classic
 *     test realtime | hpa | pgbfs-rounds
 *     test pgbfs-thread | pgbfs-openmp | pgbfs-fastflow
 *     test bb | bb-shared-thread | bb-dist-thread | bb-steal-thread
 *     test bb-mq-thread | bb-mq-openmp
 *     test bb-shared-openmp | bb-dist-openmp | bb-steal-openmp
 *     test bb-shared-fastflow | bb-dist-fastflow | bb-steal-fastflow
 *     test bb-farm-thread | bb-farm-openmp | bb-farm-fastflow
 *     test bb-det-thread | bb-det-openmp | bb-det-fastflow
 *     test bb-split-thread | bb-split-openmp | bb-split-fastflow
 *     test bb-scheme | bb-cutoff | bb-values
 *     test bb-skel-thread | bb-skel-openmp | bb-skel-fastflow
 *     test bb-ivm-thread | bb-ivm-openmp | bb-ivm-fastflow
 *     test bb-modes | bb-race | bb-checkpoint | bb-domains
 *     test moa | moa-thread | moa-openmp | moa-fastflow
 *     test sym-bfs | sym-dijkstra | sym-astar | sym-sas | sym-pdb
 *     test sym-variants | sym-topk
 *     test sym-par-thread | sym-par-openmp | sym-par-fastflow (Sylvan)
 *     test mpi | bb-mpi | pida-mpi | bb-hub-mpi   (under mpiexec)
 *     test registry <tests.txt>
 *
 * The cases are listed in test/tests.txt, which CMake and the makefile
 * read; the case registry checks that the list and the cases of the tester
 * are the same.
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

#include "test.h"

/*--------------------------------------------------------------------------*/
/*---------------------------- GLOBAL VARIABLES ----------------------------*/
/*--------------------------------------------------------------------------*/

int failures = 0;                     // the checks failed so far
bool quiet = false;                   // print only the failures

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// records the outcome of a check

void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 else
  if( ! quiet )
   std::cout << "ok: " << what << std::endl;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the skeletons, the modes, the ramp-up, the checkpoints, the
// intervals, the hubs and the domains of the B&B [see bnb.cpp]

void test_bb_skeleton_thread( const std::string & );
void test_bb_interval_thread( const std::string & );
void test_bb_modes( const std::string & );
void test_bb_race( const std::string & );
void test_bb_checkpoint( const std::string & );
void test_bb_domains( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 void test_bb_skeleton_openmp( const std::string & );
 void test_bb_interval_openmp( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 void test_bb_skeleton_fastflow( const std::string & );
 void test_bb_interval_fastflow( const std::string & );
#endif
#ifdef PBEFS_HAS_MPI
 void test_bb_hub_mpi( const std::string & );
#endif

// the multi-objective searches [see moa.cpp]
void test_moa( const std::string & );
void test_moa_thread( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 void test_moa_openmp( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 void test_moa_fastflow( const std::string & );
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the symbolic searches, with Sylvan [see symbolic.cpp]

#ifdef PBEFS_HAS_SYLVAN
 void test_sym_bfs( const std::string & );
 void test_sym_dijkstra( const std::string & );
 void test_sym_astar( const std::string & );
 void test_sym_sas( const std::string & );
 void test_sym_pdb( const std::string & );
 void test_sym_variants( const std::string & );
 void test_sym_topk( const std::string & );
 void test_sym_par_thread( const std::string & );
 #ifdef PBEFS_HAS_OPENMP
  void test_sym_par_openmp( const std::string & );
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  void test_sym_par_fastflow( const std::string & );
 #endif
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the tester, each with the runtime it needs ("-" if none):
// test/tests.txt lists the same ones [see test_registry()]

struct Case {
 const char * name;                          // the name on the command line
 const char * needs;                         // the runtime it needs
 void ( * run )( const std::string & );      // the test, given the data
 };

static const std::vector< Case > cases = {
 { "astar" , "-" , test_astar } ,
 { "open" , "-" , []( const std::string & ) { test_open(); } } ,
 { "grid" , "-" , []( const std::string & ) { test_grid(); } } ,
 { "idastar" , "-" , test_idastar } ,
 { "bfs" , "-" , test_bfs } ,
 { "batch" , "-" , test_batch } ,
 { "ddd" , "-" , test_ddd } ,
 { "idd" , "-" , test_idd } ,
 { "bidir" , "-" , test_bidirectional } ,
 { "bidir-thread" , "-" , test_bidirectional_parallel< ADS::ThreadPool > } ,
 { "pem-thread" , "-" , test_pem< ADS::ThreadPool > } ,
 { "tiles24" , "-" , test_tiles24 } ,
 { "graph" , "-" , []( const std::string & ) { test_graph(); } } ,
 { "limits" , "-" , test_limits } ,
 { "reopen" , "-" , []( const std::string & ) { test_reopen(); } } ,
 { "msa" , "-" , test_msa } ,
 { "board" , "-" , test_board } ,
 { "files" , "-" , test_files } ,
 { "shift" , "-" , test_shift } ,
 { "weighted-wrapper" , "-" , test_weighted_wrapper } ,
 { "domains" , "-" , []( const std::string & ) { test_domains(); } } ,
 { "manual" , "-" , test_manual } ,
 { "pdb" , "-" , test_pdb } ,
 { "pdb-search" , "-" , test_pdb_search } ,
 { "weighted" , "-" , test_weighted } ,
 { "anytime" , "-" , test_anytime } ,
 { "focal" , "-" , test_focal } ,
 { "lpa" , "-" , test_lpa } ,
 { "hdlpa-thread" , "-" , test_hdlpa< ADS::ThreadPool > } ,
 { "bounded-thread" , "-" , test_bounded< ADS::ThreadPool > } ,
 { "distributions" , "-" , test_distributions } ,
 { "hda-thread" , "-" , test_hda< ADS::ThreadPool > } ,
 { "hda-owners" , "-" , test_hda_owners } ,
 { "spa-thread" , "-" , test_spa< ADS::ThreadPool > } ,
 { "mqa-thread" , "-" , test_mqa< ADS::ThreadPool > } ,
 { "pbnf-thread" , "-" , test_pbnf< ADS::ThreadPool > } ,
 { "psdd-thread" , "-" , test_psdd< ADS::ThreadPool > } ,
 { "multi-thread" , "-" , test_multi< ADS::ThreadPool > } ,
 { "tds-thread" , "-" , test_tds< ADS::ThreadPool > } ,
 { "pbnf-options" , "-" , test_pbnf_options< ADS::ThreadPool > } ,
 { "pase-thread" , "-" , test_pase< ADS::ThreadPool > } ,
 { "pida-thread" , "-" , test_pida< ADS::ThreadPool > } ,
 { "permutation" , "-" , test_permutation } ,
 { "abstraction" , "-" , test_abstraction } ,
 { "pancake" , "-" , test_pancake } ,
 { "topspin" , "-" , test_topspin } ,
 { "rubik" , "-" , test_rubik } ,
 { "rubik-pdb" , "-" , test_rubik_pdb } ,
 { "hanoi" , "-" , test_hanoi } ,
 { "dimacs" , "-" , test_dimacs } ,
 { "edge" , "-" , test_edge } ,
 { "edge-thread" , "-" , test_edge_parallel< ADS::ThreadPool > } ,
 { "planning" , "-" , test_planning } ,
 { "planning-thread" , "-" , test_planning_parallel< ADS::ThreadPool > } ,
 { "bb" , "-" , test_bb } ,
 { "bb-shared-thread" , "-" , test_bb_parallel< 0 , ADS::ThreadPool > } ,
 { "bb-dist-thread" , "-" , test_bb_parallel< 1 , ADS::ThreadPool > } ,
 { "bb-steal-thread" , "-" , test_bb_parallel< 3 , ADS::ThreadPool > } ,
 { "bb-mq-thread" , "-" , test_bb_parallel< 4 , ADS::ThreadPool > } ,
 { "bb-det-thread" , "-" , test_bb_det< ADS::ThreadPool > } ,
 { "bb-split-thread" , "-" , test_bb_split< ADS::ThreadPool > } ,
 { "bb-farm-thread" , "-" , test_bb_parallel< 2 , ADS::ThreadPool > } ,
 { "farm-thread" , "-" , test_farm< ADS::ThreadPool > } ,
 { "gbfs" , "-" , test_gbfs } ,
 { "gbfs-queues" , "-" , test_gbfs_queues } ,
 { "ehc" , "-" , test_ehc } ,
 { "iterated" , "-" , test_iterated } ,
 { "width" , "-" , test_width } ,
 { "classic" , "-" , test_classic } ,
 { "realtime" , "-" , test_realtime } ,
 { "hpa" , "-" , test_hpa } ,
 { "pgbfs-rounds" , "-" , test_pgbfs_rounds } ,
 { "pgbfs-thread" , "-" , test_pgbfs< ADS::ThreadPool > } ,
 { "bb-scheme" , "-" , test_bb_scheme } ,
 { "bb-cutoff" , "-" , test_bb_cutoff } ,
 { "bb-values" , "-" , test_bb_values } ,
 { "bb-skel-thread" , "-" , test_bb_skeleton_thread } ,
 { "bb-ivm-thread" , "-" , test_bb_interval_thread } ,
 { "bb-modes" , "-" , test_bb_modes } ,
 { "bb-race" , "-" , test_bb_race } ,
 { "bb-checkpoint" , "-" , test_bb_checkpoint } ,
 { "bb-domains" , "-" , test_bb_domains } ,
 { "moa" , "-" , test_moa } ,
 { "moa-thread" , "-" , test_moa_thread } ,
 #ifdef PBEFS_HAS_OPENMP
  { "hda-openmp" , "openmp" , test_hda< ADS::OpenMPExecutor > } ,
  { "hdlpa-openmp" , "openmp" , test_hdlpa< ADS::OpenMPExecutor > } ,
  { "spa-openmp" , "openmp" , test_spa< ADS::OpenMPExecutor > } ,
  { "mqa-openmp" , "openmp" , test_mqa< ADS::OpenMPExecutor > } ,
  { "pbnf-openmp" , "openmp" , test_pbnf< ADS::OpenMPExecutor > } ,
  { "psdd-openmp" , "openmp" , test_psdd< ADS::OpenMPExecutor > } ,
  { "multi-openmp" , "openmp" , test_multi< ADS::OpenMPExecutor > } ,
  { "tds-openmp" , "openmp" , test_tds< ADS::OpenMPExecutor > } ,
  { "bb-shared-openmp" , "openmp" ,
    test_bb_parallel< 0 , ADS::OpenMPExecutor > } ,
  { "bb-dist-openmp" , "openmp" ,
    test_bb_parallel< 1 , ADS::OpenMPExecutor > } ,
  { "bb-steal-openmp" , "openmp" ,
    test_bb_parallel< 3 , ADS::OpenMPExecutor > } ,
  { "bb-mq-openmp" , "openmp" ,
    test_bb_parallel< 4 , ADS::OpenMPExecutor > } ,
  { "bb-det-openmp" , "openmp" , test_bb_det< ADS::OpenMPExecutor > } ,
  { "bb-split-openmp" , "openmp" , test_bb_split< ADS::OpenMPExecutor > } ,
  { "bb-skel-openmp" , "openmp" , test_bb_skeleton_openmp } ,
  { "bb-ivm-openmp" , "openmp" , test_bb_interval_openmp } ,
  { "bb-farm-openmp" , "openmp" ,
    test_bb_parallel< 2 , ADS::OpenMPExecutor > } ,
  { "farm-openmp" , "openmp" , test_farm< ADS::OpenMPExecutor > } ,
  { "bounded-openmp" , "openmp" , test_bounded< ADS::OpenMPExecutor > } ,
  { "pida-openmp" , "openmp" , test_pida< ADS::OpenMPExecutor > } ,
  { "bidir-openmp" , "openmp" ,
    test_bidirectional_parallel< ADS::OpenMPExecutor > } ,
  { "pem-openmp" , "openmp" , test_pem< ADS::OpenMPExecutor > } ,
  { "pgbfs-openmp" , "openmp" , test_pgbfs< ADS::OpenMPExecutor > } ,
  { "pase-openmp" , "openmp" , test_pase< ADS::OpenMPExecutor > } ,
  { "edge-openmp" , "openmp" , test_edge_parallel< ADS::OpenMPExecutor > } ,
  { "moa-openmp" , "openmp" , test_moa_openmp } ,
  { "planning-openmp" , "openmp" ,
    test_planning_parallel< ADS::OpenMPExecutor > } ,
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
  { "hda-fastflow" , "fastflow" , test_hda< ADS::FastFlowExecutor > } ,
  { "hdlpa-fastflow" , "fastflow" , test_hdlpa< ADS::FastFlowExecutor > } ,
  { "farm-fastflow" , "fastflow" , test_farm< ADS::FastFlowExecutor > } ,
  { "pbnf-fastflow" , "fastflow" , test_pbnf< ADS::FastFlowExecutor > } ,
  { "psdd-fastflow" , "fastflow" , test_psdd< ADS::FastFlowExecutor > } ,
  { "multi-fastflow" , "fastflow" , test_multi< ADS::FastFlowExecutor > } ,
  { "tds-fastflow" , "fastflow" , test_tds< ADS::FastFlowExecutor > } ,
  { "bounded-fastflow" , "fastflow" ,
    test_bounded< ADS::FastFlowExecutor , false > } ,
  { "pida-fastflow" , "fastflow" , test_pida< ADS::FastFlowExecutor > } ,
  { "bidir-fastflow" , "fastflow" ,
    test_bidirectional_parallel< ADS::FastFlowExecutor > } ,
  { "pem-fastflow" , "fastflow" , test_pem< ADS::FastFlowExecutor > } ,
  { "pgbfs-fastflow" , "fastflow" ,
    test_pgbfs< ADS::FastFlowExecutor > } ,
  { "pase-fastflow" , "fastflow" , test_pase< ADS::FastFlowExecutor > } ,
  { "edge-fastflow" , "fastflow" ,
    test_edge_parallel< ADS::FastFlowExecutor > } ,
  { "planning-fastflow" , "fastflow" ,
    test_planning_parallel< ADS::FastFlowExecutor , false > } ,
  { "bb-shared-fastflow" , "fastflow" ,
    test_bb_parallel< 0 , ADS::FastFlowExecutor > } ,
  { "bb-dist-fastflow" , "fastflow" ,
    test_bb_parallel< 1 , ADS::FastFlowExecutor > } ,
  { "bb-steal-fastflow" , "fastflow" ,
    test_bb_parallel< 3 , ADS::FastFlowExecutor > } ,
  { "bb-farm-fastflow" , "fastflow" ,
    test_bb_parallel< 2 , ADS::FastFlowExecutor > } ,
  { "bb-det-fastflow" , "fastflow" , test_bb_det< ADS::FastFlowExecutor > } ,
  { "bb-split-fastflow" , "fastflow" ,
    test_bb_split< ADS::FastFlowExecutor > } ,
  { "bb-skel-fastflow" , "fastflow" , test_bb_skeleton_fastflow } ,
  { "bb-ivm-fastflow" , "fastflow" , test_bb_interval_fastflow } ,
  { "moa-fastflow" , "fastflow" , test_moa_fastflow } ,
 #endif
 #ifdef PBEFS_HAS_SYLVAN
  { "sym-bfs" , "sylvan" , test_sym_bfs } ,
  { "sym-dijkstra" , "sylvan" , test_sym_dijkstra } ,
  { "sym-astar" , "sylvan" , test_sym_astar } ,
  { "sym-sas" , "sylvan" , test_sym_sas } ,
  { "sym-pdb" , "sylvan" , test_sym_pdb } ,
  { "sym-variants" , "sylvan" , test_sym_variants } ,
  { "sym-topk" , "sylvan" , test_sym_topk } ,
  { "sym-par-thread" , "sylvan" , test_sym_par_thread } ,
  #ifdef PBEFS_HAS_OPENMP
   { "sym-par-openmp" , "sylvan+openmp" , test_sym_par_openmp } ,
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
   { "sym-par-fastflow" , "sylvan+fastflow" , test_sym_par_fastflow } ,
  #endif
 #endif
 #ifdef PBEFS_HAS_MPI
  { "mpi" , "mpi" , test_mpi } ,
  { "bb-mpi" , "mpi" , test_bb_mpi } ,
  { "pida-mpi" , "mpi" , test_pida_mpi } ,
  { "bb-hub-mpi" , "mpi" , test_bb_hub_mpi } ,
 #endif
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the cases of the list (test/tests.txt) whose runtime is built are those
// of the tester, and vice versa

static void test_registry( const std::string & list )
{
 // a runtime of the list, or several joined by +
 const auto one = []( const std::string & needs ) {
  if( needs == "-" )
   return( true );
  #ifdef PBEFS_HAS_OPENMP
   if( needs == "openmp" )
    return( true );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
   if( needs == "fastflow" )
    return( true );
  #endif
  #ifdef PBEFS_HAS_MPI
   if( needs == "mpi" )
    return( true );
  #endif
  #ifdef PBEFS_HAS_SYLVAN
   if( needs == "sylvan" )
    return( true );
  #endif
  return( false );
  };
 const auto built = [ & ]( const std::string & needs ) {
  std::istringstream is( needs );
  bool ok = true;
  for( std::string r ; std::getline( is , r , '+' ) ; )
   ok = ok && one( r );
  return( ok );
  };

 std::ifstream in( list );
 check( bool( in ) , "the list " + list + " is read" );
 std::set< std::string > listed;
 for( std::string line ; std::getline( in , line ) ; ) {
  std::istringstream is( line );
  std::string name , needs;
  if( ( is >> name >> needs ) && ( name[ 0 ] != '#' ) && built( needs ) )
   listed.insert( name );
  }
 std::set< std::string > tester;
 for( const Case & c : cases )
  tester.insert( c.name );
 for( const auto & n : listed )
  check( tester.count( n ) == 1 , "the case " + n + " is in the tester" );
 for( const auto & n : tester )
  check( listed.count( n ) == 1 , "the case " + n + " is in the list" );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 if( argc < 2 ) {
  std::cerr << "usage: " << argv[ 0 ] << " <case> [<data dir>]" << std::endl;
  return( 2 );
  }
 const std::string what = argv[ 1 ];
 const std::string data = argc > 2 ? argv[ 2 ] : "../data";
 const auto c = std::find_if( cases.begin() , cases.end() ,
                              [ & ]( const Case & x ) {
                               return( what == x.name );
                               } );
 if( ( what != "registry" ) && ( c == cases.end() ) ) {
  std::cerr << "unknown case " << what << std::endl;
  return( 2 );
  }
 [[maybe_unused]] const bool mpi = ( c != cases.end() ) &&
                                   ( std::string( c->needs ) == "mpi" );

 #ifdef PBEFS_HAS_MPI
  int rank = 0;
  if( mpi ) {
   MPI_Init( & argc , & argv );
   MPI_Comm_rank( MPI_COMM_WORLD , & rank );
   quiet = rank != 0;
   }
 #endif

 try {
  if( what == "registry" )
   test_registry( data );
  else
   c->run( data );
  }
 catch( std::exception & e ) {
  std::cout << "FAILED: exception " << e.what() << std::endl;
  ++failures;
  }

 #ifdef PBEFS_HAS_MPI
  if( mpi ) {
   int all = 0;
   MPI_Allreduce( & failures , & all , 1 , MPI_INT , MPI_SUM ,
                  MPI_COMM_WORLD );
   failures = all;
   MPI_Finalize();
   }
 #endif

 if( ! quiet )
  std::cout << ( failures ? "FAILED: " : "passed: " ) << what << " ("
            << failures << " failures)" << std::endl;
 return( failures ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File test.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
