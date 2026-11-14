/*--------------------------------------------------------------------------*/
/*-------------------------- File AStarBench.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Benchmark of the searches of PBeFS: solves instances of the
 * 15-puzzle, of the 8-puzzle, of the 24-puzzle, of grid pathfinding, of
 * multiple sequence alignment or of the problems of the examples with a
 * chosen search, runtime, work distribution and number of workers, and
 * prints one line of comma-separated values per run.
 *
 * Compiled with PBEFS_BENCH_MPI defined, it is the benchmark of HDAStarMPI,
 * to be started under mpiexec; only the process of rank 0 prints.
 *
 * This file reads the instances and runs them; the searches are compiled by
 * AStarBenchPart.cpp, once for each Domain of AStarBenchDomains.h and each
 * of its variants [see AStarBench.h and AStarBenchSolve.h].
 *
 *     AStarBench [options]
 *
 *     --domain  tiles15 | tiles8 | tiles24 | grid4 | grid8 | board | msa |
 *               nav2d | queens | vacuum | missionaries | jugs | tsp |
 *               romania | pancake | burned | topspin | rubik | hanoi | road
 *                                                          (tiles15)
 *     --file    <instances or map>     korf100.txt, instances50.txt, a
 *                                      Moving AI map, a board, the
 *                                      sequences, the cities, korf10.txt,
 *                                      toh20.txt, a road graph (the
 *                                      name of its .gr without .gr), or
 *                                      the map of nav2d
 *                                      (below)
 *     --ids     <a-b,c,...>            instances of the file (all)
 *     --random  <W>x<H>:<p>:<seed>[:<c>]  a random grid instead of a map,
 *                                      its free cells costing 1 to c (1),
 *               or <seed>[:<n>[:<total>]]  random sequences (below)
 *     --size    <size>                 the size of an example (below)
 *     --points  <sx>,<sy>,<gx>,<gy>    start and goal (bottom corners)
 *     --origin  bottom | top           the row 0 of a board   (bottom)
 *     --matrix  <scores>               the substitution matrix of msa
 *                                                          (PAM250)
 *     --gap     <penalty>              a linear gap penalty of msa (none),
 *               or the x of GAP-x of the pancakes          (0)
 *     --plane   <d> | 1/<k> | auto       the thickness of the hyperplanes
 *                                      of --dist HWD              (auto)
 *     --metric  euclidean | manhattan  the distances of tsp   (euclidean)
 *     --algo    astar | batch | ida | pida | ddd | bfs | dfs | spa | mqa |
 *               multi | hda | pra | tds | pbnf | psdd | bfpsdd | idpsdd |
 *               dynpsdd | pase | epase | pwa | mplp | farm | kbfs | awa |
 *               ara | optimistic | ios | focal | dps | lpa | ra | bidij |
 *               mm | mme | fmm | bae | nbs | dvcbs | gbfhs | bibfs | pnba |
 *               pdij | pmm | pbae | pema | pemra | pemm | pembae | idd |
 *               gbfs | ehc | lama | iw | siw | bfws | pgbfs | hdgbfs |
 *               pea | epea | fringe | rbfs | beam | lrta | hpa   (astar)
 *     --weight  <weight w of the heuristic, or bound; inf for GBFS; that
 *               of batch and of the A* of dynpsdd>
 *                                       (1; inf for the greedy ones)
 *     --dd      0 | 1   duplicate dropping                 (0)
 *     --anytime 0 | 1   anytime search (also of epase)     (0)
 *     --schedule <w0>,<w1>,...   the weights of ARA* (w, 1) and of
 *               ARPBNF with --algo pbnf               (none: wPBNF)
 *     --optimism <weight of optimistic and ios>           (2 w - 1)
 *     --focal   d | h | opt | ees  the order of focal      (d)
 *     --paths   0 | 1   print the paths of an anytime search on stderr (0)
 *     --frontier <units of pida>                           (128 per worker)
 *     --rule    astar | breadth | depth  (the frontier of pida)   (astar)
 *               kpgbfs | puhf | puhf2 | puhf3 | puhf4 | spuhf | obat
 *               (pgbfs)                                    (kpgbfs)
 *     --units   steal | shared           (the units of pida)      (steal)
 *     --async   1 for the asynchronous scheme of pida      (0)
 *     --bounds  0 to search every unit at every bound      (1)
 *     --split   <least slack of a node given away by pida> (0: none)
 *     --table   <levels of pida looked up in the frontier> (0)
 *     --det     1 for the deterministic path of pida       (0)
 *     --pdb     the tables of rubik (korf), pancake, burned and topspin
 *               (<k>: the first k objects, then the next k, ...), hanoi
 *               (<k>+<k>+...: the ranges of disks from the smallest)
 *     --cache   <directory of the tables>                  (none)
 *     --reverse 1 to search hanoi from the goal to the start  (0)
 *     --parts   <files of a bucket of ddd and pem*>        (64)
 *     --dir     <directory of the files of ddd, pem*, idd>
 *                                                   (the temporary one)
 *     --edge    <least edge cost of the bidirectional bounds>
 *                                                   (that of the domain)
 *     --fraction <p of fmm and of the split of gbfhs>      (0.5)
 *     --side    auto | priority | alternate | cardinality  (auto)
 *     --allbounds 0 | 1   the largest of the valid bounds  (0)
 *     --eager   0 | 1   nbs and dvcbs ready at f = LB      (0)
 *     --disk    <bytes of the files of pem* and idd>       (no cap)
 *     --slots   <places of the table of idd>               (2^20)
 *     --eval    eager | lazy (gbfs)                        (eager)
 *     --queues  <queues of gbfs, see below>                (best:0)
 *     --tie     fifo | lifo | lowg | highg | highglifo | random
 *                                                          (fifo)
 *     --boost   <boost of the preferred queues of gbfs>    (0)
 *     --epsilon <probability of epsilon-greedy>            (0.2)
 *     --reopen  0 | 1 (gbfs)                               (0)
 *     --seed    <seed of the random choices>               (1)
 *     --pref    ignore | prune | first (ehc)               (ignore)
 *     --width   <width of iw, siw and bfws>                (2)
 *     --sge     0 | 1 (pgbfs)                              (0)
 *     --rounds  0 | 1 (pgbfs, synchronous rounds)          (0)
 *     --cutoff  <cutoff C of pea>                          (0)
 *     --beam    <width of beam>                            (100)
 *     --lookahead <expansions of lrta, 1 for LRTA*>        (1)
 *     --trials  <trials of lrta, until convergence>        (1)
 *     --cluster <side of the clusters of hpa>              (10)
 *     --runtime thread | openmp | fastflow                 (thread)
 *     --dist    Z | AZ | A | P | H | R | HWD               (Z)
 *     --threads <p>                                        (1)
 *     --batch   <nodes per message, or per cycle of batch> (64)
 *     --flush   <expansions between flushes>               (64)
 *     --open    buckets | heap | lowg (the heap with the ties of f to the
 *               smaller g)                                 (buckets)
 *     --block   <side of the blocks of AZH and AHDA*>   (100, msa 4)
 *     --minexp  <expansions before PBNF or BFPSDD switches nblock> (32)
 *     --safe    0 | 1: plain PBNF or SafePBNF               (1)
 *     --mult    <nblocks of a layer of BFPSDD per worker>   (4)
 *     --bound   <the bound on f of psdd and dfs>            (none)
 *     --sync    0 | 1: the synchronous sends of HDA*        (0)
 *     --recv    lockfree | locked: the inboxes of HDA*      (lockfree)
 *     --lend    <f gap at which HDA* lends a node>          (0, never)
 *     --growth  <growth of the bounds of IDA*-CR, ida, idpsdd, tds> (0)
 *     --optimal 0 | 1: the iteration of the goal completed     (1)
 *     --nblock  <side of the nblocks of PBNF>        (20, msa 16), or for a
 *               grid <w>x<h> (w x h cells) or coarse:<c>x<r> (the map cut
 *               into c x r)
 *     --eps     <inflation of the independence of PA*SE>   (1)
 *     --scan    <states of OPEN examined by pase, epase>   (64)
 *     --gen     <1: GePA*SE>                               (0)
 *     --dispatch <0: epase without its dispatcher>         (1)
 *     --step    <decrease of the weight of anytime epase>  (0.5)
 *     --time    <seconds of anytime epase>                 (0: none)
 *     --slow    <iterations of busy work per expansion>    (0)
 *     --edges   <cheap>:<expensive>:<share>[:<missing>[:<seed>]]
 *                                      edges made slow instead [see
 *                                      SlowEdgeDomain]     (none)
 *     --scale   <scale of the map of nav2d>                (5)
 *     --factors <seed of the cost factors of nav2d>        (0: none)
 *     --diag    <repetitions of the check of the diagonal
 *               motions of nav2d>                          (1)
 *     --limit   <nodes per worker>                         (none)
 *     --expect  <nodes expected per worker>               (unknown)
 *     --repeat  <runs of each instance>                    (1)
 *     --header  print the header line
 *
 * A file of the tiles domains holds one board per line, or a single
 * instance given by the cells of the tiles [see
 * SlidingTile::read_positions()]. The domain board reads a grid in the
 * format of Burns et al. (JAIR 2010) [see GridBoard], with its moves and
 * costs, the row 0 of its start and goal being the last line (bottom) or
 * the first one (top). The domain msa aligns the sequences of a FASTA file,
 * or the random ones of --random <seed>[:<n>[:<total>]] [see
 * SequenceAlignment::random_sequences()], scored by PAM250 or by the 625
 * scores of --matrix [see SequenceAlignment::read_scores()], with the gap
 * scores of the matrix or the linear gap penalty of --gap; --dist HWD is
 * its hyperplane work distribution [see HyperplaneDistribution], of the
 * thickness that --plane gives or that the rule of Kobayashi et al. gives
 * for --threads.
 *
 * The queues of gbfs are a list of kind:h[+h'] separated by commas, kind
 * being best, pref (best on the preferred successors only), type, eps or
 * pareto and h, h' the indices of the heuristics (0 the heuristic of the
 * domain, 1 the second one of the sliding-tile puzzle, the misplaced
 * tiles), e.g. best:0,pref:0 for the alternation of LAMA.
 *
 * The problems of the examples are sized by --size, and the instance of
 * the line is the size, or the size and the id after a slash: queens
 * <n>; vacuum <n>, with all the cells dirty and the robot in each cell
 * of --ids (the middle one by default, the cells counted from 0);
 * missionaries <m>:<c>:<b>; jugs <goal>:<capacity>:...; tsp <n>, with the
 * random cities of each seed of --ids (1 by default), or with the first
 * <n> cities of --file (all of them if --size is not given) [see
 * TravelingSalesman::read_cities()]; romania <from>:<to> (Arad:Bucharest
 * by default). The permutation puzzles are sized by --size too, among the
 * sizes compiled: pancake 12, 16, 20, 32, 64 and burned 8, 10, 12, 16 (a
 * random stack for each seed of --ids), topspin 12, 14, 16, 17 with a
 * turnstile of 4 (a random walk of 1000 moves for each seed); rubik takes
 * the scrambles of --file, or random walks of --size moves for each seed of
 * --ids, and hanoi the lines of --file (12, 16 or 20 disks, by --size, 20
 * by default). A road graph is read from the files <file>.gr and <file>.co,
 * and searched between the nodes of --points <s>,<t>, or between random
 * nodes for each seed of --ids, with the heuristic of the coordinates
 * (--open heap suits its large costs). The instances of nav2d are the pairs
 * of points of --ids (all by default) in the files nav2d_starts.txt and
 * nav2d_goals.txt next to the map of --file. The line ends with the weight,
 * the lower bound on the optimal cost, the number of paths found, the
 * evaluations and the iterations.
 *
 * The search pra is PRA*, i.e., hda with --sync 1 and --recv locked; with
 * --dist A it is APRA*. The search tds takes --dist as hda does; tds and
 * --recv locked run on the threads only, and only hda and farm (kbfs) on
 * FastFlow.
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

#include "AStarBench.h"

/*--------------------------------------------------------------------------*/
/*------------------------ THE RUNS OF THE DOMAINS -------------------------*/
/*--------------------------------------------------------------------------*/
// runs the sliding-tile instances of the file

template< class D >
static void run_tiles( void )
{
 const std::string text = read_file();
 if( text.find( "starting positions" ) != std::string::npos ) {
  std::istringstream is( text );
  const auto [ start , goal ] = D::read_positions( is );
  run_domain( opt[ "file" ] , D( start , goal ) );
  return;
  }
 std::istringstream in( text );
 std::vector< typename D::board_type > boards;
 std::string line;
 while( std::getline( in , line ) ) {
  if( line.empty() || ( line[ 0 ] == '#' ) )
   continue;
  boards.push_back( D::parse( line ) );
  }
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int id : parse_ids( opt[ "ids" ] , int( boards.size() ) ) ) {
  const D d( boards.at( std::size_t( id - 1 ) ) );
  for( int r = 0 ; r < repeat ; ++r )
   print( std::to_string( id ) , solve_maybe_slow( d ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the nblocks of --nblock on the map m

static GridNBlocks grid_nblocks( const GridMap & m )
{
 const std::string & v = opt[ "nblock" ];
 if( v.empty() )
  return( GridNBlocks() );
 int a , b;
 if( std::sscanf( v.c_str() , "coarse:%dx%d" , & a , & b ) == 2 )
  return( GridNBlocks::coarse( m.width() , m.height() , a , b ) );
 if( std::sscanf( v.c_str() , "%dx%d" , & a , & b ) == 2 )
  return( GridNBlocks( a , b ) );
 return( GridNBlocks( std::stoi( v ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a grid instance, from a Moving AI map or random

template< int Moves >
static void run_grid( void )
{
 GridMap m;
 std::string name;
 if( ! opt[ "random" ].empty() ) {
  int w , h , c = 1;
  double p;
  unsigned long long seed;
  if( std::sscanf( opt[ "random" ].c_str() , "%dx%d:%lf:%llu:%d" , & w , & h ,
                   & p , & seed , & c ) < 4 )
   throw( std::invalid_argument( "--random wants WxH:p:seed[:c]" ) );
  m = GridMap::random( w , h , p , seed , c );
  name = opt[ "random" ];
  }
 else {
  std::ifstream in( opt[ "file" ] );
  if( ! in )
   throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
  m = GridMap::read_movingai( in );
  name = opt[ "file" ];
  }
 int sx = 0 , sy = m.height() - 1 , gx = m.width() - 1 , gy = m.height() - 1;
 if( ! opt[ "points" ].empty() &&
     ( std::sscanf( opt[ "points" ].c_str() , "%d,%d,%d,%d" , & sx , & sy ,
                    & gx , & gy ) != 4 ) )
  throw( std::invalid_argument( "--points wants sx,sy,gx,gy" ) );
 const Grid< Moves > g( m , sx , sy , gx , gy , option( "block" , 100 ) ,
                        grid_nblocks( m ) );
 const int repeat = std::stoi( opt[ "repeat" ] );
 #ifndef PBEFS_BENCH_MPI
 if( opt[ "algo" ] == "hpa" ) {
  HPAOptions o;
  o.cluster = std::stoi( opt[ "cluster" ] );
  for( int r = 0 ; r < repeat ; ++r )
   print( name , HPAStar< Moves >( m , o ).solve( sx , sy , gx , gy ) );
  return;
  }
 #endif
 for( int r = 0 ; r < repeat ; ++r )
  print( name , solve_maybe_slow( g ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a grid in the format of Burns et al., with its moves

static void run_board( void )
{
 const std::string & o = opt[ "origin" ];
 if( ( o != "bottom" ) && ( o != "top" ) )
  throw( std::invalid_argument( "--origin wants bottom or top" ) );
 std::istringstream in( read_file() );
 const GridBoard b = GridBoard::read( in , o == "bottom" );
 const int block = option( "block" , 100 );
 const GridNBlocks nb = grid_nblocks( b.map );
 if( b.moves == 4 )
  run_domain( opt[ "file" ] , Grid< 4 >( b , block , nb ) );
 else
  run_domain( opt[ "file" ] , Grid< 8 >( b , block , nb ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the alignment of the sequences of the file

static void run_msa( void )
{
 std::vector< std::string > seqs;
 std::string name = opt[ "file" ];
 if( ! opt[ "random" ].empty() ) {
  unsigned long long seed;
  int n = 0 , total = 5000;
  if( std::sscanf( opt[ "random" ].c_str() , "%llu:%d:%d" , & seed , & n ,
                   & total ) < 1 )
   throw( std::invalid_argument( "--random wants seed[:n[:total]]" ) );
  seqs = SequenceAlignment<>::random_sequences( seed , n , total );
  name = opt[ "random" ];
  }
 else {
  std::istringstream in( read_file() );
  seqs = SequenceAlignment<>::read_fasta( in );
  }
 ScoreMatrix scores = pam250;
 if( ! opt[ "matrix" ].empty() ) {
  std::ifstream ms( opt[ "matrix" ] );
  if( ! ms )
   throw( std::invalid_argument( "cannot read " + opt[ "matrix" ] ) );
  scores = SequenceAlignment<>::read_scores( ms );
  }
 run_domain( name , SequenceAlignment<>(
              seqs , scores , option( "block" , 4 ) ,
              option( "nblock" , 16 ) , std::stoi( opt[ "gap" ] ) ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the 2D navigation on the map of --file and its pairs of points

static void run_nav2d( void )
{
 const std::string & file = opt[ "file" ];
 const int scale = std::stoi( opt[ "scale" ] );
 std::ifstream in( file );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + file ) );
 const NavMap m = NavMap::read_movingai( in , scale );
 const auto slash = file.find_last_of( "/\\" );
 const std::string dir = slash == std::string::npos ? std::string( "" )
                                                    : file.substr( 0 ,
                                                                   slash + 1 );
 std::ifstream fs( dir + "nav2d_starts.txt" ) , fg( dir + "nav2d_goals.txt" );
 const auto pts = RobotNav2D::read_points( fs , fg , scale );
 NavOptions no;
 no.seed = std::stoull( opt[ "factors" ] );
 no.factors = no.seed != 0;
 no.repeat = unsigned( std::stoul( opt[ "diag" ] ) );
 const int repeat = std::stoi( opt[ "repeat" ] );
 for( int id : parse_ids( opt[ "ids" ] , int( pts.size() ) ) ) {
  const auto & q = pts.at( std::size_t( id - 1 ) );
  const RobotNav2D d( m , q[ 0 ] , q[ 1 ] , q[ 2 ] , q[ 3 ] , no );
  for( int r = 0 ; r < repeat ; ++r )
   print( std::to_string( id ) , solve_maybe_slow( d ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the pancakes of N, burned or not, a random stack for each seed

template< int N , bool Burned >
static void run_pancake( void )
{
 const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 );
 for( int seed : ids ) {
  const std::string name = std::to_string( N ) + "/" + std::to_string( seed );
  if constexpr( Burned ) {
   const BurnedPancakePuzzle< N > d( BurnedPancakePuzzle< N >::random(
                                      std::uint64_t( seed ) ) );
   run_permutation( name , d , N );
   }
  else {
   const PancakePuzzle< N > d( PancakePuzzle< N >::random(
                                std::uint64_t( seed ) ) ,
                               unsigned( std::stoul( opt[ "gap" ] ) ) );
   run_permutation( name , d , N );
   }
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs TopSpin( N , 4 ), a random walk of 1000 moves for each seed

template< int N >
static void run_topspin( void )
{
 const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 );
 for( int seed : ids ) {
  const TopSpin< N , 4 > d( TopSpin< N , 4 >::random( std::uint64_t( seed ) ,
                                                      1000 ) );
  run_permutation( std::to_string( N ) + "/" + std::to_string( seed ) , d ,
                   N );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs Rubik's cube: the scrambles of --file, or random walks

static void run_rubik( void )
{
 std::vector< std::pair< std::string , RubiksCube::state_type > > v;
 if( ! opt[ "file" ].empty() ) {
  std::ifstream in( opt[ "file" ] );
  if( ! in )
   throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
  std::vector< RubiksCube::state_type > all;
  for( std::string line ; std::getline( in , line ) ; )
   if( ! line.empty() && ( line[ 0 ] != '#' ) )
    all.push_back( RubiksCube::scramble( line ) );
  for( int id : parse_ids( opt[ "ids" ] , int( all.size() ) ) )
   v.emplace_back( std::to_string( id ) , all.at( std::size_t( id - 1 ) ) );
  }
 else {
  const int moves = std::stoi( size_fields()[ 0 ] );
  const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                        : parse_ids( opt[ "ids" ] , 0 );
  for( int seed : ids )
   v.emplace_back( std::to_string( moves ) + "/" + std::to_string( seed ) ,
                   RubiksCube::random_walk( moves , std::uint64_t( seed ) ) );
  }
 const RubiksCube base;
 std::unique_ptr< RubikPDBHeuristic > h;
 if( opt[ "pdb" ] == "korf" )
  h = std::make_unique< RubikPDBHeuristic >(
       base , unsigned( std::stoul( opt[ "threads" ] ) ) , opt[ "cache" ] );
 else
  if( ! opt[ "pdb" ].empty() )
   throw( std::invalid_argument( "rubik knows --pdb korf only" ) );
 for( const auto & [ name , s ] : v ) {
  const RubiksCube d( s );
  if( h )
   run_domain( name , PDBDomain< RubiksCube , RubikPDBHeuristic >( d , *h ) );
  else
   run_domain( name , d );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Towers of Hanoi of D disks on the lines of --file, with the
// additive tables of --pdb (the ranges of disks, k+k+...)

template< int D >
static void run_hanoi( void )
{
 std::ifstream in( opt[ "file" ] );
 if( ! in )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] ) );
 std::vector< std::pair< std::uint64_t , std::uint64_t > > all;
 for( std::string line ; std::getline( in , line ) ; ) {
  if( line.empty() || ( line[ 0 ] != '(' ) )
   continue;
  std::istringstream is( line );
  const auto s = TowersOfHanoi< D >::read( is );
  all.emplace_back( s , TowersOfHanoi< D >::read( is ) );
  }
 std::vector< unsigned > ks;
 {
  std::istringstream is( opt[ "pdb" ].empty() ? std::to_string( D / 2 ) +
                         "+" + std::to_string( D - D / 2 ) : opt[ "pdb" ] );
  for( std::string part ; std::getline( is , part , '+' ) ; )
   ks.push_back( unsigned( std::stoul( part ) ) );
 }
 for( int id : parse_ids( opt[ "ids" ] , int( all.size() ) ) ) {
  auto [ s , g ] = all.at( std::size_t( id - 1 ) );
  if( opt[ "reverse" ] != "0" )
   std::swap( s , g );
  const TowersOfHanoi< D > d( s , g );
  const HanoiPDB h( ks , g , unsigned( std::stoul( opt[ "threads" ] ) ) ,
                    opt[ "cache" ] );
  run_domain( std::to_string( id ) , PDBDomain< TowersOfHanoi< D > ,
                                               HanoiPDB >( d , h ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a road graph of DIMACS, between --points or random pairs

static void run_road( void )
{
 std::ifstream gr( opt[ "file" ] + ".gr" ) , co( opt[ "file" ] + ".co" );
 if( ! gr || ! co )
  throw( std::invalid_argument( "cannot read " + opt[ "file" ] +
                                ".gr and .co" ) );
 auto g = read_dimacs_graph< std::int64_t >( gr );
 const auto c = read_dimacs_coordinates( co , g.size() );
 const double kappa = dimacs_scale( g , c );
 std::vector< std::pair< std::uint32_t , std::uint32_t > > pairs;
 if( ! opt[ "points" ].empty() ) {
  unsigned long s , t;
  if( std::sscanf( opt[ "points" ].c_str() , "%lu,%lu" , & s , & t ) != 2 )
   throw( std::invalid_argument( "--points wants s,t for a road graph" ) );
  pairs.emplace_back( std::uint32_t( s - 1 ) , std::uint32_t( t - 1 ) );
  }
 else
  for( int seed : opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                       : parse_ids( opt[ "ids" ] , 0 ) ) {
   std::mt19937_64 gen{ std::uint64_t( seed ) };
   const auto s = std::uint32_t( random_below( gen , g.size() ) );
   pairs.emplace_back( s , std::uint32_t( random_below( gen , g.size() ) ) );
   }
 for( const auto & [ s , t ] : pairs ) {
  g.set_source( s );
  g.set_target( t );
  set_dimacs_heuristic( g , c , kappa );
  run_domain( std::to_string( s + 1 ) + ":" + std::to_string( t + 1 ) , g );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the pancakes, the burned pancakes or TopSpin of the size --size

static void run_permutation_sized( void )
{
 const std::string & d = opt[ "domain" ];
 const int n = std::stoi( size_fields()[ 0 ] );
 if( d == "pancake" )
  switch( n ) {
   case 12: run_pancake< 12 , false >(); return;
   case 16: run_pancake< 16 , false >(); return;
   case 20: run_pancake< 20 , false >(); return;
   case 32: run_pancake< 32 , false >(); return;
   case 64: run_pancake< 64 , false >(); return;
   default: throw( std::invalid_argument( "pancake: 12, 16, 20, 32, 64" ) );
   }
 if( d == "burned" )
  switch( n ) {
   case 8: run_pancake< 8 , true >(); return;
   case 10: run_pancake< 10 , true >(); return;
   case 12: run_pancake< 12 , true >(); return;
   case 16: run_pancake< 16 , true >(); return;
   default: throw( std::invalid_argument( "burned: 8, 10, 12, 16" ) );
   }
 switch( n ) {
  case 12: run_topspin< 12 >(); return;
  case 14: run_topspin< 14 >(); return;
  case 16: run_topspin< 16 >(); return;
  case 17: run_topspin< 17 >(); return;
  default: throw( std::invalid_argument( "topspin: 12, 14, 16, 17" ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the Towers of Hanoi of --size disks (20 by default)

static void run_hanoi_sized( void )
{
 switch( opt[ "size" ].empty() ? 20 : std::stoi( opt[ "size" ] ) ) {
  case 12: run_hanoi< 12 >(); return;
  case 16: run_hanoi< 16 >(); return;
  case 20: run_hanoi< 20 >(); return;
  default: throw( std::invalid_argument( "hanoi: 12, 16, 20 disks" ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs a problem of the examples, sized by --size

static void run_example( void )
{
 const std::string & dom = opt[ "domain" ];
 if( ( dom == "tsp" ) && ! opt[ "file" ].empty() ) {
  const std::string & m = opt[ "metric" ];
  if( ( m != "euclidean" ) && ( m != "manhattan" ) )
   throw( std::invalid_argument( "--metric wants euclidean or manhattan" ) );
  std::istringstream in( read_file() );
  run_domain( opt[ "file" ] + "/" + opt[ "size" ] , TravelingSalesman(
               TravelingSalesman::read_cities(
                in , m == "euclidean" ? TravelingSalesman::Metric::euclidean
                                     : TravelingSalesman::Metric::manhattan ,
                option( "size" , 0 ) ) ) );
  return;
  }
 if( dom == "romania" ) {
  if( opt[ "size" ].empty() )
   opt[ "size" ] = "Arad:Bucharest";
  const auto f = size_fields();
  const auto from = Romania::find( f[ 0 ] );
  const auto to = Romania::find( f.size() > 1 ? f[ 1 ] : "Bucharest" );
  if( ( from == Romania::NCities ) || ( to == Romania::NCities ) )
   throw( std::invalid_argument( "unknown city in " + opt[ "size" ] ) );
  run_domain( opt[ "size" ] , Romania( from , to ) );
  return;
  }
 const auto f = size_fields();
 const int n = std::stoi( f[ 0 ] );
 if( dom == "queens" )
  run_domain( opt[ "size" ] , NQueens( n ) );
 else
  if( dom == "vacuum" ) {
   const std::uint32_t dirty = n >= 32 ? ~std::uint32_t( 0 )
                                       : ( std::uint32_t( 1 ) << n ) - 1;
   const auto ids = opt[ "ids" ].empty()
                    ? std::vector< int >{ n / 2 } : parse_ids( opt[ "ids" ] ,
                                                               0 );
   for( int start : ids )
    run_domain( opt[ "size" ] + "/" + std::to_string( start ) ,
                VacuumWorld( n , start , dirty ) );
   }
  else
   if( dom == "missionaries" ) {
    if( f.size() != 3 )
     throw( std::invalid_argument( "--size wants <m>:<c>:<b>" ) );
    run_domain( opt[ "size" ] , MissionariesCannibals(
                 n , std::stoi( f[ 1 ] ) , std::stoi( f[ 2 ] ) ) );
    }
   else
    if( dom == "jugs" ) {
     std::vector< int > cap;
     for( std::size_t i = 1 ; i < f.size() ; ++i )
      cap.push_back( std::stoi( f[ i ] ) );
     run_domain( opt[ "size" ] , WaterJugs( cap , n ) );
     }
    else
     if( dom == "tsp" ) {
      const auto ids = opt[ "ids" ].empty() ? std::vector< int >{ 1 }
                                           : parse_ids( opt[ "ids" ] , 0 );
      for( int seed : ids )
       run_domain( opt[ "size" ] + "/" + std::to_string( seed ) ,
                   TravelingSalesman( TravelingSalesman::random_cities(
                                       n , std::uint64_t( seed ) ) ) );
      }
     else
      throw( std::invalid_argument( "unknown domain " + dom ) );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char * argv[] )
{
 #ifdef PBEFS_BENCH_MPI
  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  opt[ "algo" ] = "hda";
  opt[ "runtime" ] = "mpi";
 #endif

 int ret = 0;
 try {
  bool header = false;
  for( int i = 1 ; i < argc ; ++i ) {
   std::string a = argv[ i ];
   if( a == "--header" ) {
    header = true;
    continue;
    }
   if( ( a.rfind( "--" , 0 ) != 0 ) || ! opt.count( a.substr( 2 ) ) ||
       ( i + 1 >= argc ) )
    throw( std::invalid_argument( "bad option " + a ) );
   opt[ a.substr( 2 ) ] = argv[ ++i ];
   }
  #ifdef PBEFS_BENCH_MPI
   int np;
   MPI_Comm_size( MPI_COMM_WORLD , & np );
   opt[ "threads" ] = std::to_string( np );
  #endif
  if( header && ( rank == 0 ) )
   std::cout << "domain,instance,algo,runtime,dist,threads,batch,flush,open,"
                "status,cost,expanded,reexpanded,generated,duplicates,sent,"
                "messages,stored,LB,CO,seconds,weight,bound,paths,"
                "evaluated,iterations" << std::endl;

  const std::string & d = opt[ "domain" ];
  if( d == "tiles15" )
   run_tiles< FifteenPuzzle >();
  else
   if( d == "tiles8" )
    run_tiles< EightPuzzle >();
   else
    if( d == "tiles24" )
     run_tiles< TwentyFourPuzzle >();
    else
     if( d == "grid4" )
      run_grid< 4 >();
     else
      if( d == "grid8" )
       run_grid< 8 >();
      else
       if( d == "board" )
        run_board();
       else
        if( d == "msa" )
         run_msa();
        else
         if( ( d == "pancake" ) || ( d == "burned" ) || ( d == "topspin" ) )
          run_permutation_sized();
         else
          if( d == "rubik" )
           run_rubik();
          else
           if( d == "hanoi" )
            run_hanoi_sized();
           else
            if( d == "road" )
             run_road();
            else
             if( d == "nav2d" )
              run_nav2d();
             else
              run_example();
  }
 catch( std::exception & e ) {
  std::cerr << "AStarBench: " << e.what() << std::endl;
  ret = 2;
  }

 #ifdef PBEFS_BENCH_MPI
  MPI_Finalize();
 #endif
 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File AStarBench.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
