/*--------------------------------------------------------------------------*/
/*------------------------- File AStarBenchSolve.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The searches of AStarBench, chosen by the options, and
 * solve_maybe_slow(), which runs the one of the command line on a Domain;
 * AStarBenchPart.cpp instantiates them for a Domain [see AStarBench.h].
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

#ifndef __PBeFS_AStarBenchSolve
 #define __PBeFS_AStarBenchSolve
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AStarBench.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
// the parallel options from the command line

inline ParallelOptions parallel_options( void )
{
 ParallelOptions o;
 o.threads = unsigned( std::stoul( opt[ "threads" ] ) );
 o.batch = std::stoul( opt[ "batch" ] );
 o.flush_period = unsigned( std::stoul( opt[ "flush" ] ) );
 o.node_limit = std::stoul( opt[ "limit" ] );
 o.expected_nodes = std::stoul( opt[ "expect" ] );
 o.goal_at_generation = opt[ "gengoal" ].empty() ? opt[ "domain" ] ==
                        "planning" : opt[ "gengoal" ] == "1";
 o.weight = std::stod( weight_option() );
 o.drop_duplicates = opt[ "dd" ] == "1";
 o.anytime = opt[ "anytime" ] == "1";
 return( o );
 }

#ifndef PBEFS_BENCH_MPI
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the breaking of the ties from the command line

inline TieBreak tie_break( void )
{
 const std::string & t = opt[ "tie" ];
 if( t == "fifo" )
  return( TieBreak::fifo );
 if( t == "lifo" )
  return( TieBreak::lifo );
 if( t == "lowg" )
  return( TieBreak::low_g );
 if( t == "highg" )
  return( TieBreak::high_g );
 if( t == "highglifo" )
  return( TieBreak::high_g_lifo );
 if( t == "random" )
  return( TieBreak::random );
 throw( std::invalid_argument( "unknown --tie " + t ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of gbfs from the command line

inline GreedyOptions greedy_options( const ParallelOptions & p )
{
 GreedyOptions o;
 o.node_limit = p.node_limit;
 o.expected_nodes = p.expected_nodes;
 o.evaluation = opt[ "eval" ] == "lazy" ? GreedyEvaluation::lazy
                                        : GreedyEvaluation::eager;
 o.weight = std::stod( weight_option() );
 o.tie = tie_break();
 o.boost = std::stoi( opt[ "boost" ] );
 o.epsilon = std::stod( opt[ "epsilon" ] );
 o.reopen = opt[ "reopen" ] == "1";
 o.seed = std::stoull( opt[ "seed" ] );
 o.queues.clear();
 std::istringstream in( opt[ "queues" ] );
 for( std::string q ; std::getline( in , q , ',' ) ; ) {
  const auto colon = q.find( ':' );
  const std::string kind = q.substr( 0 , colon );
  QueueSpec x;
  x.heuristics.clear();
  std::istringstream hs( colon == std::string::npos ? "0"
                                                    : q.substr( colon + 1 ) );
  for( std::string h ; std::getline( hs , h , '+' ) ; )
   x.heuristics.push_back( std::stoul( h ) );
  if( kind == "best" )
   x.kind = QueueKind::best;
  else
   if( kind == "pref" ) {
    x.kind = QueueKind::best;
    x.preferred_only = true;
    }
   else
    if( kind == "type" )
     x.kind = QueueKind::type;
    else
     if( kind == "eps" )
      x.kind = QueueKind::epsilon;
     else
      if( kind == "pareto" )
       x.kind = QueueKind::pareto;
      else
       throw( std::invalid_argument( "unknown queue " + q ) );
  o.queues.push_back( x );
  }
 return( o );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the rule of pgbfs from the command line

inline GreedyRule greedy_rule( void )
{
 const char * n[] = { "kpgbfs" , "puhf" , "puhf2" , "puhf3" , "puhf4" ,
                      "spuhf" , "obat" };
 const std::string & r = opt[ "rule" ];
 if( r.empty() )
  return( GreedyRule( 0 ) );
 for( int i = 0 ; i < 7 ; ++i )
  if( r == n[ i ] )
   return( GreedyRule( i ) );
 throw( std::invalid_argument( "unknown --rule " + opt[ "rule" ] ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with ParallelGreedy on the chosen runtime

template< class D >
inline Solution< state_t< D > , cost_t< D > > run_pgbfs(
 const D & d , const ParallelOptions & p )
{
 ParallelGreedyOptions o;
 static_cast< ParallelOptions & >( o ) = p;
 o.rule = greedy_rule();
 o.sge = opt[ "sge" ] == "1";
 o.rounds = opt[ "rounds" ] == "1";
 o.tie = tie_break();
 const std::string & r = opt[ "runtime" ];
 if( r == "thread" )
  return( ParallelGreedy< D , ADS::ThreadPool >( d , o ).solve() );
 #ifdef PBEFS_HAS_OPENMP
 if( r == "openmp" )
  return( ParallelGreedy< D , ADS::OpenMPExecutor >( d , o ).solve() );
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
 if( r == "fastflow" )
  return( ParallelGreedy< D , ADS::FastFlowExecutor >( d , o ).solve() );
 #endif
 throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with HDGBFS and the distribution Dist on the chosen runtime

template< class D , class Dist >
inline Solution< state_t< D > , cost_t< D > > run_hdgbfs(
 const D & d , const Dist & dist , const ParallelOptions & p )
{
 ParallelGreedyOptions o;
 static_cast< ParallelOptions & >( o ) = p;
 o.tie = tie_break();
 const std::string & r = opt[ "runtime" ];
 if( r == "thread" )
  return( HDGBFS< D , Dist , ADS::ThreadPool >( d , dist , o ).solve() );
 #ifdef PBEFS_HAS_OPENMP
 if( r == "openmp" )
  return( HDGBFS< D , Dist , ADS::OpenMPExecutor >( d , dist , o ).solve() );
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
 if( r == "fastflow" )
  return( HDGBFS< D , Dist , ADS::FastFlowExecutor >( d , dist , o )
          .solve() );
 #endif
 throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with the greedy, width-based, partial-expansion, real-time and
// other searches of the command line; false in done if algo is none of
// them

template< class D >
inline Solution< state_t< D > , cost_t< D > > solve_other(
 const D & d , const ParallelOptions & p , bool & done )
{
 const std::string & algo = opt[ "algo" ];
 done = true;
 if( algo == "gbfs" )
  return( gbfs( d , greedy_options( p ) ) );
 if( algo == "ehc" ) {
  EHCOptions o;
  o.node_limit = p.node_limit;
  o.preferred = opt[ "pref" ] == "prune" ? PreferredUsage::prune
              : opt[ "pref" ] == "first" ? PreferredUsage::first
                                         : PreferredUsage::ignore;
  o.fallback = true;
  o.greedy = greedy_options( p );
  return( ehc( d , o ) );
  }
 if( algo == "lama" ) {
  IteratedOptions o = lama_phases( heuristics_of( d ) ,
                                   PreferredDomain< D > );
  for( GreedyOptions & x : o.phases )
   x.node_limit = p.node_limit;
  return( IteratedSearch< D >( d , o ).solve() );
  }
 if( ( algo == "iw" ) || ( algo == "siw" ) || ( algo == "bfws" ) ) {
  WidthOptions o;
  o.node_limit = p.node_limit;
  o.width = unsigned( std::stoul( opt[ "width" ] ) );
  if constexpr( FeatureDomain< D > ) {
   if( algo == "iw" )
    return( iw( d , o ) );
   if( algo == "bfws" )
    return( bfws( d , o ) );
   if constexpr( GoalFeatureDomain< D > )
    return( siw( d , o ) );
   throw( std::invalid_argument( "no goal features for this domain" ) );
   }
  throw( std::invalid_argument( "no features for this domain" ) );
  }
 if( algo == "pgbfs" )
  return( run_pgbfs( d , p ) );
 if( algo == "hdgbfs" ) {
  const std::string & dist = opt[ "dist" ];
  if constexpr( FeatureDomain< D > )
   if( dist == "Z" )
    return( run_hdgbfs( d , ZobristDistribution< D >( d ) , p ) );
  if constexpr( AbstractFeatureDomain< D > )
   if( dist == "AZ" )
    return( run_hdgbfs( d , AbstractZobristDistribution< D >( d ) , p ) );
  if( dist == "H" )
   return( run_hdgbfs( d , StateHashDistribution< D >( d ) , p ) );
  throw( std::invalid_argument( "distribution " + dist +
                                " not available for hdgbfs" ) );
  }
 if( algo == "pea" ) {
  PartialOptions o;
  o.node_limit = p.node_limit;
  o.cutoff = std::stod( opt[ "cutoff" ] );
  return( pea_star( d , o ) );
  }
 if( algo == "epea" ) {
  if constexpr( PartialExpansionDomain< D > )
   return( epea_star( d , p ) );
  throw( std::invalid_argument( "no operator selection for this domain" ) );
  }
 if( algo == "fringe" )
  return( fringe_search( d , p ) );
 if( algo == "rbfs" )
  return( rbfs( d , p ) );
 if( algo == "beam" ) {
  BeamOptions o;
  o.node_limit = p.node_limit;
  o.width = std::stoul( opt[ "beam" ] );
  return( beam_search( d , o ) );
  }
 if( algo == "lrta" ) {
  RealTimeOptions o;
  o.lookahead = std::stoul( opt[ "lookahead" ] );
  o.trials = std::stoul( opt[ "trials" ] );
  o.until_convergence = true;
  return( lss_lrta_star( d , o ) );
  }
 done = false;
 return( Solution< state_t< D > , cost_t< D > >() );
 }
#endif
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of ParallelIDAStar from the command line

inline PIDAOptions pida_options( void )
{
 PIDAOptions o;
 static_cast< ParallelOptions & >( o ) = parallel_options();
 o.frontier = std::stoul( opt[ "frontier" ] );
 const std::string r = opt[ "rule" ].empty() ? "astar" : opt[ "rule" ];
 if( ( r != "astar" ) && ( r != "breadth" ) && ( r != "depth" ) )
  throw( std::invalid_argument( "unknown frontier rule " + r ) );
 o.rule = r == "astar" ? FrontierRule::astar : r == "breadth"
                       ? FrontierRule::breadth : FrontierRule::depth;
 o.schedule = opt[ "units" ] == "shared" ? UnitSchedule::shared
                                            : UnitSchedule::steal;
 o.asynchronous = opt[ "async" ] != "0";
 o.unit_bounds = opt[ "bounds" ] != "0";
 o.split = std::stod( opt[ "split" ] );
 o.table_depth = unsigned( std::stoul( opt[ "table" ] ) );
 o.deterministic = opt[ "det" ] != "0";
 return( o );
 }

#ifndef PBEFS_BENCH_MPI
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the options of the bidirectional searches from the command line

inline BidirectionalOptions bidirectional_options( void )
{
 BidirectionalOptions o;
 static_cast< SearchOptions & >( o ) = parallel_options();
 o.epsilon = std::stod( opt[ "edge" ] );
 o.fraction = std::stod( opt[ "fraction" ] );
 o.all_bounds = opt[ "allbounds" ] == "1";
 o.eager = opt[ "eager" ] == "1";
 const std::string & side = opt[ "side" ];
 o.side = side == "priority" ? SideChoice::priority :
          side == "alternate" ? SideChoice::alternate :
          side == "cardinality" ? SideChoice::cardinality :
                                  SideChoice::automatic;
 return( o );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with a bidirectional or an external search, if algo is one of
// them, setting done

template< class D , ADS::WorkerExecutor Executor >
inline Solution< state_t< D > , cost_t< D > > solve_bidirectional(
 const D & d , const std::string & algo , bool & done )
{
 done = true;
 const BidirectionalOptions b = bidirectional_options();
 const unsigned p = unsigned( std::stoul( opt[ "threads" ] ) );
 if constexpr( DiskDomain< D > ) {
  PEMOptions x;
  x.partitions = std::stoul( opt[ "parts" ] );
  x.directory = opt[ "dir" ];
  x.disk_limit = std::stoull( opt[ "disk" ] );
  if( algo == "pema" )
   return( PEMSearch< D , PEMKind::astar , Executor >( d , p , b , x )
           .solve() );
  if( algo == "idd" ) {
   IDDOptions y;
   y.directory = opt[ "dir" ];
   y.disk_limit = std::stoull( opt[ "disk" ] );
   y.slots = std::stoul( opt[ "slots" ] );
   return( idd_astar( d , b , y ) );
   }
  if constexpr( BidirectionalDomain< D > ) {
   if( algo == "pemm" )
    return( PEMSearch< D , PEMKind::mm , Executor >( d , p , b , x )
            .solve() );
   if( algo == "pembae" )
    return( PEMSearch< D , PEMKind::bae , Executor >( d , p , b , x )
            .solve() );
   if( algo == "pemra" ) {
    const ReverseDomain< D > r( d );
    auto s = PEMSearch< ReverseDomain< D > , PEMKind::astar , Executor >(
              r , p , b , x ).solve();
    std::reverse( s.path.begin() , s.path.end() );
    return( s );
    }
   }
  }
 if constexpr( BidirectionalDomain< D > ) {
  if( algo == "ra" ) {
   auto s = astar( ReverseDomain< D >( d ) , OpenKind::automatic , b );
   std::reverse( s.path.begin() , s.path.end() );
   return( s );
   }
  if( algo == "bidij" )
   return( bidijkstra( d , b ) );
  if( algo == "mm" )
   return( mm( d , b ) );
  if( algo == "mme" )
   return( mme( d , b ) );
  if( algo == "fmm" )
   return( fmm( d , b ) );
  if( algo == "bae" )
   return( bae( d , b ) );
  if( algo == "nbs" )
   return( nbs( d , b ) );
  if( algo == "dvcbs" )
   return( dvcbs( d , b ) );
  if constexpr( std::is_integral_v< cost_t< D > > ) {
   if( algo == "gbfhs" )
    return( gbfhs( d , b ) );
   if( algo == "bibfs" )
    return( bibfs( d , b ) );
   }
  if( algo == "pnba" )
   return( ParallelBidirectional< D , ParallelBiKind::nba , Executor >(
            d , b ).solve() );
  if( algo == "pdij" )
   return( ParallelBidirectional< D , ParallelBiKind::dijkstra ,
                                  Executor >( d , b ).solve() );
  if( algo == "pmm" )
   return( ParallelBidirectional< D , ParallelBiKind::mm , Executor >(
            d , b ).solve() );
  if( algo == "pbae" )
   return( ParallelBidirectional< D , ParallelBiKind::bae , Executor >(
            d , b ).solve() );
  }
 done = false;
 return( Solution< state_t< D > , cost_t< D > >() );
 }

#endif  /* PBEFS_BENCH_MPI */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// f.template operator()< E >() for the Executor E of the chosen runtime,
// FastFlow only if All, so that the searches that are not measured on it
// are not compiled for it

#ifndef PBEFS_BENCH_MPI
template< bool All = false , class F >
inline auto on_runtime( F && f )
{
 const std::string & r = opt[ "runtime" ];
 if( r == "thread" )
  return( f.template operator()< ADS::ThreadPool >() );
 #ifdef PBEFS_HAS_OPENMP
 if( r == "openmp" )
  return( f.template operator()< ADS::OpenMPExecutor >() );
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
 if constexpr( All )
  if( r == "fastflow" )
   return( f.template operator()< ADS::FastFlowExecutor >() );
 #endif
 throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 }
#endif
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with HDA* (or TDS) and the distribution Dist on the chosen
// runtime, with the inboxes and the sends of the command line

template< class D , class Dist >
inline auto run_hda( const D & d , const Dist & dist ,
                     const ParallelOptions & o )
{
 #ifdef PBEFS_BENCH_MPI
  return( HDAStarMPI< D , Dist >( d , dist , o ).solve() );
 #else
  HDAOptions ho( o );
  ho.sync_send = opt[ "sync" ] != "0";
  ho.outsource = std::stod( opt[ "lend" ] );
  using NodeT = Node< state_t< D > , cost_t< D > >;
  using Open = std::conditional_t< std::is_integral_v< cost_t< D > > ,
                                   BucketOpenList< NodeT > ,
                                   HeapOpenList< NodeT > >;
  const bool locked = opt[ "recv" ] == "locked";
  if( ! locked && ( opt[ "recv" ] != "lockfree" ) )
   throw( std::invalid_argument( "--recv wants lockfree or locked" ) );
  if( ( opt[ "algo" ] == "tds" ) || locked ) {        // on the threads only
   if( opt[ "runtime" ] != "thread" )
    throw( std::invalid_argument( "tds and --recv locked run on the threads"
                                  " only" ) );
   using E = ADS::ThreadPool;
   if( locked )
    return( HDAStar< D , Dist , E , Open , LockedChannelsOf< E > >(
             d , dist , ho ).solve() );
   TDSOptions to( o );
   to.growth = std::stod( opt[ "growth" ] );
   to.optimal = opt[ "optimal" ] != "0";
   return( TDS< D , Dist , E >( d , dist , to ).solve() );
   }
  return( on_runtime< true >( [ & ]< class E >() {
   return( HDAStar< D , Dist , E , Open >( d , dist , ho ).solve() );
   } ) );
 #endif
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the tiles of a list "a-b,c,..." of --dist tile:<regions>:<tiles>

inline std::vector< int > parse_tiles( const std::string & s )
{
 return( s.empty() ? std::vector< int >() : parse_ids( s , 0 ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the regions of a name of --dist tile:<regions>

inline TileRegions tile_regions( const std::string & s )
{
 static const std::map< std::string , TileRegions > r = {
  { "single" , TileRegions::single } , { "pair" , TileRegions::pair } ,
  { "line" , TileRegions::line } , { "block" , TileRegions::block } ,
  { "half" , TileRegions::half } , { "quad" , TileRegions::quad } ,
  { "cross" , TileRegions::cross } , { "odd" , TileRegions::odd } };
 const auto it = r.find( s );
 if( it == r.end() )
  throw( std::invalid_argument( "unknown regions " + s ) );
 return( it->second );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true for the domains on which the distributions built on a sample are
// offered: those with a feature projection or with feature groups, i.e.,
// the puzzles and the grids compared by Jinnai and Fukunaga, and not
// slowed down by --slow or --edges (the other cases would only add
// instantiations of HDA* to the build)

template< class D >
inline constexpr bool sampled_domain = AbstractFeatureDomain< D > ||
                                       FeatureGroupDomain< D >;

template< class D >
inline constexpr bool sampled_domain< SlowDomain< D > > = false;

template< class D >
inline constexpr bool sampled_domain< SlowEdgeDomain< D > > = false;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// run_hda(), with the distribution estimated on a sample of d first if
// --estimate is 1 and d is one of the sampled_domain

template< class D , class Dist >
inline auto run_estimated( const D & d , const Dist & dist ,
                           const ParallelOptions & o )
{
 if constexpr( sampled_domain< D > )
  if( opt[ "estimate" ] == "1" ) {
   const FeatureModel< D > m( d , std::stoul( opt[ "sample" ] ) );
   extra.estimated = true;
   extra.est = m.estimate( dist , o.threads , std::stod( opt[ "cost" ] ) );
   }
 return( run_hda( d , dist , o ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with HDA* and the distribution of --dist built on a sample of d,
// estimating it on the sample; false if --dist is not one of them

template< class D >
inline bool solve_sampled( const D & d , const ParallelOptions & o ,
                           Solution< state_t< D > , cost_t< D > > & sol )
{
 const std::string & dist = opt[ "dist" ];
 const std::string base = dist.substr( 0 , dist.find( ':' ) );
 const std::vector< std::string > built = { "GAZ" , "GRAZ" , "GRAZS" ,
                                            "GRAZE" , "FAZ" , "OZ" , "OZP" ,
                                            "DAHDA" , "AHDA" , "AUTO" ,
                                            "tile" };
 if( std::find( built.begin() , built.end() , base ) == built.end() )
  return( false );
 if constexpr( sampled_domain< D > ) {
  const FeatureModel< D > m( d , std::stoul( opt[ "sample" ] ) );
  const double c = std::stod( opt[ "cost" ] );
  const auto run = [ & ]( const TableZobristDistribution< D > & x ) {
   extra.estimated = true;
   extra.est = m.estimate( x , o.threads , c );
   #ifdef PBEFS_BENCH_MPI
    sol = run_hda( d , x , o );
   #else
    // on the C++ threads only, which keeps the instantiations of HDA*,
    // and the memory of the build, small
    if( opt[ "runtime" ] != "thread" )
     throw( std::invalid_argument( "the distributions built on a sample "
                                   "run on --runtime thread only" ) );
    sol = HDAStar< D , TableZobristDistribution< D > , ADS::ThreadPool >(
           d , x , o ).solve();
   #endif
   };
  if( base == "AUTO" ) {
   const auto x = auto_distribution( m , o.threads , c );
   extra.chosen = x.name();
   run( x );
   return( true );
   }
  if( base != "tile" ) {
   run( make_distribution( m , dist , c ) );
   return( true );
   }
  if constexpr( requires { D::tile_abstraction( TileRegions::single ); } ) {
   const auto c1 = dist.find( ':' ) , c2 = dist.find( ':' , c1 + 1 );
   const auto a = FeatureAbstraction::of( D::tile_abstraction(
    tile_regions( dist.substr( c1 + 1 , c2 - c1 - 1 ) ) ,
    parse_tiles( c2 == std::string::npos ? "" : dist.substr( c2 + 1 ) ) ) );
   run( TableZobristDistribution< D >( d , a , 1 , dist ) );
   return( true );
   }
  }
 throw( std::invalid_argument( "distribution " + dist +
                               " not available for this domain" ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the bound of --bound, the infinity if none

template< class Cost >
inline Cost bound_of( void )
{
 if( opt[ "bound" ].empty() )
  return( infinite_cost< Cost >() );
 return( Cost( std::stod( opt[ "bound" ] ) ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// runs the search S< D , Executor > on the runtime of the command line

#ifndef PBEFS_BENCH_MPI
template< template< class , class > class S , class D , class O >
inline auto run_on( const D & d , const O & o )
{
 const std::string & r = opt[ "runtime" ];
 if( r == "thread" )
  return( S< D , ADS::ThreadPool >( d , o ).solve() );
 #ifdef PBEFS_HAS_OPENMP
 if( r == "openmp" )
  return( S< D , ADS::OpenMPExecutor >( d , o ).solve() );
 #endif
 #ifdef PBEFS_HAS_FASTFLOW
 if( r == "fastflow" )
  return( S< D , ADS::FastFlowExecutor >( d , o ).solve() );
 #endif
 throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
 }

template< class D , class E > using EPASEOf = EPASE< D , E >;
template< class D , class E > using PWAStarOf = PWAStar< D , E >;
template< class D , class E > using MPLPOf = MPLP< D , E >;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with a search by edges

template< class D >
inline Solution< state_t< D > , cost_t< D > > solve_edges(
                                 const D & d , const ParallelOptions & o )
{
 const std::string & algo = opt[ "algo" ];
 if constexpr( EdgeDomain< D > ) {
  if( algo == "pwa" ) {
   PWAStarOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   return( run_on< PWAStarOf >( d , po ) );
   }
  if constexpr( PairwiseDomain< D > )
   if( algo == "epase" ) {
    EPASEOptions eo;
    static_cast< ParallelOptions & >( eo ) = o;
    eo.eps = std::stod( opt[ "eps" ] );
    eo.candidates = std::stoul( opt[ "scan" ] );
    eo.generalized = opt[ "gen" ] != "0";
    eo.dispatcher = opt[ "dispatch" ] != "0";
    eo.weight_step = std::stod( opt[ "step" ] );
    eo.time_limit = std::stod( opt[ "time" ] );
    return( run_on< EPASEOf >( d , eo ) );
    }
  if constexpr( LazyEdgeDomain< D > )
   if( algo == "mplp" ) {
    MPLPOptions mo;
    static_cast< ParallelOptions & >( mo ) = o;
    return( run_on< MPLPOf >( d , mo ) );
    }
  }
 throw( std::invalid_argument( algo + " not available for this domain "
                               "(see --edges)" ) );
 }
#endif
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d with the search, runtime and distribution of the command line;
// not inline, so that a unit can leave a variant of a Domain to another
// with an explicit instantiation declaration [see AStarBenchPart.cpp]

template< class D >
Solution< state_t< D > , cost_t< D > > solve( const D & d )
{
 const std::string & algo = opt[ "algo" ];
 const ParallelOptions o = parallel_options();
 extra = {};
 if( ( opt[ "so" ] == "1" ) && ( algo != "astar" ) )
  extra.sequential = astar( d , opt[ "open" ] == "heap" ? OpenKind::heap
                                                       : OpenKind::buckets ,
                            o ).stats.expanded;

 #ifdef PBEFS_BENCH_MPI
 if( algo == "pida" ) {
  if constexpr( std::is_trivially_copyable_v< state_t< D > > )
   return( ParallelIDAStarMPI< D >( d , pida_options() ).solve() );
  else
   throw( std::invalid_argument( "pida over MPI needs states that are "
                                 "trivially copyable" ) );
  }
 #else
 {
  bool done = false;
  const std::string & r = opt[ "runtime" ];
  Solution< state_t< D > , cost_t< D > > s;
  if( r == "thread" )
   s = solve_bidirectional< D , ADS::ThreadPool >( d , algo , done );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   s = solve_bidirectional< D , ADS::OpenMPExecutor >( d , algo , done );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   s = solve_bidirectional< D , ADS::FastFlowExecutor >( d , algo , done );
  #endif
  if( done )
   return( s );
  }
 {
  bool done = false;
  auto x = solve_other( d , o , done );
  if( done )
   return( x );
  }
 if( algo == "ida" ) {
  IDAOptions io( o );
  io.growth = std::stod( opt[ "growth" ] );
  io.optimal = opt[ "optimal" ] != "0";
  return( idastar( d , io ) );
  }
 if( algo == "pida" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( ParallelIDAStar< D , ADS::ThreadPool >( d , pida_options() )
           .solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( ParallelIDAStar< D , ADS::OpenMPExecutor >( d , pida_options() )
           .solve() );
  #endif
  #ifdef PBEFS_HAS_FASTFLOW
  if( r == "fastflow" )
   return( ParallelIDAStar< D , ADS::FastFlowExecutor >( d , pida_options() )
           .solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "ddd" ) {
  if constexpr( DiskDomain< D > ) {
   DDDOptions x;
   x.partitions = std::stoul( opt[ "parts" ] );
   x.directory = opt[ "dir" ];
   return( ddd_astar( d , o , x ) );
   }
  else
   throw( std::invalid_argument( "ddd needs integer costs and states that"
                                 " are trivially copyable" ) );
  }
 const OpenKind ok = opt[ "open" ] == "heap" ? OpenKind::heap
                     : opt[ "open" ] == "lowg" ? OpenKind::heap_low_g
                                               : OpenKind::buckets;
 if( algo == "astar" )
  return( astar( d , ok , o ) );
 if( algo == "batch" ) {
  BatchOptions bo;
  static_cast< SearchOptions & >( bo ) = o;
  bo.batch = o.batch;
  return( BatchAStar< D >( d , bo ).solve() );
  }
 if( algo == "bfs" )
  return( breadth_first( d , o ) );
 if( algo == "dfs" )
  return( cost_bounded_dfs( d , bound_of< cost_t< D > >() , o ) );
 if( algo == "multi" )
  return( on_runtime( [ & ]< class E >() {
   return( MultiAStar< D , E >( d , o ).solve() );
   } ) );
 if( ( algo == "awa" ) || ( algo == "ara" ) ) {
  AnytimeOptions x;
  static_cast< SearchOptions & >( x ) = o;
  if( algo == "ara" ) {
   std::istringstream in( opt[ "schedule" ] );
   std::string w;
   while( std::getline( in , w , ',' ) )
    x.schedule.push_back( std::stod( w ) );
   if( x.schedule.empty() )
    x.schedule = { o.weight };
   }
  return( anytime_astar( d , x , ok ) );
  }
 if( ( algo == "optimistic" ) || ( algo == "ios" ) ) {
  OptimisticOptions x;
  static_cast< SearchOptions & >( x ) = o;
  x.optimism = std::stod( opt[ "optimism" ] );
  if( algo == "ios" )
   return( ImprovedOptimisticSearch< D >( d , x ).solve() );
  return( optimistic_search( d , x , ok ) );
  }
 if( algo == "focal" ) {
  FocalOptions x;
  static_cast< SearchOptions & >( x ) = o;
  const std::string & f = opt[ "focal" ];
  x.order = f == "h" ? FocalOrder::heuristic
            : f == "opt" ? FocalOrder::optimistic
            : f == "ees" ? FocalOrder::explicit_estimation
            : FocalOrder::distance;
  return( FocalSearch< D >( d , x ).solve() );
  }
 if( algo == "dps" )
  return( PotentialSearch< D >( d , o ).solve() );
 if( algo == "lpa" ) {
  if constexpr( PredecessorDomain< D > )
   return( LPAStar< D >( d , o ).solve() );
  throw( std::invalid_argument( "no predecessors for this domain" ) );
  }
 if( algo == "spa" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( SPAStar< D , ADS::ThreadPool >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( SPAStar< D , ADS::OpenMPExecutor >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "mqa" ) {
  const std::string & r = opt[ "runtime" ];
  if( r == "thread" )
   return( MQAStar< D , ADS::ThreadPool >( d , o ).solve() );
  #ifdef PBEFS_HAS_OPENMP
  if( r == "openmp" )
   return( MQAStar< D , ADS::OpenMPExecutor >( d , o ).solve() );
  #endif
  throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
  }
 if( algo == "pbnf" ) {
  if constexpr( NBlockDomain< D > ) {
   PBNFOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.min_expansions = unsigned( std::stoul( opt[ "minexp" ] ) );
   po.safe = opt[ "safe" ] != "0";
   std::istringstream in( opt[ "schedule" ] );
   std::string sw;
   while( std::getline( in , sw , ',' ) )
    po.schedule.push_back( std::stod( sw ) );
   return( on_runtime( [ & ]< class E >() {
    return( SafePBNF< D , E >( d , po ).solve() );
    } ) );
   }
  throw( std::invalid_argument( "no nblocks for this domain" ) );
  }
 if( ( algo == "psdd" ) || ( algo == "bfpsdd" ) || ( algo == "idpsdd" ) ||
     ( algo == "dynpsdd" ) ) {
  if constexpr( NBlockDomain< D > ) {
   PSDDOptions po( o );
   po.min_expansions = unsigned( std::stoul( opt[ "minexp" ] ) );
   po.multiplier = std::stod( opt[ "mult" ] );
   return( on_runtime( [ & ]< class E >() {
    if( algo == "psdd" )
     return( PSDD< D , E >( d , po ).solve( bound_of< cost_t< D > >() ) );
    if( algo == "bfpsdd" )
     return( BFPSDD< D , E >( d , po ).solve() );
    if( algo == "idpsdd" ) {
     PSDDOptions io( o );
     io.growth = std::stod( opt[ "growth" ] );
     return( IDPSDD< D , E >( d , io ).solve() );
     }
    ParallelOptions x = o;                    // the weight is the A*'s
    x.weight = 1;
    return( DynamicBoundedPSDD< D , E >( d , o.weight , x ).solve() );
    } ) );
   }
  throw( std::invalid_argument( "no nblocks for this domain" ) );
  }
 if( ( algo == "farm" ) || ( algo == "kbfs" ) ) {
  FarmOptions fo( o );
  fo.synchronous = algo == "kbfs";
  return( on_runtime< true >( [ & ]< class E >() {
   return( FarmAStar< D , E >( d , fo ).solve() );
   } ) );
  }
 if( algo == "pase" ) {
  if constexpr( PairwiseDomain< D > ) {
   PASEOptions po;
   static_cast< ParallelOptions & >( po ) = o;
   po.eps = std::stod( opt[ "eps" ] );
   po.candidates = std::stoul( opt[ "scan" ] );
   const std::string & r = opt[ "runtime" ];
   if( r == "thread" )
    return( PASE< D , ADS::ThreadPool >( d , po ).solve() );
   #ifdef PBEFS_HAS_OPENMP
   if( r == "openmp" )
    return( PASE< D , ADS::OpenMPExecutor >( d , po ).solve() );
   #endif
   throw( std::invalid_argument( "unknown or unavailable runtime " + r ) );
   }
  throw( std::invalid_argument( "no pairwise heuristic for this domain" ) );
  }
 if( ( algo == "epase" ) || ( algo == "pwa" ) || ( algo == "mplp" ) )
  return( solve_edges( d , o ) );
 #endif
 if( algo == "pra" ) {
  opt[ "sync" ] = "1";
  opt[ "recv" ] = "locked";
  }
 else
  if( ( algo != "hda" ) && ( algo != "tds" ) )
   throw( std::invalid_argument( "unknown search " + algo ) );

 Solution< state_t< D > , cost_t< D > > sampled;
 if( solve_sampled( d , o , sampled ) )
  return( sampled );
 const std::string & dist = opt[ "dist" ];
 if constexpr( FeatureDomain< D > )
  if( dist == "Z" )
   return( run_estimated( d , ZobristDistribution< D >( d ) , o ) );
 if constexpr( AbstractFeatureDomain< D > )
  if( dist == "AZ" )
   return( run_estimated( d , AbstractZobristDistribution< D >( d ) , o ) );
 if constexpr( AbstractStateDomain< D > )
  if( dist == "A" )
   return( run_estimated( d , AbstractStateDistribution< D >( d ) , o ) );
 if( dist == "H" )
  return( run_estimated( d , StateHashDistribution< D >( d ) , o ) );
 if constexpr( LevelDomain< D > )
  if( dist == "HWD" ) {
   auto [ planes , split ] = HyperplaneDistribution< D >::thickness(
                              d , o.threads );
   const std::string & t = opt[ "plane" ];
   if( t.rfind( "1/" , 0 ) == 0 ) {
    planes = 1;
    split = unsigned( std::stoul( t.substr( 2 ) ) );
    }
   else
    if( t != "auto" ) {
     planes = unsigned( std::stoul( t ) );
     split = 1;
     }
   return( run_hda( d , HyperplaneDistribution< D >( d , 1 , planes ,
                                                       split ) , o ) );
   }
 if constexpr( PerfectHashDomain< D > )
  if( dist == "P" )
   return( run_estimated( d , PerfectHashDistribution< D >( d ) , o ) );
 #ifndef PBEFS_BENCH_MPI
 if( dist == "R" )
  return( run_estimated( d , RandomDistribution< D >( d ) , o ) );
 #endif
 throw( std::invalid_argument( "distribution " + dist +
                               " not available for this domain" ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// solves d, slowed down by --slow iterations per expansion if not 0, or
// with the edges of --edges

template< class D >
Solution< state_t< D > , cost_t< D > > solve_maybe_slow( const D & d )
{
 const unsigned long work = std::stoul( opt[ "slow" ] );
 if( ! opt[ "edges" ].empty() ) {
  EdgeWork ew;
  unsigned long long seed = 1;
  if( std::sscanf( opt[ "edges" ].c_str() , "%lu:%lu:%lf:%lf:%llu" ,
                   & ew.cheap , & ew.expensive , & ew.expensive_share ,
                   & ew.missing_share , & seed ) < 3 )
   throw( std::invalid_argument( "--edges wants cheap:expensive:share"
                                 "[:missing[:seed]]" ) );
  ew.seed = seed;
  return( solve( SlowEdgeDomain< D >( d , ew ) ) );
  }
 if( work == 0 )
  return( solve( d ) );
 return( solve( SlowDomain< D >( d , work ) ) );
 }

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_AStarBenchSolve */

/*--------------------------------------------------------------------------*/
/*----------------------- End File AStarBenchSolve.h -----------------------*/
/*--------------------------------------------------------------------------*/
