/*--------------------------------------------------------------------------*/
/*--------------------------- File planning.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for classical planning:
 * the reader of the SAS+ tasks, their heuristics and the searches on the
 * tasks of IPC of data/planning.
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
/*--------------------------------- CASES ----------------------------------*/
/*--------------------------------------------------------------------------*/
// a task in the format output.sas, from its lines

static SASTask sas( std::initializer_list< const char * > lines )
{
 std::ostringstream os;
 for( const char * l : lines )
  os << l << '\n';
 std::istringstream is( os.str() );
 return( SASTask::read( is ) );
 }

// the header of a task with one variable of three values and metric m

static std::vector< const char * > sas_head( const char * m )
{
 return( std::vector< const char * >{
          "begin_version" , "3" , "end_version" , "begin_metric" , m ,
          "end_metric" , "1" , "begin_variable" , "x" , "-1" , "3" , "a" ,
          "b" , "c" , "end_variable" , "0" , "begin_state" , "0" ,
          "end_state" , "begin_goal" , "1" , "0 2" , "end_goal" } );
 }

// the task of sas_head( m ) and the lines of rest

static SASTask sas_with( const char * m ,
                         std::initializer_list< const char * > rest )
{
 std::ostringstream os;
 for( const char * l : sas_head( m ) )
  os << l << '\n';
 for( const char * l : rest )
  os << l << '\n';
 std::istringstream is( os.str() );
 return( SASTask::read( is ) );
 }

// true if reading the task of sas_head( m ) and rest throws

static bool sas_refused( const char * m ,
                         std::initializer_list< const char * > rest )
{
 try {
  ( void ) sas_with( m , rest );
  }
 catch( std::invalid_argument & ) {
  return( true );
  }
 return( false );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the reader, the heuristics and the sequential searches on the planning
// tasks, against their optimal costs

void test_planning( const std::string & data )
{
 SearchOptions gen;                       // the goal tested when generated
 gen.goal_at_generation = true;

 // the reader: the costs, the effects with no value before, the refusals
 const auto t3 = sas_with( "1" , {
  "3" ,
  "begin_operator" , "ab" , "0" , "1" , "0 0 0 1" , "5" , "end_operator" ,
  "begin_operator" , "bc" , "0" , "1" , "0 0 1 2" , "1" , "end_operator" ,
  "begin_operator" , "any-c" , "0" , "1" , "0 0 -1 2" , "7" ,
  "end_operator" , "0" } );
 check( ( t3.num_vars() == 1 ) && ( t3.ops.size() == 3 ) && t3.use_costs &&
        ( t3.ops[ 2 ].effects[ 0 ].pre == -1 ) , "a task with costs is read" );
 for( auto h : { PlanningHeuristic::blind , PlanningHeuristic::hmax ,
                 PlanningHeuristic::lmcut , PlanningHeuristic::pdb ,
                 PlanningHeuristic::mas } ) {
  PlanningOptions o;
  o.heuristic = h;
  const PlanningD d( t3 , o );
  check_plan( d , astar( d ) , 6 , "the costs 5 + 1 against 7" );
  }
 const auto u3 = sas_with( "0" , {
  "1" , "begin_operator" , "ab" , "0" , "1" , "0 0 0 1" , "5" ,
  "end_operator" , "0" } );
 check( ! u3.use_costs && ( u3.ops[ 0 ].cost == 1 ) ,
        "with no metric every operator costs 1" );
 {
  PlanningOptions o;
  o.heuristic = PlanningHeuristic::hmax;
  const PlanningD d( u3 , o );
  const auto sol = astar( d );
  check( ( d.heuristic_or_dead( d.initial_state() ) == dead_end ) &&
         ( d.heuristic( d.initial_state() ) == 0 ) &&
         ( sol.status == Status::unsolvable ) && ( sol.stats.expanded == 1 ) ,
         "a goal that no operator reaches is a dead end" );
  }
 check( sas_refused( "1" , { "1" , "begin_operator" , "c" , "0" , "1" ,
                             "1 0 1 0 0 2" , "1" , "end_operator" , "0" } ) ,
        "a conditional effect is refused" );
 check( sas_refused( "1" , { "0" , "1" } ) , "an axiom is refused" );
 check( sas_refused( "1" , { "1" , "begin_operator" , "n" , "0" , "1" ,
                             "0 0 0 2" , "-1" , "end_operator" , "0" } ) ,
        "a negative cost is refused" );
 check( sas_refused( "1" , { "1" , "begin_operator" , "n" , "0" , "1" ,
                             "0 0 0 3" , "1" , "end_operator" , "0" } ) ,
        "a value out of range is refused" );
 bool axiom_var = false , version = false;
 try {
  ( void ) sas( { "begin_version" , "3" , "end_version" , "begin_metric" ,
                  "0" , "end_metric" , "1" , "begin_variable" , "x" , "0" ,
                  "2" , "a" , "b" , "end_variable" } );
  }
 catch( std::invalid_argument & ) { axiom_var = true; }
 try {
  ( void ) sas( { "begin_version" , "2" , "end_version" } );
  }
 catch( std::invalid_argument & ) { version = true; }
 check( axiom_var && version ,
        "a derived variable and the version 2 are refused" );
 {
  // 65 binary variables need 2 words; the goal is reached by one operator
  SASTask big;
  for( int v = 0 ; v < 65 ; ++v ) {
   big.var_names.push_back( "v" + std::to_string( v ) );
   big.range.push_back( 2 );
   big.values.push_back( { "0" , "1" } );
   big.init.push_back( 0 );
   }
  big.goal = { { 64 , 1 } };
  big.ops.push_back( { "set" , {} , { { 64 , 0 , 1 } } , 1 } );
  bool too_big = false;
  try {
   ( void ) Planning< 1 >( big );
   }
  catch( std::invalid_argument & ) { too_big = true; }
  const Planning< 2 > d2( big );
  check( too_big && ( PlanningBase::words_for( big ) == 2 ) &&
         ( astar( d2 ).cost == 1 ) && ( d2.value( d2.initial_state() , 64 )
                                       == 0 ) ,
         "a state that does not fit in the words is refused" );
  }

 // the optimal costs with every admissible heuristic and the plans of the
 // inadmissible ones, but where A* expands more than a few thousand nodes
 // (thousands of evaluations of LM-cut or h^FF take minutes in the build
 // of the coverage, which is not optimized); the expansions of all the
 // pairs are in the logs of the comparisons with the reference codes
 const std::map< PlanningHeuristic , std::set< std::string > > heavy = {
  { PlanningHeuristic::blind , { "gripper-03.sas" , "gripper-04.sas" ,
    "blocks-7-0.sas" , "blocks-7-1.sas" , "blocks-7-2.sas" ,
    "driverlog-02.sas" , "driverlog-03.sas" , "logistics-4-0.sas" ,
    "logistics-4-1.sas" , "logistics-4-2.sas" , "logistics-5-0.sas" ,
    "logistics-5-1.sas" , "logistics-6-0.sas" , "logistics-6-1.sas" ,
    "logistics-6-2.sas" } } ,
  { PlanningHeuristic::hmax , { "gripper-03.sas" , "gripper-04.sas" ,
    "blocks-7-0.sas" , "blocks-7-1.sas" , "blocks-7-2.sas" ,
    "driverlog-02.sas" , "logistics-4-0.sas" , "logistics-4-1.sas" ,
    "logistics-4-2.sas" , "logistics-5-0.sas" , "logistics-5-1.sas" ,
    "logistics-6-0.sas" , "logistics-6-2.sas" } } ,
  { PlanningHeuristic::lmcut , { "gripper-03.sas" , "gripper-04.sas" ,
    "driverlog-02.sas" } } ,
  { PlanningHeuristic::pdb , { "gripper-04.sas" } } ,
  { PlanningHeuristic::mas , { "gripper-04.sas" } } };
 const std::set< std::string > light_inadmissible = {
  "blocks-4-0.sas" , "blocks-5-1.sas" , "blocks-6-0.sas" , "gripper-01.sas" ,
  "driverlog-01.sas" , "driverlog-03.sas" , "logistics-5-2.sas" ,
  "logistics-6-1.sas" };
 const auto tasks = planning_tasks( data );
 check( tasks.size() == 28 , "the 28 planning tasks are listed" );
 for( const auto & [ name , opt ] : tasks ) {
  for( const auto & [ h , skip ] : heavy ) {
   if( skip.count( name ) )
    continue;
   const auto d = planning( data , name , h );
   check_plan( d , astar( d , OpenKind::automatic , gen ) , opt ,
               "A* " + name );
   }
  if( ! light_inadmissible.count( name ) )
   continue;
  // the inadmissible heuristics: a valid plan, not cheaper than the optimum
  for( auto h : { PlanningHeuristic::hadd , PlanningHeuristic::hff ,
                  PlanningHeuristic::goal_count } ) {
   const auto d = planning( data , name , h );
   const auto sol = astar( d );
   check( sol.solved() && ( sol.cost >= opt ) &&
          ( simulate( d.task() , d.plan( sol.path ) ) == sol.cost ) ,
          "A* with an inadmissible heuristic, " + name );
   }
  }

 // the expansions of A* with LM-cut, the same in two runs, and with the
 // OPEN list of ADS those of its order of the ties, with the goal tested at
 // the selection (the default) and at the generation
 for( const auto & [ name , exp , exp_gen ] :
      std::vector< std::tuple< std::string , std::size_t , std::size_t > >{
       { "blocks-6-2.sas" , 284 , 283 } , { "logistics-4-0.sas" , 77 , 76 } ,
       { "driverlog-03.sas" , 23 , 22 } , { "gripper-01.sas" , 83 , 82 } } ) {
  const auto d = planning( data , name , PlanningHeuristic::lmcut );
  const std::size_t e = astar( d ).stats.expanded;
  const std::size_t eg = astar( d , OpenKind::automatic , gen ).stats.expanded;
  bool ok = ( e == astar( d ).stats.expanded ) &&
            ( eg == astar( d , OpenKind::automatic , gen ).stats.expanded ) &&
            ( eg <= e );
  #ifdef PBEFS_STRUCTURES_ADS
   ok = ok && ( e == exp ) && ( eg == exp_gen );
  #else
   ( void ) exp;
   ( void ) exp_gen;
  #endif
  check( ok , "the expansions of A* with LM-cut on " + name );
  }

 // the goal test at generation: with the blind heuristic A* expands
 // exactly the states closer than the optimum (counted by a breadth-first
 // search, the costs being 1), and without it at least as many
 for( const char * name : { "blocks-6-2.sas" , "driverlog-01.sas" ,
                            "logistics-5-2.sas" , "gripper-02.sas" } ) {
  const auto d = planning( data , name , PlanningHeuristic::blind );
  std::unordered_map< std::uint64_t , std::vector< PlanningD::state_type > >
   seen;
  const auto known = [ & ]( const PlanningD::state_type & s ) {
   auto & b = seen[ d.hash( s ) ];
   if( std::find( b.begin() , b.end() , s ) != b.end() )
    return( true );
   b.push_back( s );
   return( false );
   };
  std::vector< PlanningD::state_type > layer{ d.initial_state() } , next;
  ( void ) known( layer[ 0 ] );
  std::size_t closer = 0;
  for( bool goal = false ; ! goal ; ) {
   for( const auto & s : layer )
    goal = goal || d.is_goal( s );
   if( goal )
    break;
   closer += layer.size();
   next.clear();
   for( const auto & s : layer )
    d.for_each_successor( s , 0 , [ & ]( const PlanningD::state_type & t ,
                                         int , int ) {
     if( ! known( t ) )
      next.push_back( t );
     } );
   layer.swap( next );
   }
  const auto early = astar( d , OpenKind::automatic , gen );
  const auto sel = astar( d );
  check( early.solved() && sel.solved() && ( early.cost == sel.cost ) &&
         ( early.stats.expanded == closer ) &&
         ( early.stats.goals_generated >= 1 ) &&
         ( sel.stats.goals_generated == 0 ) &&
         ( sel.stats.expanded >= closer ) ,
         std::string( "the goal test at generation, " ) + name + " (" +
         std::to_string( early.stats.expanded ) + " and " +
         std::to_string( sel.stats.expanded ) + " expansions, " +
         std::to_string( closer ) + " states closer than the optimum)" );
  check_plan( d , early , sel.cost , std::string( "A* " ) + name +
                                     " with the goal at generation" );
  }

 // along an optimal plan: the admissible heuristics below the cost to go,
 // LM-cut above h^max, the PDB consistent; the features, the nblocks and
 // the inverse operators on every edge out of the plan
 for( const char * name : { "blocks-6-2.sas" , "logistics-5-0.sas" ,
                            "driverlog-03.sas" , "gripper-01.sas" } ) {
  const auto hm = planning( data , name , PlanningHeuristic::hmax );
  const auto lm = planning( data , name , PlanningHeuristic::lmcut );
  const auto pd = planning( data , name , PlanningHeuristic::pdb );
  const auto ms = planning( data , name , PlanningHeuristic::mas );
  const auto sol = astar( lm );
  bool adm = sol.solved() , dom = true , cons = true , feat = true ,
       nb = true , inv = true , abs = true , mcons = true , dtg = true;
  std::set< std::pair< std::size_t , std::size_t > > arcs;
  lm.for_each_feature_transition( [ & ]( std::size_t i , std::size_t j ) {
   dtg = dtg && ( lm.feature_group( i ) == lm.feature_group( j ) );
   arcs.emplace( i , j );
   } );
  for( std::size_t i = 0 ; i < sol.path.size() ; ++i ) {
   const auto & s = sol.path[ i ];
   const int togo = int( sol.path.size() - 1 - i );  // unit costs
   const int a = hm.heuristic( s ) , b = lm.heuristic( s ) ,
             c = pd.heuristic( s ) , m = ms.heuristic( s );
   adm = adm && ( a <= togo ) && ( b <= togo ) && ( c <= togo ) &&
         ( m <= togo );
   dom = dom && ( a <= b );
   std::set< std::size_t > fs;
   lm.for_each_feature( s , [ & ]( std::size_t f ) {
    fs.insert( f );
    abs = abs && ( lm.abstract_feature( f ) < lm.num_abstract_features() );
    } );
   const std::size_t bs = lm.nblock( s );
   std::set< std::size_t > nbs;
   lm.for_each_nblock_successor( bs , [ & ]( std::size_t x ) {
    nbs.insert( x );
    bool back = false;
    lm.for_each_nblock_successor( x , [ & ]( std::size_t y ) {
     back = back || ( y == bs );
     } );
    nb = nb && back;
    } );
   std::vector< std::pair< PlanningD::state_type , std::uint8_t > > all;
   lm.for_each_successor( s , 0 , no_op , [ & ]( const auto & t , int ,
                                                 int , std::uint8_t op ) {
    all.emplace_back( t , op );
    } );
   for( const auto & [ t , op ] : all ) {
    cons = cons && ( c <= 1 + pd.heuristic( t ) );
    mcons = mcons && ( m <= 1 + ms.heuristic( t ) );
    for( int v = 0 ; v < int( lm.task().num_vars() ) ; ++v )
     if( lm.value( s , v ) != lm.value( t , v ) ) {
      const auto o = std::size_t( lm.strips().offset[ std::size_t( v ) ] );
      dtg = dtg && arcs.count( { o + std::size_t( lm.value( s , v ) ) ,
                                 o + std::size_t( lm.value( t , v ) ) } );
      }
    std::set< std::size_t > ft , diff;
    lm.for_each_feature( t , [ & ]( std::size_t f ) { ft.insert( f ); } );
    std::set_symmetric_difference( fs.begin() , fs.end() , ft.begin() ,
                                   ft.end() , std::inserter( diff ,
                                                             diff.end() ) );
    std::set< std::size_t > inc;
    lm.for_each_feature_change( s , t , [ & ]( std::size_t f ) {
     inc.insert( f );
     } );
    feat = feat && ( inc == diff );
    const std::size_t bt = lm.nblock( t );
    nb = nb && ( ( bt == bs ) || nbs.count( bt ) );
    // from t, generated by op, the successors skipped are only s
    if( op != no_op ) {
     std::size_t kept = 0 , total = 0;
     bool back = false;
     lm.for_each_successor( t , 0 , op , [ & ]( const auto & , int , int ,
                                                std::uint8_t ) { ++kept; } );
     lm.for_each_successor( t , 0 , no_op , [ & ]( const auto & x , int ,
                                                   int , std::uint8_t ) {
      ++total;
      back = back || ( x == s );
      } );
     inv = inv && back && ( ( kept == total ) || ( kept + 1 == total ) );
     }
    }
   }
  const std::string tag = std::string( ", " ) + name;
  check( adm , "h^max, LM-cut, the PDB and merge-and-shrink are "
               "admissible" + tag );
  check( mcons , "merge-and-shrink is consistent" + tag );
  check( dtg , "every move is an arc of a domain transition graph" + tag );
  check( dom , "LM-cut is not below h^max" + tag );
  check( cons , "the PDB is consistent" + tag );
  check( feat , "the features changed by each move" + tag );
  check( abs , "the abstract features are in range" + tag );
  check( nb , "the nblock graph holds every move and is undirected" + tag );
  check( inv , "a coded operator skips its parent only" + tag );
  }

 // IDA* and the delayed duplicate detection on the disk
 // (IDA* on the tasks with few transpositions within the optimal cost)
 for( const char * name : { "blocks-5-1.sas" , "gripper-02.sas" ,
                            "logistics-5-2.sas" , "driverlog-01.sas" } ) {
  const auto d = planning( data , name , PlanningHeuristic::lmcut );
  const int opt = astar( d ).cost;
  if( std::string( name ) != "gripper-02.sas" )
   check_plan( d , idastar( d ) , opt , std::string( "IDA* " ) + name );
  check_plan( d , ddd_astar( d ) , opt , std::string( "DDD A* " ) + name );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel searches on the planning tasks, 1 to 8 workers; SafePBNF,
// SPA* and MQA* only if All

template< class Executor , bool All >
void test_planning_parallel( const std::string & data )
{
 const auto tasks = planning_tasks( data );
 const unsigned cores = std::max( 1u , std::thread::hardware_concurrency() );
 // merge-and-shrink, whose evaluation is a few lookups, and at most 512
 // nblocks, so that the build of the coverage stays within its time
 for( const char * name : { "blocks-6-2.sas" , "logistics-4-0.sas" ,
                            "gripper-01.sas" , "driverlog-03.sas" } ) {
  const int opt = std::find_if( tasks.begin() , tasks.end() ,
                                [ & ]( const auto & x ) {
                                 return( x.first == name );
                                 } )->second;
  const auto d = planning( data , name , PlanningHeuristic::mas , 512 );
  PlanningOptions po1;
  po1.heuristic = PlanningHeuristic::mas;
  po1.max_nblocks = 512;
  const Planning< 1 > d1( SASTask::read_file( data + "/planning/" + name ) ,
                          po1 );
  for( unsigned p : { 1u , 2u , 3u , 8u } ) {
   if( p > cores )
    break;
   ParallelOptions o;
   o.threads = p;
   const std::string tag = std::string( " p " ) + std::to_string( p ) +
                           ", " + name;
   const ZobristDistribution< PlanningD > z( d );
   check_plan( d , HDAStar< PlanningD , decltype( z ) , Executor >( d , z , o )
                   .solve() , opt , "HDA*[Z]" + tag );
   const AbstractZobristDistribution< PlanningD > az( d );
   check_plan( d , HDAStar< PlanningD , decltype( az ) , Executor >( d , az ,
                    o ).solve() , opt , "HDA*[AZ]" + tag );
   const AbstractStateDistribution< PlanningD > a( d );
   check_plan( d , HDAStar< PlanningD , decltype( a ) , Executor >( d , a , o )
                   .solve() , opt , "HDA*[A]" + tag );
   const StateHashDistribution< PlanningD > h( d );
   check_plan( d , HDAStar< PlanningD , decltype( h ) , Executor >( d , h , o )
                   .solve() , opt , "HDA*[H]" + tag );
   check_plan( d , FarmAStar< PlanningD , Executor >( d , o ).solve() , opt ,
               "FarmAStar" + tag );
   const PerfectHashDistribution< Planning< 1 > > ph( d1 );
   const auto r1 = HDAStar< Planning< 1 > , decltype( ph ) , Executor >(
                    d1 , ph , o ).solve();
   check( r1.solved() && ( r1.cost == opt ) &&
          ( simulate( d1.task() , d1.plan( r1.path ) ) == opt ) ,
          "HDA*[P] on one word" + tag );
   if constexpr( All ) {
    PBNFOptions po;
    po.threads = p;
    check_plan( d , SafePBNF< PlanningD , Executor >( d , po ).solve() , opt ,
                "SafePBNF" + tag );
    check_plan( d , SPAStar< PlanningD , Executor >( d , o ).solve() , opt ,
                "SPA*" + tag );
    check_plan( d , MQAStar< PlanningD , Executor >( d , o ).solve() , opt ,
                "MQA*" + tag );
    }
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_planning_parallel< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_planning_parallel< ADS::OpenMPExecutor >
   ( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_planning_parallel< ADS::FastFlowExecutor , false >
   ( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File planning.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
