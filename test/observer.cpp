/*--------------------------------------------------------------------------*/
/*--------------------------- File observer.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the observers of the
 * searches: a search with a TraceRecorder makes the same search as with
 * no observer, and the events recorded agree with its counters.
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

#include "PBeFS/core/TraceRecorder.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- HELPERS ---------------------------------*/
/*--------------------------------------------------------------------------*/
// the id of a state for a TraceRecorder: its hash in the Domain

template< class D >
struct HashId {
 const D * d;
 std::uint64_t operator()( const state_t< D > & s ) const {
  return( d->hash( s ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the id of a node of the knapsack: the items fixed and the profit

struct KnapsackId {
 std::uint64_t operator()( const Knapsack::node_type & n ) const {
  return( ( std::uint64_t( n.k ) << 40 ) ^ std::uint64_t( n.profit ) );
  }
 };

template< class D >
using Recorder = TraceRecorder< HashId< D > >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the counters of two searches are the same, but the time

static bool same_stats( const SearchStats & a , const SearchStats & b )
{
 return( ( a.expanded == b.expanded ) && ( a.reexpanded == b.reexpanded ) &&
         ( a.generated == b.generated ) && ( a.duplicates == b.duplicates ) &&
         ( a.stored == b.stored ) && ( a.sent == b.sent ) &&
         ( a.messages == b.messages ) && ( a.expanded_by == b.expanded_by ) &&
         ( a.stored_by == b.stored_by ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// checks that the events of r agree with the counters st: the expansions,
// also per worker if per_worker, the generations and the duplicates, a
// solution if the search solved, and as many pops as pushes at most

template< class R >
static void check_events( const R & r , const SearchStats & st ,
                          bool solved , bool per_worker ,
                          const std::string & name )
{
 const unsigned all = r.workers();
 bool ok = ( r.count( TraceEvent::expand , all ) == st.expanded ) &&
           ( r.count( TraceEvent::generate , all ) == st.generated ) &&
           ( r.count( TraceEvent::duplicate , all ) == st.duplicates ) &&
           ( r.count( TraceEvent::open_pop , all ) <=
             r.count( TraceEvent::open_push , all ) ) &&
           ( ( r.count( TraceEvent::solution , all ) > 0 ) == solved ) &&
           ( r.dropped() == 0 );
 if( per_worker ) {
  ok = ok && ( st.expanded_by.size() == all );
  for( unsigned i = 0 ; ok && ( i < all ) ; ++i )
   ok = ( r.count( TraceEvent::expand , i ) == st.expanded_by[ i ] );
  }
 std::ostringstream os;
 os << name << ": the events agree with the counters (" << r.size()
    << " events, expanded " << st.expanded << ")";
 check( ok , os.str() );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the trace written by r, as lines

template< class R >
static std::vector< std::string > trace_lines( const R & r )
{
 std::ostringstream os;
 r.write( os );
 std::vector< std::string > v;
 std::istringstream is( os.str() );
 for( std::string line ; std::getline( is , line ) ; )
  v.push_back( line );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// A* on d with no observer and with a recorder: the same search, events
// that agree with its counters, the steps 0 , 1 , ... , the same trace at
// every run, and the cost of the last solution that of the path

template< class D >
static void check_astar( const D & d , const std::string & name )
{
 using A = AStar< D , default_open_t< D > , Recorder< D > >;
 const auto plain = AStar< D >( d ).solve();
 Recorder< D > rec( HashId< D >{ & d } , 0 );
 A a( d );
 a.set_observer( rec );
 const auto obs = a.solve();
 check( ( plain.status == obs.status ) && ( plain.cost == obs.cost ) &&
        ( plain.path == obs.path ) && same_stats( plain.stats , obs.stats ) ,
        "A* " + name + ": the same search with the recorder" );
 check_events( rec , obs.stats , obs.solved() , true , "A* " + name );
 const auto & ev = rec.events( 0 );
 bool steps = true;
 for( std::size_t i = 0 ; i < ev.size() ; ++i )
  steps = steps && ( ev[ i ].step == i ) && ( ev[ i ].worker == 0 );
 const auto last = std::find_if( ev.rbegin() , ev.rend() ,
                                 []( const TraceRecord & x ) {
                                  return( x.type == TraceEvent::solution );
                                  } );
 check( steps && ( ! obs.solved() ||
                   ( ( last != ev.rend() ) &&
                     ( last->g == double( obs.cost ) ) ) ) ,
        "A* " + name + ": the steps and the solution of the trace" );
 const auto first = trace_lines( rec );
 const auto again = a.solve();
 check( ( first == trace_lines( rec ) ) && ( first.size() == ev.size() + 1 )
        && same_stats( again.stats , obs.stats ) ,
        "A* " + name + ": the same trace at a second run" );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- CASES ----------------------------------*/
/*--------------------------------------------------------------------------*/
// the sequential searches with a TraceRecorder: A*, weighted A*, IDA* and
// the B&B, the cap, the mask, the sampling and the trace written

void test_observer( const std::string & data )
{
 const auto inst = korf( data );
 const EightPuzzle e8( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 } ) );
 check_astar( e8 , "8-puzzle" );
 const FifteenPuzzle f15( inst[ easy_korf[ 0 ] - 1 ].board );
 check_astar( f15 , "Korf #" + std::to_string( easy_korf[ 0 ] ) );
 const auto m = GridMap::random( 120 , 80 , 0.25 , 2 );
 check_astar( Grid< 8 >( m , 0 , 79 , 119 , 79 ) , "grid 8" );
 check_astar( romania( Arad ) , "Romania" );
 GridMap wall( 30 , 30 );
 for( int y = 0 ; y < 30 ; ++y )
  wall.set_cost( 15 , y , 0 );
 check_astar( Grid< 4 >( wall , 0 , 0 , 29 , 29 ) , "unsolvable grid" );

 {  // weighted A*, with its reopenings dropped
  SearchOptions o;
  o.weight = 2;
  o.drop_duplicates = true;
  const auto plain = AStar< FifteenPuzzle >( f15 , o ).solve();
  Recorder< FifteenPuzzle > rec( HashId< FifteenPuzzle >{ & f15 } , 0 );
  AStar< FifteenPuzzle , default_open_t< FifteenPuzzle > ,
         Recorder< FifteenPuzzle > > a( f15 , o );
  a.set_observer( rec );
  const auto obs = a.solve();
  check( ( plain.cost == obs.cost ) && same_stats( plain.stats , obs.stats ) ,
         "weighted A* w 2: the same search with the recorder" );
  check_events( rec , obs.stats , true , true , "weighted A* w 2" );
  bool keys = true;
  for( const auto & x : rec.events( 0 ) )
   if( x.type == TraceEvent::expand )
    keys = keys && ( x.f == x.g + 2 * x.h );
  check( keys , "weighted A* w 2: f = g + 2 h in the trace" );
  }

 {  // IDA*: the thresholds increase, the last one the cost
  const auto plain = IDAStar< FifteenPuzzle >( f15 ).solve();
  Recorder< FifteenPuzzle > rec( HashId< FifteenPuzzle >{ & f15 } , 0 );
  IDAStar< FifteenPuzzle , Recorder< FifteenPuzzle > > a( f15 );
  a.set_observer( rec );
  const auto obs = a.solve();
  check( ( plain.cost == obs.cost ) && ( plain.path == obs.path ) &&
         same_stats( plain.stats , obs.stats ) ,
         "IDA*: the same search with the recorder" );
  check_events( rec , obs.stats , true , true , "IDA*" );
  std::vector< double > th;
  for( const auto & x : rec.events( 0 ) )
   if( x.type == TraceEvent::threshold )
    th.push_back( x.g );
  check( ! th.empty() && std::is_sorted( th.begin() , th.end() ) &&
         ( std::adjacent_find( th.begin() , th.end() ) == th.end() ) &&
         ( th.back() == double( obs.cost ) ) &&
         ( th.front() == double( f15.heuristic( f15.initial_state() ) ) ) ,
         "IDA*: " + std::to_string( th.size() ) + " increasing thresholds "
         "from h of the root to the cost" );
  }

 {  // the cap, the mask and the sampling
  Recorder< EightPuzzle > all( HashId< EightPuzzle >{ & e8 } , 0 );
  Recorder< EightPuzzle > cap( HashId< EightPuzzle >{ & e8 } , 100 );
  Recorder< EightPuzzle > smp( HashId< EightPuzzle >{ & e8 } , 0 ,
                               trace_bit( TraceEvent::expand ) , 7 );
  using A = AStar< EightPuzzle , default_open_t< EightPuzzle > ,
                   Recorder< EightPuzzle > >;
  A a( e8 );
  a.set_observer( all );
  const auto s0 = a.solve();
  a.set_observer( cap );
  const auto s1 = a.solve();
  a.set_observer( smp );
  const auto s2 = a.solve();
  const auto & ea = all.events( 0 );
  check( same_stats( s0.stats , s1.stats ) &&
         same_stats( s0.stats , s2.stats ) && ( cap.size() == 100 ) &&
         ( cap.dropped() == all.size() - 100 ) &&
         std::equal( cap.events( 0 ).begin() , cap.events( 0 ).end() ,
                     ea.begin() , []( const TraceRecord & x ,
                                      const TraceRecord & y ) {
                      return( ( x.step == y.step ) && ( x.id == y.id ) &&
                              ( x.type == y.type ) );
                      } ) ,
         "the cap of 100 events keeps the first 100" );
  check( ( smp.size() == ( s2.stats.expanded + 6 ) / 7 ) &&
         ( smp.count( TraceEvent::expand , 1 ) == smp.size() ) &&
         ( trace_lines( cap ).front().find( "\"dropped\":" +
           std::to_string( cap.dropped() ) ) != std::string::npos ) ,
         "the mask of the expansions, one in 7, and the header" );
  }

 {  // the observer is required by a recorder
  AStar< EightPuzzle , default_open_t< EightPuzzle > ,
         Recorder< EightPuzzle > > a( e8 );
  bool thrown = false;
  try {
   static_cast< void >( a.solve() );
   }
  catch( const std::logic_error & ) {
   thrown = true;
   }
  check( thrown , "a search without its observer throws" );
  }

 // the B&B of the knapsack: evaluations, branchings, incumbents
 for( int cls : { 0 , 1 , 2 } ) {
  const auto kp = Knapsack::random( cls , 40 , 1000 , 7 );
  BBOptions o;
  const auto plain = BranchAndBoundSearch< Knapsack >( kp , o ).solve();
  TraceRecorder< KnapsackId > rec( KnapsackId{} , 0 );
  BranchAndBoundSearch< Knapsack , TraceRecorder< KnapsackId > > b( kp , o );
  b.set_observer( rec );
  const auto obs = b.solve();
  const auto & a = plain.stats;
  const auto & s = obs.stats;
  check( ( plain.value == obs.value ) && ( obs.value == kp.optimum() ) &&
         ( a.evaluated == s.evaluated ) && ( a.branched == s.branched ) &&
         ( a.generated == s.generated ) &&
         ( a.pruned_generated == s.pruned_generated ) &&
         ( a.pruned_extracted == s.pruned_extracted ) &&
         ( a.pruned_evaluated == s.pruned_evaluated ) ,
         "B&B knapsack class " + std::to_string( cls ) +
         ": the same search with the recorder" );
  double best = -1;
  bool better = true;
  for( const auto & x : rec.events( 0 ) )
   if( x.type == TraceEvent::solution ) {
    better = better && ( x.g > best );
    best = x.g;
    }
  const std::uint64_t pruned = rec.count( TraceEvent::prune , 1 );
  check( ( rec.count( TraceEvent::evaluate , 1 ) == s.evaluated ) &&
         ( rec.count( TraceEvent::branch , 1 ) == s.branched ) &&
         ( rec.count( TraceEvent::solution , 1 ) == s.incumbents ) &&
         ( pruned >= s.pruned_generated + s.pruned_evaluated ) &&
         ( pruned <= s.pruned_generated + s.pruned_evaluated +
                     s.pruned_extracted ) &&
         better && ( best == obs.value ) ,
         "B&B knapsack class " + std::to_string( cls ) + ": " +
         std::to_string( rec.size() ) + " events agree with the counters" );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the parallel searches with a TraceRecorder, 1 to 4 workers on the given
// runtime (SPA* and SafePBNF only if all, as in their own cases): the same
// cost as with no observer, events that agree with the counters (per
// worker, but for the farm, whose master expands), a message for each node
// sent, and the buffers of the workers apart

template< class Executor , bool all >
void test_observer_parallel( const std::string & data )
{
 const auto inst = korf( data );
 const FifteenPuzzle f15( inst[ easy_korf[ 0 ] - 1 ].board );
 const int c15 = inst[ easy_korf[ 0 ] - 1 ].cost;
 const EightPuzzle e8( walk8( { 1 , 2 , 5 , 4 , 3 , 6 , 7 , 8 , 5 , 2 , 1 ,
                                0 , 3 , 4 , 7 , 6 } ) );
 const int c8 = astar( e8 ).cost;
 const std::string rt = executor_name< Executor >();
 using F = FifteenPuzzle;
 using E = EightPuzzle;
 using R15 = Recorder< F >;
 using R8 = Recorder< E >;
 for( unsigned p : { 1u , 2u , 4u } ) {
  const std::string tag = ", " + rt + " p " + std::to_string( p );
  ParallelOptions o;
  o.threads = p;

  if constexpr( all ) {
   R15 rec( HashId< F >{ & f15 } , 0 );
   SPAStar< F , Executor , default_open_t< F > , R15 > s( f15 , o );
   s.set_observer( rec );
   const auto sol = s.solve();
   check_solution( f15 , sol , c15 , "SPA* with the recorder" + tag );
   check_events( rec , sol.stats , true , true , "SPA*" + tag );
   }

  {
   const ZobristDistribution< F > z( f15 );
   R15 rec( HashId< F >{ & f15 } , 0 );
   HDAStar< F , ZobristDistribution< F > , Executor , default_open_t< F > ,
            R15 > h( f15 , z , o );
   h.set_observer( rec );
   const auto sol = h.solve();
   check_solution( f15 , sol , c15 , "HDA* with the recorder" + tag );
   check_events( rec , sol.stats , true , true , "HDA*" + tag );
   std::uint64_t sent = 0;
   bool apart = true;
   for( unsigned i = 0 ; i < rec.workers() ; ++i )
    for( const auto & x : rec.events( i ) ) {
     apart = apart && ( x.worker == i );
     if( x.type == TraceEvent::message ) {
      sent += x.size;
      apart = apart && ( x.other != i ) && ( x.other < p );
      }
     }
   check( apart && ( sent == sol.stats.sent ) ,
          "HDA*: a message per node sent, each worker in its buffer" + tag );
   }

  if constexpr( all ) {
   PBNFOptions po;
   po.threads = p;
   R15 rec( HashId< F >{ & f15 } , 0 );
   SafePBNF< F , Executor , default_open_t< F > , R15 > s( f15 , po );
   s.set_observer( rec );
   const auto sol = s.solve();
   check_solution( f15 , sol , c15 , "SafePBNF with the recorder" + tag );
   check_events( rec , sol.stats , true , true , "SafePBNF" + tag );
   }

  {
   PASEOptions po;
   po.threads = p;
   R8 rec( HashId< E >{ & e8 } , 0 );
   PASE< E , Executor , R8 > s( e8 , po );
   s.set_observer( rec );
   const auto sol = s.solve();
   check_solution( e8 , sol , c8 , "PA*SE with the recorder" + tag );
   check_events( rec , sol.stats , true , true , "PA*SE" + tag );
   }

  {
   R15 rec( HashId< F >{ & f15 } , 0 );
   FarmAStar< F , Executor , default_open_t< F > , R15 > s( f15 , o );
   s.set_observer( rec );
   const auto sol = s.solve();
   check_solution( f15 , sol , c15 , "the farm with the recorder" + tag );
   check_events( rec , sol.stats , true , false , "the farm" + tag );
   std::uint64_t back = 0;
   bool from = rec.workers() == p + 1;
   for( unsigned i = 0 ; i < rec.workers() ; ++i )
    for( const auto & x : rec.events( i ) )
     if( x.type == TraceEvent::message ) {
      back += x.size;
      from = from && ( i > 0 ) && ( x.other == 0 );
      }
   check( from && ( back == sol.stats.generated ) &&
          ( rec.count( TraceEvent::expand , 0 ) == sol.stats.expanded ) ,
          "the farm: the master expands, the workers send back" + tag );
   }

  for( int dist : { 0 , 1 } ) {
   const auto kp = Knapsack::random( 1 , 40 , 1000 , 7 );
   ParallelBBOptions bo;
   bo.threads = p;
   TraceRecorder< KnapsackId > rec( KnapsackId{} , 0 );
   const auto run = [ & ]( auto && b ) {
    b.set_observer( rec );
    return( b.solve() );
    };
   const auto sol = dist ?
    run( DistributedBranchAndBound< Knapsack , Executor ,
                                    TraceRecorder< KnapsackId > >( kp , bo ) )
    : run( SharedBranchAndBound< Knapsack , Executor ,
                                 TraceRecorder< KnapsackId > >( kp , bo ) );
   const unsigned nw = rec.workers();
   bool per = ( nw == p ) && ( sol.stats.evaluated_by.size() == p );
   for( unsigned i = 0 ; per && ( i < p ) ; ++i )
    per = ( rec.count( TraceEvent::evaluate , i ) ==
            sol.stats.evaluated_by[ i ] );
   check( ( sol.value == kp.optimum() ) && per &&
          ( rec.count( TraceEvent::branch , nw ) == sol.stats.branched ) &&
          ( rec.count( TraceEvent::solution , nw ) ==
            sol.stats.incumbents ) &&
          ( ! dist || ( rec.count( TraceEvent::message , nw ) ==
                        sol.stats.sent ) ) ,
          std::string( dist ? "distributed" : "shared" ) +
          " B&B with the recorder" + tag );
   }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template void test_observer_parallel< ADS::ThreadPool >(
                                                       const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_observer_parallel< ADS::OpenMPExecutor >(
                                                       const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_observer_parallel< ADS::FastFlowExecutor , false >(
                                                       const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- End File observer.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
