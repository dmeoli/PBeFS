/*--------------------------------------------------------------------------*/
/*------------------------------ File bb.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the Branch-and-Bound
 * (B&B) of the knapsack: its frontiers, the parallel schemes, the values
 * and the cutoffs.
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
// the knapsack as a minimization of the opposite profit, which exercises
// the other sense

struct MinKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::minimize;

 const Knapsack & kp;

 node_type root( void ) const { return( kp.root() ); }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e = kp.evaluate( n );
  if( e.infeasible )
   return( e );
  e.bound = - e.bound;
  e.feasible = - *e.feasible;
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  kp.branch( n , [ & ]( node_type && c , double e ) {
   f( std::move( c ) , - e );
   } );
  }

 void serialize( const node_type & n , std::vector< char > & out ) const {
  kp.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const {
  return( kp.deserialize( in ) );
  }

 };  // end( struct( MinKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack whose workers take the incumbent by cutoff( z ), and find
// infeasible a node whose bound cannot exceed z, as a relaxation stopped
// by its cutoff does; a z above the optimum is recorded in wrong

struct CutoffKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::maximize;

 struct Worker {
  const Knapsack & kp;
  std::atomic< bool > & wrong;
  double opt;
  double z = - std::numeric_limits< double >::infinity();

  void cutoff( double v ) {
   if( v > opt )
    wrong = true;
   z = v;
   }

  Evaluation evaluate( node_type & n ) const {
   Evaluation e = kp.evaluate( n );
   if( ( ! e.infeasible ) && ( e.bound <= z ) ) {
    e = Evaluation();
    e.infeasible = true;
    }
   return( e );
   }

  template< class F >
  void branch( const node_type & n , F && f ) const {
   kp.branch( n , std::forward< F >( f ) );
   }
  };

 const Knapsack & kp;
 std::atomic< bool > & wrong;
 double opt;                      // the optimum of kp

 node_type root( void ) const { return( kp.root() ); }

 Worker make_worker( unsigned ) const {
  return( Worker{ kp , wrong , opt } );
  }

 };  // end( struct( CutoffKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack whose worker i counts its evaluations in count[ i ]

struct CountKnapsack {

 using node_type = Knapsack::node_type;

 static constexpr Sense sense = Sense::maximize;

 struct Worker {
  const Knapsack & kp;
  std::atomic< std::uint64_t > & count;

  Evaluation evaluate( node_type & n ) const {
   ++count;
   return( kp.evaluate( n ) );
   }

  template< class F >
  void branch( const node_type & n , F && f ) const {
   kp.branch( n , std::forward< F >( f ) );
   }
  };

 const Knapsack & kp;
 std::vector< std::atomic< std::uint64_t > > & count;

 node_type root( void ) const { return( kp.root() ); }

 Worker make_worker( unsigned i ) const {
  return( Worker{ kp , count.at( i ) } );
  }

 };  // end( struct( CountKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the instances of the knapsack of the tests: 2 of each class of Pisinger

static std::vector< Knapsack > knapsacks( void )
{
 std::vector< Knapsack > v;
 for( int cls : { 0 , 1 , 2 } )
  for( std::uint64_t seed : { 1 , 2 } )
   v.push_back( Knapsack::random( cls , cls == 2 ? 30 : 50 , 1000 ,
                                  seed ) );
 return( v );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the leaves of the B&B on the knapsack kp, solved by solve( q , o , seeds )
// for the problem q (Knapsack or MinKnapsack), the options o and the nodes
// to start from (none for the root): the leaves of a search carried out to
// the end bound the optimum and cover the tree, so that a search started
// from them solves the knapsack of a smaller capacity, in both senses; a
// search started from the open leaves of one stopped by a node limit, with
// its incumbent, ends it; and a search started from no node ends at once.
// With explicit children, the leaves of the smaller knapsack, the
// infeasible ones among them, cover the tree of the larger one, which a
// search started from them solves again. With MPI, each process has its
// own leaves, which may be none.

using Leaves = std::vector< BBLeaf< Knapsack::node_type > >;

template< class O , class Solve >
static void check_frontier_of( const Knapsack & kp , O o , Solve && solve ,
                               const std::string & name , bool local )
{
 constexpr double inf = std::numeric_limits< double >::infinity();
 const double opt = double( kp.optimum() );
 Knapsack kq = kp;
 kq.cap = kp.cap * 3 / 4;
 const double optq = double( kq.optimum() );
 o.keep_frontier = true;

 const auto s = solve( kp , o , std::optional< Leaves >() );
 bool ok = ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           ( local || ! s.frontier.empty() );
 Leaves seeds;
 for( const auto & l : s.frontier ) {
  ok = ok && ( l.kind != BBLeafKind::open ) && ( l.bound <= opt + 1e-9 );
  seeds.push_back( l );
  seeds.back().bound = inf;            // unknown, in the maximization
  }
 check( ok , name + " leaves of the whole tree" );
 const auto r = solve( kq , o , std::optional< Leaves >( seeds ) );
 check( ( r.status == BBStatus::optimal ) && ( r.value == optq ) &&
        r.node && kq.check( *r.node , r.value ) ,
        name + " from the leaves, a smaller capacity" );
 if( kp.explicit_children ) {
  seeds = r.frontier;
  ok = true;
  for( auto & x : seeds ) {
   ok = ok && ( ( x.kind == BBLeafKind::infeasible ) == ( x.bound == - inf ) );
   x.bound = inf;
   }
  const auto b = solve( kp , o , std::optional< Leaves >( seeds ) );
  check( ok && ( b.status == BBStatus::optimal ) && ( b.value == opt ) &&
         b.node && kp.check( *b.node , b.value ) ,
         name + " from the leaves, the capacity back" );
  }

 const auto m = solve( MinKnapsack{ kp } , o , std::optional< Leaves >() );
 ok = m.status == BBStatus::optimal;
 seeds.clear();
 for( const auto & l : m.frontier ) {
  ok = ok && ( l.bound >= - opt - 1e-9 );
  seeds.push_back( l );
  seeds.back().bound = - inf;          // unknown, in the minimization
  }
 const auto mr = solve( MinKnapsack{ kq } , o ,
                        std::optional< Leaves >( seeds ) );
 check( ok && ( mr.status == BBStatus::optimal ) && ( mr.value == - optq ) &&
        mr.node && kq.check( *mr.node , - mr.value ) ,
        name + " from the leaves, a smaller capacity, minimization" );

 O ol = o;
 ol.node_limit = 5;
 const auto l = solve( kp , ol , std::optional< Leaves >() );
 seeds.clear();
 for( const auto & x : l.frontier )
  if( x.kind == BBLeafKind::open )
   seeds.push_back( x );
 O oc = o;
 if( l.feasible() )
  oc.incumbent = l.value;
 const auto c = solve( kp , oc , std::optional< Leaves >( seeds ) );
 check( ( l.status == BBStatus::node_limit ) &&
        ( local || ! seeds.empty() ) && ( c.status == BBStatus::optimal ) &&
        ( c.value == opt ) && ( c.gap <= 1e-9 ) ,
        name + " from the open leaves of a node limit" );

 const auto e = solve( kp , o , std::optional< Leaves >( Leaves() ) );
 check( ( e.status == BBStatus::infeasible ) && ( e.stats.evaluated == 0 ) &&
        e.frontier.empty() , name + " from no node" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the same, with the children of the knapsack getting the bound of their
// parent and with explicit children, which prunes some of them when they
// are generated and gives infeasible leaves

template< class O , class Solve >
static void check_frontier( const Knapsack & kp , const O & o ,
                            Solve && solve , const std::string & name ,
                            bool local = false )
{
 Knapsack kb = kp;
 kb.explicit_children = true;
 check_frontier_of( kp , o , solve , name , local );
 check_frontier_of( kb , o , solve , name + " explicit children" , local );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a tree of a root and two children without children, both of bound 0:
// the first with the feasible value 1, the second with none, so that
// neither is pruned and both are leaves, and the optimum is 1

struct Fork {

 using node_type = int;

 static constexpr Sense sense = Sense::minimize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & n ) const {
  Evaluation e;
  if( n == 1 )
   e.feasible = 1;
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  if( n == 0 ) {
   f( 1 , 0.0 );
   f( 2 , 0.0 );
   }
  }

 };  // end( struct( Fork ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a root that is infeasible

struct NoRoot {

 using node_type = int;

 static constexpr Sense sense = Sense::minimize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & ) const {
  Evaluation e;
  e.infeasible = true;
  return( e );
  }

 template< class F >
 void branch( const node_type & , F && ) const {}

 };  // end( struct( NoRoot ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a root solved without branching: its bound is its feasible value, 5

struct Leaf {

 using node_type = int;

 static constexpr Sense sense = Sense::maximize;

 node_type root( void ) const { return( 0 ); }

 Evaluation evaluate( node_type & ) const {
  Evaluation e;
  e.bound = 5;
  e.feasible = 5;
  return( e );
  }

 template< class F >
 void branch( const node_type & , F && ) const {}

 };  // end( struct( Leaf ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the limit cases of a B&B solved by run( q , o ) for the problem q and the
// options o: a root that is infeasible, a root solved without branching,
// the tree of Fork, and, on the knapsack kp, an initial incumbent better
// than the optimum (which nothing improves), an absolute tolerance, a time
// limit that is reached at once, a node limit of 1, and plunging with the
// dive going always and never

template< class O , class Run >
static void check_bb_edges( const Knapsack & kp , const O & base , Run && run ,
                            const std::string & name )
{
 const double opt = double( kp.optimum() );

 const auto n = run( NoRoot{} , base );
 check( ( n.status == BBStatus::infeasible ) && ! n.node && ! n.feasible() ,
        name + " infeasible root" );

 const auto l = run( Leaf{} , base );
 check( ( l.status == BBStatus::optimal ) && ( l.value == 5 ) &&
        ( l.bound == 5 ) && ( l.gap == 0 ) && ( l.stats.evaluated == 1 ) ,
        name + " root solved without branching" );

 const auto f = run( Fork{} , base );
 check( ( f.status == BBStatus::optimal ) && ( f.value == 1 ) &&
        ( f.bound == 1 ) , name + " fork" );

 O o = base;
 o.incumbent = opt + 1;
 const auto z = run( kp , o );
 check( ( z.status == BBStatus::optimal ) && ( z.value == opt + 1 ) &&
        ! z.node && ( z.bound >= opt + 1 - 1e-9 ) ,
        name + " initial incumbent better than the optimum" );

 o = base;
 o.abs_gap = 0.05 * opt;
 const auto a = run( kp , o );
 check( ( a.status == BBStatus::optimal ) &&
        ( a.value >= opt - o.abs_gap - 1e-9 ) && ( a.bound >= opt - 1e-9 ) &&
        ( a.gap <= o.abs_gap + 1e-9 ) , name + " absolute tolerance" );

 o = base;
 o.time_limit = 1e-12;
 const auto t = run( kp , o );
 check( ( t.status == BBStatus::time_limit ) && ( t.bound >= opt - 1e-9 ) &&
        ( ! t.feasible() || ( t.value <= opt ) ) , name + " time limit" );

 o = base;
 o.node_limit = 1;
 const auto one = run( kp , o );
 check( ( one.status == BBStatus::node_limit ) &&
        ( one.stats.evaluated >= 1 ) && ( one.bound >= opt - 1e-9 ) ,
        name + " node limit of 1" );

 for( double g : { 0.0 , 1.0 } ) {
  o = base;
  o.selection = Selection::plunging;
  o.plunge_gap = g;
  const auto p = run( kp , o );
  check( ( p.status == BBStatus::optimal ) && ( p.value == opt ) ,
         name + " plunging with plunge_gap " + std::to_string( g ) );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the sequential B&B on the knapsack against the dynamic program, with
// every selection, in both senses, with a relative tolerance, with an
// initial incumbent and with a node limit, its leaves, and its limit cases

void test_bb( const std::string & )
{
 const auto solve = []( const auto & q , const BBOptions & o ,
                        std::optional< Leaves > seeds ) {
  BranchAndBoundSearch< std::decay_t< decltype( q ) > > b( q , o );
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   BBOptions o;
   o.selection = sel;
   const auto s = branch_and_bound( kp , o );
   const std::string t = tag + " selection " +
                         std::to_string( int( sel ) );
   check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
          s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
          "B&B optimum" + t );
   const auto m = branch_and_bound( MinKnapsack{ kp } , o );
   check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
          m.node && kp.check( *m.node , - m.value ) ,
          "B&B minimization" + t );
   check_frontier( kp , o , solve , "B&B" + t );
   }

  BBOptions o;
  o.rel_gap = 0.05;
  const auto r = branch_and_bound( kp , o );
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) &&
         ( r.gap <= 0.05 * r.value + 1e-9 ) ,
         "B&B within a relative gap of 5%" + tag );

  o = BBOptions();
  o.incumbent = opt;
  double seen = 0;
  o.on_incumbent = [ & ]( double v ) { seen = v; };
  const auto z = branch_and_bound( kp , o );
  check( ( z.status == BBStatus::optimal ) && ( z.value == opt ) &&
         ! z.node && ( seen == 0 ) , "B&B with the optimum given" + tag );

  // the optimum read from outside prunes as the one given, but it is not a
  // value of the search; a worse one leaves the search its own optimum
  o = BBOptions();
  o.shared_incumbent = [ & ] { return( opt ); };
  const auto x = branch_and_bound( kp , o );
  check( ( x.status == BBStatus::infeasible ) && ! x.node &&
         ( x.bound <= opt + 1e-9 ) &&
         ( x.stats.evaluated == z.stats.evaluated ) ,
         "B&B with the optimum shared" + tag );
  o.shared_incumbent = [ & ] { return( opt / 2 ); };
  const auto y = branch_and_bound( kp , o );
  check( ( y.status == BBStatus::optimal ) && ( y.value == opt ) &&
         y.node && kp.check( *y.node , y.value ) ,
         "B&B with a worse value shared" + tag );

  o = BBOptions();
  o.node_limit = 3;
  const auto l = branch_and_bound( kp , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) &&
         ( l.stats.evaluated == 3 ) , "B&B node limit" + tag );
  }

 BBOptions o;
 o.keep_frontier = true;
 const auto f = branch_and_bound( Fork{} , o );
 std::set< std::pair< int , int > > leaves;
 for( const auto & x : f.frontier )
  leaves.insert( { x.node , int( x.kind ) } );
 check( ( f.status == BBStatus::optimal ) && ( f.value == 1 ) &&
        ( f.bound == 1 ) &&
        ( leaves == std::set< std::pair< int , int > >{
                     { 1 , int( BBLeafKind::feasible ) } ,
                     { 2 , int( BBLeafKind::pruned ) } } ) ,
        "B&B leaves without children" );

 check_bb_edges( knapsacks().front() , BBOptions() ,
                 []( const auto & q , const BBOptions & x ) {
                  return( branch_and_bound( q , x ) );
                  } , "B&B" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a parallel B&B of the problem q: K = 0 the shared queue, K = 1 the
// distributed queues, K = 2 the farm, K = 3 the work stealing, K = 4 the
// MultiQueue

template< int K , class Executor , class Q >
static auto bb_parallel( const Q & q , const ParallelBBOptions & o ,
                         std::optional< std::vector< BBLeaf<
                          typename Q::node_type > > > seeds = std::nullopt )
{
 const auto solve = [ & ]( const auto & b ) {
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 if constexpr( K == 0 )
  return( solve( SharedBranchAndBound< Q , Executor >( q , o ) ) );
 else
  if constexpr( K == 1 )
   return( solve( DistributedBranchAndBound< Q , Executor >( q , o ) ) );
  else
   if constexpr( K == 3 )
    return( solve( StealingBranchAndBound< Q , Executor >( q , o ) ) );
   else
    if constexpr( K == 4 )
     return( solve( MultiQueueBranchAndBound< Q , Executor >( q , o ) ) );
    else
     return( solve( FarmBranchAndBound< Q , Executor >( q , o ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a parallel B&B on the knapsack against the dynamic program, with 1 to 8
// workers and every selection, in both senses (the minimization has no
// hash( n ), and hence random owners), with a relative tolerance and with
// a node limit, its leaves, and its limit cases with 8 workers (more than
// the nodes of the small trees) and with the parameters at 0, which mean 1

template< int K , class Executor >
void test_bb_parallel( const std::string & )
{
 const std::string name = std::string( K == 0 ? "shared" :
                                       K == 1 ? "distributed" :
                                       K == 2 ? "farm" :
                                       K == 3 ? "stealing" : "multiqueue" ) +
                          " B&B on " + executor_name< Executor >();
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( unsigned p : { 1u , 2u , 3u , 8u } )
   for( Selection sel : { Selection::best_first , Selection::depth_first ,
                          Selection::breadth_first ,
                          Selection::best_estimate , Selection::plunging } ) {
    ParallelBBOptions o;
    o.threads = p;
    o.selection = sel;
    const auto s = bb_parallel< K , Executor >( kp , o );
    check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
           name + tag + " p " + std::to_string( p ) + " selection " +
           std::to_string( int( sel ) ) );
    }

  ParallelBBOptions o;
  o.threads = 4;
  const auto m = bb_parallel< K , Executor >( MinKnapsack{ kp } , o );
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         name + " minimization" + tag );

  o.rel_gap = 0.05;
  const auto r = bb_parallel< K , Executor >( kp , o );
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) ,
         name + " within a relative gap of 5%" + tag );

  o = ParallelBBOptions();
  o.threads = 4;
  o.node_limit = 3;
  const auto l = bb_parallel< K , Executor >( kp , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) ,
         name + " node limit" + tag );

  for( Selection sel : { Selection::best_first , Selection::plunging } ) {
   o = ParallelBBOptions();
   o.threads = 3;
   o.selection = sel;
   check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                                std::optional< Leaves > seeds ) {
    return( bb_parallel< K , Executor >( q , x , std::move( seeds ) ) );
    } , name + tag + " selection " + std::to_string( int( sel ) ) );
   }
  }

 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( bb_parallel< K , Executor >( q , x ) );
  };
 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.batch = 0;
 o.flush_period = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.batch = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " batches of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the deterministic B&B on the knapsack against the dynamic program, with
// every selection: with 4 nodes per round its outcome (the value, the node,
// the bound and the counters) is the same with 1 to 8 workers and at every
// repetition, also when a node limit stops it, and with one node per worker
// it is the optimum in both senses; its leaves and its limit cases

template< class Executor >
void test_bb_det( const std::string & )
{
 const std::string name = "deterministic B&B on " +
                          executor_name< Executor >();
 using Out = std::tuple< int , double , double , std::uint64_t ,
                         std::uint64_t , std::uint64_t , std::uint64_t ,
                         std::uint64_t , std::vector< char > >;
 const auto outcome = []( const auto & s ) {
  return( Out{ int( s.status ) , s.value , s.bound , s.stats.evaluated ,
               s.stats.generated , s.stats.pruned_generated ,
               s.stats.pruned_extracted , s.stats.pruned_evaluated ,
               s.node ? s.node->x : std::vector< char >() } );
  };
 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( DeterministicBranchAndBound< std::decay_t< decltype( q ) > ,
                                       Executor >( q , x ).solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( std::size_t limit : { std::size_t( 0 ) , std::size_t( 7 ) } ) {
    ParallelBBOptions o;
    o.selection = sel;
    o.round = 4;
    o.node_limit = limit;
    std::optional< Out > ref;
    bool ok = true , same = true;
    for( unsigned p : { 1u , 2u , 3u , 8u } )
     for( int rep = 0 ; rep < 2 ; ++rep ) {
      o.threads = p;
      const auto s = run( kp , o );
      if( limit )
       ok = ok && ( s.status == BBStatus::node_limit ) &&
            ( s.stats.evaluated <= limit ) && ( s.bound >= opt - 1e-9 );
      else
       ok = ok && ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
            s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 );
      const Out out = outcome( s );
      if( ! ref )
       ref = out;
      else
       same = same && ( out == *ref );
      }
    const std::string t = name + tag + " selection " +
                          std::to_string( int( sel ) ) +
                          ( limit ? " node limit" : "" );
    check( ok , t );
    check( same , t + ", the same with 1 to 8 workers" );
    }

  ParallelBBOptions o;
  o.threads = 3;
  const auto m = run( MinKnapsack{ kp } , o );
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         name + " minimization, one node per worker" + tag );
  check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                               std::optional< Leaves > seeds ) {
   DeterministicBranchAndBound< std::decay_t< decltype( q ) > , Executor >
    b( q , x );
   return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
   } , name + tag );
  }

 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.round = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.round = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " rounds of 1" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the SelfSplit B&B on the knapsack against the dynamic program, with 1 to
// 8 workers, every selection, a split of 1 and 4 nodes per worker, and the
// sampling with and without its dive: its outcome is the same at every
// repetition, also when a node limit stops it, each worker evaluates its
// nodes by its own evaluator, as many as the counters say, and the values
// given to on_incumbent improve each time and end at the optimum; in both
// senses, its leaves and its limit cases; with the incumbent shared, the
// optimum with every selection and number of workers, and its limit cases;
// with 1 worker, the nodes of the sequential search

template< class Executor >
void test_bb_split( const std::string & )
{
 const std::string name = "SelfSplit B&B on " + executor_name< Executor >();
 using Out = std::tuple< int , double , double , std::uint64_t ,
                         std::uint64_t , std::vector< std::uint64_t > ,
                         std::vector< char > >;
 const auto outcome = []( const auto & s ) {
  return( Out{ int( s.status ) , s.value , s.bound , s.stats.evaluated ,
               s.stats.generated , s.stats.evaluated_by ,
               s.node ? s.node->x : std::vector< char >() } );
  };
 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( SelfSplitBranchAndBound< std::decay_t< decltype( q ) > ,
                                   Executor >( q , x ).solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( std::size_t split : { std::size_t( 1 ) , std::size_t( 4 ) } )
    for( std::size_t limit : { std::size_t( 0 ) , std::size_t( 7 ) } )
     for( bool dive : { true , false } ) {
      bool ok = true , same = true , own = true;
      for( unsigned p : { 1u , 2u , 3u , 8u } ) {
       ParallelBBOptions o;
       o.selection = sel;
       o.split = split;
       o.split_dive = dive;
       o.node_limit = limit;
       o.threads = p;
       std::optional< Out > ref;
       for( int rep = 0 ; rep < 2 ; ++rep ) {
        std::vector< std::atomic< std::uint64_t > > count( p );
        const auto s = run( CountKnapsack{ kp , count } , o );
        if( limit )
         ok = ok && ( s.status == BBStatus::node_limit ) &&
              ( s.stats.evaluated <= limit ) && ( s.bound >= opt - 1e-9 );
        else
         ok = ok && ( s.status == BBStatus::optimal ) &&
              ( s.value == opt ) && s.node && kp.check( *s.node , s.value ) &&
              ( s.gap <= 1e-9 );
        own = own && ( s.stats.evaluated_by.size() == p );
        for( unsigned j = 0 ; own && ( j < p ) ; ++j )
         own = s.stats.evaluated_by[ j ] == count[ j ].load();
        const Out out = outcome( s );
        if( ! ref )
         ref = out;
        else
         same = same && ( out == *ref );
        }
       }
      const std::string t = name + tag + " selection " +
                            std::to_string( int( sel ) ) + " split " +
                            std::to_string( split ) +
                            ( limit ? " node limit" : "" ) +
                            ( dive ? "" : " without the dive" );
      check( ok , t );
      check( same , t + ", the same at every repetition" );
      check( own , t + ", the evaluations of each worker" );
      }

  ParallelBBOptions o;
  o.threads = 3;
  std::vector< double > seen;
  std::mutex m;
  o.on_incumbent = [ & ]( double v ) {
   std::lock_guard< std::mutex > guard( m );
   seen.push_back( v );
   };
  const auto s = run( kp , o );
  bool up = ( ! seen.empty() ) && ( seen.back() == opt ) &&
            ( s.value == opt );
  for( std::size_t j = 1 ; j < seen.size() ; ++j )
   up = up && ( seen[ j ] > seen[ j - 1 ] );
  check( up , name + " incumbents reported improving" + tag );

  bool alone = true;
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   ParallelBBOptions a;
   a.selection = sel;
   a.threads = 1;
   BBOptions b;
   b.selection = sel;
   const auto x = run( kp , a );
   const auto y = BranchAndBoundSearch< Knapsack >( kp , b ).solve();
   alone = alone && ( x.value == y.value ) &&
           ( x.stats.evaluated == y.stats.evaluated ) &&
           ( x.stats.generated == y.stats.generated );
   }
  check( alone , name + " with 1 worker evaluates the nodes of the "
                        "sequential search" + tag );

  o = ParallelBBOptions();
  o.threads = 3;
  const auto mn = run( MinKnapsack{ kp } , o );
  check( ( mn.status == BBStatus::optimal ) && ( mn.value == - opt ) &&
         mn.node && kp.check( *mn.node , - mn.value ) ,
         name + " minimization" + tag );
  check_frontier( kp , o , []( const auto & q , const ParallelBBOptions & x ,
                               std::optional< Leaves > seeds ) {
   SelfSplitBranchAndBound< std::decay_t< decltype( q ) > , Executor >
    b( q , x );
   return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
   } , name + tag );

  bool shared = true;
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } )
   for( unsigned p : { 1u , 2u , 3u , 8u } ) {
    o = ParallelBBOptions();
    o.selection = sel;
    o.threads = p;
    o.split_shared = true;
    const auto s = run( kp , o );
    const auto m = run( MinKnapsack{ kp } , o );
    shared = shared && ( s.status == BBStatus::optimal ) &&
             ( s.value == opt ) && s.node && kp.check( *s.node , s.value ) &&
             ( s.gap <= 1e-9 ) && ( m.status == BBStatus::optimal ) &&
             ( m.value == - opt ) && m.node &&
             kp.check( *m.node , - m.value );
    }
  check( shared , name + " incumbent shared" + tag );
  }

 ParallelBBOptions o;
 o.threads = 8;
 check_bb_edges( knapsacks().front() , o , run , name + " 8 workers" );
 o.threads = 0;
 o.split = 0;
 check_bb_edges( knapsacks().front() , o , run , name + " parameters 0" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.split = 1;
 check_bb_edges( knapsacks().front() , o , run , name + " split of 1" );
 o.split_shared = true;
 check_bb_edges( knapsacks().front() , o , run , name + " incumbent shared" );
 o = ParallelBBOptions();
 o.threads = 3;
 o.split_dive = false;
 check_bb_edges( knapsacks().front() , o , run , name + " without the dive" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the B&B chosen by its scheme: each name gives back its scheme, an unknown
// name is refused, and every scheme solves the knapsack in both senses, from
// the root and from the leaves of a search, with 1 and 3 workers

void test_bb_scheme( const std::string & )
{
 bool names = true;
 for( BBScheme s : { BBScheme::sequential , BBScheme::shared ,
                     BBScheme::distributed , BBScheme::deterministic ,
                     BBScheme::selfsplit , BBScheme::stealing ,
                     BBScheme::multiqueue } )
  names = names && ( bb_scheme( bb_scheme_name( s ) ) == s );
 bool refused = false;
 try {
  ( void ) bb_scheme( "nonesuch" );
  }
 catch( const std::invalid_argument & ) {
  refused = true;
  }
 check( names && refused , "B&B schemes by name" );

 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( std::string_view name : { "sequential" , "shared" , "distributed" ,
                                 "deterministic" , "selfsplit" ,
                                 "stealing" , "multiqueue" } )
   for( unsigned p : { 1u , 3u } ) {
    ParallelBBOptions o;
    o.threads = p;
    const BBScheme s = bb_scheme( name );
    const auto r = parallel_branch_and_bound( kp , s , o );
    const auto m = parallel_branch_and_bound( MinKnapsack{ kp } , s , o );
    const std::string t = "B&B scheme " + std::string( name ) + tag +
                          " p " + std::to_string( p );
    check( ( r.status == BBStatus::optimal ) && ( r.value == opt ) &&
           r.node && kp.check( *r.node , r.value ) &&
           ( m.status == BBStatus::optimal ) && ( m.value == - opt ) , t );
    check_frontier( kp , o , [ & ]( const auto & q ,
                                    const ParallelBBOptions & x ,
                                    std::optional< Leaves > seeds ) {
     return( parallel_branch_and_bound( q , s , x , std::move( seeds ) ) );
     } , t );
    }
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a rational number of an exact arithmetic, as a B&B value: num / den with
// den > 0 and no common factor, or an infinity, den = 0 and num = +-1,
// whose products fit a long long for the values of the tests; it is
// trivially copyable and 16 bytes wide, so its incumbent goes under a lock
// where a std::atomic of that width is not lock-free (e.g., on x86-64)

struct Rational {

 long long num = 0;                             // the numerator
 long long den = 1;                             // the denominator

 Rational( void ) = default;

 Rational( long long n , long long d = 1 ) : num( n ) , den( d ) {
  if( den < 0 ) {
   num = - num;
   den = - den;
   }
  const long long g = std::gcd( num , den );
  if( g > 1 ) {
   num /= g;
   den /= g;
   }
  }

 static Rational infinite( int sign ) {
  Rational r;
  r.num = sign;
  r.den = 0;
  return( r );
  }

 bool operator==( const Rational & o ) const {
  return( ( num == o.num ) && ( den == o.den ) );
  }

 auto operator<=>( const Rational & o ) const {
  if( ! den || ! o.den )                       // an infinity on either side
   return( ( den ? 0 : num ) <=> ( o.den ? 0 : o.num ) );
  return( ( num * o.den ) <=> ( o.num * den ) );
  }

 Rational operator-( void ) const {
  Rational r = *this;
  r.num = - r.num;
  return( r );
  }

 Rational operator-( const Rational & o ) const {
  return( Rational( num * o.den - o.num * den , den * o.den ) );
  }

 };  // end( struct( Rational ) )

// the traits of Rational: an exact tolerance, the relative one taken as a
// fraction with denominator 10^6, and the test of a dive in long double

template<>
struct PBeFS_di_unipi_it::BBValueTraits< Rational > {
 static Rational infinity( void ) { return( Rational::infinite( 1 ) ); }
 static Rational none( void ) { return( infinity() ); }
 static bool finite( const Rational & v ) { return( v.den != 0 ); }
 static Rational tolerance( const Rational & z , const Rational & a ,
                            double r ) {
  const Rational t( ( z.num < 0 ? - z.num : z.num ) *
                    (long long)( r * 1e6 ) , z.den * 1000000 );
  return( std::max( a , t ) );
  }
 static bool within( const Rational & b , const Rational & l ,
                     const Rational & z , double g ) {
  if( ! finite( l ) )
   return( false );
  using L = long double;
  const auto x = []( const Rational & v ) { return( L( v.num ) / v.den ); };
  return( x( b ) <= x( l ) + L( g ) * ( x( z ) - x( l ) ) );
  }
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack with exact values: the fractional bound of Dantzig as a
// Rational, and the profits as integers

struct ExactKnapsack {

 using node_type = Knapsack::node_type;
 using value_type = Rational;

 static constexpr Sense sense = Sense::maximize;

 const Knapsack & kp;

 node_type root( void ) const { return( kp.root() ); }

 BasicEvaluation< Rational > evaluate( node_type & n ) const {
  BasicEvaluation< Rational > e;
  if( n.weight > kp.cap ) {
   e.infeasible = true;
   return( e );
   }
  long r = kp.cap - n.weight , b = n.profit;
  std::size_t j = n.k;
  for( ; ( j < kp.p.size() ) && ( kp.w[ j ] <= r ) ; ++j ) {
   r -= kp.w[ j ];
   b += kp.p[ j ];
   }
  e.bound = j < kp.p.size() ? Rational( b * kp.w[ j ] + kp.p[ j ] * r ,
                                        kp.w[ j ] )
                            : Rational( b );
  e.feasible = Rational( kp.greedy( n ) );
  return( e );
  }

 template< class F >
 void branch( const node_type & n , F && f ) const {
  kp.branch( n , [ & ]( node_type && c , double ) {
   f( std::move( c ) , BBValueTraits< Rational >::infinity() );
   } );
  }

 std::uint64_t hash( const node_type & n ) const { return( kp.hash( n ) ); }

 void serialize( const node_type & n , std::vector< char > & out ) const {
  kp.serialize( n , out );
  }

 node_type deserialize( const char * & in ) const {
  return( kp.deserialize( in ) );
  }

 };  // end( struct( ExactKnapsack ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the knapsack kp with values of type long

static BasicKnapsack< long > integer_knapsack( const Knapsack & kp )
{
 BasicKnapsack< long > ki;
 ki.p = kp.p;
 ki.w = kp.w;
 ki.cap = kp.cap;
 ki.explicit_children = kp.explicit_children;
 return( ki );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the B&B with values other than double, by every scheme of the core with
// 1 and 3 workers: the knapsack with integer values (bounds rounded down),
// in both selections that order by the bound and with a node limit and an
// initial incumbent, and with the exact values of Rational (by the farm
// too); the optima are those of the dynamic program,
// and the integer bounds evaluate no more nodes than the double ones in
// the sequential best-first search, summed over the instances

void test_bb_values( const std::string & )
{
 static_assert( detail::bb_lock_free< long >::value );
 std::uint64_t ev_double = 0 , ev_long = 0;
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const long opt = kp.optimum();
  const std::string tag = " knapsack " + std::to_string( i++ );
  const BasicKnapsack< long > ki = integer_knapsack( kp );
  const ExactKnapsack ke{ kp };

  ev_double += branch_and_bound( kp ).stats.evaluated;
  const auto s = branch_and_bound( ki );
  ev_long += s.stats.evaluated;
  check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
         s.node && ki.check( *s.node , s.value ) && ( s.gap == 0 ) &&
         ( s.bound == opt ) , "B&B long" + tag );

  for( Selection sel : { Selection::depth_first , Selection::plunging ,
                         Selection::best_estimate } ) {
   BasicBBOptions< long > o;
   o.selection = sel;
   const auto t = branch_and_bound( ki , o );
   check( ( t.status == BBStatus::optimal ) && ( t.value == opt ) ,
          "B&B long" + tag + " selection " + std::to_string( int( sel ) ) );
   }

  BasicBBOptions< long > o;
  o.node_limit = 3;
  const auto l = branch_and_bound( ki , o );
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt ) &&
         ( l.stats.evaluated == 3 ) , "B&B long node limit" + tag );

  o = BasicBBOptions< long >();
  o.incumbent = opt;
  long seen = 0;
  o.on_incumbent = [ & ]( long v ) { seen = v; };
  const auto z = branch_and_bound( ki , o );
  check( ( z.status == BBStatus::optimal ) && ( z.value == opt ) &&
         ! z.node && ( seen == 0 ) , "B&B long with the optimum given" +
         tag );

  const auto x = branch_and_bound( ke );
  check( ( x.status == BBStatus::optimal ) && ( x.value == Rational( opt ) ) &&
         x.node && kp.check( *x.node , double( opt ) ) &&
         ( x.gap == Rational( 0 ) ) , "B&B Rational" + tag );

  for( std::string_view name : { "shared" , "distributed" ,
                                 "deterministic" , "selfsplit" ,
                                 "stealing" } )
   for( unsigned p : { 1u , 3u } ) {
    const BBScheme sch = bb_scheme( name );
    const std::string t = std::string( name ) + tag + " p " +
                          std::to_string( p );
    BasicParallelBBOptions< long > ol;
    ol.threads = p;
    const auto r = parallel_branch_and_bound( ki , sch , ol );
    check( ( r.status == BBStatus::optimal ) && ( r.value == opt ) &&
           r.node && ki.check( *r.node , r.value ) ,
           "B&B long scheme " + t );
    BasicParallelBBOptions< Rational > oq;
    oq.threads = p;
    const auto q = parallel_branch_and_bound( ke , sch , oq );
    check( ( q.status == BBStatus::optimal ) &&
           ( q.value == Rational( opt ) ) && q.node &&
           kp.check( *q.node , double( opt ) ) ,
           "B&B Rational scheme " + t );
    const auto ml = parallel_branch_and_bound( ki , BBScheme::multiqueue ,
                                               ol );
    bool refused = false;
    try {
     ( void ) parallel_branch_and_bound( ke , BBScheme::multiqueue , oq );
     }
    catch( const std::invalid_argument & ) {
     refused = true;
     }
    check( ( ml.status == BBStatus::optimal ) && ( ml.value == opt ) &&
           refused , "B&B long multiqueue, Rational refused " + t );
    const auto fl = FarmBranchAndBound< BasicKnapsack< long > >( ki , ol )
                    .solve();
    const auto fq = FarmBranchAndBound< ExactKnapsack >( ke , oq ).solve();
    check( ( fl.status == BBStatus::optimal ) && ( fl.value == opt ) &&
           ( fq.status == BBStatus::optimal ) &&
           ( fq.value == Rational( opt ) ) , "B&B farm " + t );
    }
  }
 check( ev_long <= ev_double , "B&B long evaluates no more nodes (" +
        std::to_string( ev_long ) + " against " +
        std::to_string( ev_double ) + ")" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the incumbent given to the workers by cutoff( z ): the sequential, the
// shared, the distributed, the stealing and the deterministic B&B find the
// optimum with workers that fence the nodes by it, which they do (some
// evaluated node is infeasible, which none is in the knapsack alone), z
// never exceeds the optimum, and the outcome of the deterministic B&B is
// the same with 1 to 8 workers

void test_bb_cutoff( const std::string & )
{
 std::atomic< bool > wrong = false;
 bool fenced = false;
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  const CutoffKnapsack ck{ kp , wrong , opt };
  const auto good = [ & ]( const auto & s ) {
   fenced = fenced || ( s.stats.infeasible > 0 );
   return( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
           s.node && kp.check( *s.node , s.value ) );
   };
  BBOptions so;
  check( good( branch_and_bound( ck , so ) ) , "B&B by cutoff" + tag );
  ParallelBBOptions o;
  o.threads = 4;
  check( good( SharedBranchAndBound< CutoffKnapsack >( ck , o ).solve() ) ,
         "shared B&B by cutoff" + tag );
  check( good( DistributedBranchAndBound< CutoffKnapsack >( ck , o )
               .solve() ) , "distributed B&B by cutoff" + tag );
  check( good( StealingBranchAndBound< CutoffKnapsack >( ck , o ).solve() ) ,
         "stealing B&B by cutoff" + tag );
  o.round = 4;
  std::optional< std::tuple< double , double , std::uint64_t ,
                             std::uint64_t > > ref;
  bool ok = true , same = true;
  for( unsigned p : { 1u , 2u , 3u , 8u } ) {
   o.threads = p;
   const auto s = DeterministicBranchAndBound< CutoffKnapsack >( ck , o )
                   .solve();
   ok = ok && good( s );
   const auto out = std::make_tuple( s.value , s.bound , s.stats.evaluated ,
                                     s.stats.infeasible );
   if( ! ref )
    ref = out;
   else
    same = same && ( out == *ref );
   }
  check( ok , "deterministic B&B by cutoff" + tag );
  check( same , "deterministic B&B by cutoff" + tag +
                ", the same with 1 to 8 workers" );
  }
 check( fenced , "B&B by cutoff: some node fenced by the incumbent" );
 check( ! wrong , "B&B by cutoff: no incumbent above the optimum" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef PBEFS_HAS_MPI

// the B&B over MPI on the knapsack against the dynamic program, with every
// selection, in both senses (the minimization has no hash( n ), and hence
// random owners), with a relative tolerance and with a node limit, its
// leaves, its limit cases, and with integer and rational values, which
// travel as their bytes

void test_bb_mpi( const std::string & )
{
 const auto solve = []( const auto & q , const ParallelBBOptions & o ,
                        std::optional< Leaves > seeds ) {
  MPIBranchAndBoundSearch< std::decay_t< decltype( q ) > > b( q , o );
  return( seeds ? b.solve( std::move( *seeds ) ) : b.solve() );
  };
 int i = 0;
 for( const Knapsack & kp : knapsacks() ) {
  const double opt = double( kp.optimum() );
  const std::string tag = " knapsack " + std::to_string( i++ );
  for( Selection sel : { Selection::best_first , Selection::depth_first ,
                         Selection::breadth_first ,
                         Selection::best_estimate , Selection::plunging } ) {
   ParallelBBOptions o;
   o.selection = sel;
   o.batch = 4;
   const auto s = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
   check( ( s.status == BBStatus::optimal ) && ( s.value == opt ) &&
          s.node && kp.check( *s.node , s.value ) && ( s.gap <= 1e-9 ) ,
          "MPI B&B" + tag + " selection " + std::to_string( int( sel ) ) );
   }

  ParallelBBOptions o;
  const MinKnapsack mk{ kp };
  const auto m = MPIBranchAndBoundSearch< MinKnapsack >( mk , o ).solve();
  check( ( m.status == BBStatus::optimal ) && ( m.value == - opt ) &&
         m.node && kp.check( *m.node , - m.value ) ,
         "MPI B&B minimization" + tag );

  o.rel_gap = 0.05;
  const auto r = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
  check( ( r.status == BBStatus::optimal ) &&
         ( r.value >= opt * 0.95 - 1e-9 ) && ( r.bound >= opt - 1e-9 ) ,
         "MPI B&B within a relative gap of 5%" + tag );

  o = ParallelBBOptions();
  o.node_limit = 3;
  const auto l = MPIBranchAndBoundSearch< Knapsack >( kp , o ).solve();
  check( ( l.status == BBStatus::node_limit ) && ( l.bound >= opt - 1e-9 ) ,
         "MPI B&B node limit" + tag );

  o = ParallelBBOptions();
  o.batch = 4;
  check_frontier( kp , o , solve , "MPI B&B" + tag , true );

  const BasicKnapsack< long > ki = integer_knapsack( kp );
  BasicParallelBBOptions< long > ol;
  ol.batch = 4;
  const auto si = MPIBranchAndBoundSearch< BasicKnapsack< long > >( ki , ol )
                  .solve();
  check( ( si.status == BBStatus::optimal ) && ( si.value == long( opt ) ) &&
         si.node && ki.check( *si.node , si.value ) && ( si.gap == 0 ) ,
         "MPI B&B long" + tag );
  const ExactKnapsack ke{ kp };
  BasicParallelBBOptions< Rational > oq;
  oq.batch = 4;
  const auto sq = MPIBranchAndBoundSearch< ExactKnapsack >( ke , oq ).solve();
  check( ( sq.status == BBStatus::optimal ) &&
         ( sq.value == Rational( (long long)( opt ) ) ) && sq.node &&
         kp.check( *sq.node , opt ) , "MPI B&B Rational" + tag );
  }

 const auto run = []( const auto & q , const ParallelBBOptions & x ) {
  return( MPIBranchAndBoundSearch< std::decay_t< decltype( q ) > >( q , x )
          .solve() );
  };
 ParallelBBOptions o;
 check_bb_edges( knapsacks().front() , o , run , "MPI B&B" );
 o.batch = 0;
 o.flush_period = 0;
 check_bb_edges( knapsacks().front() , o , run , "MPI B&B parameters 0" );
 }

#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_bb_parallel< 0 , ADS::ThreadPool >( const std::string & );
template void test_bb_parallel< 1 , ADS::ThreadPool >( const std::string & );
template void test_bb_parallel< 3 , ADS::ThreadPool >( const std::string & );
template void test_bb_parallel< 4 , ADS::ThreadPool >( const std::string & );
template void test_bb_parallel< 2 , ADS::ThreadPool >( const std::string & );
template void test_bb_det< ADS::ThreadPool >( const std::string & );
template void test_bb_split< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_bb_parallel< 0 , ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_bb_parallel< 1 , ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_bb_parallel< 3 , ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_bb_parallel< 4 , ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_bb_parallel< 2 , ADS::OpenMPExecutor >
   ( const std::string & );
 template void test_bb_det< ADS::OpenMPExecutor >( const std::string & );
 template void test_bb_split< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_bb_parallel< 0 , ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_bb_parallel< 1 , ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_bb_parallel< 3 , ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_bb_parallel< 2 , ADS::FastFlowExecutor >
   ( const std::string & );
 template void test_bb_det< ADS::FastFlowExecutor >( const std::string & );
 template void test_bb_split< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- End File bb.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
