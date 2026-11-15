/*--------------------------------------------------------------------------*/
/*-------------------------- File incremental.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the incremental
 * searches: LPA* and D* Lite on grids, gridworlds and digraphs whose costs
 * change, and HD-LPA* and HD-D* Lite on the workers of HDA*.
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
// true if the path of sol is one of d from the start to a goal of the cost
// of sol and of the cost of A* from scratch, or if both find none

template< class D , class S >
static bool same_as_astar( const D & d , const S & sol )
{
 const auto a = astar( d );
 if( a.solved() != sol.solved() )
  return( false );
 if( ! a.solved() )
  return( true );
 if( ! same_cost( sol.cost , a.cost ) || sol.path.empty() ||
     ! ( sol.path.front() == d.initial_state() ) ||
     ! d.is_goal( sol.path.back() ) )
  return( false );
 cost_t< D > sum = 0;
 for( std::size_t i = 1 ; i < sol.path.size() ; ++i ) {
  cost_t< D > best = infinite_cost< cost_t< D > >();
  d.for_each_successor( sol.path[ i - 1 ] , d.heuristic( sol.path[ i - 1 ] ) ,
                        [ & ]( const state_t< D > & t , cost_t< D > c ,
                               cost_t< D > ) {
   if( ( t == sol.path[ i ] ) && ( c < best ) )
    best = c;
   } );
  if( best == infinite_cost< cost_t< D > >() )
   return( false );
  sum += best;
  }
 return( same_cost( sum , sol.cost ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// a random digraph of 6 to 25 nodes, out-degrees 1 to 3, costs 1 to 5, 1 to
// 3 goals among the candidates cands (the goals and two more nodes), whose
// hop distance is h, or h = 0 in a quarter of them, as in the simulator of
// the design of HD-LPA*

static DynamicGraph< int > random_digraph( std::mt19937_64 & r ,
                                           std::vector< std::uint32_t > &
                                           cands )
{
 const auto n = std::uint32_t( 6 + random_below( r , 20 ) );
 DynamicGraph< int > g( n , 0 );
 for( std::uint32_t u = 0 ; u < n ; ++u ) {
  const auto d = 1 + random_below( r , 3 );
  for( std::uint64_t k = 0 ; k < d ; ++k ) {
   const auto v = std::uint32_t( random_below( r , n ) );
   if( v != u )
    g.add_edge( u , v , 1 + int( random_below( r , 5 ) ) );
   }
  }
 cands.clear();
 const auto goals = 1 + random_below( r , 3 );
 for( std::uint64_t k = 0 ; k < goals + 2 ; ++k ) {
  const auto t = std::uint32_t( 1 + random_below( r , n - 1 ) );
  if( std::find( cands.begin() , cands.end() , t ) == cands.end() )
   cands.push_back( t );
  }
 for( std::size_t k = 0 ; ( k < goals ) && ( k < cands.size() ) ; ++k )
  g.set_goal( cands[ k ] , true );
 // the hop distance to the nearest candidate, over all the arcs
 std::vector< int > hop( n , -1 );
 std::vector< std::uint32_t > q( cands );
 for( const auto t : cands )
  hop[ t ] = 0;
 std::vector< std::vector< std::uint32_t > > in( n );
 for( std::size_t e = 0 ; e < g.edges() ; ++e )
  in[ g.head( e ) ].push_back( g.tail( e ) );
 for( std::size_t i = 0 ; i < q.size() ; ++i )
  for( const auto u : in[ q[ i ] ] )
   if( hop[ u ] < 0 ) {
    hop[ u ] = hop[ q[ i ] ] + 1;
    q.push_back( u );
    }
 const bool zero = random_below( r , 4 ) == 0;
 for( std::uint32_t v = 0 ; v < n ; ++v )
  g.set_heuristic( v , zero || ( hop[ v ] < 0 ) ? 0 : hop[ v ] );
 g.drain_events( []( const auto & ) {} );
 return( g );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the changes of an episode of a random digraph: 1 to 4 arcs, a fifth of
// them removed and the others given a cost of 1 to 5, and, if goals, a
// candidate that becomes a goal or stops being one, keeping one at least

static void digraph_change( DynamicGraph< int > & g ,
                            const std::vector< std::uint32_t > & cands ,
                            bool goals , std::mt19937_64 & r )
{
 const auto k = 1 + random_below( r , 4 );
 for( std::uint64_t j = 0 ; j < k ; ++j ) {
  const auto e = random_below( r , g.edges() );
  g.set_cost( e , random_below( r , 5 ) == 0 ? infinite_cost< int >()
                  : 1 + int( random_below( r , 5 ) ) );
  }
 if( goals ) {
  const auto t = cands[ random_below( r , cands.size() ) ];
  if( ! g.is_goal( t ) || ( g.goals().size() > 1 ) )
   g.set_goal( t , ! g.is_goal( t ) );
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// LPA* on grids whose costs change, told by the events of the cells or by
// update() of the 3 x 3 cells around each, against A* from scratch after
// each change

template< int Moves >
static void check_lpa_grid( std::uint64_t seed , bool typed )
{
 using G = Grid< Moves >;
 auto m = GridMap::random( 60 , 40 , 0.2 , seed , 5 );
 G g( m , 0 , 39 , 59 , 39 );
 LPAStar< G > lpa( g );
 const std::string t = " grid " + std::to_string( Moves ) + ", seed " +
                       std::to_string( seed ) + ( typed ? ", typed"
                                                        : ", update()" );
 std::uint64_t replan = 0 , scratch = 0;
 int bad = 0;
 std::mt19937_64 rng( seed );
 for( int round = 0 ; round <= 10 ; ++round ) {
  if( round > 0 ) {
   std::vector< std::pair< int , int > > changed;
   for( int k = 0 ; k < 15 ; ++k ) {
    const int x = int( rng() % 60 ) , y = int( rng() % 40 );
    if( ( ( x == 0 ) && ( y == 39 ) ) || ( ( x == 59 ) && ( y == 39 ) ) )
     continue;
    const int old = m.cost( x , y );
    m.set_cost( x , y , int( rng() % 6 ) );   // 0, a wall, to 5
    if( typed )
     g.cell_changed( x , y , old );
    changed.emplace_back( x , y );
    }
   if( typed )
    lpa.take_events( g );
   else {
    g.refresh();
    for( auto [ x , y ] : changed )
     for( int dy = -1 ; dy <= 1 ; ++dy )
      for( int dx = -1 ; dx <= 1 ; ++dx )
       if( ( x + dx >= 0 ) && ( x + dx < 60 ) && ( y + dy >= 0 ) &&
           ( y + dy < 40 ) )
        lpa.update( g.index( x + dx , y + dy ) );
    }
   }
  const auto s = lpa.solve();
  const auto a = astar( g );
  if( round > 0 ) {
   replan += s.stats.expanded;
   scratch += a.stats.expanded;
   }
  if( ! same_as_astar( g , s ) )
   ++bad;
  }
 check( bad == 0 , "LPA*" + t + ", 11 rounds against A*" );
 check( replan < scratch , "LPA* replanning on" + t + " (" +
        std::to_string( replan ) + " expansions, A* from scratch " +
        std::to_string( scratch ) + ")" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// LPA* in the gridworlds of figures 8 and 9 of Koenig, Likhachev and Furcy
// (AIJ 2004): 51 x 51 cells, 0.6 percent of the arcs given a cost of 1 or
// 2, or 8 cells blocked and 8 freed in a world of 20 percent of them
// blocked, after each episode; every episode against A* from scratch, the
// expansions of LPA* below an eighth of those of A* with the ties to the
// smaller g (version 1 of the paper, as LPA*), which the paper finds 23.71
// against 307.93 and 15.56 against 241.77, and the probability that the
// cost of the path changes near the 7.3 percent of its figure 11

static void check_lpa_paper( bool blocks )
{
 const std::string t = blocks ? "figure 9" : "figure 8";
 std::uint64_t lpa = 0 , a1 = 0 , episodes = 0 , moved = 0;
 int bad = 0;
 for( int world = 0 ; world < 10 ; ++world ) {
  std::mt19937_64 rng( 1000 + std::uint64_t( world ) );
  auto g = DynamicGraph< int >::gridworld( 51 , blocks ? 0.2 : 0 , rng );
  LPAStar< DynamicGraph< int > > s( g );
  int last = s.solve().cost;
  for( int e = 0 ; e < 50 ; ++e ) {
   g.perturb( 0.006 , blocks ? 8 : 0 , rng );
   s.take_events( g );
   const auto x = s.solve();
   const auto y = astar( g , OpenKind::heap_low_g );
   bad += ( x.solved() != y.solved() ) ||
          ( x.solved() && ( x.cost != y.cost ) );
   lpa += x.stats.expanded;
   a1 += y.stats.expanded;
   moved += x.cost != last;
   last = x.cost;
   ++episodes;
   }
  }
 check( bad == 0 , "LPA* " + t + ", 500 episodes against A*" );
 check( 8 * lpa < a1 , "LPA* " + t + ": " + std::to_string( double( lpa ) /
        double( episodes ) ) + " expansions per episode, A* " +
        std::to_string( double( a1 ) / double( episodes ) ) );
 if( ! blocks )
  check( ( 100 * moved > 4 * episodes ) && ( 100 * moved < 11 * episodes ) ,
         "LPA* figure 8: the cost changes in " + std::to_string( moved ) +
         " of " + std::to_string( episodes ) + " episodes" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// D* Lite: a robot that walks its path on a gridworld of LPA* whose costs
// change after each step, against A* from the robot

static void check_dstar_lite( void )
{
 int bad = 0 , steps = 0;
 for( int world = 0 ; world < 6 ; ++world ) {
  std::mt19937_64 rng( 77 + std::uint64_t( world ) );
  auto g = DynamicGraph< int >::gridworld( 31 , world % 2 ? 0.2 : 0 ,
                                          rng );
  const auto goal = g.goals().front();
  const auto start = g.initial_state();
  DStarLite< DynamicGraph< int > > d( g , goal , start );
  for( int step = 0 ; step < 40 ; ++step ) {
   const auto x = d.solve();
   RobotDomain< DynamicGraph< int > > r( g , goal , d.robot() );
   const auto a = astar( r );
   ++steps;
   if( ( a.solved() != x.solved() ) || ( a.solved() && ( a.cost != x.cost ) )
       || ( x.solved() && ! ( x.path.front() == d.robot() ) ) ) {
    ++bad;
    break;
    }
   if( ! x.solved() || ( x.path.size() < 2 ) )
    break;
   d.move_to( x.path[ 1 ] );
   g.perturb( 0.01 , world % 2 ? 4 : 0 , rng );
   d.take_events();
   }
  }
 check( bad == 0 , "D* Lite, " + std::to_string( steps ) +
                   " steps of a robot against A*" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// LPA* on random digraphs whose arcs and goals change, against A*

static void check_lpa_digraphs( void )
{
 int bad = 0;
 std::vector< std::uint32_t > cands;
 for( std::uint64_t t = 0 ; t < 300 ; ++t ) {
  std::mt19937_64 r( 5 * 1000003 + t );
  auto g = random_digraph( r , cands );
  g.set_heuristic_goals( cands );
  LPAStar< DynamicGraph< int > > s( g );
  bad += ! same_as_astar( g , s.solve() );
  for( int e = 0 ; e < 8 ; ++e ) {
   digraph_change( g , cands , true , r );
   s.take_events( g );
   bad += ! same_as_astar( g , s.solve() );
   }
  }
 check( bad == 0 , "LPA* on 300 random digraphs, 9 episodes each" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
void test_lpa( const std::string & )
{
 for( std::uint64_t seed = 1 ; seed <= 3 ; ++seed ) {
  check_lpa_grid< 4 >( seed , true );
  check_lpa_grid< 8 >( seed , true );
  check_lpa_grid< 4 >( seed , false );
  }
 check_lpa_paper( false );
 check_lpa_paper( true );
 check_lpa_digraphs();
 check_dstar_lite();

 const auto dist = eight_distances( EightPuzzle( EightPuzzle::korf_goal() ) );
 for( const auto & [ s , dst ] : eight_sample( dist ) ) {
  const EightPuzzle e( EightPuzzle::unpack( s ) );
  LPAStar< EightPuzzle > lpa( e );
  const auto x = lpa.solve();
  check_solution( e , x , dst , "LPA* 8-puzzle" );
  const auto y = lpa.solve();
  check( ( y.cost == x.cost ) && ( y.stats.expanded == 0 ) ,
         "LPA* with no change expands nothing" );
  }

 // a directed graph: 0 -> 1 -> 3 costs 2, 0 -> 2 -> 3 costs 4
 ExplicitGraph< int > g( 4 , 0 , 3 );
 g.add_edge( 0 , 1 , 1 , true );
 g.add_edge( 1 , 3 , 1 , true );
 g.add_edge( 0 , 2 , 2 , true );
 g.add_edge( 2 , 3 , 2 , true );
 g.add_edge( 3 , 0 , 1 , true );
 LPAStar< ExplicitGraph< int > > lg( g );
 check_solution( g , lg.solve() , 2 , "LPA* directed graph" );
 auto rg = romania( Arad );
 check_solution( rg , LPAStar< ExplicitGraph< int > >( rg ).solve() , 418 ,
                 "LPA* Romania" );
 check_solution( rg , LPAStar< ExplicitGraph< int > , LPABuckets >( rg )
                      .solve() , 418 , "LPA* Romania, buckets" );
 GridMap m( 30 , 30 );
 for( int y = 0 ; y < 30 ; ++y )
  m.set_cost( 15 , y , 0 );
 const Grid< 4 > wall( m , 0 , 0 , 29 , 29 );
 check( LPAStar< Grid< 4 > >( wall ).solve().status == Status::unsolvable ,
        "LPA* proves the wall unsolvable" );
 const EightPuzzle goal( EightPuzzle::korf_goal() );
 const auto z = LPAStar< EightPuzzle >( goal ).solve();
 check( z.solved() && ( z.cost == 0 ) && ( z.path.size() == 1 ) ,
        "LPA* initial state that is the goal" );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the owner of a node of a DynamicGraph given by a table, so that a test
// places the nodes on the workers

template< class D >
struct TableDistribution {
 std::vector< unsigned > o;                   // the owner of each node
 [[nodiscard]] unsigned owner( std::uint32_t s , unsigned p ) const {
  return( o[ s ] % p );
  }
 [[nodiscard]] static std::string name( void ) { return( "T" ); }
 };
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// HD-LPA* on gridworlds whose costs change and on random digraphs whose
// arcs and goals change, every episode against A* from scratch with the
// invariants audited at each barrier, with 1 to 8 workers, both rules of
// phase 1, the three policies, batches of 1 and 64, with and without the
// window on the least key; with Freeze::all no
// state expanded underconsistent twice in an episode (theorem 3), with
// one worker the expansions of LPA* within a tenth; the episode of the
// count to infinity of the design, which ends; HD-D* Lite with a robot
// that walks, and the restarts when a repair goes over its budget

template< class Executor >
void test_hdlpa( const std::string & )
{
 using G = DynamicGraph< int >;
 using P = IncrementalOptions::Policy;
 using F = IncrementalOptions::Freeze;
 const std::string rt = ", " + executor_name< Executor >();
 for( unsigned p : { 1u , 2u , 3u , 8u } ) {
  Executor ex( p );
  for( F fz : { F::all , F::raised_parent } )
   for( P pol : { P::repair , P::adaptive , P::restart } ) {
    if( ( pol != P::repair ) && ( fz == F::raised_parent ) )
     continue;
    int bad = 0;
    std::uint32_t raises = 0;
    std::uint64_t hx = 0 , lx = 0 , restarts = 0;
    for( int world = 0 ; world < 6 ; ++world ) {
     std::mt19937_64 rng( 300 + std::uint64_t( world ) );
     auto g = G::gridworld( world < 3 ? 21 : 31 , world % 2 ? 0.2 : 0 ,
                            rng );
     const StateHashDistribution< G > dist( g , 3 );
     HDLPAOptions o;
     o.threads = p;
     o.audit = true;
     o.freeze = fz;
     o.policy = pol;
     o.batch = world % 3 ? 64 : 1;
     o.flush_period = world % 3 ? 64 : 1;
     o.window = world % 2 ? 2 : -1;
     HDLPAStar< G , StateHashDistribution< G > , Executor > h( g , dist , o );
     h.set_executor( ex );
     LPAStar< G > l( g );
     for( int e = 0 ; e <= 20 ; ++e ) {
      if( e > 0 ) {
       g.perturb( 0.02 , world % 2 ? 4 : 0 , rng );
       g.drain_events( [ & ]( const auto & x ) {
        if constexpr( requires { x.from; } ) {
         h.edge_changed( x.u , x.v , x.from , x.to );
         l.edge_changed( x.u , x.v , x.from , x.to );
         }
        } );
       }
      const auto x = h.solve();
      const auto y = l.solve();
      bad += ! same_as_astar( g , x );
      raises = std::max( raises , h.max_raises() );
      if( e > 0 ) {
       hx += x.stats.expanded;
       lx += y.stats.expanded;
       restarts += x.stats.restarts;
       }
      }
     }
    const std::string t = " p " + std::to_string( p ) + ( fz == F::all ?
                          ", freeze all" : ", freeze raised parent" ) +
                          ( pol == P::repair ? ", repair" : pol ==
                            P::restart ? ", restart" : ", adaptive" ) + rt;
    check( bad == 0 , "HD-LPA* gridworlds, 126 episodes against A*" + t );
    if( fz == F::all )
     check( raises <= 1 , "HD-LPA* no state raised twice in an episode" + t );
    if( pol == P::restart )
     check( restarts == 120 , "HD-LPA* restarts every episode" + t );
    if( ( p == 1 ) && ( pol == P::repair ) )
     check( 10 * hx <= 11 * lx , "HD-LPA* one worker: " + std::to_string( hx )
            + " expansions, LPA* " + std::to_string( lx ) + t );
    }
  }

 // random digraphs whose arcs and goals change, batches of 1
 for( unsigned p : { 2u , 4u , 6u } ) {
  Executor ex( p );
  int bad = 0;
  std::vector< std::uint32_t > cands;
  for( std::uint64_t t = 0 ; t < 150 ; ++t ) {
   std::mt19937_64 r( 7 * 1000003 + t + 1000 * p );
   auto g = random_digraph( r , cands );
   g.set_heuristic_goals( cands );
   const StateHashDistribution< G > dist( g , t );
   HDLPAOptions o;
   o.threads = p;
   o.audit = true;
   o.batch = 1;
   o.flush_period = 1;
   o.freeze = t % 2 ? F::all : F::raised_parent;
   o.window = t % 3 == 0 ? -1 : double( t % 3 - 1 );
   HDLPAStar< G , StateHashDistribution< G > , Executor > h( g , dist , o );
   h.set_executor( ex );
   bad += ! same_as_astar( g , h.solve() );
   for( int e = 0 ; e < 8 ; ++e ) {
    digraph_change( g , cands , true , r );
    h.take_events( g );
    bad += ! same_as_astar( g , h.solve() );
    }
   }
  check( bad == 0 , "HD-LPA* 150 random digraphs, 9 episodes, p " +
                    std::to_string( p ) + rt );
  }

 // s0 -> u (1), u <-> v (1), v -> t (1), u and v on two workers: the
 // episode that removes s0 -> u, after which t is unreachable, ends
 Executor ex2( 2 ) , ex3( 3 ) , ex4( 4 );
 for( int k = 0 ; k < 50 ; ++k ) {
  G g( 4 , 0 );
  const auto e = g.add_edge( 0 , 1 , 1 );
  g.add_edge( 1 , 2 , 1 );
  g.add_edge( 2 , 1 , 1 );
  g.add_edge( 2 , 3 , 1 );
  g.set_goal( 3 , true );
  g.drain_events( []( const auto & ) {} );
  const TableDistribution< G > dist{ { 0 , 0 , 1 , unsigned( k % 2 ) } };
  HDLPAOptions o;
  o.threads = 2;
  o.audit = true;
  o.batch = 1;
  o.flush_period = 1;
  o.policy = P::repair;
  o.freeze = k % 4 < 2 ? F::all : F::raised_parent;
  HDLPAStar< G , TableDistribution< G > , Executor > h( g , dist , o );
  h.set_executor( ex2 );
  const bool first = h.solve().cost == 3;
  g.set_cost( e , infinite_cost< int >() );
  h.take_events( g );
  const auto x = h.solve();
  if( ! first || x.solved() ) {
   check( false , "HD-LPA* the count to infinity of the design" + rt );
   break;
   }
  if( k == 49 )
   check( true , "HD-LPA* the count to infinity of the design ends, 50 runs"
                 + rt );
  }

 // grids of cells with 4 and 8 moves, told by the events of the cells
 {
  auto m = GridMap::random( 60 , 40 , 0.2 , 4 , 5 );
  Grid< 8 > g( m , 0 , 39 , 59 , 39 );
  const AbstractStateDistribution< Grid< 8 > > dist( g );
  HDLPAOptions o;
  o.threads = 3;
  o.audit = true;
  HDLPAStar< Grid< 8 > , AbstractStateDistribution< Grid< 8 > > , Executor >
   h( g , dist , o );
  h.set_executor( ex3 );
  std::mt19937_64 rng( 4 );
  int bad = ! same_as_astar( g , h.solve() );
  for( int round = 0 ; round < 10 ; ++round ) {
   for( int k = 0 ; k < 15 ; ++k ) {
    const int x = int( rng() % 60 ) , y = int( rng() % 40 );
    if( ( ( x == 0 ) && ( y == 39 ) ) || ( ( x == 59 ) && ( y == 39 ) ) )
     continue;
    const int old = m.cost( x , y );
    m.set_cost( x , y , int( rng() % 6 ) );
    g.cell_changed( x , y , old );
    }
   h.take_events( g );
   bad += ! same_as_astar( g , h.solve() );
   }
  check( bad == 0 , "HD-LPA* grid 8, 11 rounds against A*, p 3" + rt );
  }

 // HD-D* Lite: a robot that walks its path, with 3 workers
 {
  int bad = 0 , steps = 0;
  for( int world = 0 ; world < 4 ; ++world ) {
   std::mt19937_64 rng( 91 + std::uint64_t( world ) );
   auto g = G::gridworld( 31 , world % 2 ? 0.2 : 0 , rng );
   const auto goal = g.goals().front();
   const StateHashDistribution< G > dist( g , 5 );
   HDLPAOptions o;
   o.threads = 3;
   o.audit = true;
   HDDStarLite< G , StateHashDistribution< G > , Executor > d(
    g , dist , goal , g.initial_state() , o );
   d.set_executor( ex3 );
   for( int step = 0 ; step < 30 ; ++step ) {
    const auto x = d.solve();
    RobotDomain< G > r( g , goal , d.robot() );
    const auto a = astar( r );
    ++steps;
    if( ( a.solved() != x.solved() ) ||
        ( a.solved() && ( a.cost != x.cost ) ) ) {
     ++bad;
     break;
     }
    if( ! x.solved() || ( x.path.size() < 2 ) )
     break;
    d.move_to( x.path[ 1 ] );
    g.perturb( 0.01 , world % 2 ? 4 : 0 , rng );
    d.take_events();
    }
   }
  check( bad == 0 , "HD-D* Lite, " + std::to_string( steps ) +
                    " steps of a robot against A*, p 3" + rt );
  }

 // a budget too small for any repair: every repair stopped and restarted
 {
  std::mt19937_64 rng( 17 );
  auto g = G::gridworld( 31 , 0 , rng );
  const StateHashDistribution< G > dist( g , 3 );
  HDLPAOptions o;
  o.threads = 4;
  o.audit = true;
  o.beta = 0.01;
  o.rho = 1e9;
  o.gamma = 1e9;
  HDLPAStar< G , StateHashDistribution< G > , Executor > h( g , dist , o );
  h.set_executor( ex4 );
  int bad = ! same_as_astar( g , h.solve() );
  std::uint64_t restarts = 0;
  for( int e = 0 ; e < 20 ; ++e ) {
   g.perturb( 0.05 , 0 , rng );
   h.take_events( g );
   const auto x = h.solve();
   bad += ! same_as_astar( g , x );
   restarts += x.stats.restarts;
   }
  check( ( bad == 0 ) && ( restarts > 0 ) , "HD-LPA* repairs over budget "
         "restarted (" + std::to_string( restarts ) + " of 20), p 4" + rt );
  }

 // the 8-puzzle: a search from scratch, then one with no change
 const EightPuzzle e( EightPuzzle::unpack( eight_sample( eight_distances(
                       EightPuzzle( EightPuzzle::korf_goal() ) ) ).back()
                       .first ) );
 const ZobristDistribution< EightPuzzle > z( e );
 HDLPAOptions o;
 o.threads = 4;
 HDLPAStar< EightPuzzle , ZobristDistribution< EightPuzzle > , Executor >
  h8( e , z , o );
 h8.set_executor( ex4 );
 const auto x = h8.solve();
 check( same_as_astar( e , x ) , "HD-LPA* 8-puzzle, p 4" + rt );
 const auto y = h8.solve();
 check( y.solved() && ( y.cost == x.cost ) && ( y.stats.expanded == 0 ) ,
        "HD-LPA* with no change expands nothing" + rt );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- INSTANTIATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template void test_hdlpa< ADS::ThreadPool >( const std::string & );
#ifdef PBEFS_HAS_OPENMP
 template void test_hdlpa< ADS::OpenMPExecutor >( const std::string & );
#endif
#ifdef PBEFS_HAS_FASTFLOW
 template void test_hdlpa< ADS::FastFlowExecutor >( const std::string & );
#endif

/*--------------------------------------------------------------------------*/
/*------------------------ End File incremental.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
