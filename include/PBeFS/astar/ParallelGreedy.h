/*--------------------------------------------------------------------------*/
/*------------------------- File ParallelGreedy.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template ParallelGreedy, the parallel greedy
 * best-first searches with one OPEN and one CLOSED list shared by the
 * workers: KPGBFS, PUHF, PUHF2, PUHF3, PUHF4, SPUHF and OBAT, each with or
 * without the separate generation and evaluation of the successors (SGE),
 * on an executor of ADS or in synchronous rounds.
 *
 * All of them follow the schema of Constrained Parallel GBFS of Shimoda
 * and Fukunaga (arXiv:2408.05682, algorithm 1): a worker removes the node
 * of least h from OPEN only if it satisfies the constraint of the rule,
 * generates and evaluates its successors, and inserts those not in CLOSED
 * into OPEN all together, so that they are ordered as GBFS orders them.
 * With \f$ s_j \f$ the node that the worker j is expanding
 * (\f$ h( s_j ) = \infty \f$ for an idle one) and \f$ \underline{h} =
 * \min_j h( s_j ) \f$, the constraints are:
 *
 * - KPGBFS: none, the K-parallel GBFS of Vidal, Bordeaux and Hamadi
 *   (SoCS 2010);
 *
 * - PUHF: the node is marked certain, which a successor \f$ s' \f$ of a
 *   node s is if \f$ h( s' ) = \min_{ t \in \{ s \} \cup succ( s ) } h( t )
 *   \f$, and the top of OPEN is when no worker is expanding [see Kuroiwa
 *   and Fukunaga, ICAPS 2020];
 *
 * - PUHF2: the node is certain, or \f$ h \leq \underline{h} \f$, which
 *   makes it certain; PUHF3 marks certain, in addition, all the nodes of
 *   OPEN with \f$ h = \underline{h} \f$ after each removal; PUHF4 adds the
 *   criterion of PUHF to PUHF3 and prefers the certain nodes among those
 *   of equal h [see Shimoda and Fukunaga, SoCS 2023];
 *
 * - SPUHF: PUHF, and a worker that has no certain node expands a node of
 *   a speculative OPEN of its own, whose states share only the evaluations
 *   with the main search; a goal found by the speculation ends the search
 *   [see Kuroiwa and Fukunaga, ICAPS 2020];
 *
 * - OBAT, One Bench At a Time: \f$ h \leq \underline{h} \f$, and a node
 *   with a successor of smaller h (a b-state, possibly a progress state)
 *   enters a queue DEFERRED with its successors, which enter OPEN only
 *   when the node leaves DEFERRED, i.e., when \f$ h \f$ of the top of
 *   DEFERRED is at most that of the top of OPEN and at most
 *   \f$ \underline{h} \f$ [see Shimoda and Fukunaga, AAAI 2025,
 *   algorithm 2].
 *
 * The guarantees, relative to the sequential GBFS, are those of the
 * papers, the bench transition system (BTS) being the states that GBFS
 * expands with some breaking of the ties [see Heusner, Keller and
 * Helmert, SoCS 2017, theorem 1]: KPGBFS is pathological, i.e., for every
 * t it expands more than t times the states of GBFS on some instance [see
 * Kuroiwa and Fukunaga, ICAPS 2020]; PUHF to PUHF4 expand only states of
 * the BTS (and so does the main search of SPUHF), but they are
 * pathological too [see Shimoda and Fukunaga, AAAI 2025, section 3];
 * OBAT expands at most \f$ N + k | p | \f$ states with k workers, N
 * being the states that GBFS expands with some breaking of the ties and p
 * the path that it finds (theorem 4 there), since the nodes leave OPEN or
 * DEFERRED for the last time in an order in which GBFS may expand them
 * (theorem 2) and at most k nodes of each h of the path remain in
 * DEFERRED (theorem 3). SGE changes the order of the work but not these
 * guarantees (section 6 there).
 *
 * In synchronous rounds (ParallelGreedyOptions::rounds), the k workers
 * are emulated in the calling thread: in each round every idle worker in
 * turn takes its work under the rules, then all of them do it, then they
 * insert their results in turn, i.e., every expansion (or, with SGE, every
 * generation and every evaluation) takes one round. The search is then
 * replicable, and its counts do not depend on the machine, which the
 * tests of the bounds above need.
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

#ifndef __PBeFS_ParallelGreedy
 #define __PBeFS_ParallelGreedy
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <atomic>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <map>

#include <memory>

#include <mutex>

#include <stdexcept>

#include <vector>

#include "PBeFS/astar/GreedySearch.h"

#include "PBeFS/astar/HDAStar.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/// the rule by which a parallel greedy search removes a node from OPEN

enum class GreedyRule {
 kpgbfs ,     ///< no constraint
 puhf ,       ///< the certain nodes of PUHF
 puhf2 ,      ///< certain, or h at most that of every node in expansion
 puhf3 ,      ///< PUHF2, the nodes of the least h in expansion certain
 puhf4 ,      ///< PUHF3 with the criterion of PUHF, certain ones first
 spuhf ,      ///< PUHF, with a speculative search for the idle workers
 obat         ///< One Bench At a Time
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the options of ParallelGreedy

struct ParallelGreedyOptions : public ParallelOptions {

 GreedyRule rule = GreedyRule::kpgbfs;   ///< the rule

 /// the successors of a node are generated by the worker that removes it
 /// and evaluated by any worker (SGE)
 bool sge = false;

 /// the order among the nodes of equal h: fifo or lifo (the others are
 /// taken as fifo)
 TieBreak tie = TieBreak::fifo;

 /// the certain nodes come first among those of equal h (always with
 /// PUHF4)
 bool certain_first = false;

 /// the workers are emulated in synchronous rounds in the calling thread
 bool rounds = false;

 };  // end( struct( ParallelGreedyOptions ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// what a parallel greedy search did, besides its SearchStats

template< class State >
struct GreedyTrace {

 /// the states removed from OPEN by the main search, in the order of
 /// their removals
 std::vector< State > expanded;

 /// the states whose successors entered OPEN, in that order (with OBAT, a
 /// b-state when it leaves DEFERRED; with the other rules, a state when
 /// its successors are inserted)
 std::vector< State > committed;

 /// the states left in DEFERRED at the end (OBAT)
 std::vector< State > deferred;

 std::uint64_t speculative = 0;      ///< the speculative expansions
 bool goal_by_speculation = false;   ///< the goal was found by them
 std::uint64_t rounds = 0;           ///< the rounds, if emulated

 };  // end( struct( GreedyTrace ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ParallelGreedy --------------------------*/
/*--------------------------------------------------------------------------*/
/// the parallel greedy best-first searches with a shared OPEN list
/** The searches of the file comment. OPEN, DEFERRED and the h of the nodes
 * in expansion are under one lock, as in the papers; CLOSED is a table of
 * the states generated, split in shards under a lock each, each state
 * evaluated once (the table is also the cache of the evaluations shared
 * with the speculative search of SPUHF), and a flag tells whether a state
 * is in CLOSED in the sense of the papers, i.e., has entered OPEN. The
 * goal test is done when a node is removed. The counters: expanded are
 * the removals of the main search (plus the speculative expansions),
 * evaluated the states evaluated, sent 0. trace(), after solve(), gives
 * what GreedyTrace records, the states only if set_record( true ). */

template< Domain D , ADS::WorkerExecutor Executor = ADS::ThreadPool >
class ParallelGreedy {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search)

 explicit ParallelGreedy( const D & domain ,
                          const ParallelGreedyOptions & options = {} )
  : f_domain( domain ) , f_options( options ) {
  f_options.normalize();
  if( f_options.rule == GreedyRule::puhf4 )
   f_options.certain_first = true;
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers, which must outlive the
 /// searches; without it, each search builds one of as many threads as
 /// there are workers

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

 /// records the states in the trace if true

 void set_record( bool r ) noexcept { f_record = r; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) {
  Stopwatch clock;
  Shared sh( *this );
  f_trace = GreedyTrace< State >();
  const unsigned p = f_options.threads;
  std::vector< Worker > w;
  w.reserve( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w.emplace_back( i );

  // the root, evaluated, in CLOSED and in OPEN
  const State s0 = f_domain.initial_state();
  const Cost h0 = Cost( f_domain.heuristic( s0 ) );
  ++sh.evaluated;
  if( h0 != infinite_cost< Cost >() ) {
   bool fresh;
   NodeT * root = sh.lookup( s0 , f_domain.hash( s0 ) , h0 , fresh );
   root->g = Cost( 0 );
   root->listed = true;
   root->certain = true;
   sh.push_open( root );
   if( spuhf() )
    sh.push_spec( root , nullptr , Cost( 0 ) , 0 );
   }

  if( f_options.rounds ) {
   while( ! sh.done ) {
    ++f_trace.rounds;
    bool any = false;
    for( Worker & x : w )
     if( x.idle() ) {
      std::lock_guard< std::mutex > lock( sh.m );
      any = acquire( sh , x ) || any;
      }
    for( Worker & x : w )
     if( ! x.idle() ) {
      any = true;
      work( sh , x );
      }
    for( Worker & x : w )
     if( ! x.idle() ) {
      std::lock_guard< std::mutex > lock( sh.m );
      commit( sh , x );
      }
    if( ! any )
     sh.done = true;                          // nothing more to do
    }
   }
  else
   run_workers( f_ex , p , [ & ]( unsigned i ) { run( sh , w[ i ] ); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( const Worker & x : w ) {
   st.expanded += x.st.expanded;
   st.generated += x.st.generated;
   st.expanded_by.push_back( x.st.expanded );
   }
  st.duplicates = sh.duplicates;
  st.evaluated = sh.evaluated.load();
  st.stored = sh.stored();
  st.stored_by = { st.stored };
  f_trace.speculative = sh.speculative;
  if( f_record && obat() )
   for( auto & [ h , q ] : sh.deferred )
    for( const Deferred & d : q )
     f_trace.deferred.push_back( d.e->node->state );
  if( sh.goal ) {
   sol.status = Status::solved;
   sol.path = path( sh.goal , sh.goal_spec );
   sol.cost = path_cost( sol.path );
   f_trace.goal_by_speculation = sh.goal_spec;
   }
  else
   if( sh.limit_hit )
    sol.status = Status::node_limit;
  st.seconds = clock.seconds();
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// what the last search did

 [[nodiscard]] const GreedyTrace< State > & trace( void ) const noexcept {
  return( f_trace );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a node: its state, its h, and the fields of the main search (under
 /// the lock of OPEN) and of the speculative one (under its own lock)

 struct NodeT {
  State state;
  Cost h;
  Cost g = Cost( 0 );
  NodeT * parent = nullptr;
  std::uint8_t op = no_op;  // the operator from the parent
  bool listed = false;      // in CLOSED in the sense of the papers
  bool certain = false;     // certain (PUHF)
  bool removed = false;     // removed from OPEN
  NodeT * sparent = nullptr;  // the speculative parent
  Cost sg = Cost( 0 );      // the speculative g
  bool sseen = false;       // pushed in the speculative OPEN
  bool sclosed = false;     // expanded by the speculation
  };

 /// a successor of a node being expanded
 struct Kid {
  State s;
  Cost c;
  Cost h;                   // known from the Domain, if not with SGE
  std::uint8_t op;
  NodeT * node = nullptr;   // its node, once evaluated
  };

 /// an expansion: the node, its successors, those left to evaluate
 struct Expansion {
  NodeT * node = nullptr;
  bool speculative = false;
  bool busy = false;        // its h is among those in expansion
  std::vector< Kid > kids;
  std::size_t left = 0;
  };

 /// an entry of OPEN or DEFERRED
 struct Entry {
  NodeT * node;
  std::uint64_t seq;
  };

 /// a node in DEFERRED, with its successors
 struct Deferred {
  std::shared_ptr< Expansion > e;
  std::uint64_t seq;
  };

 /// the nodes of OPEN of one h, the certain apart
 struct Bucket {
  std::deque< Entry > certain;
  std::deque< Entry > other;
  [[nodiscard]] bool empty( void ) const noexcept {
   return( certain.empty() && other.empty() );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a worker: what it holds, and its counters

 struct Worker {
  explicit Worker( unsigned i ) : id( i ) {}
  unsigned id;
  std::shared_ptr< Expansion > exp;    // the node it expands, if any
  std::size_t kid = 0;                 // the successor it evaluates (SGE)
  bool evaluating = false;             // it evaluates a successor
  SearchStats st;
  [[nodiscard]] bool idle( void ) const noexcept { return( ! exp ); }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the data shared by the workers

 struct Shard {
  std::mutex m;
  NodePool< NodeT > pool;
  NodeTable< NodeT > table;
  };

 static constexpr std::size_t c_shards = 64;

 struct Shared {

  explicit Shared( const ParallelGreedy & a ) : f_a( a ) ,
   shards( c_shards ) {}

  const ParallelGreedy & f_a;
  std::mutex m;                                   // OPEN, DEFERRED, ...
  std::map< Cost , Bucket > open;                 // OPEN by h
  std::size_t open_n = 0;                         // its nodes
  std::map< Cost , std::deque< Deferred > > deferred;  // DEFERRED by h
  std::map< Cost , std::size_t > busy;            // h of the nodes in
                                                  // expansion
  std::size_t busy_n = 0;                         // their number
  std::deque< std::pair< std::shared_ptr< Expansion > , std::size_t > >
   unevaluated;                                   // the successors to
                                                  // evaluate (SGE)
  std::uint64_t seq = 0;
  std::uint64_t duplicates = 0;                   // successors in CLOSED
  bool done = false;
  NodeT * goal = nullptr;
  bool goal_spec = false;
  bool limit_hit = false;
  std::vector< Shard > shards;
  std::atomic< std::uint64_t > evaluated{ 0 };
  std::atomic< std::size_t > nodes{ 0 };
  std::atomic< bool > over{ false };              // done, read freely

  std::mutex sm;                                  // the speculation
  std::map< std::pair< int , Cost > , std::deque< Entry > > sopen;
  std::uint64_t speculative = 0;
  std::size_t spec_busy = 0;

  [[nodiscard]] bool open_empty( void ) const noexcept {
   return( open_n == 0 );
   }

  // the node of s, of hash hash, evaluated with h if new
  NodeT * lookup( const State & s , std::uint64_t hash , Cost h ,
                  bool & fresh ) {
   Shard & x = shards[ hash % c_shards ];
   std::lock_guard< std::mutex > lock( x.m );
   fresh = false;
   NodeT * n = x.table.find_or_insert( s , hash , [ & ] {
    fresh = true;
    NodeT * y = x.pool.allocate();
    *y = NodeT{ s , h };
    return( y );
    } ).first;
   if( fresh )
    nodes.fetch_add( 1 , std::memory_order_relaxed );
   return( n );
   }

  // the node of s, of hash hash, if generated
  NodeT * find( const State & s , std::uint64_t hash ) {
   Shard & x = shards[ hash % c_shards ];
   std::lock_guard< std::mutex > lock( x.m );
   return( x.table.find( s , hash ) );
   }

  [[nodiscard]] std::size_t stored( void ) const {
   return( nodes.load() );
   }

  // pushes n in OPEN
  void push_open( NodeT * n ) {
   Bucket & b = open[ n->h ];
   ( n->certain && f_a.splits_certain() ? b.certain : b.other ).push_back(
    Entry{ n , seq++ } );
   ++open_n;
   }

  // the bucket of the least h of OPEN, OPEN not empty
  typename std::map< Cost , Bucket >::iterator top_bucket( void ) {
   auto it = open.begin();
   while( it->second.empty() )
    it = open.erase( it );
   return( it );
   }

  // the deque and the side of the top entry of the bucket b
  std::deque< Entry > & top_side( Bucket & b , bool & back ) const {
   back = f_a.f_options.tie == TieBreak::lifo;
   if( b.certain.empty() )
    return( b.other );
   if( b.other.empty() || f_a.f_options.certain_first )
    return( b.certain );
   const Entry & c = back ? b.certain.back() : b.certain.front();
   const Entry & o = back ? b.other.back() : b.other.front();
   return( ( back ? ( c.seq > o.seq ) : ( c.seq < o.seq ) ) ? b.certain
                                                            : b.other );
   }

  // the top node of OPEN, OPEN not empty
  NodeT * top( void ) {
   bool back;
   auto & d = top_side( top_bucket()->second , back );
   return( back ? d.back().node : d.front().node );
   }

  // removes the top node of OPEN, OPEN not empty
  NodeT * pop( void ) {
   auto it = top_bucket();
   bool back;
   auto & d = top_side( it->second , back );
   NodeT * n = back ? d.back().node : d.front().node;
   if( back )
    d.pop_back();
   else
    d.pop_front();
   --open_n;
   if( it->second.empty() )
    open.erase( it );
   return( n );
   }

  // the least h of the nodes in expansion, infinite if none
  [[nodiscard]] Cost least_busy( void ) const {
   return( busy.empty() ? infinite_cost< Cost >() : busy.begin()->first );
   }

  void add_busy( Cost h ) {
   ++busy[ h ];
   ++busy_n;
   }

  void drop_busy( Expansion & e ) {
   if( ! e.busy )
    return;
   auto it = busy.find( e.node->h );
   if( --it->second == 0 )
    busy.erase( it );
   --busy_n;
   e.busy = false;
   }

  // marks certain all the nodes of OPEN with h equal to v (PUHF3)
  void certify( Cost v ) {
   auto it = open.find( v );
   if( it == open.end() )
    return;
   Bucket & b = it->second;
   for( Entry & e : b.other )
    e.node->certain = true;
   if( ! f_a.splits_certain() )
    return;
   std::deque< Entry > all;
   std::merge( b.certain.begin() , b.certain.end() , b.other.begin() ,
               b.other.end() , std::back_inserter( all ) ,
               []( const Entry & x , const Entry & y ) {
                return( x.seq < y.seq ); } );
   b.certain.swap( all );
   b.other.clear();
   }

  // pushes n, generated by parent with g, in the speculative OPEN
  void push_spec( NodeT * n , NodeT * parent , Cost g , int level ) {
   std::lock_guard< std::mutex > lock( sm );
   if( n->sseen )
    return;
   n->sseen = true;
   n->sparent = parent;
   n->sg = g;
   sopen[ { level , n->h } ].push_back( Entry{ n , 0 } );
   }

  // removes a node of the speculative OPEN, or nullptr
  NodeT * pop_spec( void ) {
   std::lock_guard< std::mutex > lock( sm );
   while( ! sopen.empty() ) {
    auto it = sopen.begin();
    if( it->second.empty() ) {
     sopen.erase( it );
     continue;
     }
    NodeT * n = it->second.front().node;
    it->second.pop_front();
    if( n->sclosed )
     continue;
    n->sclosed = true;
    ++spec_busy;
    return( n );
    }
   return( nullptr );
   }

  };  // end( struct( Shared ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] bool obat( void ) const noexcept {
  return( f_options.rule == GreedyRule::obat );
  }

 [[nodiscard]] bool spuhf( void ) const noexcept {
  return( f_options.rule == GreedyRule::spuhf );
  }

 /// the rules that mark the certain nodes, and keep them apart in OPEN
 [[nodiscard]] bool splits_certain( void ) const noexcept {
  switch( f_options.rule ) {
   case( GreedyRule::kpgbfs ):
   case( GreedyRule::obat ): return( false );
   default: return( true );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: the worker w takes a successor to evaluate, a node of
 /// DEFERRED to commit, or a node of OPEN to expand, if the rule allows
 /// it; leaves w idle otherwise; returns true if it did one of them

 bool acquire( Shared & sh , Worker & w ) {
  if( sh.done )
   return( false );
  if( f_options.sge && ! sh.unevaluated.empty() ) {
   auto [ e , i ] = sh.unevaluated.front();
   sh.unevaluated.pop_front();
   w.exp = e;
   w.kid = i;
   w.evaluating = true;
   return( true );
   }
  const Cost low = sh.least_busy();
  if( obat() && ! sh.deferred.empty() ) {
   const Cost hd = sh.deferred.begin()->first;
   if( sh.open_empty() || ( hd <= sh.top()->h ) ) {
    if( hd <= low ) {
     auto it = sh.deferred.begin();
     auto e = it->second.front().e;
     it->second.pop_front();
     if( it->second.empty() )
      sh.deferred.erase( it );
     list( sh , *e );                       // completely expanded
     return( true );
     }
    return( false );
    }
   }
  if( sh.open_empty() ) {
   if( ( sh.busy_n == 0 ) && sh.deferred.empty() &&
       sh.unevaluated.empty() ) {
    sh.done = true;                           // the main search is over
    sh.over.store( true );
    }
   else
    if( spuhf() )
     return( speculate( sh , w ) );
   return( false );
   }
  NodeT * t = sh.top();
  bool ok = true;
  switch( f_options.rule ) {
   case( GreedyRule::kpgbfs ):
    break;
   case( GreedyRule::puhf ):
   case( GreedyRule::spuhf ):
    if( ( ! t->certain ) && ( sh.busy_n == 0 ) )
     t->certain = true;
    ok = t->certain;
    break;
   case( GreedyRule::puhf2 ):
   case( GreedyRule::puhf3 ):
   case( GreedyRule::puhf4 ):
    if( ( ! t->certain ) && ( t->h <= low ) )
     t->certain = true;
    ok = t->certain;
    break;
   case( GreedyRule::obat ):
    ok = t->h <= low;
    break;
   }
  if( ! ok ) {
   if( spuhf() )
    return( speculate( sh , w ) );
   return( false );
   }
  sh.pop();
  t->removed = true;
  sh.add_busy( t->h );
  if( ( f_options.rule == GreedyRule::puhf3 ) ||
      ( f_options.rule == GreedyRule::puhf4 ) )
   if( t->h == sh.least_busy() )
    sh.certify( t->h );
  ++w.st.expanded;
  if( f_record )
   f_trace.expanded.push_back( t->state );
  w.exp = std::make_shared< Expansion >();
  w.exp->node = t;
  w.exp->busy = true;
  w.evaluating = false;
  if( f_domain.is_goal( t->state ) ) {
   sh.goal = t;
   sh.done = true;
   sh.over.store( true );
   sh.drop_busy( *w.exp );
   w.exp.reset();
   }
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: the worker w, which has no node to expand, takes one of
 /// the speculative OPEN (SPUHF)

 bool speculate( Shared & sh , Worker & w ) {
  NodeT * n = sh.pop_spec();
  if( ! n )
   return( false );
  w.exp = std::make_shared< Expansion >();
  w.exp->node = n;
  w.exp->speculative = true;
  w.evaluating = false;
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// without the lock: the work of w, i.e., the generation of the
 /// successors of its node, evaluated unless with SGE, or the evaluation
 /// of one successor

 void work( Shared & sh , Worker & w ) {
  Expansion & e = *w.exp;
  if( w.evaluating ) {
   evaluate( sh , e.kids[ w.kid ] );
   return;
   }
  NodeT * n = e.node;
  if( e.speculative ) {
   if( f_domain.is_goal( n->state ) )
    return;                                    // found at the commit
   ++w.st.expanded;
   }
  for_each_child( f_domain , n->state , n->h ,
                  e.speculative ? no_op : n->op ,
                  [ & ]( const State & t , Cost c , Cost ht ,
                         std::uint8_t op ) {
   e.kids.push_back( Kid{ t , c , ht , op , nullptr } );
   } );
  w.st.generated += e.kids.size();
  e.left = e.kids.size();
  if( ! ( f_options.sge && ! e.speculative ) )
   for( Kid & k : e.kids )
    evaluate( sh , k );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the node of the successor k, evaluated if new

 void evaluate( Shared & sh , Kid & k ) {
  const std::uint64_t hash = f_domain.hash( k.s );
  NodeT * n = sh.find( k.s , hash );
  if( ! n ) {
   if( f_options.sge )
    k.h = Cost( f_domain.heuristic( k.s ) );
   bool fresh;
   n = sh.lookup( k.s , hash , k.h , fresh );
   if( fresh )
    sh.evaluated.fetch_add( 1 , std::memory_order_relaxed );
   }
  k.node = n;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: the results of the work of w enter the shared lists

 void commit( Shared & sh , Worker & w ) {
  auto e = w.exp;
  w.exp.reset();
  if( e->speculative ) {
   commit_speculative( sh , *e );
   return;
   }
  if( w.evaluating ) {                        // one successor evaluated
   w.evaluating = false;
   if( --e->left == 0 )
    finish( sh , e );
   return;
   }
  if( f_options.sge && ( e->left > 0 ) ) {    // to the evaluators
   for( std::size_t i = 0 ; i < e->kids.size() ; ++i )
    sh.unevaluated.emplace_back( e , i );
   return;
   }
  finish( sh , e );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: the successors of the node of e are all evaluated

 void finish( Shared & sh , const std::shared_ptr< Expansion > & e ) {
  NodeT * n = e->node;
  Cost low = n->h;
  for( const Kid & k : e->kids )
   low = std::min( low , k.node->h );
  if( obat() && ( low < n->h ) ) {            // a b-state: deferred
   sh.deferred[ n->h ].push_back( Deferred{ e , sh.seq++ } );
   sh.drop_busy( *e );
   return;
   }
  list( sh , *e );
  sh.drop_busy( *e );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: the successors of the node of e not in CLOSED enter
 /// CLOSED and OPEN

 void list( Shared & sh , Expansion & e ) {
  NodeT * n = e.node;
  if( f_record )
   f_trace.committed.push_back( n->state );
  Cost low = n->h;
  for( const Kid & k : e.kids )
   low = std::min( low , k.node->h );
  for( const Kid & k : e.kids ) {
   NodeT * m = k.node;
   if( m->listed || ( m->h == infinite_cost< Cost >() ) ) {
    ++sh.duplicates;
    continue;
    }
   m->listed = true;
   m->g = n->g + k.c;
   m->parent = n;
   m->op = k.op;
   if( ( ( f_options.rule == GreedyRule::puhf ) ||
         ( f_options.rule == GreedyRule::spuhf ) ||
         ( f_options.rule == GreedyRule::puhf4 ) ) && ( m->h == low ) )
    m->certain = true;
   sh.push_open( m );
   if( spuhf() && ! m->certain )
    sh.push_spec( m , n , m->g , 0 );
   }
  if( f_options.node_limit && ( sh.stored() >= f_options.node_limit ) ) {
   sh.limit_hit = true;
   sh.done = true;
   sh.over.store( true );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// under the lock: a speculative expansion ends (SPUHF)

 void commit_speculative( Shared & sh , Expansion & e ) {
  NodeT * n = e.node;
  {
   std::lock_guard< std::mutex > lock( sh.sm );
   --sh.spec_busy;
   ++sh.speculative;
   }
  if( f_domain.is_goal( n->state ) ) {
   sh.goal = n;
   sh.goal_spec = true;
   sh.done = true;
   sh.over.store( true );
   return;
   }
  const Cost g = n->listed ? n->g : n->sg;
  for( const Kid & k : e.kids )
   if( k.node->h != infinite_cost< Cost >() )
    sh.push_spec( k.node , n , g + k.c , 1 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the loop of a worker on the executor

 void run( Shared & sh , Worker & w ) {
  ADS::Backoff backoff;
  while( ! sh.over.load( std::memory_order_acquire ) ) {
   {
    std::lock_guard< std::mutex > lock( sh.m );
    acquire( sh , w );
    }
   if( w.idle() ) {
    backoff.pause();
    continue;
    }
   backoff.reset();
   work( sh , w );
   std::lock_guard< std::mutex > lock( sh.m );
   commit( sh , w );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the path to the goal n, through the speculative parents if spec

 std::vector< State > path( NodeT * n , bool spec ) const {
  std::vector< State > v;
  for( const NodeT * x = n ; x ; ) {
   v.push_back( x->state );
   if( spec && ! x->listed )
    x = x->sparent;
   else
    x = x->parent;
   }
  std::reverse( v.begin() , v.end() );
  return( v );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the cost of a path, the cheapest edge between two consecutive states

 Cost path_cost( const std::vector< State > & v ) const {
  Cost sum( 0 );
  for( std::size_t i = 1 ; i < v.size() ; ++i ) {
   Cost best = infinite_cost< Cost >();
   f_domain.for_each_successor( v[ i - 1 ] , Cost( 0 ) ,
                                [ & ]( const State & t , Cost c , Cost ) {
    if( ( t == v[ i ] ) && ( c < best ) )
     best = c;
    } );
   sum = sum + best;
   }
  return( sum );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 ParallelGreedyOptions f_options;                    ///< the options
 Executor * f_ex = nullptr;                          ///< given, or nullptr
 bool f_record = false;                              ///< record the states
 GreedyTrace< State > f_trace;                       ///< the last trace

 };  // end( class( ParallelGreedy ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS HDGBFS ------------------------------*/
/*--------------------------------------------------------------------------*/
/// Hash Distributed GBFS on a shared-memory machine
/** HDAStar with the order of GBFS: every worker runs GBFS on its own OPEN
 * and CLOSED lists, ordered by h and then by the order of arrival (or the
 * opposite, with TieBreak::lifo), the states being assigned to the workers
 * by the distribution Dist, and the search ends when a worker selects a
 * goal [see Kuroiwa and Fukunaga, ICAPS 2019]. Hash distribution keeps
 * the duplicate detection exact, but not the global order: a worker
 * expands the best node it owns, however bad, and HDGBFS is pathological
 * relative to GBFS in the sense of Kuroiwa and Fukunaga (ICAPS 2020). A
 * state reached again through a cheaper path is not reopened. */

template< Domain D , class Dist = StateHashDistribution< D > ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
class HDGBFS {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = Node< State , Cost >;               ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

 /// a node travelling to its owner
 struct Msg {
  State s;
  Cost g;
  Cost h;
  NodeT * parent;
  [[no_unique_address]] OperatorOf< D > o;        // the operator of s
  std::uint32_t key;                               // the key of s
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain and the distribution (which must
 /// outlive the search)

 HDGBFS( const D & domain , const Dist & dist ,
         const ParallelGreedyOptions & options = {} )
  : f_domain( domain ) , f_dist( dist ) , f_options( options ) {
  f_options.normalize();
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the Executor that runs the workers [see HDAStar::set_executor()]

 void set_executor( Executor & ex ) noexcept { f_ex = &ex; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const
 {
  const unsigned p = f_options.threads;
  Stopwatch clock;
  Shared sh( p , f_options.batch );
  const State s0 = f_domain.initial_state();
  sh.root = Msg{ s0 , Cost( 0 ) , Cost( f_domain.heuristic( s0 ) ) ,
                 nullptr , no_op , distribution_key< D >( f_dist , s0 ) };
  sh.root_owner = f_dist.owner( s0 , p );

  std::vector< std::unique_ptr< Worker > > w( p );
  for( unsigned i = 0 ; i < p ; ++i )
   w[ i ] = std::make_unique< Worker >( *this , sh , i );
  run_workers( f_ex , p , [ & ]( unsigned i ) { w[ i ]->run(); } );

  SolutionT sol;
  SearchStats & st = sol.stats;
  for( auto & x : w ) {
   st.expanded += x->f_st.expanded;
   st.generated += x->f_st.generated;
   st.duplicates += x->f_st.duplicates;
   st.sent += x->f_st.sent;
   st.messages += sh.mail.batches( x->f_id );
   st.stored += x->f_pool.size();
   st.expanded_by.push_back( x->f_st.expanded );
   st.stored_by.push_back( x->f_pool.size() );
   }
  st.evaluated = st.stored;
  if( sh.goal ) {
   sol.status = Status::solved;
   sol.cost = sh.goal->g;
   sol.path = trace_path( sh.goal );
   }
  else
   if( sh.limit_hit.load() )
    sol.status = Status::node_limit;
  st.seconds = clock.seconds();
  return( sol );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 using Mailbox = ADS::Mailbox< Msg , ADS::DefaultChannels< Msg , Executor > >;

 // the data shared by the workers
 struct Shared {
  Shared( unsigned np , std::size_t batch ) : p( np ) , mail( np , batch )
  {}
  unsigned p;
  Mailbox mail;
  Msg root;
  unsigned root_owner = 0;
  std::atomic< bool > limit_hit{ false };
  std::mutex goal_mutex;
  NodeT * goal = nullptr;
  };

 // an entry of OPEN
 struct Item {
  Cost h;
  std::uint64_t seq;
  NodeT * node;
  };

 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   return( ( a.h < b.h ) || ( ( a.h == b.h ) && ( a.seq < b.seq ) ) );
   }
  };

 // a worker
 struct alignas( cache_line ) Worker {

  Worker( const HDGBFS & a , Shared & s , unsigned i )
   : f_a( a ) , f_sh( s ) , f_id( i ) , f_dist( a.f_dist ) {}

  void run( void ) {
   if( f_id == f_sh.root_owner )
    insert( f_sh.root , f_a.f_domain.hash( f_sh.root.s ) );
   unsigned since_flush = 0;
   ADS::Backoff backoff;
   while( ! f_sh.mail.stopped() ) {
    f_sh.mail.poll( f_id , [ & ]( Msg && m ) {
     insert( m , f_a.f_domain.hash( m.s ) );
     } );
    if( ! f_open.empty() ) {
     backoff.reset();
     expand( f_open.pop().node );
     if( ++since_flush >= f_a.f_options.flush_period ) {
      since_flush = 0;
      f_sh.mail.flush( f_id );
      }
     if( f_a.f_options.node_limit &&
         ( f_pool.size() >= f_a.f_options.node_limit ) ) {
      f_sh.limit_hit.store( true );
      f_sh.mail.stop();
      }
     continue;
     }
    if( f_sh.mail.quiescent( f_id ) )
     return;
    backoff.pause();
    }
   }

  void expand( NodeT * n ) {
   n->closed = true;
   n->expanded = true;
   ++f_st.expanded;
   if( f_a.f_domain.is_goal( n->state ) ) {
    std::lock_guard< std::mutex > lock( f_sh.goal_mutex );
    if( ! f_sh.goal )
     f_sh.goal = n;
    f_sh.mail.stop();
    return;
    }
   const auto owner = successor_owners< D >( f_dist , n->state , n->key ,
                                             f_sh.p );
   for_each_child( f_a.f_domain , n->state , n->h , n->op ,
                   [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t ot ) {
    ++f_st.generated;
    if( ht == infinite_cost< Cost >() )
     return;
    const auto [ dst , key ] = owner( t );
    const Msg m{ t , n->g + c , ht , n , ot , key };
    if( dst == f_id ) {
     insert( m , f_a.f_domain.hash( t ) );
     return;
     }
    f_sh.mail.send( f_id , dst , m );
    ++f_st.sent;
    } );
   }

  void insert( const Msg & m , std::uint64_t h ) {
   auto [ x , fresh ] = f_table.find_or_insert( m.s , h , [ & ] {
    NodeT * y = f_pool.allocate();
    *y = NodeT{ m.s , m.g , m.h , m.parent , false , false , m.o.op ,
                m.key };
    return( y );
    } );
   if( fresh ) {
    const bool lifo = f_a.f_options.tie == TieBreak::lifo;
    f_open.push( Item{ m.h , lifo ? ~f_seq : f_seq , x } );
    ++f_seq;
    }
   else
    ++f_st.duplicates;
   }

  const HDGBFS & f_a;
  Shared & f_sh;
  unsigned f_id;
  Dist f_dist;
  NodePool< NodeT > f_pool;
  NodeTable< NodeT > f_table;
  ADS::DaryHeap< Item , Better > f_open;
  std::uint64_t f_seq = 0;
  SearchStats f_st;
  };

 const D & f_domain;                                   ///< the Domain
 const Dist & f_dist;                                  ///< the distribution
 ParallelGreedyOptions f_options;                      ///< the options
 Executor * f_ex = nullptr;                            ///< given, or nullptr

 };  // end( class( HDGBFS ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_ParallelGreedy */

/*--------------------------------------------------------------------------*/
/*----------------------- End File ParallelGreedy.h ------------------------*/
/*--------------------------------------------------------------------------*/
