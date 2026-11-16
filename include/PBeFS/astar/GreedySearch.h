/*--------------------------------------------------------------------------*/
/*-------------------------- File GreedySearch.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class template GreedySearch, the greedy best-first
 * search (GBFS) of Doran and Michie (1966) with the OPEN lists of the
 * satisficing planners as options that compose, for the class template
 * GreedyOpen, those lists, and for the class template HeuristicCache, the
 * values of the heuristics kept across several searches.
 *
 * GBFS selects the node of least \f$ h \f$, and is the limit
 * \f$ w \to \infty \f$ of the search ordered by \f$ g + w h \f$, which the
 * same class performs with a finite weight w. The options of GreedyOptions
 * are those of the planners of the last two decades, implemented from the
 * papers:
 *
 * - deferred (lazy) evaluation, in which a successor enters OPEN with the
 *   heuristic value of its parent and is generated and evaluated only when
 *   it is selected [see Richter and Helmert, ICAPS 2009];
 *
 * - several queues alternated in round robin, each ordered by one of the
 *   heuristics of a MultiHeuristicDomain, its ties broken by another
 *   [see Röger and Helmert, ICAPS 2010];
 *
 * - queues of the preferred successors only [see PreferredDomain], which
 *   are boosted, i.e., selected the next boost times in a row, whenever a
 *   heuristic reaches a value lower than all those seen so far [see
 *   Richter and Helmert, ICAPS 2009];
 *
 * - type-based exploration, a queue that draws a type \f$ ( h , g ) \f$ at
 *   random among those in OPEN, and then a node of that type at random
 *   [see Xie et al., AAAI 2014];
 *
 * - epsilon-greedy selection, which takes a node of OPEN at random with
 *   probability \f$ \epsilon \f$ and the best one otherwise [see Valenzano
 *   et al., ICAPS 2014];
 *
 * - Pareto queues, which draw a node at random among the vectors of the
 *   heuristic values that no other vector in OPEN dominates [see Röger and
 *   Helmert, ICAPS 2010];
 *
 * - the breaking of the ties of the key by the order of insertion (FIFO or
 *   LIFO), by the smaller or the larger g, or at random.
 *
 * The search with eager evaluation is the following, with the notation of
 * the page "Best-first search and A*" of the documentation and
 * \f$ Q_1 , \ldots , Q_q \f$ the queues of OPEN, each with a priority
 * \f$ \pi_j \f$:
 *
 * \code{.unparsed}
 *  1  procedure GBFS(s₀)
 *  2     CLOSED ← {s₀}, with g(s₀) ← 0, parent(s₀) ← nil, and the values
 *          h_i(s₀) of the heuristics;  best_i ← h_i(s₀) for each i
 *  3     push s₀ in each Q_j that is not for the preferred successors only
 *  4     loop
 *  5        j ← a nonempty queue of least π_j;  π_j ← π_j + 1
 *  6        if there is none then return unsolvable
 *  7        u ← the first node of Q_j, removed
 *  8        if u is closed or its g changed since it was pushed then
 *              continue
 *  9        close u;  if u is a goal then return the path to u
 * 10        P ← the preferred successors of u
 * 11        for each successor v of u do
 * 12           g' ← g(u) + c(u, v);  if g' ≥ bound then continue
 * 13           if v ∉ CLOSED then
 * 14              add v to CLOSED, with g(v) ← g', parent(v) ← u and the
 *                   values h_i(v);  if h(v) = ∞ then continue  ▷ dead end
 * 15              if h_i(v) < best_i for some i then
 * 16                 best_i ← h_i(v);  π_j ← π_j − boost for each queue of
 *                      preferred successors
 * 17              push v in each Q_j, those of the preferred successors
 *                   only if v ∈ P
 * 18           else if g' < g(v) and (v is not closed or reopening) then
 * 19              g(v) ← g';  parent(v) ← u;  v is not closed;  push v
 *                   again as at line 17
 * 20     ▷ with lazy evaluation, the entries are the pairs (u, v) pushed
 *          with the keys of u, and v is looked up, evaluated and expanded
 *          when its entry is selected, if it is not in CLOSED
 * \endcode
 *
 * The first goal selected is returned, so that the path is not optimal in
 * general; with a finite weight w and an admissible h, reopening the
 * closed nodes, its cost is at most \f$ w C^* \f$ [see Pohl, AIJ 1970].
 * The search is complete on a finite graph, since every state enters
 * CLOSED once and is expanded at most once without reopening (a node is
 * expanded again only by line 18 with reopening, a finite number of times
 * since \f$ g \f$ decreases each time and the paths without cycles are
 * finitely many).
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

#ifndef __PBeFS_GreedySearch
 #define __PBeFS_GreedySearch
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <cmath>

#include <cstddef>

#include <cstdint>

#include <deque>

#include <limits>

#include <map>

#include <random>

#include <stdexcept>

#include <unordered_map>

#include <vector>

#include "PBeFS/astar/Domain.h"

#include "PBeFS/astar/Node.h"

#include "PBeFS/astar/NodeTable.h"

#include "PBeFS/astar/Solution.h"

#include "PBeFS/core/SearchUtils.h"

#include "ADS/priority/DaryHeap.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- NAMESPACE --------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Parallel Best-First Search library
namespace PBeFS_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------------- TYPES ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// when the heuristic of a successor is computed

//! [best-first-sig-greedy-kinds]
enum class GreedyEvaluation {
 eager ,     ///< when it is generated, the key being its own h
 lazy        ///< when it is selected, the key being the h of its parent
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// how the ties of the key of a queue are broken

enum class TieBreak {
 fifo ,      ///< the first pushed
 lifo ,      ///< the last pushed
 low_g ,     ///< the smaller g, then the first pushed
 high_g ,    ///< the larger g, then the first pushed
 high_g_lifo ,  ///< the larger g, then the last pushed
 random      ///< at random
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the kind of a queue of OPEN

enum class QueueKind {
 best ,      ///< the least key, the ties broken by the next heuristic of
             ///< the queue and then by the TieBreak
 type ,      ///< a type at random, then a node of that type at random
 epsilon ,   ///< a node at random with probability epsilon, else the best
 pareto      ///< a vector at random among those not dominated
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a queue of OPEN: its kind, the heuristics it looks at, and whether it
/// holds only the preferred successors
/** A queue of kind best orders by the first heuristic in heuristics, the
 * ties broken by the second, if any; one of kind type takes the type of a
 * node from the first heuristic (and from g, if
 * GreedyOptions::type_with_g); one of kind epsilon takes the best node by
 * the first heuristic; one of kind pareto looks at the vector of all of
 * them. */

struct QueueSpec {

 QueueKind kind = QueueKind::best;           ///< the kind of the queue
 std::vector< std::size_t > heuristics = { 0 };  ///< the heuristics [see
                                             ///< MultiHeuristicDomain]
 bool preferred_only = false;                ///< only the preferred
                                             ///< successors [see
                                             ///< PreferredDomain]

 };  // end( struct( QueueSpec ) )
 //! [best-first-sig-greedy-kinds]

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS GreedyOptions ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the options of GreedySearch, besides those of every search; with eager
/// evaluation, SearchOptions::goal_at_generation returns a goal when it is
/// generated rather than when it is selected

struct GreedyOptions : public SearchOptions {

 /// constructor: the weight w of SearchOptions is infinite, which is GBFS;
 /// the key of a node is \f$ g + w h \f$ for a finite \f$ w \geq 0 \f$,
 /// and h alone for an infinite one

 GreedyOptions( void ) {
  weight = std::numeric_limits< double >::infinity();
  }

 /// eager or lazy evaluation
 GreedyEvaluation evaluation = GreedyEvaluation::eager;

 /// the queues of OPEN, alternated; one best queue on h by default
 std::vector< QueueSpec > queues = { QueueSpec{} };

 /// the breaking of the ties of the keys
 TieBreak tie = TieBreak::fifo;

 /// what the priorities of the queues of the preferred successors
 /// decrease by when a heuristic reaches a new least value (0: no boost)
 int boost = 0;

 /// the probability of a random node in a queue of kind epsilon
 double epsilon = 0.2;

 /// the type of a node in a queue of kind type is \f$ ( h , g ) \f$ if
 /// true, h alone otherwise
 bool type_with_g = true;

 /// a closed node reached through a cheaper path is reopened
 bool reopen = false;

 /// a node in OPEN reached through a cheaper path takes it, and is pushed
 /// again with its new g (and the place of a new entry among its ties),
 /// as in Fast Downward; if false it keeps the first path found, as in
 /// the parallel searches of ParallelGreedy
 bool improve_open = true;

 /// the successors with \f$ g \geq \mathrm{bound} \f$ are pruned, or
 /// those with \f$ g + h \geq \mathrm{bound} \f$ if prune_by_f
 double bound = std::numeric_limits< double >::infinity();

 /// the bound prunes by \f$ g + h \f$ instead of g (for an admissible h)
 bool prune_by_f = false;

 /// with lazy evaluation, the preferred successors enter OPEN before the
 /// others, which matters for the ties broken by the order
 bool preferred_first = false;

 /// the successors enter OPEN in a random order
 bool randomize_successors = false;

 /// the seed of the random choices
 std::uint64_t seed = 1;

 };  // end( struct( GreedyOptions ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS GreedyNode ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a node of the greedy searches
/** As Node, with the index of the node, by which the values of the
 * heuristics beyond the first are found [see MultiHeuristicDomain]; a dead
 * end has \f$ h = \infty \f$ and stays closed. */

template< class State , class Cost >
struct GreedyNode {

 State state;               ///< the state
 Cost g;                    ///< the cost of the best path found
 Cost h;                    ///< the first heuristic
 GreedyNode * parent;       ///< the predecessor on that path
 bool closed;               ///< true if expanded with the current g
 bool expanded;             ///< true if ever expanded
 std::uint8_t op = no_op;   ///< the operator that generated the state
 std::uint32_t index = 0;   ///< the index of the node

 };  // end( struct( GreedyNode ) )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS GreedyOpen ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the OPEN list of the greedy searches: queues alternated by priority
/** Each entry, of type P, is pushed with the values of the heuristics that
 * make its keys, its g and whether it is preferred, into every queue that
 * takes it [see QueueSpec]; pop() selects the nonempty queue of least
 * priority, the first one among the ties, adds 1 to its priority and
 * removes its first entry, and boost() lowers the priorities of the queues
 * of the preferred successors, as in the alternation of Röger and Helmert
 * (ICAPS 2010). The same entry may lie in several queues: the search
 * discards the stale ones as they surface. With the costs of a Domain, the
 * keys are kept as double, exact for integers below \f$ 2^{53} \f$. */

template< class P >
class GreedyOpen {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the options; throws std::invalid_argument if they
 /// name no queue, a heuristic out of [ 0 , m ), or a queue of kind
 /// pareto with no heuristic

 GreedyOpen( const GreedyOptions & o , std::size_t m )
  : f_o( o ) , f_rng( o.seed ) {
  if( o.queues.empty() )
   throw( std::invalid_argument( "GreedyOpen: no queue" ) );
  for( const QueueSpec & q : o.queues ) {
   if( q.heuristics.empty() )
    throw( std::invalid_argument( "GreedyOpen: a queue with no heuristic" ) );
   for( std::size_t i : q.heuristics )
    if( i >= m )
     throw( std::invalid_argument( "GreedyOpen: heuristic out of range" ) );
   v_q.emplace_back( q );
   }
  }

/*--------------------------------- QUEUES ---------------------------------*/
 /// pushes p, whose heuristics are hv[ 0 ] , ... , hv[ m - 1 ], with g
 /// and preferred or not

 template< class Cost >
 void push( const P & p , const Cost * hv , Cost g , bool preferred ) {
  const std::uint64_t seq = f_seq++;
  for( Queue & q : v_q ) {
   if( q.spec.preferred_only && ! preferred )
    continue;
   q.push( *this , p , hv , double( g ) , seq );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// removes the first entry of the nonempty queue of least priority into
 /// p; false if every queue is empty

 bool pop( P & p ) {
  Queue * best = nullptr;
  for( Queue & q : v_q )
   if( ! q.empty() && ( ! best || ( q.priority < best->priority ) ) )
    best = &q;
  if( ! best )
   return( false );
  ++best->priority;
  p = best->pop( *this );
  return( true );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// lowers by boost the priorities of the queues of the preferred
 /// successors

 void boost( void ) {
  if( f_o.boost == 0 )
   return;
  for( Queue & q : v_q )
   if( q.spec.preferred_only )
    q.priority -= f_o.boost;
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if every queue is empty

 [[nodiscard]] bool empty( void ) const {
  return( std::all_of( v_q.begin() , v_q.end() ,
                       []( const Queue & q ) { return( q.empty() ); } ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the entries in all the queues, the stale ones included

 [[nodiscard]] std::size_t size( void ) const {
  std::size_t n = 0;
  for( const Queue & q : v_q )
   n += q.size();
  return( n );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 // an entry of a heap: the keys, the order of the ties, the unique id of
 // the insertion (for the epsilon queues) and the payload
 struct Item {
  double k1;
  double k2;
  double k3;
  std::uint64_t tie;
  std::uint64_t uid;
  P p;
  };

 // the lexicographic order of the keys
 struct Better {
  bool operator()( const Item & a , const Item & b ) const noexcept {
   if( a.k1 != b.k1 )
    return( a.k1 < b.k1 );
   if( a.k2 != b.k2 )
    return( a.k2 < b.k2 );
   if( a.k3 != b.k3 )
    return( a.k3 < b.k3 );
   return( a.tie < b.tie );
   }
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the key of a node from its heuristic h and its g

 [[nodiscard]] double key( double h , double g ) const noexcept {
  if( std::isinf( f_o.weight ) )
   return( h );
  return( g + f_o.weight * h );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // the third key and the order of the ties, by the TieBreak

 void ties( double g , std::uint64_t seq , double & k3 ,
            std::uint64_t & tie ) {
  k3 = 0;
  switch( f_o.tie ) {
   case( TieBreak::fifo ): tie = seq; break;
   case( TieBreak::lifo ): tie = ~seq; break;
   case( TieBreak::low_g ): k3 = g; tie = seq; break;
   case( TieBreak::high_g ): k3 = -g; tie = seq; break;
   case( TieBreak::high_g_lifo ): k3 = -g; tie = ~seq; break;
   case( TieBreak::random ): tie = std::uint64_t( f_rng() ); break;
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // a queue

 struct Queue {

  explicit Queue( const QueueSpec & s ) : spec( s ) {}

  QueueSpec spec;                         // what it is
  long priority = 0;                      // its priority
  ADS::DaryHeap< Item , Better > heap;    // best, epsilon
  std::vector< Item > pool;               // epsilon: the entries in it
  std::vector< char > taken;              // epsilon: by uid, if removed
  std::size_t live = 0;                   // epsilon: entries not taken
  std::unordered_map< std::uint64_t , std::size_t > type_of;  // type ->
                                          // its bucket
  std::vector< std::vector< P > > buckets;  // type: the nodes of a type
  std::vector< std::size_t > nonempty;    // type: the nonempty buckets
  std::vector< std::size_t > where;       // type: place in nonempty
  std::map< std::vector< double > , std::deque< P > > front;  // pareto
  std::vector< const std::vector< double > * > undominated;   // pareto
  bool dirty = false;                     // pareto: undominated is old
  std::size_t count = 0;                  // the entries

  [[nodiscard]] bool empty( void ) const noexcept { return( count == 0 ); }
  [[nodiscard]] std::size_t size( void ) const noexcept { return( count ); }

  // pushes p
  template< class Cost >
  void push( GreedyOpen & o , const P & p , const Cost * hv , double g ,
             std::uint64_t seq ) {
   ++count;
   const auto & hs = spec.heuristics;
   switch( spec.kind ) {
    case( QueueKind::best ):
    case( QueueKind::epsilon ): {
     Item x;
     x.k1 = o.key( double( hv[ hs[ 0 ] ] ) , g );
     x.k2 = hs.size() > 1 ? double( hv[ hs[ 1 ] ] ) : 0.0;
     o.ties( g , seq , x.k3 , x.tie );
     x.uid = taken.size();
     x.p = p;
     heap.push( x );
     if( spec.kind == QueueKind::epsilon ) {
      taken.push_back( 0 );
      pool.push_back( x );
      ++live;
      }
     break;
     }
    case( QueueKind::type ): {
     const double h = double( hv[ hs[ 0 ] ] );
     const double gt = o.f_o.type_with_g ? g : 0.0;
     const std::uint64_t t = ADS::splitmix64(
      std::bit_cast< std::uint64_t >( h ) ^
      ADS::splitmix64( std::bit_cast< std::uint64_t >( gt ) ) );
     auto [ it , fresh ] = type_of.emplace( t , buckets.size() );
     if( fresh ) {
      buckets.emplace_back();
      where.push_back( 0 );
      }
     const std::size_t b = it->second;
     if( buckets[ b ].empty() ) {
      where[ b ] = nonempty.size();
      nonempty.push_back( b );
      }
     buckets[ b ].push_back( p );
     break;
     }
    case( QueueKind::pareto ): {
     std::vector< double > v( hs.size() );
     for( std::size_t i = 0 ; i < hs.size() ; ++i )
      v[ i ] = double( hv[ hs[ i ] ] );
     auto & d = front[ v ];
     if( d.empty() )
      dirty = true;
     d.push_back( p );
     break;
     }
    }
   }

  // removes the first entry
  P pop( GreedyOpen & o ) {
   --count;
   switch( spec.kind ) {
    case( QueueKind::best ):
     return( heap.pop().p );
    case( QueueKind::epsilon ):
     if( random_unit( o.f_rng ) < o.f_o.epsilon ) {
      for( ; ; ) {                     // a live entry of pool at random
       const std::size_t i = std::size_t( random_below( o.f_rng ,
                                                        pool.size() ) );
       const Item x = pool[ i ];
       pool[ i ] = pool.back();
       pool.pop_back();
       if( ! taken[ x.uid ] ) {
        taken[ x.uid ] = 1;
        --live;
        compact();
        return( x.p );
        }
       }
      }
     for( ; ; ) {                      // the best live entry of the heap
      const Item x = heap.pop();
      if( ! taken[ x.uid ] ) {
       taken[ x.uid ] = 1;
       --live;
       compact();
       return( x.p );
       }
      }
    case( QueueKind::type ): {
     const std::size_t k = std::size_t( random_below( o.f_rng ,
                                                      nonempty.size() ) );
     const std::size_t b = nonempty[ k ];
     auto & v = buckets[ b ];
     const std::size_t i = std::size_t( random_below( o.f_rng , v.size() ) );
     const P p = v[ i ];
     v[ i ] = v.back();
     v.pop_back();
     if( v.empty() ) {
      nonempty[ k ] = nonempty.back();
      where[ nonempty[ k ] ] = k;
      nonempty.pop_back();
      }
     return( p );
     }
    case( QueueKind::pareto ): {
     if( dirty )
      refresh();
     const std::size_t k = std::size_t( random_below( o.f_rng ,
                                                      undominated.size() ) );
     auto it = front.find( *undominated[ k ] );
     P p;
     if( o.f_o.tie == TieBreak::lifo ) {
      p = it->second.back();
      it->second.pop_back();
      }
     else {
      p = it->second.front();
      it->second.pop_front();
      }
     if( it->second.empty() ) {
      front.erase( it );
      dirty = true;
      }
     return( p );
     }
    }
   return( P() );
   }

  // drops from pool and heap the entries taken, when they are the most
  void compact( void ) {
   if( pool.size() < 2 * live + 64 )
    return;
   std::erase_if( pool , [ & ]( const Item & x ) { return( taken[ x.uid ] );
                                                   } );
   }

  // the vectors of the front that no other one dominates, i.e., none is
  // at most as large in every component and smaller in one
  void refresh( void ) {
   undominated.clear();
   for( const auto & [ v , d ] : front ) {
    bool dominated = false;
    for( const auto & [ u , e ] : front ) {
     if( &u == &v )
      continue;
     bool le = true , lt = false;
     for( std::size_t i = 0 ; le && ( i < v.size() ) ; ++i ) {
      le = u[ i ] <= v[ i ];
      lt = lt || ( u[ i ] < v[ i ] );
      }
     if( le && lt ) {
      dominated = true;
      break;
      }
     }
    if( ! dominated )
     undominated.push_back( &v );
    }
   dirty = false;
   }

  };  // end( struct( Queue ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const GreedyOptions & f_o;                   ///< the options
 std::mt19937_64 f_rng;                       ///< the random choices
 std::uint64_t f_seq = 0;                     ///< the insertions so far
 std::vector< Queue > v_q;                    ///< the queues

 };  // end( class( GreedyOpen ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS HeuristicCache --------------------------*/
/*--------------------------------------------------------------------------*/
/// the values of the heuristics of the states evaluated by a search, kept
/// for the next ones
/** A restarting search, such as the iterated searches of the planners
 * [see IteratedSearch], evaluates again many of the states of the search
 * before; the cache, given to each search by set_cache(), keeps the m
 * values of each state evaluated, so that a state is evaluated once over
 * all the searches [see Richter, Thayer and Ruml, ICAPS 2010]. */

template< Domain D >
class HeuristicCache {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain, whose heuristics it keeps

 explicit HeuristicCache( const D & domain )
  : f_m( heuristics_of( domain ) ) {}

/*-------------------------------- LOOKUPS ---------------------------------*/
 /// the m values of the heuristics of s, of hash hash, in a slot of the
 /// cache; true in fresh if s was not there, the slot then to be filled

 Cost * find_or_insert( const State & s , std::uint64_t hash , bool & fresh )
 {
  auto [ e , f ] = f_table.find_or_insert( s , hash , [ & ] {
   Entry * x = f_pool.allocate();
   *x = Entry{ s , std::uint32_t( f_values.size() / f_m ) };
   f_values.resize( f_values.size() + f_m );
   return( x );
   } );
  fresh = f;
  if( f )
   ++f_misses;
  else
   ++f_hits;
  return( f_values.data() + std::size_t( e->index ) * f_m );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the lookups that found the state, and those that did not

 [[nodiscard]] std::uint64_t hits( void ) const noexcept { return( f_hits ); }

 [[nodiscard]] std::uint64_t misses( void ) const noexcept {
  return( f_misses );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// a state and the index of its values
 struct Entry {
  State state;
  std::uint32_t index;
  };

 std::size_t f_m;                        ///< the number of heuristics
 NodePool< Entry > f_pool;               ///< the entries
 NodeTable< Entry > f_table;             ///< the entries by state
 std::vector< Cost > f_values;           ///< the values, m per entry
 std::uint64_t f_hits = 0;               ///< the lookups that found it
 std::uint64_t f_misses = 0;             ///< those that did not

 };  // end( class( HeuristicCache ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS GreedySearch ---------------------------*/
/*--------------------------------------------------------------------------*/
/// greedy best-first search, and the search ordered by \f$ g + w h \f$,
/// with the OPEN lists of the satisficing planners
/** The search of the file comment, with the options of GreedyOptions. The
 * heuristics beyond the first are evaluated with the first, by
 * heuristic_of(); the first one comes from for_each_successor(), or, with
 * lazy evaluation on a DeferredDomain, from heuristic() when the state is
 * selected. A heuristic value equal to infinite_cost() marks a dead end,
 * which is never expanded. The counters are those of SearchStats, and
 * evaluated is the number of states whose heuristics were computed (the
 * cached ones excluded [see HeuristicCache]); without a DeferredDomain,
 * lazy evaluation counts as evaluated only the states it selects, although
 * the Domain computes h for all the successors it generates. */

template< Domain D >
class GreedySearch {

/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using State = state_t< D >;                       ///< the state type
 using Cost = cost_t< D >;                         ///< the cost type
 using NodeT = GreedyNode< State , Cost >;         ///< the node type
 using SolutionT = Solution< State , Cost >;       ///< the solution type

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor, taking the Domain (which must outlive the search); throws
 /// std::invalid_argument if the options are not valid [see GreedyOpen]

 explicit GreedySearch( const D & domain , const GreedyOptions & options = {} )
  : f_domain( domain ) , f_options( options ) ,
    f_m( heuristics_of( domain ) ) {
  GreedyOpen< NodeT * > check( f_options , f_m );
  if( ! ( f_options.weight >= 0 ) )
   throw( std::invalid_argument( "GreedySearch: a negative weight" ) );
  f_preferred = f_options.preferred_first;
  for( const QueueSpec & q : f_options.queues )
   f_preferred = f_preferred || q.preferred_only;
  }

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the cache of the heuristic values, which must outlive the
 /// searches (nullptr: none)

 void set_cache( HeuristicCache< D > * cache ) noexcept { f_cache = cache; }

/*-------------------------------- SOLVING ---------------------------------*/
 /// searches the Domain, returning the Solution

 [[nodiscard]] SolutionT solve( void ) const {
  if( f_options.evaluation == GreedyEvaluation::lazy )
   return( solve_lazy() );
  return( solve_eager() );
  }

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

 /// an entry of the eager search: the node and its g when pushed
 struct Eager {
  NodeT * node = nullptr;
  Cost g{};
  };

 /// the record of a successor of the lazy search, not yet evaluated
 struct Lazy {
  NodeT * parent;
  State state;
  Cost g;
  Cost h;              // h( state ) if known
  std::uint8_t op;
  bool known;          // true if h is known
  };

 /// a successor generated, before it enters OPEN
 struct Kid {
  State s;
  Cost c;
  Cost h;
  std::uint8_t op;
  bool preferred;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the work space of a search

 struct Run {

  explicit Run( const GreedySearch & a )
   : open_e( a.f_options , a.f_m ) , open_l( a.f_options , a.f_m ) ,
     table( a.f_options.expected_nodes
            ? NodeTable< NodeT >::slots_for( a.f_options.expected_nodes )
            : 1024 ) , best( a.f_m , infinite_cost< Cost >() ) ,
     hv( a.f_m ) , rng( a.f_options.seed ^ 0x5DEECE66DULL ) {}

  GreedyOpen< Eager > open_e;              // OPEN of the eager search
  GreedyOpen< std::size_t > open_l;        // OPEN of the lazy search
  NodePool< NodeT > pool;                  // the nodes
  NodeTable< NodeT > table;                // CLOSED
  std::vector< Cost > extra;               // the heuristics beyond h
  std::vector< Lazy > records;             // the lazy records
  std::vector< Kid > kids;                 // the successors of a node
  std::vector< Cost > best;                // the least value of each h
  std::vector< Cost > hv;                  // the values of a node
  std::mt19937_64 rng;                     // the order of the successors
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// a new node of s, of hash hash, with g, h, parent and op

 NodeT * make( Run & r , const State & s , std::uint64_t hash , Cost g ,
               Cost h , NodeT * parent , std::uint8_t op ) const {
  return( r.table.find_or_insert( s , hash , [ & ] {
   NodeT * x = r.pool.allocate();
   *x = NodeT{ s , g , h , parent , false , false , op ,
               std::uint32_t( r.pool.size() - 1 ) };
   return( x );
   } ).first );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the values of the heuristics of the node n into r.hv, h being the
 /// first: they are computed, or taken from the cache, and stored with n;
 /// returns false if n is a dead end

 bool evaluate( Run & r , NodeT * n , Cost h , std::uint64_t hash ,
                SearchStats & st ) const {
  Cost * c = nullptr;
  bool fresh = true;
  if( f_cache )
   c = f_cache->find_or_insert( n->state , hash , fresh );
  if( fresh ) {
   r.hv[ 0 ] = h;
   for( std::size_t i = 1 ; i < f_m ; ++i )
    r.hv[ i ] = heuristic_of( f_domain , n->state , i );
   ++st.evaluated;
   if( c )
    std::copy( r.hv.begin() , r.hv.end() , c );
   }
  else
   std::copy( c , c + f_m , r.hv.begin() );
  n->h = r.hv[ 0 ];
  if( f_m > 1 ) {
   r.extra.resize( std::max( r.extra.size() ,
                             ( std::size_t( n->index ) + 1 ) * ( f_m - 1 ) ) );
   std::copy( r.hv.begin() + 1 , r.hv.end() ,
              r.extra.begin() + std::size_t( n->index ) * ( f_m - 1 ) );
   }
  return( r.hv[ 0 ] != infinite_cost< Cost >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the values of the heuristics of the evaluated node n into r.hv

 void values( Run & r , const NodeT * n ) const {
  r.hv[ 0 ] = n->h;
  for( std::size_t i = 1 ; i < f_m ; ++i )
   r.hv[ i ] = r.extra[ std::size_t( n->index ) * ( f_m - 1 ) + i - 1 ];
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if some value of r.hv is below the least seen so far, which then
 /// becomes the least

 static bool progress( Run & r ) {
  bool p = false;
  for( std::size_t i = 0 ; i < r.hv.size() ; ++i )
   if( r.hv[ i ] < r.best[ i ] ) {
    p = p || ( r.best[ i ] != infinite_cost< Cost >() );
    r.best[ i ] = r.hv[ i ];
    }
  return( p );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// true if the path of cost g to a state of heuristic h is pruned by the
 /// bound

 bool pruned( Cost g , Cost h ) const noexcept {
  if( std::isinf( f_options.bound ) )
   return( false );
  const double v = double( g ) + ( f_options.prune_by_f ? double( h ) : 0 );
  return( v >= f_options.bound );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the successors of n into r.kids, with their preference

 void children( Run & r , const NodeT * n , SearchStats & st ) const {
  r.kids.clear();
  const auto add = [ & ]( const State & t , Cost c , Cost ht ,
                          std::uint8_t op ) {
   r.kids.push_back( Kid{ t , c , ht , op , false } );
   };
  for_each_child( f_domain , n->state , n->h , n->op , add );
  st.generated += r.kids.size();
  arrange( r , n );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// marks the preferred successors of n in r.kids, if the queues need
 /// them, and puts the successors in the order of the options

 void arrange( Run & r , const NodeT * n ) const {
  if constexpr( PreferredDomain< D > )
   if( f_preferred ) {
    const auto p = f_domain.preferences( n->state );
    for( Kid & k : r.kids )
     k.preferred = bool( p( k.s ) );
    }
  if( f_options.randomize_successors )
   std::shuffle( r.kids.begin() , r.kids.end() , r.rng );
  if( f_options.preferred_first )
   std::stable_partition( r.kids.begin() , r.kids.end() ,
                          []( const Kid & k ) { return( k.preferred ); } );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// fills the solution with the path to n, and the statistics

 void finish( Run & r , SolutionT & sol , const NodeT * goal ,
              const Stopwatch & clock ) const {
  SearchStats & st = sol.stats;
  if( goal ) {
   sol.status = Status::solved;
   sol.path = trace_path( goal );
   // a node of the path may have got a cheaper g after the goal
   sol.cost = std::min( goal->g , Cost( path_cost( f_domain , sol.path ) ) );
   }
  st.stored = r.pool.size();
  st.expanded_by = { st.expanded };
  st.stored_by = { st.stored };
  st.seconds = clock.seconds();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the search with eager evaluation

 SolutionT solve_eager( void ) const {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  Run r( *this );

  const State s0 = f_domain.initial_state();
  const std::uint64_t h0 = f_domain.hash( s0 );
  NodeT * root = make( r , s0 , h0 , Cost( 0 ) , Cost( 0 ) , nullptr ,
                       no_op );
  if( ( ! evaluate( r , root , Cost( f_domain.heuristic( s0 ) ) , h0 , st ) )
      || pruned( Cost( 0 ) , root->h ) ) {
   root->closed = true;
   finish( r , sol , nullptr , clock );
   return( sol );
   }
  progress( r );
  if( f_options.goal_at_generation && f_domain.is_goal( s0 ) ) {
   finish( r , sol , root , clock );
   return( sol );
   }
  r.open_e.push( Eager{ root , root->g } , r.hv.data() , root->g , false );

  const NodeT * goal = nullptr;
  Eager e;
  while( r.open_e.pop( e ) ) {
   NodeT * n = e.node;
   if( n->closed || ( n->g != e.g ) )
    continue;
   n->closed = true;
   ++st.expanded;
   if( n->expanded )
    ++st.reexpanded;
   n->expanded = true;
   if( ( ! f_options.goal_at_generation ) && f_domain.is_goal( n->state ) ) {
    goal = n;
    break;
    }

   children( r , n , st );
   for( const Kid & k : r.kids ) {
    const Cost g = n->g + k.c;
    if( pruned( g , k.h ) )
     continue;
    const std::uint64_t hash = f_domain.hash( k.s );
    bool fresh = false;
    NodeT * m = r.table.find_or_insert( k.s , hash , [ & ] {
     fresh = true;
     NodeT * x = r.pool.allocate();
     *x = NodeT{ k.s , g , k.h , n , false , false , k.op ,
                 std::uint32_t( r.pool.size() - 1 ) };
     return( x );
     } ).first;
    if( fresh ) {
     if( ! evaluate( r , m , k.h , hash , st ) ) {
      m->closed = true;                       // a dead end
      continue;
      }
     if( f_options.goal_at_generation && f_domain.is_goal( k.s ) ) {
      goal = m;
      break;
      }
     if( progress( r ) )
      r.open_e.boost();
     r.open_e.push( Eager{ m , g } , r.hv.data() , g , k.preferred );
     continue;
     }
    if( ( g < m->g ) && ( m->h != infinite_cost< Cost >() ) &&
        ( m->closed ? f_options.reopen : f_options.improve_open ) ) {
     m->g = g;
     m->parent = n;
     m->op = k.op;
     m->closed = false;
     values( r , m );
     r.open_e.push( Eager{ m , g } , r.hv.data() , g , k.preferred );
     }
    else
     ++st.duplicates;
    }
   if( goal )
    break;

   if( f_options.node_limit && ( r.pool.size() >= f_options.node_limit ) ) {
    sol.status = Status::node_limit;
    break;
    }
   }

  finish( r , sol , goal , clock );
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the search with lazy evaluation

 SolutionT solve_lazy( void ) const {
  Stopwatch clock;
  SolutionT sol;
  SearchStats & st = sol.stats;
  Run r( *this );

  // the current entry: the root first, then the records selected
  Lazy cur{ nullptr , f_domain.initial_state() , Cost( 0 ) , Cost( 0 ) ,
            no_op , false };
  const NodeT * goal = nullptr;
  for( ; ; ) {
   const std::uint64_t hash = f_domain.hash( cur.state );
   bool fresh = false;
   NodeT * n = r.table.find_or_insert( cur.state , hash , [ & ] {
    fresh = true;
    NodeT * x = r.pool.allocate();
    *x = NodeT{ cur.state , cur.g , Cost( 0 ) , cur.parent , false , false ,
                cur.op , std::uint32_t( r.pool.size() - 1 ) };
    return( x );
    } ).first;
   bool expand = false;
   if( fresh ) {
    const Cost h = cur.known ? cur.h
                             : Cost( f_domain.heuristic( cur.state ) );
    if( evaluate( r , n , h , hash , st ) && ! pruned( n->g , n->h ) )
     expand = true;
    else
     n->closed = true;                    // a dead end, or pruned by f
    }
   else
    if( ( cur.g < n->g ) && ( n->h != infinite_cost< Cost >() ) &&
        f_options.reopen ) {
     n->g = cur.g;
     n->parent = cur.parent;
     n->op = cur.op;
     values( r , n );
     expand = true;
     }
    else
     ++st.duplicates;

   if( expand ) {
    n->closed = true;
    ++st.expanded;
    if( n->expanded )
     ++st.reexpanded;
    n->expanded = true;
    if( progress( r ) )
     r.open_l.boost();
    if( f_domain.is_goal( n->state ) ) {
     goal = n;
     break;
     }
    lazy_children( r , n , st );
    for( const Kid & k : r.kids ) {
     const Cost g = n->g + k.c;
     if( pruned( g , Cost( 0 ) ) )     // by f when it is evaluated
      continue;
     r.records.push_back( Lazy{ n , k.s , g , k.h , k.op ,
                                ! DeferredDomain< D > } );
     r.open_l.push( r.records.size() - 1 , r.hv.data() , g , k.preferred );
     }
    if( f_options.node_limit && ( r.pool.size() >= f_options.node_limit ) ) {
     sol.status = Status::node_limit;
     break;
     }
    }

   std::size_t next;
   if( ! r.open_l.pop( next ) )
    break;
   cur = r.records[ next ];
   }

  finish( r , sol , goal , clock );
  return( sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the successors of n into r.kids for the lazy search: without their h
 /// on a DeferredDomain

 void lazy_children( Run & r , const NodeT * n , SearchStats & st ) const {
  if constexpr( DeferredDomain< D > ) {
   r.kids.clear();
   f_domain.for_each_successor_state( n->state , [ & ]( const State & t ,
                                                        Cost c ) {
    r.kids.push_back( Kid{ t , c , Cost( 0 ) , no_op , false } );
    } );
   st.generated += r.kids.size();
   arrange( r , n );
   }
  else
   children( r , n , st );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const D & f_domain;                                 ///< the Domain
 GreedyOptions f_options;                            ///< the options
 std::size_t f_m;                                    ///< the heuristics
 HeuristicCache< D > * f_cache = nullptr;            ///< the cache, if any
 bool f_preferred = false;                           ///< the preferred
                                                     ///< successors matter

 };  // end( class( GreedySearch ) )

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/// searches a Domain with GreedySearch and the given options

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > gbfs(
 const D & domain , const GreedyOptions & options = {} )
{
 return( GreedySearch< D >( domain , options ).solve() );
 }

/*--------------------------------------------------------------------------*/

}  // end( namespace PBeFS_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_GreedySearch */

/*--------------------------------------------------------------------------*/
/*------------------------ End File GreedySearch.h -------------------------*/
/*--------------------------------------------------------------------------*/
