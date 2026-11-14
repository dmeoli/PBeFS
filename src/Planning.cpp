/*--------------------------------------------------------------------------*/
/*---------------------------- File Planning.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the class PlanningBase: the packing of the variables,
 * the successor generator, the codes of the operators and the
 * abstractions of a Planning domain.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <bit>

#include <stdexcept>

#include <utility>

#include "PBeFS/astar/Planning.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace PBeFS_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

PlanningHeuristic PBeFS_di_unipi_it::planning_heuristic(
                                                    const std::string & n )
{
 if( n == "blind" )
  return( PlanningHeuristic::blind );
 if( n == "goalcount" )
  return( PlanningHeuristic::goal_count );
 if( n == "hmax" )
  return( PlanningHeuristic::hmax );
 if( n == "hadd" )
  return( PlanningHeuristic::hadd );
 if( n == "hff" )
  return( PlanningHeuristic::hff );
 if( n == "pdb" )
  return( PlanningHeuristic::pdb );
 if( n == "lmcut" )
  return( PlanningHeuristic::lmcut );
 throw( std::invalid_argument( "planning_heuristic: unknown heuristic " +
                               n ) );
 }

/*--------------------------------------------------------------------------*/

namespace {

using Fact = SASTask::Fact;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the node of the generator for the operators of items, each with the
// index of its next precondition to test, appended to gen

struct GenBuilder {

 const std::vector< std::vector< Fact > > & pre;
 const std::vector< int > & range;
 std::vector< int > var;
 std::vector< std::vector< std::uint32_t > > ops;
 std::vector< std::vector< int > > child;
 std::vector< int > dont_care;

 int build( std::vector< std::pair< std::uint32_t , std::size_t > > items ) {
  const int n = int( var.size() );
  var.push_back( -1 );
  ops.emplace_back();
  child.emplace_back();
  dont_care.push_back( -1 );
  int v = -1;
  for( const auto & [ o , k ] : items )
   if( k == pre[ o ].size() )
    ops[ std::size_t( n ) ].push_back( o );
   else
    if( ( v < 0 ) || ( pre[ o ][ k ].var < v ) )
     v = pre[ o ][ k ].var;
  if( v < 0 )
   return( n );
  var[ std::size_t( n ) ] = v;
  std::vector< std::vector< std::pair< std::uint32_t , std::size_t > > >
   by( std::size_t( range[ std::size_t( v ) ] ) );
  std::vector< std::pair< std::uint32_t , std::size_t > > rest;
  for( const auto & [ o , k ] : items )
   if( k < pre[ o ].size() ) {
    if( pre[ o ][ k ].var == v )
     by[ std::size_t( pre[ o ][ k ].val ) ].emplace_back( o , k + 1 );
    else
     rest.emplace_back( o , k );
    }
  std::vector< int > c( by.size() , -1 );
  for( std::size_t d = 0 ; d < by.size() ; ++d )
   if( ! by[ d ].empty() )
    c[ d ] = build( std::move( by[ d ] ) );
  child[ std::size_t( n ) ] = std::move( c );
  if( ! rest.empty() )
   dont_care[ std::size_t( n ) ] = build( std::move( rest ) );
  return( n );
  }

 };  // end( struct( GenBuilder ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the edges of the projection of the task onto vars, in both directions,
// as a list of neighbours per abstract state; mult are the multipliers of
// vars; if max_out is not nullptr, the largest number of successors of an
// abstract state, the edges taken in their direction, is written there

std::vector< std::vector< std::uint32_t > > abstract_graph(
                        const SASTask & t ,
                        const std::vector< std::vector< Fact > > & pre ,
                        const std::vector< int > & vars ,
                        const std::vector< std::size_t > & mult ,
                        std::size_t size , std::size_t * max_out )
{
 std::vector< int > where( t.num_vars() , -1 );
 for( std::size_t k = 0 ; k < vars.size() ; ++k )
  where[ std::size_t( vars[ k ] ) ] = int( k );
 std::vector< std::vector< std::uint32_t > > g( size );
 std::vector< int > fixed( vars.size() ) , x( vars.size() );
 for( std::size_t o = 0 ; o < t.ops.size() ; ++o ) {
  bool moves = false;
  for( const auto & e : t.ops[ o ].effects )
   moves = moves || ( where[ std::size_t( e.var ) ] >= 0 );
  if( ! moves )
   continue;
  std::fill( fixed.begin() , fixed.end() , -1 );
  for( const auto & p : pre[ o ] )
   if( where[ std::size_t( p.var ) ] >= 0 )
    fixed[ std::size_t( where[ std::size_t( p.var ) ] ) ] = p.val;
  // every abstract state matching the preconditions, as an odometer
  for( std::size_t k = 0 ; k < vars.size() ; ++k )
   x[ k ] = fixed[ k ] >= 0 ? fixed[ k ] : 0;
  for( ; ; ) {
   std::size_t b = 0 , c = 0;
   for( std::size_t k = 0 ; k < vars.size() ; ++k )
    b += std::size_t( x[ k ] ) * mult[ k ];
   c = b;
   for( const auto & e : t.ops[ o ].effects ) {
    const int k = where[ std::size_t( e.var ) ];
    if( k >= 0 )
     c = c - std::size_t( x[ std::size_t( k ) ] ) * mult[ std::size_t( k ) ]
         + std::size_t( e.post ) * mult[ std::size_t( k ) ];
    }
   if( c != b )
    g[ b ].push_back( std::uint32_t( c ) );
   std::size_t k = 0;
   for( ; k < vars.size() ; ++k ) {
    if( fixed[ k ] >= 0 )
     continue;
    if( ++x[ k ] < t.range[ std::size_t( vars[ k ] ) ] )
     break;
    x[ k ] = 0;
    }
   if( k == vars.size() )
    break;
   }
  }
 for( auto & x : g ) {
  std::sort( x.begin() , x.end() );
  x.erase( std::unique( x.begin() , x.end() ) , x.end() );
  }
 if( max_out ) {
  *max_out = 0;
  for( const auto & x : g )
   *max_out = std::max( *max_out , x.size() );
  return( g );
  }
 // undirected, the predecessors of a node being among its successors
 for( std::size_t b = 0 ; b < size ; ++b )
  for( std::size_t k = 0 , m = g[ b ].size() ; k < m ; ++k )
   g[ g[ b ][ k ] ].push_back( std::uint32_t( b ) );
 for( auto & s : g ) {
  std::sort( s.begin() , s.end() );
  s.erase( std::unique( s.begin() , s.end() ) , s.end() );
  }
 return( g );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the word, the first bit and the bits of the field of each variable, the
// fields filling the words in the order of the variables, none across two
// words; returns the number of words

std::size_t fields( const SASTask & t , std::vector< std::uint32_t > & word ,
                    std::vector< std::uint32_t > & shift ,
                    std::vector< std::uint32_t > & bits )
{
 std::uint32_t w = 0 , used = 0;
 for( std::size_t v = 0 ; v < t.num_vars() ; ++v ) {
  const auto b = std::uint32_t( std::max( 1 , int( std::bit_width(
                    unsigned( t.range[ v ] - 1 ) ) ) ) );
  if( used + b > 64 ) {
   ++w;
   used = 0;
   }
  word.push_back( w );
  shift.push_back( used );
  bits.push_back( b );
  used += b;
  }
 return( std::size_t( w ) + 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// true if the operator b is an inverse of the operator a

bool inverse( const SASTask & t , const std::vector< std::vector< Fact > > &
              pre , std::size_t a , std::size_t b )
{
 const auto & ea = t.ops[ a ].effects;
 const auto & eb = t.ops[ b ].effects;
 if( ea.size() != eb.size() )
  return( false );
 for( const auto & x : ea ) {
  if( x.pre < 0 )
   return( false );
  const auto y = std::find_if( eb.begin() , eb.end() ,
                               [ & ]( const auto & e ) {
                                return( e.var == x.var );
                                } );
  if( ( y == eb.end() ) || ( y->post != x.pre ) )
   return( false );
  }
 // the preconditions of b hold after a
 for( const auto & p : pre[ b ] ) {
  const auto x = std::find_if( ea.begin() , ea.end() ,
                               [ & ]( const auto & e ) {
                                return( e.var == p.var );
                                } );
  if( x != ea.end() ) {
   if( x->post != p.val )
    return( false );
   continue;
   }
  const bool known = std::any_of( pre[ a ].begin() , pre[ a ].end() ,
                                  [ & ]( const Fact & q ) {
                                   return( ( q.var == p.var ) &&
                                           ( q.val == p.val ) );
                                   } );
  if( ! known )
   return( false );
  }
 return( true );
 }

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS PlanningBase --------------------------*/
/*--------------------------------------------------------------------------*/

std::size_t PlanningBase::words_for( const SASTask & t )
{
 std::vector< std::uint32_t > w , s , b;
 return( fields( t , w , s , b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

PlanningBase::PlanningBase( SASTask task , const PlanningOptions & o ,
                            std::size_t max_words )
 : p_task( std::make_shared< const SASTask >( std::move( task ) ) ) ,
   f_heuristic( o.heuristic )
{
 const SASTask & t = *p_task;
 p_strips = std::make_shared< const StripsTask >( t );
 const std::size_t n = t.num_vars();
 if( t.ops.size() >= ( std::size_t( 1 ) << 32 ) )
  throw( std::invalid_argument( "PlanningBase: too many operators" ) );

 // the fields of the variables, none across two words
 std::vector< std::uint32_t > bits;
 f_words = fields( t , v_word , v_shift , bits );
 v_word_vars.resize( f_words );
 for( std::size_t v = 0 ; v < n ; ++v ) {
  v_mask.push_back( ( std::uint64_t( 1 ) << bits[ v ] ) - 1 );
  v_word_vars[ v_word[ v ] ].push_back( int( v ) );
  }
 if( f_words > max_words )
  throw( std::invalid_argument( "PlanningBase: the state needs " +
                                std::to_string( f_words ) +
                                " words, more than " +
                                std::to_string( max_words ) ) );
 v_init.assign( f_words , 0 );
 for( std::size_t v = 0 ; v < n ; ++v )
  v_init[ v_word[ v ] ] |= std::uint64_t( t.init[ v ] ) << v_shift[ v ];

 // the operators: preconditions, effects, costs
 std::vector< std::vector< Fact > > pre( t.ops.size() );
 v_eff_first.push_back( 0 );
 for( std::size_t a = 0 ; a < t.ops.size() ; ++a ) {
  pre[ a ] = t.preconditions( a );
  for( std::size_t k = 1 ; k < pre[ a ].size() ; ++k )
   if( pre[ a ][ k ].var == pre[ a ][ k - 1 ].var )
    throw( std::invalid_argument( "PlanningBase: the operator " +
                                  t.ops[ a ].name + " asks two values of "
                                  "a variable" ) );
  for( const auto & e : t.ops[ a ].effects ) {
   v_eff_var.push_back( std::uint32_t( e.var ) );
   v_eff_post.push_back( std::uint64_t( e.post ) );
   }
  v_eff_first.push_back( std::uint32_t( v_eff_var.size() ) );
  v_cost.push_back( t.ops[ a ].cost );
  }

 // the successor generator
 {
  GenBuilder b{ pre , t.range , {} , {} , {} , {} };
  std::vector< std::pair< std::uint32_t , std::size_t > > all;
  for( std::size_t a = 0 ; a < t.ops.size() ; ++a )
   all.emplace_back( std::uint32_t( a ) , 0 );
  b.build( std::move( all ) );
  for( std::size_t x = 0 ; x < b.var.size() ; ++x ) {
   GenNode g;
   g.var = b.var[ x ];
   g.ops_begin = std::uint32_t( v_gen_ops.size() );
   v_gen_ops.insert( v_gen_ops.end() , b.ops[ x ].begin() ,
                     b.ops[ x ].end() );
   g.ops_end = std::uint32_t( v_gen_ops.size() );
   g.child_begin = std::uint32_t( v_gen_child.size() );
   v_gen_child.insert( v_gen_child.end() , b.child[ x ].begin() ,
                       b.child[ x ].end() );
   g.dont_care = b.dont_care[ x ];
   v_gen.push_back( g );
   }
  }

 // the codes of the operators with an inverse, and their inverses
 {
  std::vector< std::vector< std::uint32_t > > inv;
  v_code.assign( t.ops.size() , no_op );
  // the operators by the sorted variables of their effects
  std::vector< std::pair< std::vector< int > , std::uint32_t > > key;
  for( std::size_t a = 0 ; a < t.ops.size() ; ++a ) {
   std::vector< int > vs;
   for( const auto & e : t.ops[ a ].effects )
    vs.push_back( e.var );
   std::sort( vs.begin() , vs.end() );
   key.emplace_back( std::move( vs ) , std::uint32_t( a ) );
   }
  std::sort( key.begin() , key.end() );
  for( std::size_t a = 0 ; ( a < t.ops.size() ) && ( inv.size() < no_op ) ;
       ++a ) {
   std::vector< int > vs;
   for( const auto & e : t.ops[ a ].effects )
    vs.push_back( e.var );
   std::sort( vs.begin() , vs.end() );
   auto it = std::lower_bound( key.begin() , key.end() ,
                               std::pair( vs , std::uint32_t( 0 ) ) );
   std::vector< std::uint32_t > l;
   for( ; ( it != key.end() ) && ( it->first == vs ) ; ++it )
    if( inverse( t , pre , a , it->second ) )
     l.push_back( it->second );
   if( l.empty() )
    continue;
   v_code[ a ] = std::uint8_t( inv.size() );
   inv.push_back( std::move( l ) );
   }
  v_inv_first.push_back( 0 );
  for( const auto & l : inv ) {
   v_inv.insert( v_inv.end() , l.begin() , l.end() );
   v_inv_first.push_back( std::uint32_t( v_inv.size() ) );
   }
  }

 // the abstract features: GreedyAFG on the variables but the fraction
 // ignored of largest fluency
 {
  std::vector< std::size_t > changes( n , 0 );
  for( const auto & op : t.ops )
   for( const auto & e : op.effects )
    if( e.pre != e.post )
     ++changes[ std::size_t( e.var ) ];
  std::vector< std::size_t > order( n );
  for( std::size_t v = 0 ; v < n ; ++v )
   order[ v ] = v;
  std::stable_sort( order.begin() , order.end() ,
                    [ & ]( std::size_t a , std::size_t b ) {
                     return( changes[ a ] > changes[ b ] );
                     } );
  std::vector< char > ignored( n , 0 );
  const auto ni = std::size_t( std::clamp( o.ignored , 0.0 , 1.0 ) *
                               double( n ) );
  for( std::size_t k = 0 ; k < ni ; ++k )
   ignored[ order[ k ] ] = 1;
  v_abstract.assign( std::size_t( p_strips->num_facts ) , 0 );
  f_num_abstract = 1;                     // 0 for the ignored variables
  for( std::size_t v = 0 ; v < n ; ++v ) {
   if( ignored[ v ] )
    continue;
   const auto r = std::size_t( t.range[ v ] );
   const auto off = std::size_t( p_strips->offset[ v ] );
   // the undirected domain transition graph, without repeated edges
   std::vector< std::vector< char > > edge( r , std::vector< char >( r , 0 ) );
   for( const auto & op : t.ops )
    for( const auto & e : op.effects )
     if( e.var == int( v ) )
      for( std::size_t d = 0 ; d < r ; ++d )
       if( ( ( e.pre < 0 ) || ( std::size_t( e.pre ) == d ) ) &&
           ( d != std::size_t( e.post ) ) )
        edge[ d ][ std::size_t( e.post ) ] =
         edge[ std::size_t( e.post ) ][ d ] = 1;
   std::vector< std::size_t > deg( r , 0 );
   for( std::size_t d = 0 ; d < r ; ++d )
    for( std::size_t x = 0 ; x < r ; ++x )
     deg[ d ] += std::size_t( edge[ d ][ x ] );
   // S1: the value of least degree, then the value with the most edges
   // into S1, until S1 holds half of the values; S2: the others
   std::vector< char > in1( r , 0 );
   std::size_t first = 0;
   for( std::size_t d = 1 ; d < r ; ++d )
    if( deg[ d ] < deg[ first ] )
     first = d;
   in1[ first ] = 1;
   for( std::size_t size1 = 1 ; 2 * size1 < r ; ++size1 ) {
    std::size_t best = r , best_links = 0;
    for( std::size_t d = 0 ; d < r ; ++d ) {
     if( in1[ d ] )
      continue;
     std::size_t links = 0;
     for( std::size_t x = 0 ; x < r ; ++x )
      links += std::size_t( in1[ x ] && edge[ d ][ x ] );
     if( ( best == r ) || ( links > best_links ) ) {
      best = d;
      best_links = links;
      }
     }
    in1[ best ] = 1;
    }
   for( std::size_t d = 0 ; d < r ; ++d )
    v_abstract[ off + d ] = f_num_abstract + ( in1[ d ] ? 0 : 1 );
   f_num_abstract += 2;
   }
  }

 // the nblocks, built at their first use
 p_nb = std::make_shared< NBlocks >();
 p_nb->max_nblocks = o.max_nblocks;
 p_nb->max_facts = o.nblock_facts;

 // the heuristic
 switch( f_heuristic ) {
  case( PlanningHeuristic::hmax ):
   p_relaxed = std::make_shared< const RelaxedHeuristic >(
                *p_strips , RelaxedHeuristic::Kind::hmax );
   break;
  case( PlanningHeuristic::hadd ):
   p_relaxed = std::make_shared< const RelaxedHeuristic >(
                *p_strips , RelaxedHeuristic::Kind::hadd );
   break;
  case( PlanningHeuristic::hff ):
   p_relaxed = std::make_shared< const RelaxedHeuristic >(
                *p_strips , RelaxedHeuristic::Kind::hff );
   break;
  case( PlanningHeuristic::lmcut ):
   p_lmcut = std::make_shared< const LMCut >( *p_strips );
   break;
  case( PlanningHeuristic::pdb ):
   p_pdb = std::make_shared< const SASPatternDatabase >( t ,
            o.pattern.empty() ? SASPatternDatabase::greedy_pattern(
                                 t , o.pdb_size ) : o.pattern );
   break;
  default:
   break;
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void PlanningBase::build_nblocks( void ) const
{
 const SASTask & t = *p_task;
 const std::size_t n = t.num_vars();
 NBlocks & x = *p_nb;
 std::vector< std::vector< Fact > > pre( t.ops.size() );
 for( std::size_t a = 0 ; a < t.ops.size() ; ++a )
  pre[ a ] = t.preconditions( a );
 const double facts = double( p_strips->num_facts );
 std::vector< int > vars;
 std::vector< std::size_t > mult;
 std::size_t size = 1 , held = 0;
 for( ; ; ) {
  int best = -1;
  std::size_t best_deg = 0 , best_size = 0;
  for( std::size_t v = 0 ; v < n ; ++v ) {
   if( std::find( vars.begin() , vars.end() , int( v ) ) != vars.end() )
    continue;
   const auto r = std::size_t( t.range[ v ] );
   const auto sz = size * r;
   if( ( r < 2 ) || ( sz > x.max_nblocks ) || ( ( x.max_facts > 0 ) &&
       ( double( held + r ) > x.max_facts * facts ) ) )
    continue;
   auto vs = vars;
   auto ms = mult;
   vs.push_back( int( v ) );
   ms.push_back( size );
   std::size_t deg = 0;
   ( void ) abstract_graph( t , pre , vs , ms , sz , & deg );
   if( ( best < 0 ) || ( deg < best_deg ) ||
       ( ( deg == best_deg ) && ( sz > best_size ) ) ) {
    best = int( v );
    best_deg = deg;
    best_size = sz;
    }
   }
  if( best < 0 )
   break;
  vars.push_back( best );
  mult.push_back( size );
  held += std::size_t( t.range[ std::size_t( best ) ] );
  size = best_size;
  }
 const auto g = abstract_graph( t , pre , vars , mult , size , nullptr );
 x.vars = std::move( vars );
 x.mult = std::move( mult );
 x.first.assign( 1 , 0 );
 for( const auto & s : g ) {
  x.succ.insert( x.succ.end() , s.begin() , s.end() );
  x.first.push_back( x.succ.size() );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

int PlanningBase::evaluate( const std::vector< int > & s ) const
{
 switch( f_heuristic ) {
  case( PlanningHeuristic::blind ):
   return( 0 );
  case( PlanningHeuristic::goal_count ): {
   int h = 0;
   for( const auto & g : p_task->goal )
    h += s[ std::size_t( g.var ) ] != g.val;
   return( h );
   }
  case( PlanningHeuristic::lmcut ):
   return( ( *p_lmcut )( s ) );
  case( PlanningHeuristic::pdb ):
   return( ( *p_pdb )( s ) );
  default:
   return( ( *p_relaxed )( s ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File Planning.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
