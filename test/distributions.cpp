/*--------------------------------------------------------------------------*/
/*------------------------- File distributions.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The cases of the tester of PBeFS [see test.cpp] for the work
 * distributions: those built automatically on a sample and by hand, and
 * the owners of HDA* on the 15-puzzle.
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
// the 15-puzzle without its feature groups, which a FeatureModel then
// infers from the moves

struct UngroupedFifteen : public FifteenPuzzle {
 using FifteenPuzzle::FifteenPuzzle;
 void num_feature_groups( void ) const = delete;
 };
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// reduce() against the high 64 bits of h p, computed by the long
// multiplication of the 16-bit digits of h and p, on edge and random values

static void check_reduce( void )
{
 const auto reference = []( std::uint64_t h , unsigned p ) {
  std::uint64_t d[ 6 ] = {};                   // h p, 16 bits per digit
  for( int i = 0 ; i < 4 ; ++i )
   for( int j = 0 ; j < 2 ; ++j )
    d[ i + j ] += ( ( h >> ( 16 * i ) ) & 0xFFFF ) *
                  ( ( std::uint64_t( p ) >> ( 16 * j ) ) & 0xFFFF );
  for( int i = 0 ; i < 5 ; ++i ) {
   d[ i + 1 ] += d[ i ] >> 16;
   d[ i ] &= 0xFFFF;
   }
  return( unsigned( d[ 4 ] | ( d[ 5 ] << 16 ) ) );  // the digits 4 and 5
  };
 const std::vector< std::uint64_t > hs = { 0 , 1 , 0xFFFFFFFFu ,
  std::uint64_t( 1 ) << 32 , ( std::uint64_t( 1 ) << 32 ) + 0xFFFFFFFFu ,
  std::uint64_t( 1 ) << 63 , ~std::uint64_t( 0 ) , ~std::uint64_t( 0 ) - 1 ,
  0x0123456789ABCDEFu , 0xFFFFFFFF00000000u };
 const std::vector< unsigned > ps = { 1 , 2 , 3 , 7 , 64 , 65535 , 65536 ,
                                      0x7FFFFFFFu , 0x80000000u ,
                                      0xFFFFFFFEu , 0xFFFFFFFFu };
 std::size_t bad = 0 , n = 0;
 const auto one = [ & ]( std::uint64_t h , unsigned p ) {
  ++n;
  const unsigned r = reduce( h , p );
  if( ( r != reference( h , p ) ) || ( r >= p ) )
   ++bad;
  };
 for( auto h : hs )
  for( auto p : ps )
   one( h , p );
 std::mt19937_64 gen( 12345 );
 for( int k = 0 ; k < 100000 ; ++k ) {
  const unsigned p = unsigned( gen() >> ( 32 + gen() % 32 ) );
  one( gen() , p ? p : 1 );
  }
 check( bad == 0 , "reduce( h , p ) is the high word of h p on " +
                   std::to_string( n ) + " pairs" );
 check( ( reduce( ~std::uint64_t( 0 ) , 10 ) == 9 ) &&
        ( reduce( std::uint64_t( 1 ) << 63 , 10 ) == 5 ) &&
        ( key_owner( 0xFFFFFFFFu , 8 ) == 7 ) ,
        "reduce( h , p ) at the ends of the range" );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the number of distinct abstract features of the tile t in a, -1 if a
// ignores some cell of t and not all

template< class D >
static int tile_regions_of( const std::vector< std::uint32_t > & a , int t )
{
 std::set< std::uint32_t > r;
 int ignored = 0;
 for( int c = 0 ; c < D::cells ; ++c ) {
  const auto x = a[ std::size_t( t * D::cells + c ) ];
  if( x == no_abstract_feature )
   ++ignored;
  else
   r.insert( x );
  }
 if( ignored == D::cells )
  return( 0 );
 return( ignored ? -int( r.size() ) : int( r.size() ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the value of the cut side of the DTG of the group g of m by the
// objective of sparsest_cut_abstraction() with
// CutObjective::sampled_sparsity

template< class D >
static double cut_sparsity( const FeatureModel< D > & m , std::size_t g ,
                            const std::vector< unsigned char > & side )
{
 const auto & v = m.members( g );
 std::map< std::uint32_t , std::size_t > at;
 for( std::size_t k = 0 ; k < v.size() ; ++k )
  at[ v[ k ] ] = k;
 double m0 = 0 , m1 = 0 , w = 0 , total = 0;
 for( std::size_t k = 0 ; k < v.size() ; ++k ) {
  ( side[ k ] ? m1 : m0 ) += m.mass( v[ k ] );
  for( const auto & e : m.arcs( v[ k ] ) )
   if( at.count( e.to ) && ( e.to > v[ k ] ) ) {
    total += e.w;
    if( side[ at[ e.to ] ] != side[ k ] )
     w += e.w;
    }
  }
 return( m0 * m1 / std::max( w , total * 1e-9 ) );
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the work distributions built automatically and by hand: the hand-made
// projections of the sliding-tile puzzle, the sample and the groups of a
// FeatureModel, the cuts of its DTGs against enumeration, the fractions
// of FAZHDA*, DAHDA* and OZHDA*, the estimates, and HDA* with them on the
// 24-puzzle

void test_distributions( const std::string & data )
{
 check_reduce();

 // the hand-made projections
 using F = FifteenPuzzle;
 using T = TwentyFourPuzzle;
 const auto regions15 = [ & ]( TileRegions r , int t ,
                               const std::vector< int > & tiles = {} ) {
  return( tile_regions_of< F >( F::tile_abstraction( r , tiles ) , t ) );
  };
 check( ( regions15( TileRegions::single , 1 ) == 16 ) &&
        ( regions15( TileRegions::pair , 5 ) == 8 ) &&
        ( regions15( TileRegions::line , 9 ) == 4 ) &&
        ( regions15( TileRegions::block , 15 ) == 4 ) &&
        ( regions15( TileRegions::half , 2 ) == 2 ) &&
        ( regions15( TileRegions::block , 4 , { 1 , 2 , 3 , 4 } ) == 4 ) &&
        ( regions15( TileRegions::block , 5 , { 1 , 2 , 3 , 4 } ) == 0 ) &&
        ( regions15( TileRegions::single , 7 , { 1 , 7 , 12 } ) == 16 ) &&
        ( regions15( TileRegions::single , 8 , { 1 , 7 , 12 } ) == 0 ) ,
        "the regions of the hand-made projections of the 15-puzzle" );
 const auto b = F::tile_abstraction( TileRegions::block );
 check( ( b[ 16 + 0 ] == b[ 16 + 5 ] ) && ( b[ 16 + 2 ] == b[ 16 + 7 ] ) &&
        ( b[ 16 + 0 ] != b[ 16 + 2 ] ) && ( b[ 16 + 10 ] == b[ 16 + 15 ] ) &&
        ( b[ 16 + 0 ] != b[ 32 + 0 ] ) ,
        "the blocks of 2 x 2 cells of the 15-puzzle" );
 bool thrown = false;
 try {
  (void) F::tile_abstraction( TileRegions::quad );
  }
 catch( const std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "the regions of the 24-puzzle refused on the 15-puzzle" );
 const auto q = T::tile_abstraction( TileRegions::quad );
 const auto cr = T::tile_abstraction( TileRegions::cross );
 const auto od = T::tile_abstraction( TileRegions::odd );
 std::map< std::uint32_t , int > nq , nc , no;
 for( int c = 0 ; c < T::cells ; ++c ) {
  ++nq[ q[ std::size_t( 25 + c ) ] ];
  ++nc[ cr[ std::size_t( 25 + c ) ] ];
  ++no[ od[ std::size_t( 25 + c ) ] ];
  }
 bool fives = ( nc.size() == 5 ) && ( no.size() == 5 );
 for( const auto & [ x , k ] : nc )
  fives = fives && ( k == 5 );
 for( const auto & [ x , k ] : no )
  fives = fives && ( k == 5 );
 check( ( nq.size() == 5 ) && ( nq[ no_abstract_feature ] == 1 ) &&
        ( q[ 25 + 12 ] == no_abstract_feature ) && fives &&
        ( od[ 25 + 2 ] == od[ 25 + 22 ] ) &&
        ( cr[ 25 + 12 ] == cr[ 25 + 7 ] ) ,
        "the regions quad, cross and odd of the 24-puzzle" );

 // the sample, the groups and the estimates
 const auto inst = korf( data );
 const F d( inst[ easy_korf[ 0 ] - 1 ].board );
 const FeatureModel< F > m( d , 8000 );
 bool grid = true;
 for( std::size_t i = 0 ; i < m.num_features() ; ++i )
  for( const auto & e : m.arcs( i ) ) {
   const int a = int( i % 16 ) , c = int( e.to % 16 );
   grid = grid && ( i / 16 == e.to / 16 ) &&
          ( std::abs( a / 4 - c / 4 ) + std::abs( a % 4 - c % 4 ) == 1 );
   }
 check( ( m.num_groups() == 15 ) && ( m.num_expanded() > 0 ) &&
        ( m.num_expanded() <= 8000 ) && ( m.edges().size() > 0 ) && grid ,
        "the groups and the DTGs of the 15-puzzle, the grid of each tile" );
 const UngroupedFifteen u( inst[ easy_korf[ 0 ] - 1 ].board );
 static_assert( ! FeatureGroupDomain< UngroupedFifteen > );
 const FeatureModel< UngroupedFifteen > mu( u , 8000 );
 bool same = mu.num_groups() == m.num_groups();
 for( std::size_t i = 0 ; same && ( i < m.num_features() ) ; ++i )
  for( std::size_t j = 0 ; same && ( j < m.num_features() ) ; ++j )
   if( ( m.group( i ) != FeatureModel< F >::no_group ) &&
       ( m.group( j ) != FeatureModel< F >::no_group ) )
    same = ( m.group( i ) == m.group( j ) ) ==
           ( mu.group( i ) == mu.group( j ) );
 check( same , "the groups inferred from the moves are the tiles" );

 const ZobristDistribution< F > z( d );
 const auto tz = make_distribution( m , "Z" );
 bool zeq = true;
 for( const auto & s : m.states() )
  zeq = zeq && ( z.owner( s , 7 ) == tz.owner( s , 7 ) );
 const auto e1 = m.estimate( z , 1 );
 const auto e8 = m.estimate( z , 8 );
 check( zeq && ( e1.co == 0 ) && ( e1.lb == 1 ) && ( e8.co > 0.8 ) &&
        ( e8.lb >= 1 ) && ( e8.lb < 1.2 ) &&
        ( std::abs( e8.so - ( e8.lb - 1 ) ) < 1e-12 ) &&
        ( std::abs( e8.efficiency - 1 / ( ( 1 + e8.co ) * e8.lb ) ) < 1e-12 ) ,
        "the table Z is Zobrist hashing, and its estimates on 1 and 8 "
        "workers" );

 // GAZHDA*: balanced halves; GRAZHDA*: the sparsest cut of each DTG
 const auto ga = greedy_abstraction( m );
 const auto gr = sparsest_cut_abstraction( m ,
                                           CutObjective::sampled_sparsity );
 bool halves = true , sparsest = true;
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g ) {
  const auto & v = m.members( g );
  std::map< std::uint32_t , int > n;
  for( auto i : v )
   ++n[ ga.abstract[ i ] ];
  halves = halves && ( n.size() == 2 ) &&
           ( std::abs( n.begin()->second - n.rbegin()->second ) <= 1 );
  std::vector< unsigned char > side( v.size() );
  for( std::size_t k = 0 ; k < v.size() ; ++k )
   side[ k ] = gr.abstract[ v[ k ] ] != gr.abstract[ v[ 0 ] ];
  const double best = cut_sparsity( m , g , side );
  if( g < 3 )                                    // 2^15 cuts of each
   for( std::uint32_t c = 1 ; c < ( 1u << ( v.size() - 1 ) ) ; ++c ) {
    std::vector< unsigned char > other( v.size() , 0 );
    for( std::size_t k = 1 ; k < v.size() ; ++k )
     other[ k ] = ( c >> ( k - 1 ) ) & 1u;
    sparsest = sparsest &&
               ( cut_sparsity( m , g , other ) <= best * ( 1 + 1e-9 ) );
    }
  }
 check( halves , "GAZHDA* cuts each DTG in halves" );
 check( sparsest , "GRAZHDA* with the weights of the sample cuts each DTG "
                   "at its sparsest cut" );

 // GRAZHDA* with the weights of the paper: the two halves of the board, of
 // two lines each, on the 15-puzzle, two lines against three on the
 // 24-puzzle [see Jinnai and Fukunaga, JAIR 2017, figure 19]
 const auto lines = []( const auto & a , int w , int tile ) {
  std::map< std::uint32_t , std::set< int > > rows , cols;
  for( int c = 0 ; c < w * w ; ++c ) {
   const auto x = a.abstract[ std::size_t( tile * w * w + c ) ];
   rows[ x ].insert( c / w );
   cols[ x ].insert( c % w );
   }
  if( rows.size() != 2 )
   return( 0 );
  bool byrow = true , bycol = true;
  for( const auto & [ x , r ] : rows )
   for( const auto & [ y , q ] : rows )
    if( x != y )
     for( int i : r )
      byrow = byrow && ! q.count( i );
  for( const auto & [ x , r ] : cols )
   for( const auto & [ y , q ] : cols )
    if( x != y )
     for( int i : r )
      bycol = bycol && ! q.count( i );
  if( ! ( byrow || bycol ) )
   return( 0 );
  const auto & side = byrow ? rows : cols;
  return( int( std::min( side.begin()->second.size() ,
                         side.rbegin()->second.size() ) ) );
  };
 const auto gp = sparsest_cut_abstraction( m );
 bool halves15 = true;
 for( int t = 1 ; t < 16 ; ++t )
  halves15 = halves15 && ( lines( gp , 4 , t ) == 2 );
 check( halves15 , "GRAZHDA* cuts the 15-puzzle in two halves" );
 const T t24( random_walk< T >( 40 , 1 ) );
 const FeatureModel< T > m24( t24 , 64 );
 const auto gp24 = sparsest_cut_abstraction( m24 );
 bool rows24 = true;
 for( int t = 1 ; t < 25 ; ++t )
  rows24 = rows24 && ( lines( gp24 , 5 , t ) == 2 );
 check( rows24 , "GRAZHDA* cuts the 24-puzzle in two lines and three" );

 // FAZHDA* ignores the most fluent groups, DAHDA* keeps 30% of the
 // features, OZHDA* keeps the key on its most frequent operator
 const auto fa = fluency_abstraction( m , 0.3 );
 std::vector< std::pair< double , std::size_t > > fl;
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g )
  fl.emplace_back( m.fluency( g ) , g );
 std::sort( fl.rbegin() , fl.rend() );
 int ign = 0;
 bool most = true;
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g )
  if( fa.abstract[ m.members( g ).front() ] == no_abstract_feature ) {
   ++ign;
   most = most && ( m.fluency( g ) >= fl[ 3 ].first );
   }
 check( ( ign == 4 ) && most , "FAZHDA* ignores the 4 most fluent tiles" );
 const auto da = dynamic_state_abstraction( m , 0.3 );
 std::size_t kept = 0 , all = 0;
 for( std::size_t g = 0 ; g < m.num_groups() ; ++g ) {
  all += m.members( g ).size();
  if( da.abstract[ m.members( g ).front() ] != no_abstract_feature )
   kept += m.members( g ).size();
  }
 check( ( kept > 0 ) && ( kept * 10 <= all * 3 ) &&
        ( ( kept + 16 ) * 10 > all * 3 ) ,
        "DAHDA* keeps the groups of at most 30% of the features" );
 // AHDA* with at most 200 abstract states on the sample
 const auto ah = TableZobristDistribution< F >( d , state_abstraction( m ,
                                                                    200 ) );
 std::set< std::uint64_t > abs_states;
 for( const auto & s : m.states() )
  abs_states.insert( ah.hash( s ) );
 check( ( abs_states.size() <= 200 ) && ( abs_states.size() > 1 ) ,
        "AHDA* keeps at most 200 abstract states on the sample" );
 const auto oz = operator_table( m , 0.2 );
 std::uint64_t x = 0;
 for( auto i : m.operators().front().features )
  x ^= oz[ i ];
 const auto ez = m.estimate( TableZobristDistribution< F >( d , oz ) , 8 );
 check( ( ( x >> 32 ) == 0 ) && ( ez.co < e8.co ) ,
        "OZHDA* keeps the key on the most frequent operator" );

 // KZHDA*: the key vanishes on every operator of the chosen subspace, and
 // the table estimates as the domain distribution built on it; a basis
 // reduces the dependent vectors to nothing
 Gf2Basis gb( 8 );
 const bool dep = gb.insert( { 0 , 1 } ) && gb.insert( { 1 , 2 } ) &&
                  ! gb.insert( { 0 , 2 } ) && gb.insert( { 5 } ) &&
                  ( gb.rank() == 3 );
 const auto kt = kernel_table( 8 , gb , 3 );
 check( dep && ( ( ( kt[ 0 ] ^ kt[ 1 ] ) >> 32 ) == 0 ) &&
        ( ( ( kt[ 1 ] ^ kt[ 2 ] ) >> 32 ) == 0 ) && ( ( kt[ 5 ] >> 32 ) == 0 )
        && ( ( kt[ 3 ] >> 32 ) != 0 ) ,
        "the key of kernel_table() vanishes on the subspace" );
 Gf2Basis topb( m.num_features() );
 std::size_t kz_kept = 0 , tried = 0;
 const auto & kops = m.operators();
 for( std::size_t j = 0 ; j < std::min< std::size_t >( 40 , kops.size() ) ;
      ++j )
  topb.insert( kops[ j ].features );
 const auto ktop = kernel_table( m.num_features() , topb , 7 );
 for( std::size_t j = 0 ; j < std::min< std::size_t >( 40 , kops.size() ) ;
      ++j ) {
  std::uint64_t x = 0;
  for( auto i : kops[ j ].features )
   x ^= ktop[ i ];
  kz_kept += ( x >> 32 ) == 0;
  ++tried;
  }
 const auto kz = make_distribution( m , "KZ" , 1 , 1 , 8 );
 const auto ekz = m.estimate( kz , 8 );
 const auto etab = m.estimate_table( kz.table() , 8 );
 check( ( kz_kept == tried ) &&
        ( std::abs( ekz.efficiency - etab.efficiency ) < 1e-12 ) &&
        ( ekz.efficiency + 1e-12 >= m.estimate_table(
           kernel_table( m.num_features() , Gf2Basis( m.num_features() ) ,
                         1 ) , 8 ).efficiency ) ,
        "KZHDA* keeps the key on its operators and is no worse than Z on "
        "the sample" );

 // the automatic choice is the best of the estimates
 DistributionEstimate best;
 const auto au = auto_distribution( m , 8 , 1 , {} , &best );
 bool top = au.name().rfind( "auto:" , 0 ) == 0;
 for( const auto & n : auto_distribution_names() )
  top = top && ( m.estimate( make_distribution( m , n , 1 , 1 , 8 ) , 8 )
                 .efficiency <=
                 best.efficiency );
 check( top , "the automatic choice " + au.name() +
              " has the largest estimated efficiency" );

 // a grid with 8 moves, whose single group is the cell, of 100 by 60
 // cells in a light run [see light_run()]
 const int gw = light_run() ? 100 : 300 , gh = light_run() ? 60 : 200;
 const auto gmap = GridMap::random( gw , gh , 0.25 , 2 );
 const Grid< 8 > g8( gmap , 0 , gh - 1 , gw - 1 , gh - 1 , 20 );
 const auto gref = astar( Blind< Grid< 8 > >{ g8 } );
 const FeatureModel< Grid< 8 > > gm( g8 , 4096 );
 for( const char * n : { "GAZ" , "GRAZ" , "GRAZS" , "OZ" , "DAHDA" } )
  for( unsigned p : { 1u , 3u } ) {
   ParallelOptions o;
   o.threads = p;
   const auto x = make_distribution( gm , n );
   check_solution( g8 , HDAStar< Grid< 8 > , TableZobristDistribution<
                                 Grid< 8 > > >( g8 , x , o ).solve() ,
                   gref.cost ,
                   "HDA*[" + x.name() + "] grid 8, p " + std::to_string( p ) );
   }

 // the 24-puzzle, whose DTGs of 25 nodes are cut by local search, solved
 // by HDA* with the distributions built on it and the hand-made ones, on
 // one walk only in a light run
 for( std::uint64_t seed = 1 ; seed <= ( light_run() ? 1u : 2u ) ; ++seed ) {
  const T t( random_walk< T >( 40 , seed ) );
  const auto ref = astar( t );
  const FeatureModel< T > mt( t , 4096 );
  ParallelOptions o;
  o.threads = 2;
  std::vector< TableZobristDistribution< T > > tt;
  for( const char * n : { "GAZ" , "GRAZ" , "GRAZS" , "FAZ" , "OZ" ,
                          "OZP" , "DAHDA" , "AHDA" , "KZ" } )
   tt.push_back( make_distribution( mt , n ) );
  const std::pair< TileRegions , const char * > regions[] = {
   { TileRegions::quad , "tile:quad" } ,
   { TileRegions::cross , "tile:cross" } ,
   { TileRegions::odd , "tile:odd" } , { TileRegions::line , "tile:line" } };
  for( const auto & [ r , name ] : regions )
   tt.emplace_back( t , FeatureAbstraction::of( T::tile_abstraction( r ) ) ,
                    1 , name );
  tt.emplace_back( t , FeatureAbstraction::of( T::tile_abstraction(
                        TileRegions::half , { 1 , 2 , 3 , 4 , 5 } ) ) , 1 ,
                   "tile:half:1-5" );
  for( const auto & y : tt )
   check_solution( t , HDAStar< T , TableZobristDistribution< T > >(
                        t , y , o ).solve() , ref.cost ,
                   "HDA*[" + y.name() + "] 24-puzzle, walk " +
                   std::to_string( seed ) );
  }

 // the planning tasks, whose groups are their variables and whose arcs are
 // those of their domain transition graphs [see FeatureGroupDomain and
 // FeatureTransitionDomain]: every distribution on blocks-6-2 and the
 // automatic one on logistics-4-0, optimal plans with 1 and 3 workers
 static_assert( FeatureGroupDomain< PlanningD > &&
                FeatureTransitionDomain< PlanningD > );
 std::map< std::string , int > optimal;
 for( const auto & [ name , cost ] : planning_tasks( data ) )
  optimal[ name ] = cost;
 for( const char * name : { "blocks-6-2.sas" , "logistics-4-0.sas" } ) {
  const auto pd = planning( data , name , PlanningHeuristic::lmcut );
  const FeatureModel< PlanningD > pm( pd , 4096 );
  bool groups = ( pm.num_groups() > 0 ) &&
                ( pm.num_groups() <= pd.num_feature_groups() );
  std::map< std::uint32_t , std::size_t > var;  // the variable of a group
  for( std::size_t i = 0 ; i < pm.num_features() ; ++i ) {
   if( pm.group( i ) != FeatureGraph::no_group )
    groups = groups && ( var.try_emplace( pm.group( i ) ,
                                          pd.feature_group( i ) ).first->
                         second == pd.feature_group( i ) );
   for( const auto & e : pm.static_arcs( i ) )
    groups = groups && ( pd.feature_group( e.to ) == pd.feature_group( i ) );
   }
  check( groups , std::string( "the groups and the DTGs of " ) + name );
  std::vector< std::string > names = { "AUTO" };
  if( std::string( name ) == "blocks-6-2.sas" )
   names.insert( names.end() , auto_distribution_names().begin() ,
                 auto_distribution_names().end() );
  for( unsigned p : { 1u , 3u } ) {
   ParallelOptions o;
   o.threads = p;
   for( const auto & n : names ) {
    const auto x = n == "AUTO" ? auto_distribution( pm , p )
                               : make_distribution( pm , n , 1 , 1 , p );
    check_plan( pd , HDAStar< PlanningD , TableZobristDistribution<
                               PlanningD > >( pd , x , o ).solve() ,
                optimal.at( name ) , "HDA*[" + x.name() + "] " + name +
                ", p " + std::to_string( p ) );
    }
   }
  }
 }
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// the owners of HDA* on the 15-puzzle, by counting: the 48 pairs of a cell
// of the blank and a move, 8 of which change the abstract feature of AZ;
// the perfect hash modulo p, a function of the relative order of the cells
// of the last m tiles, m the least with p dividing m!, on 200000 random
// boards for each p; and the keys of each owner by key_owner, which differ
// by at most one

void test_hda_owners( const std::string & )
{
 using D = FifteenPuzzle;
 constexpr int n = D::cells;
 int pairs = 0;
 int across = 0;
 for( int b = 0 ; b < n ; ++b )
  for( int c : { b - 4 , b - 1 , b + 1 , b + 4 } ) {
   if( ( c < 0 ) || ( c >= n ) ||
       ( ( ( c == b - 1 ) || ( c == b + 1 ) ) && ( c / 4 != b / 4 ) ) )
    continue;
   ++pairs;
   if( D::abstract_feature( 5 * n + c ) != D::abstract_feature( 5 * n + b ) )
    ++across;
   }
 check( ( pairs == 48 ) && ( across == 8 ) ,
        "AZ: 8 of the 48 moves of the blank change its abstract feature" );

 const D d( D::korf_goal() );
 std::mt19937_64 gen( 1 );
 for( unsigned p : { 3u , 5u , 7u , 8u , 16u } ) {
  int m = 1;
  unsigned long long f = 1;
  while( f % p ) {
   ++m;
   f *= unsigned( m );
   }
  std::map< std::vector< int > , unsigned > sig;
  int bad = 0;
  for( int it = 0 ; it < 200000 ; ++it ) {
   D::board_type b = D::korf_goal();
   std::shuffle( b.begin() , b.end() , gen );
   std::vector< int > pos( n );
   for( int c = 0 ; c < n ; ++c )
    pos[ b[ c ] ] = c;
   std::vector< int > key;             // the relative order of the last m
   for( int k = n - m ; k < n ; ++k )
    for( int j = k + 1 ; j < n ; ++j )
     key.push_back( pos[ j ] < pos[ k ] );
   const unsigned r = unsigned( d.perfect_hash( D::pack( b ) ) % p );
   const auto [ at , fresh ] = sig.emplace( key , r );
   if( ( ! fresh ) && ( at->second != r ) )
    ++bad;
   }
  check( bad == 0 , "P: the perfect hash modulo " + std::to_string( p ) +
                    " depends on the last tiles only" );
  }

 for( unsigned p : { 3u , 7u , 64u , 96u } ) {
  unsigned long long lo = ~0ull;
  unsigned long long hi = 0;
  bool ok = true;
  for( unsigned i = 0 ; i < p ; ++i ) {
   const unsigned long long a = ( unsigned long long )( i ) << 32;
   const unsigned long long z = ( unsigned long long )( i + 1 ) << 32;
   const unsigned long long first = ( a + p - 1 ) / p;
   const unsigned long long end = ( z + p - 1 ) / p;
   ok = ok && ( key_owner( std::uint32_t( first ) , p ) == i ) &&
        ( ( first == 0 ) ||
          ( key_owner( std::uint32_t( first - 1 ) , p ) != i ) );
   lo = std::min( lo , end - first );
   hi = std::max( hi , end - first );
   }
  check( ok && ( hi - lo <= 1 ) , "key_owner: the keys of the " +
         std::to_string( p ) + " owners differ by at most one" );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File distributions.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
