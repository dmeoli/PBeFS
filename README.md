# Parallel Best-First Search

[![CI](https://github.com/dmeoli/PBeFS/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/dmeoli/PBeFS/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/dmeoli/PBeFS)](https://github.com/dmeoli/PBeFS/releases)
[![Coverage Status](https://coveralls.io/repos/github/dmeoli/PBeFS/badge.svg?branch=main)](https://coveralls.io/github/dmeoli/PBeFS?branch=main)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://dmeoli.github.io/PBeFS/)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue)](https://en.cppreference.com/w/cpp/20)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)](https://github.com/dmeoli/PBeFS/actions/workflows/ci.yml)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

PBeFS is a C++20 library for the sequential and parallel best-first search
of two kinds of problems: the A* search [[9]](#ref-9) of an implicit state-space
graph, i.e., a graph given by an initial state, a goal predicate and a
procedure that generates the successors of a state (cf. [[10]](#ref-10)), and the
Branch-and-Bound (B&B) of a tree whose nodes have a bound and possibly a
feasible value. In fact, A* is the special case of the general B&B in which
the bound of a node is $f = g + h$ and a goal gives a feasible solution of
value g (cf. [[17]](#ref-17)), and the two searches share the parallel schemes and
the runtimes; they differ in that a graph has duplicates, which A* detects
with a closed list, while a tree has none but has an incumbent. The
parallel A* searches follow the taxonomy of the survey of Fukunaga et al.
[[1]](#ref-1), and the parallel B&B searches that of the survey of Gendron and
Crainic [[16]](#ref-16). They run on C++ threads, OpenMP, FastFlow [[12]](#ref-12) and MPI;
hence, one can compare the same algorithm across runtimes, and the same
runtime across algorithms. The headers are in three layers:
`include/PBeFS/core` holds the runtimes, the messages over MPI and the
helpers that do not depend on the kind of search, `include/PBeFS/astar`
the graph searches with their domains, and `include/PBeFS/bnb` the tree
searches; `core.h`, `astar.h`, `bnb.h` and `PBeFS.h` include a layer, or
all of them, save the parts that need OpenMP, FastFlow or MPI.

## The search problem

We are concerned with finding a path of minimum cost from the initial state
to a goal state, in a graph with nonnegative edge costs, when h is
admissible (it never overestimates the cost of reaching a goal). For the
searches that do not reopen the closed nodes we also assume that h is
consistent, i.e., $h(u) \leq c(u, v) + h(v)$ for each edge $(u, v)$. Under
these conditions the searches of the library return an optimal path (cf. [[2]](#ref-2) for
HDA* and [[3]](#ref-3) for Safe PBNF). Without consistency the sequential A* is still
exact, since it reopens a closed node reached through a cheaper path; the
parallel searches, however, may pay several re-expansions for the same node.

All the algorithms are templates over the problem, which is described by a
type satisfying the concept `Domain` of `include/PBeFS/astar/Domain.h`.
Such a type gives `state_type` and an arithmetic `cost_type`, together
with the methods
`initial_state()`, `is_goal( s )`, `heuristic( s )`, `hash( s )` and
`for_each_successor( s , h_s , f )`; the latter calls `f( t , c , h_t )` for
each successor t of s. Note that h(s) is passed in, so that h(t) can be
computed incrementally when the heuristic allows it (as the Manhattan
distance of the sliding-tile puzzle does). Since all calls are resolved at
compile time, no virtual function is paid in the inner loop, which is the
first design choice of the fast implementation of Burns et al. [[4]](#ref-4).

Some algorithms need more than a `Domain`, and each extra requirement is a
refinement of the concept. Thus, `FeatureDomain` describes a state as a set
of propositions (as Zobrist hashing requires), and `AbstractFeatureDomain`
adds a projection of the features, for Abstract Zobrist hashing [[5]](#ref-5).
`AbstractStateDomain` and `NBlockDomain` give an abstraction of the state
space (for AHDA* and PBNF [[3]](#ref-3)), `PairwiseDomain` a heuristic between any two
states (for PA*SE [[6]](#ref-6)), and `PerfectHashDomain` an injective map of the
states to the integers. As examples, the library provides the sliding-tile
puzzle of any size with the Manhattan distance (`SlidingTile.h`, which
satisfies all of the concepts) and pathfinding with 4 or 8 moves on a grid
with cell costs (`Grid.h`), read from a Moving AI map [[15]](#ref-15) or generated at
random. Also, `ExplicitGraph.h` holds a weighted graph given by its arcs,
which serves for the route-finding problem on the map of Romania of Russell
and Norvig [[13]](#ref-13).

## The algorithms

We start from the sequential search (`astar( d )` in `AStar.h`), which is A*
with the goal test at the expansion and a node reopened when a cheaper path
to it is found. Ties on f are broken in favour of the larger g, which on
domains with many states in the last f layer usually reaches the goal much
earlier. Next to it, `idastar( d )` (`IDAStar.h`) is the iterative-deepening
A* of Korf [[7]](#ref-7), a depth-first search bounded by f and repeated with
a growing bound, which keeps only the current path and hence solves the
problems whose CLOSED would not fit in memory, such as the 24-puzzle, at the
price of expanding again the states reached by more than one path. The
parallel searches are the following, each one with the
reference for its correctness and its typical behavior.

- `SPAStar` (`SPAStar.h`) is the Simple Parallel A* of Algorithm 2 of [[1]](#ref-1),
  in which the workers share one OPEN and one CLOSED list, protected by
  locks. It is the baseline of the survey, and it is not expected to scale
  beyond a few threads, since each expansion goes through the shared lists.
- `HDAStar` (`HDAStar.h`) is Hash Distributed A* [[2]](#ref-2): a distribution
  function assigns each state to a worker, which keeps OPEN and CLOSED for
  the states it owns and receives, packed in batches, those generated by
  the others. The search ends when no worker has a node with f smaller than
  the incumbent and no batch is travelling. The distribution is a template
  parameter (`Distribution.h`), and the library provides those compared by
  Jinnai and Fukunaga [[5]](#ref-5), i.e., Zobrist hashing (HDA*[Z]), Abstract
  Zobrist hashing (AZHDA*), Zobrist hashing of an abstract state (AHDA*
  [[3]](#ref-3)) and the perfect hash modulo the number of workers (HDA*[P], the
  "bad" case of [[3]](#ref-3)). Furthermore, it provides the hash of the closed list
  and, as the baseline without duplicate detection across workers, a random
  assignment.
- `HDAStarMPI` (`HDAStarMPI.h`) is the same search with one MPI process per
  worker, a nonblocking exchange of the batches and the termination
  detection of Mattern [[11]](#ref-11) (four counters reduced with `MPI_Iallreduce`).
  At the end, the path is rebuilt state by state, each one sent by the
  process that holds it.
- `SafePBNF` (`SafePBNF.h`) is the livelock-free Parallel Best-NBlock-First
  search of Burns et al. [[3]](#ref-3), in which an abstraction partitions the state
  space into nblocks. A worker expands for a while the nodes of an nblock
  whose duplicate detection scope (its abstract successors) is not in use
  by another worker; hence, no lock is taken on the single node.
- `PASE` (`PASE.h`) is the Parallel A* for Slow Expansions of Phillips et
  al. [[6]](#ref-6), which expands in parallel only the states that the pairwise
  heuristic proves not to be reachable through a cheaper path from a state
  still in OPEN. It may pay off only when an expansion costs much more than
  the bookkeeping; `SlowDomain.h` adds a given amount of artificial work per
  expansion to study this regime.
- `FarmAStar` (`FarmAStar.h`) is a centralized parallel A* written as a
  FastFlow farm with feedback, i.e., the master-worker scheme of [[1]](#ref-1): a
  master holds OPEN and CLOSED, and the workers generate the successors of
  the nodes it sends them.

The shared-memory searches are templates also over the runtime that starts
the workers, i.e., `ThreadRuntime`, `OpenMPRuntime` (`OpenMPRuntime.h`) or
`FastFlowRuntime` (`FastFlowRuntime.h`); the latter runs HDA* as a mesh of
FastFlow single-producer single-consumer queues. Their options are in
`ParallelOptions`: the number of workers, the size of a batch, the number of
expansions between two flushes of a partially filled batch and the seed of
the distribution.

## Branch-and-bound

The same machinery runs a best-first search with an incumbent, i.e., a
Branch-and-Bound (B&B), on a problem satisfying the concept
`BranchAndBound` of `BranchAndBound.h`. Such a problem gives a copyable
`node_type`, a `sense` (minimization or maximization), the root node, and
two methods: `evaluate( n )` returns the bound of the node n, whether it is
infeasible and possibly a feasible value, while `branch( n , sink )` gives
the children of n, each with a bound that is cheap to compute.
Alternatively, `make_worker( i )` gives each worker an evaluator of its
own (e.g., with a copy of a model and of a solver), which may also take
the incumbent value by `cutoff( z )` before each evaluation (e.g., to stop
a relaxation that cannot improve it). Note that a B&B tree has no duplicates;
hence, there is no closed list. A node is evaluated only when it is
extracted, which pays off when the bound is "expensive" (say, a Linear
Program), and it is pruned when it is generated, when it is extracted and
after its evaluation. Also, the nodes can be extracted best-first,
depth-first, breadth-first, by an estimate of the problem, or by plunging (a dive into
the best child, then best-first). The search ends when no open node can
improve the incumbent beyond an absolute or relative tolerance (usually
well before the queue is empty); it reports the bound and the gap also
when a node or time limit stops it. On request, the `Solution` also holds
the leaves of the tree explored, and a search can start from any set of
nodes instead of the root, which allows one to reoptimize after a change of
the data that leaves the branching valid, or to resume a search stopped by
a limit.

The search is sequential (`branch_and_bound( p , options )`), parallel with a
queue shared by the workers as in SPA* (`SharedBranchAndBound`) or with a queue
per worker and the children sent to their owners as in HDA*
(`DistributedBranchAndBound`), both on the three runtimes, master-worker on a
FastFlow farm (`FarmBranchAndBound`), with one MPI process per worker
(`MPIBranchAndBoundSearch`, whose nodes travel as bytes), or by synchronous
rounds (`DeterministicBranchAndBound`), whose outcome depends neither on the
timing of the workers nor on their number, as the deterministic mode of
FiberSCIP [[18]](#ref-18), so that the experiments of a paper can be repeated
exactly, or by the SelfSplit of Fischetti, Monaci and Salvagnin [[19]](#ref-19)
(`SelfSplitBranchAndBound`), whose workers search apart, and without exchanging
anything, the subtrees of a first sampling of the tree, the outcome depending
then only on the number of the workers, or by work stealing [[20]](#ref-20)
(`StealingBranchAndBound`), whose workers keep the children they generate and
take a node from the queue of another worker only when theirs is empty, so that
each worker stays in its own part of the tree, which pays when moving the
evaluator of a worker to a distant node is costly. The searches that need no
external library can also be chosen by a value, or by its name read from a file
of parameters (`parallel_branch_and_bound( p , bb_scheme( "selfsplit" ) ,
options )`), so that a program built on the library gets every new scheme with
no change of its own. In the taxonomy of Gendron and Crainic [[16]](#ref-16),
all of them build one tree in parallel (their type 2), asynchronously but for
the rounds, with a single pool of open nodes (the shared queue, the farm and
the rounds) or with a pool per worker (the distributed queues, the work
stealing, the SelfSplit and MPI); the parallelism inside a node (their type 1)
is left to the evaluation of the problem, and the search of several trees at
once (their type 3) is not provided. SafePBNF and PA*SE have no counterpart,
since they exist to handle duplicates and re-expansions, which a tree does not
have; for the same reason, the owner of a node only balances the load. Its
design is discussed in [docs/branch-and-bound.md](docs/branch-and-bound.md),
and the tests check the variants on the 0-1 knapsack against the dynamic
program.

## The data structures

OPEN is either a 4-ary heap, for any cost type, or a two-level bucket queue
(on f first, then on g) for integer costs, which [[4]](#ref-4) reports to be somewhat
faster than a heap on the sliding-tile puzzles. Neither one supports
decrease-key: an entry made stale by a cheaper path is skipped when it
surfaces. CLOSED (`NodeTable.h`) is an open-addressing hash table that
stores the hash of a state next to its node, so that most probes do not
compare the states, and spreads the slots by Fibonacci hashing; the nodes
are allocated from an arena. Both the table and the arena are requested
from the operating system as transparent huge pages, where available. In
the sequential search the slots of the successors are also prefetched
before being probed, since on large instances most probes are cache misses
(the parallel searches do the same with the batches they receive). Each
search returns a `Solution` holding the status (solved, unsolvable or
stopped at the node limit), the cost, the path and the counters of the
search; for the parallel ones, these include the load balance and the
communication overhead of [[2]](#ref-2).

These structures can be replaced, at compile time, by those of the
standard library (`std::deque`, `std::unordered_map` and
`std::priority_queue`), of Boost (`boost::object_pool`,
`boost::unordered_flat_map` and `boost::heap::d_ary_heap`) or of
[ADS](https://github.com/dmeoli/ADS) (`Arena`, `PointerTable`, `DaryHeap`
and `BucketQueue`), with the same interface (`Structures.h`); the
searches do not change, and neither does the order of their expansions,
but for the ties of f and g in OPEN, which the sets break differently.
Neither the standard library nor Boost has a bucket queue, and hence with
them both OPEN lists are a heap. The domains on a graph take its adjacency
from the set too: the grids from the `GridGraph` of ADS, which computes the
moves from a cell and their costs from a vector of the costs of the cells,
a move paying the cell it leaves as in the domain (the other sets keep
the adjacency of the library, as neither the standard library nor Boost
has a graph of a grid with blocked cells and diagonal moves, and a graph in
compressed rows of a map of 20000 by 20000 cells would take tens of GB),
and the explicit graphs from the `ArrayGraph` of ADS or the
`compressed_sparse_row_graph` of Boost, built from the arcs before the
first expansion. The sets exist to measure the searches on each of them
(`tools/bench/structures.sh`), and the B&B, whose queue is its own, does
not depend on them.

## Pattern databases

`PatternDatabase.h` gives the sliding-tile puzzles the additive pattern
databases of Korf and Felner [[22]](#ref-22) and of Felner, Korf and Hanan
[[23]](#ref-23) in place of the Manhattan distance. A table holds, for each
placement of the tiles of a pattern (ranked as a k-permutation of the
cells), its distance from the goal in an abstraction that keeps only the
tiles of the pattern, with the blank either ignored (the table is then
consistent) or minimized over the cells it can reach for free (a stronger
but inconsistent table, as in [[22]](#ref-22)); it is built by a
breadth-first search by layers on the workers of any runtime. `AdditivePDB`
sums the tables of a partition of the tiles and takes the maximum with the
reflection about the main diagonal, and `SlidingTilePDB` is `SlidingTile`
with that heuristic, so that every search of the library runs on it
unchanged:

```cpp
#include "PBeFS/PBeFS.h"
using namespace PBeFS_di_unipi_it;

PDBOptions o;
o.blank = BlankRule::minimized;
o.cache = "pdb-cache";                  // built once, then loaded
const AdditivePDB< 4 , 4 > pdb( partition_15_7_8() ,
                                FifteenPuzzle::korf_goal() , o );
const SlidingTilePDB< 4 , 4 > d( FifteenPuzzle::parse( "14 13 15 7 11 12 "
                                 "9 5 6 0 2 1 4 8 10 3" ) , pdb );
const auto s = AStar< SlidingTilePDB< 4 , 4 > >( d ).solve();  // cost 57
```

The entries are stored in a byte, in 4 bits as the excess over the
Manhattan distance of the tiles of the pattern (the default), modulo 3 in
1.6 bits [[25]](#ref-25) or in 1 bit [[27]](#ref-27) (the heuristic of a
successor then follows from that of its parent), by the least excess of
groups of consecutive ranks (the max-compression of [[24]](#ref-24)), or,
with the structures of ADS, by rank and select over the nonzero excesses;
the tables are saved to files with a versioned header and a checksum, and a
cache directory keeps them. The theory (the admissibility of the sum and of
the reflection, the consistency, the recovery of the exact value from the
residues) is on the page "Pattern databases" of the documentation. The
nodes expanded by A* on the 100 instances of Korf, summed over the
instances, which do not depend on the machine, are the following
(`tools/bench/pdb.sh`; all the costs are the optimal ones):

| heuristic | blank | entries | memory | expanded | ratio |
|---|---|---|---|---|---|
| Manhattan | | | | 1557138490 | 1 |
| 5-5-5, max-compressed by 16 | ignored | min | 49 KB | 681766045 | 2.3 |
| 5-5-5, no reflection | ignored | 1 bit | 197 KB | 369972440 | 4.2 |
| 5-5-5 | ignored | 4 bits | 786 KB | 263234565 | 5.9 |
| 6-6-3 | ignored | 4 bits | 5.8 MB | 167667391 | 9.3 |
| 7-8, no reflection | ignored | 1 bit | 72 MB | 137386424 | 11.3 |
| 7-8, max-compressed by 8 | ignored | min | 36 MB | 119574805 | 13.0 |
| 7-8 | ignored | 4 bits | 288 MB | 43587167 | 35.7 |
| 5-5-5 | minimized | 4 bits | 786 KB | 59968561 | 26.0 |
| 7-8 | minimized | 4 bits | 288 MB | 2019506 | 771 |

The rule of the blank weighs more than the size of the patterns: with the
blank minimized the partition 7-8 expands about 20000 nodes per instance,
i.e., 21 times fewer than with the blank ignored, at the price of about
6000 reexpansions over the 100 instances, while the encodings in 1 bit and
modulo 3 need the blank ignored and no reflection. The tables of the
24-puzzle (6-6-6-6, 255 MB in 4 bits) are built and looked up as well, and
`SlidingTile` packs the boards of more than 16 cells, 5 bits per cell, into
an integer of 128 bits (`Word128`), so that `SlidingTilePDB< 5 , 5 >` with
IDA* is the search of Korf and Felner [[22]](#ref-22) on that puzzle;
`data/24-puzzle/instances50.txt` holds the 50 instances distributed with
the code of HDA* of Jinnai and Fukunaga.

## Building and using the library

The library needs only a C++20 compiler, together with CMake 3.21 or GNU
make. OpenMP, MPI and FastFlow are optional, and each of them is used if it
is found. FastFlow is header-only: CMake looks for it in `FASTFLOW_ROOT`
(`/opt/fastflow` by default) and on the include path, and if it is not
there downloads the commit `PBEFS_FASTFLOW_COMMIT` (unless
`PBEFS_FETCH_FASTFLOW` is `OFF`), while the makefiles take its headers from
`FF_ROOT` (`/opt/fastflow` by default on Linux and macOS). With CMake one
has

```sh
git clone https://github.com/dmeoli/PBeFS.git
cmake -S PBeFS -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
cmake --install build --prefix <prefix>
```

where `PBEFS_WITH_OPENMP`, `PBEFS_WITH_MPI` and `PBEFS_WITH_FASTFLOW` (all
`ON` by default) can exclude a runtime, `PBEFS_STRUCTURES` chooses the data
structures of A* (`own`, the default, `std`, `boost`, which needs Boost
1.81 or later, or `ads`, which CMake finds by `find_package(ADS)` or else
downloads at the commit `PBEFS_ADS_COMMIT`; the makefiles take it from
`PBeFSSTRUCT`), and `PBEFS_BUILD_TESTS`,
`PBEFS_BUILD_EXAMPLES` and `PBEFS_BUILD_TOOLS` (`ON` only when PBeFS
is the top-level project) the tests, the examples and the benchmark. A CMake project then
uses the library with

```cmake
find_package(PBeFS REQUIRED)
target_link_libraries(myTarget PRIVATE PBeFS::PBeFS)
```

adding `PBeFS::OpenMP`, `PBeFS::MPI` or `PBeFS::FastFlow` for the
corresponding runtime; these are the components of the package, each
present if the library was built with its runtime and the runtime is found
again by the project that uses it; hence, `find_package(PBeFS REQUIRED
COMPONENTS OpenMP)` fails without OpenMP, while `find_package(PBeFS)`
finds the package with no `PBeFS::OpenMP`. Without CMake, `pkg-config --cflags --libs
pbefs` gives the flags of the library alone. The packages of vcpkg,
conda-forge, Homebrew, Debian and Ubuntu (`libpbefs-dev`), Conan, Spack,
Nix, the AUR, MSYS2, xmake and Meson WrapDB install the same files from the
recipes of [packaging](packaging/README.md), with OpenMP where the compiler
has it, MPI only where the channel lets one choose it, and without FastFlow,
which no channel packages. Without installing anything, a CMake project can
also take PBeFS at configure time, with FetchContent (CMake 3.24 or later),

```cmake
include(FetchContent)
FetchContent_Declare(PBeFS
    URL https://github.com/dmeoli/PBeFS/releases/download/v0.1.0/PBeFS-0.1.0.tar.gz
    FIND_PACKAGE_ARGS 0.1)               # an installed PBeFS, if there is one
FetchContent_MakeAvailable(PBeFS)
target_link_libraries(<target> PRIVATE PBeFS::PBeFS)
```

or with [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake),

```cmake
CPMAddPackage("gh:dmeoli/PBeFS@0.1.0")
target_link_libraries(<target> PRIVATE PBeFS::PBeFS)
```

and as a subproject PBeFS builds neither its tests nor its examples nor its
tools. With Meson, `meson wrap install pbefs` and `dependency('pbefs')`,
from WrapDB, in C++20. Since the searches are templates over the
`Domain`, they are in the headers, and the library (`src`) holds only what
does not depend on it, i.e., the grid maps, the random tables of Zobrist
hashing, the measures of a search and the parts that depend on the
operating system. Without CMake, the makefiles in the root give the library
in three ways to the makefile of a program that uses it, which needs only
to include one of them: `makefile` builds the objects and exports them
together with the include directive and the list of the headers,
`makefile-c` adds the external libraries of the runtimes (the threads,
OpenMP and FastFlow), and `makefile-inc` gives instead the library already
built, `lib/libPBeFS.a`, which `make -f makefile-lib` builds and `make -f
makefile-lib install PREFIX=<prefix>` installs. `makefile-OS` adapts them
to the operating system and to the compiler (`CC=g++` or `CC=clang++`),
and the places of the external libraries are those of
`extlib/makefile-default-paths-linux`, `-macos` and `-win` (FastFlow in
`FF_ROOT`, ADS in `ADSSDR`, Boost in `BOOST_ROOT` and, for Apple clang, the
OpenMP runtime of LLVM in `OMP_ROOT`), unless they are given on the command
line or in `extlib/makefile-paths`, which git ignores. In `examples`,
`test` and `tools` a makefile builds the corresponding programs, in `test`
the targets `check` and `check_mpi` also run the tests, and in the root
`make`, `make test`, `make test_mpi` and `make install` do the same for the
whole repository. Both CMake and the makefiles are tested on Linux, macOS
and Windows (with MSYS2, and CMake also with MSVC), with g++ and with
clang++; on Windows FastFlow, which relies on POSIX, and MPI are left out.

As an example, a minimal program that solves the instance 1 of Korf [[7]](#ref-7)
sequentially, and then with HDA* and Abstract Zobrist hashing on 8 threads,
reads

```cpp
#include "PBeFS/PBeFS.h"
using namespace PBeFS_di_unipi_it;

int main( void )
{
 const FifteenPuzzle d( FifteenPuzzle::parse(
                         "14 13 15 7 11 12 9 5 6 0 2 1 4 8 10 3" ) );
 const auto s = astar( d );
 ParallelOptions o;
 o.threads = 8;
 const AbstractZobristDistribution< FifteenPuzzle > az( d );
 const auto p = HDAStar< FifteenPuzzle , decltype( az ) >( d , az , o )
                 .solve();
 return( ( s.cost == 57 ) && ( p.cost == 57 ) ? 0 : 1 );
 }
```

The documentation is at
[dmeoli.github.io/PBeFS](https://dmeoli.github.io/PBeFS/): a
site of Sphinx with the theme of Read the Docs (light or dark), made of the
sections of this README, of the pages of the theory in `docs/theory` (A*
with its invariants and the cost of its operations, the parallel A*
searches, branch-and-bound and its parallel schemes, with what each one
guarantees), and of the reference of the classes and of the functions that
Breathe renders from the comments of the headers (through Doxygen), where
the page of each search gives its pseudocode, arranged by layer, by group
of headers (the searches, the domains, the runtimes, and so on), by header
and by class. Locally, with
Doxygen installed, `pip install -r docs/requirements.txt`, then
`python3 docs/pages.py && sphinx-build -b html docs build-docs/html`.

## The examples and the benchmark

The programs in `examples` are, for the most part, the search problems of
chapter 3 of Russell and Norvig [[13]](#ref-13), each one with its
`Domain` in a header of `include/PBeFS/astar` (e.g., `NQueens.h` or
`Romania.h`), which the benchmark shares. `8Puzzle` solves the instance of figure 3.28
of the book, or an instance at a given distance from the goal, with the
number of misplaced tiles (h1), the Manhattan distance (h2) or the sequence
score of Nilsson [[8]](#ref-8), which is not admissible; with `-t` it
repeats the comparison of h1 and h2 of figure 3.29 (nodes generated and
effective branching factor). `15Puzzle` solves an instance of Korf,
sequentially or with HDA* (`-t <threads>`); `RouteFinding` finds a road
between two cities of the map of Romania of the book; `FindPath` finds a
path on a weighted grid, random or read from a Moving AI map, and draws it.
`VacuumWorld` is the vacuum world of section 3.2.1 on a row of n cells, and
`NQueens` places n queens with the incremental formulation of the same
section, checking the placement apart from the search. The last three come
from the exercises of the book: `MissionariesCannibals` carries m
missionaries and c cannibals over a river with a boat of b seats (exercise
3.9), `WaterJugs` measures a quantity of water with jugs (exercise 3.6(d)),
and `TravelingSalesman` solves the traveling salesperson problem with the
minimum spanning tree heuristic (exercise 3.30), checking the tour against
the dynamic program of Held and Karp [[14]](#ref-14). Each run with a known
answer (e.g., the 26 moves of figure 3.28, the 418 km from Arad to Bucharest
and the 11 crossings of the missionaries) is also a test of `ctest`,
labelled `example`. In `tools`, `AStarBench` runs any search, runtime and
distribution of the library on the 100 instances of Korf
(`data/15-puzzle/korf100.txt`), on the 8-puzzle, on the grids (with unit
or weighted cells) and on the other problems of the examples, sized by
`--size` (e.g., the number of queens or of cities), and it prints a line of
CSV per run with the cost, the counters and the time;
`BranchAndBoundBench` does the same for every scheme of the
branch-and-bound, on random 0-1 knapsacks (`include/PBeFS/bnb/Knapsack.h`)
and on the traveling salesperson problem (the tree of the partial tours,
bounded by the minimum spanning tree), and `AStarBenchMPI` and
`BranchAndBoundBenchMPI` are their MPI counterparts, to be started with
`mpirun`. If the
environment variable `PBEFS_PIN` is set, the i-th worker is pinned to the
i-th CPU of the affinity mask of the process (as given, e.g., by `numactl
--physcpubind` or `taskset`). On machines with several sockets this is
usually needed for the timings to be repeatable. Finally,
`tools/bench/structures.sh` runs sequential A*, and HDA* and SPA* on 1 to
16 threads, of four builds, one per set of data structures, on instances of
growing size (the 15-puzzle in bands of the instances of Korf, and random
grids of growing side), and
`tools/bench/structures.py` makes the tables of the comparison, with the
noise of two runs of the same build. In both campaigns (`campaign.sh` and
`structures.sh`) each run is a process of its own, pinned with its workers
to given cores of one NUMA node and with its memory on that node, and the
runs that the tables compare are made twice, in two processes, so that the
gap between the two (the A/A noise, which `tables.py` and `structures.py`
report with each table) says how far two timings can differ by the machine
alone; a campaign holds a lock on each of its cores, so that a second one
started on the same cores waits for it, and `tools/bench/launch.sh` starts
it detached, saving its PID. On x86, `AStarBench` is compiled with the
branches kept off the 32-byte boundaries, by
`-Wa,-mbranches-within-32B-boundaries` where the assembler takes it,
against the erratum JCC of Skylake, which slows a loop by its alignment
alone, and since the data
structures compared with those of PBeFS are compiled into the same program
they get it alike. The timings are meant to be taken on a machine at rest;
a virtual machine whose memory is given and taken back while it runs (such
as a VM of Hyper-V with dynamic memory, which may stay at about 2.4 GB for
hours and swap on any large instance) gives timings that are not
comparable with each other, and is not to be used for them.

## Scaling with the number of cores

We are concerned here with how each parallel search of the library scales
with the number of cores on one and the same problem, for each problem of
the examples. The runs are those of table F of `tools/bench/campaign.sh`,
and the tables below are written from its output by `tools/bench/tables.py`.
Each problem is taken at a size where the sequential search takes about 1
to 10 s on one core. For the 15-puzzle these are the instances 1, 10 and 14
of Korf, and for pathfinding three random grids of $5000 \times 5000$ cells (25%
of them blocked, the others costing 1 to 5) with 4 and with 8 moves. Also,
the vacuum world has 22 cells and three starting cells, the queens are 28,
29 and 31, the missionaries and the cannibals are 65535 each (with boats of
100 and 200 seats), the jugs are four, of 193 to 211 units, and the
salesperson visits three sets of 32 random cities. However, the 8-puzzle (the instance of
figure 3.28 of the book and one of the hardest) and the road from Arad to
Bucharest cannot grow, and they are solved in microseconds to
milliseconds; hence, on them the tables mostly measure the cost of
starting and stopping the workers. The branch-and-bound is run on three
strongly correlated knapsacks of Pisinger with 90 items, which are usually
the hardest for the bound of Dantzig, and on three traveling salesperson
problems of 24 cities, on the tree of the partial tours bounded by the
minimum spanning tree.

A search is run with the options that its `Domain` supports. HDA* takes
each work distribution that the `Domain` gives (Z, AZ, A, P, H and R in the
notation of Jinnai and Fukunaga; the abstractions of AZ and A are given
only by the sliding-tile puzzles, the grids and the vacuum world) on the
threads, OpenMP and FastFlow. SPA* and SafePBNF (the latter where the
`Domain` has nblocks) run on the threads and on OpenMP, the farm of
FastFlow has a master and one worker less than the cores, and HDA* over MPI
has one process per core, with each distribution but R (which needs a
generator per worker). The branch-and-bound runs each of its schemes on
each runtime, its farm and its search over MPI. PA*SE pays only when the
expansions are slow (cf. the table of slow expansions); thus, it is run
with HDA*, SPA* and the farm on smaller instances whose expansions carry
the artificial work of `SlowDomain`. Each run that the tables compare is
made twice, in two processes. Furthermore, the cost it finds is checked
against the one of the sequential search, and a difference is written to
`F-mismatch.csv` and marked in the tables.

The timings were taken on MACHINE (CPU, memory, compiler and flags), on the
cores CORES of the NUMA node NODE, with the memory bound to that node, the
workers pinned one per core (`PBEFS_PIN`) and no other job on the machine.
The time of a run is the wall-clock time of `solve()`, i.e., from the
creation of the workers to the solution, the tracing of its path included;
it does not include the start of the process, the reading or the
generation of the instance, the construction of the `Domain` and of the
tables of the distribution, nor the release of the memory of the search.
In each table, a row is a search (the algorithm, the distribution of HDA*
in brackets and the runtime, i.e., thread, openmp, fastflow or mpi), and a
column is the number of cores given to it. A cell is the median over the
instances of the speedup over the sequential search on 1 core, i.e., A*
with the OPEN list on buckets or the branch-and-bound with the best-first
selection. The time of a run is the median of its two repetitions, and the
smallest and the largest speedup over the instances follow in
parentheses. A cell is empty when some run did
not solve its instance (e.g., HDA*[R] beyond its limit of nodes), and `!`
marks a cell with a run whose cost differs from the sequential one. The
A/A noise under each table is the median and the largest relative gap
between the two repetitions of a run.

RESULTS: the tables of `tables.py` (section "Scaling on the problems of
the examples"), one per problem, followed by the comment table by table,
i.e., which search scales best on which problem and by how much (as a
range), where the ranking changes with the number of cores, and the
exceptions.

## Testing

The tester `test/test.cpp` checks every search against references that do
not come from the library, i.e., a breadth-first search of the whole state
space of the 8-puzzle, the optimal costs of the 100 instances of Korf
[[7]](#ref-7), Dijkstra's algorithm on random grids and the distances on the
map of Romania; for each solution it checks both the cost and the path,
which must start at the initial state, end at a goal and be made of edges
of the `Domain` whose costs add up to the one reported. Its cases (the
sequential A*, then each parallel search on each runtime with 1 to 8
workers, and HDA* over MPI with 1, 2 and 4 processes) are listed in
`test/tests.txt`, which both CMake and the makefiles read, and the case
`registry` checks that the list and the cases of the tester are the same;
`ctest` runs them together with the examples of known answer. The CI runs
them on Linux with gcc and clang, on macOS and on Windows (without
FastFlow and MPI there), with the makefiles, and once more under
AddressSanitizer with UndefinedBehaviorSanitizer and under ThreadSanitizer,
which checks the synchronization of the searches in shared memory; with
each of the other sets of data structures of A* it runs them in Release
and under AddressSanitizer with UndefinedBehaviorSanitizer (with the
preconditions of ADS checked for its set), and on the three systems it builds
`packaging/test` against the installed package. The line coverage of the
tests is on Coveralls.

## Citing

If PBeFS is useful in a publication, one can cite it as

```bibtex
@software{PBeFS,
  author = {Meoli, Donato},
  title  = {{PBeFS}: Parallel Best-First Search},
  url    = {https://github.com/dmeoli/PBeFS},
  year   = {2026}
}
```

## License

This software is released under the [GNU Lesser General Public License
version 3.0](https://opensource.org/licenses/lgpl-3.0.html) (see the file
[LICENSE](LICENSE)). Its sequential version was first written for the
course of Knowledge Engineering of prof.
[Floriana Esposito](https://www.uniba.it/it/docenti/esposito-floriana) at
the University of Bari, and the parallel ones for the course of Parallel and
Distributed Systems of prof.
[Marco Danelutto](https://www.unipi.it/en/about-us/organisation/people/marco-danelutto-4278-en/)
and prof.
[Massimo Torquati](https://www.unipi.it/en/about-us/organisation/people/massimo-torquati-78079-en/)
at the University of Pisa.

## References

<a id="ref-1"></a>[1] A. Fukunaga, A. Botea, Y. Jinnai and A. Kishimoto, [*A Survey of Parallel A\**](https://arxiv.org/abs/1708.05296), arXiv:1708.05296, 2017.

<a id="ref-2"></a>[2] A. Kishimoto, A. Fukunaga and A. Botea, [*Evaluation of a simple, scalable, parallel best-first search strategy*](https://doi.org/10.1016/j.artint.2012.10.007), Artificial Intelligence 195, 222-248, 2013.

<a id="ref-3"></a>[3] E. Burns, S. Lemons, W. Ruml and R. Zhou, [*Best-first heuristic search for multicore machines*](https://doi.org/10.1613/jair.3094), Journal of Artificial Intelligence Research 39, 689-743, 2010.

<a id="ref-4"></a>[4] E. Burns, M. Hatem, M. J. Leighton and W. Ruml, [*Implementing fast heuristic search code*](https://doi.org/10.1609/socs.v3i1.18245), Proceedings of the Symposium on Combinatorial Search (SoCS) 3, 25-32, 2012.

<a id="ref-5"></a>[5] Y. Jinnai and A. Fukunaga, [*On hash-based work distribution methods for parallel best-first search*](https://doi.org/10.1613/jair.5225), Journal of Artificial Intelligence Research 60, 491-548, 2017.

<a id="ref-6"></a>[6] M. Phillips, M. Likhachev and S. Koenig, [*PA\*SE: Parallel A\* for slow expansions*](https://doi.org/10.1609/icaps.v24i1.13652), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 24, 208-216, 2014.

<a id="ref-7"></a>[7] R. E. Korf, [*Depth-first iterative-deepening: an optimal admissible tree search*](https://doi.org/10.1016/0004-3702(85)90084-0), Artificial Intelligence 27, 97-109, 1985.

<a id="ref-8"></a>[8] N. J. Nilsson, *Principles of Artificial Intelligence*, Tioga, 1980.

<a id="ref-9"></a>[9] P. E. Hart, N. J. Nilsson and B. Raphael, [*A formal basis for the heuristic determination of minimum cost paths*](https://doi.org/10.1109/TSSC.1968.300136), IEEE Transactions on Systems Science and Cybernetics 4(2), 100-107, 1968.

<a id="ref-10"></a>[10] S. Edelkamp and S. Schroedl, [*Heuristic Search: Theory and Applications*](https://doi.org/10.1016/C2009-0-16511-X), Morgan Kaufmann, 2012.

<a id="ref-11"></a>[11] F. Mattern, [*Algorithms for distributed termination detection*](https://doi.org/10.1007/BF01782776), Distributed Computing 2(3), 161-175, 1987.

<a id="ref-12"></a>[12] M. Aldinucci, M. Danelutto, P. Kilpatrick and M. Torquati, [*FastFlow: high-level and efficient streaming on multicore*](https://doi.org/10.1002/9781119332015.ch13), in Programming Multi-core and Many-core Computing Systems, Wiley, 261-280, 2017.

<a id="ref-13"></a>[13] S. Russell and P. Norvig, [*Artificial Intelligence: A Modern Approach*](https://aima.cs.berkeley.edu/), 3rd edition, Prentice Hall, 2010.

<a id="ref-14"></a>[14] M. Held and R. M. Karp, [*A dynamic programming approach to sequencing problems*](https://doi.org/10.1137/0110015), Journal of the Society for Industrial and Applied Mathematics 10(1), 196-210, 1962.

<a id="ref-15"></a>[15] N. R. Sturtevant, [*Benchmarks for grid-based pathfinding*](https://doi.org/10.1109/TCIAIG.2012.2197681), IEEE Transactions on Computational Intelligence and AI in Games 4(2), 144-148, 2012.

<a id="ref-16"></a>[16] B. Gendron and T. G. Crainic, [*Parallel branch-and-bound algorithms: survey and synthesis*](https://doi.org/10.1287/opre.42.6.1042), Operations Research 42(6), 1042-1066, 1994.

<a id="ref-17"></a>[17] D. S. Nau, V. Kumar and L. Kanal, [*General branch and bound, and its relation to A\* and AO\**](https://doi.org/10.1016/0004-3702(84)90004-3), Artificial Intelligence 23(1), 29-58, 1984.

<a id="ref-18"></a>[18] Y. Shinano, S. Heinz, S. Vigerske and M. Winkler, [*FiberSCIP: a shared memory parallelization of SCIP*](https://doi.org/10.1287/ijoc.2017.0762), INFORMS Journal on Computing 30(1), 11-30, 2018.

<a id="ref-19"></a>[19] M. Fischetti, M. Monaci and D. Salvagnin, [*Self-splitting of workload in parallel computation*](https://doi.org/10.1007/978-3-319-07046-9_28), Integration of AI and OR Techniques in Constraint Programming (CPAIOR 2014), Lecture Notes in Computer Science 8451, 394-404, 2014.

<a id="ref-20"></a>[20] R. D. Blumofe and C. E. Leiserson, [*Scheduling multithreaded computations by work stealing*](https://doi.org/10.1145/324133.324234), Journal of the ACM 46(5), 720-748, 1999.

<a id="ref-21"></a>[21] J. C. Culberson and J. Schaeffer, [*Pattern databases*](https://doi.org/10.1111/0824-7935.00065), Computational Intelligence 14(3), 318-334, 1998.

<a id="ref-22"></a>[22] R. E. Korf and A. Felner, [*Disjoint pattern database heuristics*](https://doi.org/10.1016/S0004-3702(01)00092-3), Artificial Intelligence 134(1-2), 9-22, 2002.

<a id="ref-23"></a>[23] A. Felner, R. E. Korf and S. Hanan, [*Additive pattern database heuristics*](https://doi.org/10.1613/jair.1480), Journal of Artificial Intelligence Research 22, 279-318, 2004.

<a id="ref-24"></a>[24] A. Felner, R. E. Korf, R. Meshulam and R. C. Holte, [*Compressed pattern databases*](https://doi.org/10.1613/jair.2241), Journal of Artificial Intelligence Research 30, 213-247, 2007.

<a id="ref-25"></a>[25] T. Breyer and R. E. Korf, [*1.6-bit pattern databases*](https://doi.org/10.1609/aaai.v24i1.7558), Proceedings of the AAAI Conference on Artificial Intelligence 24, 39-44, 2010.

<a id="ref-26"></a>[26] R. E. Korf, [*Minimizing disk I/O in two-bit breadth-first search*](https://mlanthology.org/aaai/2008/korf2008aaai-minimizing), Proceedings of the AAAI Conference on Artificial Intelligence 23, 317-324, 2008.

<a id="ref-27"></a>[27] R. Clausecker and A. Reinefeld, [*Zero-aware pattern databases with 1-bit compression for sliding tile puzzles*](https://doi.org/10.1609/socs.v10i1.18500), Proceedings of the Symposium on Combinatorial Search (SoCS) 10, 35-43, 2019.

<a id="ref-28"></a>[28] F. Yang, J. Culberson, R. Holte, U. Zahavi and A. Felner, [*A general theory of additive state space abstractions*](https://doi.org/10.1613/jair.2486), Journal of Artificial Intelligence Research 32, 631-662, 2008.

<a id="ref-29"></a>[29] U. Zahavi, A. Felner, R. C. Holte and J. Schaeffer, [*Duality in permutation state spaces and the dual search algorithm*](https://doi.org/10.1016/j.artint.2007.10.019), Artificial Intelligence 172(4-5), 514-540, 2008.

<a id="ref-30"></a>[30] R. C. Holte, A. Felner, J. Newton, R. Meshulam and D. Furcy, [*Maximizing over multiple pattern databases speeds up heuristic search*](https://doi.org/10.1016/j.artint.2006.09.002), Artificial Intelligence 170(16-17), 1123-1136, 2006.

<a id="ref-31"></a>[31] A. Felner, U. Zahavi, R. Holte, J. Schaeffer, N. Sturtevant and Z. Zhang, [*Inconsistent heuristics in theory and practice*](https://doi.org/10.1016/j.artint.2011.02.001), Artificial Intelligence 175(9-10), 1570-1603, 2011.

<a id="ref-32"></a>[32] T. Nagahashi and D. Takahashi, [*Construction of large zero-aware pattern databases for sliding puzzles on distributed memory machines*](https://doi.org/10.1007/978-3-031-97000-9_17), Computational Science and Its Applications (ICCSA 2025), Lecture Notes in Computer Science, 272-284, 2025.

<a id="ref-33"></a>[33] M. Samadi, M. Siabani, A. Felner and R. Holte, [*Compressing pattern databases with learning*](https://doi.org/10.3233/978-1-58603-891-5-495), Proceedings of the European Conference on Artificial Intelligence (ECAI 2008), Frontiers in Artificial Intelligence and Applications 178, 2008.

<a id="ref-34"></a>[34] W. Myrvold and F. Ruskey, [*Ranking and unranking permutations in linear time*](https://doi.org/10.1016/S0020-0190(01)00141-7), Information Processing Letters 79(6), 281-284, 2001.
