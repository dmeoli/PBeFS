# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the sequential A* over any Domain, i.e., an initial state, a goal
  predicate and a procedure that generates the successors of a state with
  the cost of their edges and their heuristic estimates, with the OPEN list
  on two-level buckets for integer costs (whose layers become heaps by g
  once a g reaches 256, so that a layer costs its entries and not the range
  of g) and on a 4-ary heap otherwise, the closed list in an open
  addressing table with the hash stored next to the pointer, the nodes in
  an arena on huge pages, the slots of the successors prefetched before
  their lookups, and the generation of the parent pruned for the domains
  that give the operator of a move
- the parallel searches of the survey of Fukunaga et al. (2017): SPA*, HDA*
  with the work distributions of Zobrist, abstract Zobrist, abstraction,
  perfect hashing, plain hashing and random hashing, SafePBNF, PA*SE and a
  master-worker farm with feedback, on C++ threads, OpenMP, FastFlow and
  MPI (HDAStarMPI, with the termination of Mattern)
- the domains of the sliding-tile puzzles (with the Manhattan and the
  pairwise heuristics), of the grids with 4 and 8 moves (random maps and the
  maps of Moving AI), of explicit graphs, and a domain with artificial work
  per expansion
- the best-first search with an incumbent, i.e., branch-and-bound, of a
  problem satisfying the concept BranchAndBound (root(), evaluate( n ) and
  branch( n , sink ), possibly a worker of its own per thread): the nodes
  evaluated when they are extracted, pruned when they are generated, when
  they are extracted and after their evaluation, the selections
  best-first, depth-first, breadth-first, best-estimate and plunging, the
  end when no open node can improve the incumbent beyond an absolute or
  relative tolerance,
  minimization and maximization, an initial incumbent, a callback on each
  new one, and node and time limits; sequential (BranchAndBoundSearch) and
  parallel with a shared queue (SharedBranchAndBound) or with distributed
  queues (DistributedBranchAndBound) on C++ threads, OpenMP and FastFlow,
  master-worker on a farm of FastFlow (FarmBranchAndBound), and with one
  MPI process per worker (MPIBranchAndBoundSearch), tested on
  the 0-1 knapsack against the dynamic program; the design is in
  docs/branch-and-bound.md
- the values of the branch-and-bound (bounds, feasible values, incumbent and
  absolute tolerance) of the value_type of the problem, double if it gives
  none: any signed type with a total order and its BBValueTraits (given for
  the arithmetic types, to be specialized for others, e.g., the rationals
  of an exact arithmetic), with BasicBBOptions, BasicParallelBBOptions,
  BasicEvaluation and BBSolution< Node , V >, of which BBOptions,
  ParallelBBOptions, Evaluation and BBSolution< Node > are the case of
  double; the incumbent shared under a lock for a type without a lock-free
  std::atomic, and the values sent by MPI as their bytes. BasicKnapsack
  rounds its bound down for an integer type, which prunes the nodes whose
  bound rounds to the incumbent (--values long in BranchAndBoundBench), and
  the case bb-values checks every scheme with integer and rational values
- the deterministic parallel branch-and-bound (DeterministicBranchAndBound),
  by synchronous rounds of the `round` best open nodes, evaluated by the
  workers against the incumbent of the start of the round and merged in
  their order, whose outcome depends neither on the timing of the workers
  nor on their number, on the three runtimes; the tests check that it is
  the same with 1 to 8 workers and at every repetition
- StealingBranchAndBound, the parallel branch-and-bound by work stealing:
  each worker keeps the children it generates in its own queue, and takes
  the best node of the queue of another worker only when its own is empty,
  so that it stays in its part of the tree, on the three runtimes; its
  scheme is named stealing (BBScheme::stealing), and the tests check it as
  the shared and the distributed searches, by cutoff and by its name
- the data structures of A* chosen at compile time (PBEFS_STRUCTURES in
  CMake, PBeFSSTRUCT in the makefiles, Structures.h): those of the
  library, or those of the standard library, of Boost or of ADS, with the
  same interface; the CI runs the tests with each of them, and
  tools/bench/structures.sh with structures.py compares A* and HDA* on
  them, on instances of growing size
- the choice of the branch-and-bound by a value (BranchAndBoundScheme.h):
  BBScheme names the searches of the core, bb_scheme() reads it from its
  name and parallel_branch_and_bound() runs the search of the scheme, so
  that a client gets every new scheme with no change of its own
- SelfSplitBranchAndBound, the SelfSplit of Fischetti, Monaci and
  Salvagnin (2014): the tree is sampled sequentially until there are
  ParallelBBOptions::split open nodes per worker, which are then dealt to
  the workers in the order of the selection, and each worker searches its
  subtrees apart, with its own evaluator and incumbent, so that the
  outcome depends only on the split and on the number of the workers, on
  the three runtimes; BranchAndBoundSearch takes the worker whose
  evaluator it uses
- the incumbent given to the evaluator of a worker by cutoff( z ), if it
  has one, before each evaluation and branching, in every branch-and-bound
  (in the deterministic one, the incumbent of the start of the round), so
  that a relaxation can stop as soon as it cannot improve it
- the reoptimization of a branch-and-bound: with BBOptions::keep_frontier
  the Solution holds the leaves of the tree explored (pruned, feasible,
  infeasible, or left open by a limit), and solve( seeds ) starts any of
  the searches from a given set of nodes instead of the root; over MPI each
  process keeps its own leaves and starts from its own seeds
- SearchOptions::expected_nodes, which sizes the closed list at the start
- the examples after Russell and Norvig, whose results are tests: the
  8-puzzle (with the heuristics h1, h2 and that of Nilsson), the 15-puzzle,
  pathfinding on a weighted grid, route finding in Romania, the vacuum
  world, the n queens (with the placement checked apart from the search),
  the missionaries and cannibals, the water jugs and the travelling
  salesman with the bound of the minimum spanning tree
- the tester, with the instances of Korf, the grids and the explicit
  graphs, under AddressSanitizer, UndefinedBehaviorSanitizer and
  ThreadSanitizer in the CI
- tools/AStarBench and tools/AStarBenchMPI, one line of CSV per run, and
  tools/bench/campaign.sh and tools/bench/tables.py, the runs and the
  tables of the benchmarks
- table F of tools/bench/campaign.sh, the scaling of every parallel search
  with the number of cores on every problem of the examples, each sized so
  that the sequential search takes about 1 to 10 s: HDA* with every
  distribution the Domain gives on the threads, OpenMP and FastFlow, SPA*,
  SafePBNF, the farm and HDAStarMPI, PA*SE with slow expansions, and on
  the knapsack and on the traveling salesperson problem every scheme of
  the branch-and-bound on every runtime, its farm and its MPI search, with
  the cost of every run checked against the sequential one; tables.py
  gives for each problem the median speedup over the sequential search and
  its range over the instances, and its A/A noise.
  tools/BranchAndBoundBench and tools/BranchAndBoundBenchMPI run the
  branch-and-bound with the columns of AStarBench, and AStarBench solves
  the problems of the examples, sized by --size, and the grids with
  weighted cells
- the build with CMake (installed package found by
  find_package(PBeFS)) and with the makefiles; the version from the
  git tag
- the documentation on the pages of the repository, a site of Sphinx with
  the theme of Read the Docs (and a switch to its dark version, the content
  1300 pixels wide, the sidebar fixed and scrolling on its own) written
  from the README by docs/pages.py (MyST), with the API reference that
  docs/apiref.py writes from the XML of Doxygen as a tree of layers (core,
  astar, bnb), groups of headers, headers and classes, rendered by
  Breathe; a warning of Doxygen fails the build of the site. The line
  coverage is on Coveralls, and the release workflow attaches to each
  release the source tarball PBeFS-X.Y.Z.tar.gz with its sha256
- the site of Doxygen next to that of Sphinx, in its directory doxygen:
  docs/doxygen/Doxyfile with doxygen-awesome-css 2.5.0 (MIT, in
  docs/doxygen/awesome) in its sidebar-only version, light or dark as the
  system is, the README as main page, the classes without the namespace,
  and Topics that hold the tree of the site of Sphinx (PBeFS.h, then the
  layers, their groups of headers, the headers with their classes and
  concepts), which docs/apiref.py writes into build-docs/modules.dox with
  no change to the headers; the design note goes under the layer bnb by
  the input filter docs/doxygen/markdown.py, which also sends the links of
  the Markdown files that Doxygen cannot follow to the site or to GitHub.
  The overloads in the outline of a page show the types of their
  parameters (docs/doxygen/outline-overloads.js), and custom.css puts the
  long types of the members, the titles of the groups of members and the
  entries of the outline on more lines. The pages of Sphinx link their
  headers and classes to the site of Doxygen, the headers list their
  concepts, the signatures of the dark theme are on a dark background, the
  cells of the tables wrap, the titles and the captions of the menu are
  in sentence case, and a warning of Doxygen or of Sphinx fails the build
- the pages of the theory in docs/theory, under the caption Theory of the
  site of Sphinx and in the Documentation of the layers astar and bnb of
  the site of Doxygen, each in a group of its own: best-first search and A*
  (the notation, admissibility and consistency, the invariant of OPEN and
  the optimality, the reopening, the ties, and the cost of each operation
  on the heap, the bucket queue and the hash table), the parallel A*
  searches (the overheads, the condition of termination, SPA*, the farm,
  HDA* with its distributions and the counters of its end, Safe PBNF and
  PA*SE, and a comparison), branch-and-bound as best-first search (the
  bounds, the pruning, the invariant and the reported bound, the selection
  rules and A* as a special case) and the parallel branch-and-bound (the
  anomalies, and the correctness, termination and determinism of each
  scheme), which link the design note and cite the references of the
  README by section; the formulae are written between dollars, which MyST
  renders and docs/doxygen/markdown.py turns into those of Doxygen, whose
  site has MathJax 3, and pages.py sends the links of the pages to the
  headers to their pages of the API reference and the citations of the
  design note to the references
- the pseudocode of every search, with numbered lines and the conventions
  of the page of the theory on A*, in the comment of its header (AStar.h,
  SPAStar.h, HDAStar.h, HDAStarMPI.h, SafePBNF.h, PASE.h, FarmAStar.h,
  BranchAndBound.h, ParallelBranchAndBound.h,
  DeterministicBranchAndBound.h, SelfSplitBranchAndBound.h,
  StealingBranchAndBound.h, FarmBranchAndBound.h, MPIBranchAndBound.h) and
  of the detection of the end in ThreadMailbox.h and MPIMailbox.h, each
  with the section of the paper it follows and with why it is correct, and
  of the dynamic program of Held and Karp in TravelingSalesman.h, with its
  recurrence and its cost; the formulae of the comments of the headers, of
  the README and of the design note are in LaTeX (\f$ \f$ in the comments,
  dollars in the Markdown), and in the site of Doxygen the mentions of the
  design note in the comments of the headers link its page
- the recipes of packaging/ for vcpkg, conda-forge, Homebrew, Debian and
  Ubuntu (libpbefs-dev, with deb.sh for the binary package and for a PPA),
  Conan Center, Spack, Nix (with flake.nix at the root), the AUR, MSYS2,
  xmake-repo and Meson WrapDB, all on the tarball of the release, with
  OpenMP where the compiler has it, MPI only where the channel lets one
  choose it, and without FastFlow, which no channel packages; bump.sh
  writes the version and the checksums of a release into all of them, and
  packaging/test is the project that each of them builds against the
  installed package; the workflow packaging.yml builds and tests the
  formula of Homebrew on macOS, the port of vcpkg and the package of MSYS2
  on Windows, the package of the AUR on Arch, the Debian package and the
  port of WrapDB, and the README gives the use with FetchContent and
  CPM.cmake
- the pkg-config file pbefs.pc, installed with the library
- tools/bench/launch.sh, which starts a campaign detached by nohup and
  setsid, with its PID (the group of its processes) saved
- the pattern databases of the sliding-tile puzzles (PatternDatabase.h):
  TilePattern ranks the placements of a pattern, PatternTable holds the
  distances of its abstraction with the blank ignored or minimized, built
  by a breadth-first search by layers in parallel on any runtime (two bits
  per state with the blank minimized), AdditivePDB sums the tables of a
  partition (7-8, 5-5-5 and 6-6-3 of the 15-puzzle, 6-6-6-6 of the
  24-puzzle) and takes the maximum with the reflection about the main
  diagonal, and SlidingTilePDB is SlidingTile with that heuristic; the
  entries in bytes, in 4 bits as the excess over the Manhattan distance,
  modulo 3 in 1.6 bits and in 1 bit (the heuristic of a successor from
  that of its parent), compressed by the least excess of groups of ranks,
  and, with the structures of ADS, by rank and select over the nonzero
  excesses; the tables saved to and loaded from files with a versioned
  header and a checksum, and kept in a cache directory; the tests pdb and
  pdb-search, tools/PDBBench and tools/bench/pdb.sh, and the theory in
  docs/theory/pattern-databases.md

### Changed

- the domains of the examples are headers of the library (NQueens,
  VacuumWorld, MissionariesCannibals, WaterJugs, TravelingSalesman and
  Romania in include/PBeFS/astar, the knapsack of the tests in
  include/PBeFS/bnb/Knapsack.h), so that the examples, the tests and the
  benchmarks share them; NQueens takes up to 32 queens and
  TravelingSalesman up to 32 cities (Held and Karp checks the example up
  to 20), the domains give the features of Zobrist hashing and, where the
  state fits a word, a perfect hash, VacuumWorld also the abstractions,
  the nblocks and the pairwise heuristic, and Grid and ExplicitGraph a
  perfect hash; GridMap::random draws the costs of the cells up to a given
  maximum, as FindPath does, and MissionariesCannibals no longer scans the
  crossings larger than the boat

- the licence is the GNU Lesser General Public License version 3.0,
  instead of the MIT License

- the library is compiled (src/ holds what does not depend on the Domain)
  and written in C++20, instead of the single header of C++11 with linear
  searches in the OPEN and closed lists
- FastFlow is found in /opt/fastflow or on the include path, and CMake
  downloads a fixed commit when it is not there, instead of being a git
  submodule
- the library is PBeFS (Parallel Best-First Search), instead of
  AStarSearch, since it holds the branch-and-bound beside A*: the CMake
  package, the library, the directory of the headers, the include guards,
  the namespace (PBeFS_di_unipi_it) and the macros (PBEFS_*) are renamed,
  and the headers are in three layers, core/ (the runtimes, the messages
  over MPI with the termination of Mattern written once, and the helpers),
  astar/ (the graph searches and their domains) and bnb/ (the tree
  searches), each with its header (core.h, astar.h, bnb.h) and PBeFS.h for
  all of them; the messages between the workers (ThreadMailbox, in shared
  memory, and MPIMailbox), with the detection of the end, the parameters of
  the parallel searches (ParallelParams) and the start of the master-worker
  farms are written once for A* and branch-and-bound
- the version is read first from the VERSION.txt of a tarball, which git
  archive fills in with the tag, and from git only in a clone whose top
  level is the source directory, so that a tarball unpacked inside another
  git repository (the build trees of vcpkg, a FetchContent of the tarball)
  does not take the tags of that repository
- the targets of the runtimes (PBeFS::OpenMP, PBeFS::MPI, PBeFS::FastFlow)
  are installed each in an import file of its own, read by the installed
  PBeFSConfig.cmake only if the runtime is found again, and are the
  components of the package: a runtime missing on the machine of the
  consumer removes only its target, and the package fails to be found
  without it only if the runtime is among its required components
- the campaigns of tools/bench run each repetition in a process of its own
  (sequential A* with buckets twice too, the A/A of its table) and hold a
  lock (flock) on each of their cores; tables.py takes the median of the
  runs of each instance in the table of sequential A* and gives its A/A,
  and structures.py the A/A of each row; AStarBench is compiled with the
  branches kept off the 32-byte boundaries (the erratum JCC of Skylake)
  where the assembler takes the option, by CMake and by the makefile
- the CI runs macOS and Windows on their latest runners, builds
  packaging/test against the installed package on Linux, macOS and
  Windows, and runs the tests of each of the other sets of data structures
  once more under AddressSanitizer and UndefinedBehaviorSanitizer, with the
  preconditions of ADS checked
- the tests cover the limit cases of every search, on every runtime and
  over MPI: an initial state that is the goal, one move, a node limit of 5,
  an unsolvable instance, more workers than states, and the parameters at
  0 or 1 for A*; an infeasible root, a root solved without branching, an
  initial incumbent better than the optimum, an absolute tolerance, a time
  limit, a node limit of 1, plunging at the two extremes, more workers
  than nodes, and the parameters at 0 or 1 for branch-and-bound

[Unreleased]: https://github.com/dmeoli/PBeFS/commits/main
