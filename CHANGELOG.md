# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- SequenceAlignment, the optimal alignment of up to 16 protein sequences as
  a shortest path in the lattice of their positions, with the sum-of-pairs
  cost 17 - S of the PAM250 scores S (or of a matrix read by read_scores),
  the heuristic of the optimal pairwise alignments of the suffixes, the
  features of Zobrist hashing (the position in each sequence, updated by
  the move), boxes of the lattice as the abstraction of AZH and AHDA* and
  as the nblocks of PBNF, and read_fasta for the files of BAliBASE;
  a linear gap penalty in place of the gap scores, random_sequences for
  the random instances of Jinnai and Fukunaga (JAIR 2017), level() and
  max_level(); AStarBench --domain msa, --matrix, --gap and --random, the
  test msa against the dynamic program on the whole lattice and against
  Needleman and Wunsch, and in data/msa the instances of BAliBASE of the
  benchmark sets of Jinnai and the 15 of BAliBASE 3.0 of the experiments
  of Kobayashi et al. (IJCAI 2011) and of Jinnai and Fukunaga
- the concept LevelDomain and HyperplaneDistribution, the hyperplane work
  distribution of Kobayashi, Kishimoto and Watanabe (IJCAI 2011) for the
  levelled graphs, with integer and fractional thickness and their rule
  for it; AStarBench --dist HWD and --plane, and the tests of its bound on
  the workers that the successors of a worker reach and of HDA* with it
- GridBoard, the boards of Burns et al. (unit or life costs, 4 or 8 moves,
  the corners cut), read with the row 0 at the bottom or at the top, and
  Grid built from a board: GridCost::life with its heuristic, the exact
  cost on the board without blocked cells, and the cutting of the corners
  in both adjacencies; AStarBench --domain board and --origin, the test
  board against Dijkstra, and representative boards in data/grid
- SlidingTile::read_positions for the instances given by the cell of each
  tile, TravelingSalesman::read_cities for the files of cities (Euclidean
  or Manhattan distances, the first n cities), AStarBench --metric and
  --file for tsp, the 100 and the 8 instances of the 24-puzzle of the
  benchmark sets of Jinnai, korf_1, korf_12 and korf_55, ten instances of
  40 cities, and the test files
- the tests of HDA* (Z, AZ, A), SafePBNF, SPA*, MQA*, the farm and PA*SE
  run also on alignments and on boards with life costs, and a check of A*
  against the expansions that every tie-breaking makes with a consistent
  heuristic (all the states of f* < C*, none of f* > C*)
- the bounded-suboptimal searches: a weight w >= 1 in SearchOptions (the
  class Weight, with exact fractions for integer costs) orders OPEN by
  g + w h in A*, IDA*, SPA*, MQA*, HDA* (threads and MPI), SafePBNF (then
  wPBNF) and FarmAStar, which return a path of cost at most w C* with the
  lower bound on C* they prove (Solution::bound); the parallel ones prune
  by w f and stop when no node of g + w h is below the incumbent, which
  keeps the bound in any order of expansion; drop_duplicates reopens a
  closed node only when the saving exceeds (w - 1)(h(u) - h(v)), which
  keeps the bound with any admissible heuristic; anytime makes the
  parallel searches go on until the incumbent is optimal, recording each
  path in Solution::improvements; with PBNFOptions::schedule SafePBNF is
  the Anytime Restarting PBNF of Burns et al. (2010); with no weight the
  searches run with UnitWeight, the weight 1 known at compile time, and
  their instructions per expansion stay those of before (on Korf #55 with
  one worker, within 3% for A*, HDA* and SafePBNF)
- AnytimeAStar, with Anytime Weighted A* (Hansen and Zhou, 2007) and ARA*
  (Likhachev, Gordon and Thrun, 2003) with its list INCONS; the optimistic
  search (Thayer and Ruml, 2008) and the Improved Optimistic Search (Chen
  et al., 2019); FocalSearch, with the focal search A*eps (Pearl and Kim,
  1982) and Explicit Estimation Search (Thayer and Ruml, 2011);
  PotentialSearch, the Dynamic Potential Search (Gilon, Felner and Stern,
  2016); LPAStar, the Lifelong Planning A* (Koenig, Likhachev and Furcy,
  2004), which repairs its search after the caller changes some costs
- LowGOpenList, a heap whose ties of f (or of the key) go to the smaller
  g, chosen by OpenKind::heap_low_g in astar(), anytime_astar() and
  optimistic_search() and by --open lowg in AStarBench: with a weight it
  breaks the ties of g + w h by the smaller f, which on the grids of Burns
  et al. expands 1.1 to 1.5 times fewer nodes than the larger g and finds
  cheaper paths, while the larger g stays better on the 15-puzzle
- the concepts DistanceDomain and PredecessorDomain, satisfied by
  SlidingTile and Grid (and the second by ExplicitGraph), Grid::refresh(),
  which reads again the costs of the map, and the keyed push( node , g ,
  key ) of the OPEN lists
- AStarBench runs the new searches (--algo awa, ara, optimistic, ios,
  focal, dps, lpa) with --weight, --dd, --anytime, --schedule, --optimism,
  --focal and --paths, and ends each line with the weight, the bound and
  the number of paths found
- ParallelIDAStar, the parallel IDA* below a frontier of units built by A*,
  by uniform-cost search or by the tree of fixed depth, each unit with its
  own threshold so that a unit whose subtree has no node of f below the
  bound sits the iteration out; the synchronous scheme (a barrier per
  iteration, the first goal found optimal) and the asynchronous one of
  AIDA* (the units by increasing threshold, until all reach the
  incumbent), the units scheduled by work stealing or by a shared cursor,
  longest first, the splitting of the stacks of Rao and Kumar, the lookups
  of the first levels below a unit in the table of the frontier, and the
  deterministic path of the first goal in the order of the units; on the
  three executors, and ParallelIDAStarMPI on the processes of MPI; --algo
  pida of AStarBench, the section of the theory, and the tests pida-thread,
  pida-openmp, pida-fastflow and pida-mpi
- the concept TreeOperatorDomain, the moves that a tree search skips
  besides the inverse of the last one (two turns of a face of Rubik's
  cube, two moves of a disk of the Towers of Hanoi, two commuting moves of
  TopSpin), which IDAStar and ParallelIDAStar use
- Permutation.h, the lexicographic ranking of the k-permutations by
  population counts, by the binary tree of Bonet, and the ranking of
  Myrvold and Ruskey; AbstractionPDB.h, the tables of any ranked
  abstraction with integer costs, built in parallel by the algorithm of
  Dial on the table, stored in 4 or 8 bits, saved and cached;
  PermutationPDB.h, the pattern databases of the permutation puzzles with
  orientations (PermutationDomain), their maximum or sum, and PDBDomain
- the domains PancakePuzzle and BurnedPancakePuzzle with the gap heuristic
  (and GAP-x), TopSpin, RubiksCube with the three tables of Korf (1997)
  and his ten instances in data/rubik, and the Towers of Hanoi with four
  pegs, any start and goal and additive tables of ranges of disks, with the
  20 instances of Siag et al. (2024) in data/hanoi; the tests permutation,
  abstraction, pancake, topspin, rubik, rubik-pdb and hanoi, and the
  domains of AStarBench
- DimacsGraph.h, the reader of the road graphs of the 9th DIMACS Challenge
  into an ExplicitGraph (which gets set_target()) with the consistent
  heuristic of the coordinates, a synthetic sample in data/dimacs,
  tools/bench/fetch_dimacs.sh for the real graphs, --domain road of
  AStarBench, and the test dimacs
- the lazy children of the B&B: an evaluator may give children( n ), a
  cursor that gives the children one at a time, in place of or besides
  branch( n , sink ), and child( n , k ), the child k without evaluating
  n; every search takes either, and prune_level ends the children of a
  node at the first one pruned (on brock200_1, 693899 children built
  instead of 7445660, with the same 524730 nodes evaluated)
- the modes of the B&B (BBMode): decide, which seeks a value within the
  target of the options and stops at the first one, and enumerate, which
  visits the whole tree and counts its feasible nodes in BBStats::
  solutions, both in every scheme through the test of pruning
- SkeletonBranchAndBound, the depth-first B&B of the search skeletons of
  Archibald et al. (2018, 2020): depth-bounded spawning, budget spawning,
  stack stealing (one child or all of a level) and ordered tasks by
  discrepancy, with the tasks as nodes or (indexed) as the paths of their
  choices, rebuilt by child( n , k ) or by evaluating the path
  (BBStats::recomputed); the schemes depthbounded, budget, stackstealing
  and ordered of parallel_branch_and_bound()
- racing_branch_and_bound(), the racing ramp-up of the Ubiquity
  Generator: the settings of BasicParallelBBOptions::race search the root
  at once with a budget of nodes, sharing their incumbents, and the scheme
  goes on from the open leaves of the best bound; parallel_branch_and_bound()
  runs it when race is not empty
- BBCheckpoint.h: the file of a checkpoint (open nodes, incumbent and its
  node, bound of the closed region, counters), written atomically, and
  checkpointed_branch_and_bound(), which runs any search in segments,
  writing a checkpoint after each, and resumes from the file of an earlier
  run; BBOptions::keep_open_only and BBSolution::closed_bound, which keep
  only the open leaves and the bound of the others; the serialization of
  the nodes moves here from MPIBranchAndBound.h
  (SerializableBranchAndBound)
- MPIHubBranchAndBoundSearch, the B&B over MPI of a master, hubs and
  workers as in ALPS: each process searches its own subtrees a unit of
  work at a time, and the hubs pair donors and receivers by the quantity
  and by the quality (the mean bound of the elite nodes) of their work,
  the master across the clusters (BBHubParams)
- IntervalBranchAndBound, the work stealing of intervals of leaves in a
  mixed radix (the Integer-Vector-Matrix of Mezmaz et al. and Gmys et
  al.) for the trees whose nodes of a depth have the same number of
  children (RadixBranchAndBound), each node evaluated once, so that with
  the optimum given the nodes evaluated do not depend on the steals
- the example domains of the B&B: UTS (with the trees T1, T2 and T3, of
  4130071, 4117769 and 4112897 nodes), MaxClique (bitset colouring bound,
  DIMACS reader, data/dimacs/brock200_1.clq), MaxCommonSubgraph (McSplit),
  SubgraphIsomorphism (forward checking, LAD reader), QAP (Gilmore-Lawler
  bound, QAPLIB reader, 15 instances in data/qaplib) and Flowshop
  (one-machine and two-machine bounds, forward, backward, alternate and
  dynamic branching, the instances of Taillard from their seeds)
- BranchAndBoundBench: the domains uts, clique, qap, flowshop, sip and
  mcs, the algorithms of the skeletons and ivm (and hub with MPI), the
  options of the race, of the modes, of the checkpoints and of the
  domains, and the column recomputed
- the tests bb-skel, bb-ivm (on the three runtimes), bb-modes, bb-race,
  bb-checkpoint, bb-domains and bb-hub-mpi, in test/bnb.cpp
- the concepts BidirectionalDomain (the goal state, the predecessors and a
  backward heuristic), BackwardOperatorDomain and MinCostDomain (the least
  cost of an edge), satisfied by SlidingTile, Grid and ExplicitGraph, whose
  arcs are kept reversed too, and the adapters ReverseDomain, for the
  reverse searches, and BlindDomain, without heuristics;
- the front-to-end bidirectional searches over one table of the nodes of
  both directions: the bidirectional Dijkstra's algorithm, MM, MMe, fMM( p )
  and BAE* (BiBestFirst.h), NBS (NBS.h), DVCBS (DVCBS.h), GBFHS and the
  bidirectional breadth-first search (GBFHS.h), with the bounds of each one
  or the largest of all the valid ones, and the page "Bidirectional search"
  of the documentation with their proofs;
- ParallelBidirectional, the two directions on two workers of an executor,
  PNBA* and the parallel versions of the bidirectional Dijkstra's
  algorithm, MM and BAE*;
- PEMSearch, the parallel external-memory searches PEM-A*, PEM-rA*, PEMM and
  PEM-BAE* with a cap on the disk, and IDDAStar, A* with immediate
  duplicate detection on the disk (A*-IDD), with the status disk_limit;
- the searches above in AStarBench and the test cases bidir, idd,
  bidir-thread, bidir-openmp, bidir-fastflow, pem-thread, pem-openmp and
  pem-fastflow.
- GreedySearch, the greedy best-first search and the search ordered by
  g + w h (w the weight of SearchOptions, infinite by default in
  GreedyOptions), with the OPEN lists of the satisficing planners as options that
  compose (GreedyOptions): eager or lazy evaluation, queues alternated by
  priority, each on a heuristic of a MultiHeuristicDomain with its ties
  broken by another, queues of the preferred successors of a
  PreferredDomain boosted at each new least value of a heuristic,
  type-based exploration, epsilon-greedy and Pareto queues, the ties broken
  FIFO, LIFO, by g or at random, reopening, a bound on g or g + h, and
  HeuristicCache, the evaluations kept across searches;
  EnforcedHillClimbing, with the preferred successors alone or first and a
  fallback on GreedySearch; IteratedSearch, a sequence of such searches
  each bounded by the best path before it, and lama_phases(); the concepts
  MultiHeuristicDomain, PreferredDomain, DeferredDomain, GoalFeatureDomain,
  RelevantFeatureDomain and PartialExpansionDomain, satisfied by
  SlidingTile (Manhattan distance and misplaced tiles, the moves that lower
  the Manhattan distance, the tiles in their goal cells, the operator
  selection function of EPEA*) and in part by Grid
- the width-based searches of WidthSearch.h on a FeatureDomain:
  NoveltyTable, IW(k), SIW and BFWS with the novelty computed for each
  value of the goal features to achieve (or of h) and of the relevant
  features of a relaxed plan
- ParallelGreedy, the parallel greedy searches with a shared OPEN list:
  KPGBFS, PUHF, PUHF2, PUHF3, PUHF4, SPUHF and OBAT, each with or without
  the separate generation and evaluation of the successors, on the three
  executors or in synchronous rounds in one thread, and GreedyTrace, what
  a search expanded and committed; HDGBFS, GBFS on a hash distribution;
  the test pgbfs-rounds computes the bench transition system of random
  graphs and of the 8-puzzle and checks that the constrained rules expand
  only its states, that OBAT commits in an order of GBFS and leaves at most
  k states of each h in DEFERRED, and the pathologies of KPGBFS and PUHF
- PEAStar and EPEAStar (PartialExpansion.h), FringeSearch, RBFS,
  BeamSearch, RealTimeSearch (LRTA* and LSS-LRTA*) and HPAStar on the
  grids, the options of all these searches in AStarBench, the evaluations
  as a column of its CSV (SearchStats::evaluated), and the tests
  gbfs, gbfs-queues, ehc, iterated, width, classic, realtime, hpa and
  pgbfs on each executor
- the searches by edges for slow evaluations: EPASE, the edge-based PA*SE
  of Mukherjee, Aine and Likhachev (SoCS 2022), w-ePA*SE with a weight,
  GePA*SE (SoCS 2023) with generalized, which evaluates the cheap edges of
  a state together, and the anytime A-ePA*SE and A-GePA*SE (SoCS 2023) with
  SearchOptions::anytime, which lower the weight as ARA* by weight_step,
  record each path in Solution::improvements and evaluate each edge at most
  once over the iterations; one worker selects and expands, waking the
  others one at a time, or with dispatcher off every worker does all;
  PWAStar, weighted A* whose workers evaluate the edges of one expansion;
  MPLP, the lazy search of Mukherjee, Aine and Likhachev (RA-L 2022), which
  checks a path when its last edge is evaluated and repeats the search only
  after a cost has changed; the concepts EdgeDomain,
  HeterogeneousEdgeDomain and LazyEdgeDomain; SlowEdgeDomain, which makes
  the edges of any domain cheap, expensive or missing; RobotNav2D, the 2D
  navigation of their experiments, with its maps and pairs of points in
  data/robot-nav-2d; the counter iterations of SearchStats, and evaluated
  counting the edges evaluated;
  the proofs on the page of the parallel searches; --algo epase, pwa and
  mplp, --domain nav2d, --edges and --step of AStarBench, whose lines end
  with the two counters; the tests edge and edge on the three runtimes
- wPA*SE: PASE with the weight of SearchOptions on the heuristic in the
  order of OPEN (with eps raised to it), as PWAStar and MPLP take it, and
  the states of OPEN examined for a safe one as an option (candidates, 64
  by default, 0 for all); the test pase on OpenMP and FastFlow
- the concept IncrementalFeatureDomain, whose for_each_feature_change( s ,
  t , f ) gives the features that the move from s to its successor t
  changes, satisfied by SlidingTile; the Zobrist and Abstract Zobrist
  distributions then give the key of a successor, the high 32 bits of its
  hash, which decide its owner (key_owner), from the key of its parent
  (successor_owners), and HDA* and HDAStarMPI keep the key in the padding
  of each node, so that they hash each successor in two lookups instead of
  15 on the 15-puzzle and the state that they expand in none: on Korf #1
  with one worker, with the exclusive or of SlidingTile below, 1369 to 1151
  instructions per expansion with Z and 1339 to 1143 with AZ
- DDDAStar, the A* with delayed duplicate detection on the disk, for the
  domains whose CLOSED would not fit in memory and whose states are
  reached by many paths: the External A* of Edelkamp, Jabbar and Schrödl
  (2004) with the hash-based detection of Korf (2008), i.e., buckets of
  equal g and h in files on the disk (RecordFile of ADS), split by the hash
  of the state, expanded by increasing f and then g, each file freed of
  its copies in a table in memory and of the states of the closed files of
  the same h, and the path traced back from the goal by the hash of the
  parent that each record carries; the concept DiskDomain (integer costs,
  trivially copyable states), DDDOptions (partitions, directory, buffers),
  --algo ddd of AStarBench, the section of the theory, and the test ddd

- the continuous integration checks, on a pull request, the lines that it
  changes in the style of the repository code-style (its variant spaces,
  which .code-format chooses, with no tabs and a space after a cast), whose
  action posts the changes that put them there as suggestions; the README
  says how to run its formatter before a commit
- MQAStar, the MultiQueue A*: SPA* with the OPEN list on the MultiQueue of
  ADS, a relaxed priority queue of 2p queues under a lock each, and the end
  of the search detected by a counter of the nodes in OPEN or being
  expanded, and with integer costs f and the complement of g packed in the
  key, so that the ties of f go to the deeper node as in AStar; AStarBench
  --algo mqa, and the tests mqa-thread and mqa-openmp
- ConcurrentClosed, the CLOSED list shared by SPAStar and MQAStar: with the
  data structures of ADS a set of pointers to the nodes in the
  ConcurrentFlatTable of ADS, read and filled with no lock, and 256 locks
  chosen by the address of a node, taken only when a state already there
  is reached with a smaller g (a duplicate takes none); with those of std
  and Boost the 256 shards of a NodeTable under a lock each that the two
  searches had before; SPAStar and MQAStar prefetch the places in CLOSED of
  all the successors of a node before looking any of them up, and MQAStar
  puts them into one queue of the MultiQueue under one lock and updates its
  counter of termination by one atomic step per expansion instead of two
- MultiQueueBranchAndBound, the parallel B&B on the MultiQueue of ADS, as
  MQAStar: the open nodes in 2p queues under a lock each, the children of
  a node pushed into one of them under one lock, and the end detected by a
  counter of the open nodes updated by one atomic step per node; the
  scheme multiqueue of parallel_branch_and_bound() (values atomic with no
  lock only), BranchAndBoundBench --algo multiqueue, and the tests
  bb-mq-thread and bb-mq-openmp
- PASE computes the hashes of the successors of a node without its lock
  and, under it, prefetches the places in CLOSED of all of them before any
  lookup, so that the lock is held for about one miss instead of one per
  successor

- the sequential A* over any Domain, i.e., an initial state, a goal
  predicate and a procedure that generates the successors of a state with
  the cost of their edges and their heuristic estimates, with the OPEN list
  on the two-level bucket queue of ADS for integer costs and on its 4-ary
  heap otherwise, the closed list in the open-addressing table of ADS with
  the hash stored next to the pointer, the nodes in its arena on huge
  pages, the slots of the successors prefetched before
  their lookups, and the generation of the parent pruned for the domains
  that give the operator of a move
- IDA*, the iterative-deepening A* of Korf (1985), for the domains whose
  CLOSED would not fit in memory: the depth-first search bounded by f and
  repeated with the least f that exceeded the bound, the parent pruned by
  its operator or by its state, and the expansions bounded by the node
  limit (IDAStar.h, idastar(), and --algo ida in AStarBench)
- the parallel searches of the survey of Fukunaga et al. (2017): SPA*, HDA*
  with the work distributions of Zobrist, abstract Zobrist, abstraction,
  perfect hashing, plain hashing and random hashing, SafePBNF, PA*SE and a
  master-worker farm with feedback, on the executors of ADS (its pool of
  C++ threads, OpenMP and FastFlow) with its channels, its mailbox that
  detects the end and its farm, and over MPI (HDAStarMPI, on the mailbox of
  the processes of ADS with the termination of Mattern)
- the domains of the sliding-tile puzzles up to the 24-puzzle (the boards
  of more than 16 cells packed, 5 bits per cell, into an integer of 128
  bits, Word128; the 50 instances of the 24-puzzle in data/24-puzzle),
  with the Manhattan and the pairwise heuristics, of the grids with 4 and 8
  moves (random maps and the
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
  queues (DistributedBranchAndBound) and master-worker on a farm with
  feedback (FarmBranchAndBound) on the executors of ADS, and with one
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
  each worker keeps the children it generates in its own queue, and steals
  about half of the queue of another worker only when its own is empty, so
  that it stays in its part of the tree, on the three runtimes; its
  scheme is named stealing (BBScheme::stealing), and the tests check it as
  the shared and the distributed searches, by cutoff and by its name
- the data structures of A* chosen at compile time (PBEFS_STRUCTURES in
  CMake, PBeFSSTRUCT in the makefiles, Structures.h): those of ADS (the
  default), of the standard library or of Boost, with the same interface,
  the adjacency of the grids and of the explicit graphs included (the
  GridGraph of ADS, which gives the costs of the moves too, and the
  ArrayGraph of ADS or the graph in compressed rows of Boost); the CI runs
  the tests with each of them, and tools/bench/structures.sh with
  structures.py compares on them A*, and HDA* and SPA* by the number of
  threads, on instances of growing size
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
  the three runtimes; the sampling ends with a dive to a leaf, whose
  incumbent the workers start from (ParallelBBOptions::split_dive), and
  with ParallelBBOptions::split_shared the workers prune also by the best
  incumbent of all of them; BranchAndBoundSearch takes the worker whose
  evaluator it uses, and prunes also by BBOptions::shared_incumbent, a
  value read from outside at each extraction
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
  find_package(PBeFS)) and with the makefiles, which adapt to the
  operating system and to the compiler (makefile-OS), take the places of
  the external libraries from extlib/makefile-default-paths-linux, -macos
  and -win unless the command line or extlib/makefile-paths give them, and
  give a client the objects (makefile, makefile-c) or the library already
  built (makefile-inc, with makefile-lib to build and install it); the
  version from the git tag; the CI builds and tests with CMake and with the
  makefiles, with g++ and clang++, on Linux, macOS and Windows (MSYS2, and
  MSVC with CMake)
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
  StealingBranchAndBound.h, FarmBranchAndBound.h, MPIBranchAndBound.h),
  each with the section of the paper it follows and with why it is correct, and
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
  distances of its abstraction with the blank ignored or minimized, built by
  a breadth-first search by layers in parallel on any runtime (two bits per
  state with the blank minimized), AdditivePDB sums the tables of a
  partition (7-8, 5-5-5 and 6-6-3 of the 15-puzzle, 6-6-6-6 of the
  24-puzzle) and takes the maximum with the reflection about the main
  diagonal, and SlidingTilePDB is SlidingTile with that heuristic; the
  entries in bytes, in 4 bits as the excess over the Manhattan distance,
  modulo 3 in 1.6 bits and in 1 bit (the heuristic of a successor from that
  of its parent), compressed by the least excess of groups of ranks, and by
  rank and select over the nonzero excesses, all of them on the packed
  arrays and on the rank and select of ADS; the tables saved to and loaded
  from files with a versioned header and a checksum, and kept in a cache
  directory; the tests pdb and pdb-search, tools/PDBBench and
  tools/bench/pdb.sh, and the theory in docs/theory/pattern-databases.md

### Changed

- AStarBench is compiled in several translation units: the options and the
  runs in AStarBench.h and AStarBench.cpp, the searches in
  AStarBenchSolve.h, instantiated by AStarBenchPart.cpp for one Domain of
  AStarBenchDomains.h and one of its variants (plain, with slow expansions,
  with slow edges) in each unit that CMake and the makefile compile from it,
  since one unit with every search on every Domain takes more than 64 GB to
  compile, and one with the three variants of a few Domains more than 12 GB;
  the tester likewise, in a file for each family of cases with the parts
  they share in test/test.h, since test.cpp alone took 14 GB (18 GB under
  ThreadSanitizer)
- a Solution carries the lower bound on the optimal cost that its search
  proves (bound, the cost itself for the optimal searches, the parallel
  IDA*, the bidirectional, external and partial-expansion searches, Fringe
  Search and RBFS among them) and the paths of
  an anytime search (improvements); Status::solved means a path optimal or
  within the bound of the weight, and an anytime search that hits the node
  limit keeps its best path
- ConcurrentClosed::update() takes the Weight whose rule decides whether a
  cheaper path reopens a node; DDDAStar throws std::invalid_argument if the
  options give a weight other than 1, as do
  ParallelIDAStar, the bidirectional searches, PEMSearch, IDDAStar, PEA*,
  EPEA*, FringeSearch and RBFS, and IDAStar if they ask for an anytime
  search
- PASE checks the safety of a state only against the states of BE of
  smaller key, as those of OPEN, the others being independent of it by
  Theorem 2 of Phillips, Likhachev and Koenig (2014), and keeps BE as an
  ordered set
- SlidingTile::for_each_feature_change( s , t ) finds the two cells of the
  move as the nonzero fields of s xor t, instead of the blanks of s and of
  t, and Word128 has the exclusive or: 30 instructions fewer per expansion
  of HDA* with Z on Korf #1 with one worker
- PASE with integer costs orders OPEN by f and g packed in one word and by
  the address (packed_fg, as MQAStar), instead of a tuple: the same
  expansions, and 3465 to 2864 instructions per expansion on Korf #1 with
  one worker
- SelfSplitBranchAndBound with one worker gives it the root with neither
  the sampling nor the dive, and evaluates the nodes of the sequential
  search instead of up to 98 more on the knapsacks of
  BranchAndBoundBench; the page on branch-and-bound says how the batch of
  the farm trades the messages for the nodes evaluated
- the parallel searches (SPA*, MQA*, HDA* also over MPI, SafePBNF, PA*SE
  and the farm) skip the move back to the parent in an OperatorDomain, as
  the sequential A* does, the operator kept in the node and carried by the
  messages of HDA* (a member of no size when the Domain has no operators),
  so that on Korf #1 a search no longer looks up 14 million duplicates;
  for_each_child() calls the successors in either kind of Domain, and
  SlowDomain forwards the operators of the Domain it slows down

- the makefiles include the headers of ADS as those of the system
  (-isystem), as the CMake target ADS::ADS does, so that a warning of a
  header of ADS (e.g., -Wpedantic on __int128) is not one of PBeFS

- ADS is a dependency of the library, found by find_package(ADS) or
  downloaded at a fixed commit (PBEFS_ADS_COMMIT, ADSSDR for the
  makefiles): the data structures of the searches, the executors that run
  their workers, the channels and the mailboxes between them with the
  detection of the end, in shared memory and over MPI, the farm with
  feedback, the hash functions and the placement of the workers on the
  CPUs are those of ADS, instead of being written in the library; the farms
  run on every executor, and the tester has the cases farm-thread,
  farm-openmp, bb-farm-thread and bb-farm-openmp
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
  submodule; that commit and the one of ADS follow the version of PBeFS in
  a build directory already configured, unless they were given by hand
- the library is PBeFS (Parallel Best-First Search), instead of
  AStarSearch, since it holds the branch-and-bound beside A*: the CMake
  package, the library, the directory of the headers, the include guards,
  the namespace (PBeFS_di_unipi_it) and the macros (PBEFS_*) are renamed,
  and the headers are in three layers, core/ (the helpers and the
  parameters of the parallel searches),
  astar/ (the graph searches and their domains) and bnb/ (the tree
  searches), each with its header (core.h, astar.h, bnb.h) and PBeFS.h for
  all of them; the parameters of the parallel searches (ParallelParams) are
  written once for A* and branch-and-bound
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

### Fixed

- PASE returns a goal only if, besides being safe, no state of OPEN or BE
  has g + eps h smaller than its g: with several goals, the safety of a goal
  bounds its g by the cost of a cheapest path to that goal, not by the
  optimum, and on 50 instances of the 2D navigation with 5 workers and no
  weight PA*SE returned on 8 of them paths up to 4.1% longer than the
  optimum
- PASE and EPASE expand no state whose g + eps h is not smaller than the
  least g of a goal in OPEN, which the workers expanded while that goal
  waited for the test above, and EPASE evaluates no edge of such a state of
  BE (the option prune, on by default)
- the test of HDA* runs the 15-puzzle only with a core for each worker:
  with more workers than cores, a worker without a core may hold the path
  to the goal while the others expand above the optimum, and one run with
  8 workers on Windows expanded 36980409 nodes on Korf #55, which one
  worker solves in 151996, and went past the 900 seconds of the test

[Unreleased]: https://github.com/dmeoli/PBeFS/commits/main
