# Parallel Best-First Search

[![CI](https://github.com/dmeoli/PBeFS/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/dmeoli/PBeFS/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/dmeoli/PBeFS)](https://github.com/dmeoli/PBeFS/releases)
[![Coverage Status](https://coveralls.io/repos/github/dmeoli/PBeFS/badge.svg?branch=main)](https://coveralls.io/github/dmeoli/PBeFS?branch=main)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://dmeoli.github.io/PBeFS/)
[![Manual](https://img.shields.io/badge/docs-Manual-blue)](https://dmeoli.github.io/PBeFS/manual/)
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
Crainic [[16]](#ref-16). They run on C++ threads, OpenMP, FastFlow [[12]](#ref-12) and MPI,
through the executors and the mailboxes of ADS; hence, one can compare the
same algorithm across runtimes, and the same runtime across algorithms.
The headers are in three layers: `include/PBeFS/core` holds the helpers
and the parameters of the parallel searches, which do not depend on the
kind of search, `include/PBeFS/astar` the graph searches with their
domains, and `include/PBeFS/bnb` the tree searches; `core.h`, `astar.h`,
`bnb.h` and `PBeFS.h` include a layer, or all of them, save the searches
over MPI.

## Manual

The algorithms of the library, with their invariants, their proofs, their
pseudocode and their code, are the subject of a manual whose Markdown sources
are in [`manual/`](manual/), built into a PDF,
[`manual/PBeFS-manual.pdf`](manual/PBeFS-manual.pdf), which is also attached to
each release, and into the web site <https://dmeoli.github.io/PBeFS/manual/>.

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
compile time, no virtual function is paid in the inner loop, which is one
of the design choices of the fast implementation of Burns et al.
[[4]](#ref-4).

Some algorithms need more than a `Domain`, and each extra requirement is a
refinement of the concept. Thus, `FeatureDomain` describes a state as a set of
propositions (as Zobrist hashing requires), `IncrementalFeatureDomain` gives
the features that a move changes, so that the Zobrist hash of a successor is
computed from that of its parent, and `AbstractFeatureDomain` adds a projection
of the features, for Abstract Zobrist hashing [[5]](#ref-5).
`AbstractStateDomain` and `NBlockDomain` give an abstraction of the state space
(for AHDA* and PBNF [[3]](#ref-3)), `PairwiseDomain` a heuristic between any
two states (for PA*SE [[6]](#ref-6)), and `PerfectHashDomain` an injective map
of the states to the integers. `EdgeDomain` gives the actions of a state, which
cost nothing to list, and evaluates the edge of one action at a time, which may
cost much (a collision check, a simulation), as the searches by edges need
[[121]](#ref-121); `HeterogeneousEdgeDomain` tells the cheap edges from the
expensive ones [[122]](#ref-122), and `LazyEdgeDomain` gives the successor of
an action and a lower bound on the cost of its edge before the edge is
evaluated [[124]](#ref-124). As examples, the library provides the sliding-tile
puzzle of any size with the Manhattan distance (`SlidingTile.h`, which
satisfies all of the concepts) and pathfinding with 4 or 8 moves on a grid with
cell costs (`Grid.h`), read from a Moving AI map [[15]](#ref-15) or generated
at random. The same header reads the boards of Burns et al. [[3]](#ref-3),
whose moves to a corner may cut a blocked cell and whose costs are either unit
or life costs (leaving a cell of the row $y$ costs $y$); for the latter the
heuristic is the exact cost of a path on the board without blocked cells (with
the moves to a corner costing as the straight ones when there are 8 moves),
which we derive in closed form in `Grid.h` and which is therefore consistent.
`SequenceAlignment.h` is the optimal alignment of several protein sequences (up
to 10 by default) as a shortest path in the $n$-dimensional lattice of the
positions [[41]](#ref-41), [[42]](#ref-42), with the sum-of-pairs cost $17 - S$
of the PAM250 scores $S$ [[44]](#ref-44) (25 for a residue against a gap, 0 for
two gaps, which amounts to maximizing the score with a linear gap penalty of
16.5), or with any linear gap penalty $P$ (the costs $2(17 - S)$ and $17 + 2P$,
whose minimum is the maximum of the score with the penalty $P$, e.g., the
penalty 8 of [[5]](#ref-5) and [[43]](#ref-43)), and the heuristic of the
optimal pairwise alignments of the suffixes, computed once by dynamic
programming; its features are the positions in each sequence, and its
abstraction, for both Abstract Zobrist hashing and AHDA*, cuts the lattice into
boxes of $b^n$ points, which keeps the moves inside a box on one worker. The
lattice is a levelled graph (a move raises $\sum_i p_i$), and
`HyperplaneDistribution` of `Distribution.h`, for any `LevelDomain`, is the
hyperplane work distribution of Kobayashi et al. [[43]](#ref-43), which gives
the planes $\sum_i p_i = c$, or slabs of thickness $d$ of them, to the workers
in turn, modulo their number, so that the successors of the states of a worker
go to about $n/d$ workers and fewer states are reopened; `thickness()` is their
rule for $d$. Random instances as those of [[5]](#ref-5) (4 to 9 sequences of
about $5000/n$ residues) come from `SequenceAlignment::random_sequences()`.
Also, `ExplicitGraph.h` holds a weighted graph given by its arcs, which serves
for the route-finding problem on the map of Romania of Russell and Norvig
[[13]](#ref-13), and `RobotNav2D.h` is the navigation of a square robot by long
straight motions, each checked for collisions cell by cell, on Moving AI maps
scaled up, the domain on which the searches by edges were compared
[[122]](#ref-122), [[123]](#ref-123); its maps and pairs of points are in
`data/robot-nav-2d`, under the license of the code of those experiments. The
bidirectional searches need a `BidirectionalDomain`, which adds the goal state,
`for_each_predecessor( s , h_s , f )` along the edges reversed and a backward
heuristic, estimating the cost from the initial state, and they take the least
cost of an edge from a `MinCostDomain`; the sliding-tile puzzles, the grids and
the explicit graphs satisfy both concepts (and `PredecessorDomain`, the
predecessors without the backward heuristic, which LPA* needs), the
predecessors of a tile board being its successors, those of a grid cell its
free neighbours at the cost of leaving them, and those of a graph node its arcs
reversed, kept apart.
Also, `DynamicGraph.h` holds a directed graph whose arc costs, blocked nodes
and goals change, with the gridworlds of the experiments of Koenig, Likhachev
and Furcy [[54]](#ref-54) and their changes, and an `EventDomain` (as `Grid`
and `DynamicGraph`) tells each change of the cost of an edge with the cost
before and after it, coalesced per edge, the way of Ramalingam and Reps
[[144]](#ref-144) of telling insertions and deletions as changes from and to
$\infty$.

## The algorithms

We start from the sequential search (`astar( d )` in `AStar.h`), which is A*
with the goal test at the expansion and a node reopened when a cheaper path
to it is found. Ties on f are broken in favour of the larger g, which on
domains with many states in the last f layer usually reaches the goal much
earlier; optionally, the search keeps the cheapest goal generated so far
and stops as soon as no entry of OPEN has a smaller f, which with an
admissible heuristic is still optimal and, when a goal is generated from a
node of f below the optimum (as with the blind heuristic), expands no node
of the last f layer at all (`SearchOptions::goal_at_generation`, off by
default). Next to it, `idastar( d )` (`IDAStar.h`) is the iterative-deepening
A* of Korf [[7]](#ref-7), a depth-first search bounded by f and repeated with
a growing bound, which keeps only the current path and hence solves the
problems whose CLOSED would not fit in memory, such as the 24-puzzle, at the
price of expanding again the states reached by more than one path. Between
the two, `ddd_astar( d )` (`DDDAStar.h`) keeps OPEN and CLOSED in files on
the disk, in buckets of equal g and h split by the hash of the state, which
it expands by increasing f and then g, the duplicates of each file being
removed in memory against the closed files of the same h, as in the
External A* of Edelkamp, Jabbar and Schrödl [[39]](#ref-39) with the
hash-based delayed duplicate detection of Korf [[40]](#ref-40); its memory
is that of the largest file of a bucket rather than that of all the states,
at the price of the time of the disk, and it needs integer costs, a
consistent heuristic and states that are trivially copyable. Two baselines
complete the sequential side: `breadth_first( d )` (`BreadthFirst.h`), which
ignores costs and heuristic and returns a path of fewest edges, optimal only
when the edges cost the same, and `cost_bounded_dfs( d , bound )`, the single
iteration of IDA* with a given bound, whose path, if any, has no node of $f$
above the bound but need not be optimal. The batch weighted A* of Agostinelli
et al. [[126]](#ref-126) (`BatchAStar.h`) takes the $K$ best nodes of OPEN at
each cycle, orders the nodes by $g + w h$ and computes the heuristic of the new
successors of a cycle in one call, `heuristic_batch()`, when the domain
provides it (the concept `BatchHeuristicDomain`), as a learned heuristic on an
accelerator wants; a goal does not stop the search before every node of
smaller $f$ has been expanded, so that with $w = 1$ the path is optimal for
every $K$, and with $K = 1$ the search is A* with the heuristic $w h$. The
parallel searches are the following, each one with the reference for its
correctness and its typical behavior.

- `SPAStar` (`SPAStar.h`) is the Simple Parallel A* of Algorithm 2 of [[1]](#ref-1),
  in which the workers share one OPEN and one CLOSED list, protected by
  locks. It is the baseline of the survey, and it is not expected to scale
  beyond a few threads, since each expansion goes through the shared lists.
- `MultiAStar` (`MultiAStar.h`) is a portfolio: every worker runs A* alone, each with the ties broken in an order of its own,
  and the first one to end gives the path; it needs $p$ times the memory of
  A*, and it pays only when the expansions of the last $f$ layer depend much
  on the ties.
- `MQAStar` (`MQAStar.h`) is SPA* with a relaxed OPEN list, the MultiQueue
  of [[37]](#ref-37) from ADS: $2p$ queues under a lock each, a node pushed
  into one drawn at random and popped from the better of two, so that the
  workers rarely wait for one another, at the price of expanding nodes a
  little out of the order of f; the search ends when a counter of the nodes
  in OPEN or being expanded falls to 0.
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
  assignment. The channels are a template parameter as well, and with
  `LockedChannelsOf` (inboxes under a lock, which the receiver waits for) and
  the option `sync_send` (each node pushed alone when it is generated) the
  search is PRA*, the Parallel Retracting A* of Evett et al. [[127]](#ref-127)
  without its retraction, and APRA* with the abstract distribution, i.e.,
  the synchronous forms that Burns et al. [[3]](#ref-3) compare with the
  asynchronous one; with the option `outsource` a worker whose next node is
  at least that much better than the next node of another worker lends it
  to that worker, which expands it at once while the node stays in the
  CLOSED list of its owner.
- The distributions of `AutoDistribution.h` are built from the structure of a
  `FeatureDomain`, as Jinnai and Fukunaga do for planning [[5]](#ref-5),
  [[141]](#ref-141), [[138]](#ref-138): a `FeatureModel` expands a sample of the
  instance by A*, up to a given number of expansions or to the goal (so that it
  holds the states of f < f* that every distribution has to expand), takes the
  groups of the features from a `FeatureGroupDomain` or infers them from the
  moves, and weighs each edge of the transition graph of a group (its DTG) both
  by the ground operators that realize it and by the moves of the sample that
  cross it. On it, GAZHDA* cuts each DTG in two balanced halves; GRAZHDA* cuts
  it where the sparsity (the product of the sizes of the two sides over the
  weight of the edges cut) is largest, with the weights of [[5]](#ref-5) or, as
  a variant of ours, with those of the sample (the expanded states on each side
  and the moves across), by enumeration up to 20 nodes and by the passes of
  Fiduccia and Mattheyses [[142]](#ref-142) above, so that it finds the two
  halves of the 15-puzzle and the two lines against three of the 24-puzzle
  reported in [[5]](#ref-5); FAZHDA* ignores the groups that the moves change
  most often; DAHDA* keeps the groups chosen by the greedy abstraction of Zhou
  and Hansen [[139]](#ref-139) until they hold a fraction of the features, and
  ignores the others, as the AHDA* of [[3]](#ref-3) does with a bound on the
  abstract states; and OZHDA* chooses the Zobrist table so that the most
  frequent operators leave the owner unchanged. Since the key of Zobrist
  hashing is a linear map over GF(2) of the vector of the features, a move
  keeps its owner when its change vector is in the kernel of that map, and all
  these distributions are choices of a subspace of that kernel; KZHDA* chooses
  it directly, inserting the change vectors of the most frequent operators of
  the sample into a basis by Gaussian elimination, as many as give the largest
  estimated efficiency, and draws the table uniformly among those that vanish
  on it. `auto_distribution()` estimates the communication overhead, the load
  balance and the efficiency of each of them (and of those given by the caller)
  on the sample, and returns the one of largest estimated efficiency, which
  also tunes on the instance the fractions of the last three methods. The
  hand-made projections of the sliding-tile puzzle compared in [[5]](#ref-5)
  and [[141]](#ref-141) are given by `SlidingTile::tile_abstraction()`, and the
  statistics of a parallel search give the communication overhead, the load
  balance and, from the expansions of A*, the search overhead of [[5]](#ref-5),
  which `AStarBench` prints next to their estimates.
- `HDAStarMPI` (`HDAStarMPI.h`) is the same search with one MPI process per
  worker, a nonblocking exchange of the batches and the termination
  detection of Mattern [[11]](#ref-11) (four counters reduced with `MPI_Iallreduce`).
  At the end, the path is rebuilt state by state, each one sent by the
  process that holds it.
- `TDS` (`TDS.h`) is the Transposition-table Driven Scheduling of Romein
  et al. [[131]](#ref-131), the parallel IDA* from which HDA* takes its
  asynchronous scheme: the transposition table is split among the workers
  by the distribution of HDA*, each node goes to the owner of its entry,
  which pushes it on its stack only if its state was not already pushed in
  the iteration with a $g$ not larger, so that no subtree is searched twice
  in an iteration, and each iteration ends when the mailbox detects that
  no worker has work, the next bound being the least $f$ pruned over all
  the workers. The table keeps every state generated, as the CLOSED list of
  A* (there is no replacement of the entries). IDA*, IDPSDD and TDS take
  the bounds of IDA*-CR of Sarkar et al. [[132]](#ref-132) as an option
  (`DeepeningOptions` in `Deepening.h`): the next bound is the largest $f$
  of the first buckets of a histogram of the $f$ pruned by an iteration
  that hold a given multiple of its expansions, or the bound plus twice its
  last step when fewer nodes were pruned, instead of the least $f$ pruned,
  and the iteration that finds the first goal above that least $f$ is
  completed as a depth-first branch-and-bound, so that the path stays
  optimal; without the completion, its cost exceeds the optimum by at most
  the jump of the bound.
- `HDLPAStar` (`HDLPAStar.h`) is the incremental search of `LPAStar` run by
  the workers of HDA*, each of which owns the states that a distribution
  assigns to it and repairs them: an overconsistent expansion sends to the
  owner of each successor the lookahead $g + c$ that it offers, an
  underconsistent one sends an invalidation, which makes the owner scan the
  predecessors again (reading their $g$ in the tables of the other workers)
  only where it names the parent, and the owner of the virtual goal wakes
  the others when its key rises. An episode, between two barriers, runs
  the underconsistent expansions in a first phase, in which a state raised
  is not lowered again (when its new parent was raised too), and then a
  monotone second phase; without that rule two workers may count to
  infinity on a cycle that the change has cut off, while with it every
  episode ends with an optimal path whatever the interleaving (see
  section 6.9 of the [manual](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-9)).
  A first search, or a restart, is the second phase alone, and an episode
  either repairs or restarts, by the number of the changes and a budget past
  which a repair is stopped and restarted (a prediction learned from the
  past repairs may be added).
  A window on the least key of all the workers, as the buckets of
  $\Delta$-stepping [[145]](#ref-145), keeps them near the order of LPA*.
  `HDDStarLite`, in the same header, is D* Lite run in the same way.
- `SafePBNF` (`SafePBNF.h`) is the livelock-free Parallel Best-NBlock-First
  search of Burns et al. [[3]](#ref-3), in which an abstraction partitions the state
  space into nblocks. A worker expands for a while the nodes of an nblock
  whose duplicate detection scope (its abstract successors) is not in use
  by another worker; hence, no lock is taken on the single node. With the
  option `safe` false it is plain PBNF, in which a worker leaves its nblock
  only for a better free one, and which may livelock on an infinite state
  space.
- `PSDD`, `BFPSDD`, `IDPSDD` and `DynamicBoundedPSDD` (`PSDD.h`) acquire
  the nblocks as SafePBNF does, but by layers: PSDD is the parallel
  structured duplicate detection of Zhou and Hansen [[128]](#ref-128), whose
  layers are those of a breadth-first search, with the nodes of $f$ above a
  bound pruned as in their breadth-first heuristic search
  [[129]](#ref-129), and BFPSDD the best-first version of Burns et al.
  [[3]](#ref-3), whose layer holds the nblocks of least $f$. In both, the
  goal is tested at the generation, but the search ends only when no node
  below the incumbent is left, so that the path is optimal among those
  within the bound whatever the costs of the edges, and not only when they
  are all equal. IDPSDD raises the bound of PSDD as IDA* does, and
  DynamicBoundedPSDD takes it from the cost of a weighted A*.
  The nblocks of `Grid` are rectangles (`GridNBlocks`), squares by default,
  so that the coarse abstraction of [[3]](#ref-3), a map cut into a given
  number of columns and rows, is available as well.
- `PASE` (`PASE.h`) is the Parallel A* for Slow Expansions of Phillips et
  al. [[6]](#ref-6), which expands in parallel only the states that the pairwise
  heuristic proves not to be reachable through a cheaper path from a state
  still in OPEN, and with a weight on the heuristic their wPA*SE, whose cost is
  within the inflation of the test of the optimum. It may pay off only when an
  expansion costs much more than the bookkeeping; `SlowDomain.h` adds a given
  amount of artificial work per expansion to study this regime.
- `EPASE` (`EPASE.h`) is the Edge-based Parallel A* for Slow Evaluations of
  Mukherjee, Aine and Likhachev [[121]](#ref-121) (ePA*SE, and w-ePA*SE with a
  weight), which applies the test of PA*SE to the edges: a worker evaluates
  one edge at a time, the edges of a state are evaluated in parallel, and an
  edge whose state stops being the best is not evaluated at all, so that it
  evaluates fewer edges than PA*SE and keeps busy more workers than the
  branching factor; one worker selects and hands out the edges, waking the
  others one at a time, as in the paper. With an option it is GePA*SE
  [[122]](#ref-122), which evaluates the cheap edges of a state together and the
  expensive ones apart, and with another the anytime A-ePA*SE (or A-GePA*SE)
  [[123]](#ref-123), which lowers the weight after each path as ARA*
  [[48]](#ref-48) does, evaluating each edge at most once over all the
  iterations. `SlowEdgeDomain` (in `SlowDomain.h`) makes the edges of any
  domain cheap or expensive to evaluate, and some of them missing. On the 2D
  navigation of [[122]](#ref-122), with as many workers as its code, EPASE
  evaluates as many edges as that code, within the variation between two runs,
  or fewer, while PASE, which applies the successors of an expansion together
  where that code applies each one as soon as its edge is evaluated, evaluates
  up to about a quarter more; that code, however, returned paths longer than
  the optimum on up to 11 of 50 instances with no weight, which the test of the
  end of PASE and EPASE avoids (cf. the page of the parallel searches).
- `PWAStar` (`PWAStar.h`) is the Parallel Weighted A* that is the baseline of
  [[121]](#ref-121) and [[124]](#ref-124): weighted A* whose workers evaluate the
  edges of one expansion at a time, so that its speedup is bounded by the
  branching factor, while its expansions are those of the sequential search.
- `MPLP` (`MPLP.h`) is the Massively Parallelized Lazy Planning of Mukherjee,
  Aine and Likhachev [[124]](#ref-124): a worker repeats a weighted A* that uses
  the cost of the edges evaluated and a lower bound for the others, while the
  other workers evaluate the edges, first those of the paths found; a path
  whose edges are all evaluated is accepted when it costs no more than any path
  found, hence within the weight of the optimum. With one worker it is the lazy
  shortest path of Dellin and Srinivasa [[125]](#ref-125).
- `FarmAStar` (`FarmAStar.h`) is a centralized parallel A* written as a
  farm with feedback, i.e., the master-worker scheme of [[1]](#ref-1): a
  master holds OPEN and CLOSED, and the workers generate the successors of
  the nodes it sends them. With the option `synchronous` it is the
  K-best-first search of Felner, Kraus and Korf [[130]](#ref-130): the master
  waits for all the tasks of a round and inserts their successors in the
  order in which it sent them, so that its expansions are the same at every
  run.
- `ParallelIDAStar` (`ParallelIDAStar.h`) is the parallel IDA* below a
  frontier: A* (or uniform-cost search, or a tree of fixed depth) runs
  until OPEN holds a given number of units, as in the AIDA* of Reinefeld
  and Schnecke [[55]](#ref-55) and the A*+IDA* of Bu and Korf
  [[56]](#ref-56), and the iterations of IDA* then run below the units on
  the workers, each unit with its own threshold, so that a unit whose
  subtree has no node of f between its last bound and the current one sits
  the iteration out. The scheme is synchronous (a barrier per iteration,
  the first goal found being optimal) or asynchronous (the units taken by
  increasing threshold with no barrier, until every threshold reaches the
  incumbent); the units are scheduled by work stealing or from a shared
  cursor, longest first, an idle worker may take half of the successors
  not yet tried at the bottom of the stack of a busy one, as in the
  splitting of Rao and Kumar [[58]](#ref-58), and the first levels below a
  unit may be looked up in the table of the frontier, which prunes many of
  the duplicates that IDA* expands again. With one unit and one worker it
  expands exactly the nodes of `IDAStar`, and with the depth rule and
  every unit at every bound it is the parallel IDA* of the code of Siag et
  al. [[57]](#ref-57). `ParallelIDAStarMPI` (`ParallelIDAStarMPI.h`) runs
  the synchronous scheme with one MPI process per worker, each process
  building the same frontier and the units given to the processes, at each
  iteration, by the longest processing time first of Graham
  [[66]](#ref-66) on their last expansions.

The shared-memory searches are templates also over the executor of
[ADS](https://github.com/dmeoli/ADS) that runs the workers, i.e., its pool
of C++ threads (`ADS::ThreadPool`, the default), `ADS::OpenMPExecutor` or
`ADS::FastFlowExecutor`, and they take their channels, the mailbox that
detects the end of HDA* and of the distributed B&B, and the farm with
feedback from ADS as well (on FastFlow the farm is a farm of FastFlow, and
the channels a mesh of single-producer single-consumer queues); a caller
that runs many searches may give them an executor of its own by
`set_executor()`, which saves the start of the threads at each search.
Their options are in `ParallelOptions`: the number of workers, the size of
a batch, the number of expansions between two flushes of a partially
filled batch and the seed of the distribution.

When an optimal path costs too much to find, a weight $w \geq 1$ in the
options (`SearchOptions::weight`) makes A*, IDA*, SPA*, MQA*, HDA* (on
threads and on MPI), SafePBNF and the farm order OPEN by $f' = g + w h$ in
place of $f$, as in the weighted A* of Pohl [[46]](#ref-46), so that the
path they return costs at most $w C^*$, $C^*$ being the optimal cost, and
`Solution::bound` holds the lower bound on $C^*$ that the search proves
(`Weight.h`); the OPEN lists break the ties of $f'$ in favour of the larger
$g$, and `OpenKind::heap_low_g` (`LowGOpenList`) in favour of the smaller
one, i.e., of the smaller $f$, which on grids usually expands fewer nodes
and finds cheaper paths. In the parallel searches the incumbent $z$ prunes
a node when $w f \geq z$, the rule of wPBNF [[3]](#ref-3), and the search
ends when no
node of $f' < z$ is left, which keeps the bound whatever the order in which
the workers expand. With `drop_duplicates` a closed node reached through a
cheaper path is reopened only when the saving exceeds
$(w - 1)(h(u) - h(v))$, $u$ being the new parent; with a consistent
heuristic this keeps the bound in any order of expansion, while it reopens
fewer nodes than the weak duplicate dropping of [[3]](#ref-3), and with an
inconsistent one the bound may fail (Example 6.5 of the manual). With
`anytime` the parallel searches go on after the first path, the incumbent
pruning by $f$ only, until they prove it optimal, as the anytime weighted
HDA* and PBNF of [[3]](#ref-3) do, and `Solution::improvements` records each
path found with its time; with a schedule of weights `SafePBNF` is the
Anytime Restarting PBNF of [[3]](#ref-3), which lowers the weight at the end
of each iteration and orders OPEN again, as ARA* does; the sequential anytime searches are in
`AnytimeAStar.h`, i.e., Anytime Weighted A* (AWA*) of Hansen and Zhou
[[47]](#ref-47), which reports with each path the least $f$ of OPEN as a
lower bound, and Anytime Repairing A* (ARA*) of Likhachev, Gordon and
Thrun [[48]](#ref-48), which runs weighted A* with a decreasing schedule of
weights, expanding a node at most once per iteration and keeping the
nodes improved after their expansion in the list INCONS for the next one.
The other bounded-suboptimal searches are sequential: the optimistic
search of Thayer and Ruml [[49]](#ref-49) (`OptimisticSearch.h`) searches
with the weight $2 w - 1$ and then expands by $f$ only as much as is needed
to prove the path within $w$; the Improved Optimistic Search of Chen et al.
[[50]](#ref-50) (`ImprovedOptimisticSearch.h`) proves or improves the first
path with an A* restarted from the root that reopens each node at most
once; the focal search $A^*_\epsilon$ of Pearl and Kim [[51]](#ref-51)
(`FocalSearch.h`) expands, among the nodes of $f \leq w f_{\min}$, the
closest to a goal by the distance $d$ of a `DistanceDomain` (or by $h$, or
by $g + (2w - 1) h$), and Explicit Estimation Search of Thayer and Ruml
[[52]](#ref-52), in the same header, orders its focal list by an estimate
of $d$ corrected online and chooses among three lists the node to expand;
and the Dynamic Potential Search of Gilon, Felner and Stern
[[53]](#ref-53) (`PotentialSearch.h`) expands the node of largest
$(w f_{\min} - g) / h$. Last, `LPAStar` (`LPAStar.h`) is the Lifelong
Planning A* of Koenig, Likhachev and Furcy [[54]](#ref-54): after the
caller has changed the costs of some edges in the domain and has told each
change with its old and new cost (`edge_changed()`, or `take_events()` of
an `EventDomain`; `update( v )` when the old cost is not known), `solve()`
repairs the previous search rather than starting anew, a decrease costing
$O(1)$ and an increase a scan of the predecessors only where the edge was
the one of the parent, and the test of the end costing $O(1)$ whatever the
number of goals (`Grid::cell_changed()` patches the moves around a cell
whose cost changed and records their changes); `DStarLite`, in the same
header, is the D* Lite of Koenig and Likhachev [[143]](#ref-143), i.e., LPA*
on the reversed domain from the goal to a robot that moves, whose keys grow
by the estimate of each move.

The front-to-end bidirectional searches run a forward search from the
initial state and a backward one from the goal over one table of the nodes,
each node with its g in both directions, so that a state reached by one
direction with a finite g in the other lowers the incumbent by the same
lookup that detects the duplicates (the immediate solution detection of
[[94]](#ref-94)), and they stop when the incumbent is not larger than a lower
bound on the optimal cost drawn from their open lists; the bounds, with
their proofs, and the guarantees of each search are in the page "Bidirectional
search" of the documentation. `BiBestFirst.h` expands one node at a time,
the one of least priority, and gives the bidirectional Dijkstra's algorithm,
MM of Holte et al. [[61]](#ref-61), [[81]](#ref-81), which never expands a node
beyond half of the optimal cost, MMe [[82]](#ref-82), the fractional MM of
Shaham et al. [[83]](#ref-83), which meets at a given fraction of it, and BAE*
[[88]](#ref-88), [[89]](#ref-89), which adds to f the error of the heuristic
of the other direction and needs consistent heuristics. `NBS.h` is the
Near-Optimal Bidirectional Search of Chen et al. [[85]](#ref-85), which
expands both ends of a pair of least lower bound and hence at most twice the
minimum vertex cover of the must-expand graph of Eckerle et al.
[[84]](#ref-84); `DVCBS.h` is the search of Shperberg et al. [[86]](#ref-86),
which expands one cluster of nodes of equal g of a minimum vertex cover of
the dynamic must-expand graph, computed in time linear in the clusters; and
`GBFHS.h` is the Generalized Breadth-First Heuristic Search of Barley et al.
[[87]](#ref-87), which raises a bound on the cost by one unit at a time and
splits it into two limits on the depths, and which, without heuristics, is
the bidirectional breadth-first search (`bibfs()`). With consistent
heuristics a search does not open a node that the other direction has
closed, as in BS* [[90]](#ref-90), since every path through it is known
already; `BidirectionalOptions` holds this assumption, the least cost of an
edge, the choice of the direction (by priority, in turn, or by the size of
the open lists [[91]](#ref-91)), the fraction of fMM and of the split of
GBFHS, the eager variant of NBS and DVCBS, and the request of the largest of
all the valid bounds. `ReverseDomain` swaps the initial and the goal state of
a `BidirectionalDomain`, so that A* on it is the reverse A*, and
`BlindDomain` drops both heuristics.

`ParallelBidirectional.h` runs the two directions at once on two workers of
an executor, sharing the table and the incumbent: each worker writes only
the fields of its direction, publishes the least keys of its open list, from
which either worker computes the bound of the algorithm, and stops both when
the incumbent reaches it. With f as the priority and the two rejections of
the NBA* of Pijls and Post [[92]](#ref-92) this is PNBA* of Rios and
Chaimowicz [[93]](#ref-93); with the priorities of the bidirectional
Dijkstra's algorithm, of MM and of BAE* it is the parallel version of each.

The external searches go beyond the memory. `PEMSearch.h` is the framework
PEM-BiHS of Siag et al. [[57]](#ref-57), whose open and closed lists are
files of buckets of equal g and heuristics, split by the hash of the state;
the partitions of a bucket are loaded, checked for duplicates against the
closed buckets of the same heuristics and expanded by the workers of an
executor, and the meetings of the two directions are found when a copy is
expanded, the delayed solution detection of Sturtevant and Chen
[[94]](#ref-94). It gives PEM-A*, the reverse PEM-rA* (on a `ReverseDomain`),
PEMM and PEM-BAE*, with a cap on the bytes on the disk. `IDDAStar.h` is A*
with immediate duplicate detection in external memory, A*-IDD of Lin and
Fukunaga [[95]](#ref-95), whose closed list is a file of records read at
random, found through a table in memory of positions packed in as few bits
as they need, and written by blocks of one partition of the states, so that
a probe reads the disk only for the blocks of the partition of the state; it
suits a solid state drive, where DDDAStar suits sequential access.

## Greedy, width-based and other searches

When a path of least cost is not required, the greedy best-first search (GBFS)
of Doran and Michie [[96]](#ref-96), which selects the node of least h, usually
reaches a goal after a small fraction of the expansions of A*. `GreedySearch`
(`GreedySearch.h`) performs it, and the search ordered by $g + w h$ for a
finite weight w (the weight of `SearchOptions`, which `GreedyOptions` makes
infinite by default), whose path costs at most $w C^*$ when h is admissible and
the closed nodes are reopened [[46]](#ref-46), with the OPEN lists of the
satisficing planners as options that compose: the lazy evaluation of Richter
and Helmert [[97]](#ref-97), in which a successor enters OPEN with the
heuristic of its parent and is evaluated only when selected (on a
`DeferredDomain` its heuristic is then never computed before); several queues
alternated in round robin, one per heuristic of a `MultiHeuristicDomain`, a
queue breaking the ties of one heuristic by another, and Pareto queues on the
vectors of the heuristics [[98]](#ref-98); queues of the preferred successors
of a `PreferredDomain`, boosted whenever a heuristic reaches a new least value
[[97]](#ref-97); the type-based exploration of Xie et al. [[99]](#ref-99) and
the epsilon-greedy selection of Valenzano et al. [[100]](#ref-100); and the
breaking of the remaining ties by the order of insertion, by g or at random.
`EnforcedHillClimbing` is the local search of FF [[101]](#ref-101), a
breadth-first search from the current state for a better one, with the
preferred successors alone or first, and a fallback on GBFS; `IteratedSearch`
runs a sequence of such searches, each bounded by the cost of the best path
found before it and all sharing the evaluations (`HeuristicCache`), which with
`lama_phases()` is the schedule of LAMA [[102]](#ref-102) and, in general, the
restarting search of Richter, Thayer and Ruml [[103]](#ref-103). The
width-based searches of Lipovetzky and Geffner (`WidthSearch.h`) need the
features of a state: `iw()` is IW(k), the breadth-first search that prunes the
states that make no tuple of at most k features true for the first time,
`siw()` achieves the goal features of a `GoalFeatureDomain` one at a time with
IW [[104]](#ref-104), and `bfws()` is the best-first width search ordered by
novelty, then by the goal features still to achieve (or by h), then by g, the
novelty being computed for each value of these and of the relevant features of
a relaxed plan of a `RelevantFeatureDomain` [[105]](#ref-105).

The parallel greedy searches (`ParallelGreedy.h`) share one OPEN and one
CLOSED list, and differ in the condition under which a worker may take the
best node of OPEN, which is what decides how far they stray from the
sequential GBFS: KPGBFS takes it always [[106]](#ref-106) and may expand
arbitrarily more states than GBFS [[108]](#ref-108); PUHF [[108]](#ref-108)
and PUHF2, PUHF3 and PUHF4 [[109]](#ref-109) take only the nodes that some
breaking of the ties lets GBFS expand, i.e., those of the bench transition
system of Heusner, Keller and Helmert [[112]](#ref-112), and SPUHF keeps the
idle workers on a speculative search of their own; OBAT [[110]](#ref-110)
explores one bench at a time, by deferring the nodes that may leave the
bench, and expands at most $N + k |p|$ states with k workers, N being the
states that GBFS expands with some breaking of the ties and p its path.
Each of them may also separate the generation of the successors from their
evaluation, which any idle worker does [[111]](#ref-111), and may run in
synchronous rounds in one thread, which makes its counts replicable; the
tests check on random graphs and on the 8-puzzle, whose bench transition
system they compute, that the constrained searches expand only its states
and that OBAT keeps the order and the bound above. `HDGBFS` is GBFS on the
hash distribution of HDA* [[107]](#ref-107).

Last, the library has the searches that the parallel ones build on or are
compared with. Partial Expansion A* [[42]](#ref-42), [[113]](#ref-113) puts in
OPEN only the successors whose f does not exceed the stored value of their
parent, which goes back to OPEN with the least f of the others but those
already stored with a g not larger (which could only be generated again as
duplicates), and Enhanced Partial Expansion A* [[114]](#ref-114) generates only
those, on the operator selection function of a `PartialExpansionDomain`
(`PartialExpansion.h`, satisfied by `SlidingTile` and by `Grid` with 4 and 8
moves and weighted cells), so that on the 15-puzzle it stores and generates
about half of the nodes of A*; `FringeSearch` [[115]](#ref-115) keeps the
frontier of IDA* between its iterations; `RBFS` is the linear-space best-first
search of Korf
[[116]](#ref-116); `BeamSearch` keeps the best B nodes of each layer of a
breadth-first search [[117]](#ref-117); `RealTimeSearch` is LRTA*
[[118]](#ref-118) and, with a lookahead of more than one expansion, LSS-LRTA*
[[119]](#ref-119), which interleave a bounded search with the moves of an agent
and learn the heuristic; and `HPAStar` (`HPAStar.h`) is the hierarchical
pathfinding of Botea, Müller and Schaeffer [[120]](#ref-120) on the grids,
which searches a graph of the entrances between clusters of cells and refines
its path within the clusters.
These searches have been checked against those of HOG2 [[165]](#ref-165) on
counts that do not depend on the machine (the cost of the path, the
expansions, the generated nodes and the trials of LSS-LRTA* up to
convergence), on instances of Korf and on random grids of 4 and 8 moves: A*
and LSS-LRTA* agree with it up to the breaking of the ties, PEA* selects no
more nodes on the grids and as many on the 15-puzzle but on one instance,
where it selects more, and EPEA* runs also on the grids of 8 moves, where the
one of HOG2 stops on an assertion; Fringe Search, RBFS and HPA* could not be
compared, since those of HOG2 do not build in its current tree.

## Multi-objective search

A `MultiObjectiveDomain` (`MultiObjective.h`) is a domain whose edges cost a
vector of $K \geq 2$ components, with a heuristic vector consistent in each
component, and its searches return a `ParetoSolution`, one path for each vector
of the Pareto front of the paths from the start to a goal, in lexicographic
order. `MOAStar` (`MOAStar.h`) is the framework of Ren et al.
[[182]](#ref-182), in which the labels (the paths) are extracted by the
lexicographic order of f and checked, with the first component dropped, against
the frontier of the paths expanded at their state and against the solutions
[[185]](#ref-185): with the minimum of the second components it is BOA* of
Hernández et al. [[181]](#ref-181), whose checks take constant time, with a
balanced tree (the AVL tree of ADS) it is EMOA* [[182]](#ref-182), whose check
for three objectives takes logarithmic time, and with an array it is LTMOA*
[[183]](#ref-183); the option `lazy` drops the checks at the generation, which
gives LazyLTMOA*, and the option `order` permutes the objectives of the
lexicographic order; OPEN is a heap, or `MOBucketQueue`, a bucket for each
value of $f_1$ as in [[184]](#ref-184), for integer costs. `NAMOAStarDR` is NAMOA* [[180]](#ref-180) with the
dimensionality reduction of Pulido, Mandow and Pérez de la Cruz
[[185]](#ref-185), which also checks every new label against the labels of its
state in OPEN and removes those it dominates. The parallel searches are in
`ParallelMOAStar.h`: `MOAStarPar` is the parallel MOA* of Ahmadi et al.
[[184]](#ref-184), whose workers search on different orders of the objectives
and share upper bounds that, as the solutions of the other workers advance,
compare fewer and fewer components, down to one; it takes any of the three
frontiers (LTMOA*-par with the array, BOA* on the two orders for two
objectives) and up to $K!$ workers, the cyclic orders first. `HDMOAStar` is our
hash-distributed search, in which each worker owns the states of a hash, as in
HDA*, keeps the frontiers of its states with the full vectors (in a Fenwick tree
of truncated frontiers over layers of $f_1$ [[189]](#ref-189), the indexed
pruning of [[187]](#ref-187)) and sends each solution to the others; since the workers do not extract the labels in one
global order, it is a label-correcting search, which ends with the Pareto front
when no worker has work and no message travels. The proofs of all of them are
on the page "Multi-objective search" of the documentation.

The domains are `MultiObjectiveGraph` (`MultiObjectiveGraph.h`), a graph given
by its arcs, whose heuristic is the ideal point, i.e., the distance to the
target by each objective alone, and which `read_dimacs_objectives` builds from
the road graphs of DIMACS [[64]](#ref-64), the distances and the times
(`tools/bench/fetch_dimacs.sh <dir> NY --times` downloads both), with the
further objectives of the papers computed from the graph (the number of arcs
[[188]](#ref-188), the degree of the arcs [[182]](#ref-182) [[184]](#ref-184)
and a random integer in [1, 100] [[183]](#ref-183)); and `MultiObjectiveGrid`
(`MultiObjectiveGrid.h`), the grids of [[182]](#ref-182), whose moves cost
random integers. `tools/MOABench` runs them, e.g., `MOABench --graph
USA-road-d.NY.gr,USA-road-t.NY.gr --extra degree --queries pairs.txt --algo
par --threads 3 --header`, and prints the size of the front and the counters
of each run.

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
(`DistributedBranchAndBound`), master-worker on a farm with feedback
(`FarmBranchAndBound`), all of them on the three runtimes, with one MPI process
per worker (`MPIBranchAndBoundSearch`, whose nodes travel as bytes), or by
synchronous rounds (`DeterministicBranchAndBound`), whose outcome depends
neither on the timing of the workers nor on their number, as the deterministic
mode of FiberSCIP [[18]](#ref-18), so that the experiments of a paper can be
repeated exactly, or by the SelfSplit of Fischetti, Monaci and Salvagnin
[[19]](#ref-19) (`SelfSplitBranchAndBound`), whose workers search apart, and
without exchanging anything, the subtrees of a first sampling of the tree,
which ends with a dive to a leaf so that they start from an incumbent, the
outcome depending then only on the number of the workers (or, sharing their
incumbent, pruning more at the price of that determinism), or by work stealing
[[20]](#ref-20) (`StealingBranchAndBound`), whose workers keep the children
they generate and take about half of the queue of another worker
[[36]](#ref-36) only when theirs is empty, so that each worker stays in its own
part of the tree, which pays when moving the evaluator of a worker to a distant
node is costly, or on the MultiQueue of ADS (`MultiQueueBranchAndBound`), the
queue of MQA*, whose workers share the open nodes with no lock common to all
of them, at the price of a relaxed order of extraction. The searches that need no external library can also be chosen
by a value, or by its name read from a file of parameters
(`parallel_branch_and_bound( p , bb_scheme( "selfsplit" ) , options )`), so
that a program built on the library gets every new scheme with no change of its
own. In the taxonomy of Gendron and Crainic [[16]](#ref-16), all of them build
one tree in parallel (their type 2), asynchronously but for the rounds, with a
single pool of open nodes (the shared queue, the MultiQueue, the farm and the
rounds) or with a
pool per worker (the distributed queues, the work stealing, the SelfSplit and
MPI); the parallelism inside a node (their type 1) is left to the evaluation of
the problem, and the search of several trees at once (their type 3) is
provided only at the start, by the racing ramp-up below. SafePBNF and PA*SE
have no counterpart, since they exist to handle duplicates and re-expansions,
which a tree does not have; for the same reason, the owner of a node only
balances the load.

A problem may also give the children of a node one at a time, by a cursor,
so that a depth-first search builds only those it visits (with
`prune_level`, the children coming by their estimates, a child pruned ends
its siblings), and the child k of a node without its evaluation, so that a
node can be rebuilt from the path of its choices. On these rest the
depth-first searches of the skeletons of Archibald et al. [[67]](#ref-67)
[[68]](#ref-68) (`SkeletonBranchAndBound`), whose workers share the work by
spawning as tasks the children of the nodes above a depth, by spawning the
children left on their stack after a budget of backtracks, by letting an idle
worker steal from the stack of a busy one, or by taking the nodes of a depth
in a fixed order of discrepancy, which makes the runtimes repeatable; a task
is a node or, with `indexed`, the path that leads to it, a few bytes whatever
the node. The search may seek an optimal value, decide whether a value within
a target exists, stopping at the first one, or enumerate the tree, counting
its feasible nodes, and every scheme runs in the three modes. Before any of
the schemes, several settings may race on the root, sharing their
incumbents, the search going on from the tree of the one with the best
bound, as in the racing ramp-up of the Ubiquity Generator [[69]](#ref-69);
any of them may run in segments that write a checkpoint of their open nodes,
from which a run stopped at any time is resumed with any scheme and any
number of workers (`checkpointed_branch_and_bound()`). Over MPI, besides
the distributed search, the processes may be organized as a master, hubs and
workers, each searching its own subtrees and receiving work only when its
hub pairs it, by the quantity and by the quality (the mean bound of its best
nodes) of the work of the cluster, with a donor, as in ALPS [[70]](#ref-70)
(`MPIHubBranchAndBoundSearch`). For the trees whose nodes of a depth have
the same number of children, as those of the permutation problems, the work
may be an interval of leaves in a mixed radix, an idle worker taking the
second half of the interval of another one, as in the Integer-Vector-Matrix
of Mezmaz et al. and Gmys et al. [[71]](#ref-71) [[72]](#ref-72)
(`IntervalBranchAndBound`), each node being evaluated by one worker only, so
that with the optimum given the nodes evaluated do not depend on the steals.
The examples of the layer, besides the knapsack, are the trees of the
Unbalanced Tree Search of Olivier et al. [[73]](#ref-73) (`UTS.h`, with the
trees T1 to T3 of the benchmark), the maximum clique by the bound of a greedy
colouring on bitsets, as in the algorithms of Tomita and Kameda and of San
Segundo et al. [[74]](#ref-74) [[75]](#ref-75), with a reader of the graphs
of the DIMACS challenges (`MaxClique.h`), the maximum common induced
subgraph by the label classes of McSplit [[76]](#ref-76)
(`MaxCommonSubgraph.h`), the subgraph isomorphism by forward checking
(`SubgraphIsomorphism.h`), the quadratic assignment by the bound of Gilmore
and Lawler [[77]](#ref-77) with a reader of QAPLIB [[78]](#ref-78)
(`QAP.h`), and the permutation flowshop by the one-machine and two-machine
bounds of Lageweg, Lenstra and Rinnooy Kan [[79]](#ref-79), with the
branching of the head, of the tail, alternated or chosen at each node, and
the instances of Taillard [[80]](#ref-80) generated from their seeds
(`Flowshop.h`). The design is discussed in
[docs/branch-and-bound.md](docs/branch-and-bound.md), and the tests check the
variants on the 0-1 knapsack against the dynamic program, and the domains
against references outside the library (e.g., the size of the trees of UTS
by a plain recursion, the cliques by the algorithm of Bron and Kerbosch and
the common subgraphs by the cliques of the modular product of the graphs).

## The data structures

The data structures of the searches are those of ADS. OPEN is either a
4-ary heap (`ADS::DaryHeap`), for any cost type, or a two-level bucket
queue (`ADS::BucketQueue`, on f first, then on g) for integer costs, which
[[4]](#ref-4) reports to be somewhat faster than a heap on the sliding-tile
puzzles; neither one supports decrease-key, and an entry made stale by a
cheaper path is skipped when it surfaces. CLOSED (`NodeTable.h`) is an
open-addressing hash table (`ADS::PointerTable`) that stores the hash of a
state next to its node, so that most probes do not compare the states, and
spreads the slots by Fibonacci hashing; the nodes are allocated from an
arena (`ADS::Arena`), and both the table and the arena are requested from
the operating system as transparent huge pages, where available. In the
sequential search the slots of the successors are also prefetched before
being probed, since on large instances most probes are cache misses (the
parallel searches do the same with the batches they receive). The queues
inside the other searches come from ADS as well: the free list of SafePBNF
and the OPEN of PA*SE are B+-trees (`ADS::BTreeMap`), which give both the
least key and the scan in order that their tests of interference need, and
the queue of the B&B is a pair of indexed heaps on the same nodes
(`ADS::IndexedHeapPair`), one by the key of the selection and one by the
bound, which the test of termination needs. Each search returns a
`Solution` holding the status (solved, unsolvable or stopped at the node
limit), the cost, the path and the counters of the search; for the
parallel ones, these include the load balance and the communication
overhead of [[2]](#ref-2).

The structures of A* can be replaced, at compile time, by those of the
standard library (`std::deque`, `std::unordered_map` and
`std::priority_queue`) or of Boost (`boost::object_pool`,
`boost::unordered_flat_map` and `boost::heap::d_ary_heap`), with the same
interface (`Structures.h`); the searches do not change, and neither does
the order of their expansions, but for the ties of f and g in OPEN, which
the sets break differently. Neither the standard library nor Boost has a
bucket queue, and hence with them both OPEN lists are a heap. With the
structures of ADS the domains on a graph take their adjacency from ADS too:
the grids from `ADS::GridGraph`, which computes the moves from a cell and
their costs from a vector of the costs of the cells (a move paying the cell
it leaves, as in the domain), and the explicit graphs from
`ADS::ArrayGraph`, built from the arcs before the first expansion. The
other two sets compute the moves of a grid from its map, as neither the
standard library nor Boost has a graph of a grid with blocked cells and
diagonal moves (and a graph in compressed rows of a map of 20000 by 20000
cells would take tens of GB), and keep the arcs of an explicit graph in a
list per node or in the `compressed_sparse_row_graph` of Boost. The
standard library and Boost are there to measure the searches against them
(`tools/bench/structures.sh`), and the B&B does not depend on the set.

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
groups of consecutive ranks (the max-compression of [[24]](#ref-24)), or
by rank and select over the nonzero excesses, all of them in the packed
arrays and in the rank and select of ADS;
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
modulo 3 need the blank ignored and no reflection. The default of
`PDBOptions` is nevertheless the blank ignored, which every encoding
accepts and which keeps the heuristic consistent, so that the sequential
A* never expands a node twice; the blank minimized is chosen as in the
example above. The tables of the
24-puzzle (6-6-6-6, 255 MB in 4 bits) are built and looked up as well, and
`SlidingTile` packs the boards of more than 16 cells, 5 bits per cell, into
an integer of 128 bits (`Word128`), so that `SlidingTilePDB< 5 , 5 >` with
IDA* is the search of Korf and Felner [[22]](#ref-22) on that puzzle.
The 50 instances distributed with the code of HDA* of Jinnai and Fukunaga
(`data/24-puzzle/instances50.txt`) are beyond the reach of a benchmark, as
IDA* solves only one of them within $10^9$ expansions (cost 82, with
$6.5 \cdot 10^8$ expansions); `tools/bench/pdb.sh` therefore solves the 50
instances of `PDBBench --walk 150`, i.e., the ends of random walks of 150
moves of the blank from the goal (the $i$-th seeded by $i$), with at most
$10^{10}$ expansions each. 49 of them are solved, with optimal costs from
50 to 84 and $4.7 \cdot 10^{10}$ expansions in all (the median being
$3.5 \cdot 10^7$), while instance 36 exceeds the limit; grouped by the
optimal cost, the expansions are the following:

| optimal cost | instances | least expanded | mean expanded | most expanded |
|---|---|---|---|---|
| 50-60 | 3 | $3.9 \cdot 10^5$ | $7.5 \cdot 10^5$ | $1.3 \cdot 10^6$ |
| 62-66 | 12 | $2.6 \cdot 10^5$ | $1.3 \cdot 10^7$ | $7.2 \cdot 10^7$ |
| 68-70 | 9 | $3.3 \cdot 10^6$ | $4.0 \cdot 10^7$ | $2.2 \cdot 10^8$ |
| 72-74 | 11 | $2.8 \cdot 10^6$ | $4.9 \cdot 10^8$ | $2.2 \cdot 10^9$ |
| 76-78 | 9 | $2.1 \cdot 10^8$ | $2.9 \cdot 10^9$ | $9.5 \cdot 10^9$ |
| 80-84 | 5 | $2.3 \cdot 10^8$ | $3.1 \cdot 10^9$ | $5.5 \cdot 10^9$ |

The mean grows by about one order of magnitude every 6 moves of the
optimal cost up to 78, and within a group of costs the expansions of two
instances may differ by up to 3 orders of magnitude.

### The permutation puzzles and their pattern databases

`PermutationPDB.h` builds the pattern databases of any permutation puzzle,
i.e., a puzzle whose moves permute groups of objects among their positions
and possibly turn them (the concept `PermutationDomain`): the abstract
state of a pattern is the positions and the orientations of its objects,
ranked by one of the rankings of `Permutation.h` (lexicographic by
population counts, lexicographic by the binary tree of Bonet
[[62]](#ref-62), or the linear-time ranking of Myrvold and Ruskey
[[34]](#ref-34)), and its table, of 4 or 8 bits per entry, is built by
`abstract_distances()` (`AbstractionPDB.h`), the algorithm of Dial
[[65]](#ref-65) by layers on the table itself, in parallel on any runtime,
for any abstraction whose states are ranked and whose edges have small
integer costs. `PermutationHeuristic` takes the maximum (or the sum) of
several tables, and `PDBDomain` gives any Domain that heuristic besides its
own. The puzzles are the pancakes of Dweighter [[63]](#ref-63) and their
burned variant (`Pancake.h`), with the gap heuristic of Helmert
[[60]](#ref-60) and its weakened GAP-x of Holte et al. [[61]](#ref-61),
TopSpin (`TopSpin.h`), with the breakpoints, as in [[30]](#ref-30), and
Rubik's cube (`RubiksCube.h`), whose moves are those of the cubie model of
Kociemba, with the pruning of the redundant pairs of moves of Korf
[[59]](#ref-59) (`TreeOperatorDomain`, which `IDAStar` and
`ParallelIDAStar` use) and his three tables, the corners and two sets of
six edges; the tests check the size of his search tree at each depth and the
size, largest value and mean of each table against his paper, and
`data/rubik` holds his ten instances. The Towers of Hanoi with four pegs
(`TowersOfHanoi.h`) take the additive tables of ranges of disks of Felner,
Korf and Hanan [[23]](#ref-23) on the same construction, and `data/hanoi`
holds the 20 instances of 20 disks of Siag et al. [[57]](#ref-57). The
sliding-tile puzzles keep the tables of `PatternDatabase.h`, whose
encodings relative to the Manhattan distance and whose blank minimized need
the structure of the puzzle. Finally, `DimacsGraph.h` reads the road graphs
of the 9th DIMACS Implementation Challenge [[64]](#ref-64) into an
`ExplicitGraph`, with the consistent heuristic of the great-circle distance
scaled by the least ratio of the length of an arc to that of its ends;
`data/dimacs` holds a small synthetic graph in their format, and
`tools/bench/fetch_dimacs.sh` downloads the real ones.

## Classical planning

`Planning.h` makes a domain of a classical planning task in the SAS+
formalism of Bäckström and Nebel [[133]](#ref-133), read by `SASTask.h` from
the file `output.sas` that the translator of Fast Downward [[134]](#ref-134)
writes (the version 3 of the format, without axioms and conditional
effects, which are refused with an exception), so that every search and
every work distribution of the library runs on the planning tasks as on
the other problems. A state packs the values of the variables into an
array of W words of 64 bits (`Planning< W >`, each variable in a field of
the bits its range needs and none across two words), the successors are
enumerated by a decision tree over the variables, as the successor
generator of Fast Downward does, and the costs are those of the operators,
or 1 for a task without metric. The heuristics, in `SASHeuristics.h`, are
the admissible $h^{\max}$ of Bonet and Geffner [[135]](#ref-135), the
landmark-cut heuristic of Helmert and Domshlak [[136]](#ref-136), a pattern
database of a projection of the task onto a set of variables (by default
the goal variables and then their parents in the causal graph, as long as
the abstract space has at most $2^{16}$ states), computed once by a
Dijkstra search backward from the abstract goal [[137]](#ref-137), the
merge-and-shrink abstraction of Helmert, Haslum, Hoffmann and Nissim
[[140]](#ref-140), merged linearly in the order of the causal graph and of
the goal and shrunk by greedy bisimulation to at most 1000 states (the
setting of the experiments of [[5]](#ref-5)), and the blind heuristic, together with $h^{\rm add}$, $h^{\rm FF}$ [[101]](#ref-101)
and the number of unsatisfied goals, which are not admissible. The task
is also an `IncrementalFeatureDomain`, the features being the facts
$v = d$ of its STRIPS view, as in the HDA* for planning of Kishimoto,
Fukunaga and Botea [[2]](#ref-2), and a move changing only the facts of
the variables that its operator sets. The abstraction of the features is
that of FAZHDA* [[138]](#ref-138), i.e., the 30% of the variables changed by
most operators are left out, and the values of each other variable are
split in two by the greedy bisection of its domain transition graph; the
nblocks of SafePBNF, which are also the abstract states of AHDA*, are the
states of the projection chosen by the greedy abstraction of structured
duplicate detection of Zhou and Hansen [[139]](#ref-139), with at most 10000
nodes, and its abstract graph is taken undirected, since SafePBNF counts
the successors of an nblock among its predecessors. Finally, an operator
$o'$ is an inverse of $o$ if it sets back each variable that $o$ changes
to the value that $o$ asks before and its preconditions hold after $o$, so
that $o'(o(s)) = s$; the first 255 operators with an inverse are given an
operator code, and a state generated by one of them does not generate its
parent (`OperatorDomain`).

`data/planning` holds 28 small tasks of the International Planning
Competitions (gripper, blocks, logistics00 and driverlog), whose list,
`optimal.txt`, gives the cost of an optimal plan of each. They were
written by `tools/pddl2sas.py` from the PDDL files of the collection of
Fast Downward: the actions are grounded by relaxed reachability, and the
atoms are grouped into variables by the invariants "at most one atom of
the group holds", which are proved on the ground actions. With
`--domain planning --file <output.sas> --heur <heuristic>`, `AStarBench`
solves a task with any search of the library, the state in 2 or 16
words.

## Building and using the library

The library needs a C++20 compiler, CMake 3.21 or GNU make, and ADS, whose
headers CMake finds by `find_package(ADS)` (with `ADS_ROOT` as a hint) or else
downloads at the commit `PBEFS_ADS_COMMIT`, and the makefiles take from
`ADSSDR`. OpenMP, MPI and FastFlow are optional, and each of them is used if it
is found. FastFlow is header-only: CMake looks for it in `FASTFLOW_ROOT`
(`/opt/fastflow` by default) and on the include path, and if it is not there
downloads the commit `PBEFS_FASTFLOW_COMMIT` (unless `PBEFS_FETCH_FASTFLOW` is
`OFF`), while the makefiles take its headers from `FF_ROOT` (`/opt/fastflow` by
default on Linux and macOS, none on Windows). With CMake one has

```sh
git clone https://github.com/dmeoli/PBeFS.git
cmake -S PBeFS -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
cmake --install build --prefix <prefix>
```

where `PBEFS_WITH_OPENMP`, `PBEFS_WITH_MPI` and `PBEFS_WITH_FASTFLOW` (all `ON`
by default) can exclude a runtime, `PBEFS_STRUCTURES` chooses the data
structures of A* (`ads`, the default, `std`, or `boost`, which needs Boost 1.81
or later; the makefiles take it from `PBeFSSTRUCT`), and `PBEFS_BUILD_TESTS`,
`PBEFS_BUILD_EXAMPLES` and `PBEFS_BUILD_TOOLS` (`ON` only when PBeFS is the
top-level project) the tests, the examples and the benchmark. The presets of
`CMakePresets.json` give the options of the compiler for the tests and the
tools, `--preset gnu` the warnings of g++ and clang++ and `--preset msvc` those
of MSVC with `/bigobj`, without which the tester does not compile with MSVC
(e.g., `cmake --preset msvc -S PBeFS -B build`); they set `PBEFS_CXX_OPTIONS`,
which may also be given directly. A CMake project then uses the library with

```cmake
find_package(PBeFS REQUIRED)
target_link_libraries(myTarget PRIVATE PBeFS::PBeFS)
```

adding `PBeFS::OpenMP`, `PBeFS::MPI` or `PBeFS::FastFlow` for the corresponding
runtime; these are the components of the package, each present if the library
was built with its runtime and the runtime is found again by the project that
uses it; hence, `find_package(PBeFS REQUIRED COMPONENTS OpenMP)` fails without
OpenMP, while `find_package(PBeFS)` finds the package with no `PBeFS::OpenMP`.
Without CMake, `pkg-config --cflags --libs pbefs` gives the flags of the
library alone. The packages of vcpkg, conda-forge, Homebrew, Debian and Ubuntu
(`libpbefs-dev`), Conan, Spack, Nix, the AUR, MSYS2, xmake and Meson WrapDB
install the same files from the recipes of [packaging](packaging/README.md), on
the package of ADS of the same channel, with OpenMP where the compiler has it,
MPI only where the channel lets one choose it, and without FastFlow, which no
channel packages. Without installing anything, a CMake project can also take
PBeFS at configure time, with FetchContent (CMake 3.24 or later),

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
tools. With Meson, `meson wrap install pbefs` and `dependency('pbefs')`, from
WrapDB, in C++20. Since the searches are templates over the `Domain`, they are
in the headers, and the library (`src`) holds only what does not depend on it,
i.e., the grid maps and the measures of a search. Without CMake, the makefiles
in the root give the library in three ways to the makefile of a program that
uses it, which needs only to include one of them: `makefile` builds the objects
and exports them together with the include directive and the list of the
headers, `makefile-c` adds the external libraries of the runtimes (the threads,
OpenMP and FastFlow), and `makefile-inc` gives instead the library already
built, `lib/libPBeFS.a`, which `make -f makefile-lib` builds and `make -f
makefile-lib install PREFIX=<prefix>` installs. `CC` chooses the compiler
(`CC=g++` or `CC=clang++`), and the places of the external libraries and the
switches of OpenMP are those of `extlib/makefile-default-paths-linux`, or of
`-macos` and `-win` with `PBeFSPLATFORM=macos` or `win` (FastFlow in
`FF_ROOT`, ADS in `ADSSDR`, Boost in `BOOST_ROOT`, the OpenMP runtime of LLVM
in `OMP_ROOT` and the switches of OpenMP in `OMP_SW` and `OMP_LIB`, those of
Apple clang on macOS, where the g++ of Homebrew takes `OMP_SW=-fopenmp
OMP_LIB=-fopenmp`), unless they are given on the command line or in
`extlib/makefile-paths`, which git ignores. In `examples`, `test`
and `tools` a makefile builds the corresponding programs, in `test` the targets
`check` and `check_mpi` also run the tests, and in the root `make`, `make
test`, `make test_mpi` and `make install` do the same for the whole repository.
Both CMake and the makefiles are tested on Linux, macOS and Windows (with
MSYS2, and CMake also with MSVC), with g++ and with clang++; on Windows
FastFlow, which relies on POSIX, and MPI are left out.

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
chapter 3 of Russell and Norvig [[13]](#ref-13), each one with its `Domain` in
a header of `include/PBeFS/astar` (e.g., `NQueens.h` or `Romania.h`), which the
benchmark shares. `8Puzzle` solves the instance of figure 3.28 of the book, or
an instance at a given distance from the goal, with the number of misplaced
tiles (h1), the Manhattan distance (h2) or the sequence score of Nilsson
[[8]](#ref-8), which is not admissible; with `-t` it repeats the comparison of
h1 and h2 of figure 3.29 (nodes generated and effective branching factor).
`15Puzzle` solves an instance of Korf, sequentially or with HDA* (`-t
<threads>`); `RouteFinding` finds a road between two cities of the map of
Romania of the book; `FindPath` finds a path on a weighted grid, random or read
from a Moving AI map, and draws it. `VacuumWorld` is the vacuum world of
section 3.2.1 on a row of n cells, and `NQueens` places n queens with the
incremental formulation of the same section, checking the placement apart from
the search. The last three come from the exercises of the book:
`MissionariesCannibals` carries m missionaries and c cannibals over a river
with a boat of b seats (exercise 3.9), `WaterJugs` measures a quantity of water
with jugs (exercise 3.6(d)), and `TravelingSalesman` solves the traveling
salesperson problem with the minimum spanning tree heuristic (exercise 3.30),
checking the tour against the dynamic program of Held and Karp [[14]](#ref-14).
Each run with a known answer (e.g., the 26 moves of figure 3.28, the 418 km
from Arad to Bucharest and the 11 crossings of the missionaries) is also a test
of `ctest`, labelled `example`. In `tools`, `AStarBench` runs any search,
runtime and distribution of the library on the 100 instances of Korf
(`data/15-puzzle/korf100.txt`), on the 8-puzzle, on the 24-puzzle, on the grids
(with unit or weighted cells), on the boards of Burns et al. (`--domain board`,
with `--origin` telling whether the row 0 of the start and of the goal is the
last line of the board, as in the files of the authors, or the first one, as in
those of later generators), on multiple sequence alignment (`--domain msa`,
from a file or random, with another substitution matrix by `--matrix`, a gap
penalty by `--gap`, and the hyperplanes by `--dist HWD`) and on the other
problems of the examples, sized by `--size` (e.g., the number of queens, or of
the cities of a random instance or of a file, with the Euclidean or the
Manhattan distance of `--metric`), on the 2D navigation (`--domain nav2d`),
on the gridworlds of `DynamicGraph.h`, or on a Moving AI map whose cells
`--scale` multiplies, changed by `--change`, `--flips` or `--square` for
`--episodes` episodes (`--domain lpaworld`, a line per episode with the
counters of the incremental searches, for A* and HDA* from scratch, LPA*,
HD-LPA*, D* Lite and HD-D* Lite),
with the edges of any domain made slow (`--edges`), on the permutation puzzles,
the Towers of Hanoi and the road graphs (e.g., `--domain hanoi --file
data/hanoi/toh20.txt --pdb 8+12 --algo pida --threads 8 --table 6`, or
`--domain rubik --file data/rubik/korf10.txt --pdb korf --cache pdb`), and it
prints a line of CSV per run with the cost, the counters, the time, the
evaluations of the heuristic (or of the edges) and the iterations (the greedy
and the other searches of the section above have their options there, e.g.,
`--algo gbfs --eval lazy --queues best:0,pref:0 --boost 1000` for the first
phase of LAMA, or `--algo pgbfs --rule obat --sge 1`); `BranchAndBoundBench`
does the same for every scheme of the branch-and-bound, on random 0-1 knapsacks
(`include/PBeFS/bnb/Knapsack.h`), on the traveling salesperson problem (the
tree of the partial tours, bounded by the minimum spanning tree) and on the
domains of the layer (the trees of UTS, the DIMACS graphs in `data/dimacs`, the
instances of QAPLIB in `data/qaplib`, those of Taillard and random subgraph
problems), and `AStarBenchMPI` and
`BranchAndBoundBenchMPI` are their MPI counterparts, to be started with
`mpirun`. If the environment variable `PBEFS_PIN` is set, the i-th worker is
pinned to the i-th CPU of the affinity mask of the process (as given, e.g., by
`numactl --physcpubind` or `taskset`). On machines with several sockets this is
usually needed for the timings to be repeatable. Finally,
`tools/bench/structures.sh` runs sequential A*, and HDA* and SPA* on 1 to 16
threads, of four builds, one per set of data structures, on instances of
growing size (the 15-puzzle in bands of the instances of Korf, and random grids
of growing side), and `tools/bench/structures.py` makes the tables of the
comparison, with the noise of two runs of the same build. In both campaigns
(`campaign.sh` and `structures.sh`) each run is a process of its own, pinned
with its workers to given cores of one NUMA node and with its memory on that
node, and the runs that the tables compare are made twice, in two processes, so
that the gap between the two (the A/A noise, which `tables.py` and
`structures.py` report with each table) says how far two timings can differ by
the machine alone; a campaign holds a lock on each of its cores, so that a
second one started on the same cores waits for it, and `tools/bench/launch.sh`
starts it detached, saving its PID. The timings are meant to be taken on a
machine at rest; a virtual machine whose memory is given and taken back while
it runs (such as a VM of Hyper-V with dynamic memory, which may stay at about
2.4 GB for hours and swap on any large instance) gives timings that are not
comparable with each other, and is not to be used for them.

The directory `data` holds the instances of the benchmark, each set with its
origin and license in the header of its files or in a file `NOTICE` beside
them: besides the 100 instances of Korf, three of them (`korf_1`, `korf_12`,
`korf_55`) in the format of Burns et al., which gives the cell of each tile;
the 50 instances of the 24-puzzle of the code of HDA* of Jinnai and Fukunaga
[[5]](#ref-5) (`instances50.txt`) and the 100 of their benchmark set
(`instances100.txt`, eight of which are also in `fromkorf8.txt`); four boards
of Burns et al. and one of the random boards of $500 \times 500$ cells, each
free with probability 0.7, of the same benchmark set (`data/grid`); ten of its
instances of 40 random cities of the unit square (`data/tsp`), of which
`--size` takes the first $n$; and the sequences of BAliBASE 3.0 [[45]](#ref-45)
of the same set and of the code of HDA*, and the 15 instances of BAliBASE 3.0
of the experiments of [[43]](#ref-43) and of [[5]](#ref-5) (`data/msa`). The
larger boards of that set ($1000 \times 1000$ to $5000 \times 5000$ cells) are
left out: they are made by its generator `grid-gen/grid_gen.cc`, which draws
each cell free with the given probability and the start and the goal in the
first and in the last third of each side, and `AStarBench --domain board
--origin top` reads them as they are written.

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
[[7]](#ref-7), Dijkstra's algorithm on random grids and on the boards, the
dynamic program on the whole lattice of an alignment, the one of Held and
Karp on the cities, the distances on the map of Romania, and the Pareto
fronts of small graphs and grids by the enumeration of their simple paths
(`test/moa.cpp`); for each
solution it checks both the cost and the path,
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

## Code style

The code follows the style of the repository
[code-style](https://github.com/dmeoli/code-style), i.e., its `.clang-format`
with what its `format.py` adds to it (e.g., `return( x );` and a space before
each comma), which clang-format has no option for, in its variant `spaces`,
which `.code-format` chooses: the indentation and the continuations by spaces
alone, and a space after a cast (`( void ) x`). On a pull request, the
continuous integration checks the lines that the request changes, and posts
the changes that put them in the style as suggestions on them; with that
repository cloned in `~/code-style`,

```sh
~/code-style/code-format --since origin/main FILE...
```

does the same before a commit (clang-format 22 is needed, e.g., `pip
install clang-format==22.1.5`).

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

<a id="ref-8"></a>[8] N. J. Nilsson, *Principles of Artificial Intelligence*, Tioga, 1980; reprinted by [Springer, 1982](https://doi.org/10.1007/978-3-662-09438-9).

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

<a id="ref-35"></a>[35] E. W. Dijkstra and C. S. Scholten, [*Termination detection for diffusing computations*](https://doi.org/10.1016/0020-0190(80)90021-6), Information Processing Letters 11(1), 1-4, 1980.

<a id="ref-36"></a>[36] D. Hendler and N. Shavit, [*Non-blocking steal-half work queues*](https://doi.org/10.1145/571825.571876), Proceedings of the ACM Symposium on Principles of Distributed Computing (PODC), 280-289, 2002.

<a id="ref-37"></a>[37] H. Rihani, P. Sanders and R. Dementiev, [*MultiQueues: simple relaxed concurrent priority queues*](https://doi.org/10.1145/2755573.2755616), Proceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), 80-82, 2015.

<a id="ref-38"></a>[38] A. Postnikova, N. Koval, G. Nadiradze and D. Alistarh, [*Multi-queues can be state-of-the-art priority schedulers*](https://doi.org/10.1145/3503221.3508432), Proceedings of the ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP), 353-367, 2022.

<a id="ref-39"></a>[39] S. Edelkamp, S. Jabbar and S. Schrödl, [*External A\**](https://doi.org/10.1007/978-3-540-30221-6_18), KI 2004: Advances in Artificial Intelligence, Lecture Notes in Computer Science 3238, 226-240, 2004.

<a id="ref-40"></a>[40] R. E. Korf, [*Linear-time disk-based implicit graph search*](https://doi.org/10.1145/1455248.1455250), Journal of the ACM 55(6), 26:1-26:40, 2008.

<a id="ref-41"></a>[41] T. Ikeda and H. Imai, [*Enhanced A\* algorithms for multiple alignments: optimal alignments for several sequences and k-opt approximate alignments for large cases*](https://doi.org/10.1016/S0304-3975(98)00093-0), Theoretical Computer Science 210(2), 341-374, 1999.

<a id="ref-42"></a>[42] T. Yoshizumi, T. Miura and T. Ishida, [*A\* with partial expansion for large branching factor problems*](https://mlanthology.org/aaai/2000/yoshizumi2000aaai-partial), Proceedings of the AAAI Conference on Artificial Intelligence 17, 923-929, 2000.

<a id="ref-43"></a>[43] Y. Kobayashi, A. Kishimoto and O. Watanabe, [*Evaluations of Hash Distributed A\* in optimal sequence alignment*](https://doi.org/10.5591/978-1-57735-516-8/IJCAI11-105), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 22, 584-589, 2011.

<a id="ref-44"></a>[44] M. O. Dayhoff, R. M. Schwartz and B. C. Orcutt, *A model of evolutionary change in proteins*, Atlas of Protein Sequence and Structure 5, supplement 3, 345-352, 1978.

<a id="ref-45"></a>[45] J. D. Thompson, P. Koehl, R. Ripp and O. Poch, [*BAliBASE 3.0: latest developments of the multiple sequence alignment benchmark*](https://doi.org/10.1002/prot.20527), Proteins 61(1), 127-136, 2005.

<a id="ref-46"></a>[46] I. Pohl, [*Heuristic search viewed as path finding in a graph*](https://doi.org/10.1016/0004-3702(70)90007-X), Artificial Intelligence 1(3-4), 193-204, 1970.

<a id="ref-47"></a>[47] E. A. Hansen and R. Zhou, [*Anytime heuristic search*](https://doi.org/10.1613/jair.2096), Journal of Artificial Intelligence Research 28, 267-297, 2007.

<a id="ref-48"></a>[48] M. Likhachev, G. Gordon and S. Thrun, [*ARA\*: Anytime A\* with provable bounds on sub-optimality*](https://papers.nips.cc/paper/2003/file/ee8fe9093fbbb687bef15a38facc44d2-Paper.pdf), Advances in Neural Information Processing Systems (NIPS) 16, 767-774, 2003.

<a id="ref-49"></a>[49] J. T. Thayer and W. Ruml, [*Faster than weighted A\*: an optimistic approach to bounded suboptimal search*](https://www.cs.unh.edu/~ruml/papers/optimistic-icaps-08.pdf), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 18, 355-362, 2008.

<a id="ref-50"></a>[50] J. Chen, N. R. Sturtevant, W. Doyle and W. Ruml, [*Revisiting suboptimal search*](https://doi.org/10.1609/socs.v10i1.18498), Proceedings of the Symposium on Combinatorial Search (SoCS) 10, 18-25, 2019.

<a id="ref-51"></a>[51] J. Pearl and J. H. Kim, [*Studies in semi-admissible heuristics*](https://doi.org/10.1109/TPAMI.1982.4767270), IEEE Transactions on Pattern Analysis and Machine Intelligence 4(4), 392-399, 1982.

<a id="ref-52"></a>[52] J. T. Thayer and W. Ruml, [*Bounded suboptimal search: a direct approach using inadmissible estimates*](https://doi.org/10.5591/978-1-57735-516-8/IJCAI11-119), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 22, 674-679, 2011.

<a id="ref-53"></a>[53] D. Gilon, A. Felner and R. Stern, [*Dynamic potential search: a new bounded suboptimal search algorithm*](https://doi.org/10.1609/socs.v7i1.18392), Proceedings of the Symposium on Combinatorial Search (SoCS) 7, 36-44, 2016.

<a id="ref-54"></a>[54] S. Koenig, M. Likhachev and D. Furcy, [*Lifelong Planning A\**](https://doi.org/10.1016/j.artint.2003.12.001), Artificial Intelligence 155(1-2), 93-146, 2004.

<a id="ref-55"></a>[55] A. Reinefeld and V. Schnecke, *AIDA\*: asynchronous parallel IDA\**, Proceedings of the Biennial Conference of the Canadian Society for Computational Studies of Intelligence (AI 94), 295-302, 1994.

<a id="ref-56"></a>[56] Z. Bu and R. E. Korf, [*A\*+IDA\*: a simple hybrid search algorithm*](https://doi.org/10.24963/ijcai.2019/168), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 28, 1206-1212, 2019.

<a id="ref-57"></a>[57] L. Siag, S. S. Shperberg, A. Felner and N. R. Sturtevant, [*On parallel external-memory bidirectional search*](https://doi.org/10.3233/FAIA240991), Proceedings of the European Conference on Artificial Intelligence (ECAI 2024), Frontiers in Artificial Intelligence and Applications 392, 4190-4197, 2024.

<a id="ref-58"></a>[58] V. N. Rao and V. Kumar, [*Parallel depth first search. Part I. Implementation*](https://doi.org/10.1007/BF01389000), International Journal of Parallel Programming 16(6), 479-499, 1987.

<a id="ref-59"></a>[59] R. E. Korf, [*Finding optimal solutions to Rubik's Cube using pattern databases*](https://cdn.aaai.org/AAAI/1997/AAAI97-109.pdf), Proceedings of the AAAI Conference on Artificial Intelligence 14, 700-705, 1997.

<a id="ref-60"></a>[60] M. Helmert, [*Landmark heuristics for the pancake problem*](https://doi.org/10.1609/socs.v1i1.18176), Proceedings of the Symposium on Combinatorial Search (SoCS) 1, 109-110, 2010.

<a id="ref-61"></a>[61] R. C. Holte, A. Felner, G. Sharon and N. R. Sturtevant, [*Bidirectional search that is guaranteed to meet in the middle*](https://doi.org/10.1609/aaai.v30i1.10436), Proceedings of the AAAI Conference on Artificial Intelligence 30, 2016.

<a id="ref-62"></a>[62] B. Bonet, *Efficient algorithms to rank and unrank permutations in lexicographic order*, AAAI Workshop on Search in Artificial Intelligence and Robotics, Technical Report WS-08-10, 18-23, 2008.

<a id="ref-63"></a>[63] H. Dweighter, [*Problem E2569*](https://doi.org/10.2307/2318260), in *Elementary Problems: E2564-E2569*, The American Mathematical Monthly 82(10), 1009, 1975.

<a id="ref-64"></a>[64] C. Demetrescu, A. V. Goldberg and D. S. Johnson (editors), [*The Shortest Path Problem: Ninth DIMACS Implementation Challenge*](https://doi.org/10.1090/dimacs/074), DIMACS Series in Discrete Mathematics and Theoretical Computer Science 74, AMS, 2009.

<a id="ref-65"></a>[65] R. B. Dial, [*Algorithm 360: shortest-path forest with topological ordering*](https://doi.org/10.1145/363269.363610), Communications of the ACM 12(11), 632-633, 1969.

<a id="ref-66"></a>[66] R. L. Graham, [*Bounds on multiprocessing timing anomalies*](https://doi.org/10.1137/0117039), SIAM Journal on Applied Mathematics 17(2), 416-429, 1969.

<a id="ref-67"></a>[67] B. Archibald, P. Maier, R. Stewart and P. Trinder, [*Replicable parallel branch and bound search*](https://doi.org/10.1016/j.jpdc.2017.10.010), Journal of Parallel and Distributed Computing 113, 92-114, 2018.

<a id="ref-68"></a>[68] B. Archibald, P. Maier, R. Stewart and P. Trinder, [*YewPar: skeletons for exact combinatorial search*](https://doi.org/10.1145/3332466.3374537), Proceedings of the ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP) 25, 292-307, 2020.

<a id="ref-69"></a>[69] Y. Shinano, T. Achterberg, T. Berthold, S. Heinz, T. Koch and M. Winkler, [*Solving open MIP instances with ParaSCIP on supercomputers using up to 80,000 cores*](https://doi.org/10.1109/IPDPS.2016.56), Proceedings of the IEEE International Parallel and Distributed Processing Symposium (IPDPS) 30, 770-779, 2016.

<a id="ref-70"></a>[70] Y. Xu, T. K. Ralphs, L. Ladányi and M. J. Saltzman, [*ALPS: a framework for implementing parallel tree search algorithms*](https://doi.org/10.1007/0-387-23529-9_21), The Next Wave in Computing, Optimization, and Decision Technologies, 319-334, Springer, 2005.

<a id="ref-71"></a>[71] M. Mezmaz, R. Leroy, N. Melab and D. Tuyttens, [*A multi-core parallel branch-and-bound algorithm using factorial number system*](https://doi.org/10.1109/IPDPS.2014.124), Proceedings of the IEEE International Parallel and Distributed Processing Symposium (IPDPS) 28, 1203-1212, 2014.

<a id="ref-72"></a>[72] J. Gmys, R. Leroy, M. Mezmaz, N. Melab and D. Tuyttens, [*Work stealing with private integer-vector-matrix data structure for multi-core branch-and-bound algorithms*](https://doi.org/10.1002/cpe.3771), Concurrency and Computation: Practice and Experience 28(18), 4463-4484, 2016.

<a id="ref-73"></a>[73] S. Olivier, J. Huan, J. Liu, J. Prins, J. Dinan, P. Sadayappan and C.-W. Tseng, [*UTS: an unbalanced tree search benchmark*](https://doi.org/10.1007/978-3-540-72521-3_18), Languages and Compilers for Parallel Computing (LCPC) 19, Lecture Notes in Computer Science 4382, 235-250, 2007.

<a id="ref-74"></a>[74] E. Tomita and T. Kameda, [*An efficient branch-and-bound algorithm for finding a maximum clique with computational experiments*](https://doi.org/10.1007/s10898-006-9039-7), Journal of Global Optimization 37(1), 95-111, 2007.

<a id="ref-75"></a>[75] P. San Segundo, D. Rodríguez-Losada and A. Jiménez, [*An exact bit-parallel algorithm for the maximum clique problem*](https://doi.org/10.1016/j.cor.2010.07.019), Computers & Operations Research 38(2), 571-581, 2011.

<a id="ref-76"></a>[76] C. McCreesh, P. Prosser and J. Trimble, [*A partitioning algorithm for maximum common subgraph problems*](https://doi.org/10.24963/ijcai.2017/99), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 26, 712-719, 2017.

<a id="ref-77"></a>[77] P. C. Gilmore, [*Optimal and suboptimal algorithms for the quadratic assignment problem*](https://doi.org/10.1137/0110022), Journal of the Society for Industrial and Applied Mathematics 10(2), 305-313, 1962.

<a id="ref-78"></a>[78] R. E. Burkard, S. E. Karisch and F. Rendl, [*QAPLIB - a quadratic assignment problem library*](https://doi.org/10.1023/A:1008293323270), Journal of Global Optimization 10(4), 391-403, 1997.

<a id="ref-79"></a>[79] B. J. Lageweg, J. K. Lenstra and A. H. G. Rinnooy Kan, [*A general bounding scheme for the permutation flow-shop problem*](https://doi.org/10.1287/opre.26.1.53), Operations Research 26(1), 53-67, 1978.

<a id="ref-80"></a>[80] E. Taillard, [*Benchmarks for basic scheduling problems*](https://doi.org/10.1016/0377-2217(93)90182-M), European Journal of Operational Research 64(2), 278-285, 1993.

<a id="ref-81"></a>[81] R. C. Holte, A. Felner, G. Sharon, N. R. Sturtevant and J. Chen, [*MM: A bidirectional search algorithm that is guaranteed to meet in the middle*](https://doi.org/10.1016/j.artint.2017.05.004), Artificial Intelligence 252, 232-266, 2017.

<a id="ref-82"></a>[82] G. Sharon, R. C. Holte, A. Felner and N. R. Sturtevant, [*Extended abstract: an improved priority function for bidirectional heuristic search*](https://doi.org/10.1609/socs.v7i1.18399), Proceedings of the Symposium on Combinatorial Search (SoCS) 7, 139-140, 2016.

<a id="ref-83"></a>[83] E. Shaham, A. Felner, J. Chen and N. R. Sturtevant, [*The minimal set of states that must be expanded in a front-to-end bidirectional search*](https://doi.org/10.1609/socs.v8i1.18426), Proceedings of the Symposium on Combinatorial Search (SoCS) 8, 82-90, 2017.

<a id="ref-84"></a>[84] J. Eckerle, J. Chen, N. R. Sturtevant, S. Zilles and R. C. Holte, [*Sufficient conditions for node expansion in bidirectional heuristic search*](https://doi.org/10.1609/icaps.v27i1.13825), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 27, 79-87, 2017.

<a id="ref-85"></a>[85] J. Chen, R. C. Holte, S. Zilles and N. R. Sturtevant, [*Front-to-end bidirectional heuristic search with near-optimal node expansions*](https://doi.org/10.24963/ijcai.2017/69), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI), 489-495, 2017.

<a id="ref-86"></a>[86] S. S. Shperberg, A. Felner, N. R. Sturtevant, S. E. Shimony and A. Hayoun, [*Enriching non-parametric bidirectional search algorithms*](https://doi.org/10.1609/aaai.v33i01.33012379), Proceedings of the AAAI Conference on Artificial Intelligence 33, 2379-2386, 2019.

<a id="ref-87"></a>[87] M. Barley, P. Riddle, C. Linares López, S. Dobson and I. Pohl, [*GBFHS: A generalized breadth-first heuristic search algorithm*](https://doi.org/10.1609/socs.v9i1.18452), Proceedings of the Symposium on Combinatorial Search (SoCS) 9, 28-36, 2018.

<a id="ref-88"></a>[88] S. K. Sadhukhan, *Bidirectional heuristic search based on error estimate*, CSI Journal of Computing 2(1-2), S1, 2013.

<a id="ref-89"></a>[89] V. Alcázar, P. J. Riddle and M. Barley, [*A unifying view on individual bounds and heuristic inaccuracies in bidirectional search*](https://doi.org/10.1609/aaai.v34i03.5611), Proceedings of the AAAI Conference on Artificial Intelligence 34, 2327-2334, 2020.

<a id="ref-90"></a>[90] J. B. H. Kwa, [*BS\*: An admissible bidirectional staged heuristic search algorithm*](https://doi.org/10.1016/0004-3702(89)90069-6), Artificial Intelligence 38(1), 95-109, 1989.

<a id="ref-91"></a>[91] I. Pohl, *Bi-directional search*, Machine Intelligence 6, 127-140, 1971.

<a id="ref-92"></a>[92] W. Pijls and H. Post, [*A new bidirectional search algorithm with shortened postprocessing*](https://doi.org/10.1016/j.ejor.2008.09.032), European Journal of Operational Research 198(2), 363-369, 2009.

<a id="ref-93"></a>[93] L. H. O. Rios and L. Chaimowicz, *PNBA\*: A parallel bidirectional heuristic search algorithm*, Proceedings of the Encontro Nacional de Inteligência Artificial (ENIA) 8, 2011.

<a id="ref-94"></a>[94] N. R. Sturtevant and J. Chen, [*External memory bidirectional search*](https://www.ijcai.org/Abstract/16/102), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI), 676-682, 2016.

<a id="ref-95"></a>[95] S. Lin and A. Fukunaga, [*Revisiting immediate duplicate detection in external memory search*](https://doi.org/10.1609/aaai.v32i1.11526), Proceedings of the AAAI Conference on Artificial Intelligence 32, 2018.

<a id="ref-96"></a>[96] J. E. Doran and D. Michie, [*Experiments with the Graph Traverser program*](https://doi.org/10.1098/rspa.1966.0205), Proceedings of the Royal Society A 294, 235-259, 1966.

<a id="ref-97"></a>[97] S. Richter and M. Helmert, [*Preferred operators and deferred evaluation in satisficing planning*](https://doi.org/10.1609/icaps.v19i1.13345), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 19, 273-280, 2009.

<a id="ref-98"></a>[98] G. Röger and M. Helmert, [*The more, the merrier: combining heuristic estimators for satisficing planning*](https://doi.org/10.1609/icaps.v20i1.13435), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 20, 246-249, 2010.

<a id="ref-99"></a>[99] F. Xie, M. Müller, R. Holte and T. Imai, [*Type-based exploration with multiple search queues for satisficing planning*](https://doi.org/10.1609/aaai.v28i1.9036), Proceedings of the AAAI Conference on Artificial Intelligence 28, 2014.

<a id="ref-100"></a>[100] R. A. Valenzano, N. R. Sturtevant, J. Schaeffer and F. Xie, [*A comparison of knowledge-based GBFS enhancements and knowledge-free exploration*](https://doi.org/10.1609/icaps.v24i1.13681), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 24, 375-379, 2014.

<a id="ref-101"></a>[101] J. Hoffmann and B. Nebel, [*The FF planning system: fast plan generation through heuristic search*](https://doi.org/10.1613/jair.855), Journal of Artificial Intelligence Research 14, 253-302, 2001.

<a id="ref-102"></a>[102] S. Richter and M. Westphal, [*The LAMA planner: guiding cost-based anytime planning with landmarks*](https://doi.org/10.1613/jair.2972), Journal of Artificial Intelligence Research 39, 127-177, 2010.

<a id="ref-103"></a>[103] S. Richter, J. T. Thayer and W. Ruml, [*The joy of forgetting: faster anytime search via restarting*](https://doi.org/10.1609/icaps.v20i1.13412), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 20, 137-144, 2010.

<a id="ref-104"></a>[104] N. Lipovetzky and H. Geffner, [*Width and serialization of classical planning problems*](https://doi.org/10.3233/978-1-61499-098-7-540), Proceedings of the European Conference on Artificial Intelligence (ECAI 2012), Frontiers in Artificial Intelligence and Applications 242, 540-545, 2012.

<a id="ref-105"></a>[105] N. Lipovetzky and H. Geffner, [*Best-first width search: exploration and exploitation in classical planning*](https://doi.org/10.1609/aaai.v31i1.11027), Proceedings of the AAAI Conference on Artificial Intelligence 31, 2017.

<a id="ref-106"></a>[106] V. Vidal, L. Bordeaux and Y. Hamadi, [*Adaptive K-parallel best-first search: a simple but efficient algorithm for multi-core domain-independent planning*](https://doi.org/10.1609/socs.v1i1.18165), Proceedings of the Symposium on Combinatorial Search (SoCS) 1, 100-107, 2010.

<a id="ref-107"></a>[107] R. Kuroiwa and A. Fukunaga, [*On the pathological search behavior of distributed greedy best-first search*](https://doi.org/10.1609/icaps.v29i1.3485), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 29, 255-263, 2019.

<a id="ref-108"></a>[108] R. Kuroiwa and A. Fukunaga, [*Analyzing and avoiding pathological behavior in parallel best-first search*](https://doi.org/10.1609/icaps.v30i1.6659), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 30, 175-183, 2020.

<a id="ref-109"></a>[109] T. Shimoda and A. Fukunaga, [*Improved exploration of the bench transition system in parallel greedy best first search*](https://doi.org/10.1609/socs.v16i1.27285), Proceedings of the Symposium on Combinatorial Search (SoCS) 16, 74-82, 2023.

<a id="ref-110"></a>[110] T. Shimoda and A. Fukunaga, [*Parallel greedy best-first search with a bound on expansions relative to sequential search*](https://doi.org/10.1609/aaai.v39i25.34869), Proceedings of the AAAI Conference on Artificial Intelligence 39, 26668-26677, 2025.

<a id="ref-111"></a>[111] T. Shimoda and A. Fukunaga, [*Decoupling generation and evaluation for parallel greedy best-first search*](https://arxiv.org/abs/2408.05682), arXiv:2408.05682, 2024.

<a id="ref-112"></a>[112] M. Heusner, T. Keller and M. Helmert, [*Understanding the search behaviour of greedy best-first search*](https://doi.org/10.1609/socs.v8i1.18425), Proceedings of the Symposium on Combinatorial Search (SoCS) 8, 47-55, 2017.

<a id="ref-113"></a>[113] A. Felner, M. Goldenberg, G. Sharon, R. Stern, T. Beja, N. R. Sturtevant, J. Schaeffer and R. C. Holte, [*Partial-expansion A\* with selective node generation*](https://doi.org/10.1609/aaai.v26i1.8137), Proceedings of the AAAI Conference on Artificial Intelligence 26, 471-477, 2012.

<a id="ref-114"></a>[114] M. Goldenberg, A. Felner, R. Stern, G. Sharon, N. Sturtevant, R. C. Holte and J. Schaeffer, [*Enhanced partial expansion A\**](https://doi.org/10.1613/jair.4171), Journal of Artificial Intelligence Research 50, 141-187, 2014.

<a id="ref-115"></a>[115] Y. Björnsson, M. Enzenberger, R. C. Holte and J. Schaeffer, *Fringe search: beating A\* at pathfinding on game maps*, Proceedings of the IEEE Symposium on Computational Intelligence and Games (CIG 2005), 125-132, 2005.

<a id="ref-116"></a>[116] R. E. Korf, [*Linear-space best-first search*](https://doi.org/10.1016/0004-3702(93)90045-D), Artificial Intelligence 62(1), 41-78, 1993.

<a id="ref-117"></a>[117] R. Zhou and E. A. Hansen, *Beam-stack search: integrating backtracking with beam search*, Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 15, 90-98, 2005.

<a id="ref-118"></a>[118] R. E. Korf, [*Real-time heuristic search*](https://doi.org/10.1016/0004-3702(90)90054-4), Artificial Intelligence 42(2-3), 189-211, 1990.

<a id="ref-119"></a>[119] S. Koenig and X. Sun, [*Comparing real-time and incremental heuristic search for real-time situated agents*](https://doi.org/10.1007/s10458-008-9061-x), Autonomous Agents and Multi-Agent Systems 18(3), 313-341, 2009.

<a id="ref-120"></a>[120] A. Botea, M. Müller and J. Schaeffer, *Near optimal hierarchical path-finding*, Journal of Game Development 1(1), 7-28, 2004.

<a id="ref-121"></a>[121] S. Mukherjee, S. Aine and M. Likhachev, [*ePA\*SE: Edge-based parallel A\* for slow evaluations*](https://doi.org/10.1609/socs.v15i1.21761), Proceedings of the International Symposium on Combinatorial Search (SoCS) 15, 136-144, 2022.

<a id="ref-122"></a>[122] S. Mukherjee and M. Likhachev, [*GePA\*SE: Generalized edge-based parallel A\* for slow evaluations*](https://doi.org/10.1609/socs.v16i1.27295), Proceedings of the International Symposium on Combinatorial Search (SoCS) 16, 153-157, 2023.

<a id="ref-123"></a>[123] H. Yang, S. Mukherjee and M. Likhachev, [*A-ePA\*SE: Anytime edge-based parallel A\* for slow evaluations*](https://doi.org/10.1609/socs.v16i1.27297), Proceedings of the International Symposium on Combinatorial Search (SoCS) 16, 163-167, 2023.

<a id="ref-124"></a>[124] S. Mukherjee, S. Aine and M. Likhachev, [*MPLP: Massively parallelized lazy planning*](https://doi.org/10.1109/LRA.2022.3157544), IEEE Robotics and Automation Letters 7(3), 6067-6074, 2022.

<a id="ref-125"></a>[125] C. M. Dellin and S. S. Srinivasa, [*A unifying formalism for shortest path problems with expensive edge evaluations via lazy best-first search over paths with edge selectors*](https://doi.org/10.1609/icaps.v26i1.13788), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 26, 459-467, 2016.

<a id="ref-126"></a>[126] F. Agostinelli, S. McAleer, A. Shmakov and P. Baldi, [*Solving the Rubik's cube with deep reinforcement learning and search*](https://doi.org/10.1038/s42256-019-0070-z), Nature Machine Intelligence 1, 356-363, 2019.

<a id="ref-127"></a>[127] M. Evett, J. Hendler, A. Mahanti and D. Nau, [*PRA\*: massively parallel heuristic search*](https://doi.org/10.1006/jpdc.1995.1036), Journal of Parallel and Distributed Computing 25(2), 133-143, 1995.

<a id="ref-128"></a>[128] R. Zhou and E. A. Hansen, [*Parallel structured duplicate detection*](https://mlanthology.org/aaai/2007/zhou2007aaai-parallel), Proceedings of the AAAI Conference on Artificial Intelligence 22, 1217-1223, 2007.

<a id="ref-129"></a>[129] R. Zhou and E. A. Hansen, [*Breadth-first heuristic search*](https://doi.org/10.1016/j.artint.2005.12.002), Artificial Intelligence 170(4-5), 385-408, 2006.

<a id="ref-130"></a>[130] A. Felner, S. Kraus and R. E. Korf, [*KBFS: K-best-first search*](https://doi.org/10.1023/A:1024452529781), Annals of Mathematics and Artificial Intelligence 39(1-2), 19-39, 2003.

<a id="ref-131"></a>[131] J. W. Romein, A. Plaat, H. E. Bal and J. Schaeffer, [*Transposition table driven work scheduling in distributed search*](https://mlanthology.org/aaai/1999/romein1999aaai-transposition), Proceedings of the AAAI Conference on Artificial Intelligence 16, 725-731, 1999.

<a id="ref-132"></a>[132] U. K. Sarkar, P. P. Chakrabarti, S. Ghose and S. C. De Sarkar, [*Reducing reexpansions in iterative-deepening search by controlling cutoff bounds*](https://doi.org/10.1016/0004-3702(91)90100-X), Artificial Intelligence 50(2), 207-221, 1991.

<a id="ref-133"></a>[133] C. Bäckström and B. Nebel, [*Complexity results for SAS+ planning*](https://doi.org/10.1111/j.1467-8640.1995.tb00052.x), Computational Intelligence 11(4), 625-655, 1995.

<a id="ref-134"></a>[134] M. Helmert, [*The Fast Downward planning system*](https://doi.org/10.1613/jair.1705), Journal of Artificial Intelligence Research 26, 191-246, 2006.

<a id="ref-135"></a>[135] B. Bonet and H. Geffner, [*Planning as heuristic search*](https://doi.org/10.1016/S0004-3702(01)00108-4), Artificial Intelligence 129(1-2), 5-33, 2001.

<a id="ref-136"></a>[136] M. Helmert and C. Domshlak, [*Landmarks, critical paths and abstractions: what's the difference anyway?*](https://doi.org/10.1609/icaps.v19i1.13370), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 19, 162-169, 2009.

<a id="ref-137"></a>[137] S. Edelkamp, *Planning with pattern databases*, Proceedings of the European Conference on Planning (ECP), 13-24, 2001.

<a id="ref-138"></a>[138] Y. Jinnai and A. Fukunaga, [*Automated creation of efficient work distribution functions for parallel best-first search*](https://doi.org/10.1609/icaps.v26i1.13747), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 26, 184-192, 2016.

<a id="ref-139"></a>[139] R. Zhou and E. A. Hansen, [*Domain-independent structured duplicate detection*](https://cdn.aaai.org/AAAI/2006/AAAI06-170.pdf), Proceedings of the AAAI Conference on Artificial Intelligence 21, 1082-1087, 2006.

<a id="ref-140"></a>[140] M. Helmert, P. Haslum, J. Hoffmann and R. Nissim, [*Merge-and-shrink abstraction: a method for generating lower bounds in factored state spaces*](https://doi.org/10.1145/2559951), Journal of the ACM 61(3), 16:1-16:63, 2014.

<a id="ref-141"></a>[141] Y. Jinnai and A. Fukunaga, [*Abstract Zobrist hashing: an efficient work distribution method for parallel best-first search*](https://doi.org/10.1609/aaai.v30i1.10065), Proceedings of the AAAI Conference on Artificial Intelligence 30, 717-723, 2016.

<a id="ref-142"></a>[142] C. M. Fiduccia and R. M. Mattheyses, [*A linear-time heuristic for improving network partitions*](https://doi.org/10.1109/DAC.1982.1585498), Proceedings of the 19th Design Automation Conference, 175-181, 1982.

<a id="ref-143"></a>[143] S. Koenig and M. Likhachev, [*D\* Lite*](https://aaai.org/papers/00476-aaai02-072-d-lite/), Proceedings of the AAAI Conference on Artificial Intelligence (AAAI) 18, 476-483, 2002; and [*Fast replanning for navigation in unknown terrain*](https://doi.org/10.1109/TRO.2004.838026), IEEE Transactions on Robotics 21(3), 354-363, 2005.

<a id="ref-144"></a>[144] G. Ramalingam and T. Reps, [*An incremental algorithm for a generalization of the shortest-path problem*](https://doi.org/10.1006/jagm.1996.0046), Journal of Algorithms 21(2), 267-305, 1996.

<a id="ref-145"></a>[145] U. Meyer and P. Sanders, [*Δ-stepping: a parallelizable shortest path algorithm*](https://doi.org/10.1016/S0196-6774(03)00076-2), Journal of Algorithms 49(1), 114-152, 2003.

<a id="ref-146"></a>[146] A. H. Land and A. G. Doig, [*An automatic method of solving discrete programming problems*](https://doi.org/10.2307/1910129), Econometrica 28(3), 497-520, 1960.

<a id="ref-147"></a>[147] G. B. Dantzig, [*Discrete-variable extremum problems*](https://doi.org/10.1287/opre.5.2.266), Operations Research 5(2), 266-288, 1957.

<a id="ref-148"></a>[148] D. Pisinger, [*Where are the hard knapsack problems?*](https://doi.org/10.1016/j.cor.2004.03.002), Computers & Operations Research 32(9), 2271-2284, 2005.

<a id="ref-149"></a>[149] T.-H. Lai and S. Sahni, [*Anomalies in parallel branch-and-bound algorithms*](https://doi.org/10.1145/358080.358103), Communications of the ACM 27(6), 594-602, 1984.

<a id="ref-150"></a>[150] W. W. Johnson and W. E. Story, [*Notes on the "15" Puzzle*](https://doi.org/10.2307/2369492), American Journal of Mathematics 2(4), 397-404, 1879.

<a id="ref-151"></a>[151] D. Lemire, [*Fast random integer generation in an interval*](https://doi.org/10.1145/3230636), ACM Transactions on Modeling and Computer Simulation 29(1), 1-12, 2019.

<a id="ref-152"></a>[152] A. L. Zobrist, *A new hashing method with application for game playing*, Technical Report 88, Computer Sciences Department, University of Wisconsin, 1970; reprinted in [ICGA Journal 13(2), 69-73, 1990](https://doi.org/10.3233/ICG-1990-13203).

<a id="ref-153"></a>[153] V. Kumar, K. Ramesh and V. N. Rao, [*Parallel best-first search of state-space graphs: a summary of results*](https://cdn.aaai.org/AAAI/1988/AAAI88-022.pdf), Proceedings of the AAAI Conference on Artificial Intelligence 7, 122-127, 1988.

<a id="ref-154"></a>[154] R. M. Karp and Y. Zhang, [*Randomized parallel algorithms for backtrack search and branch-and-bound computation*](https://doi.org/10.1145/174130.174145), Journal of the ACM 40(3), 765-789, 1993.

<a id="ref-155"></a>[155] K. B. Irani and Y. Shih, *Parallel A\* and AO\* algorithms: an optimality criterion and performance evaluation*, Proceedings of the International Conference on Parallel Processing (ICPP), 274-277, 1986.

<a id="ref-156"></a>[156] R. E. Korf and P. Schultze, [*Large-scale parallel breadth-first search*](https://cdn.aaai.org/AAAI/2005/AAAI05-219.pdf), Proceedings of the AAAI Conference on Artificial Intelligence 20, 1380-1385, 2005.

<a id="ref-157"></a>[157] D. E. Knuth, [*The Art of Computer Programming, Volume 3: Sorting and Searching*](https://www-cs-faculty.stanford.edu/~knuth/taocp.html), 2nd edition, Addison-Wesley, 1998.

<a id="ref-158"></a>[158] D. B. Johnson, [*Priority queues with update and finding minimum spanning trees*](https://doi.org/10.1016/0020-0190(75)90001-0), Information Processing Letters 4(3), 53-57, 1975.

<a id="ref-159"></a>[159] D. R. Hanson, [*Fast allocation and deallocation of memory based on object lifetimes*](https://doi.org/10.1002/spe.4380200104), Software: Practice and Experience 20(1), 5-12, 1990.

<a id="ref-160"></a>[160] R. K. Ahuja, K. Mehlhorn, J. B. Orlin and R. E. Tarjan, [*Faster algorithms for the shortest path problem*](https://doi.org/10.1145/77600.77615), Journal of the ACM 37(2), 213-223, 1990.

<a id="ref-161"></a>[161] A. Reinefeld and T. A. Marsland, [*Enhanced iterative-deepening search*](https://doi.org/10.1109/34.297950), IEEE Transactions on Pattern Analysis and Machine Intelligence 16(7), 701-710, 1994.

<a id="ref-162"></a>[162] T. Miura and T. Ishida, [*Stochastic node caching for memory-bounded search*](https://doi.org/10.1002/(SICI)1520-684X(199902)30:2<57::AID-SCJ7>3.0.CO;2-N), Systems and Computers in Japan 30(2), 57-65, 1999.

<a id="ref-163"></a>[163] R. E. Bryant, [*Graph-based algorithms for Boolean function manipulation*](https://doi.org/10.1109/TC.1986.1676819), IEEE Transactions on Computers C-35(8), 677-691, 1986.

<a id="ref-164"></a>[164] Y. Jinnai, *Introduction to heuristic search* (in Japanese), lecture notes, 175 pages, version of 9 December 2023, available from <https://jinnaiyuu.github.io/>.

<a id="ref-165"></a>[165] N. R. Sturtevant, [*HOG2: Hierarchical Open Graph*](https://github.com/MovingAILab/hog2), software, Moving AI Lab, University of Denver.

<a id="ref-166"></a>[166] E. L. Lawler, [*The quadratic assignment problem*](https://doi.org/10.1287/mnsc.9.4.586), Management Science 9(4), 586-599, 1963.

<a id="ref-167"></a>[167] A. Felner, U. Zahavi, J. Schaeffer and R. C. Holte, *Dual lookups in pattern databases*, Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 19, 103-108, 2005.

<a id="ref-168"></a>[168] R. Nissim, J. Hoffmann and M. Helmert, *Computing perfect heuristics in polynomial time: on bisimulation and merge-and-shrink abstraction in optimal planning*, Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 22, 1983-1990, 2011.

<a id="ref-169"></a>[169] A. R. Karlin, M. S. Manasse, L. Rudolph and D. D. Sleator, [*Competitive snoopy caching*](https://doi.org/10.1007/BF01762111), Algorithmica 3(1-4), 79-119, 1988.

<a id="ref-170"></a>[170] N. J. Higham, [*Accuracy and Stability of Numerical Algorithms*](https://doi.org/10.1137/1.9780898718027), 2nd edition, SIAM, 2002.

<a id="ref-171"></a>[171] B. T. Lowerre, *The HARPY speech recognition system*, PhD thesis, Carnegie Mellon University, 1976.

<a id="ref-172"></a>[172] H. Kociemba, [*kociemba.org/cube.htm*](https://kociemba.org/cube.htm), web pages on the two-phase algorithm for Rubik's cube and its cubie model, consulted in November 2026.

<a id="ref-173"></a>[173] E. Keyder and H. Geffner, [*Heuristics for planning with action costs revisited*](https://doi.org/10.3233/978-1-58603-891-5-588), Proceedings of the European Conference on Artificial Intelligence (ECAI), Frontiers in Artificial Intelligence and Applications 178, 588-592, 2008.

<a id="ref-174"></a>[174] D. S. Johnson and M. A. Trick (editors), [*Cliques, Coloring, and Satisfiability: Second DIMACS Implementation Challenge*](https://doi.org/10.1090/dimacs/026), DIMACS Series in Discrete Mathematics and Theoretical Computer Science 26, American Mathematical Society, 1996.

<a id="ref-175"></a>[175] C. Bron and J. Kerbosch, [*Algorithm 457: finding all cliques of an undirected graph*](https://doi.org/10.1145/362342.362367), Communications of the ACM 16(9), 575-577, 1973.

<a id="ref-176"></a>[176] T. C. Koopmans and M. Beckmann, [*Assignment problems and the location of economic activities*](https://doi.org/10.2307/1907742), Econometrica 25(1), 53-76, 1957.

<a id="ref-177"></a>[177] H. W. Kuhn, [*The Hungarian method for the assignment problem*](https://doi.org/10.1002/nav.3800020109), Naval Research Logistics Quarterly 2(1-2), 83-97, 1955.

<a id="ref-178"></a>[178] S. M. Johnson, [*Optimal two- and three-stage production schedules with setup times included*](https://doi.org/10.1002/nav.3800010110), Naval Research Logistics Quarterly 1(1), 61-68, 1954.

<a id="ref-179"></a>[179] National Institute of Standards and Technology, [*Secure Hash Standard (SHS)*](https://doi.org/10.6028/NIST.FIPS.180-4), Federal Information Processing Standards Publication 180-4, 2015.

<a id="ref-180"></a>[180] L. Mandow and J. L. Pérez de la Cruz, [*Multiobjective A\* search with consistent heuristics*](https://doi.org/10.1145/1754399.1754400), Journal of the ACM 57(5), 1-25, 2010.

<a id="ref-181"></a>[181] C. Hernández, W. Yeoh, J. A. Baier, H. Zhang, L. Suazo, S. Koenig and O. Salzman, [*Simple and efficient bi-objective search algorithms via fast dominance checks*](https://doi.org/10.1016/j.artint.2022.103807), Artificial Intelligence 314, 103807, 2023.

<a id="ref-182"></a>[182] Z. Ren, R. Zhan, S. Rathinam, M. Likhachev and H. Choset, [*Enhanced multi-objective A\* using balanced binary search trees*](https://doi.org/10.1609/socs.v15i1.21764), Proceedings of the International Symposium on Combinatorial Search (SoCS) 15, 162-170, 2022.

<a id="ref-183"></a>[183] C. Hernández, W. Yeoh, J. A. Baier, A. Felner, O. Salzman, H. Zhang, S.-H. Chan and S. Koenig, [*Multi-objective search via lazy and efficient dominance checks*](https://doi.org/10.24963/ijcai.2023/850), Proceedings of the International Joint Conference on Artificial Intelligence (IJCAI) 32, 7223-7230, 2023.

<a id="ref-184"></a>[184] S. Ahmadi, N. R. Sturtevant, A. Raith, D. Harabor and M. Jalili, [*Parallelizing multi-objective A\* search*](https://doi.org/10.1609/icaps.v35i1.36109), Proceedings of the International Conference on Automated Planning and Scheduling (ICAPS) 35, 131-139, 2025.

<a id="ref-185"></a>[185] F.-J. Pulido, L. Mandow and J. L. Pérez de la Cruz, [*Dimensionality reduction in multiobjective shortest path search*](https://doi.org/10.1016/j.cor.2015.05.007), Computers & Operations Research 64, 60-70, 2015.

<a id="ref-186"></a>[186] P. Sanders and L. Mandow, [*Parallel label-setting multi-objective shortest path search*](https://doi.org/10.1109/IPDPS.2013.89), Proceedings of the IEEE International Symposium on Parallel and Distributed Processing (IPDPS) 27, 215-224, 2013.

<a id="ref-187"></a>[187] X. Xu, N. Zhang, L. Qian, X. Liu, J. Wu and H. Tang, [*Asynchronous parallel search for exact multi-objective shortest paths with versioned frontier snapshots and indexed dominance pruning*](https://arxiv.org/abs/2609.11944), arXiv:2609.11944, 2026.

<a id="ref-188"></a>[188] P. Maristany de las Casas, L. Kraus, A. Sedeño-Noda and R. Borndörfer, [*Targeted multiobjective Dijkstra algorithm*](https://doi.org/10.1002/net.22174), Networks 82(3), 277-298, 2023.

<a id="ref-189"></a>[189] P. M. Fenwick, [*A new data structure for cumulative frequency tables*](https://doi.org/10.1002/spe.4380240306), Software: Practice and Experience 24(3), 327-336, 1994.
