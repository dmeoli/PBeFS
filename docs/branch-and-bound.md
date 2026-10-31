# Best-first search with an incumbent: design note

We are concerned with extending the searches of the library to
Branch-and-Bound (B&B), i.e., to the minimization of a function over a
finite (but possibly huge) set that is recursively partitioned, each part
(a *node*) carrying a lower bound on the values it contains and, possibly,
a feasible solution. In fact, A* is a special case of this
scheme, where the bound of a node is $f = g + h$ and a goal gives a
feasible solution of value $g$.
However, the general case differs in five respects that the `Domain`
concept does not cover. First, the tree of a B&B has no duplicates, and
therefore no closed list is needed. Second, the search does not stop at the
first feasible solution, but when no open node can improve the best one
found so far (the *incumbent*), or when the gap between the two is within a
tolerance. Third, the bound of a node is usually "expensive" (e.g., the
optimal value of a Linear Program), and hence one wants to compute it when
the node is extracted, the node inheriting the bound of its parent until
then. Fourth, the order of extraction is a policy (best-first, depth-first, breadth-first,
best-estimate, or a depth-first dive from the best node). Fifth, the
objective may be maximized, and the bounds are real numbers compared with
tolerances. The structure of this note is as follows: we first give the
concept of a B&B problem, then the sequential search, then how the parallel
algorithms of the library map onto it (and which ones do not), and finally
the interface that a B&B solver built on a modeling system would use.

## The concept

A type P is a `BranchAndBound` problem if it provides the copyable type
`node_type`, a constant `sense` (`Sense::minimize` or `Sense::maximize`),
and three methods. The method
`root()` returns the node of the whole set. The method `evaluate( n )`
computes the bound of the node n, possibly storing in n what the branching
needs (say, the fractional solution of a relaxation), and returns an
`Evaluation` made of the bound, a flag for an infeasible node and an
optional feasible value, of which n is the solution. Also,
`branch( n , sink )` calls `sink( child , estimate )` for the children of
the evaluated node n, where `estimate` is a valid bound of the child that
is "cheap" to compute, typically the bound of n; the search uses it to
order the child and to prune it before its evaluation.

Optionally, P may provide `hash( n )` for the distribution of the nodes
among the workers, `depth( n )` for the depth-first policies (otherwise the
search counts the levels), and `make_worker( i )`. The latter returns, for
the worker i of a parallel search, an object with the same `evaluate` and
`branch`: this is needed when the evaluation keeps a state (e.g., a model
and a solver that are modified at each node), since then each worker has
its copy and P is only read. Such an object may also give `cutoff( z )`,
which the searches call before each `evaluate` and `branch` with the
incumbent value z (in the sense of P, infinite when there is none), so
that, e.g., the relaxation of a node can stop as soon as it is proven not
to improve z, the node being then infeasible to the search. The worker
should take z from here rather than from a value shared with the other
workers, since the deterministic search gives the incumbent of the start
of the round, and an evaluation against a value that the others change
meanwhile would make its outcome depend on their timing. A maximization
problem is solved by negating the values internally; thus, the
comparisons below are those of a minimization.

The values (the bounds, the feasible values, the incumbent and the absolute
tolerance) are of the type `value_type` that P may give, `double` if it gives
none, and the options, the evaluation and the solution are then
`BasicBBOptions< V >`, `BasicEvaluation< V >` and `BBSolution< Node , V >`, of
which `BBOptions`, `Evaluation` and `BBSolution< Node >` are the case of
`double`. Any signed type with a total order serves, as long as
`BBValueTraits` says what its infinity is, how large the tolerance of a gap
is and whether a child is close enough to the best bound to dive into it:
the primary template does so for the arithmetic types (the infinity of an
integer type being its largest value), and a caller with another type, e.g.,
the rationals of an exact arithmetic, specializes it. An integer type pays
off whenever the values of the solutions are integer, as in the 0-1 knapsack
with integer profits, since the bound of a node can then be rounded (down for
a maximization, up for a minimization) and a node whose rounded bound equals
the incumbent is pruned, while with doubles it is branched unless an
absolute tolerance just below 1 is given. We remark that with a type whose
`std::atomic` is not lock-free the parallel searches read the shared
incumbent under a lock, and that the search over MPI needs trivially
copyable values, which travel as their bytes.

Note that the `Domain` concept and the searches that use it are left
untouched: the B&B layer is a new concept with its searches (`BranchAndBound.h` and
`ParallelBranchAndBound.h`), which share with A* the runtimes, the
placement of the workers (cf. `pin_worker()`) and the detection of the
termination. Yet, they do not share the node table or the OPEN lists of A*,
which are specialized for integer costs and for the lazy deletion of
duplicates. We remark that A* is not rewritten on top of the new layer,
since this would add an indirection to its inner loop for features that it
does not use.

## The sequential search

We keep the open nodes in a pool with a free list (a node is
destroyed when it is branched, while in A* a generated node stays in the
closed list) and a priority queue of pairs (key, index). Its key depends on the policy: the bound for best-first, minus the depth for depth-first, the depth for breadth-first (the oldest node first among those of the same depth), and
a problem-given estimate of the best feasible value in the subtree for
best-estimate. In the dive (*plunging*), the children of the node just
branched are extracted depth-first as long as their bound is within a
given fraction of the gap from the best bound, after which the search goes
back to best-first. At an extraction, a node whose inherited bound cannot
improve the incumbent by more than the tolerance is dropped without being
evaluated. Otherwise it is evaluated, and it is dropped if it is infeasible
or if its bound cannot improve the incumbent; its feasible value, if any,
updates the incumbent and calls the `on_incumbent` callback of the options
(which gives the anytime behavior). Then the node is branched, and a child
whose estimate cannot improve the incumbent is dropped at once.

The search ends when the queue is empty, or when the best bound among the
open nodes is within the tolerance of the incumbent, which we take as
$\max \{ \tau_a, \tau_r |z| \}$ for the incumbent value $z$, with
$\tau_a$ and $\tau_r$ the tolerances `abs_gap` and `rel_gap`. With a policy
other than best-first the best open bound is not at the top of the queue, and
therefore a second heap of the bounds (with lazy deletion) keeps it; in
this way the termination test and the final gap are exact for all the
policies. The `Solution` reports the incumbent and its node, the best
bound, the gap, the nodes evaluated, branched and pruned (at generation and
at extraction), and the status, i.e., optimal within the tolerance,
infeasible, or stopped by a node or time limit (where the gap is still
valid). An initial incumbent value, e.g., from a heuristic, can be given
in the options.

## The parallel searches

We now turn to the parallel searches, which differ in how they share the
open nodes. Since a B&B tree has no duplicates, the algorithms whose purpose is
the balance of the load carry over, while those whose purpose is the
detection of duplicates or the avoidance of re-expansions do not.

- The shared queue of SPA* gives the classical centralized parallel B&B:
  the workers extract from one queue under a lock and evaluate the nodes
  outside of it. The contention is rather low when the evaluation dominates (as
  is usually the case in B&B), and the order of extraction is close to the
  sequential one.
- The distribution of HDA* gives a decentralized B&B with work sharing.
  Each worker has a queue, the children are sent in batches to the owner
  given by the hash of the node (or chosen at random when P gives no hash),
  and the termination is that of HDA*, i.e., a counter of the nodes in
  flight on the threads and the algorithm of Mattern over MPI. The
  incumbent is an atomic variable (a reduction over MPI). The distribution
  has no duplicates to detect here, and its only role is to balance the
  load; hence, Zobrist and random hashing are of use, while the
  abstraction-based distributions are not.
- The farm with feedback (FastFlow) gives the master-worker B&B: the
  master holds the queue and the incumbent, and the workers evaluate and
  branch the nodes it sends, returning the children. This scheme is
  probably the best suited to an expensive evaluation, since the order of
  extraction is the one of the policy and the master is idle most of the
  time.
- The synchronous rounds give a deterministic B&B, as the deterministic mode
  of FiberSCIP [[18]](../README.md#ref-18). At each round the `round` best
  open nodes are extracted, the worker $i$ evaluates the nodes
  $i, i + p, \ldots$ of
  the round against the incumbent of its start, and the outcomes are merged in
  the order of the nodes; hence, the nodes evaluated, the incumbent and the
  bound depend only on `round`, and neither on the timing of the workers nor
  on their number, provided that the evaluation of a node by a worker depends
  only on the node and on the nodes that worker evaluated before. The price is
  the wait of each round for its slowest node, and an incumbent that may be
  one round old; the node limit stays exact, while the time limit, checked
  between the rounds, does not.
- The SelfSplit of Fischetti, Monaci and Salvagnin
  [[19]](../README.md#ref-19) gives a B&B whose workers communicate only at
  the end. The tree is first explored sequentially (the *sampling*) until
  the open nodes are `split` per worker; then the k-th of them, in the
  order of the selection, goes to the worker $k \bmod p$, which searches its
  subtrees sequentially with its own evaluator, its own incumbent and its
  share of the node limit, and the incumbents, the bounds and the leaves of
  the workers are merged at the end. Hence, the outcome depends only on
  `split` and on p, under the same condition on the evaluation as above.
  We do the sampling once, by the worker 0, rather than in each worker as
  in the original scheme, since the workers here share the memory and can
  read its open nodes; over MPI, where a node would otherwise travel as
  bytes, each process would repeat it. The price is that the load is
  balanced only as far as the subtrees of the split are, and that an
  incumbent found by a worker does not prune the nodes of the others.
- Work stealing [[20]](../README.md#ref-20) gives a decentralized B&B whose
  workers move nodes only when one of them is idle. Each worker pushes the
  children it generates in its own queue and extracts its best node; a
  worker whose queue is empty takes the best node of the queue of another
  one, the victims being tried from a random one on, and the search ends
  when a counter of the nodes open or held by a worker is zero. Hence, a
  worker stays in its own part of the tree, which pays when the evaluator
  of a worker follows the path from the root to its node (e.g., by applying
  and undoing the changes of the branchings), since a node sent elsewhere,
  as in the distributed queues, forces the receiver to move far. The price
  is a selection that is best-first within each queue only, and hence an
  order of extraction possibly farther from the sequential one than that of
  the shared queue.
- SafePBNF and PA*SE do not apply. The nblocks of PBNF make the detection
  of duplicates local, and the independence test of PA*SE avoids the
  re-expansions; a tree has neither, and PA*SE without duplicates would be
  a shared queue with a costly test at each extraction.

Over MPI, the distributed search has to send the nodes between processes,
and it therefore requires P to provide `serialize( n , out )` and
`deserialize( in )` (or a trivially copyable node). Each process then has a
queue, the incumbent value is sent to all the processes when it improves,
the node of the incumbent stays on the process that found it until the end
(when it is sent to all), and the end is detected with the method of
Mattern, as in HDA* over MPI; the node limit is split evenly among the
processes.

## Reoptimization

When the data of the problem change in a way that leaves the branching
valid (e.g., a restriction of the feasible region, or a change of the
objective that the nodes do not store), the tree explored by a previous
search need not be built again: its leaves cover the subtrees of the nodes
that the search started from, and a new search can start from them. With
`keep_frontier` in the options, the `Solution` holds these leaves, each
with its bound, its depth and the reason why it is a leaf, i.e., pruned by
its bound, evaluated with a feasible value and not branched, found
infeasible, or left open by a limit; the overload `solve( seeds )` (and
`branch_and_bound( p , seeds , options )`) then starts the search from any
set of such nodes instead of the root, each of them being an open node that
is evaluated only when it is extracted. Which leaves are given back, and
with which bounds, is left to the caller, since only the caller knows what
has changed: after a restriction of the region the bounds of the pruned
leaves remain valid and the infeasible leaves need not be explored again,
after an enlargement the infeasible leaves have to be given back with an
unknown bound, and after a change of the objective all of them have to be
given back with an unknown bound, unless the old ones can be corrected
(say, by a bound on the change of the objective over the node). In the same
way, a search stopped by a limit is resumed by giving back its open leaves
together with its incumbent. The bound of the `Solution` is that of the
subtrees of the seeds, and the caller combines it with the bounds of the
leaves that it leaves out.

In the parallel searches each worker collects its leaves, which are merged
at the end, and the seeds go to their owners as the children do; over MPI,
instead, each process keeps the leaves that it owns and starts from its own
seeds, which avoids sending them across the processes. We remark that in a
tree where every inner node has two children the leaves are one more than
the inner nodes, and therefore keeping them may take much more memory than
the open nodes alone; hence, `keep_frontier` is off by default.

## The interface for a B&B solver

A B&B solver built on a modeling system would provide a problem whose
`node_type` is the list of the branching decisions from the root (e.g.,
the bounds changed on the variables, which is "small" and hence cheap to
copy and to send). Its `make_worker( i )` would return an object holding a
copy of the model and of the solver of the relaxation. The `evaluate` of a
worker applies the decisions of the node to its model (possibly undoing
only those that differ from the previous node, which makes a dive cheap),
solves the relaxation, and returns its value as the bound, its
infeasibility, and its value as a feasible one if the solution is
integral. Its `branch` chooses a fractional variable in the stored
solution and gives the two children, with the bound of the parent as their
estimate; also, if its `cutoff( z )` keeps the incumbent value and its
`evaluate` returns as feasible the value of the solution found by the
node itself (keeping that solution in the node), the outcome is that of
the scheme it runs, whichever it is. Such a solver then calls

    auto s = parallel_branch_and_bound< P , Executor >( problem ,
                                                        bb_scheme( name ) ,
                                                        options );

with the name of the scheme (`sequential`, `shared`, `distributed`,
`deterministic`, `selfsplit` or `stealing`) read, e.g., from its own
parameters, so that it needs no change when the library gets a new
scheme, or else one of the classes, e.g.,

    auto s = SharedBranchAndBound< P , Executor >( problem , options ).solve();
    auto s = FarmBranchAndBound< P , Executor >( problem , options ).solve();
    auto s = MPIBranchAndBoundSearch< P >( problem , options , comm ).solve();

and reads the incumbent, its node, the bound and the gap from the
`Solution`. In the options one gives the policy, the tolerances, the
initial incumbent and the callback, and for the parallel searches the
number of workers, the batch size, the nodes of a round and those of the
split; the executor of ADS (its pool of C++ threads, OpenMP or FastFlow)
is the second parameter, and `set_executor()` gives a search an executor
that outlives the search. We expect the choice among the
schemes to depend somewhat on the cost of an evaluation w.r.t. that of
moving a node, which only a measurement on the instances of interest can
settle.

## The test

As a test of the extension we use the 0-1 knapsack problem with the bound
of the fractional relaxation of Dantzig. Its items are sorted by decreasing
ratio of profit to weight, and a node fixes the first k of them; its bound
is that of the greedy fractional completion, the greedy integer completion
gives a feasible value, and the branching fixes the critical item to 1 and
to 0. Hence, one can check the optimum against the dynamic program on the capacity,
on random instances of the uncorrelated, weakly correlated and strongly
correlated classes of Pisinger with up to a few hundred items. We run the policies sequentially and the parallel schemes with 1 to 8
workers, as a maximization (which exercises the change of sense), and with
both a zero and a positive relative tolerance; in the latter case the value
found has to be within the tolerance of the optimum. The leaves are checked
as well, in both senses: those of a search carried out to the end have to
bound the optimum, and a search started from them has to solve the knapsack
of three quarters of the capacity (and, when the children that do not fit
are generated and found infeasible, the original knapsack again from the
leaves of the smaller one), while a search started from the open leaves of
one stopped by a node limit, with its incumbent, has to end it. For the
deterministic search, with four nodes per round, the value, the node, the
bound and the counters have to be the same with 1, 2, 3 and 8 workers and
at every repetition, also when a node limit stops it. For the SelfSplit,
with 1 and 4 open nodes per worker, they have to be the same at every
repetition, each worker has to evaluate its nodes by its own evaluator and
as many as the counters say, and the values given to the callback of the
incumbent have to improve each time.
