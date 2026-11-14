# Best-first search with an incumbent: design note

The extension of the searches of the library to Branch-and-Bound (B&B), i.e.,
to the minimization of a function over a finite (but possibly huge) set that
is recursively partitioned into nodes with a lower bound and possibly a
feasible solution, is treated in chapter 13 of the
[manual](https://dmeoli.github.io/PBeFS/manual/index.html): the concept
`BranchAndBound` of [BranchAndBound.h](../include/PBeFS/bnb/BranchAndBound.h)
and its requirements in [section 13.1](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-1), the nodes, the bounds and
the incumbent in [section 13.2](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-2), the sequential search and its
correctness in [sections 13.3 and 13.4](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-3), the selection rules in
[section 13.5](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-5), the parallel searches and what each one
guarantees in [sections 13.7 to 13.10](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-7), the reoptimization in
[section 13.11](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-11) and the interface that a B&B solver built on a
modeling system would use in [section 13.12](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-12). The test of the
searches is the case `bb` and the cases `bb-*` of
[test/tests.txt](../test/tests.txt).

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

The evaluator may give the children one at a time instead, and the child
k of a node without evaluating it [see the section on the lazy children
below]. Optionally, P may provide `hash( n )` for the distribution of the nodes
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
  abstraction-based distributions are not. Since a worker extracts the
  best node of its own queue only, it may evaluate nodes far worse than
  the best open one before an incumbent prunes them, and how many depends
  on when the incumbent is found, i.e., on the timing of the workers;
  plunging, whose dives reach a leaf early, keeps them close to those of
  the sequential search.
- The farm with feedback (FastFlow) gives the master-worker B&B: the
  master holds the queue and the incumbent, and the workers evaluate and
  branch the nodes it sends, returning the children. This scheme is
  probably the best suited to an expensive evaluation, since the order of
  extraction is the one of the policy and the master is idle most of the
  time. The `batch` nodes of a task are extracted together, and the
  children of the first of them cannot overtake the others; hence, the
  batch trades the messages for nodes evaluated beyond those of the
  sequential search, and a batch of 1 is the one for an expensive
  evaluation.
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
  incumbent found by a worker does not prune the nodes of the others. We
  end the sampling with a dive from a copy of its first open node to a leaf
  (`split_dive`, on by default), since a best-first sampling seldom reaches
  one and a worker that starts without an incumbent, on subtrees that do
  not hold the optimum, evaluates every node of them whose bound beats
  their best value; and `split_shared` lets the workers prune by the best
  incumbent of all of them, giving up the determinism of the nodes
  evaluated (not of the value). With one worker there is nothing to split:
  the root goes to it with neither the sampling nor the dive, and the
  nodes evaluated are those of the sequential search.
- Work stealing [[20]](../README.md#ref-20) gives a decentralized B&B whose
  workers move nodes only when one of them is idle. Each worker pushes the
  children it generates in its own queue and extracts its best node; a
  worker whose queue is empty steals about half of the queue of another
  one, drawn at random, and the search ends when no worker holds a node nor
  works on one. Hence, a
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

## Lazy children and the modes

The children given all at once by `branch( n , sink )` suit the searches
that keep a queue, since each child is pushed anyway; a depth-first search,
instead, visits the children of a node one at a time and often abandons the
rest of them, e.g., when the children come by decreasing bound and the first
one pruned makes the others prunable too. Hence, the evaluator may give
`children( n )` in place of (or besides) `branch( n , sink )`, a cursor
whose `next()` gives the next child with its estimate, or nothing, so that
the children that the search never asks for are never built; the searches
that take all the children at once read the cursor to its end, and those
that take them one at a time wrap `branch( n , sink )` in a cursor over a
vector, so that every problem runs on every search. The option `prune_level`
states that the children come by their estimates, the best first, and then
a child pruned when it is generated ends its siblings, which with a cursor
are never built: on the clique of a DIMACS graph, whose children come by the
colours of a greedy colouring, this builds about a tenth of the children
(693899 against 7445660 on brock200_1, with the same 524730 nodes
evaluated). The evaluator may also give `child( n , k )`, the child k of n
built without evaluating n, which the searches that rebuild a node from the
path of its choices, or from its interval of leaves, use instead of
evaluating every node of the path.

The mode of the options says what the search looks for. In the mode
`optimize` it seeks an optimal value, as above. In the mode `decide` it
seeks a feasible value within the `target` of the options (any feasible
value if none is given), pruning the nodes whose bound is worse than the
target, and it stops as soon as one is found; the status is then optimal if
one was found and infeasible otherwise, which is the decision problem of,
e.g., the subgraph isomorphism. In the mode `enumerate` it visits the whole
tree, pruning only the infeasible nodes, and the counter `solutions` of the
statistics holds the feasible values found, i.e., the feasible nodes of the
tree, which is how the nodes of a tree of UTS or the embeddings of a
subgraph are counted. Both modes are implemented in the test that prunes a
node (`bb_prunable`), which in the mode `decide` returns true for every node
once the target is met, and therefore every scheme stops there without a
signal of its own: a queue whose best node is prunable is dropped, as it is
when no node can improve the incumbent. The only exceptions are the
SelfSplit, whose workers share their incumbent in the mode `decide` even
without `split_shared`, so that one of them that meets the target stops the
others, and its dive, which is skipped in the mode `enumerate`, since its
nodes would be counted twice. A richer enumeration than the count (e.g., the
number of the nodes of each depth, or any other monoid) is left to the
evaluator: the evaluators that `make_worker( i )` gives may accumulate into
storage of the problem indexed by i, which the caller combines at the end.

## The skeletons

The parallel searches above share a queue, or keep one per worker, and
order it by the selection. The search skeletons of Archibald et al.
[[67]](../README.md#ref-67) [[68]](../README.md#ref-68) are depth-first
instead: a worker searches a task (a node and its subtree) with a stack of
frames, each holding an evaluated node and the cursor of its children, and
the work is shared by one of four rules, which `SkeletonBranchAndBound`
implements on the executors of ADS, each worker keeping its tasks in a pool
by depth, the deepest taken by its owner and the shallowest by a thief.

- With `depth_bounded`, the children of the nodes of depth less than
  `spawn_depth` become tasks, and those below are searched in place; the
  tasks are as many as the nodes of that depth, whatever the shape of the
  tree.
- With `budget`, a worker that has backtracked `budget` times turns into
  tasks the children left in all the frames of its stack but the top, and
  starts counting again, so that the tasks adapt to the tree: a worker in a
  small subtree finishes before its budget, and one in a large subtree
  shares it.
- With `stack_stealing`, an idle worker asks a busy one, which answers at
  its next step (i.e., after at most an evaluation) with the next child of
  the lowest frame of its stack that has one, the top excepted, or with all
  the children of that frame (`steal_all`), which are the largest pieces of
  work it holds.
- With `ordered`, the nodes of depth `spawn_depth` are the tasks, made at
  the start by the worker 0 (which evaluates and branches the nodes above
  them) and ordered by their discrepancy, i.e., the sum of the indices of
  the children on their path (or left to right, `discrepancy` off); all the
  workers take them in this order, so that the tasks done by the time a task
  is started are always those before it, which is what makes the runtimes
  of the Ordered skeleton repeatable and never worse than those of the
  sequential search in the order of the tasks [[67]](../README.md#ref-67).

The tasks are counted when they are made and when their subtree is done, so
that a worker that finds no task anywhere ends when the count is 0, which
is exact. A steal from a stack follows a protocol of three states on the
victim (free, asked by the worker j, claimed): the thief asks, the victim
claims the request before it builds the answer and releases it after, and a
thief that sees the count fall to 0 or a limit withdraws its request only if
it is not claimed yet, so that an answer is never lost; a worker that waits
for an answer answers the requests it gets with nothing, and a worker that
ends does so too. With `indexed`, a task is not a node but the path of the
choices that lead to it from the node it comes from (the root, or a seed),
and the worker that takes it rebuilds the node, by `child( n , k )` if the
evaluator has it and otherwise by evaluating and branching each node of the
path, the evaluations counted apart (`recomputed`); a task then takes four
bytes per level, whatever the node, which is what the messages of a search
over many processes need. After a limit, the children left on the stacks
and the tasks left in the pools are open leaves, the latter rebuilt if they
are paths, and the search can be resumed from them as any other. The order
is depth-first in the order of the children, and the selection of the
options is not used. In the taxonomy of Gendron and Crainic
[[16]](../README.md#ref-16) these are again searches of one tree (their type
2), asynchronous and with a pool per worker.

## The racing ramp-up, the checkpoints and the hubs

At the start of a search there are fewer open nodes than workers, and the
schemes above can only wait for the tree to grow. The racing ramp-up of the
Ubiquity Generator [[69]](../README.md#ref-69) uses this time to search the
tree from the root with several settings at once, keeping the tree of the
best one, and `racing_branch_and_bound()` does so before any of the schemes
(`parallel_branch_and_bound()` calls it when `race` lists the settings):
the racer k runs the sequential search with the evaluator of the worker k,
the selection, the dive and `prune_level` of its setting and at most
`race_nodes` evaluations, all the racers pruning by the best incumbent of
all of them, and the scheme goes on from the open leaves of the racer with
the best bound (the fewest open leaves on ties), with the best incumbent of
all, unless a racer ended the search. The search is correct since the open
and the closed leaves of the winner cover the root, and each closed leaf was
pruned by an incumbent not better than the final one; the bound of the
closed leaves of the winner enters the bound of the search through
`closed_bound`. The settings that differ in the evaluator (e.g., in the
branching rule, or in the cuts of a relaxation) are obtained by giving the
racer k an evaluator by `make_worker( k )`, which the problem builds as it
likes.

A search that may run longer than the time granted to it (say, by a batch
system) needs a checkpoint, i.e., a file from which it can be resumed.
Since a search already starts from any set of nodes, what is missing is the
file: `checkpointed_branch_and_bound()` runs any search in segments of a
given number of evaluations, each from the open nodes of the previous one
with the incumbent found so far and with `keep_open_only`, which keeps only
the leaves left open by the limit and the least bound of the others in
`closed_bound`, and writes after each segment the open nodes with their
bounds and depths, the incumbent and its node, the bound of the closed
region and the counters, as the bytes that the search over MPI sends
(`serialize( n , out )` and `deserialize( in )`, or the bytes of a
trivially copyable node). The file is written to a temporary one that then
takes its name, so that a run stopped while writing leaves the previous
checkpoint whole, and a run given the name of an existing checkpoint starts
from it, with any scheme and any number of workers. The open nodes of the
last checkpoint and the closed regions of all the segments cover the tree,
by induction on the segments, and the bound of the search is the least of
their bounds and of the incumbent; the price of a segment is the
synchronization of a search that ends and one that starts.

Over many processes, a search that sends every child to its owner moves
much more than it needs. In the master-hub-worker organization of ALPS
[[70]](../README.md#ref-70), each process searches its own subtrees, i.e.,
the open nodes it holds, a unit of work at a time, and reports its workload
to the hub of its cluster: its quantity (the number of its open nodes) and
its quality (the mean bound of its best nodes, the elite). A hub pairs the
workers of its cluster whose quantity is at most `need_work` with those that
have more, or, when none is short of work, those whose quality is worse than
the mean by more than a fraction with those better than it by more than
another one, and tells each donor to give half of its nodes (the second,
the fourth and so on in the order of their bounds, keeping the best) to its
receiver; a receiver with no donor in its cluster is sent to the master,
which passes it to the hub with the largest quantity. `MPIHubBranchAndBound`
implements this on the mailbox of ADS, whose four-counter method of Mattern
detects the end as in `MPIBranchAndBoundSearch`: every message is counted, a
process with no node sends only the report that it is empty, and the hubs
and the master answer only reports, so that the messages stop when every
queue is empty. A receiver paired with a donor is not paired again until
the donor reports, which it does after each order, and the balance of the
quality runs once every as many reports as the workers of the cluster, so
that a donor is not asked again before the others have spoken.

## Intervals of leaves

In a tree whose nodes of depth d all have the same number $r_d$ of children,
as those of a permutation problem ($r_d = n - d$), a leaf is the vector of
the indices of the children on its path, an integer of
$[ 0 , \prod_d r_d )$ in a mixed radix (the factorial number system for the
permutations), and the leaves of a node are an interval. The
Integer-Vector-Matrix of Mezmaz et al. [[71]](../README.md#ref-71) and
Gmys et al. [[72]](../README.md#ref-72) makes such an interval the unit of
work: a worker explores depth-first the leaves of its interval
$[ a , b )$, holding the path to $a$, and an idle worker takes the second
half of the interval of another one, which moves two vectors of digits,
whatever the node. `IntervalBranchAndBound` implements it for any problem
that gives `arity( d )` and an evaluator with `child( n , k )`: the steal
takes the lock of the victim, halves $b - a$ in the mixed radix and moves
$b$ of the victim to the middle, and the thief builds the path to its new
$a$ by `child( n , k )`. A node is evaluated as such only by the worker
whose interval holds its first leaf; the ancestors of the leaf at which a
worker took its interval are evaluated again, for their bound, and counted
as recomputed. Hence, each node is evaluated once, and with the incumbent
fixed (e.g., given as the optimum, as the experiments on the flowshop
usually do) the nodes evaluated are those of the sequential depth-first
search, whatever the steals, which the test checks with 1 to 8 workers.
After a limit, the leaves left in each interval are covered by the fewest
nodes, each the shallowest one whose leaves begin at $a$ and end before $b$,
which become open leaves with the estimate along their path as their bound.
The dynamic branching of the flowshop, whose direction depends on the
incumbent, gives a tree that changes during the search, and it is therefore
refused.

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
`deterministic`, `selfsplit`, `stealing`, `multiqueue`, `depthbounded`,
`budget`, `stackstealing` or `ordered`) read, e.g., from its own parameters,
so that it needs no change when the library gets a new scheme (the settings
of a race in `race` putting the racing ramp-up before any of them), or else
one of the classes, e.g.,

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
with 1 and 4 open nodes per worker and with and without the dive, they have
to be the same at every repetition, each worker has to evaluate its nodes by
its own evaluator and as many as the counters say, and the values given to
the callback of the incumbent have to improve each time; with the incumbent
shared, the value has to be the optimum with every selection and with 1, 2,
3 and 8 workers.

The extensions of the last sections are checked on more problems, each
against a reference that does not come from the library. The four
skeletons, with nodes and with paths as tasks and with 1 to 8 workers on
the three executors, have to give the optimum of the knapsack and its node,
a bound not worse than the optimum when a node limit stops them, and the
optimum again from the open leaves of the stopped search with its
incumbent; on a tree of UTS they have to count its nodes as a plain
recursion does, each node being evaluated once; and on the cliques of a
random graph they have to give the same clique with and without
`prune_level`, the latter building fewer children. The modes are checked
with every scheme on small knapsacks: `decide` has to find a value within a
target below the optimum, none beyond it, and some value without a target,
and `enumerate` has to count the feasible nodes of the tree as a recursion
on `evaluate` and `branch` does. The race, of three selections and budgets
of 5, 50 and 100000 nodes, has to give the optimum, its node and the bound
before every scheme. The checkpoints have to give the optimum in segments
of 7 evaluations, to leave an empty checkpoint at the end, and, for a run
stopped after 20 evaluations, a checkpoint from which a second run ends the
search, with the counters going on. The intervals have to give the optimal
makespan of small flowshops and the optimal cost of small quadratic
assignments (symmetric and not), against the enumeration of all the
permutations, with every bound and the static branchings, the same nodes
evaluated with 1 to 8 workers when the optimum is given, and the optimum
again from the leaves left by a node limit. Finally, the domains have to
match their references: the cliques of random graphs those of the algorithm
of Bron and Kerbosch, and that of brock200_1 its known size 21; the common
subgraphs the largest clique of the modular product of the two graphs,
which is a different algorithm on a different graph; the subgraph
isomorphisms, induced or not, the count and the existence of the
embeddings found by enumerating all the injections; the instances chr12a,
had12 and nug12 of QAPLIB their published optima; the small flowshops the
enumeration with both bounds and the four branchings, and ta001 its optimum
1278; and three trees of UTS the plain recursion.
