# Parallel branch-and-bound

The parallel B&B searches of the library (the shared queue, the distributed
queues, the work stealing, the MultiQueue, the deterministic rounds, the
SelfSplit, the farm and the processes of MPI), with what each one guarantees,
how it ends, whether its outcome is deterministic, and the anomalies in the
number of evaluations of a parallel B&B, are treated in
[chapter 13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of the [manual](https://dmeoli.github.io/PBeFS/manual/index.html);
the pseudocode of each scheme is also on the page of its header.

We now turn to the parallel B&B searches of the library, all of which build one
tree with $p$ workers; they differ in how they share the open nodes and the
incumbent. In the classification of Gendron and Crainic
[[16]](../../README.md#ref-16) they are of the second type, i.e., the
parallelism is in the exploration of the tree. Hence, the parallelism inside
the evaluation of a node (the first type) is left to the problem, and the
search of several trees at once (the third type) is provided only at the
start, by the racing ramp-up. Within the second type the searches are
asynchronous (but for the rounds of the deterministic search), and they keep
one pool of open nodes (the shared queue, the farm and the rounds) or one pool
per worker (the distributed queues, the work stealing, the SelfSplit, the
skeletons, the intervals of leaves and the processes of MPI). This page states what a
parallel B&B has to guarantee and the anomalies in the number of its
evaluations; then, for each scheme, it discusses the correctness, the
termination and whether the outcome is deterministic. The notation is that of
[Branch-and-bound as best-first search](branch-and-bound-search.md), the
interface and the design choices are discussed in the
[design note](../branch-and-bound.md), and the pseudocode of each scheme is on
the page of its header.

## What a parallel search has to guarantee

The invariant of the sequential search holds for a parallel one if the open
nodes are taken in a wider sense, i.e., as the nodes in any queue, those held
by a worker (being evaluated, or the next node of a dive) and those traveling
in a message. Each element of $X$ below the threshold
$z - \max \{ \tau_a, \tau_r |z| \}$ then belongs to the subset of one of them.
Each worker preserves the invariant by the same steps as the sequential search,
provided that it prunes with an incumbent that is never smaller than a value
actually found; this holds, since $z$ only decreases and each value stored in
it is feasible, and a worker that prunes with an older (larger) value of $z$
only prunes less. Hence, a parallel search is correct if it stops only when no
node is open in this wider sense, or when the least bound among them is within
the tolerance of $z$. Its reported bound is the least bound over the leaves of
all the workers, which each worker collects during the search and which are
merged at the end.

In the library, the shared incumbent (`BBIncumbent` in
[ParallelBranchAndBound.h](../../include/PBeFS/bnb/ParallelBranchAndBound.h))
keeps its value in an atomic variable, which the workers read without a lock
when they prune. An improvement takes a lock, checks the value again, and
stores the node and calls the callback under the lock.

Of course, the evaluations of a parallel search are not those of the sequential
one, and they may be more or fewer (often, but not always, more). On the one
hand, a worker may evaluate a node that the sequential search prunes with an
incumbent which, in the parallel search, another worker finds only later (a
*detrimental* anomaly, which is the search overhead of a B&B). On the other
hand, another worker may find an incumbent earlier than the sequential search
and prune nodes that the latter evaluates (an *acceleration* anomaly, which
gives a speedup larger than $p$). Both are discussed in
[[16]](../../README.md#ref-16) and in the references therein. With the
best-first selection the detrimental part is somewhat limited, since the nodes
extracted by the workers are among the best open ones; it grows with the
distance between the order of extraction and the sequential one, and this
distance distinguishes the schemes below.

## One pool of open nodes

In the shared queue (`SharedBranchAndBound`) the workers extract the nodes from
one queue under a lock, evaluate and branch them outside of it, and push the
children back under it; with depth-first or plunging, the child that dives
stays with the worker and is not pushed. In this scheme the order of extraction
is the sequential one but for the $p$ nodes evaluated at once, and the lock is
rather cheap when the evaluation dominates, which is the usual case of a B&B.
The search ends when the queue is empty and no worker holds a node, which is
detected by a counter of the busy workers. A worker raises the counter under
the lock when it extracts a node, and lowers it only after having pushed the
children of the node; thus, an empty queue and a null counter, read under the
lock, mean that no node is open.

The farm (`FarmBranchAndBound`) is the master-worker form of the same scheme,
written with FastFlow [[12]](../../README.md#ref-12). Here the master holds the
queue, extracts the nodes by the selection and sends them, in batches, to the
free workers. A worker evaluates and branches each node (going on with the
child of a dive in the same task), and it sends back the children that survive
the pruning. The order of extraction is that of the selection, and the master
only moves nodes; the search ends when the queue holds no node that can improve
the incumbent and no task is out at a worker, which the master alone can see.
Hence, the farm suits an evaluation that is "expensive" w.r.t. the transfer of
a node.

The rounds (`DeterministicBranchAndBound`) make the outcome independent of the
timing, as in the deterministic mode of FiberSCIP
[[18]](../../README.md#ref-18). At each round the $r$ best open nodes are
extracted from one queue ($r = p$ by default), and the worker $i$ evaluates and
branches the nodes $i, i + p, i + 2p, \ldots$ of the round against the
incumbent of the start of the round. When all of them are done, one worker
merges the outcomes in the order of the nodes of the round (the feasible values
update the incumbent, the children that cannot improve it are pruned and the
others are pushed) and extracts the next round, while the others wait at a
barrier. The nodes evaluated, the incumbent and the bound then depend only on
$r$, and neither on the timing nor on the number of the workers. This holds
provided that the evaluation of a node by the worker $i$ depends only on the
node and on the nodes evaluated by $i$ before, which a separate evaluator per
worker (`make_worker( i )`) allows. The node limit is exact, the last round
being shortened to it, while the time limit, checked between the rounds, makes
the end depend on the timing. The price is the wait of each round for its
slowest node, which is small when the evaluation dominates and its time varies
little among the nodes, and an incumbent that may be one round old; also, no
dive is taken, and the plunging selection extracts the nodes by their bound.

## One pool per worker

In the distributed queues (`DistributedBranchAndBound`) each worker has a
queue, and it sends each child, in batches, to the worker that owns it, i.e.,
to $\mathrm{hash}(n) \bmod p$ if the problem gives a hash of its nodes, and to
a random worker otherwise, as HDA* does with the states of A*. Since a tree has
no duplicates, the owners only balance the load, and a uniform hash is the
suitable distribution (the abstractions of HDA* have nothing to gain here).
Their termination is detected as in HDA*, by the deficit counting of the
mailbox of ADS (see [Parallel A*](parallel-astar.md)); a worker is passive when
its queue holds no node that can improve the incumbent and it holds no node of
a dive.

Over MPI (`MPIBranchAndBoundSearch`) the same scheme runs with one process per
worker, and the nodes travel as bytes. The incumbent value is sent to all the
processes when it improves (its node stays with the process that found it until
the end), and the node limit is split evenly among the processes. Again, the
end is detected by the four-counter method of Mattern
[[11]](../../README.md#ref-11), as in HDA* over MPI.

In the work stealing (`StealingBranchAndBound`, after Blumofe and Leiserson
[[20]](../../README.md#ref-20)) each worker pushes the children that it
generates in its queue, and extracts its best node; the queues are those of
the frontier of work stealing of ADS, from which the owner extracts its nodes
in the exact order of the selection and with no lock. A worker whose queue is
empty steals about half of the queue of another one, drawn at random, the half
that keeps out the best node of the victim (cf. the steal-half queues of
Hendler and Shavit [[36]](../../README.md#ref-36)). A worker then stays
in its part of the tree, which pays off when moving its evaluator to a distant
node costs (e.g., when the evaluator applies and undoes the branchings along
the path from the root), while the nodes of the distributed queues land on any
worker. The price is a selection that is best-first within each queue only, and
hence an order of extraction that may be farther from the sequential one than
that of the shared queue. The search ends when no worker holds a node nor
works on one, which the frontier tells by a counter of the workers that hold
or work on a node plus the thieves that are trying a steal (a thief enters it
before it claims the half of a victim and leaves it if it fails); hence, the
counter is never null while a node exists.

The SelfSplit of Fischetti, Monaci and Salvagnin [[19]](../../README.md#ref-19)
(`SelfSplitBranchAndBound`) gives a B&B whose workers communicate only at the
end. In this scheme the tree is first explored sequentially (the *sampling*),
by the evaluator of the worker 0 and without plunging, until the open nodes
are $\sigma p$ ($\sigma$ being `split`, 4 by default) or the tree is solved,
and then (with `split_dive`, the default) the same evaluator dives from a copy
of the first open node to a leaf, extracting at each level the child of least
estimate (or of least bound, if the problem gives no estimate). Then the
$k$-th open node, in the order of the selection, goes to the worker
$k \bmod p$, which searches its subtrees sequentially with its own evaluator,
its own incumbent (starting from that of the sampling) and its share of the
node limit; the incumbents, the bounds and the leaves of the workers are merged
in the order of the workers. The search is correct, since the subsets of the
nodes of the split cover what the sampling left open, the nodes of the dive
included. Its outcome depends only on $\sigma$ and on $p$, under the same
condition on the evaluation as the rounds. In the library the sampling is done
once, rather than by each worker as in the original scheme, since the workers
share the memory. The price is a load that is balanced only as far as the
subtrees of the split are, and an incumbent found by one worker that does not
prune the nodes of the others. The dive limits the second: a best-first
sampling seldom reaches a leaf, and without an incumbent a worker whose
subtrees do not hold the optimum evaluates each of their nodes whose bound is
better than the best value in them (on 20 weakly correlated knapsacks of 40
items, with $\sigma = 4$ and 2, 4 and 8 workers, up to 53, 850 and 118 times
the nodes of the sequential search, against at most 1.4, 1.3 and 2.8 times
with the dive), at the cost of the few nodes of the dive, which the worker
whose subtree holds them evaluates again. With `split_shared` the workers also
prune by the best incumbent of all of them, which each reads at every
extraction (at most 1.5 times the nodes of the sequential search on the same
knapsacks); the value is still the optimum, while the nodes evaluated, and the
node of the optimum among equal values, depend on the timing.

## Depth-first skeletons, intervals of leaves and hubs

The skeletons [see `SkeletonBranchAndBound.h`] keep the open nodes in a
wider sense still: besides the tasks in the pools, they are the children
that the cursors of the frames of the stacks have not given yet, and the
tasks of an answer to a steal that the thief has not taken yet. A child
leaves a cursor only by the worker that owns the stack (when it visits the
child, or when it turns it into a task, for the budget, the depth bound or
a thief), and a task leaves a pool under its lock; hence, no open node is
lost or taken twice, and the invariant of the sequential search holds by
the same steps. Let $\pi$ be the number of the tasks made and not finished,
a task being finished when the stack of its depth-first search is empty;
$\pi$ is raised before a task is published (in a pool, or in an answer) and
lowered after its subtree is done. A worker that finds no task in its pool
nor in the others, and $\pi = 0$, ends: then no task is in a pool, in an
answer or under search, and since every open node belongs to a task, none
is left. The steal of a stack needs that an answer is never lost: the
request of the worker $j$ at the victim $v$ goes from free to "asked by
$j$", the victim moves it to "claimed" before it builds the answer and back
to free after it has published the answer (and raised $\pi$), and the thief
withdraws its request (when $\pi = 0$ or after a limit) only by moving it
from "asked by $j$" to free; hence, either the victim claims it, and the
thief waits for the answer, which comes, or the thief withdraws it, and no
answer is built. Since a thief that waits answers its own requests with
nothing, two workers that steal from each other do not wait for each other.
The rebuild of a task given as a path needs that the child $k$ of a node be
the same at every call, i.e., that the children of a node come in an order
that depends on the node alone; with this, the node rebuilt is the node that
was spawned, and its bound, carried by the task, is valid. With the ordered
tasks, the outcome is not deterministic in the nodes evaluated, since a
task is pruned by the incumbent found by the tasks before it, which may be
still running; however, the tasks done when a task starts are always those
before it in the fixed order, and the worker that takes the next task is
never idle while one is left, so that, as shown by Archibald et al.
[[67]](../../README.md#ref-67), the time is never worse than that of one
worker searching the tasks in the same order and repeated runs take times
close to each other.

The intervals of leaves [see `IntervalBranchAndBound.h`] keep as open the
leaves of $[ a_i , b_i )$ for each worker $i$, the leaves before $a_i$ in
the interval of $i$ having been explored or pruned. A steal takes the lock
of the victim, cuts its interval at the middle and releases the lock,
after which the thief, under its own lock only (so that no lock is held
while another is taken), makes the second half its interval; the intervals
stay disjoint, and their union loses nothing, since in between the second
half belongs to the thief, whose interval is empty and is not taken by
anyone; a worker moves $a_i$ forward only past a leaf
it has explored or past the subtree of a node it has pruned, which holds no
solution better than the incumbent. The search ends when the number of the
workers holding a non-empty interval is 0; a thief raises it under the lock
of a victim that holds one (and is therefore counted), so that it never
drops to 0 while some leaf is open. A node is evaluated as such only by the
worker whose interval holds its first leaf, when $a_i$ is that leaf: the
leaves of an interval are visited in order and $a_i$ never goes back, so
that each node is evaluated once, by one worker; the ancestors of the leaf
at which a worker took its interval are evaluated again by it, for their
bound, but counted apart. With the incumbent fixed, whether a node is pruned
depends on the node alone, and therefore the set of the nodes evaluated is
that of the sequential depth-first search, whatever the steals and the
number of the workers.

The master-hub-worker search over MPI [see `MPIHubBranchAndBound.h`] moves
nodes only in messages (a donation), and the incumbent as in the distributed
search; the open nodes are those of the queues and of the messages in
transit, and its end is detected by the method of Mattern, as for the
distributed search, every message being counted. The reports and the orders
of the balance carry no node, and they only make the messages stop later:
a process with no node reports once, and the hubs and the master answer
only reports, so that when every queue is empty the messages stop and a
wave of Mattern ends the search.

Finally, the racing ramp-up and the checkpoints are correct by the covering
argument of the reoptimization: the open and the closed leaves of the
winner of a race, or of a segment, cover the region that the race, or the
segment, started from, and a closed leaf has a bound not better than the
incumbent by which it was pruned, which is never better than the final one.

## The schemes compared

The table gathers, for each scheme, where the open nodes are, how the end is
detected and on what its outcome depends. Here "timing" means that two runs
with the same data may evaluate different nodes, while the incumbent value is
optimal within the tolerance in each run that is not stopped by a limit.

The MultiQueue B&B (`MultiQueueBranchAndBound`) is the shared queue with the
queue of MQA*: the open nodes are in the MultiQueue of ADS, whose pop takes the
better of the first nodes of two of its $2p$ queues drawn at random, and the
children of a node go into one of them under one lock, so that no lock is
common to all the workers. The termination is a counter $m$ of the nodes in the
queue or held by a worker, which counts the children and drops the node
processed in one atomic step before pushing them, so that $m = 0$ means that no
node is open, as in MQA*. Since the least bound of the queue is not known, a
node is pruned when popped rather than with all the others, and the dive
compares the children with the bound of the node processed. The order of
extraction is relaxed, i.e., a node popped has rank less than
$\frac{5}{6} \cdot 2p$ in expectation among the open ones, which matters when an
evaluation is so fast that the lock of the shared queue is the bottleneck; the
optimum is found within the tolerance, as with the other schemes.

| Scheme | Open nodes | End | Outcome depends on |
|---|---|---|---|
| sequential | one queue | queue empty, or within the tolerance | nothing |
| shared | one queue, under a lock | queue empty and no busy worker | the timing |
| farm | the queue of the master | at the master | the timing |
| deterministic | one queue, by rounds | between two rounds | $r$ |
| distributed | a queue per worker | deficit counting | the timing |
| MPI | a queue per process | four counters of Mattern | the timing |
| stealing | a queue per worker, half stolen | a counter of the busy workers | the timing |
| multiqueue | $2p$ queues, each under a lock | a counter of the open nodes | the timing |
| selfsplit | a queue per worker, no exchange | at each worker | $\sigma$ and $p$ (the timing with `split_shared`) |

We expect the choice among the schemes to depend on the cost of an evaluation
w.r.t. that of moving a node, and on the variance of that cost among the nodes,
which only a measurement on the instances of interest can settle. Note that the
schemes that need no external library can be chosen at run time by their name
(`bb_scheme()` in
[BranchAndBoundScheme.h](../../include/PBeFS/bnb/BranchAndBoundScheme.h)), and
therefore one can compare them with no change of the program.
