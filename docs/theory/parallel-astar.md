# Parallel A*

We now turn to the parallel searches of the library, which run A* with $p$
workers on a shared-memory machine (on the threads of C++, of OpenMP or of
FastFlow) or with one MPI process per worker. They follow the taxonomy of the
survey of Fukunaga et al. [[1]](../../README.md#ref-1), and the notation and
the invariants of the sequential search are those of
[Best-first search and A*](best-first-search.md). This page states the
overheads by which a parallel search is measured and the condition under which
it may stop with an optimal path. Then, for each search, it states what the
search guarantees, how it ends and which overhead it trades for which; the
pseudocode of each search is on the page of its header.

## The overheads

Let $E^{\mathrm{seq}}$ be the expansions of the sequential A*, and $E_i$ those
of the worker $i$ of a parallel search. A parallel search pays three overheads
(cf. section 3.1 of [[1]](../../README.md#ref-1), section 2.1 of
[[5]](../../README.md#ref-5) and section 4.2 of [[2]](../../README.md#ref-2)).
The *search overhead*

$$
\mathrm{SO} = \frac{\sum_{i=0}^{p-1} E_i}{E^{\mathrm{seq}}} - 1
$$

counts the expansions that the sequential search does not make, i.e., the nodes
with $f > C^*$, the extra nodes of the last layer and the re-expansions. A
parallel search pays the latter even with a consistent heuristic, since a
worker may expand a node before another worker has found the cheaper path to
it. The *communication overhead*

$$
\mathrm{CO} = \frac{\text{nodes sent to another worker}}
                   {\text{nodes generated}}
$$

counts the nodes that change worker, each of which costs a copy, a message or a
lock. The *load balance*

$$
\mathrm{LB} = \frac{\max_i E_i}{\frac{1}{p} \sum_{i=0}^{p-1} E_i} \geq 1
$$

measures how far the busiest worker is from the mean (1 is a perfect balance).
Finally, the *coordination overhead* is the time that a worker spends waiting
for a lock, for a message or for the others. In the library, the `SearchStats`
of a search give $\mathrm{CO}$ and $\mathrm{LB}$ (by `communication_overhead()`
and `load_balance()`), while $\mathrm{SO}$ is obtained against a sequential
run.

These overheads are not independent. A distribution of the states that keeps
the successors of a node on its worker lowers $\mathrm{CO}$, but it tends to
gather the nodes of least $f$ on a few workers; this raises $\mathrm{LB}$ and,
through it, $\mathrm{SO}$, since the other workers then expand nodes that the
sequential search would not reach (the band effect of section 3.1.1 of
[[5]](../../README.md#ref-5)). Conversely, a uniform distribution balances the
load at the price of sending nearly all the nodes, a fraction $1 - \frac{1}{p}$
of them.

## The termination of a parallel search

A goal selected by one worker is not necessarily optimal, since another worker
may hold a node of smaller $f$, or such a node may be traveling towards its
owner. Hence, the cost $z$ of the best goal found (the *incumbent*) is shared,
and each node with $f \geq z$ is pruned. Therefore, the search may stop only
when (i) no worker has a node with $f < z$ in its OPEN list, (ii) no worker is
expanding a node and (iii) no node is in transit.

The invariant of the sequential search then holds over the union of the OPEN
lists, of the nodes being expanded and of the messages. In fact, if a goal is
reachable, some node $u$ of an optimal path is in one of them with
$g(u) = g^*(u)$ as long as $z > C^*$, and $f(u) \leq C^*$ if the heuristic is
admissible; thus, the three conditions cannot hold while $z > C^*$, and the
search stops with $z = C^*$. At the end, the path is rebuilt from the goal of
the incumbent through the parents, after all the workers have finished. All the
parallel searches of the library select the goals at the expansion and reopen
the closed nodes, which keeps this argument valid. They differ in how they test
the three conditions, which is the delicate part of a decentralized search,
since no worker sees the others at the same instant (cf. section 3.3.1 of
[[1]](../../README.md#ref-1)).

## Centralized searches

In the Simple Parallel A* (SPA*, Algorithm 2 of [[1]](../../README.md#ref-1);
cf. also section 2.2 of [[3]](../../README.md#ref-3)) the workers share one
OPEN and one CLOSED list. Each worker takes the best node of OPEN, expands it
and puts back its successors. The order of expansion is close to the one of A*,
but $p$ nodes are expanded at once, and a node may be expanded before the node
that would give it a cheaper path. This gives some search overhead and some
re-expansions, but no communication overhead, and the load is balanced by
construction. The price is the coordination: OPEN is protected by one lock, and
CLOSED is split into 256 shards by the hash of the state, each with a lock
(which also protects the nodes of its states). Each extraction and each
insertion goes through the lock of OPEN; therefore, SPA* is not expected to
scale beyond a few workers when an expansion is as fast as on the sliding-tile
puzzle, while it may pay off when the expansions are slow. Its termination is
simple, since OPEN and a counter of the busy workers are read under the same
lock: the search ends when OPEN has no node with $f < z$ and no worker is busy.

The farm with feedback (`FarmAStar`) is the master-worker form of the same
scheme (cf. section 3.2 of [[1]](../../README.md#ref-1)), written with FastFlow
[[12]](../../README.md#ref-12). A master owns OPEN, CLOSED and the nodes, and
hence it needs no lock; it sends the best nodes of OPEN, in batches, to the
free workers, which generate the successors with their costs and heuristic
values and send them back. Since the master does all the bookkeeping (i.e., the
detection of the duplicates and the updates of OPEN), it is a bottleneck unless
an expansion costs much more than the bookkeeping of its successors. Its
termination is local to the master, which ends the search when no node of OPEN
has $f < z$ and no batch is out at a worker.

## Hash Distributed A*

In Hash Distributed A* (HDA*, section 3 of [[2]](../../README.md#ref-2) and
section 4.1 of [[1]](../../README.md#ref-1)) each worker runs A* on its own
OPEN and CLOSED lists. A function $\mathrm{owner}(u)$ with values in
$\{0, \ldots, p - 1\}$ assigns each state to one worker: a successor generated
by a worker is sent to its owner, which checks it against its CLOSED list and
possibly puts it in its OPEN list. Thus, the duplicates are detected exactly,
since all the copies of a state reach the same worker, and the workers share no
list and take no lock on a node. The nodes travel in batches (of 64 nodes by
default); an outbox that is not full is sent after a given number of expansions
and before the worker waits, since packing many nodes in a message lowers the
cost per node (section 5.2 of [[2]](../../README.md#ref-2)). Within each worker
the order of expansion is that of A*, but each worker selects its best node
only among those that it owns. This gives the search overhead of the load
imbalance and of the initial burst, when the OPEN lists are still empty
(section 3.1 of [[5]](../../README.md#ref-5)).

The behavior of HDA* is decided by the distribution, which is a template
parameter (cf. [Distribution.h](../../include/PBeFS/astar/Distribution.h)). In
the notation $\mathrm{HDA}^*[\mathit{hash}, \mathit{abstraction}]$ of section
2.6 of [[5]](../../README.md#ref-5), the library gives the following ones.

- Zobrist hashing, HDA*[Z] (section 6.2 of [[1]](../../README.md#ref-1)): with
  the state described by the features $x_0, \ldots, x_n$ that hold in it and a
  table $R$ of random 64-bit strings, one per feature,
  $Z(s) = R[x_0] \oplus \cdots \oplus R[x_n]$, with $\oplus$ the exclusive or.
  It spreads the states uniformly over the workers (it is 3-independent) and
  can be updated incrementally by a move; the load is balanced, but nearly each
  successor changes worker, i.e., $\mathrm{CO} \approx 1 - \frac{1}{p}$.
- Abstract Zobrist hashing, AZHDA* or HDA*[Z, Afeature] (section 4 of
  [[5]](../../README.md#ref-5)): with a projection $A$ of the features onto
  fewer abstract features, $AZ(s) = R[A(x_0)] \oplus \cdots \oplus R[A(x_n)]$.
  A move that changes no abstract feature keeps the successor on the worker of
  its parent, which lowers $\mathrm{CO}$, while the hashing of all the abstract
  features keeps the load balanced; on the domains of section 4.1 of
  [[5]](../../README.md#ref-5), the sliding-tile puzzles among them, it is the
  best of these distributions.
- Abstraction, AHDA* or HDA*[Z, Astate] (cf. section 6.4 of
  [[1]](../../README.md#ref-1)): all the states with the same abstract state go
  to the same worker, which keeps $\mathrm{CO}$ low, but nothing balances the
  load, and the search overhead grows.
- Perfect hashing, HDA*[P]: the rank of the state modulo $p$, which on the
  sliding-tile puzzles depends on the last few tiles only and gives both a high
  search and a high communication overhead (sections 3.2 and 4.1 of
  [[5]](../../README.md#ref-5)); it is kept to reproduce that comparison.
- The hash of the closed list, remixed with a seed, which behaves as Zobrist
  hashing for any domain but cannot be updated incrementally.
- A random worker for each generated state, as in the randomized strategy of
  Kumar, Ramesh and Rao (cf. section 6 of [[2]](../../README.md#ref-2)). A
  duplicate is then detected only when its copies happen to reach the same
  worker, and this is the baseline without detection of the duplicates across
  the workers.

On shared memory, the end of HDA* is detected by three shared counters (cf.
`ThreadMailbox` in
[ThreadMailbox.h](../../include/PBeFS/core/ThreadMailbox.h)). Here *pending*
holds the nodes sent and not yet processed (a worker raises it before a batch
becomes visible, and lowers it after having processed the batch), *idle* holds
the workers waiting, and *epoch* is raised by a worker that stops waiting,
after it has left *idle*. A waiting worker declares the end when it reads
$\mathit{idle} = p$ and $\mathit{pending} = 0$ while $\mathit{epoch}$ does not
change. At the instant at which *idle* was read nobody had work, and no node
could be in transit, since its receiver would have had to stop waiting, i.e.,
to change *epoch*.

Over MPI (`HDAStarMPI`) the processes share no memory, and the end is detected
by the four-counter method of Mattern [[11]](../../README.md#ref-11) (cf.
section 3.3.1 of [[1]](../../README.md#ref-1)). A process that has no node with
$f < z$ and no batch to send joins a *wave*, i.e., a nonblocking reduction of
the numbers $S$ and $R$ of the messages sent and received by all the processes;
the search ends when two consecutive waves give the same totals with $S = R$.
Since the counters of a process only grow, equal totals in two waves mean that
no message was sent or received between them, and $S = R$ that none is in
transit. All the processes see the same totals, and therefore they stop at the
same wave. The incumbent travels as a message to all the processes, and the
path is rebuilt at the end state by state, each state being sent by its owner;
this is why the distribution of `HDAStarMPI` has to be deterministic.

## Safe Parallel Best-NBlock-First

Safe PBNF (section 3 and appendix A of [[3]](../../README.md#ref-3)) removes
the locks from the detection of the duplicates by another partition. An
abstraction maps the states to the nodes of an abstract graph, and the states
with the same abstract state form an *nblock*, with its own OPEN and CLOSED
lists. The children of a node of the nblock $b$ fall in $b$ or in its abstract
successors, which form the *duplicate detection scope* of $b$. A worker
acquires an nblock whose scope is disjoint from the scopes of the nblocks in
use, and then it expands its nodes and checks their children for duplicates
without any lock (as in the parallel structured duplicate detection of Zhou and
Hansen; cf. section 2.4 of [[3]](../../README.md#ref-3)).

In Safe PBNF, the nblocks that can be acquired form the *free list*, ordered by
the least $f$ of their OPEN lists. A worker keeps its nblock while it is not
worse than the best free one and than the best nblock that its use keeps from
being free, after a minimum number of expansions (32 by default); otherwise, it
releases it and acquires the best free one. Therefore, the order of expansion
is best-first only approximately, which gives some search overhead, while the
communication is limited to the counters of the abstract graph, under one lock.

Plain PBNF may livelock on an infinite state space, since the nblock that holds
a goal may never become free (section 3.1 of [[3]](../../README.md#ref-3)).
Safe PBNF prevents this by the *hot* nblocks. The *interference scope* of $b$
holds the nblocks that, if in use, keep $b$ from being free (the predecessors
of the successors of $b$); $\sigma(b)$ counts the nblocks in use whose
interference scope contains $b$, $\sigma_h(b)$ counts the hot nblocks with $b$
in their interference scope, and an nblock is free only if
$\sigma(b) = \sigma_h(b) = 0$. A worker whose nblock interferes with a better
one marks the latter hot; the nblocks that interfere with a hot one leave the
free list, and their workers release them at their next check, so that each hot
nblock eventually becomes free. Burns et al. checked the absence of livelock on
a model of the algorithm (section 3.2 of [[3]](../../README.md#ref-3)).

The search ends when no nblock is in use and the free list is empty, the
nblocks without nodes below $z$ having left the free list. Since each node of
an nblock is in its OPEN list or held by the worker of the nblock, this is the
condition stated above. In the library, the search follows the pseudocode of
appendix A of [[3]](../../README.md#ref-3), but the best $f$ of an nblock in
use is a value that its worker publishes after each expansion; the decision to
switch may then read a value slightly out of date, which changes the order of
expansion but not the cost of the path found.

## Parallel A* for slow expansions

PA*SE [[6]](../../README.md#ref-6) expands in parallel only the states whose
$g$ is already optimal, at the price of a test before each expansion; hence, no
state is ever expanded twice. Let $h(u, v)$ be a pairwise heuristic that never
overestimates the cost of a path from $u$ to $v$ and satisfies
$h(u, w) \leq h(u, v) + h(v, w)$ (forward-backward consistency). A state $s$ is
*independent* of a state $s'$ if

$$
g(s) - g(s') \leq \varepsilon_p \, h(s', s) , \qquad \varepsilon_p \geq 1 ,
$$

since then no path through $s'$ can lower $g(s)$ below $g(s) / \varepsilon_p$.
A state $s$ is *safe* if it is independent of the states being expanded (the
set BE) and of the states of OPEN with a smaller $f$, which by Theorem 2 of
[[6]](../../README.md#ref-6) are the only states of OPEN to be checked.

Each worker takes from OPEN the safe state of least $f$, adds it to BE, expands
it outside the lock and puts back the successors whose $g$ has improved. With
$\varepsilon_p = 1$ each state is expanded once and with its optimal $g$, and
the path returned is optimal; with $\varepsilon_p > 1$ its cost is within a
factor $\varepsilon_p$ of the optimum (Theorems 1 and 3 of
[[6]](../../README.md#ref-6)). The search ends when the goal is safe, or when
OPEN and BE are both empty. A safe state costs a scan of BE and of the states
of OPEN before it, all under one lock; therefore, PA*SE pays off only when an
expansion is much slower than these scans (say, in motion planning with
collision checks). In the library OPEN is an ordered set, of which the search
scans at most the first 64 states for a safe one, which bounds the time under
the lock. Also, `SlowDomain` adds a given amount of work to each expansion,
with which one can study this regime on any domain.

## The searches compared

The following table gathers what each search shares, how it detects the
duplicates and the end, and in which regime it may pay off. The overheads are
those expected from the scheme, and only the measurements on a given machine
and domain (cf. `AStarBench`) can decide among the searches.

| Search | Shared among the workers | Duplicates | End | It may pay off when |
|---|---|---|---|---|
| SPA* | OPEN and CLOSED, under locks | in CLOSED | OPEN and busy workers, under the lock | the expansions are slow |
| FarmAStar | nothing (the master owns all) | at the master | at the master | the expansions are slow |
| HDA* | the incumbent and the batches | at the owner | three counters | the expansions are fast and the distribution balances the load |
| HDAStarMPI | the incumbent, by messages | at the owner | four counters of Mattern | the memory of one machine is not enough |
| SafePBNF | the abstract graph, under one lock | in the scope of an nblock | free list empty, no nblock in use | the abstraction has many nblocks with small scopes |
| PA*SE | OPEN, BE and CLOSED, under one lock | in CLOSED, never re-expanded | goal safe, or OPEN and BE empty | the expansions are much slower than the scans |

In the comparison of section 3.2 of [[5]](../../README.md#ref-5), on the
15-puzzle and on the 24-puzzle, HDA*[Z] is faster than SafePBNF on the
instances that take more than about 10 seconds, while SafePBNF is faster on the
easier ones, where the initial burst of HDA* weighs more; the outcome on other
domains and machines is not necessarily the same. Note that the shared-memory
searches are templates over the runtime that starts the workers, and that they
take the same `ParallelOptions`. Hence, one can measure the same search on the
threads, on OpenMP and on FastFlow, and the same runtime across the searches.
