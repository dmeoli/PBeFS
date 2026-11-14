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
the closed nodes, which keeps this argument valid; as the sequential A*, they
also skip the move back to the parent in an `OperatorDomain` (the operator
travels with the node, and in HDA* with its message), which never gives a
cheaper path and on the sliding-tile puzzles is about one successor in three,
so that its lookup in CLOSED (and in HDA* its message) is saved. They differ
in how they test the three conditions, which is the delicate part of a
decentralized search, since no worker sees the others at the same instant (cf.
section 3.3.1 of [[1]](../../README.md#ref-1)).

## Centralized searches

In the Simple Parallel A* (SPA*, Algorithm 2 of [[1]](../../README.md#ref-1);
cf. also section 2.2 of [[3]](../../README.md#ref-3)) the workers share one
OPEN and one CLOSED list. Each worker takes the best node of OPEN, expands it
and puts back its successors. The order of expansion is close to the one of A*,
but $p$ nodes are expanded at once, and a node may be expanded before the node
that would give it a cheaper path. This gives some search overhead and some
re-expansions, but no communication overhead, and the load is balanced by
construction. The price is the coordination: OPEN is protected by one lock.
With the data structures of ADS, CLOSED is a set of pointers to the nodes,
compared by their states, in a hash table that is read and filled with no lock
(the `ConcurrentFlatTable` of ADS), and the fields of a node that change are
protected by one of 256 locks, taken only when a state already in CLOSED is
reached with a smaller $g$: since $g$ only decreases, a $g$ read without the
lock that is not larger than the new one means a duplicate, which is the most
frequent case, and takes no lock. With the other data structures, CLOSED is
split into 256 shards by the hash of the state, each with a lock (which also
protects the nodes of its states). Each extraction and each insertion goes
through the lock of OPEN; therefore, SPA* is not expected to
scale beyond a few workers when an expansion is as fast as on the sliding-tile
puzzle, while it may pay off when the expansions are slow. Its termination is
simple, since OPEN and a counter of the busy workers are read under the same
lock: the search ends when OPEN has no node with $f < z$ and no worker is busy.

The MultiQueue A* (`MQAStar`) removes the lock of OPEN by replacing the
priority queue with a relaxed one, the MultiQueue of Rihani, Sanders and
Dementiev [[37]](../../README.md#ref-37), which ADS provides: $2p$ sequential
queues, each under a lock of its own, a push into a queue drawn at random, and
a pop from the better of two queues drawn at random. A pop therefore does not
return the best node of OPEN, but one whose rank among the nodes of OPEN is
less than $\frac{5}{6} \cdot 2p$ in expectation; since SPA* expands $p$ nodes at
once anyway, the order of expansion departs from the one of A* by little more,
while two workers wait for each other only when they draw the same queue (cf.
[[38]](../../README.md#ref-38), where the MultiQueue schedules parallel searches
on graphs). A node popped with $f \ge z$ is dropped instead of ending the
search, since a better node may still be in another queue, and the termination
becomes a counter $m$ of the nodes that are in OPEN or being expanded: a worker
adds the number of its successors to $m$, minus 1 for the node it expanded, in
one atomic step before pushing them, so that $m = 0$ means that OPEN is empty
and that no worker is busy (one step per expansion, rather than one before and
one after the pushes, halves the contention on the line of $m$). The successors
of a node go into one queue of the MultiQueue under one lock, and their places
in CLOSED are prefetched before any of them is looked up, so that the misses of
their lookups overlap; SPA* does the same with CLOSED. With integer costs the
key of a node packs $f$ and the complement of $g$ in one word, $f$ in the high
half, so that the ties of $f$ go to the deeper node as in A*: in the last layer
of $f$, which holds most of the nodes of the 15-puzzle, a goal is then reached
after a part of the layer, while with $f$ alone the pops of the layer come in
any order and most of it is expanded (on the instance 14 of Korf, 68 million
expansions against the 26 million of A*, which MQA* matches with the packed
key).

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
  successor changes worker, i.e., $\mathrm{CO} \approx 1 - \frac{1}{p}$. When
  the Domain gives the features that a move changes
  (`IncrementalFeatureDomain`), HDA* keeps with each node the high 32 bits of
  its hash, which decide its owner and fit in the padding of the node, and the
  key of each successor is that of its parent with the entries of those
  features added, i.e., two lookups on the sliding-tile puzzles instead of one
  per tile, and none for the state that it expands.
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

On shared memory, the end of HDA* is detected by the mailbox of ADS through
which the workers send their nodes, i.e., by the deficit counting of Dijkstra
and Scholten [[35]](../../README.md#ref-35) on a single shared counter $C$,
which holds the number of the workers that are *active* plus the number of
the nodes pushed and not yet processed. A worker adds a batch to $C$ before
the batch can be received, adds 1 when it finds nodes while passive (before
processing the first of them), subtracts the nodes it has processed, and
subtracts its own activity when it has no node with $f < z$ and all its
outboxes are sent; a node is therefore counted from before its push to after
its processing, and a worker that sends is active, so that $C = 0$ only when
no worker has work and no node is in transit, a state from which nothing
raises $C$ again. A batch costs one atomic addition at the sender and one at
the receiver, and a passive worker only reads $C$.

Over MPI (`HDAStarMPI`) the processes share no memory, and the mailbox of the
processes of ADS detects the end by the four-counter method of Mattern
[[11]](../../README.md#ref-11) (cf. section 3.3.1 of
[[1]](../../README.md#ref-1)). A process that has no node with $f < z$ and no
batch to send joins a *wave*, i.e., a nonblocking reduction of the numbers $S$
and $R$ of the messages sent and received by all the processes, and a process
starts a wave only after the previous one has completed; the search ends at the
wave $k$ whose sends $S_k$ equal the receipts $R_{k-1}$ of the wave before.
Since the counters of a process only grow, this equality means that no message
was in transit when the wave $k - 1$ completed, and all the processes, having
the same sums, stop at the same wave. The incumbent travels as a message to all
the processes, and the path is rebuilt at the end state by state, each state
being sent by its owner; this is why the distribution of `HDAStarMPI` has to be
deterministic.

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
scans by default at most the first 64 states for a safe one (the option
`candidates`, 0 for all of them), which bounds the time under the lock and
changes the order of the expansions, not the bound on the cost. Also,
`SlowDomain` adds a given amount of work to each expansion, with which one can
study this regime on any domain.

With a weight $w \geq 1$ OPEN is ordered by $f' = g + w h$, which is the wPA*SE
of section 4 of [[6]](../../README.md#ref-6); the bound $\varepsilon_p$ still
holds if $w \leq \varepsilon_p$, and the library raises $\varepsilon_p$ to $w$
when it is smaller. With several goals, a goal is returned only if, besides
being safe, no state of OPEN or BE has an $f'$ smaller than its $g$, for the
reason given in the next section. The test at the expansion checks only the
states of smaller
$f'$, those of BE included, by the following lemma, which the searches by edges
below use as well.

**Lemma 1.** Let $w \leq \varepsilon$ and $h(t) \leq h(t, s) + h(s)$ for all
$s$ and $t$ (which holds when $h(t)$ is the pairwise heuristic to the nearest
goal). If $g(s) + w h(s) \leq g(t) + w h(t)$, then $s$ is independent of $t$.

*Proof.* $g(s) - g(t) \leq w (h(t) - h(s)) \leq w \, h(t, s) \leq
\varepsilon \, h(t, s)$. $\square$

## Parallel A* for slow evaluations of the edges

When the cost of an expansion is that of evaluating its edges (a collision
check, a simulation), ePA*SE [[121]](../../README.md#ref-121) makes the edge,
not the state, the unit of the parallel work. The library implements it in
`EPASE.h`, on a domain that evaluates one edge at a time (`EdgeDomain`). A
state of OPEN stands for all its edges (the dummy edge of the paper); when it is
safe it is *expanded*, i.e., moved to BE with its edges put in READY, which
costs no evaluation; an edge of READY is evaluated by a worker outside the lock,
and its successor relaxed under it; when all its edges are evaluated, a state
leaves BE. The keys of OPEN, READY and BE are $f' = g + w h$ of the state, with
the $g$ at its expansion for BE and READY, and the search takes the first safe
state of OPEN or the first edge of READY, whichever has the smaller key. In
the terms of Equations 2 and 3 of [[121]](../../README.md#ref-121) an edge is
safe when its state is independent of BE and of the states of OPEN of smaller
key; the edges of a state in BE are always safe, since the state was safe when
it was expanded and its $g$ has not changed. Let $g^*(s)$ be the cost of a
cheapest path to $s$, $v(s)$ the $g$ of $s$ at its expansion, and let the
assumptions of Lemma 1 hold, with an admissible and forward-backward consistent
pairwise heuristic.

**Theorem 2 (bounded expansions).** When a state $s$ is expanded,
$g(s) \leq \varepsilon \, g^*(s)$.

*Proof.* Suppose not, and let $s$ be the first state expanded with
$g(s) > \varepsilon \, g^*(s)$. Let $\pi$ be a cheapest path from $s_0$ to $s$
and $s_m$ the first state of $\pi$ that is in OPEN or in BE (there is one, $s$
itself being in OPEN). The states of $\pi$ before $s_m$ are closed: each was
expanded before $s$, hence with $v \leq \varepsilon g^*$, and its edges were
evaluated, each relaxing its successor with a $g$ not larger than $v$. Thus,
if $s_m = s_0$ then $g(s_m) = 0$, and otherwise, $s_{m-1}$ being its
predecessor on $\pi$,
$g(s_m) \leq v(s_{m-1}) + c(s_{m-1}, s_m) \leq \varepsilon g^*(s_{m-1}) +
c(s_{m-1}, s_m) \leq \varepsilon g^*(s_m)$; if $s_m$ is in BE, also
$v(s_m) \leq \varepsilon g^*(s_m)$ directly. Let $\bar g$ be $g(s_m)$ if
$s_m$ is in OPEN and $v(s_m)$ if it is in BE: $s$ is independent of $s_m$ with
$\bar g$, either by the test at the expansion, if the key of $s_m$ is smaller,
or by Lemma 1. Hence $g(s) \leq \bar g + \varepsilon h(s_m, s) \leq
\varepsilon g^*(s_m) + \varepsilon g^*(s_m, s) = \varepsilon g^*(s)$, where
$g^*(s_m, s)$ is the cost of the part of $\pi$ from $s_m$ to $s$, a
contradiction. $\square$

The proof is that of Theorem 1 of [[121]](../../README.md#ref-121), written for
the test that skips the states of larger key in BE (the improvement of
GePA*SE [[122]](../../README.md#ref-122)) and with $v$ for the states of BE. It
gives the properties of the search:

- a state is expanded at most once, since a state in BE or CLOSED is never put
  back into OPEN, and an edge is evaluated at most once;
- a goal $s$ is returned when it is safe and no state $t$ of OPEN or BE has
  $\bar g(t) + \varepsilon h(t) < g(s)$ (which holds for the states of key not
  smaller than $g(s)$, since $w \leq \varepsilon$); then
  $g(s) \leq \varepsilon C^*$, and $g(s) = C^*$ when $\varepsilon = 1$.
  Indeed, a cheapest path to a goal has a first state $s_m$ in OPEN or BE, with
  $\bar g(s_m) + \varepsilon h(s_m) \leq \varepsilon g^*(s_m) +
  \varepsilon h(s_m) \leq \varepsilon C^*$ by Theorem 2 and the admissibility
  of $h$. Safety alone gives $g(s) \leq \varepsilon g^*(s)$
  (Theorem 3 of [[121]](../../README.md#ref-121)), which is $\varepsilon C^*$
  only if $s$ is the goal of a cheapest path: with a goal region, as in the
  navigation below, a safe goal may be returned while a goal of smaller $f$ is
  still waiting: on the 50 instances of den520d with $w = \varepsilon = 1$ and
  5 workers the test of safety alone returned paths up to 4% longer than the
  optimum on 8 instances with PA*SE, and on den520d and hrt201n with 5 to 20
  threads the code of [[122]](../../README.md#ref-122) returned paths up to 24%
  longer on up to 11 instances out of 50, while with the test above every path
  was optimal;
- while a goal waits for this test, the workers would expand any safe state;
  the library expands no state $x$ with $g(x) + \varepsilon h(x) \geq z$, $z$
  being the least $g$ of a goal in OPEN, which keeps the bound: if the first
  state $s_m$ of a cheapest path is such a state, then $z \leq \varepsilon
  C^*$ by the inequality above, and the goal of $g = z$ is already within the
  bound;
- the search ends: each step either expands a state, which happens at most
  $|V|$ times, or hands out an edge, at most $|E|$ times, or waits; it waits
  only while an edge is being evaluated or no worker is idle, since with no edge
  in READY and none being evaluated BE is empty (a state stays in BE only while
  one of its edges is in READY or being evaluated), and then the first state of
  OPEN is safe; and if a path exists, OPEN or BE holds a state of a cheapest
  path to a goal until a goal is returned, so that the search does not end
  without one.

The thread management of [[121]](../../README.md#ref-121) has one worker, the
dispatcher, select and expand, and gives each edge to the idle worker of least
index, so that the workers beyond those that the independent edges keep busy
are never woken; without it, every worker selects, expands and evaluates under
the lock as in PA*SE. With one worker the search is the sequential eA* of the
paper, which evaluates an edge only when its state is the best: when the
successor reached by the first edge evaluated has a smaller key than its
parent, the other edges of the parent wait, and with an informative heuristic
many of them are never evaluated. GePA*SE [[122]](../../README.md#ref-122)
evaluates the cheap edges of a state (`HeterogeneousEdgeDomain`) as one item
of READY, given to one worker, and each expensive edge apart; the proof above
does not depend on how the edges of a state in BE are grouped.

The anytime search A-ePA*SE [[123]](../../README.md#ref-123) (and A-GePA*SE)
repairs the path with the weights $w_0 > w_1 > \dots > w_k = 1$,
$w_{i+1} = \max\{1, w_i - \delta\}$ and $\varepsilon = w_i$, as ARA*
[[125]](../../README.md#ref-125) does. A state of BE or CLOSED whose $g$ improves
goes to INCON if it was expanded in the current iteration, and back to OPEN
otherwise; a goal selected as safe ends the iteration, which waits for the
edges being evaluated, publishes the best path with the bound $w_i$, and moves
INCON and BE back into OPEN with the keys of $w_{i+1}$. The edges already
evaluated stay with their state and are relaxed again, without a new
evaluation, when the state is expanded again; thus, each edge is evaluated at
most once over all the iterations.

**Theorem 3 (anytime bound).** At the end of the iteration $i$, the best path
costs at most $w_i C^*$, and $C^*$ at the end of the last one.

*Proof.* At the start of an iteration every state whose successors may have
been relaxed with a $g$ larger than its present one is in OPEN (INCON and BE
are moved there), so that every other state $u$ not expanded in the iteration
is *consistent*: its successors $x$ have $g(x) \leq g(u) + c(u, x)$. The proof
of Theorem 2 then holds within the iteration, with $\varepsilon = w_i$, if the
states of $\pi$ before $s_m$ may also be consistent: for such a state
$g(s_j) \leq g(s_{j-1}) + c(s_{j-1}, s_j)$, and by induction from $s_0$, or
from the last state expanded in the iteration, $g(s_j) \leq w_i g^*(s_j)$. A
state of INCON is not in OPEN, but it was expanded in the iteration with
$v \leq w_i g^*$, and the induction uses $v$ for it. Hence the goal selected
at the end of the iteration, safe and with no state of OPEN or BE of key
smaller than its $g$, has $g \leq w_i C^*$, as above (Theorem 1 of
[[123]](../../README.md#ref-123)). $\square$

Note that the iteration ends when a goal is safe, not when its $g$ is below the
least key of OPEN as in ARA*: the latter would also need the least key of BE,
whose edges being evaluated may still lower the $g$ of the goal.

Two baselines complete the family. The Parallel Weighted A* (PwA*,
`PWAStar.h`) of [[121]](../../README.md#ref-121) and
[[124]](../../README.md#ref-124) is the weighted A* whose workers evaluate the
edges of one expansion at a time: its expansions and edges evaluated are those
of the sequential search, whatever the number of workers, and its speedup is
bounded by the branching factor. MPLP [[124]](../../README.md#ref-124)
(`MPLP.h`) is a lazy search: on a domain whose successors are known before
their edges are evaluated (`LazyEdgeDomain`), with a lower bound on the cost of
each edge, the worker 0 repeats a weighted A* that uses the cost of the edges
evaluated and the bound of the others, and records each path it finds, while
the other workers evaluate the edges, those of the paths found first. Since the
bounds do not exceed the costs, each path found costs at most $w C^*$ (Lemma 3
of [[124]](../../README.md#ref-124)), and so does the path accepted, which is one
whose edges are all evaluated and whose cost does not exceed that of any path
found (Theorem 5); if a search finds no path, none exists. The library checks a
path when its last edge is evaluated, instead of on a thread of its own, and
repeats the search only after an evaluation has changed a cost, since otherwise
it would find the same path; with one worker it evaluates the edges of each path
before searching again, which is the lazy shortest path of
[[126]](../../README.md#ref-126).

`SlowEdgeDomain` (in `SlowDomain.h`) makes the edges of any domain cheap or
expensive to evaluate, and some of them missing, so that the regimes of these
papers can be reproduced on the sliding-tile puzzle or on a grid, and
`RobotNav2D.h` is the 2D navigation of their experiments (a square robot of 32
cells a side, 8 motions of 24 cells checked cell by cell at 25 poses, on Moving
AI maps scaled by 5 or 10), whose maps and pairs of points are in
`data/robot-nav-2d`.

## Parallel iterative deepening

IDA* [[7]](../../README.md#ref-7) stores only the path, and parallelizing it
is a matter of distributing its depth-first iterations. `ParallelIDAStar`
first builds a *frontier*: A* (or uniform-cost search, or a tree of fixed
depth without duplicate detection) runs until OPEN holds $N$ live nodes, the
*units*, as in the AIDA* of Reinefeld and Schnecke
[[55]](../../README.md#ref-55) and the A*+IDA* of Bu and Korf
[[56]](../../README.md#ref-56); the iterations then run below the units on
the $p$ workers. Each unit $u$ keeps a threshold $t_u$, initially
$f(u) = g(u) + h(u)$, and its cost-bounded search with a bound
$\beta \geq t_u$ returns the least $f$ of the nodes that it cuts, which
becomes $t_u$. Two schemes share this machinery. In the synchronous one,
each iteration takes $\beta = \min_u t_u$, searches the units with
$t_u \leq \beta$, and ends at a barrier; in the asynchronous one, there is
no barrier, and the workers take the units by increasing threshold, each
searched with its own, while a goal of cost $z$ prunes every node with
$f \geq z$. The rest of this section proves both optimal with an admissible
$h$ and positive costs.

Let $C^*$ be the optimal cost, and let the frontier be built by a best-first
search with reopening (whatever its order) that has not expanded a goal. We
use the invariant of A* recalled in
[Best-first search and A*](best-first-search.md): a closed node was expanded
with its current $g$, and when a node of optimal $g$ is expanded, each of its
successors gets a $g$ not larger than the cost through it. Fix an optimal path
$P$ to a goal, and let $z$ be a node of OPEN with $g(z) = g^*(z)$ that lies on
some optimal path, chosen with $g^*(z)$ largest; one exists, since the first
node of $P$ that is not closed with its optimal $g$ is open with it. Consider
an optimal path $S$ from $z$ to a goal. If a node $w$ of $S$ after $z$ were
in the table of the frontier with $g(w) = g^*(w)$, then either $w$ would be
open, against the choice of $z$, or $w$ would be closed, and the first node of
$S$ after it that is not closed with its optimal $g$ would be open with it,
again against the choice of $z$. Hence:

- the lookups in the table (cf. `table_depth`), which prune a node found
  there with a $g$ not larger, never prune $S$, on which every $g$ is
  optimal and every node of the table has a larger $g$;
- the moves that a tree search skips (cf. `TreeOperatorDomain`) leave a path
  of the same cost and end, whose nodes have $f \leq C^*$ since $h$ is
  admissible, and the search of $z$ starts from the operator stored with
  $z$ only for its inverse, which leads to a node of the table with a
  smaller $g$;
- the threshold of $z$ never exceeds $C^*$ before the goal of $S$ is found:
  initially $t_z = f(z) \leq C^*$, and a search of $z$ with a bound below
  $C^*$ cuts the first node of $S$ whose $f$ exceeds it, whose $f$ is at most
  $C^*$.

In the synchronous scheme, therefore, $\beta \leq t_z \leq C^*$ at every
iteration, a goal found has a cost at most $\beta$, and the first goal found is
optimal; every worker can stop as soon as one is found. A unit with
$t_u > \beta$ need not be searched with $\beta$ (the *unit bounds*): no node
below it has $f$ in $(\beta', t_u)$, $\beta'$ being the bound of its last
search, so that the search would expand the same nodes as the last one and
find nothing. Without unit bounds the scheme is the classical one, which
searches every unit at every bound, as in the parallel IDA* of Sturtevant
and of Siag et al. [[57]](../../README.md#ref-57), and expands at least as
many nodes. In the asynchronous scheme, the threshold of a unit is a lower
bound on the cost of every goal below it not found yet, so that $z$ is
optimal when every threshold is at least $z$, which is when the search ends;
the work that it saves is the idle time at the barriers, and the work that it
may waste is that of the units searched with a threshold above $C^*$ before
the goal is found.

The units are scheduled by work stealing (a range of units per worker, of
which an idle worker steals the upper half) or by one shared cursor, in the
decreasing order of the expansions of their last search, i.e., the longest
processing time first, since the expansions of a unit grow by a similar
factor at each bound. A single unit can still hold most of an iteration, and
the splitting of the stacks of Rao and Kumar
[[58]](../../README.md#ref-58) balances it: a worker that sees an idle one
gives it half of the successors not yet tried at the shallowest level of its
stack, which then reports its least cut $f$ to the unit. `ParallelIDAStarMPI`
runs the synchronous scheme on the processes of MPI: every process builds the
same frontier, an iteration gives the units to the processes by the longest
processing time first, and two reductions of the vectors of the thresholds
and of the expansions end it.

## The searches compared

The following table gathers what each search shares, how it detects the
duplicates and the end, and in which regime it may pay off. The overheads are
those expected from the scheme, and only the measurements on a given machine
and domain (cf. `AStarBench`) can decide among the searches.

| Search | Shared among the workers | Duplicates | End | It may pay off when |
|---|---|---|---|---|
| SPA* | OPEN and CLOSED, under locks | in CLOSED | OPEN and busy workers, under the lock | the expansions are slow |
| MQAStar | OPEN in $2p$ queues and CLOSED, under locks | in CLOSED | a counter of the nodes in OPEN or being expanded | the expansions are too fast for the lock of SPA* |
| FarmAStar | nothing (the master owns all) | at the master | at the master | the expansions are slow |
| HDA* | the incumbent and the batches | at the owner | deficit counting | the expansions are fast and the distribution balances the load |
| HDAStarMPI | the incumbent, by messages | at the owner | four counters of Mattern | the memory of one machine is not enough |
| SafePBNF | the abstract graph, under one lock | in the scope of an nblock | free list empty, no nblock in use | the abstraction has many nblocks with small scopes |
| PA*SE | OPEN, BE and CLOSED, under one lock | in CLOSED, never re-expanded | goal safe, or OPEN and BE empty | the expansions are much slower than the scans |
| ParallelIDAStar | the units and their thresholds | in the frontier, and in its table near it | a barrier per iteration, or no unit below the incumbent | the memory of A* is not enough, and the graph has few short cycles |
| ePA*SE, GePA*SE | OPEN, READY, BE and CLOSED, under one lock | in CLOSED, never re-expanded | goal safe, or OPEN, READY and BE empty | the edges are slow to evaluate and fewer than the workers |
| PwA* | nothing (one expansion at a time) | in CLOSED | as A* | the edges are slow and the branching factor is about the workers |
| MPLP | the queue of the edges and the paths, under one lock | in CLOSED of each lazy search | a path evaluated and short enough, or none in a lazy search | the successors are known before their edges are evaluated |

In the comparison of section 3.2 of [[5]](../../README.md#ref-5), on the
15-puzzle and on the 24-puzzle, HDA*[Z] is faster than SafePBNF on the
instances that take more than about 10 seconds, while SafePBNF is faster on the
easier ones, where the initial burst of HDA* weighs more; the outcome on other
domains and machines is not necessarily the same. Note that the shared-memory
searches are templates over the executor of ADS that runs the workers, and that
they take the same `ParallelOptions`. Hence, one can measure the same search on
the threads, on OpenMP and on FastFlow, and the same runtime across the
searches.
