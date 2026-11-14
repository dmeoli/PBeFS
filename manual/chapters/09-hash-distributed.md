# 9. Hash-distributed search {#ch-9}

[Source: `include/PBeFS/astar/HDAStar.h` (`HDAStar`,
`ParallelOptions`), `include/PBeFS/astar/HDAStarMPI.h` (`HDAStarMPI`),
`include/PBeFS/astar/Distribution.h` (the work distributions),
`include/PBeFS/astar/Domain.h` (`FeatureDomain`,
`IncrementalFeatureDomain`, `AbstractFeatureDomain`,
`AbstractStateDomain`, `PerfectHashDomain`),
`include/PBeFS/core/SearchUtils.h` (`ParallelParams`)]

We are concerned in this chapter with Hash Distributed A$^*$ (HDA$^*$)
of Kishimoto, Fukunaga and Botea [[2]](Z-bibliography.md#ref-2), the
parallel search in which each state has an *owner* among the $p$
workers, fixed by a function of the state alone. Only the owner stores
the state; hence, the workers share no OPEN and no CLOSED list, and yet
the duplicates are detected as in the sequential search. In the
searches of [Chapter 8](08-shared-memory.md#ch-8) (from the Parallel
A$^*$ with one shared OPEN list onwards) each node costs a lock, or an
atomic operation, on a structure that all the workers touch. HDA$^*$
replaces them by messages between pairs of workers, at the price of the
nodes that change worker and of an order of expansion that is only (and
often rather loosely) that of A$^*$. We use the notation of
[Chapter 1](01-introduction.md#ch-1) and the invariants of A$^*$ of
[Chapter 2](02-best-first.md#ch-2); the search overhead (SO), the
communication overhead (CO) and the load balance (LB) of a parallel
search, and the general question of when a parallel search may stop,
are those of [Chapter 8](08-shared-memory.md#ch-8). The structure
of the chapter is as follows. [§9.1](#sec-9-1) states the interface
(the search, its options, the concept of a distribution and the
refinements of the `Domain` that the distributions need),
[§9.2](#sec-9-2) gives the algorithm with its invariants and the proof
of its optimality, and [§9.3](#sec-9-3) the messages, their batches and
the detection of the end on a shared memory. [§9.4](#sec-9-4) and
[§9.5](#sec-9-5) are devoted to the distributions and to the
incremental Zobrist hash, [§9.6](#sec-9-6) to the search over MPI with
the termination detection of Mattern, and [§9.7](#sec-9-7) to the
engineering of the implementation; the chapter ends with the benchmark,
the notes and the exercises.

## 9.1 The search and its interface {#sec-9-1}

A search is an object of the class template

```text
HDAStar< D , Dist = StateHashDistribution< D > ,
         Executor = ThreadPool ,
         Open = buckets if the costs of D are integers, a heap otherwise >
```

(`HDAStar.h:138-144`), built from a `Domain` `D`, a distribution `Dist`
and the options, and run by `solve()`, which returns the `Solution` (the
status, the cost, the path and the counters of [Chapter
2](02-best-first.md#ch-2)). Its workers run on an `Executor`, i.e., on a
pool of C++ threads (the default), on OpenMP or on FastFlow
[[12]](Z-bibliography.md#ref-12). One can give the search an executor
that outlives it by `set_executor()` (`HDAStar.h:183`), which saves the
start of the threads at each search; otherwise, `solve()` builds one of
$p$ threads. The options are a `ParallelOptions`, i.e., the options of
the sequential search plus the parameters of the parallel ones
(`SearchUtils.h:94-111`): the number $p$ of workers (`threads`, 1 by
default), the number of nodes packed in a message (`batch`, 64), the
number of expansions between two flushes of the outboxes that are not
full (`flush_period`, 64) and a `seed`, which `HDAStar` does not read (a
distribution takes its seed in its constructor, 1 by default). A
parameter equal to 0 is taken as 1 (`normalize()`), and therefore no
setting of the options can stop the exchange of the nodes. Also, the
same search with one MPI process per worker is `HDAStarMPI< D , Dist >`
(`HDAStarMPI.h:111-112`), described in [§9.6](#sec-9-6).

### The distributions

A *distribution* is a class with a method `owner( s , p )` that returns
the worker in $\{0, \ldots, p - 1\}$ which owns the state $s$, and with
a static method `name()` that gives its name in the notation
$\mathrm{HDA}^*[\mathit{hash}, \mathit{abstraction}]$ of Jinnai and
Fukunaga (section 2.6 of [[5]](Z-bibliography.md#ref-5)). We remark that
the search on a shared memory requires no more of a distribution; in
particular, the owner need not be deterministic there (the random
distribution of [§9.4](#sec-9-4) is not). `HDAStarMPI`, on the
contrary, needs a deterministic owner, since it rebuilds the path by
asking the owner of each state (cf. [§9.6](#sec-9-6)). A distribution
that can follow the moves models the concept `IncrementalDistribution`
(`Distribution.h:391`): it also gives a 32-bit *key* `key( s )`, whose
owner among $p$ workers is `key_owner( key( s ) , p )`, and the key
`key( k , s , t )` of a successor $t$ of $s$ computed from the key $k$
of $s$, which must equal `key( t )`. HDA$^*$ keeps the key of each
node in the node, and the function `successor_owners()`
(`Distribution.h:420`) returns, for the node being expanded, a function
that maps each successor to its owner and its key, computed from the key
of the parent when the distribution is incremental and from scratch otherwise.

### What the distributions need from the Domain

Each distribution requires of the `Domain` the least that it uses, through
a refinement of the concept `Domain` of
[Chapter 1](01-introduction.md#ch-1) (all in `Domain.h`).

* `FeatureDomain` (`Domain.h:129`): a state is a set of propositions,
  the *features*, taken from a universe of `num_features()` of them,
  and `for_each_feature( s , f )` calls `f( i )` for each feature $i$
  that holds in $s$. Zobrist hashing needs it.
* `IncrementalFeatureDomain` (`Domain.h:145`): moreover,
  `for_each_feature_change( s , t , f )` calls `f( i )` for each feature
  that holds in exactly one of $s$ and of its successor $t$, i.e., for
  the features of the symmetric difference of the two sets.
* `AbstractFeatureDomain` (`Domain.h:157`): moreover, a projection
  `abstract_feature( i )` of the features onto
  `num_abstract_features()` abstract features, which is the projection
  $A$ of Abstract Zobrist hashing.
* `AbstractStateDomain` (`Domain.h:170`): a map `abstract_state( s )`
  of the states to the integers, i.e., a state abstraction.
* `PerfectHashDomain` (`Domain.h:234`): an injective map
  `perfect_hash( s )` of the states to the integers.

On the sliding-tile puzzles (`SlidingTile.h:309-376`) a feature is a
pair (tile, cell), i.e., the index $t \cdot n + c$ for the tile
$t \ne 0$ in the cell $c$ of a board of $n$ cells. An abstract feature
is a pair (tile, half of the board), the half being the top rows or the
bottom rows, and the abstract state of a state is the triple of the
cells of the tiles 1, 2 and 3, which is the abstraction that Jinnai and
Fukunaga chose among the candidates of footnote 4 of section 3.1.4 of
[[5]](Z-bibliography.md#ref-5). Furthermore, the grids and the vacuum
world give features and abstractions as well (cf. [Chapter
14](14-domains.md#ch-14)), and any `Domain` gives the hash of its closed
list, which is all that the default distribution needs (one can
therefore run HDA$^*$ on any `Domain`).

## 9.2 Hash Distributed A$^*$ {#sec-9-2}

We now describe the search. Each worker $i$ runs A$^*$ on its lists
$\mathrm{OPEN}_i$ and $\mathrm{CLOSED}_i$, which hold only the states $v$
with $\mathrm{owner}(v) = i$. A successor generated by a worker is sent to
its owner, which checks it against its CLOSED list and possibly
inserts it in its OPEN list, as the sequential A$^*$ does with a state
that it generates. The *incumbent* $z$, i.e., the cost of the best goal
found so far, is shared by the workers, and each node with $f \ge z$
is pruned. We write a message as $(v, g', u)$: it carries a state $v$ with
the cost $g'$
of the path that reached it and its predecessor $u$ on that path (and,
in the implementation, the operator that generated $v$, its heuristic
and its key: the `Msg` of `HDAStar.h:158-165`). Its pseudocode is that
of the header (`HDAStar.h:18-46`).

**Algorithm 9.1** (HDA$^*$).

```text
 1  procedure HDA*(s₀, p, owner)
 2     z ← ∞;  INSERT(s₀, 0, nil) on owner(s₀)
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        while the search is not over do
 5           for each message (v, g', u) received do INSERT(v, g', u)
 6           if OPEN_i has a node of f < z then
 7              u ← the first node of OPEN_i, removed (stale entries
                   dropped);  close u
 8              if u is a goal then
 9                 z ← min(z, g(u)), recording u if g(u) < z
10              else for each successor v of u do
11                 g' ← g(u) + c(u, v)
12                 if g' + h(v) ≥ z then continue
13                 if owner(v) = i then INSERT(v, g', u)
14                 else send (v, g', u) to owner(v)
15              every k expansions, send the outboxes not yet full
16           else
17              send all the outboxes;  wait for mail or for the end
18     return the path to the goal of cost z, if z < ∞
19  procedure INSERT(v, g', u)                       ▷ on the owner i of v
20     if g' + h(v) ≥ z then return
21     if v ∉ CLOSED_i then
22        add v to CLOSED_i, with g(v) ← g' and parent(v) ← u
23     else if g' < g(v) then
24        g(v) ← g';  parent(v) ← u;  v is not closed
25     else return                                         ▷ a duplicate
26     push (g' + h(v), g', v) in OPEN_i
```

We comment on the lines that are not those of A$^*$.

* Line 2: the root is inserted by its owner before the first expansion
  (`HDAStar.h:194-199` and `:267-268`), and the other workers start with
  empty lists and wait for mail; the first expansions are therefore
  those of a single worker, which is the initial burst of the search
  (section 3.1.2 of [[5]](Z-bibliography.md#ref-5)).
* Line 5: the messages are processed before each expansion; hence, a
  cheaper path that has arrived is used at once; the implementation
  drains the whole inbox of the worker and inserts its nodes in one
  pass (cf. [§9.7](#sec-9-7)).
* Line 7: CLOSED holds all the nodes that a worker has stored, closed or
  not, and OPEN is a priority queue of pairs (priority, node) with lazy
  deletion: when the $g$ of a node decreases (line 24) a new entry is
  pushed, and the old one becomes *stale* and is dropped when it reaches
  the top (`HDAStar.h:314-328`), as in the sequential search of
  [Chapter 2](02-best-first.md#ch-2).
* Lines 8 and 9: the goal test is made at the expansion, not at the
  generation, and therefore a goal is recorded only when it is the best node
  of its owner; the incumbent is an atomic variable, and the goal is
  recorded under a mutex (`HDAStar.h:338-345`).
* Lines 12 and 20: a node with $f \ge z$ is pruned twice, by the worker
  that generates it (which saves the message) and by its owner (since
  $z$ may have decreased while the node was travelling).
* Lines 13 and 14: a successor owned by the expanding worker does not go
  through the mailbox, while one owned by another worker is appended to
  the outbox of that worker, which leaves as one message when it holds
  `batch` nodes (cf. [§9.3](#sec-9-3)).
* Lines 15 and 17: an outbox that is not full leaves after each block of
`flush_period` expansions, and all the outboxes leave before the worker
  waits; the wait ends when mail arrives or when the mailbox detects the
  end of the search, i.e., when no worker has a node of $f < z$ and no
  node is travelling ([§9.3](#sec-9-3)).
* Line 18: the parents cross the workers (a node records a pointer to
  the node of its predecessor, which may belong to another worker), and
  they are followed only after all the workers have returned
  (`HDAStar.h:205-227`).

As the sequential search, HDA$^*$ skips the move back to the parent in
an `OperatorDomain`: the operator that generated a state travels in its
message, and the successor that the inverse operator leads to is never
generated, since it is never cheaper than the parent (cf.
[Chapter 2](02-best-first.md#ch-2)); on the sliding-tile puzzles this
typically saves about one successor out of three and, in HDA$^*$, its
message as well.

![Figure 9.1. One expansion of HDA$^*$ with three workers. Before
(left), worker 0 takes $u$, the first node of its OPEN list; after
(right), $u$ is closed, its successor $v_1$, which worker 0 owns, is in
its OPEN list, and $v_2$ and $v_3$ wait in the outboxes of the workers
1 and 2 (the colour of a state is its owner, and the entries added by the
step are in the darker rust
colour).](../figures/hash-distributed-step.svg){#fig-9-1}

[Figure 9.1](#fig-9-1) shows an expansion, and [Figure 9.2](#fig-9-2)
what the owner of a state does with the copies of it that the other
workers send. All the copies of $v$ reach the same worker, which keeps in
its CLOSED list the cheapest of those that it has processed. Hence, the
duplicate detection is that of the sequential search. The difference
is that the copies arrive in an order that depends on the
timing of the workers, and a copy that arrives after its state has been
expanded with a larger $g$ reopens it. This is why HDA$^*$ may reopen a
state even when the heuristic is consistent (one can build a small
example, cf. [Exercise 9.5](#sec-9-10)) (footnote 6 of section 3 of
[[2]](Z-bibliography.md#ref-2)), and the reopenings are part of its
search overhead.

![Figure 9.2. The copies of a state $v$, owned by worker 2, sent by
three other workers. (a) The first copy creates the node of $v$, with
$g = 5$ and parent $u_0$, and pushes it with $f = 9$; (b) a cheaper copy
lowers $g$ and changes the parent, and the old entry of OPEN becomes
stale; (c) a copy that is not cheaper is dropped as a
duplicate.](../figures/hash-distributed-duplicate.svg){#fig-9-2}

In the header, the insertion is the following excerpt; `find_or_insert()`
of the closed list (cf. [Chapter 3](03-data-structures.md#ch-3))
returns the node of the state and whether it has just been created, the
node being built only in the latter case.

```cpp
--8<-- "snippets/hash-distributed-insert.cpp"
```

### Invariants

We now state the invariants of Algorithm 9.1, from which its optimality
follows. We assume that the costs are nonnegative and the heuristic is
admissible, as
in [Chapter 2](02-best-first.md#ch-2), and $C^*$ is the cost of an
optimal path from $s_0$ to a goal. The first invariant is the reason why the
workers need no lock on the nodes.

<a id="lem-9-1"></a>**Lemma 9.1** (Ownership). At any time, the node
of a state $v$ is created, stored and written only by the worker
$\mathrm{owner}(v)$; if the distribution is deterministic, at most one
node of $v$ exists, and each copy of $v$ that is generated is compared with it.

*Proof.* A node is created only by INSERT (lines 21-22), which a worker
calls only on a state that it owns: either on the root, by its owner
(line 2), or on a message, which was sent to the owner of its state
(line 14), or on a successor that it owns (line 13). A node is written
by INSERT (lines 22-24) and by the expansion (line 7), and a worker
expands only the nodes of its OPEN list, which INSERT fills. The
other workers only hold the address of a node, as the parent in a
message, and they do not dereference it before the end (line 18). If
the owner is a function of the state, all the copies of $v$ reach the
same worker, whose CLOSED list holds at most one node of $v$ (line 21).
$\square$

Without a deterministic owner (the random distribution of
[§9.4](#sec-9-4)) the lemma holds for each node, but a state may have up
to $p$ nodes, one per worker that has received a copy of it; the
duplicates are then detected only among the copies that happen to reach
the same worker, a weakness that the hashing of the state was introduced
to remove (cf. section 6 of [[2]](Z-bibliography.md#ref-2)).

<a id="lem-9-2"></a>**Lemma 9.2** (Values and parents). At any time,
(i) $z$ is nonincreasing; (ii) for each stored node $v \ne s_0$ with
parent $u$, one has $g(v) \ge g(u) + c(u, v)$, and $g(v)$ only
decreases; (iii) the parent pointers form a tree rooted at $s_0$; and
(iv) if $z < \infty$, the chain of the parents of the recorded goal $t$
is a path from $s_0$ to $t$ of cost at most $g(t) = z$.

*Proof.* (i) $z$ changes only at line 9, to a smaller value, under the
mutex. (ii) When INSERT sets $g(v) \gets g'$ and $\mathrm{parent}(v)
\gets u$, one has $g' = g(u) + c(u, v)$ with the value of $g(u)$ at the
expansion of $u$ (line 11). Afterwards $g(u)$ may only decrease (line 24
is reached only with $g' < g(v)$), which keeps the inequality, and the
parent of $v$ changes only together with a smaller $g(v)$. (iii)
Suppose that setting $\mathrm{parent}(v) \gets u$ closes a cycle, i.e.,
that the chain of the parents of $u$ leads to $v$. Summing (ii) along
that chain gives $g(u) \ge g(v)$ for the value of $g(v)$ before the
update, since the costs are nonnegative; but then $g' = g(u) + c(u, v)
\ge g(v)$, and line 24 is not reached, a contradiction. (iv) The goal
recorded at line 9 has $g(t) = z$ at that moment; summing (ii) along the
chain of its parents, which ends at $s_0$ by (iii), with $g(s_0) = 0$,
gives a path of cost at most $g(t)$. $\square$

Our next invariant is the parallel form of the one that makes A$^*$
optimal, i.e., that some node of an optimal path is always in OPEN with
its optimal $g$ (cf. [Chapter 2](02-best-first.md#ch-2)); here the node
may also be in the hands of a worker or in a message, as
[Figure 9.3](#fig-9-3) shows.

<a id="lem-9-3"></a>**Lemma 9.3** (A node of an optimal path is always at
hand). Let $P = (s_0 = u_0, u_1, \ldots, u_k)$ be an
optimal path to a goal. While $z > C^*$, some $u_j$ of $P$ with $g(u_j)
= g^*(u_j)$ is (a) in the OPEN list of its owner, in an entry that is not
stale, or (b) being expanded, or (c) in a message (in an outbox, in
transit, or in an inbox not yet processed) that carries $g' = g^*(u_j)$.

*Proof.* Note first that each $u_j$ has
$f(u_j) = g^*(u_j) + h(u_j) \le C^* < z$, by admissibility; hence, no
copy of it with its optimal $g$ is ever pruned at line 12 or at line 20.
At the start, $s_0$ is in case (a), or about to be inserted by its owner
(line 2). Afterwards, let $u_j$ be the node of $P$ of largest index that
is closed with $g(u_j) = g^*(u_j)$ (if there is none, $s_0$ has not been
expanded yet and is in case (a) or (b)). If the expansion of $u_j$ with
that $g$ is still in progress, $u_j$ is in case (b). Otherwise, since
$u_j$ is not the goal (otherwise line 9 would have set $z \le C^*$), its
expansion generated $u_{j+1}$ with
$g' = g^*(u_j) + c(u_j, u_{j+1}) = g^*(u_{j+1})$ (line 11), and this
copy was either given to INSERT at line 13 or sent at line 14. In the
first case, INSERT either stored $u_{j+1}$ with $g^*(u_{j+1})$ and
pushed it, or found $g(u_{j+1}) \le g'$, i.e.,
$g(u_{j+1}) = g^*(u_{j+1})$ already; in both cases $u_{j+1}$ is open
with its optimal $g$, or being expanded, or closed with its optimal $g$,
and the last contradicts the choice of $j$. In the second case the copy
is in a message, case (c), until its owner processes it at line 5, after
which the argument of the first case applies. $\square$

![Figure 9.3. The three places of the node $u_3$ of an optimal path
in [Lemma 9.3](#lem-9-3), when its predecessors are closed with their
optimal $g$ (darker tint) and $z > C^*$: (a) in the OPEN list of its owner, (b)
being expanded, (c) in a message of the worker that expanded $u_2$ to
the owner of $u_3$.](../figures/hash-distributed-invariant.svg){#fig-9-3}

The end of the search is detected by the mailbox ([§9.3](#sec-9-3)). For
the moment we take its correctness as a hypothesis, i.e., we assume that
the workers stop only when no worker has a node of $f < z$ in its OPEN
list, no worker is expanding a node and no message exists (the three
conditions of the end of a parallel search of
[Chapter 8](08-shared-memory.md#ch-8)), and that they do stop once these
hold; [Lemma 9.6](#lem-9-6) proves it on a shared memory and
[Theorem 9.13](#thm-9-13) over MPI.

<a id="thm-9-4"></a>**Theorem 9.4** (Optimality of HDA$^*$). Let the
graph be finite, the costs nonnegative, the heuristic admissible and the
distribution deterministic, and let the detection of the end be
correct. Then HDA$^*$ terminates; if a goal is reachable, it returns a
path of cost $C^*$, and otherwise it reports that there is none.

*Proof.* Termination. Each insertion of a state with a given $g$ is
followed by at most one expansion of it with that $g$, and the $g$ of a
state only decreases ([Lemma 9.2](#lem-9-2)(ii)). Note that each value of
$g$ is the cost of a walk from $s_0$ that is not larger than the first value
assigned. On a finite graph with nonnegative costs these costs form a
finite set, since a walk of cost at most $B$ has at most $B / \varepsilon$
edges of positive cost, $\varepsilon$ being the least positive cost, and
its cost is a sum of at most that many positive costs. Hence the
expansions are finitely many, and so are the messages, each of which is
generated by an expansion; after the last of them is processed no worker
has work and no message is in transit, and the detection of the end stops
the workers.

Optimality. At the end the three conditions hold; if a goal is
reachable and $z > C^*$, [Lemma 9.3](#lem-9-3) gives a node with $f \le
C^* < z$ in an OPEN list, being expanded or in a message, which
contradicts them. Hence $z \le C^*$, and $z \ge C^*$ since $z$ is the
cost of a path ([Lemma 9.2](#lem-9-2)(iv)). Note that the recorded goal
$t$ still has $g(t) = C^*$ at the end, since a later decrease of $g(t)$
would have pushed $t$ with $f < z$ and recorded it again at its expansion.
The chain of its parents is a path of cost at most $g(t)$
([Lemma 9.2](#lem-9-2)(iv)), i.e., an optimal one. If no goal is
reachable, $z$ stays $\infty$, and the search reports the instance
unsolvable. $\square$

We now bound the scope of the theorem by two remarks. First, on an
infinite (locally finite) graph with costs bounded below by some
$\varepsilon > 0$ the argument of termination also needs that each
worker makes progress. Indeed, a worker that is never scheduled (e.g.,
with more workers than cores) may hold the only node of an optimal path
while the others expand nodes above the optimum that no incumbent prunes
yet (on a finite graph the same situation costs time and memory rather
than termination, and the tester runs the 15-puzzle with at most one
worker per core for this reason, cf. `test_hda` in `test/test.cpp`).
Second, the theorem does not need a consistent heuristic, since the
reopenings are allowed; with a consistent one the sequential A$^*$ never
reopens (cf. [Chapter 2](02-best-first.md#ch-2)), while HDA$^*$ still
may, and its search overhead counts the reopenings. When the option
`node_limit` is set, a worker that stores that many nodes stops the
mailbox, all the workers return, and the status is `node_limit`, with no
claim on the cost.

## 9.3 The messages and the end of the search {#sec-9-3}

### Batches and flushes

Sending a node usually costs a push on a queue that another worker reads (on a
shared memory) or a message of MPI (e.g., over a network); in both cases a
rather large fixed part of the cost does not depend on the size of what is
sent, and the nodes for the
same worker are therefore (rather naturally) packed into one message, an
idea that HDA$^*$
takes from the transposition-table-driven scheduling of Romein et al.
(cf. section 3 of [[2]](Z-bibliography.md#ref-2)). In the library, each
worker has one outbox per other worker; a node is appended to the outbox
of its owner, an outbox that holds `batch` nodes leaves at once as one
message, and the outboxes that are not full leave after each block of
`flush_period`
expansions and before the worker waits ([Figure 9.4](#fig-9-4)). We remark
that the flush also matters for the search overhead: a node
of an optimal path that waits in an outbox for others to fill it is a
node that its owner cannot expand, while the owner expands worse nodes
and sends their successors around. Kishimoto et al. observed that packing
"too many" nodes in a message raises the expansions, while packing "few"
raises the cost of the communication (section 5.2 of
[[2]](Z-bibliography.md#ref-2)); hence, `batch` and `flush_period`
trade one overhead for the other.

![Figure 9.4. The outboxes of worker 0 with batches of 3 nodes. (a)
Two nodes wait for worker 1 and one for worker 2; (b) the third node for
worker 1 fills its outbox, which leaves as one message; (c) the flush,
after $k$ expansions or before the worker waits, sends the outbox of
worker 2, which is not full.](../figures/hash-distributed-batch.svg){#fig-9-4}

On a bounded family of channels (e.g., the mesh of single-producer
single-consumer queues used on FastFlow) the queue of a pair of workers
may be full. The part of a batch that does not fit then stays in the
outbox, and it is tried again when a further batch for that worker is
ready. A worker whose flush has left something behind keeps polling its inbox
(which empties the queues that the others are waiting on); thus, no
worker waits for another.

### The end on a shared memory

A worker that finds no node of $f < z$ in its OPEN list cannot stop,
since a node may be travelling to it, or another worker may still
generate one. In fact, the mailbox detects the end by counting, with the
deficit
counting of Dijkstra and Scholten [[35]](Z-bibliography.md#ref-35) on a
single counter $C$ in the shared memory instead of an acknowledgment of
each message. A worker is *active* from the start until it calls
`quiescent()` with its outboxes sent and no node of $f < z$, and it is
*passive* from then until it drains a message.

<a id="def-9-5"></a>**Definition 9.5** (The deficit counter). The
counter $C$ starts at $p$ (all the workers are active). A worker adds $k$
to $C$ before it pushes a batch of $k$ messages, and a passive worker that
drains messages adds 1 (it becomes active) before it processes the first
of them. A worker subtracts the $k$ messages that it has drained after it
has processed them, and a worker that becomes passive subtracts 1. A
passive worker only reads $C$, and it stops when it reads $C = 0$.

<a id="lem-9-6"></a>**Lemma 9.6** (Detection of the end on a shared
memory). Let each worker send only while active and become
active only by draining a message. Then (i) at any time $C$ equals the
number of the active workers plus the number of the messages that have
been counted by their sender and are still to be subtracted by their receiver;
(ii) $C = 0$ holds only when no worker is active and no message exists,
and then it holds forever; (iii) if at some time no worker is active and
no message exists, then $C = 0$ from then on, and all the workers stop.

*Proof.* (i) By induction on the operations of
[Definition 9.5](#def-9-5): each changes $C$ by the change of the
quantity on the right, since a message is counted before it is pushed
(hence before it can be drained) and subtracted after it is processed,
and the activity of a worker is added when it leaves the passive state
and subtracted when it enters it. (ii) Both terms of (i) are
nonnegative, so $C = 0$ forces both to be 0. From such a state no
operation raises $C$: an addition for a batch needs an active sender,
and an addition for an activity needs a message to drain. (iii) If no
worker is active and no message exists, (i) gives $C = 0$, and by (ii)
it stays 0; each passive worker reads it at its next call of
`quiescent()`. $\square$

These hypotheses hold for HDA$^*$. Indeed, a worker sends only during an
expansion (line 14) or a flush, both of which it makes while it has
work. A passive worker has no node of $f < z$, and since $z$ never
increases ([Lemma 9.2](#lem-9-2)(i)) no node of its OPEN list becomes
expandable again: only a message can give it work. Hence the end that
[Theorem 9.4](#thm-9-4) assumes is detected correctly: the workers stop
when, and only when, no worker has a node of $f < z$, none is expanding
and no node is in transit. A batch costs one atomic addition at the
sender and one at the receiver, and a passive worker only reads the
counter; [Figure 9.5](#fig-9-5) follows $C$ on a run of two workers.

![Figure 9.5. The counter $C$ on a run of two workers (the thick
segments are the active periods). Worker 0 pushes a batch of 3 nodes
($C$ from 2 to 5); worker 1 becomes passive (4), drains the batch (5)
and processes its 3 nodes (2); the two workers become passive in turn,
and the end is the moment at which $C$ reaches
0.](../figures/hash-distributed-counter.svg){#fig-9-5}

Counting cannot be dropped. A scheme sometimes used in "simple"
implementations is made of two rounds of idle flags without any count of
the messages: a worker sets its first flag when it has no work, sets its
second flag when it sees all the first flags set, and stops when it sees
all the second flags set. This scheme is not sound, since a flag carries
no information on
a message that was pushed after its owner looked at its inbox.

<a id="prop-9-7"></a>**Proposition 9.7** (Idle flags are not enough).
With two rounds of idle flags and no count of the
messages, two workers may both stop while a node is in the inbox of one
of them.

*Proof.* [Figure 9.6](#fig-9-6) gives such a run. Worker 1 has set its
first flag in an earlier round, while it had no work, and it looks at its
inbox and finds it empty (event 1). Worker 0 then expands its last node
and sends its successor $x$ to worker 1 (event 2); having no work left, it
sets its first flag, sees all the first flags set (that of
worker 1 is still set), sets its second flag, and goes on, since the
second flag of worker 1 is not set (event 3). Worker 1, acting on what
it saw at event 1, has no work: it sees all the first flags and then all
the second flags set, and stops (event 4); worker 0, at its next round,
sees the same and stops (event 5). Thus, the node $x$ is never processed; if
it is the only node of an optimal path at hand in the sense of
[Lemma 9.3](#lem-9-3), the search returns a suboptimal cost, or reports
an unsolvable instance. $\square$

![Figure 9.6. A run in which two rounds of idle flags stop both
workers while the node $x$ is still in the inbox of worker 1; the table
gives the flags after each event, those that change in the rust colour.
Counting the
messages, as in [Lemma 9.6](#lem-9-6) or in [Theorem 9.13](#thm-9-13),
excludes such a run.](../figures/hash-distributed-flags.svg){#fig-9-6}

## 9.4 The work distributions {#sec-9-4}

Much of the behavior of HDA$^*$ is decided by its distribution, which balances
two overheads against each other: a distribution that keeps the
successors of a node on the worker of the node lowers the communication,
but it tends to gather the nodes of least $f$ on a few workers, while
one that spreads the states uniformly balances the load at the price of
sending almost all the nodes. The *band effect* of Jinnai and Fukunaga
(section 3.1.1 of [[5]](Z-bibliography.md#ref-5)) explains the (somewhat
indirect) link between the load balance and the search overhead: if, with
two workers,
a worker owns a fraction $q$ of the nodes of each value of $f$, its
$n$-th expansion is roughly the $(n / q)$-th of A$^*$, and the farther
$q$ is from $1 / 2$ the more the workers expand out of the order of
A$^*$, i.e., nodes that A$^*$ would not reach. We now describe the six
distributions of `Distribution.h`, in the notation
$\mathrm{HDA}^*[\mathit{hash}, \mathit{abstraction}]$ of section 2.6 of
[[5]](Z-bibliography.md#ref-5) and by their names in the library;
[Figure 9.7](#fig-9-7) colours the 12 states of the $2 \times 2$ puzzle
by their owner under three of them.

![Figure 9.7. The 12 states of the $2 \times 2$ puzzle that are
reachable from the goal, which form a cycle, coloured by their owner
among $p = 3$ workers under Zobrist hashing (left), Abstract Zobrist
hashing with the projection of a tile onto its row (middle) and the
abstraction by the cell of tile 1 (right); the edges between two workers
are thick and in the rust colour. The random tables were drawn for the
figure, and on so small
a space the counts change from one draw to
another.](../figures/hash-distributed-owners.svg){#fig-9-7}

### From a hash to a worker

A distribution that computes a 64-bit hash $h$ maps it to a worker by
the multiply-high reduction of Lemire [[151]](Z-bibliography.md#ref-151),
$$\mathrm{reduce}(h, p) = \left\lfloor \frac{h \, p}{2^{64}}
\right\rfloor,$$
which uses the high bits of $h$ and costs one multiplication
(`Distribution.h:86-96`), and a 32-bit key $k$ is mapped by
`key_owner`, i.e., to $\lfloor k p / 2^{32} \rfloor$, the key standing
for the high 32 bits of the hash.

```cpp
--8<-- "snippets/hash-distributed-reduce.cpp"
```

<a id="prop-9-8"></a>**Proposition 9.8** (The reduction is balanced).
For $0 \le i < p \le 2^{64}$, the number of the values $h \in \{0,
\ldots, 2^{64} - 1\}$ with $\mathrm{reduce}(h, p) = i$ is either
$\lfloor 2^{64} / p \rfloor$ or $\lceil 2^{64} / p \rceil$; the same
holds for `key_owner` with $2^{32}$ in place of $2^{64}$.

*Proof.* $\mathrm{reduce}(h, p) = i$ if and only if $i \le h p / 2^{64}
< i + 1$, i.e., $h \in [\, i \, 2^{64} / p, (i + 1) \, 2^{64} / p)$, a
half-open interval of length $2^{64} / p$, which contains either
$\lfloor 2^{64} / p \rfloor$ or $\lceil 2^{64} / p \rceil$ integers. $\square$

Hence, for a "good" (uniform) hash, the share of the states of each worker
differs from $1 / p$ by less than $2^{-64}$ (or $2^{-32}$ for a key),
i.e., the reduction adds (almost) no imbalance; the modulo $h \bmod
p$, which the library uses only where the comparison of
[[5]](Z-bibliography.md#ref-5) requires it (below), has the same
property but costs a division.

### Zobrist hashing, HDA$^*$[Z]

With the state $s$ described by the features $x_0, \ldots, x_n$ that
hold in it and a table $R$ of random 64-bit strings, one per feature,
the Zobrist hash [[152]](Z-bibliography.md#ref-152) is
$$Z(s) = R[x_0] \oplus R[x_1] \oplus \cdots \oplus R[x_n],$$
with $\oplus$ the exclusive or of the bits (equation (1) of
[[5]](Z-bibliography.md#ref-5)); `ZobristDistribution`
(`Distribution.h:124`, for a `FeatureDomain`) draws $R$ from the seed,
takes as the key the high 32 bits of $Z(s)$ and as the owner that of the
key. This is the original distribution of HDA$^*$ (section 3 of
[[2]](Z-bibliography.md#ref-2)). When a state has exactly one feature
per position, as on the sliding-tile puzzles (one cell per tile),
Zobrist hashing is simple tabulation hashing, which is known to be
3-independent; with sets of features of varying size it is only pairwise
independent in general (e.g., $Z(\{x, y\}) = Z(\{x\}) \oplus Z(\{y\})$).
It spreads the states uniformly over the workers, regardless of the
structure of the graph, and the load is balanced; however, a successor
lands on the worker of its parent only with probability about $1 / p$,
and hence $\mathrm{CO} \approx 1 - 1 / p$ (section 2.4 of
[[5]](Z-bibliography.md#ref-5), where this is noted for an ideal uniform
assignment). Its hash can be updated by a move in the time of the
features that the move changes ([§9.5](#sec-9-5)).

### Abstract Zobrist hashing, HDA$^*$[Z, Afeature]

Jinnai and Fukunaga combine the hashing with a projection $A$ of the
features onto fewer *abstract* features (section 4 of
[[5]](Z-bibliography.md#ref-5), equation (2)):
$$\mathit{AZ}(s) = R'[A(x_0)] \oplus R'[A(x_1)] \oplus \cdots \oplus
R'[A(x_n)],$$
with a table $R'$ of random strings indexed by the abstract features. A
move that changes no abstract feature keeps the successor on the worker of
its parent, which lowers the communication. Also, the hashing of all the
abstract features keeps the load (rather) balanced; on the domains of
section 4.1 of [[5]](Z-bibliography.md#ref-5), the 15-puzzle and the 24-puzzle
among them, it was the best of the distributions compared there.
`AbstractZobristDistribution` (`Distribution.h:195`, for an
`AbstractFeatureDomain`) lays its table out by raw feature, i.e., it
stores $R'[A(x)]$ at the index of each raw feature $x$ (Algorithm 3 of
[[5]](Z-bibliography.md#ref-5)); thus, the hash costs as much
as $Z$ and is updated in the same way:

```cpp
--8<-- "snippets/hash-distributed-az-table.cpp"
```

On the sliding-tile puzzles one can take as the abstract feature of a tile
the half
of the board that holds it, as the library does, and then a move changes
it only when it slides a tile vertically across the middle line; on the $2
\times 2$ puzzle of
[Figure 9.7](#fig-9-7) (middle) the halves are the rows, a
horizontal move always keeps the owner, and at most half of the edges are cut.

We remark that a state in which two raw features have the same abstract
feature has the two entries cancelled in $\mathit{AZ}$, since $R'[a]
\oplus R'[a] = 0$. The hash is still a function of the state (which is
all that the correctness of [§9.2](#sec-9-2) requires), but its values
are less uniform. On the sliding-tile puzzles this never happens, since
each tile is in exactly one cell.

### Abstraction, HDA$^*$[Z, Astate]

`AbstractStateDistribution` (`Distribution.h:264`, for an
`AbstractStateDomain`) sends all the states with the same abstract state
to the same worker, hashing the abstract state with a seeded mix (a
64-bit mixing function of the abstract state, not a Zobrist hash of its
features). This is AHDA$^*$ of Burns et al.
[[3]](Z-bibliography.md#ref-3) with a uniform hash in place of their
perfect hash (hence the notation HDA$^*$[Z, Astate] rather than
HDA$^*$[P, Astate], cf. section 2.5 of [[5]](Z-bibliography.md#ref-5)).
Its communication is the least of these distributions, since a successor
changes worker only when it changes the abstract state, but no mechanism
balances the load: the abstract states are usually "few", and those that
contain the nodes of least $f$ may belong to few workers, which raises
the search overhead. On the 15-puzzle the abstraction is the triple of
the cells of tiles 1, 2 and 3, i.e., $16 \cdot 15 \cdot 14 = 3360$
abstract states.

### Perfect hashing, HDA$^*$[P]

`PerfectHashDistribution` (`Distribution.h:300`, for a
`PerfectHashDomain`) takes the rank of the state modulo $p$. This is the
distribution called "HDA$^*$" in [[3]](Z-bibliography.md#ref-3), and
Jinnai and Fukunaga found that on the sliding-tile puzzles it has both a
high search and a high communication overhead (sections 3.2 and 4.1 of
[[5]](Z-bibliography.md#ref-5)); the library keeps it to reproduce that
comparison. The reason, at least on the puzzles, lies in the following
proposition, which states for
the rank of the library what footnote 5 of section 3.1 of
[[5]](Z-bibliography.md#ref-5) observes for $p = 8$.

<a id="prop-9-9"></a>**Proposition 9.9** (The perfect hash modulo $p$).
Let the board have $n$ cells, let
$\mathrm{pos}(k)$ be the cell of tile $k$ (the blank being tile 0), and
let the rank be the one of `SlidingTile.h:441-454`, i.e.,
$$\mathrm{rank}(s) = \sum_{k=0}^{n-1} d_k \, (n - 1 - k)!, \qquad d_k =
|\{\, j > k : \mathrm{pos}(j) < \mathrm{pos}(k) \,\}|.$$
If $m$ is the least integer with $p \mid m!$, then $\mathrm{rank}(s)
\bmod p$ depends only on the relative order of the cells of the tiles
$n - m, \ldots, n - 1$.

*Proof.* The code computes $\mathrm{rank} \gets \mathrm{rank} \cdot (n -
k) + d_k$ for $k = 0, \ldots, n - 1$, which expands to the sum above.
Also, $d_k$ is the cell of tile $k$ minus the number of the cells below it
that the tiles $0, \ldots, k - 1$ occupy, i.e., the number of the tiles
$j > k$ below it. For $k < n - m$ one has $n - 1 - k \ge m$, hence $p$
divides $(n - 1 - k)!$ and the term vanishes modulo $p$; the remaining
$d_k$, $k \ge n - m$, compare the cells of the tiles $n - m, \ldots, n -
1$ only. $\square$

On the 15-puzzle with $p = 8$ one has $m = 4$, and the owner depends on
the relative order of the cells of the tiles 12, 13, 14 and 15 alone: a
move that changes this order may change the owner, any other move
keeps it, and the classes of states that share an owner follow neither
the load nor the locality of the search.

### The hash of the closed list, and a random owner

`StateHashDistribution` (`Distribution.h:328`, the default of both
searches) remixes the hash that any `Domain` gives to its closed list
with a seed; it behaves as Zobrist hashing for any `Domain`, with no
feature to define, but it cannot be updated by a move.
`RandomDistribution` (`Distribution.h:358`) sends each generated state
to a worker drawn at random, as in the randomized strategies of Kumar,
Ramesh and Rao and of Karp and Zhang (cf. section 6 of
[[2]](Z-bibliography.md#ref-2)). Hence, a state has no owner, a duplicate is
detected only when its copies happen to reach the same worker, and each
worker needs a separate generator (which the class declares by
`per_worker`). It
is the baseline without detection of the duplicates across the workers;
Kishimoto et al. report that on their planning benchmarks it ran out of
memory where HDA$^*$ did not, and that with more memory it expanded more
than an order of magnitude more nodes than HDA$^*$ (section 6 of
[[2]](Z-bibliography.md#ref-2)).

### A summary

The following table gathers what each distribution needs and what it
trades; the overheads are those that the scheme leads one to expect,
and only measurements (cf. [§9.8](#sec-9-8)) decide among them on a
given domain.

| Name | Needs | Incr. | Det. | CO | Load |
|:--|:--|:-:|:-:|:--|:--|
| Z | `FeatureDomain` | yes | yes | about $1 - 1/p$ | balanced |
| AZ | `AbstractFeatureDomain` | yes | yes | below Z | balanced |
| A | `AbstractStateDomain` | no | yes | lowest | not controlled |
| P | `PerfectHashDomain` | no | yes | high | uneven |
| H | `Domain` | no | yes | about $1 - 1/p$ | balanced |
| R | `Domain` | no | no | about $1 - 1/p$ | balanced |

The column "Incr." gives whether the distribution is an
`IncrementalDistribution`, "Det." whether the owner is a function of the
state, and the load of AZ is balanced only if the abstract features are
many with respect to $p$. All the distributions above are made "by hand": the
features, the projection or the abstraction come from the `Domain`. The
choice of the distribution from the problem (by partitioning a
graph of the features, by the operators, or dynamically) is the subject
of [Chapter 10](10-auto-distributions.md#ch-10).

## 9.5 Incremental Zobrist hashing {#sec-9-5}

A hash computed from scratch reads all the features of the state, i.e.,
$n - 1$ entries of the table on a puzzle of $n$ cells; Zobrist hashing can
instead be updated by a move, from the hash of the parent and the
features that the move changes (section 3 of
[[2]](Z-bibliography.md#ref-2) and section 2.4 of
[[5]](Z-bibliography.md#ref-5)). We now state the identity on which the
update rests, for the sets of features $F(s)$ and $F(t)$ of a state and
of its successor, whose symmetric difference $F(s) \triangle F(t)$ is
enumerated by `for_each_feature_change()`.

<a id="prop-9-10"></a>**Proposition 9.10** (Incremental Zobrist hash).
For any two states $s$ and $t$,
$$Z(t) = Z(s) \oplus \bigoplus_{x \in F(s) \triangle F(t)} R[x],$$
and the same holds for $\mathit{AZ}$ with the table laid out by raw
feature, and for the keys, i.e., for the high 32 bits of $Z$ and of
$\mathit{AZ}$.

*Proof.* The set $\{0, 1\}^{64}$ with $\oplus$ is an abelian group in
which each element is its own inverse. Hence
$$Z(s) \oplus Z(t) = \bigoplus_{x \in F(s)} R[x] \oplus \bigoplus_{x
\in F(t)} R[x] = \bigoplus_{x \in F(s) \triangle F(t)} R[x],$$
since each feature of $F(s) \cap F(t)$ appears twice and cancels, and
$\oplus$ with $Z(s)$ on both sides gives the claim. For $\mathit{AZ}$
the same computation holds with the entry $R'[A(x)]$ stored at the
index of $x$, which is how `AbstractZobristDistribution` computes both
the hash and its update. Finally, taking the high 32 bits, i.e., $y
\mapsto \lfloor y / 2^{32} \rfloor$, commutes with $\oplus$, since the
exclusive or acts bit by bit; hence the identity holds for the keys as
well. $\square$

<a id="cor-9-11"></a>**Corollary 9.11** (Owners from the parent). For a
distribution that is an `IncrementalDistribution`, the function that
`successor_owners()` returns maps each successor $t$ of a node of state
$s$ and key $k$ to `( owner( t , p ) , key( t ) )` (the key of a state
for any distribution is `distribution_key()`, which is 0 when the
distribution is not an `IncrementalDistribution`).

*Proof.* `key( k , s , t )` applies the identity of
[Proposition 9.10](#prop-9-10) to the key $k$ of $s$ (the excerpt
below), and the owner is `key_owner()` of the key, as for
`owner( t , p )`. $\square$

Our tester checks the corollary on random walks of the 8-, 15- and
24-puzzles, for $p \in \{1, 2, 3, 7, 64\}$, against the owners and the
keys computed from scratch (`check_successor_owners` in
`test/test.cpp`).

```cpp
--8<-- "snippets/hash-distributed-zobrist-key.cpp"
```

On the sliding-tile puzzles a move slides one tile $k$ from a cell $c$
to the cell $b$ of the blank, and the symmetric difference is $\{(k, c),
(k, b)\}$: two entries of the table instead of 15 on the 15-puzzle (and
none for the state being expanded, whose key is stored in its node).
Both cells are found without locating the blanks: the packed states
$s$ and $t$ differ only in the fields of the cells $b$ and $c$, both
of which hold $k$ in $s \oplus t$, and therefore two counts of trailing zeros
of $s \oplus t$ give them (`SlidingTile.h:335-346`).
[Figure 9.8](#fig-9-8) shows the update for a move inside a half of the
board and for one across the middle line.

![Figure 9.8. The updates of the Zobrist hash and of the abstract
Zobrist hash by two moves from the state $s$ (middle) of the 15-puzzle.
Sliding 7 to the right (left) changes the features $(7, 4)$ and $(7, 5)$,
in the rust colour, but not the abstract feature (7, top half), so
that $t_1$ stays on the worker of $s$ under AZ; sliding 3 up across the
middle line (right) also changes its abstract feature, and $t_2$ may
change worker.](../figures/hash-distributed-zobrist.svg){#fig-9-8}

## 9.6 HDA$^*$ over MPI {#sec-9-6}

`HDAStarMPI` (`HDAStarMPI.h`) is the search of [§9.2](#sec-9-2) with one
MPI process per worker, as in the original implementation of Kishimoto
et al., where the workers are MPI processes even on a single multicore
machine (section 3 of [[2]](Z-bibliography.md#ref-2)). Since the processes
share no memory, a node records the *state* of its parent rather
than a pointer (`HDAStarMPI.h:158-168`), the incumbent travels as a
message to all the other processes, and the state and the cost types must be
trivially copyable, since the nodes travel as bytes
(`HDAStarMPI.h:114`). MPI is initialized by the caller, all the
processes of the communicator call `solve()`, and all of them return the
same `Solution`.

**Algorithm 9.2** (HDA$^*$ over MPI, on the process $i$).

```text
 1  procedure HDA*-MPI(s₀) on the process i
 2     z ← ∞;  if owner(s₀) = i then INSERT(s₀, 0, nil)
 3     loop
 4        for each message m received do
 5           if m is a batch of nodes then INSERT each of them
 6           else if m is an incumbent z' then z ← min(z, z')
 7           else if m is a stop then stop expanding
 8        if OPEN_i has a node of f < z then
 9           expand it as at lines 7 to 15 of HDA*, a goal that lowers z
               sending z to all the other processes
10        else send all the outboxes
11        if OPEN_i has no node of f < z and the mailbox detects the end
12           then break
13     the sends in progress complete
14     the process of the least goal broadcasts it, and the path is rebuilt
         from it, each state broadcast with its parent by its owner
```

In the library, the batches are sent with nonblocking sends (`MPI_Isend`),
whose buffers
are kept until MPI has completed them, and received with `MPI_Iprobe`
and `MPI_Recv` from any process; the mailbox duplicates the
communicator, and its messages never meet those of the program. The
messages are a variant of a node and an incumbent
(`HDAStarMPI.h:181-187`), and the reception sorts them, then hashes and
prefetches all the nodes before inserting them:

```cpp
--8<-- "snippets/hash-distributed-mpi-receive.cpp"
```

We let each process keep an incumbent $z_i$, the least cost that it has
found or received; since $z_i$ is never smaller than the least cost $z$
found by any process, a process prunes less than a shared $z$ would,
which may cost expansions (how many depends on the latency of the
incumbent messages) but not correctness.

### The termination detection of Mattern

Over MPI, the end is detected by the four counters of Mattern
[[11]](Z-bibliography.md#ref-11), in the form described in section 3 of
[[2]](Z-bibliography.md#ref-2) and in section 3.3.1 of
[[1]](Z-bibliography.md#ref-1). Kishimoto et al. used the time
algorithm of Mattern, a variation that needs one wave only, while the
library uses the basic method with two waves. Each process $i$ counts
the messages that it has sent, $s_i$, and those that it has received,
$r_i$, both nondecreasing. We say that a process is *passive* when it has
no node of
$f < z_i$ and has posted all its outboxes, and it becomes active only by
receiving a message (an incumbent never gives work, and a node usually
may). Note that a process with work never contributes to a wave.

<a id="def-9-12"></a>**Definition 9.12** (Waves). A *wave* is a
nonblocking reduction (`MPI_Iallreduce`) of the pairs $(s_i, r_i)$ of
all the processes into their sums $(S, R)$; a process contributes to a
wave only when it is passive, i.e., from `quiescent()`, at an instant
$t_i$, and it keeps receiving and working while the wave completes. A
process starts the wave $k$ only after the wave $k - 1$ has completed
for it. The end is declared at the first wave $k \ge 2$ whose sum of the
sends $S_k$ equals the sum of the receipts $R_{k-1}$ of the wave before.

Since a reduction completes for a process only after all the processes
have contributed to it, each contribution to the wave $k$ comes after all the
contributions to the wave $k - 1$: if $t_i^{k}$ is the instant of the
contribution of $i$ to the wave $k$, then $t_j^{k} > \max_i t_i^{k-1}$
for all $j$. A single wave is not enough, as the left half of [Figure
9.9](#fig-9-9) shows. In fact, the counters of a wave are read at
different instants. A message sent before the cut and received
after it is counted by $S$ but not by $R$, while one sent after the cut
and received before it is counted by $R$ only. The two errors may cancel
while a message is in flight.

![Figure 9.9. (a) One wave, with the contributions $(s_i, r_i)$ of
the three processes on a dashed cut: $m_1$ crosses the cut forwards and
$m_2$ backwards, so that $S = R = 1$ although $m_3$ is in flight and
will make $P_1$ active. (b) Two waves at the end of a computation: $T$
is the last contribution to the wave $k - 1$, and the equality of $S_k$
and $R_{k-1}$ forces each inequality of the chain to be an
equality.](../figures/hash-distributed-waves.svg){#fig-9-9}

<a id="thm-9-13"></a>**Theorem 9.13** (Mattern's four counters detect the
end exactly). Let each process send only while active and become
active only by receiving, and let $T = \max_i t_i^{k-1}$. (i) If $S_k =
R_{k-1}$, then at the instant $T$ all the processes were passive and no
message was in transit; hence the computation had ended at $T$, and the
end is never declared before it happens. (ii) If the computation ends at
some instant $T^*$ and each process keeps calling `quiescent()`, then
the end is declared, at the latest at the second wave all of whose
contributions follow $T^*$, and all the processes declare it at the same wave.

*Proof.* (i) Let $S(t) = \sum_i s_i(t)$ and $R(t) = \sum_i r_i(t)$. A
message is received only after it is sent; hence, $R(t) \le S(t)$ at all
times, and the counters are nondecreasing. Since $t_i^{k-1} \le T$
for all $i$, one has $R_{k-1} = \sum_i r_i(t_i^{k-1}) \le R(T)$, and
since $t_i^{k} > T$ for all $i$, one has $S(T) \le \sum_i s_i(t_i^{k})
= S_k$. Therefore
$$R_{k-1} \le R(T) \le S(T) \le S_k = R_{k-1},$$
and all are equalities. From $R(T) = R_{k-1}$ and the monotonicity of
each $r_i$ it follows that $r_i(T) = r_i(t_i^{k-1})$ for each $i$, i.e.,
no process received a message between its contribution to the wave $k - 1$
and $T$. Since it was passive at its contribution and it can become active
only by receiving, it was passive at $T$. From $S(T) = R(T)$,
no message was in transit at $T$. A computation in which all the processes
are passive and no message is in transit has ended, and it stays so,
since nobody can send.

(ii) After $T^*$ no message is sent or received; thus, $S$ and $R$
are constant and equal, say to $N$. Every process is passive and calls
`quiescent()`, and therefore each wave completes in finite time; let $W$ be
the first wave all of whose contributions follow $T^*$ (it exists, since
the contributions to a wave follow those to the wave before). Then $R_W
= N$ and $S_{W+1} = N$, and the end is declared at the wave $W + 1$ at
the latest. The sums of a reduction are the same on all the processes, and
so is the decision. $\square$

The end that [Theorem 9.13](#thm-9-13) detects is the end assumed by
[Theorem 9.4](#thm-9-4). At that instant all the incumbents sent have been
received; hence, each $z_i$ equals the least cost found by any
process; no process has a node of $f < z$, none is expanding, and no
node is in transit. Hence, the proof of [Theorem 9.4](#thm-9-4) applies
word for word, with $z_i$ in the pruning of the process $i$ (which, as
remarked above, prunes no more than $z$ would). One can expect a wave to
cost $O(\log p)$ steps of latency on $p$ processes, and it runs while the
processes
keep working; a process never waits for a wave.

A process that reaches the node limit sends a note (an empty batch) to
each other process; a process that sees it stops expanding and
contributes to the waves as stopped, and the computation ends at the
first wave to which all have contributed stopped. At the end, the
messages still in flight are received once, by a collective exchange of
the number of the messages that each process has sent to each other one
(`MPI_Reduce_scatter_block`), and no send is left pending when the
mailbox is destroyed.

### The path

We rebuild the path after the end. First, the counters of the processes
are gathered (`MPI_Allgather`), and the least goal and its process are
found by a reduction with `MPI_MINLOC` (`HDAStarMPI.h:380-385`). Then each
state of the path is broadcast with its parent by its owner, starting
from the goal ([Figure 9.10](#fig-9-10)); each process knows the next
parent, hence its owner, and the owner looks it up in its closed list.
This is why the distribution of `HDAStarMPI` must be deterministic: with
the random owner, the process queried for a state need not be the one that
stores it. By [Lemma 9.2](#lem-9-2) the chain of the parents is a path
of cost at most that of the goal, and the rebuilding takes one broadcast
per state of the path.

![Figure 9.10. The path rebuilt over MPI: the process of the least
goal broadcasts $(t, b)$, i.e., the goal and its parent; the owner of
$b$ broadcasts $(b, a)$; and so on to $s_0$, whose parent is nil.
The colour of a state is its
owner.](../figures/hash-distributed-path.svg){#fig-9-10}

```cpp
--8<-- "snippets/hash-distributed-mpi-path.cpp"
```

## 9.7 The engineering of the implementation {#sec-9-7}

The algorithm of [§9.2](#sec-9-2) leaves much freedom to the
implementation, and we use it as follows.

* The inbox is drained as a whole, and the hashes of all its nodes are
  computed and their slots in the closed list prefetched before the
  first lookup, so that the cache misses of the lookups overlap instead
  of following one another (the excerpt below; the same is done over
  MPI). A message from another worker typically carries a state that the
  receiver has not touched recently, and its lookup is then likely to be a
  "cold" miss.
* The incumbent is an atomic variable on a separate cache line, and
  each worker is aligned to a cache line (`HDAStar.h:250` and `:260`),
  and the counters and the lists of one worker do not share a line
  with those of another. Also, the pruning tests read the incumbent with a
  relaxed load: a stale value is larger than the current one (the
  incumbent only decreases), and therefore it prunes less, never more.
* Each worker allocates its nodes from a private pool and keeps its
  closed list in a private table (cf.
  [Chapter 3](03-data-structures.md#ch-3)); the memory of a node is
  therefore written by one thread only, which is the content of
  [Lemma 9.1](#lem-9-1).
* The key of a node, 32 bits, is stored in the node, where it fits in
  the padding of a node whose state is one or two words
  (`Node.h:69-83`); the key of the state being expanded is never recomputed.
* With integer costs the OPEN list is a list of buckets, as in the
  sequential search; otherwise it is a binary heap.
* A batch costs one atomic addition at the sender and one at the
  receiver on the counter of [Definition 9.5](#def-9-5), and a passive
  worker only reads that counter.

```cpp
--8<-- "snippets/hash-distributed-batch.cpp"
```

The expansion, with the owners of the successors computed from the key
of the parent, the pruning by the incumbent and the choice between a
local insertion and a message, is the following excerpt.

```cpp
--8<-- "snippets/hash-distributed-expand.cpp"
```

## 9.8 Benchmark {#sec-9-8}

**Table 9.1.** HDA$^*$ with the distributions Z, AZ and P on the
threads, and with Z and AZ on OpenMP, on FastFlow and over MPI, on the
15-puzzle.

--8<-- "tables/hash-distributed-tiles.md"

The experiment runs on 12 of the 100 instances of Korf
[[7]](Z-bibliography.md#ref-7) (the instances 1, 3, 10, 14, 15, 22, 26,
27, 32, 33, 37 and 52, as listed in `tools/bench/campaign.sh`), with the
Manhattan distance, the default batches of 64 nodes and flushes every 64
expansions, on 1, 2, 4, 8, 16, 32, 48, 64 and 96 cores. A row is a
search, i.e., HDA$^*$ with its distribution in brackets and its runtime
(thread, openmp, fastflow or mpi), or the sequential A$^*$ with the OPEN
list on buckets, against which the speedups are computed; a column is
the number of cores given to the search (one worker per core, the
workers pinned, and over MPI one process per core). A cell gives the
seconds summed over the instances, each instance counting with the
median of its two runs, and in parentheses the speedup over the
sequential A$^*$ on 1 core; a cell is empty when the search did not solve
all the instances. The machine, the compiler and the A/A noise (the
gap between the two runs of the same search, which are made in two
processes) are stated in the note under the table. All costs are
integers, and the cost found by each run is checked against that of the
sequential search. We measure the wall-clock time of `solve()`, i.e., from
the creation of the workers to the solution, the tracing of the path
included; since the benchmark gives no executor to the search, the start
of the threads is included as well. The time does not include the
start of the process, the reading of the instance, the construction of
the `Domain` and of the tables of the distribution, nor the release of
the memory of the search. We remark that the table is a placeholder until
the campaign is run on an idle machine.

## 9.9 Notes and references {#sec-9-9}

The ownership of a state by a hash comes from PRA$^*$ of Evett et al.,
and the asynchronous exchange of the nodes, packed per destination, from
the transposition-table-driven scheduling of Romein et al. HDA$^*$ of
Kishimoto, Fukunaga and Botea [[2]](Z-bibliography.md#ref-2) combined
the two with Zobrist hashing and showed that the scheme scales on
clusters and on multicore machines (we cite PRA$^*$ and TDS through
section 2.3 of [[5]](Z-bibliography.md#ref-5) and section 3 of
[[2]](Z-bibliography.md#ref-2), where they are summarized). Also, the survey
of Fukunaga et al. [[1]](Z-bibliography.md#ref-1) places HDA$^*$ among
the decentralized parallel searches (section 4.1) and describes the
termination detection of Mattern [[11]](Z-bibliography.md#ref-11)
(section 3.3.1); the deficit counting is that of Dijkstra and Scholten
[[35]](Z-bibliography.md#ref-35). Burns et al.
[[3]](Z-bibliography.md#ref-3) compared HDA$^*$ with the perfect hash
and AHDA$^*$ with SafePBNF (cf. [Chapter 11](11-safepbnf.md#ch-11)).
Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) showed that the
perfect hash was a poor representative of HDA$^*$, introduced the
notation $\mathrm{HDA}^*[\mathit{hash}, \mathit{abstraction}]$, Abstract
Zobrist hashing and the automatic distributions of
[Chapter 10](10-auto-distributions.md#ch-10), and explained the search
overhead of HDA$^*$ by the band and burst effects. In their comparison
on the 15-puzzle and the 24-puzzle, HDA$^*$[Z] was faster than SafePBNF
on the instances that take more than about 10 seconds, and SafePBNF on
the easier ones, where the initial burst weighs more (section 3.2 of
[[5]](Z-bibliography.md#ref-5)); the outcome on other domains and machines
is not necessarily the same, and on grid pathfinding the abstraction-based
searches were faster than HDA$^*$[Z], whose communication overhead is high
there (sections 3.3 and 3.4 of [[5]](Z-bibliography.md#ref-5)). Finally,
the Zobrist hash comes from the
game-playing programs [[152]](Z-bibliography.md#ref-152), and the
reduction of a hash to an interval without a division from Lemire
[[151]](Z-bibliography.md#ref-151).

As to where the library stands, the following facts are stated from the
code and the papers, not from measurements. The library offers all the
distributions of the comparison of Jinnai and Fukunaga that need no
automatic construction (Z, AZ, A and P), together with the hash of the
closed list and the random owner, on four runtimes and under the same
search. The original implementation of [[2]](Z-bibliography.md#ref-2)
ran on MPI only, and that of [[5]](Z-bibliography.md#ref-5) used
Pthreads with a try-lock on each message queue (section 4.1 of
[[5]](Z-bibliography.md#ref-5)), while the shared-memory search of the
library exchanges batches through a mailbox whose end detection costs
one atomic addition per batch at each side. Furthermore, the incremental
key, the
prefetching of whole inboxes and the skipping of the inverse move are
engineering choices that these papers do not discuss (the last one is
in [[4]](Z-bibliography.md#ref-4), for the sequential search). On the
other hand, the operator-based Zobrist hashing (OZHDA$^*$) and the
dynamic abstraction (DAHDA$^*$) of Jinnai and Fukunaga, and the
automatic feature abstractions of section 4.2 of
[[5]](Z-bibliography.md#ref-5), are not in the library yet (cf.
[Chapter 10](10-auto-distributions.md#ch-10)), and neither is the time
algorithm of Mattern, which detects the end with one wave instead of two
at the price of a time stamp on each message.

## 9.10 Exercises {#sec-9-10}

**Exercise 9.1** (`ctest -R PBeFS.hda-thread`). Prove that the key of
`ZobristDistribution` computed by `key( k , s , t )` equals `key( t )` for
each successor $t$ of $s$, using [Proposition 9.10](#prop-9-10),
and that the same holds for `AbstractZobristDistribution`; the test
checks it on random walks of the 8-, 15- and 24-puzzles for five values of $p$.

**Exercise 9.2** (`ctest -R PBeFS.hda-thread`). On the 15-puzzle with
the abstract features of `SlidingTile.h` (a tile and the half of the
board that holds it), count, over the 16 cells of the blank, the pairs
(cell of the blank, move) and those whose move changes an abstract
feature. Show that the latter are 8 of 48, and discuss why this is not
the fraction of the generated nodes that change worker under AZ (hint:
the cells of the blank are not equally frequent among the generated
states, and a node whose abstract features change still keeps its owner
with probability about $1 / p$).

**Exercise 9.3** (`ctest -R PBeFS.hda-thread`). Show that HDA$^*$ with
`RandomDistribution` still returns an optimal path (which parts of the
proof of [Theorem 9.4](#thm-9-4) need the deterministic owner, and which
do not?), and that it stores at most $p$ nodes of each state. Its test
runs it on 21 states of the 8-puzzle, whose $9!/2$ reachable states
bound the memory whatever the timing.

**Exercise 9.4** (`ctest -R PBeFS.hda-thread`). Show that the proofs of
[Theorem 9.4](#thm-9-4) and of [Lemma 9.6](#lem-9-6) do not depend on
the values of `batch` and `flush_period`, and that with both equal to 1
the search sends each node as soon as it is generated. The test runs the
search with batches of 1 and with the parameters set to 0 (which
`normalize()` turns into 1).

**Exercise 9.5** (`ctest -R PBeFS.hda-thread`). Build a graph of four
states with a consistent heuristic, a distribution on two workers and an
order of the messages in which HDA$^*$ expands a state twice. Which
counter of the `SearchStats` records the second expansion?

**Exercise 9.6** (`ctest -R PBeFS.mpi`). Generalize the run of
[Figure 9.9](#fig-9-9)(a): show that for every $p \ge 3$ there are runs
in which a single wave has $S = R$ while a message is in transit, and
find the least number of messages that such a run needs. Does such a
run exist with $p = 2$ processes?

**Exercise 9.7** (`ctest -R PBeFS.mpi`). Bound the number of waves that
complete after the end $T^*$ of the computation before the end is
declared ([Theorem 9.13](#thm-9-13)(ii)), counting the wave that may be
in flight at $T^*$, and show that the bound is attained.

**Exercise 9.8** (`ctest -R PBeFS.mpi`). Show that two processes that
build `StateHashDistribution` with the same seed assign each state to
the same owner, and describe how the rebuilding of the path of
[§9.6](#sec-9-6) would fail with `RandomDistribution`. The test runs the
search over MPI with each deterministic distribution, with 1, 2 and 4
processes.

**Exercise 9.9** (`ctest -R PBeFS.mpi`). Show that, after a process
reaches the node limit, all the processes return from `solve()` with the
status `node_limit`, even if some process found a goal in the meantime,
and that no message is left pending (the test stops the search over MPI
at 2000 nodes).

**Exercise 9.10** (`ctest -R PBeFS.hda-thread`). Using
[Proposition 9.9](#prop-9-9), find the tiles on which the owner of a
state of the 15-puzzle depends under `PerfectHashDistribution` for $p =
3$, $p = 5$, $p = 7$ and $p = 16$, and argue for which of these values
the distribution is closest to a uniform one.

**Exercise 9.11** (`ctest -R PBeFS.hda-thread`). Turn the two rounds of
idle flags of [Proposition 9.7](#prop-9-7) into a sound detection by
adding counters of the messages, and show that the result is the counter
of [Definition 9.5](#def-9-5) or the waves of
[Definition 9.12](#def-9-12), depending on whether the workers share memory.

**Exercise 9.12** (discussion, no test). Which problems are likely to
parallelize well with HDA$^*$, and which badly? Consider the cost of an
expansion against that of a message, the number of states with $f <
C^*$ against $p$, the branching factor, the number of states with the
same $f$, and whether a good abstraction of the features exists; then
place the domains of [Chapter 14](14-domains.md#ch-14) in your
classification.
