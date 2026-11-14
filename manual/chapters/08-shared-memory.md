# 8. Parallel A* in shared memory {#ch-8}

[Source: `include/PBeFS/astar/SPAStar.h` (`SPAStar`),
`include/PBeFS/astar/MQAStar.h` (`MQAStar`),
`include/PBeFS/astar/FarmAStar.h` (`FarmAStar`),
`include/PBeFS/astar/PASE.h` (`PASE`, `PASEOptions`),
`include/PBeFS/astar/ConcurrentClosed.h` and its three implementations
in `include/PBeFS/astar/structures/`, `include/PBeFS/core/SearchUtils.h`
(`ParallelParams`, `run_workers`), `include/PBeFS/astar/Solution.h`
(`SearchStats`), `include/PBeFS/astar/SlowDomain.h`]

We are concerned in this chapter with the parallel searches of the
library in which the $p$ workers share one OPEN list, or one owner of
it, i.e., Simple Parallel A$^*$ (SPA$^*$), MultiQueue A$^*$ (MQA$^*$),
the master-worker farm and Parallel A$^*$ for Slow Expansions
(PA$^*$SE). Since this is the first chapter on parallel search, we also
state in it by which overheads the parallel searches of the library are
measured ([§8.2](08-shared-memory.md#sec-8-2)) and what they must
guarantee, i.e., an optimal path with an admissible heuristic and a
sound decision of when to stop ([§8.3](08-shared-memory.md#sec-8-3)).
The decentralized searches, Hash Distributed A$^*$ and Safe PBNF, are
the subject of [Chapter 9](09-hash-distributed.md#ch-9) and [Chapter
11](11-safepbnf.md#ch-11), and the last sections of this chapter compare
all of them. We use the notation of [Chapter
1](01-introduction.md#ch-1), and we refer to [Chapter
2](02-best-first.md#ch-2) for the sequential A$^*$ and the lemmas on
which the proofs below rest.

## 8.1 The interface of the parallel searches {#sec-8-1}

Each search of this chapter is a class template over the `Domain` $D$ of
the problem (see [Chapter 1](01-introduction.md#ch-1)) and over the
executor that runs its workers; SPA$^*$ and the farm also take the OPEN
list as a template parameter, which is the bucket queue when the costs
are integers and the binary heap otherwise ([Chapter
3](03-data-structures.md#ch-3)). In turn, the executor is one of the
three runtimes of the library (a pool of C++ threads, the default, a
team of OpenMP and the threads of FastFlow
[[12]](Z-bibliography.md#ref-12)), and the same code of a search serves
each of them (the tests run SPA$^*$, MQA$^*$ and PA$^*$SE on the C++
threads, SPA$^*$ and MQA$^*$ also on OpenMP, and the farm on all three).
A search is built from the `Domain` and from a `ParallelOptions`, and
`solve()` returns a `Solution`, as the sequential A$^*$ does:

```cpp
SPAStar< FifteenPuzzle > spa( domain , options );
const auto sol = spa.solve();   // sol.status, sol.cost, sol.path
```

`ParallelOptions` joins the options of the sequential search (the node
limit, among others) with the `ParallelParams` of
`include/PBeFS/core/SearchUtils.h:94`: the number of workers `threads`
($p$, 1 by default), the nodes per message `batch` (64), the expansions
between two flushes of an outbox `flush_period` (64) and the `seed` of
the random choices (1); `normalize()` sets to 1 any of the first three
that is 0 (`SearchUtils.h:102`). Of these, the searches of this chapter
read `threads`, the farm also `batch` (the nodes per task) and MQA$^*$
the `seed` of its draws; `PASEOptions` adds the inflation
$\varepsilon_p \ge 1$ of PA$^*$SE, `eps`, which is raised to 1 when
smaller (`PASE.h:126` and `PASE.h:174`). PA$^*$SE needs a
`PairwiseDomain` (`Domain.h:198`), i.e., a `Domain` with a function
`pairwise_heuristic(s, t)` that estimates the cost of a path from $s$ to
$t$, whose requirements are stated in
[§8.8](08-shared-memory.md#sec-8-8).

One either gives the executor with `set_executor()`, and it is then
reused by all the searches, or lets `solve()` build one with as many
threads as workers (one more for the farm, whose master runs on a
separate thread). In both cases `run_workers()` (`SearchUtils.h:121`)
calls the body of worker $i$ on a separate thread, pinned to the core
$i$ of the affinity mask when the environment variable `PBEFS_PIN` is
set to a value other than 0 (`SearchUtils.h:77`; the thread of the
caller, when an executor runs a worker on it, is not pinned). As in the
sequential search, a `Solution` holds the status (`solved`, `unsolvable`
or `node_limit`), the cost and the path, which is rebuilt from the goal
through the parent pointers after all the workers have returned, and the
counters of `SearchStats` (`Solution.h:64`): the expansions
(re-expansions included), the re-expansions, the successors generated,
the duplicates, the nodes sent to another worker, the batches sent, the
nodes stored, the seconds of the wall clock, and the expansions and the
stored nodes of each worker. In a parallel search the node limit bounds
the nodes stored by each worker, and the search that reaches it returns
`Status::node_limit` with no path.

## 8.2 The overheads of a parallel search {#sec-8-2}

A parallel search with $p$ workers is rarely $p$ times faster than the
sequential one, and the gap is usually split, following the survey of
Fukunaga et al. ([[1]](Z-bibliography.md#ref-1), section 3.1) and the
analysis of Jinnai and Fukunaga ([[5]](Z-bibliography.md#ref-5),
section 2.1), into four overheads. We give them in the form in which
the counters of the library compute them; $E_i$ denotes the expansions
of the worker $i$ (the entry $i$ of `expanded_by`, re-expansions
included) and $E^{\mathrm{seq}}$ those of the sequential A$^*$ on the
same instance, with the same heuristic and the same tie-breaking.

<a id="def-8-1"></a>**Definition 8.1** (Search overhead). The search
overhead of a parallel run is

$$
\mathrm{SO} = \frac{\sum_{i=0}^{p-1} E_i}{E^{\mathrm{seq}}} - 1 .
\tag{8.1}
$$

It counts the expansions that the sequential search does not make, i.e.,
the nodes with $f > C^*$, the extra nodes of the last layer $f = C^*$
and the re-expansions, which a parallel search pays even with a
consistent heuristic, since a worker may expand a node before another
worker has found the cheaper path to it ([Figure
8.3](08-shared-memory.md#fig-8-3)). A negative $\mathrm{SO}$ is possible
(and gives a superlinear speedup), since a parallel search may reach a
goal of the last layer after fewer of its nodes than the sequential one.
Since it needs a second, sequential run, the library does not compute
$\mathrm{SO}$; one obtains it from the field `expanded` of the two
`SearchStats`, typically from a run of `AStar` on the same instance.

<a id="def-8-2"></a>**Definition 8.2** (Communication overhead). The
communication overhead is the fraction of the generated nodes that
change worker,

$$
\mathrm{CO} = \frac{\mathrm{sent}}{\mathrm{generated}} ,
\tag{8.2}
$$

and 0 when no node is generated, as `communication_overhead()` computes
it (`src/Solution.cpp:51`). Each such node costs a copy, a message or a
lock. Note that the numerator counts nodes, as opposed to messages, as in
section 2.1 of [[5]](Z-bibliography.md#ref-5); the number of the
messages, i.e., of the batches, is the separate counter `messages`.

<a id="def-8-3"></a>**Definition 8.3** (Load balance). The load balance
is the ratio of the busiest worker to the mean,

$$
\mathrm{LB} = \frac{\max_i E_i}{\frac{1}{p} \sum_{i=0}^{p-1} E_i}
\ \ge 1 ,
\tag{8.3}
$$

1 being a perfect balance; `load_balance()` (`src/Solution.cpp:36`)
returns 1 when no worker has expanded a node.

<a id="def-8-4"></a>**Definition 8.4** (Coordination overhead). The
coordination (or synchronization) overhead is the fraction of the time
of the workers spent waiting, for a lock, for a message or for the
other workers, instead of expanding.

No counter of the library measures the coordination overhead, which
depends on the machine (and, in shared memory, also on the contention
for the memory bus, cf. section 2.1 of [[5]](Z-bibliography.md#ref-5)).
It shows in the wall clock and in a profiler. Note also that the
shared-memory searches of this chapter report $\mathrm{sent} = 0$, and
hence $\mathrm{CO} = 0$ by (8.2), although each node of SPA$^*$ and of
MQA$^*$ goes through a shared structure and each node of the farm
travels to a worker and back. In fact, (8.2) measures the nodes that
change owner in a distributed search ([Chapter
9](09-hash-distributed.md#ch-9)), and the searches of this chapter have
no owners. For the same reason the farm reports the expansions of the
master only, and its $\mathrm{LB}$ is 1 by construction.

These overheads are not independent ([[5]](Z-bibliography.md#ref-5),
section 2.1). A distribution of the states that keeps the successors of
a node on its worker lowers $\mathrm{CO}$, but it tends to gather the
nodes of least $f$ on a few workers. This raises $\mathrm{LB}$ and,
through it, $\mathrm{SO}$, since the other workers then expand nodes
that the sequential search would not reach. This is the band effect of
section 3.1.1 of [[5]](Z-bibliography.md#ref-5), which is still visible
when the expansions are made artificially slow, and is therefore usually
attributed to the load imbalance. Conversely, a uniform random
distribution balances the load at the price of sending a fraction $1 -
1/p$ of the nodes. The centralized searches of this chapter sit at one
end of this trade-off, since they have no owners, hence no $\mathrm{CO}$
in the sense of (8.2), and a balanced load, and they pay with
coordination, since each expansion goes through a shared structure or
through a master.

## 8.3 Correctness and termination of a parallel search {#sec-8-3}

The sequential A$^*$ stops at the first goal it selects, since by
[Lemma 2.8](02-best-first.md#lem-2-8) every node it selects has
$f \le C^*$. A parallel search cannot, since, when a worker selects a goal of
cost $g$, another worker may hold a node of smaller $f$, or be
expanding one, or such a node may be travelling towards the worker
that will expand it ([[1]](Z-bibliography.md#ref-1), sections 3.1 and
3.3.1). Hence, the parallel searches of the library share the cost $z$
of the best goal found so far (the *incumbent*, $z = \infty$ before the
first goal), prune each node with $f \ge z$, test the goals when a node
is selected for expansion, and stop only when no node that could
improve $z$ is left in any worker or in transit, a condition that we
now make precise.

An *entry* is a pair $(u, \gamma)$ of a node and a value of its $g$, put
in OPEN when $g(u)$ takes the value $\gamma$. The entry is *stale* when
$u$ has been closed or $g(u) < \gamma$ since then (the lazy deletion of
[Chapter 3](03-data-structures.md#ch-3), `OpenList.h:80`), and its key
is $f = \gamma + h(u)$.

<a id="def-8-5"></a>**Definition 8.5** (Entries in flight). At an
instant $\tau$ of a parallel search, an entry is *in flight* if it is
(i) in an OPEN list, (ii) held by a worker that has selected it and has
not yet finished processing it, the successors it generates being part
of what it holds until they are pushed, or (iii) in a message (or task)
that has been sent and is still to be received.

We say that a parallel search is *well-formed* if the following three
conditions hold. (a) Each entry is in flight from its creation until it
is selected. (b) A worker that selects a non-stale entry $(u, g(u))$
with $f < z$ expands $u$, i.e., for each successor $v$ it either prunes
$v$, because $g(u) + c(u, v) + h(v) \ge z'$ for a value $z'$ that $z$
had at some earlier instant, or lowers $g(v)$ to $g(u) + c(u, v)$ when
this is smaller (or $v$ is new) and creates the entry of the new value,
while still holding $u$. (c) A goal selected with $g(u) < z$ sets
$z \leftarrow g(u)$ and records $u$. All the searches of this chapter
are well-formed, as the proofs of their sections show, and so are those
of [Chapter 9](09-hash-distributed.md#ch-9) and
[Chapter 11](11-safepbnf.md#ch-11). We remark that pruning against a value
$z'$ of the past is sound because $z$ never increases, and hence
$z' \ge z$.

<a id="thm-8-6"></a>**Theorem 8.6** (The stopping criterion). Let $h$
be admissible and the search well-formed, and let $\tau$ be an instant
at which no entry in flight that is not stale has $f < z$. Then
$z = C^*$ at $\tau$ (with $C^* = \infty$ if no goal is reachable), and
the goal recorded with $z$ is reached by an optimal path.

*Proof.* The argument is that of [Lemma 2.7](02-best-first.md#lem-2-7),
with the entries in flight in the place of OPEN. If no goal is
reachable, $z = \infty = C^*$, since $z$ is only set to the $g$ of a
goal. Otherwise, suppose $z > C^*$ at $\tau$, and let $P = (s_0, \ldots,
s_k)$ be an optimal path to a goal $s_k$, with $g^*(s_k) = C^*$. Call
the expansion of $s_j$ *done* if a worker has selected an entry $(s_j,
g^*(s_j))$ and finished processing it before $\tau$. Clearly, the
expansion of $s_k$ is not done, since a goal selected with $g = C^* < z$
would have set $z \le C^*$ by (c); let $j$ be the least index whose
expansion is not done. If $j = 0$, the entry $(s_0, 0)$ created at the
start is in flight by (a), or held by a worker that has selected it. If
$j > 0$, the expansion of $s_{j-1}$ with $g^*(s_{j-1})$ is done. By (b)
it did not prune $s_j$, since $g^*(s_{j-1}) + c(s_{j-1}, s_j) + h(s_j) =
g^*(s_j) + h(s_j) \le C^* < z \le z'$ ($h$ being admissible). Hence, it
left $g(s_j) \le g^*(s_j)$, creating the entry $(s_j, g^*(s_j))$ if
$g(s_j)$ was larger. In both cases $g(s_j) = g^*(s_j)$ from then on,
since $g(s_j)$ is the cost of a path (cf. [Lemma
2.6](02-best-first.md#lem-2-6), whose induction holds in the parallel
searches because each assignment of $g$ is made under the lock of its
node, or by its single owner). Therefore, the entry $(s_j, g^*(s_j))$ is
in flight, or it has been selected; in the second case its processing is
not finished (by the choice of $j$), and it is held. It is not stale,
since $g(s_j)$ cannot decrease below $g^*(s_j)$ and $s_j$ can only have
been closed by the selection of an entry with the same $g$, i.e., by a
processing that is not finished. Its key is $g^*(s_j) + h(s_j) \le
g^*(s_j) + h^*(s_j) = C^* < z$, which contradicts the hypothesis. Hence
$z \le C^*$; since $z$ is the $g$ of a goal, i.e., the cost of a path to
it, $z = C^*$. $\square$

We draw the theorem in [Figure 8.1](08-shared-memory.md#fig-8-1): a
search may stop when the least $f$ in flight is at least the incumbent,
and only then. In fact, the searches differ only in how they decide that
this condition holds, which is easy when one lock or one master sees all
the entries (SPA$^*$ and the farm), needs a counter when the OPEN list
is spread over many queues (MQA$^*$), and is the delicate part of a
decentralized search, in which no worker sees the others at the same
instant ([[1]](Z-bibliography.md#ref-1), section 3.3.1). Once the search
has stopped, the path is rebuilt from the recorded goal through the
parent pointers, after all the workers have returned. Its cost is at
most $z$ by the inequality of [Lemma 2.10](02-best-first.md#lem-2-10),
which holds in the parallel searches at the end for the same reason as
$g$ above, and therefore exactly $C^*$.

![Figure 8.1. The stopping criterion of Theorem 8.6. In (a) an
entry with $f = 6$ is in a message and one with $f = 8$ is being
expanded, both below the incumbent $z = 11$, and the search may not stop,
even if all the OPEN lists held only the entries of $f \ge 11$. In (b) a
goal of cost 8 has been found; all the entries in flight have
$f \ge z$ (and are pruned), and $z = 8$ is
optimal.](../figures/shared-memory-termination.svg){#fig-8-1}

<a id="prop-8-7"></a>**Proposition 8.7** (Finiteness). On a finite graph
with nonnegative costs, a well-formed search creates finitely many
entries. Hence, if each of its workers selects an entry, or waits, in a
finite time, and a worker without an entry to select waits only for a
step of another worker, the search stops.

*Proof.* An entry is created only when $g(v)$ takes a new value, which
is strictly smaller than the previous one and is the cost of a path from
$s_0$ to $v$. The argument (i) of [Theorem
2.11](02-best-first.md#thm-2-11) shows that these values are finitely
many for each of the finitely many states. Each selection removes an
entry from the flight set; hence, after finitely many selections no
entry is left and the stopping condition of the search (whose form is
proved for each search below) is met. $\square$

On an infinite graph a parallel search needs, as A$^*$ does, costs at
least some $\varepsilon > 0$, a finite branching and a reachable goal
(argument (ii) of [Theorem 2.11](02-best-first.md#thm-2-11), applied to
the entries with $f < C^*$). With costs that only avoid nonpositive
cycles the search may run forever (the example after Theorem 2.11).

<a id="rem-8-8"></a>**Remark 8.8** (Idle flags are not enough). A test
of termination that only polls each worker, possibly twice in a row,
whether it is idle is not sound when the entries travel in messages:
worker 0 may send an entry of small $f$ to worker 1 and become idle,
while worker 1, idle, has not received it yet. Both answer "idle" in
both rounds if the message is slow, the search stops, and the entry,
which may be on the only optimal path, is lost. A sound test must also
count the entries in transit, as condition (iii) of Definition 8.5 and
Theorem 8.6 require. In this chapter, the centralized searches count
them, as we have seen, by construction (a counter under the lock of
OPEN, an atomic counter, the tasks out at the master). The decentralized
ones use the counting of Dijkstra and Scholten
[[35]](Z-bibliography.md#ref-35) in shared memory and the four counters
of Mattern [[11]](Z-bibliography.md#ref-11) over MPI, which [Chapter
9](09-hash-distributed.md#ch-9) treats.

Finally, all the parallel searches of the library skip the move back
to the parent in an `OperatorDomain` (`Domain.h:220`), as the
sequential A$^*$ does ([§2.9](02-best-first.md#sec-2-9) and
[[4]](Z-bibliography.md#ref-4)): the inverse of a move costs as much
as the move, and therefore the parent is never reached more cheaply through
its child, and the lookup of about one successor in three of the
sliding-tile puzzles is saved. Since the operator travels with the
node, this changes neither the well-formedness nor the proofs above,
by [Proposition 2.16](02-best-first.md#prop-2-16).

## 8.4 Simple Parallel A* {#sec-8-4}

The simplest parallel A$^*$ lets the $p$ workers share the two lists of
the sequential search, i.e., each worker takes the best entry of the
shared OPEN, expands its node and puts the successors back
([[1]](Z-bibliography.md#ref-1), Algorithm 2, and section 2.2 of
[[3]](Z-bibliography.md#ref-3)). Thus, the order of the expansions is
usually close to that of A$^*$, but $p$ nodes are expanded at once, and
a node may be expanded before the node that gives it a cheaper path.
There are no owners, hence no communication in the sense of (8.2), and
the load is balanced by construction. Its price is the coordination,
since each selection and each insertion goes through the lock of OPEN.

### The algorithm

**Algorithm 8.1** (SPA$^*$). Here, the lock $L$ protects OPEN and the
counter
busy of the workers that hold an entry; the fields of the nodes that
change are protected by the locks of CLOSED
([§8.5](08-shared-memory.md#sec-8-5)), and $z$ is atomic.

```text
 1  procedure SPA*(s₀, p)
 2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {(f(s₀), g(s₀), s₀)}
 3     z ← ∞;  busy ← 0;  done ← false
 4     in parallel, on each worker i = 0, ..., p − 1 do
 5        while not done do
 6           lock L
 7           if OPEN has an entry of f < z then
 8              (f, g, u) ← the first entry of OPEN, removed
 9              busy ← busy + 1
10           else if busy = 0 then done ← true
11           unlock L
12           if no entry was removed then wait a little;  continue
13           F ← ∅;  lock u
14           if u is closed or g ≠ g(u) then u is not expanded  ▷ stale
15           else close u
16           unlock u
17           if u is expanded and u is a goal then
18              z ← min(z, g(u)), recording u if g(u) < z
19           else if u is expanded then
20              for each successor v of u with g(u) + c(u,v) + h(v) < z do
21                 g' ← g(u) + c(u, v);  lock v
22                 if v ∉ CLOSED or g' < g(v) then
23                    g(v) ← g';  parent(v) ← u;  v is not closed
24                    F ← F ∪ {(g' + h(v), g', v)}
25                 unlock v
26           lock L;  push the entries of F in OPEN;  busy ← busy − 1
27           unlock L
28     return the path to the goal of cost z, if z < ∞
```

This is the pseudocode of `SPAStar.h:16`. A worker takes the lock of
OPEN only twice per expansion, to select an entry (lines 6 to 11) and to
push all the successors of its node at once (lines 26 and 27). The
expansion itself, at lines 13 to 25, runs without it, and this is where
the parallelism is. [Figure 8.2](08-shared-memory.md#fig-8-2) draws one
selection and one insertion. At line 14, the test of staleness is the
lazy deletion of the sequential search, done under the lock of the node
since another worker may be lowering its $g$. The goal test is at the
selection (line 17), and a goal does not end the search, which goes on
until line 10 finds no entry of $f < z$ and no busy worker.

![Figure 8.2. A step of SPA$^*$ with three workers. (a) Worker
$w_2$ is expanding a node and busy $= 1$. (b) Worker $w_0$ takes $L$,
removes the best entry $a$ of OPEN and raises busy, while $w_1$ waits
for $L$. (c) After the expansion of $a$, outside $L$, $w_0$ takes $L$
again, pushes its successors $c$ and $v$ and lowers busy in the same
critical section.](../figures/shared-memory-spa-step.svg){#fig-8-2}

### Invariants and correctness

We now prove that SPA$^*$ is correct. The proof rests on one invariant,
which ties the counter busy to the entries that are not in OPEN.

<a id="lem-8-9"></a>**Lemma 8.9** (The invariant of busy). Whenever $L$
is free, every entry in flight is in OPEN or is held by one of the
busy workers, i.e., by a worker between lines 9 and 26, and busy is the
number of these workers.

*Proof.* Since SPA$^*$ has no messages, an entry in flight is in OPEN
or held (Definition 8.5). An entry leaves OPEN only at line 8,
where the worker raises busy in the same critical section, and it is
held (with the entries of $F$ that it creates) until line 26, where the
worker pushes $F$ and lowers busy in the same critical section. Since
both changes of busy are made under $L$ together with the change of
OPEN, the claim holds whenever $L$ is released. $\square$

<a id="thm-8-10"></a>**Theorem 8.10** (Correctness of SPA$^*$). If $h$
is admissible and SPA$^*$ ends, it returns an optimal path (or none,
when no goal is reachable); on a finite graph with nonnegative costs it
ends.

*Proof.* SPA$^*$ is well-formed: (a) and (c) are lines 8, 24, 26 and 18,
and (b) is lines 20 to 25, where the update of $v$ is made under its
lock, and where $z$ is read without a lock. Hence, the pruning at line
20 uses a value that $z$ had at some instant of the expansion. When a
worker sets done at line 10, it holds $L$ and finds no entry of $f < z$
in OPEN and busy $= 0$. By Lemma 8.9 no entry in flight is held, and
hence none has $f < z$, and Theorem 8.6 gives $z = C^*$. A worker that
finds no entry to select at line 12 waits only for the busy workers,
which end their expansion without waiting for anything but locks held
for a bounded number of steps; hence Proposition 8.7 gives the
termination. $\square$

### Re-expansions

Unlike A$^*$ (cf. [Theorem 2.13](02-best-first.md#thm-2-13)), SPA$^*$
may expand a node more than once even with a consistent heuristic
([[1]](Z-bibliography.md#ref-1), section 3.2): in [Figure
8.3](08-shared-memory.md#fig-8-3), two workers expand $a$ and $b$ at
once, the second best entry $b$ generates $v$ through an expensive edge,
and $v$ is selected and closed by $w_1$ before $w_0$, through $c$, finds
the cheaper path to it. As one can check on the figure, lines 22 and 23
then reopen $v$ (`v is not closed`), and it is expanded again; the
counter `reexpanded` counts these expansions, which belong to the search
overhead (8.1). With $p = 1$, SPA$^*$ is A$^*$, and the re-expansions
are those of an inconsistent heuristic
([§2.7](02-best-first.md#sec-2-7)).

![Figure 8.3. A re-expansion in SPA$^*$ with a consistent
heuristic (values $g/h$ next to the states). (a) $w_0$ and $w_1$
expand $a$ and $b$ at once. (b) $w_1$ selects $v$, of $f = 5$, with
$g(v) = 4$, while $w_0$ expands $c$. (c) The expansion of $c$ gives
$g(v) = 3$: $v$ is reopened and will be expanded again, now with its
optimal $g$.](../figures/shared-memory-reopen.svg){#fig-8-3}

### The code

In `SPAStar.h`, the loop of a worker is the following excerpt; the
backoff of an idle worker spins for a growing number of iterations
before it yields its core, so that a worker that finds OPEN empty only
for a moment does not sleep.

```cpp
--8<-- "snippets/shared-memory-spa-run.cpp"
```

Within an expansion, the worker computes the hashes of all the
successors first and prefetches their places in CLOSED before it updates
any of them, and the cache misses of their lookups therefore overlap.
The successors with $f \ge z$ are pruned at the generation, with a
relaxed read of $z$:

```cpp
--8<-- "snippets/shared-memory-spa-successors.cpp"
```

Furthermore, the implementation makes the following choices. The lock of
OPEN and the counter busy sit on a separate cache line
(`SPAStar.h:212`), and so does each worker with its counters
(`SPAStar.h:227`); thus, the counters of two workers never share a line.
All the successors are pushed in one critical section, which takes $L$
once per expansion instead of once per successor, and the goal is
recorded under a separate mutex (`SPAStar.h:290`), which is taken only
when a goal is selected. Nevertheless, the lock of OPEN is contended by
all the workers twice per expansion, and the survey reports that, when
the expansions are as fast as on the sliding-tile puzzle, SPA$^*$ is
usually slower than the sequential A$^*$, even with lock-free data
structures ([[1]](Z-bibliography.md#ref-1), section 3.2, after
[[3]](Z-bibliography.md#ref-3)). It may pay off when the expansions are
slow, since the time under $L$ is then a small part of an expansion.
MQA$^*$ ([§8.6](08-shared-memory.md#sec-8-6)) removes this lock.

## 8.5 The shared closed list {#sec-8-5}

SPA$^*$ and MQA$^*$ share one CLOSED list, the class
`ConcurrentClosed< D , Lock >` of `ConcurrentClosed.h`, which maps each
state to its canonical node and protects the fields of a node that
change during the search ($g$, the parent, the operator and the flags
closed and expanded); $h$ never changes and is read freely. Its
operation `update(t, hash, g, h, parent, op, ...)` records the path to
$t$ of cost $g$: it returns the node of $t$ when $t$ is new or $g$ is
smaller than its $g$, which is then set with the parent and the
operator and makes the node not closed, and returns a null pointer
otherwise,
counting a duplicate. Its implementation is chosen at build time with
the set of data structures (`Structures.h`).

With the default data structures, CLOSED is a lock-free hash table of
pointers to the nodes, compared by their states, one word of 8 bytes per
node, which a lookup reads and an insertion fills by a compare-and-swap
without a lock. The fields that change are protected by 256 locks, the
lock of a node being chosen by a multiplicative hash of its address
(line 116 of the default implementation). Such a lock is taken only
when the state is already in the table and $g$ is smaller than its $g$
read without the lock ([Figure 8.4](08-shared-memory.md#fig-8-4)). With
the data structures of the standard library or of Boost, CLOSED is split
into 256 shards by the top byte of the hash of the state, each a
sequential table with a lock, which also protects the nodes of its
states (`structures/std/ConcurrentClosed.h:128`).

We now show that the shared CLOSED loses no improvement of $g$, although
most of its operations take no lock.

<a id="lem-8-11"></a>**Lemma 8.11** (The lock-free update). With the
lock-free table, (i) each state has at most one node; (ii) the $g$ of a
node never increases; and (iii) if `update` returns a null pointer, the
$g$ of the node was at most the new value at some instant of the call,
and if it returns the node, it has set $g$ to the new value, which was
smaller than $g$ at that instant. Hence no improvement of $g$ is lost,
and each improvement creates exactly one entry.

*Proof.* (i) An insertion fills an empty slot by a compare-and-swap, and
an insertion that loses the race finds the slot taken by a node of the
same state and does not insert. (ii) The only write of $g$ is made under
the lock of the node, after a test that the new value is smaller. (iii)
The read without the lock (an atomic read, line 138 of the default
implementation) returns a value that $g$ had at some instant; since $g$
never increases, the current value is not larger. If the value read is
at most the new $g'$, the current one is too, and the call is a
duplicate. Otherwise the call takes the lock and tests again, and the
test and the write under the lock are exact, since all the writes are
made under it. $\square$

The duplicates, which are often the most frequent case of a search with
many transpositions, thus take no lock and write no shared line. A
worker fills a spare node with the state before the insertion and keeps
it for the next call when the state was already there. Hence, a node is
allocated once per state inserted; the nodes are never moved, and a
pointer to a node stays valid until the end of the search.

![Figure 8.4. The three cases of `update` on the lock-free CLOSED.
(a) A new state: a compare-and-swap fills the empty slot with the
spare node of the worker. (b) A duplicate: the $g$ read without a lock
is not larger than $g'$, and no lock is taken. (c) A better path: the
worker takes the stripe of locks of the node, tests again, writes $g$,
the parent and the operator and reopens the
node.](../figures/shared-memory-closed.svg){#fig-8-4}

## 8.6 MultiQueue A* {#sec-8-6}

MQA$^*$ keeps the shared CLOSED of SPA$^*$ and replaces its OPEN, a
priority queue under one lock, by a relaxed priority queue, the
MultiQueue of Rihani, Sanders and Dementiev
[[37]](Z-bibliography.md#ref-37): $n = 2p$ sequential priority queues,
each with a separate lock, a push into a queue drawn at random, and a
pop from the better of two queues drawn at random ([Figure
8.5](08-shared-memory.md#fig-8-5)). Two workers then contend only when
they happen to draw the same queue (which is rather rare with $2p$
queues), and even then the one that finds the lock taken draws again
instead of waiting, while the order of the pops typically stays close to
that of the keys.

![Figure 8.5. The OPEN list of SPA$^*$ and that of MQA$^*$ with
$p = 4$ workers. (a) In SPA$^*$ the four workers contend for one lock,
and three of them wait. (b) In MQA$^*$ each worker reads the least keys
of two queues drawn at random (the pairs are adjacent here only for the
sake of the drawing) and pops from the better one, under the lock of
that queue only.](../figures/shared-memory-contention.svg){#fig-8-5}

### The MultiQueue

Each queue keeps the key of its minimum in an atomic word, which the
other workers read without its lock. A pop reads the minima of two
distinct queues drawn at random, locks the one with the smaller key, if
its lock is free (and draws two other queues otherwise), and removes its
minimum. A push puts the pair in a random queue whose lock it takes at
the first try. Hence, the pop returns a pair that is not always the
least of all ([Figure 8.6](08-shared-memory.md#fig-8-6)): its *rank
error* is the number of the pairs present with a smaller key. In the
sequential model of the analysis (the pairs spread at random over the
$n$ queues and removed by the better of two random queues), the expected
rank error is linear in $n$, and its maximum over a polynomial number of
pops is $O(n \log n)$ with high probability, while a pop from a single
random queue would make it grow without bound (the bounds are those
recalled in [[38]](Z-bibliography.md#ref-38)). Since SPA$^*$ expands $p$
nodes at once anyway, a rank error of order $n = 2p$ departs from the
order of A$^*$ by somewhat more than SPA$^*$ already does, and arguably
not by much. Also, a pop returns nothing only after a scan of all the
queues, each locked in turn, has found them empty. With pushes running
at the same time this does not prove that the structure is empty, and
MQA$^*$ therefore does not use it for its termination.

![Figure 8.6. A pop from a MultiQueue of six queues (their least
keys on top). (a) $w_0$ draws $q_1$ and $q_4$ and compares their least
keys, 3 and 5. (b) It removes 3 from $q_1$, whose next key becomes its
top; the rank error of this pop is 1, since the key 2 of $q_3$ is
smaller.](../figures/shared-memory-mq-pop.svg){#fig-8-6}

Two refinements of the MultiQueue save cache misses. A worker may keep
the two queues it drew for $s$ operations, the *stickiness*, which keeps
their lines in its cache, and draw again after $s$ operations or when a
lock is held ([Figure 8.7](08-shared-memory.md#fig-8-7)). The rank
errors then grow (about linearly in $s$, according to the measures
available, since there is no theory yet), and MQA$^*$ uses $s = 1$,
i.e., each operation draws its two queues (`MQAStar.h:239`), which keeps
the bounds above. Also, a worker that has several pairs to push at once,
as the successors of a node, may push them all into one queue under one
lock (`pushAll`, `MQAStar.h:288`). This takes one lock per expansion
instead of one per successor, and relaxes the order as a stickiness of
$|F|$ would, since the successors of a node then sit in the same queue.
MQA$^*$ does so.

![Figure 8.7. The two devices of the MultiQueue that save cache
misses. (a) With stickiness $s = 3$ a worker keeps the pair of queues
it drew for three operations and then draws a new pair. (b) A pushAll
puts all the successors $F$ of a node into one queue of the pair, under
one lock, where they take their places among the keys of that
queue.](../figures/shared-memory-stickiness.svg){#fig-8-7}

### The algorithm and its termination

**Algorithm 8.2** (MQA$^*$). Here, the atomic counter $m$ of the
entries that
are in OPEN or held by a worker replaces busy and the lock $L$ of
SPA$^*$.

```text
 1  procedure MQA*(s₀, p)
 2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {(f(s₀), g(s₀), s₀)}
 3     z ← ∞;  m ← 1;  done ← false
 4     in parallel, on each worker i = 0, ..., p − 1 do
 5        while not done do
 6           e ← a pop of OPEN, possibly nothing
 7           if no entry was popped then
 8              if m = 0 then done ← true
 9              else wait a little
10              continue
11           (f, g, u) ← e;  F ← ∅
12           if f < z then expand u as at lines 13 to 25 of SPA*
13           m ← m + |F| − 1;  push the entries of F in OPEN
14     return the path to the goal of cost z, if z < ∞
```

This is the pseudocode of `MQAStar.h:19`. Since the entries of OPEN are
not ordered across the queues, an entry popped with $f \ge z$ does not
mean that OPEN has no better entry, and it is dropped (line 12) instead
of ending the search. The end is decided by $m$ only ([Figure
8.8](08-shared-memory.md#fig-8-8)).

We now show that $m$ counts all the entries in flight, whatever the
interleaving of the workers.

<a id="lem-8-12"></a>**Lemma 8.12** (The invariant of $m$). At every
instant, $m = |\mathrm{OPEN}| + H + Q$, where $|\mathrm{OPEN}|$ is the
number of entries in the queues, $H$ the number of the workers that hold
a popped entry and have not yet executed the addition of line 13, and
$Q$ the number of entries that have been counted at line 13 and are still
to be pushed. Hence $m = 0$ if and only if OPEN is empty and no worker holds
an entry.

*Proof.* At line 3, $m = 1$ and OPEN holds the root. A pop (line 6)
moves one entry from OPEN to a holding worker and leaves $m$ unchanged.
In turn, the atomic addition of line 13 ends the holding of one worker
($H$ falls by 1) and counts the $|F|$ entries to be pushed ($Q$ rises by
$|F|$), and it adds $|F| - 1$ to $m$. Each push of line 13 moves an
entry from $Q$ to OPEN. The three terms are nonnegative, which gives the
last claim. $\square$

The order of line 13 matters. Indeed, the entries of $F$ are counted *before*
they are pushed, and therefore $m$ never misses an entry. Had they been
counted after the pushes, another worker could pop one of them, expand
it and subtract its 1 before it was added, and read $m = 0$ while
entries were still in flight. Also, the code makes one atomic addition
per expansion (`MQAStar.h:283`), and none when $|F| = 1$, instead of one
before and one after the pushes, which halves the traffic on the cache
line of $m$. All the workers share that line, which holds no other data
(`MQAStar.h:243`).

<a id="thm-8-13"></a>**Theorem 8.13** (Correctness of MQA$^*$). If $h$
is admissible and MQA$^*$ ends, it returns an optimal path (or none,
when no goal is reachable); on a finite graph with nonnegative costs it
ends.

*Proof.* MQA$^*$ is well-formed for the reasons given for SPA$^*$ in
Theorem 8.10, the expansion being the same. A worker sets done only when
it reads $m = 0$, and then, by Lemma 8.12, no entry is in OPEN and none
is held, i.e., none is in flight, and Theorem 8.6 gives $z = C^*$.
Whatever the order of the pops, each entry with $f < z$ that was created
has been expanded. A worker whose pop fails while $m > 0$ waits for the
workers that hold entries or are about to push them, and these end their
step without waiting; Proposition 8.7 gives the termination. $\square$

The order of the pops changes only how many nodes are expanded before
$z$ is optimal and how many are reopened, never the cost of the path.

![Figure 8.8. The counter $m$ of MQA$^*$ in a run with two workers
(the entries with their $f$). The expansion of $s_0$ adds $|F| - 1 = 2$
before its three successors are pushed; the goal $t$ sets $z = 8$ and
adds $-1$; the expansion of $a$ prunes its successors ($F = \emptyset$);
$d$, popped with $f = 9 \ge z$, is dropped; $m = 0$ then shows to
$w_0$, whose pop fails, that no entry is
left.](../figures/shared-memory-mq-counter.svg){#fig-8-8}

### The packed key

With integer costs the key of an entry of MQA$^*$ (and of PA$^*$SE) is
one 64-bit word that packs $f$ and $g$, built by `packed_fg` (`Node.h:119`):

$$
\mathrm{key}(u) = f(u) \cdot 2^{32} + \bigl(2^{32} - 1 - g(u)\bigr) ,
\tag{8.4}
$$

i.e., $f$ in the high half and the complement of $g$ in the low half
([Figure 8.9](08-shared-memory.md#fig-8-9)), and `f_of` takes $f$ back
by a shift.

<a id="lem-8-14"></a>**Lemma 8.14** (The order of the packed key). If
$g(u), g(v) \ge 0$, $h(u), h(v) \ge 0$ and $f(u), f(v) < 2^{32}$, then
$\mathrm{key}(u) < \mathrm{key}(v)$ if and only if $f(u) < f(v)$, or
$f(u) = f(v)$ and $g(u) > g(v)$.

*Proof.* Both halves of (8.4) lie in $[0, 2^{32})$, since
$0 \le g \le f < 2^{32}$ (as $h \ge 0$), and the keys therefore compare
as the pairs $(f, 2^{32} - 1 - g)$ in the lexicographic order, and the
second component decreases as $g$ grows. $\square$

The ties of $f$ thus go to the deeper node, as in the sequential A$^*$
(cf. [§2.8](02-best-first.md#sec-2-8)), and the comparison of two keys
is one comparison of integers. This matters for a relaxed queue rather
more than for an exact one, as we now explain. The last layer $f = C^*$
often holds a large part of the nodes (on the sliding-tile puzzles, for
instance), and a search that prefers the deeper nodes reaches a goal
after a part of that layer. With $f$ alone, instead, the pops of the
layer come in an arbitrary order, and typically much more of it is
expanded before a goal is found. Since a MultiQueue scatters the entries
of one layer over $2p$ queues, only the key can carry the preference for
depth into each of them. Note that the lemma requires $f < 2^{32}$,
which the integer costs of the domains of the library satisfy by a wide
margin. With real costs the key is $f$ only (`MQAStar.h:216`), and the
ties are broken by the queues.

```cpp
--8<-- "snippets/shared-memory-mqa-key.cpp"
```

![Figure 8.9. The packed key (8.4) of MQA$^*$ and PA$^*$SE: $f$ in
bits 63 to 32, the complement $\bar g$ of $g$ in bits 31 to 0. The
node $w$ comes first by its smaller $f$; $u$ comes before $v$, of the
same $f$, since its larger $g$ gives a smaller
$\bar g$.](../figures/shared-memory-packed-key.svg){#fig-8-9}

### The code

In `MQAStar.h`, the loop of a worker is the following; the addition to
$m$ precedes the `pushAll` of the successors, which builds their keys
on the fly.

```cpp
--8<-- "snippets/shared-memory-mqa-run.cpp"
```

As for the expansion and CLOSED, the engineering is that of SPA$^*$
(prefetches, cache lines, lock-free duplicates). Finally, the test
`mqa-thread` runs MQA$^*$ with 1, 2, 3 and 8 workers on two instances of
Korf of the 15-puzzle, on the map of Romania and on the corner cases of
`check_astar_edges` (e.g., an initial state that is a goal, one move
and a node limit), and `mqa-openmp` does the same on OpenMP.

## 8.7 The master-worker farm {#sec-8-7}

The farm with feedback is a master-worker form of the centralized
approach of section 3.2 of [[1]](Z-bibliography.md#ref-1), written with
the farm of the executor of the library (a farm of FastFlow
[[12]](Z-bibliography.md#ref-12) on that runtime). A master owns OPEN,
CLOSED, all the nodes and the incumbent, and therefore needs no lock. It
sends the best nodes of OPEN, in tasks of up to `batch` nodes, to the
workers that are free (the farm schedules them on demand). A worker
generates the successors of the nodes of its task, with their costs and
their heuristic values, and sends them back on the feedback channel
([Figure 8.10](08-shared-memory.md#fig-8-10)). The master then checks
them against CLOSED and puts the new and the improved ones in OPEN. In
other words, the scheme is that of SPA$^*$ with the master in the place
of the locks.

![Figure 8.10. The farm with feedback: only the master holds OPEN,
CLOSED, the incumbent $z$ and the number out of the tasks at the
workers; the tasks (solid) carry up to `batch` nodes to a free worker,
and the feedback channel (dashed) brings back their successors with
their costs and heuristic
values.](../figures/shared-memory-farm.svg){#fig-8-10}

**Algorithm 8.3** (the farm). INSERT is the insertion of the sequential
A$^*$ on the lists of the master (a new state goes in OPEN; a state
reached with a smaller $g$ is updated, reopened and put in OPEN again;
otherwise it is a duplicate).

```text
 1  procedure MASTER(s₀, p)
 2     z ← ∞;  out ← 0;  INSERT(s₀, 0, nil)
 3     loop
 4        DISPATCH()
 5        if out = 0 and OPEN has no node of f < z then return
 6        receive a task T back from a worker;  out ← out − 1
 7        for each successor v of a node u of T, with c(u, v) and h(v) do
 8           if g(u) + c(u, v) + h(v) < z then INSERT(v, g(u) + c(u, v), u)
 9  procedure DISPATCH()
10     while out < p and OPEN has a node of f < z do
11        T ← ∅
12        while |T| < batch and OPEN has a node of f < z do
13           u ← the first node of OPEN, removed;  close u
14           if u is a goal then z ← min(z, g(u)), recording u if g(u) < z
15           else T ← T ∪ {u}
16        if T = ∅ then return
17        send T to a free worker;  out ← out + 1
18  procedure WORKER(T)                       ▷ on each of the p workers
19     for each u ∈ T, add to T the successors v of u, with c(u, v), h(v)
20     send T back to the master
```

This is the pseudocode of `FarmAStar.h:17`; the goal test is at the
selection (line 14), which is the expansion point of the master, and
[Figure 8.11](08-shared-memory.md#fig-8-11) draws one return of a task.

As for SPA$^*$, one can now prove the correctness of the farm from Theorem 8.6.

<a id="thm-8-15"></a>**Theorem 8.15** (Correctness of the farm). If
$h$ is admissible and the farm ends, it returns an optimal path (or
none, when no goal is reachable); on a finite graph with nonnegative
costs it ends.

*Proof.* Here, the entries in flight are those of OPEN and the nodes of
the tasks that are out, together with their successors. The counter out
counts the tasks; only the master changes OPEN, $z$ and out, and it
handles one returned task at a time. Thus, the search is well-formed: a
node selected at line 13 is either a goal (line 14) or carried in a
task, whose successors are inserted at line 8, before the next test of
line 5. Note that $g(u)$ at line 8 is the current $g$ of the parent,
which may be smaller than when the task left, if $u$ was improved
meanwhile. It is still the cost of a path to $u$, and $g(u) + c(u, v)$
is therefore the cost of a path to $v$, and $u$, reopened, is selected
again. Hence, the test of line 5 is exact: out $= 0$ and no node of $f <
z$ in OPEN mean that no entry with $f < z$ is in flight, and Theorem 8.6
gives $z = C^*$. A worker never waits for anything but a task, and the
master waits at line 6 only when out $> 0$; Proposition 8.7 gives the
termination. $\square$

![Figure 8.11. A return of a task. (a) The task of $w_1$, the node
$b$, comes back with its successors $v$ (with $g = 3$ and $h = 3$) and
$a$ (with $g = 2$). (b) The master inserts $v$ in OPEN, counts $a$,
already closed with $g = 1$, as a duplicate, and dispatches the best
node $d$ of OPEN to the free worker $w_1$; out is 3
again.](../figures/shared-memory-farm-step.svg){#fig-8-11}

A task that returns is handled by the following method of the master,
which computes the hashes of all the successors of the task and
prefetches their places in CLOSED before inserting any of them, as the
sequential search does with the successors of a node:

```cpp
--8<-- "snippets/shared-memory-farm-back.cpp"
```

The tasks are reused through a list of free tasks, and the farm
allocates no memory for them in the steady state. The farm runs
`threads` + 1 threads, and with `PBEFS_PIN` the master takes the first
core of the affinity mask and the worker $i$ the core $i + 1$
(`FarmAStar.h:165`). Since the master does all the bookkeeping (the
lookups in CLOSED, the updates of OPEN and the selections), it is
typically a bottleneck unless the generation of the successors, the only
part done by the workers in parallel, costs much more than the
bookkeeping of the same successors. Also, the batch trades the number of
the tasks against the order of the expansions; with tasks of $b$ nodes,
up to $p\,b$ nodes are out at once, and the nodes of a task are the $b$
best of OPEN at the time of the dispatch, not at the time of their
expansion. Finally, the test `farm-thread` (and `farm-openmp`,
`farm-fastflow`) runs the farm with 1 to 8 workers and tasks of 1 and of
4 nodes, on the 15-puzzle, on a random grid and on the map of Romania,
proves a grid split by a wall unsolvable, and checks the node limit.

## 8.8 Parallel A* for slow expansions {#sec-8-8}

The searches above expand in parallel nodes whose $g$ may not be
optimal yet, and pay for it with re-expansions. PA$^*$SE of Phillips,
Likhachev and Koenig [[6]](Z-bibliography.md#ref-6) expands in parallel
only states whose $g$ is already optimal (or within a factor of it),
proved so by a pairwise heuristic, at the price of a test before each
expansion; hence, no state is expanded twice.

### Independence and safety

Let $h(t, s)$ be the pairwise heuristic of the domain, an estimate of
the cost $d(t, s)$ of a path from $t$ to $s$, and let $h(s)$ be the
heuristic of the search, an estimate of the cost $h^*(s)$ from $s$ to
the nearest goal. PA$^*$SE needs the following properties, which we call
(H): (H1) $h(t, s) \le d(t, s)$ for all $t$ and $s$ (the pairwise
heuristic is admissible), (H2) $h$ is admissible, and (H3)
$h(t) \le h(t, s) + h(s)$ for all $t$ and $s$ (the forward-backward
consistency of [[6]](Z-bibliography.md#ref-6), in the form that the
proofs below use). The library leaves them to the `Domain`
(`Domain.h:192`); the pairwise heuristics of the library are the
distance between two cells of a grid, the sum of the displacements of
the tiles of a sliding-tile puzzle and the cells that are dirty in $t$
and clean in $s$ of the vacuum world, each with the matching heuristic
$h$.

<a id="def-8-16"></a>**Definition 8.16** (Independence and safety).
For $\varepsilon_p \ge 1$, a state $s$ is *independent* of a state $t$
if

$$
g(s) - g(t) \le \varepsilon_p \, h(t, s) .
\tag{8.5}
$$

A state $s$ of OPEN is *safe* if it is independent of all the states
being expanded (the set BE) and of each state $t$ of OPEN with
$f(t) < f(s)$.

Intuitively, a path through $t$ reaches $s$ with a cost at least $g(t) +
d(t, s) \ge g(t) + h(t, s)$; hence, if $s$ is independent of $t$, such a
path cannot lower $g(s)$ below $g(s) / \varepsilon_p$ (for $g(t)$
optimal). The states of OPEN with $f(t) \ge f(s)$ need no check, by
Theorem 2 of [[6]](Z-bibliography.md#ref-6), which is also what the
proof below uses. [Figure 8.12](08-shared-memory.md#fig-8-12) draws
(8.5) on a line of states: the states that depend on $t$ are those above
a cone with its apex in $t$ and slopes $\pm \varepsilon_p$, and a larger
$\varepsilon_p$ makes the cone steeper and more states independent.

![Figure 8.12. The independence (8.5) on a line of states, with
$h(t, s) = |x_t - x_s|$ and $g$ on the vertical axis. The state $t$,
being expanded, casts a cone of slope $\pm 1$ for $\varepsilon_p = 1$
(solid, shaded above) and $\pm 2$ for $\varepsilon_p = 2$ (dashed); a
state above the cone of $t$ may still get a smaller $g$ through $t$
and is not safe. $s_3$ depends on $t$ for both values, $s_1$ only for
$\varepsilon_p = 1$, and $s_2$ is independent of
$t$.](../figures/shared-memory-pase-cones.svg){#fig-8-12}

### The algorithm

**Algorithm 8.4** (PA$^*$SE). One lock $L$ protects OPEN, BE and
CLOSED; OPEN is an ordered set by increasing $f$ and, for equal $f$, by
decreasing $g$.

```text
 1  procedure PA*SE(s₀, p)
 2     CLOSED ← {s₀}, with g(s₀) ← 0;  OPEN ← {s₀};  BE ← ∅
 3     in parallel, on each worker do
 4        lock L
 5        while not done do
 6           s ← the first state of OPEN, among its first 64, that is SAFE
 7           if there is none then
 8              if OPEN = ∅ and BE = ∅ then done ← true;  break
 9              wait on L until OPEN or BE change;  continue
10           remove s from OPEN
11           if s is a goal then record s;  done ← true;  break
12           add s to BE;  close s;  unlock L
13           S ← the successors of s, with their costs   ▷ the slow part
14           lock L
15           for each (v, c(s, v)) ∈ S do
16              g' ← g(s) + c(s, v)
17              if v ∉ CLOSED then
18                 add v to CLOSED and OPEN, with g(v) ← g', parent(v) ← s
19              else if v is not closed and g' < g(v) then
20                 g(v) ← g';  parent(v) ← s;  reorder v in OPEN
21           remove s from BE;  wake the waiting workers
22        unlock L
23  procedure SAFE(s)
24     for each t ∈ BE do
25        if g(s) − g(t) > ε_p h(t, s) then return false
26     for each t ∈ OPEN with f(t) < f(s) do
27        if g(s) − g(t) > ε_p h(t, s) then return false
28     return true
```

This is the pseudocode of `PASE.h:17`, i.e., Algorithm 1 of
[[6]](Z-bibliography.md#ref-6) with the scan of line 6 bounded. [Figure
8.13](08-shared-memory.md#fig-8-13) draws one selection, in which the
first state of OPEN depends on the state being expanded and is skipped,
and the second is safe and moves to BE. Note that, by line 19, a closed
state is never improved and never reopened, which is sound only because
of the next lemma.

![Figure 8.13. A selection of PA$^*$SE. (a) $w_0$ expands $t$;
under $L$, $w_1$ scans OPEN: $q$ depends on $t$, while $r$ is
independent of $t$ and of $q$, the only state of OPEN with a smaller
$f$. (b) $r$ leaves OPEN and joins BE, and $w_1$ releases $L$ and
expands it.](../figures/shared-memory-pase-step.svg){#fig-8-13}

### Correctness

We first show that a state found safe has a nearly optimal $g$; this
allows line 19 to leave the closed states alone.

<a id="lem-8-17"></a>**Lemma 8.17** (Safe states are nearly optimal).
Under (H), each state $s$ that PA$^*$SE finds safe at line 6 has
$g(s) \le \varepsilon_p \, g^*(s)$; with $\varepsilon_p = 1$, each
state is expanded with its optimal $g$, and hence at most once.

*Proof.* By induction on the selections (cf. the proof of Theorem 8.6,
of which this is the analogue for BE). Suppose that the states selected
before have been selected with $g \le \varepsilon_p g^*$, and that $s$
is found safe with $g(s) > \varepsilon_p g^*(s)$. Let $P$ be an optimal
path from $s_0$ to $s$, and call a state of $P$ *done* if it is closed
and outside BE, i.e., its expansion is complete. Let $t$ be the first
state of $P$ that is not done; it exists, since $s$ is in OPEN. If $t =
s_0$, then $t$ is in OPEN or in BE with $g(t) = 0$. Otherwise the
predecessor $t'$ of $t$ on $P$ is done, and it was selected with $g(t')
\le \varepsilon_p g^*(t')$. Its expansion left $g(t) \le g(t') + c(t',
t) \le \varepsilon_p (g^*(t') + c(t', t)) = \varepsilon_p g^*(t)$, using
$\varepsilon_p \ge 1$, unless $t$ was already closed, and then $t$,
which is not done, is in BE (line 19 never updates a closed state, and a
closed state that is not done is being expanded). In all cases $t$ is in
OPEN or in BE, with $g(t) \le \varepsilon_p g^*(t)$, and $t \ne s$,
since, if $t = s$, the bound just proved contradicts the assumption on
$s$. Since $t$ precedes $s$ on an optimal path, $g^*(s) = g^*(t) + d(t,
s)$, and

$$
g(s) - g(t) > \varepsilon_p g^*(s) - \varepsilon_p g^*(t)
 = \varepsilon_p \, d(t, s) \ge \varepsilon_p \, h(t, s)
$$

by (H1), and $s$ depends on $t$. If $t \in$ BE, $s$ is not safe by
line 25. If $t \in$ OPEN, then by (H3), (H1) and $\varepsilon_p \ge 1$

$$
f(t) = g(t) + h(t) \le \varepsilon_p g^*(t) + h(t, s) + h(s)
\le \varepsilon_p \bigl(g^*(t) + d(t, s)\bigr) + h(s)
= \varepsilon_p g^*(s) + h(s) < g(s) + h(s) = f(s) ,
$$

and therefore $t$ is among the states checked at line 26, and $s$ is
not safe. Both cases contradict the safety of $s$. With $\varepsilon_p = 1$
the bound is $g(s) = g^*(s)$, and a state, once closed, is never put
back in OPEN (line 19). $\square$

For the goal returned at line 11 we need one more property, which holds
trivially when there is one goal state: (G) $h(t, s) \le h^*(t)$ for all
the states $t$ and goals $s$, i.e., the pairwise heuristic towards any
goal does not exceed the distance to the nearest goal. With one goal
$s$, (G) is (H1); in the vacuum world, whose goals are all the states
without dirt (with the robot anywhere), $h(t, s)$ is the number of the
dirty cells of $t$ (`VacuumWorld.h:203`), which is at most $h^*(t)$,
since each of them needs a separate move. The grids and the puzzles of
the library have one goal state.

<a id="thm-8-18"></a>**Theorem 8.18** (Bounded suboptimality of
PA$^*$SE). Under (H) and (G), if PA$^*$SE returns a goal $s$, then
$g(s) \le \varepsilon_p C^*$; with $\varepsilon_p = 1$ the path returned
is optimal.

*Proof.* Suppose $g(s) > \varepsilon_p C^*$, and let $P$ be an optimal
path to a goal $s'$, with $g^*(s') = C^*$. As in the proof of Lemma
8.17, the first state $t$ of $P$ that is not done is in OPEN or in BE
with $g(t) \le \varepsilon_p g^*(t)$, and $t \ne s$, since $g(t) \le
\varepsilon_p g^*(t) \le \varepsilon_p C^*$ for all the states of $P$.
Then $f(t) \le \varepsilon_p g^*(t) + h^*(t) \le \varepsilon_p (g^*(t) +
h^*(t)) = \varepsilon_p C^* < g(s) = f(s)$ by (H2). Hence, $t$ is
checked by SAFE whether it is in BE or in OPEN, and the safety of $s$
gives $g(s) - g(t) \le \varepsilon_p h(t, s) \le \varepsilon_p h^*(t)$
by (G). Hence $g(s) \le \varepsilon_p g^*(t) + \varepsilon_p h^*(t) =
\varepsilon_p C^*$, a contradiction. Finally, the path rebuilt through
the parents has cost at most $g(s)$, as in [Lemma
2.10](02-best-first.md#lem-2-10). $\square$

These are Theorems 1 and 3 of [[6]](Z-bibliography.md#ref-6), stated
there for one goal state. The proofs above also make explicit where the
forward-backward consistency (H3) enters, i.e., only to bring the states
of OPEN with $f(t) < f(s)$ into the check. Without (G) a domain with
several goal states may let PA$^*$SE return a goal that is not optimal,
since the cheaper goal may be reached through a state of which the
returned one is independent; Exercise 8.6 is about an example.

Finally, we show that PA$^*$SE neither waits forever nor runs forever on
a finite graph.

<a id="prop-8-19"></a>**Proposition 8.19** (Progress and termination of
PA$^*$SE). If OPEN is not empty and BE is empty, the first state of OPEN
is safe. Hence no worker waits forever while another could proceed, and
on a finite graph PA$^*$SE ends, with a goal or with OPEN and BE empty,
in which case no goal is reachable.

*Proof.* With BE empty and $s$ the first state of OPEN, no state of OPEN
has $f < f(s)$, and SAFE returns true; the first state is among the
first 64 that line 6 scans. A worker waits at line 9 only when BE is not
empty, i.e., while another worker expands outside $L$. That worker then
takes $L$, changes BE at line 21 and wakes all the waiting workers, and
no wake-up is lost, since the wait releases $L$ atomically and all the
changes are made under $L$. By line 19 each state is closed at most once
(one can see that a closed state never re-enters OPEN), and on a finite
graph the selections are therefore finitely many. When OPEN and BE are
both empty, all the states reachable from $s_0$ have been closed, since
the successors of a closed state were all added at line 18, and no goal
was selected, so that none is reachable. $\square$

### The code and its engineering

The test of safety is the following excerpt; with real costs it allows
a relative error of $10^{-9}$ (`PASE.h:292`), since otherwise the
rounding might make the goal dependent on the states of equal $f$
(hence, with real costs, Lemma 8.17 and Theorem 8.18 hold up to this
tolerance, and exactly with integer costs):

```cpp
--8<-- "snippets/shared-memory-pase-safe.cpp"
```

and the selection, with the scan of at most `c_candidates` $= 64$
states (`PASE.h:281`) and the wait on a condition variable, is the
following:

```cpp
--8<-- "snippets/shared-memory-pase-select.cpp"
```

We now turn to the implementation, which we have kept close to the
pseudocode. OPEN is an ordered set (a B-tree), whose key with integer
costs is the packed word (8.4) followed by the address of the node
(`PASE.h:249`). Thus, the states of equal $f$ come by decreasing $g$ and
the scan of SAFE stops at the first state of OPEN with $f(t) \ge f(s)$
(`PASE.h:310`). Since the states of equal $f$ are not checked, as the
definition allows, a layer of equal $f$ costs a scan of BE per state
instead of a scan of the layer. Also, the scan of line 6 is bounded by
64 states, which bounds the time under $L$ (a state further down would
have to be checked against a longer and longer prefix of OPEN, and by
Proposition 8.19 the bound costs no progress. Finally, the hashes of the
successors are computed outside $L$, and under $L$ their places in
CLOSED are prefetched before any lookup (`PASE.h:377`), and the lock is
therefore held for about one cache miss instead of one per successor.

As the pseudocode shows, PA$^*$SE has a price, since a safe state costs
a scan of BE and of the states of OPEN before it, all under one lock.
PA$^*$SE may therefore pay off only when an expansion is much slower
than these scans (perhaps by orders of magnitude), as in motion planning
with collision checks, which is the setting of
[[6]](Z-bibliography.md#ref-6). To study this regime on any domain,
`SlowDomain< D >` (`SlowDomain.h`) wraps a domain and adds to each
expansion a given number of iterations of
$j \leftarrow 11 i \bmod 9999943$, $i \leftarrow j$, which keep a core
busy without touching the memory, as in the experiment of section 3.1.1
of [[5]](Z-bibliography.md#ref-5). The test `pase-thread` checks, with
1, 2, 3 and 8 workers, that PA$^*$SE finds the distances of a
breadth-first search on 21 states spread over the state space of the
8-puzzle, and that it solves a slowed grid optimally. The test also
checks that the counter `reexpanded`, which `PASE` increments when it
selects a state already expanded, is 0 on the slowed grid, i.e., that no
closed state is reopened (the code never reopens one, as Lemma 8.17
requires). Finally, it checks that
$\varepsilon_p = 1.5$ gives a cost within the bound of Theorem 8.18 on
the 15-puzzle, and that a grid split by a wall is unsolvable.

## 8.9 The searches compared {#sec-8-9}

The following table gathers, for each parallel search of the library,
what the workers share, where the duplicates are detected, how the end
is detected and in which regime the search may pay off. Its overheads
are those that one may expect from the scheme; only measurements on a given
machine and domain (e.g., by `AStarBench`,
[§8.10](08-shared-memory.md#sec-8-10)) can decide among the searches.

| Search | Shared | Duplicates | End | Regime |
|:--|:--|:--|:--|:--|
| SPA$^*$ | OPEN, CLOSED | in CLOSED | busy, under $L$ | slow expansions |
| MQA$^*$ | $2p$ queues, CLOSED | in CLOSED | counter $m$ | contended OPEN |
| farm | nothing | at the master | at the master | slow generation |
| PA$^*$SE | OPEN, BE, CLOSED | in CLOSED | safe goal | very slow expansions |
| HDA$^*$ | batches | at the owner | deficit counting | balanced owners |
| HDA$^*$, MPI | messages | at the owner | Mattern | cluster memory |
| Safe PBNF | abstract graph | in an nblock | no nblock left | small scopes |

In the table, the end of PA$^*$SE is a safe goal, or OPEN and BE both
empty; that of HDA$^*$ on MPI is the test of the four counters of
Mattern, and that of Safe PBNF is the state in which no nblock is free
and none is in use. The last column is the regime in which the search
may pay off: expansions slow with respect to the time under the shared
lock (SPA$^*$), a lock of OPEN that is the bottleneck (MQA$^*$), a
generation of successors much slower than their bookkeeping (the farm),
expansions much slower than the scans of the test of safety (PA$^*$SE),
a distribution that balances the load (HDA$^*$), more states than the
memory of one machine holds (HDA$^*$ on MPI), and an abstraction with
many nblocks of small scope (Safe PBNF).

HDA$^*$ is treated in [Chapter 9](09-hash-distributed.md#ch-9) and Safe
PBNF in [Chapter 11](11-safepbnf.md#ch-11). In the comparison of section
3.2 of [[5]](Z-bibliography.md#ref-5), on the 15-puzzle and on the
24-puzzle, HDA$^*$ with Zobrist hashing is faster than Safe PBNF on the
"hard" instances, while Safe PBNF is faster on the "easy" ones, where
the initial burst of HDA$^*$, i.e., the expansions made while the OPEN
lists of the workers are still nearly empty, weighs more. The outcome on
other domains and machines is not necessarily the same. Since the
shared-memory searches are templates over the executor and take the same
`ParallelOptions`, one can measure the same search on each runtime, and
all the searches on the same runtime.

## 8.10 Benchmarks {#sec-8-10}

The tables of this section are written by `manual/build/tables.py` from
the runs of `tools/bench/campaign.sh`, whose machine and compiler are
stated in the note of each table. Each run is a separate process,
started by `tools/AStarBench` with the search, the runtime and the
number of workers as options, and pinned by `numactl` to cores of one
NUMA node, with its memory on that node and one worker per core
(`PBEFS_PIN`). Each parallel run is made twice, in two processes, and
the gap between the two runs of a pair is the A/A noise reported in the
note. On the 15-puzzle the costs are integers and exact, and each run
reports the cost it found with its status; only the solved runs enter
the tables. The clock is the `Stopwatch` started at the beginning of
`solve()` and read at its end. It includes the allocation of the shared
structures, the creation of the threads (the benchmark gives no
executor, and each `solve()` builds one), the search and the rebuilding
of the path. It does not include the reading of the file of the
instances, the construction of the domain and of its tables, nor the
printing of the results.

**Table 8.1.** The centralized searches on the 15-puzzle (campaign B):
SPA$^*$ and MQA$^*$ with the C++ threads, and the farm on FastFlow, on
the 6 instances of Korf on which campaign B runs SPA$^*$ and MQA$^*$ (1,
3, 10, 14, 15 and 26; the farm also runs on 6 more, which the table
leaves out), with the Manhattan distance, for each number of cores of
the campaign. Each cell holds the seconds summed over the 6 instances,
each the median of its two runs, and in parentheses the speedup over the
sequential A$^*$ with the bucket queue on 1 core, summed over the same
instances. The farm runs its master and one worker less than the cores
(one worker on 1 core), and it therefore uses as many threads as the
others but on 1 core. An empty cell is a search that did not solve every
instance.

--8<-- "tables/shared-memory-tiles.md"

**Table 8.2.** Slow expansions (campaign D): the instance 79 of Korf, an
"easy" one, in two parts, one with 1000 and one with 10000 iterations of
the work of `SlowDomain` added to each expansion. Each part gives
SPA$^*$, MQA$^*$, the farm (on FastFlow, with one worker less than the
cores), PA$^*$SE and HDA$^*$ with Abstract Zobrist hashing on the C++
threads, for each number of cores, and the reference is the sequential
A$^*$ with the bucket queue on the same slowed domain on 1 core; the
columns are those of Table 8.1, on one instance.

--8<-- "tables/shared-memory-slow.md"

Both tables are placeholders until the campaign is run on an idle
machine. In the meantime, we remark that the experiment of Table 8.1
measures the regime of "fast" expansions, in which the time under the
shared lock (SPA$^*$), the queues (MQA$^*$) or the master (the farm) is
often a large part of an expansion. That of Table 8.2 measures the
regime of "slow" expansions, in which it is a small part and in which
PA$^*$SE is meant to work ("slow" meaning, say, a collision check per
successor), the sequential A$^*$ on the same slowed domain being the
reference.

## 8.11 Notes and references {#sec-8-11}

The overheads of [§8.2](08-shared-memory.md#sec-8-2) and their trade-off
are those of the survey of Fukunaga et al.
[[1]](Z-bibliography.md#ref-1), section 3.1, with the precise
definitions of Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5),
section 2.1, who also identified the band effect (load imbalance) and
the burst effect (the start of HDA$^*$) as the sources of the search
overhead of the decentralized searches (their section 3.1); Kishimoto,
Fukunaga and Botea [[2]](Z-bibliography.md#ref-2) measured the same
overheads for HDA$^*$ (their section 4.2). The stopping criterion of
Theorem 8.6 is the condition stated informally in section 3.3.1 of
[[1]](Z-bibliography.md#ref-1) for the decentralized searches. We have
written it once, for all the searches of the library, as a property of
the entries in flight; each search then only has to prove that its test
implies it. Counting the entries in transit is classical in distributed
termination detection (Dijkstra and Scholten
[[35]](Z-bibliography.md#ref-35), Mattern
[[11]](Z-bibliography.md#ref-11)).

SPA$^*$ is Algorithm 2 of [[1]](Z-bibliography.md#ref-1), the baseline
of that survey; Burns et al. [[3]](Z-bibliography.md#ref-3) found that
the shared OPEN list makes it slower than the sequential A$^*$ when the
expansions are fast, even with lock-free data structures, which is the
motivation of the decentralized searches. Our SPA$^*$ differs from
Algorithm 2 of the survey in its termination (busy counted under the
lock of OPEN, instead of a separate procedure of termination detection)
and in its CLOSED, whose duplicates take no lock. The MultiQueue is due
to Rihani, Sanders and Dementiev [[37]](Z-bibliography.md#ref-37);
Postnikova et al. [[38]](Z-bibliography.md#ref-38) studied it as the
priority scheduler of parallel algorithms on graphs. We are not aware of
a published A$^*$ on a MultiQueue with the termination of Lemma 8.12;
MQA$^*$ should therefore be read as an engineering choice of the
library, whose expansions are those of SPA$^*$ up to the order. Also,
the farm is a master-worker form of the centralized approach of section
3.2 of [[1]](Z-bibliography.md#ref-1), which goes back to Kumar, Ramesh
and Rao [[153]](Z-bibliography.md#ref-153), written on the farm of
FastFlow [[12]](Z-bibliography.md#ref-12) on that runtime. PA$^*$SE is
due to Phillips, Likhachev and Koenig [[6]](Z-bibliography.md#ref-6),
who proved its bounded suboptimality for one goal state and measured it
on motion planning with expensive collision checks. Condition (G) and
the bound of 64 on the scan are ours.

Further parallel schemes of the literature (e.g., the retracting
parallel A$^*$ for machines with little memory, the transposition-driven
scheduling of a distributed table of states, and A$^*$ on graphics
processors) are outside this chapter;
[Chapter 12](12-more-parallel-schemes.md#ch-12) collects the schemes
that the library does not implement.

## 8.12 Exercises {#sec-8-12}

**Exercise 8.1** (`ctest -R PBeFS.spa-thread`). Remove the counter busy
from Algorithm 8.1, so that a worker sets done as soon as it finds no
entry of $f < z$ in OPEN. On the graph of
[Figure 8.3](08-shared-memory.md#fig-8-3), whose optimal cost is 4,
give an execution with two workers in which $w_0$ is slow on $a$ while
$w_1$ expands $b$, $v$ and selects the goal $t$ with $g = 5$, and the
search stops with $z = 5$. Then check that the case `spa-thread`, which
compares the cost of SPA$^*$ with the known costs of two instances of
Korf with 1, 2, 3 and 8 workers, passes with the counter.

**Exercise 8.2** (`ctest -R PBeFS.spa-thread`). Prove that SPA$^*$ with
$p = 1$ makes the expansions of A$^*$ with the same OPEN list, and that
with $p$ workers and a consistent heuristic each re-expansion of a
node $v$ happens because a node $u$ on an optimal path to $v$ had not
been expanded with its optimal $g$ when $v$ was selected. Bound the
number of the expansions of $v$ by the number of its distinct values
of $g$.

**Exercise 8.3** (`ctest -R PBeFS.mqa-thread`). In Algorithm 8.2,
replace line 13 by "push the entries of F in OPEN; m ← m + |F| − 1"
(the addition after the pushes). Show an interleaving of two workers in
which a worker reads $m = 0$ while an entry is in OPEN, and explain why
the version of the library, with the addition first, never does
(Lemma 8.12). The case `mqa-thread` runs the library's version on 1 to
8 workers.

**Exercise 8.4** (`ctest -R PBeFS.mqa-thread`). Prove Lemma 8.14 and
show that the order it gives is that of the sequential A$^*$ with the
ties of $f$ broken towards the larger $g$. What goes wrong in (8.4) if
$f \ge 2^{32}$, or if $g < 0$? For which domains of
[Chapter 14](14-domains.md#ch-14) is the packed key used?

**Exercise 8.5** (`ctest -R PBeFS.farm-thread`). With tasks of $b$
nodes and $p$ workers, show that at most $p\,b$ nodes are out at the
workers at once, and that the farm with $p = 1$ and $b = 1$ expands
the nodes in the order of A$^*$. Explain why the master computes $g$
at line 8 of Algorithm 8.3 from the current $g$ of the parent rather
than from the $g$ it had when the task left, and why both choices would
be correct. The case `farm-thread` runs the farm with tasks of 1 and 4
nodes.

**Exercise 8.6** (`ctest -R PBeFS.pase-thread`). Build a small graph
with two goal states and a pairwise heuristic that satisfies (H) but
not (G) on which PA$^*$SE with $\varepsilon_p = 1$ returns a goal that
is not optimal. Then show that the vacuum world, the grids and the
sliding-tile puzzles of the library satisfy (G). The case `pase-thread`
checks the optimality on the 8-puzzle and on a grid and the bound of
Theorem 8.18 with $\varepsilon_p = 1.5$.

**Exercise 8.7** (`ctest -R PBeFS.pase-thread`). Prove that the bound of
64 on the scan of line 6 of Algorithm 8.4 changes neither the
correctness nor the progress of PA$^*$SE, and that it may change the
number of the states expanded in parallel. The case `pase-thread` checks
the optimal cost on a slowed grid with 1, 2, 3 and 8 workers.

**Exercise 8.8** (`ctest -R PBeFS.pase-thread`). Which problems may
parallelize well with the searches of this chapter? Relate the answer to
the ratio between the cost of an expansion and the time spent under the
shared lock (or at the master), and to the number of distinct values of
$f$ of a problem. Then use `SlowDomain`, as the case `pase-thread` does,
to make a domain of the library slow, and say from which work per
expansion one may expect SPA$^*$, the farm and PA$^*$SE to overtake the
sequential A$^*$ on a given machine.

**Exercise 8.9** (`ctest -R PBeFS.farm-fastflow`). Run the same search
on the three runtimes (the cases `spa-openmp`, `mqa-openmp` and
`farm-fastflow` do so for SPA$^*$, MQA$^*$ and the farm) and check that
the costs agree. Explain why the expansions may not, and which of the
counters of `SearchStats` are deterministic for SPA$^*$ with one worker.
