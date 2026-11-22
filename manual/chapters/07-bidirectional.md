# 7. Bidirectional and external-memory search {#ch-7}

[Source: [`include/PBeFS/astar/Bidirectional.h`][Bidirectional.h]
(`BidirectionalOptions`, `BiNode`, `BiQueue`, `BiFronts`, `ReverseDomain`,
`BlindDomain`), [`BiBestFirst.h`][BiBestFirst.h] (`BiBestFirst`, `bidijkstra`,
`mm`, `mme`, `fmm`, `bae`), [`NBS.h`][NBS.h] (`NBS`, `nbs`),
[`DVCBS.h`][DVCBS.h] (`DVCBS`, `dvcbs`), [`GBFHS.h`][GBFHS.h] (`GBFHS`,
`gbfhs`, `bibfs`), [`ParallelBidirectional.h`][ParallelBidirectional.h]
(`ParallelBidirectional`, `pnba`), [`PEMSearch.h`][PEMSearch.h] (`PEMSearch`,
`PEMOptions`), [`IDDAStar.h`][IDDAStar.h] (`IDDAStar`, `IDDOptions`,
`idd_astar`), [`include/PBeFS/astar/Domain.h`][Domain.h]
(`BidirectionalDomain`, `BackwardOperatorDomain`, `MinCostDomain`,
`for_each_parent`); the test `bidir` of
[`test/bidirectional.cpp`][bidirectional.cpp], the test `bidir-thread` of
[`test/pbidirectional.cpp`][pbidirectional.cpp], and the tests `pem-thread` and
`idd` of [`test/pem.cpp`][pem.cpp]]

We are concerned in this chapter with two ways of reducing what a search must
hold to find a path of least cost. First, one can search from both ends: a
*front-to-end bidirectional search* grows a forward search from the initial
state $s$ and a backward search from the goal $t$, along the reversed edges,
each guided by a heuristic toward the root of the other. It stops when the two
searches have met and a lower bound on the optimal cost $C^*$, drawn from their
open lists, proves that no cheaper path joins them. Since the work of a search
grows, typically, exponentially with the depth it reaches, two searches of half
the depth may cost far less than one of the whole depth (cf. the last row of
Table 1.1 in [§1.3](01-introduction.md#sec-1-3)). Here the difficulty is to
stop as soon as the two searches meet, without giving up optimality and without
letting either of them go beyond the middle; it took some decades to settle.
Second, one can keep the lists on the disk, as the delayed duplicate detection
of [§4.6](04-memory-bounded.md#sec-4-6) does for A$^*$. Here we treat the
parallel external-memory bidirectional searches of Siag et al.
[[57]](Z-bibliography.md#ref-57), which combine the two ideas, and the A$^*$
with *immediate* duplicate detection on the disk of Lin and Fukunaga
[[95]](Z-bibliography.md#ref-95), whose closed list is a file read at random
(which, as we discuss in section 7.10, suits a solid state drive).

The structure of the chapter is as follows. Section 7.1 states the interface,
i.e., the concepts that a domain has to satisfy to be searched backward, the
options of the searches and the two adapters of a domain. In section 7.2 we
describe the table of the nodes shared by the two directions, the meeting of
the two fronts and the two lemmas on which the proofs of the chapter rest;
section 7.3 gives the lower bounds on $C^*$ by which the searches stop.
Sections 7.4 to 7.7 give the searches that hold their lists in memory: those
that expand one node at a time (the bidirectional Dijkstra's algorithm, MM,
MMe, fMM and BAE$^*$), NBS with the must-expand graph, DVCBS with its dynamic
vertex cover, and GBFHS. Section 7.8 runs the two directions at once on two
workers (PNBA$^*$ and the parallel versions of the others), section 7.9 gives
the parallel external-memory searches of the framework PEM-BiHS, and section
7.10 the A$^*$ with immediate duplicate detection on the disk; sections 7.11 to
7.13 give the benchmark, the notes and references, and the exercises. We use,
without repeating them, the notation of [Chapter 1](01-introduction.md#ch-1),
the results on A$^*$ of [Chapter 2](02-best-first.md#ch-2), the lazy priority
queues and the tables of [Chapter 3](03-data-structures.md#ch-3), and the
buckets and files of the external search of
[§4.6](04-memory-bounded.md#sec-4-6).

## 7.1 The problem and the interface {#sec-7-1}

We consider the problem of [Chapter 1](01-introduction.md#ch-1) with a single
goal state $t$, for which a bidirectional search needs, besides a `Domain`, the
goal (from which the backward search starts), the predecessors of a state
(along which that search moves) and a heuristic toward $s$ (by which it is
guided). These are the requirements of the concept `BidirectionalDomain`
([`Domain.h:476`][Domain.h:476]): a `Domain` that also has `goal_state()`, the
one state for which `is_goal()` holds, `backward_heuristic( s )`, an estimate
of the cost of a path from the initial state to `s`, and `for_each_predecessor(
s , hb_s , f )`. This function calls `f( t , c , hb_t )` for each state `t`
with an edge $(t, s)$ of cost `c`, where `hb_t` is the backward heuristic of
`t`, computed from that of `s` (incrementally, if the domain can do it). The
refinement `BackwardOperatorDomain` ([`Domain.h:496`][Domain.h:496]) has,
moreover, `for_each_predecessor( s , hb , op , f )`, which skips the
predecessor reached by the inverse of `op`, the operator by which the backward
search generated `s`; this is the backward form of the skip by which an
`OperatorDomain` avoids the move back to the parent of a successor ([Lemma
1.8](01-introduction.md#lem-1-8)). The function template `for_each_parent`
([`Domain.h:773`][Domain.h:773]) is the backward counterpart of
`for_each_child` (Algorithm 1.1): it calls the form with the operator when the
domain has it, and the plain one otherwise. Finally, `MinCostDomain`
([`Domain.h:512`][Domain.h:512]) is a `Domain` with `min_edge_cost()`, a lower
bound $\varepsilon \ge 0$ on the cost of an edge, which several of the bounds
below use. We remark that these concepts are satisfied by the sliding-tile
puzzles (whose predecessors are their successors, since the moves are
reversible at the same cost), by the grids (whose predecessors are the reversed
moves, at the cost of the cell left) and by the explicit graphs, and on random
grids with weighted cells the test `bidir` checks that the edges seen backward
are those seen forward, with their costs.

<a id="def-7-1"></a>**Definition 7.1** (Notation of the bidirectional
searches). We write $d(u, v)$ for the cost of a cheapest path from $u$ to $v$,
$C^* = d(s, t)$, $F$ and $B$ for the forward and the backward direction, and
$\bar D$ for the direction other than $D$. For a state $n$, $g_F(n)$ is the
cost of the best path found from $s$ to $n$ and $g_B(n)$ that of the best path
found from $n$ to $t$, each $\infty$ until its direction reaches $n$, and we
write $d_F(n) = d(s, n)$ and $d_B(n) = d(n, t)$ for their optimal values. The
heuristic $h_F(n)$ estimates $d(n, t)$, $h_B(n)$ estimates $d(s, n)$, and $f_D
= g_D + h_D$; the heuristic $h_F$ is *admissible* if $h_F(n) \le d(n, t)$ and
*consistent* if $h_F(u) \le c(u, v) + h_F(v)$ for each edge $(u, v)$
(Definitions [2.1](02-best-first.md#def-2-1) and
[2.2](02-best-first.md#def-2-2)), and symmetrically $h_B$ is admissible if
$h_B(n) \le d(s, n)$ and consistent if $h_B(v) \le c(u, v) + h_B(u)$. Finally,
$\varepsilon \ge 0$ is a lower bound on the cost of each edge, and we call
*incumbent* the cost $U$ of the best path from $s$ to $t$ found so far
($\infty$ if none); this is the $z$ of Table 1.2, which we call $U$ as the
headers and the literature on bidirectional search do.

We provide two adapters that change a `BidirectionalDomain` without copying it.
First, `ReverseDomain< D >` ([`Bidirectional.h:640`][Bidirectional.h:640])
swaps the initial and the goal state: its successors are the predecessors of
`D`, its heuristic is the backward heuristic of `D`, and vice versa. Thus,
A$^*$ on it is the *reverse A$^*$* (rA$^*$) of Siag et al.
[[57]](Z-bibliography.md#ref-57), and a path found on it, reversed, is a path
of `D`; also, the operators are kept when `D` is both an `OperatorDomain` and a
`BackwardOperatorDomain`, and `min_edge_cost()` when `D` is a `MinCostDomain`.
Second, `BlindDomain< D >` ([`Bidirectional.h:709`][Bidirectional.h:709])
replaces both heuristics by 0: A$^*$ on it is Dijkstra's algorithm, and MM on
it is the brute-force bidirectional search MM$_0$ of Holte et al.
[[61]](Z-bibliography.md#ref-61). Both adapters hold a reference to `D` (which
must outlive them).

The in-memory searches take the options

```cpp
--8<-- "snippets/bidirectional-sig-options.cpp"
```

i.e., those of all the searches and the following ones; of the `SearchOptions`,
`node_limit` bounds the nodes stored by the two directions together,
`expected_nodes` sizes the table of the nodes at the start, and a `weight`
other than 1 is refused with `std::invalid_argument`. Among the others,
`epsilon` is the $\varepsilon$ of Definition 7.1: a negative value (the
default) takes `min_edge_cost()` of a `MinCostDomain`, and 0 for any other
domain. With `consistent` (true by default) one asserts that both heuristics
are consistent, which lets a search drop a node that the other direction has
closed ([Lemma 7.4](#lem-7-4)), and which BAE$^*$ and the parallel searches
require; it is an assumption on the domain, which the in-memory searches do not
check: a domain whose heuristics are only admissible must set it to false, or
the cost returned may be wrong. With `all_bounds`, a search stops at the
largest of the lower bounds of [§7.3](#sec-7-3) that hold under the heuristics,
instead of the bound of the algorithm, and `side` chooses the direction of the
next expansion, among the values of

```cpp
--8<-- "snippets/bidirectional-sig-sidechoice.cpp"
```

where `automatic` is the choice of each algorithm as published. In turn, the
directions are indices into the arrays of a node,

```cpp
--8<-- "snippets/bidirectional-sig-sides.cpp"
```

Finally, `fraction` is the $p$ of fMM($p$) and the share of the forward search
in the split of GBFHS, and `eager` selects the variant of NBS and DVCBS that
covers all the optimal paths ([§7.5](#sec-7-5)). Each search of sections 7.4 to
7.8 is a class template whose constructor takes the domain (which must outlive
the search) and the options, with a method `solve()` that returns a `Solution`
(the path, its cost, the status and the counters of
[§2.10](02-best-first.md#sec-2-10)); a function of one line, named after the
search, builds it and calls `solve()`, e.g., `mm( domain , options )`.

## 7.2 The two fronts {#sec-7-2}

We now describe the table of the nodes that the searches of sections 7.4 to 7.7
share, built by the class template `BiFronts`
([`Bidirectional.h:331`][Bidirectional.h:331]), since these searches differ
only in the order in which they expand its nodes and in the bound at which they
stop.

### The node, the table and the queues

A node of the table is a `BiNode`
([`Bidirectional.h:191`][Bidirectional.h:191]), which holds the state and, for
each direction $D$ (the index 0 for $F$ and 1 for $B$), $g_D$, $h_D$, the
parent on the path of $D$ and the operator by which $D$ generated the node,
with a byte of flags. For each direction, these flags record whether the node
is closed, whether it has ever been expanded (to count the reexpansions),
whether it is open, and whether it is *ready* in the sense of NBS and DVCBS. A
node is created by the first direction that reaches its state, with $g$ of the
other direction infinite; $h$ of the other direction is computed when that
direction reaches it (or at once, for the searches that need both heuristics of
each node, e.g., BAE$^*$). Therefore, the two directions look up their states
in one table (the open-addressing table of
[§3.6](03-data-structures.md#sec-3-6)), and the lookup by which a direction
detects its duplicates also shows whether the other direction has reached the
state: this is the *immediate solution detection* of Sturtevant and Chen
[[94]](Z-bibliography.md#ref-94), which costs no more than the duplicate
detection.

The open lists are lazy, as that of A$^*$ ([Lemma
3.1](03-data-structures.md#lem-3-1)): an open list holds entries $(n,
\mathit{key}, g)$, not the nodes. An entry is *stale* when its node has been
closed in that direction or reached again with a smaller $g$ since the entry
was pushed, or, with `consistent`, when its node has been closed in the other
direction. A stale entry is dropped when it reaches the top
([`Bidirectional.h:476`][Bidirectional.h:476]), and the entries that are not
stale are *fresh*. Such a list is the class template `BiQueue`
([`Bidirectional.h:239`][Bidirectional.h:239]), ordered by the key and then by
the larger $g$. It is the two-level bucket queue of
[§3.3](03-data-structures.md#sec-3-3) when the key and the cost are the same
integer type (the keys being nonnegative), and the $d$-ary heap of
[§3.2](03-data-structures.md#sec-3-2) otherwise, e.g., for the real keys of fMM
(neither of them decreases a key, and a node whose key changes is pushed
again).

### The expansion and the meeting

At the start, `BiFronts::start` creates the two roots and opens $s$ forward and
$t$ backward; when they are the same state, the incumbent is $U = 0$ at once.
Then `BiFronts::expand( d , n , push )`
([`Bidirectional.h:391`][Bidirectional.h:391]) closes $n$ in the direction $d$
and enumerates its successors (forward) or its predecessors (backward) by
`for_each_child` or `for_each_parent`, prefetching the place of each in the
table ([§3.8](03-data-structures.md#sec-3-8)), and then, for each state $m$ so
generated, with $g = g_d(n) + c$ and its heuristic $h$ in the direction $d$, it
does the following.

1. If $g \ge U$ or $g + h \ge U$, it drops $m$, since no path through it
   in this direction is cheaper than the incumbent (with an admissible
   heuristic).
2. Otherwise it looks $m$ up, creating its node if absent; if the node has
   $g_d(m) \le g$ it counts a duplicate and drops it, and otherwise it sets
   $g_d(m) = g$ and the parent and the operator of $d$, clears the flag
   closed of $d$ (a reopening, if it was set), and, if $g_{\bar d}(m)$ is
   finite and $g + g_{\bar d}(m) < U$, sets $U = g + g_{\bar d}(m)$ and
   records $m$ as the meeting node.
3. If `consistent` holds and $m$ is closed in the direction $\bar d$, it
   does not open $m$, since every path through it is known ([Lemma
   7.4](#lem-7-4)); otherwise it opens $m$ and calls `push`, by which the
   search inserts $m$ in its own lists.

When the search ends, the path is the chain of the forward parents from the
meeting node back to $s$, reversed, followed by the chain of the backward
parents from it to $t$ ([`Bidirectional.h:497`][Bidirectional.h:497]). Figure
7.1 shows the two fronts of the bidirectional Dijkstra's algorithm on a grid,
before and after they meet.

![Figure 7.1. The two fronts of a bidirectional search on a grid with 4
moves of unit cost, from $s$ and from $t$, $C^* = 6$; each cell carries its
$g$ in the direction that reached it. (a) After the expansion of the states
within distance 1 of each root: no state is reached by both, $U = \infty$,
and the bound $g_{\min F} + g_{\min B} + \varepsilon$ of the bidirectional
Dijkstra's algorithm is $2 + 2 + 1 = 5$. (b) After those within distance 2:
the state $m$ is reached by both, at $g_F(m) = g_B(m) = 3$, which makes $U =
6$, while the bound is now $3 + 3 + 1 = 7 \ge U$, and the search
stops.](../figures/bidirectional-fronts.svg){#fig-7-1}

### The invariants

We now state the facts on the table that the searches of the chapter use:
first, $U$ is always the cost of a path; second (and this is the main fact), an
optimal path always lies between the two fronts while $U > C^*$. This is the
argument of Holte et al. (section 4 of [[81]](Z-bibliography.md#ref-81)), which
we write for any search that detects the meetings at once and reopens a node
reached by a cheaper path. In what follows, a search *closes a node unexpanded*
when it pops a node whose $f$ in the direction of the pop is not smaller than
$U$, which all the searches do (e.g.,
[`BiBestFirst.h:342`][BiBestFirst.h:342]); a *step* of a search is an expansion
or such a closure.

<a id="lem-7-2"></a>**Lemma 7.2** (Values of $g$ and of $U$). At any time, each
finite $g_F(n)$ is the cost of a path from $s$ to $n$ and each finite $g_B(n)$
that of a path from $n$ to $t$, whence $g_D(n) \ge d_D(n)$ (cf. Definition
7.1). The incumbent $U$ is $\infty$ or the value $g_F(m) + g_B(m)$ at some time
of a node $m$, hence $U \ge C^*$, and the parent chains of the meeting node
give a path from $s$ to $t$ of cost at most $U$.

*Proof.* The value $g_D$ is set only to $0$ at the root of $D$ and to $g_D(n) +
c$ along an edge of $D$ from a node $n$, which by induction is the cost of a
path. $U$ is lowered only to $g_F(m) + g_B(m)$ for the node $m$ recorded with
it, and a later decrease of $g_F(m)$ or of $g_B(m)$ lowers $U$ again through
the same test (cf. step 2). Also, the parent of $m$ in $D$ is set together with
$g_D(m) = g_D(n) + c$, and $g_D(n)$ can only decrease afterwards; hence, as for
A$^*$ ([Lemma 2.10](02-best-first.md#lem-2-10)), the parent pointers of $D$
form a tree in which the cost of the chain of a node is at most its $g_D$.
$\square$

Thus, when a search can prove $U \le C^*$, one can return $U$ with the path of
the lemma as an optimal solution. Let now $P = (s = p_0, p_1, \dots, p_k = t)$
be an optimal path; we say that $p_l$ is *settled forward* if it is closed
forward with $g_F(p_l) = d_F(p_l)$, and *settled backward* if it is closed
backward with $g_B(p_l) = d_B(p_l)$. For the flag `consistent` we need the
following hypothesis, which [Lemma 7.4](#lem-7-4) proves for the searches of
the chapter.

(H) While $U > C^*$, each node of an optimal path from $s$ to $t$ that is
closed in a direction $D$ has its optimal $g_D$.

<a id="lem-7-3"></a>**Lemma 7.3** (An optimal path between the fronts). Let a
search expand by the procedure above and, if it uses the flag `consistent`, let
(H) hold. Then, between two steps, if $U > C^*$ there are indices $i < j$ such
that $p_i$ is open forward with $g_F(p_i) = d_F(p_i)$ and a fresh entry in the
forward open list, and $p_j$ is open backward with $g_B(p_j) = d_B(p_j)$ and a
fresh entry in the backward open list (in particular, both open lists have a
fresh entry).

*Proof.* We first remark that no node $p_l$ of $P$ ever has both $g_F(p_l) =
d_F(p_l)$ and $g_B(p_l) = d_B(p_l)$ while $U > C^*$: in fact, the later of the
two assignments finds the other $g$ finite and lowers $U$ to at most
$d_F(p_l) + d_B(p_l) = C^*$ (step 2, or the start if $s = t$), and $U$ never
increases.

Let $i$ be the least index such that $p_i$ is not settled forward; it exists,
since otherwise $t = p_k$ would have $g_F(t) = C^*$ and $g_B(t) = 0$. If $i =
0$, $p_0 = s$ has $g_F(s) = 0$; if $i > 0$, $p_{i-1}$ is settled forward; note
that it was closed with its present $g_F$, since a later decrease would have
reopened it, and it was not closed unexpanded, since $f_F(p_{i-1}) \le
d_F(p_{i-1}) + d_B(p_{i-1}) = C^* < U$ by admissibility ($U$ was larger than
$C^*$ then, as it is now). Hence it was expanded with $g_F(p_{i-1}) =
d_F(p_{i-1})$, and it generated $p_i$ with $g = d_F(p_i)$ and $g + h_F(p_i) \le
C^* < U$, which step 1 did not drop, and after that expansion $g_F(p_i) =
d_F(p_i)$, which by Lemma 7.2 cannot decrease. In both cases $p_i$ has its
optimal $g_F$, and it is not closed forward, since it is not settled; when
$g_F(p_i)$ took its present value, $p_i$ was either opened and pushed with it
or, under `consistent`, found closed backward; in the latter case it is settled
backward by (H), which by the first remark contradicts $U > C^*$. Also, the
entry pushed then is fresh: $p_i$ has not been closed forward since (it would
be settled forward), nor backward under `consistent` (it would be settled
backward by (H)). A symmetric argument gives the largest index $j$ such that
$p_j$ is not settled backward: $p_j$ is open backward with $g_B(p_j) =
d_B(p_j)$ and a fresh entry.

Finally, if $j \le i$, the node $p_j$ would have both optimal values of $g$:
$g_B(p_j)$ by the choice of $j$, and $g_F(p_j)$ because $p_j$ is settled
forward if $j < i$, and by the argument above if $j = i$, which by the first
remark contradicts $U > C^*$, and hence $i < j$. $\square$

We remark that Lemma 7.3 requires nothing of the order of the expansions, and
that without the flag `consistent` it holds with heuristics that are only
admissible, since a node reached by a cheaper path is then always reopened. The
flag is the "nipping" of the BS$^*$ of Kwa [[90]](Z-bibliography.md#ref-90): a
direction that reaches a node closed by the other one does not open it, which
is sound when the node has its optimal $g$ in the other direction, which the
order of the expansions guarantees under consistent heuristics.

<a id="lem-7-4"></a>**Lemma 7.4** (Closed with the optimal $g$). Let both
heuristics be consistent, and let each direction $D$ of a search pop its nodes
by the least key $\mathrm{pr}_D$, ties broken by the larger $g_D$, where
$\mathrm{pr}_D$ is nondecreasing along an edge of $D$ and strictly increasing
in $g_D$ at a fixed state. (a) Without the flag `consistent`, each node is
expanded in each direction with its optimal $g$, and none is reexpanded; (b)
with the flag, (H) holds. The keys $f$, $\max\{f, 2g\}$, $\max\{f, 2g +
\varepsilon\}$, $\max\{f, g/p\}$, $g$ and $b_D = f_D + g_D - h_{\bar D}$ have
these properties, and (a) and (b) hold also for the orders of NBS, DVCBS and
GBFHS.

*Proof.* The keys are monotone along an edge of cost $c$: $f$ by consistency
([Lemma 2.12](02-best-first.md#lem-2-12)), $g$ and $2g$ trivially, and $b_D$
because it changes by $2c + \Delta h_D - \Delta h_{\bar D} \ge 2c - c - c = 0$,
both differences of the heuristics being at least $-c$ along the edge by
consistency, and each of them is strictly increasing in $g$ at a fixed state.

(a) Let $n$ be the first node expanded in $D$ with $g_D(n) > d_D(n)$, $Q$ an
optimal path of $D$ from the root of $D$ to $n$, and $x$ the first node of $Q$
that is not closed in $D$ with its optimal $g_D$. Note that the nodes of $Q$
before $x$ were expanded, not closed unexpanded: if one of them, $q$, had been
closed unexpanded, then $f_D(q) \ge U$ at that time, and $f_D(n) > d_D(n) +
h_D(n) \ge d_D(q) + h_D(q) = f_D(q)$ by consistency along $Q$; since $U$ has
not increased, $n$ too would be closed unexpanded. Hence the predecessor of $x$
on $Q$ (if any) generated $x$ with $g_D(x) = d_D(x)$ and $f_D(x) \le d_D(n) +
h_D(n) < f_D(n) < U$, and thus $x$ is open in $D$ with a fresh entry (without
the flag, no node is dropped as closed by the other direction), and $x \ne n$.
Along $Q$, $\mathrm{pr}_D(x) \le \mathrm{pr}_D(n \text{ at } d_D(n)) <
\mathrm{pr}_D(n \text{ at } g_D(n))$, and therefore $x$ is popped before $n$, a
contradiction. Finally, a node expanded with its optimal $g$ is never reached
again with a smaller one, and it is not reexpanded.

(b) Let $p = p_l$ be the first node, over the run, of an optimal path $P$ from
$s$ to $t$ that is closed in a direction while $U > C^*$ with a $g$ larger than
its optimal one; say $B$, the forward case being symmetric. As in the proof of
Lemma 7.3, let $a > l$ be the largest index such that $p_a$ is not settled
backward, which exists since otherwise $p_{l+1}$ would have been expanded with
its optimal $g_B$ (it is not closed unexpanded, its $f_B$ being at most $C^* <
U$), and would have given $p_l$ its optimal $g_B$, which cannot then increase.
Thus $p_a$ has its optimal $g_B$ and is not closed backward. Suppose that it
was not opened backward, or that its entry is stale, because it is closed
forward; then, by the choice of $p$, it is closed forward with its optimal
$g_F$, and the first remark of the proof of Lemma 7.3 contradicts $U > C^*$.
Hence $p_a$ is open backward with a fresh entry and, as in (a),
$\mathrm{pr}_B(p_a) < \mathrm{pr}_B(p_l)$; it is therefore popped first, a
contradiction.

For NBS and DVCBS, a node is expanded from the ready nodes, by least $g_D$, and
each open node with $f_D$ below the bound ($C_{lb}$ or LB) is ready when a pair
or a cluster is chosen ([Lemma 7.10](#lem-7-10)). The node $x$ (or $p_a$) has
$f_D$ smaller than that of the node popped, which is at most the bound, and
hence it is ready; since its $g_D$ is smaller, it precedes the node popped (and
it cannot be generated by the cluster being expanded, whose children have a
$g_D$ at least that of the cluster). For GBFHS the nodes of a level are popped
by $f$, and $x$ (or $p_a$) has $f_D$ and $g_D$ smaller than those of the node
popped; hence it is within the limits of its level, and it is popped first.
$\square$

Hence, under the flag `consistent`, a direction that generates a node of an
optimal path closed by the other direction does not open it: its two values of
$g$ give the path through it, and the meeting is recorded at the generation
(step 2). We remark that (b) does not extend to the nodes off the optimal
paths: such a node may be expanded with a $g$ larger than its optimal one when
the cheaper path to it goes through a node that the flag did not open; this may
cost some search, but it does not affect the correctness. With this flag, the
test `bidir` runs all the searches on the 8-puzzle, the 15-puzzle, the grids
and the graphs whose heuristics are consistent, and without the flag on a graph
with an admissible heuristic that is not consistent (that of the test of the
reopenings of [§2.7](02-best-first.md#sec-2-7)), where all the searches but
BAE$^*$ return the optimal cost.

## 7.3 The lower bounds {#sec-7-3}

A front-to-end search stops when $U$ is not larger than a lower bound
$\mathrm{LB}$ on $C^*$, or when one of its open lists is empty; all the bounds
below are computed from the least keys of the open lists, and all of them rest
on [Lemma 7.3](#lem-7-3). They need not be lower bounds on $C^*$ at all times:
it suffices, for the stopping rule, that they are while $U > C^*$. We write
$f_{\min D}$, $g_{\min D}$, $b_{\min D}$ and $d_{\min D}$ for the least $f_D$,
$g_D$, $b_D$ and $\delta_D = g_D - h_{\bar D}$ over the fresh entries of the
open list of $D$.

<a id="thm-7-5"></a>**Theorem 7.5** (The lower bounds). Under the
hypotheses of [Lemma 7.3](#lem-7-3), whenever $U > C^*$ each of the
following is at most $C^*$:

1. $f_{\min F}$ and $f_{\min B}$, the bounds of A$^*$;
2. $g_{\min F} + g_{\min B} + \varepsilon$, the bound of the bidirectional
   Dijkstra's algorithm;
3. for $p \in (0, 1)$, the least key of fMM($p$),
   $C_p = \min\{\min_F \max\{f_F, g_F / p\},
   \min_B \max\{f_B, g_B / (1 - p)\}\}$, which for $p = 1/2$ is the least
   key $C$ of MM, and the least key of MMe,
   $\min_D \min \max\{f_D, 2 g_D + \varepsilon\}$;
4. if both heuristics are consistent, $(b_{\min F} + b_{\min B}) / 2$,
   the bound of BAE$^*$;
5. if both heuristics are consistent, $f_{\min F} + d_{\min B}$ and
   $d_{\min F} + f_{\min B}$.

With integer costs each of them may be rounded up.

*Proof.* Let $n = p_i$ and $m = p_j$ be the nodes of Lemma 7.3, whence $C^* =
g_F(n) + d(n, m) + g_B(m)$; note that $d(n, m) \ge \varepsilon$, since $i < j$.
(1) $f_F(n) = g_F(n) + h_F(n) \le d(s, n) + d(n, t) = C^*$ by admissibility,
and the same holds for $f_B(m)$. (2) $g_F(n) + g_B(m) + \varepsilon \le C^*$.
(3) Since $g_F(n) + g_B(m) \le C^*$, either $g_F(n) \le p\, C^*$ or $g_B(m) \le
(1 - p)\, C^*$; in the first case $\max\{f_F(n), g_F(n)/p\} \le C^*$, by (1),
and in the second the same holds for $m$ in $B$, and for MMe, if $g_F(n) \le
g_B(m)$, then $2 g_F(n) + \varepsilon \le g_F(n) + g_B(m) + \varepsilon \le
C^*$, and symmetrically otherwise. (4) By consistency along the subpath from
$n$ to $m$ ([Lemma 2.3](02-best-first.md#lem-2-3)), $h_F(n) \le d(n, m) +
h_F(m)$ and $h_B(m) \le d(n, m) + h_B(n)$, whence
$$b_F(n) + b_B(m) = 2 g_F(n) + 2 g_B(m) + (h_F(n) - h_F(m)) + (h_B(m) - h_B(n)) \le 2 (g_F(n) + d(n, m) + g_B(m)) = 2 C^*.$$
(5) $f_F(n) + \delta_B(m) = g_F(n) + g_B(m) + h_F(n) - h_F(m) \le g_F(n) +
g_B(m) + d(n, m) = C^*$, and symmetrically. Each least key over the fresh
entries is at most the key of $n$ or of $m$, which have fresh entries; note
that with integer costs $C^*$ is an integer, and a lower bound on it may be
rounded up. $\square$

<a id="thm-7-6"></a>**Theorem 7.6** (The stopping rule). Let a search satisfy
the hypotheses of [Lemma 7.3](#lem-7-3), and stop when an open list has no
fresh entry, or when $U < \infty$ and $U \le \mathrm{LB}$ for a value
$\mathrm{LB}$ that is at most $C^*$ whenever $U > C^*$. Then, when it stops
(and the node limit has not stopped it first), $U = C^*$ if $t$ is reachable
from $s$, and $U = \infty$ otherwise.

*Proof.* If $t$ is not reachable, $U = \infty$ by Lemma 7.2. Otherwise, suppose
that $U > C^*$ when the search stops; then, in fact, both open lists have a
fresh entry by Lemma 7.3, and $\mathrm{LB} \le C^* < U$; thus neither condition
holds, a contradiction, and hence $U \le C^*$, while $U \ge C^*$ by Lemma 7.2.
$\square$

By default each search uses the bound of its algorithm, and with `all_bounds`
the largest of those that hold under its assumptions: the bounds 1 and 2
always, the bound 3 for MMe (whose key is never smaller than that of MM) and,
with `consistent`, the bounds 4 and 5. Of course, all of them are dominated by
the bound of the pairs of NBS of [§7.5](#sec-7-5), which, however, needs a
search over the pairs of open nodes. We remark that the bounds of Theorem 7.5
are lower bounds on $C^*$ only while $U > C^*$: once $U = C^*$ they may exceed
$C^*$; this is harmless, since they then only stop the search.

## 7.4 One node at a time: MM, fMM, BAE$^*$ and the bidirectional Dijkstra's algorithm {#sec-7-4}

At each step, the class template `BiBestFirst< D , P >`
([`BiBestFirst.h`][BiBestFirst.h]) expands the node of least priority in the
direction chosen, with a priority among

```cpp
--8<-- "snippets/bidirectional-sig-bipriority.cpp"
```

i.e., for a node $n$ open in $D$, respectively,
$$g_D(n), \quad \max\{f_D(n), 2 g_D(n)\}, \quad \max\{f_D(n), 2 g_D(n) + \varepsilon\}, \quad \max\{f_D(n), g_D(n) / p_D\}, \quad b_D(n) = f_D(n) + g_D(n) - h_{\bar D}(n),$$
for the bidirectional Dijkstra's algorithm, MM of Holte et al.
[[61]](Z-bibliography.md#ref-61), [[81]](Z-bibliography.md#ref-81), its variant
MMe of Sharon et al. [[82]](Z-bibliography.md#ref-82), the fractional MM of
Shaham et al. [[83]](Z-bibliography.md#ref-83), with $p_F = p$ and $p_B = 1 -
p$, and BAE$^*$ of Sadhukhan [[88]](Z-bibliography.md#ref-88) as defined by
Alcázar et al. [[89]](Z-bibliography.md#ref-89). The priority of BAE$^*$ adds
to $f$ the error $g_D - h_{\bar D}$ of the heuristic of the other direction at
the node. Its declarations, without the bodies, are

```cpp
template< BidirectionalDomain D , BiPriority P >
class BiBestFirst {
 public:
 using Key = std::conditional_t< P == BiPriority::fmm , double , Cost >;
 explicit BiBestFirst( const D & domain ,
                       const BidirectionalOptions & options = {} );
 [[nodiscard]] SolutionT solve( void ) const;
 };

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > mm(
 const D & domain , const BidirectionalOptions & options = {} );
```

and the same for `mme`, `fmm` (whose $p$ is `options.fraction`), `bae` and
`bidijkstra`. In particular, `bidijkstra` is `BiBestFirst< BlindDomain< D > ,
BiPriority::dijkstra >`, i.e., the MM$_0$ of [[61]](Z-bibliography.md#ref-61),
which ignores the heuristics. The constructor throws `std::invalid_argument`
for BAE$^*$ when `consistent` is false, and for fMM when the fraction is not in
$(0, 1)$; the key is real for fMM, whose $g / p$ is not an integer, and of the
cost type otherwise.

**Algorithm 7.1** (The searches one node at a time, the pseudocode of the
header [`BiBestFirst.h:29-40`][BiBestFirst.h:29-40]).

```text
 1  procedure BiBestFirst(s, t)
 2     U ← ∞; OPEN_F ← {s}; OPEN_B ← {t}
 3     while OPEN_F ≠ ∅ and OPEN_B ≠ ∅ do
 4        if U ≤ LB then return U
 5        D ← the direction chosen
 6        n ← the node of least pr_D in OPEN_D, of largest g_D among those
 7        if f_D(n) < U then expand n in the direction D, and
             U ← min(U, g_F(m) + g_B(m)) for each m generated
 8     return U
```

Here OPEN$_D$ is the open list of $D$ by $\mathrm{pr}_D$, whose stale entries
are dropped at the top, and the expansion of line 7 is that of [§7.2](#sec-7-2)
(we write $g_F(m) + g_B(m)$ where the header writes $g_D(m) + g_{\bar D}(m)$);
a node of $f_D(n) \ge U$ is closed unexpanded. The bound of line 4 is
$$g_{\min F} + g_{\min B} + \varepsilon \ \text{(Dijkstra)}, \qquad \max\{C, f_{\min F}, f_{\min B}, g_{\min F} + g_{\min B} + \varepsilon\} \ \text{(MM, MMe, fMM)}, \qquad \lceil (b_{\min F} + b_{\min B}) / 2 \rceil \ \text{(BAE}^*\text{)},$$
where $C$ is the least priority of the two open lists (section 4 of
[[81]](Z-bibliography.md#ref-81)), and the ceiling applies to integer costs
only, and with `all_bounds` it is the largest of the bounds of Theorem 7.5 that
hold. For `SideChoice::automatic`, the direction of line 5 is, for MM, MMe, fMM
and Dijkstra's algorithm, the one of least priority (the forward one on a tie,
unless the top of the backward list has a smaller $g$), which is
`SideChoice::priority`, and for BAE$^*$ it is the alternation of the two, as in
its definition and in PEM-BAE$^*$ [[57]](Z-bibliography.md#ref-57). The choices
`alternate` and `cardinality` (the direction with fewer open nodes, the
criterion of Pohl [[91]](Z-bibliography.md#ref-91)) are available to all the
priorities.

<a id="thm-7-7"></a>**Theorem 7.7** (Correctness, and the meeting in the
middle). Algorithm 7.1 returns $C^*$ (or reports $t$ unreachable) for each
priority, with admissible heuristics and `consistent` false, or with consistent
heuristics; for BAE$^*$ the heuristics must be consistent. Moreover, with the
choice of the direction by priority, MM expands no node with $g_D(n) > C^*/2$,
MMe none with $2 g_D(n) + \varepsilon > C^*$, and fMM($p$) none with $g_F(n) >
p\, C^*$ or $g_B(n) > (1 - p)\, C^*$.

*Proof.* The priorities of the five searches satisfy [Lemma 7.4](#lem-7-4);
hence (H) holds under the flag `consistent`, and with it Lemma 7.3. Each bound
is one of those of Theorem 7.5 (the bound $C$ of MM, MMe and fMM is the bound
3, and that of BAE$^*$ the bound 4, which needs consistency), and Theorem 7.6
gives the first claim. For the second, let $n$ be expanded in $D$: if $U > C^*$
at that time, with the direction of least priority one has $\mathrm{pr}_D(n) =
C \le C^*$ by Theorem 7.5, and if $U \le C^*$, then $U = C^*$, and since the
search did not stop at line 4, $C \le \mathrm{LB} < U = C^*$. In both cases
$\mathrm{pr}_D(n) \le C^*$, i.e., $2 g_D(n) \le C^*$ for MM, $2 g_D(n) +
\varepsilon \le C^*$ for MMe, and $g_D(n) \le p_D\, C^*$ for fMM. $\square$

The meeting in the middle depends on the choice of the direction: with
`alternate` or `cardinality`, the node popped in the direction chosen may have
a priority larger than $C$, and then larger than $C^*$; thus the second claim
of the theorem holds only for `SideChoice::priority` (the default of these
searches). We check it in the test `bidir`, on 10 states of the 8-puzzle, for
MM and for fMM with $p = 0.25$ and $p = 0.8$. The test records the states that
each direction expands, and it compares their distances from $s$ and to $t$,
computed by a breadth-first search, with $p\, C^*$ and $(1 - p)\, C^*$; with
the same distances, it also checks that MM, NBS, DVCBS and GBFHS expand no
state whose optimal $g$ in the direction plus its heuristic exceeds $C^*$.

<a id="prop-7-8"></a>**Proposition 7.8** (Termination). On a finite graph
with nonnegative costs, Algorithm 7.1 terminates, and so do the searches of
sections 7.5 to 7.7.

*Proof.* Each step pops an entry, and an entry is pushed only after a strict
decrease of some $g_D(n)$, which is then the cost of a path (Lemma 7.2). If
$\delta$ is the least positive cost of an edge of the finite graph, a path of
cost at most a given value $K$ has at most $K / \delta$ edges of positive cost,
and therefore the distinct costs of such paths are finitely many. Hence each
$g_D(n)$ decreases finitely many times, the entries are finitely many, and so
are the steps, as in [Theorem 2.11](02-best-first.md#thm-2-11). In the
procedure that chooses the pair of NBS (or the cluster of DVCBS), each raise of
the bound is followed by a return or by the move of an entry from waiting to
ready; indeed, the new bound is the least of the quantities that the two tests
compare with it, and an entry is moved at most once. Finally, GBFHS raises its
bound by one unit at a time, and when the bound and the limits exceed the
largest $f$ and $g$ of the finitely many paths without repeated states, all the
open nodes are expandable. $\square$

### The implementation

For each direction, the search keeps the open list by priority and, for the
bounds, up to five further lazy queues: by $f$, by $g$, by the priority of MMe,
by $b$ and by $\delta = g - h_{\bar D}$. Each of them is kept only when the
bound of the algorithm (or, with `all_bounds`, one that holds) needs it; e.g.,
MM keeps the queues by $f$ and by $g$, and the bidirectional Dijkstra's
algorithm none. A node is pushed in each queue kept when its $g$ changes, and
the least key of a queue is read after dropping its stale entries. Hence each
bound costs as much as the open list, i.e., an amortized $O(1)$ per entry on
the bucket queue ([Proposition 3.5](03-data-structures.md#prop-3-5)) and
$O(\log n)$ on the heap. Both heuristics of a node are computed at its creation
only for BAE$^*$ and for the bounds that need $h_{\bar D}$; the other searches
compute $h_{\bar D}$ only when $\bar D$ reaches the node.

## 7.5 Expanding pairs: NBS and the must-expand graph {#sec-7-5}

For a node $u$ open forward and a node $v$ open backward, $u \ne v$, the value
$$\mathrm{lb}(u, v) = \max\{f_F(u),\ f_B(v),\ g_F(u) + g_B(v) + \varepsilon\}$$
is a lower bound on the cost of a path from $s$ through $u$ and then $v$ to $t$
(the bounds 1 and 2 of Theorem 7.5 for one pair). The Near-Optimal
Bidirectional Search (NBS) of Chen, Holte, Zilles and Sturtevant
[[85]](Z-bibliography.md#ref-85) expands *both* nodes of a pair of least lb,
and it is the only search of the chapter with a guarantee on the number of its
expansions relative to all the other front-to-end searches, which we now state.

<a id="def-7-9"></a>**Definition 7.9** (Must-expand graph). Let both heuristics
be consistent. The *must-expand graph* $G_{MX}$ is the bipartite graph with a
forward copy and a backward copy of each state, in which the forward copy of
$u$ and the backward copy of $v$ are joined if
$$\max\{d_F(u) + h_F(u),\ d_B(v) + h_B(v),\ d_F(u) + d_B(v) + \varepsilon\} < C^*,$$
i.e., if the pair has $\mathrm{lb} < C^*$ when both nodes have their optimal
$g$. We write $\mathrm{VC}$ for the size of a minimum vertex cover of $G_{MX}$.

Theorem 1 of Eckerle et al. [[84]](Z-bibliography.md#ref-84), as reported by
the headers and by Chen et al. [[85]](Z-bibliography.md#ref-85), concerns the
admissible front-to-end bidirectional searches, i.e., those that know only the
values of $g$, $h_F$ and $h_B$ of their nodes and $\varepsilon$. Such a search
must expand $u$ forward or $v$ backward for each edge $(u, v)$ of $G_{MX}$,
since otherwise a path of cost $\mathrm{lb} < C^*$ joining them could exist on
an instance indistinguishable from the given one, and hence one can conclude
that VC is a lower bound on the expansions of all such searches.

**Algorithm 7.2** (The choice of the pair of NBS, the pseudocode of the
header [`NBS.h:35-45`][NBS.h:35-45]).

```text
 1  procedure PrepareBest
 2     move from waiting_D to ready_D every node with f_D < C_lb
 3     loop
 4        if ready_D ∪ waiting_D = ∅ for some D then return none
 5        if ready_F.g + ready_B.g + ε ≤ C_lb then return the two tops
 6        if waiting_D.f ≤ C_lb for some D then
             move the top of waiting_D to ready_D (for each such D)
 7        else C_lb ← min(waiting_F.f, waiting_B.f,
                          ready_F.g + ready_B.g + ε)
```

We find the pair as in section 7 of [[85]](Z-bibliography.md#ref-85), with a
lower bound $C_{lb}$ on the least lb, raised until a pair reaches it: in each
direction $D$ there is a queue waiting$_D$ of the open nodes by $f_D$, the
larger $g_D$ first on a tie, and a queue ready$_D$ of those with $f_D \le
C_{lb}$, by $g_D$ and then by the smaller $f_D$. In the pseudocode, `ready_D.g`
is the least $g$ of ready$_D$ and `waiting_D.f` the least $f$ of waiting$_D$,
both over the fresh entries. At each step the search calls PrepareBest, and it
stops if PrepareBest returns none or if $U \le C_{lb}$. Otherwise it expands
the top of ready$_F$ and then that of ready$_B$, each closed unexpanded if its
$f$ is not below $U$ ([`NBS.h:234-246`][NBS.h:234-246]). With `eager`, line 2
moves the nodes with $f_D \le C_{lb}$, and line 6 is then never taken. This is
the NBS$_A$ of Shperberg et al. [[86]](Z-bibliography.md#ref-86), whose pairs
cover those with $\mathrm{lb} \le C^*$, as needed to find all the optimal
paths. With $\varepsilon = 0$ the bound is that of Chen et al. without
$\varepsilon$ (the base case of Shperberg et al.). Its declarations are those
of `BiBestFirst` without the parameter `P`, i.e., `NBS< D >( domain , options
)` with `solve()`, and `nbs( domain , options )`.

<a id="lem-7-10"></a>**Lemma 7.10** (The bound of the pairs). (a) Within a call
of PrepareBest, after line 2, each node of waiting$_D$ has $f_D \ge C_{lb}$
(with `eager`, $f_D > C_{lb}$); hence each open node with $f_D < C_{lb}$ is
ready, and $C_{lb}$ strictly increases at each execution of line 7. (b) Each
time $C_{lb}$ is set at line 7, it is at most the least lb of the pairs of open
nodes, and hence at most $C^*$ if $U > C^*$. (c) A pair returned at line 5 has
$\mathrm{lb}(u, v) \le C_{lb}$.

*Proof.* (a) Line 2 leaves in waiting$_D$ the nodes with $f_D \ge C_{lb}$, and
line 7 is reached only when the tops of waiting have $f > C_{lb}$ (line 6) and
the ready tops have $g_F + g_B + \varepsilon > C_{lb}$ (line 5); hence the new
value is larger, and it is again at most the least $f$ of waiting (no node
enters waiting within the call). (b) Let $u$ be open forward and $v$ open
backward, each with a fresh entry in ready or in waiting. If $u$ or $v$ is
waiting, $\mathrm{lb}(u, v)$ is at least its $f$, which is at least the least
$f$ of its waiting queue, and if both are ready, $\mathrm{lb}(u, v) \ge
g_F(u) + g_B(v) + \varepsilon$, which is at least the sum of the least values
of $g$ of the two ready queues plus $\varepsilon$. The new $C_{lb}$ is the
minimum of these quantities. If $U > C^*$, the pair $(p_i, p_j)$ of [Lemma
7.3](#lem-7-3) has $\mathrm{lb} \le C^*$ (Theorem 7.5, bounds 1 and 2), whence
$C_{lb} \le C^*$. (c) A node enters ready$_D$ with $f_D \le C_{lb}$, $C_{lb}$
only grows, and the $f$ of a ready node does not change while its entry is
fresh; the condition of line 5 bounds the third term of lb. $\square$

<a id="thm-7-11"></a>**Theorem 7.11** (Correctness and near-optimality of NBS).
NBS returns $C^*$ (or reports $t$ unreachable) under the hypotheses of Theorem
7.7. Let both heuristics be consistent and no node be reexpanded (which holds
with `consistent` false, by [Lemma 7.4](#lem-7-4)(a)). Then the pairs that NBS
expands with $\mathrm{lb} < C^*$ form a matching of $G_{MX}$, and it expands at
most $2\, \mathrm{VC}$ nodes in them; the other pairs expanded have
$\mathrm{lb} = C^*$.

*Proof.* NBS satisfies Lemma 7.4, and therefore (H) and Lemma 7.3 hold. By
Lemma 7.10(b), the value of $C_{lb}$ at the test $U \le C_{lb}$ is at most
$C^*$ whenever $U > C^*$, since $U$ only decreases after $C_{lb}$ is set;
Theorem 7.6 gives the correctness. Let now the heuristics be consistent, and
$(u, v)$ be a pair expanded: since $g_F(u) \ge d_F(u)$ and $g_B(v) \ge d_B(v)$,
the quantity of Definition 7.9 is at most $\mathrm{lb}(u, v)$, and $(u, v)$ is
an edge of $G_{MX}$ if $\mathrm{lb}(u, v) < C^*$. If no node is reexpanded,
each state is expanded at most once in each direction, and two pairs expanded
share no forward copy and no backward copy: they form a matching $M$. A vertex
cover of $G_{MX}$ contains an end of each edge of $M$, and the edges of $M$ are
disjoint; hence $|M| \le \mathrm{VC}$, and the pairs of $M$ make at most $2 |M|
\le 2\, \mathrm{VC}$ expansions. Finally, by Lemma 7.10(c), a pair is expanded
only if $U > C_{lb} \ge \mathrm{lb}(u, v)$; if $U = C^*$ this gives
$\mathrm{lb} < C^*$, and if $U > C^*$ it gives $\mathrm{lb} \le C_{lb} \le
C^*$. $\square$

This is theorem 3 of [[85]](Z-bibliography.md#ref-85), and their theorem 5
shows that no admissible front-to-end search can guarantee less than $2\,
\mathrm{VC}$ in the worst case. With the flag `consistent`, a node off the
optimal paths may be expanded with a $g$ larger than its optimal one (see the
remark after Lemma 7.4) and later reexpanded; such cases are reported by the
counter `reexpanded` of the `Solution`, and the bound of the theorem holds for
the runs in which it is 0. We also remark that the two expansions of a step are
those of the two tops returned by PrepareBest. Indeed, the forward expansion
cannot make the backward top stale, since it writes only the forward fields and
$u \ne v$ (if $u = v$, the node would have both values of $g$, and $U \le
g_F(u) + g_B(u) \le C_{lb}$ would have stopped the search).

## 7.6 Expanding clusters: DVCBS {#sec-7-6}

NBS expands both ends of one edge of $G_{MX}$ at a time, which is optimal only
up to the factor 2. Instead, the Dynamic Vertex Cover Bidirectional Search
(DVCBS) of Shperberg, Felner, Sturtevant, Shimony and Hayoun
[[86]](Z-bibliography.md#ref-86) expands, at each step, one vertex of a minimum
vertex cover of the part of the must-expand graph that it knows. The class
template `DVCBS< D >` ([`DVCBS.h`][DVCBS.h]), with the function `dvcbs( domain
, options )`, keeps the high level of NBS: a lower bound LB on $C^*$, raised
until $U \le \mathrm{LB}$, the queues waiting$_D$ by $f_D$, and the ready
nodes, i.e., those with $f_D < \mathrm{LB}$ ($\le$ with `eager`). The nodes
move from waiting$_D$ to ready$_D$ one value of $f_D$ at a time, and the ready
nodes are grouped in *clusters*, one per value of $g_D$. These clusters are the
vertices of the *dynamic must-expand graph* DGMX, in which a forward cluster of
$g_F = a$ and a backward one of $g_B = b$ are joined if $a + b + \varepsilon
\le \mathrm{LB}$; then each pair of their nodes has $\mathrm{lb} \le
\mathrm{LB}$. DVCBS expands one cluster of a minimum vertex cover of DGMX,
weighted by the sizes of the clusters (algorithm 3 of
[[86]](Z-bibliography.md#ref-86)). The header has no pseudocode of its own: the
clusters are chosen as in Algorithm 7.2, with the moves of whole values of $f$
and the test of line 5 on the least $g$ of the clusters
([`DVCBS.h:235-272`][DVCBS.h:235-272]), and the cover is computed as follows.

Let $a_0 < a_1 < \dots < a_{m-1}$ be the values of $g$ of the forward clusters
with an edge, and $b_0 < \dots < b_{n-1}$ those of the backward ones, with
weights $w(a_i)$ and $w(b_j)$, and we write $N(a_i) = \{b_j : a_i + b_j +
\varepsilon \le \mathrm{LB}\}$, a prefix of the backward clusters.

<a id="lem-7-12"></a>**Lemma 7.12** (The covers of DGMX). For $0 \le i \le m$,
the set $K_i = \{a_0, \dots, a_{i-1}\} \cup N(a_i)$ (with $N(a_m) = \emptyset$)
is a vertex cover of DGMX, and each vertex cover $K$ contains some $K_i$. Hence
a minimum cover is a $K_i$ of least weight $\sum_{l < i} w(a_l) + \sum_{b \in
N(a_i)} w(b)$.

*Proof.* The neighborhoods are nested: $N(a_{i'}) \subseteq N(a_i)$ for $i' >
i$, since $b \le \mathrm{LB} - \varepsilon - a_{i'} \le \mathrm{LB} -
\varepsilon - a_i$. An edge $(a_l, b)$ is covered by $K_i$ through $a_l$ if $l
< i$, and through $b \in N(a_l) \subseteq N(a_i)$ if $l \ge i$; hence $K_i$ is
a cover. Given a cover $K$, let $i$ be the least index with $a_i \notin K$ ($i
= m$ if there is none). Then $a_0, \dots, a_{i-1} \in K$, and each edge of
$a_i$ is covered by its backward end, i.e., $N(a_i) \subseteq K$, and thus $K_i
\subseteq K$, and the weight of $K_i$ is at most that of $K$. $\square$

This structure of the covers was shown by Shaham et al.
[[83]](Z-bibliography.md#ref-83) for the must-expand graph, and it makes the
minimum cover a pass over $i = 0, \dots, m$, with a pointer to the end of
$N(a_i)$ that only moves down and the prefix sums of the backward weights, in
time linear in the number of clusters ([`DVCBS.h:199-232`][DVCBS.h:199-232]).
Among the minimum covers, the cluster expanded is the forward $a_0$ or the
backward $b_0$ that lies in one of them: the one with fewer nodes, and the
backward one on a tie, which is the variant reported by Shperberg et al.
(section 5.1 of [[86]](Z-bibliography.md#ref-86)). Its nodes are expanded in
the order in which they became ready, each closed unexpanded if its $f$ is not
below $U$, and the expansion of the cluster stops as soon as $U \le
\mathrm{LB}$. Figure 7.2 shows the clusters, the edges of DGMX and a minimum
cover.

![Figure 7.2. The pairs of NBS and the clusters of DVCBS. The forward
clusters of $g_F = a$ (rows) and the backward ones of $g_B = b$ (columns)
carry their sizes, and each cell the bound $a + b + \varepsilon$ of the
pairs of its two clusters, with $\varepsilon = 1$; with $\mathrm{LB} = 5$
the cells within it (darker) are the edges of DGMX, a staircase since
$N(a_i)$ is a prefix that shrinks as $i$ grows. Of the covers $K_0, \dots,
K_4$ of Lemma 7.12, of weights 16, 18, 13, 14 and 20, the least is $K_2 =
\{a_0, a_1, b_0, b_1, b_2\}$; both $a_0$ and $b_0$ lie in it, and DVCBS
expands the smaller, $b_0$.](../figures/bidirectional-pairs.svg){#fig-7-2}

<a id="prop-7-13"></a>**Proposition 7.13** (Correctness of DVCBS). DVCBS
returns $C^*$ (or reports $t$ unreachable) under the hypotheses of Theorem
7.7.

*Proof.* One can check that Lemma 7.10 holds for LB with the clusters in place
of the ready queues, with the same proof, since the moves of a whole value of
$f$ keep the nodes of waiting at $f \ge \mathrm{LB}$. Also, the cluster
expanded is a cluster of least $g$ of its direction, and hence Lemma 7.4 holds,
and with it Lemma 7.3; since the test $U \le \mathrm{LB}$ is made after each
node of a cluster, Theorem 7.6 applies. $\square$

DVCBS has no bound of the kind of Theorem 7.11 on its expansions (section 5.2
of [[86]](Z-bibliography.md#ref-86)), but it expanded fewer nodes than NBS on
all the domains of their table 1. As for NBS, the test `bidir` runs it with
each variant of the options (all the bounds, eager, $\varepsilon = 0$, and the
choices of the direction, which these two searches ignore), and it checks that
DVCBS never expands a node of $f$ above $C^*$.

### The implementation

We keep a cluster as a vector of nodes in the order in which they became ready,
with a count of its live nodes, and the clusters of a direction are a map by
$g$ ([`DVCBS.h:123-134`][DVCBS.h:123-134]). A ready node whose $g$ decreases
leaves its cluster (the count is decremented, and the cluster erased when it
reaches 0) and goes back to waiting. Under `consistent`, a node expanded in $D$
also leaves its cluster of $\bar D$. A node may also die in its cluster by a
closure in the other direction, which is not counted as a leaving; for this
reason, the cluster expanded is erased as a whole at the end, and the weights
used by the cover are the live counts, i.e., the sizes of the clusters up to
these dead nodes.

## 7.7 GBFHS and the bidirectional breadth-first search {#sec-7-7}

The Generalized Breadth-First Heuristic Search (GBFHS) of Barley et al.
[[87]](Z-bibliography.md#ref-87) works with integer costs, and it raises a
bound $f_{Lim}$ on $C^*$ one unit at a time. At each value it splits $f_{Lim} -
\varepsilon + 1$ into two limits $g_{Lim F} + g_{Lim B}$ on the depths of the
two searches, and in each direction $D$ it expands all the open nodes with $f_D
\le f_{Lim}$ and $g_D < g_{Lim D}$.

**Algorithm 7.3** (GBFHS, the pseudocode of the header
[`GBFHS.h:16-27`][GBFHS.h:16-27]).

```text
 1  procedure GBFHS(s, t)
 2     U ← ∞; f_Lim ← max(h_F(s), h_B(t), ε); g_LimF ← g_LimB ← 0
 3     while OPEN_F ≠ ∅ and OPEN_B ≠ ∅ do
 4        if U ≤ f_Lim then return U
 5        raise g_LimF and g_LimB so that g_LimF + g_LimB = f_Lim − ε + 1
 6        while some n in OPEN_D has f_D(n) ≤ f_Lim and g_D(n) < g_LimD do
 7           expand the one of least f_D, of largest g_D among those,
             of the direction of the least such f
 8           if U ≤ f_Lim then return U
 9        f_Lim ← f_Lim + 1
```

The class template `GBFHS< D >` requires `std::is_integral_v< cost_t< D > >`,
as do `gbfhs( domain , options )` and `bibfs( domain , options )`, where
`bibfs` is GBFHS on the `BlindDomain` with each unit of depth given to the
direction with fewer open nodes, i.e., the bidirectional breadth-first search
with the criterion of Pohl [[91]](Z-bibliography.md#ref-91). The split of line
5 is chosen by `side`. With `alternate`, and also with `automatic` (the
default), each new unit goes to the smaller limit (the backward one on a tie)
after a first split in halves; with `cardinality` it goes to the direction with
fewer open nodes, and with `priority`, as in fMM($p$), the forward limit gets
the share `fraction` of the sum, $\lfloor p\,(f_{Lim} - \varepsilon + 1)
\rfloor$, clamped so that no limit decreases
([`GBFHS.h:177-206`][GBFHS.h:177-206]).

<a id="lem-7-14"></a>**Lemma 7.14** (A level finds its paths). Let the costs be
integers, $\varepsilon$ a lower bound on them, and the hypotheses of Lemma 7.3
hold. When the loop of lines 6 to 8 ends at the level $f_{Lim}$, each path $P$
from $s$ to $t$ of cost $C \le f_{Lim}$ has $U \le C$.

*Proof.* Let $P = (p_0, \dots, p_k)$, suppose $U > C$, and let $a$ be the least
index such that $p_a$ is not closed forward with $g_F(p_a)$ at most the cost of
the prefix of $P$ up to $p_a$, and $b$ the largest index such that $p_b$ is not
closed backward with $g_B(p_b)$ at most the cost of the suffix of $P$ from
$p_b$. We argue as in the proof of Lemma 7.3, with the costs of the prefixes in
place of the optimal values of $g$, and $C < U$ in place of $C^* < U$: thus
$p_a$ has $g_F(p_a)$ at most the cost of its prefix, it is not closed forward,
and it is open forward with a fresh entry; symmetrically for $p_b$. If $a \ge
b$, the node $p_b$ has both values of $g$ within the costs of its prefix and of
its suffix, whose sum, at most $C$, was offered to $U$; hence $a < b$. Then
$f_F(p_a) \le g_F(p_a) + d(p_a, t) \le C \le f_{Lim}$ by admissibility, and
since the level ended, $g_F(p_a) \ge g_{Lim F}$; symmetrically, $g_B(p_b) \ge
g_{Lim B}$. Since $P$ has at least one edge between $p_a$ and $p_b$, one has
$$C \ge g_F(p_a) + \varepsilon + g_B(p_b) \ge g_{Lim F} + g_{Lim B} + \varepsilon = f_{Lim} + 1,$$
a contradiction. $\square$

<a id="thm-7-15"></a>**Theorem 7.15** (Correctness of GBFHS). Under the
hypotheses of Lemma 7.14, GBFHS returns $C^*$ (or reports $t$ unreachable), and
$f_{Lim} \le C^*$ at each expansion; hence it expands no node with $f_D > C^*$.

*Proof.* If $s = t$, $U = 0$ stops it at once. Otherwise $f_{Lim} \le C^*$ at
the start, since $h_F(s) \le C^*$, $h_B(t) \le C^*$ and $\varepsilon \le C^*$.
While $f_{Lim} < C^*$, one has $U \ge C^* > f_{Lim}$ (Lemma 7.2), and the
search does not stop at lines 4 or 8. If the level $f_{Lim} = C^*$ is reached,
at its end Lemma 7.14 gives $U \le C^*$, i.e., $U = C^*$, and the search then
returns it at line 4 of the next iteration, or earlier, within the level $C^*$,
with $U \le f_{Lim} = C^*$. Hence each expansion happens at a level $f_{Lim}
\le C^*$, of a node with $f_D \le f_{Lim}$. An empty open list stops the search
with $U = C^*$ or $U = \infty$, by Lemma 7.3 and Theorem 7.6. $\square$

We remark that the order of line 7 is ours: Barley et al. leave free the order
of the expandable nodes of a level, and the order of A$^*$ (least $f$, then
largest $g$) typically reaches the meeting point within the last level much
earlier than the order of $g$; the header reports 2.5 times fewer expansions on
the 15-puzzle, a measure that we have not repeated here. Under consistent
heuristics this order satisfies Lemma 7.4, which makes the flag `consistent`
sound.

### The implementation

In each direction we keep three lazy queues: waiting, with the nodes of $f >
f_{Lim}$ by $f$; ready, with those of $f \le f_{Lim}$ by $f$ and then by the
larger $g$; and deep, with those of ready whose $g$ reached the limit, by $g$.
A node at the top of ready with $g_D \ge g_{Lim D}$ goes to deep, and it comes
back to ready when the limit passes it ([`GBFHS.h:163-174`][GBFHS.h:163-174]).
Thus a level is a sequence of moves between the queues and of expansions (cf.
Algorithm 7.3), and each move is paid once per entry and per change of the
limits. We remark that $f_{Lim}$ grows by one unit even when no open node has
its $f$ near it: on a domain with large costs, a run may then pass through many
empty levels, each of a constant cost.

## 7.8 The two directions in parallel {#sec-7-8}

The class template `ParallelBidirectional< D , K , Executor >`
([`ParallelBidirectional.h`][ParallelBidirectional.h]) runs the two directions
at once, each on a worker of an executor (the executors of
[§1.6](01-introduction.md#sec-1-6)), and the two workers share the table of the
nodes and the incumbent. A worker runs one of the searches

```cpp
--8<-- "snippets/bidirectional-sig-parallelkind.cpp"
```

where `nba` is the NBA$^*$ of Pijls and Post [[92]](Z-bibliography.md#ref-92),
which expands by $f$ and rejects some nodes as described below; with it, the
whole scheme is PNBA$^*$ of Rios and Chaimowicz
[[93]](Z-bibliography.md#ref-93). The other three kinds are the parallel
versions of the bidirectional Dijkstra's algorithm, of MM and of BAE$^*$. Its
declarations are

```cpp
template< BidirectionalDomain D , ParallelBiKind K ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
class ParallelBidirectional {
 public:
 explicit ParallelBidirectional( const D & domain ,
                                 const BidirectionalOptions & options = {} );
 void set_executor( Executor & ex ) noexcept;
 [[nodiscard]] SolutionT solve( void ) const;
 };

template< BidirectionalDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > pnba(
 const D & domain , const BidirectionalOptions & options = {} );
```

The constructor refuses `consistent` false and a weight other than 1 with
`std::invalid_argument`, since the scheme needs consistent heuristics, and
without `set_executor` each search builds an executor of two threads; the
counters `expanded_by` and `stored_by` have one entry per direction. The node
limit is checked every 64 expansions: a worker stops both workers when twice
the nodes it has created reach it.

Each worker runs the best-first search of its direction on a separate open
list, writes only the fields of its direction in the nodes ($g_D$, the parent,
the operator and the flag closed$_D$), and reads those of the other direction.
The values of $g$ and the flags closed are accessed through atomic references,
with sequentially consistent operations, and both workers lower $U$, with its
meeting node, under a lock
([`ParallelBidirectional.h:264`][ParallelBidirectional.h:264]). After each
expansion, a worker publishes the least key of its open list and, for MM, its
least $f$ and its least $g$, each as the maximum of the value published before
and the present one
([`ParallelBidirectional.h:392`][ParallelBidirectional.h:392]). It stops both
workers when its open list is empty, or when $U$ is not larger than the bound
of the algorithm computed from the two published records
([`ParallelBidirectional.h:275`][ParallelBidirectional.h:275]); for Dijkstra's
algorithm, MM and BAE$^*$ this is the bound of `BiBestFirst`, and for NBA$^*$
it is $\max\{f_{\min F}, f_{\min B}\}$. Three facts, which we now prove, make
the scheme exact.

<a id="lem-7-16"></a>**Lemma 7.16** (No meeting is missed). If both directions
set $g_F(m)$ and $g_B(m)$ of a node $m$ to finite values, then at least one of
the two workers reads the value set by the other, and offers their sum to $U$.

*Proof.* When the worker of $D$ sets $g_D(m)$ it then reads $g_{\bar D}(m)$,
and the worker of $\bar D$ does the converse
([`ParallelBidirectional.h:445-449`][ParallelBidirectional.h:445-449]). With
sequentially consistent accesses, there is a total order of the four operations
compatible with the order of each worker, and the worker whose write comes
second in that order reads after both writes, and hence reads the value of the
other. A node created by one worker with its $g$ is published by its insertion
in the concurrent table, and the other worker that finds it sees its fields;
this is the argument of the mutual exclusion of Dekker. $\square$

<a id="lem-7-17"></a>**Lemma 7.17** (The published keys). With consistent
heuristics, the least key of the nodes open or being expanded in a direction
never decreases. Hence each value published by a worker, however old, is at
most the least key of the nodes open in its direction at any later time.

*Proof.* Each worker is a sequential best-first search of its direction, whose
key satisfies Lemma 7.4 (e.g., $\max\{f, 2g\}$ for MM). A node popped has the
least key of the open list, and its children have keys not smaller than its
key, by the monotonicity along an edge; the nodes that the other worker closes,
or that $U$ prunes, only leave the list. Hence the least key over the open
nodes and the node being expanded never decreases, and a published value, which
is a least key at the time of its publication (or an earlier one, by the
maximum), stays below it. $\square$

<a id="lem-7-18"></a>**Lemma 7.18** (The rejections of NBA$^*$). Let the
heuristics be consistent, $F_{\bar D}$ be a value published by the other
direction, and $x$ be popped in $D$ with $g_D(x) + F_{\bar D} - h_{\bar D}(x)
\ge U$. Then no path from the root of $D$ to $x$ of cost $g_D(x)$, continued
through a node open in $\bar D$, costs less than $U$.

*Proof.* Let $y$ be open in $\bar D$, and $d(x, y)$ the cost of the subpath
between them. By consistency $h_{\bar D}(y) \le d(x, y) + h_{\bar D}(x)$, and a
path through $x$ and $y$ costs at least
$$g_D(x) + d(x, y) + g_{\bar D}(y) \ge g_D(x) + f_{\bar D}(y) - h_{\bar D}(x) \ge g_D(x) + F_{\bar D} - h_{\bar D}(x) \ge U,$$
since $f_{\bar D}(y) \ge F_{\bar D}$ by Lemma 7.17. $\square$

NBA$^*$ closes unexpanded the node $x$ popped in $D$ when $f_D(x) \ge U$, or
when $g_D(x) + F_{\bar D} - h_{\bar D}(x) \ge U$ (section 3 of
[[92]](Z-bibliography.md#ref-92)). The second test is the *rejection* of Lemma
7.18 ([`ParallelBidirectional.h:472-476`][ParallelBidirectional.h:472-476]).

<a id="thm-7-19"></a>**Theorem 7.19** (Correctness of the parallel
searches). With consistent heuristics, `ParallelBidirectional` returns $C^*$
for each kind, or reports $t$ unreachable, unless the node limit stops it.

*Proof.* By Lemma 7.16, each node that gets both values of $g$ lowers $U$ to at
most their sum, which is the first remark of the proof of Lemma 7.3. Each
direction pops by a key of Lemma 7.4, and hence one can see that (H) holds,
with the same proof for the two workers interleaved. Also, a node of an optimal
path is neither pruned nor rejected while $U > C^*$: one has $f_D \le C^* < U$
and, by Lemma 7.18 with the node of the other direction given by Lemma 7.3,
$g_D(x) + F_{\bar D} - h_{\bar D}(x) \le C^* < U$. Hence Lemma 7.3 holds at any
time at which a worker reads its open list: if $U > C^*$, the open list of each
direction has a node of $P$ with its optimal $g$, and a worker that finds its
list empty thus has $U \le C^*$. Consider now a worker that finds $U \le
\mathrm{LB}$, with LB computed from the published values: by Lemma 7.17, its LB
is at most the bound of Theorem 7.5 computed from the present open lists (each
bound is a nondecreasing function of the least keys), and hence at most $C^*$
if $U > C^*$; Theorem 7.6 concludes. $\square$

The test `bidir-thread` (and `bidir-openmp`, `bidir-fastflow` on the other
executors) runs the four kinds on 50 states of the 8-puzzle, on two instances
of the 15-puzzle of Korf, on random grids with 8 moves and on random directed
graphs, and half of the runs on the 8-puzzle build a new executor, and the
others are given one. Of course, a parallel bidirectional search has two
workers, and therefore its speedup is at most 2; also, it pays the overheads of
[§8.2](08-shared-memory.md#sec-8-2), i.e., the search overhead of a bound read
late (a worker may expand nodes that a sequential search would have pruned by a
bound already reached) and the coordination overhead of the atomic accesses to
the shared nodes.

## 7.9 Parallel external-memory bidirectional search {#sec-7-9}

The class template `PEMSearch< D , K , Executor >`
([`PEMSearch.h`][PEMSearch.h]) is the framework PEM-BiHS of Siag, Shperberg,
Felner and Sturtevant [[57]](Z-bibliography.md#ref-57). Its open and closed
lists are files on the disk, and the nodes of one bucket are loaded, checked
for duplicates and expanded by the workers of an executor. Each search of the
class is one of

```cpp
--8<-- "snippets/bidirectional-sig-pemkind.cpp"
```

i.e., PEM-A$^*$ (forward only), PEMM, the external MM of Sturtevant and Chen
[[94]](Z-bibliography.md#ref-94), and PEM-BAE$^*$; the reverse PEM-rA$^*$ is
PEM-A$^*$ on a `ReverseDomain`. The domain must be a `DiskDomain` of
[§4.6](04-memory-bounded.md#sec-4-6) (integer costs and trivially copyable
states), and also a `BidirectionalDomain` for PEMM and PEM-BAE$^*$. Its
declarations are

```cpp
template< DiskDomain D , PEMKind K ,
          ADS::WorkerExecutor Executor = ADS::ThreadPool >
 requires( ( K == PEMKind::astar ) || BidirectionalDomain< D > )
class PEMSearch {
 public:
 PEMSearch( const D & domain , unsigned threads = 1 ,
            const SearchOptions & options = {} , PEMOptions pem = {} );
 void set_executor( Executor & ex ) noexcept;
 [[nodiscard]] SolutionT solve( void ) const;
 };
```

with the options of the disk

```cpp
--8<-- "snippets/bidirectional-sig-pemoptions.cpp"
```

i.e., the number $P$ of partitions of a bucket (16 by default, 0 being
refused), the directory of the files, the bytes of the buffer of each file that
is written, and a cap on the bytes held on the disk, at which the search stops
with `Status::disk_limit`. In the directory given (or in the temporary
directory of the system) the search makes a directory for its files, which it
removes at the end, also when it stops at a limit or throws. There are
`threads` workers, and of the `SearchOptions`, `node_limit` bounds the states
written to the closed files, and a weight other than 1 is refused. Its counters
are `expanded` (the states expanded, and `expanded_by` by direction),
`generated`, `duplicates` (the copies dropped), `stored` (the states written to
the closed files), `sent` (the bytes written to the disk) and `messages` (here,
the most bytes held on the disk at once).

### The buckets and the meetings

A node of the direction $D$ is kept in the bucket $(D, g_D, h_F, h_B)$ (with
$h_B = 0$ for PEM-A$^*$), which carries all that the priorities and the bounds
need. Each bucket is split into $P$ partitions by the hash of the state, and a
partition is a file per worker that writes to it. Two copies of a state, in
either direction, have the same $h_F$ and $h_B$ and the same partition; hence
the duplicates of a partition are in the same partition of the buckets with the
same heuristics, and so are the copies of the state in the other direction.

**Algorithm 7.4** (PEM-BiHS, the pseudocode of the header
[`PEMSearch.h:23-39`][PEMSearch.h:23-39]).

```text
 1  procedure PEM-BiHS(s, t)
 2     U ← ∞; OPEN_F ← {s}; OPEN_B ← {t}
 3     while OPEN_F ≠ ∅ ∧ OPEN_B ≠ ∅ ∧ U > LB do
 4        D ← the direction chosen; b ← the best bucket of OPEN_D
 5        for each partition p of b, in parallel, do
 6           T ← the states of the files of p, without copies
 7           T ← T ∖ the partition p of the closed buckets of D with the
                same h_F and h_B and a smaller g
 8           write T to the partition p of CLOSED_D(b)
 9           U ← min(U, g + g') for each state of T found, with g', in
                the partition p of the buckets of the other direction
                with the same h_F and h_B
10           for each u in T and each child v of u in D do
                append v to the partition of v of its bucket of OPEN_D
11     return U
```

We remark that lines 6 to 8 are the processing of a file of DDD-A$^*$ (cf.
Algorithm 4.6), here for a partition made of the files of all the workers. Line
9 is the *delayed solution detection* of Sturtevant and Chen
[[94]](Z-bibliography.md#ref-94): a meeting is found when one of its two copies
is expanded, by a scan of the open and the closed buckets of the other
direction with the same heuristics, in the same partition (Figure 7.3), while
PEM-A$^*$ tests the goal when a node is generated, as Siag et al. do. The best
bucket of a direction is the one of least priority ($f$ for PEM-A$^*$,
$\max\{f, 2g\}$ for PEMM, $b$ for PEM-BAE$^*$), of least $g$ among those, and
then of least $h$ of the direction; this is the lower-$g$-first order of Siag
et al. (section 3.3 of [[57]](Z-bibliography.md#ref-57)). PEMM takes the
direction of least priority, the forward one on a tie, and PEM-BAE$^*$
alternates them. Here the workers take the partitions of a bucket from a shared
counter, and each worker writes the children through a buffer per file of the
buckets it fills, and the new buckets are added to the open ones when the
bucket is done ([`PEMSearch.h:519-684`][PEMSearch.h:519-684]). At the end, the
path is traced back through the hashes of the parents, as by Algorithm 4.7,
from each of the two records of the meeting.

![Figure 7.3. The buckets $(D, g, h_F, h_B)$ of PEM-BiHS on the disk, each
split into $P = 4$ partitions by the hash of the state. The partition 2 of
the forward bucket $(F, 5, 4, 3)$ is loaded and cleared of its copies and of
the states of the partition 2 of the closed forward buckets of the same
heuristics and a smaller $g$ (here $(F, 3, 4, 3)$); its states are looked up
in the partition 2 of the backward buckets of the same heuristics, open and
closed, to find the meetings; and each child is appended to the partition of
its state in the bucket of its $g$ and heuristics, by the worker that
generated it.](../figures/bidirectional-pem.svg){#fig-7-3}

Note that, since a meeting is found only when one of its copies is expanded,
two open copies of a state on an optimal path may be unknown to each other.
Thus the bound $g_{\min F} + g_{\min B} + \varepsilon$ is not valid (section
3.3 of [[94]](Z-bibliography.md#ref-94)), and the bound of line 3 is $f_{\min
F}$ for PEM-A$^*$,
$$\max\{C, f_{\min F}, f_{\min B}, g_{\min F} + g_{\min B}\}$$
for PEMM and, for PEM-BAE$^*$, the largest of $\lceil (b_{\min F} + b_{\min B})
/ 2 \rceil$, $g_{\min F} + g_{\min B}$, $f_{\min F} + d_{\min B}$ and $d_{\min
F} + f_{\min B}$ (the bounds of Alcázar et al.
[[89]](Z-bibliography.md#ref-89), which the PEM-BAE$^*$ of Siag et al.
computes), all of them computed from the keys of the open buckets
([`PEMSearch.h:380-406`][PEMSearch.h:380-406]).

<a id="lem-7-20"></a>**Lemma 7.20** (The order of the buckets). Let the costs
be positive, and let the key of a bucket be nondecreasing along each edge
generated, which the search checks. Then, in each direction, no bucket gets a
record after its expansion, and the copy of a state with its optimal $g$ in
that direction, if generated, is expanded before all the other copies of it,
and the later ones are dropped at line 7.

*Proof.* Let the bucket $b$ of $D$ be expanded: the children of $b$ and of the
later buckets of $D$ have keys not smaller than their parents and, for equal
keys, a larger $g$ (since the costs are positive). Since the best bucket is the
one of least key and then of least $g$, the pairs $(\mathrm{key}, g)$ of the
buckets expanded in $D$ never decrease. A child of a later bucket $b'$ has
either a key larger than that of $b$, or the same key and $g > g(b') \ge g(b)$;
in both cases it is not in $b$. For the second claim, the copies of a state $v$
have the same heuristics, and in its buckets the key is a strictly increasing
function of $g$ (as in Lemma 7.4). Indeed, the copy with the optimal $g$ is in
the bucket of least key among them. By the induction along an optimal path of
[Lemma 4.14](04-memory-bounded.md#lem-4-14), whose argument applies to any key
monotone along the edges, it is there when that bucket is expanded. A later
copy, of larger $g$, is in a later bucket of the same heuristics and the same
partition, and line 7 finds $v$ in the closed partition. $\square$

The search throws `std::logic_error` on an edge of a cost that is not positive,
and on a child whose bucket has a smaller key than its parent
([`PEMSearch.h:638`][PEMSearch.h:638]). For PEM-A$^*$, the second case is,
e.g., a heuristic that is not consistent, and the test `pem-thread` checks the
refusal on the graph of [§2.7](02-best-first.md#sec-2-7). For PEMM and
PEM-BAE$^*$ the check is on their key, $\max\{f, 2g\}$ or $b$, and a heuristic
that is not consistent but keeps the key monotone on the edges generated is not
refused; for PEM-BAE$^*$ its bounds then need not hold, and the consistency of
both heuristics is an assumption of the search, as stated in the header.

<a id="lem-7-21"></a>**Lemma 7.21** (Delayed meetings). Let the hypotheses of
Lemma 7.20 hold, $P = (p_0, \dots, p_k)$ be an optimal path, and $U > C^*$
between two buckets; then there are $i \le j$ such that $p_i$ has a copy in an
open forward bucket with $g_F = d_F(p_i)$, and $p_j$ one in an open backward
bucket with $g_B = d_B(p_j)$.

*Proof.* We say that $p_l$ is settled in $D$ if its copy with the optimal $g$
has been expanded in $D$. If $p_l$ is settled in both directions, the later of
the two expansions scanned the closed partition that holds the earlier one, and
set $U \le C^*$. As in Lemma 7.3, the first node $p_i$ not settled forward has
its optimal copy in an open forward bucket (generated by $p_{i-1}$ with $g +
h_F \le C^* < U$, which line 7 does not drop, by Lemma 7.20). Similarly, the
last node $p_j$ not settled backward has its optimal copy in an open backward
bucket. Suppose $j < i$: then $p_j$ is settled forward and, if $j + 1 < i$,
$p_{j+1}$ is settled in both directions, which is excluded, and thus $j = i -
1$, $p_{i-1}$ is settled forward and open backward, and $p_i$ is settled
backward and open forward. The forward expansion of $p_{i-1}$ generated the
forward copy of $p_i$, and the backward expansion of $p_i$ generated the
backward copy of $p_{i-1}$. Whichever of the two expansions came later found
the copy generated by the earlier one (the forward copy of $p_i$, or the
backward copy of $p_{i-1}$, both still open) in an open bucket of the other
direction with the same heuristics and partition, and it set $U \le C^*$. Hence
$i \le j$. $\square$

<a id="thm-7-22"></a>**Theorem 7.22** (Correctness of PEM-BiHS). Under the
hypotheses of Lemma 7.20, and with consistent heuristics for PEM-BAE$^*$, the
searches of `PEMSearch` return $C^*$, or report $t$ unreachable, unless a limit
stops them first. Each state is expanded at most once in each direction.

*Proof.* With $n = p_i$ and $m = p_j$ of Lemma 7.21, one has $g_F(n) + g_B(m)
\le C^*$ (now possibly with $n = m$), which proves the bound $g_{\min F} +
g_{\min B}$. The bounds of $f$ and $C$ follow as in Theorem 7.5 (1 and 3),
whose proofs use only $g_F(n) + g_B(m) \le C^*$, and those of BAE$^*$ as in
Theorem 7.5 (4 and 5), whose proofs use $d(n, m) \ge 0$. The key of an open
bucket is the key of the copies in it, and hence the least keys of the buckets
are at most those of the copies of $n$ and $m$. Theorem 7.6, with Lemma 7.21 in
place of Lemma 7.3, gives the result for PEMM and PEM-BAE$^*$, while PEM-A$^*$
is DDD-A$^*$ with the goal tested at the generation and the stop at $U \le
f_{\min F}$, which is correct by the same argument with $j = k$. By Lemma 7.20
a state is written to a closed file of $D$ only once, and it is expanded only
from there. $\square$

### The implementation

We index the closed buckets by $(D, h_F, h_B)$, each with the list of the
values of $g$ that have a closed file, and thus the scan of line 7 reads only
the closed partitions of the same heuristics and of a smaller $g$, and that of
line 9 those of the other direction with the same heuristics. Also, the open
buckets are a set, scanned linearly to choose the best and to compute the
bound, at a cost proportional to their number, which is usually rather small
relative to the records. For each worker, the memory holds the hash table of
the partition it is processing and a buffer of `buffer_bytes` for each file it
writes, i.e., at most $P$ for each distinct bucket of the children of the
current one. As for DDD-A$^*$ ([Proposition
4.16](04-memory-bounded.md#prop-4-16)), the disk holds all the closed files to
the end, since a duplicate may be in any closed bucket of the same heuristics
and the path is traced through them. The test `pem-thread` (and `pem-openmp`,
`pem-fastflow`) runs PEM-A$^*$, PEMM, PEM-BAE$^*$ and PEM-rA$^*$ against the
cost of A$^*$: with 1, 2 and 4 workers and 5 partitions on 3 states of the
8-puzzle, on random grids with 4 moves and on the map of Romania, and with 2
workers on random walks of 40 moves of the 15-puzzle. It checks the stop at a
cap of 4096 bytes on the disk and at a node limit of 1000, the refusal of an
inconsistent heuristic, and that no file is left when the searches end, stop or
throw.

## 7.10 A$^*$ with immediate duplicate detection on the disk {#sec-7-10}

The external searches above, and DDD-A$^*$ of
[§4.6](04-memory-bounded.md#sec-4-6), detect the duplicates in batches by
scanning files, which suits a disk read sequentially. On a solid state drive a
random read usually costs much less, and Lin and Fukunaga
[[95]](Z-bibliography.md#ref-95) revisited the *immediate* duplicate detection
on the disk: their A$^*$-IDD keeps the closed list in a file read at random,
and in memory only a "small" table of the positions of its records. The class
template `IDDAStar< D >` ([`IDDAStar.h`][IDDAStar.h]) implements it for a
`DiskDomain`, with the declarations

```cpp
template< DiskDomain D >
class IDDAStar {
 public:
 explicit IDDAStar( const D & domain , const SearchOptions & options = {} ,
                    IDDOptions idd = {} );
 [[nodiscard]] SolutionT solve( void ) const;
 };

template< DiskDomain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > idd_astar(
 const D & domain , const SearchOptions & options = {} ,
 const IDDOptions & idd = {} );
```

and the options

```cpp
--8<-- "snippets/bidirectional-sig-iddoptions.cpp"
```

i.e., the directory of the files (as for `PEMSearch`), the number $S$ of places
of the table of the positions, which bounds the closed states to $7S/8$, the
number $P$ of partitions, the number $B$ of states of a block of the closed
file, the number of records of the top of an open bucket kept in memory, and
the cap on the bytes on the disk. A zero option, or $S < 8B$, is refused with
`std::invalid_argument`, as is a weight other than 1.

### The lists

We keep OPEN as a set of buckets $(f, g)$, ordered by least $f$ and then by
largest $g$, each a stack of records on the disk with its top in memory. A push
adds to the top, and when the top holds twice the records of `open_block`, its
lower half goes to the file of the bucket; a pop takes from the top, and when
the top is empty the last block of the file is read back
([`IDDAStar.h:350-398`][IDDAStar.h:350-398]). A record holds the state, its
$g$, the state of its parent and, in an `OperatorDomain`, its operator. CLOSED
has three parts. First, the records are in a file, written only at its end (or
rewritten in place). Second, there is a table in memory of $S$ places of
$\lceil \log_2(S + 1) \rceil$ bits each, packed in words; a place holds the
position plus one of a record in the file (0 for an empty place), and the
places are probed linearly from the hash of the state. Third, there are $P$
buffers in memory, one per partition of the states by a second hash. A buffer
collects the states closed until it holds $B$ of them, and it is then written
to the file as one block, whose positions are put in the table. Since each
block holds the states of one partition, a vector of one partition number per
block shows whether the record that a place points to may be the state looked
for; the file is read only when it may (section 4 of
[[95]](Z-bibliography.md#ref-95)).

**Algorithm 7.5** (A$^*$-IDD, the pseudocode of the header
[`IDDAStar.h:25-39`][IDDAStar.h:25-39]).

```text
 1  procedure A*-IDD(s₀)
 2     push (s₀, g = 0) in OPEN
 3     while OPEN ≠ ∅ do
 4        n ← pop the top of the bucket of least f and largest g
 5        if n is in the buffer of its partition or, probing the table
             and reading only the records of blocks of its partition,
             in the file, with a g not larger, then continue
 6        put n in the buffer of its partition (or update its record)
 7        if the buffer has B states, append them to the file as a
             block, and put their positions in the table
 8        if n is a goal then return the path, traced by the parents
 9        push each successor of n, with its g and parent, in OPEN
10     return unsolvable
```

The duplicates are thus detected when a node is popped, before its expansion,
instead of when it is generated. OPEN holds all the copies generated, as the
lazy list of [§3.1](03-data-structures.md#sec-3-1), and a state popped with a
smaller $g$ than its closed record is reopened: its record is rewritten in the
buffer or in place in the file. At the end, the path is traced from the goal by
looking up the parent of each record among the closed states, down to the root,
whose recorded parent is the root.

<a id="lem-7-23"></a>**Lemma 7.23** (The lookup is exact). At each pop, (a)
each closed state is either in the buffer of its partition or in the file, but
not in both, and in the file it is at a position $a$ such that the place
holding $a + 1$ is reached by the linear probing from its hash before an empty
place; (b) the block $\lfloor a / B \rfloor$ of the file holds only states of
one partition, the one recorded for it; (c) the table has an empty place. Hence
the probe of line 5 terminates, and it finds the closed record of $n$ if there
is one, reading from the disk only records of the partition of $n$.

*Proof.* (b) A buffer is written only when it holds exactly $B$ states, all of
its partition; hence each block is full and of one partition, and its partition
is appended to the vector at the same time
([`IDDAStar.h:327-348`][IDDAStar.h:327-348]). (a) A state is put in the buffer
only when it is found neither there nor in the file. It leaves the buffer only
when the buffer is written, and then each record of the block gets a place
found by the linear probing from its hash. The places are never cleared; hence
the probing from the hash of a state reaches its place before any empty one
([Lemma 3.9](03-data-structures.md#lem-3-9)), and a reopening rewrites the
record at the same position. (c) The search stops when $8\, \mathit{used} \ge
7S$ after a step; before a flush, then, fewer than $7S/8$ places are taken, and
a flush adds $B \le S/8$ of them (the constructor requires $S \ge 8B$), and
hence fewer than $S$ places are ever taken, the stop coming right after the
flush that reaches $7S/8$. Thus the probe meets an empty place, and a place
whose block has another partition is skipped without reading the record, which
by (b) cannot be that of $n$. $\square$

<a id="thm-7-24"></a>**Theorem 7.24** (Correctness of A$^*$-IDD). With an
admissible heuristic and positive costs, A$^*$-IDD returns a path of cost
$C^*$, or reports the instance unsolvable on a finite graph, unless a limit
stops it. With a consistent heuristic each state is expanded once, and
$\mathit{expanded} = \mathit{stored}$.

*Proof.* By Lemma 7.23 the test of line 5 is exact, and hence the search is
A$^*$ with the duplicates of OPEN dropped at the selection instead of at the
generation, the reopening of a closed state popped with a smaller $g$, and the
goal tested at the selection. A state that is popped and kept has a $g$ smaller
than that of all the earlier closures of it. Thus the invariant of OPEN of
[Lemma 2.7](02-best-first.md#lem-2-7) holds (an optimal path has a copy of one
of its nodes in OPEN with its optimal $g$, all its predecessors on the path
being closed with the optimal one), and the argument of [Theorem
2.9](02-best-first.md#thm-2-9) gives that the first goal popped has $g = C^*$.
The path traced through the records of the parents has a cost at most $g$,
since each parent was closed with a $g$ that has only decreased since, and it
ends at the root, since $g$ strictly decreases along it. With a consistent
heuristic the order by $f$ closes each state with its optimal $g$ ([Theorem
2.13](02-best-first.md#thm-2-13)); then no reopening occurs, and each closed
state is stored and expanded once. On a finite graph the closures are finitely
many, as in [Theorem 2.11](02-best-first.md#thm-2-11). $\square$

<a id="prop-7-25"></a>**Proposition 7.25** (Memory and accesses of A$^*$-IDD).
The memory holds the table, of $S \lceil \log_2(S+1) \rceil / 8$ bytes, the $P$
buffers, with fewer than $P B$ records, the vector of the partitions of the
blocks, and the tops of the open buckets, with fewer than $2\,
\mathit{open\_block}$ records each, while on the disk there are the closed
records and the bottoms of the open buckets. A lookup in the file reads at most
one record for each place of its partition that the probe visits.

*Proof.* The bounds follow from the construction above and from Lemma 7.23.
$\square$

Hence, for $S = 2^{20}$ (the default), a place takes 21 bits (i.e., $\lceil
\log_2(2^{20} + 1) \rceil$) and the table about 2.75 MB, for up to 917 504
closed states. Here `messages` counts the records read at random from the
closed file, i.e., the probes that reach the disk, and `sent` the bytes
written. The search stops with `Status::node_limit` when the table is $7/8$
full or when the closed states reach `node_limit`, and with
`Status::disk_limit` when the files exceed the cap. The test `idd` runs
A$^*$-IDD with $S = 2^{19}$, $B = 16$, 7 partitions and open blocks of 8
records, so that most of the closed list is on the disk. It runs it on 50
states of the 8-puzzle against their exact distances, on the instance 55 of the
15-puzzle of Korf (where it checks that $\mathit{expanded} = \mathit{stored}$
and that the closed file is read), on random grids and on the map of Romania.
It checks the reopening on the graph with an inconsistent heuristic of
[§2.7](02-best-first.md#sec-2-7), the stop when the table of 1024 places is
full and at a cap of 4096 bytes on the disk, the removal of the files, and the
refusal of a table smaller than 8 blocks.

We remark that a duplicate found in the file is read twice, once by the probe
and once to compare its $g$ ([`IDDAStar.h:429-435`][IDDAStar.h:429-435]), and
`messages` counts two reads for it. Also, the table is never cleared: a state
closed once keeps its place for the rest of the search, which is why the closed
states are bounded by $7S/8$ rather than by the disk.

## 7.11 Benchmark {#sec-7-11}

We compare the bidirectional searches with A$^*$ and with the reverse A$^*$
(A$^*$ on the `ReverseDomain`) on two domains whose heuristics are consistent
and whose moves are reversible: the instances of the 15-puzzle of Korf, with
the Manhattan distance in both directions, and random grids with 8 moves and
the octile heuristic. Also, we compare the external searches with DDD-A$^*$ of
[§4.6](04-memory-bounded.md#sec-4-6) on the same instances of the 15-puzzle.
The runs are made by `AStarBench`, whose option `--algo` names the searches of
this chapter (`ra`, `bidij`, `mm`, `mme`, `fmm`, `bae`, `nbs`, `dvcbs`,
`gbfhs`, `bibfs`, `pnba`, `pdij`, `pmm`, `pbae`, `pema`, `pemra`, `pemm`,
`pembae` and `idd`). Each run is made twice, in separate processes on an idle
machine, and the gap between the two runs is the A/A noise given in the note of
the tables.

**Table 7.1.** The in-memory bidirectional searches. For each set of instances
and each search, the columns give the expansions in the forward and in the
backward direction, the nodes stored, and the seconds summed over the instances
(each the median of its two runs). The machine and the compiler are given in
the note of the table.

--8<-- "tables/bidirectional-memory.md"

**Table 7.2.** The external searches on the 15-puzzle. For each search and
number of workers, the columns give the states expanded, the states stored in
the closed files, the bytes written to the disk, the most bytes held on the
disk at once, and the seconds summed over the instances. The directory of the
files is on the local drive of the machine named in the note.

--8<-- "tables/bidirectional-external.md"

We measure the time by the timer of the `Solution`. It includes the creation of
the table, the start of the threads of the executor when none is lent, the
search, the writing and the reading of the files and the rebuilding of the
path; it does not include the construction of the domain nor the removal of the
directory of the files. We leave the tables as placeholders until the campaign
is run on an idle machine.

## 7.12 Notes and references {#sec-7-12}

Bidirectional search is due to Pohl [[91]](Z-bibliography.md#ref-91), who also
proposed to expand in the direction with fewer open nodes. The "nipping" of the
nodes closed by the other direction is from the BS$^*$ of Kwa
[[90]](Z-bibliography.md#ref-90). The front-to-end searches of this chapter
come from a line of work started by MM of Holte, Felner, Sharon and Sturtevant
[[61]](Z-bibliography.md#ref-61). Its journal version
[[81]](Z-bibliography.md#ref-81) gives the bounds and the proof of the meeting
in the middle that we followed in sections 7.2 to 7.4. MMe, with its improved
priority, is from Sharon et al. [[82]](Z-bibliography.md#ref-82), and fMM and
the structure of the minimum covers of the must-expand graph are from Shaham et
al. [[83]](Z-bibliography.md#ref-83). The sufficient conditions for node
expansion and the must-expand graph are due to Eckerle et al.
[[84]](Z-bibliography.md#ref-84); NBS is due to Chen et al.
[[85]](Z-bibliography.md#ref-85), whose section 7 gives the queues waiting and
ready that we use, DVCBS, NBS$_A$ and the eager variants to Shperberg et al.
[[86]](Z-bibliography.md#ref-86), and GBFHS to Barley et al.
[[87]](Z-bibliography.md#ref-87). BAE$^*$ was introduced by Sadhukhan
[[88]](Z-bibliography.md#ref-88), whose paper we could not consult. Alcázar,
Riddle and Barley [[89]](Z-bibliography.md#ref-89) derived its bound, and the
individual bounds used with `all_bounds`, from a unified view of the errors of
the heuristics, and we followed their definition. The immediate and the delayed
solution detection, and the external MM, are due to Sturtevant and Chen
[[94]](Z-bibliography.md#ref-94), and the framework PEM-BiHS, with PEM-BAE$^*$
and the reverse A$^*$, to Siag et al. [[57]](Z-bibliography.md#ref-57). NBA$^*$
is due to Pijls and Post [[92]](Z-bibliography.md#ref-92), and its parallel
version PNBA$^*$ to Rios and Chaimowicz [[93]](Z-bibliography.md#ref-93), whose
paper we know only from its abstract. Lemma 7.16, by which no meeting is
missed, uses the argument of the mutual exclusion of Dekker; we give it for the
atomic accesses of the library, instead of for the scheme of their paper.
A$^*$-IDD is due to Lin and Fukunaga [[95]](Z-bibliography.md#ref-95), and the
hash-based delayed duplicate detection on which the external searches rest is
due to Korf (see [§4.6](04-memory-bounded.md#sec-4-6)).

We implement only front-to-end searches, i.e., searches whose heuristics
estimate the distance to the root of the other direction; the front-to-front
searches, which estimate the distance to the nodes of the other front, are not
treated here. Also, the searches of sections 7.4 to 7.7 are sequential. Over
the two directions, the only parallelism is that of [§7.8](#sec-7-8), with two
workers, while the parallelism within a direction is that of the external
searches of [§7.9](#sec-7-9), over the partitions of a bucket. A comparison of
the expansions of MM, NBS, DVCBS, GBFHS and BAE$^*$ across domains is in
[[86]](Z-bibliography.md#ref-86) (table 1) and in
[[89]](Z-bibliography.md#ref-89); Table 7.1 is to give our own counts, on our
instances.

## 7.13 Exercises {#sec-7-13}

**Exercise 7.1** (`ctest -R PBeFS.bidir`). Show that the predecessors of a cell
of a grid with weighted cells, as defined by the domain, are the reversed moves
at the cost of the cell left. Deduce that the forward and the backward graph
are the same graph with the edges reversed. The test checks it on 4 random
grids of $11 \times 7$ cells with 8 moves.

**Exercise 7.2** (`ctest -R PBeFS.bidir`). Prove the bound of MMe, $\min_D
\max\{f_D, 2 g_D + \varepsilon\} \le C^*$ while $U > C^*$, directly from
Lemma 7.3, and show that it is never smaller than the bound $C$ of MM. Then
build a graph of three nodes on which MM expands a node that MMe does not.
The test runs both on 20 states of the 8-puzzle with 6 variants of the
options.

**Exercise 7.3** (`ctest -R PBeFS.bidir`). Theorem 7.7 assumes the choice
of the direction by priority. Build a small instance on which MM with
`SideChoice::alternate` expands a node with $g_D > C^*/2$, and explain why
the cost returned is still $C^*$. The test checks the meeting in the middle
only with the default choice.

**Exercise 7.4** (`ctest -R PBeFS.bidir`). On the graph of five nodes of
the test, with the admissible but inconsistent heuristic $h = (0, 3, 0, 1,
0)$, run MM by hand with `consistent` false, and show where a node is
reopened. Then find the step of the proof of Lemma 7.4 that fails with this
heuristic, and explain why `bae` refuses `consistent` false.

**Exercise 7.5** (`ctest -R PBeFS.bidir`). Prove that, with `eager`, every
pair with $\mathrm{lb} \le C^*$ has an end expanded by NBS before the search
stops. Which line of Algorithm 7.2 is then never executed, and why?

**Exercise 7.6** (`ctest -R PBeFS.bidir`). In Figure 7.2, change the size
of the cluster $a_0$ to 1 and recompute the weights of the covers $K_0,
\dots, K_4$ of Lemma 7.12; which cluster does DVCBS expand? Show that NBS,
on the same clusters, may expand nodes of the four clusters $a_0$, $a_1$,
$b_0$ and $b_1$ before it expands any of $b_2$.

**Exercise 7.7** (`ctest -R PBeFS.bidir`). Show that GBFHS with $\varepsilon =
0$ is still correct, and that each level then raises the sum of the two limits
by one. The test runs GBFHS and the bidirectional breadth-first search on
random directed graphs with arcs of cost 0 to 9, and with the option
$\varepsilon = 0$ on the 8-puzzle and the grids. Which bound of the other
searches becomes weaker with $\varepsilon = 0$?

**Exercise 7.8** (`ctest -R PBeFS.bidir-thread`). In Lemma 7.16, replace the
sequentially consistent accesses by a store with release order and a load with
acquire order, and give an execution of the two workers on a path of three
nodes in which neither reads the $g$ written by the other. Then show that the
published keys of Lemma 7.17 may be read late without harm, whereas a published
key that could decrease would break Theorem 7.19.

**Exercise 7.9** (`ctest -R PBeFS.pem-thread`). Lemma 7.21 allows $i = j$,
which makes the bound $g_{\min F} + g_{\min B} + \varepsilon$ invalid for the
delayed detection. Build an instance with two open copies of a state on an
optimal path, one per direction, on which that bound exceeds $C^*$ while $U >
C^*$, and check that PEMM uses the bound without $\varepsilon$. The test runs
PEMM, PEM-BAE$^*$, PEM-A$^*$ and PEM-rA$^*$ with 1, 2 and 4 workers.

**Exercise 7.10** (`ctest -R PBeFS.pem-thread`). Show that the bytes on the
disk of `PEMSearch` never fall below the bytes of the closed files, and
estimate them, for the 15-puzzle, from the number of states stored and the
size of a record (the state, the hash of the parent and the operator). The
test checks the stop at a cap of 4096 bytes.

**Exercise 7.11** (`ctest -R PBeFS.idd`). With $S = 2^{19}$ places, how
many bits does a place of the table of A$^*$-IDD take, and how many closed
states can the search hold? Show that, with 1024 places and blocks of 64
states as in the test, the search stops before the probing of Lemma 7.23
could find the table full, and count the places taken at the stop.

**Exercise 7.12** (`ctest -R PBeFS.idd`). On the graph of five nodes with the
inconsistent heuristic, follow Algorithm 7.5, show which state is reopened, and
whether its record is rewritten in the buffer or in the file with the options
of the test ($B = 16$, 7 partitions). Then explain why the test checks
$\mathit{expanded} = \mathit{stored}$ on the 15-puzzle, whose heuristic is
consistent, but not on this graph.

<!-- the links to the sources, written by srclinks.py -->

[BiBestFirst.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/BiBestFirst.h
[BiBestFirst.h:29-40]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/BiBestFirst.h#L29-L40
[BiBestFirst.h:342]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/BiBestFirst.h#L342
[bidirectional.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/bidirectional.cpp
[Bidirectional.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h
[Bidirectional.h:191]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L191
[Bidirectional.h:239]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L239
[Bidirectional.h:331]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L331
[Bidirectional.h:391]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L391
[Bidirectional.h:476]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L476
[Bidirectional.h:497]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L497
[Bidirectional.h:640]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L640
[Bidirectional.h:709]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Bidirectional.h#L709
[Domain.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h
[Domain.h:476]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h#L476
[Domain.h:496]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h#L496
[Domain.h:512]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h#L512
[Domain.h:773]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h#L773
[DVCBS.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/DVCBS.h
[DVCBS.h:123-134]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/DVCBS.h#L123-L134
[DVCBS.h:199-232]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/DVCBS.h#L199-L232
[DVCBS.h:235-272]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/DVCBS.h#L235-L272
[GBFHS.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/GBFHS.h
[GBFHS.h:16-27]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/GBFHS.h#L16-L27
[GBFHS.h:163-174]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/GBFHS.h#L163-L174
[GBFHS.h:177-206]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/GBFHS.h#L177-L206
[IDDAStar.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/IDDAStar.h
[IDDAStar.h:25-39]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/IDDAStar.h#L25-L39
[IDDAStar.h:327-348]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/IDDAStar.h#L327-L348
[IDDAStar.h:350-398]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/IDDAStar.h#L350-L398
[IDDAStar.h:429-435]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/IDDAStar.h#L429-L435
[NBS.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/NBS.h
[NBS.h:234-246]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/NBS.h#L234-L246
[NBS.h:35-45]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/NBS.h#L35-L45
[ParallelBidirectional.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h
[ParallelBidirectional.h:264]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h#L264
[ParallelBidirectional.h:275]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h#L275
[ParallelBidirectional.h:392]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h#L392
[ParallelBidirectional.h:445-449]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h#L445-L449
[ParallelBidirectional.h:472-476]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelBidirectional.h#L472-L476
[pbidirectional.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/pbidirectional.cpp
[pem.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/pem.cpp
[PEMSearch.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/PEMSearch.h
[PEMSearch.h:23-39]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/PEMSearch.h#L23-L39
[PEMSearch.h:380-406]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/PEMSearch.h#L380-L406
[PEMSearch.h:519-684]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/PEMSearch.h#L519-L684
[PEMSearch.h:638]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/PEMSearch.h#L638
