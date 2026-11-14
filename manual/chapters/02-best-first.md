# 2. Best-first search and A* {#ch-2}

[Source: `include/PBeFS/astar/AStar.h` (`AStar`, `astar`),
`include/PBeFS/astar/Node.h` (`Node`, `trace_path`),
`include/PBeFS/astar/Solution.h` (`Solution`, `SearchStats`,
`SearchOptions`, `Status`), `include/PBeFS/astar/Domain.h`
(`OperatorDomain`, `for_each_child`),
`include/PBeFS/astar/ExplicitGraph.h`; the tests `astar`, `grid`,
`graph`, `limits` and `reopen` of `test/test.cpp`]

We are concerned in this chapter with the sequential best-first search
of the graph of a `Domain` (cf. [Chapter 1](01-introduction.md#ch-1),
whose notation we use throughout), and in particular with A$^*$
[[9]](Z-bibliography.md#ref-9), on which the other searches of the
library are built. We first describe the generic best-first search, the
difference between a search of the tree of the paths and a search of the
graph, and the priorities that give breadth-first search, depth-first
search, the algorithm of Dijkstra and A$^*$. Then we state A$^*$ as the
library implements it, and we prove its correctness under an admissible
heuristic and its stronger properties under a consistent one. The last
sections discuss the reopening of the closed nodes, the breaking of the
ties, the pruning of the move back to the parent and the counters by
which a search is measured. The data structures of OPEN and CLOSED, with
the cost of their operations, are the subject of [Chapter
3](03-data-structures.md#ch-3), and we only point to them here. Also,
the searches that hold less than CLOSED in memory (IDA$^*$ and the
search with delayed duplicate detection on the disk) are those of
[Chapter 4](04-memory-bounded.md#ch-4), and the heuristics themselves
are those of [Chapter 5](05-heuristics.md#ch-5).

Throughout the chapter we assume that the costs of the edges are
nonnegative, $c(u, v) \ge 0$, and that the heuristic is nonnegative,
$h(u) \ge 0$, as it is in all the domains of the library. The second
assumption is used only to have $h(t) \le h^*(t) = 0$ imply $h(t) = 0$
at a goal $t$ (a heuristic that is positive at a goal may make A$^*$
return a path that is not optimal, cf. [Exercise
2.6](02-best-first.md#sec-2-13)).

## 2.1 Searching a tree and searching a graph {#sec-2-1}

A best-first search keeps the nodes it knows in two sets: OPEN, the
*frontier*, i.e., the nodes generated but still unexpanded, and CLOSED,
the nodes already expanded. At each iteration it removes from OPEN a
node of least priority $f(u)$, tests whether it is a goal and, if it is
not, *expands* it, i.e., generates its successors by the domain (which
are then *generated*) and puts them in OPEN. A state is thus, at any
time, in one of three conditions: not generated yet, generated and
waiting in OPEN, or expanded and in CLOSED (the faint, light and darker
states of the figures of this chapter).

In its simplest form, a search of this kind does not look at what it has
already seen: each successor is put in OPEN as a new node, even if its
state has been generated before. Such a search explores the *tree of the
paths* from $s_0$, in which a state typically appears many times, once
for each path that reaches it, and it is called a *tree search*. [Figure
2.1](02-best-first.md#fig-2-1) shows a graph in which two paths reach
$c$, and the tree that a tree search explores on it: $c$ and $t$ appear
twice, and so does everything below them. On a graph with many short
cycles, e.g., a grid, the paths of length $d$ are exponentially many in
$d$, while the states within distance $d$ are only polynomially many,
and a tree search may even fail to terminate on a finite graph, since a
cycle gives paths of any length (cf. section 3.1 of
[[10]](Z-bibliography.md#ref-10)).

![Figure 2.1. A graph in which two paths reach $c$ (left), and the
tree of its paths from $s_0$, which a tree search explores (right): the
second copies of $c$ and $t$, in rust, are the *duplicates* that a
graph search detects in CLOSED and drops, since their $g$ is not
smaller than that of the first
copies.](../figures/best-first-tree-graph.svg){#fig-2-1}

A *graph search* keeps, instead, one node per state: it looks each
successor up among the states already generated (the *duplicate
detection*), and it drops the successor if its state is known with a $g$
that is not larger. This needs CLOSED to be a dictionary of the
generated states, and the library keeps one (the `NodeTable` of [Chapter
3](03-data-structures.md#ch-3)), but it bounds the nodes by the states.
One case remains: a successor $v$ whose state is known, but which has
now been reached through a cheaper path. If $v$ is still in OPEN, its
$g$ is lowered and its priority with it; if $v$ has already been
expanded, it is *reopened*, i.e., moved from CLOSED back to OPEN with
the new $g$, so that the cheaper path is propagated to its successors
when $v$ is expanded again. Whether a search must reopen depends on the
priority and, for A$^*$, on the heuristic (see
[§2.6](02-best-first.md#sec-2-6) and [§2.7](02-best-first.md#sec-2-7)).
Dropping a duplicate without comparing the values of $g$ is safe only
when the first copy of a state is known to be reached through a path of
least cost, e.g., in a breadth-first search with unit costs, or when an
optimal path is not required.

**Algorithm 2.1** (best-first graph search with a priority $f$).

```text
 1  procedure BEST-FIRST(s₀, f)
 2     CLOSED ← {s₀}, with g(s₀) ← 0 and parent(s₀) ← nil
 3     OPEN ← {s₀}
 4     while OPEN ≠ ∅ do
 5        u ← a node of OPEN of least f(u); remove it from OPEN
 6        if u is a goal then return the path to u, of cost g(u)
 7        mark u as expanded
 8        for each successor v of u do
 9           g' ← g(u) + c(u, v)
10           if v ∉ CLOSED then         ▷ a new state
11              add v to CLOSED, with g(v) ← g' and parent(v) ← u
12              put v in OPEN
13           else if g' < g(v) then     ▷ a cheaper path to v
14              g(v) ← g';  parent(v) ← u
15              put v in OPEN, whether it was expanded or not
16     return unsolvable
```

As in the library, CLOSED holds all the generated states, each with a
flag that records whether it has been expanded (line 7); the closed list
of the textbooks is the subset of CLOSED with the flag set, and OPEN is
a queue of the others. At line 6 the path is rebuilt by following the
pointers `parent` from the goal back to $s_0$ (the function `trace_path`
of `Node.h:126`). Note that the goal test is done when a node is
*selected* at line 5, rather than when it is generated at line 10. In
fact, a goal generated through a path of cost $g$ is not necessarily
reached optimally, since a node of OPEN of smaller priority may still
lead to a cheaper one; only the test at the selection gives the
guarantee of [Theorem 2.9](02-best-first.md#thm-2-9).

## 2.2 The priority {#sec-2-2}

Algorithm 2.1 is a family of searches, one for each priority $f$; the
classical ones are the following, with $d(u)$ the number of edges of
the path found to $u$ (cf. chapter 3 of
[[13]](Z-bibliography.md#ref-13) and chapter 2 of
[[10]](Z-bibliography.md#ref-10)).

* *Breadth-first search* (BrFS), $f(u) = d(u)$, with the ties broken
  first in, first out: the nodes are expanded by layers of depth, and
  the first goal selected is one with the fewest edges, which is
  optimal when all the costs are equal.
* *Depth-first search* (DFS), $f(u) = -d(u)$, with the ties broken last
  in, first out: the search follows one path as deep as it goes and
  backtracks only at a dead end; OPEN then behaves as a stack, and the
  recursive formulation does without it altogether, keeping only the
  current path (which is the basis of the memory-bounded searches of
  [Chapter 4](04-memory-bounded.md#ch-4)).
* *The algorithm of Dijkstra*, also called uniform-cost search,
  $f(u) = g(u)$: the nodes are expanded by nondecreasing cost from
  $s_0$, and the first goal selected is optimal for any nonnegative
  costs.
* *A$^*$*, $f(u) = g(u) + h(u)$, with $h$ the heuristic of the domain:
  $f(u)$ estimates the cost of a solution through $u$, and with an
  admissible $h$ the first goal selected is optimal
  ([§2.5](02-best-first.md#sec-2-5)). With $h = 0$, A$^*$ is the
  algorithm of Dijkstra.
* *Greedy best-first search*, $f(u) = h(u)$, which goes towards the
  goal that looks closest and is not optimal in general, and *weighted
  A$^*$*, $f(u) = g(u) + w\,h(u)$ with $w \ge 1$, which returns a path
  of cost at most $w\,C^*$ when $h$ is admissible.

For a graph with branching factor $b$ (the largest number of successors
of a state), an optimal solution of $d$ edges and cost $C^*$ (the
optimal cost of Table 1.2 of [§1.4](01-introduction.md#sec-1-4)), costs
of at least $\varepsilon > 0$, and largest depth $m$ of a state, the
first four have the following properties, the time and the space being
measured in nodes generated.

| search | complete | optimal | time and space |
|---|---|---|---|
| breadth-first | finite $b$ | equal costs | $O(b^{d+1})$ |
| depth-first, graph | finite graphs | no | $O(b^{m})$ |
| Dijkstra | finite $b$, $c\ge\varepsilon$ | yes | $O(b^{1+D})$ |
| A$^*$ | as Dijkstra | $h$ admissible | see the text |

Here $D = \lfloor C^* / \varepsilon \rfloor$, and the cost of A$^*$ is
[Proposition 2.19](02-best-first.md#prop-2-19). The rows of breadth-first
search and of the algorithm of Dijkstra are
[Proposition 2.20](02-best-first.md#prop-2-20), which derives them from
the theorems on A$^*$ with $h = 0$; the bound of breadth-first search
has the exponent $d + 1$ because the goal is tested when it is selected
(with the test at the generation, as in section 3.4 of
[[13]](Z-bibliography.md#ref-13), it would be $d$). The row of the
depth-first search of a graph, which the library does not implement, is
taken from section 3.4 of [[13]](Z-bibliography.md#ref-13): it
terminates on a finite graph since it expands each state at most once,
and it is not optimal since it returns the first goal that it reaches,
whatever the cost of its path.

With a consistent heuristic, the nodes expanded by A$^*$ are a subset of
those expanded by the algorithm of Dijkstra, but for the ties of the
last layer: by [Theorem 2.14](02-best-first.md#thm-2-14), A$^*$ expands
only nodes with $g^*(u) + h(u) \le C^*$, hence with $g^*(u) \le C^*$,
while the algorithm of Dijkstra expands all those with $g^*(u) < C^*$
([Theorem 2.15](02-best-first.md#thm-2-15) with $h_1 = 0$ gives the same
conclusion when $h$ is positive on the states that are not goals). Note
that the depth-first search of a tree needs only $O(b\,m)$ memory, the
current path and the siblings of its nodes, which is the reason for its
use in [Chapter 4](04-memory-bounded.md#ch-4); the depth-first search of
a graph keeps CLOSED and needs as much memory as the states it visits.
Consider the graph of five states $s_0, a, b, c, t$ drawn in [Figure
2.2](02-best-first.md#fig-2-2), whose seven edges $(s_0, a)$, $(s_0,
b)$, $(a, b)$, $(a, c)$, $(b, c)$, $(b, t)$ and $(c, t)$ carry their
costs and whose optimal path $s_0, a, b, c, t$ costs 7 (the figure is read in
[§2.3](02-best-first.md#sec-2-3)). Breadth-first search expands $s_0,
a, b, c$ and $t$ and returns the path $s_0, b, t$ of two edges and cost
10; depth-first search (with the successors pushed in the order $a, b$)
expands $s_0, b, t$ and returns the same path; the algorithm of Dijkstra
and A$^*$ both expand $s_0, a, b, c, t$ and return the optimal path.
That A$^*$ expands no fewer nodes than Dijkstra there is due to the size
of the example. On the grid of [Figure 2.6](02-best-first.md#fig-2-6), a
map of $9 \times 7$ cells with 4 moves of unit cost and a wall between
$s_0$ and $t$ (the figure is read in [§2.6](02-best-first.md#sec-2-6)),
a run of Algorithm 2.2 with $h = 0$ and the ties broken first in, first
out expands 56 cells, while A$^*$ with the Manhattan distance and the
ties broken in favor of the larger $g$, and last in, first out among
equal $f$ and $g$ (the rule of the bucket queue of the library),
expands 26 of them. These counts come from a separate simulation of
Algorithm 2.2 with these rules of the ties, not from a run of the
library, whose heap breaks the ties among equal $f$ and $g$ in no fixed
order.

Of these searches, the library implements A$^*$ (`AStar.h`), and
gets the algorithm of Dijkstra from it with a domain whose heuristic is
zero. A wrapper of a `Domain` that returns 0 from `heuristic` and passes
0 as the heuristic of each successor is all that one needs, and the
tests use one (`Blind`, `test/test.cpp:272`) to check A$^*$ on the grids
against the costs of Dijkstra (test `grid`). On a domain with unit costs
the same wrapper gives a breadth-first search, up to the order of the
ties within a layer.

> **Status — planned.** The library has no depth-first search of a
> graph, no greedy best-first search and no weighted A$^*$ as searches
> of their own. The depth-first searches of a tree are those of
> [Chapter 4](04-memory-bounded.md#ch-4); the bounded-suboptimal
> searches, weighted A$^*$ among them, are the subject of
> [Chapter 6](06-bounded-suboptimal.md#ch-6). A weighted A$^*$ with an
> integer weight is obtained today by a wrapper of the domain whose
> heuristic is $w\,h$, on which `AStar` runs unchanged (cf.
> [Exercise 2.7](02-best-first.md#sec-2-13)).

## 2.3 A* in the library {#sec-2-3}

### The interface

A search of the library takes a type `D` that models the concept
`Domain` of `Domain.h:103`. Such a type has the type `state_type` (a
copyable value compared with `==`), the arithmetic type `cost_type` of
$g$, $h$ and the costs, and the members `initial_state()`, `is_goal(s)`,
`heuristic(s)`, `hash(s)` (a well spread 64-bit hash, used by CLOSED)
and `for_each_successor(s, h_s, f)`. The last one calls `f(t, c, h_t)`
for each successor $t$ of $s$, with $c = c(s, t)$ and $h_t = h(t)$; the
value $h(s)$ is passed in, and hence the domain can compute $h(t)$
incrementally, as the sliding-tile puzzles do with the Manhattan
distance (cf. [Chapter 1](01-introduction.md#ch-1) and [Chapter
14](14-domains.md#ch-14)).

The class template `AStar<D, Open>` (`AStar.h:104`) is constructed from
the domain, which must outlive it, and from the `SearchOptions` of
`Solution.h:120`. These are `node_limit`, the number of stored nodes at
which the search gives up and returns `Status::node_limit` (0, the
default, for no limit; [Proposition 2.17](02-best-first.md#prop-2-17)
states when the test is made), and `expected_nodes`, the number of nodes
that CLOSED is sized for at the start, instead of growing by doublings
as it fills (0, the default, for a small initial table). Then `solve()`
runs the search and returns a `Solution` (`Solution.h:101`), i.e., the
`Status` (`solved`, `unsolvable` or `node_limit`), the `cost` of the
path, the `path` itself as the sequence of its states from $s_0$ to the
goal (empty unless the status is `solved`), and the counters `stats` of
[§2.10](02-best-first.md#sec-2-10); `solved()` tells whether the status
is `solved`. The OPEN list `Open` defaults to the bucket queue
(`BucketOpenList`) for integer costs and to the heap (`HeapOpenList`)
otherwise. The function `astar(domain, open, options)` (`AStar.h:240`)
chooses it at run time by an `OpenKind`, i.e., `automatic` (the default
of the class), `heap` or `buckets`; for a domain whose costs are not
integers it runs the heap whatever the `OpenKind`, since the bucket
queue needs integer keys. Both lists are described in [Chapter
3](03-data-structures.md#ch-3). Their declarations are the following
(quoted by hand, since in the header they come with their bodies).

```cpp
template< Domain D , class Open = /* buckets if cost_t< D > is integral,
                                     else the heap */ >
class AStar {
 public:
 explicit AStar( const D & domain , const SearchOptions & options = {} );
 [[nodiscard]] Solution< state_t< D > , cost_t< D > > solve( void ) const;
 };

template< Domain D >
[[nodiscard]] Solution< state_t< D > , cost_t< D > > astar(
 const D & domain , OpenKind open = OpenKind::automatic ,
 const SearchOptions & options = {} );
```

The following example, from the test case `manual` of the library
(`ctest -R PBeFS.manual`), solves an instance of the 8-puzzle with the
Manhattan distance by both lists, and then with a node limit.

```cpp
--8<-- "snippets/best-first-usage.cpp"
```

The examples of the repository (`examples/`) do the same on all the
domains of the library.

### The algorithm

In the library, the search is Algorithm 1 of the survey of Fukunaga et
al. [[1]](Z-bibliography.md#ref-1), i.e., Algorithm 2.1 with $f(u) =
g(u) + h(u)$. Two choices on the data structures change its form, while
its outcome stays the same. First, OPEN holds *entries* $(f, g, u)$
rather than nodes, ordered by increasing $f$ and, for equal $f$, by
decreasing $g$ ([§2.8](02-best-first.md#sec-2-8)). Second, neither OPEN
list of the library can lower the key of an entry: a node reached
through a cheaper path is pushed again with its new $g$, and the old
entry, which records the $g$ of its push, becomes *stale*, i.e., its
node is closed or has a different $g$. A stale entry is dropped when it
reaches the top of OPEN. We give the pseudocode of the Doxygen comment
of `AStar.h`.

**Algorithm 2.2** (A$^*$).

```text
 1  procedure A*(s₀)
 2     CLOSED ← {s₀}, with g(s₀) ← 0 and parent(s₀) ← nil
 3     OPEN ← {(f(s₀), g(s₀), s₀)}
 4     while OPEN ≠ ∅ do
 5        (f, g, u) ← the first entry of OPEN; remove it from OPEN
 6        if u is closed or g ≠ g(u) then continue      ▷ a stale entry
 7        close u
 8        if u is a goal then return the path to u, of cost g(u)
 9        for each successor v of u do
10           g' ← g(u) + c(u, v)
11           if v ∉ CLOSED then
12              add v to CLOSED, with g(v) ← g' and parent(v) ← u
13              push (g' + h(v), g', v) in OPEN
14           else if g' < g(v) then           ▷ a cheaper path: reopen v
15              g(v) ← g';  parent(v) ← u;  v is not closed
16              push (g' + h(v), g', v) in OPEN
17     return unsolvable
```

We now comment on the lines that are not obvious. Line 6 drops the stale
entries. An entry whose node is closed was pushed before the node was
expanded with a smaller (or equal) $g$, and an entry whose $g$ differs
from that of the node was overtaken by a cheaper path, whose entry is in
OPEN with a smaller $f$. In fact, the flag *closed* of line 7 records
that the node has been expanded *with its current $g$*; line 15 clears
it when a cheaper path is found, whether the node was waiting in OPEN or
already expanded, which is the reopening of
[§2.7](02-best-first.md#sec-2-7). Also, the node keeps a second flag,
*expanded*, which is never cleared and by which the expansions of a node
already expanded are counted (`Node.h:81`). With a consistent heuristic
the branch at lines 14 to 16 is taken only for nodes still in OPEN
([Theorem 2.13](02-best-first.md#thm-2-13)).

Algorithm 2.2 leaves out, as all the pseudocode of the book
([§1.4](01-introduction.md#sec-1-4)), the counters, the node limit and
the prefetches. Since the claims of [§2.10](02-best-first.md#sec-2-10)
are about the counters, we also give the loop of `AStar::solve()`
(`AStar.h:127`) with all of them. Its expansion is split in two passes:
the first asks the domain for all the successors of $u$, computes the
hash of each one and *prefetches*, i.e., asks the processor to start
loading into its cache, the slot of CLOSED where the lookup of that
successor will start; the second does the lookups, whose cache misses
thus overlap instead of being paid one after the other
([§3.8](03-data-structures.md#sec-3-8)).

**Algorithm 2.3** (A$^*$ as `AStar::solve()` runs it, with a node
limit $L$, $L = 0$ meaning none).

```text
 1  procedure SOLVE(s₀, L)
 2     CLOSED ← {s₀}, with g(s₀) ← 0, h(s₀) and parent(s₀) ← nil
 3     OPEN ← {(h(s₀), 0, s₀)}
 4     while OPEN ≠ ∅ do
 5        (f, g, u) ← the first entry of OPEN; remove it from OPEN
 6        if u is closed or g ≠ g(u) then continue      ▷ a stale entry
 7        close u;  expanded ← expanded + 1
 8        if u was expanded before then reexpanded ← reexpanded + 1
 9        mark u as expanded
10        if u is a goal then return (solved, g(u), the path to u)
11        K ← ∅                              ▷ pass 1: the successors
12        for each (v, c(u, v), h(v), o) of CHILDREN(u, h(u), op(u)) do
13           compute the hash of v and prefetch its slot of CLOSED
14           append (v, c(u, v), h(v), o) to K
15        generated ← generated + |K|
16        for each (v, c, h(v), o) of K, in order do  ▷ pass 2: lookups
17           g' ← g(u) + c
18           if v ∉ CLOSED then
19              add v to CLOSED, with g(v) ← g', parent(v) ← u, op(v) ← o
20              push (g' + h(v), g', v) in OPEN
21           else if g' < g(v) then
22              g(v) ← g';  parent(v) ← u;  op(v) ← o;  v is not closed
23              push (g' + h(v), g', v) in OPEN
24           else duplicates ← duplicates + 1
25        if L > 0 and |CLOSED| ≥ L then return (node_limit)
26     return (unsolvable)
```

Here CHILDREN is Algorithm 1.1 of [§1.4](01-introduction.md#sec-1-4),
which passes the operator $o$ of each successor and skips the move back
to the parent in a domain that knows its operators
([§2.9](02-best-first.md#sec-2-9)); `op(u)` is the operator stored in
the node, `no_op` for $s_0$. Line 19 is the lookup and the insertion of
CLOSED (`find_or_insert` of `NodeTable`, which builds the node only when
the state is new), and `stored` is $|$CLOSED$|$ when the search ends,
since CLOSED keeps one node per state. Also, line 7 sets the flag
*closed* and line 9 the flag *expanded* of `Node.h:80` and `Node.h:81`;
the path of line 10 is rebuilt by `trace_path` (`Node.h:126`).

[Proposition 2.17](02-best-first.md#prop-2-17) proves that Algorithm
2.3 does what Algorithm 2.2 does, and what the node limit changes.

### A trace

[Figure 2.2](02-best-first.md#fig-2-2) and [Figure
2.3](02-best-first.md#fig-2-3) follow Algorithm 2.2 on a graph of five
states with a consistent heuristic (the values $g/h$ are written next to
each state, in rust where they change, and the entries of OPEN are
listed in the order of the queue, the new ones in rust and the stale
ones struck out). In step (3) the expansion of $a$ finds a path to $b$
of cost 3, cheaper than the edge of cost 4 from $s_0$: $b$ is still in
OPEN, and the entry $(8, 4, b)$ becomes stale. In step (4) the goal $t$
is generated with $g = 9$, but it is not tested, and in step (5) the
expansion of $c$ lowers its $g$ to 7; in step (6) $t$ is selected and
the path $s_0, a, b, c, t$ of cost 7 is returned, while three stale
entries are still in OPEN and are never popped.

![Figure 2.2. Algorithm 2.2 on a graph of five states, steps (1) to (3),
the initial state $s_0$ being labelled $s$:
the node being expanded is violet, the closed ones darker, those not
generated yet faint; $g/h$ is written next to each state, in rust where
it changes. The entries $(f, g, u)$ of OPEN are in the order of the
queue, the new ones in rust, the stale ones struck out.
](../figures/best-first-trace-a.svg){#fig-2-2}

![Figure 2.3. Algorithm 2.2 on the graph of Figure 2.2, steps (4)
to (6): the goal is generated in step (4) with $g = 9$, improved in
step (5) to $g = 7$, and returned when it is selected in step (6); the
optimal path is drawn in violet.](../figures/best-first-trace-b.svg){#fig-2-3}

## 2.4 Admissible and consistent heuristics {#sec-2-4}

Much of the behavior of A$^*$ depends on two properties of $h$, which we now
define (the heuristics of the library and how they are obtained are
discussed in [Chapter 5](05-heuristics.md#ch-5)). Let $h^*(u)$ be the
least cost of a path from $u$ to a goal ($+\infty$ if there is none),
and $k(u, v)$ the least cost of a path from $u$ to $v$ ($+\infty$ if
there is none).

<a id="def-2-1"></a>**Definition 2.1** (Admissible heuristic). A
heuristic $h$ is *admissible* if $h(u) \le h^*(u)$ for each state $u$.

<a id="def-2-2"></a>**Definition 2.2** (Consistent heuristic). A
heuristic $h$ is *consistent* if $h(t) = 0$ for each goal $t$ and
$h(u) \le c(u, v) + h(v)$ for each edge $(u, v)$.

The inequality of Definition 2.2 is the triangle inequality of $h$
along the edges ([Figure 2.4](02-best-first.md#fig-2-4)): the estimate
from $u$ is never larger than the cost of an edge plus the estimate
from its head. Some textbooks define consistency over all the pairs of
states, with $k(u, v)$ in place of $c(u, v)$, and call the edgewise
property *monotonicity*; the two are equivalent, as the next lemma
shows, and the property of $f$ along a path that is sometimes also
called monotonicity is a consequence of either
([Lemma 2.12](02-best-first.md#lem-2-12)).

<a id="lem-2-3"></a>**Lemma 2.3** (Edgewise and pairwise consistency).
A heuristic satisfies $h(u) \le c(u, v) + h(v)$ for each edge
$(u, v)$ if and only if it satisfies $h(u) \le k(u, v) + h(v)$ for each
pair of states $u, v$.

*Proof.* Since $k(u, v) \le c(u, v)$ for each edge, the pairwise
inequality implies the edgewise one. Conversely, if $v$ is not
reachable from $u$, then $k(u, v) = +\infty$ and there is nothing to
prove; otherwise, let $u = w_0, w_1, \ldots, w_q = v$ be a path of
least cost. Summing the edgewise inequalities
$h(w_{i-1}) - h(w_i) \le c(w_{i-1}, w_i)$ for $i = 1, \ldots, q$, the
left-hand side telescopes to $h(u) - h(v)$ and the right-hand side is
$k(u, v)$. $\square$

<a id="lem-2-4"></a>**Lemma 2.4** (Consistency implies
admissibility). A consistent heuristic is admissible.

*Proof.* If no goal is reachable from $u$, then $h^*(u) = +\infty$.
Otherwise, let $t$ be a goal with $k(u, t) = h^*(u)$; by
[Lemma 2.3](02-best-first.md#lem-2-3),
$h(u) \le k(u, t) + h(t) = h^*(u)$, since $h(t) = 0$. $\square$

Of course, the converse does not hold. On the graph of
[Figure 2.7](02-best-first.md#fig-2-7), whose states $v_0, \ldots, v_4$
are joined by the arcs $(v_0, v_1)$, $(v_1, v_2)$, $(v_2, v_3)$ and
$(v_3, v_4)$ of cost 1 and by the arc $(v_0, v_2)$ of cost 3, each
directed from the state of the smaller index to that of the larger one
(the figure draws them as lines), with $s_0 = v_0$ and the goal $v_4$ (the
figure is read in [§2.7](02-best-first.md#sec-2-7)), the heuristic
$h = (0, 3, 0, 1, 0)$ is admissible, since $h^* = (4, 3, 2, 1, 0)$, but
$h(v_1) = 3 > c(v_1, v_2) + h(v_2) = 1$. The sequence score of Nilsson
(`SequenceScore.h`, used by the example `8Puzzle`) is a heuristic of
the sliding-tile puzzles that is not even admissible: A$^*$ then
usually expands fewer nodes, but it may return a path that is not
optimal (`SequenceScore.h:50`).

![Figure 2.4. Consistency (left): the estimate $h(u)$ is at most
the cost of the edge $(u, v)$ plus the estimate $h(v)$, i.e., $h$
satisfies the triangle inequality along the edges. Its consequence
(right): $f = g + h$ does not decrease along a path, which is
Lemma 2.12.](../figures/best-first-consistency.svg){#fig-2-4}

Both properties are preserved by the maximum of heuristics, which is
used by
the pattern databases of [Chapter 5](05-heuristics.md#ch-5) (the
maximum with the reflection of a table) and in general to combine
several heuristics of a domain.

<a id="lem-2-5"></a>**Lemma 2.5** (The maximum of heuristics). Let
$h_1, \ldots, h_m$ be heuristics and $h(u) = \max_i h_i(u)$. If all the
$h_i$ are admissible, so is $h$; if all are consistent, so is
$h$.

*Proof.* The first claim is immediate, since $h_i(u) \le h^*(u)$ for
each $i$. For the second, $h(t) = \max_i h_i(t) = 0$ at a goal; for an
edge $(u, v)$, let $i$ be an index with $h(u) = h_i(u)$: then
$h(u) = h_i(u) \le c(u, v) + h_i(v) \le c(u, v) + h(v)$. $\square$

We say that a heuristic $h_2$ is *more informed* than $h_1$ if $h_2(u) >
h_1(u)$ for each state $u$ that is not a goal; [Theorem
2.15](02-best-first.md#thm-2-15) below shows that, between two
consistent heuristics, the more informed one never makes A$^*$ expand a
node that the other does not expand.

## 2.5 The correctness of A* {#sec-2-5}

We now prove that A$^*$ with an admissible heuristic returns an optimal
path, and that it terminates. Our argument rests on one invariant on
OPEN (cf. Lemma 1 of [[9]](Z-bibliography.md#ref-9) and section 2.1 of
[[1]](Z-bibliography.md#ref-1)), which we state for any state $w$ (a
goal being a special case), since [§2.6](02-best-first.md#sec-2-6) needs
it in that form. Let $g^*(w)$ be the least cost of a path from $s_0$ to
$w$, and recall that $C^* = h^*(s_0)$.

<a id="lem-2-6"></a>**Lemma 2.6** (Values of $g$). At any time of
Algorithm 2.2, $g(u)$ is the cost of a path from $s_0$ to $u$, and
hence $g(u) \ge g^*(u)$, for each node $u$ in CLOSED.

*Proof.* By induction on the assignments to $g$: $g(s_0) = 0$ is the
cost of the empty path, and lines 12 and 15 set $g(v) = g(u) + c(u, v)$,
the cost of the path to $u$ of value $g(u)$ followed by the edge
$(u, v)$. $\square$

<a id="lem-2-7"></a>**Lemma 2.7** (The invariant of OPEN). Let $w$ be a
state reachable from $s_0$ that, at the start of an iteration of the
loop at line 4, is not closed with $g(w) = g^*(w)$. Then each path of
least cost $P$ from $s_0$ to $w$ has a node $u$ that is not closed, has
$g(u) = g^*(u)$ and has in OPEN the entry $(g^*(u) + h(u), g^*(u), u)$,
which is not stale.

*Proof.* Let $u$ be the first node of $P$ that is not closed with its
optimal $g$; it exists, since $w$ is such a node. If $u = s_0$, then
$s_0$ has not been expanded yet, and it is in OPEN with $g = 0$ by line
3. Otherwise the predecessor $p$ of $u$ on $P$ is closed with $g(p) =
g^*(p)$; since the prefix of $P$ to $u$ is a path of least cost, $g^*(u)
= g^*(p) + c(p, u)$. When $p$ was last expanded, with $g(p) = g^*(p)$,
line 10 computed $g' = g^*(u)$ for $u$. Then either $u$ was new, or it
had $g(u) > g'$, or it already had $g(u) \le g' = g^*(u)$. In the first
two cases lines 12 to 16 set $g(u) = g^*(u)$ and pushed the entry
$(g^*(u) + h(u), g^*(u), u)$, while in the third case $g(u) = g^*(u)$ by
[Lemma 2.6](02-best-first.md#lem-2-6) and such an entry had been pushed
when $g(u)$ took that value. Since $g(u)$ cannot fall below $g^*(u)$, it
has not changed since then, and by the choice of $u$ the node has not
been closed with it; hence the entry is still in OPEN, and it is valid.
$\square$

In other words ([Figure 2.5](02-best-first.md#fig-2-5)), the expanded
prefix of an optimal path is followed by one node of OPEN that carries
its optimal $g$. Its expansion generates the next node of the path with
its optimal $g$, which takes its place. Note that the proof needs the
reopening of lines 14 to 16: the next node might have been closed before
with a larger $g$, and without the reopening it would never be put back
in OPEN.

![Figure 2.5. The invariant of OPEN (Lemma 2.7): on an optimal path
from $s_0$ to a goal, the nodes before $u$ are closed with their
optimal $g$, and $u$ is in OPEN with $g(u) = g^*(u)$; if $h$ is
admissible, $f(u) \le C^*$.](../figures/best-first-invariant.svg){#fig-2-5}

<a id="lem-2-8"></a>**Lemma 2.8** (Every selected node has
$f \le C^*$). If $h$ is admissible and a goal is reachable from $s_0$,
each entry $(f, g, u)$ that is selected at line 5 and is not stale,
before a goal is returned, has $f \le C^*$, and therefore
$g^*(u) + h(u) \le C^*$.

*Proof.* Let $t^*$ be a goal with $g^*(t^*) = C^*$. Before the search
returns, $t^*$ is not closed with $g(t^*) = C^*$ (closing it would end
the search at line 8), and by [Lemma 2.7](02-best-first.md#lem-2-7) an
optimal path to $t^*$ has a node $u'$ with an entry of key
$f' = g^*(u') + h(u') \le g^*(u') + h^*(u') = C^*$ in OPEN, the
inequality being the admissibility of $h$ and the equality the
optimality of the path through $u'$. The entry selected at line 5 has
the least key in OPEN, and hence $f \le f' \le C^*$; finally,
$g^*(u) + h(u) \le g(u) + h(u) = f$ by
[Lemma 2.6](02-best-first.md#lem-2-6), since a selected entry that is
not stale has $g = g(u)$. $\square$

<a id="thm-2-9"></a>**Theorem 2.9** (Admissibility of A$^*$). If $h$
is admissible and A$^*$ returns a path, the path is optimal, i.e., its
cost is $C^*$.

*Proof.* Let the search return at line 8 with the goal $t$ selected
from a non-stale entry $(f, g(t), t)$. Since $h(t) \le h^*(t) = 0$ and
$h \ge 0$, $h(t) = 0$ and $f = g(t)$; by
[Lemma 2.8](02-best-first.md#lem-2-8), $g(t) \le C^*$, and by
[Lemma 2.6](02-best-first.md#lem-2-6), $g(t)$ is the cost of a path
from $s_0$ to $t$, and therefore $g(t) \ge C^*$. $\square$

It remains to prove that A$^*$ terminates. On a finite graph this
holds for any heuristic, admissible or not; on an infinite graph it
needs a positive lower bound on the costs, and the absence of cycles of
nonpositive cost is not enough.

<a id="lem-2-10"></a>**Lemma 2.10** (The parent pointers form a tree).
At any time of Algorithm 2.2, the pointers `parent` of the nodes of
CLOSED form a tree rooted at $s_0$, and
$g(\mathrm{parent}(v)) + c(\mathrm{parent}(v), v) \le g(v)$ for each
node $v \ne s_0$.

*Proof.* The inequality holds when `parent(v)` is set, with equality,
and is preserved afterwards, since $g$ never increases and a later
change of `parent(v)` sets it again. Suppose that the pointers form a
cycle $v_1, \ldots, v_q$, with $v_{i-1}$ the parent of $v_i$ (indices
modulo $q$), and let $v_j$ be the node of the cycle whose pointer was
set last, at line 15, to $u = v_{j-1}$. Just before that assignment the
other pointers of the cycle were already in place, and the inequality
along the path $v_j, v_{j+1}, \ldots, v_{j-1} = u$ gives
$g(u) \ge g(v_j)$, the costs being nonnegative; the new value
$g(u) + c(u, v_j)$ is then not smaller than the old $g(v_j)$, while
line 14 requires it to be strictly smaller. Hence there is no cycle,
and since each node but $s_0$ has a parent, the pointers form a tree
rooted at $s_0$. $\square$

<a id="thm-2-11"></a>**Theorem 2.11** (Termination of A$^*$). Algorithm
2.2 terminates (i) on any finite graph, and (ii) on any graph with
finite branching whose costs are at least some $\varepsilon > 0$,
provided that $h$ is admissible and that a goal is reachable from $s_0$.

*Proof.* Each iteration either drops a stale entry or closes a node, and
the entries are pushed only by lines 3, 13 and 16, once per value taken
by $g$ of a node. Hence, it suffices to show that the pairs $(v, g(v))$
taken by the nodes are finitely many. (i) By [Lemma
2.6](02-best-first.md#lem-2-6), each value of $g(v)$ is the cost of a
path from $s_0$ to $v$, and the successive values of $g(v)$ are strictly
decreasing. A path decomposes into a simple path and a multiset of
simple cycles, whose costs are nonnegative. Since the simple paths and
the simple cycles of a finite graph are finitely many, there is an $\eta
> 0$ not larger than the cost of each simple cycle of positive cost, and
a path whose cost is at most the first value of $g(v)$ contains at most
that value over $\eta$ cycles of positive cost. Thus the costs of such
paths are finitely many, and so are the values of $g(v)$. (ii) By [Lemma
2.8](02-best-first.md#lem-2-8), each node $u$ expanded before the end
has $g(u) \le f \le C^*$; by [Lemma 2.10](02-best-first.md#lem-2-10),
$g(u)$ is at least the cost of the path of the tree from $s_0$ to $u$,
which therefore has at most $C^* / \varepsilon$ edges. Under a finite
branching, the paths of at most $C^* / \varepsilon$ edges are finitely
many, and so are the nodes that may be expanded and their successors.
The values of $g$ are then finitely many by (i), applied to the finite
subgraph that these nodes span. $\square$

We remark that the condition of (ii) cannot be weakened to the absence
of cycles of nonpositive cost. On an infinite chain $s_0, w_1, w_2,
\ldots$ whose $i$-th edge costs $2^{-i}$, with a goal reached from $s_0$
by one edge of cost 1 and $h = 0$, each node of the chain has $g < 1$
and is selected before the goal, and the search never ends. On a finite
graph without a reachable goal, A$^*$ ends with OPEN empty and returns
`Status::unsolvable`, which the test `limits` checks on a grid split by
a wall; on an infinite graph it may end only by the node limit.

## 2.6 Consistent heuristics {#sec-2-6}

With a consistent heuristic A$^*$ has stronger properties, which we now
prove. The values of $f$ do not decrease along a path, and the nodes are
selected by nondecreasing $f$; each node is closed with its optimal $g$
and never reopened; and the set of the expanded nodes is determined by
$C^*$, up to the ties on the last value of $f$.

<a id="lem-2-12"></a>**Lemma 2.12** (Monotonicity of $f$). If $h$ is
consistent and $v$ is a successor of $u$ with
$g(v) = g(u) + c(u, v)$, then $f(v) \ge f(u)$.

*Proof.* $f(v) = g(u) + c(u, v) + h(v) \ge g(u) + h(u) = f(u)$, by
Definition 2.2. $\square$

<a id="thm-2-13"></a>**Theorem 2.13** (No reopening). If $h$ is
consistent, each node is closed by Algorithm 2.2 with
$g(u) = g^*(u)$; hence no closed node is ever reopened, each state is
expanded at most once, and the keys of the selected non-stale entries
are nondecreasing.

*Proof.* Suppose that a node $u$ is selected from a non-stale entry with
$g(u) > g^*(u)$. By [Lemma 2.7](02-best-first.md#lem-2-7) applied to $w
= u$, an optimal path to $u$ has a node $w'$ with an entry of key
$g^*(w') + h(w')$ in OPEN, and $w' \ne u$, since $g(u) \ne g^*(u)$. Now,
by [Lemma 2.3](02-best-first.md#lem-2-3), $h(w') \le k(w', u) + h(u)$,
and $g^*(w') + k(w', u) = g^*(u)$ since $w'$ is on an optimal path to
$u$; hence $g^*(w') + h(w') \le g^*(u) + h(u) < g(u) + h(u)$. Thus the
entry of $w'$ has a smaller key than that of $u$, and $u$ would not have
been selected. Each node is therefore closed with its optimal $g$, which
no path can improve, and line 14 never holds for a closed node; a state
is then expanded only when its node is closed, i.e., once. Finally, when
$u$ is expanded with key $f(u)$, all the keys in OPEN are at least
$f(u)$, and by [Lemma 2.12](02-best-first.md#lem-2-12) so are the keys
pushed by the expansion; hence the next selected key is not smaller.
$\square$

Thus, the keys of the selected entries grow by layers of equal $f$, and
$C^*$ is the last layer. Equivalently, A$^*$ is the algorithm of
Dijkstra on the reduced costs $c'(u, v) = c(u, v) - h(u) + h(v)$, which
consistency makes nonnegative. Indeed, along any path
$s_0 = w_0, \ldots, w_q = u$ the reduced costs telescope to
$\sum_i c(w_{i-1}, w_i) + h(u) - h(s_0)$: the reduced cost of a path is
its cost plus $h(u) - h(s_0)$, which depends only on its ends, the least
reduced cost from $s_0$ to $u$ is $g^*(u) + h(u) - h(s_0)$, and the key
$g(u) + h(u)$ of A$^*$ is the key of the algorithm of Dijkstra on $c'$
plus the constant $h(s_0)$, which orders the entries in the same way.
The next two results characterize the nodes
that A$^*$ expands ([Figure 2.6](02-best-first.md#fig-2-6)); the first
part holds for any admissible heuristic and is [Lemma
2.8](02-best-first.md#lem-2-8) again.

<a id="thm-2-14"></a>**Theorem 2.14** (Necessary and sufficient
conditions of expansion). Let $h$ be admissible and a goal be
reachable. (i) Every node $u$ expanded by Algorithm 2.2 before the goal
it returns has $g^*(u) + h(u) \le C^*$. (ii) If $h$ is consistent,
each node $u$ with $g^*(u) + h(u) < C^*$ is expanded before the search
returns.

*Proof.* (i) is [Lemma 2.8](02-best-first.md#lem-2-8). (ii) Suppose that
the search returns, selecting a goal $t$ with key $C^*$ (by [Theorem
2.9](02-best-first.md#thm-2-9)), while $u$ has not been expanded. By
[Lemma 2.7](02-best-first.md#lem-2-7), an optimal path to $u$ then has a
node $w$ with an entry of key $g^*(w) + h(w)$ in OPEN. As in the proof
of [Theorem 2.13](02-best-first.md#thm-2-13), this key is at most
$g^*(u) + h(u) < C^*$, and $t$, of key $C^*$, would not have been
selected before it. $\square$

A$^*$ with a consistent heuristic thus expands all the nodes with
$f^*(u) = g^*(u) + h(u) < C^*$, none with $f^*(u) > C^*$, and some of
those with $f^*(u) = C^*$, which ones depending on the order of the
ties ([§2.8](02-best-first.md#sec-2-8)). The comparison of two
heuristics follows.

<a id="thm-2-15"></a>**Theorem 2.15** (A more informed heuristic
expands no more nodes). Let $h_1$ and $h_2$ be consistent, and $h_2$
more informed than $h_1$. Every node expanded by A$^*$ with $h_2$, but
the goal it returns, is expanded by A$^*$ with $h_1$.

*Proof.* Let $u$ be a node that is not a goal and is expanded with $h_2$. By
[Theorem 2.14](02-best-first.md#thm-2-14) (i),
$g^*(u) + h_2(u) \le C^*$, and since $h_1(u) < h_2(u)$,
$g^*(u) + h_1(u) < C^*$; by [Theorem 2.14](02-best-first.md#thm-2-14)
(ii), $u$ is expanded with $h_1$. $\square$

If $h_2$ only dominates $h_1$, i.e., $h_2(u) \ge h_1(u)$, the same
argument covers the nodes with $g^*(u) + h_2(u) < C^*$, and the nodes
of the last layer may differ by the ties (cf. section 2.1 of
[[1]](Z-bibliography.md#ref-1)). Theorem 2 of
[[9]](Z-bibliography.md#ref-9) also shows that A$^*$ is optimal, in the
number of expanded nodes, among the admissible searches that use the
same information; the precise sense of that claim, which depends on
the ties, is discussed in chapter 2 of
[[10]](Z-bibliography.md#ref-10).

![Figure 2.6. The cells expanded on a grid with 4 moves of unit cost,
from $s_0$ to $t$ around a wall, by (a) the algorithm of Dijkstra (A$^*$
with $h = 0$, the ties broken first in, first out), and by A$^*$ with
the Manhattan distance, breaking the ties of $f$ (b) in favor of the
larger $g$ and (c) in favor of the smaller one. The colours classify the
expanded cells
by $f^* = g^* + h$, with $h$ the Manhattan distance ($C^* = 12$): both
runs of A$^*$ expand the 18 cells with $f^* < C^*$ and none with $f^* >
C^*$, but (b) expands 8 of the 31 cells with $f^* = C^*$, and (c) all of
them; the path returned is drawn as a dark line.
](../figures/best-first-layers.svg){#fig-2-6}

## 2.7 Reopening with an inconsistent heuristic {#sec-2-7}

When $h$ is admissible but not consistent, [Theorem
2.13](02-best-first.md#thm-2-13) fails: a node may be closed with a $g$
larger than $g^*$, and it is then reopened when a cheaper path to it is
found, as line 14 of Algorithm 2.2 allows. [Figure
2.7](02-best-first.md#fig-2-7) shows the graph of the test `reopen`
(`test/test.cpp:1006`), in which $h(v_1) = 3$ makes $v_1$ look expensive.
After $v_0$, A$^*$ expands $v_2$ (key 3) with $g(v_2) = 3$, through the
edge of cost 3, and then $v_1$ (key 4), which finds the path of cost 2
to $v_2$ and reopens it. Then $v_2$ is expanded a second time with key
2, smaller than the key 4 of the previous expansion (the keys of the
selected entries are thus not monotone), and the search returns the path
$v_0, v_1, v_2, v_3, v_4$ of cost 4, which is optimal, as [Theorem
2.9](02-best-first.md#thm-2-9) guarantees. The test checks the cost and
that the counter `reexpanded` is positive.

![Figure 2.7. A$^*$ with the admissible but inconsistent heuristic
$h = (0, 3, 0, 1, 0)$ on the graph of the test `reopen`, whose arcs go from
left to right. (a) $v_2$ is closed with $g = 3$; (b) the expansion of
$v_1$ finds $g(v_2) = 2$ and reopens $v_2$ (rust); (c) $v_2$ is
expanded again and lowers $g(v_3)$, making the old entry of $v_3$
stale.](../figures/best-first-reopen.svg){#fig-2-7}

Reopening keeps A$^*$ optimal, but it may cost a great deal: with an
inconsistent heuristic the number of expansions may grow much faster
than the number of states; a family of graphs on which it grows
exponentially in the number of states is discussed in chapter 2 of
[[10]](Z-bibliography.md#ref-10), and we do not reproduce it here. The
analysis of Felner et al.
[[31]](Z-bibliography.md#ref-31) discusses when inconsistent heuristics
are nevertheless useful. In the library a domain does not declare
whether its heuristic is consistent, and `AStar` always reopens. The
counter `reexpanded` of the `SearchStats` reports how many expansions
were of nodes already expanded; it is zero when the heuristic is
consistent (by [Theorem 2.13](02-best-first.md#thm-2-13)), and it is the
first quantity to look at when a search is slower than expected (cf.
[Chapter 5](05-heuristics.md#ch-5)). The pattern databases with the
blank minimized, for instance, are admissible but not consistent, and
A$^*$ reopens nodes on them. Note that a search that cannot reopen,
such as the search with delayed duplicate detection of [Chapter
4](04-memory-bounded.md#ch-4), refuses a successor whose $f$ is smaller
than that of its parent.

## 2.8 The ties within the last layer {#sec-2-8}

[Theorem 2.14](02-best-first.md#thm-2-14) leaves free the order of the
nodes with the same $f$. This order matters on the domains with few
distinct values of $f$, where the last layer $f = C^*$ may hold a large
part of the nodes generated. Both OPEN lists of the library break the
ties in favor of the larger $g$, i.e., of the node with the smaller $h$,
which is closer to a goal if the heuristic is informative. That on the
domains with unit costs the search then usually reaches a goal after
expanding a small part of the last layer is an empirical observation,
made in the experiments of Burns et al. [[4]](Z-bibliography.md#ref-4),
and not a theorem: an adversarial order of the successors can make
either rule expand the whole layer. On the grid of [Figure
2.6](02-best-first.md#fig-2-6), the simulation of
[§2.2](02-best-first.md#sec-2-2) expands 26 cells with the ties in favor
of the larger $g$, against 49 with the ties in favor of the smaller one;
the difference is entirely in the layer $f^* = C^*$ (8 of its 31 cells
against all of them), as [Theorem 2.14](02-best-first.md#thm-2-14)
requires. Among the entries with the same $f$ and $g$, the bucket queue
takes the last one pushed, which tends to follow one path down the
layer, while the heap takes any. Also, the library can be built with one
of three *sets of data structures*, i.e., three implementations of OPEN,
CLOSED and the arena of the nodes chosen at compile time (`Structures.h`,
[§3.12](03-data-structures.md#sec-3-12)); they break these last ties
differently, and hence they may expand different nodes of the last
layer, while by [Theorem 2.14](02-best-first.md#thm-2-14) they expand
the same nodes below it when $h$ is consistent.

## 2.9 Pruning the move back to the parent {#sec-2-9}

On a domain whose moves are reversible at the same cost, such as the
sliding-tile puzzles and the grids with symmetric costs, each
expansion generates the parent of the node it expands. Such a
successor is a duplicate that CLOSED would drop, since it is reached with a $g$
larger than that of the parent; but detecting it costs a hash and a
lookup, i.e., one of the about $b$ lookups of each expansion. A domain
that models the concept `OperatorDomain` (`Domain.h:217`) avoids it:
its `for_each_successor(s, h_s, op, f)` receives the operator `op`
that generated $s$ (`no_op` at the root), skips the successor that the
inverse of `op` leads to, and passes to `f` the operator of each
successor it generates, which the node stores in its field `op`
([Figure 2.8](02-best-first.md#fig-2-8), left). The function
`for_each_child` (`Domain.h:276`) calls this form when the domain has
it and the plain one otherwise, and all the searches of the library go
through it.

<a id="prop-2-16"></a>**Proposition 2.16** (Parent pruning changes no
outcome). If the inverse of each operator has the same cost, skipping
the successor that the inverse of the operator of $u$ leads to changes
neither the path returned by Algorithm 2.2 nor the order of its
expansions; it changes only the counters `generated` and `duplicates`
of [§2.10](02-best-first.md#sec-2-10), each by the number of the
successors skipped.

*Proof.* Let $p = \mathrm{parent}(u)$ be the state skipped, so that
$g(u) \ge g(p) + c(p, u)$ by [Lemma 2.10](02-best-first.md#lem-2-10).
The skipped successor would have
$g' = g(u) + c(u, p) \ge g(p) + 2\,c(p, u) \ge g(p)$, and $p$ is in
CLOSED; hence line 14 would not hold, and the successor would only be
counted as a duplicate, without any push. $\square$

A search that keeps no CLOSED, such as IDA$^*$, cannot detect the
duplicate at all, and without the pruning it would walk back and forth
along each edge; on a domain without operators it compares each
successor with the state two places back on its path
([Figure 2.8](02-best-first.md#fig-2-8), right, and
[Chapter 4](04-memory-bounded.md#ch-4)). The pruning by operator and
its use in the fast implementations of A$^*$ are discussed by Burns et
al. [[4]](Z-bibliography.md#ref-4).

![Figure 2.8. Parent pruning. (a) In an `OperatorDomain`, $u$
records the operator that generated it (here the blank moved down),
and its expansion does not generate the inverse move, so that $p$ is
never looked up. (b) Without operators, a search that keeps only the
current path skips the successor equal to the state before $u$ on the
path.](../figures/best-first-pruning.svg){#fig-2-8}

## 2.10 The cost of a search and its counters {#sec-2-10}

The `SearchStats` of a `Solution` (`Solution.h:65`) report what a search
did. For `AStar` the counters are the following, written by the lines of
Algorithm 2.3 given in parentheses.

* `expanded`, the expansions, the reexpansions included (line 7; an
  iteration that drops a stale entry is not an expansion);
* `reexpanded`, the expansions of nodes already expanded before
  (line 8, [§2.7](02-best-first.md#sec-2-7));
* `generated`, the successors generated, after the parent pruning
  (line 15);
* `duplicates`, the successors whose state was already in CLOSED with
  a $g$ not larger, and which were dropped (line 24);
* `stored`, the distinct nodes stored at the end, i.e., the states
  generated and the root, which measure the memory;
* `seconds`, the wall-clock time of `solve()`, from its first line to
  its return ([§2.11](02-best-first.md#sec-2-11)).

The other fields of `SearchStats` are for the parallel searches: `sent`
and `messages` count the nodes and the batches of nodes sent to another
worker, and `expanded_by` and `stored_by` hold the expansions and the
nodes of each worker, from which `load_balance()` and
`communication_overhead()` compute the load balance and the
communication overhead ([Chapter 8](08-shared-memory.md#ch-8) and
[Chapter 9](09-hash-distributed.md#ch-9)). For `AStar`, `sent` and
`messages` are 0, and the two vectors have the single entries
`expanded` and `stored`. We first prove that Algorithm 2.3, from which
these counters are read, is Algorithm 2.2, and then the relations
between them.

<a id="prop-2-17"></a>**Proposition 2.17** (Algorithm 2.3 is Algorithm
2.2). Let both algorithms run with the same OPEN list and the domain
enumerate the successors of each state in the same order. (i) With
$L = 0$, Algorithm 2.3 pops the same entries, closes the same nodes in
the same order and returns the same status, cost and path as Algorithm
2.2. (ii) With $L > 0$, it does the same until the first expansion after
which $|\mathrm{CLOSED}| \ge L$, and then returns `node_limit`; at that
moment $L \le |\mathrm{CLOSED}| \le \max(L - 1, 1) + b$, with $b$ the
largest number of successors of a state, and a goal that is selected is
always returned, since the test of line 10 precedes that of line 25.

*Proof.* (i) The counters and the hashes are written but never read by
the search, and a prefetch changes the content of the cache, not that
of CLOSED or of OPEN. Pass 1 calls only the domain, which reads neither
list; hence $K$ is the sequence of the successors that Algorithm 2.2
visits at its line 9, in the same order. During pass 2, $g(u)$ does not
change: an improvement of $u$ through one of its own successors would
need $g(u) + c(u, u) < g(u)$, i.e., a negative cost. Thus each
iteration of lines 16 to 24 does, on the same CLOSED, what lines 10 to
16 of Algorithm 2.2 do for the same successor (a state that occurs twice
in $K$ is found at its second occurrence where the first left it, in
both algorithms). Lines 19 and 22 also store the operator, which only
CHILDREN reads, and it reads it to skip a successor that Algorithm 2.2
would treat as a duplicate ([Proposition 2.16](02-best-first.md#prop-2-16)).
By induction on the iterations, the two algorithms go through the same
states of CLOSED and OPEN. (ii) The test of line 25 is made only after an
expansion, and before it the two algorithms agree by (i). When the
search stops, either the last expansion was the first one, made with
$|\mathrm{CLOSED}| = 1$, or the test failed after the expansion before
it, with $|\mathrm{CLOSED}| \le L - 1$, the iterations that drop a stale
entry adding no node; and an expansion adds at most $b$ nodes.
$\square$

Note that a search that stops at the limit may have found no goal only
because it stopped, and also that a search whose OPEN would have become
empty at the same iteration returns `node_limit` rather than
`unsolvable`.

<a id="prop-2-18"></a>**Proposition 2.18** (The counters of A$^*$). Let
$I$ be the number of the improvements of $g$ at line 22 of Algorithm
2.3. When `solve()` returns, (i) `generated` $=$ (`stored` $- 1$) $+ I
+$ `duplicates`; (ii) the entries pushed in OPEN are `stored` $+ I$, at
most $I$ of them are dropped as stale, and the iterations of the loop at
line 4 are at most `stored` $+ I$; (iii) `expanded` $-$ `reexpanded` is
the number of the distinct states expanded, which is at most `stored`;
(iv) if $h$ is consistent, `reexpanded` $= 0$ and the $I$ improvements
are all of nodes in OPEN.

*Proof.* (i) Each element of $K$ goes through exactly one of the
branches of lines 18, 21 and 24; the first adds a node to CLOSED, which
starts with $s_0$ and never loses a node, the second is an improvement
and the third a duplicate. (ii) The entries are pushed by lines 3, 20
and 23, i.e., $1 + (\mathtt{stored} - 1) + I$ times. Let $v$ be a node
improved $I_v$ times, so that it has $1 + I_v$ entries, pushed with
strictly decreasing values of $g$ (line 21); call *latest* the entry
pushed last at a given time, the only one whose $g$ equals $g(v)$. The
flag *closed* of $v$ is set only by the pop of its latest entry (line
7, after the test of line 6), and it is cleared by each later push
(line 22); hence the latest entry is never popped while $v$ is closed,
and its pop is not stale. A stale pop of $v$ is thus the pop of an entry
that is no longer the latest, and there are at most $I_v$ of these,
since the last entry pushed is always the latest. Summing over the
nodes, at most $I$ entries are stale; each iteration pops one entry, and
an entry is popped once. (iii) Line 7 counts each expansion and line 8
those of the nodes whose flag *expanded*, never cleared, was already
set; their difference is the number of the nodes on which the flag is
set, i.e., of the states expanded at least once, which are in CLOSED.
(iv) By [Theorem 2.13](02-best-first.md#thm-2-13), no closed node is
reopened, i.e., line 21 holds only for nodes that are not closed, and
each state is expanded at most once. $\square$

For instance, a search whose initial state is a goal returns at line
10 of its first iteration, with `expanded` $= 1$, `generated` $= 0$
and `stored` $= 1$, a path of one state and the cost 0 (the edge cases
checked by `check_astar_edges`, `test/test.cpp:400`). From the counters
one can derive the *effective branching factor* `generated / expanded`,
the fraction of the successors that are duplicates and the rate of
expansions per second; the number $I$ of the improvements is not
counted on its own, but it follows from (i).

The time of a best-first search is thus the number of its expansions
times the time of an expansion, which is made of (i) the extraction of
an entry from OPEN, (ii) the generation of the successors by the
domain, (iii) the evaluation of their heuristic, (iv) the duplicate
detection in CLOSED and (v) the insertions in CLOSED and OPEN. We state
the cost of `AStar` with the costs of the data structures of [Chapter
3](03-data-structures.md#ch-3) as parameters.

<a id="prop-2-19"></a>**Proposition 2.19** (The cost of A$^*$). Let
$X$, $\Gamma$ and $N$ be the values of `expanded`, `generated` and
`stored` at the end of a search, $I$ the improvements, $b$ the largest
number of successors of a state and $t_D$ a bound on the time that the
domain takes to generate one successor, with its heuristic and its hash.
Then Algorithm 2.3 takes time
$$
O\bigl(X + \Gamma\, t_D\bigr) + \Gamma\, t_{\mathrm{look}} +
(N + I)\,(t_{\mathrm{push}} + t_{\mathrm{pop}}) + O(N),
$$
with $t_{\mathrm{look}}$, $t_{\mathrm{push}}$ and $t_{\mathrm{pop}}$
the times of a lookup or insertion of CLOSED and of a push and a pop of
OPEN, and the last term the rebuilding of the path; it holds $N$ nodes,
at most $N + I$ entries of OPEN, the slots of CLOSED and $b$ successors
at a time. In particular: (i) with the heap of $d = 4$ children
([Proposition 3.3](03-data-structures.md#prop-3-3)) the term of OPEN is
$O((N + I) \log (N + I))$; (ii) with the bucket queue, integer costs and
a consistent $h$ ([Proposition 3.5](03-data-structures.md#prop-3-5),
under its hypotheses) it is $O(N + I + f_{\max})$, with $f_{\max}$ the
largest key pushed; (iii) with uniform hashing and the load factor of
CLOSED at most $7/10$ ([Proposition 3.11](03-data-structures.md#prop-3-11))
the expected value of $t_{\mathrm{look}}$ is $O(1)$, amortized over the
doublings; (iv) if $h$ is consistent and a goal is reachable, $I \le
\Gamma$, $X$ is at most the number of the states $u$ with $g^*(u) + h(u)
\le C^*$, and $N \le 1 + \Gamma \le 1 + b\,X$.

*Proof.* Each iteration of the loop pops one entry, and those that are
not stale run lines 7 to 25 once. Lines 11 to 15 call the domain once
per successor, i.e., $\Gamma$ times overall, besides an $O(1)$ work per
expansion; lines 16 to 24 make one lookup per successor and one push per
new node or improvement; the pops are at most the pushes, which are $N +
I$ by [Proposition 2.18](02-best-first.md#prop-2-18) (ii); and the path
of line 10 has at most $N$ states, since the pointers `parent` form a
tree on CLOSED ([Lemma 2.10](02-best-first.md#lem-2-10)). The memory
follows in the same way, the entries in OPEN at any time being at most
those pushed. Then (i), (ii) and (iii) are the propositions quoted. For
(iv), $I \le \Gamma$ by [Proposition 2.18](02-best-first.md#prop-2-18)
(i); by [Theorem 2.13](02-best-first.md#thm-2-13) each state is expanded
at most once, and by [Theorem 2.14](02-best-first.md#thm-2-14) (i) only
if $g^*(u) + h(u) \le C^*$; finally, $N - 1 \le \Gamma \le b\,X$ by
[Proposition 2.18](02-best-first.md#prop-2-18) (i). $\square$

On the sliding-tile puzzles and on the grids, $t_D$ is a few
instructions (the Manhattan distance of a successor is that of its
parent plus or minus one), and the time goes mostly to the data
structures, i.e., to $t_{\mathrm{look}}$, $t_{\mathrm{push}}$ and
$t_{\mathrm{pop}}$; with a costly heuristic, such as a pattern database
whose tables do not fit in the cache, $t_D$ may dominate instead
([Chapter 5](05-heuristics.md#ch-5)). Also, A$^*$ stores every state
that it generates, i.e., $N \ge X$; whether a search runs out of memory
before it runs out of time depends on the domain and on the machine,
and it is not a property that the bounds above decide.

The uninformed searches of [§2.2](02-best-first.md#sec-2-2) are A$^*$
with $h = 0$, and their properties follow from the results of this
chapter.

<a id="prop-2-20"></a>**Proposition 2.20** (The uninformed searches).
Let $h = 0$, so that Algorithm 2.2 is the algorithm of Dijkstra, and the
costs be nonnegative. (i) If it returns a path, the path is optimal; it
terminates on a finite graph, and on a graph of finite branching with
costs at least $\varepsilon > 0$ and a reachable goal. (ii) In the
latter case it expands each state at most once, only states with $g^*(u)
\le C^*$, and at most $\sum_{i = 0}^{D} b^i$ of them, $D = \lfloor C^* /
\varepsilon \rfloor$; it thus generates and stores
$O(b^{1 + D})$ nodes. (iii) With unit costs, the key $g(u)$ of a node
selected from a non-stale entry is the least number of edges of a path
from $s_0$ to $u$, the search is a breadth-first search up to the order
of the ties among the nodes of equal depth, and it generates and stores
$O(b^{d + 1})$ nodes, with $d = C^*$ the number of edges of an optimal
path.

*Proof.* The heuristic $h = 0$ is consistent (Definition 2.2: $0 \le
c(u, v) + 0$, and 0 at the goals), hence admissible by [Lemma
2.4](02-best-first.md#lem-2-4). (i) is [Theorem
2.9](02-best-first.md#thm-2-9) and [Theorem
2.11](02-best-first.md#thm-2-11). (ii) By [Theorem
2.13](02-best-first.md#thm-2-13) each state is expanded at most once,
with $g(u) = g^*(u)$, and by [Theorem 2.14](02-best-first.md#thm-2-14)
(i) only if $g^*(u) \le C^*$. By [Lemma 2.10](02-best-first.md#lem-2-10)
the path of the tree of the pointers `parent` from $s_0$ to $u$ costs at
most $g(u) \le C^*$, and since each of its edges costs at least
$\varepsilon$ it has at most $D$ edges. A node has at most $b$ children
in that tree, so that the nodes at depth $i$ of the tree are at most
$b^i$; the expanded ones are thus at most $\sum_{i \le D} b^i$, each
generates at most $b$ successors, and the stored nodes are at most one
more than the successors generated ([Proposition
2.18](02-best-first.md#prop-2-18) (i)). (iii) With unit costs, $g^*(u)$
is the least number of edges of a path to $u$, which is the key of $u$
when it is selected (by [Theorem 2.13](02-best-first.md#thm-2-13)), and
the keys selected are nondecreasing: the nodes are expanded by layers of
depth, which is breadth-first search up to the order within a layer.
The bound is (ii) with $\varepsilon = 1$ and $C^* = d$. $\square$

## 2.11 Benchmark {#sec-2-11}

Our benchmark of this chapter is the sequential A$^*$ on the 100
instances of the 15-puzzle of Korf [[7]](Z-bibliography.md#ref-7), with
the Manhattan distance and the goal with the blank in the first cell,
run by `tools/AStarBench` from `tools/bench/campaign.sh` (its table A).
Each instance is solved in a process of its own, pinned to one core and
to the memory of its node, first with the bucket queue, then with the
heap, and then with the bucket queue again; the second run of the bucket
queue is an A/A comparison, which measures the noise of the machine.
When the reference solver of Burns et al. [[4]](Z-bibliography.md#ref-4)
is available to the campaign, it is run on the same instances and its
row is added to the table. Costs are integers and the search is exact:
each instance is solved to its optimal cost, which is known from
[[7]](Z-bibliography.md#ref-7). The machine and the compiler are those
stated in the note of the table. The clock is the `Stopwatch` started at
the first line of `AStar::solve()` (`AStar.h:129`) and read before it
returns (`AStar.h:209`). Hence it does not include the start of the
process, the reading of the instance and the construction of the domain;
nor does it include the release of CLOSED and of the arena of the nodes,
which happens when the local objects of `solve()` are destroyed, after
the clock has been read, while it includes their allocation.

**Table 2.1.** Sequential A$^*$ on the 100 instances of Korf: for each
search (the OPEN list, and the reference solver when present), the
instances solved, the seconds summed over the instances (each instance
the median of its repetitions), the nodes expanded summed over the
instances, and the millions of expansions per second, i.e., the ratio
of the last two.

--8<-- "tables/best-first-korf.md"

## 2.12 Notes and references {#sec-2-12}

A$^*$ and the proof of its admissibility are due to Hart, Nilsson and
Raphael [[9]](Z-bibliography.md#ref-9), who also stated the consistency
assumption and, under it, the optimality of A$^*$ in the number of
expansions (their Theorem 2). The form of the algorithm with the
reopening of the closed nodes, which the library follows, is Algorithm 1
of the survey of Fukunaga et al. [[1]](Z-bibliography.md#ref-1). The
textbook of Edelkamp and Schrödl [[10]](Z-bibliography.md#ref-10) treats
in its chapters 2 and 3 the general best-first search, the reopening and
the worst cases of inconsistent heuristics, and chapter 3 of Russell and
Norvig [[13]](Z-bibliography.md#ref-13) compares the uninformed
searches. The textbook of Nilsson [[8]](Z-bibliography.md#ref-8) treats
the graph search and the heuristic estimates in the same framework. The
analysis of inconsistent heuristics by Felner et al.
[[31]](Z-bibliography.md#ref-31) shows that they can be useful in
practice, mostly in depth-first searches with the bidirectional pathmax,
which the library does not implement. The engineering of a fast A$^*$
(the bucket queue, the ties broken in favor of the larger $g$, the
parent pruning by operator and the states stored in place) is studied by
Burns et al. [[4]](Z-bibliography.md#ref-4), whose implementation is the
reference of the benchmark. The library makes the same choices, adds the
prefetch of the lookups of CLOSED and the huge pages ([Chapter
3](03-data-structures.md#ch-3)), and keeps one node per state in an
arena for the whole search. The parallel searches of [Chapter
8](08-shared-memory.md#ch-8) and [Chapter
9](09-hash-distributed.md#ch-9) keep the structure of Algorithm 2.2 in
each worker.

## 2.13 Exercises {#sec-2-13}

**Exercise 2.1** (`ctest -R PBeFS.astar`). Run Algorithm 2.2 by hand
on the graph of [Figure 2.2](02-best-first.md#fig-2-2) with the ties
broken in favor of the smaller $g$, and then in an arbitrary order;
show that the cost returned is the same, as
[Theorem 2.9](02-best-first.md#thm-2-9) requires. The test checks that
both OPEN lists, which break the ties differently, return the optimal
cost on 200 states of the 8-puzzle and on instances of the 15-puzzle.

**Exercise 2.2** (`ctest -R PBeFS.reopen`). On the graph of
[Figure 2.7](02-best-first.md#fig-2-7), list the entries of OPEN after
each iteration of Algorithm 2.2 and count the expansions and the
reexpansions. Which inequality of Definition 2.2 does
$h = (0, 3, 0, 1, 0)$ violate? Find the largest consistent heuristic
that is not larger than $h^*$ on that graph, and show that with it no
node is reopened.

**Exercise 2.3** (`ctest -R PBeFS.grid`). Show that A$^*$ on a domain
whose heuristic is zero is the algorithm of Dijkstra, and that the
wrapper `Blind` of `test/test.cpp:272` gives it from any domain. On a
grid with 4 moves of unit cost, compare the cells expanded with $h = 0$
and with the Manhattan distance (cf.
[Figure 2.6](02-best-first.md#fig-2-6)), and explain by
[Theorem 2.15](02-best-first.md#thm-2-15) why the second never expands
a cell that the first does not, but for the ties of the last layer.

**Exercise 2.4** (`ctest -R PBeFS.graph`). On the map of Romania of
Russell and Norvig [[13]](Z-bibliography.md#ref-13), with the
straight-line distance to Bucharest as the heuristic, run A$^*$ by hand
from Arad, and check that the path returned passes through Sibiu and
costs 418. Prove that the straight-line distance is consistent if the
length of each road is at least the distance between its ends.

**Exercise 2.5** (`ctest -R PBeFS.limits`). Show that a tree search
(Algorithm 2.1 without the lookups of line 10) may not terminate on a
finite graph with a reachable goal if some costs are zero, and that the
graph search of Algorithm 2.2 terminates on any finite graph
([Theorem 2.11](02-best-first.md#thm-2-11)). The test checks that
A$^*$ proves a grid split by a wall unsolvable, and that it stops at
the node limit.

**Exercise 2.6** (`ctest -R PBeFS.shift`). Let $h$ be consistent and
$\delta > 0$. (i) Is $h + \delta$ admissible? Does it satisfy the
inequality of Definition 2.2 on the edges? (ii) Show that A$^*$ with
$h' = h + \delta$ on all the states, goals included, selects the same
entries in the same order as with $h$, and returns an optimal path,
although $h'$ is not admissible. (iii) Show by an example that if $h$
is shifted on some states only, the path returned may not be optimal.
The test checks (ii) on 200 states of the 8-puzzle, comparing the
costs and the expansions with those of the Manhattan distance.

**Exercise 2.7** (`ctest -R PBeFS.weighted`). Let $h$ be admissible and
$w \ge 1$. Prove that A$^*$ with the priority $g + w\,h$ returns a path
of cost at most $w\,C^*$, following the proof of
[Lemma 2.8](02-best-first.md#lem-2-8) with the node $u'$ of the
invariant. The test runs `AStar` on a wrapper of the 8-puzzle whose
heuristic is twice the Manhattan distance, and checks the bound on 200
states.

**Exercise 2.8** (`ctest -R PBeFS.astar`). Give a graph and an
admissible but inconsistent heuristic on which $I > 0$ and `reexpanded`
$> 0$, and compute the four counters of [Proposition
2.18](02-best-first.md#prop-2-18) by hand; check (i) of that
proposition on them. Then show that a node limit $L = 5$ on a start at
least six moves from the goal of the 8-puzzle makes the search return
`node_limit` with between 5 and $\max(4, 1) + 4 = 8$ stored nodes
([Proposition 2.17](02-best-first.md#prop-2-17)). The test checks the
edge cases of a start that is a goal, a start one move away and the node
limit of 5 (`check_astar_edges`, `test/test.cpp:400`).

**Exercise 2.9** (`ctest -R PBeFS.reopen`). Prove that the pointers
`parent` of A$^*$ never form a cycle even when nodes are reopened
([Lemma 2.10](02-best-first.md#lem-2-10)), and that the path returned
on the graph of [Figure 2.7](02-best-first.md#fig-2-7) is the path of
the tree from $v_0$ to $v_4$. What would break in the proof if an edge
had a negative cost?
