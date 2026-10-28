# Best-first search and A*

We are concerned with finding a path of minimum cost from an initial state to a
goal state of a graph that is too large to be stored. Such a graph is given
implicitly, i.e., by its initial state, by a predicate that recognizes the
goals and by a procedure that generates the successors of a state with the
costs of the edges that lead to them (cf. chapter 2 of
[[10]](../../README.md#ref-10)). This page gives the notation shared by all the
pages of the theory, the best-first search and A* [[9]](../../README.md#ref-9),
the invariant on which the correctness of A* rests and the role of the
properties of the heuristic. Then it discusses the reopening of the closed
nodes and the breaking of the ties, and finally the cost of each operation on
the data structures of the library; the parallel searches of
[Parallel A*](parallel-astar.md) and the branch-and-bound of
[Branch-and-bound as best-first search](branch-and-bound-search.md) build upon
all of this. The pseudocode of each search is on the page of its header in the
API reference, with the conventions given at the end of this page.

## The problem and the notation

The search space is a directed graph $G = (V, E)$ with a cost $c(u, v) \geq 0$
on each edge $(u, v) \in E$, an initial state $s_0 \in V$ and a set of goals
$T \subseteq V$. The set $V$ may be infinite, but each state has finitely many
successors, which a type satisfying the concept `Domain` of
[Domain.h](../../include/PBeFS/astar/Domain.h) enumerates by
`for_each_successor`. The cost of a path is the sum of the costs of its edges;
we denote by $g^*(u)$ the least cost of a path from $s_0$ to $u$, by $h^*(u)$
the least cost of a path from $u$ to a goal ($+\infty$ if there is none), and
by $C^* = h^*(s_0)$ the optimal cost, which is the value sought.

A search generates *nodes*. A node records a state $u$, the cost $g(u)$ of the
best path from $s_0$ to $u$ found so far, the estimate $h(u)$ of $h^*(u)$ given
by the domain (the *heuristic*) and the predecessor of $u$ on that path, from
which the path is rebuilt at the end. The library keeps one node per generated
state (cf. `Node` in [Node.h](../../include/PBeFS/astar/Node.h)), and hence we
use node and state interchangeably. The value

$$
f(u) = g(u) + h(u)
$$

estimates the cost of a solution through $u$; note that $g(u) \geq g^*(u)$ at
any time, since $g(u)$ is the cost of an actual path.

The heuristic is *admissible* if $h(u) \leq h^*(u)$ for each state $u$, i.e.,
if it never overestimates the cost of reaching a goal. It is *consistent* if
$h(t) = 0$ for each goal $t$ and

$$
h(u) \leq c(u, v) + h(v) \qquad \text{for each edge } (u, v) \in E ,
$$

which is the triangle inequality of $h$ along the edges. A consistent heuristic
is admissible, as one sees by summing the inequality along an optimal path from
$u$ to a goal; the converse does not hold (e.g., the sequence score of Nilsson
used by the example `8Puzzle` is neither). On the sliding-tile puzzle the
Manhattan distance is consistent, since a move changes the distance of one tile
by one unit at the cost of one unit. On the grids the Manhattan distance (4
moves) or the octile distance (8 moves) times the least cost of a cell is
consistent for the same reason; these are the heuristics of
[SlidingTile.h](../../include/PBeFS/astar/SlidingTile.h) and of the grids. For
the traveling salesperson of
[TravelingSalesman.h](../../include/PBeFS/astar/TravelingSalesman.h), the
heuristic of a partial tour is the cost of a minimum spanning tree of the
unvisited cities together with the current and the first one (cf. exercise 3.30
of [[13]](../../README.md#ref-13)). It is admissible, since the rest of the
tour is a path through all these cities, and a path costs at least as much as a
spanning tree of its cities; the branch-and-bound of the benchmark uses the
same bound on the tree of the partial tours.

## Best-first search

A best-first search keeps two sets of nodes: OPEN, which holds the nodes
generated and not yet expanded (the frontier), and CLOSED, which holds the
nodes already expanded. At each iteration it removes from OPEN a node of least
priority and *expands* it, i.e., generates its successors; then it puts in OPEN
each successor that is new, or that is now reached through a cheaper path. With
the priority $f(u)$ the search is A*, with $g(u)$ (i.e., with $h = 0$) it is
the algorithm of Dijkstra, and with $h(u)$ it is the greedy best-first search,
which is not optimal.

In the library the table of the generated nodes plays the role of CLOSED (cf.
`NodeTable` in [NodeTable.h](../../include/PBeFS/astar/NodeTable.h)). Each
generated state has its node there, together with a flag that records whether
the node is *closed*, i.e., expanded with its current $g$; OPEN is a priority
queue of entries $(f, g, u)$ (cf.
[OpenList.h](../../include/PBeFS/astar/OpenList.h)). In A* the goal test is
done when a node is selected for expansion, rather than when it is generated.
In fact, a goal generated through a path of cost $g$ is not necessarily
optimal, since a node of OPEN with $f < g$ may still lead to a cheaper one; the
test at the selection makes the first goal selected an optimal one, as shown
below.

## Optimality and the invariants of A*

The correctness of A* rests on one invariant on OPEN (cf. Lemma 1 of
[[9]](../../README.md#ref-9) and section 2.1 of [[1]](../../README.md#ref-1)):
*before the search ends, if a goal is reachable from $s_0$, then each optimal
path $P$ from $s_0$ to a goal has a node $u$ in OPEN with $g(u) = g^*(u)$*.
Clearly, the invariant holds after the initialization, since $s_0$ is in OPEN
with $g(s_0) = 0$. Also, let $u$ be the first node of $P$ that is in OPEN with
its optimal $g$ (all its predecessors on $P$ having been expanded with their
optimal $g$). When $u$ is expanded, its successor $v$ on $P$ is generated with
$g(v) = g^*(u) + c(u, v) = g^*(v)$, and $v$ is put in OPEN unless it is already
there with that $g$; hence, $v$ takes the place of $u$. Note that the argument
needs the reopening of a closed node reached through a cheaper path, since $v$
might have been closed before with $g(v) > g^*(v)$.

If $h$ is admissible, the node $u$ of the invariant has
$f(u) = g^*(u) + h(u) \leq g^*(u) + h^*(u) = C^*$. Thus, when a goal $t$ is
selected, $g(t) = f(t) \leq f(u) \leq C^*$, and since $g(t)$ is the cost of an
actual path, $g(t) = C^*$: A* with an admissible heuristic returns an optimal
path as soon as it selects a goal (Theorem 1 of [[9]](../../README.md#ref-9)).
It also terminates whenever a goal is reachable, provided that the costs of the
edges are bounded below by some $\varepsilon > 0$, since each node selected has
$f \leq C^*$ and only finitely many paths have a cost of at most $C^*$. On a
finite graph it terminates in any case: a node is reopened only with a strictly
smaller $g$, which is the cost of one of finitely many simple paths, and the
search ends with an empty OPEN when no goal is reachable.

If $h$ is also consistent, then $f$ does not decrease along a path, since
$f(v) = g(u) + c(u, v) + h(v) \geq g(u) + h(u) = f(u)$ for a successor $v$ of
$u$. Therefore, the nodes are selected by nondecreasing $f$, and each node is
selected with $g(u) = g^*(u)$; hence, no closed node is ever reopened, and each
state is expanded at most once. Furthermore, A* expands all the nodes with
$f(u) < C^*$, some of those with $f(u) = C^*$ and none with $f(u) > C^*$ (cf.
section 2.1 of [[1]](../../README.md#ref-1)), and in this sense it is optimally
efficient among the admissible searches that use the same heuristic. With an
admissible but inconsistent heuristic the search remains optimal, but a state
may be expanded several times. Their number may then grow much faster than the
number of states (cf. chapter 2 of [[10]](../../README.md#ref-10)), and the
counter `reexpanded` of the `SearchStats` of a search reports them.

## Iterative deepening

The memory of A* grows with the nodes it generates, which on the 24-puzzle
exceeds that of any machine long before an optimal path is found; IDA*
([IDAStar.h](../../include/PBeFS/astar/IDAStar.h)), the iterative-deepening A*
of [[7]](../../README.md#ref-7), keeps only the current path. Each iteration is
a depth-first search from $s_0$ that cuts every node with $g(u) + h(u) >
\beta$, the bound $\beta$ being $h(s_0)$ in the first iteration and, in each of
the others, the least $f$ that exceeded the bound of the previous one; the
search ends when an iteration reaches a goal, or with no goal reachable when no
node was cut. With $h$ admissible the bound never exceeds $C^*$, since the node
$u$ of the invariant of the previous section, on an optimal path, has $f(u)
\leq C^*$ and is cut only if $\beta < C^*$; hence, the goal reached first is
reached with $g = \beta \leq C^*$, i.e., at an optimal cost (cf.
[[7]](../../README.md#ref-7)). Of course, a depth-first search detects no
duplicates: a state reached by several paths is expanded once per path in each
iteration, and the library prunes only the move back to the parent, by its
operator in an `OperatorDomain` or else by comparing the successor with the
state two places back on the path. Thus, IDA* suits the domains whose graph has
few short cycles (such as the sliding-tile puzzles, where with a pattern
database it is the method of [[22]](../../README.md#ref-22) for the 24-puzzle),
while on the grids, whose cells are reached by very many paths, its expansions
grow exponentially with the length of the path. We remark that the costs must
be positive on every cycle, or an iteration might not end, that without a
reachable goal the search ends only on an acyclic graph (on one with cycles the
bound grows without end, and only the node limit stops it), and that with
integer costs the number of iterations is at most $C^* - h(s_0) + 1$, since the
bound grows by at least 1 at each iteration and never exceeds $C^*$, which on
the sliding-tile puzzles, where $f$ changes by 2 at each move, is $(C^* -
h(s_0)) / 2 + 1$.

## The breaking of the ties

The invariant leaves free the order of the nodes with the same $f$. This order
matters on the domains with few distinct values of $f$, where the last layer
$f = C^*$ may hold a large part of the nodes. The library breaks the ties in
favor of the larger $g$, i.e., of the node with the smaller $h$, which is
closer to a goal if the heuristic is informative. On the domains with unit
costs the search then usually reaches a goal after expanding a small part of
the last layer (cf. [[4]](../../README.md#ref-4)). Among the nodes with the
same $f$ and $g$ the bucket queue takes the last one pushed (i.e., it behaves
as a stack, which tends to follow one path down), while the heap takes any of
them. Note that the sets of data structures of
[Structures.h](../../include/PBeFS/astar/Structures.h) break these last ties
differently; hence, they may expand different nodes of the last layer, while
they expand the same nodes below it.

## The cost of the operations

An iteration of A* costs one extraction from OPEN, the expansion of the node
(which the domain pays) and, for each of the $b$ successors, one evaluation of
the hash of the state, one lookup in CLOSED and possibly one insertion in OPEN.
Neither OPEN list supports the decrease of a key: a node reached through a
cheaper path is pushed again with its new $g$. The old entry, which records the
$g$ of its push, becomes *stale*, i.e., its node is closed or has a different
$g$, and it is dropped when it reaches the top (line 6 of the pseudocode of A*
in [AStar.h](../../include/PBeFS/astar/AStar.h)). Thus, the entries are as many
as the pushes, i.e., the nodes plus the improvements of $g$, which with a
consistent heuristic are at most the generated nodes. This costs less than
keeping the position of each node in the heap, which each swap of the heap
would have to update.

For any cost type, OPEN is a $d$-ary heap with $d = 4$ (`HeapOpenList`), whose
items keep $f$ and $g$ next to the pointer to the node (so that a comparison
does not touch the node). With $N$ entries, a push costs $O(\log_d N)$
comparisons and a pop $O(d \log_d N)$. With $d = 4$ the heap has half the depth
of a binary heap, and the four children of an item lie in one cache line for
integer costs (four items of 16 bytes).

For integer costs, the default OPEN is a two-level bucket queue
(`BucketOpenList`) in the manner of the bucket queue of Dial: an array of
layers indexed by $f$, each one an array of stacks indexed by $g$, which
[[4]](../../README.md#ref-4) reports to be faster than a binary heap on the
sliding-tile puzzles. A push costs $O(1)$ amortized time (the arrays grow by
doubling), and a pop costs $O(1)$ plus the empty buckets that it skips. When
the least $f$ in OPEN never decreases (as with a consistent heuristic), the
skipped layers are at most $f_{\max}$ over the whole search, with $f_{\max}$
the largest $f$ pushed, and the skipped stacks of a layer are at most its
largest $g$. Of course, the memory of the buckets grows with the range of the
values of $f$ and $g$ rather than with the entries; this is why the bucket
queue is used only for integer costs, and why the stacks of a layer are
released when the layer is emptied. The range of $g$ is "small" on the
sliding-tile puzzles ($g \le 80$ on the 15-puzzle), but on a grid $g$ grows to
thousands, while a layer holds few entries; in a worker of HDA*, whose least
layer empties and fills again each time the nodes sent by the other workers
arrive, rebuilding the array of a layer then costs $O(g)$ each time, and it
was measured to take most of the time of the search. Hence, the queue keeps
the arrays only as long as every $g$ pushed is below `DenseG` $= 256$: at the
first larger $g$, every layer becomes a binary heap ordered by $g$ and, for
equal $g$, by the order of the pushes, so that the order of the pops does not
change, a push and a pop cost $O(\log k)$ time with $k$ the entries of the
layer, and the memory of a layer is proportional to its entries.

CLOSED is a hash table with open addressing and linear probing. Its slots hold
the pointer to a node together with the 64-bit hash of its state, so that a
probe compares two states only when their hashes are equal, and the table grows
without evaluating any hash again. The slot of a hash value $x$ in a table of
$2^k$ slots is

$$
\left\lfloor \frac{(x \cdot \lfloor 2^{64} / \phi \rfloor) \bmod 2^{64}}
                   {2^{64 - k}} \right\rfloor ,
\qquad \phi = \frac{1 + \sqrt{5}}{2} ,
$$

i.e., the top $k$ bits of the product of $x$ by the integer part of
$2^{64} / \phi$ (Fibonacci hashing). Thus, the low bits of $x$, which a
distribution such as $x \bmod p$ makes equal on the states of a worker, do not
decide the slot. In the library the load factor $\alpha$ is kept below $7/10$.
By the classical analysis of linear probing under uniform hashing, a successful
lookup inspects about $\frac{1}{2}(1 + \frac{1}{1 - \alpha})$ slots and an
unsuccessful one about $\frac{1}{2}(1 + \frac{1}{(1 - \alpha)^2})$, i.e., at
most about 2.2 and 6.1 slots, which are consecutive in memory; the growth by
doubling costs $O(1)$ amortized time per insertion.

Hence, besides the expansion, an iteration costs $O(b + d \log_d N)$ expected
time with the heap, and $O(b)$ amortized time with the bucket queue (plus the
skipped buckets, and with a factor $\log k$ once its layers are heaps). In practice the time of a large search is dominated by the
cache misses of the lookups, since the slots of the successors lie at random
addresses of a table much larger than the cache. For this reason, the
sequential search computes the hashes of all the successors first and
prefetches their slots, so that the misses overlap, and the table and the arena
of the nodes are allocated on huge pages, which reduce the misses of the
translation lookaside buffer. Also, the arena of the nodes is never compacted,
and a pointer to a node (in OPEN, in CLOSED, or as the parent of another node)
stays valid until the end of the search. Overall, the memory of the search is
$O(N)$ nodes and entries, which is usually (but not always) its limit, since A*
typically stores each node that it generates.

## The conventions of the pseudocode

The pseudocode on the page of each header follows the code of the library, with
the following conventions (the same for A* and for branch-and-bound).

- The lines are numbered, a block is given by its indentation, and the keywords
  are `procedure`, `while ... do`, `for each ... do`, `if ... then`, `else`,
  `loop`, `return`, `continue` and `break`; `▷` starts a comment.
- `x ← e` assigns, `∞` is the infinite cost and `∅` the empty set, and the
  relations are `=`, `≠`, `<`, `≤`, `≥` and `∈`.
- OPEN and CLOSED are the lists of the search and `s₀` is the initial state;
  `g(u)`, `h(u)`, `f(u) = g(u) + h(u)` and `parent(u)` are the fields of the
  node of the state `u`, and `c(u, v)` is the cost of the edge $(u, v)$. Also,
  `z` is the incumbent, i.e., the cost of the best goal found so far ($+\infty$
  if none), and `p` is the number of workers, indexed by `i` from 0 to $p - 1$.
- `lock L` and `unlock L` take and release the lock `L`, and `send m to j` puts
  the message `m` in the outbox of the worker `j`, which may deliver it later,
  in a batch.
- The node limits, the counters of the statistics and the prefetches are left
  out, since they change neither the outcome nor the order of the search.
