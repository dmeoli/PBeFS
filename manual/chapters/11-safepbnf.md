# 11. Parallel best-nblock-first search {#ch-11}

[Source: `include/PBeFS/astar/SafePBNF.h` (`SafePBNF`, `PBNFOptions`),
`include/PBeFS/astar/Domain.h` (`NBlockDomain`), and the nblocks of
`include/PBeFS/astar/Grid.h`, `SlidingTile.h` and `VacuumWorld.h`]

We are concerned in this chapter with a parallel A$^*$ in which the
workers never lock the lists of the search while they expand nodes: the
Parallel Best-NBlock-First search (PBNF) of Burns, Lemons, Ruml and Zhou
[[3]](Z-bibliography.md#ref-3), in the livelock-free form that the
authors call Safe PBNF, which is the one of the library. Whereas the
searches of [Chapter 8](08-shared-memory.md#ch-8) share one OPEN list
under a lock, and HDA$^*$ ([Chapter 9](09-hash-distributed.md#ch-9))
sends each generated node to the worker that owns it, PBNF partitions
the state space by an abstraction into *nblocks*, each with an OPEN and
a CLOSED list, and lets a worker expand the nodes of an nblock only when
no other worker can touch the lists that these expansions may write (we
say that the worker holds the nblock). Hence, the synchronization is
confined to a graph of counters (the abstract graph of the nblocks),
which is usually (but not always) much smaller than the state space and
which the workers lock only to choose the next nblock. The price is
twofold: the order of expansion is best-first only approximately, which
often costs some search overhead (cf. [Chapter
8](08-shared-memory.md#ch-8)), and the `Domain` has to supply an
abstraction of the problem.

The idea comes from the Structured Duplicate Detection (SDD) of Zhou and
Hansen, developed for the search on external memory and turned into the
Parallel Structured Duplicate Detection (PSDD) for breadth-first
heuristic search (both summarized in section 2.4 of
[[3]](Z-bibliography.md#ref-3), pp. 694-696). We describe PSDD, with its
synchronization by layers, in [Chapter
12](12-more-parallel-schemes.md#ch-12). Here we use the notation of
[Chapter 1](01-introduction.md#ch-1). A$^*$ and its invariants are in
[Chapter 2](02-best-first.md#ch-2), the OPEN and CLOSED lists of a
single nblock are those of [Chapter 3](03-data-structures.md#ch-3), and
the overheads by which a parallel search is measured, as well as the
general condition under which a parallel search may stop, are in
[Chapter 8](08-shared-memory.md#ch-8).

The structure of the chapter is as follows. Section 11.1 states the
interface, i.e., the concept `NBlockDomain`, the abstractions of the
domains of the library and the options of the search. Section 11.2 gives
the duplicate detection scope and the interference between nblocks, on
which the search is built, and section 11.3 the search PBNF. Section
11.4 shows how plain PBNF may livelock and how the hot nblocks of Safe
PBNF prevent it. In section 11.5 we give the pseudocode of the library
line by line, and in section 11.6 we prove its properties (exclusive
access, absence of deadlock, optimality, termination, what can be
proved of the livelock freedom, and the detection of the end). Section
11.7 describes the implementation (the creation of the nblocks, the free
list and the expansion, in pseudocode) with the cost of its operations
and its memory, section 11.8 the benchmark, and sections 11.9 and 11.10
give the references and the exercises.

## 11.1 The problem and the interface {#sec-11-1}

We consider the problem of [Chapter 1](01-introduction.md#ch-1): a
`Domain` gives the initial state $s_0$, the successors of a state $u$
with the costs $c(u, v)$ of the edges, the goal test and an admissible
heuristic $h$, and the search returns a path of least cost $C^*$ from
$s_0$ to a goal, or proves that none exists. Throughout the chapter the
costs are positive, and a statement that needs them bounded away from
zero states it. As in the other parallel searches
([§8.3](08-shared-memory.md#sec-8-3)), the workers share the
*incumbent* $z$, i.e., the least $g$ of a goal selected so far
($z = +\infty$ before the first one), and they prune the nodes with
$f \ge z$. PBNF needs, besides this, an abstraction of the state space,
which we now define.

<a id="def-11-1"></a>**Definition 11.1** (Abstraction and nblocks). An
*abstraction* is a map $\phi$ from the states to a finite set
$B = \{0, \ldots, N - 1\}$ of *abstract states*, together with a
relation $\mathrm{succ}(b) \subseteq B \setminus \{b\}$, the *abstract
successors* of $b$, such that (i) for each edge $(u, v)$ of the state
space, either $\phi(v) = \phi(u)$ or
$\phi(v) \in \mathrm{succ}(\phi(u))$, and (ii) the relation is
symmetric, i.e., $a \in \mathrm{succ}(b)$ if and only if
$b \in \mathrm{succ}(a)$. The *nblock* $b$ is the set of the states $u$
with $\phi(u) = b$ (and, by the usual abuse, the abstract state $b$
itself), and the *abstract graph* has the nodes $B$ and the edges
$\{a, b\}$ with $a \in \mathrm{succ}(b)$.

Condition (i) means that the abstraction is a homomorphism of the state
space onto the abstract graph with its self-loops, i.e., that an edge
never jumps over an nblock. Burns et al. state it as "if two states are
successors in the state space, then their images are successors in the
abstract graph" ([[3]](Z-bibliography.md#ref-3), section 2.4, p. 694),
with a self-loop on an nblock when one of its states has a successor in
it. Condition (ii) is not required by Burns et al., whose appendix A
also covers directed abstract graphs. We assume it in the library, since
`SafePBNF` builds the interference scope of an nblock from its
successors alone (the comment "undirected abstract graphs: the
predecessors are the successors", `SafePBNF.h:333`). Both conditions
hold for the domains of the library, and we check neither of them; a
`Domain` that breaks them gives a search with data races (see Exercise
11.9).

![Figure 11.1. An abstraction of a $6 \times 6$ grid of states into
nine nblocks of $2 \times 2$ states (left, the shaded squares) and the
abstract graph of the nine nblocks (right). The edge of the state space
from the
state in the primary colour to the state in the accent colour crosses
from the nblock A1 to the nblock B1, and it is mapped onto the abstract
edge between them, as condition (i) of Definition 11.1
requires.](../figures/safepbnf-abstraction.svg){#fig-11-1}

In the library, a `Domain` provides an abstraction by modelling the
concept `NBlockDomain` (`Domain.h:182`), i.e., by three members:
`num_nblocks()`, the number $N$ of nblocks; `nblock(s)`, the nblock
$\phi(s) \in \{0, \ldots, N - 1\}$ of a state; and
`for_each_nblock_successor(b, f)`, which calls `f(a)` for each
$a \in \mathrm{succ}(b)$, without $b$ (the search adds $b$ to the scope,
see below). Figure 11.1 gives the picture on a grid; in the library,
three domains model the concept, as follows (one can of course add the
three members to any other `Domain` that has a "natural" abstraction).

* `Grid<Moves>` (`Grid.h:376-400`): the nblocks are the blocks of
  $s \times s$ cells, $s$ being the argument `nblock_side` of the
  constructor (20 by default, and 20 in `AStarBench`), and the abstract
  successors of a block are the 4 blocks next to it with 4 moves, and
  the 8 blocks around it with 8 moves, since a diagonal move may cross a
  corner of a block.
* `SlidingTile<W, H>` (`SlidingTile.h:374-412`): the nblock of a board
  is the triple of the cells of the blank and of the tiles 1 and 2,
  which is the abstraction used by Burns et al. for the 15-puzzle
  ([[3]](Z-bibliography.md#ref-3), section 4.4, p. 712); its ids are
  $(p_0 \cdot WH + p_1) \cdot WH + p_2$, i.e., $16^3 = 4096$ ids for the
  15-puzzle, of which only the $16 \cdot 15 \cdot 14 = 3360$ with three
  distinct cells are abstract states (the others have no successors and
  are never reached). A move of the blank to a cell that holds the tile
  1 or 2 swaps two cells of the triple, and any other move changes only
  the cell of the blank.
* `VacuumWorld` (`VacuumWorld.h:176-196`): a row of $n$ cells with a
  robot that moves left or right or sucks the dirt of its cell. The
  nblock of a state is the cell of the robot, and the abstract graph is
  the path $0 - 1 - \cdots - (n - 1)$, which is the abstract graph of
the figures of this chapter from Figure 11.3 on (where, for
readability, the nblocks are numbered from 1).

In the library, the search is the class template `SafePBNF`
(`SafePBNF.h:200-205`), whose parameters are the `Domain` (a model of
`NBlockDomain`), the executor that runs the workers (the default thread
pool unless another one is given; the executors are described in
[§1.6](01-introduction.md#sec-1-6) and
[§8.1](08-shared-memory.md#sec-8-1)) and the OPEN list of the nblocks
(the bucket queue for integer costs and the binary heap otherwise,
[Chapter 3](03-data-structures.md#ch-3)). Its constructor takes the
`Domain` (which must outlive the search) and a `PBNFOptions`, and its
`solve()` returns the `Solution` of [Chapter 2](02-best-first.md#ch-2),
with the status `solved`, `unsolvable` or `node_limit` and the counters
of the search. Each of the $p$ workers runs on a separate thread of the
executor (e.g., on the threads of C++ or on OpenMP, the two executors
that the tests and the benchmark use). One can lend an executor that
already exists by `set_executor()`; without it, `solve()` builds one of
$p$ threads. The members that a user calls are declared as follows (the
declarations of `SafePBNF.h`, without their bodies; `SolutionT` is
`Solution< state_t< D > , cost_t< D > >`):

```cpp
explicit SafePBNF( const D & domain , const PBNFOptions & options = {} );
void set_executor( Executor & ex ) noexcept;
[[nodiscard]] SolutionT solve( void ) const;
```

`PBNFOptions` (`SafePBNF.h:155-160`) derives from the `ParallelOptions`
of the other parallel searches ([§8.1](08-shared-memory.md#sec-8-1)) and
adds one field,

```cpp
--8<-- "snippets/safepbnf-sig-options.cpp"
```

whose meaning, with that of the inherited fields, is the following:

* `min_expansions`, the number $m$ of expansions that a worker makes in
  an nblock before it may switch to another one (32 by default; 0 is
  taken as 1, `SafePBNF.h:225-226`);
* of the inherited fields, `threads` is the number $p$ of workers (0 is
  taken as 1), and `node_limit` stops the search with the status
  `node_limit` when a worker has allocated that many nodes; the fields
  `batch`, `flush_period` and `seed`, which serve the messages and the
  owners of HDA$^*$, and `expected_nodes` are not used, since the closed
  lists are per nblock and grow from a "small" size (16 slots).

The following example, from the case `manual` of the tester, which the
build compiles and runs (`ctest -R PBeFS.manual`), solves a random grid
of $300 \times 200$ cells with 4 moves, whose nblocks are squares of
$10 \times 10$ cells (the last argument of the constructor of the
`Grid`), with $p = 4$ workers and $m = 32$ (the default, set explicitly
in the options `o`); the case then checks the `Solution` `r` against the
cost of the sequential A$^*$ on the same grid.

```cpp
--8<-- "snippets/safepbnf-usage.cpp"
```

The declaration needs only the `Domain` type, since the executor and the
OPEN list have defaults, and the `Solution` is that of every search of
the library.

## 11.2 Structured duplicate detection {#sec-11-2}

When a worker expands a node $u$ of the nblock $b$, the children may
fall only in $b$ or in an abstract successor of $b$, by condition (i) of
Definition 11.1. Hence, a worker that has exclusive access to the lists
of these nblocks can check each child against the CLOSED list of its
nblock and insert it in the OPEN list of its nblock with no lock at all.
This is the observation behind SDD, which the following definitions make
precise.

<a id="def-11-2"></a>**Definition 11.2** (Scopes and interference). The
*duplicate detection scope* of the nblock $b$ is
$S(b) = \{b\} \cup \mathrm{succ}(b)$. Two distinct nblocks $a$ and $b$
*interfere* if $S(a) \cap S(b) \neq \emptyset$, and the *interference
scope* of $b$ is the set $I(b) = \{a \in B \setminus \{b\} : S(a) \cap
S(b) \neq \emptyset\}$ of the nblocks that interfere with $b$.

Note that the definition of Burns et al. reads $\mathrm{IntScope}(b) =
\{l : \exists x \in \mathrm{Successors}(b),\ l \in
\mathrm{Predecessors}(x)\}$ ([[3]](Z-bibliography.md#ref-3), section
3.1, p. 698), with the self-loops of the abstract graph among the
successors. On an undirected abstract graph it is the set above plus
$b$, which plays no role. The following two facts are immediate, and the
rest of the chapter uses them throughout.

<a id="lem-11-3"></a>**Lemma 11.3** (The children stay in the scope). If
$u$ is a state of the nblock $b$ and $v$ is a successor of $u$, then
$\phi(v) \in S(b)$.

*Proof.* By condition (i) of Definition 11.1, $\phi(v) = b$ or
$\phi(v) \in \mathrm{succ}(b)$, i.e., $\phi(v) \in S(b)$. $\square$

<a id="lem-11-4"></a>**Lemma 11.4** (Form and symmetry of the interference
scope). For each nblock $b$,
$$I(b) = \Big(\mathrm{succ}(b) \cup \bigcup_{a \in \mathrm{succ}(b)}
\mathrm{succ}(a)\Big) \setminus \{b\},$$
and $a \in I(b)$ if and only if $b \in I(a)$. In particular,
$\mathrm{succ}(b) \subseteq I(b)$.

*Proof.* We obtain the symmetry from that of the condition
$S(a) \cap S(b) \neq \emptyset$. For the form, let $a \neq b$. If
$a \in \mathrm{succ}(b)$, then $a \in S(a) \cap S(b)$; if
$a \in \mathrm{succ}(x)$ for some $x \in \mathrm{succ}(b)$, then
$x \in S(a) \cap S(b)$. Conversely, let $x \in S(a) \cap S(b)$: if
$x = b$, then $b \in \mathrm{succ}(a)$, and hence
$a \in \mathrm{succ}(b)$ by condition (ii); if $x = a$, then
$a \in \mathrm{succ}(b)$; otherwise, $x \in \mathrm{succ}(b)$ and
$x \in \mathrm{succ}(a)$, and therefore $a \in \mathrm{succ}(x)$ by
condition (ii). $\square$

We remark that this is how `SafePBNF` builds the interference scope of
an nblock the first time that the nblock is touched: the successors of
$b$ and their successors, sorted, without repetitions and without $b$
(`SafePBNF.h:322-350`, Algorithm 11.4 in section 11.7). Figure 11.2
gives the scopes on a grid of nblocks. With 4 moves, $S(b)$ is the cross
of the nblock and its four neighbours and $I(b)$ the diamond of the 12
nblocks within Manhattan distance 2. With 8 moves, $S(b)$ is the
$3 \times 3$ square and $I(b)$ the 24 other nblocks of the $5 \times 5$
square around $b$ (fewer on the border of the grid);
[Lemma 11.19](#lem-11-19) proves these counts and the general bound
$|I(b)| \le d^2$ for an abstract graph of largest degree $d$. Yet, the
interference scope is usually a tiny part of the abstract graph, e.g.,
24 nblocks out of the 62 500 of the benchmark of section 11.8.

![Figure 11.2. The scopes of an nblock $b$ (dark) in the middle of a
grid of nblocks, with 4 moves (left) and with 8 moves (right): the other
nblocks of the duplicate detection scope $S(b)$ (light violet), which
belong to the interference scope $I(b)$ as well, and the rest of $I(b)$
(accent colour). While $b$ is in use, none of the coloured nblocks can
be acquired.](../figures/safepbnf-scopes.svg){#fig-11-2}

Figure 11.3 gives the interference on the path of the vacuum world. The
scopes of the nblocks 3 and 5 share the nblock 4, and a worker expanding
the nodes of 3 and another expanding those of 5 could both write the
lists of 4. On the contrary, the scopes of 3 and 6 are disjoint, and the
two nblocks can be searched in parallel.

![Figure 11.3. Interference on the abstract graph of the vacuum
world, a path. (a) The scopes of the nblocks 3 and 5 share the nblock 4:
the two interfere, and at most one of them can be in use. (b) The scopes
of the nblocks 3 and 6 are disjoint: two workers may hold them at the
same time, and they expand their nodes with no
lock.](../figures/safepbnf-interference.svg){#fig-11-3}

A worker that holds the nblock $b$ thus expands a node $u$ of $b$ much
as the sequential A$^*$ does (cf. [Chapter 2](02-best-first.md#ch-2)),
except that each child $v$ is looked up in the CLOSED list of $\phi(v)$
and inserted in the OPEN list of $\phi(v)$, both of which lie in $S(b)$
by Lemma 11.3. In the expansion of Figure 11.4, one of the three
children of $u$ falls in the nblock 2, one in the nblock 3 of $u$, and
the third is a state already closed in the nblock 4 with a path that is
not longer, which is dropped as a duplicate.

![Figure 11.4. An expansion by the worker that holds the nblock 3
(in the colour of worker 0). (a) The node $u$, the first of OPEN(3), is
selected. (b) Its children $c_1$ and $c_2$ go into OPEN(2) and OPEN(3),
and the third child, the closed state $w$ of the nblock 4, is found in
CLOSED(4) with a value of $g$ that is not larger, and it is dropped. All
three lists lie in $S(3)$, and the worker writes them with no
lock.](../figures/safepbnf-expand.svg){#fig-11-4}

PSDD (section 2.4 of [[3]](Z-bibliography.md#ref-3)) uses these scopes
to expand in parallel the nodes of one layer of a breadth-first
heuristic search. In each layer, the workers acquire nblocks whose
scopes are disjoint, expand their nodes of the current depth, and
synchronize all together before the next layer. PSDD and the variants of
Burns et al. that order the layers by $f$ (BFPSDD) or deepen them
iteratively (IDPSDD) are described in [Chapter
12](12-more-parallel-schemes.md#ch-12). PBNF keeps the scopes and drops
the layers, whose synchronization is often "costly" when they are small.
Indeed, Burns et al. observe that BFPSDD "may incur excessive
synchronization overhead if there are few nodes in each f layer"
(section 2.4.1, p. 696).

## 11.3 Parallel best-nblock-first search {#sec-11-3}

PBNF approximates the order of A$^*$ by a heap of nblocks. Each worker
expands the nodes of the nblock with the best OPEN list among those that
can be taken, and keeps it while it is not worse than the others
(section 3 of [[3]](Z-bibliography.md#ref-3), pp. 696-697, which borrows
the heap of sets of nodes from the Localized A$^*$ of Edelkamp and
Schrödl). Besides the lists of the nblocks, we call the state of the
search the following.

<a id="def-11-5"></a>**Definition 11.5** (The nblock graph). At each
instant, an nblock is *in use* if a worker holds it. For each nblock
$b$, $f(b)$ is the least $f$ of the nodes in OPEN($b$) ($+\infty$ if
none), $\sigma(b)$ is the number of the nblocks in use whose
interference scope contains $b$, and, in Safe PBNF, $\sigma_h(b)$ is the
number of the *hot* nblocks (section 11.4) whose interference scope
contains $b$ (in plain PBNF, $\sigma_h \equiv 0$). The nblock $b$ is
*free* if it is not in use, $\sigma(b) = \sigma_h(b) = 0$, and
$f(b) < z$, where $z$ is the incumbent. The *free list* holds the free
nblocks ordered by $f(b)$.

The counters, the hot flags and the free list form the *nblock graph*,
which one lock protects; its critical sections are rather short, since
each of them typically touches only the interference scopes of a few
nblocks. A worker that holds the nblock $b$ searches it as follows. It
takes the best node of OPEN($b$), stops if its $f$ is not smaller than
$z$, records it as the incumbent if it is a goal, and otherwise expands
it as in Figure 11.4. After each $m$ expansions, it *checks* whether it
should switch, i.e., whether the best free nblock or the best nblock of
$I(b)$ is better than $b$. If it should, it locks the graph, *releases*
$b$ (it decrements $\sigma$ on $I(b)$ and puts in the free list the
nblocks that become free) and *acquires* the first nblock $b'$ of the
free list (it takes $b'$ out of the free list and increments $\sigma$ on
$I(b')$, whose nblocks leave the free list). Figures 11.5 and 11.6 give
the two operations on the path of the vacuum world, with $f(b)$ above
each nblock and $\sigma$ and $\sigma_h$ below it.

![Figure 11.5. ACQUIRE on the path of the vacuum world, with 10
nblocks. (a) Worker 0 holds the nblock 3, and hence $\sigma = 1$ on
$I(3) = \{1, 2, 4, 5\}$, and the free list holds the nblocks 6 to 10,
ordered by $f$ (each cell gives the nblock and its $f$). (b) Worker 1
acquires 7, the first of the free list: $\sigma$ grows by one on
$I(7) = \{5, 6, 8, 9\}$, and 6, 8 and 9 leave the free list. The nblock
5, which interferes with both, has $\sigma(5) = 2$, and only 10 stays
free.](../figures/safepbnf-acquire.svg){#fig-11-5}

![Figure 11.6. RELEASE on the state reached in Figure 11.5, after
some expansions in the nblock 3. (a) Workers 0 and 1 hold 3 and 7. (b)
Worker 0 releases 3, and $\sigma$ drops by one on $I(3)$: 1, 2 and 4
become free, and so does 3, while 5 stays out, since $\sigma(5) = 1$
because of the nblock 7. The new entries of the free list (in the accent
colour) keep the order of $f$, ties broken by the smaller
nblock.](../figures/safepbnf-release.svg){#fig-11-6}

Three choices of Burns et al. reduce the overhead of this scheme, and we
keep all of them. First, the minimum number $m$ of expansions between
two checks bounds the frequency of the lock of the graph, at the price
of "speculative" expansions, i.e., of nodes that may be worse than the
best node of another nblock (section 3, p. 697). Their tuning on twenty
$5000 \times 5000$ grids with 62 500 nblocks found that the time spent
on the lock decreased as $m$ grew from 1 to 64, and that the benefit
levelled off at about $m = 8$, while the number of expanded nodes grew
with $m$ and with the number of threads (section 4.1 and Figure 4, pp.
702-703). In their implementation without instruments, somewhat larger
values of $m$ gave the best times, and $m = 32$, the default of the
library, is the value of their experiment on the size of the abstraction
(Figure 5, p. 704). Second, a worker whose nblock is not empty does not
wait for the lock when it would switch: it tries the lock and, if
another worker holds it, it usually goes on with its nblock (line 25 of
Algorithm 11.2, in section 11.5). Third, the first goal found is not
necessarily optimal
(with several workers it may well not be), and the search goes on until
no OPEN list holds a node better than the incumbent. It ends when no
nblock is in use and the free list is empty (line 32), which section
11.6 proves to be the right condition.

## 11.4 Livelock and the hot nblocks {#sec-11-4}

In plain PBNF a worker chooses its next nblock greedily among the free
ones, and nothing forces the nblocks of $I(b)$ to be out of use at the
same time. Hence, the nblock $b$ that holds the best nodes may never
become free. Burns et al. found such an execution with a model checker
([[3]](Z-bibliography.md#ref-3), section 3.2.3, p. 699), and the
following proposition gives one explicitly.

<a id="prop-11-6"></a>**Proposition 11.6** (Plain PBNF may livelock).
There are an infinite state space with an admissible and consistent
heuristic, an abstraction, and an execution of plain PBNF with two
workers in which a node of least $f$ is never expanded and the search
never ends, although a goal is reachable.

*Proof.* We take the state space of Figure 11.7 with $h \equiv 0$ (and
hence $f = g$), the abstract graph being the path of the nblocks
$1, \ldots, 9$: the root $s_0$ in the nblock 4, the states $v$, $t$ and
the goal $\gamma$ in 5, $x$ in 6, and two infinite chains
$u_1, u_2, \ldots$ in 3 and $y_1, y_2, \ldots$ in 7, with the costs of
the figure. Condition (i) of Definition 11.1 holds, since each edge
stays in an nblock or joins two adjacent ones. We let $m = 1$ and
consider the following execution. Worker 1 acquires 4 (the only free
nblock) and expands $s_0$, which puts $u_1$ into OPEN(3) and $v$ into
OPEN(5). It releases 4 and acquires 5, the better of the free 3 and 5,
while worker 0 finds no free nblock (3 is in $I(5)$) and waits. Worker 1
expands $v$, which puts $x$ into OPEN(6), releases 5 and acquires 6, and
worker 0, woken by the release, acquires 3, whose scope $\{2, 3, 4\}$ is
disjoint from $S(6) = \{5, 6, 7\}$. Worker 1 expands $x$, which puts $t$
(with $f = 5$) into OPEN(5) and $y_1$ into OPEN(7). It then releases 6
and acquires 7, since 5 is not free ($\sigma(5) = 1$ because of the
nblock 3) and 7 is. From now on $f(5) = 5$, while $f(3), f(7) \ge 12$
grow by one at each expansion. At each of its checks, worker 0 sees the
nblock 5 of $I(3)$ better than its own, releases 3, which brings
$\sigma(5)$ from 2 to 1, and, in the same critical section, acquires the
first free nblock, which is 3 again (the nblocks 1, 2 and 4 are empty,
and 5 is not free). Worker 1 does the same with 7 (Figure 11.8). Since
each release is followed by an acquisition under the same lock,
$\sigma(5) \ge 1$ at every instant, the nblock 5 is never free, $t$ is
never expanded, and the chains are expanded forever. $\square$

![Figure 11.7. The state space of Proposition 11.6, drawn in the
columns of its nblocks (the path $1 - 2 - \cdots - 9$; the empty nblocks
1 and 9 are not drawn). With $h \equiv 0$, $f = g$ is written next to
each state. The node $t$ (in the primary colour) has the least $f$ among
the nodes generated after $x$, and the goal $\gamma$ is reached only
through it. The chains of the nblocks 3 and 7 are
infinite.](../figures/safepbnf-livelock-graph.svg){#fig-11-7}

![Figure 11.8. The livelock of plain PBNF on the state space of
Figure 11.7, from the instant at which worker 1 acquires 7 ("--" is an
empty OPEN list). (a) Worker 0 holds 3 and worker 1 holds 7, and
$\sigma(5) = 2$. (b) Worker 0 switches, since $f(5) < f(3)$, but the
release of 3 and the acquisition of the first free nblock, which is 3
again, happen under one lock, and $\sigma(5)$ never reaches 0. (c)
Worker 1 does the same with 7, and the state is that of (a) with larger
values of $f$.](../figures/safepbnf-livelock.svg){#fig-11-8}

Note that one worker alone cannot livelock in this way, since when it
releases its nblock $\sigma$ drops to 0 on the whole interference scope.
Also, the execution above needs a particular interleaving of the two
workers (which a scheduler may, or may not, produce). We remark that a
finite state space cannot livelock (Theorem 11.13 below), although it
may well pay some search overhead for the same reason, since the "best"
nblock may stay blocked for rather long.

Safe PBNF (section 3.1 of [[3]](Z-bibliography.md#ref-3), pp. 697-699)
prevents the livelock by what they call the hot nblocks. At its check, a
worker that holds $b$ looks also at the best nblock $b'$ of $I(b)$. If
$b'$ is better than $b$ and better than the best free nblock, the worker
sets $b'$ *hot*, i.e., it increments $\sigma_h$ on $I(b')$, and no
nblock of $I(b')$ is free any more and none can be acquired. Then it
releases $b$, and each other worker that holds an nblock of $I(b')$
releases it at its next check. Hence $\sigma(b')$ can only decrease, and
it eventually reaches 0. When it does, $b'$ is set cold and freed, with
the least key of the free list. To keep this from blocking too much of
the graph, at most one nblock of an interference scope is hot. Indeed,
an nblock is not set hot if a better hot nblock interferes with it, and
setting it hot cools the worse hot nblocks that interfere with it. These
are the three cases of Burns et al. for undirected abstract graphs
(section 3.1, p. 698), whose two further cases concern directed graphs,
which the library does not support. Also, a worker that keeps its nblock
at a check, because no nblock of $I(b)$ is better than $b$, cools the
hot nblocks of $I(b)$, which would otherwise wait for an nblock that is
not worse than them.

![Figure 11.9. SETHOT on the state of Figure 11.8 (a). (a) Worker 0,
at its check, finds $f(5) < f(3)$, and the free list is empty. (b) It
sets 5 hot: $\sigma_h = 1$ on $I(5) = \{3, 4, 6, 7\}$, and none of these
nblocks can be acquired, even when it is not in use, until 5 is cold
again.](../figures/safepbnf-sethot.svg){#fig-11-9}

![Figure 11.10. Safe PBNF after Figure 11.9 (b). (a) Worker 0
releases 3, which is not free ($\sigma_h(3) = 1$). The free list is
empty, and worker 0 waits. (b) At its check, worker 1 releases 7:
$\sigma(5)$ drops to 0, 5 is set cold, $\sigma_h$ drops to 0 on $I(5)$,
and 5, 3 and 7 enter the free list, 5 first. (c) Worker 1 acquires 5,
and the nodes $t$ and $\gamma$ of Figure 11.7 are expanded
next.](../figures/safepbnf-safe.svg){#fig-11-10}

In Figures 11.9 and 11.10 we run Safe PBNF from the state of Figure 11.8
(a). The first check of worker 0 sets 5 hot, and two releases later
worker 1 acquires 5. Section 11.6 proves what holds in general.

## 11.5 The pseudocode {#sec-11-5}

The pseudocode of the library (the Doxygen comment of `SafePBNF.h`)
follows appendix A of [[3]](Z-bibliography.md#ref-3) (pp. 737-739), with
the conventions of [§1.4](01-introduction.md#sec-1-4). We split it here
into three parts, which keep the numbers of its lines. The nblock graph
is protected by one lock, $m$ is `min_expansions`, and $f(b)$ of an
nblock in use is the value that its worker publishes after each
expansion.

**Algorithm 11.1** (Safe PBNF: the worker).

```text
 1  procedure THREADSEARCH()                        ▷ on each worker
 2     b ← nil
 3     loop
 4        b ← NEXTNBLOCK(b);  if b = nil then return
 5        k ← 0
 6        while not SHOULDSWITCH(b, k) do
 7           u ← the first node of OPEN(b), removed (stale entries dropped)
 8           if f(u) ≥ z then OPEN(b) ← ∅;  break
 9           close u
10           if u is a goal then z ← min(z, g(u)), recording u if g(u) < z
11           else for each successor v of u with g(u) + c(u,v) + h(v) < z do
12              insert v in the lists of its nblock, as A* does
13           k ← k + 1;  publish f(b)
```

Before the first call, the root $s_0$ is inserted in OPEN($\phi(s_0)$)
and its nblock is put in the free list (`SafePBNF.h:248-255`). Line 7
drops the stale entries of OPEN, i.e., those whose node has been closed
or has received a smaller $g$ since the entry was pushed (the lazy
deletion of [§3.1](03-data-structures.md#sec-3-1)). Line 8 empties
OPEN($b$) when its first
node is not better than the incumbent, since then no node of the list is
(e.g., after a goal of cost close to the optimum has been found, for an
nblock that holds only worse nodes). We make the goal test at the
selection (line 10), as in the other searches of the library, under a
separate lock that protects the incumbent and its node. Line 11 prunes
at the generation the children that cannot be better than $z$. Line 12
is the insertion of A$^*$ (a lookup in the CLOSED list of $\phi(v)$, a
new node or a reopened one if the path is shorter, and a push on
OPEN($\phi(v)$)), made with no lock, since $\phi(v) \in S(b)$ by Lemma
11.3. Line 13 publishes $f(b)$, an atomic value that the other workers
read without the lock at lines 18 and 28. The insertions into an nblock
$a \neq b$ lower its published value as well.

**Algorithm 11.2** (Safe PBNF: the switch).

```text
14  procedure SHOULDSWITCH(b, k)
15     if OPEN(b) = ∅ then return true
16     if k < m then return false
17     k ← 0;  lock the graph
18     b' ← the nblock of I(b) of least published f
19     if f(free list) < f(b) or f(b') < f(b) then
20        if f(b') < f(free list) then SETHOT(b')
21        return true
22     SETCOLD each hot nblock of I(b);  return false
23  procedure NEXTNBLOCK(b)
24     if b ≠ nil, OPEN(b) ≠ ∅ and b was not just set hot then
25        if the graph is locked by another worker then return b
26     lock the graph
27     if b ≠ nil then
28        if f(b) < z, f(b) ≤ f(b'), f(b) ≤ f(free list) and b was not
             just set hot then return b
29        RELEASE(b)
30     drop from the free list the nblocks with no node of f < z
31     while the free list is empty do
32        if no nblock is in use then done ← true;  wake all;  return nil
33        sleep until woken
34     b ← the first nblock of the free list;  ACQUIRE(b);  return b
```

The check (lines 17 to 22) is made after $m$ expansions at least, and
always when OPEN($b$) is empty (line 15). Line 19 is the switching rule
of PBNF, and line 20 is the addition of Safe PBNF to it. The best
interfering nblock $b'$ is set hot if no free nblock is as "good", since
otherwise the worker can switch to a free nblock that is at least as
good, since a free nblock can be acquired at once, while a hot one waits
for the releases of the others. "Just set hot" in lines 24 and 28 means
that the check before the call has called SETHOT (the code sets the flag
even when the call returns at once because $b'$ was already hot). Then
the worker always takes the lock and releases $b$, and this release
eventually frees the hot nblock. Line 25 is the try-lock of PBNF, and line 28
tests again, under the lock, the condition of line 19 on the values of
that moment: a worker keeps $b$ if, in the meantime, it has become the
best again. Line 30 drops from the free list the nblocks whose key is
not below the incumbent (they were free when $z$ was larger), which is a
scan of the entries with a key not below $z$ only. Lines 31 to 33 are
the end and the wait: a worker that finds no free nblock sleeps on a
condition variable until FREE wakes it, and the last worker to find no
nblock in use declares the end.

**Algorithm 11.3** (Safe PBNF: the operations on the nblock graph, all
under its lock).

```text
35  procedure ACQUIRE(b)
36     b leaves the free list and is in use
37     for each n ∈ I(b) do σ(n) ← σ(n) + 1;  n leaves the free list
38  procedure RELEASE(b)
39     b is not in use
40     for each n ∈ I(b) do
41        σ(n) ← σ(n) − 1
42        if σ(n) = 0 then
43           if σ_h(n) = 0 and n is hot then SETCOLD(n)
44           FREE(n)
45     FREE(b)
46  procedure SETHOT(b)
47     if b is hot, σ(b) = 0, or a hot nblock of I(b) has a smaller f
          then return
48     b is hot
49     for each n ∈ I(b) do
50        if n is hot then SETCOLD(n)
51        n leaves the free list;  σ_h(n) ← σ_h(n) + 1
52  procedure SETCOLD(b)
53     b is not hot
54     for each n ∈ I(b) do
55        σ_h(n) ← σ_h(n) − 1
56        if σ(n) = 0 and σ_h(n) = 0 then
57           if n is hot then SETCOLD(n)
58           FREE(n)
59  procedure FREE(n)
60     if n is not in use, σ(n) = σ_h(n) = 0 and f(n) < z then
61        put n in the free list, with the key f(n);  wake all
```

Note that FREE is the only way into the free list, and it tests exactly
the condition of Definition 11.5, computing $f(n)$ from OPEN($n$) (no
worker can touch OPEN($n$) at that moment, Theorem 11.8). SETHOT refuses
an nblock with $\sigma(b) = 0$ (which no nblock in use keeps from being
free) and one that interferes with a better hot nblock (case 2 of Burns
et al.). Otherwise, it cools the worse hot nblocks of $I(b)$ (case 3)
and blocks the whole of $I(b)$. Line 43 cools a hot nblock when the last
nblock in use of its interference scope is released, which is the event
that the hot flag waits for, and line 57 is a guard that Lemma 11.7
shows never to fire.

Our pseudocode differs from appendix A of [[3]](Z-bibliography.md#ref-3)
in six details, none of which changes the arguments of section 11.6.
First, appendix A prunes a single node with $f$ not better than the
incumbent, while line 8 empties the whole list, which is equivalent
since OPEN is ordered by $f$. Second, an nblock enters the free list
only if $f(n) < z$, and line 30 removes those that the incumbent has
made useless, while appendix A admits any nblock that is not empty.
Third, line 28 keeps $b$ also on ties ($\le$ instead of $<$). Fourth,
ACQUIRE takes the nblocks of $I(b)$ out of the free list explicitly.
Fifth, $I(b)$ does not contain $b$ (in appendix A, whose interference
scope is built with the self-loops, it does): being in use is a flag of
its own, and RELEASE puts $b$ itself back with FREE at line 45. Finally,
the end is detected by the number of nblocks in use instead of the
condition $\sigma \equiv 0$ of appendix A; the former implies the
latter by Lemma 11.7 (a), and it is the right one in the library, since
an nblock in use with an empty interference scope (e.g., the only
nblock of a trivial abstraction) leaves $\sigma \equiv 0$.

## 11.6 Analysis {#sec-11-6}

We now prove the properties of the search. Our arguments are about the
states of the nblock graph between two critical sections, i.e., at the
instants at which no worker holds its lock. The operations under the
lock (Algorithms 11.2 and 11.3, and the insertion of the root) are
atomic with respect to each other, while the expansions (lines 7 to 13)
run in parallel with the rest.

<a id="lem-11-7"></a>**Lemma 11.7** (Invariants of the nblock graph). At
each instant at which the lock of the graph is free, the following hold.

* (a) $\sigma(n)$ is the number of the nblocks in use in $I(n)$.
* (b) $\sigma_h(n)$ is the number of the hot nblocks in $I(n)$.
* (c) Each nblock $n$ of the free list is not in use and has
  $\sigma(n) = \sigma_h(n) = 0$, and its key is the least $f$ of
  OPEN($n$).
* (d) No two nblocks in use interfere, i.e., their scopes are disjoint.
* (e) Each hot nblock $b$ has $\sigma(b) \ge 1$ and $\sigma_h(b) = 0$,
  i.e., it interferes with an nblock in use and with no hot nblock.
* (f) Each nblock $n$ that is not in use, has
  $\sigma(n) = \sigma_h(n) = 0$ and has a node of $f < z$ in OPEN($n$)
  is in the free list.

*Proof.* We proceed by induction on the critical sections. At the start
no nblock is in use or hot, all the counters are 0, and the only nblock
with a node is that of $s_0$, which is put in the free list; hence (a)
to (f) hold.

(a) and (b). ACQUIRE($b$) and RELEASE($b$) change $\sigma$ by one on
$I(b)$ exactly when $b$ starts and stops being in use, and SETHOT($b$)
and SETCOLD($b$) change $\sigma_h$ by one on $I(b)$ exactly when $b$
becomes hot and cold (SETHOT returns at line 47 when $b$ is already hot,
and SETCOLD is called only on a hot nblock, at lines 22, 43, 50 and 57).
Hence, we have that $\sigma(n)$ counts the nblocks $c$ in use with
$n \in I(c)$, which by Lemma 11.4 are the nblocks in use of $I(n)$. The
same holds for $\sigma_h$.

(c). An nblock enters the free list only through FREE, which checks the
condition and takes the key from OPEN($n$). It leaves the list at its
own ACQUIRE, at the ACQUIRE of an nblock of $I(n)$ (line 37, when
$\sigma(n)$ becomes positive), at the SETHOT of an nblock of $I(n)$
(line 51, when $\sigma_h(n)$ becomes positive), and at line 30. While it
is in the list, $\sigma(n) = 0$ and $n$ is not in use, and therefore, by
(a), no nblock in use has $n$ in its scope (if $n \in S(c)$ with
$c \neq n$, then $n \in \mathrm{succ}(c) \subseteq I(c)$ by Lemma 11.4);
hence, no worker writes OPEN($n$), and the key stays exact.

(d). ACQUIRE takes the first nblock $b$ of the free list, for which
$\sigma(b) = 0$ by (c), i.e., no nblock of $I(b)$ is in use, by (a).

(e). SETHOT($b$) makes $b$ hot only if $\sigma(b) \ge 1$ (line 47), and
it cools the hot nblocks of $I(b)$ (line 50), after which
$\sigma_h(b) = 0$ by (b). A hot $b$ keeps $\sigma_h(b) = 0$, since an
nblock $x \in I(b)$ becomes hot only through SETHOT($x$), which cools
$b$ first at line 50, unless it returns at line 47. Also, $\sigma(b)$
decreases only at line 41 of a RELEASE, and when it reaches 0 the
nblock, hot with $\sigma_h(b) = 0$, is cooled at line 43. By (e), the
condition of line 57 (a hot nblock with $\sigma = 0$) never holds at the
start of a critical section, and SETHOT and SETCOLD do not create it;
hence, line 57 is never executed.

(f). The nodes of OPEN($n$) are written only by the worker that holds
$n$ or an nblock $c$ with $n \in S(c)$ (lines 7 to 13), i.e., while $n$
is in use or $\sigma(n) \ge 1$, by (a) and Lemma 11.4. Each event that
brings $n$ to "not in use and $\sigma(n) = \sigma_h(n) = 0$" calls
FREE($n$): the release of $n$ (line 45), the release of the last nblock
in use of $I(n)$ (line 44), and the cooling of the last hot nblock of
$I(n)$ (line 58). FREE puts $n$ in the list if $f(n) < z$. If not, $n$
has no node of $f < z$, and it gets none until it is in a scope again,
since $z$ only decreases. Finally, $n$ leaves the list only when one of
the three conditions fails, or at line 30 when its key, which is exact
by (c), is not below $z$. $\square$

Indeed, Lemma 11.7 (e) is Lemma 1 of Burns et al. ("If an nblock n is
hot, there is at least one other nblock in its interference scope that
is in use. Also, n is not interfering with any other hot nblocks",
[[3]](Z-bibliography.md#ref-3), p. 700), proved here on the pseudocode
of the library. As a first consequence, we can show that the lists of an
nblock have one writer at a time, which is the purpose of the whole
construction.

<a id="thm-11-8"></a>**Theorem 11.8** (Exclusive access). At each
instant, the OPEN and CLOSED lists of an nblock $n$ are read or written
by at most one worker outside the lock of the graph, namely the worker
that holds the nblock $c$ with $n \in S(c)$, if there is one. Under the
lock, they are read only by that worker or when no such worker exists.

*Proof.* We first observe that a worker that holds $c$ reads and writes,
in lines 7 to 13, the lists of $c$ (lines 7 and 8) and those of the
nblocks $\phi(v)$ of the children of the nodes of $c$, which lie in
$S(c)$ by Lemma 11.3. If two workers held $c \neq c'$ with
$n \in S(c) \cap S(c')$, then $c$ and $c'$ would interfere, against
Lemma 11.7 (d); hence, at most one worker touches the lists of $n$.
Under the lock, OPEN($n$) is read by FREE (line 60), which requires $n$
not in use and $\sigma(n) = 0$, and therefore no nblock in use with $n$
in its scope, by Lemma 11.7 (a), and at line 28 by the worker that holds
$n$. $\square$

<a id="cor-11-9"></a>**Corollary 11.9** (Duplicate detection). Each
state has at most one node in the search, and the test of line 12 that
compares a child with the CLOSED list of its nblock is the test of the
sequential A$^*$: a child is a duplicate if and only if a node of its
state exists, and it is reopened if and only if its new path is shorter.

*Proof.* We store the node of a state $v$ only in the CLOSED list (the
table of the nodes) of $\phi(v)$, which by Theorem 11.8 has one writer
at a time; the lookup and the insertion of line 12 are therefore those
of a sequential table. $\square$

One can therefore say that the duplicate detection is "perfect" in the
sense of SDD, and no node is sent from a worker to another: the only
shared writes outside the lock of the graph are the published values
$f(b)$ (atomic), the incumbent (under a separate lock) and the creation
of an nblock the first time that it is touched (under a third lock,
section 11.7). Since the search has more than one lock, the absence of
deadlock needs a (short) argument.

<a id="prop-11-10"></a>**Proposition 11.10** (No deadlock). No worker
waits forever for a lock, and a worker that sleeps at line 33 is woken
when an nblock enters the free list or the search ends.

*Proof.* We recall that the three locks of the search are those of the
graph, of the incumbent and of the creation of the nblocks. The second
is taken only in the goal test of line 10, which holds no other lock and
takes no other; the third is taken in the creation of an nblock,
possibly under the lock of the graph, and nothing is locked under it.
Hence, we always take the locks in the order graph, creation, and a
cycle of waits is impossible. The condition variable of line 33 is
waited on under the lock of the graph and notified under it by FREE and
at line 32, so that no notification is lost. Also, a worker that sleeps
has released its nblock at line 29, or had none. $\square$

This is section 3.2.2 of Burns et al. for the one lock of their
pseudocode, which "the thread that currently holds it never attempts to
acquire a second time" (p. 699). We now come to the optimality, for
which one has to read the argument of A$^*$ on the union of all the OPEN
lists.

<a id="lem-11-11"></a>**Lemma 11.11** (An optimal path is always
represented). Let $h$ be admissible, let
$P = (s_0 = w_0, w_1, \ldots, w_k)$ be an optimal path to a goal, and
let $z > C^*$ at some instant. Then some $w_j$ has a node with
$g(w_j) = g^*(w_j)$ that is either in OPEN($\phi(w_j)$) or being
expanded, and therefore $f(w_j) \le C^* < z$.

*Proof.* We let $j$ be the largest index such that $w_j$ has a node with
$g(w_j) = g^*(w_j)$; it exists, since $s_0$ has $g = 0$ from the start,
and $g$ of a node never increases. Suppose that the node of $w_j$ is
neither in OPEN nor being expanded: then its expansion with
$g = g^*(w_j)$ has completed. If $j = k$, the goal $w_k$ was selected
with $g = C^*$ at line 10, and $z \le C^*$. If $j < k$, the expansion
generated $w_{j+1}$ with $g = g^*(w_{j+1})$, which line 11 did not
prune, since $g^*(w_{j+1}) + h(w_{j+1}) \le C^* < z$ by admissibility;
then the node of $w_{j+1}$ received $g = g^*(w_{j+1})$ at line 12, or
had it already (Corollary 11.9), against the choice of $j$. Hence, we
conclude that the node of $w_j$ is in OPEN or being expanded, and
$f(w_j) = g^*(w_j) + h(w_j) \le C^*$. Line 8 does not remove it, since
it empties OPEN only when its first node, and hence each node, has
$f \ge z$. $\square$

<a id="thm-11-12"></a>**Theorem 11.12** (Optimality). Let $h$ be
admissible. If `SafePBNF` ends with the status `solved`, the path that
it returns is optimal. If it ends with the status `unsolvable`, no goal
is reachable.

*Proof.* We recall that the search ends at line 32, when no nblock is in
use and the free list is empty (line 31, after line 30). Then
$\sigma \equiv 0$ by Lemma 11.7 (a), no nblock is hot by (e), and
therefore $\sigma_h \equiv 0$ by (b); by (f), each nblock with a node of
$f < z$ in OPEN is in the free list, which is empty. No node is being
expanded, since an expansion happens only in an nblock in use. Hence, if
a goal were reachable and $z > C^*$, the node of Lemma 11.11 would be in
some OPEN with $f < z$, a contradiction. Therefore, $z \le C^*$ if a
goal is reachable, and $z = +\infty$ (the status `unsolvable`) only if
none is. The incumbent $z$ is the value $g$ of a goal node, and the path
is rebuilt from that node through the parents after all the workers have
returned. Along this chain
$g(v) \ge g(\mathrm{parent}(v)) + c(\mathrm{parent}(v), v)$, since the
parent's $g$ can only have decreased after it was used; hence $g$
decreases strictly along the chain (the costs are positive), the chain
is acyclic, and its cost is at most $z \le C^*$, i.e., it is optimal.
$\square$

This is Theorem 3 of Burns et al. (p. 701), in whose short proof one
reads that "the search continues to expand nodes until the minimum f
value among all frontier nodes is greater than or equal to the
incumbent". The part that needs the invariants is that the end condition
of line 32 implies it, i.e., implies the condition of the end of a
parallel search of [Theorem 8.6](08-shared-memory.md#thm-8-6): no entry
in flight that is not stale has $f < z$, i.e., no such entry in an OPEN
list, no worker expanding a node, and no node in transit, the last being
trivial here, since a node is never in transit.

<a id="thm-11-13"></a>**Theorem 11.13** (Termination on a finite state
space). Let the state space reachable from $s_0$ be finite, with
positive costs. Then `SafePBNF` ends, with any number of workers and
any $m \ge 1$, under any scheduling of the threads in which a worker
that can proceed eventually does (the assumption (A1) below).

*Proof.* A node is expanded once for each value that its $g$ takes, and
$g$ takes strictly decreasing values, each the cost of the path of its
parents at the time, which is simple since $g$ grows strictly along it;
as the graph is finite, the simple paths are finitely many, and so are
the expansions and the insertions into OPEN. Between two acquisitions, a
worker expands at least one node or empties a list: after ACQUIRE,
$k = 0 < m$ and OPEN($b$) holds a node of $f < z$ by Lemma 11.7 (c)
(line 30 has just removed the keys not below $z$), and line 6 lets it
expand (or, if a goal has lowered $z$ meanwhile, which happens under the
lock of the incumbent and not under that of the graph, line 8 empties
the list). Also, a worker that keeps its nblock at line 25 or 28
returns to line 6 with $k = 0$. Each emptying of a list at line 8
removes at least one entry, and the entries are finitely many. Hence,
we have finitely many critical sections as well, and each of them ends.
A worker that sleeps at line 33 holds no nblock, and if no other worker
holds one, the last one to look at the free list finds it empty with no
nblock in use and ends the search (line 32), and by Proposition 11.10 no
worker is left asleep while an nblock is free. $\square$

Burns et al. prove completeness for an admissible heuristic or a finite
space (Corollary 1, p. 701), relying for the former on their Theorem 2.
On an infinite space the argument above fails at its first step (and a
node limit is then the only "safeguard" that one may rely on), and
termination needs that the node of Lemma 11.11 be expanded in finite
time, which is the progress property that the hot nblocks are meant to
give. It also needs, as usual, costs bounded below by some
$\varepsilon > 0$, as the sequential A$^*$ does ([Chapter
2](02-best-first.md#ch-2)). In the following theorem we state what the
pseudocode of the library guarantees for a hot nblock. It needs three
assumptions on the execution, which do not follow from the pseudocode
(although they are, arguably, rather mild):

* (A1) each worker is scheduled, and each expansion takes a finite time,
  so that a worker that holds an nblock reaches its next check (line 14
  with $k \ge m$, or line 15 with an empty list) in finite time;
* (A2) the lock is fair, and the try-lock of line 25 does not fail at
  infinitely many consecutive attempts of a worker;
* (A3) at each check, the published value of each nblock is the least
  $f$ of its OPEN list, and the test of line 28 gives the same answer as
  the test of line 19 of the same check.

The first two concern the threads and the lock, and they usually hold
for the threads of C++ and of OpenMP (although neither standard promises
a fair lock). Instead, (A3) closes the window between the check, whose
lock is released at its end, and NEXTNBLOCK, which takes the lock again
(the same window is in appendix A, whose SHOULDSWITCH reads the values
without the lock).

<a id="thm-11-14"></a>**Theorem 11.14** (A hot nblock is freed or
overtaken). Assume (A1) to (A3), and let the nblock $b$ be hot at some
instant. Then, after finitely many checks, one of the following happens.

* (i) $\sigma(b)$ reaches 0, and $b$ is set cold and freed (lines 42 to
  44), entering the free list if it has a node of $f < z$.
* (ii) $b$ is set cold at line 22 by a worker whose nblock $c \in I(b)$
  has $f(c) \le f(b)$.
* (iii) $b$ is set cold at line 50 by SETHOT($x$) for some $x \in I(b)$
  with $f(x) \le f(b)$.

*Proof.* We first show that, while $b$ is hot, $\sigma(b)$ does not
increase: each $n \in I(b)$ has $\sigma_h(n) \ge 1$ by Lemma 11.7 (b),
and hence it is not in the free list, by (c); and an nblock is acquired
only from the free list; hence no nblock of $I(b)$ starts being in use.
Let $c \in I(b)$ be in use, held by the worker $i$. By (A1), $i$ reaches
a check in finite time. If OPEN($c$) is empty, $i$ calls NEXTNBLOCK,
which takes the lock (lines 24 and 26) and releases $c$, since
$f(c) = +\infty$ fails line 28. Otherwise, let $b'$ be the best
nblock of $I(c)$; since
$b \in I(c)$ by Lemma 11.4, $f(b') \le f(b)$. If the test of line 19
fails, then $f(c) \le f(b') \le f(b)$, and line 22 cools $b$, which is
case (ii). If it succeeds and $f(b') < f(\text{free list})$, line 20
calls SETHOT($b'$), which either returns at once ($b' = b$, or $b'$
already hot, or a better hot nblock interferes with $b'$) or makes $b'$
hot, cooling $b$ at line 50 if $b \in I(b')$ (case (iii), since
$f(b') \le f(b)$). In all the other sub-cases the flag "just set hot"
makes NEXTNBLOCK take the lock and release $c$. If it succeeds and
$f(\text{free list}) \le f(b')$, NEXTNBLOCK obtains the lock by (A2),
and line 28 fails by (A3); hence $c$ is released. Each release of an
nblock of $I(b)$ decreases $\sigma(b)$ by one, and
$\sigma(b) \le p$, since at most $p$ nblocks are in use; hence,
unless (ii) or (iii) happens first, after at most $\sigma(b)$ checks of
these workers $\sigma(b) = 0$, at which point line 43 cools $b$ (which
has $\sigma_h(b) = 0$ by Lemma 11.7 (e)) and line 44 frees it, which is
case (i). Lines 43 and 57 are the only other places where $b$ could be
cooled, and the second is never executed, by Lemma 11.7. $\square$

Note that in cases (ii) and (iii) the hot flag passes to, or is
overridden by, an nblock that is not worse than $b$ and is in use or
hot. When $b$ holds the unique best nodes of the search, neither can
happen, and the theorem gives the progress for which Safe PBNF is
designed.

<a id="cor-11-15"></a>**Corollary 11.15** (Progress with a unique best
nblock). Assume (A1) to (A3) and a consistent $h$. Let $\varphi < z$ be
the least $f$ of the nodes in the OPEN lists at some instant, and let
all the nodes with $f = \varphi$ lie in one nblock $B$, not in use. Then
the following hold.

* (a) If $B$ is free, a worker acquires $B$ at the next acquisition, and
  one happens after finitely many checks.
* (b) If $B$ is hot, a worker acquires $B$ after finitely many checks.
* (c) If $\sigma(B) \ge 1$, then after finitely many checks $B$ is hot,
  or it is free, or $\sigma(B) = 0 < \sigma_h(B)$.

*Proof.* We first note that, since $h$ is consistent, a child has an $f$
not smaller than that of its parent; hence, until $B$ is acquired, no
node of $f \le \varphi$ is generated outside $B$, and the other nblocks
have $f > \varphi$. (a) While $B$ is free it is the first of the free
list, and it can leave the list only when it is acquired: an nblock of
$I(B)$ cannot be acquired before it (it is behind $B$ in the list), and
SETHOT($x$) with $B \in I(x)$ needs
$f(x) < f(\text{free list}) = \varphi$ at line 20, which no $x \neq B$
has. A worker that sleeps was woken when $B$ entered the list, and a
worker that holds an nblock $c$ finds at its next check, by (A1),
$f(\text{free list}) = \varphi < f(c)$, does not set anything hot
($f(b') > \varphi$, since $\sigma(B) = 0$ means that no nblock in use
interferes with $B$), and releases $c$ by (A2) and (A3), reaching
line 34. (b) By Theorem 11.14, since cases (ii) and (iii) need an nblock
other than $B$ with $f \le \varphi$, $B$ is freed in case (i); the
worker that releases the last nblock in use of $I(B)$ goes on, under the
same lock, to line 34, where $B$ is the first of the free list. (c) Let
$c \in I(B)$ be in use. At its check, the best of $I(c)$ is $B$ (by
Lemma 11.4, $B \in I(c)$), with $\varphi < f(c)$, and
$f(\text{free list}) > \varphi$ since $B$ is not free; hence line 20
calls SETHOT($B$), which makes $B$ hot unless $\sigma(B)$ has dropped to
0 in the meantime (no hot nblock of $I(B)$ has a smaller $f$). In the
latter case $B$ is free if $\sigma_h(B) = 0$, and the third alternative
holds otherwise. $\square$

<a id="rem-11-16"></a>**Remark 11.16** (What is taken from Burns et
al.). The case that we leave open in Corollary 11.15 is that of an
nblock $B$ with the best nodes, no nblock in use around it, and a worse
hot nblock $x$ in $I(B)$: by Theorem 11.14, $x$ is freed or overtaken in
finite time, which decreases $\sigma_h(B)$, but new hot nblocks may
appear in $I(B)$ meanwhile, and we have no bound on such a sequence.
Burns et al. state the progress property as their Theorem 2 ("A node n
with minimum f value will eventually be expanded", pp. 700-701), whose
third case (an nblock kept from being free by an nblock in use) is
reduced to their Theorem 1 ("If an nblock n becomes hot, it will
eventually be added to the free list and will no longer be hot", p.
700). The latter is proved, in TLA$^+$, on a model in which a search
step only chooses an nblock to set hot and then releases the nblock that
it holds (section 3.2.3, p. 699, and appendix B, pp. 740-743). In the
pseudocode of appendix A, as in ours (line 22), a worker that keeps its
nblock may also cool a hot nblock, which is why Theorem 11.14 has the
cases (ii) and (iii). Burns et al. also checked with the model checker
TLC that their model of Safe PBNF has no livelock with up to three
threads and 12 nblocks on an undirected ring, and up to three threads
and eight nblocks on a directed graph (p. 699). We take the livelock
freedom of Safe PBNF beyond Corollary 11.15 from these results, within
the limits of their model, and we do not claim more. In particular, no
test of the library can check it, since a livelock needs an infinite
state space.

It remains to state how the workers learn that the search is over. The
flag `done` of the nblock graph is the only signal: it is set at line 32,
or by a worker that reaches the node limit (under the lock of the graph,
`SafePBNF.h:529-536`), and a worker returns only when it has read it, or
at line 32. The following proposition shows that the detection is
both sound (the flag is never set too early) and live (once it is set,
every worker returns).

<a id="prop-11-17"></a>**Proposition 11.17** (The detection of the
end). (i) When `done` is set at line 32, no worker holds an nblock or
expands a node, and no OPEN list holds an entry that is not stale with
$f < z$; hence the incumbent is final, and Theorem 11.12 applies.
(ii) Under (A1) and (A2) of section 11.6, once `done` is set (at line
32 or by the node limit),
every worker returns in finite time. (iii) A worker that sleeps at
line 33 when the search is not over is woken by the next FREE that puts
an nblock in the free list, or by the setting of `done`.

*Proof.* (i) Line 32 runs under the lock of the graph, with the free
list empty after line 30 and no nblock in use. No worker expands a node
outside an nblock in use, and a worker inserts nodes only during an
expansion; hence, no OPEN list changes until a worker acquires an
nblock, which needs a free one. By Lemma 11.7 (f), as in the proof of
Theorem 11.12, no nblock has a node of $f < z$, and the entries that are
not stale are the nodes of the lists. (ii) A worker is in one of three
places: asleep at line 33, in a check or a call of NEXTNBLOCK, or in
the expansions of lines 7 to 13. The flag is set under the lock and
followed by a notification of all the sleepers under the same lock
(`wake all`); a sleeper wakes, takes the lock and returns at its next
test of `done` (`SafePBNF.h:591-592`), and a worker that waits for the
lock does the same once it has it. A worker in its expansions reaches a
check in finite time by (A1); there, either OPEN($b$) is empty and
NEXTNBLOCK takes the lock, or SHOULDSWITCH takes it after $m$
expansions and returns true on reading `done` (`SafePBNF.h:553-554`).
Then NEXTNBLOCK takes the lock (a try-lock that fails only returns to
the expansions, which reach the next check, and by (A2) it does not fail
forever), does not keep $b$ (its
test reads `done`), releases it and returns at the test of `done`.
(iii) The waits and the notifications are under the same lock, and the
notification of FREE (line 61) follows the insertion into the free
list; a worker tests the free list and the flag under the lock before
each wait, and a spurious wakeup only repeats the test, so that no
notification is lost. $\square$

Note that the node limit leaves the worker that reaches it holding its
nblock (it returns at once), so that the counter of the nblocks in use
does not reach 0 afterwards; this is harmless, since every other worker
returns on reading `done`, and the status is then `node_limit` whatever
the incumbent (`SafePBNF.h:270-271`), as Exercise 11.5 asks to show.

## 11.7 The implementation {#sec-11-7}

`SafePBNF::solve()` (`SafePBNF.h:239-280`) builds the nblock graph and
the workers, inserts the root, runs the workers on the executor, and
then collects their counters and rebuilds the path. The timer of the
`Solution` starts at the beginning of `solve()` and stops at its end. We
now go through the parts of the code that are specific to PBNF, i.e.,
the creation of the nblocks, the free list and the expansion, with the
cost and the memory of each; the OPEN lists, the tables of the nodes and
the node arena are those of the sequential A$^*$ ([Chapter
2](02-best-first.md#ch-2) and [Chapter 3](03-data-structures.md#ch-3)).
The pseudocode below continues the numbering of the lines of
Algorithms 11.1 to 11.3.

### The nblocks

An nblock (`SafePBNF.h:291-303`), which typically groups many states
(400 cells on the grids of the benchmark), holds its OPEN list (the
bucket queue for integer costs and the binary heap otherwise), its table
of nodes (the CLOSED list, which starts with 16 slots and grows by
doubling), its scope and its interference scope as vectors of ids, the
counters $\sigma$ and $\sigma_h$, the flags hot, in use and in the free
list, its key in the free list, and the published value $f(b)$, an
atomic read and written with relaxed order (i.e., with no ordering
w.r.t. the other memory operations). The nblocks are created the first
time that a worker touches them, so that only the visited part of the
abstract graph costs memory (Burns et al. do the same, section 4, p.
702). This matters somewhat for the 15-puzzle, whose 4096 ids include
736 that are not abstract states, and much more for "large" grids, of
which the search often reaches only a part. The search keeps an array
$A$ of $N$ pointers, all nil at the start, and creates an nblock by a
double-checked lookup under a separate lock, as follows.

**Algorithm 11.4** (Safe PBNF: the creation of an nblock).

```text
62  function BLOCK(b)                   ▷ the nblock b, built at its first use
63     x ← A[b], read with acquire order;  if x ≠ nil then return x
64     lock the creation
65     x ← A[b];  if x ≠ nil then unlock the creation;  return x
66     T ← succ(b);  J ← T
67     for each a ∈ T do J ← J ∪ succ(a)      ▷ predecessors = successors
68     I(b) ← J ∖ {b}, sorted, each nblock once;  S(b) ← T ∪ {b}
69     OPEN(b) ← ∅;  CLOSED(b) ← a table of 16 slots;  σ(b) ← σ_h(b) ← 0
70     A[b] ← the new nblock, written with release order;  unlock the creation
71     return A[b]
```

Line 67 is Lemma 11.4, which holds only on an undirected abstract graph
(condition (ii) of Definition 11.1); sorting $I(b)$ once makes the loops
of Algorithm 11.3 visit the nblocks in a fixed order, which (arguably)
makes an execution somewhat easier to follow in a debugger.

<a id="lem-11-18"></a>**Lemma 11.18** (The creation of the nblocks).
For each $b$, Algorithm 11.4 creates at most one nblock, every call
returns that nblock, and a worker that obtains it at line 63 sees its
scopes, counters and empty lists as line 69 wrote them.

*Proof.* We first note that the creation happens only between lines 64
and 70, under the lock of the creation, and that line 65 tests $A[b]$
again under that lock; hence, a second creation of $b$ would need a
first one that has not reached line 70, which the lock excludes. A call
that does not create returns the pointer that it read, which, once not
nil, never changes (the nblocks are destroyed only when `solve()`
returns). Finally, line 70 is a store with release order that follows
all the writes of lines 66 to 69, and line 63 is a load with acquire
order: by the memory model of C++, a load with acquire order that reads
the value of a store with release order sees all the writes that
precede the store in its thread. A reader at line 65 sees them through
the lock of the creation. $\square$

After the creation, the scopes of an nblock are only read, and its lists
are written only by the worker of Theorem 11.8; the creation lock is
thus taken at most twice for each nblock that a worker reaches first.
The size of the scopes, on which the cost of every operation of the
graph depends, is bounded as follows.

<a id="lem-11-19"></a>**Lemma 11.19** (The size of the scopes). Let $d$
be the largest number of abstract successors of an nblock. Then
$|S(b)| \le d + 1$ and $|I(b)| \le d^2$ for each $b$. On a grid of
nblocks with 4 moves $|I(b)| \le 12$, with 8 moves $|I(b)| \le 24$, both
with equality for an nblock at distance at least 2 from the border, and
in the vacuum world $|I(b)| \le 4$.

*Proof.* By definition $|S(b)| = 1 + |\mathrm{succ}(b)| \le 1 + d$. By
Lemma 11.4, $I(b)$ is made of the at most $d$ successors of $b$ and, for
each successor $a$, of the successors of $a$ other than $b$, which are
at most $d - 1$ since $b \in \mathrm{succ}(a)$ by condition (ii); hence
$|I(b)| \le d + d(d - 1) = d^2$. On a grid of nblocks with 4 moves the
successors of $b$ are the nblocks at Manhattan distance 1, and their
successors those at distance 0 or 2 from $b$; hence $I(b)$ is the set of
the nblocks at distance 1 or 2, i.e., at most 4 + 8 = 12, all of them
when none falls outside the grid. With 8 moves the successors are the
nblocks at Chebyshev distance 1, and the nblocks reached by two of them
are those at Chebyshev distance at most 2; hence $I(b)$ is the
$5 \times 5$ square around $b$ without $b$, i.e., at most 24 nblocks.
In the vacuum world the abstract graph is a path, $d = 2$, and $I(b)$ is
made of the nblocks at distance 1 or 2 along it. $\square$

For the sliding-tile puzzles with the abstraction of section 11.1 one
has $d \le 4$ (one abstract successor for each move of the blank), and
hence $|I(b)| \le 16$, a bound that we have not checked to be tight.

### The free list and the graph operations

ACQUIRE and RELEASE are lines 35 to 45 of the pseudocode, and SETHOT and
SETCOLD are lines 46 to 58, with the same names in lowercase in the code
(`SafePBNF.h:406-470`); SETHOT reads the published values of the hot
nblocks of $I(b)$ to apply case 2 of Burns et al. We keep the free list
$F$ in an ordered map (a B-tree of the default data structures) keyed by
the pairs $(f(b), b)$, so that the first entry is the best nblock, the
ties go to the smaller id, the removal of a given nblock is exact, and
line 30 is a single scan from the first key not below $z$. An nblock
knows whether it is in $F$ and with which key, and the operations of the
pseudocode that touch $F$ are the following (`make_free`, `make_busy`,
`prune_free`, `best_free` and `best_scope`, `SafePBNF.h:366-404` and
`472-479`).

**Algorithm 11.5** (Safe PBNF: the free list, under the lock of the
graph).

```text
72  procedure FREE(n)                          ▷ line 59 of Algorithm 11.3
73     if n ∈ F, n is in use, σ(n) ≠ 0 or σ_h(n) ≠ 0 then return
74     drop the stale entries at the top of OPEN(n);  φ ← their least f, or ∞
75     if φ ≥ z then return
76     key(n) ← φ;  publish f(n) ← φ;  F ← F ∪ {(φ, n)};  wake all
77  procedure BUSY(n)                          ▷ "n leaves the free list"
78     if n ∈ F then F ← F ∖ {(key(n), n)}
79  procedure PRUNEFREE()                      ▷ line 30
80     for each (φ, n) ∈ F with φ ≥ z, in increasing order do F ← F ∖ {(φ, n)}
81  function BESTFREE()                        ▷ f(free list)
82     return the key of the first pair of F, or ∞ if F = ∅
83  function BESTSCOPE(b)                      ▷ line 18
84     return the least published f(a) over the built nblocks a ∈ I(b),
          with its nblock, or (∞, nil)
```

Line 74 computes $f(n)$ from OPEN($n$), not from the published value,
as the proof of Lemma 11.7 (c) requires; line 76 publishes it, so
that the value read by the others at line 18 is exact for a free
nblock. BESTSCOPE skips the nblocks not yet built, whose OPEN lists are
empty. We can now bound the cost of the critical sections.

<a id="prop-11-20"></a>**Proposition 11.20** (The cost of the graph
operations). Let $F$ also denote the number of the entries of the free
list ($F \le N$), let $\Delta$ be the largest $|I(b)|$ of an nblock
built, and let each insertion, removal and lower bound of the B-tree
cost $O(\log F)$. Then, under the lock of the graph, (i) BUSY and
BESTFREE take $O(\log F)$ time, and FREE takes $O(\log F)$ time plus one
extraction for each stale entry that it drops; (ii) ACQUIRE and SETCOLD
take $O(\Delta \log F)$ time, RELEASE $O(\Delta^2 \log F)$, and SETHOT
$O(\Delta^2 \log F)$, the stale entries dropped by FREE apart;
(iii) PRUNEFREE takes $O((r + 1) \log F)$ time, $r$ being the number of
the pairs that it removes; (iv) the check of lines 17 to 22 takes
$O(\Delta)$ time besides its calls of SETHOT and SETCOLD. Over the whole
search, each entry of an OPEN list is dropped at most once, and each
pair removed by PRUNEFREE was inserted by an earlier FREE.

*Proof.* (i) is the cost of one operation of the B-tree, plus the
constant work of the tests. (ii) ACQUIRE makes one BUSY and constant
work for each $n \in I(b)$. SETCOLD makes, for each $n \in I(b)$, one
decrement and possibly one FREE, and no recursive call, since line 57 is
never executed (Lemma 11.7). RELEASE makes, for each $n \in I(b)$, one
decrement, possibly one SETCOLD (line 43) and one FREE, and one more
FREE for $b$. SETHOT reads the published values of the hot nblocks of
$I(b)$ in $O(\Delta)$, and then makes, for each $n \in I(b)$, possibly
one SETCOLD and one BUSY. (iii) Each step of the scan removes one pair,
and the scan starts with one lower bound. (iv) BESTSCOPE reads one
atomic value for each nblock of $I(b)$, BESTFREE is (i), and line 22
tests the flag of each nblock of $I(b)$. Finally, an entry dropped by
FREE or at line 7 leaves its list, and a pair enters $F$ only at line
76. $\square$

By Lemma 11.19, $\Delta \le 24$ on the grids and $\Delta \le 4$ in the
vacuum world; hence the critical sections cost $O(\log N)$ on the
domains of the library. The number of the critical sections is bounded
as well. A worker takes the lock of the graph at a check only after $m$
expansions in its nblock, or when its OPEN list is empty: in the first
case at most twice (lines 17 and 26), in the second once (line 26), and
the counter $k$ of line 5 starts again from 0 after each check. Hence,
if $E$ is the number of the expansions of the search and $Y$ the number
of the times that a worker finds its OPEN list empty (or empties it at
line 8), the lock is taken
at most $2 E / m + Y + p$ times ($p$ for the first call of NEXTNBLOCK of
each worker), the waits at line 33 and the creations of the nblocks
(whose lock is another one) apart. Note that this bound says nothing on
the time that a worker spends waiting for the lock, which depends on the
scheduling; the try-lock of line 25 is there to keep that time "small".

### The worker

The check and the choice of the next nblock follow Algorithm 11.2
(`should_switch` and `next_nblock`, `SafePBNF.h:541-606`). The lock of
the graph is the one of the executor (`Executor::lock_type`), so that
the same code runs on the threads of C++ and on OpenMP, and a worker
waits on a `std::condition_variable_any`, which works with any lock. The
expansion is that of the sequential A$^*$ of the library, in two phases:
the successors are generated first (without the move back to the parent,
in a `Domain` that knows its operators), each with its nblock and the
slot of its table computed and prefetched, i.e., requested from memory
without waiting for it, and the lookups follow, so that the cache misses
of the tables overlap (the design of [[4]](Z-bibliography.md#ref-4);
cf. [Chapter 3](03-data-structures.md#ch-3)). Lines 9 to 12 of
Algorithm 11.1 are thus the following (`expand` and `insert`,
`SafePBNF.h:608-665`).

**Algorithm 11.6** (Safe PBNF: the expansion of $u$, in the nblock $b$
held by the worker, with no lock but that of the incumbent).

```text
85  procedure EXPAND(u)
86     close u;  count an expansion, and a re-expansion if u was expanded
87     if u is a goal then
88        lock the incumbent;  if g(u) < z then z ← g(u);  record u
89        unlock the incumbent;  return
90     K ← ()                                     ▷ phase 1: generate, prefetch
91     for each successor v of u, but its parent if the Domain knows it, do
92        γ ← g(u) + c(u, v);  if γ + h(v) ≥ z then continue
93        a ← BLOCK(nblock(v));  hv ← hash(v)
          prefetch the slot hv of CLOSED(a)
94        append (v, γ, hv, a) to K
95     for each (v, γ, hv, a) ∈ K do               ▷ phase 2: the lookups
96        if INSERT(a, v, γ, hv) then publish f(a) ← min(f(a), γ + h(v))
97  function INSERT(a, v, γ, hv)
98     x ← the node of v in CLOSED(a), found from hv, or a new node of v
          with g = γ and parent u, taken from the arena of the worker
99     if x is not new then
100       if γ ≥ g(x) then count a duplicate;  return false
101       g(x) ← γ;  parent(x) ← u;  x is not closed      ▷ a reopening
102    push the entry (x, γ) on OPEN(a);  return true
```

By Lemma 11.3, $a \in S(b)$ at line 93, and therefore lines 98 to 102
write only lists that Theorem 11.8 reserves to the worker. The read of
$z$ at line 92 is relaxed, which may only prune less (cf.
[§8.3](08-shared-memory.md#sec-8-3), where pruning by a past value of
$z$ is shown to be sound). The insertion into a neighbouring nblock
lowers its published value with a relaxed atomic store at line 96,
which is safe since only the worker that holds the scope writes it (a
reader may see a "stale" value, as discussed below). Each worker
allocates its nodes from a private arena (`f_pool`, `SafePBNF.h:683`),
whose nodes are never moved. Hence, the parent pointers that cross the
workers stay valid until the end, and they are followed only after all
the workers have returned. The worker object is aligned to a cache line
(`SafePBNF.h:499`), and the counters of two workers never share a line.
After each expansion, a worker publishes the least $f$ of its OPEN list
(line 13) and compares the size of its arena with `node_limit`
(`SafePBNF.h:526-536`).

The cost of an expansion is thus that of the sequential A$^*$ with the
same lists: the generation of the successors, one call of BLOCK per
child (one atomic load once the nblock exists), one lookup of the table
and one push of OPEN per child that is not pruned, i.e., the costs of
[Chapter 3](03-data-structures.md#ch-3) for the table with open
addressing and for the bucket queue or the binary heap, the latter on
the list of one nblock, which is shorter than the single OPEN list of
A$^*$ (cf. section 11.9). The memory is the following.

<a id="prop-11-21"></a>**Proposition 11.21** (The memory of the
search). Let $B$ be the number of the nblocks created, $V$ the number of
the states reached and $J$ the number of the pushes on the OPEN lists.
Then the search uses $O(N + B \Delta + V + J)$ memory, $\Delta$ being as
in Proposition 11.20, and `solve()` takes $\Omega(N)$ time to initialize
the array $A$.

*Proof.* We note that $A$ has $N$ pointers, set to nil one by one
(`SafePBNF.h:310-314`). Each nblock created holds $S(b)$ and $I(b)$,
i.e., $O(\Delta)$ ids by Lemma 11.19, and constant data besides its
lists. By Corollary 11.9 each state has at most one node, allocated once
in the arena of a worker; the tables of the nodes hold one pointer per
node, with the load factor of [Chapter 3](03-data-structures.md#ch-3),
and an OPEN list holds at most one entry per push, the stale ones
included, until they are dropped. The free list has at most $B$
entries. $\square$

The term $N$ is the reason why the abstraction should not be made
"huge" without need: for the grids of the benchmark it is 62 500
pointers, for the 15-puzzle 4096, while an abstraction into, say,
$10^8$ ids would cost about 800 MB before the first expansion.

### Relaxed values

In fact, we read the published values $f(b)$ without the lock and with
relaxed order (lines 18 and 28, and line 47 for the hot nblocks), and
therefore a check may see a value slightly out of date (typically by a
few expansions). This changes the order of expansion and the choice of
the hot nblock, but none of Lemma 11.7, Theorem 11.8, Theorem 11.12,
Theorem 11.13 and Proposition 11.17, whose proofs read only the
counters, the flags and the lists under the lock (FREE computes $f(n)$
from OPEN($n$), line 74, not from the published value). It is one of
the reasons for the assumption (A3) of Theorem 11.14. Note also that
the value published at line 13 is the key of the first entry of
OPEN($b$), which may be a stale entry; since that key is the least of
all the entries of the list, the live ones included, the published value
is a lower bound on $f(b)$, never an overestimate of it.

## 11.8 The benchmark {#sec-11-8}

We compare Safe PBNF with HDA$^*$ on grids of the kind used by Burns et
al., which have a "natural" abstraction into nblocks. The instances are
three random grids of $5000 \times 5000$ cells, each cell blocked with
probability 0.35 and the free cells of cost 1, generated by
`GridMap::random` with the seeds 1, 2 and 5 (whose corners are
connected); the start is the bottom-left corner and the goal the
bottom-right one (the defaults of `AStarBench`). The moves are 8, a
diagonal move costing $\sqrt{2}$ times a straight one, with the octile
heuristic. Note that Burns et al. blocked a cell with probability 0.45
on their grids with 8 moves ([[3]](Z-bibliography.md#ref-3), section
4.3, p. 708), and the instances are not theirs. The nblocks are blocks
of $20 \times 20$ cells (62 500 nblocks, as in section 4.1 of
[[3]](Z-bibliography.md#ref-3)) and $m = 32$, the defaults of
`AStarBench` (`--nblock 20 --minexp 32`).

We run each search on 1, 2, 4, 8, 16, 32, 48, 64 and 96 cores, one
worker per core, pinned to the cores of one NUMA node with the memory on
that node, each run twice in separate processes, the gap between the two
runs being the A/A noise given in the note of the table. The rows are
Safe PBNF on the threads of C++ and HDA$^*$ with the distributions of
[Chapter 9](09-hash-distributed.md#ch-9) (Z, AZ and P on the threads, Z
and AZ on OpenMP, on FastFlow and over MPI), as `tools/bench/campaign.sh`
runs them. Since the costs are not integers, the searches keep their
OPEN lists in binary heaps, and the reference labelled "A$^*$ buckets",
the sequential A$^*$ of `AStarBench` on 1 core, runs on the heap as well
(`astar()` takes the bucket queue only for integer costs). The costs are
those of the domain, in double precision.

**Table 11.1.** Safe PBNF and HDA$^*$ on three random grids of
$5000 \times 5000$ cells with 8 moves. Each cell gives the seconds
summed over the three instances (each the median of its two runs) and,
in parentheses, the speedup over sequential A$^*$ on 1 core; an empty
cell is a search that did not solve every instance. The machine and the
compiler are given in the note of the table.

--8<-- "tables/safepbnf-grids.md"

The clock is the timer of the `Solution`: it starts at the beginning of
`solve()` and stops at its end, and therefore it includes the creation
of the nblock graph (the nblocks themselves are created during the
search), the start of the threads of the executor (`AStarBench` does not
lend one), the search and the rebuilding of the path. It does not
include the generation of the grid, the construction of the `Domain`,
the start of the process and the printing of the results. We leave the
table as a placeholder until the campaign is run on an idle machine.

## 11.9 Notes and references {#sec-11-9}

PBNF and Safe PBNF are due to Burns, Lemons, Ruml and Zhou
[[3]](Z-bibliography.md#ref-3), who also brought the abstraction into
the hash-based distribution of PRA$^*$ and HDA$^*$ (AHDA$^*$, treated in
[Chapter 9](09-hash-distributed.md#ch-9)) and compared the two families
on grids, sliding-tile puzzles and STRIPS planning; their appendix A is
the pseudocode followed by the library, and their appendix B the TLA$^+$
model of the hot nblocks. The scopes come from the SDD and the PSDD of
Zhou and Hansen, which Burns et al. summarize in their section 2.4, and
the ordering of the nblocks by their best node from the Localized A$^*$
of Edelkamp and Schrödl, which grouped the nodes by the page of memory
where they reside (both cited there; neither is in the bibliography of
this manual, since we could not verify a doi for them). We refer to the
survey of Fukunaga et al. [[1]](Z-bibliography.md#ref-1) for the place
of PBNF among the parallel best-first searches.

Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) revisited the
comparison of Burns et al. with HDA$^*$ on Zobrist hashing, on "easy"
and "hard" instances, and with both the heap and the bucket queue as
OPEN list (section 3.2, pp. 507-510). They observed that HDA$^*$[Z] is
faster than Safe PBNF on the instances of the 15-puzzle and of the
24-puzzle that need "much" search, and Safe PBNF on those solved
quickly, because of the start-up burst of HDA$^*$ (Observation 3 and
Figure 5, p. 509). On grid pathfinding, instead, Safe PBNF outperformed
HDA$^*$[Z] although it expanded more nodes, because of the communication
overhead of the latter (Observation 4, p. 512). They also remarked that
PBNF gains less than HDA$^*$ from a bucket queue in place of a binary
heap, since each nblock has its own queue: the difference is $O(\log (N
/ B))$ against $O(1)$ per insertion, with $N$ states and $B$ nblocks (p.
509). Finally, they attributed the smaller burst of Safe PBNF to the
exclusive access to the scopes (section 3.1.2, p. 505). We implement
both searches with the same OPEN lists and the same tables, on the same
executors, so that the comparison of Table 11.1 is not biased by the
data structures in the way that Jinnai and Fukunaga pointed out.

Two limits of the method are stated by Burns et al. in their discussion
(section 9.1, p. 735): the abstraction is chosen by hand, and "it seems
likely that abstraction size should change when number of threads
changes or perhaps even based on features of the domain or problem
instance". In the library the abstraction is given by the `Domain`, and
its size by its parameters (e.g., `nblock_side` for the grids), which
one has to choose by hand as well. The automatic construction of
abstractions for HDA$^*$ is the subject of [Chapter
10](10-auto-distributions.md#ch-10), and it does not cover the nblocks.
Also, the library supports only undirected abstract graphs (condition
(ii) of Definition 11.1), while Burns et al. treat directed ones with
two further cases of SETHOT, and it offers neither the
bounded-suboptimal nor the anytime variants of PBNF of their sections 5
and 7, which belong with [Chapter 6](06-bounded-suboptimal.md#ch-6).

## 11.10 Exercises {#sec-11-10}

**Exercise 11.1** (`ctest -R PBeFS.pbnf-options`). In the vacuum world
of $n$ cells, write $S(b)$ and $I(b)$ for each nblock $b$, and show that
$|I(b)| \le 4$ and that at most $\lceil n / 3 \rceil$ workers can hold
nblocks at the same time. Check that `SafePBNF` returns the cost of
A$^*$ on the vacuum world with 1, 3 and 8 workers and each value of
$m$ of the test.

**Exercise 11.2** (`ctest -R PBeFS.pbnf-thread`). On a grid of
$W \times H$ nblocks with 4 moves, show that the nblocks whose two
coordinates are multiples of 3 have pairwise disjoint scopes, and hence
that $\lceil W / 3 \rceil \lceil H / 3 \rceil$ workers can be busy at
once, and that the same holds with 8 moves. Then count $|I(b)|$ on the
border and in a corner of the grid. The test runs `SafePBNF` on grids of
$300 \times 200$ cells with nblocks of $10 \times 10$ cells and 1, 2, 3
and 8 workers.

**Exercise 11.3** (`ctest -R PBeFS.pbnf-thread`). Prove directly from
Algorithm 11.3 that no reachable state of the nblock graph has two hot
nblocks that interfere, and that line 57 is never executed. Which lines
would have to change for the library to accept a directed abstract
graph?

**Exercise 11.4** (`ctest -R PBeFS.pbnf-thread`). The test runs
`SafePBNF` with 4 workers on a $60 \times 60$ grid split by a wall of
blocked cells, with nblocks of $7 \times 7$ cells, and expects the
status `unsolvable`. Follow the proof of Theorem 11.12 on this instance:
which nblocks are ever acquired, and why does the search end at line 32
rather than sleep forever at line 33?

**Exercise 11.5** (`ctest -R PBeFS.pbnf-thread`). The node limit of
`SafePBNF` is per worker, and the worker that reaches it sets the flag
of the end under the lock of the graph and wakes all the others. Show
that all the workers then return, whatever it was doing, and that the
status is `node_limit` even if a goal was found before. The test checks
a limit of 5 nodes on an instance of the 8-puzzle that needs at least 6
moves.

**Exercise 11.6** (`ctest -R PBeFS.pbnf-options`). Show that the cost
returned does not depend on $m$ (Theorem 11.12), and explain why
`SafePBNF` takes $m = 0$ as $m = 1$: which step of the proof of Theorem
11.13 fails if a worker may switch before its first expansion in an
nblock? Discuss the trade-off of $m$ between the time spent on the lock
and the "speculative" expansions, as measured by Burns et al. (section
4.1 of [[3]](Z-bibliography.md#ref-3)).

**Exercise 11.7** (`ctest -R PBeFS.pbnf-thread`). Show that plain PBNF
with one worker never livelocks (which is rather easy), and build an
execution of plain PBNF with three workers on a cycle of nblocks in
which the nblock with the best node is never free. Then run Safe PBNF on
your execution by hand, as in Figures 11.9 and 11.10.

**Exercise 11.8** (`ctest -R PBeFS.pbnf-openmp`). List the parts of
`SafePBNF` that depend on the executor, and show that the proofs of
section 11.6 hold for any executor whose lock provides mutual exclusion
and a try-lock. The test runs the cases of `pbnf-thread` on OpenMP.

**Exercise 11.9** (`ctest -R PBeFS.pbnf-thread`). Give a `Domain` whose
abstraction violates condition (ii) of Definition 11.1 (a directed
abstract graph), and two nblocks that interfere although the
interference scope built by `SafePBNF` (Lemma 11.4) does not record it.
Show an interleaving in which two workers write the same table at the
same time, and explain why the tests, which use only undirected
abstractions, cannot detect the defect.
