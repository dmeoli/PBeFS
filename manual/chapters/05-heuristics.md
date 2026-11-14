# 5. Heuristics and pattern databases {#ch-5}

[Source: `include/PBeFS/astar/PatternDatabase.h` (`BlankRule`,
`PDBOptions`, `TilePattern`, the encodings, `PatternTable`,
`AdditivePDB`, `SlidingTilePDB`), `include/PBeFS/astar/SlidingTile.h`,
`include/PBeFS/astar/Grid.h`, `include/PBeFS/astar/TravelingSalesman.h`,
`include/PBeFS/astar/SequenceScore.h`; the tests `pdb`, `pdb-search`,
`tiles24` and `grid` of `test/test.cpp`; `tools/PDBBench.cpp` and
`tools/bench/pdb.sh`]

We are concerned in this chapter with the heuristic $h$ that a domain
gives to the searches, i.e., with how admissible and consistent
heuristics are obtained (the definitions and their consequences on A$^*$
are those of [§2.4](02-best-first.md#sec-2-4) to
[§2.7](02-best-first.md#sec-2-7)). Most of the chapter is devoted to the
pattern databases of the sliding-tile puzzles, which store the exact
distances of an abstraction of the puzzle in a table: the databases of
Culberson and Schaeffer [[21]](Z-bibliography.md#ref-21), summed over a
partition of the tiles as in the disjoint databases of Korf and Felner
[[22]](Z-bibliography.md#ref-22) and of Felner, Korf and Hanan
[[23]](Z-bibliography.md#ref-23). We first derive the heuristics of the
domains of the library from relaxations and abstractions. Then we define
the abstraction of a pattern under the two rules of the blank of the
library, rank its states and build its table by a backward breadth-first
search. We prove that the sum of the tables of a partition and its
maximum with the reflection about the main diagonal are admissible, and
we establish when they are consistent; then we give the encodings of the
entries, from one byte down to one bit, with the argument by which the
relative encodings recover the exact value. We end the chapter with the
files of the tables, their sizes, a guide to the searches whose
heuristic does not help, and the benchmark.

## 5.1 Relaxations and abstractions {#sec-5-1}

A heuristic is usually the exact cost of a simpler problem. In a
*relaxation*, the problem is given more edges or cheaper ones, and each
path of the original graph is still a path, of no larger cost; hence the
optimal cost of the relaxation from $u$ is a lower bound on $h^*(u)$.
For instance, the Manhattan distance of the sliding-tile puzzles is the
optimal cost of the relaxation in which a tile may move onto an occupied
cell. An *abstraction* maps the states to fewer abstract states and
keeps the edges between their images; the distances of the abstract
graph, which may be small enough to be computed once and stored, are
again lower bounds. We state the property in general, since the pattern
databases of the next sections are abstractions of this kind.

<a id="def-5-1"></a>**Definition 5.1** (Abstraction). Let
$G = (V, E)$ be the graph of a domain, with costs $c$ and goals
$T \subseteq V$. An *abstraction* is a map $\psi : V \to V'$ onto the
states of a graph $G' = (V', E')$ with costs $c'$ and goals $T'$, such
that (i) $\psi(t) \in T'$ for each $t \in T$, and (ii) for each edge
$(u, v) \in E$, either $\psi(u) = \psi(v)$, or $(\psi(u), \psi(v)) \in
E'$ with $c'(\psi(u), \psi(v)) \le c(u, v)$. The *abstraction
heuristic* is $h_\psi(u) = d'(\psi(u))$, the least cost in $G'$ of a
path from $\psi(u)$ to a goal of $T'$.

<a id="prop-5-2"></a>**Proposition 5.2** (Abstraction heuristics are
consistent). The heuristic $h_\psi$ of an abstraction is consistent,
and hence admissible.

*Proof.* At a goal $t$, $\psi(t) \in T'$ and $h_\psi(t) = 0$. For an
edge $(u, v)$ with $\psi(u) = \psi(v)$, $h_\psi(u) = h_\psi(v) \le
c(u, v) + h_\psi(v)$; otherwise, a path of $G'$ from $\psi(v)$ to a
goal, preceded by the edge $(\psi(u), \psi(v))$, gives $h_\psi(u) \le
c'(\psi(u), \psi(v)) + h_\psi(v) \le c(u, v) + h_\psi(v)$. Admissibility
follows from [Lemma 2.4](02-best-first.md#lem-2-4). $\square$

Of course, the more states $\psi$ merges, the smaller $G'$ and the
(typically) weaker the bound, and the art of an abstraction is to keep the
part of the state that makes the problem hard. Two abstractions combine
by the maximum, which keeps both properties ([Lemma
2.5](02-best-first.md#lem-2-5)); they combine by the sum only when they
do not count the same cost twice, which is the subject of
[§5.6](05-heuristics.md#sec-5-6).

## 5.2 The heuristics of the domains of the library {#sec-5-2}

Each domain of [Chapter 14](14-domains.md#ch-14) comes with one of the
following heuristics, each of which is the exact cost of a relaxation.

* *Sliding-tile puzzles* (`SlidingTile.h`): the Manhattan distance, the
  sum over the tiles of the distance, in moves of a rook of one cell,
  between the cell of the tile and its goal cell. A move changes the
  cell of one tile by one unit, and hence the distance by exactly one,
  at the cost of one. Thus the heuristic is consistent, it has the
  parity of the true distance, and the domain computes it incrementally,
  from the value of the parent passed to `for_each_successor`. The
  number of the tiles out of their goal cell (the heuristic $h_1$ of the
  example `8Puzzle`) is also consistent, and weaker.
* *Grids* (`Grid.h`): with 4 moves, the Manhattan distance between the
  cell and the goal; with 8 moves, the *octile* distance $\max(\Delta x,
  \Delta y) + (\sqrt{2} - 1)\min(\Delta x, \Delta y)$, the cost of the
  path with as many diagonal moves as possible. Both are multiplied by
  the least cost of a cell (`Grid.h:444`), since a move pays the cost of
  the cell it leaves, times $\sqrt{2}$ on a diagonal. Both are
  consistent, being the exact distances of the grid without walls and
  with each cell at the least cost.
* *Traveling salesperson* (`TravelingSalesman.h`): the cost of a minimum
  spanning tree of the cities still to be visited, together with the
  current one and the first (cf. exercise 3.30 of
  [[13]](Z-bibliography.md#ref-13)). The rest of the tour is a path
  through all these cities, and a path costs at least as much as a
  spanning tree of its cities; hence the heuristic is admissible.
* *Explicit graphs* (`ExplicitGraph.h`): the heuristic given node by
  node by `set_heuristic`, e.g., the straight-line distance of the map
  of Romania [[13]](Z-bibliography.md#ref-13), which is consistent when
  each road is at least as long as the segment between its ends.

On a grid with 8 moves the Manhattan distance is *not* admissible,
whatever the cost of a diagonal move between 1 and $\sqrt{2}$. For
instance, three diagonal moves of unit cells cost $3\sqrt{2} \approx
4.24$, while the Manhattan distance of their ends is 6 ([Figure
5.1](05-heuristics.md#fig-5-1)); the Euclidean distance is admissible
there, and consistent, but it is weaker than the octile distance, which
is the distance of the relaxation without walls. Finally, the sequence
score of Nilsson (`SequenceScore.h`) adds to the Manhattan distance a
term that counts the tiles out of sequence; it is not admissible, and
A$^*$ then usually expands fewer nodes, but it may return a path that is
not optimal ([§2.4](02-best-first.md#sec-2-4)).

![Figure 5.1. On a grid with 8 moves, a path of three diagonal moves
(violet) costs $3\sqrt{2} \approx 4.24$, which the octile distance gives,
while the Manhattan distance of its ends, the cost of the path
with 4 moves (dashed), is 6: the Manhattan distance overestimates,
and it is not admissible with 8
moves.](../figures/heuristics-octile.svg){#fig-5-1}

## 5.3 The abstraction of a pattern {#sec-5-3}

In a sliding-tile puzzle, the board has $n = W H$ cells, and a state
$s$ places the tiles $1, \ldots, n - 1$ and the blank on distinct
cells; a move exchanges the blank with a tile in an adjacent cell, at
cost 1. A *pattern* is a set $\pi = \{t_1, \ldots, t_k\}$ of tiles,
and the *placement* of $\pi$ in $s$ is the vector
$p = \psi_\pi(s) = (q_{t_1}, \ldots, q_{t_k})$ of the cells of its
tiles, the other tiles being "don't care"
([Figure 5.2](05-heuristics.md#fig-5-2)).
Thus, the placements are the $k$-permutations of the $n$ cells, i.e.,
$N = n! / (n - k)!$ of them, and the *goal placement* is the placement
of the goal. In the library, `TilePattern` (`PatternDatabase.h:233`) takes
a pattern as the list of its tiles, in a chosen order, and the
board of the goal; its constructor refuses a pattern with the blank,
an empty one, or one with a tile twice (`PatternDatabase.h:282` and `:287`).

![Figure 5.2. The abstraction of the pattern $\pi = \{1, 2, 3\}$ of
the 8-puzzle: the states $s$ and $s'$, which differ in the other tiles
and in the blank, have the same placement $p = (1, 4, 2)$ of the tiles
of $\pi$ (the cells numbered row by row from 0); the table of $\pi$
holds the distance from each placement to the goal
placement.](../figures/heuristics-abstraction.svg){#fig-5-2}

The *abstract space* of $\pi$ is a graph on the placements whose edges
move one tile of $\pi$ to an adjacent cell, at cost 1. As [Figure
5.3](05-heuristics.md#fig-5-3) shows, the library builds it under two
rules of the blank (`BlankRule`, `PatternDatabase.h:196`).

* With the blank *ignored*, an abstract state is a placement, and each
  move of a tile of $\pi$ to an adjacent cell free of the other tiles
  of $\pi$ is allowed, as if the blank were wherever needed. The table
  holds the distance $d_\pi(p)$ from $p$ to the goal placement in this
  graph, which is undirected and has unit costs.
* With the blank *minimized*, an abstract state is a placement together
  with the cell $b$ of the blank, a move of the blank into a cell that
  no tile of $\pi$ holds costs 0 (it moves a tile out of the pattern),
  and a move of a tile of $\pi$ into the cell of the blank costs 1. The
  blank thus reaches at no cost each cell of its *region*, i.e., of
  its connected component among the cells free of $\pi$, and the table
  holds the least distance over the cells of the blank,
  $d_\pi(p) = \min_b d_\pi(p, b)$; this is the rule of the disjoint
  databases of [[22]](Z-bibliography.md#ref-22).

![Figure 5.3. The moves of the abstraction of
$\pi = \{1, 2, 3\}$. (a) With the blank ignored, a tile of $\pi$ moves
to any adjacent cell free of the other tiles of $\pi$. (b) With the
blank minimized, the blank (top left) reaches its region (shaded) at no
cost, and a tile of $\pi$ moves only into a cell of that region; the
moves in rust would need the blank in the other region, which it cannot
reach without moving a tile of
$\pi$.](../figures/heuristics-blank.svg){#fig-5-3}

In both cases a move of the puzzle that moves a tile of $\pi$ is an edge
of cost 1 of the abstract space. A move of any other tile changes no
placement (with the blank ignored) or is an edge of cost 0 (with the
blank minimized); hence $\psi_\pi$, with the cell of the blank added in
the second case, is an abstraction in the sense of Definition 5.1. We
now compare the two rules.

<a id="lem-5-3"></a>**Lemma 5.3** (The minimized blank dominates).
For each placement $p$, $d_\pi(p)$ with the blank minimized is at
least $d_\pi(p)$ with the blank ignored.

*Proof.* A path of the minimized space from $(p, b)$ to the goal,
deprived of its moves of the blank, is a sequence of moves of the tiles
of $\pi$, each into a cell free of the other tiles of $\pi$, i.e., a
path of the space with the blank ignored, of the same cost. Hence the
distance with the blank ignored is at most $d_\pi(p, b)$ for each $b$,
and at most their minimum. $\square$

Note that the test `pdb` checks Lemma 5.3, and the values of the tables against
a breadth-first search of the abstraction written apart from the
library, on each placement of the 8-puzzle, with two goals and both
rules.

## 5.4 Ranking the placements {#sec-5-4}

A table is an array indexed by the placements, which need a *perfect
ranking*, i.e., a bijection $\mathrm{rank}$ from the placements onto
$[0, N)$, computed at each lookup. We rank a placement as a
$k$-permutation in lexicographic order: the digit of the $i$-th tile
is the number $d_i$ of the cells below $p_i$ that the tiles before it
do not hold, and

$$
\mathrm{rank}(p) = \sum_{i=1}^{k} d_i \prod_{j=i+1}^{k} (n - j + 1),
$$

which is the number $d_1 d_2 \cdots d_k$ written in the mixed radix
$(n, n - 1, \ldots, n - k + 1)$, the $i$-th digit $d_i$ being in
$[0, n - i + 1)$ since $p_i$ has at most $n - i + 1$ free cells below
it, itself excluded ([Figure 5.4](05-heuristics.md#fig-5-4)).

<a id="lem-5-4"></a>**Lemma 5.4** (The rank is a bijection). The map
$\mathrm{rank}$ is a bijection from the $k$-permutations of $n$ cells
onto $[0, n! / (n - k)!)$.

*Proof.* The mixed radix representation is a bijection from the
vectors of digits with $d_i \in [0, n - i + 1)$ onto
$[0, \prod_i (n - i + 1)) = [0, n! / (n - k)!)$. The map from a
placement to its digits is a bijection too: given the cells
$p_1, \ldots, p_{i-1}$, the cell $p_i$ is the $(d_i + 1)$-th of the
$n - i + 1$ cells that they leave free, in increasing order, so that
$p$ is rebuilt from its digits one cell after the other. $\square$

In `TilePattern::rank` and `unrank`, the code keeps the cells already
used in a word of 32 bits. Thus the digit is a population count of the
used cells below $p_i$, and the inverse takes a division and the
selection of the $d_i$-th free bit for each tile:

```cpp
--8<-- "snippets/heuristics-rank.cpp"
```

Two placements whose ranks differ by less than $n - k + 1$ and lie in
the same block of that size differ only in the cell of the last tile,
which the max-compression of [§5.9](05-heuristics.md#sec-5-9) exploits.
The linear-time ranking of Myrvold and Ruskey
[[34]](Z-bibliography.md#ref-34) would do as well, and it does without
the population counts. We chose the lexicographic one because it keeps
close placements in close entries, and its cost of $k$ population counts
is small against the cache miss that the lookup of a large table pays
anyway.

![Figure 5.4. The rank of the placement $p = (4, 0, 7)$ of three
tiles on 9 cells: each row is the board before the tile $i$ is placed,
with the cells of the tiles before it shaded; the digit $d_i$ is the
number of the free cells below $p_i$, and the rank is accumulated in
the mixed radix $(9, 8, 7)$, giving
$(4 \cdot 8 + 0) \cdot 7 + 5 = 229$.](../figures/heuristics-rank.svg){#fig-5-4}

## 5.5 Building a table {#sec-5-5}

Each table is filled by a breadth-first search of the
abstract space, backward from the goal placement and by layers: at the
layer $d$ each entry of value $d$ is expanded, and its neighbours not
seen yet get the value $d + 1$; the search ends when a layer adds no
entry. Since the abstract space with the blank ignored is undirected,
the backward distances are the forward ones. We give the pseudocode of
the Doxygen comment of `PatternDatabase.h`.

**Algorithm 5.1** (BUILD-IGNORED, the table of a pattern with the blank
ignored).

```text
 1  procedure BUILD-IGNORED(π)               ▷ the blank ignored
 2     for each r ∈ [0, N) do h[r] ← ∞          ▷ N = n! / (n − k)!
 3     h[rank(goal(π))] ← 0;  d ← 0
 4     loop
 5        new ← 0
 6        for each r with h[r] = d do        ▷ by blocks, in parallel
 7           p ← unrank(r)
 8           for each tile i of π, each cell c next to p_i, c ∉ p do
 9              r' ← rank(p with p_i ← c)
10              if h[r'] = ∞ then          ▷ by compare-and-swap
11                 h[r'] ← d + 1;  new ← new + 1
12        if new = 0 then return h
13        d ← d + 1
```

[Figure 5.5](05-heuristics.md#fig-5-5) follows the algorithm on the
board of $3 \times 2$ cells with the pattern $\{1, 2\}$, whose 30
placements are laid out by rank as a $6 \times 5$ matrix: the layer 0
is the goal placement, the layer 1 its three neighbours, and the table
is complete after the layer 5.

![Figure 5.5. Algorithm 5.1 on the board of $3 \times 2$ cells with $\pi
= \{1, 2\}$: the 30 entries, in the order of their ranks (the row is the
cell of the tile 1, the column the digit of the tile 2), are $\infty$
until reached. (a) The layer 0 (violet) is expanded and its neighbours
get 1 (rust); (b) the layer 1 is expanded and gives the layer 2; (c) the
table at the end, whose largest distance is 5.
](../figures/heuristics-bfs.svg){#fig-5-5}

Each layer is scanned in blocks of $2^{14}$ consecutive ranks, which the
workers of a runtime take from a shared counter. An entry is claimed by
an atomic compare-and-swap from $\infty$ to $d + 1$ (an entry of one
byte, with the value 255 standing for $\infty$):

```cpp
--8<-- "snippets/heuristics-build.cpp"
```

<a id="lem-5-5"></a>**Lemma 5.5** (The parallel construction is
deterministic). Algorithm 5.1, with the layers scanned by any number of
workers in any interleaving, gives each entry the distance $d_\pi$ of
its placement from the goal placement.

*Proof.* We show by induction on $d$ that, at the end of the layer $d -
1$, the entries of value at most $d$ are the placements at distance at
most $d$, with their distances. During the layer $d$, the workers read
the entries equal to $d$ and write only the value $d + 1$ into entries
equal to $\infty$. Hence the set of the entries equal to $d$ does not
change during the layer; each neighbour of such an entry that was
$\infty$ at the start of the layer is set to $d + 1$ (by one
compare-and-swap, the others failing); and no other entry changes. These
are the placements at distance $d + 1$, by the induction hypothesis and
the definition of the breadth-first layers. Finally, the search ends at
the first layer that adds no entry, i.e., when each placement reachable
from the goal has its distance. With the blank ignored each placement is
reachable, which the code checks, and a distance beyond 254 is refused
(`PatternDatabase.h:1490`). $\square$

Indeed, the test `pdb` checks that the tables built with 1 and with 3
workers are equal. With the blank minimized the states are the pairs
$(p, b)$, with $b$ coded by its index among the $n - k$ cells free of
the pattern. Each state takes 2 bits (unseen, current, next, old), as in
the two-bit breadth-first search of Korf
[[26]](Z-bibliography.md#ref-26). A state is expanded once per region,
from the least cell of the region, and its successors are marked
together with all the cells of the region of the blank after the move,
which realizes the moves of cost 0 without a second queue.

**Algorithm 5.2** (BUILD-MINIMIZED, the table of a pattern with the
blank minimized).

```text
 1  procedure BUILD-MINIMIZED(π)             ▷ the blank minimized
 2     S[r, b] ← unseen for each r and each cell b ∉ unrank(r)
 3     S[rank(goal(π)), b] ← current for each b ∉ goal(π)
 4     h[r] ← ∞ for each r;  h[rank(goal(π))] ← 0;  d ← 0
 5     loop
 6        for each (r, b) with S[r, b] = current do    ▷ in parallel
 7           p ← unrank(r);  R ← region(p, b)  ▷ the cells b reaches
 8           if b ≠ min R then continue      ▷ one cell per region
 9           for each c ∈ R, each tile i of π with p_i next to c do
10              p' ← p with p_i ← c;  r' ← rank(p')
11              for each b' ∈ region(p', old p_i) do
12                 if S[r', b'] = unseen then S[r', b'] ← next
13              if h[r'] = ∞ then h[r'] ← d + 1
14        if no state is next then return h
15        S ← S with current → old and next → current;  d ← d + 1
```

Line 3 starts from the goal placement with the blank anywhere, since the
backward search looks for the states from which the goal placement is
reached, with the blank wherever it ends. Line 11 marks the whole region
of the blank after the move, which contains the cell left by the tile,
since the states of a region are at the same distance (the blank moves
within it at no cost). Note that the first value written in $h[r']$ is
the least distance over the cells of the blank, since the layers are
visited by increasing distance. Overall, the construction takes $O(N k)$
time with the blank ignored and $O(N (n - k) k)$ with the blank
minimized; it takes $N$ bytes, plus $N (n - k) / 4$ bytes for the states
of the second rule. Also, the patterns of a partition are built one
after the other, each with all the workers (`PDBOptions::threads`, all
the cores if 0).

## 5.6 Additive pattern databases {#sec-5-6}

Clearly, the table of one pattern is a rather weak heuristic, since it
ignores the other tiles. Taking the maximum of the tables of several
patterns is admissible ([Lemma 2.5](02-best-first.md#lem-2-5)), but it
is not much stronger than the best of them; the sum is much stronger,
and it is admissible when the patterns are disjoint and each table
counts only the moves of its own tiles. Let $\pi_1, \ldots, \pi_m$ be
disjoint patterns (a partition of the tiles, or of a subset of them;
[Figure 5.6](05-heuristics.md#fig-5-6)), and let

$$
h(s) = \sum_{j=1}^{m} d_{\pi_j}(\psi_{\pi_j}(s)).
$$

<a id="thm-5-6"></a>**Theorem 5.6** (Admissibility of the additive
heuristic). If the patterns are disjoint and none holds the blank, $h$
is admissible, with either rule of the blank, and $h = 0$ at the goal.

*Proof.* Consider a path $P$ from $s$ to the goal, and let $c_j(P)$ be
the number of its moves that move a tile of $\pi_j$. The projection of
$P$ on $\pi_j$ (with the cell of the blank, under the minimized rule) is
a path of the abstract space from $\psi_{\pi_j}(s)$ to the goal
placement whose cost is $c_j(P)$, since the moves of the other tiles
cost nothing there. Hence $d_{\pi_j}(\psi_{\pi_j}(s)) \le c_j(P)$ (with
the minimized rule, $d_{\pi_j}$ is the minimum over the cells of the
blank, and the projection starts from one of them). Each move of the
puzzle moves one tile, which belongs to at most one pattern; thus
$\sum_j c_j(P) \le c(P)$, and summing over $j$ gives $h(s) \le c(P)$ for
each path, i.e., $h(s) \le h^*(s)$. At the goal each placement is the
goal placement, of distance 0. $\square$

This is the argument of Korf and Felner
[[22]](Z-bibliography.md#ref-22), which Yang et al.
[[28]](Z-bibliography.md#ref-28) cast for any family of abstractions
whose costs split the cost of each edge. This argument needs the
patterns to be disjoint and the abstract spaces to charge only the moves
of their tiles. For this reason the blank belongs to no pattern: a table
that charged the moves of the blank would count each move twice, once
for the blank and once for the tile. `AdditivePDB`
(`PatternDatabase.h:1681`) refuses a tile in two patterns, and
`TilePattern` refuses the blank. In fact, the Manhattan distance is the
additive database of the $n - 1$ patterns of one tile each, whose
abstract space with the blank ignored is the board. A partition into
larger patterns dominates it, since merging two patterns can only
increase the sum (the distance of two tiles together is at least the sum
of their separate distances, by the same projection argument). The
partitions of [Figure 5.6](#fig-5-6) are given by the functions
`partition_15_7_8()`, `partition_15_5_5_5()`, `partition_15_6_6_3()` and
`partition_24_6_6_6_6()` of `PatternDatabase.h`.

![Figure 5.6. The partitions of the tiles of the library, each
pattern in a colour: 7-8 of the 15-puzzle
[[22]](Z-bibliography.md#ref-22), 5-5-5 and 6-6-3 of the 15-puzzle
[[23]](Z-bibliography.md#ref-23), and 6-6-6-6 of the 24-puzzle
[[22]](Z-bibliography.md#ref-22); the blank, in the first cell of the
goal, belongs to no pattern.](../figures/heuristics-partition.svg){#fig-5-6}

The lookup of a state gathers, for each pattern, the cells of its tiles
from the array $q$ of the cells of all the tiles, and it sums the values
of the tables:

```cpp
--8<-- "snippets/heuristics-lookup.cpp"
```

and `SlidingTilePDB` (`PatternDatabase.h:1919`) is the domain
`SlidingTile` whose heuristic is that of an `AdditivePDB`, so that each
search of the library runs on it unchanged. Its `for_each_successor`
computes the heuristic of a successor by moving one tile in the array of
the cells of the parent and looking it up again (or, with a relative
encoding, by the change of one table, [§5.9](05-heuristics.md#sec-5-9)):

```cpp
--8<-- "snippets/heuristics-child.cpp"
```

Here $q$ holds the cells of the tiles of the parent, the tile $x$ moves
into the cell $b$ of the blank of the parent, and the blank (tile 0)
moves to the cell that $x$ leaves. The following use comes from the
README of the repository. Its options of the construction (`PDBOptions`,
`PatternDatabase.h:209`) are the rule of the blank, the reflection, the
compression, the workers and the directory of the cached tables.

```cpp
PDBOptions o;
o.blank = BlankRule::minimized;
o.cache = "pdb-cache";                  // built once, then loaded
const AdditivePDB< 4 , 4 > pdb( partition_15_7_8() ,
                                FifteenPuzzle::korf_goal() , o );
const SlidingTilePDB< 4 , 4 > d( FifteenPuzzle::parse( "14 13 15 7 11 12 "
                                 "9 5 6 0 2 1 4 8 10 3" ) , pdb );
const auto s = AStar< SlidingTilePDB< 4 , 4 > >( d ).solve();  // cost 57
```

## 5.7 Consistency {#sec-5-7}

With the blank ignored the additive heuristic is consistent, and the
proof gives a finer property of the tables, on which the encodings of
[§5.9](05-heuristics.md#sec-5-9) rest. Let $m_\pi(p)$ be the Manhattan
distance of the tiles of $\pi$ in the placement $p$.

<a id="lem-5-7"></a>**Lemma 5.7** (Parity and excess). With either
rule of the blank, $d_\pi(p) \ge m_\pi(p)$ and
$d_\pi(p) \equiv m_\pi(p) \pmod 2$ for each placement $p$. With the
blank ignored, two adjacent placements $p$ and $p'$ have
$|d_\pi(p) - d_\pi(p')| = 1$.

*Proof.* A move of cost 1 moves one tile of $\pi$ by one cell, and
hence changes $m_\pi$ by exactly one unit, while a move of cost 0 (of
the blank, with the minimized rule) does not change it. Along a path of
least cost from $p$ to the goal placement, where $m_\pi = 0$, the value
of $m_\pi$ thus goes to 0 by $d_\pi(p)$ steps of $\pm 1$, which gives
both claims of the first sentence. With the blank ignored, the abstract
space is undirected with unit costs, so that
$|d_\pi(p) - d_\pi(p')| \le 1$ for adjacent placements, and the parity
of $d_\pi$, which is that of $m_\pi$, differs between them, since
$m_\pi(p)$ and $m_\pi(p')$ differ by one; hence the difference is not
zero. $\square$

<a id="thm-5-8"></a>**Theorem 5.8** (Consistency with the blank
ignored). With the blank ignored, the additive heuristic of disjoint
patterns is consistent.

*Proof.* A move of the puzzle moves one tile; if it belongs to no
pattern, no placement changes and $h$ does not change; otherwise only
the placement of its pattern changes, to an adjacent one, and by
[Lemma 5.7](05-heuristics.md#lem-5-7) $h$ changes by one. In either
case $h(s) \le 1 + h(s')$ for each move $(s, s')$, and $h = 0$ at the
goal by [Theorem 5.6](05-heuristics.md#thm-5-6). $\square$

With the blank minimized, the same move may change the value of a table
by more than one. In fact, the minimum over the cells of the blank
forgets in which region the blank actually is, and the region that the
true blank occupies after the move may be far from the one that realized
the minimum before it. Hence, the sum is admissible but usually not
consistent. On the 8-puzzle with the partition 4-4 it is inconsistent
along about 37000 of its moves (the count printed by the test `pdb`),
and A$^*$ stays optimal by reopening the closed nodes
([§2.7](02-best-first.md#sec-2-7) and [[31]](Z-bibliography.md#ref-31)).
In exchange, the heuristic is larger: by [Lemma
5.3](05-heuristics.md#lem-5-3) each table with the blank minimized is at
least the one with the blank ignored, which in turn is at least the
Manhattan distance of its tiles ([Lemma 5.7](05-heuristics.md#lem-5-7));
`PDBBench --mode table` reports the mean of each heuristic on $10^6$
random states (column `mean_h`), against that of the Manhattan distance
(column `mean_manhattan`). By default, `PDBOptions` has the blank
ignored, which all the encodings accept and which keeps the heuristic
consistent; the sequential A$^*$ then never expands a node twice. One
chooses the blank minimized explicitly, as in the example of
[§5.6](05-heuristics.md#sec-5-6).

## 5.8 The reflection about the main diagonal {#sec-5-8}

Let the board be square, and let the goal have the blank on the main
diagonal, as the goal of Korf (the blank in the first cell) does. The
transposition $T$ of the cells, $T(\text{row } i, \text{column } j) =
(\text{row } j, \text{column } i)$, maps adjacent cells to adjacent
cells. Let $\sigma(t)$ be the tile that the goal holds in the cell
$T(c^*_t)$, with $c^*_t$ the goal cell of $t$, and let $\phi(s)$ be the
state that holds the tile $\sigma(t)$ in the cell $T(q_t)$, for each
tile $t$ of $s$ (the blank included,
[Figure 5.7](05-heuristics.md#fig-5-7)).

<a id="lem-5-9"></a>**Lemma 5.9** (The reflection is an automorphism).
The map $\phi$ is a bijection of the states, maps the moves of the
puzzle to moves, and fixes the goal; hence $h^*(\phi(s)) = h^*(s)$ for
each state $s$.

*Proof.* Since $T$ and $\sigma$ are bijections, so is $\phi$. The blank
is fixed by $\sigma$, since $T(c^*_0) = c^*_0$ when the goal cell of
the blank is on the diagonal. A move of $s$ exchanges the blank, in
the cell $b$, with the tile $x$ in an adjacent cell $a$; in $\phi(s)$
the blank is in $T(b)$ and $\sigma(x)$ in $T(a)$, which is adjacent to
$T(b)$, and exchanging them gives $\phi$ of the successor. The goal
holds $t$ in $c^*_t$, and $\phi$ puts in $T(c^*_t)$ the tile
$\sigma(t)$, which is by definition the tile of the goal in
$T(c^*_t)$, i.e., it gives the goal back. Hence $\phi$
maps the paths from $s$ to the goal onto the paths from $\phi(s)$ to
the goal, of the same costs. $\square$

Thus, the sum of the same tables on $\phi(s)$ is another admissible
heuristic of $s$, since $h(\phi(s)) \le h^*(\phi(s)) = h^*(s)$, and it
is consistent if $h$ is, since $\phi$ maps moves to moves; the library
takes

$$
h_\phi(s) = \max \{ h(s), h(\phi(s)) \},
$$

which is the reflection of [[22]](Z-bibliography.md#ref-22), admissible
and, with the blank ignored, consistent by
[Lemma 2.5](02-best-first.md#lem-2-5). The reflection costs a second
lookup per table; it is used only when it is valid, i.e., on a square
board whose goal has the blank on the diagonal (`AdditivePDB::reflected()`),
and the constructor precomputes $T$ and $\sigma$:

```cpp
--8<-- "snippets/heuristics-reflect.cpp"
```

![Figure 5.7. The reflection about the main diagonal (the shaded cells)
on the 15-puzzle with the goal of Korf (right): the state $s$
(left), two moves from the goal, has the tiles 1 and 5 out of place;
its reflection $\phi(s)$ (middle) holds $\sigma(1) = 4$ in the cell
transposed of that of 1 and $\sigma(5) = 5$ in the cell transposed of
that of 5, and it is also two moves from the
goal.](../figures/heuristics-reflection.svg){#fig-5-7}

The *dual lookup* of Zahavi et al. [[29]](Z-bibliography.md#ref-29),
which looks up the placement of the inverse permutation of the state,
is another admissible lookup of the same tables; it is not consistent,
and it pays off with the bidirectional pathmax of IDA$^*$, which the
library does not have.

## 5.9 The encodings of the entries {#sec-5-9}

A table of the partition 7-8 of the 15-puzzle has more than half a
billion entries ([§5.11](05-heuristics.md#sec-5-11)), and the memory of
the entries decides both whether the tables fit and how often a lookup
misses the cache. Let $e_\pi(p) = (d_\pi(p) - m_\pi(p)) / 2$ be the
*excess* of the distance over the Manhattan distance of the tiles of the
pattern; by [Lemma 5.7](05-heuristics.md#lem-5-7) it is a nonnegative
integer, typically much smaller than the distance (cf. [Figure
5.8](05-heuristics.md#fig-5-8)). At the lookup, the Manhattan distance
of the tiles is computed from a table of $k n$ bytes, and the entries
store the excess, or less. We list the encodings of the library below
([Figure 5.8](05-heuristics.md#fig-5-8)).

* `ByteEntries` (`PatternDatabase.h:513`) stores $d_\pi$ in 8 bits per
  entry.
* `NibbleEntries` (`PatternDatabase.h:578`, the default) stores
  $e_\pi$ in 4 bits per entry, saturated at 15; the value is
  $m_\pi + 2 e_\pi$. A saturated entry gives a smaller value, which is
  still admissible, and `clamped()` counts them (the test `pdb` checks
  that the tables of the 8-puzzle have none).
* `Mod3Entries` (`PatternDatabase.h:666`) stores $d_\pi \bmod 3$, five
  residues in a byte as the digits of a number in base 3
  ($3^5 = 243 \le 256$), i.e., $8 / 5 = 1.6$ bits per entry, after
  Breyer and Korf [[25]](Z-bibliography.md#ref-25).
* `ParityBitEntries` (`PatternDatabase.h:745`) stores
  $\lfloor d_\pi / 2 \rfloor \bmod 2$, one bit per entry: with the
  parity of $d_\pi$, which is that of $m_\pi$, it gives
  $d_\pi \bmod 4$, as in the 1-bit encoding of Clausecker and
  Reinefeld [[27]](Z-bibliography.md#ref-27).
* `MinNibbleEntries` (`PatternDatabase.h:823`) stores, for each group
  of $c$ consecutive ranks ($c$ a power of 2, `PDBOptions::compress`),
  the least excess of the group in 4 bits, i.e., $4 / c$ bits per
  entry: the lossy max-compression of Felner et al.
  [[24]](Z-bibliography.md#ref-24), which here keeps the Manhattan
  distance of each placement exact.
* `SparseEntries` (`PatternDatabase.h:925`) stores a bit per entry that
  records whether its excess is nonzero, with a directory of rank and
  select of about $0.4$ bits per entry, and the nonzero excesses in 4
  bits each, in the order of their ranks; it is smaller than
  `NibbleEntries` when fewer than about 64% of the excesses are
  nonzero, at the price of a rank per lookup.

![Figure 5.8. The encodings of the entries of ranks 0 to 11 of the
table of [Figure 5.5](05-heuristics.md#fig-5-5): the distance $d$, the
Manhattan distance $m$ of the tiles of the pattern and the excess $e$;
the entries stored in a byte, in a nibble, modulo 3 and in one bit; the
least excess of each group of $c = 2$ ranks and the value it gives, in
rust where it is smaller than $d$ (the rank 5, of distance 4, gets
2).](../figures/heuristics-encodings.svg){#fig-5-8}

The values of `ByteEntries`, of `NibbleEntries` when no entry saturates,
and of `SparseEntries` (which saturates in the same way) when no entry
saturates, are the distances themselves, so that they are admissible and
consistent as in [§5.6](05-heuristics.md#sec-5-6) and
[§5.7](05-heuristics.md#sec-5-7). Each of the other three needs an
argument.

### The relative encodings

In `Mod3Entries` and `ParityBitEntries`, the residues do not give the
value of a placement, but they give the change of the value along a
move.

<a id="lem-5-10"></a>**Lemma 5.10** (The change from the residues).
With the blank ignored, let $p$ and $p'$ be adjacent placements and
$r$, $r'$ the residues of their distances modulo $M \in \{3, 4\}$.
Then

$$
d_\pi(p') - d_\pi(p) = \bigl( (r' - r + 1) \bmod M \bigr) - 1 .
$$

*Proof.* By [Lemma 5.7](05-heuristics.md#lem-5-7) the difference
$\delta$ is $+1$ or $-1$. If $\delta = 1$, then $r' - r \equiv 1$ and
the right-hand side is $2 - 1 = 1$; if $\delta = -1$, then
$r' - r \equiv -1$, $(r' - r + 1) \bmod M = 0$, and the right-hand side
is $-1$. Both cases are told apart since $1 \not\equiv -1 \pmod M$ for
$M \ge 3$. $\square$

A move of the puzzle moves the tile of one pattern only, so that the
heuristic of a successor is that of its parent plus the change of one
table (`AdditivePDB::delta()`, and `SlidingTilePDB` calls it for each
successor):

```cpp
--8<-- "snippets/heuristics-delta.cpp"
```

and A$^*$ never needs the absolute value but for the initial state.
That one is found from the residues only, by a descent to the goal
([Figure 5.9](05-heuristics.md#fig-5-9)).

<a id="thm-5-11"></a>**Theorem 5.11** (The value from the residues).
With the blank ignored, starting from a placement $p$ and moving each
time to a neighbour whose residue is that of the current placement
minus one modulo $M$, one reaches the goal placement in exactly
$d_\pi(p)$ steps; such a neighbour exists at each step before the
goal.

*Proof.* Let $q$ be the current placement, not the goal. The abstract
space is undirected with unit costs, and the next placement on a path
of least cost from $q$ to the goal is a neighbour of distance
$d_\pi(q) - 1$, whose residue is $r - 1 \bmod M$. Each neighbour has
distance $d_\pi(q) \pm 1$ ([Lemma 5.7](05-heuristics.md#lem-5-7)), and
those of distance $d_\pi(q) + 1$ have residue $r + 1 \bmod M \ne
r - 1 \bmod M$; hence each neighbour with the residue sought is one
step closer, and the descent reaches the goal placement, the only one
of distance 0, after exactly $d_\pi(p)$ steps. $\square$

Each step tries at most $4 k$ neighbours, i.e., $4 k$ ranks and lookups;
the value of the initial state thus costs $O(k\, d_\pi(p))$ lookups per
table (`PatternTable::exact()`):

```cpp
--8<-- "snippets/heuristics-exact.cpp"
```

We remark that the argument needs adjacent placements to differ by
exactly one, i.e., the blank ignored (the constructor of `PatternTable`
refuses a relative encoding with the blank minimized,
`PatternDatabase.h:1096`), and the heuristic to be the sum of the
tables, since from the maximum of two sums the change of one table
cannot be told: this is why `SlidingTilePDB` refuses a relative encoding
with the reflection (`PatternDatabase.h:1943`).

![Figure 5.9. The value of a placement from its residues modulo 3,
on the table of [Figure 5.5](05-heuristics.md#fig-5-5). (a) From
$p = (2, 0)$, each step goes to a neighbour whose residue is one less,
and the goal placement is reached in 5 steps, i.e., $d = 5$. (b) The
neighbours of $(4, 0)$, of residue 0 and distance 3: four of them have
residue 1, i.e., distance 4, and one has residue 2, i.e., distance 2,
which is the one taken.](../figures/heuristics-descent.svg){#fig-5-9}

### The max-compression

In `MinNibbleEntries`, the value is $m_\pi(p) + 2 \min_{p' \in G(p)}
e_\pi(p') \le m_\pi(p) + 2 e_\pi(p) = d_\pi(p)$, with $G(p)$ the group
of the rank of $p$; hence each table is a lower bound on its distance,
and the projection argument of [Theorem 5.6](05-heuristics.md#thm-5-6)
holds for the sum. It is not consistent in general, since two adjacent
placements may fall in groups with very different minima, and A$^*$ may
then reopen nodes. Felner et al. [[24]](Z-bibliography.md#ref-24) store
the least distance of the group; storing the least excess loses less,
because the placements of a group differ mostly in the cell of the last
tile (the groups are aligned to $c$ and the blocks of equal prefix to $n
- k + 1$, so that a group with $c \le n - k + 1$ meets at most two
blocks), and the Manhattan distance, which accounts for most of their
differences, is kept exact.

## 5.10 Files and the cache of the tables {#sec-5-10}

Building the largest tables (the 24-puzzle, the partition 7-8 of the
15-puzzle) is far more expensive than a search on them; hence the tables
are built once and saved. `PatternTable::save()`
(`PatternDatabase.h:1272`) writes a table with a versioned header (the
magic, the version of the format, the board, the rule of the blank, the
encoding and its parameter, the tiles and the goal), the payload and a
checksum of everything, to a temporary file in the same directory, which
is then renamed; a reader never sees it half written. Conversely,
`load()` rejects a file that is truncated, corrupted, of another version
or of another table. With a cache directory (`PDBOptions::cache`),
`PatternTable::cached()` (`PatternDatabase.h:1362`) loads a table if it
is there; otherwise it encodes the table from the table in bytes of the
same pattern, if that one is there, or it builds it, and the tables
built are saved. Also, the test `pdb` checks the round trip, the
rejection of a truncated, a corrupted and a foreign file, and the reuse
of the cache.

## 5.11 Sizes {#sec-5-11}

The table of a pattern of $k$ tiles has $N = n! / (n - k)!$ entries. The
partitions of the library (7-8 and 6-6-6-6 of
[[22]](Z-bibliography.md#ref-22), 6-6-3 and 5-5-5 of
[[23]](Z-bibliography.md#ref-23), on the 15-puzzle or the 24-puzzle as
the number in parentheses indicates) give the following sizes, computed
from $N$ and from the bits per entry of each encoding. In the last
column we give the max-compression by $c = 8$ for 7-8 and by $c = 16$
for 5-5-5, the groups of the benchmark of
[§5.13](05-heuristics.md#sec-5-13). We remark that the entries are
summed over the patterns, e.g., $57657600 + 518918400$ for the patterns
of 7 and 8 tiles of the partition 7-8.

| partition | entries | 8 bits | 4 bits | mod 3 | 1 bit | min |
|---|---|---|---|---|---|---|
| 7-8 (15) | 576576000 | 576 MB | 288 MB | 115 MB | 72 MB | 36 MB |
| 6-6-3 (15) | 11534880 | 11.5 MB | 5.8 MB | 2.3 MB | 1.4 MB | |
| 5-5-5 (15) | 1572480 | 1.6 MB | 786 KB | 314 KB | 197 KB | 49 KB |
| 6-6-6-6 (24) | 510048000 | 510 MB | 255 MB | 102 MB | 64 MB | |

A table with the blank minimized has the same size (cf.
[§5.5](05-heuristics.md#sec-5-5)), but its construction needs $N (n - k)
/ 4$ more bytes for the states of Algorithm 5.2, i.e., about 1 GB for
the pattern of 8 tiles of the 15-puzzle.

## 5.12 When the heuristic does not help {#sec-5-12}

Recall that the time of a search is the number of its expansions times
the time of an expansion ([§2.10](02-best-first.md#sec-2-10)). A search
that is too slow is thus too slow for one of these two reasons, and the
counters of the `SearchStats` usually indicate which one, and what one may try.

* *Reexpansions.* If `reexpanded` is positive, the heuristic is not
  consistent ([Theorem 2.13](02-best-first.md#thm-2-13)); a heuristic
  meant to be consistent then has an error, and an inconsistent one
  (e.g., the blank minimized) may cost more in reexpansions than it
  saves in expansions. On a domain whose states can all be enumerated,
  the tests of the library check the admissibility and the consistency
  of a heuristic on each state and each move (the 8-puzzle has $9!/2 =
  181440$ states, `check_heuristic8` in `test/test.cpp`). The same check
  on a small instance of a new domain is the quickest diagnosis.
* *Expansions.* The number of expansions is fixed by the heuristic up to
  the ties ([Theorem 2.14](02-best-first.md#thm-2-14)); a stronger
  admissible heuristic (larger patterns, the reflection, the blank
  minimized) reduces it, and the mean heuristic over a sample of states,
  which `PDBBench --mode table` reports, is a first estimate of the
  gain. If even a heuristic with a large weight does not solve the
  instance (a weighted A$^*$, which the library obtains by a wrapper of
  the domain, [§2.2](02-best-first.md#sec-2-2)), the instance is
  probably very hard for any optimal search. If an instance seems never
  to be solved, one should check that it is solvable at all:
  `SlidingTile` refuses a start whose permutation has the wrong parity
  (`SlidingTile.h:212`), while on a grid A$^*$ proves the absence of a
  path only after exhausting the component of the start.
* *Time per expansion.* The rate of expansions (`expanded / seconds`)
  falls when the heuristic is costly: a lookup of a table much larger
  than the cache misses on most expansions, and a smaller encoding
  ([§5.9](05-heuristics.md#sec-5-9)) may then be faster even with more
  expansions. The data structures of
  [Chapter 3](03-data-structures.md#ch-3) are the other cost.
* *Memory.* If the search stops with `Status::node_limit`, or the memory
  runs out, `stored` measures the nodes held. Furthermore, the searches
  of [Chapter 4](04-memory-bounded.md#ch-4) trade time for memory, and
  IDA$^*$ with a pattern database is the method of
  [[22]](Z-bibliography.md#ref-22) for the 24-puzzle (the test `tiles24`
  runs it with patterns of 3 tiles).

## 5.13 Benchmark {#sec-5-13}

Our benchmark of the heuristics is the sequential A$^*$ on the 100
instances of the 15-puzzle of Korf [[7]](Z-bibliography.md#ref-7), with
the Manhattan distance and with the pattern databases, run by
`tools/PDBBench` from `tools/bench/pdb.sh`. First, the script builds (or
loads from its cache directory) the tables of the partitions 5-5-5,
6-6-3 and 7-8, with the blank ignored and, for 5-5-5 and 7-8, minimized;
the searches use the encodings in a nibble, in one bit (without the
reflection), max-compressed (by 16 for 5-5-5, by 8 for 7-8), and by rank
and select when it is available. Then it solves each instance with each
heuristic, twice, running up to `JOBS` searches at the same time (4 by
default). The costs are integers and exact, i.e., each heuristic is
admissible and each instance is solved to its optimal cost. The machine
and the compiler are those stated in the note of the table. The clock is
the `Stopwatch` of `AStar::solve()`, as in
[§2.11](02-best-first.md#sec-2-11): it does not include the construction
or the loading of the tables, which precede the search, nor the start of
the process and the reading of the instances.

**Table 5.1.** A$^*$ with the Manhattan distance and with the pattern
databases on the 100 instances of Korf: for each heuristic (the
partition and its encoding), the instances solved, the nodes expanded
summed over the instances, and the seconds summed over the instances,
each instance the median of its repetitions.

--8<-- "tables/heuristics-pdb.md"

The numbers of expansions do not depend on the machine, and those of
the README of the repository (section "Pattern databases", from the
same script) are the following, summed over the 100 instances, with
the ratio to the Manhattan distance.

| heuristic | blank | entries | memory | expanded | ratio |
|---|---|---|---|---|---|
| Manhattan | | | | 1557138490 | 1 |
| 5-5-5, max-compressed by 16 | ignored | min | 49 KB | 681766045 | 2.3 |
| 5-5-5, no reflection | ignored | 1 bit | 197 KB | 369972440 | 4.2 |
| 5-5-5 | ignored | 4 bits | 786 KB | 263234565 | 5.9 |
| 6-6-3 | ignored | 4 bits | 5.8 MB | 167667391 | 9.3 |
| 7-8, no reflection | ignored | 1 bit | 72 MB | 137386424 | 11.3 |
| 7-8, max-compressed by 8 | ignored | min | 36 MB | 119574805 | 13.0 |
| 7-8 | ignored | 4 bits | 288 MB | 43587167 | 35.7 |
| 5-5-5 | minimized | 4 bits | 786 KB | 59968561 | 26.0 |
| 7-8 | minimized | 4 bits | 288 MB | 2019506 | 771 |

The rule of the blank weighs more than the size of the patterns there:
with the blank minimized the partition 7-8 expands about 20000 nodes
per instance, 21 times fewer than with the blank ignored, at the price
of about 6000 reexpansions over the 100 instances. On the 24-puzzle the
same script solves with IDA$^*$ and the partition 6-6-6-6 the ends of
50 random walks of 150 moves of the blank from the goal, with at most
$10^{10}$ expansions each; the README reports 49 of them solved, with
optimal costs from 50 to 84, while the 50 instances used by Jinnai
and Fukunaga [[5]](Z-bibliography.md#ref-5) are beyond such a limit.

## 5.14 Notes and references {#sec-5-14}

Pattern databases were introduced by Culberson and Schaeffer
[[21]](Z-bibliography.md#ref-21) for the 15-puzzle, with patterns that
count each move. Korf and Felner [[22]](Z-bibliography.md#ref-22) made
them additive by counting only the moves of the tiles of the pattern,
with the partitions 7-8 and 6-6-6-6 and the reflection, and Felner, Korf
and Hanan [[23]](Z-bibliography.md#ref-23) studied the static and
dynamic partitions, among them 5-5-5 and 6-6-3. Yang et al.
[[28]](Z-bibliography.md#ref-28) gave the general theory of the additive
abstractions, and Holte et al. [[30]](Z-bibliography.md#ref-30) showed
that the maximum over several smaller databases can beat one larger
database. Zahavi et al. [[29]](Z-bibliography.md#ref-29) introduced the
dual lookups, which Felner et al. [[31]](Z-bibliography.md#ref-31)
exploit with the bidirectional pathmax. Felner et al.
[[24]](Z-bibliography.md#ref-24) compressed the tables by the maximum
(lossy) and without loss when the merged entries are known to differ by
little. Breyer and Korf [[25]](Z-bibliography.md#ref-25) stored the
distances modulo 3, and Clausecker and Reinefeld
[[27]](Z-bibliography.md#ref-27) in one bit with the zero-aware
encoding. Samadi et al. [[33]](Z-bibliography.md#ref-33) learned a
smaller model of the entries, and Nagahashi and Takahashi
[[32]](Z-bibliography.md#ref-32) built the largest tables on clusters,
by partitioning the ranks among the processes as the blocks of the
library are partitioned among the threads. The ranking of Myrvold and
Ruskey [[34]](Z-bibliography.md#ref-34) is the linear-time alternative
to the lexicographic one, and the two-bit breadth-first search of Korf
[[26]](Z-bibliography.md#ref-26) is the model of the construction with
the blank minimized. The library differs from these works in storing the
excess over the Manhattan distance rather than the distance. According
to the documentation of the repository, this makes 4 bits per entry
enough for the tables of the library, and it keeps the Manhattan
distance exact in the max-compression. Yet, the library has
neither the dual lookups nor the lossless compression, and the searches
that use the tables best on the 24-puzzle, i.e., IDA$^*$ with the
bidirectional pathmax [[31]](Z-bibliography.md#ref-31), are beyond its
IDA$^*$, which has no pathmax ([Chapter 4](04-memory-bounded.md#ch-4)).

## 5.15 Exercises {#sec-5-15}

**Exercise 5.1** (`ctest -R PBeFS.grid`). On a grid with 8 moves, where
a move pays the cost of the cell it leaves and $\sqrt{2}$ times that on
a diagonal, prove that the octile distance times the least cost of a
cell is consistent, that the Euclidean distance times the same cost is
consistent and dominated by it, and that the Manhattan distance is not
admissible ([Figure 5.1](05-heuristics.md#fig-5-1)). The test checks
A$^*$ with the octile distance against the costs of Dijkstra.

**Exercise 5.2** (`ctest -R PBeFS.pdb`). Compute by hand the rank of
the placement $(2, 7, 1)$ of three tiles on 9 cells, and unrank 229 back
to $(4, 0, 7)$ ([Figure 5.4](05-heuristics.md#fig-5-4)). Show that the
ranks $r$ and $r + 1$ in the same block of $n - k + 1$ consecutive
ranks differ only in the cell of the last tile.

**Exercise 5.3** (`ctest -R PBeFS.pdb`). Prove that with $k = n - 1$
tiles (all the tiles in the pattern) both rules of the blank give the
true distance of the puzzle, and give a placement of a pattern of 2
tiles of the 8-puzzle for which the two rules give different
distances. The test checks
[Lemma 5.3](05-heuristics.md#lem-5-3) on each state of the 8-puzzle.

**Exercise 5.4** (`ctest -R PBeFS.pdb`). Find a move of the 8-puzzle
along which the sum of the tables of the partition 4-4 with the blank
minimized decreases by more than one, and explain it by the regions of
the blank. Why does the proof of
[Theorem 5.8](05-heuristics.md#thm-5-8) not apply? The test counts these
moves.

**Exercise 5.5** (`ctest -R PBeFS.pdb`). Prove that the reflection of
[Lemma 5.9](05-heuristics.md#lem-5-9) is not an automorphism when the
goal has the blank off the diagonal, by an example on the 8-puzzle, and
check that the library then does not use it. Show that
$\phi(\phi(s)) = s$.

**Exercise 5.6** (`ctest -R PBeFS.pdb`). Prove
[Lemma 5.10](05-heuristics.md#lem-5-10) for $M = 4$ with the bit of
`ParityBitEntries` and the parity of $m_\pi$, and show by an example
that the residues modulo 2 would not suffice. Why does a relative
encoding with the reflection fail? The test checks that the value of
each state of the 8-puzzle is recovered from the residues only, and
that of each successor from its parent.

**Exercise 5.7** (`ctest -R PBeFS.pdb`). Show by an example that the sum
of the tables of two patterns that share a tile, or of a pattern that
holds the blank, may exceed $h^*$; the test checks that the library
refuses both.

**Exercise 5.8** (`ctest -R PBeFS.pdb-search`). Prove that the
Manhattan distance is the additive database of the patterns of one
tile with the blank ignored, and that merging two disjoint patterns
never decreases the heuristic. By
[Theorem 2.15](02-best-first.md#thm-2-15), what can be said of the
expansions of A$^*$ with the partition 4-4 of the 8-puzzle against
those with the Manhattan distance? The test checks that they are fewer.

**Exercise 5.9** (`ctest -R PBeFS.tiles24`). The test builds the
additive database of the 24-puzzle with the 8 patterns of 3 consecutive
tiles, checks that it dominates the Manhattan distance on random walks
of 40 moves, and solves them with IDA$^*$. Compute the number of the
entries of these tables, and compare it with that of the partition
6-6-6-6 ([§5.11](05-heuristics.md#sec-5-11)).

**Exercise 5.10** (`ctest -R PBeFS.pdb`). Prove that the value of
`MinNibbleEntries` is admissible, and give two adjacent placements of
the table of [Figure 5.5](05-heuristics.md#fig-5-5) whose values with
$c = 2$ differ by 3. The test checks that the max-compressed tables
never exceed the exact ones, that `SparseEntries` gives the same values
as `NibbleEntries`, and that a corrupted file is rejected.
