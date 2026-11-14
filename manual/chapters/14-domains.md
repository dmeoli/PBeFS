# 14. The example domains {#ch-14}

[Source: `include/PBeFS/astar/SlidingTile.h`, `SequenceScore.h`,
`Grid.h` (with `src/Grid.cpp`), `ExplicitGraph.h`, `Romania.h`,
`NQueens.h`, `VacuumWorld.h`, `WaterJugs.h`, `MissionariesCannibals.h`,
`TravelingSalesman.h`, `SlowDomain.h`; `examples/*.cpp`,
`tools/AStarBench.cpp`, `data/`]

This chapter describes the problems that the library ships as models of
the concept `Domain` of [Chapter 1](01-introduction.md#ch-1), on which
the searches of the other chapters are run, tested and measured. For
each one we give the state and its encoding (with the trade-offs we made
in it), the operators and their costs, the heuristic with the proof of
its admissibility or consistency (or a counterexample where it fails),
the refinements of `Domain` it models, its instances, and how the
examples and the benchmark `AStarBench` run it. The sliding-tile puzzles
([§14.2](#sec-14-2) and [§14.3](#sec-14-3)) and the grids
([§14.4](#sec-14-4)) are the usual benchmarks of the parallel searches
in the literature; the puzzles are the only domain that models all the
refinements, and the grids model all of them but the incremental
features and the operators. The explicit graphs ([§14.5](#sec-14-5)),
the problems of the textbook of Russell and Norvig ([§14.6](#sec-14-6),
cf. chapter 3 of [[13]](Z-bibliography.md#ref-13)) and the traveling
salesperson problem ([§14.7](#sec-14-7)) cover the other shapes of a
state space that a search may (and usually does) meet (a tree with dead
ends, an irreversible problem, real costs, a heuristic that is
admissible but not consistent), and `SlowDomain` ([§14.8](#sec-14-8))
makes any of them slow to expand. We leave the knapsack of the
branch-and-bound searches to [Chapter 13](13-branch-and-bound.md#ch-13),
and the multiple sequence alignment in [Chapter
15](15-msa-instances.md#ch-15).

## 14.1 The domains and their concepts {#sec-14-1}

Table 14.1 gives, for each domain, the refinements of `Domain` it models;
they decide which searches and which work distributions of HDA$^*$ can be
run on it ([Figure 1.5](01-introduction.md#fig-1-5)). A domain built on
another one inherits its refinements where the code forwards them
(`SequenceScore` forwards none, `SlowDomain` all but
`IncrementalFeatureDomain`), and the last column is the concept
`DiskDomain` of the search with delayed duplicate detection of
[Chapter 4](04-memory-bounded.md#ch-4), i.e., integer costs and trivially
copyable states.

**Table 14.1.** The refinements of `Domain` modeled by each domain: F
feature, IF incremental feature, AF abstract feature, AS abstract state,
NB nblock, PW pairwise, OP operator, PH perfect hash, DK disk (x
modeled, (a) only up to 20 cells, (b) only with 4 moves and with integer
costs; the cost type of the grids is `int` with 4 moves and `double`
with 8).

| domain | cost | F | IF | AF | AS | NB | PW | OP | PH | DK |
|---|---|---|---|---|---|---|---|---|---|---|
| `SlidingTile< W , H >` | `int` | x | x | x | x | x | x | x | (a) | x |
| `SequenceScore< W , H >` | `int` | | | | | | | | | x |
| `Grid< 4 >`, `Grid< 8 >` | both | x | | x | x | x | x | | x | (b) |
| `ExplicitGraph< Cost >` | `Cost` | x | | | | | | | x | (b) |
| `Romania` | `double` | x | | | | | | | x | |
| `NQueens` | `int` | x | | | | | | | | x |
| `VacuumWorld` | `int` | x | | x | x | x | x | | x | x |
| `WaterJugs` | `int` | x | | | | | | | x | x |
| `MissionariesCannibals` | `int` | x | | | | | | | x | x |
| `TravelingSalesman` | `double` | x | | | | | | | x | |

Every domain packs its state into an integer, or into a small trivially
copyable record (the board of `NQueens`). A state is thus compared and
copied as a word, and its hash is the finalizer of SplitMix64 applied to
it; this is the representation that Burns et al.
[[4]](Z-bibliography.md#ref-4) recommend, since the state is then at the
same time its key and the input of its hash.

## 14.2 The sliding-tile puzzles {#sec-14-2}

The sliding-tile puzzle on a board of $W \times H$ cells has $WH - 1$
tiles, numbered from 1, and one empty cell, the *blank*; a move slides a
tile adjacent to the blank into it, at cost 1, and the goal is a given
arrangement of the tiles. In the library, the class template
`SlidingTile< W , H >` (`SlidingTile.h:170`) models it for $WH \le 25$,
with the aliases `EightPuzzle`, `FifteenPuzzle` and `TwentyFourPuzzle`
for the boards of $3 \times 3$, $4 \times 4$ and $5 \times 5$ cells; the
goal is by default that of Korf [[7]](Z-bibliography.md#ref-7), with the
blank in the first cell and the tiles $1, 2, \ldots$ in row-major order.

### The state and the moves

A state is packed into one integer, in which the field of $b$ bits of
the cell $i = rW + k$ holds the tile in row $r$ and column $k$ (0 for
the blank). With $b = 4$ the boards of up to 16 cells fit a word of 64 bits,
and with $b = 5$ the larger ones fit a `Word128`, an unsigned integer of
128 bits written on two words of 64 in order to be the same with all
compilers (`SlidingTile.h:71`). This integer is at the same time the
state, its equality test and the key of its hash. A move exchanges the
blank, in the cell $b$, with the tile $t$ in a cell $c$ adjacent to it:
the field of $c$ is cleared and $t$ is written in the field of $b$, two
operations on the word ([Figure 14.1](#fig-14-1)).

![Figure 14.1. A move of the 8-puzzle, on the instance of figure 3.28
of [[13]](Z-bibliography.md#ref-13) ($u$, whose optimal solution has 26
moves, towards the goal on the right). Sliding tile 2 down into the blank
moves it from distance 1 to distance 2 from its goal cell, and the
Manhattan distance of the successor follows from that of $u$ by the change
of that tile alone.](../figures/domains-tiles-move.svg){#fig-14-1}

The cell of the blank is not stored: it is found from the word by a few
bitwise operations ([Figure 14.2](#fig-14-2)). In fact, the bits of each
field are or-ed into its lowest bit, $s \mid (s \gg 1) \mid \cdots \mid
(s \gg (b - 1))$. Hence, the lowest bit of a field is 0 if and only if
the whole field is 0, i.e., in the field of the blank; the complement of
the result, masked by the lowest bit of each field, has a single bit
set, and the number of its trailing zeros divided by $b$ is the cell of
the blank. Note that the shifts let the bits of the next field into the
top of a field, but never into its lowest bit, which is the only one
that is read.

```cpp
--8<-- "snippets/domains-tiles-blank.cpp"
```

![Figure 14.2. The cell of the blank of the instance of
[Figure 14.1](#fig-14-1), from its packed word (four bits per cell): the
lowest bit of $s \mid s{\gg}1 \mid s{\gg}2 \mid s{\gg}3$ is 0 only in the
field of the blank, and the trailing zeros of the complement (masked by
the lowest bits) give the cell 4.
](../figures/domains-tiles-packing.svg){#fig-14-2}

The moves of each cell are precomputed at construction, in the order up,
left, right and down (`SlidingTile.h:220`), and the expansion reads them
from that table (at most 25 entries of 5 integers, hence usually in the
cache). Also, the domain is an `OperatorDomain` whose operator is the
cell of the blank in the parent: in the expansion with an operator, the
move that sends the blank back to that cell is skipped ([Lemma
1.5](01-introduction.md#lem-1-5)).

```cpp
--8<-- "snippets/domains-tiles-successors.cpp"
```

### The reachable states

Only half of the arrangements of the tiles can be reached from a given
one, and the constructor refuses a start and a goal that are not in the
same half (`SlidingTile.h:212`). We now state the invariant, after
fixing the notation. For a board $\sigma$, let $N(\sigma)$ be the number
of *inversions* of the tiles read in row-major order without the blank
(the pairs of tiles $t > t'$ with $t$ before $t'$), and $r(\sigma)$ the
row of the blank.

**Proposition 14.1** (Parity of the sliding-tile puzzles).
<a id="prop-14-1"></a> If $W$ is odd, $N(\sigma) \bmod 2$ is invariant
under the moves; if $W$ is even, $(N(\sigma) + r(\sigma)) \bmod 2$ is.
Hence, at most half of the $(WH)!$ boards are reachable from a given one.

*Proof.* A horizontal move does not change the order of the tiles read
in row-major order without the blank, nor the row of the blank. A
vertical move takes a tile $t$ over the $W - 1$ tiles that lie between
its cell and that of the blank in row-major order, and changes the row
of the blank by 1. Each of these $W - 1$ tiles forms an inversion with
$t$ before the move but not after it, or conversely, and therefore $N$
changes by $W - 1 - 2j$ for some $j$, i.e., by an even number if $W$ is
odd and by an odd one if $W$ is even. In the first case $N \bmod 2$ is
preserved; in the second, $N$ and $r$ change parity together, and their
sum keeps it. Finally, the boards on which the invariant takes the other
value are half of the $(WH)!$ boards. In fact, exchanging two tiles
changes $N$ by an odd number and leaves $r$ unchanged (and an exchange
is a bijection between the two halves). $\square$

The converse, that every board on which the invariant takes the value of
the start is reachable, is the theorem of Johnson and Story
[[150]](Z-bibliography.md#ref-150) for the 15-puzzle, and it holds for all
the boards of at least $2 \times 2$ cells; hence, the 8-puzzle has $9!/2
= 181440$ reachable states, which the test `astar` checks by a
breadth-first search of the whole space, the 15-puzzle $16!/2 \approx
1.05 \cdot 10^{13}$ and the 24-puzzle $25!/2 \approx 7.8 \cdot 10^{24}$.
In the excerpt below, the function `parity` is the invariant of
[Proposition 14.1](#prop-14-1), the row of the blank being added on the
boards of even width.

```cpp
--8<-- "snippets/domains-tiles-parity.cpp"
```

On the puzzles, the branching factor of a state is the number of the
neighbours of its blank, 2 in a corner, 3 on a side and 4 inside the
board, and its mean over the reachable states has a closed form ([Figure
14.3](#fig-14-3)).

**Proposition 14.2** (Branching factor of the puzzles).
<a id="prop-14-2"></a> Over the reachable states of the puzzle of $W
\times H$ cells, the mean number of moves is $2(2WH - W - H)/(WH)$, i.e.,
$8/3$ for the 8-puzzle, 3 for the 15-puzzle and $16/5$ for the 24-puzzle;
without the move back to the parent it is 1 less, but at $s_0$.

*Proof.* Fix a cell of the blank. The reachable boards with the blank in
that cell are as many as the arrangements of the $WH - 1$ tiles on the
other cells that satisfy the invariant of [Proposition
14.1](#prop-14-1), i.e., $(WH - 1)!/2$ (by the theorem of Johnson and
Story and by the argument of the exchange of two tiles). Hence, the
blank is uniformly distributed over the cells, and the mean number of
moves is the mean degree of a cell in the grid graph of the board,
$2|E|/(WH)$, with $|E| = H(W - 1) + W(H - 1) = 2WH - W - H$ edges. The
three values follow with $|E| = 12, 24, 40$ (cf. [Figure
14.3](#fig-14-3)). Skipping the parent removes one successor of each
state but $s_0$. $\square$

![Figure 14.3. The number of moves of the blank in each cell of the
boards of $3 \times 3$ and $4 \times 4$ cells. Since the blank is
uniformly distributed over the cells of the reachable states, the mean
branching factor is the mean of these numbers
([Proposition 14.2](#prop-14-2)).
](../figures/domains-tiles-blank.svg){#fig-14-3}

### The Manhattan distance

The heuristic is the Manhattan distance: the sum, over the tiles, of the
distance between the cell of the tile and its goal cell in moves of the
grid, $h(\sigma) = \sum_{t \ne 0} \delta_t(c_\sigma(t))$, with
$\delta_t(c)$ the Manhattan distance of the cell $c$ from the goal cell
of $t$; the table `v_md[ t ][ c ]` holds $\delta_t(c)$ for each tile and
cell (`SlidingTile.h:217`).

**Proposition 14.3** (The Manhattan distance is consistent).
<a id="prop-14-3"></a> On the sliding-tile puzzle, for each move from
$\sigma$ to $\sigma'$ of a tile $t$ from the cell $c$ to the cell $b$ of
the blank,
$h(\sigma') = h(\sigma) + \delta_t(b) - \delta_t(c)$ and
$|h(\sigma') - h(\sigma)| = 1$; hence, $h$ is consistent, and therefore
admissible.

*Proof.* Only the tile $t$ changes its cell, which gives the first
identity. The cells $b$ and $c$ are adjacent. Therefore, their distances
from the goal cell of $t$ differ by 1 (one coordinate changes by 1, and
the parity of the distance changes). Hence $h(\sigma) \le 1 + h(\sigma')
= c(\sigma, \sigma') + h(\sigma')$, and $h$ is 0 at the goal;
consistency implies admissibility ([Chapter 2](02-best-first.md#ch-2)).
$\square$

The first identity is how the expansion computes $h$: it receives
$h(\sigma)$ and adds `v_md[ t ][ b ] - v_md[ t ][ c ]`, two lookups
instead of a sum over the tiles, which is the incremental heuristic of
[[4]](Z-bibliography.md#ref-4). Furthermore, the second identity has a
consequence on the shape of the search: $f$ changes by $0$ or $2$ at
each move, and hence on the puzzles the values of $f$ met by A$^*$ all
have the parity of $h(s_0)$, and iterative deepening raises its bound by
2 at each iteration ([Chapter 4](04-memory-bounded.md#ch-4)). Of course,
the pattern databases of [Chapter 5](05-heuristics.md#ch-5) replace this
heuristic in `SlidingTilePDB`, which keeps the rest of the domain.

### The refinements

The puzzles model all the refinements of `Domain` ([Figure
14.4](#fig-14-4)). Its feature $i = tWH + c$ is "tile $t$ is in cell
$c$", for $t \ne 0$. A move changes two features, and
`for_each_feature_change` finds them without knowing the moved tile: the
fields of $s \oplus s'$ are 0 except in the two cells $b$ and $c$, where
they both hold $t$: the lowest set bit gives the first cell and the
tile, and the lowest set bit of the rest the second cell.

```cpp
--8<-- "snippets/domains-tiles-feature-change.cpp"
```

Next, the abstract feature of AZH keeps of the cell of a tile only
whether it lies in the upper or in the lower half of the board, which is
the projection of figure 8(b) of Jinnai and Fukunaga
[[5]](Z-bibliography.md#ref-5) for the 15-puzzle; a move changes an
abstract feature only when a tile moves between the two halves of the
board, and therefore
most successors (usually) stay on the worker of their parent. As in
Burns et al. [[3]](Z-bibliography.md#ref-3), the abstract state of
AHDA$^*$ is the triple of the cells of tiles 1, 2 and 3. The nblock of
PBNF is the triple of the cells of the blank and of tiles 1 and 2; the
abstract successors of an nblock are the triples reached by a move of
the blank, which may also move tile 1 or 2 (`SlidingTile.h:399`). Also,
the pairwise heuristic of PA$^*$SE is the Manhattan distance between two
boards, consistent for the same reason as [Proposition
14.3](#prop-14-3). Finally, the perfect hash is the lexicographic rank
of the permutation of the cells of the tiles. Burns et al.
[[3]](Z-bibliography.md#ref-3) used it as a work distribution modulo
$p$, and, as Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) show,
depends then on very few tiles; it is kept as the poor distribution to
compare against, and only up to 20 cells, since $21! > 2^{64}$.

![Figure 14.4. The abstractions of the 15-puzzle on the instance 1 of
Korf: (a) the abstract feature of AZH, a tile in the upper or in the lower
half of the board; (b) the abstract state of AHDA$^*$, the cells of tiles 1,
2 and 3; (c) the nblock of PBNF, the cells of the blank and of tiles 1 and
2.](../figures/domains-tiles-abstraction.svg){#fig-14-4}

### Instances

The file `data/15-puzzle/korf100.txt` holds the 100 instances of the
15-puzzle of Korf [[7]](Z-bibliography.md#ref-7), each with its optimal
cost, which the tests and the benchmarks use as the reference, and
`data/24-puzzle/instances50.txt` the 50 instances of the 24-puzzle
distributed with the code of HDA$^*$ of Jinnai and Fukunaga; the example
`8Puzzle` solves the instance of figure 3.28 of
[[13]](Z-bibliography.md#ref-13), or the first state found at a given
distance from the goal, and `15Puzzle` an instance of Korf, sequentially
or with HDA$^*$.

## 14.3 The sequence score of Nilsson {#sec-14-3}

`SequenceScore< W , H >` (`SequenceScore.h:54`) is the same puzzle with
the heuristic of Nilsson [[8]](Z-bibliography.md#ref-8), $h = M + W S$,
where $M$ is the Manhattan distance and $S$ the *sequence score*: on the
$3 \times 3$ board, $S$ counts 1 if the centre does not hold its goal
tile, and 2 for each tile of the border not followed, clockwise, by the
tile that follows it in the goal; on the other boards it counts 1 if the
first cell does not hold its goal tile, and 2 for each tile not followed,
in row-major order, by its successor in the goal.

```cpp
--8<-- "snippets/domains-sequence-heuristic.cpp"
```

**Proposition 14.4** (The sequence score is not admissible).
<a id="prop-14-4"></a> On the 8-puzzle with the goal of Korf, the state
one move away from the goal in which tile 1 has moved into the first cell
has $h = 13$, while $h^* = 1$.

*Proof.* Only tile 1 is out of place, at distance 1, and $M = 1$. Since
the centre holds tile 4, its goal tile, it contributes 0. Along the
border ([Figure 14.5](#fig-14-5)), tile 1 in the first cell is followed
by the blank instead of tile 2, and tile 3 in the cell 3 is followed (in
the first cell) by tile 1 instead of the blank, which follows it in the
goal; each other tile of the border is followed by its goal follower.
Hence $S = 4$ and $h = 1 + 3 \cdot 4 = 13 > 1$. $\square$

![Figure 14.5. The sequence score on the 8-puzzle. Left: the goal,
with the clockwise order of the border that defines the follower of each
tile. Right: the state one move away of
[Proposition 14.4](#prop-14-4), where the two pairs in rust break the
order.](../figures/domains-sequence.svg){#fig-14-5}

Since the heuristic may overestimate, A$^*$ with it is not guaranteed (cf.
[Chapter 2](02-best-first.md#ch-2)) to return an optimal path, and it is
not: on the instance of [Figure 14.1](#fig-14-1) (optimal cost 26, with
$M = 18$ and $h = 69$), the example `8Puzzle -h seq` returns a path of
28 moves after very few expansions. This domain is there to show this
trade, which weighted A$^*$ makes with a bound ([Chapter
6](06-bounded-suboptimal.md#ch-6)). It forwards no refinement of
`SlidingTile`, and its heuristic is computed from scratch at each
successor.

## 14.4 Pathfinding on grids {#sec-14-4}

We call *grid* a rectangular map of $W \times H$ cells, each "blocked"
or free with a cost in $\{1, \ldots, 255\}$ (`GridMap`, `Grid.h:73`), on
which the agent moves from cell to cell; it is the benchmark of the
video games and of robotics [[15]](Z-bibliography.md#ref-15), and many
problems reduce to it (the multiple sequence alignment of [Chapter
15](15-msa-instances.md#ch-15) is a shortest path on a lattice of higher
dimension). `Grid< Moves >` (`Grid.h:259`) models pathfinding with 4
moves (up, left, right, down) or 8 (the four diagonals too), generated
by a `GridAdjacency< Moves >` with the std and Boost data structures and
by the grid graph of the default ones. A state is
the index $yW + x$ of the cell $(x, y)$, which is also the perfect hash.
A move *leaves* a cell at its cost, times $\sqrt2$ for a diagonal, i.e.,
$c(u, v) = \kappa(u)$ or $\sqrt2\, \kappa(u)$, with $\kappa(u)$ the cost
of the cell $u$; and a diagonal move is allowed only if both the cells
it cuts are free ([Figure 14.6](#fig-14-6)). Hence, the costs are
integers with 4 moves and real numbers with 8.

![Figure 14.6. (a) The moves from a cell $u$ of cost 2 with the cell
above it blocked: the straight moves cost 2, the diagonal ones
$2\sqrt2$, and the two diagonals that would cut the corner of the blocked
cell are not moves. (b) The octile distance between two cells: as many
diagonal moves as the smaller difference of the coordinates, then
straight moves.](../figures/domains-grid-moves.svg){#fig-14-6}

Note that the cost of an edge depends on the cell left: $c(u, v) =
\kappa(u)$ and $c(v, u) = \kappa(v)$ differ on a map with several costs:
the graph is directed, and the problem is reversible (in the sense of
[§1.2](01-introduction.md#sec-1-2)) only on the maps of unit cost. This
is why a grid is not an `OperatorDomain`, while its moves are symmetric.

### The heuristics

The heuristic is the distance on the empty grid times the least cost
$\kappa_{\min}$ of a free cell: the Manhattan distance $d_4(u, t) =
\Delta x + \Delta y$ with 4 moves, and the octile distance

$$
d_8(u, t) = \max(\Delta x, \Delta y) + (\sqrt2 - 1) \min(\Delta x,
\Delta y)
$$

with 8 moves, $\Delta x$ and $\Delta y$ being the absolute differences of
the coordinates of $u$ and of the goal $t$.

```cpp
--8<-- "snippets/domains-grid-distance.cpp"
```

Our proof of consistency rests on $d_4$ and $d_8$ being the costs of the
shortest paths on the grid with no blocked cell and unit cell costs; for
$d_4$ this is immediate, and for $d_8$ it is the following lemma.

**Lemma 14.5** (The octile distance). <a id="lem-14-5"></a> On the
infinite grid with 8 moves, straight moves of cost 1 and diagonal moves
of cost $\sqrt2$, the least cost of a path between two cells whose
coordinates differ by $\Delta x$ and $\Delta y$ is $d_8$.

*Proof.* A path with $k$ diagonal and $l$ straight moves changes the sum
of the coordinates by at most $2k + l$ and each coordinate by at most $k
+ l$; hence $2k + l \ge \Delta x + \Delta y$ and $k + l \ge \max(\Delta x,
\Delta y)$. Its cost is
$\sqrt2\, k + l = (\sqrt2 - 1)(2k + l) + (2 - \sqrt2)(k + l) \ge (\sqrt2
- 1)(\Delta x + \Delta y) + (2 - \sqrt2) \max(\Delta x, \Delta y)$,
which equals $d_8$ since $\Delta x + \Delta y = \max + \min$. The path of
$\min(\Delta x, \Delta y)$ diagonal moves followed by $|\Delta x - \Delta
y|$ straight ones attains the bound. $\square$

**Proposition 14.6** (The grid heuristics are consistent).
<a id="prop-14-6"></a> On a grid with 4 moves, $h = \kappa_{\min}\, d_4$
is consistent; on a grid with 8 moves, $h = \kappa_{\min}\, d_8$ is
consistent.

*Proof.* Both $d_4$ and $d_8$ are shortest-path distances on the empty
grid ([Lemma 14.5](#lem-14-5) for $d_8$), and therefore they satisfy the
triangle inequality $d(u, t) \le d(u, v) + d(v, t)$. For a move from $u$
to a neighbour $v$, $d(u, v)$ is 1 for a straight move and $\sqrt2$ for
a diagonal one, while $c(u, v)$ is $\kappa(u) \ge \kappa_{\min}$ times
the same factor; hence $h(u) = \kappa_{\min} d(u, t) \le \kappa_{\min}
d(u, v) + \kappa_{\min} d(v, t) \le c(u, v) + h(v)$, and $h(t) = 0$.
Note that the blocked cells and the rule on the corners only remove
edges, which keeps the inequality on the remaining ones. $\square$

We remark that the Manhattan distance is *not* admissible with 8 moves,
whether a diagonal move costs $\sqrt2$ or 1: from $(0, 0)$ to $(1, 1)$
it gives 2, while one diagonal move costs $\sqrt2$ (or 1). Also, with 8
moves the costs and the heuristic are computed in double precision, and
[Proposition 14.6](#prop-14-6) holds in exact arithmetic: two paths to
the same cell whose exact costs are equal may differ by an ulp, and a
search may then (rarely) reopen a cell, which costs time but not
optimality. As a check, the test `grid` compares A$^*$ with the algorithm of
Dijkstra, written apart from the domain, on random maps of unit costs
and of costs 1 to 5, with 4 and 8 moves.

### The refinements and the maps

The features of Zobrist hashing are the column and the row of the cell,
the abstract features of AZH the column and the row of the block of $B
\times B$ cells containing it ($B$ = `block()`, 100 by default), and the
abstract state of AHDA$^*$ that block, the two distributions coinciding
on this domain as in section 4.1 of [[5]](Z-bibliography.md#ref-5). In
turn, the nblocks of PBNF are the blocks of $S \times S$ cells ($S$ =
`nblock_side()`, 20 by default), and the abstract successors of an
nblock are its 4 or 8 neighbouring nblocks ([Figure 14.7](#fig-14-7));
the pairwise heuristic of PA$^*$SE is the distance between two cells. A
grid models neither `IncrementalFeatureDomain` (a move changes one or
two of its features, but the domain does not enumerate them) nor
`OperatorDomain` (see above).

![Figure 14.7. (a) The blocks of AZH and AHDA$^*$ on a grid, coloured
by owner: a move changes the owner only across a thick line. (b) The
nblocks of Safe PBNF: the duplicate detection scope of the nblock $b$ is
the set of its neighbouring nblocks, which no other worker may use while
$b$ is acquired ([Chapter 11](11-safepbnf.md#ch-11)).
](../figures/domains-grid-blocks.svg){#fig-14-7}

A map is read from a file of the Moving AI benchmarks of Sturtevant
[[15]](Z-bibliography.md#ref-15) by `GridMap::read_movingai` (a header
`type octile`, `height`, `width`, `map`, then the rows, in which `.`,
`G` and `S` are free cells of cost 1 and any other character is
blocked), or drawn at random by `GridMap::random( W , H , p , seed ,
c_max )`, which blocks each cell with probability $p$ except the
bottom-left and the bottom-right corners, the start and the goal of the
instances of Burns et al. [[3]](Z-bibliography.md#ref-3), and gives each
free cell the cost 1 or a cost drawn uniformly in $\{1, \ldots,
c_{\max}\}$, from the words of `std::mt19937_64` and hence a seed gives
the same map on all platforms ([§1.6](01-introduction.md#sec-1-6)). In
the examples, `FindPath` finds and draws a path on such a map.

## 14.5 Explicit graphs and the map of Romania {#sec-14-5}

`ExplicitGraph< Cost >` (`ExplicitGraph.h:138`) is a weighted graph given
by its arcs, stored in an `ExplicitAdjacency` of the chosen set of data
structures ([Chapter 3](03-data-structures.md#ch-3)), with nodes
$0, \ldots, n - 1$: `add_edge( u , v , c )` adds an
edge (or an arc, with `directed`), `set_heuristic( u , h )` gives $h$ node
by node (0 by default), and the search goes from a source to a target.
It serves the small graphs of the tests (e.g., the graph of
[Figure 1.3](01-introduction.md#fig-1-3)) and the map of Romania; since
$h$ is given, its admissibility and consistency are the responsibility of
whoever gives it. Its arcs are kept in the adjacency of the chosen set of
data structures ([Chapter 3](03-data-structures.md#ch-3)), in the order
of their insertion but with the set of Boost, and the default set drops a
loop $(u, u)$, which never shortens a path.

`Romania` (`Romania.h:61`) is the road map of figure 3.2 of Russell and
Norvig [[13]](Z-bibliography.md#ref-13), with 20 cities, 23 roads whose
cost is their length in km, and real costs. Towards Bucharest the
heuristic is the straight-line distance of figure 3.22 of the book;
towards any other city it is the Euclidean distance between the
coordinates of the cities in the code that accompanies the book
(aima-python), which on this map never exceeds the length of a road (one
can check it on the data, see below).

```cpp
--8<-- "snippets/domains-romania-heuristic.cpp"
```

**Proposition 14.7** (The heuristics of Romania are consistent).
<a id="prop-14-7"></a> On the map of Romania, the straight-line distance
to Bucharest of the book and the Euclidean distance between the
coordinates of aima-python, towards any city, are consistent.

*Proof.* For the Euclidean distance $e$, the triangle inequality gives
$e(u, t) \le e(u, v) + e(v, t)$. Also, $e(u, v)$ is at most the length
of the road $(u, v)$ for each of the 23 roads (the largest ratio of
$e(u, v)$ to the length of the road is about 0.996, which we checked on
the data of `Romania.h`); hence $h(u) \le c(u, v) + h(v)$, and $h(t) =
0$. However, the table of the book is not derived from the coordinates,
and its consistency is checked edge by edge: for each of the 23 roads
$(u, v)$, $|h(u) - h(v)| \le c(u, v)$, the least slack being 1 km, on
the road from Pitesti to Bucharest ($100 - 0 \le 101$), and the next one
4 km, from Rimnicu Vilcea to Pitesti ($193 - 100 \le 97$). Since the
roads are two-way, the inequality holds in both directions. $\square$

Of course, the proof is a finite check, which the test `domains`
performs on each edge towards every city, and the test `graph` checks
what matters to the searches: A$^*$ finds the 418 km from
Arad to Bucharest of the book, through Sibiu, Rimnicu Vilcea and Pitesti
([Figure 14.8](#fig-14-8)), and the distance to Bucharest from each
city, against a table computed apart. Also, the example `RouteFinding`
finds a road between any two cities (e.g., `RouteFinding Lugoj Neamt`).

![Figure 14.8. The map of Romania of
[[13]](Z-bibliography.md#ref-13), with the length of each road in km and,
next to each city, its straight-line distance to Bucharest; the path
returned by A$^*$ from Arad to Bucharest (violet) costs
$140 + 80 + 97 + 101 = 418$.](../figures/domains-romania.svg){#fig-14-8}

## 14.6 Four problems of the textbook {#sec-14-6}

The next four domains are problems of chapter 3 of Russell and Norvig
[[13]](Z-bibliography.md#ref-13), each generalized by a size, by which
it can be made hard enough to measure a parallel search. They are small
in code, and each shows one shape of state space that the searches must
handle.

### The n-queens problem

`NQueens` (`NQueens.h:67`) places $n \le 32$ queens on an $n \times n$
board in the incremental formulation of section 3.2.1 of the book. A
state is the number $k$ of queens placed and the row of the queen of
each column $c < k$, and a move adds a queen in the column $k$, in a row
that no queen already on the board attacks, at cost 1. The state space
is a tree whose paths from the root have length at most $n$, and the
goals are the states with $k = n$; a state of the tree may have no
successor, a dead end ([Figure 14.9](#fig-14-9)).

```cpp
--8<-- "snippets/domains-queens-successors.cpp"
```

Its heuristic is the number of the empty columns, $h = n - k$, which is
exact on the states that lead to a solution: such a state needs $n - k$
more queens, and from a dead end no goal is reachable, $h^* = +\infty$.
Hence $h$ is admissible; it is also consistent, with equality $h(u) = 1
+ h(v)$ on each edge, and hence $f = n$ on all the states. Hence, the
ties of $f$ are broken by the larger $g$ ([Chapter
3](03-data-structures.md#ch-3)), and A$^*$ expands the states in
depth-first order, backtracking at the dead ends; the counters of the
example `NQueens` measure how many it meets. The 8-queens problem has 92
solutions, and its tree 2057 states (the root included); for $n = 2$ and
$n = 3$ the problem has no solution, and A$^*$ reports it after
exhausting the tree.

![Figure 14.9. The tree of the 4-queens problem, each state labelled
by the rows of the queens placed, column by column. The dead ends (rust)
have no free row in the next column; the two goals are the two solutions,
and all the states have $f = 4$.](../figures/domains-queens.svg){#fig-14-9}

### The vacuum world

`VacuumWorld` (`VacuumWorld.h:68`) is the vacuum world of section 3.2.1
of the book on a row of $n \le 32$ cells instead of 2: a robot moves
`Left` or `Right`, or `Suck`s the dirt of its cell, each action costing
1, and the goal is a clean row. A state packs the dirty cells in the
bits 0 to 31 of a word and the cell of the robot in the bits 32 to 37,
and the word is also the perfect hash. There are $n 2^n$ states, all
reachable from the state with all the cells dirty ([Figure
14.10](#fig-14-10)).

```cpp
--8<-- "snippets/domains-vacuum-successors.cpp"
```

Its heuristic, the number of dirty cells, is consistent: `Suck` lowers
it by 1 at cost 1, a move leaves it unchanged, and it is 0 at a goal; it
is also "weak", since it ignores the moves needed to reach the dirt,
which usually cost about as much as the cleaning. Furthermore, the
domain models the abstractions of the parallel searches: the abstraction
keeps the dirt and forgets the robot. Thus, the abstract features of AZH are
the dirty cells (the features of the robot being mapped to one abstract
feature, which all the states have) and the abstract state of AHDA$^*$
is the set of the dirty cells; hence, a move keeps the owner of a state,
and only a `Suck` changes it. Its nblocks of PBNF are the cells of the
robot, and the pairwise heuristic of PA$^*$SE between $s$ and $t$ is the
number of the cells dirty in $s$ and clean in $t$, which is consistent
for the same reason as $h$.

![Figure 14.10. The state space of the vacuum world of two cells of
the book: each state is the row with the robot (R) and the dirt (dots),
the edges are the actions, and the goals (dark violet) are the clean rows. The
path returned by A$^*$ (violet) is Suck, Right,
Suck.](../figures/domains-vacuum.svg){#fig-14-10}

### The water jugs

`WaterJugs` (`WaterJugs.h:67`) is the problem of exercise 3.6(d) of the
book (the `PourProblem` of aima-python). Up to 4 jugs of capacities of
at most 65535 units are empty at the start, a move fills a jug from the
faucet, empties it on the ground, or pours it into another jug until
either the first is empty or the second is full, and the goal is a jug
holding a given quantity ([Figure 14.11](#fig-14-11)). Every move costs
1, and the moves that do not change the state are not generated. In the
state, the level of the jug $j$ is packed in the bits $16j$ to $16j +
15$, and the word is also the perfect hash.

```cpp
--8<-- "snippets/domains-jugs-successors.cpp"
```

Its heuristic is 0 at a goal and 1 elsewhere, which is consistent, since
it changes by at most 1 along an edge of cost 1; it is (perhaps) the
best estimate that does not depend on the number theory of the
capacities (a goal is reachable only if the quantity is a multiple of
the greatest common divisor of the capacities and fits a jug). Note that
the problem is not reversible: emptying a jug cannot be undone in one
move.

![Figure 14.11. Measuring 4 units with jugs of 5 and 3: the six moves
of an optimal solution, with the levels of the two jugs under each
state.](../figures/domains-jugs.svg){#fig-14-11}

### The missionaries and cannibals

`MissionariesCannibals` (`MissionariesCannibals.h:65`) is the problem of
exercise 3.9 of the book with $m$ missionaries, $c$ cannibals and a boat
of $b$ seats: everybody starts on the left bank with the boat, a crossing
carries from 1 to $b$ people, and on neither bank may the missionaries be
outnumbered by the cannibals if any missionary is there. Every crossing
costs 1. A state packs the missionaries and the cannibals on the left bank
and the bank of the boat, and the word is the perfect hash; with 3, 3 and
2 the answer is 11 crossings, and with 4, 4 and 2 there is no solution.

Its heuristic counts the crossings that would carry the $n$ people left
on the left bank over the river with the boat always full: 0 if $n = 0$,
1 if $0 < n \le b$, $2 \lceil (n - b)/(b - 1) \rceil + 1$ if $n > b$ and
$b > 1$ ($2n - 1$ if $b = 1$), and one more if the boat is on the right
bank.

```cpp
--8<-- "snippets/domains-missionaries-heuristic.cpp"
```

**Proposition 14.8** (Admissible, not consistent). <a id="prop-14-8"></a>
The heuristic of `MissionariesCannibals` is admissible. It is not
consistent: with 3 missionaries, 3 cannibals and a boat of 2 seats, the
first crossing of an optimal solution lowers it by 3.

*Proof.* Let $h_L(n)$ be the value with the boat on the left bank and
$n$ people there. With the boat on the left, a solution alternates
crossings to the right, each carrying at most $b$ people, and crossings
back, each bringing at least 1. After $k$ round trips at most $k(b - 1)$
people have left, and the last crossing carries at most $b$; hence $n
\le k(b - 1) + b$, i.e., $k \ge \lceil (n - b)/(b - 1) \rceil$ for $n >
b$ and $b > 1$, and the solution has at least $2k + 1$ crossings (for $b
= 1$ no round trip makes progress, no solution exists when $n > 1$, and
any value is a lower bound). With the boat on the right and $n > 0$, the
first crossing brings at least one person back, and the rest needs at
least $h_L(n + 1) \ge h_L(n)$ crossings, $h_L$ being nondecreasing;
hence $h^* \ge 1 + h_L(n)$, the value of $h$. For the second claim, the
start has $n = 6$ and $h = 2 \lceil 4/1 \rceil + 1 = 9$; carrying a
missionary and a cannibal leaves $n = 4$ with the boat on the right, $h
= 2 \lceil 2/1 \rceil + 1 + 1 = 6$, and $9 > 1 + 6$ ([Figure
14.12](#fig-14-12)). $\square$

The proof shows where the slack may be: on the right bank the heuristic
adds 1 to $h_L(n)$, while the argument gives $1 + h_L(n + 1)$, which is
larger by 2 for some $n \ge b$ (for each $n \ge b$ when $b = 2$). A$^*$
is optimal with an admissible heuristic, reopening the closed states
when needed ([Chapter 2](02-best-first.md#ch-2)); the searches that
assume consistency are not. For instance, the search with delayed
duplicate detection of [Chapter 4](04-memory-bounded.md#ch-4) refuses a
successor of a smaller $f$ than its parent with an exception
(`DDDAStar.h:59`), and it does so on this domain.

![Figure 14.12. The 11 crossings of an optimal solution with 3
missionaries, 3 cannibals and a boat of 2 seats; a state is the numbers
of missionaries and cannibals on the left bank and the bank of the boat.
The heuristic drops by 3 on the first crossing (rust), which costs 1
([Proposition 14.8](#prop-14-8)).
](../figures/domains-missionaries.svg){#fig-14-12}

## 14.7 The traveling salesperson problem {#sec-14-7}

`TravelingSalesman` (`TravelingSalesman.h:75`) finds a tour of least
cost through $n \le 32$ cities, starting and ending at the city 0, given
the matrix of the distances (`random_cities( n , seed )` draws $n$
points of the unit square and gives their Euclidean distances). A state
is the set of the visited cities and the current one, packed as the
visited cities in the bits 0 to 31, the current city in the bits 32 to
36 and, in the bit 37, whether the tour has gone back to the city 0; the
word is also the perfect hash. A move goes to an unvisited city $k$ at
cost $d(c, k)$, or, once all the cities are visited, back to the city 0
at cost $d(c, 0)$, which reaches the goal. Its state space has $(n - 1)
2^{n - 2} + 2$ states, i.e., the start, the closed tour, and the pairs
of a nonempty set of visited cities other than 0 and of a current city
in it. One can compare this number with the $(n - 1)!/2$ tours of a
symmetric instance: this is the state space of the dynamic program of
Held and Karp [[14]](Z-bibliography.md#ref-14), and A$^*$ on it explores
a part of that table, the part selected by the heuristic.

```cpp
--8<-- "snippets/domains-tsp-successors.cpp"
```

Here, the heuristic of a state $u$ with current city $c$ and set $U$ of
unvisited cities is the cost $\mathrm{MST}(U \cup \{c, 0\})$ of a minimum
spanning tree of the unvisited cities together with the current and the
first one, computed by the algorithm of Prim in $O(n^2)$ time (exercise
3.30 of [[13]](Z-bibliography.md#ref-13)); it is 0 at the goal.

```cpp
--8<-- "snippets/domains-tsp-heuristic.cpp"
```

**Proposition 14.9** (The MST heuristic is consistent).
<a id="prop-14-9"></a> If the distances are symmetric and nonnegative,
the heuristic of `TravelingSalesman` is consistent, and hence
admissible.

*Proof.* Let $u$ have current city $c$ and unvisited set $U \ne
\emptyset$, and let $v$ be the successor that moves to $k \in U$; thus,
$v$ has current city $k$ and unvisited set $U \setminus \{k\}$, and
$h(v) = \mathrm{MST}(U \cup \{0\})$. A minimum spanning tree of $U \cup
\{0\}$ together with the edge $(c, k)$ is a spanning tree of $U \cup
\{c, 0\}$, of cost $h(v) + d(c, k)$; hence $h(u) \le d(c, k) + h(v) =
c(u, v) + h(v)$ ([Figure 14.13](#fig-14-13)). When $U = \emptyset$ and
$c \ne 0$, the only successor is the goal, at cost $d(c, 0)$, and $h(u)
= \mathrm{MST}(\{c, 0\}) = d(c, 0) = c(u, v) + 0$; at the goal $h = 0$.
$\square$

We remark that the proof needs the symmetry, since a spanning tree is a
set of undirected edges; on an asymmetric matrix the procedure of Prim
of the excerpt still returns a number (it reads $d$ from the tree
towards the other cities), but the proposition no longer applies.
Admissibility is also seen directly: the rest of a tour from $u$ is a
path from $c$ through all the cities of $U$ to 0, i.e., a spanning tree
of $U \cup \{c, 0\}$, and it costs at least the minimum one.

![Figure 14.13. The MST heuristic on seven cities. (a) At city 5,
after the cities 0 and 2 (violet path), $h$ is the cost of the tree
(dashed) on the unvisited cities and on 5 and 0. (b) After the move to 3,
the tree on the cities left and on 3 and 0, together with the edge
$(5, 3)$, spans the set of (a); hence $h$ drops by at most the cost of
the move ([Proposition 14.9](#prop-14-9)).
](../figures/domains-tsp.svg){#fig-14-13}

The header also gives the dynamic program of Held and Karp, as
`held_karp()`, by which the example `TravelingSalesman` checks the cost
found by A$^*$ up to 20 cities. With $C(S, j)$ the least cost of a path
that starts at the city 0, visits exactly the cities of the set $S \ni 0$
and ends at $j \in S$, it rests on the recurrence $C(S \cup \{k\}, k) =
\min_{j \in S} \{C(S, j) + d(j, k)\}$ for $k \notin S$, and the optimal
tour costs $\min_{j \ne 0} \{C(V, j) + d(j, 0)\}$, with $V$ the set of all
the cities.

**Algorithm 14.1** (The dynamic program of Held and Karp).

```text
1  procedure HELDKARP(d, n)
2     C({0}, 0) ← 0;  C(S, j) ← ∞ for each other pair
3     for each set S ∋ 0, by increasing value of its bits, do
4        for each j ∈ S with C(S, j) < ∞ do
5           for each k ∉ S do
6              C(S ∪ {k}, k) ← min(C(S ∪ {k}, k), C(S, j) + d(j, k))
7     return min over j ≠ 0 of C(V, j) + d(j, 0)
```

Line 3 visits each set after all its subsets, since adding a city raises
the value of the bits; the table has $2^n n$ entries and lines 3 to 6 take
$O(2^n n^2)$ time, which is why the method serves to check the searches
on small instances only.

## 14.8 Slow expansions {#sec-14-8}

`SlowDomain< D >` (`SlowDomain.h:55`) wraps a domain `D` and adds to each
expansion a given number of iterations of the computation $i \leftarrow
11 i \bmod 9999943$, which keeps a core busy with no access to memory,
as in the experiment of section 3.1.1 of Jinnai and Fukunaga
[[5]](Z-bibliography.md#ref-5) on the effect of slow expansions. It models
the problems in which the expansion dominates the cost of a search
(planning with an expensive heuristic, or motion planning with collision
checks), the regime in which PA$^*$SE of
[Chapter 8](08-shared-memory.md#ch-8) is meant to pay.

```cpp
--8<-- "snippets/domains-slow-work.cpp"
```

Since the result of the loop is added to a thread-local variable, the
compiler cannot remove it. The wrapper forwards all the refinements of
`D` (the operator included, with the work paid in both expansions) but
`IncrementalFeatureDomain`, and therefore the Zobrist keys of the
successors of a slow domain are computed from scratch; with an expansion
that is slow by construction, this cost is negligible. In the benchmark,
the option `--slow` gives the number of iterations per expansion.

## 14.9 Problems not in the library {#sec-14-9}

> **Status — planned.** The problems of this section are not domains of
> the library; the section states how they fit the model, and what the
> searches would need.

Two problems of the textbook of Jinnai (section 2.2 of
[[164]](Z-bibliography.md#ref-164)) and of the
literature are not provided. Sokoban, PSPACE-complete, is the typical
irreversible problem: the player moves at no cost and pushes a box at a
positive cost; a domain would thus have edges of cost 0 (which the
searches of the library accept, but which the delayed duplicate
detection refuses) and many dead ends (a box pushed into a corner),
which a heuristic must detect to be useful. Also, the Rubik's cube is
the classical benchmark of the pattern databases (and arguably the
hardest of the puzzles commonly used). Its state is the position and the
orientation of the 20 movable cubies, which fit an integer of 128 bits
(like the boards of `Word128`), a move turns one of the 6 faces by a
quarter or a half turn, and its reachable states are about $4.3 \cdot
10^{19}$. A domain would model `OperatorDomain` (a face is never turned
twice in a row) and the features of Zobrist hashing (cubie, position,
orientation), and its heuristic would (typically) be a maximum of
pattern databases of corners and of edges ([Chapter
5](05-heuristics.md#ch-5)). Finally, the multiple sequence alignment,
reduced to a shortest path on a lattice of $n$ dimensions with $2^n - 1$
moves, is the subject of [Chapter 15](15-msa-instances.md#ch-15).

## 14.10 Running the domains {#sec-14-10}

Each domain is run in two ways. First, the programs of `examples` solve
one instance with A$^*$ and print the path and its cost, the counters
and the time of the search. They are `8Puzzle` (the instance of figure
3.28 of the book, or a state at a given distance, with the misplaced
tiles, the Manhattan distance or the sequence score; with `-t` it
repeats the comparison of figure 3.29 of the book), `15Puzzle` (an
instance of Korf, sequentially or with HDA$^*$), `FindPath`,
`RouteFinding`, `NQueens`, `VacuumWorld`, `WaterJugs`,
`MissionariesCannibals` and `TravelingSalesman`; each run with a known
answer is also a test of `ctest`, labelled `example`. Second, the
benchmark `AStarBench` runs any search of the library, on any runtime
and with any distribution, on the domain chosen by `--domain` (`tiles8`,
`tiles15`, `tiles24`, `grid4`, `grid8`, `queens`, `vacuum`,
`missionaries`, `jugs`, `tsp`, `romania`): the puzzles read their
instances from `--file` (the instances chosen by `--ids`), the grids
from a Moving AI map or from `--random WxH:p:seed[:c]`, and the other
domains are sized by `--size` (e.g., `queens 28`, `missionaries
65535:65535:100`, `jugs 1:12:8:3`, `tsp 16` with the seeds of `--ids`).
It prints one line of CSV per run, with the cost, the counters of
[§1.3](01-introduction.md#sec-1-3), LB, CO and the seconds.

**Table 14.2.** The problems of the examples, solved by sequential
A$^*$.

--8<-- "tables/domains-examples.md"

Table 14.2 is written by `tables.py` from the table E of the campaign of
`tools/bench/campaign.sh`. The campaign runs each example program 10
times with the arguments of the second column (no argument being the
default instance of the program, e.g., the instance of figure 3.28 for
`8Puzzle` and the instance 1 of Korf for `15Puzzle`; `FindPath -r
1000x1000:0.3:2` is a random map of $1000 \times 1000$ cells with 30% of
them blocked and the others costing 1 to 5, with 4 moves and, with `-8`,
with 8). Its columns give the cost of the path returned (in moves, in km
or in the units of the domain, exact for the integer costs and in double
precision for the real ones), the number of expansions, which does not
depend on the machine, and the median of the 10 wall-clock times in
milliseconds. Note that the time is the field `seconds` of the
`Solution`, measured by the `Stopwatch` started at the beginning of
`solve()` (`AStar.h:129`) and read after the path has been traced; it
does not include the start of the process, the construction of the
domain (the reading or the drawing of the map, the precomputed tables of
the puzzle), nor the printing of the path. The note of the table gives
the machine, the compiler and the flags. Since the problems of the
examples are meant to be solved in milliseconds, we measure (in turn)
the scaling of the parallel searches on larger instances of the same
domains in [Chapter 8](08-shared-memory.md#ch-8) and [Chapter
9](09-hash-distributed.md#ch-9).

## 14.11 Notes and references {#sec-14-11}

The sliding-tile puzzles became the benchmark of heuristic search with
the 100 instances of the 15-puzzle that Korf
[[7]](Z-bibliography.md#ref-7) solved by IDA$^*$ with the Manhattan
distance; the parity argument of [Proposition 14.1](#prop-14-1) and its
converse are due to Johnson and Story [[150]](Z-bibliography.md#ref-150).
In turn, the 24-puzzle needs pattern databases
[[22]](Z-bibliography.md#ref-22). The instances of `data/24-puzzle` come
with the code of HDA$^*$ of Jinnai and Fukunaga
[[5]](Z-bibliography.md#ref-5), whose abstractions of the puzzles and of
the grids the domains of this chapter reproduce; the abstractions of
AHDA$^*$ and of PBNF are those of Burns et al.
[[3]](Z-bibliography.md#ref-3), and the packed state, the incremental
heuristic and the pruning of the parent follow Burns et al.
[[4]](Z-bibliography.md#ref-4). Our grids follow the benchmarks of
Sturtevant [[15]](Z-bibliography.md#ref-15) and the random maps of
[[3]](Z-bibliography.md#ref-3); the map of Romania, the n-queens, the
vacuum world, the water jugs, the missionaries and cannibals and the MST
heuristic of the traveling salesperson come from chapter 3 of Russell
and Norvig [[13]](Z-bibliography.md#ref-13), the sequence score from
Nilsson [[8]](Z-bibliography.md#ref-8), and the dynamic program of the
traveling salesperson from Held and Karp
[[14]](Z-bibliography.md#ref-14).

It is fair to remark that the literature of the parallel searches uses
these domains, and therefore they allow the searches of the library to
be compared with the published ones; however, they are written by hand
for one problem each, with hand-made abstractions, and the choice of the
abstraction from the problem itself is open ([Chapter
10](10-auto-distributions.md#ch-10)). Also, the representation of the
puzzles in one or two words is that of the fastest published codes
[[4]](Z-bibliography.md#ref-4), while the other domains use a word for
convenience more than for speed; their expansions are cheap, and on them
the measures of the parallel searches are dominated by the bookkeeping
of the searches, which is the regime in which those searches are hardest
to scale.

## 14.12 Exercises {#sec-14-12}

**Exercise 14.1** (discussion, no test). Model the Rubik's cube as a
`Domain`: give the state and its packing in a word, the 18 moves (a
quarter turn either way and a half turn of each face) and their
operators, the features for Zobrist hashing, and a pair of refinements
that the cube can model; explain why the operator should also skip the
moves of the opposite face in one of the two orders.

**Exercise 14.2** (`ctest -R PBeFS.astar`). Show that the 8-puzzle has
$9!/2 = 181440$ reachable states and the 15-puzzle $16!/2$, using
[Proposition 14.1](#prop-14-1) and its converse; the test checks the first
number by a breadth-first search of the whole space from the goal.

**Exercise 14.3** (`ctest -R PBeFS.domains`). Compute the mean branching
factor of a grid of $5 \times 5$ free cells with 4 moves (it is $80/25 =
3.2$) and with 8 moves ($144/25 = 5.76$); then measure it with the
counters of `AStarBench`, as `generated / expanded` on a search that
expands all the cells (e.g., towards a blocked goal), and explain the
difference between the two numbers. The test checks the two counts of
moves.

**Exercise 14.4** (`ctest -R PBeFS.domains`). Prove
[Proposition 14.2](#prop-14-2) in detail, i.e., that the mean branching
factor of the 8-puzzle over its reachable states is exactly $8/3$, and
$5/3$ without the parent; the test checks the first value on the whole
state space.

**Exercise 14.5** (`ctest -R PBeFS.grid`). Show that on a grid with 8
moves the heuristic $\kappa_{\min}\, d_4$ (the Manhattan distance) is not
admissible, and that $\kappa_{\min} \max(\Delta x, \Delta y)$ (the
Chebyshev distance) is admissible but less informed than $\kappa_{\min}
d_8$. The test checks A$^*$ with $d_8$ against Dijkstra on random maps.

**Exercise 14.6** (`ctest -R PBeFS.domains`). Prove that the heuristic
$1 + h_L(n + 1)$, used on the right bank in place of $1 + h_L(n)$, is
still admissible for the missionaries and cannibals, and decide whether
it is consistent. The test checks that the heuristic of the library is
admissible, by the distances to the goal of all the reachable states,
but not consistent, on the instance with 3, 3 and 2.

**Exercise 14.7** (`ctest -R PBeFS.domains`). Count the states of the
traveling salesperson problem of [§14.7](#sec-14-7) with $n$ cities, and
show that A$^*$ with the MST heuristic never expands more states than
there are entries in the table of Held and Karp. The test compares the
cost of A$^*$ with `held_karp()` on random instances and checks the
consistency of the heuristic on all the edges of their state spaces.

**Exercise 14.8** (`ctest -R PBeFS.domains`). Show that A$^*$ on
`NQueens` expands the states in the order of a depth-first search that
tries the rows in decreasing order (the successors are given by
increasing row, and the ties of $f$ and $g$ are broken in favor of the
last pushed), and count the states of the tree for
$n = 8$; the test checks that the tree has 2057 states and 92 goals, and
that the problem with 3 queens is unsolvable.

**Exercise 14.9** (`ctest -R PBeFS.graph`). Check [Proposition
14.7](#prop-14-7) on the data of `Romania.h`, and explain why the table
of the book cannot serve as the heuristic towards a city other than
Bucharest, while the coordinates can. The test checks the distances to
Bucharest from each city.

