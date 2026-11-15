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
in it), the operators and their costs as pseudocode, the heuristic with
the proof of its admissibility or consistency (or a counterexample where
it fails), the refinements of `Domain` it models with the proofs of the
conditions that the searches put on them, its instances, and how the
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
and the multiple sequence alignment to [Chapter
15](15-msa-instances.md#ch-15).

## 14.1 The domains and their concepts {#sec-14-1}

Table 14.1 gives, for each domain, the refinements of `Domain` it models;
they decide which searches and which work distributions of HDA$^*$ can be
run on it ([Figure 1.5](01-introduction.md#fig-1-5), where each
refinement is introduced). A domain built on another one inherits its
refinements where the code forwards them (`SequenceScore` forwards none,
`SlowDomain` all but `IncrementalFeatureDomain`), and the last column is
the concept `DiskDomain` of the search with delayed duplicate detection
of [Chapter 4](04-memory-bounded.md#ch-4), i.e., integer costs and
trivially copyable states (`DDDAStar.h:136`).

**Table 14.1.** The refinements of `Domain` modeled by each domain: F
feature, IF incremental feature, AF abstract feature, AS abstract state,
NB nblock, PW pairwise, OP operator, PH perfect hash, DK disk (x
modeled, (a) only up to 20 cells, (b) only with 4 moves, (c) only with
an integral `Cost`; the cost type of the grids is `int` with 4 moves and
`double` with 8).

| domain | cost | F | IF | AF | AS | NB | PW | OP | PH | DK |
|---|---|---|---|---|---|---|---|---|---|---|
| `SlidingTile< W , H >` | `int` | x | x | x | x | x | x | x | (a) | x |
| `SequenceScore< W , H >` | `int` | | | | | | | | | x |
| `Grid< 4 >`, `Grid< 8 >` | both | x | | x | x | x | x | | x | (b) |
| `ExplicitGraph< Cost >` | `Cost` | x | | | | | | | x | (c) |
| `Romania` | `double` | x | | | | | | | x | |
| `NQueens` | `int` | x | | | | | | | | x |
| `VacuumWorld` | `int` | x | | x | x | x | x | | x | x |
| `WaterJugs` | `int` | x | | | | | | | x | x |
| `MissionariesCannibals` | `int` | x | | | | | | | x | x |
| `TravelingSalesman` | `double` | x | | | | | | | x | |

Every domain packs its state into an integer, or into a small trivially
copyable record (the board of `NQueens`). A state is thus compared and
copied as a word, and its hash is the mixing function of SplitMix64
applied to it ([§1.5](01-introduction.md#sec-1-5)), which is injective
on the words of 64 bits ([Lemma 1.7](01-introduction.md#lem-1-7)); this
is the representation that Burns et al. [[4]](Z-bibliography.md#ref-4)
recommend, since the state is then at the same time its key and the
input of its hash. The declarations by which a program builds the
instances are the constructors below (quoted by hand from the headers,
without their bodies; `board_type` is an array of `W * H` bytes,
`City` an enumeration of the 20 cities and `matrix_type` a vector of
vectors of `double`).

```cpp
explicit SlidingTile( const board_type & start ,
                      const board_type & goal = korf_goal() );
explicit SequenceScore( const board_type & start ,
                        const board_type & goal = Base::korf_goal() );
static GridMap GridMap::random( int width , int height , double p ,
                                std::uint64_t seed , int max_cost = 1 );
static GridMap GridMap::read_movingai( std::istream & in );
Grid( const GridMap & map , int sx , int sy , int gx , int gy ,
      int block = 100 ,
      int nblock_side = 20 );
ExplicitGraph( std::size_t n , state_type s , state_type t );
Romania( City from = Arad , City to = Bucharest );
explicit NQueens( int n );
VacuumWorld( int n , int start , std::uint32_t dirty );
WaterJugs( std::vector< int > capacity , int goal );
MissionariesCannibals( int m , int c , int b );
explicit TravelingSalesman( matrix_type dist );
SlowDomain( const D & d , unsigned long work );
```

The example below, compiled and run by the tester of the library (in
the case `manual`, `ctest -R PBeFS.manual`), builds the 8-puzzle of
figure 3.28 of [[13]](Z-bibliography.md#ref-13)
and a random grid with 4 and with 8 moves, solves them with the
sequential A$^*$ of [Chapter 2](02-best-first.md#ch-2), and checks the
costs: 26 moves for the puzzle with both OPEN lists, the cost found by
the algorithm of Dijkstra (i.e., by A$^*$ with $h = 0$) for the grid
with 4 moves, and, for the grid with 8 moves, a cost not larger than
that with 4.

```cpp
--8<-- "snippets/domains-usage.cpp"
```

## 14.2 The sliding-tile puzzles {#sec-14-2}

The sliding-tile puzzle on a board of $W \times H$ cells has $WH - 1$
tiles, numbered from 1, and one empty cell, the *blank*; a move slides a
tile adjacent to the blank into it, at cost 1, and the goal is a given
arrangement of the tiles. In the library, the class template
`SlidingTile< W , H >` (`SlidingTile.h:170`) models it for $WH \le 25$
and $W, H \ge 2$, with the aliases `EightPuzzle`, `FifteenPuzzle` and
`TwentyFourPuzzle` for the boards of $3 \times 3$, $4 \times 4$ and $5
\times 5$ cells; the goal is by default that of Korf
[[7]](Z-bibliography.md#ref-7) (`korf_goal()`), with the blank in the
first cell and the tiles $1, 2, \ldots$ in row-major order.

### The state and the moves

A state is packed into one integer, in which the field of $\beta$ bits
of the cell $i = rW + k$ holds the tile in row $r$ and column $k$ (0 for
the blank), i.e., the integer $s = \sum_{i} \sigma(i)\, 2^{\beta i}$ for
the board $\sigma$. With $\beta = 4$ the boards of up to 16 cells fit a
word of 64 bits, and with $\beta = 5$ the larger ones fit a `Word128`,
an unsigned integer of 128 bits written on two words of 64 in order to
be the same with all compilers (`SlidingTile.h:71`). This integer is at
the same time the state, its equality test and the key of its hash. A
move exchanges the blank, in the cell $b$, with the tile $t$ in a cell
$c$ adjacent to it: the field of $c$ is cleared and $t$ is written in
the field of $b$, two operations on the word ([Figure 14.1](#fig-14-1)).

![Figure 14.1. A move of the 8-puzzle, on the instance of figure 3.28
of [[13]](Z-bibliography.md#ref-13) ($u$, whose optimal solution has 26
moves, towards the goal on the right). Sliding tile 2 down into the blank
moves it from distance 1 to distance 2 from its goal cell, and the
Manhattan distance of the successor follows from that of $u$ by the change
of that tile alone.](../figures/domains-tiles-move.svg){#fig-14-1}

The cell of the blank is not stored: it is found from the word by a few
bitwise operations (Algorithm 14.1, `SlidingTile.h:465`, and
[Figure 14.2](#fig-14-2)), in which $L = \sum_i 2^{\beta i}$ is the mask
of the lowest bit of each field, $\gg$ the right shift, $\mid$, $\&$ and
$\sim$ the bitwise or, and and complement, and TZ the number of trailing
zeros of a word.

**Algorithm 14.1** (The cell of the blank).

```text
1  procedure BLANK(s)
2     x ← s
3     for k ← 1 to β − 1 do
4        x ← x | (s >> k)       ▷ or the bits of each field into its lowest
5     x ← ~x & L                ▷ one bit set: the lowest of the zero field
6     return TZ(x) / β          ▷ integer division
```

<a id="lem-14-1"></a>**Lemma 14.1** (The cell of the blank). For every
board $\sigma$ with its blank in the cell $b$, BLANK$(s)$ returns $b$,
with $\beta - 1$ shifts and ors and one count of the trailing zeros.

*Proof.* The bit $\beta i$ of $x$ after line 4 is the or of the bits $\beta
i, \beta i + 1, \ldots, \beta i + \beta - 1$ of $s$, all of which belong
to the field of the cell $i$ (the shifts bring into the higher bits of a
field the bits of the next field, but never into its lowest bit, since
$k \le \beta - 1$). Hence, the bit $\beta i$ of $x$ is 0 if and only if
the field of $i$ is 0, i.e., if and only if $i = b$, the tiles being
nonzero and the blank unique. After line 5, $x = 2^{\beta b}$, and its
trailing zeros divided by $\beta$ give $b$. $\square$

![Figure 14.2. The cell of the blank of the instance of
[Figure 14.1](#fig-14-1), from its packed word (four bits per cell): the
lowest bit of $s \mid s{\gg}1 \mid s{\gg}2 \mid s{\gg}3$ is 0 only in the
field of the blank, and the trailing zeros of the complement (masked by
the lowest bits) give the cell 4.
](../figures/domains-tiles-packing.svg){#fig-14-2}

The cells adjacent to each cell are precomputed at construction, in the
order up, left, right and down (`SlidingTile.h:220`), in a table $N(b)$
of at most 25 entries of 5 integers (hence usually in the cache), and so
is the table $\delta_t(c)$ of the Manhattan distance of the cell $c$
from the goal cell of the tile $t$ (`SlidingTile.h:217`, with $\delta_0
= 0$ for the blank). Also, the domain is an `OperatorDomain` whose
operator is the cell of the blank in the parent, and Algorithm 14.2 is
its expansion with an operator (`SlidingTile.h:292`); the expansion
without an operator is the same without line 5 (`SlidingTile.h:274`).

**Algorithm 14.2** (The successors of a board).

```text
1  procedure TILESUCCESSORS(s, h(s), o, f)   ▷ o: the cell the blank left
2     b ← BLANK(s)
3     for each cell c ∈ N(b), in the order up, left, right, down, do
4        if c = o then
5           continue                ▷ the move back to the parent
6        t ← the tile in the field of c of s
7        s' ← s with the field of c set to 0 and that of b set to t
8        f(s', 1, h(s) + δ_t(b) − δ_t(c), b)   ▷ b: the operator of s'
```

Line 7 is $s' = (s \mathbin{\&} \sim(m \ll \beta c)) \mid (t \ll \beta
b)$, with $m = 2^\beta - 1$ the mask of a field; since the field of $b$
is 0 in $s$, the or writes $t$ there, and $s'$ is the packed word of the
board with $t$ and the blank exchanged. Hence an expansion costs one
BLANK and, per successor, a constant number of operations on the word
and two lookups of $\delta$, with no loop over the cells; the value of
line 8 is the Manhattan distance of $s'$ by [Proposition
14.5](#prop-14-5) below. The inverse of the move from $s$ to $s'$ is the
move of the blank of $s'$ (now in $c$) back to $b$, i.e., to the
operator of $s'$, at the same cost; hence the hypothesis of
[Lemma 1.8](01-introduction.md#lem-1-8) holds, and line 5 skips only
duplicates.

### The reachable states

Only half of the arrangements of the tiles can be reached from a given
one, and the constructor refuses a start and a goal that are not in the
same half (`SlidingTile.h:212`), as well as a board that is not a
permutation of $0, \ldots, WH - 1$. We now state the invariant, after
fixing the notation. For a board $\sigma$, let $N(\sigma)$ be the number
of *inversions* of the tiles read in row-major order without the blank
(the pairs of tiles $t > t'$ with $t$ before $t'$), and $r(\sigma)$ the
row of the blank. The function `parity` (`SlidingTile.h:578`) computes
$N(\sigma) \bmod 2$ if $W$ is odd and $(N(\sigma) + r(\sigma)) \bmod 2$
if $W$ is even, by a double loop over the cells, in $O((WH)^2)$ time.

<a id="prop-14-2"></a>**Proposition 14.2** (Parity of the sliding-tile
puzzles). If $W$ is odd, $N(\sigma) \bmod 2$ is invariant under the
moves; if $W$ is even, $(N(\sigma) + r(\sigma)) \bmod 2$ is. Hence, at
most half of the $(WH)!$ boards are reachable from a given one.

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
changes $N$ by an odd number and leaves $r$ unchanged, and it is
therefore a bijection between the two halves. $\square$

The converse, that every board on which the invariant takes the value of
the start is reachable, is the theorem of Johnson and Story
[[150]](Z-bibliography.md#ref-150) for the 15-puzzle, and it holds for all
the boards of at least $2 \times 2$ cells; we use it without
reproducing its proof. Hence, the 8-puzzle has $9!/2 = 181440$
reachable states, which the test `domains` checks by a search of the
whole space, the 15-puzzle $16!/2 \approx 1.05 \cdot 10^{13}$ and the
24-puzzle $25!/2 \approx 7.8 \cdot 10^{24}$.

<a id="lem-14-3"></a>**Lemma 14.3** (The blank is uniform). Assume the
converse of [Proposition 14.2](#prop-14-2). Then, for each cell $b$,
exactly $(WH - 1)!/2$ reachable boards have the blank in $b$.

*Proof.* With the blank fixed in $b$, the boards are the $(WH - 1)!$
arrangements of the tiles on the other cells, and $r(\sigma)$ is the same
for all of them; hence the invariant is the parity of $N$, possibly
shifted by the constant $r$. Exchanging two tiles changes the parity of
$N$, and is a bijection between the arrangements of the two parities,
which are therefore $(WH - 1)!/2$ each; by the converse, those of the
right parity are exactly the reachable ones. $\square$

On the puzzles, the branching factor of a state is the number of the
neighbours of its blank, 2 in a corner, 3 on a side and 4 inside the
board, and by the lemma its mean over the reachable states has a closed
form ([Figure 14.3](#fig-14-3)).

<a id="prop-14-4"></a>**Proposition 14.4** (Branching factor of the
puzzles). Over the reachable states of the puzzle of $W \times H$ cells,
the mean number of moves is $2(2WH - W - H)/(WH)$, i.e., $8/3$ for the
8-puzzle, 3 for the 15-puzzle and $16/5$ for the 24-puzzle; without the
move back to the parent it is 1 less on every state but $s_0$.

*Proof.* By [Lemma 14.3](#lem-14-3) the blank is uniformly distributed
over the cells, and the mean number of moves is the mean degree of a
cell in the grid graph of the board, $2|E|/(WH)$, with $|E| = H(W - 1)
+ W(H - 1) = 2WH - W - H$ edges. The three values follow with $|E| =
12, 24, 40$ (cf. [Figure 14.3](#fig-14-3)). Skipping the parent removes
exactly one successor of each state but $s_0$, since the cell left by
the blank is adjacent to its current cell. $\square$

![Figure 14.3. The number of moves of the blank in each cell of the
boards of $3 \times 3$ and $4 \times 4$ cells. Since the blank is
uniformly distributed over the cells of the reachable states, the mean
branching factor is the mean of these numbers
([Proposition 14.4](#prop-14-4)).
](../figures/domains-tiles-blank.svg){#fig-14-3}

### The Manhattan distance

The heuristic is the Manhattan distance: the sum, over the tiles, of the
distance between the cell of the tile and its goal cell in moves of the
grid, $h(\sigma) = \sum_{t \ne 0} \delta_t(c_\sigma(t))$, with
$c_\sigma(t)$ the cell of $t$ in $\sigma$; computed from scratch, by
`heuristic( s )`, it takes one lookup per cell.

<a id="prop-14-5"></a>**Proposition 14.5** (The Manhattan distance is
consistent). On the sliding-tile puzzle, for each move from $\sigma$ to
$\sigma'$ of a tile $t$ from the cell $c$ to the cell $b$ of the blank,
$h(\sigma') = h(\sigma) + \delta_t(b) - \delta_t(c)$ and $|h(\sigma') -
h(\sigma)| = 1$; hence, $h$ is consistent, and therefore admissible.

*Proof.* Only the tile $t$ changes its cell, which gives the first
identity. The cells $b$ and $c$ are adjacent, so that one coordinate
changes by 1 and the other not; hence their distances from the goal cell
of $t$ differ by exactly 1, since the distance changes by $\pm 1$ in
that coordinate. Therefore $h(\sigma) \le 1 + h(\sigma') = c(\sigma,
\sigma') + h(\sigma')$, and $h$ is 0 at the goal; consistency implies
admissibility ([Lemma 2.4](02-best-first.md#lem-2-4)). $\square$

The first identity is how Algorithm 14.2 computes $h$, two lookups
instead of a sum over the tiles, which is the incremental heuristic of
[[4]](Z-bibliography.md#ref-4). The second identity has a consequence on
the shape of the search.

<a id="cor-14-6"></a>**Corollary 14.6** (The parity of $f$). On the
puzzles with the Manhattan distance, $f = g + h$ changes by 0 or 2 along
each edge, and hence the values of $f$ of all the states generated by a
search from $s_0$ have the parity of $h(s_0)$.

*Proof.* Along an edge, $g$ grows by 1 and $h$ changes by $\pm 1$
([Proposition 14.5](#prop-14-5)), so that $f$ changes by $0$ or $2$. A
generated state has $f = g + h$ for the cost $g$ of a path from $s_0$,
and by induction on the length of the path its parity is that of $f(s_0)
= h(s_0)$. $\square$

Hence, iterative deepening raises its bound by exactly 2 at each
iteration on the puzzles, since the least $f$ above the bound
([Chapter 4](04-memory-bounded.md#ch-4)) has the parity of the bound.
Of course, the pattern databases of [Chapter 5](05-heuristics.md#ch-5)
replace this heuristic in `SlidingTilePDB`, which keeps the rest of the
domain.

### The refinements

The puzzles model all the refinements of `Domain` ([Figure
14.4](#fig-14-4)). Their feature $i = tWH + c$ is "tile $t$ is in cell
$c$", for $t \ne 0$ (`num_features()` is $(WH)^2$, of which those with
$t = 0$ are never used). A move changes two features, and
`for_each_feature_change` (Algorithm 14.3, `SlidingTile.h:332`) finds
them without knowing the moved tile, nor the blank.

**Algorithm 14.3** (The features changed by a move).

```text
1  procedure FEATURECHANGE(s, s', f)
2     d ← s xor s'                ▷ bitwise exclusive or of the two words
3     if d = 0 then
4        return
5     i ← TZ(d) / β             ▷ the first cell where s and s' differ
6     k ← the field of i of d   ▷ the moved tile
7     j ← TZ(d with the field of i cleared) / β
8     f(k WH + i);  f(k WH + j)
```

<a id="lem-14-7"></a>**Lemma 14.7** (The features changed by a move).
If $s'$ is a successor of $s$, Algorithm 14.3 calls `f` exactly on the
two features of the symmetric difference of the features of $s$ and
$s'$, in constant time.

*Proof.* Let the move slide the tile $k$ from $c$ to the cell $b$ of the
blank of $s$. The fields of $s$ and $s'$ coincide in every other cell;
in $b$ they are $0$ and $k$, and in $c$ they are $k$ and $0$. Hence $d$
has the field $k \ne 0$ in the two cells $b$ and $c$ and 0 elsewhere,
line 5 finds the lower of the two, line 6 reads $k$ in it, and line 7
the higher one. The features that differ are "$k$ in $c$" (in $s$ only)
and "$k$ in $b$" (in $s'$ only), the blank having no feature, and these
are the two calls of line 8. $\square$

Next, the abstract feature of AZH (`abstract_feature( i )`) keeps of
the cell of a tile only whether its row is at least $\lfloor H/2
\rfloor$ (the lower half of the board) or not, which is the projection
of figure 8(b) of Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) for
the 15-puzzle; there are $2WH$ abstract features, two per tile.

<a id="prop-14-8"></a>**Proposition 14.8** (The moves that change the
abstract features). A horizontal move never changes the abstract features
of a board; a vertical move changes them (and then two of them) if and
only if it moves a tile between the rows $\lfloor H/2 \rfloor - 1$ and
$\lfloor H/2 \rfloor$.

*Proof.* A move changes the cell of one tile only. A horizontal move
keeps its row, and hence the half of its cell; a vertical move changes
the row by 1, and the half changes if and only if one of the two rows is
below $\lfloor H/2 \rfloor$ and the other is not, i.e., if and only if
the rows are $\lfloor H/2 \rfloor - 1$ and $\lfloor H/2 \rfloor$; then
the abstract feature of the tile in the old half is replaced by that in
the new one. $\square$

Hence, under HDA$^*$[AZ], a successor can have an owner other than that
of its parent only through the moves across the middle of the board (one
pair of rows out of $H - 1$), and the others stay on the worker of their
parent; how often such moves occur along a search depends on the
instance and is only measured, not bounded, by the counters of
[Chapter 9](09-hash-distributed.md#ch-9). As in Burns et al.
[[3]](Z-bibliography.md#ref-3), the abstract state of AHDA$^*$
(`abstract_state( s )`) is the triple of the cells of tiles 1, 2 and 3,
coded as an integer.

The nblock of PBNF (`nblock( s )`) is the triple $(p_0, p_1, p_2)$ of the
cells of the blank and of tiles 1 and 2, coded as $(p_0 WH + p_1) WH +
p_2$ among `num_nblocks()` $= (WH)^3$ codes, of which only those of three
distinct cells are used. Its abstract successors
(`for_each_nblock_successor`, `SlidingTile.h:395`) are the triples
reached by a move of the blank to a cell $c \in N(p_0)$, which moves tile
1 or 2 if it lies in $c$: they are $(c, p_0, p_2)$ if $c = p_1$, $(c,
p_1, p_0)$ if $c = p_2$, and $(c, p_1, p_2)$ otherwise. Safe PBNF
requires the two conditions of
[Definition 11.1](11-safepbnf.md#def-11-1).

<a id="prop-14-9"></a>**Proposition 14.9** (The nblocks of the puzzles).
The nblocks of `SlidingTile` satisfy [Definition
11.1](11-safepbnf.md#def-11-1): (i) for each move from $u$ to $v$, the
nblock of $v$ is an abstract successor of that of $u$, and (ii) the
relation of the abstract successors is symmetric and irreflexive.

*Proof.* (i) A move takes the blank from $p_0$ to an adjacent cell $c$
and the tile of $c$ to $p_0$; tile 1 moves if and only if $c = p_1$, and
tile 2 if and only if $c = p_2$, which are the three cases of the
abstract successors. (ii) The blank changes cell in each case, so that
no nblock is its own successor. The successor of $(p_0, p_1, p_2)$
through $c$ has its blank in $c$, adjacent to $p_0$, and the move of
that blank back to $p_0$ puts back in $p_0$ the tile it moved, if any;
by the three cases this abstract successor is $(p_0, p_1, p_2)$.
$\square$

Also, the pairwise heuristic of PA$^*$SE (`pairwise_heuristic( s , t )`) is
the sum over the tiles of the distance between their cells in the two boards,
whose properties are proved with those of the grids and of the vacuum world in
[Proposition 14.18](#prop-14-18). Finally, the perfect hash
(`perfect_hash( s )`) is the lexicographic rank of the permutation of the
cells of the tiles $0, 1, \ldots, WH - 1$, i.e., the rank of the pattern
databases of [Chapter 5](05-heuristics.md#ch-5) with all the $WH$ tiles; by
[Lemma 5.4](05-heuristics.md#lem-5-4) with $k = n = WH$, it is a bijection
from the boards onto $[0, (WH)!)$, and it is given only up to 20 cells, since
$20! < 2^{64} < 21!$. Burns et al. [[3]](Z-bibliography.md#ref-3) used it as a
work distribution modulo $p$, and, as Jinnai and Fukunaga
[[5]](Z-bibliography.md#ref-5) show, modulo $p$ it then depends on very few
tiles ([Proposition 9.14](09-hash-distributed.md#prop-9-14)); it is kept as
the poor distribution to compare against.

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
distributed with the code of HDA$^*$ of Jinnai and Fukunaga. In the
examples, `8Puzzle` solves the instance of figure 3.28 of
[[13]](Z-bibliography.md#ref-13), or the first state found at a given
distance from the goal, and `15Puzzle` an instance of Korf, sequentially
or with HDA$^*$.

## 14.3 The sequence score of Nilsson {#sec-14-3}

`SequenceScore< W , H >` (`SequenceScore.h:54`) is the same puzzle with
the heuristic of Nilsson [[8]](Z-bibliography.md#ref-8), $h = M + W S$
(i.e., $M + 3S$ on the 8-puzzle), where $M$ is the Manhattan distance
and $S$ the *sequence score*. Let $\mathrm{next}(c)$ be the cell after
$c$ in the *ring*, which is, on the $3 \times 3$ board, the cycle of the
eight cells of the border in clockwise order (the centre having no next
cell), and, on the other boards, the cells $0, 1, \ldots, WH - 1$ in
row-major order (the last cell having no next cell). The *follower*
$\mathrm{fol}(t)$ of a tile $t$ is the tile that the goal puts in the
next cell of the goal cell of $t$ or, when that cell has no next cell,
in the cell 0 on the $3 \times 3$ board and in the cell 1 on the
others. The centre $c_0$ is the cell 4 on the $3 \times 3$ board and the
cell 0 on the others. Then
`heuristic( s )` is Algorithm 14.4 (`SequenceScore.h:110`).

**Algorithm 14.4** (The sequence score).

```text
1  procedure SEQUENCESCORE(σ)
2     M ← 0;  S ← 0
3     if σ(c₀) ≠ the goal tile of c₀ then
4        S ← 1
5     for each cell c with a tile t = σ(c) ≠ 0 do
6        M ← M + δ_t(c)
7        if next(c) exists and σ(next(c)) ≠ fol(t) then
8           S ← S + 2
9     return M + W · S
```

It takes $O(WH)$ time per state, and the successors of `SequenceScore`
are those of the puzzle with this value computed from scratch for each
of them.

<a id="prop-14-10"></a>**Proposition 14.10** (The sequence score is
not admissible). On the 8-puzzle with the goal of Korf, the state one
move away from the goal in which tile 1 has moved into the first cell
has $h = 13$, while $h^* = 1$.

*Proof.* Only tile 1 is out of place, at distance 1, and $M = 1$. Since
the centre holds tile 4, its goal tile, line 3 adds 0. Along the
border ([Figure 14.5](#fig-14-5)), tile 1 in the first cell is followed
by the blank instead of tile 2, and tile 3 in the cell 3 is followed (in
the first cell) by tile 1 instead of the blank, which follows it in the
goal; each other tile of the border is followed by its goal follower,
and the blank adds nothing (line 5). Hence $S = 4$ and $h = 1 + 3 \cdot
4 = 13 > 1$. $\square$

![Figure 14.5. The sequence score on the 8-puzzle. Left: the goal,
with the clockwise order of the border that defines the follower of each
tile. Right: the state one move away of
[Proposition 14.10](#prop-14-10), where the two pairs in rust break the
order.](../figures/domains-sequence.svg){#fig-14-5}

Since the heuristic may overestimate, A$^*$ with it is not guaranteed (cf.
[Theorem 2.9](02-best-first.md#thm-2-9)) to return an optimal path, and
it does not: on the instance of [Figure 14.1](#fig-14-1) (optimal cost
26, with $M = 18$ and $h = 69$), the example `8Puzzle -h seq` returns a
path of 28 moves after fewer than a hundred expansions (a run, not a
bound; the test `domains` only checks the overestimate of
[Proposition 14.10](#prop-14-10)). This domain is there to show this
trade, which weighted A$^*$ makes with a bound ([Chapter
6](06-bounded-suboptimal.md#ch-6)). It forwards no refinement of
`SlidingTile`.

## 14.4 Pathfinding on grids {#sec-14-4}

We call *grid* a rectangular map of $W \times H$ cells, each "blocked" or free
with a cost $\kappa(u) \in \{1, \ldots, 255\}$ (`GridMap`, `Grid.h:73`, where
`cost( x , y )` is 0 for a blocked cell or one outside the map), on which the
agent moves from cell to cell. It is the benchmark of the video games and of
robotics [[15]](Z-bibliography.md#ref-15), and many problems reduce to it (the
multiple sequence alignment of [Chapter 15](15-msa-instances.md#ch-15) is a
shortest path on a lattice of higher dimension). In the library,
`Grid< Moves >` (`Grid.h:259`) models pathfinding with 4 moves (up, left,
right, down) or 8 (the four diagonals too), from the cell $(s_x, s_y)$ to the
cell $(g_x, g_y)$, and its constructor throws if either is blocked. A state is
the index $yW + x$ of the cell $(x, y)$, which is also the perfect hash. A
move *leaves* a cell at its cost, times $\sqrt2$ for a diagonal, i.e.,
$c(u, v) = \kappa(u)$ or $\sqrt2\, \kappa(u)$ (written cost and sqrt2 in
Algorithm 14.5), and a diagonal move is allowed only if both the cells it cuts
are free ([Figure 14.6](#fig-14-6)); hence, the costs are integers with 4
moves and real numbers with 8. The successors are generated by a
`GridAdjacency< Moves >` (`Grid.h:175`) with the std and Boost data
structures, and by the grid graph of the default ones, which generates the
same moves at the same costs, in the same order with 4 moves and by increasing
index of the cell with 8 (an order that changes only the ties of OPEN);
Algorithm 14.5 is the former.

**Algorithm 14.5** (The moves on a grid).

```text
1  procedure GRIDSUCCESSORS(x, y, f)
2     D ← (0,−1), (−1,0), (1,0), (0,1), (−1,−1), (1,−1), (−1,1), (1,1)
3     for each of the first Moves pairs (δx, δy) of D do
4        (x', y') ← (x + δx, y + δy)
5        if (x', y') is blocked or outside the map then
6           continue
7        if δx ≠ 0 and δy ≠ 0 and ((x', y) or (x, y') is blocked) then
8           continue               ▷ no corner cutting
9        c ← cost(x, y) if δx · δy = 0 else sqrt2 · cost(x, y)
10       f((x', y'), c, h(x', y'))
```

An expansion thus costs $O(\mathrm{Moves})$, the heuristic of each
successor being a formula of its coordinates (below).

![Figure 14.6. (a) The moves from a cell $u$ of cost 2 with the cell
above it blocked: the straight moves cost 2, the diagonal ones
$2\sqrt2$, and the two diagonals that would cut the corner of the blocked
cell are not moves. (b) The octile distance between two cells: as many
diagonal moves as the smaller difference of the coordinates, then
straight moves.](../figures/domains-grid-moves.svg){#fig-14-6}

Note that the cost of an edge depends on the cell left: $c(u, v) =
\kappa(u)$ and $c(v, u) = \kappa(v)$ differ on a map with several costs.
The graph is therefore directed, and the problem is reversible (in the
sense of [§1.2](01-introduction.md#sec-1-2)) only on the maps of uniform
cost; this is why a grid is not an `OperatorDomain`, while its moves are
symmetric.

### The heuristics

The heuristic is the distance on the empty grid times the least cost
$\kappa_{\min}$ of a free cell (`GridMap::min_cost()`, 1 if no cell is
free): the Manhattan distance $d_4(u, t) = \Delta x + \Delta y$ with 4
moves, and the octile distance

$$
d_8(u, t) = \max(\Delta x, \Delta y) + (\sqrt2 - 1) \min(\Delta x,
\Delta y)
$$

with 8 moves, $\Delta x$ and $\Delta y$ being the absolute differences of
the coordinates of $u$ and of the goal $t$ (`Grid.h:443`). Our proof of
consistency rests on $d_4$ and $d_8$ being the costs of the shortest
paths on the grid with no blocked cell and unit cell costs; for $d_4$
this is immediate (each move changes $\Delta x + \Delta y$ by at most 1,
and the straight path attains it), and for $d_8$ it is the following
lemma.

<a id="lem-14-11"></a>**Lemma 14.11** (The octile distance). On the
infinite grid with 8 moves, straight moves of cost 1 and diagonal moves
of cost $\sqrt2$, the least cost of a path between two cells whose
coordinates differ by $\Delta x$ and $\Delta y$ is $d_8$.

*Proof.* A path with $k$ diagonal and $l$ straight moves changes the sum
of the coordinates by at most $2k + l$ and each coordinate by at most $k
+ l$; hence $2k + l \ge \Delta x + \Delta y$ and $k + l \ge \max(\Delta x,
\Delta y)$. Its cost is
$\sqrt2\, k + l = (\sqrt2 - 1)(2k + l) + (2 - \sqrt2)(k + l) \ge (\sqrt2
- 1)(\Delta x + \Delta y) + (2 - \sqrt2) \max(\Delta x, \Delta y)$,
since both coefficients are positive, and this equals $d_8$ since
$\Delta x + \Delta y = \max + \min$. The path of $\min(\Delta x, \Delta
y)$ diagonal moves followed by $|\Delta x - \Delta y|$ straight ones
attains the bound. $\square$

<a id="prop-14-12"></a>**Proposition 14.12** (The grid heuristics are
consistent). On a grid with 4 moves, $h = \kappa_{\min}\, d_4$ is
consistent; on a grid with 8 moves, $h = \kappa_{\min}\, d_8$ is
consistent.

*Proof.* Both $d_4$ and $d_8$ are shortest-path distances on the empty
grid ([Lemma 14.11](#lem-14-11) for $d_8$), and therefore they satisfy the
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
[Proposition 14.12](#prop-14-12) holds in exact arithmetic: two paths to
the same cell whose exact costs are equal may differ by an ulp, and a
search may then (rarely) reopen a cell, which costs time but not
optimality (the search of [Chapter 2](02-best-first.md#ch-2) reopens
whenever it finds a smaller $g$). As a check, the test `grid` compares
A$^*$ with the algorithm of Dijkstra, written apart from the domain, on
random maps of unit costs and of costs 1 to 5, with 4 and 8 moves, and
the test `domains` checks the consistency on every edge of three random
maps of costs 1 to 5.

<a id="prop-14-13"></a>**Proposition 14.13** (The moves of an open
grid). On a grid of $W \times H$ free cells, there are $2(2WH - W - H)$
moves with 4 moves, and $2(2WH - W - H) + 4(W - 1)(H - 1)$ with 8.

*Proof.* The straight moves are the two orientations of the $H(W - 1)$
horizontal and of the $W(H - 1)$ vertical pairs of adjacent cells. With
no blocked cell no diagonal is cut, and each of the $(W - 1)(H - 1)$
squares of $2 \times 2$ cells has two diagonals, each in two
orientations. $\square$

On $5 \times 5$ cells this gives 80 and 144 moves, i.e., the mean
branching factors $3.2$ and $5.76$, which the test `domains` checks.

### The refinements and the maps

The features of Zobrist hashing are the column and the row of the cell
($W + H$ features), the abstract features of AZH the column and the row
of the block of $B \times B$ cells containing it ($B$ = `block()`, 100
by default), and the abstract state of AHDA$^*$ that block, the two
distributions coinciding on this domain as in section 4.1 of
[[5]](Z-bibliography.md#ref-5). In turn, the nblocks of PBNF are the
blocks of $S \times S$ cells ($S$ = `nblock_side()`, 20 by default), and
the abstract successors of an nblock are its 4 (with 4 moves) or 8
(with 8 moves) neighbouring nblocks inside the map ([Figure
14.7](#fig-14-7)).

<a id="prop-14-14"></a>**Proposition 14.14** (The nblocks of the
grids). The nblocks of `Grid` satisfy [Definition
11.1](11-safepbnf.md#def-11-1).

*Proof.* (i) A move changes each coordinate by at most 1, and hence the
coordinates of the block, $\lfloor x/S \rfloor$ and $\lfloor y/S
\rfloor$, each by at most 1; with 4 moves only one coordinate changes.
Thus the nblock of the successor is that of the parent or one of its
neighbours of the right kind, both inside the map. (ii) Being a
neighbour is symmetric and irreflexive. $\square$

A grid models neither `IncrementalFeatureDomain` (a move changes one or
two of its features, but the domain does not enumerate them) nor
`OperatorDomain` (see above); its pairwise heuristic for PA$^*$SE is
the heuristic between any two cells ([Proposition 14.18](#prop-14-18)).

![Figure 14.7. (a) The blocks of AZH and AHDA$^*$ on a grid, coloured
by owner: a move changes the owner only across a thick line. (b) The
nblocks of Safe PBNF: the duplicate detection scope of the nblock $b$ is
the set of its neighbouring nblocks, which no other worker may use while
$b$ is acquired ([Chapter 11](11-safepbnf.md#ch-11)).
](../figures/domains-grid-blocks.svg){#fig-14-7}

A map is read from a file of the Moving AI benchmarks of Sturtevant
[[15]](Z-bibliography.md#ref-15) by `GridMap::read_movingai` (a header
`type octile`, `height`, `width`, `map`, then the rows, in which `.`, `G` and
`S` are free cells of cost 1 and any other character is blocked). Otherwise,
it is drawn at random by `GridMap::random( W , H , p , seed , c_max )`
(`src/Grid.cpp:53`), which blocks each cell with probability $p$ and gives
each free cell the cost 1 (if $c_{\max} = 1$) or a cost drawn uniformly in
$\{1, \ldots, c_{\max}\}$, row by row. Then it makes the bottom-left and the
bottom-right corners free with cost 1, since they are the start and the goal
of the instances of Burns et al. [[3]](Z-bibliography.md#ref-3). The draws are
made from the words of `std::mt19937_64` ([Lemma
1.9](01-introduction.md#lem-1-9)), and hence a seed gives the same map on all
platforms. In the examples, `FindPath` finds and draws a path on such a map.

## 14.5 Explicit graphs and the map of Romania {#sec-14-5}

`ExplicitGraph< Cost >` (`ExplicitGraph.h:138`) is a weighted graph given
by its arcs, with nodes $0, \ldots, n - 1$: `add_edge( u , v , c )` adds
an edge (or an arc, with `directed`), refusing a negative cost,
`set_heuristic( u , h )` gives $h$ node by node (0 by default),
`set_source( s )` changes the source, and the search goes from the
source to the target. It serves the small graphs of the tests (e.g., the
graph of [Figure 1.3](01-introduction.md#fig-1-3)) and the map of
Romania; since $h$ is given, its admissibility and consistency are the
responsibility of whoever gives it. Its arcs are kept in an
`ExplicitAdjacency` of the chosen set of data structures
([Chapter 3](03-data-structures.md#ch-3)), which enumerates the arcs that
leave a node in the order of their insertion with the std and the default
sets, but not with that of Boost; the default set also drops a loop $(u,
u)$, which never shortens a path. Its single feature is the node, and
its perfect hash the node itself.

`Romania` (`Romania.h:61`) is the road map of figure 3.2 of Russell and
Norvig [[13]](Z-bibliography.md#ref-13), with 20 cities, 23 roads whose
cost is their length in km, and real costs. Towards Bucharest the
heuristic is the straight-line distance of figure 3.22 of the book;
towards any other city $t$ it is the Euclidean distance $e(u, t)$
between the coordinates of the cities in the code that accompanies the
book (aima-python), both written into the graph by the constructor
(`Romania.h:83`).

<a id="prop-14-15"></a>**Proposition 14.15** (The heuristics of Romania
are consistent). On the map of Romania, the straight-line distance to
Bucharest of the book and the Euclidean distance between the
coordinates of aima-python, towards any city, are consistent.

*Proof.* For the Euclidean distance $e$, the triangle inequality gives
$e(u, t) \le e(u, v) + e(v, t)$. Also, $e(u, v)$ is at most the length
of the road $(u, v)$ for each of the 23 roads (the largest ratio of
$e(u, v)$ to the length of the road is about 0.996, on the road from
Rimnicu Vilcea to Pitesti, on the data of `Romania.h`); hence $h(u) \le
c(u, v) + h(v)$, and $h(t) = 0$. However, the table of the book is not
derived from the coordinates, and its consistency is checked edge by
edge: for each of the 23 roads $(u, v)$, $|h(u) - h(v)| \le c(u, v)$,
the least slack being 1 km, on the road from Pitesti to Bucharest ($100
- 0 \le 101$), and the next one 4 km, from Rimnicu Vilcea to Pitesti
($193 - 100 \le 97$). Since the roads are two-way, the inequality holds
in both directions. $\square$

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

`NQueens` (`NQueens.h:67`) places $1 \le n \le 32$ queens on an $n
\times n$ board in the incremental formulation of section 3.2.1 of the
book. A state is the number $k$ of queens placed and the row $\mathrm{row}(c)$
of the queen of each column $c < k$ (a record of 33 bytes), and a move,
Algorithm 14.6 (`NQueens.h:124`), adds a queen in the column $k$, in a
row that no queen already on the board attacks, at cost 1.

**Algorithm 14.6** (The successors of a placement).

```text
1  procedure QUEENSUCCESSORS(k, row, h, f)
2     for r ← 0 to n − 1 do
3        if r ≠ row(c) and |r − row(c)| ≠ k − c for each column c < k then
4           f((k + 1, row with row(k) = r), 1, h − 1)
```

The test of line 3 excludes the rows and the two diagonals of the queens
placed; hence an expansion takes $O(nk)$ time. The state space is a tree
whose paths from the root have length at most $n$, since each state has
a single parent (the placement without its last queen), and the goals
are the states with $k = n$; a state of the tree may have no successor,
a dead end ([Figure 14.9](#fig-14-9)). The features of Zobrist hashing
are the pairs (column, row) of the queens placed, and the domain has no
perfect hash.

<a id="prop-14-16"></a>**Proposition 14.16** (The heuristic of the
queens). The heuristic $h = n - k$, the number of the empty columns, is
consistent, with $h(u) = c(u, v) + h(v)$ on each edge, and it equals
$h^*$ on the states from which a goal is reachable.

*Proof.* Each edge adds one queen at cost 1, so that $h$ drops by
exactly 1, and $h = 0$ at a goal. A state with $k$ queens from which a
goal is reachable needs exactly $n - k$ more moves to reach it, each
move adding one queen. $\square$

Hence $f = n$ on all the states. The ties of $f$ are broken by the larger
$g$ ([Chapter 3](03-data-structures.md#ch-3)), and A$^*$ then expands the
states in a depth-first order, backtracking at the dead ends; the
counters of the example `NQueens` measure how many it meets. The
8-queens problem has 92 solutions, and its tree 2057 states (the root
included), two counts that the test `domains` obtains by enumerating the
tree; for $n = 2$ and $n = 3$ the problem has no solution, and A$^*$
reports it after exhausting the tree.

![Figure 14.9. The tree of the 4-queens problem, each state labelled
by the rows of the queens placed, column by column. The dead ends (rust)
have no free row in the next column; the two goals are the two solutions,
and all the states have $f = 4$.](../figures/domains-queens.svg){#fig-14-9}

### The vacuum world

`VacuumWorld` (`VacuumWorld.h:68`) is the vacuum world of section 3.2.1
of the book on a row of $1 \le n \le 32$ cells instead of 2: a robot
moves `Left` or `Right`, or `Suck`s the dirt of its cell, each action
costing 1, and the goal is a clean row. A state packs the set $D$ of the
dirty cells in the bits 0 to 31 of a word and the cell $r$ of the robot
in the bits 32 to 37, and the word is also the perfect hash; Algorithm
14.7 is its expansion (`VacuumWorld.h:128`).

**Algorithm 14.7** (The successors in the vacuum world).

```text
1  procedure VACUUMSUCCESSORS(r, D, h, f)
2     if r ∈ D then
3        f((r, D ∖ {r}), 1, h − 1)      ▷ Suck
4     if r > 0 then
5        f((r − 1, D), 1, h)            ▷ Left
6     if r + 1 < n then
7        f((r + 1, D), 1, h)            ▷ Right
```

<a id="prop-14-17"></a>**Proposition 14.17** (The vacuum world). From a
start with all the cells dirty, all the $n 2^n$ states are reachable;
the heuristic $h = |D|$ is consistent; and the nblocks of PBNF, i.e.,
the cells of the robot with the neighbouring cells as abstract
successors, satisfy [Definition 11.1](11-safepbnf.md#def-11-1).

*Proof.* A state $(r, D)$ is reached from $(r_0, \{0, \ldots, n-1\})$ by
walking along the row and sucking every cell not in $D$, and then by
walking to $r$. A `Suck` lowers $h$ by 1 at cost 1 and a move leaves it
unchanged, so that $h(u) \le c(u, v) + h(v)$ on each edge, and $h = 0$
at a goal. Finally, an action keeps the cell of the robot or moves it to
a neighbour, and being neighbours is symmetric and irreflexive. $\square$

The heuristic is also "weak", since it ignores the moves needed to reach
the dirt, which usually cost about as much as the cleaning. Furthermore,
the domain models the other abstractions of the parallel searches, all
of which keep the dirt and forget the robot. Thus, the features of
Zobrist hashing are the dirty cells $i < n$ and the cell $n + r$ of the
robot, the abstract features of AZH are the dirty cells (the $n$
features of the robot being mapped to one abstract feature, which all the
states have), and the abstract state of AHDA$^*$ is $D$; hence, a move
keeps the owner of a state under these two distributions, and only a
`Suck` may change it.

The pairwise heuristic of PA$^*$SE between $s$ and $t$,
`pairwise_heuristic( s , t )`, is the number of the cells dirty in $s$
and clean in $t$. We now prove the properties (H1) and (H3) that PA$^*$SE
requires (stated in [Chapter 8](08-shared-memory.md#ch-8), where
PA$^*$SE is given) for the pairwise heuristics
of the three domains that have one, $h(t, s)$ being the pairwise
heuristic and $d(t, s)$ the least cost of a path from $t$ to $s$.

<a id="prop-14-18"></a>**Proposition 14.18** (The pairwise
heuristics). For the sliding-tile puzzles with the Manhattan distance,
the grids with $\kappa_{\min} d_4$ or $\kappa_{\min} d_8$ and the vacuum
world with $|D|$, the pairwise heuristic satisfies (H1) $h(t, s) \le
d(t, s)$ and (H3) $h(t) \le h(t, s) + h(s)$, for all states $t$ and
$s$.

*Proof.* On the puzzles, $h(t, s) = \sum_{k \ne 0} \mathrm{dist}(c_t(k),
c_s(k))$, with dist the Manhattan distance between cells; a move changes
the cell of one tile by one step, and hence $h(t, \cdot)$ by at most 1,
which gives (H1), and the triangle inequality of dist, tile by tile,
gives (H3), since $h(t)$ is the same sum with the goal in place of $s$.
On the grids, $h(t, s) = \kappa_{\min} d(t, s)$ for $d = d_4$ or $d_8$:
(H1) is the argument of [Proposition 14.12](#prop-14-12) with $s$ in
place of the goal, and (H3) is the triangle inequality of $d$. On the
vacuum world, $h(t, s) = |D_t \setminus D_s|$: each cell of $D_t
\setminus D_s$ needs a `Suck` on any path from $t$ to $s$, which gives
(H1), and $|D_t| \le |D_t \setminus D_s| + |D_s|$, which is (H3).
$\square$

![Figure 14.10. The state space of the vacuum world of two cells of
the book: each state is the row with the robot (R) and the dirt (dots),
the edges are the actions, and the goals (dark violet) are the clean rows. The
path returned by A$^*$ (violet) is Suck, Right,
Suck.](../figures/domains-vacuum.svg){#fig-14-10}

### The water jugs

`WaterJugs` (`WaterJugs.h:67`) is the problem of exercise 3.6(d) of the
book (the `PourProblem` of aima-python). Up to 4 jugs of capacities
$\gamma_j \le 65535$ units are empty at the start, a move fills a jug
from the faucet, empties it on the ground, or pours it into another jug
until either the first is empty or the second is full, and the goal is a
jug holding a given quantity $q$ ([Figure 14.11](#fig-14-11)). Every move
costs 1, and the moves that do not change the state are not generated
(Algorithm 14.8, `WaterJugs.h:133`, where $L$ is the vector of the
levels). In the state, the level $L_j$ of the jug $j$ is packed in the
bits $16j$ to $16j + 15$, and the word is also the perfect hash; the
features of Zobrist hashing are the pairs
(jug, level).

**Algorithm 14.8** (The successors of the water jugs).

```text
1  procedure JUGSUCCESSORS(L, f)
2     for each jug i do
3        MOVE(L with Lᵢ = γᵢ)                     ▷ fill i
4        MOVE(L with Lᵢ = 0)                      ▷ empty i
5        for each jug j ≠ i do
6           δ ← min(Lᵢ, γⱼ − Lⱼ)
7           MOVE(L with Lᵢ = Lᵢ − δ and Lⱼ = Lⱼ + δ)   ▷ pour i into j
8  procedure MOVE(L')
9     if L' ≠ L then
10       f(L', 1, h(L'))
```

An expansion generates at most $J^2 + J$ states with $J$ jugs, each with
the goal test of $O(J)$ that gives its heuristic, in $O(J^3)$ time.

<a id="prop-14-19"></a>**Proposition 14.19** (The water jugs). The
heuristic, 0 at a goal and 1 elsewhere, is consistent. A goal is
reachable only if $q$ is a multiple of the greatest common divisor
$\gamma$ of the capacities and $q \le \max_j \gamma_j$.

*Proof.* Along an edge of cost 1 the heuristic changes by at most 1, and
it is 0 at the goals. For the second claim, each level is a multiple of
$\gamma$ at the start, and lines 3, 4 and 7 keep it so: $\gamma_i$, 0,
and $\delta$, the minimum of two multiples of $\gamma$, are multiples of
$\gamma$; also, no level exceeds the capacity of its jug. $\square$

The heuristic is (perhaps) the best estimate that does not use the
number theory of the capacities; the converse of the condition is not
needed by the domain, which leaves an instance that violates it to be
proved unsolvable by the search. Note that the problem is not
reversible: emptying a jug cannot be undone in one move.

![Figure 14.11. Measuring 4 units with jugs of 5 and 3: the six moves
of an optimal solution, with the levels of the two jugs under each
state.](../figures/domains-jugs.svg){#fig-14-11}

### The missionaries and cannibals

`MissionariesCannibals` (`MissionariesCannibals.h:65`) is the problem of
exercise 3.9 of the book with $m \le 65535$ missionaries, $c \le 65535$
cannibals and a boat of $b \ge 1$ seats: everybody starts on the left
bank with the boat, a crossing carries from 1 to $b$ people, and on
neither bank may the missionaries be outnumbered by the cannibals if any
missionary is there. Every crossing costs 1. A state packs the numbers
$m_L$ and $c_L$ of the missionaries and the cannibals on the left bank
and the bank of the boat, and the word is the perfect hash; the features
of Zobrist hashing are $m_L$, $c_L$ and the bank. Algorithm 14.9 is the
expansion (`MissionariesCannibals.h:146`), where $(a_m, a_c)$ are the
people on the bank of the boat and *safe*$(x, y)$ means $x = 0$ or $x
\ge y$.

**Algorithm 14.9** (The crossings).

```text
1  procedure CROSSINGS(m_L, c_L, boat, f)
2     (a_m, a_c) ← the people on the bank of the boat
3     for d_m ← 0 to min(a_m, b) do
4        for d_c ← (0 if d_m > 0 else 1) to min(a_c, b − d_m) do
5           (m', c') ← the people left on the left bank after the crossing
6           if safe(m', c') and safe(m − m', c − c') then
7              f((m', c', the other bank), 1, h(m', c', the other bank))
```

The heuristic counts the crossings that would carry the $n = m_L + c_L$
people left on the left bank over the river with the boat always full,
`heuristic( s )` (`MissionariesCannibals.h:129`):

$$
h_L(n) =
\begin{cases}
0 & \text{if } n = 0, \\
1 & \text{if } 0 < n \le b, \\
2 \lceil (n - b)/(b - 1) \rceil + 1 & \text{if } n > b > 1, \\
2n - 1 & \text{if } n > b = 1,
\end{cases}
$$

with $h = h_L(n)$ if the boat is on the left bank and $h = h_L(n) + 1$
otherwise (for $n > 0$).

<a id="prop-14-20"></a>**Proposition 14.20** (Admissible, not
consistent). The heuristic of `MissionariesCannibals` is admissible. It
is not consistent: with 3 missionaries, 3 cannibals and a boat of 2
seats, the first crossing of an optimal solution lowers it by 3.

*Proof.* With the boat on the left, a solution alternates crossings to
the right, each carrying at most $b$ people, and crossings back, each
bringing at least 1. After $k$ round trips at most $k(b - 1)$ people
have left, and the last crossing carries at most $b$; hence $n \le k(b -
1) + b$, i.e., $k \ge \lceil (n - b)/(b - 1) \rceil$ for $n > b$ and $b
> 1$, and the solution has at least $2k + 1$ crossings; for $0 < n \le
b$ it has at least 1. For $b = 1$ no round trip makes progress, no
solution exists when $n > 1$, and $h^* = +\infty$ bounds any value. With
the boat on the right and $n > 0$, the first crossing brings back $j \ge
1$ people, and the rest needs at least $h_L(n + j) \ge h_L(n)$
crossings, $h_L$ being nondecreasing; hence $h^* \ge 1 + h_L(n)$, the
value of $h$. For the second claim, the start has $n = 6$ and $h = 2
\lceil 4/1 \rceil + 1 = 9$; carrying a missionary and a cannibal (the
first crossing of the solution of [Figure 14.12](#fig-14-12)) leaves $n
= 4$ with the boat on the right, $h = 2 \lceil 2/1 \rceil + 1 + 1 = 6$,
and $9 > 1 + 6$. $\square$

The proof shows where the slack may be: on the right bank the heuristic
adds 1 to $h_L(n)$, while the argument gives $1 + h_L(n + 1)$, which is
larger by 2 for some $n \ge b$ (for each $n \ge b$ when $b = 2$). With
3, 3 and 2 the answer is 11 crossings, and with 4, 4 and 2 there is no
solution, as the test `domains` checks by exhausting the state space.
A$^*$ is optimal with an admissible heuristic, reopening the closed
states when needed ([Theorem 2.9](02-best-first.md#thm-2-9)); the
searches that assume consistency are not. For instance, the search with
delayed duplicate detection of [Chapter 4](04-memory-bounded.md#ch-4)
refuses a successor of a smaller $f$ than its parent with an exception
(`DDDAStar.h:360`), and it does so on this domain.

![Figure 14.12. The 11 crossings of an optimal solution with 3
missionaries, 3 cannibals and a boat of 2 seats; a state is the numbers
of missionaries and cannibals on the left bank and the bank of the boat.
The heuristic drops by 3 on the first crossing (rust), which costs 1
([Proposition 14.20](#prop-14-20)).
](../figures/domains-missionaries.svg){#fig-14-12}

## 14.7 The traveling salesperson problem {#sec-14-7}

`TravelingSalesman` (`TravelingSalesman.h:76`) finds a tour of least
cost through $2 \le n \le 32$ cities, starting and ending at the city 0,
given the matrix $d$ of the distances (`random_cities( n , seed )` draws
$n$ points of the unit square and gives their Euclidean distances). A
state is the set $S$ of the visited cities and the current one $c$,
packed as $S$ in the bits 0 to 31, $c$ in the bits 32 to 36 and, in the
bit 37, whether the tour has gone back to the city 0; the word is also
the perfect hash, and the features of Zobrist hashing are the visited
cities, the current one and the closing of the tour. The start is $S =
\{0\}$, $c = 0$, and Algorithm 14.10 is the expansion
(`TravelingSalesman.h:185`).

**Algorithm 14.10** (The successors of a partial tour).

```text
1  procedure TOURSUCCESSORS(S, c, f)
2     if S = {0, …, n − 1} then
3        f(the closed tour, d(c, 0), 0)     ▷ the goal
4        return
5     for each city k ∉ S, by increasing k, do
6        f((S ∪ {k}, k), d(c, k), MST((V ∖ (S ∪ {k})) ∪ {k, 0}))
```

The heuristic of a state $u$ with current city $c$ and set $U = V
\setminus S$ of unvisited cities ($V$ being the set of all the cities)
is the cost $\mathrm{MST}(U \cup \{c, 0\})$ of a minimum spanning tree of
the unvisited cities together with the current and the first one
(exercise 3.30 of [[13]](Z-bibliography.md#ref-13)), and 0 at the goal;
it is computed by the algorithm of Prim on the complete graph of these
$m$ cities, Algorithm 14.11 (`TravelingSalesman.h:148`), in $O(m^2)$
time.

**Algorithm 14.11** (The algorithm of Prim).

```text
1  procedure MST(u₁, …, u_m)
2     key(u₁) ← 0;  key(u_k) ← ∞ for k > 1;  T ← ∅;  w ← 0
3     for it ← 1 to m do
4        x ← the city not in T of least key
5        T ← T ∪ {x};  w ← w + key(x)
6        for each city y ∉ T do
7           key(y) ← min(key(y), d(x, y))
8     return w
```

<a id="lem-14-21"></a>**Lemma 14.21** (The algorithm of Prim). If $d$ is
symmetric, Algorithm 14.11 returns the cost of a minimum spanning tree
of the complete graph on $u_1, \ldots, u_m$ with the costs $d$.

*Proof.* Let the city $x$ of line 4 be joined to the city of $T$ for
which its key was last lowered at line 7 (none for $u_1$); then $w$ is
the cost of the tree $F$ so built on $T$, since key$(x)$ is the least
cost of an edge from $x$ to $T$. We show by induction on $|T|$ that $F$
is contained in a minimum spanning tree. Let $F^*$ be a minimum spanning
tree containing $F$ before $x$ joins, and let $e = (x, y)$, $y \in T$,
be the edge added. If $e \notin F^*$, the path of $F^*$ from $x$ to $y$
crosses the cut between $T$ and the rest by an edge $e'$, of cost at
least key$(x)$, since $x$ has the least key among the cities outside $T$
and the key of a city is its least edge to $T$; replacing $e'$ by $e$
gives a spanning tree of cost at most that of $F^*$, which contains $F
\cup \{e\}$. After $m$ steps $F$ spans all the cities, and it is a
minimum spanning tree. $\square$

Its state space is small, as an exhaustive count shows.

<a id="prop-14-22"></a>**Proposition 14.22** (The states of the tour).
The states reachable from the start are $(n - 1)\, 2^{n-2} + 2$: the
start, the closed tour, and the pairs of a nonempty set $S' \subseteq
\{1, \ldots, n - 1\}$ and of a current city $c \in S'$, the visited set
being $S' \cup \{0\}$.

*Proof.* Line 6 adds to $S$ the city it moves to, so that every state
but the start and the closed tour has its current city in $S \setminus
\{0\}$; conversely, each pair $(S', c)$ is reached by visiting the
cities of $S' \setminus \{c\}$ in any order and then $c$. There are
$\sum_{j=1}^{n-1} \binom{n-1}{j} j = (n - 1)\, 2^{n-2}$ such pairs.
$\square$

One can compare this number with the $(n - 1)!/2$ tours of a
symmetric instance: this is the state space of the dynamic program of
Held and Karp [[14]](Z-bibliography.md#ref-14), and A$^*$ on it explores
a part of that table, the part selected by the heuristic.

<a id="prop-14-23"></a>**Proposition 14.23** (The MST heuristic is
consistent). If the distances are symmetric and nonnegative, the
heuristic of `TravelingSalesman` is consistent, and hence admissible.

*Proof.* Let $u$ have current city $c$ and unvisited set $U \ne
\emptyset$, and let $v$ be the successor that moves to $k \in U$; thus,
$v$ has current city $k$ and unvisited set $U \setminus \{k\}$, and
$h(v) = \mathrm{MST}(U \cup \{0\})$. A minimum spanning tree of $U \cup
\{0\}$ together with the edge $(c, k)$ is a spanning tree of $U \cup
\{c, 0\}$, of cost $h(v) + d(c, k)$ (if $c = 0$ the edge is not needed);
hence $h(u) \le d(c, k) + h(v) = c(u, v) + h(v)$ ([Figure
14.13](#fig-14-13)). When $U = \emptyset$ and $c \ne 0$, the only
successor is the goal, at cost $d(c, 0)$, and $h(u) = \mathrm{MST}(\{c,
0\}) = d(c, 0) = c(u, v) + 0$; at the goal $h = 0$. The values of $h$
are those of a minimum spanning tree by [Lemma 14.21](#lem-14-21).
$\square$

We remark that the proof needs the symmetry, since a spanning tree is a
set of undirected edges; on an asymmetric matrix Algorithm 14.11 still
returns a number (it reads $d$ from the tree towards the other cities),
but the proposition no longer applies. Admissibility is also seen
directly: the rest of a tour from $u$ is a path from $c$ through all the
cities of $U$ to 0, i.e., a spanning tree of $U \cup \{c, 0\}$, and it
costs at least the minimum one. An expansion costs $O(n^3)$, i.e., one
MST of $O(n^2)$ for each of the at most $n - 1$ successors; the test
`domains` checks the consistency on every edge of the state spaces of 5
to 9 random cities, and the number of states of
[Proposition 14.22](#prop-14-22).

![Figure 14.13. The MST heuristic on seven cities. (a) At city 5,
after the cities 0 and 2 (violet path), $h$ is the cost of the tree
(dashed) on the unvisited cities and on 5 and 0. (b) After the move to 3,
the tree on the cities left and on 3 and 0, together with the edge
$(5, 3)$, spans the set of (a); hence $h$ drops by at most the cost of
the move ([Proposition 14.23](#prop-14-23)).
](../figures/domains-tsp.svg){#fig-14-13}

The header also gives the dynamic program of Held and Karp
[[14]](Z-bibliography.md#ref-14), as `held_karp()`
(`TravelingSalesman.h:258`), by which the test `domains` and the example
`TravelingSalesman` check the cost found by A$^*$ (the example up to 20
cities). With $C(S, j)$ the least cost of a path that starts at the city
0, visits exactly the cities of the set $S \ni 0$ and ends at $j \in S$,
it rests on the recurrence $C(S \cup \{k\}, k) = \min_{j \in S} \{C(S,
j) + d(j, k)\}$ for $k \notin S$, and the optimal tour costs $\min_{j
\ne 0} \{C(V, j) + d(j, 0)\}$.

**Algorithm 14.12** (The dynamic program of Held and Karp).

```text
1  procedure HELDKARP(d, n)
2     C({0}, 0) ← 0;  C(S, j) ← ∞ for each other pair
3     for each set S ∋ 0, by increasing value of its bits, do
4        for each j ∈ S with C(S, j) < ∞ do
5           for each k ∉ S do
6              C(S ∪ {k}, k) ← min(C(S ∪ {k}, k), C(S, j) + d(j, k))
7     return min over j ≠ 0 of C(V, j) + d(j, 0)
```

<a id="prop-14-24"></a>**Proposition 14.24** (Held and Karp). Algorithm
14.12 returns the cost of an optimal tour, in $O(2^n n^2)$ time with a
table of $2^n n$ entries.

*Proof.* The recurrence holds since the last edge of an optimal path
that ends at $k$ is some $(j, k)$, and the rest is a path that visits
$S$ and ends at $j$, which may be replaced by an optimal one. Line 3
visits each set after all its subsets, since adding a city raises the
value of the bits; hence, when $S$ is visited, every $C(S, j)$ has
received at line 6 all the terms of its recurrence, and by induction on
$|S|$ it is final. A tour is a path through $V$ ending at some $j \ne
0$, followed by $(j, 0)$, which gives line 7. The table has $2^n$ sets
times $n$ cities, and lines 4 to 6 take $O(n^2)$ per set. $\square$

The time and the memory are why the method serves to check the searches
on small instances only.

## 14.8 Slow expansions {#sec-14-8}

`SlowDomain< D >` (`SlowDomain.h:55`) wraps a domain `D` and adds to
each expansion (both `for_each_successor`, the one with an operator
included) a given number $w$ of iterations of the computation $i
\leftarrow 11 i \bmod 9999943$, from $i = 7$ (`SlowDomain.h:143`), which
keeps a core busy with no access to memory, as in the experiment of
section 3.1.1 of Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) on
the effect of slow expansions. It models the problems in which the
expansion dominates the cost of a search (planning with an expensive
heuristic, or motion planning with collision checks), the regime in
which PA$^*$SE of [Chapter 8](08-shared-memory.md#ch-8) is meant to pay.
The result of the loop is added to a thread-local variable, which is
visible outside the function, and therefore the compiler may not remove
the loop.

<a id="prop-14-25"></a>**Proposition 14.25** (A slow domain is the same
problem). `SlowDomain< D >` has the states, the initial state, the
goals, the edges with their costs, the heuristic, the hash and the
refinements of `D` but `IncrementalFeatureDomain`; hence every guarantee
of a search on `D` holds on it, and each expansion costs $\Theta(w)$
more.

*Proof.* Each method of the wrapper calls the method of the same name of
`D` with the same arguments, the expansions after the busy loop, which
changes no state of the search. $\square$

Since `IncrementalFeatureDomain` is not forwarded, the Zobrist keys of
the successors of a slow domain are computed from scratch; with an
expansion that is slow by construction, this cost is negligible. In the
benchmark, the option `--slow` gives $w$.

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

Table 14.2 collects the cost of an expansion of each domain, which
follows from the count of the loops of the algorithm cited (with $J$
jugs, $k$ queens placed, and the $O(1)$ of the operations on a word).

**Table 14.2.** The cost of an expansion, and of the heuristic of a
successor, in each domain.

| domain | algorithm | successors | expansion | $h$ of a successor |
|---|---|---|---|---|
| `SlidingTile` | 14.1, 14.2 | at most 4 | $O(\beta)$ | $O(1)$, incremental |
| `SequenceScore` | 14.2, 14.4 | at most 4 | $O(WH)$ | $O(WH)$ |
| `Grid< Moves >` | 14.5 | at most Moves | $O(\mathrm{Moves})$ | $O(1)$ |
| `ExplicitGraph` | | out-degree | $O(\text{out-degree})$ | $O(1)$, a lookup |
| `NQueens` | 14.6 | at most $n$ | $O(nk)$ | $O(1)$ |
| `VacuumWorld` | 14.7 | at most 3 | $O(1)$ | $O(1)$ |
| `WaterJugs` | 14.8 | at most $J^2 + J$ | $O(J^3)$ | $O(J)$, the goal test |
| `MissionariesCannibals` | 14.9 | $O(b^2)$ | $O(b^2)$ | $O(1)$ |
| `TravelingSalesman` | 14.10, 14.11 | at most $n - 1$ | $O(n^3)$ | $O(n^2)$ |

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
`missionaries`, `jugs`, `tsp`, `romania`). The puzzles read their
instances from `--file` (the instances chosen by `--ids`), the grids
from a Moving AI map or from `--random WxH:p:seed[:c]`, and the other
domains are sized by `--size` (e.g., `queens 28`, `missionaries
65535:65535:100`, `jugs 1:12:8:3`, `tsp 16` with the seeds of `--ids`).
It prints one line of CSV per run, with the cost, the counters of
[§1.3](01-introduction.md#sec-1-3), LB, CO and the seconds.

**Table 14.3.** The problems of the examples, solved by sequential
A$^*$.

--8<-- "tables/domains-examples.md"

Table 14.3 is written by `tables.py` from the table E of the campaign of
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
distance; the parity argument of [Proposition 14.2](#prop-14-2) and its
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
[[14]](Z-bibliography.md#ref-14). The domains of section 2.2 of the
textbook of Jinnai [[164]](Z-bibliography.md#ref-164) (the grid, the
puzzles and the traveling salesperson) are all here.

It is fair to remark that the literature of the parallel searches uses
these domains, and therefore they allow the searches of the library to
be compared with the published ones; however, they are written by hand
for one problem each, with hand-made abstractions, and the choice of the
abstraction from the problem itself is open ([Chapter
10](10-auto-distributions.md#ch-10)). Also, the representation of the
puzzles in one or two words is that of the fastest published codes
[[4]](Z-bibliography.md#ref-4), while the other domains use a word for
convenience more than for speed; their expansions are cheap
(Table 14.2), and on them the measures of the parallel searches are
dominated by the bookkeeping of the searches, which is the regime in
which those searches are hardest to scale.

## 14.12 Exercises {#sec-14-12}

**Exercise 14.1** (discussion, no test). Model the Rubik's cube as a
`Domain`: give the state and its packing in a word, the 18 moves (a
quarter turn either way and a half turn of each face) and their
operators, the features for Zobrist hashing, and a pair of refinements
that the cube can model; explain why the operator should also skip the
moves of the opposite face in one of the two orders.

**Exercise 14.2** (`ctest -R PBeFS.astar`). Show that the 8-puzzle has
$9!/2 = 181440$ reachable states and the 15-puzzle $16!/2$, using
[Proposition 14.2](#prop-14-2) and its converse; the test checks the first
number by a breadth-first search of the whole space from the goal.

**Exercise 14.3** (`ctest -R PBeFS.domains`). Using
[Proposition 14.13](#prop-14-13), compute the mean branching factor of a
grid of $5 \times 5$ free cells with 4 moves ($80/25 = 3.2$) and with 8
moves ($144/25 = 5.76$); then measure it with the counters of
`AStarBench`, as `generated / expanded` on a search that expands all the
cells (e.g., towards a goal walled off by blocked cells), and explain
the difference between the two numbers. The test checks the two counts
of moves.

**Exercise 14.4** (`ctest -R PBeFS.domains`). Prove
[Proposition 14.4](#prop-14-4) in detail for the 8-puzzle, i.e., that
the mean branching factor over its reachable states is exactly $8/3$,
and $5/3$ without the move back to the parent; the test checks the first
value on the whole state space, as $3$ times the number of the moves
equal to $8$ times the number of the states.

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

**Exercise 14.7** (`ctest -R PBeFS.domains`). Using
[Proposition 14.22](#prop-14-22) and [Proposition 14.23](#prop-14-23),
show that A$^*$ with the MST heuristic expands each state at most once in
exact arithmetic, and hence never expands more states than there are
entries in the table of Held and Karp. The test compares the cost of
A$^*$ with `held_karp()` on random instances and checks the consistency
of the heuristic on all the edges of their state spaces.

**Exercise 14.8** (`ctest -R PBeFS.domains`). Show that A$^*$ on
`NQueens`, with its default OPEN list (the bucket queue), expands the
states in the order of a depth-first search that tries the rows in
decreasing order (the successors are given by increasing row, and the
ties of $f$ and $g$ are broken in favor of the last pushed), and count
the states of the tree for
$n = 8$; the test checks that the tree has 2057 states and 92 goals, and
that the problem with 3 queens is unsolvable.

**Exercise 14.9** (`ctest -R PBeFS.graph`). Check [Proposition
14.15](#prop-14-15) on the data of `Romania.h`, and explain why the table
of the book cannot serve as the heuristic towards a city other than
Bucharest, while the coordinates can. The test checks the distances to
Bucharest from each city.
