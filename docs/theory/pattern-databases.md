# Pattern databases

We are concerned with the heuristics of the sliding-tile puzzles that store
the exact distances of an abstraction of the puzzle in a table, i.e., the
pattern databases of Culberson and Schaeffer [[21]](../../README.md#ref-21),
summed over a partition of the tiles as in the disjoint databases of Korf and
Felner [[22]](../../README.md#ref-22) and of Felner, Korf and Hanan
[[23]](../../README.md#ref-23). This page uses the notation of
[Best-first search and A*](best-first-search.md); it defines the abstraction
of a pattern under the two rules of the blank of the library, proves that the
sum of the tables of a partition and its maximum with the reflection about the
main diagonal are admissible, and when they are consistent, and then gives the
encodings of the entries, from one byte down to one bit, with the argument by
which the relative encodings recover the exact value. The pseudocode of the
construction and of the lookup is on the page of
[PatternDatabase.h](../../include/PBeFS/astar/PatternDatabase.h), and the last
section relates the classes of the library to the literature.

## The abstraction of a pattern

The board has $n = W H$ cells, and a state $s$ places the tiles
$1, \ldots, n - 1$ and the blank on distinct cells; a move exchanges the blank
with a tile in an adjacent cell, at cost 1. A *pattern* is a set
$\pi = \{t_1, \ldots, t_k\}$ of tiles, and the *placement* of $\pi$ in $s$ is
the vector $p = \psi_\pi(s) = (q_{t_1}, \ldots, q_{t_k})$ of the cells of its
tiles, the other tiles being "don't care". The placements are the
$k$-permutations of the $n$ cells, i.e., $N = n! / (n - k)!$ of them, and the
library ranks them in lexicographic order (cf. `TilePattern`): the digit of
the tile $t_i$ is the number $d_i$ of the cells below $p_i$ that the tiles
before it do not hold, and

$$
\mathrm{rank}(p) = \sum_{i=1}^{k} d_i \prod_{j=i+1}^{k} (n - j + 1)
$$

is a bijection onto $[0, N)$, computed by $k$ population counts and inverted
by $k$ divisions and $k$ selections of a bit. Two consecutive ranks differ
only in the cell of the last tile, which the max-compression below exploits.
Of course, the linear-time ranking of Myrvold and Ruskey
[[34]](../../README.md#ref-34) would do as well; the lexicographic one is
chosen because it keeps close placements in close entries, and its cost is
small against the miss of the cache that the lookup of a large table pays
anyway.

The *abstract space* of $\pi$ is the graph whose nodes are the placements and
whose edges are the moves of one tile of $\pi$ to an adjacent cell free of
the other tiles of $\pi$, each at cost 1. The library builds it under two
rules of the blank (cf. `BlankRule`).

- With the blank *ignored*, an abstract state is a placement and every such
  move is allowed, as if the blank were wherever needed. The table holds the
  distance $d_\pi(p)$ from $p$ to the goal placement in this graph, which is
  undirected and has unit costs.
- With the blank *minimized*, an abstract state is a placement together with
  the cell of the blank, a move of the blank into a cell that no tile of
  $\pi$ holds costs 0 (it moves a tile out of the pattern), and a move of a
  tile of $\pi$ costs 1. The blank then reaches at no cost every cell of its
  *region*, i.e., of its connected component among the cells free of $\pi$,
  and the table holds the least distance over the cells of the blank,
  $d_\pi(p) = \min_b d_\pi(p, b)$; this is the rule of the disjoint databases
  of [[22]](../../README.md#ref-22).

In both cases a move of the puzzle that moves a tile of $\pi$ is an edge of
cost 1 of the abstract space, and a move of any other tile changes no
placement (with the blank ignored) or is an edge of cost 0 (with the blank
minimized). Since the minimized space only forbids the abstract moves for
which the blank cannot reach the destination, $d_\pi$ is at least as large
with the blank minimized as with the blank ignored, and the tests check it on
every state of the 8-puzzle.

## Admissibility of the additive sum

Let $\pi_1, \ldots, \pi_m$ be disjoint patterns, and let

$$
h(s) = \sum_{j=1}^{m} d_{\pi_j}(\psi_{\pi_j}(s)) .
$$

Consider a path $P$ from $s$ to the goal, and let $c_j(P)$ be the number of
its moves that move a tile of $\pi_j$. The projection of $P$ on $\pi_j$ is a
path of the abstract space from $\psi_{\pi_j}(s)$ to the goal placement whose
cost is $c_j(P)$, since the moves of the other tiles cost nothing there;
hence $d_{\pi_j}(\psi_{\pi_j}(s)) \leq c_j(P)$. Each move moves exactly one
tile, which belongs to at most one pattern, so that
$\sum_j c_j(P) \leq c(P)$, and summing over $j$ gives $h(s) \leq c(P)$ for
every path, i.e., $h(s) \leq h^*(s)$. This is the argument of
[[22]](../../README.md#ref-22), which Yang et al.
[[28]](../../README.md#ref-28) cast for any abstraction whose costs split the
cost of each edge; we remark that it needs the patterns to be disjoint and the
abstract spaces to charge only the moves of their own tiles, which is why the
blank belongs to no pattern (the constructor of `TilePattern` refuses it). It
also gives $h = 0$ at the goal, whose placements are the goal placements.

## Consistency

With the blank ignored, the abstract space is undirected with unit costs, so
that $|d_\pi(p) - d_\pi(p')| \leq 1$ for adjacent placements $p$ and $p'$; a
move of the puzzle changes the placement of at most one pattern, by one
abstract edge, and hence $|h(s) - h(t)| \leq 1 = c(s, t)$ for each move
$(s, t)$, i.e., $h$ is consistent. Furthermore, the distance has the parity
of the Manhattan distance $m_\pi(p)$ of the tiles of the pattern, since each
abstract move changes $m_\pi$ by one unit, and two adjacent placements differ
by exactly one, never by zero; this is what the relative encodings below use.

With the blank minimized, the same move may change $d_\pi$ by more than one:
the minimum over the cells of the blank forgets in which region the blank
actually is, and the region that the true blank occupies after the move may
be far from the one that realized the minimum before it. The sum is then
admissible but not consistent in general; on the 8-puzzle with the partition
4-4 it is inconsistent along about 37000 of its moves (the test `pdb` prints
the count), while A* stays optimal by reopening the closed nodes (cf.
[Best-first search and A*](best-first-search.md) and
[[31]](../../README.md#ref-31)). The rule pays this with a larger heuristic:
on the 100 instances of Korf the mean of the partition 5-5-5 is 40.5 with the
blank minimized against 38.8 with the blank ignored (37.1 for the Manhattan
distance).

## The reflection about the main diagonal

Let the board be square, and let the goal have the blank on the main
diagonal, as the goal of Korf (the blank in the first cell) does. The
transposition $T$ of the cells maps adjacent cells to adjacent cells, so that
the map $\phi$ which moves the tile in the cell $c$ to the cell $T(c)$ and
renames the tile $t$ as $\sigma(t)$, the tile that the goal holds in the cell
$T(c^*_t)$, with $c^*_t$ the goal cell of $t$, is an automorphism of the
graph of the puzzle that fixes the goal. Hence $h^*(\phi(s)) = h^*(s)$, and
the sum of the same tables on $\phi(s)$ is also admissible for $s$, since
$h(\phi(s)) \leq h^*(\phi(s)) = h^*(s)$; the library takes

$$
h_\phi(s) = \max \{ h(s), h(\phi(s)) \} ,
$$

which is the reflection of [[22]](../../README.md#ref-22). The maximum of two
consistent heuristics is consistent, since $h(s) \leq 1 + h(t) \leq
1 + h_\phi(t)$ for each move $(s, t)$ and the same holds for
$h \circ \phi$, which is consistent because $\phi$ maps moves to moves; for
the same reason the maximum over several abstractions (cf.
[[30]](../../README.md#ref-30)) keeps the properties of its terms. The
reflection costs a second lookup per table, and it is used only when it is
valid (cf. `AdditivePDB::reflected()`), i.e., never on a goal whose blank is
off the diagonal.

The *dual lookup* of Zahavi et al. [[29]](../../README.md#ref-29), which
looks up the placement of the inverse permutation of the state, is another
admissible lookup of the same tables; it is not consistent, and its use pays
off with the bidirectional pathmax of IDA*, which the library does not have
(cf. the last section).

## The encodings of the entries

Let $e_\pi(p) = (d_\pi(p) - m_\pi(p)) / 2$ be the *excess* of the distance
over the Manhattan distance of the tiles of the pattern, which is a
nonnegative integer since $m_\pi$ is a lower bound of $d_\pi$ (each abstract
move moves one tile by one cell) and the two have the same parity (with both
rules the moves of cost 1 change $m_\pi$ by one unit, and those of cost 0 do
not change it). The Manhattan distance of the tiles is computed at the lookup
from a table of $k n$ bytes, and the entries store the excess, which is much
smaller than the distance; the encodings are the following, all in
[PatternDatabase.h](../../include/PBeFS/astar/PatternDatabase.h).

- `ByteEntries` stores $d_\pi$ in 8 bits per entry.
- `NibbleEntries` stores $e_\pi$ in 4 bits per entry, saturated at 15. A
  saturated entry gives a smaller value, which is still admissible, and
  `clamped()` counts them; none of the tables of the library saturates.
- `Mod3Entries` stores $d_\pi \bmod 3$, five residues in a byte as the digits
  of a number in base 3 ($3^5 = 243 \leq 256$), i.e., $8 / 5 = 1.6$ bits per
  entry, after Breyer and Korf [[25]](../../README.md#ref-25).
- `ParityBitEntries` stores $\lfloor d_\pi / 2 \rfloor \bmod 2$, one bit per
  entry: with the parity of $d_\pi$, which is that of $m_\pi$, it gives
  $d_\pi \bmod 4$, as in the 1-bit encoding of Clausecker and Reinefeld
  [[27]](../../README.md#ref-27).
- `MinNibbleEntries` stores, for each group of $c$ consecutive ranks ($c$ a
  power of 2), the least excess of the group in 4 bits, i.e., $4 / c$ bits per
  entry: the lossy max-compression of Felner et al.
  [[24]](../../README.md#ref-24), which here keeps the Manhattan distance of
  each placement exact.
- `SparseEntries`, only with the data structures of ADS (cf.
  [Structures.h](../../include/PBeFS/astar/Structures.h)), stores a bit per
  entry that tells whether its excess is nonzero, with rank and select, and
  the nonzero excesses in 4 bits each; it is smaller than `NibbleEntries`
  when fewer than about 64% of the excesses are nonzero.

The values of `ByteEntries`, `NibbleEntries` (when no entry saturates) and
`SparseEntries` are the distances themselves, so that they are admissible and
consistent exactly as above.

**The relative encodings.** The residues of `Mod3Entries` and
`ParityBitEntries` do not give the value of a placement, but they give the
change of the value along a move. With the blank ignored, two adjacent
placements have distances that differ by exactly one, so that the change
$\delta \in \{-1, +1\}$ is told apart by the residues modulo 3,

$$
\delta = \bigl( (r' - r + 1) \bmod 3 \bigr) - 1 ,
$$

with $r$ and $r'$ the residues before and after the move, and likewise modulo
4. A move of the puzzle moves the tile of one pattern only, so that the
heuristic of a successor is that of its parent plus the change of one table
(cf. `AdditivePDB::delta()`), and A* never needs the absolute value but for
the initial state. This one is found from the residues alone: from a
placement $p$ other than the goal there is an adjacent placement whose
distance is one less (the next one on a shortest path, the graph being
undirected), its residue is $(r - 1) \bmod 3$, and no other neighbour has
that residue, since the neighbours differ from $p$ by $\pm 1$ and
$+1 \not\equiv -1 \pmod 3$. Following such neighbours reaches the goal
placement in exactly $d_\pi(p)$ steps, which count the value (cf.
`PatternTable::exact()`); each step costs at most $4 k$ lookups. The argument
needs adjacent placements to differ by exactly one, i.e., the blank ignored,
and the heuristic to be the sum of the tables, since from the maximum of two
sums the change of one table cannot be told; this is why `SlidingTilePDB`
refuses a relative encoding with the reflection.

**The max-compression.** The value of `MinNibbleEntries` is
$m_\pi(p) + 2 \min_{p' \in G(p)} e_\pi(p') \leq m_\pi(p) + 2 e_\pi(p) =
d_\pi(p)$, with $G(p)$ the group of the rank of $p$; hence it is admissible,
and summing over a partition preserves the inequality. It is not consistent
in general, since two adjacent placements may fall in groups with very
different minima, and A* may then reopen nodes. Felner et al.
[[24]](../../README.md#ref-24) store the least distance of the group; storing
the least excess instead loses less, because the placements of a group differ
in the cell of the last tile only, and the Manhattan distance, which accounts
for most of their differences, is kept exact.

## The construction

The table of a pattern is filled by a breadth-first search of the abstract
space backward from the goal placement, by layers: at the layer $d$ every
entry of value $d$ is expanded and its unseen neighbours get $d + 1$, and the
search ends when a layer adds no entry (BUILD-IGNORED in the header). The
layer is scanned in blocks of ranks shared among the workers of a runtime,
and an entry is claimed by an atomic compare-and-swap, so that the result
does not depend on the number of the workers (the test checks 1 against 3).
With the blank minimized, the states $(p, b)$ take 2 bits each (unseen,
current, next, old), as in the two-bit breadth-first search of Korf
[[26]](../../README.md#ref-26), the cell of the blank being coded by its index
among the $n - k$ cells free of the pattern; a state is expanded once per
region, from the least cell of the region, and its successors are marked
together with all the cells of the region of the blank after the move, which
realizes the moves of cost 0 without a second queue (BUILD-MINIMIZED in the
header). The construction takes $O(N k)$ time with the blank ignored and
$O(N (n - k) k)$ with the blank minimized, and $N$ bytes (plus
$N (n - k) / 4$ bytes with the blank minimized), and the partitions are
built one after the other, each with all the workers.

The tables of the 24-puzzle take minutes, those of the partition 7-8 of the
15-puzzle about one; `PatternTable::save()` writes a table to a file with a
versioned header (the magic, the version, the board, the rule, the encoding
and its parameter, the tiles and the goal), the payload and a checksum of
everything, written to a temporary file and renamed, and `load()` rejects a
file that is truncated, corrupted, of another version or of another table.
With a cache directory (`PDBOptions::cache`), a table is loaded if it is
there, and otherwise it is encoded from the table in bytes of the same
pattern, if that one is there, or built; the tables built are saved.

## Sizes and the literature

The table of a pattern of $k$ tiles has $N = n! / (n - k)!$ entries: for the
partition 7-8 of the 15-puzzle [[22]](../../README.md#ref-22), 57657600 and
518918400 entries, i.e., 576 MB in bytes, 288 MB in nibbles, 115 MB modulo 3
and 72 MB in bits; for the partition 5-5-5 [[23]](../../README.md#ref-23),
three tables of 524160 entries, i.e., 0.8 MB in nibbles; for the partition
6-6-6-6 of the 24-puzzle [[22]](../../README.md#ref-22), four tables of
127512000 entries, i.e., 255 MB in nibbles and 64 MB in bits. The benchmark
`PDBBench` and the script `tools/bench/pdb.sh` give the size, the time of
the construction, the mean heuristic and the time of a lookup of each
encoding, and the expansions of A* on the 100 instances of Korf with each
heuristic.

Beyond the encodings above, the literature compresses the tables without
loss when the merged entries are known to differ by little
[[24]](../../README.md#ref-24), or by learning a smaller model of the
entries [[33]](../../README.md#ref-33), and builds the largest tables on
clusters, by partitioning the ranks among the processes as the blocks of the
library are partitioned among the threads [[32]](../../README.md#ref-32). The
searches
that use the tables best on the 24-puzzle are depth-first, i.e., IDA* with
the bidirectional pathmax of the inconsistent lookups
[[31]](../../README.md#ref-31), whose memory does not grow with the nodes
expanded; the library has no IDA*, and the states of `SlidingTile` hold at
most 16 cells, so that the tables of the 24-puzzle are built and looked up
but not searched.
