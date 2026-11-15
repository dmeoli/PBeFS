# 15. Multiple sequence alignment and the instances {#ch-15}

[Source: [`include/PBeFS/astar/SequenceAlignment.h`][SequenceAlignment.h]
(`SequenceAlignment`, `ScoreMatrix`, `pam250`),
[`include/PBeFS/astar/Distribution.h`][Distribution.h]
(`HyperplaneDistribution`), the concept `LevelDomain` of
[`include/PBeFS/astar/Domain.h`][Domain.h], the readers of the instance files
of the domains, and the instance sets of [`data/`][data]]

> **Status — beta.** The code that this chapter describes is in `main`
> but not yet in a release, and its interface may still change.

We are concerned in this chapter with two subjects that the benchmarks of the
parallel searches share. The first is the optimal alignment of several protein
sequences (the Multiple Sequence Alignment, MSA), formulated as a shortest path
on an $n$-dimensional lattice by Ikeda and Imai
[[41]](Z-bibliography.md#ref-41) and used, after Yoshizumi, Miura and Ishida
[[42]](Z-bibliography.md#ref-42), as the test of a search with a large
branching factor ($2^n - 1$ successors of a state) and nonuniform costs; on it
Kobayashi, Kishimoto and Watanabe [[43]](Z-bibliography.md#ref-43) evaluated
HDA$^*$ and proposed the hyperplane work distribution, which we describe here
together with the concept of a levelled graph on which it rests. The second is
the collection of the instance sets of [`data/`][data] (where each set comes
from, in which format it is written, which reader of the library parses it,
and under which license it is redistributed), together with the options of
`AStarBench` that select them. The notation is that of [Chapter
1](01-introduction.md#ch-1); the heuristics and their properties
(admissibility, consistency) are those of [§2.4](02-best-first.md#sec-2-4), the
work distributions of HDA$^*$ those of [§9.4](09-hash-distributed.md#sec-9-4),
and the nblocks of Safe PBNF those of [Definition
11.1](11-safepbnf.md#def-11-1).

The structure of the chapter is as follows. Section 15.1 states the problem and
the interface of the domain, section 15.2 the costs of a column and their
relation with the score of an alignment, and section 15.3 the pairwise
heuristic, with its proof of admissibility and consistency, and the cost of an
expansion. Section 15.4 gives the features and the abstractions of the domain
for the parallel searches, and section 15.5 the levelled graphs and the
hyperplane work distribution, with the bound on the number of workers that the
successors of the states of a worker reach. Section 15.6 describes the instance
sets, section 15.7 the readers of their files, and section 15.8 the options of
`AStarBench`; sections 15.9 to 15.11 give the benchmark, the notes and the
exercises.

## 15.1 The problem and the interface {#sec-15-1}

We are given $n \ge 2$ sequences $s_0, \ldots, s_{n-1}$ over the alphabet of
the amino acids, of lengths $L_0, \ldots, L_{n-1}$, and a cost $c(a, b) \ge 0$
for each pair of symbols, the gap $-$ included, with $c(-, -) = 0$. An
alignment inserts gaps into the sequences so that they have a common length
$m$, and it is read as $m$ columns, the column $k$ holding the $k$-th symbol of
each padded sequence; no column is made of gaps only. The cost of a column
$x = (x_0, \ldots, x_{n-1})$ is the *sum of pairs*
$$\mathrm{cost}(x) = \sum_{0 \le i < j < n} c(x_i, x_j),$$
and the cost of an alignment is the sum of the costs of its columns. We want an
alignment of least cost.

<a id="def-15-1"></a>**Definition 15.1** (The lattice of an alignment). The
*lattice* of the sequences is the set of the points
$p = (p_0, \ldots, p_{n-1})$ with $0 \le p_i \le L_i$, where $p_i$ is the
number
of the residues of $s_i$ already aligned. There is an edge from $p$ to
$p + e_M$ for each nonempty set $M \subseteq \{i : p_i < L_i\}$, where $e_M$ is
the vector with $1$ on the coordinates of $M$ and $0$ elsewhere, which stands
for the column whose symbol $i$ is the residue $s_i[p_i]$ if $i \in M$ and the
gap otherwise; its cost is the cost of that column. Its initial state is the
origin and the goal is $(L_0, \ldots, L_{n-1})$.

The alignments are in one-to-one correspondence with the paths from the origin
to the goal (the $k$-th edge of the path is the $k$-th column), and the cost
of a path is the cost of its alignment; hence an optimal alignment is a
shortest path of the lattice [[41]](Z-bibliography.md#ref-41). The lattice has
$\prod_i (L_i + 1)$ points and, at a point where all the sequences are still
open,
$2^n - 1$ successors; it is acyclic, since each edge raises $\sum_i p_i$.
[Figure 15.1](#fig-15-1) shows the lattice of two short sequences and an
optimal path.

![Figure 15.1. The lattice of the sequences $s_0 = \mathtt{WAC}$ and
$s_1 = \mathtt{WC}$: the point $p_0 p_1$ has the successors that advance $s_0$
(right), $s_1$ (up) or both (diagonal). The path in the primary colour from the
origin to the goal is an optimal alignment with the costs of PAM250: the column
$\mathtt{W}/\mathtt{W}$ costs $0$, $\mathtt{A}$ against a gap $25$ and
$\mathtt{C}/\mathtt{C}$ $5$, in all $30$ (right, the gap in the accent
colour).](../figures/msa-instances-lattice.svg){#fig-15-1}

In the library, the domain is the class template `SequenceAlignment< K >`
([`SequenceAlignment.h:230`][SequenceAlignment.h:230]), whose parameter is the
largest number of sequences (10 by default, from 2 to 16 by a
`static_assert`). A state is a `std::array< std::uint16_t , K >` of positions
(those beyond the $n$-th are 0), which is why a sequence must have fewer than
65535 residues, and the cost type is `int`. Its symbols and matrices of scores
are declared at namespace scope as follows:

```cpp
--8<-- "snippets/msa-instances-sig-symbols.cpp"
```

i.e., `amino_symbols` lists the 20 amino acids, the ambiguity codes `B`, `J`,
`Z` and `X`, and the gap `*` (`amino_count` symbols, the gap at the index
`amino_gap`), and a `ScoreMatrix` gives the score $S(a, b)$ of each pair of
them. The constant `pam250` is the matrix at 250 PAMs of Dayhoff, Schwartz and
Orcutt [[44]](Z-bibliography.md#ref-44), with the scores of the ambiguity codes
and the score $-8$ of a residue against a gap. The members of the class that a
user calls are the following (quoted from the header, without their bodies):

```cpp
explicit SequenceAlignment( const std::vector< std::string > & sequences ,
                            const ScoreMatrix & scores = pam250 ,
                            int block = 4 , int nblock_side = 16 ,
                            int penalty = 0 );
static std::vector< std::string > read_fasta( std::istream & in );
static ScoreMatrix read_scores( std::istream & in );
static std::vector< std::string > random_sequences( std::uint64_t seed ,
                                                    int n = 0 ,
                                                    int total = 5000 );
static int symbol( char a ) noexcept;
int size( void ) const noexcept;
int length( int i ) const noexcept;
cost_type pair_cost( int a , int b ) const noexcept;
std::vector< std::string > alignment(
                            const std::vector< state_type > & path ) const;
int block( void ) const noexcept;
int nblock_side( void ) const noexcept;
```

The constructor takes the sequences as strings of the symbols of
`amino_symbols` but the gap, the matrix of scores, the side `block` of the
boxes of Abstract Zobrist hashing and of AHDA$^*$ (section 15.4), the side
`nblock_side` of the nblocks of Safe PBNF, and a linear gap penalty (section
15.2; 0 keeps the gap scores of the matrix). It throws `std::invalid_argument`
if there are fewer than 2 or more than `K` sequences, if a sequence is empty or
has 65535 residues or more, if it holds a symbol that is not in
`amino_symbols` or the gap, if a side is not positive, if the penalty is
negative, or if a cost of section 15.2 is negative (i.e., if a gap scores more
than the best pair of residues). It then fills the tables of the heuristic
(section 15.3), and hence its cost is that of the dynamic programs of all the
pairs. `symbol( a )` is the index of the character `a` in `amino_symbols`
($-1$ if it is not there); `size()` and `length( i )` give $n$ and $L_i$;
`pair_cost( a , b )` gives $c(a, b)$ by the indices of the symbols; and
`alignment( path )` turns a path of states, e.g., the `path` of a `Solution`,
into the $n$ rows of its alignment, `-` standing for a gap. Its three static
members read the sequences of a FASTA file, read a matrix of scores, and draw
random sequences; we describe them in section 15.7.

The domain models the concepts `IncrementalFeatureDomain`,
`AbstractFeatureDomain`, `AbstractStateDomain`, `NBlockDomain`, `DiskDomain`
and `LevelDomain` ([§1.5](01-introduction.md#sec-1-5)), which the case `msa` of
the tester checks by a `static_assert`; it is neither an `OperatorDomain` nor
a `PairwiseDomain` nor a `PerfectHashDomain`. Hence it runs under A$^*$,
IDA$^*$, the search with delayed duplicate detection of
[§4.6](04-memory-bounded.md#sec-4-6) (whose states must be trivially
copyable), HDA$^*$ with the distributions Z, AZ, A, H and the hyperplanes, Safe
PBNF, and the other searches that need no more than these concepts.

## 15.2 The costs {#sec-15-2}

A substitution matrix gives a score, which is to be maximized, whereas a search
minimizes a cost. Hence the library turns the scores into costs by a shift. Let
$S_{\max}$ be the largest score $S(a, b)$ between two residues ($a$ and $b$
ranging over the 24 symbols that are not the gap), and let $P \ge 0$ be the
penalty given to the constructor. The cost of a pair of symbols is
$$c(a, b) = S_{\max} - S(a, b), \qquad c(-, -) = 0 \qquad (P = 0),$$
$$c(a, b) = 2\,(S_{\max} - S(a, b)), \quad c(a, -) = c(-, a) = S_{\max} + 2P,
\quad c(-, -) = 0 \qquad (P > 0),$$
$a$ and $b$ being residues
([`SequenceAlignment.h:284-304`][SequenceAlignment.h:284-304]). With PAM250,
$S_{\max} = S(\mathtt{W}, \mathtt{W}) = 17$, and with $P = 0$ we have
$c(\mathtt{W}, \mathtt{W}) = 0$, $c(\mathtt{C}, \mathtt{C}) = 5$ and a residue
against a gap costs $17 + 8 = 25$, the values that the case `msa` checks. The
costs between two residues are nonnegative by the choice of $S_{\max}$; those
against a gap are nonnegative unless, with $P = 0$, the matrix gives a gap a
score larger than $S_{\max}$, which the constructor refuses.

The shift by $S_{\max}$ changes the problem only by a constant and a gap
penalty, as follows.

<a id="prop-15-2"></a>**Proposition 15.2** (Minimum cost and maximum score).
Let an alignment induce on the pair $i < j$, once its columns of two gaps are
dropped, $m_{ij}$ columns of two residues, whose scores sum to $\Sigma_{ij}$,
and $g_{ij}$ columns of a residue and a gap. If $P = 0$ and $G$ is the score
of a residue against a gap (the same for all the residues, as in PAM250), the
cost of the alignment is
$$\sum_{i < j} \Bigl( S_{\max} \frac{L_i + L_j}{2} - \bigl( \Sigma_{ij} -
\bigl(\tfrac{S_{\max}}{2} - G\bigr) g_{ij} \bigr) \Bigr),$$
and if $P > 0$ it is
$$\sum_{i < j} \Bigl( S_{\max} (L_i + L_j) - 2\,( \Sigma_{ij} - P\,g_{ij} )
\Bigr).$$
Hence an alignment of least cost is one of largest sum-of-pairs score with the
linear gap penalty $S_{\max}/2 - G$ (16.5 for PAM250) if $P = 0$, and with the
linear gap penalty $P$ if $P > 0$.

*Proof.* Each residue of $s_i$ and of $s_j$ appears in exactly one of the
columns kept, a column of two residues holding two of them and a column of a
residue and a gap one, whence $2 m_{ij} + g_{ij} = L_i + L_j$. Columns of two
gaps cost $0$. If $P = 0$, the cost of the pair is
$$\sum (S_{\max} - S) + (S_{\max} - G)\,g_{ij} = S_{\max} (m_{ij} + g_{ij}) -
\Sigma_{ij} - G\,g_{ij},$$
the first sum running over the $m_{ij}$ columns of two residues, and
$m_{ij} + g_{ij} = (L_i + L_j)/2 + g_{ij}/2$ gives the first formula. If
$P > 0$, the cost of the pair is
$2 S_{\max} m_{ij} - 2 \Sigma_{ij} + (S_{\max} + 2P)\,g_{ij} =
S_{\max}(2 m_{ij} + g_{ij}) - 2(\Sigma_{ij} - P g_{ij})$, which is the second.
In both cases the first term does not depend on the alignment, and the cost is
a decreasing function of the score. $\square$

The penalty $P = 8$ is that of the experiments of Jinnai and Fukunaga
([[5]](Z-bibliography.md#ref-5), section 4.1.3, as the header reports it), and
with $P > 0$ the costs stay integers, the doubling absorbing the half of
$S_{\max}$. We check both statements, in the case `msa`, on the pair
`HEAGAWGHEE`,
`PAWHEAE`: the optimal cost with $P = 8$ equals $17 (L_0 + L_1) - 2\sigma$,
$\sigma$ being the best score with the penalty 8 computed by the recurrence of
Needleman and Wunsch, and the optimal cost with $P = 0$ equals
$(17 (L_0 + L_1) - 2\sigma')/2$, $\sigma'$ being the best score with the
penalty 16.5. A user who wants another matrix gives it to the constructor,
e.g., as read by `read_scores` (section 15.7).

## 15.3 The pairwise heuristic {#sec-15-3}

The heuristic of the library is the classical one of the optimal pairwise
alignments [[41]](Z-bibliography.md#ref-41): for each pair $i < j$, let
$D_{ij}(q_i, q_j)$ be the least cost of an alignment of the suffixes
$s_i[q_i..]$ and $s_j[q_j..]$ with the costs $c$, and let
$$h(p) = \sum_{0 \le i < j < n} D_{ij}(p_i, p_j).$$
The tables $D_{ij}$ are filled once by the constructor, backwards from the
corner $(L_i, L_j)$, by the dynamic program of the following algorithm, which
transcribes the private member `pairwise_tables`
([`SequenceAlignment.h:642`][SequenceAlignment.h:642]; the header has no
pseudocode of it).

**Algorithm 15.1** (The tables of the pairwise heuristic).

```text
 1  procedure PAIRWISETABLES()
 2     for each pair i < j do
 3        for q_i ← L_i down to 0 do
 4           for q_j ← L_j down to 0 do
 5              D ← 0 if (q_i, q_j) = (L_i, L_j), ∞ otherwise
 6              if q_i < L_i then D ← min(D, D_ij(q_i + 1, q_j) + c(s_i[q_i], −))
 7              if q_j < L_j then D ← min(D, D_ij(q_i, q_j + 1) + c(−, s_j[q_j]))
 8              if q_i < L_i and q_j < L_j then
 9                 D ← min(D, D_ij(q_i + 1, q_j + 1) + c(s_i[q_i], s_j[q_j]))
10              D_ij(q_i, q_j) ← D
```

Then `heuristic( p )` ([`SequenceAlignment.h:442`][SequenceAlignment.h:442])
sums $n(n-1)/2$ entries of the tables.

<a id="lem-15-3"></a>**Lemma 15.3** (The tables). Algorithm 15.1 sets
$D_{ij}(q_i, q_j)$ to the least cost of an alignment of $s_i[q_i..]$ and
$s_j[q_j..]$, in time and memory $O((L_i + 1)(L_j + 1))$ for the pair; the
tables of all the pairs take $O(\sum_{i<j} (L_i + 1)(L_j + 1))$ words.

*Proof.* The lattice of the two suffixes is the lattice of Definition 15.1 for
$n = 2$, restricted to the points $\ge (q_i, q_j)$, and the first column of an
alignment of the suffixes advances $s_i$, $s_j$ or both; lines 5 to 9 are the
Bellman equation of the shortest path to the corner on this acyclic graph,
whose successors come earlier in the order of the loops of lines 3 and 4. By
induction on $(L_i - q_i) + (L_j - q_j)$ each entry is the least cost of a path
to the corner. Each entry costs $O(1)$. $\square$

<a id="thm-15-4"></a>**Theorem 15.4** (The pairwise heuristic is admissible
and consistent). For each point $p$ of the lattice, $h(p)$ is at most the
least cost of a path from $p$ to the goal, $h$ is $0$ at the goal, and
$h(p) \le c(p, p') + h(p')$ for each edge $(p, p')$.

*Proof.* A path from $p$ to the goal is an alignment $A$ of the suffixes
$s_i[p_i..]$; its cost is $\sum_{i<j} \sum_{x \in A} c(x_i, x_j)$, and for the
pair $i < j$ the columns of $A$, once those with two gaps on $i$ and $j$ are
dropped (they cost $c(-, -) = 0$), form an alignment of $s_i[p_i..]$ and
$s_j[p_j..]$, whose cost is at least $D_{ij}(p_i, p_j)$ by Lemma 15.3. Summing
over the pairs gives $h(p) \le$ the cost of $A$. At the goal each table is
read at its corner, where it is $0$. Let now the edge $(p, p')$ be the column
$x$. For the pair $i < j$, either $x_i = x_j = -$, and then $p'_i = p_i$,
$p'_j = p_j$ and $c(x_i, x_j) = 0$, or $x$ restricted to $\{i, j\}$ is an edge
of the lattice of the pair from $(p_i, p_j)$ to $(p'_i, p'_j)$, of cost
$c(x_i, x_j)$, and the Bellman equation of the proof of Lemma 15.3 gives
$D_{ij}(p_i, p_j) \le c(x_i, x_j) + D_{ij}(p'_i, p'_j)$. In both cases the
inequality holds for the pair, and summing over the pairs gives
$h(p) \le \mathrm{cost}(x) + h(p')$. $\square$

Since $h$ is consistent, A$^*$ never reopens a state of the lattice
([Theorem 2.13](02-best-first.md#thm-2-13)), which the case `msa` checks on the
instances of 4 sequences of [`data/msa/easy`][easy] against the exact
distances of all the points (`check_consistent_astar`), together with the cost
of A$^*$ against a dynamic program over the whole lattice written from the
definition of the sum of pairs (`alignment_dp` of
[`test/files.cpp`][files.cpp]). It is worth noting that $h$ is exact for
$n = 2$ and that, with more sequences, it ignores the interaction among the
pairs; the tighter heuristics of the literature (e.g., from the optimal
alignments of triples) are not in the library.

Each expansion generates the successors by `for_each_successor`
([`SequenceAlignment.h:457`][SequenceAlignment.h:457]), which collects the set
$O$ of the sequences not yet ended and the residues $s_i[p_i]$ of these, and
then enumerates the nonempty subsets $M \subseteq O$ as bit masks in
decreasing order (`m = open` first, then `m = ( m - 1 ) & open`); thus the
first successor is the column that advances all the open sequences. For each
one
it computes the cost of the column and the heuristic of the successor.

<a id="prop-15-5"></a>**Proposition 15.5** (The cost of an expansion). An
expansion of a point whose $k$ sequences are open generates $2^k - 1$
successors in time $O(2^k n^2)$, the cost of a column and the heuristic of a
successor costing $O(n^2)$ each; the goal test and the hash of a state cost
$O(K)$.

*Proof.* The loop over the masks visits each nonempty subset of $O$ once, and
for each one the two double loops over the pairs take $O(n^2)$, the table
lookups costing $O(1)$. For the goal test, the code compares two arrays of $K$
positions, and
`hash` mixes the positions four at a time, 16 bits each, by SplitMix64.
$\square$

Hence the domain is the opposite of the puzzles of [Chapter
14](14-domains.md#ch-14): the expansion is costly and the branching factor is
large (127 successors for 7 open sequences), and with HDA$^*$ the
communication is relatively cheap and the search overhead and the load balance
weigh more.

## 15.4 The features, the abstractions and the nblocks {#sec-15-4}

The parallel searches need the structure of the lattice in three forms, which
the domain provides as follows.

### The features of Zobrist hashing

The features of a point $p$ are the pairs $(i, p_i)$, one for each sequence,
numbered from the offset of the sequence $i$ in the list of all the pairs
(`num_features()` is $\sum_i (L_i + 1)$, and `for_each_feature( p , f )` calls
`f` for the $n$ features of $p$). A move changes the features of the sequences
it advances, and `for_each_feature_change( p , p' , f )` calls `f` for the old
and the new feature of each of them, whence the incremental Zobrist key of
[§9.5](09-hash-distributed.md#sec-9-5) costs $O(|M|)$ per successor; the case
`msa` checks that the keys of the successors computed from the key of the
parent are those computed from scratch, for 1, 3 and 64 workers.

### The boxes of Abstract Zobrist hashing and of AHDA$^*$

With the side $b$ (`block()`, 4 by default, as in section 4.1.3 of
[[5]](Z-bibliography.md#ref-5) as the header reports it), the abstract feature
of $(i, p_i)$ is $(i, \lfloor p_i / b \rfloor)$ (`abstract_feature( k )` maps
the feature $k$ to it), and `num_abstract_features()` is
$\sum_i (\lfloor L_i / b \rfloor + 1)$. The abstract state of AHDA$^*$,
`abstract_state( p )`, is the box
$(\lfloor p_0 / b \rfloor, \ldots, \lfloor p_{n-1} / b \rfloor)$ of $b^n$
points, written in mixed radix. Under both distributions the owner of a point
is a function of its box, and hence each move inside a box keeps the owner of
the state; the moves that cross the boundary of a box (a coordinate $p_i$
reaching a multiple of $b$) may change it. The abstraction of AHDA$^*$ of
[[5]](Z-bibliography.md#ref-5) keeps the position in the longest sequence only
(as the header reports), i.e., it cuts the lattice into slices instead of
boxes; ours is the product abstraction. On the instance `test3`, with boxes of
side 3 and 64 workers, the case `msa` checks that both AZH and AHDA$^*$ keep
on the owner of the parent more than four times as many successors as Zobrist
hashing.

### The nblocks of Safe PBNF

The nblocks are the boxes of side $\nu$ (`nblock_side()`, 16 by default),
$(\lfloor p_0 / \nu \rfloor, \ldots, \lfloor p_{n-1} / \nu \rfloor)$, with
$\lfloor L_i / \nu \rfloor + 1$ boxes along the coordinate $i$
(`num_nblocks()` is their product, and `nblock( p )` the box of $p$ in mixed
radix), and `for_each_nblock_successor( b , f )`
([`SequenceAlignment.h:573`][SequenceAlignment.h:573]) calls `f` for each box
that differs from $b$ by $+1$ in a nonempty set of coordinates, or by $-1$ in a
nonempty set of coordinates, within the lattice.

<a id="prop-15-6"></a>**Proposition 15.6** (The boxes are an abstraction). The
boxes of side $\nu$ with the relation of `for_each_nblock_successor` satisfy
the conditions (i) and (ii) of [Definition 11.1](11-safepbnf.md#def-11-1), and
each box has at most $2(2^n - 1)$ abstract successors.

*Proof.* (i) An edge $p \to p + e_M$ raises each coordinate in $M$ by one and
leaves the others. Hence each coordinate of the box rises by $0$ or $1$, all
upwards; if none rises the box is the same, and otherwise the box of $p + e_M$
differs from that of $p$ by $+1$ in a nonempty set of coordinates, which the
code lists (both boxes lie within the lattice). (ii) The code lists the boxes
that differ by $+1$ on a set $T$ and those that differ by $-1$ on a set $T$;
$b'$ is listed from $b$ with $+1$ on $T$ if and only if $b$ is listed from $b'$
with $-1$ on $T$, and $b$ itself is never listed, since $T$ is not empty. The
bound counts the nonempty subsets in the two directions. $\square$

The edges of the lattice all go upwards, and the graph of the boxes is in
fact directed; the downward successors are listed only because `SafePBNF`
takes the interference scope of an nblock from its successors alone and needs
the relation symmetric ([§11.1](11-safepbnf.md#sec-11-1)). In the case `msa`,
we check both conditions on all the moves of the instance `test3` with boxes of
side 4.

## 15.5 Levelled graphs and the hyperplane work distribution {#sec-15-5}

The level of a point $p$ is $\ell(p) = \sum_i p_i$ (`level( p )`), the number
of the residues aligned, and `max_level()` is $\sum_i L_i$, the level of the
goal. Each move raises the level by $|M| \in \{1, \ldots, n\}$, and the lattice
is therefore a *levelled graph*, the structure that the concept `LevelDomain`
of [`Domain.h`][Domain.h] states: a nonnegative integer `level( s )` that every
move raises, and the largest one, `max_level()`, that of the goals.

The hyperplane work distribution of Kobayashi, Kishimoto and Watanabe
([[43]](Z-bibliography.md#ref-43), section 4.2, as the header reports it) gives
the hyperplanes $\ell(p) = c$ of the lattice, or slabs of $d$ consecutive ones,
to the workers in turn. Its rationale is the following. In the lattice of MSA a
point is reached by many paths of different costs, and with Zobrist hashing its
copies are generated by workers spread over the whole machine (the successors
of a point go to up to $\min\{p, 2^n - 1\}$ workers); the duplicates are then
detected only after the messages have crossed, and a worker that receives a
cheaper path late has to reopen a state that it had already expanded, which
matters with the nonuniform costs of this domain. If, instead, the successors
of the states of a worker go to a few workers only, the copies of a state meet
more often on the same worker in the order of their $f$, and fewer states are
reopened.

<a id="def-15-7"></a>**Definition 15.7** (The hyperplane distribution). Let
$\ell$ be the level of a levelled graph, $p$ the number of workers, $Z$ a hash
of the states, and the thickness $d$ either a positive integer or the inverse
$1/k$ of one. The plane of a state $s$ is
$$\mathit{Plane}(s) = \lfloor \ell(s) / d \rfloor \quad (d \ge 1), \qquad
\mathit{Plane}(s) = k\,\ell(s) + (Z(s) \bmod k) \quad (d = 1/k),$$
and the owner of $s$ is $\mathit{Plane}(s) \bmod p$.

With $d \ge 1$ the slab of $d$ consecutive levels goes to one worker, and with
$d = 1/k$ each level is split by the hash among $k$ consecutive planes, i.e.,
a level is shared by $k$ workers; consecutive planes go to consecutive workers.
In the library, the distribution is the class template
`HyperplaneDistribution< D >`
([`Distribution.h:536-591`][Distribution.h:536-591]), for any `D` that models
`LevelDomain`; its members are the following (quoted from the header, without
their bodies):

```cpp
explicit HyperplaneDistribution( const D & domain , std::uint64_t seed = 1 ,
                                 unsigned planes = 1 , unsigned split = 1 );
unsigned owner( const state_t< D > & s , unsigned p ) const;
static std::pair< unsigned , unsigned > thickness(
               const D & domain , unsigned p , double lambda = 0.003 );
static std::string name( void );
```

The thickness is $d = $ `planes` if `split` is 1 and $d = 1/$`split` if
`planes` is 1; the constructor throws `std::invalid_argument` if either is 0 or
both exceed 1. `owner( s , p )` ([`Distribution.h:555`][Distribution.h:555])
computes the plane of Definition 15.7, with $Z(s)$ the hash of the `Domain`
mixed with the seed by SplitMix64 (Kobayashi et al. take the Zobrist hash, of
which only the residue modulo $k$ is used), and reduces it modulo $p$; unlike
the other distributions of [§9.4](09-hash-distributed.md#sec-9-4), it does not
use the multiply-high reduction, since the consecutive planes must go to
consecutive workers. `name()` is `"Hyperplane"`. This distribution has no
incremental key, and HDA$^*$ calls `owner` on each successor, at the cost of
`level` (a sum of $n$ positions) and, for $d < 1$, of the hash of the state.

![Figure 15.2. (a) The owners of the points of a lattice of two sequences under
the hyperplane distribution of thickness $d = 2$ with $p = 3$ workers: the
dashed lines bound the slabs of two levels, and the successors of the point in
the primary colour (level 2) go to the workers 1 and 2. (b) The boxes of side
$b = 2$ of AZH and of AHDA$^*$, each owned by one worker (the owners drawn for
the figure): each move inside a box keeps the owner of its state, while the
three successors of the point in the primary colour lie in three boxes.
](../figures/msa-instances-hyperplanes.svg){#fig-15-2}

The following theorem bounds the number of workers that the successors of the
states of a worker reach, over all the states of that worker; [Figure
15.2](#fig-15-2) shows it for $d = 2$.

<a id="thm-15-8"></a>**Theorem 15.8** (The workers reached from a worker). Let
every move raise the level by at least 1 and at most $\Delta$ ($\Delta = n$ for
the alignment of $n$ sequences), and let $w$ be a worker. Under the hyperplane
distribution, the successors of all the states owned by $w$ are owned by at
most
$$\min\{p, \lceil \Delta / d \rceil + 1\} \ \ (d \ge 2), \qquad
\min\{p, \Delta\} \ \ (d = 1), \qquad
\min\{p, (\Delta + 1) k - 1\} \ \ (d = 1/k)$$
workers, namely the workers $w + j \bmod p$ for $j$ in
$\{0, \ldots, \lceil \Delta / d \rceil\}$, in $\{1, \ldots, \Delta\}$ and in
$\{1, \ldots, (\Delta + 1) k - 1\}$ respectively.

*Proof.* Let $s$ be a state of $w$ and $t$ a successor of $s$, so that
$\ell(s) + 1 \le \ell(t) \le \ell(s) + \Delta$. Let first $d \ge 1$ be an
integer and $P = \lfloor \ell(s) / d \rfloor$, so that
$P d \le \ell(s) \le P d + d - 1$. Then $\mathit{Plane}(t) \ge P$, and
$$\mathit{Plane}(t) \le \Bigl\lfloor \frac{P d + d - 1 + \Delta}{d}
\Bigr\rfloor = P + \Bigl\lfloor \frac{\Delta + d - 1}{d} \Bigr\rfloor = P +
\lceil \Delta / d \rceil .$$
Since $P \equiv w \pmod p$, the owner of $t$ is $w + j \bmod p$ with
$0 \le j \le \lceil \Delta / d \rceil$, a set that does not depend on $s$. If
$d = 1$, $\mathit{Plane}(t) = \ell(t) > \ell(s) = P$, so that $j \ge 1$, which
leaves $\Delta$ values. Let now $d = 1/k$, with
$\mathit{Plane}(s) = k \ell(s) + r$ and $\mathit{Plane}(t) = k \ell(t) + r'$,
$0 \le r, r' < k$. Then
$$\mathit{Plane}(t) - \mathit{Plane}(s) \ge k(\ell(s) + 1) - k\ell(s) - r =
k - r \ge 1,$$
$$\mathit{Plane}(t) - \mathit{Plane}(s) \le k(\ell(s) + \Delta) + k - 1 -
k\ell(s) - r \le (\Delta + 1) k - 1,$$
and the owner of $t$ is $w + j \bmod p$ with $1 \le j \le (\Delta + 1)k - 1$.
In each case the number of the values of $j$, and the number $p$ of the
workers, bound the number of the owners. $\square$

<a id="rem-15-9"></a>**Remark 15.9** (The bound of Kobayashi et al.). The
header of `HyperplaneDistribution` reports that Remark 1 of
[[43]](Z-bibliography.md#ref-43) writes the bound $n/d + 1$ for $d \ge 1$,
which fails when $d$ does not divide $n$: with $n = 4$ and $d = 5$, a state of
level 4 (plane 0) has a successor of level 5 (plane 1), and the successors of a
worker reach two workers, more than $4/5 + 1$; Theorem 15.8 gives
$\lceil n / d \rceil + 1$ instead. For $d = 1/k$, the case `msa` checks
the bound $(n + 1)/d = (n + 1) k$, which Theorem 15.8 lowers by one. We
have not checked the paper itself. In the case `msa`, we check the bounds
$\lceil n/d \rceil + 1$ and $(n + 1) k$ with 64 workers on the instance `test3`
of 4 sequences, for $d \in \{1, 2, 5, 1/2, 1/4\}$, against the owners of all
the moves of the lattice.

Under Zobrist hashing, in contrast, a worker owns states scattered over the
whole lattice, and the successors of its states may reach all the $p$ workers;
under the boxes, the successors of a single state lie in at most $2^n - 1$
boxes other than its own, but the successors of all the states of a worker are
bounded by nothing smaller than $p$. Of course, the hyperplanes cost some load
balance: the levels are expanded roughly in the order of $f$, and if $d$ is
large, the few slabs that hold the nodes of the current $f$ belong to a few
workers. The thickness trades the two effects, and `thickness( domain , p ,
lambda )` gives the rule of section 4.3 of [[43]](Z-bibliography.md#ref-43) (as
the header reports it): with
$$x = \frac{\lambda\, L}{\log_2 p}, \qquad L = \texttt{max\_level()},$$
the thickness is $\mathrm{round}(x)$ if $x \ge 1$ and $1 / \mathrm{round}(1/x)$
otherwise, with $\lambda = 0.003$ as in that paper, returned as the pair
(`planes`, `split`) of the constructor, and $(1, 1)$ for $p < 2$. That paper
does not give the base of the logarithm (the header remarks), and the library
takes 2. For instance, on the instance `test3`, of total length $L = 73$, with
$p = 64$, $x = 0.003 \cdot 73 / 6 = 0.0365$ and the rule gives $d = 1/27$, the
value that the case `msa` checks.

## 15.6 The instance sets {#sec-15-6}

The directory [`data/`][data] holds the instances on which the tests and the
benchmarks run, each set with a `NOTICE` (or a `LICENSE`, or comments at the
head of the file) that says where it comes from and under which license it is
redistributed. Table 15.1 summarizes them; the paragraphs that follow give the
details, and section 15.7 the readers.

**Table 15.1.** The instance sets of [`data/`][data]: the directory, what it
holds, the reader that parses it (section 15.7, or the chapter that describes
the domain), the source and the license of the files.

| directory | contents | reader | source, license |
|---|---|---|---|
| `15-puzzle` | [`korf100.txt`][korf100.txt] (100 boards), `korf_1`, `korf_12`, `korf_55` | `SlidingTile::parse`, `SlidingTile::read_positions` | Korf [[7]](Z-bibliography.md#ref-7); the code of PBNF, MIT |
| `24-puzzle` | [`instances50.txt`][instances50.txt], [`instances100.txt`][instances100.txt], [`fromkorf8.txt`][fromkorf8.txt] | `SlidingTile::parse` | the code of HDA$^*$ and the benchmark set of Jinnai, MIT |
| `grid` | 5 boards of $10 \times 8$ to $500 \times 500$ cells | `GridBoard::read` | the code of PBNF and the benchmark set of Jinnai, MIT |
| `msa` | 36 files `.tfa` in `balibase`, `easy`, `papers` | `SequenceAlignment::read_fasta` | BAliBASE 3.0 [[45]](Z-bibliography.md#ref-45); the sets of Jinnai, MIT |
| `tsp` | `tsp_40_0` to `tsp_40_9`, 40 cities each | `TravelingSalesman::read_cities` | the benchmark set of Jinnai, MIT |
| `rubik` | [`korf10.txt`][korf10.txt], 10 scrambles of 100 moves | `RubiksCube::parse_moves` | Korf [[59]](Z-bibliography.md#ref-59); HOG2 [[165]](Z-bibliography.md#ref-165), MIT |
| `hanoi` | [`toh20.txt`][toh20.txt], 20 instances of 4 pegs and 20 disks | `TowersOfHanoi::read` | Siag et al. [[57]](Z-bibliography.md#ref-57), MIT |
| `dimacs` | `sample.gr`, `sample.co` (900 nodes), `brock200_1.clq` | `read_dimacs_graph`, `read_dimacs_coordinates` | written by [`tools/bench/dimacs_sample.py`][dimacs_sample.py]; the clique graph from YewPar, MIT |
| `qaplib` | 15 instances `.dat` of 12 to 16 facilities | the QAP of [Chapter 13](13-branch-and-bound.md#ch-13) | QAPLIB [[78]](Z-bibliography.md#ref-78); P3D-DFS, BSD 3-clause |
| `planning` | 28 tasks `output.sas` and [`optimal.txt`][optimal.txt] | `SASTask::read_file` | IPC, translated by [`tools/pddl2sas.py`][pddl2sas.py] ([Chapter 16](16-planning.md#ch-16)) |
| `robot-nav-2d` | 5 Moving AI maps, 50 pairs of points each | `NavMap::read_movingai`, `RobotNav2D::read_points` | the code of ePA$^*$SE [[121]](Z-bibliography.md#ref-121), MIT |

### The sliding-tile puzzles

[`data/15-puzzle/korf100.txt`][korf100.txt] holds the 100 random instances of
the 15-puzzle of Korf [[7]](Z-bibliography.md#ref-7), one per line, the 16
tiles cell by cell (0 being the blank and the goal $0, 1, \ldots, 15$) followed
by the optimal cost, as described in [§14.2](14-domains.md#sec-14-2). The files
`korf_1`, `korf_12` and `korf_55` are the instances 1, 12 and 55 in the format
of Burns et al. [[3]](Z-bibliography.md#ref-3), which gives the cell of each
tile (section 15.7), taken from the directory src/tiles of the code of PBNF under the MIT
license of its authors (`NOTICE`). In [`data/24-puzzle`][24-puzzle],
[`instances50.txt`][instances50.txt] holds the 50 instances of the 24-puzzle
distributed with the code of HDA$^*$ of Jinnai and Fukunaga
[[5]](Z-bibliography.md#ref-5), [`instances100.txt`][instances100.txt] the 100
instances of the benchmark set of Jinnai (of which the instances 1 to 7 and 11
are those of [`fromkorf8.txt`][fromkorf8.txt], and the instances 82, 83 and 84
repeat the instances 8, 9 and 10), and [`fromkorf8.txt`][fromkorf8.txt] the 8
instances of that set taken from Korf, which its list `instance_numbers_korf`
numbers 1, 2, 7, 8, 10, 12, 15 and 18; each file gives, in its comments, its
origin and the MIT license of the code it comes from, and holds one board per
line, the 25 tiles cell by cell. In the case `files`, we check that the three
files `korf_*` are the corresponding instances of [`korf100.txt`][korf100.txt],
and that the instances of [`instances100.txt`][instances100.txt] and
[`fromkorf8.txt`][fromkorf8.txt] are boards from which the goal is reachable (a
board of the wrong parity is refused by `SlidingTile::parse`,
[§14.2](14-domains.md#sec-14-2)).

### The grids

[`data/grid`][grid] holds five boards in the format of Burns et al.
[[3]](Z-bibliography.md#ref-3) (section 15.7): `ethan_unit_four-way_1` (a board
of $10 \times 8$ cells), `jordan_unit_four-way_100x100_1`,
`jordan_unit_four-way_500x500_1` and `jordan_life_four-way_500x500_1`, from
the directory src/grid of the code of PBNF, whose row 0 is the last line of the board; and
`grid_500-500-0.7_0`, a board of $500 \times 500$ cells each free with
probability 0.7, from the benchmark set of Jinnai, whose row 0 is the first
line. The larger boards of that set ($1000 \times 1000$ to $5000 \times 5000$
cells, and the boards of $500 \times 500$ cells with the probabilities 0.45 to
0.6) are not redistributed; they are written by its generator
grid-gen/grid_gen.cc, which draws each cell free with the given probability
and the start and the goal in the first and in the last third of each side,
and the reader of the library reads them as they are. Both repositories are
under the MIT license (`NOTICE`). In the case `board`, we check the optimal
costs of
the five boards (23, 143, 763, 362873 and 755) against Dijkstra's algorithm on
the board and against the costs that the code of their authors finds.

### Multiple sequence alignment

[`data/msa`][msa] holds 36 files in the FASTA format (`.tfa`) of BAliBASE 3.0
of Thompson et al. [[45]](Z-bibliography.md#ref-45), in three directories:

* `balibase`: the 15 files of the directory `msa/BAliBASE` of the benchmark set
  of Jinnai (4 to 6 sequences of 51 to 974 residues), sequences of BAliBASE 3.0
  as selected there;
* `easy`: the 6 files `test0` to `test5` of `BAliBASE/easy` of the code of
  HDA$^*$ of Jinnai and Fukunaga, short sequences cut from the same benchmark
  (3 to 7 sequences of 3 to 46 residues);
* `papers`: the 15 files of the directories RV11 and RV12 of BAliBASE 3.0 used
  in the experiments of Kobayashi et al. [[43]](Z-bibliography.md#ref-43)
  (`BB12003`, `BB12022`, `BB12023`, `BB12032`, `BB12036`, `BBS11026`,
  `BBS11037`, `BBS12005`, `BBS12010`, `BBS12023`, `BBS12032`) and of Jinnai
  and Fukunaga ([[5]](Z-bibliography.md#ref-5), section 4.1.3: `BB12021`,
  `BB12022`, `BB12036`, `BBS11010`, `BBS11026`, `BBS11035`, `BBS11037`,
  `BBS12016`, `BBS12023`, `BBS12032`), taken from the archive of BAliBASE (4
  to 9 sequences of 58 to 738 residues); `BBS11035` is also in `balibase`.

The files of `balibase` and `easy` come from repositories under the MIT
license; the archive of BAliBASE states no license, and its README names J. D.
Thompson as the contact (`NOTICE`). In the case `msa`, we check the optimal
costs of five of the instances of `easy` against the dynamic program on the
whole
lattice (those of 4 sequences) and against the costs that the code of HDA$^*$
of Jinnai and Fukunaga finds; Table 15.2 gives them, with the size of the
lattice. It also reads the 15 files of `balibase` and computes their
heuristic at the origin, solves `test2` (7 sequences) with A$^*$, and solves
`test0` with IDA$^*$ and with the search with delayed duplicate detection.

**Table 15.2.** The instances of [`data/msa/easy`][easy] with the costs of
PAM250 and no gap penalty ($P = 0$): the lengths of the sequences, the number
of the points of the lattice, and the optimal cost, as checked by the case
`msa` of [`test/files.cpp`][files.cpp].

| instance | lengths | points of the lattice | optimal cost |
|---|---|---|---|
| `test0` | 9, 5, 16, 15 | 16 320 | 1741 |
| `test1` | 31, 12, 6, 11 | 34 944 | 2665 |
| `test2` | 17, 13, 13, 9, 16, 3, 5 | 14 394 240 | 6055 |
| `test3` | 19, 19, 20, 15 | 134 400 | 2022 |
| `test4` | 46, 38, 39 | 73 320 | 2312 |
| `test5` | 17, 4, 7, 3, 6, 3, 5 | 483 840 | 3990 |

Random instances, as those of section 4.1.3 of [[5]](Z-bibliography.md#ref-5),
are drawn by `random_sequences` (section 15.7) and are not stored.

### The traveling salesperson

[`data/tsp`][tsp] holds the first 10 files of `tsp/40_cities` of the benchmark
set of Jinnai, `tsp_40_0` to `tsp_40_9`, each the number of the cities (40)
followed by the coordinates $x$ and $y$ of each city in the unit square, under
the MIT license (`NOTICE`). The domain of [§14.7](14-domains.md#sec-14-7) takes
the first $n$ cities of a file. In the case `files`, we check the optimal tours
of
the first 6, 8, 10 and 12 cities of `tsp_40_0` with the Manhattan distances
(28975, 31254, 32734 and 33861, those that the code of HDA$^*$ of its authors
finds). It also checks the tours of 6, 9 and 12 cities with both metrics
against the dynamic program of Held and Karp [[14]](Z-bibliography.md#ref-14).

### Rubik's cube and the Towers of Hanoi

[`data/rubik/korf10.txt`][korf10.txt] holds the ten instances of Korf
[[59]](Z-bibliography.md#ref-59), each a scramble of 100 moves in the notation
of Singmaster (`-` for a quarter turn counterclockwise, `2` for a half turn),
one per line, from `apps/rubik` of HOG2 [[165]](Z-bibliography.md#ref-165)
under the MIT license. [`data/hanoi/toh20.txt`][toh20.txt] holds the 20
instances of the Towers of Hanoi with four pegs and 20 disks of Siag et al.
[[57]](Z-bibliography.md#ref-57), each line giving the four pegs of the start
and then the four pegs of the goal, a peg being its number in parentheses
followed by its disks from the bottom, from the code of that paper under the
MIT license. The domains and the readers of both (`RubiksCube::parse_moves`,
`TowersOfHanoi::read`) are described with the permutation puzzles and their
pattern databases.

### The DIMACS graphs

[`data/dimacs`][dimacs] holds `sample.gr` and `sample.co`, a synthetic graph of
900 nodes and 3106 arcs (a perturbed grid of $30 \times 30$ nodes, the lengths
in meters) in the format of the road graphs of the 9th DIMACS Implementation
Challenge [[64]](Z-bibliography.md#ref-64), written by
[`tools/bench/dimacs_sample.py`][dimacs_sample.py] and small enough for the
tests; the graphs of the challenge (from NY, of 264346 nodes, to USA, of
23947347) are downloaded by [`tools/bench/fetch_dimacs.sh`][fetch_dimacs.sh].
It also holds `brock200_1.clq`, a graph of 200 vertices and 14834 edges of the
second DIMACS Implementation Challenge (by M. Brockington)
[[174]](Z-bibliography.md#ref-174),
copied from the repository YewPar of B. Archibald under the MIT license
(`LICENSE`), on which the case `bb-domains` runs the maximum clique of [Chapter
13](13-branch-and-bound.md#ch-13).

### The other sets

[`data/qaplib`][qaplib] holds 15 instances of QAPLIB
[[78]](Z-bibliography.md#ref-78) (`chr12a` to `tai12a`), copied from the
repository P3D-DFS of G. Helbecque under the BSD 3-clause license (`LICENSE`),
for the quadratic assignment of [Chapter 13](13-branch-and-bound.md#ch-13).
[`data/planning`][planning] holds the 28 planning tasks of [Chapter
16](16-planning.md#ch-16), whose list [`optimal.txt`][optimal.txt] states their
origin and the absence of a license of the tasks of the International Planning
Competitions. [`data/robot-nav-2d`][robot-nav-2d] holds five maps of the Moving
AI benchmarks [[15]](Z-bibliography.md#ref-15) (`brc203d`, `den501d`,
`den520d`, `hrt201n` and `ht_chantry`), each in a directory with the files
*nav2d_starts.txt* and *nav2d_goals.txt* of its 50 pairs of points, the
instances of the experiments of the edge-based parallel searches
[[121]](Z-bibliography.md#ref-121), [[122]](Z-bibliography.md#ref-122); the
directory has the MIT license of the code of those experiments (`LICENSE`) and
no `NOTICE`.

## 15.7 The readers of the instance files {#sec-15-7}

Each reader is a static member or a free function that takes a
`std::istream` (a test may thus read a file or a string alike) and throws
`std::invalid_argument` on a malformed input. The readers of the Moving AI maps
(`GridMap::read_movingai`) and of the boards of the puzzles
(`SlidingTile::parse`) are described in [§14.4](14-domains.md#sec-14-4) and
[§14.2](14-domains.md#sec-14-2).

### The sequences and the matrices

`SequenceAlignment::read_fasta( in )`
([`SequenceAlignment.h:336`][SequenceAlignment.h:336]) reads the FASTA format:
each line that starts with `>` (the name of a sequence) opens a new sequence,
which may span several lines; the lines before the first name are skipped, the
spaces and the characters `-` and `.` of an aligned file are dropped, and the
symbols are turned into upper case. It throws if no sequence is read. It does
not check the symbols, which the constructor of the domain does: in particular,
a file whose sequences end with `*` (as some FASTA files do) is read, and then
refused by the constructor, since `*` is the gap. In the case `msa`, we check
the names, the lines, the case, the gaps and the blank lines on a small string.

`SequenceAlignment::read_scores( in )`
([`SequenceAlignment.h:363`][SequenceAlignment.h:363]) reads a `ScoreMatrix`,
i.e., 625 integers row by row in the order of `amino_symbols`, the gap
included, and throws if fewer are read; the score of two gaps is read and
ignored, since $c(-, -)$ is always 0. In the case `msa`, we write `pam250` and
read it back.

`SequenceAlignment::random_sequences( seed , n , total )`
([`SequenceAlignment.h:383`][SequenceAlignment.h:383]) draws the random
instances of section 4.1.3 of [[5]](Z-bibliography.md#ref-5) (as the header
reports them): $n$ sequences ($n$ drawn uniformly in $\{4, \ldots, 9\}$ if it
is 0), each of a length drawn uniformly among the integers $l$ with
$0.9\,T/n < l < 1.1\,T/n$, $T$ being `total` (5000 there), and each residue
drawn uniformly among the 20 amino acids, all from `std::mt19937_64` seeded by
`seed` with the draws of [Lemma 1.9](01-introduction.md#lem-1-9), whence a
seed gives the same sequences on all platforms. It throws if $n < 2$, if $n$
exceeds `K`, or if no integer lies in the range of the lengths. For the seeds
1 to 20, the case `msa` checks the number of the sequences, the range of the
lengths, the alphabet and the repeatability, and solves a small instance
against the dynamic program.

### The puzzles by the cells of their tiles

`SlidingTile< W , H >::read_positions( in )`
([`SlidingTile.h:751`][SlidingTile.h:751]) reads the format of Burns et al.
[[3]](Z-bibliography.md#ref-3): the height and the width, the words
`starting positions for each tile:`, the cell of each tile from the blank (the
tile 0) on, the words `goal positions:` and the goal cell of each tile, the
cells counted row by row. It returns the pair of the start board and of
the goal board, and throws if the size is not $H \times W$ or if the cells are
not a permutation. A `SlidingTile` built from both boards searches toward that
goal; `AStarBench` recognizes the format by its words and otherwise reads one
board per line. We check the three files of Korf in the case `files`, and the
refusal of a board of another size and of a repeated cell.

### The boards of Burns et al.

`GridBoard::read( in , from_bottom )` ([`Grid.h:196`][Grid.h:196], the
`struct GridBoard`) reads a pathfinding instance in the format of Burns et al.
[[3]](Z-bibliography.md#ref-3): the height $H$ and the width $W$, the word
`Board:`, $H$ lines of $W$ characters (`#` a blocked cell, any other character
a free one; a short line is padded with free cells), a line `Unit` or `Life`
(the cost of leaving a cell: 1, or its row $y$ for the life costs), a line
`Four-way` or `Eight-way` (the moves; a diagonal move is allowed even if it
cuts a blocked corner), and the four numbers $x_s$, $y_s$, $x_g$, $y_g$ of the
start and of the goal. The files of Burns et al. count $y$ from the last line
of the board up, those of later generators from the first line down:
`from_bottom` tells which, and the cells are stored with $y$ counted from the
first line, as in a `GridMap`. A blocked start is made free, which changes no
shortest path, since a shortest path never comes back to its start (some
instances of the authors have one). The `GridBoard` holds the `GridMap`, the
start, the goal, the moves and the costs (`GridCost::cells` or
`GridCost::life`), and a `Grid< 4 >` or `Grid< 8 >` is built from it with the
moves it states (a `Grid< 4 >` refuses a board of 8 moves). In the case
`board`, we check the reader on both origins, its refusals (no `Board:`, other
costs or
moves, a goal outside, missing rows), the heuristics of the life costs, which
are exact on the boards without blocked cells, the cut corners, random boards
against Dijkstra's algorithm, and the five boards of the data.

### The cities

`TravelingSalesman::read_cities( in , metric , first )`
([`TravelingSalesman.h:142`][TravelingSalesman.h:142]) reads the number $n$ of
the cities and their coordinates, the format of the random instances of the
unit square of Jinnai and Fukunaga, keeps the first `first` cities if it is
positive, and returns the matrix of the distances, Euclidean
(`Metric::euclidean`, the default) or Manhattan (`Metric::manhattan`,
multiplied by $10^4$ and truncated to integer costs). It throws
on a malformed file or on one with fewer than `first` cities.

### The road graphs

The free functions of [`DimacsGraph.h`][DimacsGraph.h] read the road graphs of
the 9th DIMACS Implementation Challenge [[64]](Z-bibliography.md#ref-64).
`read_dimacs_graph< Cost >( in , source , target )` reads a file `.gr` (the
line `p sp n m` and $m$ lines `a u v w`, the arc $(u, v)$ of the integer cost
$w$, the nodes numbered from 1, each direction of a road an arc of its own)
into an `ExplicitGraph< Cost >` searched from `source` to `target`, and throws
if the line `p sp` is missing, if an arc is out of range or of negative cost,
or if fewer arcs are read than stated. `read_dimacs_coordinates( in , n )`
reads a file `.co` (the line `p aux sp co n` and $n$ lines `v v x y`, the
longitude and the latitude in millionths of a degree) into a
`DimacsCoordinates`, the two vectors `lon` and `lat`. The heuristic of a node
$u$ for the target $t$ is
$$h(u) = \lfloor \kappa\, \delta(u, t) \rfloor, \qquad \kappa = \min_{(u, v):
\ \delta(u, v) > 0} \frac{w(u, v)}{\delta(u, v)},$$
$\delta$ being the great-circle distance in meters, by the formula of the
haversine on the sphere of radius 6371008.8 m (`dimacs_distance`);
`dimacs_scale( g , c )` computes $\kappa$, lowered by the factor
$1 - 10^{-9}$ against the rounding of the distances, and
`set_dimacs_heuristic( g , c , kappa )` stores $h$ in the graph for its
current target (with $\kappa$ from `dimacs_scale` if `kappa` is negative) and
returns $\kappa$.

<a id="prop-15-10"></a>**Proposition 15.10** (The heuristic of the
coordinates). If $\delta$ satisfies the triangle inequality, $h$ is 0 at the
target and consistent, whatever the unit of the costs.

*Proof.* $h(t) = \lfloor 0 \rfloor = 0$. For an arc $(u, v)$ with
$\delta(u, v) > 0$, $\kappa\,\delta(u, v) \le w(u, v)$ by the definition of
$\kappa$, and the inequality holds trivially if $\delta(u, v) = 0$. Hence
$\kappa\,\delta(u, t) \le \kappa\,\delta(u, v) + \kappa\,\delta(v, t) \le
w(u, v) + \kappa\,\delta(v, t)$, and since $w(u, v)$ is an integer,
$h(u) \le \lfloor w(u, v) + \kappa\,\delta(v, t) \rfloor = w(u, v) + h(v)$.
$\square$

The great-circle distance is a metric, but its value in floating point is not
exactly one; the factor $1 - 10^{-9}$ guards against the rounding, and we have
not proved that it suffices on all graphs. In the case `dimacs`, we check the
consistency of $h$ on each arc of the sample for five random pairs of nodes,
and the cost of A$^*$ against the blind search.

## 15.8 Running the instances {#sec-15-8}

The benchmark `AStarBench` ([`tools/AStarBench.cpp`][AStarBench.cpp]) selects
the domain by `--domain`, the file by `--file` and the instances of a file by
`--ids` (all by default); the options that concern the sets of this chapter are
the following.

* `--domain tiles15 | tiles8 | tiles24 --file <file>`: a file of boards, one
  per line (the lines starting with `#` skipped), or a single instance in the
  format of `read_positions`, recognized by its words.
* `--domain board --file <board> --origin bottom | top`: a board of Burns et
  al., with the row 0 of its start and goal on the last line (`bottom`, the
  default, for the files of the code of PBNF) or on the first (`top`, for the
  files of the benchmark set of Jinnai); `--block` and `--nblock` give the
  blocks of AZH and AHDA$^*$ (100 by default) and the nblocks of PBNF.
* `--domain msa --file <file.tfa>` or `--random <seed>[:<n>[:<total>]]`: the
  sequences of a FASTA file, or the random ones of `random_sequences`; the
  matrix is PAM250 or the 625 scores of `--matrix <file>`, the gap penalty that
  of the matrix or the linear one of `--gap <P>` (0 by default), `--block` the
  side of the boxes (4 by default for `msa`) and `--nblock` that of the nblocks
  (16 by default for `msa`); `--dist HWD` selects the hyperplane distribution,
  of the thickness `--plane <d>`, `--plane 1/<k>` or `--plane auto` (the
  default), which takes the rule of `thickness()` for the number of workers of
  `--threads`.
* `--domain tsp --file <file> [--size <n>] [--metric euclidean | manhattan]`:
  the first $n$ cities of a file (all of them without `--size`).
* `--domain rubik --file data/rubik/korf10.txt`, `--domain hanoi --file
  data/hanoi/toh20.txt`, `--domain road --file <name>` (the files `<name>.gr`
  and `<name>.co`, between the nodes of `--points <s>,<t>` or random ones),
  `--domain nav2d --file <map>` (with the pairs of points next to the map), and
  `--domain planning --file <task.sas>` ([Chapter 16](16-planning.md#ch-16)).

For instance, `AStarBench --domain msa --file data/msa/papers/BBS11035.tfa
--algo hda --dist HWD --threads 8` aligns the five sequences of `BBS11035` with
HDA$^*$ on 8 workers and the thickness of the rule of Kobayashi et al., and
prints one line of CSV with the cost, the counters of
[§1.3](01-introduction.md#sec-1-3), LB, CO and the seconds.

## 15.9 The benchmark {#sec-15-9}

We compare the work distributions of HDA$^*$ on the alignments of
[`data/msa/papers`][papers], with the costs of PAM250 and the linear gap
penalty $P = 8$ of the experiments of Jinnai and Fukunaga: Zobrist hashing
(Z), Abstract Zobrist hashing (AZ) and AHDA$^*$ (A) with boxes of side 4, and
the hyperplanes (HWD) with the thickness of the rule of `thickness()`, each on
the threads of C++ with 8 and 32 workers, against the sequential A$^*$ on 1
core. All the searches are those of `AStarBench` (`--domain msa --gap 8 --algo
hda
--dist Z | AZ | A | HWD`); the costs are integers, and the OPEN lists are the
bucket queues.

**Table 15.3.** HDA$^*$ with four work distributions on the alignments of
[`data/msa/papers`][papers], with PAM250 and the gap penalty 8. For each
instance and each distribution the table gives the seconds (the median of two
runs), the speedup over sequential A$^*$, the search overhead SO, the
communication overhead CO and the load balance LB, as defined in
[§8.2](08-shared-memory.md#sec-8-2); the machine and the compiler are given in
the note of the table.

--8<-- "tables/msa-instances-hwd.md"

The table is written by [`manual/build/tables.py`][tables.py] from the CSV of
the campaign; until the runs are made on an idle machine, it is a placeholder
that says so. The clock is the field `seconds` of the `Solution`, which does
not include the reading of the file nor the construction of the domain (the
dynamic programs of the pairs, which take $O(\sum_{i<j} L_i L_j)$ time).

## 15.10 Notes and references {#sec-15-10}

The formulation of MSA as a shortest path on the lattice and the pairwise
heuristic are due to Ikeda and Imai [[41]](Z-bibliography.md#ref-41), and the
domain became the benchmark of the searches with a large branching factor with
the Partial Expansion A$^*$ of Yoshizumi, Miura and Ishida
[[42]](Z-bibliography.md#ref-42). Kobayashi, Kishimoto and Watanabe
[[43]](Z-bibliography.md#ref-43) evaluated HDA$^*$ on it and proposed the
hyperplane work distribution, whose bound we restate in Theorem 15.8 (see
Remark 15.9); Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5) compared the
distributions of HDA$^*$ on the same instances, with the gap penalty 8 and the
boxes of side 4 that the domain takes by default. As for the data, the matrix
PAM250 is that of
Dayhoff, Schwartz and Orcutt [[44]](Z-bibliography.md#ref-44), and the
instances are those of BAliBASE 3.0 [[45]](Z-bibliography.md#ref-45). The
claims on [[5]](Z-bibliography.md#ref-5) and [[43]](Z-bibliography.md#ref-43)
in this chapter (the sections, the penalty, the remark and the rule of the
thickness) are those of the comments of the headers, which we have not checked
against the papers.

It is worth noting that the domain is the plainest form of the alignment: the
costs are linear in the gaps (an affine gap penalty would need the state to
remember which sequences are in a gap), the heuristic is the sum of the
pairwise alignments (the heuristics from triples, and the pruning of the
lattice by an upper bound, are not in the library), and the tables of the pairs
are kept whole; the domain is therefore meant for the instances of a few
hundred
residues per sequence of the papers above, not for the alignment of long
sequences. The instance sets of [`data/`][data] are, for each domain, those of
the papers whose searches the library reproduces (a result of the
library can thus be compared with a published one), and their sources and
licenses are recorded next to them.

## 15.11 Exercises {#sec-15-11}

**Exercise 15.1** (`ctest -R PBeFS.msa`). Show that, with PAM250 and $P = 0$,
$c(\mathtt{W}, \mathtt{C}) = 25$ and $c(\mathtt{C}, \mathtt{C}) = 5$, and that
an optimal alignment of `HEAGAWGHEE` and `PAWHEAE` with $P = 8$ costs
$17 \cdot 17 - 2\sigma$, $\sigma$ being the best score with the gap penalty 8
(Proposition 15.2). Check both values against the case.

**Exercise 15.2** (`ctest -R PBeFS.msa`). Prove that the lattice of `test5` has
483840 points and that the expansion of its origin generates 127 successors
(Proposition 15.5). Why does the case compare the cost of A$^*$ with the
dynamic program on the lattice only for the instances of 4 sequences?

**Exercise 15.3** (`ctest -R PBeFS.msa`). For the instance `test3` (4
sequences, $L = 73$) and 64 workers, compute the bound of Theorem 15.8 for
$d = 1, 2, 5, 1/2, 1/4$, and give a state of the lattice of `test3` whose
successors reach $\lceil n/d \rceil + 1 = 2$ workers for $d = 5$. Check that
the case finds the thickness $1/27$ of the rule of Kobayashi et al.

**Exercise 15.4** (`ctest -R PBeFS.msa`). Show that the relation of
`for_each_nblock_successor` without its downward boxes satisfies condition (i)
of [Definition 11.1](11-safepbnf.md#def-11-1) but not condition (ii), and
count the abstract successors of an inner box for $n = 4$ (Proposition 15.6).

**Exercise 15.5** (`ctest -R PBeFS.files`). Read
[`data/15-puzzle/korf_12`][korf_12] by `read_positions` and explain why the
case finds the board of the line 12 of [`korf100.txt`][korf100.txt]; then
explain why a board of [`instances100.txt`][instances100.txt] that
`SlidingTile::parse` refused would make the case fail.

**Exercise 15.6** (`ctest -R PBeFS.board`). Show that freeing a blocked start
changes no optimal cost, and that the life heuristic of a board without blocked
cells and with 4 moves is exact. Check the costs of the five boards of
[`data/grid`][grid] given in section 15.6.

**Exercise 15.7** (`ctest -R PBeFS.dimacs`). Prove Proposition 15.10, and
say which inequalities of the proof hold exactly for the distances computed in
floating point and which may fail by a rounding error; explain why lowering
$\kappa$ by the factor $1 - 10^{-9}$ makes a failure less likely without
excluding it. Check that the case finds $h$ consistent on every arc of the
sample.

<!-- the links to the sources, written by srclinks.py -->

[24-puzzle]: https://github.com/dmeoli/PBeFS/tree/main/data/24-puzzle
[AStarBench.cpp]: https://github.com/dmeoli/PBeFS/blob/main/tools/AStarBench.cpp
[data]: https://github.com/dmeoli/PBeFS/tree/main/data
[dimacs]: https://github.com/dmeoli/PBeFS/tree/main/data/dimacs
[dimacs_sample.py]: https://github.com/dmeoli/PBeFS/blob/main/tools/bench/dimacs_sample.py
[DimacsGraph.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/DimacsGraph.h
[Distribution.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Distribution.h
[Distribution.h:536-591]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Distribution.h#L536-L591
[Distribution.h:555]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Distribution.h#L555
[Domain.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Domain.h
[easy]: https://github.com/dmeoli/PBeFS/tree/main/data/msa/easy
[fetch_dimacs.sh]: https://github.com/dmeoli/PBeFS/blob/main/tools/bench/fetch_dimacs.sh
[files.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/files.cpp
[fromkorf8.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/24-puzzle/fromkorf8.txt
[grid]: https://github.com/dmeoli/PBeFS/tree/main/data/grid
[Grid.h:196]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Grid.h#L196
[instances100.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/24-puzzle/instances100.txt
[instances50.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/24-puzzle/instances50.txt
[korf10.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/rubik/korf10.txt
[korf100.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/15-puzzle/korf100.txt
[korf_12]: https://github.com/dmeoli/PBeFS/blob/main/data/15-puzzle/korf_12
[msa]: https://github.com/dmeoli/PBeFS/tree/main/data/msa
[optimal.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/planning/optimal.txt
[papers]: https://github.com/dmeoli/PBeFS/tree/main/data/msa/papers
[pddl2sas.py]: https://github.com/dmeoli/PBeFS/blob/main/tools/pddl2sas.py
[planning]: https://github.com/dmeoli/PBeFS/tree/main/data/planning
[qaplib]: https://github.com/dmeoli/PBeFS/tree/main/data/qaplib
[robot-nav-2d]: https://github.com/dmeoli/PBeFS/tree/main/data/robot-nav-2d
[SequenceAlignment.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h
[SequenceAlignment.h:230]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L230
[SequenceAlignment.h:284-304]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L284-L304
[SequenceAlignment.h:336]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L336
[SequenceAlignment.h:363]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L363
[SequenceAlignment.h:383]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L383
[SequenceAlignment.h:442]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L442
[SequenceAlignment.h:457]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L457
[SequenceAlignment.h:573]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L573
[SequenceAlignment.h:642]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SequenceAlignment.h#L642
[SlidingTile.h:751]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/SlidingTile.h#L751
[tables.py]: https://github.com/dmeoli/PBeFS/blob/main/manual/build/tables.py
[toh20.txt]: https://github.com/dmeoli/PBeFS/blob/main/data/hanoi/toh20.txt
[TravelingSalesman.h:142]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/TravelingSalesman.h#L142
[tsp]: https://github.com/dmeoli/PBeFS/tree/main/data/tsp
