# 17. Multi-objective search {#ch-17}

[Source: [`include/PBeFS/astar/MultiObjective.h`][MultiObjective.h] (the
concept, the options, the solution, the frontiers, the bucket queue and the
indices of the full vectors), [`include/PBeFS/astar/MOStructures.h`][MOStructures.h]
(the heap and the ordered set of each set of data structures),
[`include/PBeFS/astar/MOAStar.h`][MOAStar.h] (the sequential searches),
[`include/PBeFS/astar/ParallelMOAStar.h`][ParallelMOAStar.h] (the parallel
ones), [`include/PBeFS/astar/MultiObjectiveGraph.h`][MultiObjectiveGraph.h]
and [`include/PBeFS/astar/MultiObjectiveGrid.h`][MultiObjectiveGrid.h] (the
domains); [`tools/MOABench.cpp`][MOABench.cpp]; [`test/moa.cpp`][moa.cpp]]

> **Status — beta.** The code that this chapter describes is in `main`
> but not yet in a release, and its interface may still change.

We are concerned in this chapter with the search of a graph whose edges cost
a vector of $K \geq 2$ nonnegative components (say, the length, the time and
the number of turns of a route). In general no single path is cheapest, and
the answer is therefore the set of the paths from the start to a goal that no
other path improves in all the components at once. The searches of the
previous chapters order their nodes by one number $f = g + h$; here a search
orders paths (not states) by a vector, and for each path that it generates it
has to decide whether a path already known is at least as good in each
component (in which case the new one is useless). The algorithms of the
literature differ mainly in how fast they answer this question (the
*dominance check*), from NAMOA$^*$ of Mandow and Pérez de la Cruz
[[180]](Z-bibliography.md#ref-180) to BOA$^*$ of Hernández et al.
[[181]](Z-bibliography.md#ref-181), EMOA$^*$ of Ren et al.
[[182]](Z-bibliography.md#ref-182) and LTMOA$^*$ of Hernández et al.
[[183]](Z-bibliography.md#ref-183), and the library implements all of them on
one framework. For the parallel case we give the search of Ahmadi et al.
[[184]](Z-bibliography.md#ref-184), whose workers run on different orders of
the objectives and share their bounds, and a search of ours distributed by the
hash of the states, in the manner of HDA$^*$ [[2]](Z-bibliography.md#ref-2)
([Chapter 9](09-hash-distributed.md#ch-9)).

The structure of the chapter is as follows. [§17.1](#sec-17-1) gives the
problem, the notation and the interface. [§17.2](#sec-17-2) gives the
framework of the sequential searches with its proof of correctness, which
holds whatever the structure that answers the dominance checks, and
[§17.3](#sec-17-3) these structures, i.e., the frontiers of BOA$^*$, LTMOA$^*$
and EMOA$^*$, together with the queue of OPEN. In [§17.4](#sec-17-4) we treat
NAMOA$^*$ with dimensionality reduction, whose checks are made at the
generation of a path. The two parallel searches follow: in [§17.5](#sec-17-5)
the searches on several orders of the objectives, and in [§17.6](#sec-17-6)
the hash-distributed search, with the indices that answer its checks on the
full vectors. Finally, [§17.7](#sec-17-7) describes the domains of the chapter
(the road graphs of DIMACS with several costs, and grids whose moves cost
random vectors), and [§17.8](#sec-17-8) to [§17.10](#sec-17-10) give the
benchmark, the notes and the exercises. We use the notation of the
single-objective searches of [Chapter 1](01-introduction.md#ch-1) and the
conventions of the pseudocode of [§1.4](01-introduction.md#sec-1-4).

## 17.1 The problem and its interface {#sec-17-1}

The cost of a path $\pi$ is the vector $c(\pi) \in \mathbb{R}^K_{\geq 0}$,
i.e., the sum of the vectors of its edges. We write $u \preceq v$ ($u$
*weakly dominates* $v$) if $u_i \leq v_i$ for $i = 1, \ldots, K$, $u \prec v$
($u$ *dominates* $v$) if moreover $u \neq v$, and $u <_{lex} v$ for the
lexicographic order. Clearly $u \preceq v$ implies $u \leq_{lex} v$, and $u
\prec v$ implies $u <_{lex} v$: all the searches of the chapter rest on this
property. A path to a goal is *Pareto-optimal* if no path to a goal dominates
it, and the *Pareto front* $\mathcal{F}$ is the set of the costs of the
Pareto-optimal paths. Since several paths may have the same cost, the
searches return a *cost-unique Pareto-optimal set*, i.e., one path for each
vector of $\mathcal{F}$ [[180]](Z-bibliography.md#ref-180). As the costs are
nonnegative, any path to a goal is weakly dominated by a simple one (drop its
cycles), and hence by a vector of $\mathcal{F}$, which is finite. Finally, the
*truncation* $\mathrm{Tr}(v) = (v_2, \ldots, v_K)$ drops the first component
of a vector.

A *label* $x$ is a path from $s_0$ to its state $s(x)$, of cost $g(x)$, with
$f(x) = g(x) + h(s(x))$. The heuristic $h$ is a vector too, and we call it
*consistent* if $h(t) = 0$ at the goals and $h(s) \preceq c(s, t) + h(t)$ for
each edge $(s, t)$, i.e., if each component is a consistent heuristic of the
corresponding objective ([§2.6](02-best-first.md#sec-2-6)). It is then
admissible in each component, and if $Q$ is a path from $s(x)$ to a goal one
has $f(x) \preceq g(x) + c(Q)$. A component equal to infinity marks a state
from which no goal is reachable; the searches do not generate such states,
which loses no path to a goal. Unlike the searches of the previous chapters,
which keep one node per state, these searches keep many labels per state (one
for each cost that is not dominated at that state), and their number usually
makes a multi-objective search much "harder" than its single-objective
version. In fact, already with two objectives the front may have a number of
vectors exponential in the size of the graph (although on the road graphs of
[§17.7](#sec-17-7) it rarely exceeds a few thousand).

A problem is given to these searches as a type `D` that satisfies the concept
`MultiObjectiveDomain` ([`MultiObjective.h:123`][MultiObjective.h:123]), i.e.,
that provides (i) the types `state_type` and `cost_type` (an arithmetic type
for the components) and the constant `objectives`, the number $K \geq 2$,
(ii) the methods `initial_state()`, `is_goal( s )` and `hash( s )` of a
`Domain` ([§1.5](01-introduction.md#sec-1-5)), (iii) the method
`heuristic( s )`, which returns the vector $h(s)$ as a `std::array`, and (iv)
the method `for_each_successor( s , h_s , f )`, which calls `f( t , c , h_t )`
for each successor $t$ of $s$, where $c$ is the vector of the edge and $h_t =
h(t)$. A domain that computes $h$ together with the successors may thus
compute it once. In the header the concept reads as follows:

```cpp
--8<-- "snippets/multi-objective-sig-concept.cpp"
```

The searches return a `ParetoSolution`
([`MultiObjective.h:318`][MultiObjective.h:318]). Its member `front` holds one
`Path` (the vector `cost` and the `states` along it) for each vector of the
front, in lexicographic order, with the costs in the order of the objectives
of the domain whatever the order of the search, and `costs()` gives the
vectors without the paths. Its `status` is `Status::solved` when the front is
complete (`Status::unsolvable` if no goal is reachable, the front being then
empty) and `Status::node_limit` when the search stopped at its limit. Its
member `stats` is a `MOSearchStats`
([`MultiObjective.h:294`][MultiObjective.h:294]), which adds to the counters
of [§2.10](02-best-first.md#sec-2-10) the labels `extracted` from OPEN (pruned
or not), the labels `inserted` into OPEN, the dominance `checks` (the
comparisons of two vectors, or of two components for the frontiers that keep
one component) and the `solutions` found, which a parallel search may find
more than once. The options are a `MOSearchOptions`
([`MultiObjective.h:357`][MultiObjective.h:357]):

```cpp
--8<-- "snippets/multi-objective-sig-options.cpp"
```

i.e., the limit on the labels stored (for each worker, in a parallel search),
the number of states expected (by which the table of the states is sized), the
lazy checks (cf. [§17.2](#sec-17-2)), the order of the objectives in the
lexicographic order of OPEN, and whether the paths are kept besides their
costs.

The header also gives the operations on vectors that the searches share, with
which one can check a front. Among them, `weakly_dominates( u , v , from )`
tests $u_i \leq v_i$ for the components from `from` on (all of them by
default), `dominates( u , v )` tests $u \prec v$, `mo_dead_end( v )` tests
whether a component is infinite, and `truncate( v )` returns $\mathrm{Tr}(v)$;
`pareto_front( v )` returns the vectors of a set that no other one weakly
dominates, each once and in lexicographic order (this is the front by
enumeration, on which the tests of the library rest). In all the searches, the
order of the objectives is a permutation $p$, checked by `objective_order(
order )`, which returns the identity for an empty order and throws
`std::invalid_argument` if the order is not a permutation; `permute( v , p )`
gives $(v_{p_0}, \ldots, v_{p_{K-1}})$ and `unpermute( w , p )` its inverse. A
search in the order $p$ works on the permuted vectors throughout, and
unpermutes the costs of its front at the end; hence, in the rest of the chapter
we may assume the identity order without loss of generality. Also,
[`MOAStar.h`][MOAStar.h] gives `add_costs( u , v )`, the sum by components,
`trace_labels( x )`, the states of the path of a label followed back through
its parents, and `finish_front`, which writes the front and the status of a
`ParetoSolution` from the labels of the solutions found.

## 17.2 The framework of the sequential searches {#sec-17-2}

All the sequential searches of the library but NAMOA$^*$dr are instances of the
framework of Ren et al. [[182]](Z-bibliography.md#ref-182) (its Algorithm 1),
which generalizes BOA$^*$ [[181]](Z-bibliography.md#ref-181) to any number of
objectives. OPEN holds labels, ordered by the lexicographic order of $f$; each
state $s$ has the *frontier* $G^{Tr}(s)$ of the truncated $g$ of the labels
expanded there, and $S$ is the frontier of the truncated $f$ of the solutions
found. A frontier supports two operations: the question whether one of its
vectors weakly dominates a given vector (we then write $v \in \mathrm{dom}\,
F$), and the addition of a vector, which removes those that it dominates. In
the library, the search is the class template `MOAStar< D , F , Q >`
([`MOAStar.h:206`][MOAStar.h:206]), with $F$ the frontier and $Q$ the queue of
OPEN; its pseudocode, which is that of the Doxygen comment of
[`MOAStar.h`][MOAStar.h], is the following.

**Algorithm 17.1** (The framework of the multi-objective searches).

```text
 1  procedure MOA*(s₀)
 2     OPEN ← {x₀}, with s(x₀) = s₀, g(x₀) = 0;  sols ← ∅
 3     while OPEN ≠ ∅ do
 4        x ← the first label of OPEN by the lexicographic order of f;
             remove it from OPEN
 5        if s(x) is not a goal and Tr(g(x)) ∈ dom G(s(x)) then
             continue                          ▷ the frontier check
 6        if Tr(f(x)) ∈ dom S then continue    ▷ the solution check
 7        if s(x) is a goal then
 8           add Tr(f(x)) to S, removing what it dominates;
             add x to sols;  continue
 9        add Tr(g(x)) to G(s(x)), removing what it dominates
10        for each successor t of s(x) with no infinite h(t) do
11           y ← the label of t with g(y) = g(x) + c(s(x), t), parent x
12           if not lazy and (the check of line 5 or 6 prunes y) then
                continue
13           push y in OPEN
14     return sols
```

Ties of $f$ go to the label pushed first, and a label of a goal is not
expanded, since a path through a goal is weakly dominated by its prefix that
ends at that goal. With the option `lazy` the checks of line 12 are not made,
as in the lazy searches of Hernández et al. [[183]](Z-bibliography.md#ref-183):
OPEN then receives more labels, but the search makes fewer checks, since a
label that is pruned anyway at its extraction is checked once instead of
twice. Which of the two variants is faster depends (rather unpredictably) on
the cost of a check against that of a push and a pop. Note that the
truncation at lines 5 and 6 is the
*dimensionality reduction* of Pulido, Mandow and Pérez de la Cruz
[[185]](Z-bibliography.md#ref-185), which the next lemmas justify.

<a id="lem-17-1"></a>**Lemma 17.1** (The order of the extractions). With a
consistent heuristic, Algorithm 17.1 extracts the labels by nondecreasing
$f_1$, and by a nondecreasing $f$ in the lexicographic order.

*Proof.* By consistency $f(y) \succeq f(x)$ for a successor $y$ of $x$, and
hence $f(y) \geq_{lex} f(x)$. A label pushed after the extraction of $x$
descends from labels extracted no earlier than $x$, and therefore its $f$ is
not smaller than $f(x)$; since OPEN gives the smallest of its labels, the claim
follows. $\square$

<a id="lem-17-2"></a>**Lemma 17.2** (Dimensionality reduction). If the label
$y$ was expanded at $s$ before the label $x$ of $s$ is checked, then
$\mathrm{Tr}(g(y)) \preceq \mathrm{Tr}(g(x))$ implies $g(y) \preceq g(x)$;
if the solution $z$ was found before $x$ is checked, then $\mathrm{Tr}(z)
\preceq \mathrm{Tr}(f(x))$ implies $z \preceq f(x)$.

*Proof.* If $x$ is checked when it is extracted, then $f_1(y) \leq f_1(x)$ by
[Lemma 17.1](#lem-17-1), and $g_1(y) \leq g_1(x)$, since $h(s)$ is the same.
If $x$ is checked when it is generated, its parent was extracted after $y$,
and $f_1(y)$ is not larger than the $f_1$ of the parent, which is not larger
than $f_1(x)$. For $z$ the argument is the same, since its $f$ is its cost.
$\square$

We call a label $x$ *useful* for $c \in \mathcal{F}$ if a path $Q$ from $s(x)$
to a goal has $g(x) + c(Q) = c$; then $f(x) \preceq c$, and $h(s(x))$ is
finite. Clearly, the root is useful for each $c \in \mathcal{F}$, and we write
$|Q|$ for the number of edges of the shortest such $Q$.

<a id="thm-17-3"></a>**Theorem 17.3** (Correctness of the framework). With a
consistent heuristic, Algorithm 17.1 terminates and returns a cost-unique
Pareto-optimal set, whatever its frontier and with or without the lazy
checks; when it stops at its node limit, each path that it returns is
Pareto-optimal.

*Proof.* Termination: a label expanded at $s$ is not weakly dominated by the
labels expanded there before ([Lemma 17.2](#lem-17-2), for its truncation).
If its path visited $s$ twice, the prefix that ends at the first visit would
have been expanded at $s$ before it, with a $g$ not larger; hence each label
expanded is a simple path, of which there are finitely many, and each
expansion pushes finitely many labels.

Completeness: let $c \in \mathcal{F}$ never be found, and let $x$ be a label
useful for $c$, among all the labels ever generated, with the least $|Q|$ (the
labels are finitely many). If $|Q| = 0$, then $s(x)$ is a goal and $g(x) = c$,
and $x$ either became a solution or was pruned by $S$, i.e., a solution $z
\preceq f(x) = c$ was found ([Lemma 17.2](#lem-17-2)), and $z = c$ since $c$ is
Pareto-optimal; both contradict the assumption. If $|Q| > 0$, then $x$ was (i)
expanded, and its successor on $Q$, useful with a shorter path, was generated,
against the choice of $x$, or (ii) pruned by $S$, which as above means that $c$
was found (cf. the case $|Q| = 0$), or (iii) pruned by $G^{Tr}(s(x))$. In case
(iii) a label $y$ expanded at $s(x)$ has $g(y) \preceq g(x)$ by [Lemma
17.2](#lem-17-2), and hence $g(y) + c(Q) \preceq c$, i.e., $= c$; thus $y$ is
useful with the same $Q$, and its successor on $Q$, useful with a shorter path,
was generated when $y$ was expanded, again against the choice of $x$. OPEN is
empty at the end, and no other case is left.

Soundness: let the label $z$ of a goal be extracted and not pruned. If a
solution of cost $g(z)$ had been found before, its truncation (or a vector
that weakly dominates it) would be in $S$, and $z$ would be pruned; hence no
cost is found twice. Suppose now that $c \prec g(z)$ for some $c \in
\mathcal{F}$. A label $w$ useful for $c$ has $f(w) \preceq c \prec g(z) =
f(z)$, hence $f(w) <_{lex} f(z)$, and none is in OPEN when $z$ is extracted.
Applied to the labels generated so far, the argument of the completeness
shows that $c$ has been found, and $\mathrm{Tr}(c) \preceq \mathrm{Tr}(f(z))$
makes $S$ prune $z$, a contradiction. Hence each solution is Pareto-optimal;
the argument does not need the search to have ended, and therefore it holds
at the node limit too. $\square$

The proof uses only what the frontiers have in common, i.e., a set of
truncated vectors that weakly dominates each vector added to it, and the
checks of line 12, when they are made, are those of lines 5 and 6 made
earlier; hence the theorem holds for each frontier of [§17.3](#sec-17-3), lazy
or not. Of course, it does not hold with an inconsistent heuristic, for which
[Lemma 17.1](#lem-17-1) fails, and the library does not check consistency.

A program chooses a search by `mo_search( domain , a , options )`
([`MOAStar.h:576`][MOAStar.h:576]), where `a` is a `MOAlgorithm`
([`MOAStar.h:128`][MOAStar.h:128]): `boa` (BOA$^*$, for two objectives only,
with `std::invalid_argument` otherwise), `emoa` (EMOA$^*$, i.e., TOA$^*$ for
three objectives), `ltmoa` (LTMOA$^*$, and LazyLTMOA$^*$ with `lazy`) or
`namoa_dr` (NAMOA$^*$dr, [§17.4](#sec-17-4)). Furthermore,
`mo_algorithm( name )` gives the `MOAlgorithm` of the names `boa`, `emoa`,
`ltmoa` and `namoa-dr`, e.g., of an option of a command line. The following
example, which is compiled and run by the case `manual`, finds the front of a
graph of four nodes with BOA$^*$ and with the hash-distributed search of
[§17.6](#sec-17-6) on two workers:

```cpp
--8<-- "snippets/multi-objective-usage.cpp"
```

Both fronts are $\{(2, 6), (3, 5), (4, 2)\}$. The second vector is the cost of
the path through all four nodes, which has more arcs than the other two and is
Pareto-optimal nonetheless.

## 17.3 The frontiers and the queue of OPEN {#sec-17-3}

We remark that the frontiers trade the cost of a check for the cost of their
structure, and the library has three of them, those of the three algorithms of
the framework. With $K = 2$ the truncated vectors are numbers, and the least of
them, $g_2^{\min}$, weakly dominates them all: `ScalarFrontier`
([`MultiObjective.h:409`][MultiObjective.h:409]) keeps only this number, and
both operations take one comparison, which is the check in constant time of
BOA$^*$ [[181]](Z-bibliography.md#ref-181). For more objectives,
`ArrayFrontier` ([`MultiObjective.h:448`][MultiObjective.h:448]) keeps the
truncated vectors in an array, which a check scans and from which a vector
dominated by a new one is removed by moving the last one into its place, as in
LTMOA$^*$ [[183]](Z-bibliography.md#ref-183). The array is small and
contiguous, which often (but not always) makes its linear scan faster than a
tree, at least for the "small" frontiers of most states. `TreeFrontier`
([`MultiObjective.h:510`][MultiObjective.h:510]) keeps them in a balanced tree
ordered lexicographically, as in EMOA$^*$ [[182]](Z-bibliography.md#ref-182):
since a dominating vector is lexicographically smaller, a check looks only at
the keys $\leq_{lex} \mathrm{Tr}(v)$. For $K = 3$ the keys have two components
and form a staircase (the second decreasing as the first grows), and the check
is answered in $O(\log n)$ by the largest key $\leq_{lex} \mathrm{Tr}(v)$,
which has the least second component among them (TOA$^*$ of
[[182]](Z-bibliography.md#ref-182)); the keys that $\mathrm{Tr}(v)$ dominates
follow it in a run, which an insertion removes. The aliases `AutoFrontier`
(that of BOA$^*$ for two objectives and that of LTMOA$^*$ otherwise, the
default of the searches) and `EMOAFrontier` (that of BOA$^*$ for two objectives
and the tree otherwise, as in [[182]](Z-bibliography.md#ref-182)) give the
frontier by the number of objectives.

Both the tree and the heap of OPEN come from
[`MOStructures.h`][MOStructures.h], which takes them from the set of data
structures chosen ([§3.12](03-data-structures.md#sec-3-12)). `MOHeap` is a
priority queue of plain items, and `LexSet` an ordered set of vectors whose
operations are those of the checks, i.e., (i) whether a key up to a bound
satisfies a predicate (stopping at the first one), (ii) the largest key up to
a bound, and (iii) the removal of the keys from a bound on that satisfy a
predicate, either all of them or only those of the first run. The default set
has the $d$-ary heap and the AVL tree, whose keys up to a bound are visited
from the root down as in the check of [[182]](Z-bibliography.md#ref-182),
while the sets of the standard library and of Boost visit them in order.

By default OPEN is `MOHeap`, and for integer costs it may be `MOBucketQueue`
([`MultiObjective.h:569`][MultiObjective.h:569]), i.e., a bucket for each
value of $f_1$ on a ring that doubles as needed, with each bucket a heap in
the order of OPEN. A queue of this kind is suggested by Ahmadi et al.
[[184]](Z-bibliography.md#ref-184), whose buckets, however, keep the ties of
$f_1$ in the order of their pushes. That order may expand a label that a later
one of the same $f_1$ dominates at its state, which their search repairs at
the solutions; the buckets of the library keep the order of the heap instead.

<a id="prop-17-4"></a>**Proposition 17.4** (The bucket queue). With integer
costs and a consistent heuristic, `MOBucketQueue` returns the labels in the
same order as `MOHeap`, and Algorithm 17.1 then makes the same extractions,
expansions and checks with either.

*Proof.* By [Lemma 17.1](#lem-17-1) no label pushed has an $f_1$ smaller than
that of the last label extracted; hence, the least nonempty bucket holds the
labels of least $f_1$, and its heap returns the least of them in the order of
OPEN, i.e., the lexicographic order of $f$ with the ties to the first pushed.
This is also the label that the heap returns, and since the rest of the search
depends on OPEN only through the order of its extractions, the two runs
coincide. $\square$

The case `moa` checks the proposition on random graphs, where the two queues
give the same counters, and it also checks the traces of the examples of the
papers. On the graph of figure 1 of [[181]](Z-bibliography.md#ref-181)
BOA$^*$ extracts 10 labels and prunes 2 of its 11 successors when it generates
them (tables 1 and 2 there), and on that of figure 1 of
[[183]](Z-bibliography.md#ref-183) LTMOA$^*$ extracts 7 labels and
LazyLTMOA$^*$ 9 (tables 1 and 2 there).

## 17.4 NAMOA$^*$ with dimensionality reduction {#sec-17-4}

NAMOA$^*$ [[180]](Z-bibliography.md#ref-180), the earliest of these searches,
checks a label mainly when it is generated. Besides the frontier $G_{cl}(s)$
of the labels expanded at $s$, each state keeps the set $G_{op}(s)$ of its
labels in OPEN, against which a new label is compared in full, and from which
(and from OPEN) the labels that the new one dominates are removed.
`NAMOAStarDR< D , F >` ([`MOAStar.h:385`][MOAStar.h:385]) is this search with
the dimensionality reduction of NAMOA$^*$dr [[185]](Z-bibliography.md#ref-185),
for any number of objectives; it ignores the option `lazy`, since the checks
at the generation are those that define the search.

**Algorithm 17.2** (NAMOA$^*$dr).

```text
 1  procedure NAMOA*dr(s₀)
 2     OPEN ← {x₀};  G_op(s₀) ← {x₀};  sols ← ∅
 3     while OPEN ≠ ∅ do
 4        x ← the first label of OPEN by the lexicographic order of f;
             remove it from OPEN and from G_op(s(x))
 5        if Tr(f(x)) ∈ dom S then continue    ▷ the solution check
 6        if s(x) is a goal then
 7           add Tr(f(x)) to S;  add x to sols;  continue
 8        add Tr(g(x)) to G_cl(s(x))
 9        for each successor t of s(x) with no infinite h(t) do
10           y ← the label of t with g(y) = g(x) + c(s(x), t), parent x
11           if Tr(g(y)) ∈ dom G_cl(t) or g(y) is weakly dominated by
                a label of G_op(t) or Tr(f(y)) ∈ dom S then continue
12           remove from G_op(t) and from OPEN the labels that y
                dominates
13           push y in OPEN and in G_op(t)
14     return sols
```

A label extracted is never weakly dominated by a label expanded at its state:
a label $z$ expanded there after $y$ was generated was in $G_{op}$ together
with $y$, and the later of the two was checked against the earlier one at
line 11. Hence the frontier check of line 5 of Algorithm 17.1 is not needed
at the extraction, and only the solution check is made there.

<a id="prop-17-5"></a>**Proposition 17.5** (Correctness of NAMOA$^*$dr). With
a consistent heuristic, Algorithm 17.2 terminates and returns a cost-unique
Pareto-optimal set.

*Proof.* The termination and the soundness are those of
[Theorem 17.3](#thm-17-3), whose arguments use only the order of the
extractions and the checks against $G_{cl}$ and $S$; here both checks are made
at the generation, with the same outcome by [Lemma 17.2](#lem-17-2). For the
completeness, let $x$ be, among the useful labels of least $|Q|$, the one
generated last. If it was removed from OPEN by a later label, or pruned by a
label of $G_{op}$ that is removed later by another one, that label is useful
with the same $|Q|$ and was generated later, against the choice of $x$. If it
was pruned by a label of $G_{op}$ that is expanded, or by $G_{cl}$, the
successor of that label on $Q$ is useful with a shorter path. The other cases
are those of [Theorem 17.3](#thm-17-3). $\square$

In the implementation, the sets $G_{op}$ are arrays, and a label removed from
OPEN is marked dead and skipped at its extraction instead of being taken out of
the heap; this keeps the heap free of deletions, at the price of its size.

## 17.5 Searches on several orders of the objectives {#sec-17-5}

The first parallel search, `MOAStarPar< D , F , Executor , Q >`
([`ParallelMOAStar.h:199`][ParallelMOAStar.h:199]), is the parallel MOA$^*$
of Ahmadi et al. [[184]](Z-bibliography.md#ref-184). Each worker $j$ runs
Algorithm 17.1 with its own order of the objectives, whose first one is its
*primary objective* $i_j$, and with separate OPEN, frontiers and solutions;
thus the workers search the same graph from different "sides" of the front.
`objective_orders< K >( n )` gives the first $n$ orders, i.e., the $K$ cyclic
ones first, as in the paper, and then the other permutations in lexicographic
order, up to $K!$; the workers are as many as the threads (hence at most
$K!$). The options are a `ParallelMOOptions`
([`ParallelMOAStar.h:135`][ParallelMOAStar.h:135]), i.e., the options of the
sequential searches (of which `MOAStarPar` ignores `order`) and the parameters
of the parallel ones ([§1.6](01-introduction.md#sec-1-6)), with three more
fields:

```cpp
--8<-- "snippets/multi-objective-sig-parallel-options.cpp"
```

The workers share two vectors: $F_i$ is the largest $f_i$ that a worker of
primary objective $i$ has extracted (with several workers on one primary
objective, the largest of theirs), and $U_i$ is an upper bound on the $i$-th
component of the vectors of the front that are still to be found. A worker
keeps the costs of its solutions, in its order, in the sets $V_1, \ldots,
V_{K-1}$, where $V_\lambda$ holds vectors of which only the components after
the $\lambda$-th (in the order of the worker) are compared; hence $V_1$ is the
frontier $S$ of Algorithm 17.1. At line 11 of Algorithm 17.3 the components
are numbered from 0, as in the code, and the component $\lambda$ there is
therefore the first one that $V_\lambda$ compares.

**Algorithm 17.3** (MOA$^*$ on several orders of the objectives).

```text
 1  procedure MOA*-par(s₀)
 2     F ← 0;  U ← ∞
 3     in parallel, on each worker j with its order (i_j first) do
 4        MOA*(s₀) in the order of j, where
 5           extracting x:  F_{i_j} ← max(F_{i_j}, f_{i_j}(x));
                if f_{i_j}(x) ≥ U_{i_j} then stop all the workers
 6           the solution check of x (and of a successor, if not lazy):
                f_i(x) ≥ U_i for some i, or the components from λ on
                of f(x) weakly dominated by a vector of V_λ, some λ
 7           a solution z:  add z to V₁;  UPDATE
 8     return the union of the solutions of the workers, each cost once
 9  procedure UPDATE                              ▷ on the worker j
10     for λ = 1, …, K − 2 do
11        i ← the objective of the component λ in the order of j
12        for each v ∈ V_λ with v_i ≤ F_i do
13           move v from V_λ to V_{λ+1}, unless a vector there
                dominates it, removing those that it dominates
14     if V_{K−1} ≠ ∅ then U_i ← min(U_i, the least component of
           V_{K−1}), i the last objective in the order of j
```

In the implementation, the values of $F$ and $U$ are written with release and
read with acquire; hence, what a worker did before writing a value happens
before what another one does after reading it, and the proofs below reason on
this order. We call a vector of $V_\lambda$, kept for a solution $z$,
*justified* if, for each objective $m$ whose component it no longer compares,
$z_m \leq F_m$ held for a value of $F_m$ read by the worker, i.e., some worker
of primary objective $m$ had extracted before a label $x'$ with $f_m(x') \geq
z_m$. The bound $U_m = z_m$ is the vector truncated in all the components but
$m$.

<a id="lem-17-6"></a>**Lemma 17.6** (What an extraction certifies). Let a
worker of primary objective $i$ extract a label $x$. Then each $c \in
\mathcal{F}$ with $c_i < f_i(x)$ is found by some worker by the end of the
search.

*Proof.* By induction on the order of the events, which is a strict partial
order on the finitely many events of a run. Let $c$ be such a vector, not
found, and choose the useful label $w$ for $c$ of least $|Q|$ among those that
the worker generated before extracting $x$. Since $f_i(w) \leq c_i <
f_i(x)$, $w$ is not in OPEN when $x$ is extracted ([Lemma 17.1](#lem-17-1), in
the order of the worker), and the cases of [Theorem 17.3](#thm-17-3) apply,
with one more: $w$ (or a goal label $w$ of cost $c$) was pruned by a
justified vector of a solution $z$, i.e., the compared components of $z$ are
not larger than those of $f(w) \preceq c$. If also $z_m \leq c_m$ for each
component $m$ that is no longer compared, then $z \preceq c$, and $c = z$ is
found. Otherwise $c_m < z_m \leq f_m(x')$ for some $m$ and some label $x'$
extracted by a worker of primary objective $m$ before the pruning, and hence
before the extraction of $x$; by the induction, $c$ is found. $\square$

<a id="thm-17-7"></a>**Theorem 17.7** (Correctness of the search on several
orders). With a consistent heuristic, Algorithm 17.3 returns a cost-unique
Pareto-optimal set: each vector of $\mathcal{F}$ is found by at least one
worker, each solution found is Pareto-optimal, and the copies found by
several workers are merged. With a node limit, each path that it returns is
Pareto-optimal.

*Proof.* Each worker terminates, as Algorithm 17.1 does. Let the worker $j$ end
first. If its OPEN is empty, the argument of [Lemma 17.6](#lem-17-6) applies to
all the labels it generated, and all the vectors of $\mathcal{F}$ are found; if
it stopped at a label with $f_{i_j} \geq U_{i_j}$, all the labels left in its
OPEN have $f_{i_j} \geq U_{i_j}$ and are pruned by that bound, which is
justified, and the same argument applies. For the soundness, let a worker add
the solution $z$ with $c \prec z$ for some $c \in \mathcal{F}$. As in [Theorem
17.3](#thm-17-3), no label useful for $c$ is in its OPEN, and therefore a
useful label $w$ of the worker was pruned by its frontier of a state (and then
a useful label with a shorter path exists) or by its solutions. Furthermore,
the rules of pruning are monotone, i.e., a vector weakly dominated by a pruned
one is pruned too at any later time: in fact, the sets $V_\lambda$ only grow in
power (a vector moves to a set that compares fewer components, or is replaced
by one that dominates it), and the bounds only decrease. Hence $z \succeq c
\succeq f(w)$ would have been pruned too. Within a worker no cost is found
twice, by its $V_1$, and the merge removes the copies of the costs found by
several workers. $\square$

Two choices of the implementation are ours. First, the sets $V_\lambda$ are
updated at each solution, as in [[184]](Z-bibliography.md#ref-184), and also
every `flush_period` extractions when $F$ has grown since, which justifies the
truncations earlier. Second, the check of a label against the solutions
starts from $V_{K-1}$, the smallest and strongest set, which saves
comparisons and changes nothing else. With `share_bounds` false no bound is
shared, which is the variant without upper bounds of
[[184]](Z-bibliography.md#ref-184). The frontier `F` may be any of
[§17.3](#sec-17-3): with the array (the default for $K \geq 3$) the search is
LTMOA$^*$-par, and with two objectives and the frontier of BOA$^*$ its two
workers are BOA$^*$ in the orders $(1, 2)$ and $(2, 1)$, sharing their bounds.
A worker beyond the first $K$ adds an order whose primary objective is already
taken, and it shares the $F_i$ of that objective.

## 17.6 Hash-distributed multi-objective search {#sec-17-6}

The second parallel search, `HDMOAStar< D , Executor >`
([`ParallelMOAStar.h:573`][ParallelMOAStar.h:573]), divides the states
instead of the orders. Each state belongs to the worker to which a function of
its hash assigns it, as in HDA$^*$ ([§9.2](09-hash-distributed.md#sec-9-2)),
and this owner keeps the labels of its states in its OPEN, by the
lexicographic order of $f$, together with their frontiers; the labels
generated by a worker travel in batches to the owner of their last state,
through the mailbox of [§9.3](09-hash-distributed.md#sec-9-3). A solution is
kept by the owner of its goal, which sends it to the other workers, and each
worker prunes with its copy of the solutions known.

Since the workers do not extract the labels in one global order, neither
[Lemma 17.1](#lem-17-1) nor [Lemma 17.2](#lem-17-2) holds across the workers: a
label expanded may later be dominated by a label that arrives at the same
state, and a solution may be dominated by a later one. Hence the frontier of a
state keeps the full $g$ of the labels expanded there, a label is pruned if its
$g$ is weakly dominated by that frontier or its $f$ by a known solution, and
the owner of a goal removes from its solutions those that a new one dominates.
This is a *label-correcting* search, in which a label may be expanded and
dominated later, while the searches of the previous sections are
label-setting (i.e., a label expanded is final).

**Algorithm 17.4** (The worker of the hash-distributed search).

```text
 1  procedure WORKER(i)              ▷ OPEN_i, the frontiers G, sols_i
 2     loop
 3        for each message m received from a worker j do
 4           if m is a solution z then add z to the copy of the
                solutions;  continue
 5           y ← the label of m;  received[j] ← received[j] + 1
 6           if not lazy and (g(y) is weakly dominated by G(s(y))
                or f(y) by a solution) then continue
 7           push y in OPEN_i
 8        publish received[ ];  L_i ← min(f₁ of the top of OPEN_i,
             the least f₁ sent and not yet received);  publish L_i
 9        if OPEN_i = ∅ or f₁(top) > min_k L_k + w then
10           flush the outboxes;  wait for a message or the end;
             continue
11        x ← the top of OPEN_i, removed
12        if g(x) weakly dominated by G(s(x)) or f(x) by a solution
             then continue                     ▷ checked again
13        if s(x) is a goal then
14           remove from sols_i the solutions that x dominates;
             add x to sols_i;  send x to every worker;  continue
15        add g(x) to G(s(x))
16        for each successor t of s(x) with no infinite h(t) do
17           y ← the label of t, parent x
18           if not lazy and f(y) is weakly dominated by a solution
                then continue
19           if t is owned by i then push y in OPEN_i, as at lines 6-7
20           else put y in the outbox of the owner of t;  sent ← sent + 1
```

A worker never blocks: when it has no label, or its top is outside the
window, it flushes its outboxes, pauses briefly and starts the loop again, and
hence it polls its inbox and publishes its values also while it waits. The
window decides the order of the extractions. A worker expands a label only if
its $f_1$ is at most $L + w$ (line 9), where $L = \min_k L_k$ is the least
$f_1$ of the tops of the OPEN lists and of the labels still travelling, as
published by the workers, and $w$ is the option `window` (0 by default,
negative for none); thus the workers advance together by $f_1$, close to the
order of the sequential search. With `acknowledge` (the default) a worker
counts the labels it has sent to each owner, each owner counts the labels it
has received from each worker, and the least $f_1$ of the labels sent stays in
the $L_i$ of the sender until all of them have arrived. Without it, the sender
forgets that $f_1$ as soon as its outboxes are pushed, and $L$ may then exceed
the $f_1$ of a label still travelling, which the window then fails to hold
back. The mailbox detects the end when no worker has a label and no message
travels ([§9.3](09-hash-distributed.md#sec-9-3)), and the front is then the
union of the solutions of the owners, each cost once.

<a id="lem-17-8"></a>**Lemma 17.8** (The window does not stall). Assume that
the values published by a worker become visible to the others in finite time.
Then, with $w \geq 0$, with or without `acknowledge`, as long as some label is
in an OPEN list or travelling, some worker eventually extracts a label or
receives one.

*Proof.* Suppose that from some moment on no worker extracts or receives a
label. A label that travels is in an outbox that its sender flushes at its
next wait, and its owner receives it at its next poll; hence, after a while,
no label travels, and all the labels are in OPEN lists that no longer change.
Each sender then sees (with `acknowledge`) that all its labels have arrived,
or has forgotten them at its flush (without it), and so publishes as $L_k$ the
$f_1$ of its top ($\infty$ if its OPEN is empty). Let $L^*$ be the least $f_1$
of the tops, and $k$ a worker whose top has $f_1 = L^*$. Once the values
above are visible, $k$ reads $\min_k L_k = L^*$, its top satisfies $L^* \leq
L^* + w$, and it extracts it, a contradiction. $\square$

<a id="thm-17-9"></a>**Theorem 17.9** (Correctness of the hash-distributed
search). With a consistent heuristic, without a node limit and under the
assumption of [Lemma 17.8](#lem-17-8), Algorithm 17.4 terminates and returns
a cost-unique Pareto-optimal set, whatever the order in which its workers
extract their labels.

*Proof.* Termination: a label is expanded at $s$ only if no label expanded
there before weakly dominates it, and therefore, as in [Theorem
17.3](#thm-17-3), each label expanded is a simple path. Hence the labels and
the messages are finitely many, by [Lemma 17.8](#lem-17-8) they are all
extracted or received eventually, and the mailbox ends the search when no
worker has a label and no message travels.

Completeness: let $c \in \mathcal{F}$, and choose the useful label $x$ of
least $|Q|$ among all the labels created. If it was expanded, its successor on
$Q$ was sent to its owner, which either made a useful label of shorter path,
against the choice of $x$, or pruned it. In the latter case the pruning was
either by the frontier of the state, whose label $y$ with $g(y) \preceq$ that
of the successor is useful and was expanded (and then the successor of $y$ on
$Q$ is useful with a path shorter still), or by a solution $z \preceq c$,
i.e., $z = c$. The same holds if $x$ itself was pruned; if $s(x)$ is a goal
with $g(x) = c$, then $x$ became a solution unless a solution of cost $c$ was
there already. A solution is removed only by one that weakly dominates it,
which for $c$ has cost $c$; hence, for each $c \in \mathcal{F}$, the owner of
some goal ends with a solution of cost $c$.

Soundness: the solutions of each owner are mutually non-dominated, and a
solution dominated by a vector of $\mathcal{F}$ found by another owner is
removed by the final filter of the union. $\square$

By [Theorem 17.9](#thm-17-9), a value of $L$ read late may change the order of
the expansions, but it leaves the front unchanged. With a node limit the
solutions are mutually non-dominated but need not be Pareto-optimal, as in any
label-correcting search, and the status `Status::node_limit` marks this case.

In this search the frontiers of the states and the copies of the solutions
keep full vectors,
whose checks are dearer than the truncated ones of the sequential search. The
structures that answer them rest on the following remark.

<a id="lem-17-10"></a>**Lemma 17.10** (Truncation above a threshold). Let a set
$Z$ of vectors be kept truncated in a frontier of EMOA$^*$, a truncation being
added only if no key dominates it, and let $T = \max_{z \in Z} z_1$. For a
vector $x$ with $x_1 \geq T$, a key weakly dominates $\mathrm{Tr}(x)$ if and
only if some $z \in Z$ weakly dominates $x$.

*Proof.* A key $\mathrm{Tr}(z) \preceq \mathrm{Tr}(x)$ comes from some $z \in
Z$ with $z_1 \leq T \leq x_1$, hence $z \preceq x$. Conversely, if $z \preceq
x$, then either the truncation of $z$ was added, or a key that dominates it
was there; since a key is removed only by one that dominates it, some key
dominates $\mathrm{Tr}(z) \preceq \mathrm{Tr}(x)$. $\square$

Hence we keep each set of vectors twice. Truncated, as in [Lemma
17.10](#lem-17-10), it answers in one check of the sequential search (one
comparison for two and three objectives) the vectors with $x_1 \geq T$, i.e.,
those of the labels that arrive at a state in the order of $f_1$, which the
window makes the common case (cf. line 9 of Algorithm 17.4). In full, it
answers the others. For two objectives we keep the full vectors in the tree of
EMOA$^*$ of the vectors $(0, v)$, whose truncations are the vectors $v$ and
form a staircase, and a check takes $O(\log n)$. For three, we keep them in a
`LayerIndex` ([`MultiObjective.h:655`][MultiObjective.h:655]), which keeps the
vectors by increasing first component in layers of a few of them, with a
Fenwick tree [[189]](Z-bibliography.md#ref-189) over the layers; each node of
the tree keeps the vectors of the layers it covers truncated as in [Lemma
17.10](#lem-17-10). A check of $x$ queries the $O(\log L)$ nodes that cover the
layers before that of $x$, whose vectors have a first component at most $x_1$
(and therefore each answer is exact by [Lemma 17.10](#lem-17-10)), and then
compares $x$ with the vectors of its own layer whose first component is not
larger; these are the layers of $f_1$ of the indexed pruning of Xu et al.
[[187]](Z-bibliography.md#ref-187), here on prefixes. For more objectives a
node of the Fenwick tree costs too much to keep, and a `BlockIndex`
([`MultiObjective.h:787`][MultiObjective.h:787]) keeps the vectors in blocks by
increasing first component, with a segment tree of the minima of the blocks,
component by component. A check descends into a node of this tree only if its
minimum weakly dominates $x$ (the block-minimum index of
[[187]](Z-bibliography.md#ref-187)), and an insertion compares no vector.

<a id="prop-17-11"></a>**Proposition 17.11** (The indices are exact). For
each $x$, `LayerIndex` and `BlockIndex` report that $x$ is weakly dominated
if and only if some vector added to them weakly dominates $x$.

*Proof.* Note that neither removes a vector; hence, the question is about all
the vectors added. For `LayerIndex`, the vectors of the layers before that of
$x$ have a first component at most $x_1$, and each node of the Fenwick tree
answers exactly for its layers by [Lemma 17.10](#lem-17-10), with $T$ the
largest first component of the layers it covers; the vectors of the layer of
$x$ are compared in full, and the layers after it hold only vectors of first
component larger than $x_1$, which cannot dominate $x$. For `BlockIndex`, a
block or a node of the segment tree whose minimum does not weakly dominate $x$
holds no vector that does, since all its vectors are weakly dominated by the
minimum; the search descends into all the others, and compares $x$ with the
vectors of the blocks it reaches. $\square$

Keeping the dominated vectors does not change the answers, since a vector
weakly dominated by another vector of the index dominates no vector that the
other does not dominate as well. It costs memory, though, and `LayerIndex`
rebuilds its layers when one of them grows past four times its nominal size,
which keeps the scan of a layer short.

## 17.7 The domains {#sec-17-7}

`MultiObjectiveGraph< Cost , K >`
([`MultiObjectiveGraph.h:152`][MultiObjectiveGraph.h:152]) is a directed graph
whose arcs cost a vector, given by `add_arc( u , v , c )`, and which is
searched from a source to a target. Its heuristic is set by
`set_heuristic( u , h )` or by `set_ideal_heuristic()`, which computes the
*ideal point* of each node, i.e., $h_i(u) = d_i(u, t)$, the cost of a
cheapest path from $u$ to the target by the $i$-th objective only, by one
backward Dijkstra's algorithm per objective (`distances_to`,
[`MultiObjectiveGraph.h:106`][MultiObjectiveGraph.h:106]), as in
[[181]](Z-bibliography.md#ref-181). This heuristic is consistent by the
triangle inequality of each $d_i$, and it is infinite at the nodes from which
the target is not reachable, which the searches never generate.

`read_dimacs_objectives< Cost , K >( in , extra , source , target , seed )`
([`MultiObjectiveGraph.h:310`][MultiObjectiveGraph.h:310]) builds such a graph
from the road graphs of the 9th DIMACS Implementation Challenge
[[64]](Z-bibliography.md#ref-64), one file per objective with the same arcs in
the same order (the distances of `USA-road-d` and the times of `USA-road-t`).
It appends the objectives `extra`, computed from the graph as the papers do:
`ExtraCost::arcs` (1 for each arc, as in [[188]](Z-bibliography.md#ref-188)),
`ExtraCost::degree` (2 for the arcs whose end nodes have on average 4
neighbours or more and 1 otherwise, as in [[182]](Z-bibliography.md#ref-182)),
`ExtraCost::avg_degree` (the average number of neighbours of the two end
nodes, rounded down, as in [[184]](Z-bibliography.md#ref-184)) and
`ExtraCost::random` (an integer drawn uniformly in $[1, 100]$, arc by arc, by
a generator seeded with `seed`, as in [[183]](Z-bibliography.md#ref-183));
`extra_cost( name )` gives the `ExtraCost` of the names `arcs`, `degree`,
`avg-degree` and `random`. The heuristic is left to 0, and the caller sets the
ideal point once the source and the target are known. The script
[`tools/bench/fetch_dimacs.sh`][fetch_dimacs.sh] downloads the graphs (both
files of a map with the option `--times`).

`MultiObjectiveGrid< Cost , K >`
([`MultiObjectiveGrid.h:80`][MultiObjectiveGrid.h:80]) is pathfinding on the
free cells of a `GridMap` ([Chapter 14](14-domains.md#ch-14)) with 4 or 8
moves, where a diagonal move is allowed only if both the cells beside it are
free. The $i$-th cost of the move between two cells is an integer drawn
uniformly in $[lo, hi]$ by a hash of the unordered pair of cells, of $i$ and of
a seed; hence a move costs as much as its inverse, and the grid needs no
memory for its costs. With 4 moves and costs in $[1, 10]$ these are the grids
of [[182]](Z-bibliography.md#ref-182), which uses the empty grid of $10 \times
10$ cells and maps of the benchmarks of Sturtevant
[[15]](Z-bibliography.md#ref-15). The heuristic is $h_i(u) = lo \cdot
\delta(u, t)$, where $\delta$ is the Manhattan distance with 4 moves and the
Chebyshev distance with 8: a move lowers $\delta$ by 1 at most and costs at
least $lo$ in each component, and hence $h$ is consistent;
`set_ideal_heuristic()` replaces it with the ideal point.

The benchmark [`tools/MOABench.cpp`][MOABench.cpp] runs any of the searches of
the chapter on these domains, i.e., on a road graph with its queries or on a
grid with its points, and writes one line of CSV per run with the size of the
front, the counters of `MOSearchStats` and, for the parallel searches, the
labels sent, the messages and the most labels expanded by one worker.

## 17.8 Benchmark {#sec-17-8}

We measure the searches by [`tools/MOABench.cpp`][MOABench.cpp] on the road
graph of New York of [[64]](Z-bibliography.md#ref-64) with the distances, the
times and, for three objectives and more, the extra objectives of
[§17.7](#sec-17-7), on the queries of the papers. Besides the time, the
benchmark reports counts that do not depend on the machine (the labels
extracted and expanded, the dominance checks, the most labels expanded by one
worker and the messages). By these counts we compare the searches with the
codes of the papers, whose fronts and expansions they reproduce when the order
of the ties is the same.

**Table 17.1.** The multi-objective searches on the road graph of New York. A
row is a search, with its number of workers; the columns give the number of
objectives, the queries solved by all the searches within the limit, the sum
of the sizes of their fronts, the labels expanded, the most expanded by one
worker, the dominance checks and the time in seconds, summed over the queries
solved by all.

--8<-- "tables/multi-objective-roads.md"

For now, the table is a placeholder until it is measured on an idle machine.

## 17.9 Notes and references {#sec-17-9}

NAMOA$^*$ of Mandow and Pérez de la Cruz
[[180]](Z-bibliography.md#ref-180) is the multi-objective A$^*$ with
consistent heuristics in which a label expanded is final, which is the
setting of the chapter but [§17.6](#sec-17-6). The dimensionality reduction
of Pulido, Mandow and Pérez de la Cruz [[185]](Z-bibliography.md#ref-185) is
the observation behind [Lemma 17.2](#lem-17-2), and BOA$^*$ of Hernández et
al. [[181]](Z-bibliography.md#ref-181) its use for two objectives, with the
checks in constant time. Ren et al. [[182]](Z-bibliography.md#ref-182) wrote
the framework of Algorithm 17.1 and EMOA$^*$, with the balanced tree and the
check in logarithmic time for three objectives, and Hernández et al.
[[183]](Z-bibliography.md#ref-183) the linear frontier of LTMOA$^*$ and the
lazy checks, with the measure of the checks that the counter `checks`
follows. The targeted multiobjective Dijkstra algorithm of Maristany de las
Casas et al. [[188]](Z-bibliography.md#ref-188) is a label-setting search of
the same family, from which we take the objective that counts the arcs.

For the parallel searches, Sanders and Mandow
[[186]](Z-bibliography.md#ref-186) expand at once all the labels whose $f$ is
not dominated in OPEN, which keeps the label-setting property at the price of
a global synchronization per round; the searches of this chapter have none.
Ahmadi et al. [[184]](Z-bibliography.md#ref-184) proposed the search of
[§17.5](#sec-17-5), whose argument we write in full in
[Lemma 17.6](#lem-17-6) and [Theorem 17.7](#thm-17-7). The hash-distributed
search of [§17.6](#sec-17-6) is close to SIP-MOSP of Xu et al.
[[187]](Z-bibliography.md#ref-187), whose update workers own the frontiers of
the vertices assigned to them; here the owner of a state both keeps its
frontier and expands its labels, as in HDA$^*$ [[2]](Z-bibliography.md#ref-2),
and therefore no frontier is written by two workers and no worker waits for
another. The window of $f_1$, its acknowledgement and the two indices of the
full vectors (the Fenwick tree of Fenwick [[189]](Z-bibliography.md#ref-189)
over the layers, and the minima of the blocks) are the parts of that search
that make its checks on full vectors affordable. In our runs on grids and road
graphs the search still made several times the checks of EMOA$^*$, while it
divided the expansions among the workers rather evenly.

## 17.10 Exercises {#sec-17-10}

**Exercise 17.1** (`ctest -R PBeFS.moa`). Run Algorithm 17.1 by hand, with
the frontier of BOA$^*$, on the graph of figure 1 of
[[181]](Z-bibliography.md#ref-181), and check that it extracts 10 labels and
prunes 2 successors at their generation; then run it lazily, and explain why
the number of extractions grows while that of the checks does not.

**Exercise 17.2** (`ctest -R PBeFS.moa`). Show by an example of three nodes
that [Lemma 17.2](#lem-17-2) fails with an inconsistent heuristic, and that
Algorithm 17.1 may then return a dominated path; the test checks each search
against the front computed by `pareto_front` over all the simple paths of
random graphs, with consistent heuristics only.

**Exercise 17.3** (`ctest -R PBeFS.moa`). Prove that the search in the order
$p$ returns the same front as the search in the identity order, and that its
expansions may differ; the test runs all the orders of the objectives.

**Exercise 17.4** (`ctest -R PBeFS.moa`). Prove that for $K = 3$ the keys of
a `TreeFrontier` form a staircase, i.e., that the second component decreases
as the first grows, and that the largest key $\leq_{lex} \mathrm{Tr}(v)$
therefore suffices to answer the check of $v$.

**Exercise 17.5** (`ctest -R PBeFS.moa-thread`). In Algorithm 17.3 with two
workers on two objectives, show that the worker of primary objective 1 stops
at the first label it extracts with $f_1 \geq U_1$, and that by then all the
vectors of the front have been found by one of the two workers
([Lemma 17.6](#lem-17-6)).

**Exercise 17.6** (`ctest -R PBeFS.moa-thread`). Give a run of
Algorithm 17.4 with two workers and no window in which a label is expanded and
later dominated at its state by a label that arrives from the other worker.
Then show that, if all the workers read the current values $L_k$, with $w = 0$
and `acknowledge` a label expanded could later be dominated only by a label of
the same $f_1$, and explain why values read late weaken this; the test runs
the search with and without the acknowledgement and checks the front in both
cases.

**Exercise 17.7** (`ctest -R PBeFS.moa`). Show that
[Lemma 17.10](#lem-17-10) fails for a vector with $x_1 < T$, by an example
with two vectors in $Z$, and explain why `LayerIndex` compares the vectors of
the layer of $x$ in full.

**Exercise 17.8** (`ctest -R PBeFS.manual`). Change the cost of the arc from
2 to 3 in the example of [§17.2](#sec-17-2) to $(1, 4)$, and find the new
front by hand; then check it with the example.

<!-- the links to the sources, written by srclinks.py -->

[fetch_dimacs.sh]: https://github.com/dmeoli/PBeFS/blob/main/tools/bench/fetch_dimacs.sh
[moa.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/moa.cpp
[MOABench.cpp]: https://github.com/dmeoli/PBeFS/blob/main/tools/MOABench.cpp
[MOAStar.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOAStar.h
[MOAStar.h:128]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOAStar.h#L128
[MOAStar.h:206]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOAStar.h#L206
[MOAStar.h:385]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOAStar.h#L385
[MOAStar.h:576]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOAStar.h#L576
[MOStructures.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MOStructures.h
[MultiObjective.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h
[MultiObjective.h:123]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L123
[MultiObjective.h:294]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L294
[MultiObjective.h:318]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L318
[MultiObjective.h:357]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L357
[MultiObjective.h:409]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L409
[MultiObjective.h:448]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L448
[MultiObjective.h:510]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L510
[MultiObjective.h:569]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L569
[MultiObjective.h:655]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L655
[MultiObjective.h:787]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjective.h#L787
[MultiObjectiveGraph.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGraph.h
[MultiObjectiveGraph.h:106]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGraph.h#L106
[MultiObjectiveGraph.h:152]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGraph.h#L152
[MultiObjectiveGraph.h:310]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGraph.h#L310
[MultiObjectiveGrid.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGrid.h
[MultiObjectiveGrid.h:80]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/MultiObjectiveGrid.h#L80
[ParallelMOAStar.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelMOAStar.h
[ParallelMOAStar.h:135]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelMOAStar.h#L135
[ParallelMOAStar.h:199]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelMOAStar.h#L199
[ParallelMOAStar.h:573]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/ParallelMOAStar.h#L573
