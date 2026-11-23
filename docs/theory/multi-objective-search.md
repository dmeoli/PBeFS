# Multi-objective search

We are concerned here with the searches of a `MultiObjectiveDomain`, whose
edges cost a vector of $K \geq 2$ nonnegative components and whose answer is
the Pareto front of the paths from the start $s_0$ to a goal, instead of one
path. This page fixes the notation and proves that the sequential framework of
[MOAStar.h](../../include/PBeFS/astar/MOAStar.h) returns the front with any of
its frontiers, lazy or not, and with NAMOA\*dr. It then proves the same for the
two parallel searches of
[ParallelMOAStar.h](../../include/PBeFS/astar/ParallelMOAStar.h), i.e., the
searches on several orders of the objectives of Ahmadi et al.
[[184]](../../README.md#ref-184), whose argument we write in full, and the
search distributed by the hash of the states, a label-correcting search of
ours. The pseudocode of each search is in its header, and the notation of the
single-objective searches is that of [Best-first search and
A*](best-first-search.md).

## Notation

The cost of a path $\pi$ is the vector $c(\pi) \in \mathbb{R}^K_{\geq 0}$, the
sum of the vectors of its edges. We write $u \preceq v$ ($u$ weakly dominates
$v$) if $u_i \leq v_i$ for each $i$, $u \prec v$ if moreover $u \neq v$, and
$u \lt_{lex} v$ for the lexicographic order; clearly $u \preceq v$ implies
$u \leq_{lex} v$, and $u \prec v$ implies $u \lt_{lex} v$. A path to a goal is
Pareto-optimal if no path to a goal dominates it, and the Pareto front
$\mathcal{F}$ is the set of the costs of the Pareto-optimal paths. Our searches
return a cost-unique Pareto-optimal set, i.e., one path for each vector of
$\mathcal{F}$ [[180]](../../README.md#ref-180). Since the costs are
nonnegative, any path to a goal is weakly dominated by a simple one (drop its
cycles), and hence by a vector of $\mathcal{F}$, which is finite. Finally, the
truncation $\mathrm{Tr}(v) = (v_2, \ldots, v_K)$ drops the first component and
$\mathrm{Tr}^\lambda$ the first $\lambda$ ones. A label $x$ is a path from
$s_0$ to its state $s(x)$, of cost $g(x)$, with $f(x) = g(x) + h(s(x))$. We
assume that the heuristic is consistent, i.e., $h(t) = 0$ at the goals and
$h(s) \preceq c(s, t) + h(t)$ for each edge; hence, each component is a
consistent heuristic of its objective, and therefore admissible: if $Q$ is a
path from $s(x)$ to a goal, $f(x) \preceq g(x) + c(Q)$. A component equal to
infinity marks a state from which no goal is reachable, and such states are
never generated, which loses no path. On the road graphs the heuristic is the
ideal point, $h_i(u) = d_i(u, t)$, the distance of $u$ from the target by the
$i$-th objective only, computed by $K$ backward runs of Dijkstra's algorithm,
as in [[181]](../../README.md#ref-181); it is consistent by the triangle
inequality of each $d_i$. On the grids of
[MultiObjectiveGrid.h](../../include/PBeFS/astar/MultiObjectiveGrid.h) the
default is $h_i(u) = \ell\,\delta(u, t)$, $\ell$ being the least cost of a move
and $\delta$ the Manhattan (4 moves) or the Chebyshev (8 moves) distance, which
a move lowers by 1 at most, while it costs at least $\ell$ in each component.

## The sequential framework

`MOAStar` extracts the labels by the lexicographic order of $f$ (the ties to
the label pushed first). It keeps for each state $s$ the frontier $G(s)$ of the
truncated $g$ of the labels expanded there, and the frontier $S$ of the
truncated $f$ of the solutions. A label is pruned when extracted (and, unless
the search is lazy, when generated) if $\mathrm{Tr}(g(x))$ is weakly dominated
by $G(s(x))$, $s(x)$ not being a goal, or if $\mathrm{Tr}(f(x))$ is weakly
dominated by $S$; a label not pruned is a solution if $s(x)$ is a goal, and is
expanded otherwise. This is Algorithm 1 of Ren et al.
[[182]](../../README.md#ref-182), of which BOA*
[[181]](../../README.md#ref-181), EMOA\* [[182]](../../README.md#ref-182) and
LTMOA\* [[183]](../../README.md#ref-183) differ only in how the frontier
answers its two queries (whether a vector is dominated, and the addition of a
vector, which removes those it dominates). In these three, the frontier is the
minimum of the second components for $K = 2$, a balanced tree of the truncated
vectors, and an array of them. Our proofs use only what the frontiers have in
common, i.e., a set of truncated vectors that weakly dominates each vector
added to it, and hold for each of them. The lazy searches consult neither the
frontier of the solutions nor that of the state at the generation, and the
proofs hold with and without these checks, which are then the same as those at
the extraction, only earlier.

*Lemma 1.* The labels are extracted by nondecreasing $f_1$, and by a
nondecreasing $f$ in the lexicographic order.

*Proof.* By consistency $f(y) \succeq f(x)$ for a successor $y$ of $x$, and
hence $f(y) \geq_{lex} f(x)$; a label pushed after the extraction of $x$
descends from labels extracted no earlier than $x$; thus, its $f$ is not
smaller than $f(x)$, and OPEN gives the smallest of its labels. $\square$

*Lemma 2 (dimensionality reduction).* If the label $y$ was expanded at $s$
before the label $x$ of $s$ is checked, then
$\mathrm{Tr}(g(y)) \preceq \mathrm{Tr}(g(x))$ implies $g(y) \preceq g(x)$; if
the solution $z$ was found before $x$ is checked, then
$\mathrm{Tr}(z) \preceq \mathrm{Tr}(f(x))$ implies $z \preceq f(x)$.

*Proof.* If $x$ is checked when extracted, $f_1(y) \leq f_1(x)$ by Lemma 1, and
$g_1(y) \leq g_1(x)$ since $h(s)$ is the same. If $x$ is checked when
generated, its parent was extracted after $y$, and $f_1(y)$ is not larger than
the $f_1$ of the parent, which is not larger than $f_1(x)$. The same holds for
$z$, whose $f$ is its cost. $\square$

A label $x$ is *useful* for $c \in \mathcal{F}$ if a path $Q$ from $s(x)$ to a
goal has $g(x) + c(Q) = c$; then $f(x) \preceq c$, and $h(s(x))$ is finite. The
root is useful for each $c$, and we write $|Q|$ for the number of edges of the
shortest such $Q$.

*Theorem 1.* `MOAStar` terminates and returns a cost-unique Pareto-optimal set,
whatever its frontier and with or without the lazy checks; when it stops at its
node limit, each path it returns is Pareto-optimal.

*Proof.* Termination: a label expanded at $s$ is not weakly dominated by the
labels expanded there before (by Lemma 2 for its truncation). If its path
visits $s$ twice, the prefix that ends at the first visit was expanded at $s$
before it, with a $g$ not larger. Hence each label expanded is a simple path,
of which there are finitely many, and each expansion pushes finitely many
labels.

Completeness: let $c \in \mathcal{F}$ never be found, and let $x$ be a label
useful for $c$, among all the labels ever generated, of least $|Q|$ (the labels
are finitely many). If $|Q| = 0$, $s(x)$ is a goal and $g(x) = c$; then either
$x$ was pruned by $S$, i.e., a solution $z$ with $z \preceq f(x) = c$ (Lemma 2)
was found, and $z = c$ since $c$ is Pareto-optimal, or $x$ became a solution;
both contradict the assumption. If $|Q| > 0$, $x$ was (i) expanded, and its
successor on $Q$, useful with a shorter path, was generated, against the choice
of $x$; (ii) pruned by $S$, which as above means that $c$ was found; or (iii)
pruned by $G(s(x))$. In the last case a label $y$ expanded at $s(x)$ has
$g(y) \preceq g(x)$ by Lemma 2, and hence $g(y) + c(Q) \preceq c$, i.e., $= c$;
thus, $y$ is useful with the same $Q$, and its successor on $Q$, useful with a
shorter path, was generated when $y$ was expanded, again against the choice of
$x$. OPEN is empty at the end, and no other case is left.

Soundness: let the label $z$ of a goal be extracted and not pruned. If a
solution of cost $g(z)$ was found before, its truncation, or a vector that
weakly dominates it, is in $S$, and $z$ is pruned; hence, no cost is found
twice. Suppose that $c \prec g(z)$ for some $c \in \mathcal{F}$. A label $w$
useful for $c$ has $f(w) \preceq c \prec g(z) = f(z)$, hence
$f(w) \lt_{lex} f(z)$, and none is in OPEN when $z$ is extracted. Then the
argument of the completeness, applied to the labels generated so far, shows
that $c$ has been found, and $\mathrm{Tr}(c) \preceq \mathrm{Tr}(f(z))$ makes
$S$ prune $z$, a contradiction. Hence each solution is Pareto-optimal; the
argument does not need the search to have ended, and holds at the node limit.
$\square$

The traces of the examples of the papers are those of the searches, which the
case `moa` checks. On the graph of figure 1 of
[[181]](../../README.md#ref-181), BOA\* takes 10 labels from OPEN and prunes 2
of its 11 successors when it generates them (tables 1 and 2 there); on that of
figure 1 of [[183]](../../README.md#ref-183), LTMOA\* takes 7 labels and
LazyLTMOA\* 9 (tables 1 and 2 there).

For integer costs, OPEN may also be `MOBucketQueue`, a bucket for each value of
$f_1$ on a ring that doubles as needed, as in [[184]](../../README.md#ref-184),
each bucket a heap in the order of OPEN. Its labels leave it in the same order
as from the heap, and the search, its expansions and its checks are those of
the heap. The buckets of [[184]](../../README.md#ref-184), instead, keep the
ties of $f_1$ in the order of their pushes, which may expand a label that a
later one of the same $f_1$ dominates at its state. There the solutions are
kept correct by removing, when a solution arrives, those of the same $f_1$ that
it dominates, and by a final filter of the solutions in the parallel search.

The frontiers trade the checks for their structure: with $K = 2$ the frontier
of a state is a number, and both queries take one comparison. In LTMOA\*, the
array scans its vectors, and removes a dominated one by moving the last one
into its place. In EMOA\*, the tree looks only at the keys
$\leq_{lex} \mathrm{Tr}(v)$, since a dominating vector is lexicographically
smaller. For $K = 3$ the keys have two components and form a staircase (the
second decreasing as the first grows); the tree then answers with the largest
key $\leq_{lex} \mathrm{Tr}(v)$ only, which has the least second component
among them, in $O(\log n)$ (TOA\* of [[182]](../../README.md#ref-182)), and the
keys that $\mathrm{Tr}(v)$ dominates follow it in a run. In all the searches,
the counter `checks` counts the comparisons of two vectors (or of two
components, for the frontier of two objectives), the measure of
[[183]](../../README.md#ref-183).

## NAMOA\* with dimensionality reduction

`NAMOAStarDR` keeps, besides the frontier $G_{cl}(s)$ of the truncated $g$ of
the labels expanded at $s$, the set $G_{op}(s)$ of the labels of $s$ in OPEN. A
new label $y$ of $t$ is pruned if $\mathrm{Tr}(g(y))$ is weakly dominated by
$G_{cl}(t)$, if $g(y)$ is weakly dominated by a label of $G_{op}(t)$, or if
$\mathrm{Tr}(f(y))$ is weakly dominated by $S$; otherwise it removes from
$G_{op}(t)$ and from OPEN the labels it dominates
[[185]](../../README.md#ref-185) [[180]](../../README.md#ref-180). A label
extracted is then not weakly dominated by a label expanded at its state, since
a label $z$ expanded there after $y$ was generated was in $G_{op}$ with it, and
the later of the two was checked against the earlier one. Hence only the check
of $S$ is done at the extraction, and the termination and the soundness are
those of Theorem 1. For the completeness, let $x$ be, among the useful labels
of least $|Q|$, the one generated last. If it was removed from OPEN by a later
label, or pruned by a label of $G_{op}$ that is removed later by another one,
that label is useful with the same $|Q|$ and was generated later. If it was
pruned by a label of $G_{op}$ that is expanded, or by $G_{cl}$, the successor
of that label on $Q$ is useful with a shorter path. The other cases are those
of Theorem 1.

## Searches on several orders of the objectives

`MOAStarPar` runs MOA\* on each of its workers with a different order of the
objectives, the cyclic orders first, as in [[184]](../../README.md#ref-184),
then the other permutations. A worker of order $(i_1, \ldots, i_K)$ has the
primary objective $i_1$, and it extracts its labels by nondecreasing $f_{i_1}$
(Lemma 1 in its order). Each worker publishes, for its primary objective $i$,
the largest $F_i$ of the $f_i$ it has extracted (with several workers on one
primary objective, the largest of theirs). Its solutions are kept in the sets
$V_1, \ldots, V_{K-1}$ of the header, $V_\lambda$ comparing only the components
$\lambda + 1, \ldots, K$ of its order. A vector $v$ of $V_\lambda$ whose
component $\lambda + 1$, of objective $m$, satisfies $v_m \leq F_m$ is moved to
$V_{\lambda+1}$, and the vector of $V_{K-1}$, whose last component is the only
one compared, lowers the upper bound $\bar f$ of its objective, which all the
workers read. A label $y$ is pruned by the solutions if $f_m(y) \geq \bar f_m$
for some $m$, or if the compared components of $f(y)$ are weakly dominated by a
vector of some $V_\lambda$. A worker of primary objective $i$ stops when it
extracts a label of $f_i \geq \bar f_i$, and the search ends as soon as one
worker ends. The values of $F$ and $\bar f$ are written with release and read
with acquire; hence, what a worker did before writing a value happens before
what another one does after reading it. Our proofs reason on this order, which
is a strict partial order on the finitely many events of a run.

We call *justified* a truncated vector $v$, from a solution $z$ of a worker, if
for each objective $m$ whose component it no longer compares, $z_m \leq F_m$
held for a value of $F_m$ read by the worker, i.e., some worker of primary
objective $m$ had extracted a label $x'$ with $f_m(x') \geq z_m$ before. The
bound $\bar f_m = z_m$ is the vector truncated in all the components but $m$.

*Lemma 3.* Let a worker $j$ of primary objective $i$ extract a label $x$. Then
each $c \in \mathcal{F}$ with $c_i < f_i(x)$ is found by some worker, by the
end of the search.

*Proof.* By induction on the order of the events. Let $c$ be one such vector,
not found, and choose the useful label $w$ for $c$ of least $|Q|$ among those
that $j$ generated before extracting $x$. Since $f_i(w) \leq c_i < f_i(x)$, $w$
is not in OPEN when $x$ is extracted, and the cases of Theorem 1 apply, with
one more: $w$ (or a goal label $w$ of cost $c$) was pruned by a justified
vector $v$ of a solution $z$, i.e., the compared components of $z$ are not
larger than those of $f(w) \preceq c$. If $z_m \leq c_m$ also for each
component $m$ that $v$ no longer compares, then $z \preceq c$ and $c = z$ is
found. Otherwise $c_m < z_m \leq f_m(x')$ for some $m$ and some label $x'$
extracted by a worker of primary objective $m$ before the pruning, and hence
before the extraction of $x$; by the induction, $c$ is found. $\square$

*Theorem 2.* `MOAStarPar` returns a cost-unique Pareto-optimal set: each vector
of $\mathcal{F}$ is found by at least one worker, each solution found is
Pareto-optimal, and the duplicates of the workers are merged. With a node limit
each path it returns is Pareto-optimal.

*Proof.* Each worker terminates, as MOA\* does; let the worker $j$ end first.
If its OPEN is empty, the argument of Lemma 3 applies to all the labels it
generated, and each $c \in \mathcal{F}$ is found. If it stopped at a label of
$f_i \geq \bar f_i$, each label left in its OPEN has $f_i \geq \bar f_i$ and is
pruned by that bound, which is justified, and the same argument applies. For
the soundness, let a worker add the solution $z$ with $c \prec z$,
$c \in \mathcal{F}$. As in Theorem 1, no label useful for $c$ is in its OPEN;
hence, a useful label of the worker was pruned by its frontier of a state (and
then a shorter useful label exists) or by its solutions. Each pruning rule is
monotone, i.e., a vector weakly dominated by a pruned one is pruned too, at any
later time, since the sets $V_\lambda$ only grow in power (a vector moves to a
set that compares fewer components, or is replaced by one that dominates it)
and the bounds only decrease. Thus, $z \succeq c \succeq f(w)$ would have been
pruned too. Within a worker no cost is found twice, by its $V_1$, and the merge
removes the copies of the costs found by several workers. $\square$

Two choices are ours: first, the sets $V_\lambda$ are updated, besides at each
solution as in [[184]](../../README.md#ref-184), every `flush_period`
extractions when $F$ has grown, which justifies the truncations earlier. Also,
the check of a label against the solutions starts from $V_{K-1}$, the smallest
and strongest set, which saves comparisons and changes no outcome.

## Hash-distributed multi-objective search

`HDMOAStar` assigns each state to the worker that owns it by a hash of the
state, as HDA\* [[2]](../../README.md#ref-2). Thus, the paths generated by a
worker travel to the owner of their last state, which checks them and keeps
them in its OPEN, ordered lexicographically by $f$; a solution is kept by the
owner of the goal and sent to all the workers, each of which prunes with its
copy of the known solutions. Since the workers do not extract the labels in a
global order, neither Lemma 1 nor Lemma 2 holds across the workers: a label
expanded may be dominated later by a label that arrives at the same state, and
a solution may be dominated by a later one. Hence the frontier of a state keeps
the full $g$ of the labels expanded there, and a label is pruned if its $g$ is
weakly dominated by that frontier or its $f$ by a known solution; the owner of
a goal removes from its solutions those that a new one dominates. This is a
label-correcting search, with the termination detection of the mailbox of
[ADS](https://github.com/dmeoli/ADS) (see [Parallel A*](parallel-astar.md)).

*Theorem 3.* `HDMOAStar` terminates and returns a cost-unique Pareto-optimal
set.

*Proof.* Termination: a label is expanded at $s$ only if no label expanded
there before weakly dominates it; hence, as in Theorem 1, each label expanded
is a simple path. Also, the labels and the messages are finitely many, and the
mailbox ends the search when no worker has a label and no message travels.
Completeness: let $c \in \mathcal{F}$, and choose the useful label $x$ of least
$|Q|$ among all the labels created. If it was expanded, its successor on $Q$
was sent to its owner, which (i) made a useful label of shorter path, (ii)
pruned it by the frontier of the state, whose label $y$ with $g(y) \preceq$
that of the successor is useful and expanded, and then the successor of $y$ on
$Q$ is useful with a path shorter still, or (iii) pruned it by a solution
$z \preceq c$, i.e., $z = c$. The same holds if $x$ itself was pruned, and if
$s(x)$ is a goal with $g(x) = c$, $x$ became a solution unless a solution of
cost $c$ was there already. A solution is removed only by one that dominates
it, which for $c$ cannot exist. Hence the owner of some goal ends with a
solution of cost $c$ for each $c \in \mathcal{F}$, and the solutions of the
owners, filtered to the non-dominated ones with each cost once, are
$\mathcal{F}$. $\square$

The proof does not depend on the order in which the workers extract their
labels, which only decides how many of them are expanded and later dominated.
Hence a worker expands a label only if its $f_1$ is at most $L + w$, where $w$
is the option `window` (0 by default, negative for none) and $L$ is the least
$f_1$ of the tops of the OPEN lists and of the labels still traveling, which
each worker publishes. To this aim a worker counts the labels it sends, each
owner the labels it receives from each worker, and a label counts in the $f_1$
of its sender until all those it sent have arrived. With the window, the
workers advance together by $f_1$, close to the order of the sequential search,
and the worker whose top is the least is never held back. Hence, the search
does not stall, since the least label is in an OPEN list, or traveling to an
owner that polls its inbox also while it waits. A value of $L$ read late
changes the order but not the front. The frontiers of the states and the copies
of the solutions rest on one remark.

*Lemma 4.* Let a set $Z$ of vectors be kept truncated in a frontier of EMOA\*,
a truncation being added only if no key dominates it, and let
$T = \max_{z \in Z} z_1$. For a vector $x$ with $x_1 \geq T$, a key weakly
dominates $\mathrm{Tr}(x)$ if and only if some $z \in Z$ weakly dominates $x$.

*Proof.* A key $\mathrm{Tr}(z) \preceq \mathrm{Tr}(x)$ comes from some
$z \in Z$ with $z_1 \leq T \leq x_1$, hence $z \preceq x$. Conversely, if
$z \preceq x$, the truncation of $z$ was added, or a key that dominates it was
there, and a key is removed only by one that dominates it; hence, a key
dominates $\mathrm{Tr}(z) \preceq \mathrm{Tr}(x)$. $\square$

Hence each set keeps its vectors twice. Truncated, as in Lemma 4, it answers in
one check of the sequential search (one comparison for two and three
objectives) the vectors with $x_1 \geq T$, i.e., those of the labels that
arrive at a state in the order of $f_1$, which the window makes the common
case. In full, it answers the others: for two objectives they are kept in the
tree of EMOA\* of the vectors $(0, v)$, whose truncations are the vectors $v$
and form a staircase, in $O(\log n)$. For three they are kept in a
`LayerIndex`, i.e., the vectors by increasing first component in layers of a
few of them and a Fenwick tree over the layers
[[189]](../../README.md#ref-189), each node of which keeps the vectors of the
layers it covers truncated as in Lemma 4. A check of $x$ queries the
$O(\log L)$ nodes that cover the layers before that of $x$, whose vectors have
a first component at most $x_1$ (and therefore Lemma 4 makes each answer
exact), and then compares $x$ with the vectors of its layer whose first
component is not larger. These are the layers of $f_1$ of the indexed pruning
of Xu et al. [[187]](../../README.md#ref-187), here on prefixes. For more
objectives a node of the Fenwick tree costs too much to keep, and a
`BlockIndex` keeps the vectors in blocks by increasing first component, with a
segment tree of the minima of the blocks, component by component, into whose
nodes a check descends only if the minimum weakly dominates $x$ (the
block-minimum index of [[187]](../../README.md#ref-187)); an insertion compares
no vector. Both keep all the vectors added, also those that a later one
dominates, and the question of the search is whether some vector added weakly
dominates $x$, which keeping them does not change.

With a node limit the solutions found are mutually non-dominated, but they need
not be Pareto-optimal, as in any label-correcting search. This design is close
to that of SIP-MOSP of Xu et al. [[187]](../../README.md#ref-187), whose update
workers own the frontiers of the vertices assigned to them; here the owner of a
state both keeps its frontier and expands its labels, as in HDA\*, and hence no
frontier is written by two workers and no worker waits for another. In turn,
the parallel label-setting search of Sanders and Mandow
[[186]](../../README.md#ref-186) expands at once all the labels whose $f$ is
not dominated in OPEN, which keeps the label-setting property at the price of a
global synchronization per round; the searches here have no such
synchronization.
