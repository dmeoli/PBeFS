# Multi-objective search

We are concerned here with the searches of a `MultiObjectiveDomain`, whose
edges cost a vector of $K \geq 2$ nonnegative components and whose answer is
not one path but the Pareto front of the paths from the start $s_0$ to a goal.
This page fixes the notation, proves that the sequential framework of
[MOAStar.h](../../include/PBeFS/astar/MOAStar.h) returns the front with any of
its frontiers, lazy or not, and with NAMOA\*dr, and then proves the same for
the two parallel searches of
[ParallelMOAStar.h](../../include/PBeFS/astar/ParallelMOAStar.h): the searches
on several orders of the objectives of Ahmadi et al.
[[130]](../../README.md#ref-130), whose argument we write in full, and the
search distributed by the hash of the states, a label-correcting search of our
own. The pseudocode of each search is in its header, and the notation of the
single-objective searches is that of [Best-first search and A*](best-first-search.md).

## Notation

The cost of a path $\pi$ is the vector $c(\pi) \in \mathbb{R}^K_{\geq 0}$, the
sum of the vectors of its edges. We write $u \preceq v$ ($u$ weakly dominates
$v$) if $u_i \leq v_i$ for every $i$, $u \prec v$ if moreover $u \neq v$, and
$u <_{lex} v$ for the lexicographic order; clearly $u \preceq v$ implies
$u \leq_{lex} v$, and $u \prec v$ implies $u <_{lex} v$. A path to a goal is
Pareto-optimal if no path to a goal dominates it, and the Pareto front
$\mathcal{F}$ is the set of the costs of the Pareto-optimal paths; the searches
return a cost-unique Pareto-optimal set, i.e., one path for each vector of
$\mathcal{F}$ [[126]](../../README.md#ref-126). Since the costs are
nonnegative, every path to a goal is weakly dominated by a simple one (drop its
cycles), and hence by a vector of $\mathcal{F}$, which is finite. The
truncation $\mathrm{Tr}(v) = (v_2, \ldots, v_K)$ drops the first component and
$\mathrm{Tr}^\lambda$ the first $\lambda$ ones.

A label $x$ is a path from $s_0$ to its state $s(x)$, of cost $g(x)$, with
$f(x) = g(x) + h(s(x))$; the heuristic is consistent, i.e., $h(t) = 0$ at the
goals and $h(s) \preceq c(s, t) + h(t)$ for every edge, so that each component
is a consistent heuristic of its own objective, and hence admissible: if $Q$ is
a path from $s(x)$ to a goal, $f(x) \preceq g(x) + c(Q)$. A component equal to
infinity marks a state from which no goal is reachable; such states are never
generated, which loses nothing. On the road graphs the heuristic is the ideal
point, $h_i(u) = d_i(u, t)$, the distance of $u$ from the target by the $i$-th
objective alone, computed by $K$ backward runs of Dijkstra's algorithm, as in
[[127]](../../README.md#ref-127); it is consistent by the triangle inequality
of each $d_i$. On the grids of
[MultiObjectiveGrid.h](../../include/PBeFS/astar/MultiObjectiveGrid.h) the
default is $h_i(u) = \ell\,\delta(u, t)$, $\ell$ being the least cost of a move
and $\delta$ the Manhattan (4 moves) or the Chebyshev (8 moves) distance, which
a move lowers by 1 at most, while it costs at least $\ell$ in each component.

## The sequential framework

`MOAStar` extracts the labels by the lexicographic order of $f$ (the ties to
the label pushed first), keeps for each state $s$ the frontier $G(s)$ of the
truncated $g$ of the labels expanded there, and the frontier $S$ of the
truncated $f$ of the solutions; a label is pruned when extracted (and, unless
the search is lazy, when generated) if $\mathrm{Tr}(g(x))$ is weakly dominated
by $G(s(x))$, $s(x)$ not being a goal, or if $\mathrm{Tr}(f(x))$ is weakly
dominated by $S$; a label not pruned is a solution if $s(x)$ is a goal, and is
expanded otherwise. This is Algorithm 1 of Ren et al.
[[128]](../../README.md#ref-128), of which BOA* [[127]](../../README.md#ref-127),
EMOA\* [[128]](../../README.md#ref-128) and LTMOA\* [[129]](../../README.md#ref-129)
differ only in how the frontier answers its two questions (is a vector
dominated, and add a vector removing those it dominates): the minimum of the
second components for $K = 2$, a balanced tree of the truncated vectors, and an
array of them. The proofs below use only what the frontiers have in common, a
set of truncated vectors that weakly dominates every vector added to it, and
hold for each of them. The frontier of the solutions is not consulted at the
generation by the lazy searches, nor that of the state, and the proofs hold
with and without these checks, which are then the same as those at the
extraction, only earlier.

*Lemma 1.* The labels are extracted by nondecreasing $f_1$, and by a
nondecreasing $f$ in the lexicographic order.

*Proof.* By consistency $f(y) \succeq f(x)$ for a successor $y$ of $x$, and
hence $f(y) \geq_{lex} f(x)$; a label pushed after the extraction of $x$
descends from labels extracted no earlier than $x$, so its $f$ is not smaller
than $f(x)$, and OPEN gives the smallest of its labels. $\square$

*Lemma 2 (dimensionality reduction).* If the label $y$ was expanded at $s$
before the label $x$ of $s$ is checked, then $\mathrm{Tr}(g(y)) \preceq
\mathrm{Tr}(g(x))$ implies $g(y) \preceq g(x)$; if the solution $z$ was found
before $x$ is checked, then $\mathrm{Tr}(z) \preceq \mathrm{Tr}(f(x))$ implies
$z \preceq f(x)$.

*Proof.* If $x$ is checked when extracted, $f_1(y) \leq f_1(x)$ by Lemma 1, and
$g_1(y) \leq g_1(x)$ since $h(s)$ is the same; if $x$ is checked when
generated, its parent was extracted after $y$, and $f_1(y)$ is not larger than
the $f_1$ of the parent, which is not larger than $f_1(x)$. The same holds for
$z$, whose $f$ is its cost. $\square$

A label $x$ is *useful* for $c \in \mathcal{F}$ if a path $Q$ from $s(x)$ to a
goal has $g(x) + c(Q) = c$; then $f(x) \preceq c$, and $h(s(x))$ is finite.
The root is useful for every $c$. We write $|Q|$ for the number of edges of the
shortest such $Q$.

*Theorem 1.* `MOAStar` terminates and returns a cost-unique Pareto-optimal set,
whatever its frontier and with or without the lazy checks; when it stops at its
node limit, every path it returns is Pareto-optimal.

*Proof.* Termination: a label expanded at $s$ is not weakly dominated by the
labels expanded there before (by Lemma 2 for its truncation). If its path
visits $s$ twice, the prefix that ends at the first visit was expanded at $s$
before it, with a $g$ not larger: hence every label expanded is a simple path,
of which there are finitely many, and each expansion pushes finitely many
labels.

Completeness: let $c \in \mathcal{F}$ never be found, and let $x$ be a label
useful for $c$, among all the labels ever generated, of least $|Q|$ (the
labels are finitely many). If $|Q| = 0$, $s(x)$ is a goal and $g(x) = c$;
then $x$ was pruned by $S$, i.e., a solution $z$ with $z \preceq f(x) = c$
(Lemma 2) was found, and $z = c$ since $c$ is Pareto-optimal, or $x$ became a
solution, both against the assumption. If $|Q| > 0$, $x$ was (i) expanded, and
its successor on $Q$, useful with a shorter path, was generated, against the
choice of $x$; (ii) pruned by $S$, which as above means that $c$ was found; or
(iii) pruned by $G(s(x))$: a label $y$ expanded at $s(x)$ has $g(y) \preceq
g(x)$ by Lemma 2, so $g(y) + c(Q) \preceq c$, hence $= c$, $y$ is useful with
the same $Q$, and its successor on $Q$, useful with a shorter path, was
generated when $y$ was expanded, again against the choice of $x$. OPEN is empty
at the end, and no other case is left.

Soundness: let the label $z$ of a goal be extracted and not pruned. If a
solution of cost $g(z)$ was found before, its truncation, or a vector that
weakly dominates it, is in $S$, and $z$ is pruned; so no cost is found twice.
Suppose that $c \prec g(z)$ for some $c \in \mathcal{F}$. A label $w$ useful
for $c$ has $f(w) \preceq c \prec g(z) = f(z)$, hence $f(w) <_{lex} f(z)$, and
none is in OPEN when $z$ is extracted. Then the argument of the completeness,
applied to the labels generated so far, shows that $c$ has been found, and
$\mathrm{Tr}(c) \preceq \mathrm{Tr}(f(z))$ makes $S$ prune $z$: a
contradiction. Hence every solution is Pareto-optimal; the argument does not
need the search to have ended, and holds at the node limit.
$\square$

The traces of the examples of the papers are those of the searches: on the
graph of figure 1 of [[127]](../../README.md#ref-127) BOA\* takes 10 labels
from OPEN and prunes 2 of its 11 successors when it generates them (tables 1
and 2 there), and on that of figure 1 of
[[129]](../../README.md#ref-129) LTMOA\* takes 7 labels and LazyLTMOA\* 9
(tables 1 and 2 there), which the case `moa` checks.

The frontiers trade the checks for their structure. With $K = 2$ the
frontier of a state is a number, and both questions take one comparison. The
array of LTMOA\* scans its vectors, and removes a dominated one by moving the
last one into its place. The tree of EMOA\* looks only at the keys
$\leq_{lex} \mathrm{Tr}(v)$, since a dominating vector is lexicographically
smaller, and for $K = 3$, when the keys have two components and form a
staircase (the second decreasing as the first grows), it answers with the
largest key $\leq_{lex} \mathrm{Tr}(v)$ alone, which has the least second
component among them, in $O(\log n)$ (TOA\* of [[128]](../../README.md#ref-128));
the keys that $\mathrm{Tr}(v)$ dominates follow it in a run. The counter
`checks` of the searches counts the comparisons of two vectors (or of two
components, for the frontier of two objectives), the measure of
[[129]](../../README.md#ref-129).

## NAMOA\* with dimensionality reduction

`NAMOAStarDR` keeps, besides the frontier $G_{cl}(s)$ of the truncated $g$ of
the labels expanded at $s$, the set $G_{op}(s)$ of the labels of $s$ in OPEN; a
new label $y$ of $t$ is pruned if $\mathrm{Tr}(g(y))$ is weakly dominated by
$G_{cl}(t)$, if $g(y)$ is weakly dominated by a label of $G_{op}(t)$, or if
$\mathrm{Tr}(f(y))$ is weakly dominated by $S$, and otherwise it removes from
$G_{op}(t)$ and from OPEN the labels it dominates
[[131]](../../README.md#ref-131) [[126]](../../README.md#ref-126). A label
extracted is then not weakly dominated by a label expanded at its state: a
label $z$ expanded there after $y$ was generated was in $G_{op}$ with it, and
the later of the two was checked against the earlier one. Hence only the check
of $S$ is done at the extraction, and the termination and the soundness are
those of Theorem 1. For the completeness, let $x$ be, among the useful labels
of least $|Q|$, the one generated last. If it was removed from OPEN by a later
label, or pruned by a label of $G_{op}$ that is removed later by another one,
that label is useful with the same $|Q|$ and was generated later; if it was
pruned by a label of $G_{op}$ that is expanded, or by $G_{cl}$, the successor of
that label on $Q$ is useful with a shorter path; the other cases are those of
Theorem 1.

## Searches on several orders of the objectives

`MOAStarPar` runs MOA\* on each of its workers with its own order of the
objectives, the cyclic orders first, as in [[130]](../../README.md#ref-130),
then the other permutations; the worker of order $(i_1, \ldots, i_K)$ has the
primary objective $i_1$, and it extracts its labels by nondecreasing
$f_{i_1}$ (Lemma 1 in its order). Each worker publishes, for its primary
objective $i$, the largest $F_i$ of the $f_i$ it has extracted (with several
workers on one primary objective, the largest of theirs). Its solutions are
kept in the sets $V_1, \ldots, V_{K-1}$ of the header, $V_\lambda$ comparing
only the components $\lambda + 1, \ldots, K$ of its order; a vector $v$ of
$V_\lambda$ whose component $\lambda + 1$, of objective $m$, satisfies $v_m
\leq F_m$ is moved to $V_{\lambda+1}$, and the vector of $V_{K-1}$, whose last
component is the only one compared, lowers the upper bound $\bar f$ of its
objective, which every worker reads. A label $y$ is pruned by the solutions if
$f_m(y) \geq \bar f_m$ for some $m$, or if the compared components of $f(y)$
are weakly dominated by a vector of some $V_\lambda$; the worker of primary
objective $i$ stops when it extracts a label of $f_i \geq \bar f_i$, and the
search ends as soon as one worker ends. The values of $F$ and $\bar f$ are
written with release and read with acquire, so that what a worker did before
writing a value happens before what another one does after reading it; the
proofs reason on this order, which is a strict partial order on the finitely
many events of a run.

We call a truncated vector $v$, from a solution $z$ of a worker, *justified* if
for every objective $m$ whose component it no longer compares, $z_m \leq F_m$
held for a value of $F_m$ read by the worker, i.e., some worker of primary
objective $m$ had extracted a label $x'$ with $f_m(x') \geq z_m$ before; the
bound $\bar f_m = z_m$ is the vector truncated in every component but $m$.

*Lemma 3.* Let a worker $j$ of primary objective $i$ extract a label $x$.
Then every $c \in \mathcal{F}$ with $c_i < f_i(x)$ is found by some worker,
by the end of the search.

*Proof.* By induction on the order of the events. Let $c$ be one such vector,
not found, and choose the useful label $w$ for $c$ of least $|Q|$ among those
that $j$ generated before extracting $x$. Since $f_i(w) \leq c_i < f_i(x)$,
$w$ is not in OPEN when $x$ is extracted, and the cases of Theorem 1 apply,
with one more: $w$ (or a goal label $w$ of cost $c$) was pruned by a justified
vector $v$ of a solution $z$, i.e., the compared components of $z$ are not
larger than those of $f(w) \preceq c$. If $z_m \leq c_m$ also for every
component $m$ that $v$ no longer compares, then $z \preceq c$ and $c = z$ is
found. Otherwise $c_m < z_m \leq f_m(x')$ for some $m$ and some label $x'$
extracted by a worker of primary objective $m$ before the pruning, and hence
before the extraction of $x$; by the induction $c$ is found. $\square$

*Theorem 2.* `MOAStarPar` returns a cost-unique Pareto-optimal set: every
vector of $\mathcal{F}$ is found by at least one worker, every solution found
is Pareto-optimal, and the duplicates of the workers are merged. With a node
limit every path it returns is Pareto-optimal.

*Proof.* Every worker terminates, as MOA\* does. Let the worker $j$ end first.
If its OPEN is empty, the argument of Lemma 3 applies to all the labels it
generated, and every $c \in \mathcal{F}$ is found; if it stopped at a label of
$f_i \geq \bar f_i$, every label left in its OPEN has $f_i \geq \bar f_i$ and
is pruned by that bound, which is justified, and the same argument applies.
For the soundness, let a worker add the solution $z$ with $c \prec z$, $c \in
\mathcal{F}$: as in Theorem 1, no label useful for $c$ is in its OPEN, so that
a useful label of the worker was pruned by its frontier of a state (and then a
shorter useful label exists) or by its solutions. Every pruning rule is
monotone: a vector weakly dominated by a pruned one is pruned too, at any later
time, since the sets $V_\lambda$ only grow in power (a vector moves to a set
that compares fewer components, or is replaced by one that dominates it) and
the bounds only decrease; hence $z \succeq c \succeq f(w)$ would have been
pruned too. Within a worker no cost is found twice, by its own $V_1$; the
merge removes the copies of the costs found by several workers. $\square$

Two choices are ours. The sets $V_\lambda$ are updated, besides at every
solution as in [[130]](../../README.md#ref-130), every `flush_period`
extractions when $F$ has grown, which justifies the truncations earlier; and
the check of a label against the solutions starts from $V_{K-1}$, the smallest
and strongest set, which saves comparisons and changes nothing else.

## Hash-distributed multi-objective search

`HDMOAStar` assigns each state to the worker that owns it by a hash of the
state, as HDA\* [[2]](../../README.md#ref-2): the paths generated by a worker
travel to the owner of their last state, which checks them and keeps them in
its OPEN, ordered lexicographically by $f$; a solution is kept by the owner of
the goal and sent to every worker, which prunes with its copy of the known
solutions. Since the workers do not extract the labels in a global order,
neither Lemma 1 nor Lemma 2 holds across the workers: a label expanded may be
dominated later by a label that arrives at the same state, and a solution may be
dominated by a later one. Hence the frontier of a state keeps the full $g$ of
the labels expanded there, a label is pruned if its $g$ is weakly dominated by
that frontier or its $f$ by a known solution, and the owner of a goal removes
from its solutions those that a new one dominates. This is a label-correcting
search, with the termination detection of the mailbox of ADS
(see [Parallel A*](parallel-astar.md)).

*Theorem 3.* `HDMOAStar` terminates and returns a cost-unique Pareto-optimal
set.

*Proof.* Termination: a label is expanded at $s$ only if no label expanded
there before weakly dominates it, so that, as in Theorem 1, every label
expanded is a simple path; the labels and the messages are finitely many, and
the mailbox ends the search when no worker has a label and no message travels.
Completeness: let $c \in \mathcal{F}$ and choose the useful label $x$ of least
$|Q|$ among all the labels created. If it was expanded, its successor on $Q$ was
sent to its owner, which made a useful label of shorter path, or pruned it by
the frontier of the state, whose label $y$ with $g(y) \preceq$ that of the
successor is useful and expanded, so that its own successor on $Q$ is useful
with a path shorter still, or by a solution $z \preceq c$, i.e., $z = c$; the
same holds if $x$ itself was pruned, and if $s(x)$ is a goal with $g(x) = c$,
$x$ became a solution unless a solution of cost $c$ was there already. A
solution is removed only by one that dominates it, which for $c$ cannot exist.
Hence the owner of some goal ends with a solution of cost $c$ for every $c \in
\mathcal{F}$, and the solutions of the owners, filtered to the non-dominated
ones with each cost once, are $\mathcal{F}$. $\square$

The proof does not depend on the order in which the workers extract their
labels, which only decides how many of them are expanded and later dominated.
Hence a worker expands a label only if its $f_1$ is at most $L + w$, $L$ being
the least $f_1$ of the tops of the OPEN lists, which every worker publishes,
and $w$ the option `window` (0 by default, negative for none): the workers then
advance together by $f_1$, close to the order of the sequential search, and
the worker whose top is the least is never held back, so that the search does
not stall. A value of $L$ read late, or a label travelling to its owner with an
$f_1$ below $L$, changes the order and not the front. The frontiers of the
states and the copies of the solutions keep the full vectors in the tree of
EMOA\* (the vectors $(0, v)$, whose truncations are the vectors $v$), so that
for two objectives, where the vectors of a frontier form a staircase, a check
takes $O(\log n)$ comparisons.

With a node limit the solutions found are mutually non-dominated, but they need
not be Pareto-optimal, as in any label-correcting search. The design is close to
that of SIP-MOSP of Xu et al. [[133]](../../README.md#ref-133), whose update
workers own the frontiers of the vertices assigned to them; here the owner of a
state both keeps its frontier and expands its labels, as in HDA\*, so that no
frontier is written by two workers and no worker waits for another. The parallel
label-setting search of Sanders and Mandow [[132]](../../README.md#ref-132)
expands at once all the labels whose $f$ is not dominated in OPEN, which keeps
the label-setting property at the price of a global synchronization per round;
the searches here have none.
