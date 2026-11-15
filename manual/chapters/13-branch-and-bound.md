# 13. Branch-and-bound and its parallel schemes {#ch-13}

[Source: [`include/PBeFS/bnb/BranchAndBound.h`][BranchAndBound.h] (the concept,
the options, the solution, the queue and the sequential search),
[`include/PBeFS/bnb/ParallelBranchAndBound.h`][ParallelBranchAndBound.h] (the
shared incumbent, the shared and the distributed queues),
[`MultiQueueBranchAndBound.h`][MultiQueueBranchAndBound.h],
[`FarmBranchAndBound.h`][FarmBranchAndBound.h],
[`DeterministicBranchAndBound.h`][DeterministicBranchAndBound.h],
[`StealingBranchAndBound.h`][StealingBranchAndBound.h],
[`SelfSplitBranchAndBound.h`][SelfSplitBranchAndBound.h],
[`MPIBranchAndBound.h`][MPIBranchAndBound.h],
[`BranchAndBoundScheme.h`][BranchAndBoundScheme.h] and
[`Knapsack.h`][Knapsack.h] of the same directory;
[`tools/BranchAndBoundBench.cpp`][BranchAndBoundBench.cpp];
[`test/test.cpp`][test.cpp]]

We are concerned in this chapter with the minimization of a function
$v$ over a finite (but possibly "huge") set $X$, i.e., with finding

$$
z^* = \min \{ \, v(x) : x \in X \, \} ,
$$

when $X$ can be partitioned recursively into smaller sets, each of which
comes with a lower bound on the values that it contains; a maximization
is the minimization of $-v$, which is how the library treats it
internally. The method is branch-and-bound (B&B), which Land and Doig
[[146]](Z-bibliography.md#ref-146) introduced for discrete programming.
A$^*$ is a special case of it [[17]](Z-bibliography.md#ref-17): the
subsets are the nodes of a tree, the bound of a node plays the part of
$f = g + h$, and the best value found so far (the *incumbent*) prunes
the nodes that cannot improve on it. The layer of the library that
implements B&B shares with the graph searches of the previous chapters
the executors of the workers ([§1.6](01-introduction.md#sec-1-6)), their
placement on the cores and the detection of the termination. It does not
share the node table and the OPEN lists of A$^*$, which are specialized
for integer costs and for the duplicates of a graph.

The structure of the chapter is as follows. [§13.1](#sec-13-1) gives the
concept `BranchAndBound` and the types around it, and [§13.2](#sec-13-2)
the notions of node, bound, estimate and pruning, on the 0-1 knapsack
that serves as the running example of the chapter. [§13.3](#sec-13-3)
gives the sequential search, [§13.4](#sec-13-4) its correctness and the
bound it reports when it stops early, and [§13.5](#sec-13-5) the rules
by which it selects the next node, depth-first B&B among them, together
with the queue that implements them and the cost of the whole search.
[§13.6](#sec-13-6) shows A$^*$ as a special case. [§13.7](#sec-13-7)
states what a parallel B&B has to guarantee, with the shared incumbent,
the dive and the limits common to all the schemes. [§13.8](#sec-13-8)
and [§13.9](#sec-13-9) then give the eight parallel schemes of the
library, each with its pseudocode, the invariant that makes its result
optimal, the proof that its detection of the end is sound and that the
end comes, and its engineering: first those with one pool of open nodes
(the shared queue, the MultiQueue, the farm and the deterministic
rounds), then those with one pool per worker (the distributed queues,
their version over MPI, the work stealing and the SelfSplit).
[§13.10](#sec-13-10) compares the schemes and shows how a program
chooses one by its name. [§13.11](#sec-13-11) treats the reoptimization
from the leaves of a previous search, and [§13.12](#sec-13-12) the
interface that a B&B solver built on a modeling system would use; the
last sections give the benchmark, the notes and the exercises. The
notation is that of [Chapter 1](01-introduction.md#ch-1), and the
conventions of the pseudocode are those of
[§1.4](01-introduction.md#sec-1-4). In the pseudocode, as in the
headers, $\beta(n)$ is the bound of the node $n$, $z$ the incumbent
value, Q a queue of open nodes and $p$ the number of workers.

## 13.1 The problem and its interface {#sec-13-1}

A problem is given to the B&B layer as a type `P` that satisfies the
concept `BranchAndBound` ([`BranchAndBound.h:338`][BranchAndBound.h:338]),
i.e., that provides

* the copyable type `node_type` of its nodes;
* the constant `sense`, `Sense::minimize` or `Sense::maximize`;
* the method `root()`, which gives the node of the whole set $X$;
* optionally, the type `value_type` of its values, `double` if it gives
  none;
* either the methods `evaluate( n )` and `branch( n , sink )` of a
  `BBEvaluator` ([`BranchAndBound.h:315`][BranchAndBound.h:315]), which have to
  be safe to call concurrently in a parallel search, or the method
  `make_worker( i )`, which returns the `BBEvaluator` of the worker $i$.

In the header the concept reads as follows (the concepts `BBValue` and
`BBEvaluator` are described below):

```cpp
--8<-- "snippets/branch-and-bound-sig-concept.cpp"
```

The method `evaluate( n )` computes a bound of the node, possibly
storing in `n` what the branching needs (say, the fractional solution of
a relaxation), and returns a `BasicEvaluation< V >` made of the bound, a
flag for an infeasible node and an optional feasible value, of which
`n` is the solution. In turn, `branch( n , sink )` calls
`sink( c , e )` for each child `c` of the evaluated node `n`. Here `e`
is a bound of `c` that is "cheap" to compute, typically the bound of
`n` (or the infinity of the sense of `P`, which means the same thing);
the search uses it to order the child and to prune it before its
evaluation. An evaluator given by `make_worker( i )` exists for the problems
whose evaluation keeps a state (e.g., a model and a solver of the
relaxation that are modified at each node): each worker then has a copy
of that state, and `P` is only read. Such an evaluator may also
give `cutoff( z )`, which the searches call before each `evaluate` and
`branch` with the incumbent value $z$ (in the sense of `P`, infinite if
there is none). With it, e.g., the relaxation of a node can stop as
soon as it is proven not to improve $z$, the node being then returned
as infeasible. Note that a worker should take $z$ from here rather than
from a value it shares with the others: the deterministic search of
[§13.8](#sec-13-8) passes the incumbent of the start of its round, and
an evaluation against a value that the other workers change meanwhile
would make its outcome depend on their timing.

Three refinements of the concept are tested by the searches that use
them. A `HashedBranchAndBound` gives `hash( n )`, by which the
distributed searches send a node to its owner. An
`EstimatingBranchAndBound` gives `estimate( n )`, an estimate of the
best feasible value under `n`, by which the best-estimate selection
orders the nodes. Finally, the search over MPI needs an
`MPIBranchAndBound` ([`MPIBranchAndBound.h:107`][MPIBranchAndBound.h:107]),
i.e., trivially copyable values and, unless its nodes are trivially copyable
too, the methods `serialize( n , out )` and `deserialize( in )`, which write a
node to bytes and read it back.

The values (the bounds, the feasible values, the incumbent and the absolute
tolerance) are of the type `V = bb_value_t< P >`, which has to satisfy the
concept `BBValue` ([`BranchAndBound.h:262`][BranchAndBound.h:262]), i.e., to be
a copyable, totally ordered type, constructible from an `int`, with negation
and subtraction. In addition, the traits class `BBValueTraits< V >`
([`BranchAndBound.h:205`][BranchAndBound.h:205]) gives the infinity of the type
(that of a floating point type, the largest value of an integer one) and what
the value of a solution with no feasible value is (a quiet NaN, or the
infinity). It also gives the size of the tolerance of a gap, and decides
whether a bound is "close enough" to the least open bound to dive into it. Its
primary template serves the signed arithmetic types, and a caller with another
type, e.g., the rationals of an exact arithmetic, specializes it (the test
`bb-values` does so with a type `Rational`). An integer type pays off whenever
the values of the solutions are integer, since the bound of a node can then be
rounded, as [§13.2](#sec-13-2) explains.

Options are given as a `BasicBBOptions< V >`
([`BranchAndBound.h:386`][BranchAndBound.h:386]), whose fields and defaults are

```cpp
--8<-- "snippets/branch-and-bound-sig-options.cpp"
```

i.e.,

* `selection`, the order of extraction of the open nodes
  ([§13.5](#sec-13-5)), best-first by default;
* `abs_gap` and `rel_gap`, the absolute and relative tolerances
  $\tau_a$ and $\tau_r$ on the gap (both 0 by default);
* `plunge_gap`, the fraction $\gamma$ of the dive of the plunging
  selection (0.5 by default);
* `incumbent`, a known feasible value, if any;
* `shared_incumbent`, a function that the search reads at each
  extraction and prunes by when it is better than the incumbent of
  the search (e.g., the incumbent of other searches running at the
  same time);
* `node_limit` (the nodes evaluated, by all the workers together in a
  parallel search, [§13.7](#sec-13-7); 0 for none) and `time_limit`
  (seconds, 0 for none);
* `keep_frontier`, which requests the leaves of the tree explored
  ([§13.11](#sec-13-11));
* `on_incumbent`, a function called with each new incumbent value.

The parallel searches take a `BasicParallelBBOptions< V >`
([`ParallelBranchAndBound.h:136`][ParallelBranchAndBound.h:136]). It adds the
parameters of the parallel searches of the library, i.e., the number of workers
`threads`, the nodes `batch` of a message, the evaluations
`flush_period` between two flushes of the outboxes and the `seed` of
the random owners; by default the batches have 16 nodes and are flushed
at each evaluation, since an evaluation is usually "expensive". It also
adds the parameters `round`, `split`, `split_dive` and `split_shared`
of two of the schemes ([§13.8](#sec-13-8) and [§13.9](#sec-13-9)).

Each search returns a `BBSolution< Node , V >`
([`BranchAndBound.h:476`][BranchAndBound.h:476]). It holds the status
(`optimal` within the tolerance, `infeasible`, `node_limit` or `time_limit`),
the incumbent node, its value, a bound on the optimum and the gap between the
two, which is valid whatever the status ([Theorem 13.6](#thm-13-6)). It also
holds the counters of `BBStats`: the nodes evaluated, branched,
generated, found infeasible and pruned at each of the three tests of
[§13.3](#sec-13-3), the improvements of the incumbent, the nodes and
messages sent, the largest queue, the time and the evaluations of each
worker. With `keep_frontier`, finally, it holds the leaves of the tree
explored as `BBLeaf` records, each with its node, its bound, its depth
and the reason why it is a leaf (`open`, `pruned`, `feasible` or
`infeasible`). A maximization problem is solved by negating its values
on the way in and on the way out (`detail::bb_signed()`), so that the
code, the pseudocode and the proofs below only minimize.

## 13.2 Nodes, bounds and the incumbent {#sec-13-2}

A *node* $n$ stands for a subset $X(n) \subseteq X$, the root for the
whole $X$, and *branching* $n$ gives children $n_1, \ldots, n_k$ whose
subsets cover $X(n)$ (usually they partition it). The *evaluation* of
$n$ gives a bound $\beta'(n)$ with

$$
\beta'(n) \le \min \{ \, v(x) : x \in X(n) \, \} ,
$$

e.g., the optimal value of a relaxation of the problem restricted to
$X(n)$, together with two optional facts: that $X(n)$ is empty (the node
is *infeasible*), or a feasible value $v(x)$ for some $x \in X(n)$
(e.g., when the solution of the relaxation is feasible, or by a
heuristic). We call *incumbent* the best feasible value $z$ found so far
($+\infty$ if none is known and none is given at the start), and
*incumbent node*, if any, a node of value $z$.

Since a child is a subset of its parent, the bound of the parent is
also a bound of the child. For each open node the library keeps the
bound

$$
\beta(m) = \max \{ \, \beta(n), \; e(m) \, \}
$$

for the child $m$ of $n$, where $e(m)$ is the estimate given by
`branch()` (e.g., the bound obtained by fixing one variable in the
solution of the relaxation of the parent). After the evaluation of $m$
it replaces $\beta(m)$ by $\max \{ \beta(m), \beta'(m) \}$. Thus, the
evaluation, which is usually "expensive", is postponed until the node
is extracted, and the bound never decreases along a path of the tree;
the second property is the rule that the literature of A$^*$ calls
*pathmax* (cf. [§13.6](#sec-13-6)).

A node is *pruned*, i.e., discarded together with its subtree, when it
is infeasible or when its bound cannot improve the incumbent by more
than the tolerance, i.e., when

$$
\beta(n) \ge z - \max \{ \tau_a , \; \tau_r |z| \} ,
$$

with an absolute tolerance $\tau_a \ge 0$ and a relative one
$\tau_r \ge 0$ (`abs_gap` and `rel_gap`); when $z = +\infty$, the node
is pruned only if $\beta(n) = +\infty$. We write $\theta(z) = z - \max
\{ \tau_a, \tau_r |z| \}$ for the threshold, and PRUNABLE($\beta$) for
the test $\beta \ge \theta(z)$. With $\tau_a = \tau_r = 0$ the test is
$\beta(n) \ge z$, since a subtree whose bound equals $z$ cannot hold a
better value. In the library the test is the function
`detail::bb_prunable()`, i.e., lines 30 to 32 of Algorithm 13.3 below,
where the tolerance $\max \{ \tau_a, \tau_r |z| \}$ is computed by
`BBValueTraits::tolerance()` and is rounded down to
$\max \{ \tau_a, \lfloor \tau_r |z| \rfloor \}$ for an integer type
of values.

When the values of the solutions are integer, $\lceil \beta(n) \rceil$
is a bound too, and the test applied to it prunes each node with
$\beta(n) > z - 1$, i.e., those that could improve $z$ by a fraction
only (for an integer $z$ and zero tolerances, $\lceil \beta \rceil \ge
z$ if and only if $\beta > z - 1$). With doubles the same nodes are
branched, unless one gives an absolute tolerance just below 1. In the
library the rounding comes from the type of the values. A problem whose
`value_type` is an integer type returns integer bounds, rounded towards
the inside (down for a maximization, up for a minimization), and the
test above then prunes a node whose rounded bound equals the incumbent;
the test `bb-values` checks that, summed over its knapsacks, the integer
bounds evaluate no more nodes than the double ones.

### The running example: the 0-1 knapsack

The problem of the tests and of the benchmark is the 0-1 knapsack
`BasicKnapsack< V >` of [`Knapsack.h`][Knapsack.h], a maximization: given $n$
items of profit $p_j > 0$ and weight $w_j > 0$ and a capacity $C$, find $x \in
\{0, 1\}^n$ maximizing $\sum_j p_j x_j$ subject to $\sum_j w_j x_j \le
C$. Its items are sorted by decreasing ratio $p_j / w_j$, and a node
fixes the first $k$ of them (the field `k`, with the weight and the
profit of those taken). Its evaluation is the bound of the fractional
relaxation of Dantzig [[147]](Z-bibliography.md#ref-147). In it, the free
items are taken greedily in the order of the ratio while they fit, and
the first one that does not fit (the *critical* item) is taken in the
fraction that fills the capacity. As feasible value, the node gives that
of the greedy integer completion, which skips the items that do not fit
and goes on with the next ones. In pseudocode, with $u$ a node that
fixes the items $0, \ldots, k - 1$, $W(u)$ and $P(u)$ the weight and the
profit of the items that it takes, and the items numbered from 0, the
evaluation is the following (`evaluate()` and `greedy()` in
[`Knapsack.h`][Knapsack.h]).

**Algorithm 13.1** (The evaluation of a node of the knapsack).

```text
 1  function EVALUATE(u)               ▷ u fixes the items 0, …, k − 1
 2     if W(u) > C then return infeasible
 3     r ← C − W(u);  b ← P(u);  j ← k
 4     while j < n and w_j ≤ r do r ← r − w_j;  b ← b + p_j;  j ← j + 1
 5     if j < n then b ← b + p_j · r / w_j        ▷ j: the critical item
 6     if the values are integer then b ← floor(b)
 7     return the bound b and the feasible value GREEDY(u)
 8  function GREEDY(u)                 ▷ the greedy integer completion
 9     r ← C − W(u);  v ← P(u)
10     for j = k, …, n − 1 do
11        if w_j ≤ r then r ← r − w_j;  v ← v + p_j
12     return v
```

Both functions take $O(n - k)$ time, i.e., $O(n)$ per evaluation.

<a id="lem-13-1"></a>**Lemma 13.1** (The bound of Dantzig). Let a
node fix the items $0, \ldots, k-1$, with residual capacity $r \ge 0$
and profit $P$ of the items taken. Then the value $b$ computed by
Algorithm 13.1 is at least the profit of each completion of the node
that fits, and the greedy value is the profit of one of them; hence,
$v \le z^*(n) \le b$, where $z^*(n)$ is the best profit in $X(n)$.

*Proof.* We first remark that the best completion is at most $P$ plus
the optimum of the Linear Program $\max \{ \sum_{j \ge k} p_j y_j :
\sum_{j \ge k} w_j y_j \le r, \; 0 \le y \le 1 \}$. Let $y$ be an
optimal solution of the Linear Program, and suppose that $y_i < 1$ and
$y_j > 0$ for some $i < j$; take $i$ the least such index and $j$ the
largest one. Move the weight $\delta = \min \{ (1 - y_i)
w_i, \; y_j w_j \} > 0$ from $j$ to $i$, i.e., raise $y_i$ by $\delta /
w_i$ and lower $y_j$ by $\delta / w_j$: the solution stays feasible, and
its value changes by $\delta (p_i / w_i - p_j / w_j) \ge 0$, since $p_i
/ w_i \ge p_j / w_j$, so that it stays optimal. Each exchange either
fills $y_i$, which raises the least index, or empties $y_j$, which
lowers the largest one; hence, at most $n - k$ exchanges give an optimal
solution in which the items before some index $c$ are taken whole, the
item $c$ in a fraction and the items after it not at all. Its weight is
$r$ unless all the free items are taken whole, since otherwise raising
$y_c$ would increase the value ($p_c > 0$). These two properties define
the greedy fractional solution, whose value is $b - P$ (line 4 of
Algorithm 13.1 takes the items $k, k+1, \ldots$ while they fit, and the
critical item $j$ in the fraction that fills the residual capacity
left by them). Finally, the greedy integer completion fits by
construction, so its profit $v$ is the value of a point of $X(n)$.
$\square$

A node whose weight exceeds the capacity is infeasible. The branching
fixes the item $k$: the child that leaves it out is given first and the
child that takes it last, the latter only if it fits. With ties, a
depth-first search therefore goes on with the item taken, as the greedy
order does. Both children get as estimate the infinity of the
maximization, i.e., the bound of their parent. With `explicit_children`
they get their bound instead, and the item is taken also when it does
not fit, which gives children that are found infeasible at their
evaluation (the tests use this to produce infeasible leaves). With an
integer `V` the bound is rounded down, which is still a bound since the
profits are integer. Each node also has a hash (of `k`, the weight and
the profit, mixed by SplitMix64), an estimate (its greedy value), and
the two methods of `MPIBranchAndBound`. Finally, `optimum()` computes
the optimal value by the dynamic program on the capacity, in $O(n C)$
time, and the tests check each search against it. The function
`random()` draws an instance of one of the classes of Pisinger
[[148]](Z-bibliography.md#ref-148), i.e., uncorrelated ($p_j$ uniform in
$[1, R]$), weakly correlated ($p_j$ uniform in $[w_j - R/10,
w_j + R/10]$, at least 1) and strongly correlated ($p_j = w_j + R/10$),
with $w_j$ uniform in $[1, R]$ and the capacity half of the total
weight. The branching is the following (`branch()` in
[`Knapsack.h`][Knapsack.h]), where `sink` is the function of [§13.1](#sec-13-1)
and $\infty$ is the infinity of the maximization, i.e., an estimate that leaves
the child with the bound of its parent ([§13.2](#sec-13-2)).

**Algorithm 13.2** (The branching of a node of the knapsack).

```text
1  procedure BRANCH(u, sink)
2     if k = n then return                       ▷ every item is fixed
3     c ← u with the item k left out, and k + 1 items fixed
4     sink(c, the bound of c if explicit_children, else ∞)
5     if explicit_children or W(u) + w_k ≤ C then
6        t ← c with the item k taken:  W(t) ← W(u) + w_k;  P(t) ← P(u) + p_k
7        sink(t, the bound of t if explicit_children and t fits, else ∞)
```

The example of the chapter has the four items of
[Figure 13.1](#fig-13-1) and the capacity $C = 9$. Its optimum is 15
(the items 1 and 3, numbered from 1 in the figures), while the greedy
completion of the root takes the items 1 and 2 and gives 12. The root
bound is $5 + 7 + 10 \cdot 3 / 7 = 16.29$ (to two decimals, as in all
the figures). A node is named by the decisions it fixes, e.g., 10 for
"item 1 taken, item 2 left out".

![Figure 13.1. The tree that the best-first B&B explores on the
example (a maximization): under each node its bound $\beta$ (its
inherited bound for 1100, which is never evaluated) and the value $v$
of its greedy completion, and above it the order of the evaluations.
Violet nodes are evaluated; the dark violet node is the node of the optimum,
which its bound prunes, the rust node is pruned without evaluation,
and the faint nodes are not generated, since their item
does not fit.](../figures/branch-and-bound-tree.svg){#fig-13-1}

## 13.3 The sequential search {#sec-13-3}

The sequential search `BranchAndBoundSearch`
([`BranchAndBound.h:879`][BranchAndBound.h:879], or the function
`branch_and_bound( p , options )`, whose declarations, without their bodies,
are

```cpp
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
branch_and_bound( const P & problem ,
                  const BasicBBOptions< bb_value_t< P > > & o = {} );
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
branch_and_bound( const P & problem ,
                  std::vector<
                    BBLeaf< typename P::node_type , bb_value_t< P > > > seeds ,
                  const BasicBBOptions< bb_value_t< P > > & o = {} );
```

for a model `P` of `BranchAndBound`) keeps the open nodes
in a queue Q, ordered by the key of the selection
([§13.5](#sec-13-5)), and, with the plunging selection, a node $d$ of
the dive, which is extracted next. Its pseudocode, which is that of the
Doxygen comment of [`BranchAndBound.h`][BranchAndBound.h], is the following;
$R$ is the set of the nodes it starts from (the root, unless the caller gives
others, [§13.11](#sec-13-11)).

**Algorithm 13.3** (Sequential branch-and-bound).

```text
 1  procedure B&B(R)          ▷ the nodes to start from, {root} by default
 2     z ← the given incumbent, or ∞;  Q ← R;  d ← nil  ▷ d: the dive
 3     while Q ≠ ∅ or d ≠ nil do
 4        if the node or the time limit is reached then break
 5        βmin ← the least bound of Q and of d
 6        if PRUNABLE(βmin) then the nodes of Q and d are leaves;  break
 7        n ← d if d ≠ nil, else the first node of Q, removed;  d ← nil
 8        K ← PROCESS(n)
 9        if the selection is plunging and K ≠ ∅ then
10           k ← the node of K of least bound
11           βmin ← min(the least bound of Q, β(k))
12           if z = ∞ or β(k) ≤ βmin + γ (z − βmin) then d ← k, from K
13        push the nodes of K in Q
14     the nodes left in Q and d are leaves             ▷ open by a limit
15     return z, its node, and the bound min(z, least bound of the leaves)
16  procedure PROCESS(n)
17     if PRUNABLE(β(n)) then n is a leaf (pruned);  return ∅
18     cutoff(z);  (β', infeasible, v) ← evaluate(n)
19     if infeasible then n is a leaf (infeasible);  return ∅
20     β(n) ← max(β(n), β')
21     if v is given and v < z then z ← v;  n is the node of z
22     if PRUNABLE(β(n)) then n is a leaf (pruned or feasible);  return ∅
23     K ← ∅;  cutoff(z)
24     for each child m of n, with its estimate e, given by branch(n) do
25        β(m) ← max(β(n), e)
26        if PRUNABLE(β(m)) then m is a leaf (pruned)
27        else K ← K ∪ {m}
28     if n has no child then n is a leaf
29     return K
30  procedure PRUNABLE(β)
31     if z = ∞ then return β = ∞
32     return β ≥ z − max(τ_a, τ_r |z|)
```

A node is thus tested three times ([Figure 13.2](#fig-13-2)): (i) when
it is extracted (line 17), against an incumbent that may have improved
since its generation, (ii) after its evaluation (line 22), against the
bound it has just obtained, and (iii) when it is generated (line 26),
against its inherited bound. A node that survives the second test is
branched, also when it has a feasible value, since its subtree may hold
better values than the one found. Line 6 is the termination test: as
soon as the least bound of the open nodes cannot improve $z$, neither
can any of them, and all of them would be pruned at their extraction.
All the searches of the library but the MultiQueue, which does not
know the least bound of its queues ([§13.8](#sec-13-8)), use this test,
and they also end when Q is empty or a limit is reached; as in the
pseudocode, the emptiness of line 3 is tested before the limits of line
4, so that a search whose last evaluation reaches the node limit with no
open node left ends as optimal. With best-first the test at line 6 often
(but not always) ends the search while Q still holds many nodes, i.e.,
those generated before the incumbent was found and never extracted
since.

![Figure 13.2. The procedure PROCESS of Algorithm 13.3: the three
tests of pruning (i), (ii) and (iii), at lines 17, 22 and 26, and the
leaves that each one makes. With $z = +\infty$ a test prunes only a
bound $+\infty$.](../figures/branch-and-bound-process.svg){#fig-13-2}

The procedure PROCESS is `detail::bb_process()`
([`BranchAndBound.h:808`][BranchAndBound.h:808]). All the searches of the
library share it, each with the incumbent it sees (a function `inc()` that
returns $z$ in the minimization) and with the function `improve` by which it
offers a feasible value; the children that survive line 26 go to a function
`push` of the caller, which puts them in its queue or in its messages.
The leaves go to a record of the worker, `detail::BBFrontier`
([`BranchAndBound.h:700`][BranchAndBound.h:700]), which keeps the least bound
`low` of the leaves that lower the bound of the search and, only with
`keep_frontier`, the leaves themselves.

**Algorithm 13.4** (The record of the leaves of a worker).

```text
1  procedure LEAF(n, β, kind, lowers)          ▷ a leaf n of bound β
2     if lowers then low ← min(low, β)
3     if keep_frontier then append (n, β, depth(n), kind) to the leaves
4  procedure LEAVES(Q, kind)                   ▷ all the nodes of a queue
5     low ← min(low, the least bound of Q)
6     if keep_frontier then append each node of Q with its bound and depth
7     Q ← ∅
8  procedure MERGE(L')                         ▷ the record of another worker
9     low ← min(low, low of L');  append the leaves of L'
```

In PROCESS, an infeasible leaf (line 19) and a node with no children
(line 28) are recorded with `lowers` false: the former has no element,
and the latter is solved by Assumption (A3) of [§13.4](#sec-13-4), its
value having been offered to $z$ at line 21. The leaves of the three
tests are of the kind `pruned`, or `feasible` at line 22 when the node
gave a feasible value, and those left by a limit (line 14) of the kind
`open`. A child pruned at line 26 is copied into the record only with
`keep_frontier`; otherwise only its bound enters `low`.

[Figure 13.3](#fig-13-3) traces Algorithm 13.3 on the example. The
greedy completion of the root gives $z = 12$ at once. Then the search
goes down the right branch, whose nodes inherit the bound $16.29$ (the
ties among equal keys are broken in favor of the last node pushed,
cf. [§13.5](#sec-13-5)); the node 110 has the bound $13.71$, which is
still above $z = 12$, and is branched. The node 10 gives the optimum 15
and is pruned by its bound, and the node 0 is pruned after its
evaluation. Finally, the open node 1100, whose inherited bound is
$13.71 \le 15$, ends the search at line 6 without being evaluated; six
nodes are evaluated in all.

![Figure 13.3. The steps of Algorithm 13.3 on the example, one row
per evaluation: the node extracted, the outcome of PROCESS, the queue
after the step with the nodes pushed by it in rust (each with its bound;
the first node of the queue is on the left), and the incumbent. In this
maximization, the test of pruning reads
$\beta \le z$.](../figures/branch-and-bound-steps.svg){#fig-13-3}

## 13.4 Correctness and the reported bound {#sec-13-4}

We now prove that Algorithm 13.3 returns a value within the tolerance of
the optimum when no limit stops it, and that the bound it reports is a
bound on the optimum whatever the status. We need four assumptions on
the problem, which the concept cannot check:

* (A1) the children of a node cover its subset:
  $X(n) \subseteq X(n_1) \cup \cdots \cup X(n_k)$;
* (A2) the bound of an evaluation, the estimate of a child and the
  bound given with each node of $R$ are lower bounds:
  $\beta'(n) \le v(x)$ for every $x \in X(n)$, $e(m) \le v(x)$ for
  every $x \in X(m)$, and likewise for the nodes of $R$ (the root starts
  with $-\infty$);
* (A3) a node has no child only if it is solved, i.e., if its subset is
  empty or its evaluation gives a feasible value $v$ with
  $v \le v(x)$ for every $x \in X(n)$;
* (A4) every chain of branchings is finite, e.g., because the subsets
  shrink strictly and $X$ is finite.

The knapsack satisfies them, as we now show; in its case $X(u)$ is the
set of the $x \in \{0, 1\}^n$ that agree with the node $u$ on the items
$0, \ldots, k - 1$ and fit in the capacity.

<a id="lem-13-2"></a>**Lemma 13.2** (The knapsack satisfies the
assumptions). With or without `explicit_children`, and with values of a
floating point or of an integer type, the knapsack of Algorithms 13.1
and 13.2 satisfies (A1) to (A4), the bounds being read in the
minimization of the negated profits.

*Proof.* (A1) Let $x \in X(u)$ with $k < n$. If $x_k = 0$, then
$x \in X(c)$. If $x_k = 1$, then $W(u) + w_k \le \sum_j w_j x_j \le C$,
so that the child $t$ is generated at line 5 of Algorithm 13.2, and
$x \in X(t)$. (A2) If $W(u) > C$, then $X(u) = \emptyset$ and the node
is reported infeasible; otherwise the bound of Algorithm 13.1 is at
least the profit of each $x \in X(u)$ by [Lemma 13.1](#lem-13-1), and
so is its floor for an integer type, since the profits are integer. The
estimate $\infty$ of a child is a trivial bound, and with
`explicit_children` the estimate is the bound of the child, or $\infty$
if it does not fit. (A3) BRANCH gives no child only when $k = n$; then
$X(u)$ holds at most the vector of the node, and if it fits its profit
is $P(u)$, the greedy value of line 7 of Algorithm 13.1, while if it
does not fit the node is infeasible. (A4) Each child fixes one item
more than its parent, and therefore every chain of branchings has at
most $n$ steps. $\square$

We also need the threshold to behave well when $z$ improves.

<a id="lem-13-3"></a>**Lemma 13.3** (The threshold is monotone). If
$\tau_r \le 1$, the threshold $\theta(z) = z - \max \{ \tau_a,
\tau_r |z| \}$ is nondecreasing in $z$, and so is its version for an
integer type, $z - \max \{ \tau_a, \lfloor \tau_r |z| \rfloor \}$.

*Proof.* One has $\theta(z) = \min \{ z - \tau_a, \; z - \tau_r |z| \}$,
the minimum of two functions of $z$. The first is increasing, and the
second has slope $1 - \tau_r \ge 0$ for $z > 0$ and $1 + \tau_r$ for
$z < 0$; hence, both are nondecreasing, and so is their minimum. For
the integer version, a unit increase of $z$ changes
$\lfloor \tau_r |z| \rfloor$ by at most $\lceil \tau_r \rceil \le 1$,
so that $z - \lfloor \tau_r |z| \rfloor$ does not decrease, and the
same argument applies. $\square$

The condition $\tau_r \le 1$ is not checked by the library. With a
larger relative tolerance the threshold of a positive incumbent falls
when the incumbent rises, which one hardly wants (a relative gap of more
than 100% accepts any value of the same sign). We call *open*, at a
given moment of Algorithm 13.3, the nodes of Q and the node $d$ of the
dive.

<a id="lem-13-4"></a>**Lemma 13.4** (The invariant). Under (A1) to
(A3) and with $\tau_r \le 1$, at the start of each iteration of the
loop at line 3 every $x \in X(R)$ with $v(x) < \theta(z)$ belongs to
$X(n)$ for some open node $n$, where $X(R)$ is the union of the subsets
of the nodes of $R$.

*Proof.* We proceed by induction on the iterations. At the start the open nodes
are those of $R$, and the claim holds. Assume it at the start of an
iteration, which either breaks (and then there is no next iteration)
or removes a node $n$ from the open ones at line 7 and processes it.
First, $z$ can only decrease (line 21), and then so does $\theta(z)$ by
[Lemma 13.3](#lem-13-3); hence, the set of the $x$ to be covered can
only shrink. Second, the elements of $X(n)$ below the threshold are
covered again after PROCESS, in each of its four outcomes. If $n$ is
pruned at line 17 or 22, then $\beta(n) \ge \theta(z)$; since
$\beta(n)$ is the maximum of $\beta'(n)$ and of bounds of supersets of
$X(n)$, it is a bound of $X(n)$ by (A2), and no element of $X(n)$ is
below the threshold. If $n$ is infeasible (line 19), $X(n)$ is empty.
If $n$ has no child, each $x \in X(n)$ has
$v(x) \ge v \ge z \ge \theta(z)$ by (A3), where $v$ is its feasible
value, offered at line 21. Otherwise, by (A1), each $x \in X(n)$ is in
the subset of a child $m$, which is either pruned at line 26 (and then
$v(x) \ge \beta(m) \ge \theta(z)$ by (A2)) or put in Q or in $d$ at
lines 12 and 13. $\square$

<a id="thm-13-5"></a>**Theorem 13.5** (Correctness and finiteness).
Under (A1) to (A4) and with $\tau_r \le 1$, Algorithm 13.3 performs a
finite number of iterations. If no limit stops it, then either
$z < +\infty$ and

$$
z - z^*_R \le \max \{ \tau_a, \; \tau_r |z| \} ,
$$

where $z^*_R = \min \{ v(x) : x \in X(R) \}$, or $z = +\infty$ and
$X(R)$ is empty (when no initial incumbent is given).

*Proof.* We note that each iteration that does not break removes one
node from the open ones and adds only children of it; hence, each node
of the trees below the nodes of $R$ is processed at most once for each
node of $R$ above it. The set $R$ is finite, and by (A4) and by
König's lemma so are these trees, since each node has finitely many
children (`branch()` returns); hence, so is the number of iterations.
If no limit
stops the search, it ends at line 3 with no open node, or at line 6 with
each open node of bound at least $\theta(z)$. In both cases [Lemma
13.4](#lem-13-4) shows that no $x \in X(R)$ has $v(x) < \theta(z)$,
i.e., $z^*_R \ge \theta(z)$, which is the claim. If $z = +\infty$, then
$\theta(z) = +\infty$ in the sense of the test at line 31, and no $x$
exists. $\square$

Note that with an initial incumbent value below the optimum, the search
ends with that value and no node (the test `bb` checks this case). The
value $z$ of [Theorem 13.5](#thm-13-5) is in any case a value found or
given, and the inequality bounds it from above. When a limit stops the
search, or when the tolerances are positive, the leaves give a bound on
the optimum, which is the bound of `BBSolution`. The *leaves* of a
search are the nodes it does not branch: those pruned at the three
tests, the infeasible ones, those with no child, and those left open by
a limit (lines 6 and 14). Algorithm 13.4 records the least bound `low`
of those that lower the bound, i.e., all but the infeasible ones and
those with no child, and it keeps the leaves themselves only with
`keep_frontier`.

<a id="thm-13-6"></a>**Theorem 13.6** (The reported bound). Under (A1)
to (A3), whatever the status of the search,

$$
z^*_R \ge L = \min \{ \, z, \; \min \{ \beta(l) : l
\text{ a leaf that lowers the bound} \} \, \} ,
$$

and therefore the gap $z - L$ reported by `BBSolution` bounds
$z - z^*_R$.

*Proof.* If $X(R)$ is empty there is nothing to prove. Otherwise we let
$x^*$ be optimal over $X(R)$. By (A1), the subsets of the leaves cover
$X(R)$; the argument is that of [Lemma 13.4](#lem-13-4) with no
threshold, since each node of the tree is either branched (and then
covered by its children) or a leaf, the nodes still open when the search
ends being leaves of kind `open`. Let $l$ be a leaf with $x^* \in X(l)$.
It is not infeasible, since $X(l) \ne \emptyset$. If it has no child,
then $z^*_R = v(x^*) \ge v \ge z \ge L$ by (A3), where $v$ is its
feasible value, offered to the incumbent. Otherwise it lowers the bound,
and $\beta(l) \le v(x^*)$ by (A2), so that again $z^*_R \ge L$.
$\square$

The leaf that has no child is excluded from $L$ because its subtree is
the node, whose value is offered to $z$; the infeasible leaf is excluded
because its bound is $+\infty$. Note that the bound $\beta(l)$ of a leaf
pruned at its generation or at its extraction is an inherited one, which
is usually weaker than the bound its evaluation would give. The reported
bound is therefore valid but not necessarily the best one available; it
is exact (equal to $z$) when the search ends with zero tolerances and no
limit, since then each leaf has $\beta(l) \ge z$. We remark that the
theorem is stated for Algorithm 13.3, which prunes by its own $z$. With
`shared_incumbent` the search prunes by the least $c$ of $z$ and of the
values read, which are not offered to $z$, and the same proof gives
only $z^*_R \ge \min \{ L, c \}$: a leaf with no child whose value is
not better than $c$ is in neither $z$ nor $L$ (this is why the bound of
Exercise 13.2 is at most the optimum). On the example the
bound is $\max \{ 15, 14.14, 13.71, 15 \} = 15$ in the maximization, and
the gap is 0.

## 13.5 The selection rules {#sec-13-5}

[Lemma 13.4](#lem-13-4) holds whatever the order in which the open nodes
are extracted. Therefore the order is a policy, the `Selection` of the
options ([`BranchAndBound.h:137`][BranchAndBound.h:137]), which decides how
many nodes are evaluated and how much memory the queue takes, while the value
returned is within the tolerance anyway. In the queue we order the nodes by a
key: the bound with best-first and plunging, minus the depth with
depth-first, the depth with breadth-first, and the estimate of the
problem with best-estimate (the bound, if the problem gives no
estimate). Among equal keys the last node pushed comes first, but with
breadth-first, where the first comes first. [Figure 13.4](#fig-13-4)
shows the five rules on the example.

![Figure 13.4. The five selections on the example: in each tree
the number is the order of evaluation, the dark violet node is the node of
the optimum, a rust node is pruned without evaluation and a faint one
is not generated. In the plunging tree the thick violet edges are the
dive
from the root to 110; its child 1100 is not dived into, since its bound
13.71 is farther than $\gamma = 0.5$ of the gap from the least open
bound.](../figures/branch-and-bound-selection.svg){#fig-13-4}

### Best-first

With best-first the node of least bound is extracted first. When the
bound of a node does not depend on the order of the search, best-first
evaluates the fewest nodes, in the following sense.

<a id="prop-13-7"></a>**Proposition 13.7** (Best-first evaluates the
necessary nodes, and "few" others). Assume (A1) to (A4), zero tolerances,
no limit, no initial incumbent below $z^*$, and that the evaluation of
a node depends only on the node (no `cutoff`). Call a node *critical*
if its inherited bound and the evaluated bounds of all its ancestors
(taken as $+\infty$ for an infeasible one) are less than $z^*$. Then (i)
each selection evaluates each critical node, and (ii) best-first
evaluates only nodes whose inherited bound is at most $z^*$.

*Proof.* (i) We argue by induction along the path from the root. At the
start the root is open. A critical node $m$ that is open is never pruned
at line 6, 17 or 26, since its bound is less than $z^* \le z$; since the
search ends with no open node of bound below $z$ ([Theorem
13.5](#thm-13-5)), $m$ is extracted and evaluated. If $m$ has a critical
child, its evaluated bound is less than $z^*$ by the definition, and it
is not pruned at line 22; it is therefore branched, and its critical
children, whose inherited bounds are less than $z^* \le z$, survive line
26 and are open. (ii) While $z > z^*$, some open node $m$ holds an
optimal $x^*$ ([Lemma 13.4](#lem-13-4), since $v(x^*) = z^* < z$ and the
tolerance is 0), and $\beta(m) \le z^*$ by (A2); best-first extracts a
node of least bound, which is therefore at most $z^*$. Once $z = z^*$, a
node extracted with $\beta \ge z^*$ is pruned at line 17 before its
evaluation. $\square$

Hence, the nodes evaluated by best-first are the critical ones plus some
of those whose bound equals $z^*$, and how many of the latter depends on
the ties (the test `bb` checks the optimum of every selection, not the
number of the evaluations). This result is the counterpart of the
optimality of A$^*$ among the algorithms that use the same heuristic
([Chapter 2](02-best-first.md#ch-2)), and it fails in two ways that are
common in practice: with `cutoff`, the evaluation of a node depends on
the incumbent, and with a relaxation solved by a warm start from the
previous node it depends on the order. Its price is the memory: the
queue may grow to a "large" part of the tree before an incumbent is found,
since each open node whose bound is below the incumbent of the moment
stays in it until it is extracted.

### Depth-first branch-and-bound

With depth-first the deepest node is extracted first, which typically
finds feasible values early and keeps the queue "small". This is the
depth-first B&B (DFBnB) of the literature of heuristic search, i.e., a
depth-first search that prunes each node whose bound is not less than
the best value found so far. It is an *anytime* algorithm, since an
incumbent is available (through `on_incumbent`) long before the search
ends. Since the children of a node have the same depth and the ties are
broken in favor of the last pushed, the child that a problem gives last
to the sink is the first one extracted. A problem that wants its
children tried by increasing bound, as DFBnB with sorted successors
does, gives them by decreasing bound (the knapsack gives the child that
takes the item last, which is the child of the greedy order). The
version of DFBnB on the graph of a `Domain`, without a queue, is in
[Chapter 4](04-memory-bounded.md#ch-4). Here, instead, the search works
on any `BranchAndBound` problem, e.g., on the tree of the partial tours
of the traveling salesperson problem that the benchmark uses
([§13.13](#sec-13-13)); there the first complete tour is reached after
one evaluation per level of the tree, and from then on it prunes the
partial tours whose bound is not less than its cost.

<a id="prop-13-8"></a>**Proposition 13.8** (The memory of depth-first).
With the depth-first selection, started from the root only, Q never
holds more than $(b - 1) D + 1$ nodes, where $b$ is the largest number
of children of a node and $D$ the depth of the tree.

*Proof.* We show that, at the start of each iteration, the nodes of Q
at each depth $\delta \ge 1$ are children of one node, and that Q holds
at most $b - 1$ of them at each depth but its deepest one, where it
holds at most $b$. This holds at the start, when only the root is in Q. An
iteration extracts a node $n$ of the deepest level $\delta$ of Q, which
leaves at most $b - 1$ nodes at depth $\delta$, all siblings of $n$.
It pushes at most $b$ children of $n$ at depth $\delta + 1$, where Q
had no node. The levels above $\delta$ are unchanged, and each of them
lost one node (the ancestor of $n$ at that depth) when it was the
deepest; hence, each holds at most $b - 1$ nodes. Summing over the
depths $1, \ldots, D$ gives $(b - 1)(D - 1) + b = (b - 1) D + 1$.
$\square$

Depth-first pays in the evaluations: it may evaluate "many" nodes that
best-first would prune, i.e., those whose bound lies between $z^*$ and
the incumbent of the moment (the depth-first tree of [Figure
13.4](#fig-13-4) evaluates 1100, which best-first prunes unevaluated).

### Breadth-first, best-estimate and plunging

*Breadth-first* extracts the shallowest node, the oldest among those
of the same depth. It is rarely the right choice for a B&B whose
feasible values come only from the leaves (on the example it is the
worst of the five, with eight evaluations), but it gives the nodes of a
level together, which some uses may want. *Best-estimate* extracts the
node of least estimate of the best feasible value in its subtree, which
the problem gives by `estimate( n )` (for the knapsack, its greedy
value). It aims at finding a "good" incumbent early, while the bound
still drives the pruning and the termination test.

*Plunging* is the usual compromise: the best child of the node just
branched is extracted next (a *dive*) as long as its bound is within
the fraction $\gamma$ (`plunge_gap`) of the gap from the least open
bound $\underline{\beta}$, i.e., while

$$
\beta(k) \le \underline{\beta} + \gamma \, (z - \underline{\beta}) ,
$$

and always while there is no incumbent; when the dive stops, the search
goes back to best-first (lines 9 to 12 of Algorithm 13.3). With $\gamma
= 0$ a dive continues only into a child whose bound equals the least
open bound. With $\gamma = 1$ it continues into each child that is not
pruned, which makes the search depth-first in the order of the bound of
the children and best-first among the nodes left behind (the test `bb`
runs both extremes). One can expect a dive to be "cheap" when the
evaluator keeps the state of the last node, e.g., a relaxation solved
again from the basis of its parent. In the code
([`BranchAndBound.h:1011`][BranchAndBound.h:1011]) the child $k$ of line 10 is
the last of the children of least bound, and the test of line 12 is
`BBValueTraits::within()`, which computes the right-hand side in the type of
the values for a floating point type and in `long double` for an integer one,
which cannot overflow. The test is false when $\underline{\beta}$ is not finite
(e.g., $-\infty$ for a node given back with an unknown bound,
[§13.11](#sec-13-11)), and then no dive is taken while an incumbent exists.

### The queue

The open nodes are kept in a `BBQueue`
([`BranchAndBound.h:505`][BranchAndBound.h:505]). Its nodes are in an array of
slots, each holding a node with its bound and its depth, with a free list,
since a node is destroyed when it is branched (in A$^*$ a node stays in CLOSED,
while here no structure refers to a node after its processing). Its order is a
pair of indexed 4-ary heaps on the indices of the slots, taken from the default
data structures of the library ([Chapter 3](03-data-structures.md#ch-3)): one
by the key of the selection, with the number of the push to break the
ties, and one by the bound. With the second heap one has the least open
bound in $O(1)$, which the termination test at line 6 and the reported
bound need when the key is not the bound. Each change, the pop
included, is made on both heaps at the places of the index, and the
heaps grow by doubling with the slots. [Figure 13.5](#fig-13-5) shows a
pop with the depth-first selection.

![Figure 13.5. The queue of the open nodes before and after a pop
with the depth-first selection (the key is minus the depth, the ties
broken by the last push). Here the pop returns the node d of the slot 4, the
top of the heap by key, removes the index 4 from both heaps at its
places, and puts the slot on the free list, where the next push will
take it; the least bound, at the top of the other heap, is still 4.
The 4-ary heaps hold here at most four indices, i.e., a root and its
children.](../figures/branch-and-bound-queue.svg){#fig-13-5}

A push takes a free slot (or a new one), stores the node and pushes its
index in both heaps, the key paired with the number of the push,
complemented unless the order is first in, first out, and a pop takes
the index of least key out of both heaps and frees its slot. In
pseudocode, with $H_K$ the heap of the keys and $H_B$ that of the
bounds, both on the indices of the slots, and $s$ the number of the
pushes so far:

**Algorithm 13.5** (The queue of the open nodes).

```text
 1  procedure PUSH(n, β, δ, key)          ▷ a node, its bound, depth, key
 2     i ← a slot of the free list, removed, or a new slot
 3     if i is beyond the capacity c of the heaps then
 4        rebuild both heaps with capacity 2i + 16, inserting their indices
 5     s ← s + 1;  slot[i] ← (n, β, δ)
 6     insert i in H_K with the key (key, s) if breadth-first, else (key, −s)
 7     insert i in H_B with the key β
 8  function POP()
 9     i ← the top of H_K;  remove i from H_K and from H_B
10     put i on the free list;  return the content of slot[i], emptied
11  function BESTBOUND()
12     return the key of the top of H_B, or ∞ if the queue is empty
```

Here $-s$ stands for the complement of $s$ in 64 bits, which orders the
pushes in the opposite way. The costs are the following.

<a id="prop-13-9"></a>**Proposition 13.9** (The queue). Let $N$ be the
number of the open nodes and $M$ the largest number of open nodes at
any time. Then (i) the pairs (key, $\pm s$) are distinct, so that POP
returns a well-defined node, the least key first and, among equal keys,
the last pushed (the first pushed with breadth-first); (ii) BESTBOUND
and the first node take $O(1)$ time, POP takes $O(\log N)$ time, and
PUSH takes $O(\log N)$ time but for the rebuilding of line 4, whose
total cost over a sequence of pushes is $O(M \log M)$, i.e., $O(\log M)$
amortized per push when at least $M$ pushes are made; (iii) the queue
uses $O(M)$ memory, besides the nodes themselves.

*Proof.* (i) Each push takes a new value of $s$, and the order on the
pairs is lexicographic. (ii) In a 4-ary heap of $N$ elements, an
insertion moves an element up at most $\lceil \log_4 N \rceil$ levels
with one comparison per level, and the removal of an element whose
position is known (the heaps are indexed by the slots) moves the last
element up or down along one path of at most $\lceil \log_4 N \rceil$
levels with at most 4 comparisons per level; a push makes one insertion
in each heap and a pop one removal in each, and the free list is a
stack. The top of a heap is its first element. The rebuilding of line 4
happens when the index $i$ reaches the capacity $c$, and sets it to
$2i + 16 > 2c$; the capacities are therefore at least doubled, and since
$i < M + 1$ (a slot is new only when all the slots hold nodes), the
rebuildings insert at most $M$ indices each, with capacities
$c_1 < c_2 < \cdots$ at least doubling up to $O(M)$, i.e.,
$O(\sum_t c_t \log M) = O(M \log M)$ in all. (iii) The slots are at most
$M + 1$, and the heaps have capacity $O(M)$. $\square$

The move of one node per push and per pop is apart, and so is the
destruction of the node of a slot. These costs are usually negligible
when the evaluation
solves a relaxation, which is the case the layer is meant for. They
matter only for problems whose bound costs about as much as a "few"
comparisons; the knapsack of the tests is one, and it is therefore a
test of the overheads of the schemes rather than of their intended
use. One could also keep the bounds in a heap with lazy deletion, in
which the indices popped by key stay until they reach its top: this
saves the second walk of a pop, but the least bound must then skip the
stale entries, while the indexed pair keeps the bound of each index
exact.

We can now bound the cost of a whole sequential search in terms of the
work of the problem. In the following, an *evaluation* costs at most
$T_e$ time and a branching, i.e., a call of `branch()`, at most $T_b$
time besides the calls of `sink`.

<a id="prop-13-10"></a>**Proposition 13.10** (The cost of the sequential
search). Let Algorithm 13.3 process $E$ nodes (line 8), generate $G$
children, and let $M$ be the largest number of its open nodes. Then it
makes at most $E$ evaluations and $E$ branchings, and takes
$O(E (T_e + T_b + \log M) + G \log M)$ time besides the leaves
recorded at its end (lines 6 and 14, $O(M)$ time) and the calls of
`shared_incumbent` and of the clock (one per iteration). It uses $O(M)$
memory for the open nodes, plus the leaves kept with `keep_frontier`.

*Proof.* We note that each iteration of the loop at line 3 that does not
break extracts one node (line 7), by POP in $O(\log M)$ time
(Proposition 13.9) or from the dive, and processes it. PROCESS makes at
most one evaluation (line 18) and one branching (line 24), and for each
child it makes $O(1)$ work and either a record in $O(1)$ (Algorithm
13.4, amortized for the append) or a push, in $O(\log M)$ amortized
time; with plunging, lines 9 to 12 scan the children once. The least
bound of line 5 is BESTBOUND, in $O(1)$. The memory holds the queue
(Proposition 13.9), the node of the dive, the children of one node and
the record of the leaves. $\square$

Note that $E \le |R| + G$, since each node processed was pushed once,
and that $G \le b E$ when no node has more than $b$ children; hence the
overhead of the search is $O(b \log M)$ per evaluation, which is
negligible when $T_e$ is the time of a relaxation.

## 13.6 A$^*$ as a special case {#sec-13-6}

Nau, Kumar and Kanal [[17]](Z-bibliography.md#ref-17) formulated a
general B&B of which A$^*$ (and AO$^*$) are special cases. In our terms,
A$^*$ on a tree is the B&B in which a node is a path from $s_0$ to a
state $u$, and its subset is the set of the paths to a goal that extend
it. Its bound is $f(u) = g(u) + h(u)$, which is a lower bound on the
cost of these paths when $h$ is admissible, and a goal gives the
feasible value $g(u)$. In the benchmark of the chapter we use this
construction (`TourTree` in
[`tools/BranchAndBoundBench.cpp`][BranchAndBoundBench.cpp]) on the traveling
salesperson problem of [Chapter 14](14-domains.md#ch-14), whose search space is
a tree: `evaluate()` returns $g + h$ as the bound and $g$ as the feasible value
of a goal, and `branch()` gives each successor with its $g + h$ as estimate.

<a id="prop-13-11"></a>**Proposition 13.11** (A$^*$ as best-first B&B).
On the tree of the paths from $s_0$, with an admissible heuristic,
$h = 0$ on the goals, nonnegative edge costs, zero tolerances and the
best-first selection, the first goal that Algorithm 13.3 evaluates has
the optimal cost $C^*$, and the search ends at the next test of line 6
with no further evaluation.

*Proof.* We first observe that, while no goal of cost $C^*$ has been
evaluated, $z > C^*$ and some open node $m$ lies on an optimal path
([Lemma 13.4](#lem-13-4)); its bound is the largest $f$ along that path
up to $m$ (by the pathmax of [§13.2](#sec-13-2)), which is at most $C^*$
by admissibility. Let $u$ be the first goal evaluated; when it was
extracted it had the least bound among the open nodes, so that $\beta(u)
\le \beta(m) \le C^*$, while $\beta(u) \ge f(u) = g(u) \ge C^*$ since
$h(u) = 0$. Hence $g(u) = C^*$, the evaluation sets $z = C^*$, and the
least bound of the open nodes is at least $\beta(u) = C^*$ (best-first
extracted $u$ as a node of least bound, and the children pushed since
then have bounds of at least that of their parent). The test at line 6
then prunes all of them. $\square$

We remark that the proof shows the role of pathmax. With an inconsistent
heuristic $f$ may decrease along a path, but the inherited bound does
not, and best-first on the inherited bound extracts the nodes in
nondecreasing order of their key. Still, the two searches differ in two
respects. First, a graph has *duplicates*, i.e., several paths to the
same state, of which only the cheapest matters. A$^*$ detects them with
CLOSED and keeps one node per state, while a B&B tree has none and needs
no closed list; on a graph with many duplicates the B&B view enumerates
the paths, and A$^*$ is the right tool. Second, A$^*$ stops at the first
goal extracted, while a B&B goes on until no open node can improve the
incumbent. In a general B&B the first feasible value is often not
optimal, and a node with a feasible value may hold better values below,
as the root and the node 1 of [Figure 13.1](#fig-13-1) do (value 12, and
the optimum 15 below). This is why the two layers of the library share
the parallel schemes and the executors, but not the node table and the
OPEN lists of A$^*$. For the same reason SafePBNF and PA$^*$SE ([Chapter
11](11-safepbnf.md#ch-11) and [Chapter 8](08-shared-memory.md#ch-8))
have no counterpart here: they exist to handle the duplicates and the
re-expansions, which a tree does not have.

## 13.7 What a parallel search has to guarantee {#sec-13-7}

We now turn to the parallel B&B searches of the library, all of which
build one tree with $p$ workers and differ in how they share the open
nodes and the incumbent. In the classification of Gendron and Crainic
[[16]](Z-bibliography.md#ref-16) they are of the second type, i.e., the
parallelism is in the exploration of the tree. Any parallelism inside
the evaluation of a node (their first type) is left to the problem, and
the search of several trees at once (their third type) is not provided.
Within the second type the searches are asynchronous, but for the
deterministic rounds. They keep either one pool of open nodes (the
shared queue, the MultiQueue, the farm and the rounds,
[§13.8](#sec-13-8)) or one pool per worker (the distributed queues, the
processes of MPI, the work stealing and the SelfSplit,
[§13.9](#sec-13-9)).

Note that the invariant of the sequential search holds for a parallel
one if the open nodes are taken in a wider sense: the nodes in any
queue, those held by a worker (being processed, or the next node of a
dive), and those traveling in a message, a task, an offer of a steal or
a slot of a round. Since the incumbent is shared, each worker prunes by
the value it sees, which may be older than the current one; the
following lemma shows that this is safe. The incumbent is an object of
the class `BBIncumbent`
([`ParallelBranchAndBound.h:190`][ParallelBranchAndBound.h:190]), whose two
operations are lines 1 to 8 of the following algorithm; lines 9 to 17 give the
dive and the test of the limits, which all the parallel searches share and
which are described at the end of this section.

**Algorithm 13.6** (The shared incumbent, the dive and the limits of the
parallel searches).

```text
 1  function VALUE()                       ▷ read by the workers to prune
 2     return z, an atomic read with relaxed order if the values have a
          lock-free atomic type, else a read under the lock M
 3  function IMPROVE(v, n)                 ▷ the feasible value v of n
 4     if the values are atomic and v ≥ VALUE() then return false
 5     lock M
 6     if v ≥ z then unlock M;  return false
 7     z ← v;  the node of z ← n;  call on_incumbent(v);  unlock M
 8     return true
 9  function DIVE(K, bl, z)                ▷ bl: the least open bound known
10     if K = ∅ then return nil
11     if depth-first then return the last node of K
12     if not plunging then return nil
13     k ← the last node of K of least bound;  bl ← min(bl, β(k))
14     if z = ∞ or (bl is finite and β(k) ≤ bl + γ (z − bl)) then return k
15     return nil
16  function LIMIT(c)                      ▷ c: the nodes counted so far
17     return (node_limit > 0 and c ≥ node_limit) or the time limit has
          passed, setting the flags of the limit (and of the time) if so
```

<a id="lem-13-12"></a>**Lemma 13.12** (The shared incumbent). In
`BBIncumbent` ([`ParallelBranchAndBound.h:190`][ParallelBranchAndBound.h:190]),
(i) the value never increases, (ii) each value it takes is the initial one or a
feasible value offered by a worker, and (iii) the stored node is a node of the
stored value whenever the lock is free, and the callback `on_incumbent`
is called with strictly decreasing values. Hence, a worker that prunes
by a value it read earlier prunes a subset of what it would prune by
the current one.

*Proof.* We note that the value changes only at line 7 of Algorithm
13.6, under the lock $M$, after the test of line 6, which reads the
value under the same lock. That value is the last one stored, since all
the stores are made under $M$, which orders them; the store of the node
and the call of the callback follow under $M$. The test of line 4 only
returns early, when $v$ is not below a value that $z$ had at some
moment, which is not below the current one, since each store of line 7
lowers $z$ (line 6). This gives (i), (ii) and (iii). A worker that read
the value $z' \ge z$ earlier (the read without the lock may return an
older value) prunes $\beta \ge \theta(z')$, which implies $\beta \ge
\theta(z)$ by [Lemma 13.3](#lem-13-3). $\square$

In the library, when the type of the values has a lock-free
`std::atomic`, the value is an atomic variable on a line of cache of its
own, which the workers read without a lock and with relaxed order when
they prune. The first test of `improve()` discards without the lock the
values that are not better than the one seen, and only an improvement,
which is usually "rare", takes the lock. For another type (e.g., the
rationals of an exact arithmetic) the value is read under the lock too.

<a id="thm-13-13"></a>**Theorem 13.13** (Correctness of a parallel
search). Under (A1) to (A4) and with $\tau_r \le 1$, let a parallel
search process each node by PROCESS with an incumbent value read from
`BBIncumbent` (or, over MPI, received from it,
[§13.9](#sec-13-9)), and let it move nodes only between the places of
the wider open set. Then at any moment every $x \in X(R)$ with
$v(x) < \theta(z)$, $z$ the current incumbent, belongs to $X(n)$ for a
node $n$ of the wider open set. Hence, if the search ends only when the
wider open set is empty, or when each node in it has a bound of at
least $\theta(z)$, its value satisfies
[Theorem 13.5](#thm-13-5), and the bound computed from the leaves of
all the workers satisfies [Theorem 13.6](#thm-13-6).

*Proof.* We note that the proof of [Lemma 13.4](#lem-13-4) applies to
each step of each worker, the steps of different workers being
interleaved. A step that moves a node keeps it in the wider open set. A
step of PROCESS replaces the node by its children or makes it a leaf,
and a leaf pruned by a value $z'$ read earlier has no element below
$\theta(z') \ge \theta(z)$ ([Lemma 13.12](#lem-13-12)). An improvement of
$z$ shrinks the set to be covered. Finally, the leaves collected by the
workers cover $X(R)$ together with the wider open set, which gives the
bound as in [Theorem 13.6](#thm-13-6) once the nodes left open are added
as leaves (each search does so at its end). $\square$

All the schemes below satisfy the hypotheses, but for the SelfSplit,
whose workers prune by incumbents of their own and which
[Theorem 13.23](#thm-13-23) treats apart. What remains to prove for
each is that it ends only in one of the two situations of [Theorem
13.13](#thm-13-13), and that it does end. Of the two, the end is the
delicate part: a worker that finds no node cannot conclude that none
exists, since another worker may be about to generate some (cf. the
general discussion of the termination of parallel searches in [Chapter
8](08-shared-memory.md#ch-8) and [Chapter
9](09-hash-distributed.md#ch-9)). Also, one can stop each scheme at a
node or time limit (LIMIT, lines 16 and 17 of Algorithm 13.6,
`detail::bb_limit()`,
[`ParallelBranchAndBound.h:321`][ParallelBranchAndBound.h:321]), and the node
that a worker holds when it sees the limit becomes an open leaf. In every
search the count $c$ is that of the evaluations, i.e., of the calls of PROCESS
that reach line 18 (a node pruned at line 17 is not counted); in the shared
queue, the MultiQueue, the farm, the distributed queues and the work stealing
it is a counter shared by the workers, raised after each such call. The limit
is then exceeded by a bounded amount.

<a id="prop-13-14"></a>**Proposition 13.14** (The node limit of the
asynchronous schemes). Let a node limit $L > 0$ be given, and no time
limit. In the shared queue, the MultiQueue, the distributed queues and
the work stealing, the evaluations are at most $L + p - 1$; in the
farm, with batches of $k$ nodes, at most $L + p k$.

*Proof.* In the first four schemes a worker calls LIMIT on the counter
$c$ before it processes each node, the node of a dive included, and
raises $c$ after an evaluation. Consider the last evaluation that
starts: at its test $c \le L - 1$, i.e., at most $L - 1$ evaluations
had ended, and at most $p - 1$ others were under way, one per other
worker; hence the evaluations are at most $(L - 1) + (p - 1) + 1$. In
the farm the master tests the limit when a task comes back and
dispatches no task after the test has succeeded, while a worker tests
it only after each evaluation, to stop the dive of the node, and
completes the at most $k$ nodes of its task. The evaluations that end
after $c$ has reached $L$ belong to the tasks out at that moment or sent
by the dispatch under way, which are at most $p$ together, since the
master never has more than $p$ tasks out; each of them evaluates at most
its $k$ nodes, since no dive goes on once $c \ge L$. $\square$

The bound counts the evaluations; the nodes pruned when extracted, which
cost no evaluation, are not counted. The rounds, instead, never exceed
the limit (Theorem 13.18), and over MPI each process stops at its share
([§13.9](#sec-13-9)).

### Anomalies

Of course, the evaluations of a parallel search are not those of the
sequential one, and they may be more or fewer (often, but not always,
more). On the one hand, a worker may evaluate a node that the
sequential search prunes by an incumbent which, in the parallel search,
another worker finds only later (a *detrimental* anomaly, which is the
search overhead of a B&B). On the other hand, a worker may find an
incumbent earlier than the sequential search does and prune nodes that
the latter evaluates (an *acceleration* anomaly, which gives a speedup
larger than $p$). Lai and Sahni [[149]](Z-bibliography.md#ref-149) showed
that both occur in parallel best-first B&B, and that adding workers may
even "slow" a search down; Gendron and Crainic
[[16]](Z-bibliography.md#ref-16) discuss them with the literature that
followed. With the best-first selection the detrimental part is
somewhat limited, since the nodes extracted by the workers are among
the best open ones. It grows with the distance between the order of
extraction and the sequential one, and this distance distinguishes the
schemes below.

### The dive in a parallel search

With depth-first and plunging, a worker of a parallel search goes on
with one child of the node it processed without putting it in a queue
(DIVE, lines 9 to 15 of Algorithm 13.6, `detail::bb_dive()`,
[`ParallelBranchAndBound.h:296`][ParallelBranchAndBound.h:296]). With
depth-first it is the last child, which the queue would extract next anyway.
With plunging it is the last child of least bound, if that bound is "close
enough" to the least open bound bl that the worker knows: the least bound of
the shared queue (read under its lock), that of the queue of the worker in the
distributed queues, over MPI and in the work stealing, the least open
bound of the moment at which the task was dispatched in the farm, and
the bound of the node processed in the MultiQueue. We write DIVE($K$)
for that child, nil if there is none; it stays in the wider open set as
a node held by the worker.

## 13.8 One pool of open nodes {#sec-13-8}

### The shared queue

In the shared queue (`SharedBranchAndBound`,
[`ParallelBranchAndBound.h:350`][ParallelBranchAndBound.h:350]) the workers
extract the nodes from one queue under a lock $L$, evaluate and branch them
outside of it, and push the children back under it, as the workers of SPA$^*$
do with the states of A$^*$ ([Chapter 8](08-shared-memory.md#ch-8)); the child
of a dive stays with the worker. Its order of extraction is close to the
sequential one, from which it departs by the $p$ nodes processed at
once and by the times at which the incumbents are found. Also, the lock
is rather "cheap" when the evaluation dominates, which is the usual case
of a B&B. The pseudocode, with a counter `busy` of the workers that hold
a node, is the following.

**Algorithm 13.7** (The shared queue).

```text
 1  procedure SHARED(R)       ▷ one queue Q, under the lock L, for all
 2     Q ← R;  busy ← 0;  done ← false
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        d ← nil
 5        loop
 6           if d ≠ nil then n ← d;  d ← nil
 7           else n ← NEXT();  if n = nil then return
 8           if a limit is reached then
 9              n is a leaf (open);  busy ← busy − 1;  done ← true;  return
10           K ← PROCESS(n)
11           lock L;  d ← DIVE(K);  push the nodes of K but d in Q
12           if d = nil then busy ← busy − 1
13           unlock L
14  procedure NEXT()
15     loop
16        if done then return nil
17        if Q = ∅ then                             ▷ read without the lock
18           if busy = 0 then
19              lock L;  if Q = ∅ and busy = 0 then done ← true;  unlock L
20           else wait a little
21           continue
22        lock L
23        if PRUNABLE(the least bound of Q) then
24           the nodes of Q are leaves (pruned);  Q ← ∅;  unlock L
25           continue
26        busy ← busy + 1;  n ← the first node of Q, removed;  unlock L
27        return n
```

[Figure 13.6](#fig-13-6) shows two steps: a worker extracts a node,
raising `busy`, while another pushes the children of the node it
processed, lowering it.

![Figure 13.6. The shared queue before and after two steps of
Algorithm 13.7 (minimization, each node with its bound): worker 1
extracts the first node a under the lock (line 26, busy $+ 1$), and
worker 0, having processed x, pushes its children x$_1$ and x$_2$
under the lock and, having no dive, leaves the count (lines 11 and 12,
busy $- 1$).](../figures/branch-and-bound-shared.svg){#fig-13-6}

<a id="prop-13-15"></a>**Proposition 13.15** (The end of the shared
queue). In Algorithm 13.7, (i) whenever $L$ is free, `busy` is the
number of the workers that hold a node (extracted and not yet pushed
back, or kept for a dive); (ii) if a worker finds $Q = \emptyset$ and
`busy` $= 0$ under $L$ (line 19), no node is open in the wider sense;
(iii) if the evaluations terminate, each worker returns.

*Proof.* (i) We observe that a worker raises `busy` under $L$ in the
same critical section in which it takes a node from Q (line 26). It
keeps it raised while it processes the node and the nodes of its dive
(lines 6 and 10), and lowers it under $L$ only in the critical section
in which it has pushed the children and kept no dive (lines 11 to 13),
or when it stops at a limit (line 9). (ii) The wider open set of this
scheme is Q and the nodes held by the workers (no node travels), and
both are empty. (iii) A worker loops in NEXT only while Q is empty and
`busy` $> 0$, or while Q is not empty. In the first case some worker
holds a node, which it processes in finite time, and whose dive ends at
a leaf of the finite tree; that worker then pushes the children or
lowers `busy`. In the second case the worker takes a node or empties Q
at line 24. Since each node is processed once and the tree is finite
(A4), the loop ends by line 16 or 19. $\square$

The test of line 17 reads without the lock an atomic copy of the size
of Q. Hence, a worker that finds no node waits without taking
the lock, and with a growing pause, while the decision to stop is taken
under the lock, which makes the stale reads harmless. The counters
`busy`, the size of Q, the evaluations and the flag `done` sit on
separate lines of cache. Also, when the least bound of Q cannot improve
the incumbent, the whole queue becomes a set of leaves at once
(line 24). The push of the children (lines 11 to 13) holds the lock for
the time of the pushes of one node, i.e., $O(|K| \log N)$ by
Proposition 13.9, and the extraction (lines 22 to 27) for one POP; each
node processed thus costs two critical sections, besides those of
line 19, and the lock of the incumbent only when it improves it.

### The MultiQueue

The MultiQueue B&B (`MultiQueueBranchAndBound`,
[`MultiQueueBranchAndBound.h:129`][MultiQueueBranchAndBound.h:129]) is the
shared queue with the relaxed queue of MQA$^*$ ([Chapter
8](08-shared-memory.md#ch-8)). The open nodes are in a MultiQueue of Rihani,
Sanders and Dementiev [[37]](Z-bibliography.md#ref-37), i.e., $2p$ sequential
queues, each under a separate lock. A pop reads the least keys of two queues
drawn at random, without locks, and takes the least node of the better one; no
lock is common to all the workers. The price is a relaxed order of
extraction: indeed, a pop may return a node that is not the least of
all. In the sequential model of the MultiQueue (the nodes spread at
random, then removed by the better of two queues drawn at random) the
expected number of open nodes with a smaller key, the *rank error*,
converges, as the pops go on, to a value below
$\frac{5}{6} \cdot 2p$, i.e., it is $O(p)$
(cf. [[37]](Z-bibliography.md#ref-37) and
[[38]](Z-bibliography.md#ref-38)). This is a result of the sequential
model with two new queues at each pop (the library draws two new queues
at each pop, i.e., with a stickiness of 1), which we do not prove here;
it does not cover the buffers of the queues nor the push of all the
children of a node into one queue, and we claim no bound on the rank
error of the library. Its termination test cannot use a
counter of busy workers, since there is no lock under which Q and the
counter are read together. It uses instead a counter $m$ of the nodes
that are in the queue or held by a worker.

**Algorithm 13.8** (The MultiQueue).

```text
 1  procedure MULTIQUEUE(R)
 2     Q ← R;  m ← |R|;  done ← false
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        d ← nil
 5        while not done do
 6           if d ≠ nil then n ← d;  d ← nil
 7           else
 8              n ← a pop of Q, possibly nothing
 9              if no node was popped then
10                 if m = 0 then done ← true else wait a little
11                 continue
12              if PRUNABLE(the bound of n) then
13                 n is a leaf (pruned);  m ← m − 1;  continue
14           if a limit is reached then n is a leaf (open); done ← true; return
15           K ← PROCESS(n);  d ← DIVE(K)
16           m ← m + |K| − 1;  push the nodes of K but d in one queue of Q
```

![Figure 13.7. The MultiQueue with $p = 2$ workers and $2p = 4$
queues (minimization). (a) Worker 0 pops: it reads the least keys of
Q$_1$ and Q$_2$, drawn at random, and takes c, of bound 3, from the
better one, although a, of bound 2, is in Q$_0$; the counter $m$ counts
the six queued nodes and the node y held by worker 1. (b) Worker 1 has
processed y: it adds its two children and drops y in one atomic step,
$m \leftarrow m + 2 - 1$, before it pushes both children into one queue
under one lock.](../figures/branch-and-bound-multiqueue.svg){#fig-13-7}

<a id="prop-13-16"></a>**Proposition 13.16** (The end of the
MultiQueue). In Algorithm 13.8, $m$ is at any moment at least the
number of the nodes in Q plus those held by the workers; hence,
$m = 0$ at line 10 means that no node is open, and the search then
ends. If the evaluations terminate, each worker returns.

*Proof.* We first note that at the start $m = |R|$, and the claim holds
with equality. A pop moves a node from Q to a worker and leaves $m$
unchanged, and a node pruned at line 13 leaves both and is subtracted.
At line 16 the node processed leaves and $|K|$ children appear, of which
DIVE($K$) is held by the worker and the others enter Q. Since $m$ is
changed by $|K| - 1$ in one atomic step *before* the pushes, between
this step and the end of the pushes $m$ counts children not yet in Q,
which is an overcount and never an undercount. A worker that stops at a
limit sets `done` for all. For the second claim, note that a pop finds
no node only after a scan of all the queues found them empty. Hence, a
worker waits at line 10 only while $m > 0$ and each counted node is held
by a worker or being pushed; those nodes are processed in finite time,
and the tree is finite. $\square$

The least bound of Q is not known, since each queue knows only its least
key. Hence, a node is pruned when popped (line 12) rather than with all
the others, and DIVE compares the children with the bound of the node
processed. The children of a node go into one queue under one lock (the
operation `pushAll` of the MultiQueue), which relaxes the rank somewhat,
as a stickiness of their number would, and saves a lock per child.
Furthermore, the values must be atomic with no lock (the concept
`MultiQueueValue`,
[`MultiQueueBranchAndBound.h:105`][MultiQueueBranchAndBound.h:105]), since the
MultiQueue keeps the least key of each queue in an atomic word that the pops
read without the lock. Similarly, the counter of the evaluations, which each
worker would update, is touched only when a node limit is set. In the code $m$
is the counter `held`, and line 16 is one atomic addition of $|K| - 1$ (a
subtraction of 1 when $K = \emptyset$, in the modular arithmetic of the
unsigned counter), made before the operation `pushAll` that moves the nodes of
$K$ but $d$ into one queue under its lock; the case of one child, which
replaces the node processed, needs no atomic operation. We expect the
MultiQueue to pay off only when an evaluation is so "fast" that the lock of the
shared queue is the bottleneck; when the evaluation dominates, its relaxed
order perhaps costs more evaluations than the lock costs time.

### The master-worker farm

The farm (`FarmBranchAndBound`,
[`FarmBranchAndBound.h:104`][FarmBranchAndBound.h:104]) is the master-worker
form of the shared queue, written as a farm with feedback, the pattern of
FastFlow [[12]](Z-bibliography.md#ref-12), on the executor of the search (the
tests run it on the three executors: the threads of C++, OpenMP and FastFlow).
The master holds Q, extracts the nodes by the selection and sends them
in tasks of `batch` nodes to the free workers. A worker evaluates and
branches each node of its task (going on with the child of a dive in
the same task), and sends the task back with the children that survive
the pruning, which the master pushes in Q.

**Algorithm 13.9** (The master-worker farm).

```text
 1  procedure MASTER(R, p)
 2     Q ← R;  out ← 0
 3     loop
 4        DISPATCH()
 5        if out = 0 and (a limit is reached or Q = ∅) then return
 6        receive a task T back from a worker;  out ← out − 1
 7        push the children of T in Q
 8  procedure DISPATCH()
 9     while no limit is reached, out < p and Q ≠ ∅ do
10        if PRUNABLE(the least bound of Q) then
11           the nodes of Q are leaves (pruned);  Q ← ∅;  return
12        T ← the first batch nodes of Q, removed
13        send T to a free worker;  out ← out + 1
14  procedure WORKER(T)                       ▷ on each of the p workers
15     for each n ∈ T do
16        while n ≠ nil do                  ▷ the node, then its dive
17           K ← PROCESS(n);  n ← DIVE(K), unless the node limit is reached
18           add the nodes of K but n to the children of T
19     send T back to the master
```

![Figure 13.8. The farm with batches of two nodes. (a) The master
holds Q and one task, with the nodes a and b, is out at worker 0.
(b) The master has sent the task of c and d to the free worker 1, and
worker 0 has sent its task back with the children a$_1$ and b$_1$,
which the master has pushed in Q; out is
$1 + 1 - 1 = 1$.](../figures/branch-and-bound-farm.svg){#fig-13-8}

<a id="prop-13-17"></a>**Proposition 13.17** (The end of the farm).
In Algorithm 13.9 each open node is either in Q or in a task that is
out at a worker, and `out` is the number of such tasks. Hence, the test
at line 5 is exact: the master returns only when each task is back and
either no node is open, or a limit is reached, or (through line 11) no
open node can improve the incumbent.

*Proof.* We observe that only the master changes Q and `out`. It moves
nodes from Q to a task and raises `out` (lines 12 and 13), and it lowers
`out` when a task comes back, pushing its children (lines 6 and 7). A
worker changes only the nodes of its task, replacing each node by the
children that survive PROCESS and keeping the dive within the task. When
`out` = 0 no task is out, and Q holds all the open nodes. $\square$

The master also never waits for nothing: if the evaluations terminate,
the farm ends. Indeed, a task out comes back in finite time, since its
worker processes at most `batch` nodes and their dives, each of which
ends at a leaf by (A4); and after each DISPATCH either `out` $> 0$, or Q
is empty, or a limit is reached, since DISPATCH sends a task whenever
`out` $< p$, Q is not empty and no limit is reached, unless it empties
Q at line 11. Hence the receive of line 6 always has a task to wait
for, each node is processed once, and the tree is finite.

In the farm the order of extraction is that of the selection, and the
master only moves nodes; hence, the farm suits an evaluation that is
"expensive" w.r.t. the transfer of a node. The batch of a task trades
the messages for evaluations. The nodes of a task are extracted
together, and the children of the first of them cannot overtake the
others; a batch of 1 is therefore the one for an expensive evaluation,
while the default of the parallel options, 16 nodes, amortizes the
messages when the evaluation is "cheap". Also, the tasks are recycled
through a free list of the master, and the farm runs $p + 1$ threads,
the master on the first CPU of the affinity mask when the workers are
pinned. The dispatch is lines 8 to 13 of Algorithm 13.9, in which each
task also carries the least open bound of the moment of its dispatch,
the bound bl of the DIVE of its worker.

### The deterministic rounds

The rounds (`DeterministicBranchAndBound`,
[`DeterministicBranchAndBound.h:125`][DeterministicBranchAndBound.h:125]) make
the outcome independent of the timing of the workers, as the deterministic mode
of FiberSCIP [[18]](Z-bibliography.md#ref-18) does, so that a run can be
repeated exactly (e.g., for the experiments of a paper, or to find a bug). At
each round the $r$ best open nodes are extracted from one queue
($r$ = `round`, or $p$ if 0), and the worker $i$ processes the nodes
$i, i + p, i + 2p, \ldots$ of the round against the incumbent $z_0$ of
the start of the round. When all of them are done, one worker merges
the outcomes in the order of the nodes of the round and extracts the
next round, while the others wait at a barrier.

**Algorithm 13.10** (The deterministic rounds).

```text
 1  procedure DETERMINISTIC(R)
 2     Q ← R;  z ← the given incumbent, or ∞;  N ← ()
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        loop
 5           barrier, at which one worker runs NEXTROUND()
 6           if done then return
 7           for j = i, i + p, i + 2p, ... while j < |N| do
 8              (K_j, v_j) ← PROCESS(N_j) against z₀, v_j being the best
                   feasible value of N_j below z₀, if any
 9  procedure NEXTROUND()
10     for j = 0, ..., |N| − 1 do                  ▷ in the order of N
11        if v_j < z then z ← v_j;  N_j is the node of z
12        for each m ∈ K_j do
13           if PRUNABLE(β(m)) then m is a leaf (pruned) else push m in Q
14     if a limit is reached, Q = ∅ or PRUNABLE(the least bound of Q) then
15        done ← true;  return
16     N ← the first min(r, evaluations left) nodes of Q, removed
17     z₀ ← z
```

![Figure 13.9. A round with $r = 4$ nodes and $p = 2$ workers
(minimization). (a) The four best nodes of Q form the round, the node
$j$ going to worker $j \bmod 2$, and each node is processed against
$z_0 = 10$. (b) After the barrier the outcomes are merged in the order
of $j$: the feasible value 9 of N$_0$ becomes the incumbent and its
child k$_1$ is pushed; the child k$_2$ of N$_2$ is pushed, since
$8 < 9$; the value 8 of N$_3$ becomes the incumbent and prunes the
child k$_3$, of bound 8. The order of the merge, not the order in which
the workers finished, decides all of
this.](../figures/branch-and-bound-deterministic.svg){#fig-13-9}

<a id="thm-13-18"></a>**Theorem 13.18** (Determinism of the rounds).
Assume that no time limit is set, and that the evaluator of each worker
returns, for a node and an incumbent given by `cutoff`, an outcome that
depends only on them. Then the nodes evaluated, the incumbent value, its
node, the bound and the counters other than the evaluations of each
worker and the time are functions of the problem, of $R$, of the options and of
$r$, and depend neither on the timing nor on the number $p$ of the
workers. If the outcome of the evaluator of worker $i$ also depends on
the nodes that $i$ evaluated before, the same holds with $p$ among the
arguments. The node limit is exact.

*Proof.* We show by induction on the rounds that the state at the
barrier (the content of Q with the numbers of their pushes, $z$ and its
node, the leaves and the counters) is such a function. At the start it
is fixed by $R$ and the options. Given the state, NEXTROUND extracts
$N$ by the order of the queue, which depends only on the keys and the
numbers of the pushes, and sets $z_0$. Each outcome $(K_j, v_j)$ is
computed by PROCESS on $N_j$ with the incumbent $z_0$: the function
`inc()` that the workers pass returns the snapshot of the round, and
`improve` records the best value of the node without touching $z$.
Hence, the outcome is a function of $N_j$ and $z_0$ under the first
assumption. Under the second one it is a function of $N_j$, $z_0$, $p$
and the earlier rounds, since the nodes of worker $i$ are those of
index $i \bmod p$. The merge then runs on one thread in the order of
$j$. Finally, the last round is shortened to the evaluations left
(line 16), and a round never stops in the middle; thus the node limit
is met exactly. $\square$

We remark that the history of a worker depends on $p$, as the proof
shows; hence, the independence of $p$ needs an evaluation that depends
on the node (and on $z_0$) only, as that of the knapsack does, and a
history-dependent evaluator gives only the weaker property. Indeed, the
test `bb-det-thread` checks that the status, the value, the node, the
bound and the counters of the evaluated, generated and pruned nodes are
the same with 1, 2, 3 and 8 workers and at each repetition, also when a
node limit stops the search. The test `bb-cutoff` checks that the
value, the bound and the numbers of the evaluated and of the infeasible
nodes are the same with workers that fence the nodes by the incumbent
given by `cutoff`. The time limit, checked between the rounds,
makes the end depend on the timing.

The price is the wait of each round for its slowest node, which is
usually "small" when the evaluation dominates and its time varies little
among the nodes, and an incumbent that may be one round old (which may
cost somewhat more evaluations). Also, no dive is taken, and the
plunging selection extracts the nodes by their bound. In the code the
barrier is a `std::barrier` whose completion function is NEXTROUND; the
merge thus runs on one of the workers (the standard does not say which;
in practice the last to arrive) while the others are blocked.
Furthermore, the slots of a round, each with its node, its
children, its best value, its counters and its leaves, are reused from
round to round. An exception thrown by an evaluation is kept, and thrown
again by the search after the workers have stopped. The work of a worker
in a round is line 8 of Algorithm 13.10, with the function `inc()` of
PROCESS returning $z_0$ and the function `improve` keeping in the slot
the best value below $z_0$; the merge is lines 10 to 13, in the order of
the slots, where the counters of each slot are added to those of the
search, so that they do not depend on which worker processed the slot.
It remains to show that the rounds are correct and end.

<a id="prop-13-19"></a>**Proposition 13.19** (The end of the rounds).
Under (A1) to (A4) and with $\tau_r \le 1$, (i) at each barrier, after
NEXTROUND, every $x \in X(R)$ with $v(x) < \theta(z)$ is in the subset
of a node of Q or of the new round $N$; (ii) when NEXTROUND sets `done`
with no limit reached, the value satisfies [Theorem 13.5](#thm-13-5)
and the bound [Theorem 13.6](#thm-13-6); (iii) if the evaluations
terminate, the search ends after finitely many rounds.

*Proof.* (i) We argue by induction on the rounds. The slots of a round
are the nodes removed from Q at line 16, and each of them is processed
exactly once, by the worker $j \bmod p$. PROCESS prunes by $z_0$, a past
value of $z$, which is safe by [Lemma 13.12](#lem-13-12), and the value
found by a slot is offered to $z$ at line 11, before the end test, so
that a node with no child is covered as in [Lemma 13.4](#lem-13-4).
After the barrier each slot has been processed, and line 13 puts each
child in Q or makes it a leaf by the current $z$; hence, at line 14, the
open nodes are those of Q, and the claim follows from the proof of
[Lemma 13.4](#lem-13-4) applied to each slot. (ii) At line 14 no node
is held by a worker and none is in a slot, since the end is decided by
one thread while the others wait at the barrier; hence `done` is set
only when Q is empty or no node of Q can improve $z$, which are the
situations of [Theorem 13.13](#thm-13-13). (iii) A round processes at
least one node, since $r \ge 1$ at line 16 (with a node limit,
NEXTROUND stops at line 14 as soon as the evaluations reach it), and
each node is processed once in a finite tree; a round ends when its at
most $\lceil r / p \rceil$ processings per worker end. $\square$

## 13.9 One pool per worker {#sec-13-9}

### The distributed queues

In the distributed queues (`DistributedBranchAndBound`,
[`ParallelBranchAndBound.h:599`][ParallelBranchAndBound.h:599]) each worker has
a queue $Q_i$ and sends each child, in batches, to the worker that owns it, as
HDA$^*$ does with the states of A$^*$ ([Chapter
9](09-hash-distributed.md#ch-9)). Its owner is $\mathrm{hash}(m) \bmod p$ if
the problem gives a hash of its nodes, and a random worker otherwise. Since a
tree has no duplicates, the owner only balances the load, and a uniform hash is
perhaps the most suitable distribution; the distributions of HDA$^*$ that keep
related states on the same worker (abstraction, abstract Zobrist hashing) have
little to gain here. In the Doxygen comment of
[`ParallelBranchAndBound.h`][ParallelBranchAndBound.h] the pseudocode follows
that of the shared queue, as its lines 28 to 45; here its lines are numbered
from 1.

**Algorithm 13.11** (The distributed queues).

```text
 1  procedure DISTRIBUTED(R)  ▷ a queue Q_i per worker, nodes to the owners
 2     put each node of R in the queue of its owner
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        d ← nil
 5        while the search is not over do
 6           push the nodes received in Q_i
 7           if d ≠ nil then n ← d;  d ← nil
 8           else if Q_i ≠ ∅ and not PRUNABLE(the least bound of Q_i) then
 9              n ← the first node of Q_i, removed
10           else
11              the nodes of Q_i are leaves (pruned);  Q_i ← ∅
12              send all the outboxes;  wait for mail or the end; continue
13           if a limit is reached then n is a leaf (open); stop all; return
14           K ← PROCESS(n);  d ← DIVE(K)
15           for each m ∈ K but d do
16              j ← hash(m) mod p, or a random worker without hash(m)
17              if j = i then push m in Q_i else send m to j
18           every k evaluations, send the outboxes not yet full
```

![Figure 13.10. The distributed queues with three workers, the
color of a node being that of its owner. (a) Worker 0 extracts n and
branches it; its children are owned by the workers 0, 1 and 2.
(b) Worker 0 keeps n$_1$ in its queue (first, since the ties go to
the last node pushed) and sends n$_2$ and n$_3$, which are in flight
and counted as such until their owners receive
them.](../figures/branch-and-bound-distributed.svg){#fig-13-10}

<a id="prop-13-20"></a>**Proposition 13.20** (The end of the
distributed queues). Algorithm 13.11 ends only when no worker holds a
node and no node is in flight, with each queue empty (a queue whose
least bound cannot improve the incumbent having been emptied into
leaves at line 11); if the evaluations terminate, it ends.

*Proof.* As we said, the wider open set is made of the queues, the dives
and the nodes in the outboxes and in flight. A worker reaches the test
of the end (line 12) only with no dive and an empty queue, after sending
all its outboxes. The test is that of the mailbox (`quiescent()`), the
deficit counting of [Definition 9.6](09-hash-distributed.md#def-9-6),
whose hypotheses hold here: a worker sends only while it has work (at
line 17, or when the mailbox flushes its outboxes before it becomes
passive) and becomes active again only by draining a message (line 6).
By [Lemma 9.7](09-hash-distributed.md#lem-9-7) the mailbox reports the
end only when no worker is active and no message exists, i.e., when
each worker is at line 12 with no dive and an empty queue and each
node sent has been received; at that moment the wider open set is
empty. Conversely, while a node exists in the wider set, some worker
holds it or will receive it, and processes it in finite time; since
each node is processed once and the tree is finite, the end comes, and
by part (iii) of the lemma all the workers then read it. After a stop
at a limit the mailbox reports the end to all the workers at once, and the
nodes left in the queues, in the dives, in the outboxes and in the
channels are collected as open leaves once the workers have returned
(`leftovers()`), so that the bound covers them. $\square$

The batches are of 16 nodes by default, and the outboxes are flushed at
each evaluation (`flush_period` = 1), since an evaluation is usually
"expensive" and a node that waits in an outbox delays its owner. The
nodes a search starts from go to their owners, or to the worker $i \bmod
p$ for the $i$-th of them when the problem gives no hash. Without a
hash, the owners come from a generator per worker, seeded by `seed` plus
the index of the worker. Note that, since a worker extracts the best
node of its queue only, it may evaluate nodes "far" worse than the
best open one before an incumbent prunes them, and how many depends on
when the incumbent is found, i.e., on the timing. Plunging, whose dives
reach a leaf early, tends to keep them close to those of the sequential
search. The choice of the owner is lines 15 to 17 of Algorithm 13.11,
the random owner being drawn by a generator of 64 bits of each worker.

### The distributed queues over MPI

Over MPI (`MPIBranchAndBoundSearch`,
[`MPIBranchAndBound.h:225`][MPIBranchAndBound.h:225]) the same scheme runs with
one process per worker, and the nodes travel as bytes, written by `serialize()`
and read by `deserialize()` (a trivially copyable node needs neither). Each
process has its queue and a copy of the incumbent value. A process that
improves it sends the value to all the others, while the node of the incumbent
stays with the process that found it until the end, when the process of the
least value sends it to all.

**Algorithm 13.12** (The distributed queues over MPI).

```text
 1  procedure MPI-B&B(R_i) on the process i       ▷ R_i: its own nodes
 2     Q_i ← R_i;  d ← nil
 3     loop
 4        for each message m received do
 5           if m is a batch of nodes then push them in Q_i
 6           else if m is an incumbent z' then z ← min(z, z')
 7           else if m is a stop then stop evaluating
 8        if not stopped and (d ≠ nil or not PRUNABLE(βmin(Q_i))) then
 9           n ← d, or the first node of Q_i, removed;  d ← nil
10           K ← PROCESS(n), a value v < z sending v to all the others
11           d ← DIVE(K)
12           for each m ∈ K but d do
13              j ← hash(m) mod p, or a random process without hash(m)
14              if j = i then push m in Q_i else send m to j
15           every k evaluations, send the outboxes not yet full
16           if the share of the node limit or the time limit is reached
                then stop all the processes
17        else send all the outboxes
18        if the mailbox detects the end then break
19     the nodes still travelling, after a stop, are put in Q_i
20     the process of the least z sends its node to all
```

![Figure 13.11. The distributed queues over MPI with three
processes (minimization, the color of a node being that of its owner).
(a) Process 1 finds the value 9 and sends it to the other two, while a
batch of nodes travels from process 0 to process 2; the processes 0
and 2 still prune by $z = 12$, which is safe. (b) After the receipts
each process prunes by $z = 9$, and the nodes whose bound is not less
than 9 are pruned when they reach the top of their queue, or with
their whole queue.](../figures/branch-and-bound-mpi.svg){#fig-13-11}

<a id="prop-13-21"></a>**Proposition 13.21** (Correctness and end over
MPI). The value and the bound returned by Algorithm 13.12 satisfy
[Theorems 13.5](#thm-13-5) and [13.6](#thm-13-6) when no limit stops
it, and the bound satisfies [Theorem 13.6](#thm-13-6) in any case; the
search ends when no process has a node that can improve its incumbent
and no message is in flight.

*Proof.* We first note that the value $z$ of a process is at any moment
at least the least value found by all, which plays the part of the
shared incumbent of [Theorem 13.13](#thm-13-13). Indeed, it is the
initial value or a value received, each value sent is feasible, and a
process that prunes by a (possibly larger) value prunes less, as in
[Lemma 13.12](#lem-13-12). The end is detected by the four-counter method
of Mattern [[11]](Z-bibliography.md#ref-11) in the mailbox of MPI
(proved sound in [Theorem 9.18](09-hash-distributed.md#thm-9-18)), which a
process joins (line 18) only when it has no node to evaluate and has
sent all its outboxes. The method counts each message, the incumbents
included, and reports the end only when all the processes are passive
and each message sent has been received. At that moment each process has
turned its queue into leaves by its incumbent (line 8), which is at
least the least value found, and no node travels; the wider open set is
then empty, and [Theorem 13.13](#thm-13-13) applies. After a stop at a
limit, the nodes still in flight are received as open leaves (line 19),
and the bound covers them. $\square$

The node limit $L$ is split among the processes, each stopping at
$\lceil L / p \rceil$ evaluations, and a process that reaches its share,
or the time limit, stops the evaluations of all. Each process returns
the same value, node and bound (an all-gather of the values picks the
process of the least value, the least rank on ties, which broadcasts
its node). Each keeps, however, the leaves
it owns, and the reoptimization of [§13.11](#sec-13-11) then starts each
process from them, with no transfer. All sends are non-blocking. MPI
must be initialized by the caller, with at least `MPI_THREAD_FUNNELED`
if the problem uses threads. The broadcast of an improvement is line 10
of Algorithm 13.12: a value $v$ below the local $z$ becomes the local
incumbent, with its node, and is sent, as its bytes, to each other
process with a tag of its own.

### The work stealing

In the work stealing (`StealingBranchAndBound`,
[`StealingBranchAndBound.h:117`][StealingBranchAndBound.h:117]), after Blumofe
and Leiserson [[20]](Z-bibliography.md#ref-20), each worker pushes the children
it generates in its queue and extracts the best node of that queue. A
worker whose queue is empty steals about half of the queue of another
one, drawn at random, the half that leaves the victim its best node (cf.
the steal-half queues of Hendler and Shavit
[[36]](Z-bibliography.md#ref-36)). A worker then typically stays in its
part of the tree. This pays off when moving its evaluator to a distant
node costs (e.g., when the evaluator applies and undoes the branchings
along the path from the root), while the nodes of the distributed queues
land on any worker.

**Algorithm 13.13** (The work stealing).

```text
 1  procedure STEALING(R)
 2     the k-th node of R in Q_(k mod p)
 3     in parallel, on each worker i = 0, ..., p − 1 do
 4        d ← nil
 5        loop
 6           if d ≠ nil then n ← d;  d ← nil
 7           else n ← NEXT(i);  if n = nil then return
 8           if a limit is reached then
 9              n is a leaf (open);  stop all the workers;  return
10           K ← PROCESS(n)
11           d ← DIVE(K);  push the nodes of K but d in Q_i
12  procedure NEXT(i)
13     if Q_i ≠ ∅ and PRUNABLE(the least bound of Q_i) then
14        the nodes of Q_i are leaves (pruned);  Q_i ← ∅
15     if Q_i ≠ ∅ then return the first node of Q_i, removed
16     loop
17        if stopped, or no worker holds a node or works on one, return nil
18        v ← a worker ≠ i at random
19        if Q_v holds two nodes at least then
20           move about half of them, but the first, to Q_i
21           if some were moved, return the first node of Q_i, removed
22        wait a little
```

The queues are the stealable priority queues of the default data
structures of the library, whose owner extracts its nodes in the exact
order of the selection and with no lock. Each queue keeps its nodes in
two 4-ary heaps, $H_0$ and $H_1$, which no other thread reads. A push
goes into the smaller of the two, and hence the nodes pushed one after
the other (e.g., the children of a node) are dealt to the two heaps in
turn; a pop takes the better of their two least nodes. A thief that
finds no node to claim raises a flag of the queue of the victim. The
owner reads it at its next operation, and moves the heap whose least
node is not the least of the queue into the *offer* of the queue, an
atomic pointer, which a thief claims by an exchange with null ([Figure
13.12](#fig-13-12)). Note that the owner keeps its least node. Until a
thief claims the offer, its nodes still belong to the owner: a pop whose
least local key is not less than the least key of the offer takes the
offer back by the same exchange, and whoever of the two exchanges first
gets it all. Hence, the common case of the owner, a push or a pop with
no thief waiting, costs no read-modify-write and no fence. A steal costs
the thief two exchanges at most, and the owner the store of a pointer
(or the time to deal its nodes again to two heaps of equal size, when
they differ in size by more than a factor of 4). The price is that a
thief waits for the next operation of the victim, i.e., at most the
processing of one node.

![Figure 13.12. A steal, with the keys of the nodes in the two
heaps of the victim. (a) Worker 0 has no node, leaves the counter $C$
of the workers that hold or work on a node, and raises the flag of the
queue of worker 1. (b) At its next operation worker 1 moves H$_1$,
whose least key 3 is not the least of its queue, into the offer, by a
pointer. (c) Worker 0 enters $C$ and claims the offer by an exchange;
worker 1 keeps its best node 2, and each worker pops its nodes with
no lock.](../figures/branch-and-bound-stealing.svg){#fig-13-12}

<a id="prop-13-22"></a>**Proposition 13.22** (The end of the work
stealing). Let $C$ be the number of the workers that hold a node or
work on one, plus the thieves that are trying to claim an offer. A
worker leaves $C$ when its pop finds its queue and its offer empty, and
a thief enters $C$ before the exchange that claims an offer and leaves
it if the exchange fails. Then $C = 0$ only when no node is open, and
the test at line 17, which returns nil when $C = 0$, is exact; if the
evaluations terminate, each worker returns.

*Proof.* We claim that a node is at any moment in a heap of its owner,
in an offer (which belongs to the owner until claimed), in the hands of
a thief between its exchange and its first pop, or held by a worker. In
the first two cases the owner is counted in $C$, since it left $C$ only
after finding its heaps and its offer empty, and it re-enters only by a
steal. In the third case the thief entered $C$ before its exchange, and
in the fourth the worker is counted. Hence a node exists only if $C >
0$. A worker loops at lines 16 to 22 only while $C > 0$, i.e., while
some worker holds or works on a node. That worker processes it in finite
time, and then pushes children (which a thief can claim after the next
operation of the owner) or drops it; the tree is finite. $\square$

The price of the scheme is a selection that is best-first within each
queue only, and hence an order of extraction that may be farther from
the sequential one than that of the shared queue. Also, the nodes have
to be copyable with moves that do not throw, since a steal moves them to
the thief, and the counter of the nodes sent in `BBStats` holds the
steals. The extraction of a worker (lines 12 to 15 of Algorithm 13.13)
drops its queue at once when its least bound, which counts the offer
not yet claimed, cannot improve the incumbent.

### The SelfSplit

The SelfSplit of Fischetti, Monaci and Salvagnin
[[19]](Z-bibliography.md#ref-19) (`SelfSplitBranchAndBound`,
[`SelfSplitBranchAndBound.h:130`][SelfSplitBranchAndBound.h:130]) gives a B&B
whose workers communicate only at the end. The tree is first explored
sequentially (the *sampling*), by the evaluator of the worker 0 and without
plunging, until the open nodes are $\sigma p$ ($\sigma$ being `split`, 4 by
default) or the tree is solved. Then, with `split_dive` (the default),
the same evaluator dives from a copy of the first open node to a leaf,
taking at each level the child of least estimate (of least bound, if
the problem gives no estimate). Finally the $k$-th open node, in the
order of the selection, goes to the worker $k \bmod p$. Each worker
searches its subtrees sequentially (by `BranchAndBoundSearch`) with its
evaluator, an incumbent starting from that of the sampling, and its
share of the node limit and of the time left; the incumbents, the
bounds and the leaves of the workers are merged in the order of the
workers.

**Algorithm 13.14** (The SelfSplit).

```text
 1  procedure SELFSPLIT(R)
 2     Q ← R;  z ← the given incumbent, or ∞     ▷ the sampling, by worker 0
 3     while p > 1 and Q ≠ ∅ and |Q| < σ p do
 4        if a limit is reached then break
 5        if PRUNABLE(the least bound of Q) then
 6           the nodes of Q are leaves (pruned);  Q ← ∅;  break
 7        n ← the first node of Q, removed;  Q ← Q ∪ PROCESS(n)
 8     if p > 1, split_dive and Q ≠ ∅ then         ▷ the dive, by worker 0
 9        n ← a copy of the first node of Q
10        while n ≠ nil and no limit is reached do
11           n ← the child of PROCESS(n) of least estimate, or nil
12     for k = 0, 1, ..., while Q ≠ ∅ do           ▷ the split
13        add the first node of Q, removed, to R_(k mod p)
14     in parallel, on each worker i = 0, ..., p − 1 do
15        S_i ← B&B(R_i) on the evaluator of i, from the incumbent z, with
             its share of the node limit and the time left, pruning also
             by the incumbent of all the workers if split_shared
16     z ← min(z, the values of the S_i), the node of the least i on ties
17     return z, its node, the leaves of the sampling and of the S_i, and
          the least of their bounds
```

![Figure 13.13. The SelfSplit with $p = 2$ and $\sigma = 2$
(minimization, the numbers being bounds). (a) The sampling by worker 0
evaluates r, s and s$_1$ and stops at four open nodes; the dive goes
from a copy of the first of them, s$_{11}$, to a leaf of value 7.
(b) The open nodes, in the order of the selection, go to the workers
$k \bmod 2$ (the colors), and each worker runs a sequential B&B on its
nodes from the incumbent $z = 7$, with no exchange until the
merge.](../figures/branch-and-bound-selfsplit.svg){#fig-13-13}

<a id="thm-13-23"></a>**Theorem 13.23** (The SelfSplit). Under (A1)
to (A4) and with $\tau_r \le 1$, (i) when no limit stops it, the value
returned by Algorithm 13.14 satisfies [Theorem 13.5](#thm-13-5), and its
bound satisfies [Theorem 13.6](#thm-13-6) in any case. (ii) Without
`split_shared` and without a time limit, its outcome depends only on
the problem, $R$, the options, $\sigma$ and $p$, provided that the
evaluation of a node by a worker depends only on the node, the
incumbent given by `cutoff` and the nodes evaluated by that worker
before. (iii) If the evaluations terminate, the search ends, its end
being the return of the last worker.

*Proof.* (i) We observe that the sampling is Algorithm 13.3 stopped
early. Hence, by [Lemma 13.4](#lem-13-4), each element of $X(R)$ below
the threshold of its incumbent $z_s$ is in an open node at its end, and
its leaves bound the rest. As for the dive, it works on a copy and
records no leaf; it may only lower $z_s$, by a value found. Then the
split partitions the open nodes among the $R_i$, and each worker is a
sequential B&B from $R_i$ with the initial value $z_s$. Let $z = \min_i
z_i$; the value $z_i \le z_s$ of the worker $i$ satisfies $\min \{ v(x)
: x \in X(R_i) \} \ge \theta(z_i) \ge \theta(z)$ by [Theorem
13.5](#thm-13-5) and [Lemma 13.3](#lem-13-3). With `split_shared` the
worker $i$ prunes by the least $c_i$ of $z_i$ and of the values found
by the others, and the same argument gives the bound
$\theta(c_i) \ge \theta(z)$, since $c_i \ge z$. The elements that the
sampling left below $\theta(z_s)$ are in some $X(R_i)$, and those of its
leaves are not below $\theta(z_s) \ge \theta(z)$; hence $z^*_R \ge
\theta(z)$. The bound is the least of the bounds of the sampling and of
the workers, each valid by [Theorem 13.6](#thm-13-6) for the region it
covers, and these regions cover $X(R)$. (ii) The sampling and the dive
are sequential, and the split depends only on the order of the queue and
on $p$. Worker $i$ receives $R_i$, $z_s$ and its share of the limit, all
fixed, and runs with no exchange, with its evaluator; the merge goes in the
order of the workers. (iii) The sampling stops when its queue holds
$\sigma p$ nodes, or earlier, and processes each node once in a finite
tree; the dive ends at a leaf by (A4); each worker runs Algorithm 13.3,
which ends by [Theorem 13.5](#thm-13-5); and the search waits for all
the workers to return before the merge, so that no node is open
anywhere when it merges. $\square$

The proof also shows why the dive matters. A best-first sampling seldom
reaches a leaf, and without an incumbent a worker whose subtrees do not
hold the optimum evaluates each node of them whose bound is better than
the best value in them, since no incumbent prunes them. Yet, the dive
gives each worker the incumbent of one leaf, at the cost of the "few"
nodes of the dive, which the worker whose subtree holds them evaluates
again (both evaluations are counted). With `split_shared` the workers
also prune by the best incumbent of all of them, which each reads at
each extraction through `shared_incumbent`. Even so, the value is still the
optimum, since the shared value is a value found (cf. [Lemma
13.12](#lem-13-12)), while the nodes evaluated, and the node of the
optimum among equal values, depend on the timing. With one worker there
is nothing to split: the root goes to it with neither the sampling nor
the dive, and the nodes evaluated are those of the sequential search.
The test `bb-split-thread` checks this, and that the outcome is the same
at each repetition, with a split of 1 and 4 nodes per worker and with
and without the dive.

In the library the sampling is done once, by the worker 0, rather than
by each worker as in the original scheme, since the workers share the
memory and can read its open nodes. Over MPI, where a node would
otherwise travel as bytes, each process would repeat it, as Fischetti,
Monaci and Salvagnin do. Its price is a load that is
balanced only as far as the subtrees of the split are, and (without
`split_shared`) an incumbent found by a worker that does not prune the
nodes of the others. The split is lines 12 and 13 of Algorithm 13.14:
with $L'$ the evaluations left by the sampling, the worker $i$ receives
the share $\lfloor L' / p \rfloor + 1$ if $i < L' \bmod p$ and
$\lfloor L' / p \rfloor$ otherwise, so that the shares sum to $L'$, and
a node whose worker has a share of 0 becomes an open leaf. The dive is
lines 8 to 11, taking the first of the children of least estimate;
its leaves go to a record of Algorithm 13.4 that keeps nothing, since
the open nodes of the sampling still cover them.

## 13.10 The schemes compared and chosen by name {#sec-13-10}

[Table 13.1](#sec-13-10) gathers, for each scheme, where the open nodes
are, how the end is detected and on what the outcome depends. Here
"the timing" means that two runs with the same data may evaluate
different nodes, and return a different node among those of equal
value. The value is in any case within the tolerance of the optimum in
each run that no limit stops.

**Table 13.1.** The B&B searches of the library: where the open nodes
are, how the end is detected (the condition of the proposition of each
scheme), and on what the outcome depends besides the problem and the
options; for the rounds, $p$ matters only with an evaluator that
depends on the history of its worker ([Theorem 13.18](#thm-13-18)), and
the SelfSplit depends on the timing with `split_shared`.

| Scheme | Open nodes | End | Depends on |
|:--|:--|:--|:--|
| sequential | one queue | Q empty, or line 6 | nothing |
| shared | one queue, a lock | Q empty, busy $= 0$ | the timing |
| multiqueue | $2p$ queues, a lock each | $m = 0$ | the timing |
| farm | the queue of the master | out $= 0$ | the timing |
| deterministic | one queue, by rounds | between rounds | $r$ (and $p$) |
| distributed | a queue per worker | message counts | the timing |
| MPI | a queue per process | Mattern's counters | the timing |
| stealing | a queue per worker | $C = 0$ | the timing |
| selfsplit | a queue per worker | at each worker | $\sigma$, $p$ |

We expect the choice among the schemes to depend on the cost of an
evaluation w.r.t. that of moving a node, and on the variance of that
cost among the nodes, which only a measurement on the instances of
interest can settle. The shared queue and the farm keep the order of the
selection, at the price of a lock or of a master, and the MultiQueue
removes the lock at the price of the order. The distributed queues and
the work stealing remove the common pool, the latter keeping each worker
in its part of the tree. Finally, the rounds and the SelfSplit perhaps
give up some speed for an outcome that does not depend on the timing.
All the schemes but the farm and the search over MPI, which build their
workers and their channels themselves, can be chosen at run time by a
value of the enumeration `BBScheme`, or by its name, read, e.g., from
a file of parameters. The enumeration is

```cpp
--8<-- "snippets/branch-and-bound-sig-scheme.cpp"
```

and the functions of [`BranchAndBoundScheme.h`][BranchAndBoundScheme.h] are
declared as follows (without their bodies; the template head of
`parallel_branch_and_bound()` has the problem type `P` and, as its
second parameter, the executor):

```cpp
[[nodiscard]] inline std::string_view bb_scheme_name( BBScheme s );
[[nodiscard]] inline BBScheme bb_scheme( std::string_view name );
[[nodiscard]] BBSolution< typename P::node_type , bb_value_t< P > >
parallel_branch_and_bound( const P & problem , BBScheme s ,
                           const BasicParallelBBOptions< bb_value_t< P > > &
                             o = {} ,
                           std::optional< std::vector< BBLeaf<
                             typename P::node_type , bb_value_t< P > > > >
                             seeds = std::nullopt );
```

Thus, a program built on the library gets each new scheme with no change
of its code. The following example, from the case `manual` of the
tester, which the build compiles and runs (`ctest -R PBeFS.manual`),
solves a knapsack `kp` of 50 items by the sequential search `s`, with
the options `o`, and by the scheme `sch` read from its name, with the
parallel options `po` and 4 workers, and shows that an unknown name is
refused; the case then checks both solutions against the dynamic
program.

```cpp
--8<-- "snippets/branch-and-bound-usage.cpp"
```

The names are `sequential`, `shared`, `distributed`,
`deterministic`, `selfsplit`, `stealing` and `multiqueue` (an unknown
name throws `std::invalid_argument`, and so does `multiqueue` with
values that are not atomic with no lock). The farm and the search over
MPI are not among the named schemes, and they are used by their
classes, `FarmBranchAndBound< P , Executor >( problem , options )` and
`MPIBranchAndBoundSearch< P >( problem , options , comm )`, whose
`solve()` returns the solution. The executor, i.e., the pool of C++
threads, OpenMP or FastFlow, is the second template parameter of each
parallel search on threads. Also, `set_executor()` gives a search an
executor that outlives it (without one, each search builds a pool).

## 13.11 Reoptimization {#sec-13-11}

When the data of the problem change in a way that leaves the branching
valid (e.g., a restriction of the feasible region, or a change of the
objective that the nodes do not store), the tree explored by a previous
search need not be built again. Its leaves cover the subtrees of the
nodes that the search started from, and a new search can start from
them. With `keep_frontier`, the `BBSolution` holds these leaves, each
with its bound (in the sense of the problem), its depth and its kind.
The overload `solve( seeds )` of each search (and
`branch_and_bound( p , seeds , options )`) starts from any set of such
nodes instead of the root, each of them being an open node that is
evaluated only when it is extracted. [Figure 13.14](#fig-13-14) shows
the simplest use: a search stopped by a node limit is resumed from its
open leaves, together with its incumbent.

![Figure 13.14. Resuming a search from its leaves. (a) The
best-first search of the example stopped by a node limit of 3: the
open leaves 0, 10 and 110, with their inherited bound 16.29, give the
bound 16.29 and the gap 4.29. (b) A search started from these leaves
with the incumbent 12 evaluates 110, 10, where it finds 15, and 0, then
prunes the child 1100 of 110 with the open nodes, and ends with the gap
0, as the search of [Figure 13.1](#fig-13-1)
did.](../figures/branch-and-bound-frontier.svg){#fig-13-14}

<a id="prop-13-24"></a>**Proposition 13.24** (Reoptimization). Let
$L$ be the leaves of a search from $R$ with `keep_frontier`. (i) Every
element of $X(R)$ belongs to $X(l)$ for some $l \in L$. (ii) If the
problem changes to $X'$ and $v'$, with (A2) to (A4) holding for it,
with the children generated by the previous search still covering the
new subset of their parent (A1), with the subsets of the nodes
$X'(n) \subseteq X(n)$, and with
$v' = v$, then a search from the non-infeasible leaves of $L$, with
their bounds, solves the new problem over $X'(R)$ in the sense of
[Theorems 13.5](#thm-13-5) and [13.6](#thm-13-6). (iii) If instead
$X'(n) \supseteq X(n)$, or $v' \ne v$, the same holds when each leaf,
the infeasible ones included, is given back with an unknown bound
($-\infty$ in the minimization), or with a bound that the caller knows
to be valid for the new problem.

*Proof.* Part (i) is the covering shown in the proof of
[Theorem 13.6](#thm-13-6). (ii) The bound of a leaf was a bound of
$X(l)$ for $v$, and hence of $X'(l) \subseteq X(l)$; an infeasible
leaf has $X'(l) \subseteq X(l) = \emptyset$. Hence, the seeds cover
$X'(R)$ with valid bounds, and the proofs of
[Lemma 13.4](#lem-13-4) and of the theorems apply with the seeds in
place of $R$; the covering of $X'(R)$ by the seeds follows, as in (i),
from (A1) applied to the nodes that the previous search branched.
(iii) An unknown bound is valid for any problem, and the seeds cover
$X'(R)$ by the same argument, since (A1) for the new subsets is again
assumed of the children that the previous search generated (all the
leaves, the infeasible ones included, being given back). $\square$

Which leaves are given back, and with which bounds, is left to the
caller, since only the caller
knows what has changed. After a restriction of the region the bounds of
the pruned leaves remain valid, and the infeasible leaves need not be
explored again. After an enlargement the infeasible leaves have to be
given back with an unknown bound. After a change of the objective all of
them have to be given back with an unknown bound, unless the old ones
can be corrected (say, by a bound on the change of the objective over
the node). A search stopped by a limit is resumed by giving back its
open leaves together with its incumbent, as in [Figure
13.14](#fig-13-14). The bound of the new `BBSolution` is that of the
subtrees of the seeds, and the caller combines it with the bounds of the
leaves that it leaves out. The test `bb` and the parallel tests check
three uses on the knapsack. First, the leaves of a search carried out to
the end bound the optimum and solve the knapsack of three quarters of
the capacity. Second, with `explicit_children`, the leaves of the
smaller knapsack, the infeasible ones among them, solve the original one
again (without `explicit_children` the child that takes an item that
did not fit is never generated, and (A1) fails for the larger
capacity). Third, the open leaves of a search stopped by a node limit end
it.

In the parallel searches each worker collects its leaves, which are
merged at the end, and the seeds go to their owners as the children do.
Over MPI each process keeps the leaves that it owns and starts from its
seeds, which avoids sending them across the processes. We remark that
in a tree where every inner node has two children the leaves are one
more than the inner nodes, and therefore keeping them may take much
more memory than the open nodes; hence, `keep_frontier` is off by
default.

## 13.12 The interface for a branch-and-bound solver {#sec-13-12}

A B&B solver built on a modeling system (e.g., for Mixed-Integer Linear
Programs) would provide a problem whose `node_type` is the list of the
branching decisions from the root, e.g., the bounds changed on the
variables, which is "small" and hence "cheap" to copy and to send. Its
`make_worker( i )` would return an object holding a copy of the model
and of the solver of the relaxation. The `evaluate` of a worker applies
the decisions of the node to its model (possibly undoing only those
that differ from the previous node, which makes a dive "cheap") and
solves the relaxation. It returns the value of the relaxation as the
bound, its infeasibility, and the same value as a feasible one if the
solution is integral. Its `branch` chooses a fractional variable in the
solution stored in the node, and gives the two children with the bound
of the parent as their estimate. Its `cutoff( z )` should keep the
incumbent value (and pass it to the solver of the relaxation as an
objective limit), and its `evaluate` should return as feasible only the
values found at the node, keeping the solution in the node;
then the outcome is that of the scheme it runs, whichever it is. In
particular, the deterministic rounds and the SelfSplit then give a
solver whose runs can be repeated exactly, under the conditions of
[Theorem 13.18](#thm-13-18) and [Theorem 13.23](#thm-13-23); a warm
start from the previous node of the worker satisfies them only in their
weaker form, in which the outcome depends on $p$. Such a solver calls
`parallel_branch_and_bound()` with the name of a scheme read from its
parameters, as in [§13.10](#sec-13-10), and reads the incumbent, its
node, the bound and the gap from the `BBSolution`. In the options it
gives the selection, the tolerances, the initial incumbent and the
callback, and for the parallel searches the number of workers, the
batch, the nodes of a round and those of the split. The test
`bb-cutoff` checks this interface with workers that fence the nodes by
the incumbent they receive.

## 13.13 Benchmark {#sec-13-13}

The schemes are measured by
[`tools/BranchAndBoundBench.cpp`][BranchAndBoundBench.cpp], run by the table F
of [`tools/bench/campaign.sh`][campaign.sh], on two problems, the problems
`knapsack` and `tsp-bb` of that table, each of which gives one part of Table
13.2. The first is made of random 0-1 knapsacks of the strongly correlated
class of Pisinger with 90 items of weight at most 1000 (the seeds 3, 8 and 13),
with values of type `double`; the feasible values are exact, since the
profits are integer and their sums much smaller than $2^{53}$, while the
fractional bounds are rounded. Our second problem is the tree of the
partial tours of random traveling salesperson problems with 24 cities of
the unit square (the seeds 2, 3 and 4), whose B&B is the `TourTree` of
[§13.6](#sec-13-6), with the bound $g + h$ of the minimum spanning tree
heuristic ([Chapter 5](05-heuristics.md#ch-5)) in double precision. In
all the runs the selection is best-first. The sequential search runs on
one core, and each parallel scheme on $t$ cores, for the numbers of
cores in the header of the table, with one worker per core, the cores of
one NUMA node and the memory on that node. The shared queue, the
distributed queues, the rounds (one node per worker and round), the
SelfSplit (with and without the shared incumbent), the work stealing and
the MultiQueue run on each of the three executors; the farm runs with $t
- 1$ workers and its master on FastFlow, and the search over MPI with
one process per core. Each run is a separate process, and each is run
twice; the gap between the two is the A/A noise of the machine, which
the note of the table gives.

The time is the wall-clock time measured by the `Stopwatch` of the
search, which starts when `solve()` starts. Hence, it does not include
the generation of the instance, the sorting of the items, the dynamic
program that checks the optimum, nor, for MPI, the start of the
processes and `MPI_Init`. It does include the start of the threads of
the parallel searches, since the benchmark gives them no executor and
each search builds a pool. As for the machine and the compiler, they are those
stated in the note of the table.

**Table 13.2.** The B&B schemes on the two problems of table F
described above, in two parts, `knapsack` and `tsp-bb`. In each part a
row is a search and a column a number of cores; each cell is the time
in seconds summed over the three instances of the problem, each the
median of its two runs, and in parentheses the speedup over the row
`B&B sequential`, i.e., the sequential search on one core; an empty
cell is a search that did not solve every instance. The note under each
part gives the machine and the A/A noise.

--8<-- "tables/branch-and-bound-schemes.md"

The table is a placeholder until it is measured on an idle machine.
Our benchmark also writes, for each run, the counters by which one can
read the anomalies of [§13.7](#sec-13-7). They are the nodes evaluated
(column `expanded` of the CSV, in the format shared with the benchmark
of A$^*$), branched (`reexpanded`) and pruned (`duplicates`), the
largest queue (`stored`), the load balance LB (the largest number of
evaluations of a worker over their mean) and the communication overhead
CO (the nodes sent to another worker over those generated). With
`--slow` each evaluation first performs a fixed amount of busy work, as
a bound computed by an "expensive" relaxation would; this moves the
comparison towards the regime the layer is meant for.

## 13.14 Notes and references {#sec-13-14}

Branch-and-bound goes back to Land and Doig
[[146]](Z-bibliography.md#ref-146), who solved discrete programs by
partitioning the feasible region and bounding each part by the Linear
Programming relaxation. The bound of the knapsack used here is the
fractional relaxation of Dantzig [[147]](Z-bibliography.md#ref-147), and
the classes of the random instances are among those that Pisinger
[[148]](Z-bibliography.md#ref-148) discusses in his study of which
knapsacks are "hard". Nau, Kumar and Kanal
[[17]](Z-bibliography.md#ref-17) gave the general formulation of B&B
of which A$^*$ and AO$^*$ are special cases, which is the view of
[§13.6](#sec-13-6); depth-first B&B on a graph is treated with the
other memory-bounded searches in [Chapter 4](04-memory-bounded.md#ch-4).
The textbook of Jinnai [[164]](Z-bibliography.md#ref-164) introduces B&B
as a depth-first search pruned by the incumbent (its section 3.4) and
gives the recursive DFBnB with the successors sorted by $h$ (its
section 6.2, Algorithms 8 and 9); in its Algorithm 9 a successor $v$ of
$u$ is tested by $g + h(v) < U$, while its bound is
$g + c(u, v) + h(v)$, which prunes more.

For the parallel schemes, the survey of Gendron and Crainic
[[16]](Z-bibliography.md#ref-16) gives the classification by type used
in [§13.7](#sec-13-7) and the literature on the anomalies, which Lai
and Sahni [[149]](Z-bibliography.md#ref-149) analyzed early on. The
deterministic rounds follow the deterministic mode of FiberSCIP
[[18]](Z-bibliography.md#ref-18), and the SelfSplit is that of
Fischetti, Monaci and Salvagnin [[19]](Z-bibliography.md#ref-19). The
work stealing is that of Blumofe and Leiserson
[[20]](Z-bibliography.md#ref-20), with the steal-half of Hendler and
Shavit [[36]](Z-bibliography.md#ref-36) applied to priority queues.
The MultiQueue is that of Rihani, Sanders and Dementiev
[[37]](Z-bibliography.md#ref-37), whose use as a priority scheduler
Postnikova et al. [[38]](Z-bibliography.md#ref-38) studied. Finally,
the farm follows the farm with feedback of FastFlow
[[12]](Z-bibliography.md#ref-12), and the end of the search over MPI is
detected by the four-counter method of Mattern
[[11]](Z-bibliography.md#ref-11).

With respect to the parallel B&B solvers of Mixed-Integer Programming,
the layer is a framework of node-level schemes, and it is behind them in
what concerns the problem: it has no cuts, no presolve, no branching
rules and no primal heuristics. It leaves to the problem the parallelism
inside an evaluation (type 1 of [[16]](Z-bibliography.md#ref-16)), and
racing, i.e., the search of several trees with different settings (type
3), is not provided; also, the SelfSplit does its sampling once instead
of in each worker. On the other hand, each scheme runs on any problem
that models the concept, and on the three executors, with no change of
the program. The deterministic rounds and the SelfSplit give repeatable
runs under the stated conditions, the bound reported is valid whatever
the status, and the leaves of a search can seed the next one. These
properties are checked by the tests of the library, while we have not
measured their cost against a specific solver.

## 13.15 Exercises {#sec-13-15}

**Exercise 13.1** (`ctest -R PBeFS.bb`). Run Algorithm 13.3 by hand on
the example of [Figure 13.1](#fig-13-1) with the breadth-first and the
best-estimate selections, and check the orders of
[Figure 13.4](#fig-13-4). Then prove that the value returned is the
same for each selection, and explain why the numbers of evaluations
differ; the test checks the optimum of each selection, in both senses,
against the dynamic program.

**Exercise 13.2** (`ctest -R PBeFS.bb`). With `shared_incumbent`
returning the optimum, the sequential search returns the status
`infeasible`, no node, a bound at most the optimum, and as many
evaluations as with `incumbent` set to the optimum. Explain each of the
four facts from Algorithm 13.3 and [Theorem 13.6](#thm-13-6), and say
what changes when the value returned by `shared_incumbent` is worse
than the optimum.

**Exercise 13.3** (`ctest -R PBeFS.bb-values`). Show that with an
integer `value_type` and bounds rounded towards the inside the test
PRUNABLE prunes, besides the nodes that the double values prune, those
whose bound is within less than 1 of the incumbent. Then explain why
this does not imply that each instance needs fewer evaluations with
integer bounds (hint: the rounding creates ties among the keys of
best-first), and why the test compares the evaluations summed over its
instances.

**Exercise 13.4** (`ctest -R PBeFS.bb-shared-thread`). Suppose that in
Algorithm 13.7 a worker lowered `busy` right after PROCESS, outside the
lock, and pushed the children afterwards. Give an interleaving of two
workers that ends the search with a node still to be processed, and
show which statement of [Proposition 13.15](#prop-13-15) fails. The
test checks the algorithm as written: with 1 to 8 workers and every
selection it ends with the optimum of the dynamic program.

**Exercise 13.5** (`ctest -R PBeFS.bb-mq-thread`). In Algorithm 13.8,
suppose that $m$ were raised by $|K|$ after the pushes and lowered by 1
for the node processed before them. Show that $m = 0$ could then be read
while a node is in a queue, and that the search could end with a
suboptimal value. The test checks the algorithm as written, which ends
with the optimum with 1 to 8 workers and every selection.

**Exercise 13.6** (`ctest -R PBeFS.bb-det-thread`). The test checks
that, with rounds of 4 nodes, the outcome on the knapsack is the same
with 1, 2, 3 and 8 workers and at each repetition; explain why by
[Theorem 13.18](#thm-13-18), whose first assumption the knapsack
satisfies. Then give an evaluator whose outcome depends on the nodes
its worker evaluated before (e.g., one that reuses the bound of its
previous node when it is the parent), and show on a "small" tree that
the rounds then evaluate different nodes with $p = 1$ and $p = 2$, while
each $p$ repeats exactly; this is the weaker form of the theorem, which
the test does not cover.

**Exercise 13.7** (`ctest -R PBeFS.bb-cutoff`). Explain why a worker
of the deterministic rounds must take the incumbent from `cutoff( z )`
rather than from a value it shares with the other workers, by an
example in which the relaxation of a node stops early against a value
that another worker has just found. The test checks, with workers that
find infeasible each node whose bound cannot exceed the $z$ they
receive, that the rounds give the same value, bound and numbers of
evaluated and infeasible nodes with 1 to 8 workers.

**Exercise 13.8** (`ctest -R PBeFS.bb-split-thread`). Prove that the
SelfSplit with one worker evaluates exactly the nodes of the sequential
search, and that with $p$ workers and without `split_shared` two runs
evaluate the same nodes; the test checks both on the knapsack. Then
show, on a tree in which the optimum is in the subtree of the last
worker, how many nodes the other workers evaluate without the dive and
with it (a question that the test does not cover).

**Exercise 13.9** (`ctest -R PBeFS.bb-scheme`). Write a program that
reads the name of a scheme from its command line and solves a knapsack
by `parallel_branch_and_bound()`, first from the root and then, with
three quarters of the capacity, from the leaves of a search carried out
to the end on the original capacity; explain, by
[Proposition 13.24](#prop-13-24), which leaves it has to give back and
with which bounds. The test does both with each named scheme, giving
back every leaf with an unknown bound.

**Exercise 13.10** (`ctest -R PBeFS.bb-steal-thread`). Show that the
counter $C$ of [Proposition 13.22](#prop-13-22) must count a thief
from before its exchange, by an interleaving in which a thief that
entered $C$ after a successful exchange lets another worker see
$C = 0$ while the stolen nodes are open. The test checks the scheme as
written, which ends with the optimum with 1 to 8 workers and every
selection.

<!-- the links to the sources, written by srclinks.py -->

[BranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h
[BranchAndBound.h:1011]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L1011
[BranchAndBound.h:137]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L137
[BranchAndBound.h:205]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L205
[BranchAndBound.h:262]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L262
[BranchAndBound.h:315]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L315
[BranchAndBound.h:338]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L338
[BranchAndBound.h:386]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L386
[BranchAndBound.h:476]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L476
[BranchAndBound.h:505]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L505
[BranchAndBound.h:700]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L700
[BranchAndBound.h:808]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L808
[BranchAndBound.h:879]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBound.h#L879
[BranchAndBoundBench.cpp]: https://github.com/dmeoli/PBeFS/blob/main/tools/BranchAndBoundBench.cpp
[BranchAndBoundScheme.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/BranchAndBoundScheme.h
[campaign.sh]: https://github.com/dmeoli/PBeFS/blob/main/tools/bench/campaign.sh
[DeterministicBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/DeterministicBranchAndBound.h
[DeterministicBranchAndBound.h:125]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/DeterministicBranchAndBound.h#L125
[FarmBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/FarmBranchAndBound.h
[FarmBranchAndBound.h:104]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/FarmBranchAndBound.h#L104
[Knapsack.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/Knapsack.h
[MPIBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MPIBranchAndBound.h
[MPIBranchAndBound.h:107]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MPIBranchAndBound.h#L107
[MPIBranchAndBound.h:225]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MPIBranchAndBound.h#L225
[MultiQueueBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MultiQueueBranchAndBound.h
[MultiQueueBranchAndBound.h:105]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MultiQueueBranchAndBound.h#L105
[MultiQueueBranchAndBound.h:129]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/MultiQueueBranchAndBound.h#L129
[ParallelBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h
[ParallelBranchAndBound.h:136]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L136
[ParallelBranchAndBound.h:190]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L190
[ParallelBranchAndBound.h:296]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L296
[ParallelBranchAndBound.h:321]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L321
[ParallelBranchAndBound.h:350]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L350
[ParallelBranchAndBound.h:599]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/ParallelBranchAndBound.h#L599
[SelfSplitBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/SelfSplitBranchAndBound.h
[SelfSplitBranchAndBound.h:130]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/SelfSplitBranchAndBound.h#L130
[StealingBranchAndBound.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/StealingBranchAndBound.h
[StealingBranchAndBound.h:117]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/bnb/StealingBranchAndBound.h#L117
[test.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/test.cpp
