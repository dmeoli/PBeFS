# Branch-and-bound as best-first search

We are concerned with the minimization of a function $v$ over a finite (but
possibly huge) set $X$, i.e., with finding

$$
z^* = \min \{ v(x) : x \in X \} ,
$$

when $X$ can be partitioned recursively into smaller sets, each of which comes
with a lower bound on the values that it contains; a maximization is the
minimization of $-v$, which is how the library treats it internally. This page
gives the notions on which the branch-and-bound (B&B) of the library rests (the
nodes, the bounds, the incumbent, the pruning and the selection rules), the
invariant that makes the search correct and the bound that it reports when it
stops early. Then it shows A* as a special case of the same scheme (cf.
[[17]](../../README.md#ref-17)), and it gives the cost of the operations on the
queue of the library. The [design note](../branch-and-bound.md) discusses how
the concept `BranchAndBound` of
[BranchAndBound.h](../../include/PBeFS/bnb/BranchAndBound.h) maps these notions
onto a problem, and how a solver built on a modeling system would use it. The
parallel schemes are in
[Parallel branch-and-bound](parallel-branch-and-bound.md), and the conventions
of the pseudocode are those of
[Best-first search and A*](best-first-search.md).

## Nodes, bounds and the incumbent

A *node* $n$ stands for a subset $X(n) \subseteq X$, the root for $X$ itself,
and *branching* $n$ gives children $n_1, \ldots, n_k$ whose subsets cover
$X(n)$ (usually they partition it). Clearly, the tree of the nodes is finite as
long as the branching makes the subsets strictly smaller. The *evaluation* of
$n$ gives a *bound* $\beta(n)$ with

$$
\beta(n) \leq \min \{ v(x) : x \in X(n) \} ,
$$

e.g., the optimal value of a relaxation of the problem restricted to $X(n)$,
together with two optional facts: that $X(n)$ is empty (the node is
*infeasible*), or a feasible value $v(x)$ for some $x \in X(n)$ (e.g., when the
solution of the relaxation is feasible, or by a heuristic). The best feasible
value found so far is the *incumbent* $z$ ($+\infty$ if none is known, unless a
value is given at the start), and the incumbent node, if any, is a node of
value $z$.

Since a child is a subset of its parent, the bound of the parent is also a
bound of the child. For each open node the library keeps the bound
$\max \{ \beta(\text{parent}), e(n) \}$, where $e(n)$ is an *estimate* of the
bound of the child that the branching gives at little cost (e.g., by fixing one
variable in the solution of the relaxation of the parent). Thus, the evaluation
of the node, which is usually "expensive", is postponed until the node is
extracted, and its result is combined with the inherited bound by taking the
larger of the two.

A node is *pruned*, i.e., discarded with its subtree, when it is infeasible or
when its bound cannot improve the incumbent by more than the tolerance, i.e.,
when

$$
\beta(n) \geq z - \max \{ \tau_a , \tau_r |z| \} ,
$$

with an absolute tolerance $\tau_a \geq 0$ and a relative one $\tau_r \geq 0$
(`abs_gap` and `rel_gap` of `BBOptions`). With $\tau_a = \tau_r = 0$ the rule
is $\beta(n) \geq z$, since a subtree whose bound equals $z$ cannot hold a
better value. When the values of the solutions are integer, $\lceil \beta(n)
\rceil$ is a bound too, and the rule applied to it prunes every node with
$\beta(n) > z - 1$, i.e., those that could improve $z$ by a fraction only;
the library does so when the values are of an integer type and the problem
rounds its bounds (cf. `BasicKnapsack`, which rounds the bound of Dantzig).
In the library the rule is tested at three times: (i) when a
child is generated, against its inherited bound, (ii) when a node is extracted,
against the incumbent that may have improved since its generation, and (iii)
after its evaluation, against its own bound. A node that survives is branched,
also when it has a feasible value, since its subtree may hold better values
than the one found.

## Correctness and the reported bound

The search keeps the following invariant: *each $x \in X$ with
$v(x) < z - \max \{ \tau_a, \tau_r |z| \}$ belongs to $X(n)$ for some open node
$n$*. In fact, it holds at the start, when the root is the only open node, and
each step preserves it. Branching replaces a node by children that cover it, an
infeasible node has no element, and a node is pruned only when its bound shows
that its subset has no element below the threshold (which can only decrease
when $z$ improves). Hence, when no open node is left, the incumbent is optimal
within the tolerance, i.e., $z - z^* \leq \max \{ \tau_a, \tau_r |z| \}$, with
$z = +\infty$ meaning that $X$ is empty (if no value was given at the start).
Also, the search can stop as soon as the least bound of the open nodes is
within the tolerance of $z$, since then all of them would be pruned at their
extraction. This is the termination test of all the searches of the library,
which also end when a node or time limit is reached.

When a limit stops the search, or when the tolerances are positive, the
invariant gives a bound on the optimum, which the `BBSolution` reports:

$$
z^* \geq \min \{ z , \ \min_{n \text{ open}} \beta(n) ,
                 \ \min_{n \text{ pruned}} \beta(n) \} ,
$$

where the last term accounts for the subtrees pruned with $\beta(n) < z$, which
a positive tolerance allows. Accordingly, the reported gap is $z$ minus this
bound, and it is valid whatever the status of the search. With `keep_frontier`,
the leaves of the tree explored (pruned, feasible and not branched, infeasible,
or left open) are also kept. Their subsets cover those of the nodes the search
started from, and a new search started from them explores the same region,
which is the basis of the reoptimization discussed in the
[design note](../branch-and-bound.md).

## The selection rules

The invariant holds whatever the order in which the open nodes are extracted.
Therefore, the order is a policy (`Selection`), which decides how many nodes
are evaluated and how much memory the queue takes. With *best-first* the node
of least bound is extracted first. If the bound of a node does not depend on
the order of the search, each node whose bound is less than $z^*$ is evaluated
by any correct search that uses the same bounds; best-first evaluates only
those (and possibly some with bound $z^*$), and it is thus the search with the
fewest evaluations. Its price is a queue that may grow to a large part of the
tree before an incumbent is found.

*Depth-first* extracts the deepest node, which finds feasible values early and
keeps the queue proportional to the depth times the number of the children, but
it may evaluate many nodes that best-first would prune. *Breadth-first*
extracts the shallowest node (the oldest among those of the same depth), and
*best-estimate* the node of least estimate of the best feasible value in its
subtree, which the problem gives. *Plunging* is the usual compromise: the best
child of the node just branched is extracted next (a dive) as long as its bound
is within a fraction $\gamma$ (`plunge_gap`) of the gap from the least open
bound $\underline{\beta}$, i.e., while

$$
\beta(n) \leq \underline{\beta} + \gamma \, (z - \underline{\beta}) ,
$$

and always while there is no incumbent, after which the search goes back to
best-first. A dive is "cheap" when the evaluator keeps the state of the last
node (e.g., a relaxation solved again from the basis of its parent).

## A* as a special case

A* is the B&B in which a node is a path from $s_0$ to a state $u$ and its
subset is the set of the paths to a goal that extend it. Its bound is
$f(u) = g(u) + h(u)$, which is a lower bound on the cost of these paths if $h$
is admissible, and a goal gives the feasible value $g(u)$ (cf.
[[17]](../../README.md#ref-17)). In this view, the best-first selection with
the ties broken by larger $g$ is then the selection of A*, and the test at the
extraction of a goal is the pruning by the incumbent, since a goal extracted
with $f = g$ has the least bound of the open nodes.

The two searches differ in two respects. First, a graph has *duplicates*, i.e.,
several paths to the same state, of which only the cheapest matters: A* detects
them with CLOSED and keeps one node per state, while a B&B tree has no
duplicates and needs no closed list. Second, A* stops at the first goal
extracted, while a B&B goes on until no open node can improve the incumbent; in
fact, its bound does not necessarily grow along the tree as $f$ does with a
consistent heuristic, and a node with a feasible value may hold better values
below. This is why the two layers of the library share the parallel schemes and
the runtimes, but not the node table and the OPEN lists of A*, which are
specialized for integer costs and for the duplicates.

## The cost of the operations

The open nodes of a B&B are kept in an array of slots with a free list, since a
node is destroyed when it is branched. A binary heap of the pairs (key, slot)
gives the next node (cf. `BBQueue` in
[BranchAndBound.h](../../include/PBeFS/bnb/BranchAndBound.h)); the key is the
bound with best-first and plunging, the opposite of the depth with depth-first,
the depth with breadth-first and the estimate with best-estimate. Among equal
keys, the ties are broken by the order of the pushes (the last pushed first,
but the first pushed with breadth-first). With $N$ open nodes, a push and a pop
cost $O(\log N)$ comparisons and the move of one node. When the key is not the
bound, a second heap of the bounds gives the least open bound
$\underline{\beta}$, which the termination test and the reported bound need;
its entries are deleted lazily when they reach the top, which costs one more
push per node and the pops of the stale entries, i.e., $O(\log N)$ amortized
time per node. These costs are negligible when the evaluation solves a
relaxation, which is the case the layer is meant for, and they matter only for
problems whose bound costs about as much as a few comparisons. The pseudocode
of the sequential search is on the page of
[BranchAndBound.h](../../include/PBeFS/bnb/BranchAndBound.h).
