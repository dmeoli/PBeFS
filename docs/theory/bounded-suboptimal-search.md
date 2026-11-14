# Bounded-suboptimal, anytime and incremental search

We now turn to the searches that give up the optimality of the path for a
bound on its cost, those that improve a first path until they prove it
optimal, and the one that repairs a search after the costs have changed. The
notation is that of [Best-first search and A*](best-first-search.md): $C^*$
is the cost of an optimal path, $g^*(n)$ that of a cheapest path from the
root to $n$, $h$ is admissible unless said otherwise, $z$ is the incumbent
(the cost of the best path found so far, $\infty$ if none) and $w \geq 1$ the
weight; a search is *$w$-admissible* if the path it returns costs at most
$w C^*$. The pseudocode of each search is on the page of its header, and the
conditions under which a parallel search may stop are those of
[Parallel A*](parallel-astar.md), which this page extends to a bound.

## Weighted search in any order of expansion

Weighted A* (wA*) orders OPEN by $f' = g + w h$, which trades the optimality
of the path for fewer expansions (cf. [[46]](../../README.md#ref-46)); its
first goal selected costs at most $w C^*$, as the classical argument shows.
The parallel searches, however, do not expand in the order of $f'$, and the
first goal they find may be much worse; we therefore state the rules of
`Weight.h` in a form that does not depend on the order of the expansions.
A node $n$ is *pruned* by $z$ when it is dropped without being expanded, and
the search *ends* when OPEN has no node that the rules below keep.

**Proposition 1.** Let a search (sequential or parallel) reopen every closed
node reached through a cheaper path, prune a node when $w f(n) \geq z$, and
end when no node of OPEN has $f'(n) < z$ and no node is in transit between
the workers. If it ends with $z < \infty$, then $z \leq w C^*$.

*Proof.* Let $s_0, \dots, s_k$ be an optimal path, and let $s_i$ be its first
node that is not closed with $g(s_i) = g^*(s_i)$ (if there is none, $s_k$ has
been expanded with $g^*(s_k) = C^*$, and $z \leq C^*$). Its predecessor
$s_{i-1}$ was expanded with $g^*(s_{i-1})$, and therefore offered
$g^*(s_i)$ to $s_i$, which, since the search reopens, has
$g(s_i) = g^*(s_i)$ and is either in OPEN, in transit, or was pruned. In
transit it cannot be at the end, and in OPEN it has
$z \leq f'(s_i) \leq w f(s_i)$; pruned, it has $z \leq w f(s_i)$ directly. In
all cases $z \leq w (g^*(s_i) + h(s_i)) \leq w C^*$, by the admissibility
of $h$. $\square$

The test of the end, on $f'$, is weaker than the pruning, on $w f \geq f'$:
a node may be kept in OPEN and yet never be expanded, which is why a worker
of HDA* drops it when it reaches the top of its OPEN rather than stopping.
With $w = 1$ both rules reduce to those of the optimal parallel searches,
which is the case $w = 1$ of Proposition 1. The pruning by $w f$ is that of
wPBNF (cf. [[3]](../../README.md#ref-3)), and it prunes more than the
pruning by $f'$, since $w f \geq f'$; the lower bound that a search ending
with $z$ proves is $z / w$ (rounded up with integer costs), which
`Solution::bound` reports.

## Duplicate dropping

A weighted search may reach a closed node through a cheaper path much more
often than A* does, and reopening it costs expansions that rarely change the
path found; the *duplicate dropping* variants keep such a node closed. In a
sequential search with a consistent heuristic this keeps the bound (cf.
Theorem 1 of [[48]](../../README.md#ref-48)), but the argument relies on the
order of the expansions, and a parallel search may close a node with a $g$
much larger than $g^*$ before the cheaper path arrives. We therefore reopen
a closed node $v$, reached from $u$ with $g' < g(v)$, only when the saving
exceeds a budget that the weight allows, i.e., when

$$
g(v) - g' > (w - 1) \bigl( h(u) - h(v) \bigr) ,
$$

and we prune by $f'$ instead of $w f$; with $w = 1$ the rule reopens every
node reached through a cheaper path.

**Proposition 2.** Let a search reopen as above, prune a node when
$f'(n) \geq z$, and end as in Proposition 1. If it ends with $z < \infty$,
then $z \leq w C^*$, whatever the order of the expansions and with any
admissible heuristic.

*Proof.* Along an optimal path $s_0, \dots, s_k$, consider the last expansion
of $s_i$, with $g(s_i)$, which offers $g(s_i) + c(s_i, s_{i+1})$ to
$s_{i+1}$: after it, $g(s_{i+1}) \leq g(s_i) + c(s_i, s_{i+1}) + (w - 1)
(h(s_i) - h(s_{i+1}))$, since either $s_{i+1}$ takes the new path or the
saving it refuses is at most the budget, and $g(s_{i+1})$ only decreases
afterwards. Summing from $s_0$, where $g = 0$, every node $s_i$ that has been
expanded has $g(s_i) \leq g^*(s_i) + (w - 1)(h(s_0) - h(s_i))$, and so does
the first node $s_j$ of the path that is not closed, which was offered its
$g$ by $s_{j-1}$. Hence

$$
f'(s_j) = g(s_j) + w h(s_j) \leq g^*(s_j) + h(s_j) + (w - 1) h(s_0)
\leq C^* + (w - 1) C^* = w C^* ,
$$

and, as in Proposition 1, $s_j$ in OPEN or pruned gives $z \leq f'(s_j)$.
$\square$

With a consistent heuristic $h(u) - h(v) \leq c(u, v)$, and therefore the
budget is at most $(w - 1) c(u, v)$, the one of the weak duplicate dropping
of [[3]](../../README.md#ref-3) (section 5.2, Theorem 5), which keeps the old
path of $v$ when $g(v) \leq g(u) + w c(u, v)$: the rule above reopens fewer
nodes, it needs only the heuristic of the parent, which every node carries
(the MPI search computes it again from the state of the parent, only when a
closed node improves), and, unlike the weak criterion, it keeps the bound
with an inconsistent heuristic. In both cases the least $f$ of OPEN is no
longer a lower bound on $C^*$ (cf. Corollary 3 of
[[3]](../../README.md#ref-3)), and the bound reported is $z / w$. The wA*
that never reopens a closed node (cf. [[49]](../../README.md#ref-49), p. 357)
is the special case of an infinite budget, which keeps the bound only with a
consistent heuristic and in the order of $f'$.

## Anytime search

An anytime search returns a first path quickly and then improves it. Anytime
Weighted A* (AWA*, [[47]](../../README.md#ref-47)) is wA* that, after the first
goal, goes on with the incumbent pruning by $f \geq z$ instead of $w f$, and
ends when OPEN is empty; by Proposition 1 with $w = 1$ its last incumbent is
optimal, and every path it reports comes with the lower bound
$\min(z, \min_{n \in \mathrm{OPEN}} f(n))$, which is at most $C^*$ since a
node of an optimal path with $g = g^*$ is in OPEN (or $z = C^*$). We keep
this minimum with a second OPEN ordered by $f$ over the same nodes, whose
entries become stale with those of the first. The parallel searches are
anytime in the same way (anytime weighted HDA*, PBNF, SPA*, MQA* and farm),
except that their intermediate paths carry only the bound $h(s_0)$, since no
worker knows the least $f$ of the others; the first path of a parallel
search, unlike that of AWA*, need not be within $w C^*$.

Anytime Repairing A* (ARA*, [[48]](../../README.md#ref-48)) runs wA* with a
decreasing schedule $w_0 > w_1 > \dots \geq 1$ and reuses the work of each
iteration in the next: within an iteration no node is expanded twice, a
closed node whose $g$ decreases goes to the list INCONS, and at the end of
the iteration OPEN $\cup$ INCONS, ordered by the next weight, starts the next
one. With a consistent heuristic every iteration ends with $z \leq w_i C^*$
(Theorem 1 of [[48]](../../README.md#ref-48)); the bound
$L = \min(z, \min_{n \in \mathrm{OPEN} \cup \mathrm{INCONS}} f(n))$, whose
ratio $z / L$ is the $\epsilon'$ of the paper, holds with any admissible
heuristic, by the argument of Proposition 1 applied to the first node of an
optimal path with $g > g^*$ or not closed. Since our OPEN has no
decrease-key, the union is formed by emptying OPEN at the end of the
iteration, and the same pass pushes its live nodes again with the keys of the
next weight.

## Searches with a second order

The optimistic search ([[49]](../../README.md#ref-49)) orders OPEN by
$g + w_o h$ with an optimistic weight $w_o > w$ ($2 w - 1$ by default), and,
once it has a path, expands by $f$ until $w \min_{\mathrm{OPEN}} f \geq z$.
Since $\min_{\mathrm{OPEN}} f \leq C^*$ (with reopening), it ends with
$z \leq w C^*$; the order of $w_o$ only decides how soon. The Improved
Optimistic Search ([[50]](../../README.md#ref-50)) proves the first path
differently: it restarts an A* from the root that reopens each node of the
first search at most once, and it lowers the cost of the first path by the
largest saving found on one of its nodes; the largest $f$ selected by this A*
is a lower bound on $C^*$, which needs a consistent heuristic.

The focal searches look, among the nodes with $f \leq w f_{\min}$ (the
FOCAL list), for the one that a second order prefers: the distance $d$ to a
goal in edges in $A^*_\epsilon$ ([[51]](../../README.md#ref-51)), and an
estimate $\hat d$ corrected online, with OPEN by $\hat f$ and a third list by
$f$, in Explicit Estimation Search ([[52]](../../README.md#ref-52)). Both are
$w$-admissible because the goal they select has $g \leq w f_{\min}$ (in EES,
$g = \hat f \leq w f$ of the best node by $f$, or the goal is the best node
by $f$), and $f_{\min} \leq C^*$ with reopening. Dynamic Potential Search
([[53]](../../README.md#ref-53)) orders all of OPEN by the potential
$(w f_{\min} - g) / h$, which changes with $f_{\min}$ and is then recomputed
for all the nodes; its bound follows in the same way (see `PotentialSearch.h`
for the case $h = 0$).

## Incremental search

Lifelong Planning A* (LPA*, [[54]](../../README.md#ref-54)) keeps, besides
$g$, the one-step lookahead $\mathrm{rhs}(s) = \min_{s'} g(s') + c(s', s)$
over the predecessors of $s$, and expands only the states where the two
differ; after a change of some costs, the states whose lookahead changes are
the heads of the edges changed, and the next search repairs only what the
change has made inconsistent. With a consistent heuristic the first search
expands the states that A* expands, each at most twice, and a repair expands
the states whose $g$ has changed together with those that the heuristic
cannot rule out (cf. section 5 of [[54]](../../README.md#ref-54)). The goals
of a domain, which may be many, are joined to a virtual goal by edges of cost
0, so that the termination test compares the least key of the queue with the
least $g$ of the goals reached.
