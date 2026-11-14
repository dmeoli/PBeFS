# Bounded-suboptimal, anytime and incremental search

We now turn to three families of searches: those that give up the optimality
of the path for a bound on its cost, those that improve a first path until
they prove it optimal, and one that repairs a search after some costs have
changed. The notation is that of [Best-first search and A*](best-first-search.md):
$C^*$ is the cost of an optimal path, $g^*(n)$ that of a cheapest path from
the root to $n$, $z$ the incumbent (the cost of the best path found so far,
$\infty$ if none) and $w \geq 1$ the weight. Unless said otherwise $h$ is
admissible, and a search is *$w$-admissible* if the path it returns costs at
most $w C^*$. The pseudocode of each search is on the page of its header;
this page states the rules that keep the bounds and why they hold, and it
extends to a bound the conditions of [Parallel A*](parallel-astar.md).

## Weighted search in any order of expansion

Weighted A* (wA*) orders OPEN by $f' = g + w h$ (cf.
[[46]](../../README.md#ref-46)), which usually saves many expansions at the
price of a path that costs at most $w C^*$. However, the classical argument
relies on the order of $f'$, which the parallel searches do not follow: their
first goal may be much worse than $w C^*$. We therefore state the rules of
`Weight.h` in a form that holds in any order of expansion (a node is
*pruned* by $z$ when it is dropped without being expanded).

**Proposition 1.** Let a search, sequential or parallel, reopen each closed
node reached through a cheaper path, prune a node $n$ when
$w f(n) \geq z$, and end when no node of OPEN has $f'(n) < z$ and no node is
in transit between the workers. If it ends with $z < \infty$, then
$z \leq w C^*$.

*Proof.* Let $s_0, \dots, s_k$ be an optimal path, and $s_i$ its first node
that is not closed with $g(s_i) = g^*(s_i)$ (if there is none, $s_k$ has been
expanded with $g = C^*$, and $z \leq C^*$). Its predecessor $s_{i-1}$ was
expanded with $g^*(s_{i-1})$ and offered $g^*(s_i)$ to $s_i$, which (since
the search reopens) has $g(s_i) = g^*(s_i)$ and is in OPEN or pruned, no node
being in transit at the end. In OPEN it has $z \leq f'(s_i) \leq w f(s_i)$,
and pruned it has $z \leq w f(s_i)$ directly. Hence
$z \leq w (g^*(s_i) + h(s_i)) \leq w C^*$, by the admissibility of $h$.
$\square$

This is, in substance, Theorem 4 and Corollary 2 of
[[3]](../../README.md#ref-3) (section 5.1), which justify the pruning of
wPBNF; we state it with the end of the search, the condition that a
parallel search tests. Note that the end is tested on $f'$ and the pruning on
$w f \geq f'$: a node may stay in OPEN and yet never be expanded, which is
why a worker of HDA* drops it when it reaches the top of its OPEN (rather
than stopping there). With $w = 1$ both rules are those of the optimal
parallel searches. A search ending with $z$ proves the lower bound
$z / w$ (rounded up with integer costs), which `Solution::bound` reports.

## Duplicate dropping

A weighted search reaches a closed node through a cheaper path much more
often than A* does, and reopening it costs expansions that rarely change the
path found; the *duplicate dropping* variants keep such a node closed. In a
sequential search with a consistent heuristic this keeps the bound (cf.
Theorem 1 of [[48]](../../README.md#ref-48)). Yet, the argument relies on the
order of the expansions, and a parallel search may close a node with a $g$
much larger than $g^*$ before the cheaper path arrives. We therefore reopen
a closed node $v$, reached from $u$ with $g' < g(v)$, only when the saving
exceeds a budget that the weight allows, i.e., when

$$
g(v) - g' > (w - 1) \bigl( h(u) - h(v) \bigr) ,
$$

and we prune by $f' \geq z$ instead of $w f \geq z$; with $w = 1$ the rule
reopens each node reached through a cheaper path.

**Proposition 2.** Let a search reopen as above, prune a node $n$ when
$f'(n) \geq z$, and end as in Proposition 1. If it ends with $z < \infty$,
then $z \leq w C^*$, in any order of expansion and with any admissible
heuristic.

*Proof.* Along an optimal path $s_0, \dots, s_k$, consider the last expansion
of $s_i$, with the value $g(s_i)$, which offers $g(s_i) + c(s_i, s_{i+1})$ to
$s_{i+1}$. After it, $g(s_{i+1}) \leq g(s_i) + c(s_i, s_{i+1}) +
(w - 1)(h(s_i) - h(s_{i+1}))$, since either $s_{i+1}$ takes the new path or
the saving it refuses is at most the budget (and $g(s_{i+1})$ only decreases
afterwards). Summing from $s_0$, where $g = 0$, each expanded node $s_i$ of
the path has $g(s_i) \leq g^*(s_i) + (w - 1)(h(s_0) - h(s_i))$, and so does
the first node $s_j$ that is not closed, which $s_{j-1}$ offered its $g$.
Hence

$$
f'(s_j) = g(s_j) + w h(s_j) \leq g^*(s_j) + h(s_j) + (w - 1) h(s_0)
\leq w C^* ,
$$

and, as in Proposition 1, $s_j$ in OPEN or pruned gives $z \leq f'(s_j)$.
$\square$

With a consistent heuristic one has $h(u) - h(v) \leq c(u, v)$, and the
budget is at most $(w - 1) c(u, v)$, the one of the weak duplicate dropping
of [[3]](../../README.md#ref-3) (section 5.2, Theorem 5), which keeps the old
path of $v$ when $g(v) \leq g(u) + w c(u, v)$. Hence, the rule above
reopens fewer nodes, and (unlike the weak criterion) it keeps the bound with
an inconsistent heuristic; it needs only the heuristic of the parent, which
each node carries (the MPI search computes it again from the state of the
parent, and only when a closed node improves). In both cases the least $f$ of
OPEN is no longer a lower bound on $C^*$ (cf. Corollary 3 of
[[3]](../../README.md#ref-3)), and the bound reported is $z / w$. Furthermore, the wA*
that never reopens a closed node (cf. [[49]](../../README.md#ref-49), p. 357)
is the special case of an infinite budget, which keeps the bound only with a
consistent heuristic and in the order of $f'$.

## Anytime search

An anytime search returns a first path quickly and then improves it. Anytime
Weighted A* (AWA*, [[47]](../../README.md#ref-47)) is wA* that, after the first
goal, goes on with the incumbent pruning by $f \geq z$ (instead of
$w f \geq z$), and ends when OPEN is empty; by Proposition 1 with $w = 1$ its
last incumbent is optimal. Each path it reports comes with the lower bound
$\min(z, \min_{n \in \mathrm{OPEN}} f(n))$, which is at most $C^*$ since a
node of an optimal path with $g = g^*$ is in OPEN (or $z = C^*$); we keep this
minimum with a second OPEN ordered by $f$ over the same nodes, whose entries
become stale with those of the first. Also the parallel searches are anytime in the
same way (anytime weighted HDA*, PBNF, SPA*, MQA* and farm), except that
their intermediate paths carry only the bound $h(s_0)$, no worker knowing the
least $f$ of the others; also, the first path of a parallel search need not
be within $w C^*$.

Anytime Repairing A* (ARA*, [[48]](../../README.md#ref-48)) runs wA* with a
decreasing schedule $w_0 > w_1 > \dots \geq 1$ and reuses the work of each
iteration in the next. Within an iteration no node is expanded twice: a closed
node whose $g$ decreases goes to the list INCONS, and OPEN $\cup$ INCONS,
ordered by the next weight, starts the next iteration. With a consistent
heuristic each iteration ends with $z \leq w_i C^*$ (Theorem 1 of
[[48]](../../README.md#ref-48)). The bound
$L = \min(z, \min_{n \in \mathrm{OPEN} \cup \mathrm{INCONS}} f(n))$, whose
ratio $z / L$ is the $\epsilon'$ of the paper, holds with any admissible
heuristic, by the argument of Proposition 1 applied to the first node of an
optimal path that is not closed with $g = g^*$. Since our OPEN has no
decrease-key, the union is formed by emptying OPEN at the end of the
iteration, and the same pass pushes its live nodes again with the keys of the
next weight.

The Anytime Restarting PBNF of [[3]](../../README.md#ref-3) (section 7.1) is
ARA* on the nblocks of SafePBNF, each nblock with an INCONS list. In the
paper the lists are ordered again at each new incumbent, with all the workers
stopped by a dedicated protocol; we do it at the end of an iteration, when no worker
holds an nblock and the lock of the graph suffices. The paper reports that
this search improves its paths more smoothly than the anytime weighted PBNF,
but more slowly (cf. section 8.4).

## Searches with a second order

The optimistic search ([[49]](../../README.md#ref-49)) orders OPEN by
$g + w_o h$, with an optimistic weight $w_o > w$ ($2 w - 1$ by default), and,
once it has a path, expands by $f$ until $w \min_{\mathrm{OPEN}} f \geq z$.
Since $\min_{\mathrm{OPEN}} f \leq C^*$ (with reopening), it ends with
$z \leq w C^*$, and $w_o$ only decides how soon. The Improved Optimistic
Search ([[50]](../../README.md#ref-50)) adds two ideas, the lower bound
$\max (g / w_o + h)$ over the nodes selected by the first search (valid with a
consistent heuristic) and the saving of a cheaper path to a node of the
incumbent path, which lowers the incumbent without re-expanding that path.
The paper goes on from the OPEN of the first search without reopening; we
restart an A* from the root that sees each node anew, since without
reopening the least $f$ of OPEN is possibly above $C^*$ (see
`ImprovedOptimisticSearch.h`).

The focal searches choose, among the nodes with $f \leq w f_{\min}$ (the
list FOCAL), the one that a second order prefers: the distance $d$ to a goal
in edges for $A^*_\epsilon$ ([[51]](../../README.md#ref-51)), and an estimate
$\hat d$ corrected online for Explicit Estimation Search
([[52]](../../README.md#ref-52)), whose OPEN is ordered by an estimate
$\hat f$ of the cost and which keeps a third list by $f$. Both are
$w$-admissible because the goal they select has $g \leq w f_{\min}$ (in EES,
$g = \hat f$ at most $w$ times the least $f$, or else the goal is the node of
least $f$), and $f_{\min} \leq C^*$ with reopening. Dynamic Potential Search
([[53]](../../README.md#ref-53)) orders all of OPEN by the potential
$(w f_{\min} - g) / h$, which changes with $f_{\min}$ and is then recomputed
for all the nodes; its bound follows in the same way (see
`PotentialSearch.h` for the case $h = 0$).

## Incremental search

Lifelong Planning A* (LPA*, [[54]](../../README.md#ref-54)) keeps, besides
$g$, the one-step lookahead $\mathrm{rhs}(s) = \min_{s'} g(s') + c(s', s)$
over the predecessors of $s$, and expands only the states where the two
differ. After a change of some costs only the heads of the edges changed see
their lookahead change, and the next search repairs what the change has made
inconsistent. With a consistent heuristic a search expands each state at most
twice (Theorem 4 of [[54]](../../README.md#ref-54), section 5.1), and a
repair typically expands far fewer states than a search from scratch,
although a change near the root may (somewhat obviously) undo most of the
previous search. The goals of a domain, which may be several, are joined
to a virtual goal by edges of cost 0, so that the test of the end compares
the least key of the queue with the least $g$ of the goals reached.
