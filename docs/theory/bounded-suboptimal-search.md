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
The changes are told with the cost of the edge before and after them, so
that a decrease of $c(u, v)$ only lowers $\mathrm{rhs}(v)$ to
$g(u) + c(u, v)$ if smaller, and an increase recomputes $\mathrm{rhs}(v)$
only when $u$ was its parent (figure 7 of
[[54]](../../README.md#ref-54), lines {33'}-{44'}); a decrease into a
state never reached creates it, since a predecessor of finite $g$ now
reaches it. The goals of finite $g$ are kept in an indexed heap, so that
the test of the end costs $O(1)$ whatever their number. D* Lite
([[143]](../../README.md#ref-143)) is LPA* on the reversed graph, from the
goal to a robot that moves; the heuristic of a state is the estimate
$h(s_r, s)$ from the robot $s_r$, and a move from $a$ to $b$ adds $h(a, b)$
to a scalar $k_m$ that enters the first component of every key, so that,
by the triangle inequality, the keys stored before the move are lower
bounds of the current ones, and a state whose stored key is below its
current one is pushed again instead of being expanded.

## Hash-distributed incremental search

We now run LPA* on the workers of HDA*: the worker $i$ owns the states $v$
with $\mathrm{owner}(v) = i$ and is the only one that writes $g(v)$,
$\mathrm{rhs}(v)$, the parent $\mathrm{par}(v)$ and the place of $v$ in its
queue $U_i$, while the others read $g(v)$ (by an atomic load) only when
they recompute an rhs over all the predecessors (a *rescan*). The goals are
joined to a virtual target $t^*$ by edges of cost 0, and the owner of
$t^*$ keeps, for each goal $t$, the last value of
$m(t) = \min(g(t), \mathrm{rhs}(t))$ that its owner sent, and publishes
$B = \min_t m(t)$ in an atomic word. A state $u$ is *eligible* when it is
inconsistent and either $k(u) < [B + k_m; B]$, or $u$ is a goal and
$k(u) \leq [B + k_m; B]$ (the edge to $t^*$ costs 0, so that a goal at the
key of $t^*$ must be repaired too, or the cost returned may be below the
optimum). Three messages, all with absolute values, cross the workers:
$\mathrm{Relax}(v, x, u)$, sent after the overconsistent expansion of $u$
(which sets $g(u) := \mathrm{rhs}(u)$) to the owner of each successor $v$
with $x = g(u) + c(u, v)$, sets $\mathrm{rhs}(v) := x$ and
$\mathrm{par}(v) := u$ if $x < \mathrm{rhs}(v)$; $\mathrm{Invalidate}(v, u)$,
sent after the underconsistent expansion of $u$ (which sets
$g(u) := \infty$ and rescans $u$), makes the owner of $v$ rescan $v$ if
$\mathrm{par}(v) = u$, and is dropped otherwise; and $\mathrm{Wake}$, sent
by the owner of $t^*$ to the others when $B$ rises, since a worker of the
mailbox of ADS that has found no work becomes active again only by
receiving. An *episode* applies the changes at a barrier (each worker those
of the edges into its states, as LPA* does), then runs a first phase, if
an increase met a parent edge or the key of $t^*$ rose, and a second phase
in which only overconsistent states are expanded, each to quiescence and
separated by a barrier; the first search, and every restart, is the second
phase alone from an empty table, i.e., HDA* on the pair $(g, \mathrm{rhs})$
with no pruning at insertion. In the first phase a state expanded
underconsistent in the episode is *frozen*, i.e., not expanded
overconsistent, until the second phase: by default (`Freeze::raised_parent`)
only when its parent has been expanded underconsistent in the episode too,
with `Freeze::all` always. With one worker every message is applied at once
and the least eligible state is expanded, which is the order of LPA*, and
nothing is frozen.

We remark that the design without the freeze does not always end. Take
$s_0 \to u$, $u \to v$, $v \to u$ and $v \to t$ of cost 1, $h = 0$, $u$ and
$v$ on two workers, and the episode that removes $s_0 \to u$: $u$ is
rescanned at the barrier ($\mathrm{rhs}(u) = g(v) + 1 = 3$), expanded
underconsistent, and lowered to 3 by its own rescan, which reads the old
$g(v) = 2$; then $v$, invalidated, reads $g(u) = 3$, is raised and lowered
to 4, and so on, the two workers counting to infinity as the routers of an
asynchronous distance-vector protocol do, with $B$ rising along and the
keys of $u$ and $v$ staying below it. Sequential LPA* cannot do it, since
its keys do not decrease within a search. The test `hdlpa-*` runs this
episode, which ends.

The hypotheses are the following. (H1) The costs are positive, at least
$c_{\min} > 0$, the graph is locally finite and $h$ is consistent, with
$h(t) = 0$ on the goals, for the costs of each episode. (H2) The costs and
the goals change only at the barriers. (H3) Only $\mathrm{owner}(v)$ writes
$g(v)$, $\mathrm{rhs}(v)$ and $\mathrm{par}(v)$, so that the modification
order of $g(v)$ is the program order of its owner. (H4) Every write of
$g(u)$ is sequenced before the messages it implies, which go to every
successor of $u$, and every change of $m(t)$ of a goal before its message
to $t^*$. (H5) The messages of a worker to another arrive in the order in
which it sent them, and a send happens before the actions of the receiver
after the drain (the channels of ADS are such). (H6) The loads of $g$ are
atomic, so that a load returns neither a write older than one that happens
before it nor one older than what an earlier load of the same thread
returned. A *quiescent point* is a time at which no worker is active and no
message exists, which the counter of the mailbox detects; the barriers are
such. At a quiescent point we want:

- (R) $\mathrm{rhs}(s_0) = 0$ and, for every other state $v$ (those never
  generated included), $\mathrm{rhs}(v) = \min_w g(w) + c(w, v)$ over its
  predecessors, attained by $\mathrm{par}(v)$ when finite, and the mirror
  of each goal $t$ equal to $m(t)$;
- (Q) every inconsistent state has $k(v) \geq [B + k_m; B]$, and every
  inconsistent goal $k(v) > [B + k_m; B]$.

*Theorem 1 (exactness at quiescence).* Under (H1), at a point where (R)
and (Q) hold, $B$ is the least $g$ of the goals and the cost $d(t^*)$ of a
shortest path, and the parents from a goal $t$ with $g(t) = B$ give a
shortest path.

*Proof.* Let $\hat t$ be a goal with $m(\hat t) = B$. If $\hat t$ were
inconsistent, $k(\hat t) = [B + k_m; B]$ (since $h(\hat t) = 0$), against
(Q); hence $g(\hat t) = B$, and $g(t) \geq m(t) \geq B$ for every goal, so
that $B = \min_t g(t)$. Write $K^* = [B + k_m; B]$ (we take $k_m = 0$, the
general case shifting both sides of every comparison below). (a)
$B \geq d(t^*)$. Suppose $B < d(t^*)$, so that $g(t_1) = B < d(t_1)$ for a
goal $t_1$, consistent as above and not $s_0$. By (R) its parent $u_2$ has
$g(u_2) = g(t_1) - c(u_2, t_1) < d(u_2)$ and $g(u_2) + h(u_2) \leq B$, by
the consistency of $h$. Inductively, if $u_k$ has $g(u_k) < d(u_k)$,
$g(u_k) \leq B - (k - 1) c_{\min}$ and $g(u_k) + h(u_k) \leq B$, it is
consistent (otherwise $k(u_k) \leq [g(u_k) + h(u_k); g(u_k)] < K^*$, the
second component being below $B$, against (Q)), it is not $s_0$, and its
parent satisfies the same with $k + 1$; the values $g(u_k)$ fall by
$c_{\min}$ at each step and are nonnegative, a contradiction. (b)
$B \leq d(t^*)$. Let $s_0 = w_0, \ldots, w_n = t$ be a shortest path to a
goal of $d(t) = d(t^*) < B$ (otherwise there is nothing to prove); then
$g(t) \geq B > d(t)$, so some $w_i$ has $g(w_i) > d(w_i)$, and we take the
least $i$. Either $i = 0$ and $k(s_0) = [h(s_0); 0]$ with
$h(s_0) \leq d(t) < B$, or
$\mathrm{rhs}(w_i) \leq g(w_{i-1}) + c(w_{i-1}, w_i) \leq d(w_i) < g(w_i)$
and the first component of $k(w_i)$ is at most
$d(w_i) + h(w_i) \leq d(t) < B$; in both cases an inconsistent state lies
below $K^*$, against (Q). (c) With $B = d(t^*)$, the induction of (a)
shows that the chain of parents from $t$ goes through consistent states
with $g \leq d$, strictly decreasing, and therefore ends at $s_0$, the only
state with $\mathrm{rhs} = 0$ and no parent; its cost is $g(t) = B$. $\square$

The lemma is static: it does not care how the values were reached.
Everything that follows serves to establish (R) and (Q) at the end of each
episode, and to show that it ends.

*Theorem 2 (no update lost).* Under (H2)-(H6), if (R) holds at the start
of a phase, then (R) holds at the next quiescent point, whatever the
interleaving of the workers, the order in which they choose their states,
the delays of the messages and the values that the loads return within
(H6).

*Proof.* Fix $v \neq s_0$, let $g_f(w)$ be the value of $g(w)$ at the
quiescent point $Q$, and call an *assignment* of $v$ any action that writes
$\mathrm{rhs}(v)$ and $\mathrm{par}(v)$ (the state at the start of the
phase, a Relax that lowers rhs, a rescan). Upper bound: if $g(w)$ was not
written in the phase, (R) gave $\mathrm{rhs}(v) \leq g_f(w) + c(w, v)$ at
the start, and any later rescan loads $g_f(w)$, the only value since the
barrier; if the last write $W$ of $g(w)$ is finite, it is an
overconsistent expansion followed by
$\mathrm{Relax}(v, g_f(w) + c(w, v), w)$ (H4), processed before $Q$, after
which rhs only rises by rescans that happen after $W$ (H5) and therefore
load $W$ (H6); if it is $\infty$ the bound is trivial. Lower bound and
parent: let $A$ be the last assignment of $v$ before $Q$, which set
$\mathrm{rhs}(v) = \gamma + c(u, v)$ and $\mathrm{par}(v) = u$, $\gamma$
a value of $g(u)$ produced by a write $W_\gamma$, and suppose that
$\gamma \neq g_f(u)$, so that a first write $W'$ of $g(u)$ follows
$W_\gamma$, with its message $M'$ to $v$. (i) $M'$ is processed after $A$:
if $A$ is a Relax, it was sent before $W'$ and precedes $M'$ in the
channel (H5); if $A$ is a rescan and $M'$ were processed before it, $W'$
would happen before the load, which then could not return $W_\gamma$ (H6);
if $A$ is the state at the barrier, everything of the phase follows it.
(ii) $M'$ causes an assignment: when it is processed $\mathrm{par}(v) = u$,
and either $W'$ is underconsistent and $M'$ is an Invalidate naming the
parent, or $W'$ is overconsistent and writes a value below $\gamma$ (an
overconsistent expansion lowers $g$), so that $M'$ is a Relax that lowers
rhs. This contradicts the choice of $A$. The mirror of a goal is the last
$m(t)$ sent, by (H4) and (H5). $\square$

The proof uses neither the order of the expansions nor the monotonicity of
$g$, nor that a rescan is atomic: its loads may be spread over time and
see values of different ages. At the barrier nobody writes, and the
changes keep (R) state by state, as in LPA*: a decrease of $c(u, v)$ gives
the new minimum $\min(\mathrm{rhs}(v), g(u) + c(u, v))$, an increase leaves
the minimum and its parent unchanged unless $u$ was the parent, when the
rescan computes it, and a decrease into a state not in a table creates it
when $g(u)$ is finite (the only case in which a state never generated
acquires a predecessor of finite $g$).

*Theorem 3 (termination).* Under (H1), on a finite graph (or a locally
finite one whose states of key below $\sup B$ are finitely many), every
phase reaches a quiescent point in finitely many steps. With `Freeze::all`
each state is expanded underconsistent at most once per episode; with
`Freeze::raised_parent` at most once more than the number of its
predecessors; in the second phase no state is expanded underconsistent.

*Proof.* In the first phase, between two writes of $\infty$ the writes of
$g(u)$ are a strictly decreasing run of finite values, each of them an rhs,
i.e., a value of $g$ at the barrier plus the cost of a walk, or the cost of
a walk from $s_0$; by (H1) only finitely many of them lie below any bound,
so each run is finite. With `Freeze::all`, after its $\infty$ a state is
not lowered again in the phase, and it cannot be raised again, since
$g = \infty \geq \mathrm{rhs}$. With `Freeze::raised_parent`, a raised
state $u$ is lowered only through a parent $w$ not raised in the episode,
as its owner sees it; the store of $g(w)$ releases after the stamp of the
raise of $w$, and the load or the message by which $u$ received that value
acquires it, so that once $u$ has its rhs from a value written after the
raise of $w$, the raise is seen. Between two lowerings of $u$ an
Invalidate from its parent must have raised $\mathrm{rhs}(u)$, which means
that the parent of the first lowering has been raised since, and can no
longer serve; hence $u$ is lowered at most once through each predecessor,
and raised at most once more. The messages are one per write and
successor, the rescans one per Invalidate or raise, and $B$, the least of
the mirrors, rises finitely often, so that the Wakes, at most $p - 1$ per
rise, are finitely many. In the second phase there are only Relax
messages, $g$, rhs and $B$ only decrease, and the runs are finite as
above; an underconsistent state has the key $[g + h + k_m; g]$, which does
not change while it stays underconsistent, and at the end of the first
phase it is not eligible with the value $B_1$ of $B$ (by Theorem 4), while
$B$ only decreases from $B_1$: it never becomes eligible, and if its rhs
falls below $g$ it is overconsistent. $\square$

*Theorem 4 (the stop at quiescence).* At the end of each phase (Q) holds,
except, at the end of the first, for the frozen states, which the second
phase expands.

*Proof.* A worker becomes passive only after finding no eligible state with
the value $B_r$ of $B$ that it loaded, and stays passive until it drains a
message, so that its queue is the same from that load to the quiescent
point $Q$. If the final $B$ is not above $B_r$, (Q) holds for its states,
since eligibility only grows when $B$ rises. If it is, some write of $B$
after the one loaded rose, and its Wake was sent after it (H4); had the
Wake been drained before the load, the load would have returned that write
or a later one (H5, H6), so it was drained after, which makes the worker
active after its last passive moment, a contradiction. In the second phase
$B$ only decreases and the first case always applies. $\square$

By Theorems 1 to 4, every episode returns $B = d(t^*)$ and a shortest path;
an episode stopped by its budget restarts from an empty table, where (R)
holds trivially. The coalescing of the Relax messages to one state in an
outbox by the least value (as long as no Invalidate to the state lies
between them) keeps Theorem 2, the record dropped being dominated by the
one kept, but the outboxes of ADS do not coalesce yet. The rule
`Freeze::raised_parent` and the bound on $\min(g, \mathrm{rhs})$ of the
goals instead of $g$ are ours. With the freeze of every raised state and
the bound on $g$, one worker expanded about 4 times the states of LPA* on
gridworlds of 21 x 21 cells, since a goal raised in the first phase left
$B$ at $\infty$ until the second; with the bound on $\min(g, \mathrm{rhs})$
still about 1.7 times on those of figure 8 of
[[54]](../../README.md#ref-54), since a raised state waits for the second
phase, and about 1.3 times with `Freeze::raised_parent`; this is why one
worker freezes nothing and follows the order of LPA*.

Each episode chooses between the repair and a restart. With the policy
`adaptive` it restarts when the changes that matter (the decreases that
lower an rhs and the increases of a parent edge) exceed $\rho \hat S$,
$\hat S$ the mean of the expansions of the restarts, or when the prediction
of the repair, the sum over the heads of the changes of the expansions per
head that the past repairs cost by the place of the head between the start
and the goal, exceeds $\gamma \hat S$; and a repair is stopped, and
restarted, once it has expanded $\beta \hat S$ states, so that with
$\beta = 1$ (and $\hat S$ exact) an episode expands at most
$\min(R, 2 S)$ states, $R$ those of the repair and $S$ those of the
restart, the deterministic bound of ski rental. We remark that the
prediction is a guess, which the budget makes safe but not good.
