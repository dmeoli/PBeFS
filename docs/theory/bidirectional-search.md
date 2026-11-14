# Bidirectional search

We are concerned here with the front-to-end bidirectional searches of the
library, which look for a path of least cost from a start $s$ to a goal $t$
by a forward search from $s$ and a backward search from $t$ along the edges
reversed, each one guided by a heuristic towards the root of the other. This
page fixes the notation, states the lower bounds on the optimal cost on which
the searches stop, and proves them; then, for each search, it states what the
search guarantees and why, and it ends with the parallel and the external
searches. The pseudocode of each search is on the page of its header, and the
notation of the unidirectional search is that of
[Best-first search and A*](best-first-search.md).

## Notation

A `BidirectionalDomain` gives, besides a `Domain`, the goal state $t$, the
predecessors of a state and a backward heuristic. We write $d( u , v )$ for
the cost of a cheapest path from $u$ to $v$, $C^* = d( s , t )$, and, for a
direction $D \in \{ F , B \}$ with $\bar D$ the other one, $g_D( n )$ for the
cost of the best path found from the root of $D$ to $n$ (from $s$ for $F$, to
$t$ for $B$), $h_F( n )$ for the estimate of $d( n , t )$ and $h_B( n )$ for
that of $d( s , n )$, and $f_D = g_D + h_D$. The heuristic $h_F$ is
admissible if $h_F( n ) \leq d( n , t )$ and consistent if $h_F( u ) \leq
c( u , v ) + h_F( v )$ for every edge $( u , v )$, and symmetrically
$h_B( n ) \leq d( s , n )$ and $h_B( v ) \leq c( u , v ) + h_B( u )$; we
write $\epsilon \geq 0$ for a lower bound on the cost of an edge, which a
`MinCostDomain` gives. We remark that a domain whose moves are reversible at
the same cost, as the sliding-tile puzzle, has its predecessors equal to its
successors, while a grid whose moves cost the cell left, or a directed graph,
does not; the tests check the predecessors of every domain against its
successors.

The two searches share one table of the nodes, each node with $g_F$ and
$g_B$, $\infty$ until the direction reaches it. When a direction sets or
lowers the $g$ of a node whose $g$ in the other direction is finite, the
path through the node costs $g_F + g_B$, and the incumbent $U$, the cost of
the best path found, is lowered to it: this is the immediate solution
detection of Sturtevant and Chen [[94]](../../README.md#ref-94), done by the
same lookup that detects the duplicates. Every search stops when one of its
open lists is empty or when $U \leq \mathrm{LB}$, $\mathrm{LB}$ being a lower
bound on $C^*$ computed from the open lists; then $U = C^*$ if $U < \infty$,
and no path exists otherwise.

## The lower bounds

All the bounds rest on one lemma, which is the argument of Holte et al.
[[61]](../../README.md#ref-61) (section 4 of [[81]](../../README.md#ref-81))
written for any search with immediate solution detection that reopens a node
reached by a cheaper path.

*Lemma 1.* Let $P = ( s = p_0 , \ldots , p_k = t )$ be an optimal path and
suppose $U > C^*$. Then there are $i \leq j$ such that $p_i$ is open in the
forward direction with $g_F( p_i ) = d( s , p_i )$, $p_j$ is open in the
backward direction with $g_B( p_j ) = d( p_j , t )$, and $i < j$.

*Proof.* As in A*, the first node of $P$ not closed in the forward direction
with its optimal $g_F$ is open with it, since its predecessor on $P$ was
expanded with its optimal $g_F$ (by induction from $p_0$, the root), and
likewise backward from $p_k$. If $i \geq j$, the node $p_j$ has been reached
by both directions with its optimal $g$, and when the second of the two set
it the path $P$ was found, so $U \leq C^*$, against the assumption. $\square$

If a search does not open a node that the other direction has closed, the
lemma still holds, provided that this is done only when the closed $g$ is
optimal (see below): the first node of $P$ closed backward and reached
forward with its optimal $g_F$ gives the path $P$ at once. From Lemma 1, with
$n = p_i$ and $m = p_j$, $C^* = g_F( n ) + d( n , m ) + g_B( m )$ and
$d( n , m ) \geq \epsilon$, so the following are lower bounds on $C^*$
whenever $U > C^*$ (hence $\min\{ U , \mathrm{LB} \} \leq C^*$, which is what
the stopping rule needs):

1. $f_{\min F}$ and $f_{\min B}$, the least $f$ of the two open lists, by
   admissibility;
2. $g_{\min F} + g_{\min B} + \epsilon$ (the bound of MM);
3. for $p \in ( 0 , 1 )$, $C_p = \min\{ \min_F \max\{ f_F , g_F / p \} ,
   \min_B \max\{ f_B , g_B / ( 1 - p ) \} \}$: either $g_F( n ) \leq p C^*$
   or $g_B( m ) \leq ( 1 - p ) C^*$, since their sum is at most $C^*$, so one
   of the two nodes has priority at most $C^*$; with $p = 1 / 2$ this is the
   $C$ of MM, and $\max\{ f , 2 g + \epsilon \}$ gives a bound too, since
   $2 g_F( n ) + \epsilon \leq g_F( n ) + g_B( m ) + \epsilon \leq C^*$ when
   $g_F( n ) \leq g_B( m )$;
4. if both heuristics are consistent, $( b_{\min F} + b_{\min B} ) / 2$,
   with $b_D = f_D + g_D - h_{\bar D}$ (the bound of BAE*
   [[89]](../../README.md#ref-89)): by consistency $h_F( n ) \leq d( n , m ) +
   h_F( m )$ and $h_B( m ) \leq d( n , m ) + h_B( n )$, hence
   $b_F( n ) + b_B( m ) = 2 g_F( n ) + 2 g_B( m ) + h_F( n ) - h_F( m ) +
   h_B( m ) - h_B( n ) \leq 2 C^*$;
5. if both heuristics are consistent, $f_{\min F} + d_{\min B}$ and
   $d_{\min F} + f_{\min B}$, with $d_D = g_D - h_{\bar D}$ (two of the
   individual bounds of [[89]](../../README.md#ref-89)): $f_F( n ) + d_B( m ) =
   g_F( n ) + g_B( m ) + h_F( n ) - h_F( m ) \leq C^*$, and symmetrically.

With integer costs every bound may be rounded up, since $C^*$ is an integer.
Each search uses its own bound by default, and the largest of those that
hold under its assumptions with `BidirectionalOptions::all_bounds`. Of
course, all of them are dominated by the bound of pairs of NBS below, which,
however, needs a search over pairs.

## Not opening what the other direction has closed

With consistent heuristics, a node expanded in a direction has its optimal
$g$ there when the priority of the search does not decrease along an edge:
if $Q$ is a cheaper path to $n$ than the one by which $n$ was expanded, every
node of $Q$ has a priority not larger than that of $n$ on $Q$, which is
smaller than that of $n$ on its path, since the priority grows with $g$ at
equal heuristics; so, by induction along $Q$, the prefix of $Q$ is expanded
before $n$, which then is reached through $Q$ first. The priorities $f$,
$\max\{ f , 2 g \}$, $\max\{ f , g / p \}$ and $b$ do not decrease along an
edge under consistency, since $f$ does not and $b$ changes by $2 c +
\Delta h_D - \Delta h_{\bar D} \geq 2 c - c - c = 0$; NBS, DVCBS and GBFHS
expand by $g$ the nodes of $f$ within a bound that only grows. Hence, as in
BS* [[90]](../../README.md#ref-90), every path through a node closed in the
direction $\bar D$ is known when the direction $D$ reaches the node, and the
node is not opened in $D$ (`BidirectionalOptions::consistent`); without
consistency the searches but BAE* still return $C^*$, reopening the nodes.

## The searches one node at a time

`BiBestFirst` expands, in the direction chosen, the node of least priority.
For MM the priority is $\max\{ f_D , 2 g_D \}$ and the bound is $\max\{ C ,
f_{\min F} , f_{\min B} , g_{\min F} + g_{\min B} + \epsilon \}$; since
$C \leq C^*$ whenever $U > C^*$, every node expanded has priority at most
$C^*$, so MM never expands a node with $g_D > C^* / 2$, the property that
names it [[61]](../../README.md#ref-61). The same argument gives, for the
fractional MM of Shaham et al. [[83]](../../README.md#ref-83), no node with
$g_F > p C^*$ or $g_B > ( 1 - p ) C^*$, and the tests check both properties
on the true distances of the 8-puzzle. MMe uses $\max\{ f , 2 g + \epsilon
\}$ [[82]](../../README.md#ref-82), and the bidirectional Dijkstra's algorithm
(the MM0 of [[61]](../../README.md#ref-61)) $g$, with the bound $g_{\min F} +
g_{\min B} + \epsilon$. BAE* expands by $b_D$, alternating the directions as in
its definition [[88]](../../README.md#ref-88) and in PEM-BAE*
[[57]](../../README.md#ref-57), and stops on the bound 4 above, valid only
with consistent heuristics, which is why `bae()` refuses the others.

## NBS and the must-expand graph

For $u$ open forward and $v$ open backward, $\mathrm{lb}( u , v ) = \max\{
f_F( u ) , f_B( v ) , g_F( u ) + g_B( v ) + \epsilon \}$ bounds the cost of
the paths through $u$ and then $v$. Theorem 1 of Eckerle et al.
[[84]](../../README.md#ref-84) states that, with consistent heuristics, every
admissible front-to-end search expands $u$ or $v$ for every pair of optimal
paths to them with $\mathrm{lb} < C^*$; these pairs are the edges of the
must-expand graph $G_{MX}$, and the size VC of its minimum vertex cover bounds
the expansions from below. NBS [[85]](../../README.md#ref-85) expands both
ends of a pair of least lb, i.e., the greedy cover of Papadimitriou and
Steiglitz, so it expands at most $2 \mathrm{VC}$ states to cover the pairs
(theorem 3 of [[85]](../../README.md#ref-85)), and no admissible front-to-end
search does better in the worst case (theorem 5); the expansions of the nodes
with $\mathrm{lb} = C^*$ are outside the bound. The least lb is found by the
queues waiting (by $f$) and ready (by $g$) of section 7 of
[[85]](../../README.md#ref-85), with a lower bound $C_{lb}$ raised only when
no pair reaches it, so that $C_{lb}$ is the least lb of the pairs when a pair
is returned and the stopping rule $U \leq C_{lb}$ is that of Lemma 1. The
eager variant makes ready at once the nodes of $f = C_{lb}$, which covers the
pairs of $\mathrm{lb} \leq C^*$, those needed to find every optimal path
[[86]](../../README.md#ref-86).

## DVCBS and the dynamic vertex cover

DVCBS [[86]](../../README.md#ref-86) builds, at each step, the graph of the
clusters of the ready nodes (one per value of $g$ in each direction), a
forward cluster $a$ and a backward one $b$ being joined when $a + b +
\epsilon \leq \mathrm{LB}$, and expands one cluster of a minimum vertex cover
weighted by the sizes. The cover is found in linear time by the following
fact, which is the structure of the covers of $G_{MX}$ shown by Shaham et al.
[[83]](../../README.md#ref-83).

*Lemma 2.* Some minimum vertex cover is a prefix $a_0 < \cdots < a_{i-1}$ of
the forward clusters together with the prefix of the backward clusters joined
to $a_i$.

*Proof.* The cluster $b$ is joined to $a$ if and only if $b \leq
\mathrm{LB} - \epsilon - a$, so the neighbours of $a_i$ contain those of every
$a_k$ with $k > i$. Given a cover, let $a_i$ be its forward cluster of least
$g$ missing; all the neighbours of $a_i$ are in it, hence those of every later
cluster, and dropping the forward clusters after $a_i$ and the backward ones
not joined to $a_i$ leaves a cover of no larger weight; adding the missing
forward clusters before $a_i$ instead of those backward neighbours is
accounted for by the other values of $i$. $\square$

The pass over $i = 0 , \ldots , m$ keeps the backward prefix by a pointer that
only moves down, and the cluster expanded is the forward or the backward one
of least $g$ that lies in a minimum cover, the smaller of the two, as in the
variant reported by Shperberg et al. DVCBS has no bound of the kind of NBS
(section 5.2 of [[86]](../../README.md#ref-86)).

## GBFHS

GBFHS [[87]](../../README.md#ref-87) raises a bound $f_{Lim}$ by one unit at a
time and, for each value, expands in each direction every node with
$f_D \leq f_{Lim}$ and $g_D < g_{LimD}$, where $g_{LimF} + g_{LimB} =
f_{Lim} - \epsilon + 1$ with integer costs.

*Lemma 3.* After the level $f_{Lim}$, every path of cost $C \leq f_{Lim}$
whose nodes have their optimal $g$ has been found.

*Proof.* Let $u$ follow the last node of the path expanded forward and $v$
precede the first one expanded backward. If $u$ precedes $v$, then
$g_F( u ) \geq g_{LimF}$, $g_B( v ) \geq g_{LimB}$ and $C \geq g_F( u ) +
\epsilon + g_B( v ) \geq f_{Lim} + 1$, a contradiction; so a node of the path
is reached by both directions, and the path is found. $\square$

With consistent heuristics the nodes are expanded with their optimal $g$, since
they are expanded by increasing $g$ among those of $f$ within a bound that only
grows and $f$ does not decrease along an optimal path; hence $f_{Lim}$ stays a
lower bound on $C^*$ and GBFHS expands no node of $f_D > C^*$. The split of
each new unit between the two limits is a parameter: in turn, to the direction
with fewer open nodes (Pohl's criterion [[91]](../../README.md#ref-91)), or by
a fraction, as fMM. Without heuristics GBFHS is the bidirectional breadth-first
search.

## The two directions in parallel

`ParallelBidirectional` runs each direction on a worker, as PNBA* of Rios and
Chaimowicz [[93]](../../README.md#ref-93) does for NBA* of Pijls and Post
[[92]](../../README.md#ref-92). Two facts make it exact. First, the worker of
$D$ writes $g_D$ of a node and then reads $g_{\bar D}$, and the other worker
does the converse, with sequentially consistent accesses: in any execution
one of the two reads follows both writes, so a node reached by both
directions lowers $U$. Second, the bound is computed from the least keys that
each worker publishes after each expansion; with consistent heuristics the
least key of the nodes open or being expanded in a direction never
decreases, because a node popped has a key not smaller than the least one and
its children keys not smaller than its own, so a value read late is still a
lower bound, and the stopping rule holds as in the sequential search. NBA*
rejects a node $x$, closing it unexpanded, when $f_D( x ) \geq U$ or
$g_D( x ) + F_{\bar D} - h_{\bar D}( x ) \geq U$, with $F_{\bar D}$ the least
$f$ published by the other direction: by consistency, a path through $x$ and
a node $y$ open in $\bar D$ costs at least $g_D( x ) + h_{\bar D}( y ) -
h_{\bar D}( x ) + g_{\bar D}( y ) \geq g_D( x ) + F_{\bar D} -
h_{\bar D}( x )$ (section 3 of [[92]](../../README.md#ref-92)).

## The external searches

`PEMSearch` is the framework PEM-BiHS of Siag et al.
[[57]](../../README.md#ref-57). The nodes of a direction are kept in buckets of
equal $( g , h_F , h_B )$, split into partitions by the hash of the state, so
that the copies of a state, in either direction, lie in the same partition of
buckets of the same heuristics. A bucket is expanded with its partitions in
parallel: the copies within a partition are dropped by a hash table, those
already closed by a scan of the closed partitions of the same heuristics and
smaller $g$ (the bucket of least $g$ among those of least priority is the next
one, so that, with consistent heuristics and positive costs, a bucket gets no
node after its expansion), and the meetings are found by a scan of the
partitions of the same heuristics of the other direction, the delayed solution
detection of Sturtevant and Chen [[94]](../../README.md#ref-94). Since a
meeting is found only when one of its two copies is expanded, a pair of open
copies on an optimal path may be unknown, and the bound $g_{\min F} + g_{\min
B} + \epsilon$ is not valid (section 3.3 of [[94]](../../README.md#ref-94));
the bound without $\epsilon$ is, by the same argument as Lemma 1 with $i \leq
j$, and so are the bounds 1, 3 (with $\max\{ f , 2 g \}$), 4 and 5 computed
from the keys of the open buckets.

`IDDAStar` is A*-IDD of Lin and Fukunaga [[95]](../../README.md#ref-95). Its
invariant is that every closed state is either in the buffer of its partition
or in the file, at a position held, plus one, by a place of the table that the
linear probing from its hash reaches before an empty place, and that the block
of a position holds only states of one partition; a probe that meets a block of
another partition skips it without reading the disk, and the duplicates are
dropped when a node is popped, before its expansion, so that with a consistent
heuristic each state is expanded once, as by A*.
