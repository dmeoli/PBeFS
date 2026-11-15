# Symbolic search

We are concerned here with the searches of `SymbolicSearch.h`, which store
sets of states instead of states, each set represented by the reduced ordered
Binary Decision Diagram (BDD) [[126]](../../README.md#ref-126) of its
characteristic function. A set is expanded as a whole, by a relational product
with the transition relation, as in symbolic model checking
[[127]](../../README.md#ref-127). This page fixes the encoding of a
`SymbolicDomain` and proves that the relations and the sets of the heuristic
built by `SymbolicModel` are the intended ones; then it proves the optimality
of BDDA* and of the stopping rule of the bidirectional search, and it ends
with the parallel computation of the images, the choice of the library of
BDDs and the limits of the implementation. The notation of the explicit
searches is that of [Best-first search and A*](best-first-search.md) and of
[Bidirectional search](bidirectional-search.md).

## The encoding

A `SymbolicDomain` has the variables $v = 0 , \ldots , n - 1$, the variable
$v$ with the values $D_v = \{ 0 , \ldots , | D_v | - 1 \}$. A state $x \in D_0
\times \cdots \times D_{n-1}$ is encoded by writing each $x_v$ in binary on
$b_v = \lceil \log_2 | D_v | \rceil$ bits (the most significant first) and by
concatenating the codes in the order of the variables. The bit $k$ of the
state before a move is the variable $2 k$ of the library, and the same bit
after the move is the variable $2 k + 1$; hence, a set of states is a BDD $S(
x )$ on the even variables and a relation is a BDD $T( x , x' )$ on both. The
codes that are not values (when $| D_v |$ is not a power of two) are excluded
by the set $V$ of the valid states, those with $x_v < | D_v |$ for all $v$.
The initial state is the minterm of its code, and the goal is $G = V \cap
\bigcup_p \{ x : x_u = d \text{ for } ( u , d ) \in p \}$ over the partial
assignments $p$ given by the domain (one for the puzzles, the goal of a SAS+
task [[133]](../../README.md#ref-133) for `SymbolicSAS`).

An operator $o$ has a cost $c( o ) \geq 0$, a precondition $\mathrm{pre}_o$ (a
set of facts $x_u = d$) and effects, each of which sets a variable $w$ either
to a constant ($x'_w = d$) or to the value that a variable $u$ had before the
move ($x'_w = x_u$). A variable is written by one effect at most, and the
variables that are not written keep their value. With $R_o$ the variables that
$o$ reads or writes and $W_o \subseteq R_o$ those it writes, its relation is

$$
T_o( x , x' ) = \mathrm{pre}_o( x ) \wedge \bigwedge_{e \text{ of } o} E_e(
x , x' ) \wedge \bigwedge_{v \in R_o \setminus W_o} ( x'_v = x_v ) ,
$$

where $E_e$ is $x'_w = d$ or $x'_w = x_u$, the latter written bit by bit as
$x'_w \leftrightarrow x_u$ (with the high bits of $x'_w$ set to zero if $b_u <
b_w$). The image of $S$ under $T_o$ on the variables $R_o$ is

$$
\mathrm{img}( S , T_o ) = \{ x' : \exists x \in S , \; T_o( x , x' ) , \;
x'_v = x_v \text{ for } v \notin R_o \} ,
$$

which the relational product of the library computes when it is given the
variables of $R_o$ (in both copies), leaving the others unchanged; the
preimage $\mathrm{pre}( S , T_o )$ is the same with $x$ and $x'$ exchanged.

*Lemma 1.* For any set $S$ of valid states, $\mathrm{img}( S , T_o )$ is the
set of the states reached from $S$ by applying $o$, and $\mathrm{pre}( S , T_o
)$ the set of the states from which applying $o$ reaches $S$.

*Proof.* A pair $( x , x' )$ is in $T_o$, with $x'_v = x_v$ outside $R_o$, if
and only if $x$ satisfies $\mathrm{pre}_o$, each written variable of $x'$ has
the value given by its effect (computed on $x$, since each $E_e$ reads the
unprimed copy), and each other variable of $x'$ equals that of $x$ (inside
$R_o$ by the frame of $T_o$, outside by the relational product). This is the
definition of applying $o$ to $x$, and both sets follow by projection.
$\square$

The operators of equal cost are merged, in the order in which the domain gives
them, into partitions $T_P = \bigvee_{o \in P} T_o^{R_P}$ on the union $R_P$ of
their variables ($T_o^{R_P}$ is $T_o$ with the frames $x'_v = x_v$ of $R_P
\setminus R_o$ added). A partition is closed when the next operator would bring
its diagram above `max_relation_nodes` nodes, and with 0 each partition holds
one operator. This is the disjunctive partitioning bounded by size of the
symbolic planners [[130]](../../README.md#ref-130): the image of a disjunction
is the union of the images, while the diagram of the disjunction of all the
operators is usually much larger than those of the parts together. By Lemma 1
applied to each $T_o^{R_P}$ (whose added frames repeat what the relational
product does anyway), $\mathrm{img}( S , T_P ) = \bigcup_{o \in P}
\mathrm{img}( S , T_o )$; hence, the successors of cost $c$ of $S$ are the
union of its images under the partitions of cost $c$.

## The heuristic as a partition of the states

The searches do not evaluate a heuristic on a state; instead, they split a set
by the value of the heuristic, i.e., they intersect it with the sets $H_k = \{
x \in V : h( x ) = k \}$. A `SymbolicHeuristicDomain` gives $h$ as a sum of
terms $h( x ) = \sum_{t = 1}^m h_t( x_{U_t} )$, each a table on the values of a
few variables $U_t$ (one cell for the Manhattan distance of the sliding-tile
puzzle, two adjacent positions for the gap heuristic of the pancakes
[[60]](../../README.md#ref-60)), and a term equal to $\infty$ marks a dead end.
`SymbolicModel` computes

$$
P_0( 0 ) = \mathbf{1} , \qquad P_t( a ) = \bigvee_{b} \big( P_{t-1}( a - b
) \wedge B_t( b ) \big) , \qquad B_t( b ) = \{ x : h_t( x_{U_t} ) = b \} ,
$$

and takes $H_k = P_m( k )$, as Edelkamp and Reffel did for the Manhattan
distance in their BDDA* [[128]](../../README.md#ref-128).

*Lemma 2.* $P_m( k )$ is the set of the states whose terms are finite and add
up to $k$; in particular, the sets $H_k$ are disjoint, and their union holds
all the states that are not marked as dead ends.

*Proof.* By induction on $t$, $P_t( a )$ is the set of the states whose first
$t$ terms are finite and add up to $a$. This holds for $t = 0$, and a state
is in $P_t( a )$ if and only if its term $t$ has some finite value $b$ (only
the finite values have a set $B_t( b )$) and its first $t - 1$ terms add up
to $a - b$. $\square$

Note that the diagrams stay rather small when the terms are on disjoint
variables, ordered as the variables are: a level of the diagram of $P_t( a )$
then has at most as many nodes as there are partial sums, i.e., $O( h_{\max}
)$, and at most $O( h_{\max} | D_v | )$ when consecutive terms share a variable
$v$ (as the pairs of the pancakes do). We remark that a heuristic given as a
table on many variables (a large pattern database, say) would need a different
construction, which is not implemented (cf. the last section).

Since $h$ is never evaluated on a state, its consistency cannot be checked on
the states expanded, as an explicit search may do; it can, however, be
checked on all the valid states at once, by images.

*Lemma 3.* `SymbolicModel::consistent()` returns true if and only if (i) no
goal has $h > 0$ and (ii) for each partition $T_P$ of cost $c$ and each $k$,
$\mathrm{img}( H_k , T_P ) \cap H_j = \emptyset$ for all $j < k - c$. Under
(i) and (ii), $h$ is consistent along the moves between states that are not
dead ends, and it is admissible.

*Proof.* The first claim is the code. By Lemmas 1 and 2, (ii) says that no
move of cost $c$ leads from a state with $h = k$ to one with $h < k - c$, i.e.,
$h( x ) \leq c + h( y )$ for any move $( x , y )$ between states with finite
$h$. Summing along a cheapest path from $x$ to a goal $y$ (whose states reach
a goal, and hence are not dead ends) gives $h( x ) \leq d( x , y ) + h( y ) =
d( x , y )$ by (i). $\square$

Also, the backward heuristic (towards the initial state, from the terms of a
`SymbolicBackwardHeuristicDomain`) is checked by preimages in the same way,
with the initial state in the place of the goals. The test `sym-astar` checks
that the Manhattan distance and the gap heuristic are consistent in both
directions, and that twice the Manhattan distance is found inconsistent.

## BDDA*

A direction keeps an open list of sets $O_{f,g}$ (the states reached with
cost $g$ and with $f = g + h$), the closed set $C$ of the states expanded, the
sets $R_g$ of the states reached with cost $g$, and the sequence of the sets
expanded with their $g$. At the start the root (the initial state $s$, or the
goals backward) is split into the sets $O_{k,0}$ by the sets $H_k$ of the
direction. A step takes the nonempty set $O_{f,g}$ of least $f$ and, among
those, of least $g$, and sets $F = O_{f,g} \setminus C$; then, while $F \neq
\emptyset$, it (i) adds $F$ to $C$, (ii) stops with cost $g$ if $F$ meets the
target, (iii) computes $Y_c = \mathrm{img}( F , T_c ) \setminus C$ for each
cost $c$ ($T_c$ the partitions of cost $c$), and (iv) adds $Y_c \cap H_k$
to $O_{g+c+k,g+c}$ for each $k$, except $Y_0 \cap H_{f-g}$, which becomes the
next $F$. These are the buckets of $g$ and $h$ of BDDA*
[[128]](../../README.md#ref-128), close to the SetA* of Jensen, Veloso and
Bryant [[129]](../../README.md#ref-129). With $h = 0$ (the option `heuristic
= false`, or a domain without terms) the search is the symbolic Dijkstra's
algorithm of the planners with costs [[134]](../../README.md#ref-134), and
with all costs set to 1 (`unit_costs`, `symbolic_bfs()`) it is the
breadth-first search, whose layers are the sets $O_{g,g}$. The loop on the
moves of cost 0 is the breadth-first search within a bucket that such moves
require: by consistency they cannot lower $h$, and the successors with the
same $h$ belong to the bucket being expanded.

*Proposition 4.* Let $h$ be consistent, zero on the target and infinite only
on states from which no target is reachable. Then a state is expanded at most
once, a state $x$ from which a target is reachable with $g = d( s , x )$, and
the search stops with $g = C^*$ if a target is reachable; otherwise, the open
list empties.

*Proof.* A state is expanded at most once, since it enters $C$ when it is
expanded and any set taken from the open list is reduced by $C$. Suppose that
a state from which a target is reachable is expanded with $g$ larger than its
distance, let $x$ be the first one in time, expanded at $( f , g )$ with $g >
d( s , x )$, and let $P$ be a cheapest path from $s$ to $x$. When $x$ is
expanded, let $y$ be the first state of $P$ not in $C$ (a target is reachable
from $y$ through $x$, and hence $h( y )$ is finite). The predecessor of $y$ on
$P$ was expanded before $x$, and therefore with its distance; hence, $y$ was
inserted in $O_{f_y,g_y}$ with $g_y = d( s , y )$ and, by consistency along
the rest of $P$, $f_y = g_y + h( y ) \leq d( s , x ) + h( x ) < f$. Since $y
\notin C$, that set has not been expanded yet, which contradicts the choice of
the set of least $f$ (the successors of the set being expanded have $f' \geq
f$ by consistency). The same argument, applied to a cheapest path to a target
$t$ with $d( s , t ) < g$ when the search stops at $g$, gives a set of the
open list with $f_y \leq d( s , t ) + h( t ) = d( s , t ) < g \leq f$, again a
contradiction; hence the search stops with $g = C^*$. It ends, since $C$ grows
at each step and the states are finite. $\square$

With an inconsistent heuristic the search still returns the cost of a path
(the sets record the costs of real paths), but not necessarily an optimal one,
since a state may be closed with a $g$ larger than its distance; this is why
`consistent()` is there, to check a heuristic before it is used. The states
marked as dead ends are in no $H_k$, and they are dropped when they are
reached.

*The path.* The search records each set $F$ expanded with its $g$, in order.
If $F_i$ (with cost $g_i$) holds the state $t$ that meets the target, a
predecessor of $t$ is looked for in the sets $F_j$, $j < i$, from the last
one, with $g_j = g_i - c$ for some cost $c$ and $\mathrm{pre}( \{ t \} , T_c )
\cap F_j \neq \emptyset$; one state of that intersection is taken, and the
walk repeats from it until it reaches the root. The walk does not fail, since
each state of $F_i$ entered its bucket as a successor of a state of a set
expanded before, with $g$ lower by the cost of the move. It costs one preimage
of a single state per partition and per cost tried (the counter
`path_images`). Backward, the walk takes images, and it gives the path from
the initial state to a goal directly.

## The bidirectional search

Both directions run as above, the forward one by images from the initial state
with $h_F$, the backward one by preimages from the goals with $h_B$. Each set
$X$ inserted in the open list of a direction with cost $g$ (the roots and the
layers of cost 0 included) is intersected with the sets $R_{g'}$ of the other
direction by increasing $g'$, and the first nonempty intersection lowers the
incumbent $U$ to $g + g'$ if it is smaller. This is the immediate solution
detection of the explicit searches [[94]](../../README.md#ref-94), done on
sets. A successor with $f \geq U$ is not inserted, since no path through it
costs less than $U$. Let $f_D$ and $g_D$ be the least $f$ and the least $g$ of
the nonempty sets in the open list of the direction $D$, and $\epsilon$ the
least cost of an operator; the search stops when an open list is empty or
when

$$
U \leq \mathrm{LB} = \max \{ f_F , f_B , g_F + g_B + \epsilon \}
$$

(more precisely, the largest value of $\mathrm{LB}$ seen so far).

*Proposition 5.* Let $h_F$ and $h_B$ be consistent and admissible. When the
search stops, $U = C^*$ if a path exists, and $U = \infty$ otherwise.

*Proof.* $U$ is always the cost of a path, since $g$ and $g'$ are costs of
paths from $s$ to the state and from the state to a goal. Suppose that $U >
C^*$ at some time, and let $P = ( s = p_0 , \ldots , p_k )$ be an optimal path
to a goal. As in the proof of Proposition 4, the first state $p_i$ of $P$ not
closed forward is in the forward open list with $g = d( s , p_i )$ (it is not
pruned, since its $f$ is at most $C^* < U$), and the last state $p_j$ not
closed backward is in the backward open list with $g = d( p_j , p_k )$. If $i
\geq j$, the state $p_j$ has been inserted by both directions with its optimal
$g$, and the second of the two insertions met the set $R_{g'}$ of the other
direction that holds it, giving $U \leq C^*$, a contradiction. Hence $i < j$,
and $f_F \leq d( s , p_i ) + h_F( p_i ) \leq C^*$, $f_B \leq C^*$ likewise,
and $g_F + g_B + \epsilon \leq C^*$, since the segment from $p_i$ to $p_j$ has
at least one move. Thus $\mathrm{LB} \leq C^* < U$ whenever $U > C^*$, and
also for the largest $\mathrm{LB}$ seen, since $U$ does not increase; the stop
on $U \leq \mathrm{LB}$ therefore gives $U = C^*$. If the forward open list
empties, all the states reachable from $s$ that are not dead ends have been
closed forward, among them the goal $p_k$ of $P$, which was inserted with $g =
C^*$ and met the set $R_0$ of the backward direction; hence $U = C^*$, and
symmetrically for the backward list. With no path, $U$ stays $\infty$.
$\square$

This is the argument of Lemma 1 of the page on the bidirectional search, with
sets in the place of nodes. Without the two bounds on $f$, the bound $g_F +
g_B + \epsilon$ is the stopping rule of the bidirectional Dijkstra's
algorithm, and of the bidirectional symbolic search of SymBA*
[[130]](../../README.md#ref-130), whose heuristics are abstractions computed
during the search instead of the terms given by the domain.

*The direction of a step.* Any interleaving of the two directions satisfies
Proposition 5, and the choice is therefore free. A step expands the direction
whose next set is estimated cheaper to expand; the estimate is the nodes of
the set times the ratio between the nodes of the image and those of the set at
the last step of that direction. SymBA* chooses the direction by the estimated
time of its next step (as described in [[130]](../../README.md#ref-130)); we
use counts of nodes instead, so that the search and its counters do not
depend on the load of the machine or on the workers. The path is the forward
walk from a state of the intersection that set $U$ (with its $g$ forward),
reversed, followed by the backward walk from the same state.

## The parallel images

In Sylvan [[131]](../../README.md#ref-131) the operations are parallel
inside: each recursive call of a relational product, a conjunction or a
quantification spawns its two cofactors as tasks of the work-stealing
scheduler Lace, so that one product of a large set can use all the workers.
On top of that, `SymbolicImages` chooses how the images of the partitions of a
set are computed: (i) `sequential`, one product after the other; (ii)
`lace`, all of them in one task of Lace, which spawns the products of the two
halves of the partitions recursively and then joins them by a balanced tree of
disjunctions (`images()` and `disjoin()` in `Bdd.h`), so that the scheduler
balances products of rather different sizes; (iii) `executor`, the partitions
dealt round-robin to the workers of an executor of ADS (the threads of the
library, OpenMP or FastFlow), each of which submits its products to Lace,
which accepts them from several threads at once.

*Lemma 6.* The costs, the paths and the counters of `SymbolicStats` other than
`images_by` do not depend on the choice of `SymbolicImages`, on the workers of
the executor or on the workers of Lace.

*Proof.* A reduced ordered BDD is canonical for its function and its order of
the variables [[126]](../../README.md#ref-126), and the order is fixed by the
encoding. Each set computed by the search is a function of the sets before it,
whatever way its parts were computed; hence so are its nodes, its states, the
choice of the direction (which depends on nodes only) and the state taken from
an intersection by the path (the first in the order of the variables, by
`sylvan_sat_one`). The products are counted one per partition and per set,
whoever computes them. $\square$

Indeed, the test `sym-par-*` checks Lemma 6 on the 8-puzzle, the pancakes and
the 8-puzzle of heavy tiles with moves of cost 0, in the three directions, with
1, 2 and 4 workers of Lace and of each executor. The lemma also makes the
counters a fair measure of the work of these searches, since the nodes and the
products are the same on any machine while the time of a product is not. The
speedup itself is not measured here.

One point of the implementation needs care. The garbage collection of Sylvan
drops the nodes that no registered handle reaches, and a result received by a
thread outside Lace is not registered until it is stored in a `Bdd`. With one
thread driving the search this cannot matter, since a collection runs only
inside a product; with the executor, however, several threads submit products
at once, and `image_into()` therefore stores its result in the registered
handle of its output inside the task of Lace, before the task ends.

## The library of BDDs

PBeFS interfaces an existing library instead of implementing the diagrams,
and we compared the two that fit a C++ library distributed on three systems.

| | Sylvan 1.12.6 | CUDD 3.0 |
|---|---|---|
| licence | Apache-2.0 | BSD-3-Clause |
| build | CMake (3.23 or later), downloads Lace | autotools only |
| Linux, macOS | in its CI (gcc, clang) | configure and make |
| Windows | in its CI: MSYS2 (3 environments) and MSVC | no supported build |
| vcpkg, Homebrew, MSYS2 packages | none | none |
| parallel operations | work stealing on Lace | none, one manager per thread |

We took the table from the repositories of the two libraries (the workflows of
Sylvan, the files at the root of the mirrors of CUDD) and from the indices of
vcpkg, Homebrew and MSYS2, as of November 2026. Sylvan builds on the three
systems, with CMake as PBeFS does, and its operations are parallel, which is
the side of symbolic search that concerns a library of parallel search; we
therefore use Sylvan, while CUDD would have needed a port to
Windows and would have computed each product on one core. Sylvan is optional:
CMake finds an installed one with `find_package(sylvan)`, or else downloads
the commit pinned in `CMakeLists.txt` (unless `PBEFS_FETCH_SYLVAN` is off),
and the makefiles take it from `SYLVAN_ROOT`. Without it, the library
`PBeFS::Symbolic` and the cases `sym-*` of the tester are left out, and the
case `registry` accounts for them.

## Limits

Of course, the encodings by cells of the sliding-tile puzzle and by positions
of the pancakes are the natural ones, but not necessarily the best. The flip of
$k$ pancakes reverses $k$ positions, so that its relation keeps $\lfloor k / 2
\rfloor$ values while it reads the positions in order, and its diagram grows
roughly as $N^{k/2}$; the stacks of the tests have at most 9 pancakes, and of
the 15-puzzle the tests take the easiest instances of Korf. The operators have
no conditional effects, the SAS+ tasks have no axioms, and the heuristic is a
sum of terms on few variables; the symbolic pattern databases and the perimeter
abstractions of SymBA* [[130]](../../README.md#ref-130) and the search for the
$k$ best plans of symk [[132]](../../README.md#ref-132) are not implemented.
The partitions follow the order of the operators, and the order of the
variables is that of the domain, without reordering. Finally, a search that
fills more than `max_filled` of the largest table of nodes after a collection
stops with `Status::node_limit`, with the least $f$ of its open list as its
bound; a single product that fills the table ends the process, as Sylvan does.