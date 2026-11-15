# Symbolic search

We are concerned here with the searches of `SymbolicSearch.h`, which do not
store the states one by one but sets of states, each represented by a reduced
ordered Binary Decision Diagram (BDD) [[126]](../../README.md#ref-126) of its
characteristic function, and which expand a whole set at once by a relational
product with the transition relation, as in symbolic model checking
[[127]](../../README.md#ref-127). This page fixes the encoding of a
`SymbolicDomain`, proves that the relations and the sets of the heuristic
built by `SymbolicModel` are the intended ones, proves the optimality of
BDDA* and of the stopping rule of the bidirectional search, and ends with the
parallel computation of the images, the choice of the library of BDDs and
the limits of the implementation. The notation of the explicit searches is
that of [Best-first search and A*](best-first-search.md) and of
[Bidirectional search](bidirectional-search.md).

## The encoding

A `SymbolicDomain` has the variables $v = 0 , \ldots , n - 1$, the variable
$v$ with the values $D_v = \{ 0 , \ldots , | D_v | - 1 \}$; a state is a
vector $x \in D_0 \times \cdots \times D_{n-1}$, encoded by writing each $x_v$
in binary on $b_v = \lceil \log_2 | D_v | \rceil$ bits, the most significant
first, and concatenating the codes in the order of the variables. The bit $k$
of the state before a move is the variable $2 k$ of the library and the same
bit after the move the variable $2 k + 1$, so that a set of states $S$ is a
BDD $S( x )$ on the even variables and a relation a BDD $T( x , x' )$ on
both; the codes of $x_v$ that are not values (when $| D_v |$ is not a power of
two) are excluded by the set $V = \{ x : x_v < | D_v | \text{ for every } v
\}$ of the valid states. The initial state is the minterm of its code, and
the goal is $G = V \cap \bigcup_p \{ x : x_u = d \text{ for } ( u , d ) \in p
\}$ over the partial assignments $p$ that the domain gives (one for the
puzzles, the goal of a SAS+ task [[133]](../../README.md#ref-133) for
`SymbolicSAS`).

An operator $o$ has a cost $c( o ) \geq 0$, a precondition $\mathrm{pre}_o$
(a set of facts $x_u = d$) and effects, each of which sets a variable $w$
either to a constant ($x'_w = d$) or to the value that a variable $u$ had
before the move ($x'_w = x_u$); a variable is written by at most one effect,
and the variables not written keep their value. Let $R_o$ be the variables
that $o$ reads or writes and $W_o \subseteq R_o$ those it writes; its
relation is

$$
T_o( x , x' ) = \mathrm{pre}_o( x ) \wedge \bigwedge_{e \text{ of } o} E_e(
x , x' ) \wedge \bigwedge_{v \in R_o \setminus W_o} ( x'_v = x_v ) ,
$$

where $E_e$ is $x'_w = d$ or $x'_w = x_u$, the latter being the conjunction
over the bits of $x'_w \leftrightarrow x_u$ (the high bits of $x'_w$ set to
zero if $b_u < b_w$). The image of $S$ under $T_o$ on the variables $R_o$ is

$$
\mathrm{img}( S , T_o ) = \{ x' : \exists x \in S , \; T_o( x , x' ) , \;
x'_v = x_v \text{ for } v \notin R_o \} ,
$$

which is what the relational product of the library computes when it is given
the set of the variables of $R_o$ (in both copies) and leaves the others
unchanged; the preimage $\mathrm{pre}( S , T_o )$ is the same with the roles
of $x$ and $x'$ exchanged.

*Lemma 1.* For every set $S$ of valid states, $\mathrm{img}( S , T_o )$ is the
set of the states reached from a state of $S$ by applying $o$, and
$\mathrm{pre}( S , T_o )$ is the set of the states from which applying $o$
reaches a state of $S$.

*Proof.* A pair $( x , x' )$ is in $T_o$ with $x'_v = x_v$ outside $R_o$ if
and only if $x$ satisfies $\mathrm{pre}_o$, every written variable of $x'$
has the value its effect gives (computed on $x$, since every $E_e$ reads
the unprimed copy), and every other variable of $x'$ equals that of $x$,
inside $R_o$ by the frame of $T_o$ and outside it by the semantics of the
relational product; this is the definition of applying $o$ to $x$. Both sets
follow by projection. $\square$

The operators of equal cost are merged, in the order in which the domain
gives them, into partitions $T_P = \bigvee_{o \in P} T_o^{R_P}$ on the union
$R_P$ of their variables, $T_o^{R_P}$ being $T_o$ conjoined with the frames
$x'_v = x_v$ of the variables of $R_P \setminus R_o$; a partition is closed
when adding the next operator would bring its diagram above
`max_relation_nodes` nodes, and with 0 each operator is a partition of its
own. This is the disjunctive partitioning bounded by size used by symbolic
planners [[130]](../../README.md#ref-130), since the image of a disjunction
is the union of the images, while the BDD of the disjunction of all the
operators may be much larger than the sum of those of the parts. By Lemma 1
applied to each $T_o^{R_P}$ (whose frame on $R_P \setminus R_o$ adds nothing
to that of the relational product), $\mathrm{img}( S , T_P ) = \bigcup_{o
\in P} \mathrm{img}( S , T_o )$, and the successors of cost $c$ of $S$ are
the union of the images under the partitions of cost $c$.

## The heuristic as a partition of the states

The searches do not evaluate a heuristic on a state: they split a set of
states by the value of the heuristic, i.e., they intersect it with the sets
$H_k = \{ x \in V : h( x ) = k \}$. A `SymbolicHeuristicDomain` gives $h$ as
a sum of terms $h( x ) = \sum_{t = 1}^m h_t( x_{U_t} )$, each a table on the
values of a few variables $U_t$ (one cell for the Manhattan distance of the
sliding-tile puzzle, two adjacent positions for the gap heuristic of the
pancakes [[60]](../../README.md#ref-60)), a term equal to $\infty$ marking a
dead end. `SymbolicModel` computes

$$
P_0( 0 ) = \mathbf{1} , \qquad P_t( a ) = \bigvee_{b} \big( P_{t-1}( a - b
) \wedge B_t( b ) \big) , \qquad B_t( b ) = \{ x : h_t( x_{U_t} ) = b \} ,
$$

and takes $H_k = P_m( k )$, as Edelkamp and Reffel did for the Manhattan
distance in their BDDA* [[128]](../../README.md#ref-128).

*Lemma 2.* $P_m( k )$ is the set of the states whose terms are all finite and
add up to $k$; in particular the sets $H_k$ are pairwise disjoint and their
union holds every state that is not a dead end.

*Proof.* By induction on $t$, $P_t( a )$ is the set of the states whose
first $t$ terms are finite and add up to $a$: true for $t = 0$, and a state
is in $P_t( a )$ if and only if its term $t$ takes some finite value $b$ (it
is in $B_t( b )$, which holds only the finite values) and its first $t - 1$
terms add up to $a - b$. $\square$

The size of these diagrams is what makes the encoding by terms practical:
when the terms are on disjoint variables, ordered as the variables are, a
level of the diagram of $P_t( a )$ distinguishes at most as many nodes as
there are partial sums, i.e., $O( h_{\max} )$, and when two consecutive terms
share a variable (as the pairs of the pancakes do) at most $O( h_{\max} | D_v
| )$. We remark that a heuristic given as a table on all the states (a
pattern database over many variables) would need a different construction,
which is not implemented [see the last section].

Since the searches never evaluate $h$ on a state, its consistency cannot be
checked on the states expanded, as an explicit search may do; it can,
however, be checked on all the valid states at once, by images.

*Lemma 3.* `SymbolicModel::consistent()` returns true if and only if (i) no
goal state has $h > 0$, and (ii) for every partition $T_P$ of cost $c$ and
every $k$, $\mathrm{img}( H_k , T_P ) \cap H_j = \emptyset$ for all $j < k -
c$. Under (i) and (ii), $h$ is consistent along every move between states
that are not dead ends and admissible on every state.

*Proof.* The first claim is the code. By Lemmas 1 and 2, (ii) says that no
move of cost $c$ leads from a state with $h = k$ to a state with $h = j < k
- c$, i.e., $h( x ) \leq c + h( y )$ for every move $( x , y )$ with both
finite. Summing along a cheapest path from $x$ to a goal $y$, whose states
are not dead ends since a goal is reachable from each, gives $h( x ) \leq
d( x , y ) + h( y ) = d( x , y )$ by (i). $\square$

The backward heuristic (towards the initial state, from the terms of a
`SymbolicBackwardHeuristicDomain`) is checked by preimages in the same way,
the initial state taking the place of the goals. The test `sym-astar` checks
that the Manhattan distance and the gap heuristic are consistent in both
directions, and that twice the Manhattan distance is found inconsistent.

## BDDA*

A direction keeps an open list of sets $O_{f,g}$, the states reached with
cost $g$ and $f = g + h$, the closed set $C$ of the states expanded, the sets
$R_g$ of the states reached with cost $g$, and the sequence of the sets
expanded with their $g$. At the start $O_{h,0} = H_h \cap \{ s \}$ for the
initial state $s$ (or the goals, backward, split by the backward $h$). A step
takes the nonempty set $O_{f,g}$ of least $f$ and, among those, of least $g$,
sets $F = O_{f,g} \setminus C$, and, while $F \neq \emptyset$, (i) adds $F$
to $C$, (ii) stops with cost $g$ if $F$ meets the target, (iii) computes
$Y_c = \mathrm{img}( F , T_c ) \setminus C$ for each cost $c$ ($T_c$ the
partitions of cost $c$), (iv) adds $Y_c \cap H_k$ to $O_{g+c+k,g+c}$ for each
$k$, except the part $Y_0 \cap H_{f-g}$, which becomes the next $F$. These are
the buckets of $g$ and $h$ of BDDA* [[128]](../../README.md#ref-128), close
to the SetA* of Jensen, Veloso and Bryant [[129]](../../README.md#ref-129);
with $h = 0$ (option `heuristic = false`, or a domain without terms) the
search is the symbolic version of Dijkstra's algorithm, with buckets of $g$
only, used by symbolic planners for costs [[134]](../../README.md#ref-134),
and with every cost set to 1 (`unit_costs`, `symbolic_bfs()`) it is the
breadth-first search whose layers are the sets $O_{\cdot,g}$. The inner loop
on the moves of cost 0 is the breadth-first search within a bucket that a
cost of 0 requires: by consistency such a move cannot lower $h$, and the
successors with the same $h$ stay in the bucket being expanded.

*Proposition 4.* Let $h$ be consistent, zero on the target and infinite only
on states from which no target is reachable. Then every state is expanded at
most once, a state $x$ from which a target is reachable with $g = d( s , x
)$, and if a target is reachable the search stops with $g = C^*$; if none
is, the open list empties.

*Proof.* A state is expanded at most once, since it enters $C$ when it is
expanded and every set taken from the open list is reduced by $C$. Suppose
now that some state from which a target is reachable is expanded with a $g$
larger than its distance, let $x$ be the first one (in time), expanded at $( f , g )$ with $g > d( s , x )$,
and let $P$ be a cheapest path from $s$ to $x$. When $x$ is expanded, let $y$
be the first state of $P$ not in $C$ (a target is reachable from it, through
$x$, so that $h( y )$ is finite); then $y \neq s$, and the predecessor
of $y$ on $P$ was expanded before $x$, hence with its distance, so that $y$
was inserted in $O_{f_y,g_y}$ with $g_y = d( s , y )$ and, by consistency
along the rest of $P$, $f_y = g_y + h( y ) \leq d( s , x ) + h( x ) < f$.
Since $y \notin C$, that set has not been expanded yet; however, the search
takes the set of least $f$ (and a successor of the set being expanded has
$f' \geq f$ by consistency), which contradicts $f_y < f$. The same argument,
applied to a cheapest path to a target $t$ with $d( s , t ) < g$ when the
search stops at $g$, gives a state of the open list with $f_y \leq d( s , t
) + h( t ) = d( s , t ) < g \leq f$ (since $h( t ) = 0$), again a
contradiction; hence the search stops with $g = C^*$. It ends, since $C$
grows at each step and the states are finite. $\square$

With an inconsistent heuristic the search still returns the cost of a path
(the sets record real costs), but not necessarily an optimal one, since a
state may be closed with a $g$ that is not the least; this is why the
options do not take one, and `consistent()` lets one check a heuristic
before using it. The dead ends (an infinite term) are in no $H_k$ and are
dropped when they are reached.

*The path.* The search records each set $F$ expanded with its $g$, in order.
If $F_i$ (with cost $g_i$) holds the state $t$ that meets the target, a
predecessor of $t$ is found in the last $F_j$, $j < i$, with $g_j = g_i - c$
for a cost $c$ such that $\mathrm{pre}( \{ t \} , T_c ) \cap F_j \neq
\emptyset$; one state of that intersection is taken and the walk repeats from
it until it reaches the root. The walk always succeeds, since every state of
$F_i$ entered its bucket as a successor of a state of a set expanded before,
with $g$ lower by the cost of the move; it costs one preimage of a single
state per cost and per step, i.e., at most $| \{ c \} |$ products per move of
the path, which `path_images` counts. Backward, the walk uses images and
yields the path from the initial state to a goal directly.

## The bidirectional search

Both directions run as above, the forward one by images from the initial
state with $h_F$ and the backward one by preimages from the goals with
$h_B$; every set $X$ inserted in the open list of a direction with cost $g$
(including the roots, and the layers of cost 0 of a step) is intersected with
the sets $R_{g'}$ of the other direction, by increasing $g'$, and the first
nonempty intersection lowers the incumbent $U$ to $g + g'$ if it is smaller.
This is the immediate solution detection of the explicit searches
[[94]](../../README.md#ref-94), done on sets. Successors with $f \geq U$ are
not inserted, since no path through them costs less than $U$. Let $f_D$ and
$g_D$ be the least $f$ and the least $g$ of the nonempty sets in the open
list of the direction $D$, and $\epsilon$ the least cost of an operator; the
search stops as soon as

$$
U \leq \mathrm{LB} = \max \{ f_F , f_B , g_F + g_B + \epsilon \}
$$

(the largest value of $\mathrm{LB}$ seen so far, more precisely), or when an
open list is empty.

*Proposition 5.* Let $h_F$ and $h_B$ be consistent and admissible. When the
search stops, $U = C^*$ if a path exists, and $U = \infty$ otherwise.

*Proof.* $U$ is always the cost of a path, since both $g + g'$ are costs of
real paths from $s$ to the state and from it to a goal. Suppose $U > C^*$ at
some time, and let $P = ( s = p_0 , \ldots , p_k )$ be an optimal path to a
goal. As in the proof of Proposition 4, the first state $p_i$ of $P$ not
closed forward is in the forward open list with $g = d( s , p_i )$ (it is
not pruned, since $f \leq C^* < U$), and the last state $p_j$ not closed
backward is in the backward open list with $g = d( p_j , p_k )$. If $i \geq
j$, the state $p_j$ has been inserted by both directions with its optimal
$g$, and when the second of the two insertions happened it met the set
$R_{g'}$ of the other direction that holds it, so that $U \leq C^*$: a
contradiction. Hence $i < j$, and $f_F \leq d( s , p_i ) + h_F( p_i ) \leq
C^*$, $f_B \leq C^*$ likewise, and $g_F + g_B + \epsilon \leq d( s , p_i ) +
d( p_j , p_k ) + c( p_i , \ldots , p_j ) = C^*$, the middle segment having at
least one move. Thus $\mathrm{LB} \leq C^* < U$ whenever $U > C^*$; since $U$
does not increase, this holds also for the largest $\mathrm{LB}$ seen, and
stopping on $U \leq \mathrm{LB}$ gives $U = C^*$. If the forward open list
empties, every state reachable from $s$ that is not a dead end has been
closed forward, among them the goal $p_k$ of $P$, which was inserted with $g
= C^*$ and met $R_0$ of the backward direction, which holds it; hence $U =
C^*$, and symmetrically for the backward list. With no path, $U$ stays
$\infty$. $\square$

This is the argument of Lemma 1 of the bidirectional page for the explicit
searches, with the sets in place of the nodes; the bound $g_F + g_B +
\epsilon$ alone gives the stopping rule of the bidirectional Dijkstra's
algorithm and of the bidirectional symbolic search of SymBA*
[[130]](../../README.md#ref-130), which uses abstraction heuristics in the
place of our $h_F$ and $h_B$.

*The direction of a step.* Any interleaving of the two directions satisfies
Proposition 5, so the choice is free; a step expands the direction whose
next set is estimated cheaper, the estimate being the nodes of the set times
the ratio between the nodes of the image and those of the set at the last
step of that direction. SymBA* makes the same choice with the time of the
previous step [[130]](../../README.md#ref-130); we use counts of nodes, so
that the search, and every counter it reports, does not depend on the load of
the machine or on the number of workers. The path is the forward walk from
a state of the intersection that set $U$ (its $g$ forward recorded with it),
reversed, followed by the backward walk from the same state.

## The parallel images

The operations on BDDs of Sylvan [[131]](../../README.md#ref-131) are
parallel inside: each recursive call of a relational product, of a
conjunction or of a quantification spawns its two cofactors as tasks of the
work-stealing scheduler Lace, so that a single product of a large set uses
all the workers of Lace. On top of that, `SymbolicImages` chooses how the
images of the partitions of a set are computed: (i) `sequential`, one product
after the other; (ii) `lace`, all of them in one task of Lace that spawns the
products of the two halves of the partitions recursively, and then their
union by a balanced tree of tasks (`images()` and `disjoin()` in `Bdd.h`),
so that the scheduler balances products of very different sizes across the
workers; (iii) `executor`, the partitions dealt round-robin to the workers of
an executor of ADS (the threads of the library, OpenMP or FastFlow), each
submitting its products to Lace from outside, which Lace accepts from several
threads at once.

*Lemma 6.* The costs, the paths and every counter of `SymbolicStats` except
`images_by` are the same for every choice of `SymbolicImages`, of the workers
of the executor and of the workers of Lace.

*Proof.* A reduced ordered BDD is canonical for its function and its order of
the variables [[126]](../../README.md#ref-126), and the order is fixed by the
encoding; every set computed by the search is a function of the sets before
it, independent of how its parts were computed, and hence so are its nodes,
its states, the choices of the direction (which depend only on nodes) and the
state taken from an intersection by the path (the first in the order of the
variables, `sylvan_sat_one`). The images are counted one per partition and
per set, whoever computes them. $\square$

The test `sym-par-*` checks Lemma 6 on the 8-puzzle, the pancakes and the
puzzle of heavy tiles with moves of cost 0, in the three directions, with 1,
2 and 4 workers of Lace and of each executor. The lemma is also why the
counters of these searches are a fair measure of their work: the nodes and
the products are the same on any machine, while the time of a product is
not. We remark that the speedup itself is not measured here.

One point of the implementation needs care: the garbage collection of Sylvan
drops the nodes that no registered handle reaches, and a result that a
thread outside Lace has received but not yet stored in a `Bdd` is not
registered. With one thread driving the search this cannot happen, since the
collection runs only inside a product; with the executor several threads
submit products at once, and `image_into()` therefore stores its result in
the registered handle of its output inside the task of Lace, before the task
ends, so that no collection can find it unregistered.

## The library of BDDs

We interface PBeFS with an existing library instead of writing BDDs of our
own, and we compared the two that fit a C++ library distributed on three
systems.

| | Sylvan 1.12.6 | CUDD 3.0 |
|---|---|---|
| licence | Apache-2.0 | BSD-3-Clause |
| build | CMake (3.23 or later), downloads Lace | autotools only |
| Linux, macOS | in its CI (gcc, clang) | configure and make |
| Windows | in its CI: MSYS2 (MINGW64, UCRT64, CLANG64) and MSVC | no supported build |
| vcpkg, Homebrew, MSYS2 packages | none | none |
| parallel operations | yes, work stealing on Lace | no, one manager per thread |

The table is from the repositories of the two libraries (the workflows of
Sylvan, the files at the root of the mirrors of CUDD) and from the indices
of vcpkg, Homebrew and MSYS2, in November 2026. Sylvan builds on the three
systems, with CMake as PBeFS, and its operations are parallel, which is the
side of symbolic search a library of parallel best-first search is concerned
with; we therefore use Sylvan alone. CUDD would have needed a build of our
own on Windows and would have computed every product on one core. Sylvan is
optional: CMake finds an installed one with `find_package(sylvan)`, or else
downloads the commit pinned in `CMakeLists.txt` (PBEFS_FETCH_SYLVAN), and the
makefiles take it from `SYLVAN_ROOT`; without it the library
`PBeFS::Symbolic` and the cases `sym-*` of the tester are left out, and the
registry of the cases accounts for it.

## Limits

The encoding by cells of the sliding-tile puzzle and by positions of the
pancakes is the natural one, but not necessarily the best: the flip of $k$
pancakes reverses $k$ positions, so that its relation must remember $\lfloor
k / 2 \rfloor$ values while it reads the positions in order, and its diagram
grows roughly as $N^{k/2}$; the stacks of the tests have at most 9 pancakes.
The 15-puzzle is notoriously hard for BDDs of any encoding
[[128]](../../README.md#ref-128), and the tests take its easiest instances of
Korf. The operators have no conditional effects, the SAS+ tasks no axioms,
and the heuristic is additive in terms of few variables: the symbolic
pattern databases and the perimeter abstractions of SymBA*
[[130]](../../README.md#ref-130), and the search for the $k$ best plans of
symk [[132]](../../README.md#ref-132), are not implemented. The relations
are partitioned by size in the order of the operators, without the
conjunctive partitioning by variables of [[130]](../../README.md#ref-130),
and the order of the variables is that of the domain, without reordering.
Finally, a search that fills more than `max_filled` of the largest table of
nodes after a collection stops with `Status::node_limit`, its bound being the
least $f$ in its open list; a single product that fills the table on its own
ends the process, as Sylvan does.
