# Symbolic search

We are concerned here with the searches of `SymbolicSearch.h`, which store sets
of states instead of states, each set represented by the reduced ordered Binary
Decision Diagram (BDD) [[163]](../../README.md#ref-163) of its characteristic
function. A set is expanded as a whole, by a relational product with the
transition relation, as in symbolic model checking
[[190]](../../README.md#ref-190). We first fix the encoding of a
`SymbolicDomain` and prove that the relations and the heuristics built by
`SymbolicModel` are the intended ones. Then we prove the optimality of BDDA*
and of the stopping rule of the bidirectional search, and we end with the
parallel computation of the images, the choice of the library of BDDs and the
limits of the implementation. As for notation, that of the explicit searches is
of [Best-first search and A*](best-first-search.md) and of [Bidirectional
search](bidirectional-search.md).

## The encoding

A `SymbolicDomain` has the variables $v = 0 , \ldots , n - 1$, the variable $v$
with the values $D_v = \{ 0 , \ldots , | D_v | - 1 \}$. A state
$x \in D_0 \times \cdots \times D_{n-1}$ is encoded by writing each $x_v$ in
binary on $b_v = \lceil \log_2 | D_v | \rceil$ bits (the most significant
first) and by concatenating the codes in the order of the variables. The bit
$k$ of the state before a move is the variable $2 k$ of the library, and the
same bit after the move is the variable $2 k + 1$; hence, a set of states is a
BDD $S( x )$ on the even variables and a relation is a BDD $T( x , x' )$ on
both. Those codes that are not values (when $| D_v |$ is not a power of two)
are excluded by the set $V$ of the valid states, those with $x_v < | D_v |$ for
all $v$. The initial state is the minterm of its code, and the goal is $G = V
\cap \bigcup_p \{ x : x_u = d \text{ for } ( u , d ) \in p \}$ over the partial
assignments $p$ given by the domain (one for the puzzles, the goal of a SAS+
task [[133]](../../README.md#ref-133) for `SymbolicSAS`).

An operator $o$ has a cost $c( o ) \geq 0$, a precondition $\mathrm{pre}_o$ (a
set of facts $x_u = d$) and effects, each of which sets a variable $w$ either
to a constant ($x'_w = d$) or to the value that a variable $u$ had before the
move ($x'_w = x_u$). A variable is written by one effect at most, and the
variables that are not written keep their value. With $R_o$ the variables that
$o$ reads or writes and $W_o \subseteq R_o$ those it writes, its relation is

$$ T_o( x , x' ) = \mathrm{pre}_o( x ) \wedge \bigwedge_{e \text{ of } o} E_e(
x , x' ) \wedge \bigwedge_{v \in R_o \setminus W_o} ( x'_v = x_v ) , $$

where $E_e$ is $x'_w = d$ or $x'_w = x_u$, the latter written bit by bit as
$x'_w \leftrightarrow x_u$ (with the high bits of $x'_w$ set to zero if
$b_u < b_w$). The image of $S$ under $T_o$ on the variables $R_o$ is

$$ \mathrm{img}( S , T_o ) = \{ x' : \exists x \in S , \; T_o( x , x' ) , \;
x'_v = x_v \text{ for } v \notin R_o \} , $$

which the relational product of the library computes when it is given the
variables of $R_o$ (in both copies), leaving the others unchanged; the preimage
$\mathrm{pre}( S , T_o )$ is the same with $x$ and $x'$ exchanged.

*Lemma 1.* For any set $S$ of valid states, $\mathrm{img}( S , T_o )$ is the
set of the states reached from $S$ by applying $o$, and
$\mathrm{pre}( S , T_o )$ the set of the states from which applying $o$ reaches
$S$.

*Proof.* A pair $( x , x' )$ is in $T_o$, with $x'_v = x_v$ outside $R_o$, if
and only if $x$ satisfies $\mathrm{pre}_o$, each written variable of $x'$ has
the value given by its effect (computed on $x$, since each $E_e$ reads the
unprimed copy), and each other variable of $x'$ equals that of $x$ (inside
$R_o$ by the frame of $T_o$, outside by the relational product). This is the
definition of applying $o$ to $x$, and both sets follow by projection.
$\square$

An effect may also have a condition $C_e$, a set of facts: it happens only if
$C_e$ holds before the move, and its relation is then
$( C_e \wedge E_e ) \vee ( \neg C_e \wedge x'_w = x_w )$, the variables of
$C_e$ joining $R_o$; Lemma 1 holds with the same proof, since the relation
still describes the move fact by fact. These are the conditional effects of the
planning tasks, e.g., the objects carried by a move of the briefcase of the
test `sym-variants`.

A `SymbolicAxiomDomain` has derived variables as well, two-valued functions of
the others defined by axioms, as in the planning tasks with axioms of the
translator of Fast Downward: they take no bits of the state. The set $D_v$ of
the states where the derived variable $v$ does not have its default value is
computed once, layer by layer, as the least fixed point of
$D_v = \bigvee_{a} C_a$ over the axioms $a$ of $v$; a fact of a condition $C_a$
on a derived variable (of a lower layer, or of the same layer, positively) is
read from its set. A fact on a derived variable is then the set $D_v$ or its
complement, and an operator that reads it reads the variables of the state on
which $D_v$ depends; the projections drop such facts, which only enlarges the
projected moves. The test `sym-sas` searches a task whose goal and one
precondition are derived variables of two layers.

The operators of equal cost are merged, in the order in which the domain gives
them, into partitions $T_P = \bigvee_{o \in P} T_o^{R_P}$ on the union $R_P$ of
their variables ($T_o^{R_P}$ is $T_o$ with the frames $x'_v = x_v$ of
$R_P \setminus R_o$ added). A partition is closed when the next operator would
bring its diagram above `max_relation_nodes` nodes, and with 0 each partition
holds one operator. This is the disjunctive partitioning bounded by size of the
symbolic planners [[193]](../../README.md#ref-193): the image of a disjunction
is the union of the images, while the diagram of the disjunction of all the
operators may be much larger (and usually is) than those of the parts together.
By Lemma 1 applied to each $T_o^{R_P}$ (whose added frames repeat what the
relational product does anyway), $\mathrm{img}( S , T_P ) = \bigcup_{o \in P}
\mathrm{img}( S , T_o )$; hence, the successors of cost $c$ of $S$ are the
union of its images under the partitions of cost $c$.

## The heuristic as a partition of the states

Our searches do not evaluate a heuristic on a state; instead, they split a set
by the value of the heuristic, i.e., they intersect it with the sets
$H_k = \{ x \in V : h( x ) = k \}$. A `SymbolicHeuristicDomain` gives $h$ as a
sum of terms $h( x ) = \sum_{t = 1}^m h_t( x_{U_t} )$, each a table on the
values of a few variables $U_t$ (one cell for the Manhattan distance of the
sliding-tile puzzle, two pancakes of consecutive sizes for the gap heuristic of
the pancakes [[60]](../../README.md#ref-60)), and a term equal to $\infty$
marks a dead end. `SymbolicModel` computes

$$ P_0( 0 ) = \mathbf{1} , \qquad P_t( a ) = \bigvee_{b} \big( P_{t-1}( a - b )
\wedge B_t( b ) \big) , \qquad B_t( b ) = \{ x : h_t( x_{U_t} ) = b \} , $$

and takes $H_k = P_m( k )$, the list of the sets of each value of $h$ that
BDDA* takes as its heuristic (cf. section 3.3 of
[[193]](../../README.md#ref-193)).

*Lemma 2.* $P_m( k )$ is the set of the states whose terms are finite and add
up to $k$; in particular, the sets $H_k$ are disjoint, and their union holds
all the states that are not marked as dead ends.

*Proof.* By induction on $t$, $P_t( a )$ is the set of the states whose first
$t$ terms are finite and add up to $a$. This holds for $t = 0$, and a state is
in $P_t( a )$ if and only if its term $t$ has some finite value $b$ (only the
finite values have a set $B_t( b )$) and its first $t - 1$ terms add up to
$a - b$. $\square$

Note that the diagrams stay rather small when the terms are on disjoint
variables, ordered as the variables are: a level of the diagram of $P_t( a )$
then has at most as many nodes as there are partial sums, i.e.,
$O( h_{\max} )$, and at most $O( h_{\max} | D_v | )$ when consecutive terms
share a variable $v$ (as the pairs of the pancakes do). We remark that a
heuristic given as a table on many variables (e.g., a large pattern database)
would need a different construction, which is not implemented (cf. the last
section).

Since $h$ is never evaluated on a state, one cannot check its consistency on
the states expanded, as an explicit search may do; it can, however, be checked
on all the valid states at once, by images.

*Lemma 3.* Let $Q$ be the states reachable from the initial state. Then
`SymbolicModel::consistent()` returns true if and only if (i) no goal in $Q$
has $h > 0$ and (ii) for each partition $T_P$ of cost $c$ and each $k$,
$\mathrm{img}( H_k \cap Q , T_P ) \cap H_j = \emptyset$ for all $j < k - c$.
Under (i) and (ii), $h$ is consistent along the moves between states of $Q$
that are not dead ends, and it is admissible on $Q$.

*Proof.* The first claim is the code, which finds $Q$ first by a breadth-first
search of the sets. By Lemmas 1 and 2, condition (ii) means that no move of
cost $c$ from a state of $Q$ leads from $h = k$ to $h < k - c$, i.e.,
$h( x ) \leq c + h( y )$ for any move $( x , y )$ between states with finite
$h$, $x \in Q$ (and hence $y \in Q$). Summing along a cheapest path from
$x \in Q$ to a goal $y$ (whose states reach a goal, and hence are not dead
ends) gives $h( x ) \leq d( x , y ) + h( y ) = d( x , y )$ by (i). $\square$

The check is restricted to $Q$ because the encoding has codes that are no state
of the domain: two pancakes in one position, for instance, on which the gap
heuristic is not consistent, although no search ever meets them.

Also, the backward heuristic (towards the initial state, from the terms of a
`SymbolicBackwardHeuristicDomain`) is checked by preimages in the same way,
with the initial state in the place of the goals. In the tests, `sym-astar`
checks that the Manhattan distance and the gap heuristic are consistent in both
directions, and that twice the Manhattan distance is found inconsistent.

## Pattern databases and perimeters

The terms of a domain are not the only source of $H_k$:
`SymbolicModel::strengthen()` raises the heuristic of a direction to the
largest of itself and of heuristics computed by searches on sets, in the
original space and in projections of it, which are the symbolic Pattern
Databases (PDBs) and the perimeters of the symbolic planners
[[196]](../../README.md#ref-196), [[198]](../../README.md#ref-198). A pattern
is a set $P$ of variables, and the projection of a state keeps $x_P$. An
operator is projected by dropping its facts and its effects outside $P$ (a copy
into $P$ from a variable outside it leaves its variable free); an operator that
writes nothing in $P$ is dropped, since its projection would be a loop. By
Lemma 1 applied to the projected relations, a move $x \to y$ of cost $c$ is
projected onto a move $x_P \to y_P$ of cost $c$, or onto no move when
$x_P = y_P$. Given seeds (sets $Z_j$ of states with costs $\gamma_j$),
Dijkstra's algorithm on the projected sets computes the sets of

$$ h_P( x ) = \min_j \; \big( \gamma_j + d_P( x_P , ( Z_j )_P ) \big) , $$

where $d_P$ is the distance in the projected space, found by preimages for the
heuristic towards the goals and by images for that towards the initial state.
With the targets as the single seed, at cost 0, this is the PDB of $P$ (cf.
section 3.3 of [[193]](../../README.md#ref-193)). With `perimeter` $= m > 0$,
Dijkstra's algorithm in the original space first expands $m$ sets from the
targets. Let $C$ be the states it closes, each with its distance $d( x )$, and
$F$ the states it has reached without closing them, each with the least $g$
found; $R$ is the least of these $g$ (or $\infty$ if $F$ is empty). Then the
perimeter gives

$$ h_C( x ) = \begin{cases} d( x ) & x \in C \\ R & x \notin C \end{cases} ,
\qquad h_{C,P}( x ) = \begin{cases} d( x ) & x \in C \\ \min_{u \in F} \big( g(
u ) + d_P( x_P , u_P ) \big) & x \notin C \end{cases} , $$

the latter being the projected search seeded by the sets of $F$ with their $g$,
in the spirit of the perimeter abstractions of Felner and Ofek
[[199]](../../README.md#ref-199) and of SymBA*, whose abstract searches start
from the frontiers of the search in the original space
[[198]](../../README.md#ref-198). The states from which no seed is reachable in
the projection are in no $H_k$, i.e., they are taken as dead ends.

*Lemma 4.* The PDB of a pattern is consistent and zero on the targets, and so
are $h_C$ and $h_{C,P}$; the largest of consistent heuristics zero on the
targets is consistent and zero on the targets, and hence (cf. the proof of
Lemma 3) admissible. The dead ends of all of them are true dead ends.

*Proof.* For the PDB, a move $x \to y$ of cost $c$ is a projected move of cost
$c$ or none. The distance to a set satisfies the triangle inequality; hence
$h_P( x ) \leq c + h_P( y )$, and the targets are seeds at cost 0. For the
perimeter, Dijkstra's algorithm closes each state at a distance below $R$, and
the states of $F$ have $g \geq R$. Take a move $x \to y$ of cost $c$ (towards
the targets). If both are in $C$, $d( x ) \leq c + d( y )$. If $x \notin C$ and
$y \in C$, the expansion of $y$ reached $x$ with $g \leq c + d( y )$, and
$x \in F$; hence $h_C( x ) = R \leq g( x ) \leq c + d( y )$, and
$h_{C,P}( x ) \leq g( x ) + 0$ by taking $u = x$. If $x \in C$ and
$y \notin C$, $d( x ) < R \leq h( y )$ for both, since
$h_{C,P}( y ) \geq \min_F g \geq R$. If neither is in $C$, $h_C$ is $R$ on
both, and $h_{C,P}$ satisfies the triangle inequality of the projected
distance. Also, the targets are in $C$ with $d = 0$. Furthermore, the largest
of heuristics with $h( x ) \leq c + h( y )$ along each move has the same
property. Finally, let a target be reachable from a state $x \notin C$. On a
path from $x$ to a target, the last state outside $C$ is in $F$, since it was
reached by the expansion of the next one; the projection of the path reaches
the projection of that state, and $x$ is therefore not taken as a dead end.
$\square$

Thus, the PDBs are checked like any other heuristic, by `consistent()`. The
test `sym-pdb` does so for the PDBs of two patterns, for a perimeter alone and
for both, on the 8-puzzle, on its heavy tiles with and without moves of cost 0,
on the pancakes and on the SAS+ task of the robot, before comparing the costs
of BDDA* with them against A*. We remark that the PDBs of several patterns are
combined by their largest, which needs no condition on the patterns, while
their sum, which would be larger, is admissible only when no operator writes
into two of them (as for the additive PDBs of the sliding-tile puzzle
[[23]](../../README.md#ref-23)). With `additive`, the PDBs are added, provided
that the patterns are disjoint and that no operator writes into two of them
(`strengthen()` checks both).

*Lemma 5.* Under that condition, the sum of the PDBs of the patterns is
consistent and zero on the targets, and hence admissible.

*Proof.* Let $x \to y$ be a move of cost $c$ by an operator $o$. For each
pattern that $o$ does not write, $x_P = y_P$, and the PDB has the same value on
$x$ and on $y$. For the pattern that $o$ writes, if any, the PDB of $x$ is at
most $c$ plus that of $y$ (Lemma 4). Summing, the sum at $x$ is at most $c$
plus the sum at $y$; on a target each term is 0. $\square$

The test checks the sum on independent robots, where it is exact at the start,
and that patterns written together by one operator are refused.

*Abstractions chosen during the search.* With `abstraction_nodes` $= T$, the
bidirectional search starts blind, as SymBA* does
[[198]](../../README.md#ref-198), and it raises once the heuristic of each
direction by the perimeter of the other, when the next step of both directions
is estimated above $T$ nodes. A blind direction is Dijkstra's algorithm on
sets. Hence its closed sets hold all the states at a distance below the least
$g$ of its open sets (its frontier), and these are the $C$ and $F$ of
$h_{C,P}$, with the targets of the other direction. The pattern $P$ starts from
all the variables and loses, one at a time, the variable whose removal leaves
the smallest projection of the frontier, until the projected search from the
frontier keeps all its sets within $T$ nodes; the open sets of each direction
are then split again by its new heuristic. This greedy relaxation by size is
our choice; SymBA* relaxes the frontier with other abstraction strategies
[[198]](../../README.md#ref-198). By Lemma 4 the new heuristics are consistent
and admissible, with one remark. The bidirectional search does not insert the
states with $f \geq U$, which may leave out of the frontier a state whose path
to a target enters $C$ there; such a state has no path cheaper than $U$, while
the proofs of Propositions 6 and 8 need the consistency only along an optimal
path (whose states have $f \leq C^* < U$ as long as $U > C^*$). Changing the
heuristic between two steps keeps the proof valid too, since it uses only the
current heuristic and the fact that the closed states have their distances. We
remark that the moves of the sliding-tile puzzle write two cells, and those of
the pancakes up to $N$ positions, so that their encodings admit no additive
patterns but trivial ones; the additive PDBs of [[23]](../../README.md#ref-23)
are on the positions of the tiles, a different encoding.

## BDDA*

A direction keeps an open list of sets $O_{f,g}$ (the states reached with cost
$g$ and with $f = g + h$), the closed set $C$ of the states expanded, the sets
$R_g$ of the states reached with cost $g$, and the sequence of the sets
expanded with their $g$. At the start the root (the initial state $s$, or the
goals backward) is split into the sets $O_{k,0}$ by the sets $H_k$ of the
direction. A step takes the nonempty set $O_{f,g}$ of least $f$ and, among
those, of least $g$, and sets $F = O_{f,g} \setminus C$. Then, while
$F \neq \emptyset$, it (i) adds $F$ to $C$, (ii) stops with cost $g$ if $F$
meets the target, (iii) computes $Y_c = \mathrm{img}( F , T_c ) \setminus C$
for each cost $c$ ($T_c$ the partitions of cost $c$), and (iv) adds
$Y_c \cap H_k$ to $O_{g+c+k,g+c}$ for each $k$, except $Y_0 \cap H_{f-g}$,
which becomes the next $F$. These are the buckets of $g$ and $h$ of BDDA*
[[191]](../../README.md#ref-191), close to the SetA* of Jensen, Veloso and
Bryant [[192]](../../README.md#ref-192). With $h = 0$ (the option
`heuristic = false`, or a domain without terms) the search is the symbolic
Dijkstra's algorithm of the planners with costs
[[196]](../../README.md#ref-196), and with all costs set to 1 (`unit_costs`,
`symbolic_bfs()`) it is the breadth-first search, whose layers are the sets
$O_{g,g}$. Note that the loop on the moves of cost 0 is the breadth-first
search within a bucket that such moves require: by consistency they cannot
lower $h$, and the successors with the same $h$ belong to the bucket being
expanded.

*Proposition 6.* Let $h$ be consistent, zero on the target and infinite only on
states from which no target is reachable. Then a state is expanded at most
once, a state $x$ from which a target is reachable with $g = d( s , x )$, and
the search stops with $g = C^*$ if a target is reachable; otherwise, the open
list empties.

*Proof.* A state is expanded at most once, since it enters $C$ when it is
expanded and any set taken from the open list is reduced by $C$. Suppose that a
state from which a target is reachable is expanded with $g$ larger than its
distance, let $x$ be the first one in time, expanded at $( f , g )$ with
$g > d( s , x )$, and let $P$ be a cheapest path from $s$ to $x$. When $x$ is
expanded, let $y$ be the first state of $P$ not in $C$ (a target is reachable
from $y$ through $x$, and hence $h( y )$ is finite). The predecessor of $y$ on
$P$ was expanded before $x$, and therefore with its distance; hence, $y$ was
inserted in $O_{f_y,g_y}$ with $g_y = d( s , y )$ and, by consistency along the
rest of $P$, $f_y = g_y + h( y ) \leq d( s , x ) + h( x ) < f$. Since
$y \notin C$, that set has not been expanded yet, which contradicts the choice
of the set of least $f$ (the successors of the set being expanded have
$f' \geq f$ by consistency). The same argument, applied to a cheapest path to a
target $t$ with $d( s , t ) < g$ when the search stops at $g$, gives a set of
the open list with $f_y \leq d( s , t ) + h( t ) = d( s , t ) < g \leq f$,
again a contradiction; hence the search stops with $g = C^*$. It ends, since
$C$ grows at each step and the states are finite. $\square$

With an inconsistent heuristic the search still returns the cost of a path (the
sets record the costs of real paths), but not necessarily an optimal one, since
a state may be closed with a $g$ larger than its distance; this is why
`consistent()` is there, to check a heuristic before it is used. States marked
as dead ends are in no $H_k$, and they are dropped when they are reached.

*The path.* The search records each set $F$ expanded with its $g$, in order. If
$F_i$ (with cost $g_i$) holds the state $t$ that meets the target, a
predecessor of $t$ is looked for in the sets $F_j$, $j < i$, from the last one,
with $g_j = g_i - c$ for some cost $c$ and
$\mathrm{pre}( \{ t \} , T_c ) \cap F_j \neq \emptyset$; one state of that
intersection is taken, and the walk repeats from it until it reaches the root.
This walk cannot fail, since each state of $F_i$ entered its bucket as a
successor of a state of a set expanded before, with $g$ lower by the cost of
the move. It costs one preimage of a single state per partition and per cost
tried (the counter `path_images`). Backward, the walk takes images, and it
gives the path from the initial state to a goal directly.

*List BDDA\*.* With the option `list`, the open list keeps one set $L_g$ per
$g$, not split by $h$, as the List BDDA* of Torralba et al. (section 3.3 of
[[193]](../../README.md#ref-193)). Before a step, the least $h$ of each $L_g$
is found by intersecting it with $H_0 , H_1 , \ldots$ from a lower bound kept
for it, by increasing $g$, until $g$ reaches the least $f$ found; the set
$L_g \cap H_h$ of the least $( f , g )$ is then moved to the open sets
$O_{f,g}$ and expanded as above. The set taken is that of the least $( f , g )$
over all the states reached, as in BDDA*; hence the two expand the same sets,
and only the number of intersections with the sets $H_k$ (the counter `splits`)
differs, since the states never expanded are never split.

## The k cheapest plans

`SymbolicSearch::top_k( k )` returns the $k$ cheapest plans, i.e., paths from
the initial state to a goal, by nondecreasing cost, a path being allowed to
visit a state more than once, as in the top-$k$ planning of Speck et al.
[[195]](../../README.md#ref-195); two plans differ when their sequences of
states or their costs do. Let $W$ be the states from which a goal is reachable,
computed first by preimages from the goals. The search keeps, for each $g$, the
set $S_g$ of the states of $W$ reached from $s$ by a path of cost $g$, without
removing the states reached with a smaller $g$: $S_0 = \{ s \} \cap W$, and
taking the least $g$ still open, it records $S_g$, enumerates the plans of cost
$g$ and then adds $\mathrm{img}( S_g , T_c ) \cap W$ to $S_{g + c}$ for each
cost $c$. A plan of cost $g$ ending in $t \in S_g \cap G$ is enumerated
backward: the predecessors of a state $u$ with cost $g_u$ are the states of
$\mathrm{pre}( \{ u \} , T_c ) \cap S_{g_u - c}$ for each $c$, and the walk
ends at $s$ with cost 0.

*Proposition 7.* Let each operator cost more than 0. Then $S_g$ is the set of
the states of $W$ reached from $s$ by a path of cost $g$, all the plans of cost
$g$ are enumerated when $S_g$ is taken, each once, and the plans are returned
by nondecreasing cost.

*Proof.* The first claim is by induction on $g$. A state $y$ reached by a path
of cost $g > 0$ whose states are in $W$ (all its states reach a goal through
$y$) has a predecessor $x$ reached with cost $g - c$, $c > 0$, and $x$ is in
$S_{g - c}$, which is taken before $S_g$ and expanded into $S_g$. Conversely,
$S_g$ holds only images of such states. When $S_g$ is taken, each $S_{g'}$ with
$g' < g$ is complete, since the costs are positive and the sets are taken by
increasing $g$. A plan of cost $g$ has each of its states $u$, with the cost
$g_u$ of its prefix, in $S_{g_u}$, so that the backward walk finds it; two
branches of the walk differ in a state or in a cost, and give different plans.
$\square$

With moves of cost 0, $S_g$ would not be complete when it is taken, and a cycle
of cost 0 would give infinitely many plans of the same cost; `top_k()`
therefore refuses them. The search ends, since either $k$ plans are found or
$W$ has no cycle through the reachable states and the sets $S_g$ become empty,
or $g$ exceeds the bound `max_cost`. The test `sym-topk` compares the costs of
the plans with those of a uniform-cost search on the paths themselves, on the
8-puzzle and on its heavy tiles, forward and bidirectionally.

*Bidirectionally.* With the bidirectional direction, `top_k()` also keeps the
sets $S^B_b$ of the states from which a goal is reached by a path of cost $b$,
extending the direction whose next set has fewer nodes; let the forward sets be
complete up to $A$ and the backward ones up to $B$. Before the plans of cost
$c$ are enumerated, the steps go on until $c \leq A$ or $c \leq A + B + 1$. A
plan of cost $c$ is then split after its longest prefix of cost at most $A$:
the prefix goes forward through the sets $S^F$, and the next state (reached by
a move that brings the cost above $A$) has a rest of cost at most
$c - A - 1 \leq B$, which goes through the sets $S^B$. Each plan has exactly
one such split, and is therefore enumerated once, without the checks of
duplicates that a meeting of the two directions at any state would need.

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
least cost of an operator; the search stops when an open list is empty or when

$$ U \leq \mathrm{LB} = \max \{ f_F , f_B , g_F + g_B + \epsilon \} $$

(more precisely, the largest value of $\mathrm{LB}$ seen so far).

*Proposition 8.* Let $h_F$ and $h_B$ be consistent and admissible. When the
search stops, $U = C^*$ if a path exists, and $U = \infty$ otherwise.

*Proof.* $U$ is always the cost of a path, since $g$ and $g'$ are costs of
paths from $s$ to the state and from the state to a goal. Suppose that
$U > C^*$ at some time, and let $P = ( s = p_0 , \ldots , p_k )$ be an optimal
path to a goal. As in the proof of Proposition 6, the first state $p_i$ of $P$
not closed forward is in the forward open list with $g = d( s , p_i )$ (it is
not pruned, since its $f$ is at most $C^* < U$), and the last state $p_j$ not
closed backward is in the backward open list with $g = d( p_j , p_k )$. If
$i \geq j$, the state $p_j$ has been inserted by both directions with its
optimal $g$, and the second of the two insertions met the set $R_{g'}$ of the
other direction that holds it, giving $U \leq C^*$, a contradiction. Hence
$i < j$, and $f_F \leq d( s , p_i ) + h_F( p_i ) \leq C^*$, $f_B \leq C^*$
likewise, and $g_F + g_B + \epsilon \leq C^*$, since the segment from $p_i$ to
$p_j$ has at least one move. Thus $\mathrm{LB} \leq C^* < U$ whenever
$U > C^*$, and also for the largest $\mathrm{LB}$ seen, since $U$ does not
increase; the stop on $U \leq \mathrm{LB}$ therefore gives $U = C^*$. If the
forward open list empties, all the states reachable from $s$ that are not dead
ends have been closed forward, among them the goal $p_k$ of $P$, which was
inserted with $g = C^*$ and met the set $R_0$ of the backward direction; hence
$U = C^*$, and symmetrically for the backward list. With no path, $U$ stays
$\infty$. $\square$

This is the argument of the stopping rule of the explicit bidirectional
searches (sections 7.2 and 7.3 of the manual, cf. [Bidirectional
search](bidirectional-search.md)), with sets in the place of nodes. Without the
two bounds on $f$, the bound $g_F + g_B + \epsilon$ is the stopping rule of the
bidirectional Dijkstra's algorithm, and of the bidirectional symbolic search of
SymBA* [[193]](../../README.md#ref-193), whose heuristics are abstractions
computed during the search instead of the terms given by the domain.

*The direction of a step.* Any interleaving of the two directions satisfies
Proposition 8, and the choice is therefore free. A step expands the direction
whose next set is estimated cheaper to expand; the estimate is the nodes of the
set times the ratio between the nodes of the image and those of the set at the
last step of that direction. Torralba et al. estimate instead the time of the
next step from the time of the previous one and from the sizes of the two sets
(section 3.2 of [[193]](../../README.md#ref-193)); we use counts of nodes
instead, so that the search and its counters do not depend on the load of the
machine or on the workers. The path is the forward walk from a state of the
intersection that set $U$ (with its $g$ forward), reversed, followed by the
backward walk from the same state.

## The parallel images

In Sylvan [[194]](../../README.md#ref-194) the operations are parallel inside:
each recursive call of a relational product, a conjunction or a quantification
spawns its two cofactors as tasks of the work-stealing scheduler Lace, so that
one product of a large set can use all the workers. On top of that,
`SymbolicImages` chooses how the images of the partitions of a set are
computed. These are (i) `sequential`, one product after the other, (ii) `lace`,
all of them in one task of Lace, and (iii) `executor`, the partitions dealt
round-robin to the workers of an executor of ADS. With `lace`, the task spawns
the products of the two halves of the partitions recursively and then joins
them by a balanced tree of disjunctions (`images()` and `disjoin()` in
`Bdd.h`), and hence the scheduler balances products of rather different sizes.
With `executor`, the executor is the threads of the library, OpenMP or
FastFlow, and each of its workers submits its products to Lace, which accepts
them from several threads at once.

*Lemma 9.* The costs, the paths and the counters of `SymbolicStats` other than
`images_by` do not depend on the choice of `SymbolicImages`, on the workers of
the executor or on the workers of Lace.

*Proof.* A reduced ordered BDD is canonical for its function and its order of
the variables [[163]](../../README.md#ref-163), and the order is fixed by the
encoding. Each set computed by the search is a function of the sets before it,
whatever way its parts were computed; hence so are its nodes, its states, the
choice of the direction (which depends on nodes only) and the state taken from
an intersection by the path (the first in the order of the variables, by
`sylvan_sat_one`). Finally, the products are counted one per partition and per
set, whoever computes them. $\square$

Indeed, the test `sym-par-*` checks Lemma 9 on the 8-puzzle, the pancakes and
the 8-puzzle of heavy tiles with moves of cost 0, in the three directions, with
1, 2 and 4 workers of Lace and of each executor. Lemma 9 also makes the
counters a fair measure of the work of these searches, since the nodes and the
products are the same on any machine while the time of a product is not. We
have not measured the speedup itself.

One point of the implementation needs care. The garbage collection of Sylvan
drops the nodes that no registered handle reaches, and a result received by a
thread outside Lace is not registered until it is stored in a `Bdd`. With one
thread driving the search this cannot matter, since a collection runs only
inside a product; with the executor, however, several threads submit products
at once, and `image_into()` therefore stores its result in the registered
handle of its output inside the task of Lace, before the task ends.

## The library of BDDs

PBeFS interfaces an existing library instead of implementing the diagrams, and
we compared the two that fit a C++ library distributed on three systems. In the
following table each column is one library at the version compared, and each
row a property of it: its licence, its build system, how it builds on Linux and
macOS and on Windows, its packages on vcpkg, Homebrew and MSYS2, and how its
operations run in parallel.

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
therefore use Sylvan, while CUDD would have needed a port to Windows and would
have computed each product on one core. Sylvan is optional, and how CMake and
the makefiles find it or download it (at the commit pinned in `CMakeLists.txt`)
is described in the section on building of the [README](../../README.md).
Without it, the library `PBeFS::Symbolic` and the cases `sym-*` of the tester
are left out, and the case `registry` accounts for them.

## Limits

Of course, the encodings are not necessarily (and perhaps rarely) the best, and
their choice matters a great deal. The stack of pancakes as the state, with a
variable per position, makes the flip of $k$ pancakes copy the position $i$
into the position $k - 1 - i$: its relation reads the positions in order and
must keep up to $k$ values at once, and with 8 pancakes our first searches in
that encoding filled a table of $2^{24}$ nodes. `SymbolicPancake` takes instead
the position of each pancake as its variable, the inverse of the stack, and a
flip is then a map of each variable on itself, $x'_p = k - 1 - x_p$ for
$x_p < k$ (an effect `map_of()` of `SymbolicDomain.h`), whose relation has a
few nodes per variable; the gap heuristic is again a sum of terms, on pairs of
pancakes of consecutive sizes. The tests take at most 9 pancakes, and of the
15-puzzle four of the six easiest instances of Korf (numbers 55, 42, 79 and
12), whose largest sets have under 100000 nodes, while the numbers 16 and 71,
whose sets reach 400000 nodes, are measured only by AStarBench
(`--algo sym-bdda`). The patterns of the PDBs are given by the user, and the
abstractions chosen during the search are chosen once per direction, greedily
by size, while SymBA* refines its abstractions when a step is estimated too
hard [[198]](../../README.md#ref-198); the search for the $k$ best plans of
`top_k()` is blind, as that of symk [[195]](../../README.md#ref-195), and needs
positive costs. Our partitions follow the order of the operators. The order of
the variables is that of the domain, or one given, or one found by
`SymbolicLayout::optimize`, which moves one variable at a time to the place
that most lowers the sum, over the operators and the terms of the heuristic, of
the squared span of their variables in bits. It is a static order, chosen once,
since Sylvan has no dynamic reordering, and it is good for that sum only.
Finally, a search that fills more than `max_filled` of the largest table of
nodes after a collection stops with `Status::node_limit`, with the least $f$ of
its open list as its bound; a single product that fills the table ends the
process, as Sylvan does.