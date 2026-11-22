# 18. Symbolic search {#ch-18}

[Source: [`include/PBeFS/symbolic.h`][symbolic.h] (the header of the whole
chapter), [`include/PBeFS/symbolic/Bdd.h`][Bdd.h] and [`src/Bdd.cpp`][Bdd.cpp]
(the diagrams and the scope of
Sylvan), [`include/PBeFS/symbolic/SymbolicDomain.h`][SymbolicDomain.h] (the
concepts), [`include/PBeFS/symbolic/SymbolicModel.h`][SymbolicModel.h] (the
compilation of a
domain), [`include/PBeFS/symbolic/SymbolicSearch.h`][SymbolicSearch.h] (the
searches), [`include/PBeFS/symbolic/SymbolicPuzzles.h`][SymbolicPuzzles.h]
and [`include/PBeFS/symbolic/SymbolicSAS.h`][SymbolicSAS.h] (the
domains); [`test/symbolic.cpp`][symbolic.cpp]]

> **Status — beta.** The code that this chapter describes is in `main`
> but not yet in a release, and its interface may still change.

We are concerned in this chapter with the searches that store sets of states
instead of states, and that expand a set as a whole by a relational product
with the transition relation, as symbolic model checking
does [[190]](Z-bibliography.md#ref-190); each set is represented by the reduced
ordered Binary Decision Diagram (BDD) of
Bryant [[163]](Z-bibliography.md#ref-163) of its characteristic function. Since
the size of a diagram depends on the structure of the set rather than on its
number of states, a set of billions of states may have a diagram of a few
thousand nodes; hence, these searches may solve, within a given memory,
problems whose explicit searches would not fit. The converse also happens (and
quite often), since the diagrams of some sets grow exponentially whatever the
order of the variables. We give the breadth-first search, Dijkstra's algorithm
and BDDA$^*$ of Edelkamp and Reffel [[191]](Z-bibliography.md#ref-191), in the
forward, backward and bidirectional directions, with heuristics that are sums
of "small" tables and the abstractions of the symbolic planners of Torralba et
al. [[193]](Z-bibliography.md#ref-193), [[198]](Z-bibliography.md#ref-198). We
also give the $k$ cheapest plans in the manner of Speck et
al. [[195]](Z-bibliography.md#ref-195) and the computation of the images of a
set in parallel, all on the diagrams of
Sylvan [[194]](Z-bibliography.md#ref-194), an optional dependency of the
library (cf. [§18.11](#sec-18-11)).

The structure of the chapter is as follows. In [§18.1](#sec-18-1) we give the
diagrams and the scope in which Sylvan runs, and in [§18.2](#sec-18-2) the
concepts by which a problem is given to the searches. We show
in [§18.3](#sec-18-3) how a domain is compiled into relations and why the
images computed are the "intended" ones; [§18.4](#sec-18-4) does the same for
the heuristics, together with the check of their consistency. The pattern
databases, the perimeters and the abstractions chosen during the search are
treated in [§18.5](#sec-18-5). Then come BDDA$^*$ with its proof
([§18.6](#sec-18-6)), the bidirectional search ([§18.7](#sec-18-7)), the $k$
cheapest plans ([§18.8](#sec-18-8)), the parallel images ([§18.9](#sec-18-9))
and the domains ([§18.10](#sec-18-10)). Finally, [§18.11](#sec-18-11)
to [§18.13](#sec-18-13) give the choice of the library with the limits of the
implementation, the notes and the exercises. As for notation, that of the
explicit searches is the one of [Chapter 1](01-introduction.md#ch-1)
and [Chapter 7](07-bidirectional.md#ch-7).

## 18.1 The diagrams and their runtime {#sec-18-1}

A `Bdd` ([`Bdd.h:63`][Bdd.h:63]) is a Boolean function of the variables $0, 1,
\ldots$ of Sylvan, ordered by their index and held through a handle that the
garbage collector of Sylvan knows: each `Bdd` registers the address of its
handle when it is built, and releases it when it is destroyed. Since the
diagrams are reduced and ordered, equal functions have the same diagram, and
`==` compares two functions in constant time. Besides the constants
`Bdd::one()` and `Bdd::zero()`, one has `Bdd::variable( v )`, the function of
the variable $v$, and `Bdd::cube( vars , values )`, the conjunction of the
literals of `vars`, and `Bdd::variables( vars )`, the set of variables on which
the quantifications and the counts act. The operators `&`, `|`, `-` and `!` are
the conjunction, the disjunction, the difference and the negation; `equiv( b )`
is the equivalence and `exists( vars )` the existential quantification. A
function is queried by `empty()`, `full()`, `disjoint( b )` and `nodes()`, by
`count( vars )` (the number of assignments of `vars` that satisfy it, as a
`double`) and by `pick( vars )` (one of them, the first in the order of the
variables).

In Sylvan, the variables of a state are the even ones, $x_k = 2k$, and those of
a successor the odd ones, $x'_k = 2k + 1$: a set of states $S(x)$ is then a
diagram on the even variables, and a relation $T(x, x')$ a diagram on both. The
free functions of the header give the relational products: `image( S , T , vars
)` is $\exists x : S(x) \wedge T(x, x')$ renamed to $x$, and `preimage( S , T ,
vars )` is $\exists x' : T(x, x') \wedge S(x')$, where `vars` holds the
variables that $T$ reads or writes (the variables of a state outside `vars` are
left unchanged by $T$). Then `image_into( S , T , vars , backward , out )`
writes the image (or the preimage) into `out`, and it is the only operation
that several threads may call at once. We compute the images of $S$ under all
the relations `T[ i ]` at once, as tasks of Lace, by `images( S , T , vars ,
backward , out )`. Finally, `disjoin( v )` is the disjunction of a vector as a
balanced tree of tasks ([§18.9](#sec-18-9)), and `nodes( v )` counts the nodes
of several diagrams, those shared counted once.

Sylvan and Lace (its scheduler) are one per process. A `BddRuntime`
([`Bdd.h:255`][Bdd.h:255]) starts Lace with its workers (0 for one per core)
and Sylvan with its table of nodes and its cache within a given number of
bytes, and it stops both when it is destroyed. Each `Bdd` must be destroyed
before the `BddRuntime`, and a second `BddRuntime` alive at the same time
throws `std::logic_error`. Its static members give the workers of Lace
(`workers()`), whether a runtime is alive (`running()`) and a garbage
collection (`collect()`); they also give the fraction of the largest table of
nodes in use (`filled()`), by which a search decides to stop before the table
is full.

## 18.2 The problem and its interface {#sec-18-2}

We give a problem to the searches of the chapter as a type `D` that satisfies
the concept `SymbolicDomain` ([`SymbolicDomain.h:202`][SymbolicDomain.h:202]).
Its states are vectors of values of the variables $v = 0, \ldots, n - 1$, where
$v$ ranges over $D_v = \{0, \ldots, |D_v| - 1\}$. The domain provides the
integral type `cost_type` and the methods `num_variables()`, `range( v )` (that
is, $|D_v|$) and `initial_values()`. It also provides `for_each_goal( f )`,
which calls `f` with one or more partial assignments (a state is a goal if it
satisfies at least one of them), and `for_each_operator( f )`, which calls `f`
with each operator. A fact $x_v = d$ is a `SymbolicFact`
([`SymbolicDomain.h:65`][SymbolicDomain.h:65]), and an operator is a
`SymbolicOperator` ([`SymbolicDomain.h:112`][SymbolicDomain.h:112]), with a
cost $\geq 0$, a precondition `pre` (a vector of facts) and its effects. Each
effect is a `SymbolicEffect` ([`SymbolicDomain.h:75`][SymbolicDomain.h:75]),
built by `assign( v , d )` ($x'_v = d$), by `copy_of( v , u )` ($x'_v = x_u$,
the value that $u$ had before the move), or by `map_of( v , u , table )`
($x'_v$ is then the entry of `table` at $x_u$). An effect whose vector `cond`
of facts is not empty is a conditional effect, which happens only where its
condition holds before the move. An operator applies all its effects at once,
each reading the values before the move, and a variable is written by one
effect at most.

Four concepts refine `SymbolicDomain`. A `SymbolicHeuristicDomain`
([`SymbolicDomain.h:242`][SymbolicDomain.h:242]) gives the terms of an additive
heuristic towards the goals by `for_each_heuristic_term( f )`, while a
`SymbolicBackwardHeuristicDomain`
([`SymbolicDomain.h:255`][SymbolicDomain.h:255]) gives the terms of an estimate
of the cost from the initial state, by `for_each_backward_heuristic_term( f )`.
A term is a `SymbolicTerm` ([`SymbolicDomain.h:138`][SymbolicDomain.h:138]),
i.e., a table of values on the values of a few variables (the first variable
varying fastest); a value equal to `infinite_term< Cost >()` marks the states
from which no goal is reachable. A `SymbolicAxiomDomain`
([`SymbolicDomain.h:227`][SymbolicDomain.h:227]) has derived variables, i.e.,
two-valued functions of the other variables that no operator writes and that
are not part of the state, given by `is_derived( v )`, `default_value( v )`,
`layer( v )` and `for_each_axiom( f )`. Its axioms are `SymbolicAxiom`s
([`SymbolicDomain.h:122`][SymbolicDomain.h:122]), each of which sets its
variable to its value where its condition holds, as the axioms of the planning
tasks of the translator of Fast Downward [[134]](Z-bibliography.md#ref-134).
Finally, a `DecodableSymbolicDomain`
([`SymbolicDomain.h:268`][SymbolicDomain.h:268]) has the type `state_type` and
the method `decode( x )`, which gives the state of an explicit `Domain`
([§1.5](01-introduction.md#sec-1-5)) whose variables have the values `x`: one
can thus check (or use) on the explicit domain a path found symbolically, and
the tests rest on it.

All the searches take their options in a `SymbolicOptions`
([`SymbolicSearch.h:93`][SymbolicSearch.h:93]), i.e., the `direction`
(`forward`, `backward` or `bidirectional`), whether the `heuristic` and its
`terms` are used, the `unit_costs` of a breadth-first search, and the
`patterns`, the `perimeter` and the `additive` combination
of [§18.5](#sec-18-5). Then come the `layout` of the variables, the
`abstraction_nodes` of the abstractions chosen during the search, the open list
as `list` ([§18.6](#sec-18-6)) and the `max_relation_nodes` of the partitions.
Finally, one sets the way the `images` are computed and the `threads` of the
executor ([§18.9](#sec-18-9)), the fraction `max_filled` of the table of nodes
at which a search stops, and whether the `path` is rebuilt. A search returns a
`SymbolicSolution` ([`SymbolicSearch.h:188`][SymbolicSearch.h:188]), whose
`status` is `Status::solved`, `Status::unsolvable` or `Status::node_limit`, and
which holds the optimal `cost`, a lower `bound` on it, the `path` (the values
of the variables of each state, from the initial state to a goal) and its
counters, a `SymbolicStats` ([`SymbolicSearch.h:148`][SymbolicSearch.h:148]).
These count the sets expanded in each direction, the relational products of the
expansions and of the path, the states and the nodes expanded, and the largest
set and image. They also count the partitions and their nodes, the nodes of the
sets of the heuristic, the splits of a set by them, and the products and steps
of the abstractions. None of them depends on the number of workers ([Lemma
18.9](#lem-18-9)), and we therefore take them as the measure of these searches.

The class `SymbolicSearch< D , Executor >`
([`SymbolicSearch.h:248`][SymbolicSearch.h:248]) compiles the domain when it is
built, and it runs a search with `solve()`, the $k$ cheapest plans with `top_k(
k )` ([§18.8](#sec-18-8)) and the states reachable with `reachable( backward
)`. The latter returns a `SymbolicReachable`
([`SymbolicSearch.h:206`][SymbolicSearch.h:206]), i.e., the number of states
first reached after each number of moves, with the products made and the nodes
of the largest layer. Also, `set_executor( ex )` gives the search the executor
of [§18.9](#sec-18-9), and `model()` returns the compiled model. The free
functions `symbolic_search( d , o )`, `symbolic_bfs( d , o )` and
`symbolic_reachable( d , backward , o )` build the search and run it once. Note
that a `BddRuntime` must be alive while a search exists; the
header [`PBeFS/symbolic.h`][symbolic.h] includes all those of the chapter.

## 18.3 The encoding and the relations {#sec-18-3}

A `SymbolicModel` ([`SymbolicModel.h:114`][SymbolicModel.h:114]) compiles a
`SymbolicDomain` into diagrams. We encode a state $x \in D_0 \times \cdots
\times D_{n-1}$ by writing each $x_v$ in binary on $b_v = \lceil \log_2 |D_v|
\rceil$ bits, the most significant first, and by concatenating the codes in the
order of the variables. This order is that of the domain, unless a
`SymbolicLayout` ([`SymbolicModel.h:74`][SymbolicModel.h:74]) gives another
`order`, or asks with `optimize` for one found by a local search
([§18.11](#sec-18-11)). The bit $k$ of the state before a move is the variable
$2k$ of Sylvan, and the same bit after the move is the variable $2k + 1$. Thus,
the two copies of a bit are next to each other, which keeps "small" the
relations of the moves that leave a variable as it is. When $|D_v|$ is not a
power of two, the codes that are not values are excluded by the set $V$ of the
valid states, i.e., those with $x_v < |D_v|$ for all $v$. Then, the initial
state is the minterm of its code, and the goal is $G = V \cap \bigcup_p \{x :
x_u = d \text{ for } (u, d) \in p\}$ over the partial assignments $p$ of the
domain. In the model one can find the diagrams of these sets (`initial()`,
`goal()`, `valid()`), the set of a value of a variable (`value( v , d , next
)`) and the set of one state (`state( x )`), together with the inverse
operations (`decode( S )`, the values of the first state of a set, and `count(
S )`, its number of states).

With $R_o$ the variables that the operator $o$ reads or writes and $W_o
\subseteq R_o$ those it writes, the relation of $o$ is

$$
T_o(x, x') = \mathrm{pre}_o(x) \wedge \bigwedge_{e \text{ of } o} E_e(x, x')
\wedge \bigwedge_{v \in R_o \setminus W_o} (x'_v = x_v) ,
$$

where $E_e$ is $x'_w = d$ for an assignment and $x'_w = x_u$ for a copy. The
latter is written bit by bit as $x'_w \leftrightarrow x_u$, with the high bits
of $x'_w$ set to zero if $b_u < b_w$ (a copy into a variable of fewer bits is
refused with `std::invalid_argument`); a map is the disjunction over $d$ of
$x_u = d \wedge x'_w = \mathrm{table}(d)$. The image of $S$ under $T_o$ on the
variables $R_o$ is

$$
\mathrm{img}(S, T_o) = \{x' : \exists x \in S, \; T_o(x, x'), \; x'_v = x_v
\text{ for } v \notin R_o\} ,
$$

which `image` computes when it is given the variables of $R_o$ in both copies;
the preimage $\mathrm{pre}(S, T_o)$ is the same with $x$ and $x'$ exchanged.

<a id="lem-18-1"></a>**Lemma 18.1** (The images).
For any set $S$ of valid states, $\mathrm{img}(S, T_o)$ is the set of the
states reached from $S$ by applying $o$, and $\mathrm{pre}(S, T_o)$ is the set
of the states from which applying $o$ reaches $S$.

*Proof.* A pair $(x, x')$ is in $T_o$, with $x'_v = x_v$ outside $R_o$, if and
only if (i) $x$ satisfies $\mathrm{pre}_o$, (ii) each written variable of $x'$
has the value given by its effect (computed on $x$, since each $E_e$ reads the
unprimed copy), and (iii) each other variable of $x'$ equals that of $x$
(inside $R_o$ by the frame of $T_o$, outside by the relational product). This
is the definition of applying $o$ to $x$, and both sets follow by projection.
$\square$

A conditional effect $e$ with the condition $C_e$ has the relation $(C_e \wedge
E_e) \vee (\neg C_e \wedge x'_w = x_w)$, the variables of $C_e$ joining $R_o$,
and [Lemma 18.1](#lem-18-1) holds with the same proof, since the relation still
describes the move fact by fact. As for the derived variables of a
`SymbolicAxiomDomain`, they take no bits of the state. Instead, we compute
once, layer by layer, the set $D_v$ of the states where the derived variable
$v$ does not have its default value; it is the least fixed point of $D_v =
\bigvee_a C_a$ over the axioms $a$ of $v$; a fact of a condition on a derived
variable (of a lower layer, or of the same layer with the value that is not the
default) is read from its set. A fact on a derived variable is then the set
$D_v$ or its complement, and an operator that reads it reads the variables of
the state on which $D_v$ depends.

We merge the operators of equal cost, in the order in which the domain gives
them, into partitions $T_P = \bigvee_{o \in P} T_o^{R_P}$ on the union $R_P$ of
their variables, which `partitions()` returns (here $T_o^{R_P}$ is $T_o$ with
the frames $x'_v = x_v$ of $R_P \setminus R_o$ added). A partition is closed
when the next operator would bring its diagram above `max_relation_nodes`
nodes, and with 0 each partition holds one operator. This is the disjunctive
partitioning bounded by size of the symbolic
planners [[193]](Z-bibliography.md#ref-193): the image of a disjunction is the
union of the images, while the diagram of the disjunction of all the operators
may be much larger (and usually is) than those of its parts together. By [Lemma
18.1](#lem-18-1) applied to each $T_o^{R_P}$, whose added frames repeat what
the relational product does anyway, one has $\mathrm{img}(S, T_P) = \bigcup_{o
\in P} \mathrm{img}(S, T_o)$; hence, the successors of cost $c$ of $S$ are the
union of its images under the partitions of cost $c$.

## 18.4 The heuristic as a partition of the states {#sec-18-4}

A symbolic search does not evaluate a heuristic on a state: it splits a set by
the value of the heuristic, i.e., it intersects the set with the sets $H_k =
\{x \in V : h(x) = k\}$. The heuristic of a `SymbolicHeuristicDomain` is a sum
of terms $h(x) = \sum_{t=1}^m h_t(x_{U_t})$, each a table on the values of a
few variables $U_t$ (e.g., one cell for the Manhattan distance of the
sliding-tile puzzle, and two pancakes of consecutive sizes for the gap
heuristic [[60]](Z-bibliography.md#ref-60)). For such a sum, `SymbolicModel`
computes

$$
P_0(0) = \mathbf{1} , \qquad P_t(a) = \bigvee_b \big( P_{t-1}(a - b) \wedge
B_t(b) \big) , \qquad B_t(b) = \{x : h_t(x_{U_t}) = b\} ,
$$

and takes $H_k = P_m(k)$, which gives the list of sets (`levels( backward )`)
that BDDA$^*$ uses as its heuristic (cf. section 3.3
of [[193]](Z-bibliography.md#ref-193)).

<a id="lem-18-2"></a>**Lemma 18.2** (The sets of the heuristic).
$P_m(k)$ is the set of the states whose terms are finite and add up to $k$; in
particular, the sets $H_k$ are disjoint, and their union holds all the states
that are not marked as dead ends.

*Proof.* By induction on $t$, $P_t(a)$ is the set of the states whose first $t$
terms are finite and add up to $a$. Indeed, this holds for $t = 0$.
Furthermore, a state is in $P_t(a)$ if and only if its term $t$ has some finite
value $b$ (only the finite values have a set $B_t(b)$) and its first $t - 1$
terms add up to $a - b$. $\square$

The diagrams usually stay rather "small" when the terms are on disjoint
variables, ordered as the variables are: a level of the diagram of $P_t(a)$
then has at most as many nodes as there are partial sums, i.e., $O(h_{\max})$,
and at most $O(h_{\max} |D_v|)$ when consecutive terms share a variable $v$ (as
the pairs of the pancakes do). Of course, a heuristic given as a table on many
variables (e.g., a "large" pattern database of an explicit search) would
usually not fit this construction; the pattern databases of the symbolic
searches are computed as sets instead ([§18.5](#sec-18-5)).

Since $h$ is never evaluated on a state, its consistency cannot be checked on
the states expanded, as an explicit search may do. One can, however, check it
on all the reachable states at once, by images.

<a id="lem-18-3"></a>**Lemma 18.3** (The check of consistency).
Let $Q$ be the states reachable from the initial state. Then `consistent( false
)` returns true if and only if (i) no goal in $Q$ has $h > 0$ and (ii) for each
partition $T_P$ of cost $c$ and each $k$, $\mathrm{img}(H_k \cap Q, T_P) \cap
H_j = \emptyset$ for all $j < k - c$. Under (i) and (ii), $h$ is consistent
along the moves between states of $Q$ that are not dead ends, and it is
admissible on $Q$.

*Proof.* The first claim is the code, which finds $Q$ first by a breadth-first
search of the sets. By [Lemma 18.1](#lem-18-1) and [Lemma 18.2](#lem-18-2),
condition (ii) means that no move of cost $c$ from a state of $Q$ leads from $h
= k$ to $h < k - c$, i.e., $h(x) \leq c + h(y)$ for any move $(x, y)$ between
states with finite $h$, $x \in Q$ (and hence $y \in Q$). Take a cheapest path
from $x \in Q$ to a goal $y$, whose states reach a goal and are therefore not
dead ends; summing along it gives $h(x) \leq d(x, y) + h(y) = d(x, y)$ by (i).
$\square$

We restrict the check to $Q$ because the encoding has codes that are no state
of the domain (for instance, two pancakes in one position, on which the gap
heuristic is not consistent, although no search ever meets them). The backward
heuristic is checked in the same way by `consistent( true )`, by preimages,
with the initial state in the place of the goals.

## 18.5 Pattern databases, perimeters and abstractions {#sec-18-5}

We now raise the heuristic of a direction by heuristics computed with searches
on sets, in the original space and in projections of it, which are the symbolic
Pattern Databases (PDBs) and the perimeters of the symbolic
planners [[196]](Z-bibliography.md#ref-196), [[198]](Z-bibliography.md#ref-198)
(the explicit abstractions of the library are
in [§16.6](16-planning.md#sec-16-6)). In `SymbolicModel`, `strengthen( backward
, patterns , perimeter , additive )` replaces the sets $H_k$ of the direction
by those of the largest of its heuristic and of these. A pattern is a set $P$
of variables, and the projection of a state keeps $x_P$. We project an operator
by dropping its facts and its effects outside $P$ (a copy into $P$ from a
variable outside it leaves its variable free), and an operator that writes
nothing in $P$ is dropped, since its projection would be a loop. By [Lemma
18.1](#lem-18-1) applied to the projected relations, a move $x \to y$ of cost
$c$ is projected onto a move $x_P \to y_P$ of cost $c$, or onto no move when
$x_P = y_P$. Given seeds, i.e., sets $Z_j$ of states with costs $\gamma_j$,
Dijkstra's algorithm on the projected sets computes the sets of

$$
h_P(x) = \min_j \big( \gamma_j + d_P(x_P, (Z_j)_P) \big) ,
$$

where $d_P$ is the distance in the projected space, found by preimages for the
heuristic towards the goals and by images for that towards the initial state.
With the targets as the single seed, at cost 0, this is the PDB of $P$. With
`perimeter` $= m > 0$, Dijkstra's algorithm in the original space first expands
$m$ sets from the targets. Let $C$ be the states it closes, each with its
distance $d(x)$, and $F$ those it has reached without closing them, each with
the least $g$ found, and let $R$ be the least of these $g$ (or $\infty$ if $F$
is empty). The perimeter then gives

$$
h_C(x) = \begin{cases} d(x) & x \in C \\ R & x \notin C \end{cases} ,
\qquad
h_{C,P}(x) = \begin{cases} d(x) & x \in C \\
\min_{u \in F} \big( g(u) + d_P(x_P, u_P) \big) & x \notin C
\end{cases} ,
$$

where $h_{C,P}$ is the projected search seeded by the sets of $F$ with their
$g$, in the spirit of the perimeter abstractions of Felner and
Ofek [[199]](Z-bibliography.md#ref-199). The states from which no seed is
reachable in the projection are in no $H_k$, i.e., they are taken as dead ends.

<a id="lem-18-4"></a>**Lemma 18.4** (The abstractions are consistent).
The PDB of a pattern is consistent and zero on the targets, and so are $h_C$
and $h_{C,P}$; the largest of consistent heuristics zero on the targets is
consistent and zero on the targets, and hence admissible. The dead ends of all
of them are true dead ends.

*Proof.* For the PDB, a move $x \to y$ of cost $c$ is a projected move of cost
$c$ or none; the distance to a set satisfies the triangle inequality, and hence
$h_P(x) \leq c + h_P(y)$, the targets being seeds at cost 0. For the perimeter,
Dijkstra's algorithm closes each state at a distance below $R$, and the states
of $F$ have $g \geq R$. Take a move $x \to y$ of cost $c$ (towards the
targets). If both are in $C$, then $d(x) \leq c + d(y)$. If $x \notin C$ and $y
\in C$, the expansion of $y$ reached $x$ with $g \leq c + d(y)$, and $x \in F$;
hence $h_C(x) = R \leq g(x) \leq c + d(y)$, and $h_{C,P}(x) \leq g(x) + 0$ by
taking $u = x$. If $x \in C$ and $y \notin C$, then $d(x) < R \leq h(y)$ for
both, since $h_{C,P}(y) \geq \min_F g \geq R$. If neither is in $C$, $h_C$ is
$R$ on both, and $h_{C,P}$ satisfies the triangle inequality of the projected
distance. The targets are in $C$ with $d = 0$, and the largest of heuristics
with $h(x) \leq c + h(y)$ along each move has the same property; admissibility
follows by summing along a path, as in the proof of [Lemma 18.3](#lem-18-3).
Finally, let a target be reachable from a state $x \notin C$. On a path from
$x$ to a target, the last state outside $C$ is in $F$, since it was reached by
the expansion of the next one. The projection of the path reaches the
projection of that state, and hence $x$ is not taken as a dead end. $\square$

We generally combine the PDBs of several patterns by their largest, which needs
no condition on the patterns. Their sum would be larger, but it is admissible
only when no operator writes into two of them (as for the additive PDBs of the
sliding-tile puzzle [[23]](Z-bibliography.md#ref-23)). With `additive` the PDBs
are added, after `strengthen()` has checked that the patterns are disjoint and
that no operator writes into two of them (it throws `std::invalid_argument`
otherwise); the sum then raises the heuristic, and the perimeter, if any,
raises it separately, without seeding the patterns.

<a id="lem-18-5"></a>**Lemma 18.5** (Additive pattern databases).
Under that condition, the sum of the PDBs of the patterns is consistent and
zero on the targets, and hence admissible.

*Proof.* Let $x \to y$ be a move of cost $c$ by an operator $o$. For each
pattern that $o$ does not write, $x_P = y_P$, and the PDB has the same value on
$x$ and on $y$, while for the pattern that $o$ writes, if any, the PDB of $x$
is at most $c$ plus that of $y$ ([Lemma 18.4](#lem-18-4)). Summing, the sum at
$x$ is at most $c$ plus the sum at $y$, and on a target each term is 0.
$\square$

Note that the moves of the sliding-tile puzzle in the encoding
of [§18.10](#sec-18-10) write two cells, and those of the pancakes up to $N$
positions; hence, their encodings admit no additive patterns but "trivial"
ones. The additive PDBs of [[23]](Z-bibliography.md#ref-23) are on the
positions of the tiles, which is a different encoding.

**Abstractions chosen during the search.** With `abstraction_nodes` $= T > 0$,
the bidirectional search starts "blind", as SymBA$^*$
does [[198]](Z-bibliography.md#ref-198). (The constructor of `SymbolicSearch`
throws `std::invalid_argument` unless the search is bidirectional, with the
option `heuristic` on and no other heuristic.) When the next step of both
directions is estimated above $T$ nodes, the search raises, once, the heuristic
of each direction by the perimeter of the other, through `strengthen_from(
backward , closed , frontier , budget )` of `SymbolicModel`. A blind direction
is Dijkstra's algorithm on sets, and hence its closed sets hold all the states
at a distance below the least $g$ of its open sets: these are the $C$ and $F$
of $h_{C,P}$, with the targets of the other direction. The pattern $P$ starts
from all the variables and loses, one at a time, the variable whose removal
leaves the smallest projection of the frontier. It stops when the projected
search from the frontier keeps all its sets within $T$ nodes, and the open sets
of each direction are then split again by its new heuristic. We remark that
this greedy relaxation by size is ours, while SymBA$^*$ relaxes the frontier
with abstraction strategies of its own. By [Lemma 18.4](#lem-18-4) the new
heuristics are consistent and admissible, with one caveat: since the
bidirectional search does not insert the states with $f \geq U$, a state whose
path to a target enters $C$ there may be left out of the frontier. Such a state
has no path cheaper than $U$, while the proofs of [Proposition
18.6](#prop-18-6) and [Proposition 18.7](#prop-18-7) need the consistency only
along an optimal path (whose states have $f \leq C^* < U$ as long as $U >
C^*$). Also, changing the heuristic between two steps keeps these proofs valid,
since they use only the current heuristic and the fact that the closed states
have their distances.

## 18.6 BDDA$^*$ {#sec-18-6}

In each direction we keep an open list of sets $O_{f,g}$ (the states reached
with cost $g$ and $f = g + h$), the closed set $C$ of the states expanded, and
the sets $R_g$ of the states reached with cost $g$, as well as the sequence of
the sets expanded with their $g$. At the start, the root (the initial state
$s$, or the goals backward) is split into the sets $O_{k,0}$ by the sets $H_k$
of the direction.

**Algorithm 18.1** (A step of BDDA$^*$). Take the nonempty set $O_{f,g}$ of
least $f$ and, among those, of least $g$, and let $F = O_{f,g} \setminus C$.
While $F \neq \emptyset$: (i) add $F$ to $C$; (ii) stop with cost $g$ if $F$
meets the target; (iii) compute $Y_c = \mathrm{img}(F, T_c) \setminus C$ for
each cost $c$, $T_c$ being the partitions of cost $c$; (iv) add $Y_c \cap H_k$
to $O_{g+c+k, g+c}$ for each $k$, except $Y_0 \cap H_{f-g}$, which becomes the
next $F$.

These are the buckets of $g$ and $h$ of
BDDA$^*$ [[191]](Z-bibliography.md#ref-191), close to the SetA$^*$ of Jensen,
Veloso and Bryant [[192]](Z-bibliography.md#ref-192). With $h = 0$ (the option
`heuristic` off, or a domain without terms) the search is the symbolic
Dijkstra's algorithm of the planners with
costs [[196]](Z-bibliography.md#ref-196). With all costs set to 1
(`unit_costs`, `symbolic_bfs()`) it is the breadth-first search, whose layers
are the sets $O_{g,g}$. The loop of step (iv) on the moves of cost 0 is the
breadth-first search within a bucket that such moves require: by consistency
they cannot lower $h$, and the successors with the same $h$ belong to the
bucket being expanded.

<a id="prop-18-6"></a>**Proposition 18.6** (Optimality of BDDA$^*$).
Let $h$ be consistent, zero on the target and infinite only on states from
which no target is reachable. Then a state is expanded at most once, a state
$x$ from which a target is reachable is expanded with $g = d(s, x)$, and the
search stops with $g = C^*$ if a target is reachable; otherwise, the open list
empties.

*Proof.* A state is expanded at most once, since it enters $C$ when it is
expanded and each set taken from the open list is reduced by $C$. Suppose that
a state from which a target is reachable is expanded with $g$ larger than its
distance. Let $x$ be the first one in time, expanded at $(f, g)$ with $g > d(s,
x)$, and let $P$ be a cheapest path from $s$ to $x$. When $x$ is expanded, let
$y$ be the first state of $P$ not in $C$ (a target is reachable from $y$
through $x$, and hence $h(y)$ is finite). The predecessor of $y$ on $P$ was
expanded before $x$, and therefore with its distance; hence, $y$ was inserted
in $O_{f_y, g_y}$ with $g_y = d(s, y)$ and, by consistency along the rest of
$P$, $f_y = g_y + h(y) \leq d(s, x) + h(x) < f$. Since $y \notin C$, that set
has not been expanded yet, which contradicts the choice of the set of least $f$
(the successors of the set being expanded have $f' \geq f$ by consistency).
Apply the same argument to a cheapest path to a target $t$ with $d(s, t) < g$
when the search stops at $g$: it gives a set of the open list with $f_y \leq
d(s, t) + h(t) = d(s, t) < g \leq f$, again a contradiction; hence the search
stops with $g = C^*$. It ends, since $C$ grows at each step and the states are
finite. $\square$

With an inconsistent heuristic the search still returns the cost of a path,
since the sets record the costs of "real" paths, but not necessarily an optimal
one: a state may be closed with a $g$ larger than its distance. This is why
`consistent()` is provided, to check a heuristic, typically before it is used.
The states marked as dead ends are in no $H_k$, and they are dropped when they
are reached.

**The path.** Each set $F$ expanded is recorded with its $g$, in order. If
$F_i$ (with cost $g_i$) holds the state $t$ that meets the target, we look for
a predecessor of $t$ in the sets $F_j$, $j < i$, from the last one, with $g_j =
g_i - c$ for some cost $c$ and $\mathrm{pre}(\{t\}, T_c) \cap F_j \neq
\emptyset$. We take one state of that intersection, and the walk repeats from
it until it reaches the root. The walk cannot fail, since each state of $F_i$
entered its bucket as a successor of a state of a set expanded before, with $g$
lower by the cost of the move; it costs one preimage of a single state per
partition and per cost tried (the counter `path_images`). Backward, the walk
takes images, and it gives the path from the initial state to a goal directly.

**List BDDA$^*$.** With the option `list`, the open list keeps one set $L_g$
per $g$, not split by $h$, as the List BDDA$^*$ of Torralba et al. (section 3.3
of [[193]](Z-bibliography.md#ref-193)). Before a step, the least $h$ of each
$L_g$ is found by intersecting it with $H_0, H_1, \ldots$ from a lower bound
kept for it, by increasing $g$ until $g$ reaches the least $f$ found. The set
$L_g \cap H_h$ of the least $(f, g)$ is then moved to the open sets $O_{f,g}$
and expanded as above. Since the set taken is that of the least $(f, g)$ over
all the states reached, as in BDDA$^*$, the two expand the same sets, and only
the number of intersections with the sets $H_k$ (the counter `splits`) differs,
as the states never expanded are never split.

## 18.7 The bidirectional search {#sec-18-7}

Both directions run as in [§18.6](#sec-18-6), the forward one by images from
the initial state with $h_F$, the backward one by preimages from the goals with
$h_B$. Each set $X$ inserted in the open list of a direction with cost $g$ (the
roots and the layers of cost 0 included) is intersected with the sets $R_{g'}$
of the other direction by increasing $g'$, and the first nonempty intersection
lowers the incumbent $U$ to $g + g'$ if it is smaller. This is the immediate
detection of the solutions of the explicit searches
([§7.2](07-bidirectional.md#sec-7-2)), done on sets. We do not insert a
successor with $f \geq U$, since no path through it costs less than $U$. Let
$f_D$ and $g_D$ be the least $f$ and the least $g$ of the nonempty sets in the
open list of the direction $D$, and $\epsilon$ the least cost of an operator.
The search stops when an open list is empty or when

$$
U \leq \mathrm{LB} = \max \{f_F, f_B, g_F + g_B + \epsilon\}
$$

(more precisely, the largest value of $\mathrm{LB}$ seen so far).

<a id="prop-18-7"></a>**Proposition 18.7** (The stopping rule).
Let $h_F$ and $h_B$ be consistent and admissible. When the search stops, $U =
C^*$ if a path exists, and $U = \infty$ otherwise.

*Proof.* At any time, $U$ is the cost of a path, since $g$ and $g'$ are the
costs of paths from $s$ to the state and from the state to a goal. Suppose that
$U > C^*$ at some time, and let $P = (s = p_0, \ldots, p_k)$ be an optimal path
to a goal. As in the proof of [Proposition 18.6](#prop-18-6), the first state
$p_i$ of $P$ not closed forward is in the forward open list with $g = d(s,
p_i)$ (it is not pruned, since its $f$ is at most $C^* < U$); likewise, the
last state $p_j$ not closed backward is in the backward open list with $g =
d(p_j, p_k)$. If $i \geq j$, the state $p_j$ has been inserted by both
directions with its optimal $g$, and the second of the two insertions met the
set $R_{g'}$ of the other direction that holds it, giving $U \leq C^*$, a
contradiction. Hence $i < j$, and $f_F \leq d(s, p_i) + h_F(p_i) \leq C^*$,
$f_B \leq C^*$ likewise, and $g_F + g_B + \epsilon \leq C^*$, since the segment
from $p_i$ to $p_j$ has at least one move. Thus $\mathrm{LB} \leq C^* < U$
whenever $U > C^*$, and also for the largest $\mathrm{LB}$ seen, since $U$ does
not increase; the stop on $U \leq \mathrm{LB}$ therefore gives $U = C^*$. If
the forward open list empties, all the states reachable from $s$ that are not
dead ends have been closed forward, among them the goal $p_k$ of $P$, which was
inserted with $g = C^*$ and met the set $R_0$ of the backward direction; hence
$U = C^*$, and symmetrically for the backward list. With no path, $U$ stays
$\infty$. $\square$

This is the argument of [Theorem 7.6](07-bidirectional.md#thm-7-6), with sets
in the place of nodes. Without the two bounds on $f$, the bound $g_F + g_B +
\epsilon$ is the stopping rule of the bidirectional Dijkstra's algorithm, and
that of the bidirectional symbolic search of
SymBA$^*$ [[193]](Z-bibliography.md#ref-193), whose heuristics are abstractions
computed during the search.

**The direction of a step.** Any interleaving of the two directions
satisfies [Proposition 18.7](#prop-18-7), and one can therefore choose it
freely. A step expands the direction whose next set is estimated "cheaper" to
expand (ties going forward); its estimate is the nodes of the set times the
ratio between the nodes of the image and those of the set at the last step of
that direction. Torralba et al. estimate instead the time of the next step from
the time of the previous one and from the sizes of the two sets (section 3.2
of [[193]](Z-bibliography.md#ref-193)); we use counts of nodes, and hence the
search and its counters depend neither on the load of the machine nor on the
workers. The path is the forward walk from a state of the intersection that set
$U$ (with its $g$ forward), reversed, followed by the backward walk from the
same state.

## 18.8 The $k$ cheapest plans {#sec-18-8}

We now turn to `top_k( k )`, which returns the $k$ cheapest plans, i.e., paths
from the initial state to a goal, by nondecreasing cost. As in the top-$k$
planning of Speck et al. [[195]](Z-bibliography.md#ref-195), a path may visit a
state more than once, and two plans differ when their sequences of states or
their costs do. We first compute, by preimages from the goals, the set $W$ of
the states from which a goal is reachable. For each $g$, the search keeps the
set $S_g$ of the states of $W$ reached from $s$ by a path of cost $g$, without
removing the states reached with a smaller $g$. Starting from $S_0 = \{s\} \cap
W$, it takes the least $g$ still open, records $S_g$, enumerates the plans of
cost $g$, and adds $\mathrm{img}(S_g, T_c) \cap W$ to $S_{g+c}$ for each cost
$c$. A plan of cost $g$ ending in $t \in S_g \cap G$ is enumerated backward:
the predecessors of a state $u$ with cost $g_u$ are the states of
$\mathrm{pre}(\{u\}, T_c) \cap S_{g_u - c}$ for each $c$, and the walk ends at
$s$ with cost 0.

<a id="prop-18-8"></a>**Proposition 18.8** (The $k$ cheapest plans).
Let every operator cost more than 0. Then $S_g$ is the set of the states of $W$
reached from $s$ by a path of cost $g$, all the plans of cost $g$ are
enumerated when $S_g$ is taken, each once, and the plans are returned by
nondecreasing cost.

*Proof.* The first claim is by induction on $g$. A state $y$ reached by a path
of cost $g > 0$ whose states are in $W$ (all its states reach a goal through
$y$) has a predecessor $x$ reached with cost $g - c$, $c > 0$, and $x$ is in
$S_{g-c}$, which is taken before $S_g$ and expanded into $S_g$; conversely,
$S_g$ holds only images of such states. When $S_g$ is taken, each $S_{g'}$ with
$g' < g$ is complete, since the costs are positive and the sets are taken by
increasing $g$. A plan of cost $g$ has each of its states $u$, with the cost
$g_u$ of its prefix, in $S_{g_u}$, and the backward walk therefore finds it.
Two branches of the walk differ in a state or in a cost, and give different
plans. $\square$

With moves of cost 0, $S_g$ would not be complete when it is taken, and a cycle
of cost 0 would give infinitely many plans of the same cost; `top_k()`
therefore refuses them (with `std::invalid_argument`). The search ends, since
either $k$ plans are found, or $W$ has no cycle through the reachable states
and the sets $S_g$ become empty, or $g$ exceeds its bound `max_cost`.

**Bidirectionally.** With the bidirectional direction, `top_k()` (through
`top_k_bidirectional()`) also keeps the sets $S^B_b$ of the states from which a
goal is reached by a path of cost $b$, and it extends the direction whose next
set has fewer nodes. Let the forward sets be complete up to $A$ and the
backward ones up to $B$. Before the plans of cost $c$ are enumerated, the steps
go on until $c \leq A$ or $c \leq A + B + 1$. A plan of cost $c$ is then split
after its longest prefix of cost at most $A$: the prefix goes forward through
the sets $S^F$, and the next state (reached by a move that brings the cost
above $A$) has a rest of cost at most $c - A - 1 \leq B$, which goes through
the sets $S^B$. Note that each plan has exactly one such split, and it is
therefore enumerated once, without the checks of duplicates that a meeting of
the two directions at any state would need.

## 18.9 The parallel images {#sec-18-9}

In Sylvan the operations are parallel inside: each recursive call of a
relational product, a conjunction or a quantification spawns its two cofactors
as tasks of the work-stealing scheduler Lace, and one product of a "large" set
may thus use all the workers. On top of that, the option `images` (a
`SymbolicImages`) chooses how the images of the partitions of a set are
computed. The choices are (i) `sequential`, one product after the other, (ii)
`lace`, all of them in one task of Lace, and (iii) `executor`, the partitions
dealt round-robin to the workers of an executor of ADS (the threads of the
library, OpenMP or FastFlow, as in [Chapter 8](08-shared-memory.md#ch-8)). With
`lace`, the task spawns the products of the two halves of the partitions
recursively and then joins them by a balanced tree of disjunctions (`images()`
and `disjoin()`); the scheduler typically balances products of quite different
sizes. With `executor`, each worker submits its products to Lace by
`image_into()`, which Lace accepts from several threads at once.

<a id="lem-18-9"></a>**Lemma 18.9** (Independence of the workers).
The costs, the paths and the counters of `SymbolicStats` other than `images_by`
do not depend on the choice of `SymbolicImages`, on the workers of the executor
or on those of Lace.

*Proof.* A reduced ordered BDD is canonical for its function and its order of
the variables [[163]](Z-bibliography.md#ref-163), and the order is fixed by the
encoding. Each set computed by the search is a function of the sets before it,
in whatever way its parts were computed. Hence so are its nodes, its states,
the choice of the direction (which depends on nodes only) and the state taken
from an intersection by the path (the first in the order of the variables).
Finally, the products are counted one per partition and per set, whoever
computes them. $\square$

In fact, the lemma also makes the counters a "fair" measure of the work of
these searches, since the nodes and the products are the same on any machine
while the time of a product is not. We remark that two points of the
implementation need care. First, the garbage collection of Sylvan drops the
nodes that no registered handle reaches, and a result received by a thread
outside Lace is not registered until it is stored in a `Bdd`. Since several
threads submit products at once with the executor, `image_into()` stores its
result in the registered handle of its output inside the task of Lace, before
the task ends. Second, Lace hands a task to a thief with memory fences that the
tools which check the synchronization of threads do not model. Hence, the tasks
of the library take the data of their job, and give it back when they end,
through a flag read with acquire and written with release; ThreadSanitizer then
sees the order that the fences give.

## 18.10 The domains {#sec-18-10}

In `SymbolicSlidingTile< W , H >`
([`SymbolicPuzzles.h:74`][SymbolicPuzzles.h:74]), the sliding-tile puzzle on $W
\times H$ cells ([§14.2](14-domains.md#sec-14-2)), the variable $c$ is the cell
$c$ and its value the tile in it (0 the blank); thus, a state takes $WH \lceil
\log_2 WH \rceil$ bits, i.e., 36 for the 8-puzzle and 64 for the 15-puzzle. The
move of the blank from the cell $b$ to the next cell $c$ is the operator with
the precondition $x_b = 0$ and the effects $x'_b = x_c$ and $x'_c = 0$, which
reads and writes two variables only. The heuristic is the Manhattan distance, a
term per cell, and the backward heuristic the same distance from the cells of
the start; the boards are not checked for parity, and an unsolvable instance
can thus be given. In `SymbolicPancake< N >`
([`SymbolicPuzzles.h:223`][SymbolicPuzzles.h:223]), the pancake problem
([§14.9](14-domains.md#sec-14-16)), the variable of each pancake is its
position (the inverse of the stack). The flip of the top $k$ pancakes then maps
each variable on itself, $x'_p = k - 1 - x_p$ for $x_p < k$ (an effect
`map_of`), and the gap heuristic is a sum of terms on the pairs of pancakes of
consecutive sizes, with one more for the largest pancake against the plate.
Both are `DecodableSymbolicDomain`s, whose states decode to those of
`SlidingTile` and `PancakePuzzle`.

Finally, `SymbolicSAS< Task >` ([`SymbolicSAS.h:67`][SymbolicSAS.h:67]) is a
SAS+ task [[133]](Z-bibliography.md#ref-133) ([§16.1](16-planning.md#sec-16-1))
as a `SymbolicDomain`, whose variables are those of the task and whose goal is
the partial assignment of the task. It takes the conditional effects of the
operators of the task, and, if the task has axioms (members `axioms` and
`axiom_layer`), it is a `SymbolicAxiomDomain` whose derived variables are those
of the layers of the translator; the task must outlive the domain.

## 18.11 The library of BDDs and the limits {#sec-18-11}

PBeFS interfaces an existing package of BDDs instead of implementing the
diagrams, and we compared the two that fit a C++ library distributed on three
systems, Sylvan [[194]](Z-bibliography.md#ref-194) and
CUDD [[197]](Z-bibliography.md#ref-197). Sylvan (Apache-2.0) builds with CMake
on Linux, macOS and Windows (MSYS2 and MSVC are in its continuous integration),
and its operations are parallel on Lace; CUDD (BSD-3-Clause) builds with the
autotools only, has no supported build on Windows, and computes each product on
one core, one manager per thread. Neither is packaged by vcpkg, Homebrew or
MSYS2 (as of November 2026). We therefore use Sylvan, whose parallel operations
are the side of symbolic search that concerns a library of parallel search.
Sylvan is optional: CMake finds an installed one with `find_package(sylvan)`,
or else downloads the commit pinned in [`CMakeLists.txt`][CMakeLists.txt]
(unless `PBEFS_FETCH_SYLVAN` is off). The makefiles take it from `SYLVAN_ROOT`;
without it, the library `PBeFS::Symbolic` and the cases `sym-*` of the tester
are left out.

Of course, the encodings are not necessarily (and perhaps rarely) the best, and
their choice matters a great deal. With the stack of pancakes as the state,
i.e., a variable per position, the flip of $k$ pancakes copies the position $i$
into the position $k - 1 - i$, and its relation then reads the positions in
order, and it must keep up to $k$ values at once. With 8 pancakes, the searches
in that encoding fill a table of $2^{24}$ nodes, while the relations of the
encoding of [§18.10](#sec-18-10) have a few nodes per variable. In the
partitions we follow the order of the operators. The order of the variables is
that of the domain, or one given, or one found with `SymbolicLayout::optimize`,
which moves one variable at a time to the place that most lowers the sum, over
the operators and the terms of the heuristic, of the squared span of their
variables in bits. It is a static order, chosen once (Sylvan has no dynamic
reordering), and it is "good" for that sum only. We find the $k$ cheapest plans
by a blind search, as in symk [[195]](Z-bibliography.md#ref-195), which needs
positive costs (cf. [Proposition 18.8](#prop-18-8)). Finally, a search that
fills more than `max_filled` of the largest table of nodes after a collection
stops with `Status::node_limit`, with the least $f$ of its open list as its
bound (the largest $\mathrm{LB}$ of [§18.7](#sec-18-7) for the bidirectional
search); a single product that fills the table ends the process, as Sylvan
does.

## 18.12 Notes and references {#sec-18-12}

The reduced ordered BDDs are those of
Bryant [[163]](Z-bibliography.md#ref-163), and the relational product with a
partitioned relation is that of symbolic model checking of Burch et
al. [[190]](Z-bibliography.md#ref-190). BDDA$^*$ is due to Edelkamp and
Reffel [[191]](Z-bibliography.md#ref-191), and Jensen, Veloso and
Bryant [[192]](Z-bibliography.md#ref-192) gave SetA$^*$ and the branching on
sets of states. Edelkamp and Kissmann [[196]](Z-bibliography.md#ref-196)
extended the symbolic searches to action costs, with Dijkstra's algorithm and
BDDA$^*$ on buckets and the symbolic PDBs. The symbolic planners of Torralba et
al. [[193]](Z-bibliography.md#ref-193), [[198]](Z-bibliography.md#ref-198) gave
the bidirectional search with abstractions chosen during the search, List
BDDA$^*$ and the partitioning of the relation by size, and the choice of the
direction by the estimated time of the next step, which [§18.7](#sec-18-7)
replaces by a count of nodes. The perimeter seeding the PDBs follows Felner and
Ofek [[199]](Z-bibliography.md#ref-199), and the top-$k$ planning Speck et
al. [[195]](Z-bibliography.md#ref-195), whose planner symk is blind as ours is.
The library of diagrams is Sylvan of van Dijk and van de
Pol [[194]](Z-bibliography.md#ref-194).

## 18.13 Exercises {#sec-18-13}

**Exercise 18.1** (`ctest -R PBeFS.sym-astar`). Show by an example of two cells
that the copy $x'_w = x_u$ of [§18.3](#sec-18-3) needs the high bits of $x'_w$
set to zero when $b_u < b_w$, and that [Lemma 18.1](#lem-18-1) fails without
them; the test compares BDDA$^*$ with A$^*$ on the 8-puzzle and on four
instances of the 15-puzzle.

**Exercise 18.2** (`ctest -R PBeFS.sym-astar`). Twice the Manhattan distance is
not consistent: find a move of the 8-puzzle on which condition (ii) of [Lemma
18.3](#lem-18-3) fails, and check that `consistent()` refuses it.

**Exercise 18.3** (`ctest -R PBeFS.sym-pdb`). Give two patterns of a SAS+ task
written together by one operator, and a state at which the sum of their PDBs
exceeds the optimal cost; then explain why their largest is still admissible
([Lemma 18.4](#lem-18-4)).

**Exercise 18.4** (`ctest -R PBeFS.sym-variants`). Show that BDDA$^*$ and List
BDDA$^*$ expand the same sets, and give a case in which the second makes fewer
intersections with the sets $H_k$.

**Exercise 18.5** (`ctest -R PBeFS.sym-topk`). Prove that the split of a plan
in [§18.8](#sec-18-8) after its longest prefix of cost at most $A$ is unique,
and show by an example that splitting at any state where the two directions
meet would enumerate some plans twice.

**Exercise 18.6** (`ctest -R PBeFS.sym-par-thread`). Run the bidirectional
search of the 8-puzzle with 1 and 4 workers of Lace, and check that the
counters other than `images_by` agree ([Lemma 18.9](#lem-18-9)); explain why
the choice of the direction would not be reproducible if it were made by times.

<!-- the links to the sources, written by srclinks.py -->

[Bdd.cpp]: https://github.com/dmeoli/PBeFS/blob/main/src/Bdd.cpp
[Bdd.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/Bdd.h
[Bdd.h:255]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/Bdd.h#L255
[Bdd.h:63]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/Bdd.h#L63
[CMakeLists.txt]: https://github.com/dmeoli/PBeFS/blob/main/CMakeLists.txt
[symbolic.cpp]: https://github.com/dmeoli/PBeFS/blob/main/test/symbolic.cpp
[symbolic.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic.h
[SymbolicDomain.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h
[SymbolicDomain.h:112]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L112
[SymbolicDomain.h:122]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L122
[SymbolicDomain.h:138]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L138
[SymbolicDomain.h:202]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L202
[SymbolicDomain.h:227]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L227
[SymbolicDomain.h:242]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L242
[SymbolicDomain.h:255]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L255
[SymbolicDomain.h:268]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L268
[SymbolicDomain.h:65]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L65
[SymbolicDomain.h:75]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicDomain.h#L75
[SymbolicModel.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicModel.h
[SymbolicModel.h:114]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicModel.h#L114
[SymbolicModel.h:74]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicModel.h#L74
[SymbolicPuzzles.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicPuzzles.h
[SymbolicPuzzles.h:223]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicPuzzles.h#L223
[SymbolicPuzzles.h:74]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicPuzzles.h#L74
[SymbolicSAS.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSAS.h
[SymbolicSAS.h:67]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSAS.h#L67
[SymbolicSearch.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h
[SymbolicSearch.h:148]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h#L148
[SymbolicSearch.h:188]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h#L188
[SymbolicSearch.h:206]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h#L206
[SymbolicSearch.h:248]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h#L248
[SymbolicSearch.h:93]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/symbolic/SymbolicSearch.h#L93
