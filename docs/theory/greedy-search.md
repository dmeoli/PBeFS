# Greedy best-first search and its parallelization

We are concerned here with the searches that give up the optimality of the
path for speed: the greedy best-first search (GBFS), which selects the node
of least $h$, the width-based searches, and the parallel versions of GBFS.
The notation is that of [Best-first search and A*](best-first-search.md);
the pseudocode of each search is on the page of its header. Unlike A*, for
which a parallel search must expand every node of $f < C^*$ and no other
one is needed, GBFS fixes no set of states that a correct parallelization
has to expand, and the question of how far a parallel GBFS may stray from
the sequential one is the subject of the second half of this page.

## The breaking of the ties and the benches

GBFS with a given breaking of the ties is deterministic, but different
breakings expand very different sets of states, and the analysis of
Heusner, Keller and Helmert [[112]](../../README.md#ref-112) characterizes
the union of these sets. The *high-water mark* of a state $s$ is
$\mathrm{hwm}(s) = \min_{p} \max_{s' \in p} h(s')$ over the paths $p$ from
$s$ to a goal ($\infty$ if there is none); it is the largest $h$ that any
GBFS must accept on its way from $s$, since every state of $h <
\mathrm{hwm}(s)$ reachable from $s$ without crossing a larger $h$ leads to
no goal. A state $s$ is a *progress state* if $\mathrm{hwm}(s) >
\min_{t \in \mathrm{succ}(s)} \mathrm{hwm}(t)$ (and so is a goal); the
*bench* of a progress state $s$ has the level $\ell = \min_t
\mathrm{hwm}(t)$, the inner states reachable from $s$ through
non-progress states of $h \leq \ell$, and the exits, the progress states
of $h = \ell$ among the successors of $s$ and of its inner states. The
*bench transition system* (BTS) is the set of the benches reached from
the bench of the initial state through their non-goal exits, and theorem 1
of [[112]](../../README.md#ref-112) states that a state is expanded by GBFS
with some breaking of the ties if and only if it lies in a bench of the
BTS. Since $\mathrm{hwm}(s) = \max( h(s) , \min_t \mathrm{hwm}(t) )$, the
high-water marks of a finite space are those of a shortest-path problem in
which a path costs the largest $h$ along it, solved by Dijkstra's algorithm
backward from the goals with $\max$ in place of $+$; the tester computes
the BTS in this way on random graphs and on the 8-puzzle, and checks it
against GBFS, FIFO and LIFO, before using it on the parallel searches.

## Parallel GBFS with a shared OPEN list

All the parallel searches of `ParallelGreedy.h` take the best node of the
shared OPEN only under a condition, and insert the successors of a node all
together, so that they are ordered among themselves as GBFS orders them
(algorithm 1 of [[111]](../../README.md#ref-111)). Let $s_j$ be the node that
the worker $j$ expands, $h(s_j) = \infty$ for an idle worker, and
$\underline{h} = \min_j h(s_j)$. KPGBFS [[106]](../../README.md#ref-106) has
no condition; Kuroiwa and Fukunaga [[108]](../../README.md#ref-108) show that
it is *pathological*, i.e., that for every $t$ there is an instance on which
it expands more than $t$ times the states of GBFS, and the simplest witness
is a state outside the BTS, which a worker takes while another one expands
the state that GBFS would take: the tester builds such an instance and
checks that KPGBFS expands a state outside the BTS with two workers in
synchronous rounds.

PUHF [[108]](../../README.md#ref-108) takes only the nodes marked *certain*,
those whose membership in the BTS is proven: a successor $s'$ of a certain
node $s$ with $h(s') = \min_{t \in \{s\} \cup \mathrm{succ}(s)} h(t)$ is
certain, since GBFS with some breaking of the ties expands it right after
$s$ (theorem 2 there), and so is the best node of OPEN when no worker is
expanding, which is then exactly the choice of GBFS. PUHF2
[[109]](../../README.md#ref-109) marks certain the best node of OPEN when
$h \leq \underline{h}$: if such a node $n$ were outside the BTS, a path to
a goal of high-water mark below $h(n)$ would leave the bench of its parent,
and a state of that path, of $h < h(n)$, would lie in OPEN or be in
expansion (lemma 1 there), against the choice of $n$. PUHF3 marks certain
all the nodes of OPEN of $h = \underline{h}$ after each removal, since any
of them could have been chosen in place of the one removed, and PUHF4 adds
the criterion of PUHF and prefers the certain nodes among those of equal
$h$; these are all BTS-constrained. They are not bounded, though: section 3
of [[110]](../../README.md#ref-110) shows two benches entered at once by two
workers, both in the BTS but never both explored by one GBFS, and nests the
construction to make the ratio to GBFS unbounded; the tester builds the
first level (an initial state with two successors of equal $h$, each the
head of a line of $x$ states of a smaller $h$ leading to its own goal) and
checks that PUHF with two workers expands both lines while GBFS expands
one.

OBAT [[110]](../../README.md#ref-110) takes the best node of OPEN when $h
\leq \underline{h}$, as PUHF2, and calls *b-state* a node with a successor
of smaller $h$, the only kind that can be a progress state; a b-state does
not insert its successors but enters DEFERRED, ordered by $h$, and its
successors enter OPEN only when it leaves DEFERRED, which happens when the
best $h$ of DEFERRED is at most that of OPEN and at most $\underline{h}$.
Let a node be *committed* when it leaves OPEN or DEFERRED for the last
time. Theorem 2 of [[110]](../../README.md#ref-110) states that the committed
nodes, in the order of their commits, are an order in which GBFS may expand
them: when $s_{i+1}$ is committed, the nodes generated by $s_1 , \ldots ,
s_i$ and not committed are in OPEN, in DEFERRED, in expansion as b-states,
or successors of the a-states in expansion, which have no successor of
smaller $h$ than themselves, and the conditions above make $h(s_{i+1})$
the least among all of them. Theorem 3 bounds the nodes left in DEFERRED at
the end by $k |p|$, $k$ the workers and $p$ the path: at most $k$ of each
$h$ (two of equal $h$ from one worker would mean a removal from OPEN while
DEFERRED held a node of the same $h$), and only the values of $h$ of the
states of the path. Hence OBAT expands at most $N + k |p|$ states, $N$ the
expansions of GBFS with some breaking of the ties (theorem 4). The tester
replays the commits of OBAT as a GBFS, checking at each step that the
committed node has the least $h$ of the simulated OPEN, and counts the
nodes of each $h$ left in DEFERRED, on random graphs, on the 8-puzzle and,
under real threads, on the three executors.

The separate generation and evaluation (SGE) of
[[111]](../../README.md#ref-111) lets any idle worker evaluate a successor
generated by another one, the node staying in expansion until all its
successors are evaluated and inserted; since the conditions above compare
only the $h$ of the nodes in OPEN, in DEFERRED and in expansion, SGE keeps
all the guarantees, and the tester checks them with SGE too. In synchronous
rounds the $k$ workers are emulated in one thread, every unit of work
taking one round, which makes the counts of a parallel search independent
of the machine: the bounds above hold for every interleaving, so they hold
for this one, and the witnesses of the pathologies are replicable.

## Width-based search

The *novelty* of a state is the size of the smallest tuple of its
features that no state generated before it made true. IW(k) is a
breadth-first search that prunes the states of novelty larger than $k$;
it generates at most $O(n^k)$ states for $n$ features, since each state
kept makes a new tuple of at most $k$ features true, and it reaches a goal
of a problem of width at most $k$ by an optimal path when the costs are
unit (theorem 4 of [[104]](../../README.md#ref-104)). Since the tuples of the
parent of a state have been recorded, only the tuples with a feature that
the move adds can be new, which bounds the work of a novelty test by the
tuples of the features changed by the move. On a grid whose features are
the column and the row of the cell, a cell is a tuple of two features, so
that IW(2) never prunes and is a breadth-first search, which the tester
checks against A*.

## Partial expansion

PEA* and EPEA* keep with each node a stored value $F(n) \geq f(n)$ and put
in OPEN, at a selection of $n$, only the successors of $f = F(n)$ (or of
$f \leq F(n) + C$), sending $n$ back to OPEN with the least $f$ of the
others. The invariant of A* holds with $F$ in place of $f$: every node of
OPEN has $F(n)$ at most the $f$ of each successor of $n$ not yet put in
OPEN, so that with an admissible $h$ the least $F$ of OPEN is a lower bound
of the cost of the paths not yet explored and the first goal selected is
optimal [[114]](../../README.md#ref-114). EPEA* generates only the successors
it puts in OPEN, by an operator selection function that knows the
$\Delta f$ of each move without making it; on the 15-puzzle with the
Manhattan distance, where $\Delta f \in \{ 0 , 2 \}$, it generates and
stores about half of the nodes of A*, at the price of expanding again the
nodes with successors of both values.
