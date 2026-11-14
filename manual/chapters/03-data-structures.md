# 3. The data structures of search {#ch-3}

[Source: `include/PBeFS/astar/OpenList.h` (`OpenEntry`, `OpenKind`,
`HeapOpenList`, `BucketOpenList`), `include/PBeFS/astar/NodeTable.h`,
`include/PBeFS/astar/Node.h` (`Node`, `NodePool`, `packed_fg`),
`include/PBeFS/astar/Structures.h`, `include/PBeFS/astar/structures/`
(the sets `std`, `boost` and the default one), `include/PBeFS/astar/AStar.h`;
the tests `open`, `astar` and `ddd` of `test/test.cpp`;
`tools/bench/structures.sh`]

We are concerned in this chapter with the data structures that a
best-first search runs on: the OPEN list, i.e., a priority queue of the
entries $(f, g, u)$ of [§2.3](02-best-first.md#sec-2-3), the CLOSED
list, i.e., a dictionary from the generated states to their nodes, and
the memory of the nodes themselves. On the domains whose expansion takes
a few instructions, such as the sliding-tile puzzles and the grids,
these structures take most of the time of a search
([§2.10](02-best-first.md#sec-2-10)), and the choices made here decide
the speed of all the searches of the library, the parallel ones
included. We first list the operations that a search performs on them
and the lazy deletion that replaces the decrease of a key. Then we
describe the heaps, the two-level bucket queue and the packed keys, and
we prove what the breaking of the ties changes. In the second half of
the chapter we treat CLOSED, i.e., the hash table with open addressing,
the arena of the nodes, the prefetch of the lookups, the perfect hashes
of the domains, and the two variants of duplicate detection that the
memory-bounded searches of [Chapter 4](04-memory-bounded.md#ch-4) use.
Finally, we describe the three sets of data structures among which a
build of the library chooses.

## 3.1 The operations of a search {#sec-3-1}

An iteration of A$^*$ (Algorithm 2.2 of
[§2.3](02-best-first.md#sec-2-3)) performs the following operations on
its data structures: (i) one extraction of the first entry of OPEN, (ii)
for each of the $b$ successors of the node, one evaluation of the hash
of its state and one lookup in CLOSED, which inserts the state if it is
new, and (iii) for each successor that is new or reached through a
cheaper path, one push in OPEN. In the library, the three structures
have the following interfaces, which the three sets of
[§3.12](03-data-structures.md#sec-3-12) implement in the same way.

* *OPEN* (`OpenList.h`): `push(node)` pushes the node with its current
  $f$ and $g$, `push(node, g)` pushes it with the given $g$ (and
  $f = g + h$) without reading `node->g`, `top()` is the entry with the
  least $f$ and, among those, the largest $g$, `top_f()` is its $f$,
  `pop()` removes and returns it, and `empty()` and `size()` count the
  entries, the stale ones included.
* *CLOSED* (`NodeTable.h`): `find_or_insert(state, hash, make)` looks
  the state up and, if it is absent, inserts the node returned by
  `make()`, returning the node of the state and whether it has just
  been inserted; `find(state, hash)` returns the node or `nullptr`;
  `prefetch(hash)` brings into the cache the first place probed for
  the hash; `slots_for(n)` is the number of places that holds $n$ nodes
  without growing; `size()` and `bytes()` measure the table.
* *The arena of the nodes* (`NodePool`, `Node.h`): `allocate()` returns
  a new node, whose address never changes until the end of the search;
  the nodes are never freed one by one.

The hash is computed once, by the search, from the `hash` of the domain,
and passed to each operation of CLOSED; hence neither the lookup, nor
the growth of the table, nor the prefetch evaluates it again.

### Lazy deletion

Neither OPEN list of the library can lower the key of an entry. A node
reached through a cheaper path is pushed again with its new $g$, and the
old entry, which records the $g$ of its push, becomes *stale*. Then the
search drops it when it reaches the top (line 6 of Algorithm 2.2). An
entry is the pointer to the node and the value of $g$ at the push, and
it is stale when its node is closed or has another $g$:

```cpp
--8<-- "snippets/data-structures-entry.cpp"
```

This is the *lazy deletion* of the textbooks. It makes the entries plain
values, without the handles by which a heap would find an entry to move
it, and which each swap of the heap would have to update. The price is
in the number of entries, which the next lemma bounds.

<a id="lem-3-1"></a>**Lemma 3.1** (The entries of OPEN). The entries
pushed by A$^*$ are $1 + N + I$, with $N$ the nodes stored besides the
root and $I$ the improvements of $g$ at line 14 of Algorithm 2.2; for any
heuristic they are at most the successors generated plus one, and with a
consistent one the improved nodes are all still in OPEN.

*Proof.* Line 3 pushes the root, line 13 pushes each new node once, and
line 16 pushes once per improvement. Each successor generated is either
new, or an improvement, or a duplicate, and hence $N + I$ is at most the
number of the successors generated. With a consistent heuristic, moreover,
no closed node is improved ([Theorem
2.13](02-best-first.md#thm-2-13)), and the improvements are those of
nodes still in OPEN. $\square$

The stale entries are thus at most $I$, and they cost one extraction
each, without an expansion. On the sliding-tile puzzles with the
Manhattan distance and on the grids with unit costs, a node of OPEN is
rather rarely improved, since its first path is usually a shortest one, and
the stale entries are typically a small fraction of the pushes.

## 3.2 Heaps {#sec-3-2}

A *binary heap* keeps $N$ entries in an array, the entry in the place
$i$ being not worse than its children in the places $2i + 1$ and $2i +
2$; the first entry is in the place 0. A *$d$-ary heap* has $d$ children
per node, in the places $d i + 1, \ldots, d i + d$, and the parent of
the place $i > 0$ is the place $\lfloor (i - 1) / d \rfloor$ (cf.
Johnson [[158]](Z-bibliography.md#ref-158)). A push puts the new entry in
the place $N$ and moves it up, swapping it with its parent as long as it
is better ([Figure 3.1](03-data-structures.md#fig-3-1)). A pop takes the
entry of the place 0, moves the last entry there and moves it down,
swapping it with the best of its children as long as that child is
better ([Figure 3.2](03-data-structures.md#fig-3-2)).

<a id="prop-3-2"></a>**Proposition 3.2** (The cost of a $d$-ary heap).
On a $d$-ary heap of $N$ entries, a push takes at most
$\lceil \log_d ((d - 1) N + 1) \rceil$ comparisons and a pop at most
$d$ times as many, i.e., $O(\log_d N)$ and $O(d \log_d N)$.

*Proof.* The levels $0, \ldots, \ell$ of a $d$-ary heap hold
$(d^{\ell + 1} - 1) / (d - 1)$ places, and hence a heap of $N$
entries has at most $\lceil \log_d ((d - 1) N + 1) \rceil$ levels. A
push makes one comparison per level that the entry climbs, and a pop
makes, per level that the entry descends, $d - 1$ comparisons to find
the best child and one with the entry. $\square$

In the default set, `HeapOpenList` is a $d$-ary heap with $d = 4$. Its
items carry $f$ and $g$ next to the pointer to the node, and hence a
comparison never follows the pointer; they are ordered by the least $f$
and, for equal $f$, by the largest $g$:

```cpp
--8<-- "snippets/data-structures-heap-order.cpp"
```

With $d = 4$ the heap has half the levels of a binary heap; a push makes
half the comparisons, and a pop makes $4$ comparisons per level against
$2$, on half the levels, i.e., about as many in all. Furthermore, with
integer costs an item takes
16 bytes, and the four children of a place are 64 consecutive bytes; with
the root in the place 0 they start 16 bytes after the beginning of a line,
and hence they straddle two lines of cache (an array shifted by three
places would align them on one). The descent of a pop thus reads at most
two lines per level, on half the levels of a binary heap. Also, the
default heap of the
library is allocated, when large, on huge pages
([§3.7](03-data-structures.md#sec-3-7)), and the walk from the root to a
leaf, whose places are far apart, does not miss the translation
lookaside buffer at each level. In the set `std` the heap is the binary
heap of `std::priority_queue`, and that of the set `boost` is the 4-ary
`boost::heap::d_ary_heap`. All three order the entries in the same way,
but for the ties of $f$ and $g$.

![Figure 3.1. A push in a 4-ary heap of entries $(f, g)$, as a tree and
as an array (the children of the place $i$ are the places $4i + 1$ to
$4i + 4$). (a) The new entry $(5, 4)$ goes to the place 8, a child of
the place 1. (b) It moves up twice, past $(6, 3)$, which has a larger
$f$, and past $(5, 2)$, which has the same $f$ and a smaller
$g$.](../figures/data-structures-heap-push.svg){#fig-3-1}

![Figure 3.2. A pop from the heap of Figure 3.1. (a) The first entry
$(5, 4)$ is returned and the last one, $(6, 3)$, fills the root. (b) It
moves down past the best of the four children of the root, $(5, 2)$,
and it stops at the place 1, being better than its children; the
places read are those of two lines of
cache.](../figures/data-structures-heap-pop.svg){#fig-3-2}

## 3.3 The two-level bucket queue {#sec-3-3}

With integer costs and a bounded range of $f$, a priority queue can do
without comparisons: the queue of Dial [[65]](Z-bibliography.md#ref-65)
keeps an array of buckets indexed by the key, and it finds the least key
by a scan that moves forward. The order of A$^*$ has two keys, the least
$f$ and then the largest $g$. Accordingly, `BucketOpenList`, the default
OPEN for integer costs, nests two levels of buckets, as the nested bucket
queue of Burns et al. [[4]](Z-bibliography.md#ref-4) (which they
distinguish from the 2-level bucket queue of
[[10]](Z-bibliography.md#ref-10), meant to shrink the array when the
values of $f$ are many): an array of *layers*
indexed by $f$, each of them a set of LIFO lists indexed by $g$ ([Figure
3.3](03-data-structures.md#fig-3-3)).

* *Push* $(f, g, u)$: the entry goes on top of the list $g$ of the
  layer $f$, and the bound of the least nonempty layer is lowered to
  $f$ if $f$ is smaller.
* *Top and pop*: the scan moves forward from that bound to the first
  nonempty layer, and within it to the nonempty list of the largest
  $g$, from a bound that only moves down between two pushes of a
  larger $g$; the entry is the top of that list, i.e., the last one
  pushed with that $f$ and $g$.

The items lie in one vector of cells, each with the index of the next
cell of its list, and the cells freed by a pop form a free list that a
push takes from first; hence a push or a pop calls no allocator once the
vector has reached the largest size of the queue.

<a id="prop-3-3"></a>**Proposition 3.3** (The cost of the bucket
queue). If the least $f$ of the entries never decreases between two
pops (as with a consistent heuristic, by
[Theorem 2.13](02-best-first.md#thm-2-13)), the layers skipped by all
the pops of a search are at most $f_{\max} + 1$, with $f_{\max}$ the
largest $f$ pushed; a push costs $O(1)$ amortized time, and a pop
$O(1)$ amortized time plus the lists of its layer that it skips.

*Proof.* Under the hypothesis, the bound of the least layer only moves
forward, and it crosses each of the $f_{\max} + 1$ layers at most once.
The cost of the scan over all the pops is thus at most $f_{\max} + 1$,
plus one per pop. A push writes a cell and a head, and the vectors grow
by doubling, which costs $O(1)$ amortized per element. $\square$

Without the hypothesis, as in a worker of HDA$^*$, whose least layer
empties and fills again each time the nodes sent by the other workers
arrive, the bound moves back, and the scan may cross a layer many times.
In that case the cost is that of the layers between the new bound and the old
one. For the same reason, the lists of a layer need care. An array of
heads indexed by $g$ costs $\Theta(G)$ to rebuild, with $G$ the largest
$g$ of the layer, and on a grid $g$ grows to thousands while a layer
holds few entries. The default queue therefore keeps the heads of the
lists of $g < D = 256$ in a dense array (the first 4 in the layer), and
those of $g \ge D$ as a sorted vector of at most $R = 64$ runs $(g,
\text{head})$, the last of which is the top. A layer that receives a
65-th distinct $g \ge D$ turns dense on all its values of $g$ until it
empties. A filling of a layer with $m$ entries then costs $O(m + D)$ if
they have at most $R$ distinct values $g \ge D$, and $O(m + G)$
otherwise, where the $O(G)$ is paid by more than $R$ pushes, whereas a
dense array alone would cost $\Theta(G)$ at each filling, even of one
entry. Of course, the memory of a bucket queue grows with the range of
the keys rather than with the entries, and the queue is therefore used
only for integer costs. The sets `std` and `boost`, which have no bucket
queue, use the heap for all the cost types.

![Figure 3.3. The two-level bucket queue: a layer per value of $f$,
and in each layer a LIFO list per value of $g$ (the light places are
empty lists). (a) The push of $(5, 4, x)$ puts $x$ on top of the list of
$f = 5$ and $g = 4$. (b) The pop returns the entry of the least $f$,
then of the largest $g$, then the last pushed, i.e., $x$; the next pops
return $c$, $b$, $a$ and then $y$, the top of the layer
6.](../figures/data-structures-buckets.svg){#fig-3-3}

The test `open` checks both OPEN lists against an ordered set of the
pairs $(f, -g)$ on random sequences of pushes and pops in which the
least $f$ moves back below the $f$ of the last pop, as in a worker of
HDA$^*$, with $g$ up to 200 and up to 2000, i.e., on both sides of $D$.

> **Status — planned.** The *radix heap* of Ahuja et al.
> [[160]](Z-bibliography.md#ref-160), which needs only that no key pushed
> be smaller than the last key popped (the case of A$^*$ with a
> consistent heuristic) and costs $O(\log C)$ amortized time per
> extraction for keys in a range of $C$, is not in the library; the
> bucket queue covers the integer costs of its domains, and the heap the
> others.

## 3.4 Packed keys {#sec-3-4}

The parallel searches whose OPEN is shared or relaxed (MultiQueue A$^*$
and PA$^*$SE of [Chapter 8](08-shared-memory.md#ch-8)) compare the keys
of the entries as single machine words, and a comparison, or an atomic
read of a key, is then one instruction. For integer costs, `packed_fg`
(`Node.h`) puts $f = g + h$ in the high half of a word of 64 bits and
the complement of $g$ in the low one:

```cpp
--8<-- "snippets/data-structures-packed.cpp"
```

<a id="lem-3-4"></a>**Lemma 3.4** (The order of the packed keys). Let $0
\le g_i$, $0 \le h_i$ and $g_i + h_i < 2^{32}$ for $i = 1, 2$, and let
$k_i$ be the packed key of $(g_i, h_i)$. Then $k_1 < k_2$ if and only
if $f_1 < f_2$, or $f_1 = f_2$ and $g_1 > g_2$.

*Proof.* The key is $k_i = 2^{32} f_i + (2^{32} - 1 - g_i)$, since $0 \le
g_i \le f_i < 2^{32}$ (here $h_i \ge 0$ is used) and the complement of
$g_i$ in 32 bits is
$2^{32} - 1 - g_i$. Note that the low half lies in $[0, 2^{32})$. Hence the
comparison of the keys is that of the high halves, and for equal high
halves that of the low ones, i.e., the comparison of $f$ and then the
reverse comparison of $g$. $\square$

[Figure 3.4](03-data-structures.md#fig-3-4) shows two keys of the same
$f$: the one of the larger $g$ is the smaller word, and comes first.

![Figure 3.4. The packed keys of two entries with $f = 7$, in
hexadecimal: the high half holds $f$, the low half the complement of
$g$, so that the entry with $g = 5$ has the smaller word and comes
before the one with $g = 3$, as in the order of
A$^*$.](../figures/data-structures-packed.svg){#fig-3-4}

## 3.5 Breaking the ties {#sec-3-5}

Both OPEN lists break the ties of $f$ in favor of the larger $g$, i.e.,
of the smaller $h$, since $g = f - h$ when $f$ is fixed; ordering by
the least $h$ among the least $f$, as some textbooks put it, is the same
rule. We now make precise what the choice changes, with a consistent
heuristic and the notation of [§2.6](02-best-first.md#sec-2-6).

<a id="prop-3-5"></a>**Proposition 3.5** (The ties change only the
last layer). With a consistent heuristic and any rule for the ties,
A$^*$ expands the same nodes with $g^*(u) + h(u) < C^*$ and none with
$g^*(u) + h(u) > C^*$; the rules differ only in the nodes of the last
layer, $g^*(u) + h(u) = C^*$, that they expand before a goal.

*Proof.* The claim is [Theorem 2.14](02-best-first.md#thm-2-14), whose
proof does not depend on the order of the entries with equal keys.
$\square$

<a id="prop-3-6"></a>**Proposition 3.6** (The ties with a perfect
heuristic). Let the costs be positive and $h = h^*$. With the ties of
$f$ broken in favor of the larger $g$, A$^*$ expands the nodes of one
optimal path from $s_0$ to a goal, and no other node.

*Proof.* Since $h = h^*$, a node $v$ has $g^*(v) + h(v) = C^*$ if and
only if it lies on an optimal path, and $g^*(v) + h(v) > C^*$ otherwise;
the heuristic is consistent (it is the exact distance), and by [Theorem
2.14](02-best-first.md#thm-2-14) no node of the second kind is expanded.
We show by induction that each node $x$ selected before the goal lies on
an optimal path and that, when it is selected, its entry has the largest
$g$ among the entries of key $C^*$. This holds for $s_0$. When $x$ is
expanded, with $g(x) = g^*(x)$ by [Theorem
2.13](02-best-first.md#thm-2-13), each successor $y$ of $x$ on an
optimal path through $x$ gets an entry of key $C^*$ and $g(y) = g(x) +
c(x, y) > g(x)$, which is larger than the $g$ of all the other entries
of key $C^*$; no key is smaller than $C^*$, and hence the next selected
entry is one of these, which lies on an optimal path. The selected nodes
thus form a path from $s_0$, each the successor of the previous one on
an optimal path, and it ends at the first goal selected. $\square$

In words, with the ties in favor of the larger $g$ the search dives
along the last layer, and with an exact heuristic it expands only the
$d + 1$ nodes of the path it returns, $d$ being its number of edges; an
informative heuristic, exact on part of the last layer, gives a similar
dive on that part. Conversely, the opposite rule does the reverse: since
$h \ge 0$, a goal $t$ has the largest $g$ of the last layer, $g(t) =
C^*$, and with the ties in favor of the smaller $g$ it is selected only
when no other entry of key $C^*$ is left in OPEN, i.e., after all the
nodes of the last layer that is generated before it has been expanded.
On the grid of [Figure 2.6](02-best-first.md#fig-2-6) the two rules
expand 8 and 31 of the 31 cells of the last layer; on the domains with
unit costs and few values of $f$, such as the sliding-tile puzzles, the
last layer may hold a large part of the nodes, and the experiments of
Burns et al. [[4]](Z-bibliography.md#ref-4) report that the larger $g$
is the better rule. Among the entries with the same $f$ and $g$, the
bucket queue takes the last one pushed, which follows the most recent
path, while the heaps take any of them; this last order differs among
the three sets of [§3.12](03-data-structures.md#sec-3-12), which
therefore may expand different nodes of the last layer, and only of it.

## 3.6 CLOSED: a hash table with open addressing {#sec-3-6}

CLOSED is a dictionary from the generated states to their nodes, with no
removal, since a search never forgets a state. The default `NodeTable`
is a hash table with *open addressing* and *linear probing*: its places
are an array of $2^k$ cells, each holding the pointer to a node and the
64-bit hash of its state, and a state is looked for from its first place
onward, along consecutive places, up to its place or an empty one
([Figure 3.5](03-data-structures.md#fig-3-5)). Two choices make a probe
cheap. The hash is stored, and hence a probe compares two states only
when their hashes are equal, and a growth moves the places without
touching the nodes. Also, the places are consecutive in memory, and a
probe sequence reads one or two lines of cache instead of following the
pointers of a list.

The first place of a hash $x$ in a table of $2^k$ places is

$$
\left\lfloor \frac{(x \cdot \lfloor 2^{64} / \phi \rfloor) \bmod 2^{64}}
{2^{64 - k}} \right\rfloor ,
\qquad \phi = \frac{1 + \sqrt{5}}{2} ,
$$

i.e., the top $k$ bits of the product of $x$ by the integer part of
$2^{64} / \phi$, which is the multiplicative hashing of Knuth with the
golden ratio (*Fibonacci hashing*, cf. section 6.4 of
[[157]](Z-bibliography.md#ref-157)). However, the low bits of $x$ are
often regular: in a worker of HDA$^*$, for instance, all the states have
the same $x \bmod p$, and a table that took the low $k$ bits of $x$
would use one place in $p$. In contrast, the product spreads all the
bits of $x$ into the top bits, and the states of a worker fill its whole
table.

The table doubles when its load factor $\alpha$ would exceed $7/10$.
Thus `slots_for(n)`, i.e., $10 \lfloor n / 7 \rfloor + 16$, is the
number of places that holds $n$ nodes without growing, and
`SearchOptions::expected_nodes` sizes the table at the start of a search
([§2.3](02-best-first.md#sec-2-3)). We recall the classical analysis of
linear probing for the cost of a lookup.

<a id="prop-3-7"></a>**Proposition 3.7** (The cost of CLOSED). With hash
values independent and uniformly distributed, and a load factor $\alpha <
1$, a successful lookup
inspects about $\frac{1}{2}\bigl(1 + \frac{1}{1 - \alpha}\bigr)$
places on average and an unsuccessful one about
$\frac{1}{2}\bigl(1 + \frac{1}{(1 - \alpha)^2}\bigr)$; with
$\alpha \le 7/10$ these are at most about 2.2 and 6.1 places. The
insertions, with the doublings, cost $O(1)$ amortized time each.

*Proof.* The two averages are the classical results on linear probing
(section 6.4 of [[157]](Z-bibliography.md#ref-157)), and at $\alpha =
7/10$ they are $\frac{1}{2}(1 + 10/3) \approx 2.17$ and $\frac{1}{2}(1 +
100/9) \approx 6.06$. A doubling from $2^k$ places moves at most
$\frac{7}{10} 2^k$ nodes and follows at least $\frac{7}{20} 2^k$
insertions since the previous doubling; its cost, charged to them, is at
most 2 moves per insertion. $\square$

The unsuccessful lookups are those of the new states, i.e., about one
per node stored, and the successful ones those of the duplicates. A probe
sequence of 6 places of 16 bytes spans 96 bytes, i.e., two or three lines
of cache.
Indeed, the insertion is done by the same probe that finds the state absent
(`find_or_insert`), with one inspection only. In the set `std`, CLOSED
is a `std::unordered_map` (a table of chains) whose key is the state
together with its hash, which the search computed; the map thus does not
evaluate the hash again:

```cpp
--8<-- "snippets/data-structures-std-key.cpp"
```

```cpp
--8<-- "snippets/data-structures-find-or-insert.cpp"
```

Its `prefetch` does nothing, since a table of chains cannot prefetch a
bucket. In the set `boost`, CLOSED is `boost::unordered_flat_map`, a
table with open addressing, reserved for $7/10$ of its places as the
default one.

![Figure 3.5. Linear probing in CLOSED, whose places hold a hash
(here 3 hexadecimal digits) and a pointer to a node. (a) The insertion
of $v$, whose first place is 5: the places 5, 6 and 7 hold other hashes,
compared without reading their states, and $v$ goes into the empty place
8. (b) The lookup of $w$, of the same first place: the place 5 has the
same hash and another state, which costs one comparison of states, and
the place 7 holds $w$.](../figures/data-structures-probe.svg){#fig-3-5}

## 3.7 The arena of the nodes {#sec-3-7}

A node holds the state, $g$, $h$, the pointer to its parent, the flags
*closed* and *expanded* ([§2.3](02-best-first.md#sec-2-3)), the
operator that generated it and the key of its owner in the parallel
searches that need one:

```cpp
--8<-- "snippets/data-structures-node.cpp"
```

With a state of one word and integer costs, a node takes 32 bytes, the
operator and the key fitting in the padding after the flags. The
infinite cost of a cost type, which stands for no incumbent and for the
bound of a search not yet cut, is `infinite_cost()` (`Node.h`), i.e.,
the infinity of a floating-point type and the largest value of an
integer one. The nodes are allocated from an *arena* (`NodePool`): a
block of $2^{16}$ nodes is taken from the system when the previous one
is full, and a node is allocated by advancing a pointer in the block,
without a call to the allocator. Nodes are never freed one by one, but
all together at the end of the search ([Figure
3.6](03-data-structures.md#fig-3-6)). The arena has some consequences.
The address of a node never changes, and the pointers to it (from
CLOSED, from the entries of OPEN and from its children as their parent)
stay valid until the end; hence CLOSED and OPEN can store pointers, and
neither moves a node. Also, the nodes made one after another are next to
each other in memory, and the children of a node, which are usually made
together, often share lines of cache; furthermore, a node costs no
header of the allocator. In the set `std` the nodes are kept in a
`std::deque`, whose `emplace_back` does not move the nodes already there
either, and the set `boost` in a `boost::object_pool`:

```cpp
--8<-- "snippets/data-structures-std-pool.cpp"
```

The blocks of the arena, the places of CLOSED and the array of a large
heap of the default set are requested from the operating system as
*huge pages* where possible: a request of at least 2 MB is aligned to 2
MB and, on Linux, advised to the kernel for the transparent huge pages,
which the kernel may ignore. A structure whose places are reached in no
useful order, as those of CLOSED by hashes, then misses the translation
lookaside buffer once per 2 MB instead of once per 4 KB page, and its
memory is faulted in a huge page at a time.

![Figure 3.6. The arena of the nodes, in chunks: each node points to
its parent (the arcs above), the places of CLOSED and the entries of
OPEN point to the nodes, and no node ever moves; a chunk that is full is
followed by a new one, from which the next node ($t$) is
allocated.](../figures/data-structures-arena.svg){#fig-3-6}

## 3.8 Prefetching the lookups {#sec-3-8}

On a large search, CLOSED is much larger than the cache, and the place
of a successor lies at an address that the hash makes random; almost all
the lookups then miss the cache, and the time of an iteration is
dominated by $b$ misses. For this reason, the sequential A$^*$ splits
the expansion in two passes (the excerpt of
[§2.3](02-best-first.md#sec-2-3)): the first generates all the
successors, computes their hashes and calls `prefetch(hash)` on each of
them, which issues a prefetch of the first place of the hash into the
cache without waiting for it. The second does the lookups, by which time
the lines are on their way or already there. Thus, the misses of the $b$
lookups overlap ([Figure 3.7](03-data-structures.md#fig-3-7)), and the
iteration waits for about one miss instead of $b$. Also, the parallel
searches do the same with the batches of nodes they receive ([Chapter
9](09-hash-distributed.md#ch-9)), whose lookups are independent of each
other.

![Figure 3.7. The lookups of the four successors of a node in time:
one after the other, each waits for its miss (above); prefetched, the
four misses overlap after the four hashes, and the comparisons follow
(below).](../figures/data-structures-prefetch.svg){#fig-3-7}

## 3.9 Perfect hashing {#sec-3-9}

A *perfect hash* of the states is an injective map of the states to the
integers, and a *minimal* one maps $n$ states onto $[0, n)$. A domain
with a perfect hash models the concept `PerfectHashDomain`
(`Domain.h:234`), whose `perfect_hash(s)` the distribution HDA$^*$[P] of
Burns et al. [[3]](Z-bibliography.md#ref-3) uses to assign the states to
the workers ([Chapter 9](09-hash-distributed.md#ch-9)). On a grid of
width $W$, the cell $(x, y)$ has the index $x + W y$, which is a minimal
perfect hash of the cells, and the state of the domain `Grid` is that
index itself ([Figure 3.8](03-data-structures.md#fig-3-8)); the explicit
graphs use the index of the node, and the traveling salesperson, the
vacuum world and the missionaries and cannibals pack their states into
words that are their perfect hashes. On the permutations, the
lexicographic ranking of [§5.4](05-heuristics.md#sec-5-4) is a minimal
perfect hash of the placements of a pattern, which the pattern databases
use as the indices of their tables (the ranking of Myrvold and Ruskey
[[34]](Z-bibliography.md#ref-34) would serve as well, cf.
[§5.4](05-heuristics.md#sec-5-4)).

![Figure 3.8. The perfect hash of a grid: the cell $(x, y)$ of a grid
of width $W = 6$ has the index $x + W y$, a bijection of the cells
onto $[0, W H)$.](../figures/data-structures-perfect.svg){#fig-3-8}

With a minimal perfect hash and a known number of states, CLOSED could
be a plain array indexed by it, with no probe at all; the library keeps
the hash table for all the domains, which does not depend on the range
of the hash nor on the number of states, and leaves the perfect hashes
to the distributions and to the tables of the heuristics.

## 3.10 Transposition tables {#sec-3-10}

A *transposition table* is a CLOSED that may forget: it has a fixed
number of places, and when it is full a new state replaces an old one by
a replacement rule (e.g., the one with the smaller $g$, or the one
reached less often), or the new state is not stored. Such a table gives
the depth-first searches, which keep no CLOSED, part of the duplicate
detection of a graph search within a bounded memory: a state found in
the table with a $g$ not larger is not searched again. With the
*stochastic node caching*, a state is stored with a probability $p$ at
each visit, so that a state visited $n$ times is in the table with
probability $1 - (1 - p)^n$, and the states reached by many paths, which
are those whose detection saves the most, are stored with a high
probability.

> **Status — planned.** The library has no transposition table; its
> depth-first search, IDA$^*$, detects only the move back to the parent
> ([§2.9](02-best-first.md#sec-2-9) and
> [Chapter 4](04-memory-bounded.md#ch-4)).

## 3.11 Immediate and delayed duplicate detection {#sec-3-11}

The searches of [Chapter 2](02-best-first.md#ch-2) detect the duplicates
*immediately*: each successor is looked up in CLOSED when it is
generated, and it is dropped at once if its state is known with a $g$
that is not larger. A search with *delayed* duplicate detection instead
appends the successors to a buffer, or to a file on the disk, without
looking them up, and removes the copies later, in a batch, by sorting or
hashing the buffer and comparing it with the states already expanded
([Figure 3.9](03-data-structures.md#fig-3-9)). The batch replaces the
random lookups by sequential passes, as the disk requires, and it lets
CLOSED live outside the memory. Its price is that the copies are stored
until the batch, and that the order of the expansions is that of the
batches, not of the single nodes. The search of [Chapter
4](04-memory-bounded.md#ch-4) (`DDDAStar`) keeps its buckets in files
and detects the duplicates of a file at once, before expanding it, so that
each state is expanded at most once; the test `ddd` checks it, and checks
that the records of a search that exhausts its buckets split, with the
root, into the states stored and the copies dropped.

![Figure 3.9. Immediate and delayed duplicate detection. (a) Each
successor is looked up in CLOSED when it is generated. (b) The
successors are appended to a buffer or a file, copies included, and the
copies are removed later in a batch, together with the states already
closed.](../figures/data-structures-ddd.svg){#fig-3-9}

## 3.12 The sets of data structures {#sec-3-12}

The arena, CLOSED and the OPEN lists of A$^*$ come from one of three
sets, chosen when the library is built (`Structures.h`): the default
one; that of the standard library, chosen by the macro
`PBEFS_STRUCTURES_STD`; and that of Boost, chosen by
`PBEFS_STRUCTURES_BOOST`. CMake defines the macro from its option
`PBEFS_STRUCTURES` (its default value, `std` or `boost`) and the makefile
from
`PBeFSSTRUCT`, and defining both macros is an error. `Node.h`,
`NodeTable.h`, `OpenList.h` and `ConcurrentClosed.h` then include the
classes of the set from `structures/std/`, `structures/boost/` or the
default directory, with the same interface (the adjacency of the explicit
graphs has a version in `structures/boost/` and in the default directory,
that of the grids only in the default one, and the other sets use the
adjacency written in `ExplicitGraph.h` and `Grid.h`); the sets are the
following.

| set | `NodePool` | `NodeTable` | `HeapOpenList` | `BucketOpenList` |
|---|---|---|---|---|
| default | arena, huge pages | linear probing | 4-ary heap | two levels |
| `std` | `std::deque` | `std::unordered_map` | binary heap | binary heap |
| `boost` | `object_pool` | `unordered_flat_map` | 4-ary heap | 4-ary heap |

The searches are the same with all the sets, and so is the order of
their expansions, but for the ties of $f$ and $g$ in OPEN, which the
sets break differently ([§3.5](03-data-structures.md#sec-3-5)); the
costs returned are the same. A search is parametrized by its OPEN list
at compile time, `AStar<D, Open>` taking `BucketOpenList` by default for
integer costs and `HeapOpenList` otherwise, and the function
`astar(domain, open, options)` chooses at run time by `OpenKind`
(`automatic`, `heap` or `buckets`); the arena and CLOSED have no
parameter, being those of the set. Furthermore, the parallel searches
use the same classes in each worker, and their shared structures (e.g.,
the CLOSED of SPA$^*$) are described in [Chapter
8](08-shared-memory.md#ch-8). We keep the standard library and Boost to
measure the searches against them, with the script of the next section.

## 3.13 Benchmark {#sec-3-13}

> **Status — planned.** The table of this section is not yet produced
> by the tables of the manual.

The script `tools/bench/structures.sh` runs `tools/AStarBench` from
three builds of the library, one per set of data structures, on
instances of growing size: the 15-puzzle in four bands of the instances
of Korf [[7]](Z-bibliography.md#ref-7), by the nodes that A$^*$ expands
(about $10^5$, $10^6$, $10^7$ and above $7 \cdot 10^7$), and random
grids of 2000 to 20000 cells a side, with 35% of the cells blocked and 4
and 8 moves. Each instance is solved by the sequential A$^*$ (with both
OPEN lists on the 15-puzzle) and by the parallel searches (by default
HDA$^*$, whose structures are those of each worker, and SPA$^*$, whose
CLOSED is shared) on 1 to 16 threads, twice with each set, the sets
taking turns so that a drift of the machine falls on all of them; each
run is a process of its own, pinned to the cores and the memory of one
node, and the two runs of a set are its A/A. As in
[§2.11](02-best-first.md#sec-2-11), the clock is that of `solve()`,
which does not include the reading of the instance nor the construction
of the domain. Since the expansions are the same for the three sets but
for the last layer, the comparison is one of time and of memory per
node.

## 3.14 Notes and references {#sec-3-14}

The bucket queue for integer keys is due to Dial
[[65]](Z-bibliography.md#ref-65), and its nested use with two levels, the
ties broken in favor of the larger $g$ and the pruning of the move back to
the parent in a fast A$^*$ are studied by Burns et al.
[[4]](Z-bibliography.md#ref-4), who report the bucket queue to be faster
than a binary heap on the sliding-tile puzzles (they update the position
of a node in OPEN when its $g$ improves, where the library uses the lazy
deletion); the textbook of Edelkamp
and Schrödl [[10]](Z-bibliography.md#ref-10) treats the priority queues,
the hash tables and the duplicate detection of the searches at length.
The $d$-ary heaps go back to Johnson [[158]](Z-bibliography.md#ref-158),
the linear probing and the multiplicative hashing are analyzed in
section 6.4 of [[157]](Z-bibliography.md#ref-157), and the arenas, i.e.,
the allocation by lifetime, are those of Hanson
[[159]](Z-bibliography.md#ref-159). The delayed duplicate detection on the
disk is that of Korf [[40]](Z-bibliography.md#ref-40) and of the
External A$^*$ of Edelkamp, Jabbar and Schrödl
[[39]](Z-bibliography.md#ref-39). With respect to
[[4]](Z-bibliography.md#ref-4), the library adds the prefetch of the
lookups, the huge pages, and the sparse lists of the large values of $g$
in the bucket queue, which keep the queue usable on the grids and in the
workers of HDA$^*$; it does not reuse the nodes of the closed states,
which the arena keeps until the end, and it has neither the radix heap
nor the transposition tables.

## 3.15 Exercises {#sec-3-15}

**Exercise 3.1** (`ctest -R PBeFS.open`). Show that the order of
`BucketOpenList` (least $f$, then largest $g$, then the last pushed)
and that of `HeapOpenList` agree on the pairs $(f, g)$ that they pop,
and explain why the test compares the pops with an ordered multiset of
the pairs $(f, -g)$ rather than of the entries.

**Exercise 3.2** (`ctest -R PBeFS.open`). In a worker of HDA$^*$ the
least $f$ of OPEN may fall below the $f$ of the last pop. Show that
[Proposition 3.3](03-data-structures.md#prop-3-3) then fails, give a
sequence of $m$ pushes and pops on which the scan of the layers costs
$\Theta(m f_{\max})$, and explain what the sparse lists of the values
$g \ge 256$ save in the same situation.

**Exercise 3.3** (`ctest -R PBeFS.astar`). On the 8-puzzle, break the
ties of $f$ by the least $h$ and then by the largest $h$, and compare the
nodes expanded; show by
[Proposition 3.5](03-data-structures.md#prop-3-5) that the two rules
differ only in the last layer, and that the costs are equal, as the
test checks for the two OPEN lists of the library.

**Exercise 3.4** (`ctest -R PBeFS.mqa-thread`). Prove [Lemma
3.4](03-data-structures.md#lem-3-4) also when $g + h$ reaches $2^{32} -
1$, and show that it fails when $g + h \ge 2^{32}$. Which representation
of the key should a domain whose costs exceed that range use? (MultiQueue
A$^*$, whose test runs here, orders its queues by these keys.)

**Exercise 3.5** (`ctest -R PBeFS.hda-thread`). Prove that $x + W y$ is
a minimal perfect hash of the cells of a grid of width $W$ and height
$H$, and that the distribution HDA$^*$[P], which takes it modulo the
number $p$ of the workers, gives each worker the cells of whole columns
when $p$ divides $W$. The test runs HDA$^*$[P] on the 15-puzzle, with the
other distributions.

**Exercise 3.6** (`ctest -R PBeFS.ddd`). Compare, on the 8-puzzle, the
nodes expanded and generated by A$^*$ with immediate duplicate
detection and by the search with delayed duplicate detection of
[Chapter 4](04-memory-bounded.md#ch-4); explain why the second expands
each state at most once but may expand more states, and why, on a search
that ends without a goal (as on the grid with the wall of the test), its
records satisfy `generated` $+ 1 =$ `stored` $+$ `duplicates`.

**Exercise 3.7** (`ctest -R PBeFS.astar`). With a load factor of at
most $7/10$, compute the memory of CLOSED per node stored, and that of
the nodes in the arena, for the 15-puzzle (a state of 8 bytes and
integer costs); how much memory does a search that stores $10^8$ nodes
need, without OPEN? Show that `slots_for(n)` places hold $n$ nodes
without a doubling.
