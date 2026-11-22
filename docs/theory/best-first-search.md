# Best-first search and A*

The sequential searches of the library are treated in the
[manual](https://dmeoli.github.io/PBeFS/manual/index.html), which gives for
each of them the idea, the invariants, the pseudocode, the proofs and the code.
In it, the notation shared by all the searches and the conventions of the
pseudocode of the headers are in [chapter
1](https://dmeoli.github.io/PBeFS/manual/chapters/01-introduction.html), and
best-first search and A* in [chapter
2](https://dmeoli.github.io/PBeFS/manual/chapters/02-best-first.html). This
chapter gives its invariant, admissibility, consistency, the reopening of the
closed nodes, the pruning of the parent and the goal test at the generation of
a state, and its sections 2.12 to 2.15 treat the satisficing searches of the
planners (cf. [Greedy best-first search](greedy-search.md)). Next, the data
structures of a search are in [chapter
3](https://dmeoli.github.io/PBeFS/manual/chapters/03-data-structures.html).
Finally, [chapter
4](https://dmeoli.github.io/PBeFS/manual/chapters/04-memory-bounded.html)
treats the searches that hold less than CLOSED in memory or bound its growth,
i.e., beam search, iterative deepening with the bounds of IDA*-CR, recursive
best-first search, delayed duplicate detection on the disk, partial expansion,
fringe search, the real-time searches and the hierarchical path-finding on
grids. Also, the pseudocode of each search is on the page of its header in the
API reference.
