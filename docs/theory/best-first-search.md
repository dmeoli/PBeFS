# Best-first search and A*

The sequential searches of the library are treated in the
[manual](https://dmeoli.github.io/PBeFS/manual/index.html), which gives for each of them the idea, the
invariants, the pseudocode, the proofs and the code: the notation shared by all
the searches and the conventions of the pseudocode of the headers are in
[chapter 1](https://dmeoli.github.io/PBeFS/manual/chapters/01-introduction.html), best-first search and A* (its invariant, admissibility, consistency,
the reopening of the closed nodes, the pruning of the parent and the goal test
at the generation of a state) in [chapter 2](https://dmeoli.github.io/PBeFS/manual/chapters/02-best-first.html), together with the satisficing searches
of the planners (greedy best-first search and its OPEN lists, enforced
hill-climbing, the iterated searches of LAMA and the width-based searches) in
[sections 2.12 to 2.15](https://dmeoli.github.io/PBeFS/manual/chapters/02-best-first.html#sec-2-12), the data structures of a search in [chapter 3](https://dmeoli.github.io/PBeFS/manual/chapters/03-data-structures.html),
and the searches that hold less than CLOSED in memory or bound its growth
(beam search, iterative deepening with the bounds of IDA*-CR, recursive
best-first search, delayed duplicate detection on the disk, partial expansion,
fringe search, the real-time searches and the hierarchical path-finding on
grids) in [chapter 4](https://dmeoli.github.io/PBeFS/manual/chapters/04-memory-bounded.html). The pseudocode of each search is also on the page of its
header in the API reference.
