# 12. Further parallel schemes {#ch-12}

[Source: the further parallel searches of [`include/PBeFS/astar/`][astar], on
the branch `feature/more-parallel-schemes`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with the parallel searches that divide the
work in ways other than those of the previous chapters: Parallel
Structured Duplicate Detection (PSDD), which expands in parallel the
layers of disjoint parts of an abstract graph; K-Best-First Search
(KBFS), which expands the $k$ best nodes of OPEN at a time; Multi-A$^*$,
in which independent searches share the incumbent; and Parallel
Retracting A$^*$ (PRA$^*$), the first of the searches that distribute
the nodes by a hash of the state, with a bounded memory per worker.

The figures planned for the chapter are the following:

* Figure 12.1: the layers of PSDD on the abstract graph, with the
  disjoint duplicate detection scopes expanded in parallel;
* Figure 12.2: the $k$ nodes taken by KBFS from OPEN in one round, before
  and after the round;
* Figure 12.3: the independent searches of Multi-A$^*$ with their shared
  incumbent;
* Figure 12.4: the retraction of PRA$^*$: the memory of a worker before
  and after a retraction.

<!-- the links to the sources, written by srclinks.py -->

[astar]: https://github.com/dmeoli/PBeFS/tree/main/include/PBeFS/astar
