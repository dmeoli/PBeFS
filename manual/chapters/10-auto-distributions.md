# 10. Automatic work distribution {#ch-10}

[Source: [`include/PBeFS/astar/Distribution.h`][Distribution.h], on the branch
`feature/auto-distributions`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with the choice of the distribution of HDA$^*$
from the problem itself, without an abstraction made by hand: the
graph-partitioning distributions GAZHDA$^*$ and GRAZHDA$^*$, which build
an abstract graph and cut it so as to balance the load and keep most
edges inside a worker, the feature-based FAZHDA$^*$, the operator-based
OZHDA$^*$ and the dynamic DAHDA$^*$; and the measures by which the
distributions are compared, i.e., the search overhead (SO), the
communication overhead (CO) and the load balance (LB), as defined by
Jinnai and Fukunaga [[5]](Z-bibliography.md#ref-5).

The figures planned for the chapter are the following:

* Figure 10.1: the abstract graph of a small puzzle, its partition by
  GRAZHDA$^*$, and the edges cut by it;
* Figure 10.2: the same state space coloured by owner under Z, AZ and the
  automatic distribution;
* Figure 10.3: SO, CO and LB of each distribution as a point in the plane
  of CO against LB;
* Figure 10.4: the features selected by FAZHDA$^*$ on a state, and the
  operators of OZHDA$^*$.

<!-- the links to the sources, written by srclinks.py -->

[Distribution.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/Distribution.h
