# 15. Multiple sequence alignment and the instances {#ch-15}

[Source: the MSA domain of `include/PBeFS/astar/` and the instance sets
of `data/`, on the branch `feature/msa-and-instances`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with the Multiple Sequence Alignment (MSA)
problem as a shortest path on a lattice, with its pairwise heuristic,
which is the usual test of a parallel search with a large branching
factor and a costly heuristic, and with the sets of instances on which
the searches of this manual are compared (the 15-puzzle and the
24-puzzle, the grids, the traveling salesman and the alignments), with
their generators and the format of their files.

The figures planned for the chapter are the following:

* Figure 15.1: the lattice of the alignment of two and of three
  sequences, with an optimal path and the gaps it stands for;
* Figure 15.2: the pairwise heuristic as the sum of the projections of a
  state on the faces of the lattice;
* Figure 15.3: the distribution of the optimal costs of each instance
  set.
