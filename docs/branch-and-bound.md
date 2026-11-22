# Best-first search with an incumbent: design note

The extension of the searches of the library to Branch-and-Bound (B&B), i.e.,
to the minimization of a function over a finite (but possibly huge) set that is
recursively partitioned into nodes with a lower bound and possibly a feasible
solution, is treated in [chapter
13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of
the [manual](https://dmeoli.github.io/PBeFS/manual/index.html). This note gives
the section of each part. Thus, the concept `BranchAndBound` of
[BranchAndBound.h](../include/PBeFS/bnb/BranchAndBound.h) and its requirements
are in [section
13.1](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-1),
the nodes, the bounds and the incumbent in [section
13.2](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-2),
the sequential search and its correctness in [sections 13.3 and
13.4](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-3),
and the selection rules in [section
13.5](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-5).
Next, the parallel searches and what each one guarantees are in [sections 13.7 to
13.10](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-7),
the reoptimization in [section
13.11](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-11),
and the interface that a B&B solver built on a modeling system would use in
[section
13.12](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-12).
Section
[13.13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-16)
gives the children built one at a time and the modes that decide or enumerate,
[sections 13.14 to
13.18](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-17)
the search skeletons, the racing ramp-up, the checkpoints, the hubs over MPI
and the intervals of leaves, and [section
13.19](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-22)
the example problems. Finally, the test of the searches is the case `bb` and the cases
`bb-*` of [test/tests.txt](../test/tests.txt).
