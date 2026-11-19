# Best-first search with an incumbent: design note

The extension of the searches of the library to Branch-and-Bound (B&B), i.e.,
to the minimization of a function over a finite (but possibly huge) set that
is recursively partitioned into nodes with a lower bound and possibly a
feasible solution, is treated in [chapter 13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of the [manual](https://dmeoli.github.io/PBeFS/manual/index.html): the
concept `BranchAndBound` of [BranchAndBound.h](../include/PBeFS/bnb/BranchAndBound.h)
and its requirements in [section 13.1](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-1), the nodes, the bounds and the incumbent in
[section 13.2](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-2), the sequential search and its correctness in
[sections 13.3 and 13.4](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-3), the selection rules in [section 13.5](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-5), the
parallel searches and what each one guarantees in
[sections 13.7 to 13.10](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-7), the reoptimization in [section 13.11](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-11), the
interface that a B&B solver built on a modeling system would use in
[section 13.12](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-12), the children built one at a time and the modes that decide or
enumerate in [section 13.13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-13), the search skeletons, the racing ramp-up, the
checkpoints, the hubs over MPI and the intervals of leaves in
[sections 13.14 to 13.18](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-14), and the example problems in [section 13.19](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html#sec-13-19).
The test of the searches is the case `bb` and the cases `bb-*` of
[test/tests.txt](../test/tests.txt).
