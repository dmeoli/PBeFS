# 16. Classical planning {#ch-16}

[Source: the planning domain of [`include/PBeFS/astar/`][astar], on the branch
`feature/planning-domain`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with classical planning as a domain of the
library: a task in STRIPS form, whose states are sets of facts and whose
operators have preconditions, additions and deletions, the
domain-independent heuristics computed on it, and the distributions of
HDA$^*$ that use the structure of the task (e.g., the abstraction of a
subset of its variables) instead of that of a particular puzzle.

The figures planned for the chapter are the following:

* Figure 16.1: a small STRIPS task (e.g., a logistics problem) with its
  facts, operators and a plan;
* Figure 16.2: the relaxed planning graph of the delete relaxation, with
  the values of the heuristic;
* Figure 16.3: the abstraction of a subset of the variables used by the
  distribution, with the owners it gives.

<!-- the links to the sources, written by srclinks.py -->

[astar]: https://github.com/dmeoli/PBeFS/tree/main/include/PBeFS/astar
