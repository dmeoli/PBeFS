# 6. Bounded-suboptimal and anytime search {#ch-6}

[Source: the weighted and anytime searches of [`include/PBeFS/astar/`][astar],
on the branch `feature/bounded-suboptimal`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with the searches that give up the optimality
of the solution for a smaller number of expansions while still bounding
its cost: weighted A$^*$, whose solution costs at most $w$ times the
optimum when the heuristic is admissible and $w \ge 1$, the parallel
searches that keep this bound (e.g., weighted PBNF and weighted
HDA$^*$), and the anytime searches (e.g., anytime repairing A$^*$),
which find a first solution quickly and then improve it, together with a
lower bound on the optimum, until they prove it optimal or are stopped.

The figures planned for the chapter are the following:

* Figure 6.1: the contours of $f = g + w\,h$ against those of $f = g +
  h$ on a small grid, with the states expanded by each search;
* Figure 6.2: the proof of the bound $w$ drawn on a path: the node of
  OPEN on an optimal path and its $f$-value;
* Figure 6.3: the successive solutions of an anytime search and the
  lower bound closing on them, as a staircase over time;
* Figure 6.4: the INCONS list of anytime repairing A$^*$ before and
  after a decrease of $w$.

<!-- the links to the sources, written by srclinks.py -->

[astar]: https://github.com/dmeoli/PBeFS/tree/main/include/PBeFS/astar
