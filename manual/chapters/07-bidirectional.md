# 7. Bidirectional search {#ch-7}

[Source: the bidirectional and external searches of
`include/PBeFS/astar/`, on the branch `feature/bidirectional-external`]

> **Status — under development.** The code that this chapter describes
> is being written and is not in a release yet; the chapter is filled
> when that code is.

This chapter is concerned with the searches that grow a frontier from
the initial state and one from the goal, and stop when the two meet
with a proof that no cheaper path joins them: the front-to-end
bidirectional searches, the meet-in-the-middle condition by which a
bidirectional search with an admissible heuristic never expands a
node beyond half of the optimal cost on either side, and the stopping
rules that make the first meeting of the frontiers an optimal
solution; and with the external variants of these searches, whose
frontiers are kept on disk as in the delayed duplicate detection of
[Chapter 4](04-memory-bounded.md#ch-4).

The figures planned for the chapter are the following:

* Figure 7.1: the two frontiers on a small grid, before and after they
  meet, with the cost of the best path found and the bound that
  stops the search;
* Figure 7.2: the regions expanded by unidirectional A$^*$ and by a
  meet-in-the-middle search on the same instance, as the sets of
  states with $g$ below half of the optimal cost on each side;
* Figure 7.3: the stopping rule as the two lower bounds of the
  frontiers rising until they reach the cost of the incumbent path.
