# Branch-and-bound as best-first search

With its proofs, the branch-and-bound (B&B) of the library is treated in
[chapter
13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of the [manual](https://dmeoli.github.io/PBeFS/manual/index.html). There one finds the nodes, the bounds, the incumbent, the pruning and the selection rules, the invariant that
makes the search correct, the bound that it reports when it stops early, and A*
as a special case of the same scheme. In turn, the [design note](../branch-and-bound.md) gives the sections of the chapter that treat each part of the concept
`BranchAndBound` of
[BranchAndBound.h](../../include/PBeFS/bnb/BranchAndBound.h).
