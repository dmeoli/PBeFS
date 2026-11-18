# Branch-and-bound as best-first search

The branch-and-bound (B&B) of the library, i.e., the nodes, the bounds, the
incumbent, the pruning and the selection rules, the invariant that makes the
search correct, the bound that it reports when it stops early and A* as a
special case of the same scheme, is treated in
[chapter 13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of the [manual](https://dmeoli.github.io/PBeFS/manual/index.html),
with its proofs; the [design note](../branch-and-bound.md) points to the
sections of the chapter that treat each part of the concept `BranchAndBound`
of [BranchAndBound.h](../../include/PBeFS/bnb/BranchAndBound.h).
