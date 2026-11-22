# Parallel branch-and-bound

The parallel B&B searches of the library are treated in [chapter
13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of
the [manual](https://dmeoli.github.io/PBeFS/manual/index.html), and their
sections are listed in the [design note](../branch-and-bound.md). These are the
shared queue, the distributed queues, the work stealing, the MultiQueue, the
deterministic rounds, the SelfSplit, the farm, the processes of MPI, the search
skeletons, the racing ramp-up, the master, hubs and workers over MPI and the
work stealing of intervals of leaves. For each of them the chapter gives what
it guarantees, how it ends and whether its outcome is deterministic, and it
also treats the anomalies in the number of evaluations of a parallel B&B. Also, the pseudocode of each scheme is on the page of its header.
