# Parallel branch-and-bound

The parallel B&B searches of the library (the shared queue, the distributed
queues, the work stealing, the MultiQueue, the deterministic rounds, the
SelfSplit, the farm, the processes of MPI, the search skeletons, the racing
ramp-up, the master, hubs and workers over MPI and the work stealing of
intervals of leaves), with what each one guarantees, how it ends, whether its
outcome is deterministic, and the anomalies in the number of evaluations of a
parallel B&B, are treated in [chapter 13](https://dmeoli.github.io/PBeFS/manual/chapters/13-branch-and-bound.html) of the [manual](https://dmeoli.github.io/PBeFS/manual/index.html); the
pseudocode of each scheme is also on the page of its header.
