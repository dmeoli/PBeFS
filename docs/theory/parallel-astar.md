# Parallel A*

The parallel searches of the library, which run A* with $p$ workers on a
shared-memory machine (on the threads of C++, of OpenMP or of FastFlow) or with
one MPI process per worker, are treated in the
[manual](https://dmeoli.github.io/PBeFS/manual/index.html). [Chapter
8](https://dmeoli.github.io/PBeFS/manual/chapters/08-shared-memory.html) gives
the overheads by which a parallel search is measured, the condition under which
it may stop with an optimal path, the searches on a shared OPEN (SPA*, MQA*,
the farm and PA*SE) and the searches for slow evaluations of the edges (ePA*SE,
GePA*SE, A-ePA*SE, PwA* and MPLP). Then the hash-distributed search HDA*, its work
distributions and its termination detection on MPI are in [chapter
9](https://dmeoli.github.io/PBeFS/manual/chapters/09-hash-distributed.html),
and the distributions built from a sample of the instance and their automatic
choice in [chapter
10](https://dmeoli.github.io/PBeFS/manual/chapters/10-auto-distributions.html).
Safe PBNF and plain PBNF are in [chapter
11](https://dmeoli.github.io/PBeFS/manual/chapters/11-safepbnf.html), and the
further schemes (PSDD and its variants, KBFS, the portfolio of A*, PRA* and the
lending of nodes, TDS, the parallel IDA* below a frontier and the parallel
greedy searches of [Greedy best-first search](greedy-search.md)) in [chapter
12](https://dmeoli.github.io/PBeFS/manual/chapters/12-more-parallel-schemes.html).
Also, the pseudocode of each search is on the page of its header in the API
reference.
