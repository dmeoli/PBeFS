# PBeFS — Parallel Best-First Search

**Donato Meoli**
*Dipartimento di Informatica, Università di Pisa*
`donato.meoli@outlook.com`

*A manual of the algorithms of PBeFS, their proofs and their code*

---

## Abstract

PBeFS is a C++-20 library for the best-first search of an implicit
state-space graph, i.e., A$^*$ and its variants, and of a tree with an
incumbent, i.e., branch-and-bound, sequentially and in parallel on
C++ threads, OpenMP, FastFlow and MPI. A problem is given to the
library as a `Domain` (the initial state, the successors of a state
with the costs of the edges, the goal test and the heuristic) or as a
branch-and-bound problem (the root, the branching, the bound and the
value of a solution), and every search of the library runs on every
problem that models the concepts it needs.

This manual describes each search of the library in the same way: the
idea, the invariants it keeps, figures of the state of the search
before and after a step, the pseudocode (the same as in the Doxygen
comments of the headers), the analysis (optimality, bounded
suboptimality, termination of the parallel searches) as lemmas and
theorems with their proofs, the excerpt of the code that implements
it, taken from the headers at build time, and the engineering that
makes the implementation fast. Each chapter ends with the benchmark of
its searches and with exercises, each of which is checked by a test of
the library. The manual is a companion to the
[API reference](https://dmeoli.github.io/PBeFS/) and to the README of
the [repository](https://github.com/dmeoli/PBeFS), and it is
available both as a web site and as a
[PDF](PBeFS-manual.pdf).

## License

PBeFS is distributed under the
[GNU Lesser General Public License v3.0](https://opensource.org/licenses/lgpl-3.0.html)
(LGPL-3); this manual is distributed under the same terms.

---

## Conventions

* **Class and function names** are typeset in `monospace`, e.g.,
  `AStar`, `HDAStar::solve()`, and so are the files of the
  repository, e.g., [`include/PBeFS/astar/AStar.h`][AStar.h]; a precise pointer
  into the source has the form `file:line`.
* **Each chapter** starts with a line `[Source: ...]` that lists the
  files it describes.
* **Mathematical notation** is fixed in
  [Chapter 1](chapters/01-introduction.md#ch-1); in particular, $s_0$
  is the initial state, $c(u, v)$ the cost of the edge $(u, v)$, $g$,
  $h$ and $f = g + h$ the usual functions of A$^*$, $z$ the incumbent
  and $p$ the number of workers.
* **Pseudocode** has numbered lines and the conventions of the Doxygen
  comments of the headers, which are listed in
  [§1.4](chapters/01-introduction.md#sec-1-4).
* **Lemmas, theorems and exercises** are numbered within their
  chapter; the end of a proof is marked by $\square$. Each exercise
  names the test (`ctest -R PBeFS.<case>`) that checks its solution.
* **Benchmark tables** are written by [`manual/build/tables.py`][tables.py]
  from the results of the benchmark tools of the repository; a table whose
  results have not been measured yet on an idle machine is an empty
  table marked as a placeholder.
* **Status admonitions** flag the maturity of a feature explicitly.
  Four forms appear in the text:
    * **Status — mature.** The feature is stable, tested and in tagged
      releases. Implicit unless stated.
    * **Status — beta.** The feature is in `main` but its interface or
      semantics may change in incompatible ways.
    * **Status — under development.** The feature is partly present
      or is being written on a branch; the manual describes what
      exists today and the direction of the work.
    * **Status — planned.** The feature is documented as intended but
      is not yet written. The manual indicates explicitly that it
      cannot be used today.

When in doubt about a claim made in this manual, the authoritative
sources are, in order of decreasing specificity: the headers of the
repository at <https://github.com/dmeoli/PBeFS>, and the API reference.

---

## Table of contents (high-level)

* **Part I — Sequential search.**
  [Chapter 1](chapters/01-introduction.md#ch-1) introduces the problem,
  the notation and the concepts of the library;
  [Chapter 2](chapters/02-best-first.md#ch-2) treats best-first search
  and A$^*$;
  [Chapter 3](chapters/03-data-structures.md#ch-3) the data structures
  of a search, i.e., the open and closed lists, the bucket and radix
  queues and the hash tables;
  [Chapter 4](chapters/04-memory-bounded.md#ch-4) the searches whose
  memory is bounded or on disk, among them iterative deepening and the
  delayed duplicate detection;
  [Chapter 5](chapters/05-heuristics.md#ch-5) the heuristics and the
  pattern databases;
  [Chapter 6](chapters/06-bounded-suboptimal.md#ch-6) the
  bounded-suboptimal, anytime and incremental searches;
  [Chapter 7](chapters/07-bidirectional.md#ch-7) the bidirectional and
  external-memory searches.
* **Part II — Parallel search.**
  [Chapter 8](chapters/08-shared-memory.md#ch-8) treats the parallel
  searches on a shared memory;
  [Chapter 9](chapters/09-hash-distributed.md#ch-9) the hash-distributed
  search and its distributions;
  [Chapter 10](chapters/10-auto-distributions.md#ch-10) the automatic
  choice of the distribution;
  [Chapter 11](chapters/11-safepbnf.md#ch-11) the parallel
  best-nblock-first search;
  [Chapter 12](chapters/12-more-parallel-schemes.md#ch-12) further
  parallel schemes.
* **Part III — Branch-and-bound.**
  [Chapter 13](chapters/13-branch-and-bound.md#ch-13) treats
  branch-and-bound and its parallel schemes.
* **Part IV — Problems.**
  [Chapter 14](chapters/14-domains.md#ch-14) describes the example
  domains of the library;
  [Chapter 15](chapters/15-msa-instances.md#ch-15) the multiple
  sequence alignment and the sets of instances;
  [Chapter 16](chapters/16-planning.md#ch-16) classical planning.
* **Part V — Multi-objective search.**
  [Chapter 17](chapters/17-multi-objective.md#ch-17) treats the searches
  whose edges cost a vector and whose answer is the Pareto front, sequential
  and parallel.
* **Part VI — Symbolic search.**
  [Chapter 18](chapters/18-symbolic.md#ch-18) treats the searches on sets
  of states represented by binary decision diagrams, with their heuristics,
  abstractions and parallel images.
* **Back matter.** The [Glossary](chapters/Y-glossary.md#glossary)
  and the [Bibliography](chapters/Z-bibliography.md#bib).

<!-- the links to the sources, written by srclinks.py -->

[AStar.h]: https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/astar/AStar.h
[tables.py]: https://github.com/dmeoli/PBeFS/blob/main/manual/build/tables.py
