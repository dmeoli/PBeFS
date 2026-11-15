# PBeFS Manual — Maintenance & Handoff Notes

*Internal document for whoever maintains this manual. It is excluded
from the MkDocs site and from the PDF build.*

Read this first. It captures what the manual is, the rules it is
written by, how it is built and deployed, and the traps already hit so
they are not re-hit.

---

## 1. What this is

A manual of **PBeFS**, the C++-20 library of sequential and parallel
best-first search (A$^*$ and its variants, branch-and-bound) on C++
threads, OpenMP, FastFlow and MPI. Audience: readers who want to
understand each search of the library (its idea, its invariants, its
proofs, its code) and authors of new domains and searches. It is a
companion to the API reference (Doxygen and Sphinx, published by
`.github/workflows/docs.yml`) and to the README, not a replacement.

Architecture: Part I (sequential search), Part II (parallel search),
Part III (branch-and-bound), Part IV (problems), then a glossary and
the bibliography. Each chapter on an algorithm follows the same
template, in the manner of a textbook of data structures:

1. the problem and the concept of the library it uses;
2. for each algorithm: the idea, the invariants, figures of the state
   before and after a step, the pseudocode, the analysis (lemmas and
   theorems with proofs), the excerpt of the code, and the engineering
   of the implementation;
3. the benchmark (a table written by `build/tables.py`);
4. notes and references;
5. exercises, each checked by a test of the library.

Single source = Markdown under `manual/chapters/`, built to **both** a
PDF (Pandoc + pdflatex) and a web site (MkDocs Material). Every chapter
describes code that is in the tree; those whose code has not been through
a release yet are marked **Status — beta**.

---

## 2. Authoring conventions — DO NOT drift from these

- **Language:** scientific English, *understatement*. No marketing
  language, no first person singular, no colloquialisms; sentences
  long by subordination, the punctuation carrying the structure, no
  em dashes in the prose, no bold inside the prose except the labels
  below. Comparisons with other codes only in the benchmark sections
  and in the "Notes and references" section, with their source.
- **Only PBeFS.** No other library of the authors is named (the
  default data structures are "the default data structures", the
  mailbox of HDA$^*$ is "the mailbox"), and no design choice is
  presented as taken from someone else's code: the papers are cited.
- **Section numbering:** headings carry *manual* numbers (`# 2.`,
  `## 2.3`); LaTeX's numbering is suppressed in `build/preamble.tex`.
  Sub-sections (`###`) are not numbered and carry no id. Keep the
  sequence gap-free.
- **Status admonitions** (defined in `index.md`): **Status — mature**
  (implicit), **Status — beta**, **Status — under development**,
  **Status — planned**, as block quotes, with these exact spellings.
  An algorithm that the library does not implement is described only
  under **Status — planned**.
- **Source line:** each chapter starts with `[Source: ...]`, the files
  it describes.
- **Math delimiters — always `$…$` / `$$…$$`, never `\(…\)`**, and
  only what both KaTeX and pdflatex render.
- **Cross-references are explicit-id Markdown links**: `{#ch-N}`,
  `{#sec-N-M}`, figures `{#fig-N-M}`, lemmas and theorems by an anchor
  `<a id="lem-N-M"></a>` at the start of their paragraph, citations
  `[[N]](Z-bibliography.md#ref-N)` with the numbers of the README
  (the bibliography has the same list, the same anchors). From a
  chapter the target is the bare sibling file name; from `index.md` it
  is prefixed with `chapters/`. `build/crossref-link.lua` turns them
  into `\hyperref` for the PDF, and the `<a id>` anchors into labels.
- **Lemmas, theorems, definitions, algorithms, exercises:**
  bold-labelled paragraphs numbered within the chapter
  (`**Lemma 2.3** (Name).`), proofs as `*Proof.* ... $\square$`;
  pseudocode in ```text blocks with numbered lines and the conventions
  of the Doxygen comments of the headers (listed in §1.4), using only
  the Unicode symbols mapped in `build/preamble.tex`; each exercise
  names the test that checks it, `` (`ctest -R PBeFS.<case>`) ``.
- **Pseudocode and signatures, never implementation code.** An
  algorithm is given as numbered pseudocode with the notation of the
  theory; the only C++ of a chapter is (i) the API signatures, i.e., the
  declarations of the public classes, functions and option structs that
  a user calls or fills, and (ii) one short usage example. Both are
  taken from the sources at build time, never copied: a signature by two
  lines `//! [<topic>-sig-<what>]` around a region of a header that holds
  only declarations (or quoted by hand when the declaration is mixed with
  its body), a usage example by the markers `//! [<topic>-usage]` in the
  case `manual` of `test/manual.cpp`, so that it is compiled and run by
  `ctest -R PBeFS.manual`; the chapter includes them by a block
  ```` ```cpp ```` whose only line is `--8<-- "snippets/<name>.cpp"`.
  The Doxygen filter `docs/doxygen/markdown.py` blanks the marker lines
  of the headers, so that they are not taken as documentation.
- **Numbers are never typed** in a benchmark section: the table is the
  line `--8<-- "tables/<name>.md"`, written by `build/tables.py` from
  the CSV of `tools/bench/campaign.sh` and `pdb.sh` in
  `manual/results/` (with `machine.txt`, one line describing the
  machine). Until the runs are made on an idle machine the table is an
  empty placeholder that says so. Machine-independent counts quoted
  from the README, the tests or the papers carry their source.

---

## 3. Source of truth & verification discipline

- The headers of the repository are the authority, then the tests, then
  the README; a claim on the code cites `file:line`, a claim on a paper
  its section, page or table. What was only read in an abstract is
  marked as such.
- Every excerpt comes from the headers at build time, so it cannot
  drift; every exercise is tied to a case of `test/tests.txt`.
- After any sizeable edit, run a consistency pass over all chapters
  (cross-references, numbering, status admonitions, tone) and the
  style measure of §5.

---

## 4. Layout

```
manual/
  MAINTENANCE.md        # this file
  index.md              # front matter (title, abstract, conventions, ToC)
  PBeFS-manual.pdf      # built PDF deliverable (committed)
  chapters/             # 01-introduction.md … 16-planning.md,
                        # Y-glossary.md, Z-bibliography.md
  figures/              # NAME.tex (TikZ) + NAME.pdf (LaTeX) + NAME.svg (web)
    tikz-pbefs.tex      # the shared styles and palette of every figure
  build/
    Makefile  pandoc.yaml  preamble.tex  mkdocs.yml  katex-init.js
    include.lua  figure-ext.lua  crossref-link.lua  break-code.lua
    figcaption.js  lamalunga.sty
    snarf.py  tables.py  bibcheck.py  bibcheck-verified.tsv  coverage.py
    refcheck.py  srclinks.py  hooks.py
  generated/            # written by snarf.py and tables.py (not committed)
  results/              # the CSV of the benchmarks, when measured
```

---

## 5. Building locally

From `manual/build/`:

```sh
make pdf        # -> ../PBeFS-manual.pdf   (Pandoc + pdflatex)
make html       # -> ../../manual-site/    (MkDocs Material, --strict)
make serve      # live preview at 127.0.0.1:8000
make figures    # ../figures/NAME.pdf and .svg from NAME.tex
make generated  # the excerpts (snarf.py) and the tables (tables.py)
```

`pandoc` and `mkdocs` may be given as `make pdf PANDOC=...` and
`make html MKDOCS=...`. The PDF is committed after each change of the
text; `.github/workflows/manual.yml` rebuilds both outputs at each push
that touches `manual/` or the headers and keeps them as artifacts, and
`release.yml` attaches the PDF to each release.

The style of each chapter is measured on its pages of the PDF
(`pdftotext -f A -l B`), and every value out of band is fixed.

Pieces of the pipeline and why they exist:

- `pandoc.yaml` — the input files **in reading order**; keep it and the
  MkDocs `nav` in sync when adding a chapter. `number-sections: false`.
- `preamble.tex` — geometry, `\setcounter{secnumdepth}{-1}`, captions
  without LaTeX's label (the caption carries the manual number), the
  framed code blocks, and the `\DeclareUnicodeCharacter` mappings of
  the symbols of the pseudocode, so that pdflatex typesets them.
- `include.lua` — for the PDF, expands the `--8<--` lines that MkDocs
  expands by `pymdownx.snippets`; a missing file is an error.
- `bibcheck.py` — fails if a work cited by the README, the headers
  (by author and year, or by doi), the pages of `docs` or the chapters
  is missing from `chapters/Z-bibliography.md`, or if a DOI or arXiv
  identifier of an entry is not in `bibcheck-verified.tsv` with the
  title its registry gives; `--verify` asks Crossref and arXiv for the
  new ones. The README and the bibliography hold the same numbered
  list: a new work gets the next number in both.
- `refcheck.py` — fails if a link whose text names a numbered object
  (Lemma 3.6, Figure 9.2, ...) points to a missing anchor or to an object
  of another kind or number, which catches a renumbering that left a
  reference behind (a reference to the right number of the wrong
  statement is not caught: renumber with care).
- `coverage.py` — fails if a header, or a class, struct, concept or
  function at namespace scope of `include/PBeFS` (outside `detail` and
  `structures/`), is named in no chapter and is not listed in
  `coverage-omitted.txt` with the reason; a new public name gets its
  place in the manual in the same change.
- `srclinks.py` — makes each file of the tree named in a chapter
  (`Domain.h`, `test/test.cpp`, `AStar.h:104`, `SafePBNF.h:270-271`, a
  directory) a link to it on GitHub, a reference link whose definition
  goes at the end of the chapter; run it after writing a new name.
  `--check` fails on a name that is not a link and on a link to a path
  missing in the tree or to a line past the end of its file. It cannot
  see a line number that a change of the code moved: when the code
  changes, the lines cited by the chapters are checked by hand.
- `hooks.py` — the hook of MkDocs that sets each `$$` display apart
  from its paragraph, so that the Markdown of the site does not read
  `\{`, `*` and `_` in it (Pandoc needs no such care). A line of a
  paragraph must not begin with `>`, `+ `, `- ` or `1. `, which the site
  takes as a quote or a list: inside a formula, wrap the line elsewhere.
- `crossref-link.lua`, `break-code.lua`, `figure-ext.lua` — inter-file
  links to `\hyperref`, breakable inline code, `.svg` to `.pdf` for
  LaTeX.

---

## 6. Adding or updating a figure

1. Write `figures/NAME.tex`, a `standalone` document whose preamble is
   `\usepackage{lmodern}`, `\usepackage{amsmath,amssymb}`,
   `\usepackage{tikz}`, `\usepackage{lamalunga}` and
   `\input{tikz-pbefs}`. The grammar is flat: filled shapes, one colour
   per role, no outlines, gradients or shadows, the ink chosen by
   contrast on the fill (`\readableon`). The colours are computed by
   lamalunga from one primary, `\pbefsprimary` (violet 7030A0): the
   primary is what the step reads or follows (the node being expanded,
   the entries touched: `current`, `touched`, `edge-touched`, `path`),
   the triad colour of the harmony is what the step adds (generated
   nodes, new entries, messages: `added`, `reopened`, `msg`), and the
   golden tints of the primary are the rest (`open`, `closed`,
   `unseen`, `worker`); the four owners of the parallel searches
   (`owner=k`, `ownerbox=k`) are tints of the harmony. Edges are ink at
   45% and 1.1pt, pointers `ptr`, variables `var` (italic, primary).
   `pbefs-tikz.tex` is a shim that lets the first figures, written with
   an earlier palette, build; a figure that still inputs it is to be
   converted.
2. `make figures` builds `NAME.pdf` and `NAME.svg`, with `TEXINPUTS`
   pointing at lamalunga (the copy of its author in `~/GDrive/lamalunga`
   if there is one, else `build/lamalunga.sty`); look at the rendered
   figure: no label may touch a shape or another label, and no text
   smaller than `\scriptsize`.
3. Reference it as `![Figure N.M. Caption.](../figures/NAME.svg){#fig-N-M}`
   in a paragraph of its own (no span or bold before the image, which
   would cost the caption in the PDF); `build/figcaption.js` shows the
   alt text as the caption on the site, Pandoc in the PDF.
4. Commit `NAME.tex`, `NAME.pdf` and `NAME.svg`.

---

## 7. Deployment

The docs workflow builds the MkDocs site and copies it to `manual/` of
the site of the documentation, next to the Sphinx and Doxygen sites;
the home page of Sphinx links it. The PDF is part of the MkDocs site
(it sits in `docs_dir`) and is attached to each release.

Gotchas: `site_dir` must be outside `docs_dir`; `build/katex-init.js`
must ship (the `exclude_docs` list excludes `build/` selectively); the
figures must be committed in both forms.

---

## 8. Updating for a new version — checklist

1. When a branch with the code of a stub chapter is merged, write the
   chapter (the planned figures are listed in it), remove its
   **Status — under development**, and update §1 here.
2. Re-run `snarf.py`: an excerpt whose markers moved or disappeared is
   an error.
3. When the benchmarks are measured on an idle machine, put the CSV in
   `manual/results/` with `machine.txt`, rebuild, and write the
   comments of each benchmark section from the tables.
4. Rebuild the PDF and commit it.

*End of maintenance notes.*
