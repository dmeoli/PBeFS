# Conventions of the manual

## Files

- One chapter per file, `chapters/<topic>.tex`, labelled `ch:<topic>`; one
  picture per file, `figs/<topic>-<name>.tex`, holding the `tikzpicture` only
  and input by the chapter inside its `figure`; labels `sec:`, `fig:`, `alg:`,
  `lem:`, `thm:`, `tab:`, `ex:` followed by `<topic>-<name>`.
- Code is never copied into the manual: an excerpt is marked in its header
  by two lines `//! [<name>]`, which `snarf.py` extracts and `\snippet{<name>}`
  prints (the same markers serve `\snippet` in Doxygen).
- Numbers are never typed: each table is written by `tables.py` from the CSV
  of `tools/bench`, measured on an idle machine, and is a marked placeholder
  until then.
- Each exercise names the case of `test/tests.txt` that checks its solution,
  `\begin{exercise}[<case>]`.
- Lines of at most 79 columns; the boxes at the top and at the bottom of each
  file as in the existing ones.

## Prose

- Sentences long by subordination (mean 28 to 33 words, median 25 to 30),
  with the punctuation carrying the structure: a colon opens an explanation, a
  semicolon chains a contrast or a consequence, parentheses hold asides,
  hedges and pointers ("(cf. (2.4))", "(see \cref{...})"); no em dashes, no
  bold in the prose, italics only for a term where it is defined.
- `we` plural and active for the authors, `one` for the reader; `i.e.,` and
  `e.g.,` with the comma; never `etc.`; enumerations inline as (i), (ii),
  (iii), and `itemize` only for parallel technical objects.
- No rhetorical antithesis ("not X, but Y"), no paragraph opened for effect
  or closed on a summary, no words of hype (novel, key, crucial without a
  reason, powerful, leverage, seamless, robust in a vague sense), no tables or
  algorithms that "show", "tell" or "buy"; the hedge goes on the exact point
  of doubt ("usually", "typically", "at least", "possibly").
- Every claim of optimality, consistency, termination or bound is a Lemma or
  a Theorem, announced by the sentence before it, with a proof made of
  pointers; the intuition comes after the proof.
- An algorithm is given as pseudocode with the conventions of the headers
  (numbered lines, OPEN and CLOSED in small capitals, `g(u)`, `h(u)`,
  `f(u)`, `c(u, v)`, the incumbent `z`, `p` workers indexed by `i`), then
  commented line by line where it is not obvious.
- Each fact is said once, where it belongs, and the other places point to it.
- Computational results give the machine, the compiler, the precision, the
  measure and what the clock does not include, explain every column, and
  compare by ranges ("2 to 3 times").
- Before a chapter is committed it is measured on the PDF:

      make chapter CH=<topic> OUT=/tmp/<dir>
      pdftotext /tmp/<dir>/<topic>.pdf out.txt

  with a checker of the shape of the sentences and of the punctuation, and
  every value out of the band above is fixed.

## Figures

- TikZ, with the styles of `preamble.sty`: `open` (white) and `closed`
  (grey) states, `goal` (double border), `current` (thick border) for the
  node being expanded, `edge` for the arcs of the graph, `list` for the cells
  of OPEN, CLOSED and the queues, `worker` for the box of a worker and of
  what it owns, `msg` (dashed) for a message in flight.
- A step of an algorithm is drawn as two pictures side by side, before and
  after, with the same layout, so that only what changes moves; the values
  `g/h` or `f` are written under the state, never on its edges.
- Black and grey only, readable when printed; no label touches a line or
  another label, which is checked on the rendered page.
