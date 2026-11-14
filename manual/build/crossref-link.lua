-- crossref-link.lua
-- Single-source cross-reference handling. In the Markdown source, an
-- internal reference is written as an ordinary relative link to another
-- chapter file plus an explicit anchor, e.g. `[Chapter 7](07-physical-
-- abstract.md#ch-7)` or `[§7.6](07-physical-abstract.md#sec-7-6)`. MkDocs
-- serves these natively (one HTML page per source file, anchors from the
-- explicit `{#id}` heading attributes). For the LaTeX (PDF) build all the
-- sources are concatenated into a single document, so the file part of the
-- target is meaningless. Doxygen (https://…) and figure links are untouched.
--
-- Why \hyperref and not a plain `#anchor` link: Pandoc renders an internal
-- `[text](#id)` link as `\hyperlink{id}{text}`, whose destination is the
-- `\hypertarget{id}{…}` it wraps around the heading. For a `\chapter`, that
-- hypertarget sits *before* the `\clearpage` the chapter triggers, so the
-- destination lands at the bottom of the previous (often half-empty) page.
-- Pandoc also emits a `\label{id}` *after* the sectioning command; that
-- label's anchor is the heading's real page anchor (the same one the auto
-- table of contents uses). So we emit `\hyperref[id]{text}` instead, which
-- targets the label and lands exactly on the heading's page. (Subsection
-- links are unaffected — no page break — but stay correct under \hyperref.)
--
-- The anchors that are not headings, i.e., the entries of the bibliography
-- and the lemmas, are written `<a id="ref-N"></a>` (as in the README of
-- the repository); for the PDF build they become a \label and a
-- \hypertarget at that point, so that a citation
-- `[[N]](Z-bibliography.md#ref-N)` lands on its entry.

function RawInline(el)
  -- Pandoc reads `<a id="x"></a>` as two raw inlines, the opening tag
  -- and the closing one: the first becomes the label, the second nothing
  if FORMAT:match('latex') and el.format == 'html' then
    local id = el.text:match('^<a id="([%w%-]+)">$')
      or el.text:match('^<a id="([%w%-]+)"></a>$')
    if id then
      -- a \label for \hyperref (the links to another file) and a
      -- \hypertarget for \hyperlink (the links inside the same file)
      return pandoc.RawInline('latex', '\\phantomsection\\label{' .. id
                              .. '}\\hypertarget{' .. id .. '}{}')
    end
    if el.text == '</a>' then
      return {}
    end
  end
  return el
end
