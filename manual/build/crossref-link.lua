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
-- The anchors that are not headings, i.e., the entries of the bibliography,
-- are written `<a id="ref-N"></a>` (as in the README of the repository);
-- for the PDF build they become a \label at that point, so that a citation
-- `[[N]](Z-bibliography.md#ref-N)` lands on its entry.

function RawInline(el)
  if FORMAT:match('latex') and el.format == 'html' then
    local id = el.text:match('^<a id="([%w%-]+)"></a>$')
    if id then
      return pandoc.RawInline('latex', '\\phantomsection\\label{' .. id .. '}')
    end
  end
  return el
end

function Link(el)
  if FORMAT:match('latex') or FORMAT:match('beamer') then
    local anchor = el.target:match('^[%w%./%-]+%.md#(.+)$')
    if anchor then
      local out = { pandoc.RawInline('latex', '\\hyperref[' .. anchor .. ']{') }
      for _, inl in ipairs(el.content) do out[#out + 1] = inl end
      out[#out + 1] = pandoc.RawInline('latex', '}')
      return out
    end
  end
  return el
end
