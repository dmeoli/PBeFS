-- include.lua
-- Single-source inclusion of what the build generates into ../generated:
-- the excerpts of code that snarf.py takes from the headers and the tables
-- that tables.py writes from the results of the benchmarks. The Markdown
-- names them with the syntax of pymdownx.snippets, which the web build
-- expands by itself:
--
--     ```cpp
--     --8<-- "snippets/best-first-expand.cpp"
--     ```
--
--     --8<-- "tables/best-first-korf.md"
--
-- For the PDF build this filter does the same: a code block whose only
-- line is such a marker gets the text of the file, and a paragraph that
-- is only a marker is replaced by the blocks of the Markdown file. A file
-- that is not there is an error, not a silent omission.

local GENERATED = '../generated/'

local function target(s)
  -- the path of a marker, also after the smart punctuation of Pandoc
  -- (the dashes and the quotes are multi-byte, so they are replaced first)
  s = s:gsub('—', '---'):gsub('–', '--'):gsub('“', '"'):gsub('”', '"')
  return s:match('^%s*%-%-8<%-%-%s*"([^"]+)"%s*$')
end

local function slurp(path)
  local f = io.open(GENERATED .. path, 'r')
  if not f then
    error('include.lua: ' .. GENERATED .. path .. ' is not there; run '
          .. 'snarf.py and tables.py first (make pdf does)')
  end
  local text = f:read('a')
  f:close()
  return text
end

function CodeBlock(el)
  local path = target(el.text)
  if path then
    el.text = slurp(path):gsub('\n+$', '')
    return el
  end
end

function Para(el)
  local path = target(pandoc.utils.stringify(el))
  if path then
    return pandoc.read(slurp(path), 'markdown+pipe_tables').blocks
  end
end
