-- break-code.lua
-- For the LaTeX/PDF build, make inline code spans (\texttt) breakable so
-- long file paths, URL patterns and API signatures written as `code` do not
-- overflow the right margin. We escape LaTeX specials per character and
-- insert a zero-penalty \allowbreak after each one: because the penalty is
-- zero, TeX only ever breaks inside a code span when the line would
-- otherwise overflow, so short identifiers stay intact. (Bare links and
-- \url/\href are handled separately by the xurl package.)

local function esc_char(c)
  if c == '\\' then
    return '\\textbackslash{}'
  elseif c:match('[%%%$#&_{}]') then
    return '\\' .. c
  elseif c == '^' then
    return '\\textasciicircum{}'
  elseif c == '~' then
    return '\\textasciitilde{}'
  else
    return c
  end
end

function Code(el)
  if not FORMAT:match('latex') then
    return nil
  end
  local out = {}
  for _, cp in utf8.codes(el.text) do
    out[#out + 1] = esc_char(utf8.char(cp))
    out[#out + 1] = '\\allowbreak{}'
  end
  return pandoc.RawInline('latex', '\\texttt{' .. table.concat(out) .. '}')
end
