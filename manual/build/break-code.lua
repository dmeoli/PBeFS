-- break-code.lua
-- For the LaTeX/PDF build, make inline code spans (\texttt) breakable so
-- long file paths, URL patterns and API signatures written as `code` do not
-- overflow the right margin, without cutting a name in two. A zero-penalty
-- \allowbreak follows the characters that separate the parts of a name
-- (/ . _ , ( the second colon of :: and a - between two letters or
-- digits, as in a ctest name), and a high penalty the boundary
-- of two words of a CamelCase name, which TeX takes only when no other
-- break is left; within a word there is no break at all, and a line that
-- cannot be broken otherwise stretches (\emergencystretch in preamble.tex).
-- (Bare links and \url/\href are handled separately by the xurl package.)

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
  local chars = {}
  for _, cp in utf8.codes(el.text) do
    chars[#chars + 1] = utf8.char(cp)
  end
  local out = {}
  for i, c in ipairs(chars) do
    local nxt = chars[i + 1]
    out[#out + 1] = esc_char(c)
    if nxt then
      local prv = chars[i - 1]
      if c:match('[/._,(]') or (c == ':' and prv == ':')
         or (c == '-' and prv and prv:match('%w') and nxt:match('%w')) then
        out[#out + 1] = '\\allowbreak{}'
      elseif c:match('[%l%d]') and nxt:match('%u') then
        out[#out + 1] = '\\penalty5000{}'
      end
    end
  end
  return pandoc.RawInline('latex', '\\texttt{' .. table.concat(out) .. '}')
end
