-- figure-ext.lua
-- Single-source figure handling: the Markdown always references the web
-- form of a figure (`.svg`), which MkDocs serves natively. For the LaTeX
-- (PDF) build we rewrite the extension to `.pdf`, the vector form that
-- pdflatex/graphicx can include directly, so no rasterisation is needed.
-- Both forms live side by side under manual/figures/.

function Image(img)
  if FORMAT:match('latex') or FORMAT:match('beamer') then
    img.src = img.src:gsub('%.svg$', '.pdf')
  end
  return img
end
