// Render math delimited by $...$ and $$...$$ once each page (and each
// instant-navigation swap) has loaded, using the KaTeX auto-render
// extension. Paired with pymdownx.arithmatex in generic mode.
document$.subscribe(() => {
  renderMathInElement(document.body, {
    delimiters: [
      { left: "$$", right: "$$", display: true },
      { left: "$", right: "$", display: false },
      { left: "\\(", right: "\\)", display: false },
      { left: "\\[", right: "\\]", display: true }
    ],
    throwOnError: false
  });
});
