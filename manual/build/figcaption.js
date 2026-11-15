// Show the caption of each figure of the manual under it: a figure is an
// image with an id fig-N-M whose alt text is its caption ("Figure N.M.
// ..."), which Pandoc prints under the image in the PDF; here the image
// is wrapped in a <figure> with that text as its <figcaption>, whose
// mathematics KaTeX renders.
document$.subscribe(() => {
  document.querySelectorAll('img[id^="fig-"]').forEach((img) => {
    if (img.closest('figure')) return;
    const figure = document.createElement('figure');
    const caption = document.createElement('figcaption');
    caption.textContent = img.getAttribute('alt') || '';
    img.replaceWith(figure);
    figure.appendChild(img);
    figure.appendChild(caption);
    // the caption is added after katex-init.js has rendered the page
    renderMathInElement(caption, {
      delimiters: [
        { left: "$", right: "$", display: false },
        { left: "\\(", right: "\\)", display: false }
      ],
      throwOnError: false
    });
  });
});
