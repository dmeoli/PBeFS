// Let the code within the prose break as the PDF breaks it (break-code.lua):
// a <wbr> after the characters that separate the parts of a name, i.e.,
// / . _ , ( and the second colon of ::, and nowhere within a word, so that
// SMS++/include/Block.h or Block::get_Variable() can go to the next line
// at a separator without cutting a name in two. It also marks the
// [Source: ...] line under the title of a chapter (see extra.css).
document$.subscribe(() => {
  // the [Source: ...] line under the title of a chapter, set ragged right
  // by extra.css
  const src = document.querySelector(".md-typeset h1 + p");
  if (src && src.textContent.trim().startsWith("[Source:"))
    src.classList.add("source-line");
  const SEP = /([\/._,(]|::)(?=.)/g;
  document.querySelectorAll(".md-typeset :not(pre) > code").forEach((code) => {
    if (code.dataset.wbr) return;
    code.dataset.wbr = "1";
    const walker = document.createTreeWalker(code, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach((node) => {
      const parts = node.nodeValue.split(SEP);
      if (parts.length < 2) return;
      const frag = document.createDocumentFragment();
      parts.forEach((part, i) => {
        if (!part) return;
        frag.append(part);
        if (i % 2) frag.append(document.createElement("wbr"));
      });
      node.replaceWith(frag);
    });
  });
});
