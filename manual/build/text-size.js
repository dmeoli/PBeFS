// Two magnifiers in the header, next to the choice of the theme, that make
// the whole page smaller or larger, as the zoom of the browser does: they
// scale the root of Material (set by the width of the screen in extra.css),
// so the column grows with the text and a line keeps its number of
// characters. The choice is kept in the browser for the next pages.
(() => {
  const STEPS = [0.85, 1, 1.15, 1.3, 1.5];
  const KEY = "manual-zoom";
  let step = 1;
  try {
    const s = parseInt(localStorage.getItem(KEY), 10);
    if (s >= 0 && s < STEPS.length) step = s;
  } catch (e) {}
  const path = (sign) =>
    "M15.5,14L20.5,19L19,20.5L14,15.5V14.71L13.73,14.43C12.59,15.41 " +
    "11.11,16 9.5,16A6.5,6.5 0 0,1 3,9.5A6.5,6.5 0 0,1 9.5,3A6.5,6.5 0 0," +
    "1 16,9.5C16,11.11 15.41,12.59 14.43,13.73L14.71,14H15.5M9.5,14C12,14 " +
    "14,12 14,9.5C14,7 12,5 9.5,5C7,5 5,7 5,9.5C5,12 7,14 9.5,14M7,9H12V10" +
    "H7V9Z" + (sign > 0 ? "M9,7H10V12H9V7Z" : "");
  const root = document.documentElement;
  const apply = () => {
    root.style.fontSize = "";
    if (step !== 1)
      root.style.fontSize =
        parseFloat(getComputedStyle(root).fontSize) * STEPS[step] + "px";
    document.querySelectorAll(".text-size").forEach((b) => {
      const d = +b.dataset.d;
      b.disabled = step + d < 0 || step + d >= STEPS.length;
    });
  };
  const button = (d, title) => {
    const b = document.createElement("button");
    b.className = "md-header__button md-icon text-size";
    b.dataset.d = d;
    b.title = b.ariaLabel = title;
    b.innerHTML = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">'
      + '<path d="' + path(d) + '"/></svg>';
    b.addEventListener("click", () => {
      step = Math.min(STEPS.length - 1, Math.max(0, step + d));
      try { localStorage.setItem(KEY, step); } catch (e) {}
      apply();
    });
    return b;
  };
  apply();
  matchMedia("screen and (min-width: 95em)").addEventListener("change", apply);
  document$.subscribe(() => {
    if (document.querySelector(".text-size")) return apply();
    const at = document.querySelector(".md-header__option") ||
               document.querySelector(".md-header__inner > label[for=__search]");
    if (!at) return;
    at.before(button(-1, "Smaller text"), button(+1, "Larger text"));
    apply();
  });
})();
