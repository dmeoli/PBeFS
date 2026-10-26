/*
 * The outline of a page (the panel on the right) names the members alone,
 * so that the overloads of a method are the same line more times: the
 * entries of a name shared by more members get the types of the parameters
 * of the member they point to, e.g., solve() and
 * solve(std::vector< BBLeaf< Node > >) of BranchAndBoundSearch. The
 * workflow of the documentation appends this file to the navtree.js of
 * Doxygen [see .github/workflows/docs.yml].
 */

(function () {
  // a type without the const in front and the reference at the end
  function bare(t) {
    return t.replace(/\s+/g, ' ').replace(/^const /, '')
            .replace(/\s*&+\s*$/, '').trim();
  }

  // the parameters of a declaration "name (T a, U b = x)", split at the
  // commas outside the angle brackets and the parentheses
  function declaredTypes(text) {
    const open = text.indexOf('(');
    const close = text.lastIndexOf(')');
    if (open < 0 || close < open) return null;
    const args = [];
    let depth = 0, cur = '';
    for (const c of text.substring(open + 1, close)) {
      if ('<(['.includes(c)) ++depth;
      if ('>)]'.includes(c)) --depth;
      if (c === ',' && depth === 0) { args.push(cur); cur = ''; }
      else cur += c;
    }
    args.push(cur);
    return args.map(a => a.replace(/=.*$/, '').trim())
      .filter(a => a.length && a !== 'void')
      .map(a => bare(a.replace(/^(.*[\s*&>])\s*[A-Za-z_]\w*$/, '$1')));
  }

  // the types of the parameters of the member at the anchor id, from the
  // cells paramtype of its detailed documentation or, if it has none, from
  // its line among the declarations
  function parameterTypes(id) {
    let types = null;
    const anchor = document.getElementById(id.replace(/^r_/, ''));
    let e = anchor ? anchor.nextElementSibling : null;
    while (e && !(e.classList && e.classList.contains('memitem')))
      e = e.nextElementSibling;
    if (e)
      types = Array.from(e.querySelectorAll('td.paramtype'))
        .map(td => bare(td.textContent))
        .filter(t => t.length && t !== 'void');
    else {
      // the anchor of a template is on the line of its template<...>,
      // and its declaration is on the next one
      const row = document.getElementById('r_' + id.replace(/^r_/, ''));
      let cell = row ? row.querySelector('td.memItemRight') : null;
      if (row && !cell && row.nextElementSibling)
        cell = row.nextElementSibling.querySelector('td.memItemRight');
      if (cell) types = declaredTypes(cell.textContent);
    }
    return types ? '(' + types.join(', ') + ')' : null;
  }

  function qualify(nav) {
    const links = Array.from(nav.querySelectorAll('li > div.item > a'))
      .filter(a => a.getAttribute('href') &&
                   a.getAttribute('href').startsWith('#'));
    // the members of each name, a member being listed twice (among the
    // declarations, at #r_<id>, and among the details, at #<id>)
    const members = {};
    links.forEach(a => {
      const n = a.textContent.trim();
      (members[n] = members[n] || new Set())
        .add(a.getAttribute('href').substring(1).replace(/^r_/, ''));
    });
    links.forEach(a => {
      const n = a.textContent.trim();
      if (!members[n] || members[n].size < 2 || a.dataset.qualified) return;
      const types = parameterTypes(a.getAttribute('href').substring(1));
      if (types) {
        a.textContent = n + types;
        a.title = n + types;
      }
      a.dataset.qualified = '1';
    });
  }

  function start() {
    const nav = document.getElementById('page-nav-contents');
    if (!nav) return;
    qualify(nav);
    new MutationObserver(() => qualify(nav))
      .observe(nav, { childList: true, subtree: true });
  }

  if (document.readyState === 'complete') start();
  else window.addEventListener('load', start);
})();
