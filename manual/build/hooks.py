# --------------------------------------------------------------------------- #
#                                 hooks.py                                    #
#                                                                             #
# The hooks of MkDocs for the web site of the manual. Before the Markdown of  #
# a page is converted, (i) each displayed formula between two lines $$ gets  #
# a blank line before and after, so that pymdownx.arithmatex takes it as a   #
# block and hands it to KaTeX untouched (Pandoc, which builds the PDF, also  #
# takes a display inside a paragraph, and the Markdown of MkDocs would read  #
# \{, * and _ in it), and (ii) the name of a file or folder of the           #
# repository written as code, e.g., `lemon/dijkstra.h`, becomes a link to it #
# on GitHub, as build/source-links.lua does in the PDF.                      #
# --------------------------------------------------------------------------- #

import os
import re

REPO = 'https://github.com/dmeoli/PBeFS'  # the repository
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
SPAN = re.compile(r'(?<![\[`])`([A-Za-z0-9_.][A-Za-z0-9_./-]*)`(?![`\]])')


def source_url(path):
    """the GitHub URL of a file or folder of the repository, or None"""
    if '/' not in path or '..' in path:
        return None
    full = os.path.join(ROOT, path)
    if os.path.isfile(full):
        return f'{REPO}/blob/main/{path}'
    if os.path.isdir(full):
        return f'{REPO}/tree/main/{path.rstrip("/")}'
    return None


def link_sources(line):
    def repl(m):
        url = source_url(m.group(1))
        return f'[`{m.group(1)}`]({url})' if url else m.group(0)
    return SPAN.sub(repl, line)


def on_page_markdown(markdown, **kwargs):
    out, fence, disp = [], None, False
    for line in markdown.split('\n'):
        s = line.strip()
        m = re.match(r'(`{3,}|~{3,})', s)
        if fence:
            out.append(line)
            if m and s.startswith(fence):
                fence = None
            continue
        if m:
            fence = m.group(1)
            out.append(line)
            continue
        if not disp and len(s) > 4 and s.startswith('$$') and \
           s.endswith('$$'):
            ind = line[:len(line) - len(line.lstrip())]
            if out and out[-1].strip():
                out.append(ind)
            out += [line, ind]
            continue
        if s == '$$':
            ind = line[:len(line) - len(line.lstrip())]
            if not disp:
                if out and out[-1].strip():
                    out.append(ind)
                out.append(line)
            else:
                out.append(line)
                out.append(ind)
            disp = not disp
            continue
        out.append(line if disp else link_sources(line))
    return '\n'.join(out)
