#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    markdown.py                                                              #
#                                                                             #
#    The input filter of Doxygen for the Markdown files, which are pages of   #
#    the site of Doxygen (the README is its main page): it prints the file    #
#    with the links that Doxygen cannot follow made links of the site or of   #
#    GitHub, i.e., the citations of the README [n] that another page points   #
#    to at ../README.md#ref-n go to index.html#ref-n, and the links to the    #
#    files of the repository that are not pages go to GitHub. Given a group,  #
#    it ends the page with "\ingroup <group>", so that the page goes under    #
#    its layer in the tree of the site instead of at its top; the             #
#    FILTER_PATTERNS that give the groups are written by docs/apiref.py into  #
#    build-docs/pages.cfg, which the second run of Doxygen reads.             #
#                                                                             #
#        markdown.py [<group>] <file>                                         #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
BLOB = 'https://github.com/dmeoli/PBeFS/blob/main/'
TREE = 'https://github.com/dmeoli/PBeFS/tree/main/'

# the Markdown files that are pages of the site of Doxygen [see INPUT in
# docs/doxygen/Doxyfile]
PAGES = {'README.md', 'docs/branch-and-bound.md'}


def links(text, path):
    def fix(m):
        target = m.group(2)
        if re.match(r'[a-z]+:', target) or target.startswith('#'):
            return m.group(0)
        name, _, anchor = target.partition('#')
        full = os.path.relpath(os.path.normpath(os.path.join(
            os.path.dirname(path), name)), ROOT)
        if full == 'README.md' and anchor:
            return '%s(index.html#%s)' % (m.group(1), anchor)
        if full in PAGES:
            return m.group(0)
        base = TREE if target.endswith('/') else BLOB
        return '%s(%s%s)' % (m.group(1), base, full)
    return re.sub(r'(\]|\]\s)\(([^)\s]+)\)', fix, text)


def main():
    group = sys.argv[1] if len(sys.argv) > 2 else None
    path = os.path.abspath(sys.argv[-1])
    with open(path) as f:
        sys.stdout.write(links(f.read(), path))
    if group:
        sys.stdout.write('\n\\ingroup %s\n' % group)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
#                            End of markdown.py                               #
# --------------------------------------------------------------------------- #
