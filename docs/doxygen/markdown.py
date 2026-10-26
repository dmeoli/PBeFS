#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    markdown.py                                                              #
#                                                                             #
#    The input filter of Doxygen for the Markdown files, which are pages of   #
#    the site of Doxygen (the README is its main page): it prints the file    #
#    with the links that Doxygen cannot follow made links of the site or of   #
#    GitHub, i.e., the citations of the README [n] that another page points   #
#    to at ../README.md#ref-n go to index.html#ref-n, the links to the        #
#    headers go to their pages, and the links to the other files of the       #
#    repository that are not pages go to GitHub. The formulae, written for    #
#    MyST and GitHub between dollars, become those of Doxygen outside the     #
#    code, i.e., $...$ becomes \f$...\f$ and $$...$$ becomes \f[...\f]. Given #
#    a group, it puts "\addtogroup <group>" in place of the title, so that    #
#    the page is the documentation of its group, which goes under its layer   #
#    in the tree of the site, and the links to the page go to the group; the  #
#    FILTER_PATTERNS that give the groups are written by docs/apiref.py into  #
#    build-docs/pages.cfg, which the second run of Doxygen reads (--topics).  #
#    In that run it also filters the headers, whose mentions of the pages of  #
#    the layers (e.g., docs/branch-and-bound.md) become references to their   #
#    groups.                                                                  #
#                                                                             #
#        markdown.py [--topics] [<group>] <file>                              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))))
import apiref  # noqa: E402

# the group of each page of a layer, whose documentation the page is,
# which the links to the page point to when the groups exist, i.e., in the
# run of Doxygen with the Topics (--topics)
GROUPS = apiref.doc_groups()
TOPICS = False
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
BLOB = 'https://github.com/dmeoli/PBeFS/blob/main/'
TREE = 'https://github.com/dmeoli/PBeFS/tree/main/'

# the Markdown files that are pages of the site of Doxygen [see INPUT in
# docs/doxygen/Doxyfile]
PAGES = {'README.md', 'docs/branch-and-bound.md'} | {
    os.path.join('docs', 'theory', f)
    for f in os.listdir(os.path.join(ROOT, 'docs', 'theory'))
    if f.endswith('.md')}


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
        if full in GROUPS and TOPICS:
            return '%s(@ref %s)' % (m.group(1), GROUPS[full])
        if full in PAGES:
            return m.group(0)
        if full.startswith('include/') and full.endswith('.h'):
            return '%s(@ref %s)' % (m.group(1), full[len('include/'):])
        base = TREE if target.endswith('/') else BLOB
        return '%s(%s%s)' % (m.group(1), base, full)
    return re.sub(r'(\]|\]\s)\(([^)\s]+)\)', fix, text)


def formulae(text):
    """the formulae between dollars as those of Doxygen, but in the fenced
    code blocks and in the code spans; a formula may span more lines, but
    not a blank one"""
    parts = re.split(r'(^\s*(?:```|~~~).*?^\s*(?:```|~~~)[^\n]*$|`[^`\n]*`)',
                     text, flags=re.M | re.S)
    for i in range(0, len(parts), 2):
        t = re.sub(r'\$\$((?:(?!\n\s*\n).)+?)\$\$', r'\\f[\1\\f]',
                   parts[i], flags=re.S)
        parts[i] = re.sub(r'(?<![\\\w])\$(?=\S)((?:(?!\n\s*\n)[^$])+?)'
                          r'(?<=\S)\$', r'\\f$\1\\f$', t)
    return ''.join(parts)


def header(text):
    """the header with its mentions of the pages of a layer in the comments
    (e.g., "docs/branch-and-bound.md") made references to their groups"""
    for page, grp in GROUPS.items():
        text = re.sub(r'(?<![\w/])%s\b' % re.escape(page),
                      r'\\ref %s "%s"' % (grp, page), text)
    return text


def main():
    global TOPICS
    args = sys.argv[1:]
    TOPICS = '--topics' in args
    args = [a for a in args if a != '--topics']
    group = args[0] if len(args) > 1 else None
    path = os.path.abspath(args[-1])
    with open(path) as f:
        text = f.read()
    if path.endswith('.h'):
        sys.stdout.write(header(text))
        return
    if group:                       # the title is that of the group
        text = re.sub(r'\A# [^\n]*\n', '\\\\addtogroup %s\n' % group, text)
    sys.stdout.write(formulae(links(text, path)))


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
#                            End of markdown.py                               #
# --------------------------------------------------------------------------- #
