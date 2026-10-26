#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    pages.py                                                                 #
#                                                                             #
#    Writes the pages of the site of the documentation into docs/generated,   #
#    from the README, which stays the one source of the text: each section    #
#    of the README becomes a page, the citations [n] of the README point to   #
#    the page of the references, the links to the files of the repository     #
#    point to GitHub, the design note of the branch-and-bound becomes a page  #
#    of the user guide and the pages of the theory in docs/theory pages of    #
#    their own, their links to the headers going to the API reference, and    #
#    the CHANGELOG becomes the release notes; then docs/apiref.py writes the  #
#    pages of the API reference. The workflow .github/workflows/docs.yml      #
#    runs it before Sphinx.                                                   #
#                                                                             #
#        python3 docs/pages.py && sphinx-build -b html docs build-docs/html   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

import os
import re

import apiref

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'generated')
BLOB = 'https://github.com/dmeoli/PBeFS/blob/main/'
TREE = 'https://github.com/dmeoli/PBeFS/tree/main/'

# the sections of the README and the page each one goes to; the sections
# not listed here go to the home page
PAGES = {
    'The search problem': 'user-guide/search-problem.md',
    'The algorithms': 'user-guide/algorithms.md',
    'Branch-and-bound': 'user-guide/branch-and-bound.md',
    'The data structures': 'user-guide/data-structures.md',
    'Building and using the library': 'getting-started/installation.md',
    'The examples and the benchmark': 'user-guide/examples.md',
    'Testing': 'user-guide/testing.md',
    'References': 'user-guide/references.md',
}
REFERENCES = 'user-guide/references.md'
THEORY = os.path.join(ROOT, 'docs', 'theory')
# the Markdown files of docs and their pages
NOTES = {'docs/branch-and-bound.md': 'user-guide/branch-and-bound-design.md'}
NOTES.update(('docs/theory/' + f, 'theory/' + f) for f in os.listdir(THEORY)
             if f.endswith('.md'))


def links(text, page, source='README.md'):
    """the links of a Markdown file of the repository (the README, by
    default) for the page: the citations [n] to the targets of the
    references, the design note and the theory to their pages, the headers
    to their pages of the API reference, the other files of the repository
    on GitHub"""
    up = '../' * page.count('/')
    text = re.sub(r'\[\[(\d+)\]\]\((?:(?:\.\./)*README\.md)?#(ref-\d+)\)',
                  r'{ref}`[\1] <\2>`', text)

    def fix(m):
        target = m.group(2)
        if target.startswith('#ref-'):
            if page == REFERENCES:
                return m.group(0)
            return '%s(%s%s%s)' % (m.group(1), up, REFERENCES, target)
        if re.match(r'[a-z]+:', target) or target.startswith('#'):
            return m.group(0)
        name, _, anchor = target.partition('#')
        full = os.path.normpath(os.path.join(os.path.dirname(source), name))
        if full == 'README.md' and anchor.startswith('ref-'):
            return '%s(%s%s#%s)' % (m.group(1), up, REFERENCES, anchor)
        if full in NOTES:
            return '%s(%s%s)' % (m.group(1), up, NOTES[full])
        header = re.match(r'include/PBeFS/(\w+)/(\w+)\.h$', full)
        if header:
            return '%s(%s../api/%s/header_%s_%s_h.rst)' % (
                m.group(1), up, header.group(1), header.group(1),
                header.group(2))
        base = TREE if target.endswith('/') else BLOB
        return '%s(%s%s)' % (m.group(1), base, full)
    return re.sub(r'(\]|\]\s)\(([^)]+)\)', fix, text)


def write(path, text):
    full = os.path.join(OUT, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, 'w') as f:
        f.write(text)


def main():
    readme = open(os.path.join(ROOT, 'README.md')).read()
    parts = re.split(r'^## ', readme, flags=re.M)
    home, pages = [parts[0]], {}
    for part in parts[1:]:
        title, _, body = part.partition('\n')
        page = PAGES.get(title.strip())
        text = '# %s\n%s' % (title.strip(), body)
        if page:
            pages[page] = pages.get(page, '') + text
        else:
            home.append('## ' + part)
    write('index.md', links(''.join(home), 'index.md'))
    for page, text in pages.items():
        if page == REFERENCES:          # the anchors as targets of MyST
            text = re.sub(r'^<a id="(ref-\d+)"></a>', r'(\1)=\n', text,
                          flags=re.M)
        write(page, links(text, page))
    for source, page in NOTES.items():
        write(page, links(open(os.path.join(ROOT, source)).read(), page,
                          source))
    log = open(os.path.join(ROOT, 'CHANGELOG.md')).read()
    write('release_notes.md', links(log.replace('# Changelog',
                                                '# Release notes', 1),
                                    'release_notes.md'))
    n = apiref.main()
    print('Doxygen: %s' % ('not run' if n is None else '%d warnings' % n))


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
#                              End of pages.py                                #
# --------------------------------------------------------------------------- #
