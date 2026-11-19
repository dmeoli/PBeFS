#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    refcheck.py                                                              #
#                                                                             #
#    Checks the cross-references of the manual: every link whose text names   #
#    a numbered object, e.g., [Lemma 3.6](03-data-structures.md#lem-3-6),     #
#    must point to an anchor that exists in that file, whose number is the    #
#    one of the text, and whose labelled paragraph (or heading, or figure)    #
#    is of that kind and has that number, e.g., **Lemma 3.6** (an exercise    #
#    or a table may point to the section that holds it). A link that          #
#    fails is printed and the exit status is 1, so that a renumbering does    #
#    not leave a reference that resolves to the wrong statement.              #
#                                                                             #
#        python3 manual/build/refcheck.py                                     #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

import glob
import os
import re
import sys

BUILD = os.path.dirname(os.path.abspath(__file__))
CHAPTERS = os.path.join(os.path.dirname(BUILD), 'chapters')

KINDS = ('Lemma', 'Theorem', 'Proposition', 'Corollary', 'Definition',
         'Remark', 'Algorithm', 'Exercise', 'Figure', 'Table', 'Chapter',
         'Example', 'Observation')
LINK = re.compile(r'\[(%s)s? ([0-9A-Z]+(?:\.[0-9]+)?)\]'
                  r'\(([\w.-]*\.md)?#([\w-]+)\)' % '|'.join(KINDS))


def objects(text):
    """{ anchor : ( kind , number ) } of the labelled objects of a file"""
    out = {}
    for m in re.finditer(r'<a id="([\w-]+)"></a>\s*\*\*(%s) ([0-9A-Z.]+?)'
                         r'\.?\*\*' % '|'.join(KINDS), text):
        out[m.group(1)] = (m.group(2), m.group(3))
    for m in re.finditer(r'^#+ (?:Appendix )?([0-9A-Z]+(?:\.[0-9]+)?)\.? .*?'
                         r'\{#([\w-]+)\}', text, re.M):
        out[m.group(2)] = ('Chapter' if '.' not in m.group(1) else 'Section',
                           m.group(1))
    for m in re.finditer(r'!\[(Figure|Table) ([0-9.]+?)\.? .*?\]\([^)]*\)'
                         r'\{#([\w-]+)\}', text, re.S):
        out[m.group(3)] = (m.group(1), m.group(2))
    for m in re.finditer(r'<a id="([\w-]+)"></a>\s*\*\*(Table) ([0-9.]+?)'
                         r'\.?\*\*', text):
        out[m.group(1)] = (m.group(2), m.group(3))
    return out


def main():
    files = sorted(glob.glob(os.path.join(CHAPTERS, '*.md')))
    objs = {os.path.basename(f): objects(open(f, encoding='utf-8').read())
            for f in files}
    bad = 0
    for f in files:
        name = os.path.basename(f)
        for n, line in enumerate(open(f, encoding='utf-8'), 1):
            for m in LINK.finditer(line):
                kind, num, target, anchor = m.groups()
                where = objs.get(target or name)
                got = where.get(anchor) if where is not None else None
                if got is None:
                    print('refcheck: %s:%d: %s points to a missing anchor'
                          % (name, n, m.group(0)))
                    bad += 1
                elif got[0] == 'Section' and kind in ('Exercise', 'Table'):
                    continue          # the section that holds it
                elif got[1] != num or got[0] != kind:
                    print('refcheck: %s:%d: %s points to %s %s'
                          % (name, n, m.group(0), got[0], got[1]))
                    bad += 1
    print('refcheck: %d bad references' % bad)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
# ----------------------------- END refcheck.py ----------------------------- #
# --------------------------------------------------------------------------- #
