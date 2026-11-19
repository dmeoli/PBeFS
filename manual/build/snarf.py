#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    snarf.py                                                                 #
#                                                                             #
#    Writes into manual/generated/snippets the excerpts of code of the        #
#    manual: each excerpt between two markers //! [name] of a header or of    #
#    the tester test/manual.cpp (the examples of use, which are thus          #
#    compiled and run), the markers that the command \snippet of Doxygen      #
#    reads, goes to <name>.cpp without its common indentation, and a          #
#    chapter includes it by the line --8<-- "snippets/<name>.cpp" in a block  #
#    of code (see include.lua), so that the code is written once. A marker    #
#    opened and not closed, or a name used twice, is an error.                #
#                                                                             #
#        python3 manual/build/snarf.py                                        #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

import glob
import os
import re
import sys

MANUAL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(MANUAL)
OUT = os.path.join(MANUAL, 'generated')
MARKER = re.compile(r'^\s*//! \[([A-Za-z0-9_.-]+)\]\s*$')


def excerpts(path):
    """{ name : lines } of the markers of a file"""
    found, open_ = {}, {}
    with open(path, encoding='utf-8') as f:
        for n, line in enumerate(f, 1):
            m = MARKER.match(line)
            if not m:
                for lines in open_.values():
                    lines.append(line.rstrip('\n'))
                continue
            name = m.group(1)
            if name in open_:
                found[name] = open_.pop(name)
            elif name in found:
                sys.exit('%s:%d: the excerpt %s is there twice'
                         % (path, n, name))
            else:
                open_[name] = []
    if open_:
        sys.exit('%s: the excerpts %s are not closed'
                 % (path, ', '.join(sorted(open_))))
    return found


def dedent(lines):
    """the lines without the indentation common to the ones not empty"""
    full = [l for l in lines if l.strip()]
    cut = min((len(l) - len(l.lstrip(' ')) for l in full), default=0)
    return [l[cut:] for l in lines]



def main():
    os.makedirs(os.path.join(OUT, 'snippets'), exist_ok=True)
    seen = {}
    paths = sorted(glob.glob(os.path.join(ROOT, 'include', '**', '*.h'),
                             recursive=True))
    paths += sorted(glob.glob(os.path.join(ROOT, 'test', '*.cpp')))
    for path in paths:
        for name, lines in excerpts(path).items():
            if name in seen:
                sys.exit('the excerpt %s is in both %s and %s'
                         % (name, seen[name], path))
            seen[name] = path
            with open(os.path.join(OUT, 'snippets', name + '.cpp'), 'w',
                      encoding='utf-8') as f:
                f.write('\n'.join(dedent(lines)) + '\n')
    print('snarf.py: %d excerpts' % len(seen))


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
# ------------------------------ END snarf.py ------------------------------- #
# --------------------------------------------------------------------------- #
