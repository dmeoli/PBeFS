#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    snarf.py                                                                 #
#                                                                             #
#    Writes into manual/generated the pieces of the manual that come from     #
#    the code: each excerpt of a header between two markers //! [name],       #
#    the same markers that the command \snippet of Doxygen reads, goes to     #
#    generated/snippets/<name>.cpp, without its common indentation, so that   #
#    the code of a listing is written once, in the header; and the version    #
#    of the library, the last git tag or "unreleased", to                     #
#    generated/version.tex. A marker opened and not closed, or a name used    #
#    twice, is an error.                                                      #
#                                                                             #
#        python3 manual/snarf.py                                              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

import glob
import os
import re
import subprocess
import sys

MANUAL = os.path.dirname(os.path.abspath(__file__))
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


def version():
    """the last git tag, without its v, or unreleased"""
    try:
        return subprocess.run(['git', 'describe', '--tags', '--abbrev=0'],
                              cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip().lstrip('v')
    except (OSError, subprocess.CalledProcessError):
        return 'unreleased'


def main():
    os.makedirs(os.path.join(OUT, 'snippets'), exist_ok=True)
    seen = {}
    for path in sorted(glob.glob(os.path.join(ROOT, 'include', '**', '*.h'),
                                 recursive=True)):
        for name, lines in excerpts(path).items():
            if name in seen:
                sys.exit('the excerpt %s is in both %s and %s'
                         % (name, seen[name], path))
            seen[name] = path
            with open(os.path.join(OUT, 'snippets', name + '.cpp'), 'w',
                      encoding='utf-8') as f:
                f.write('\n'.join(dedent(lines)) + '\n')
    with open(os.path.join(OUT, 'version.tex'), 'w') as f:
        v = version()
        f.write(('Version %s' % v if v != 'unreleased' else 'Unreleased')
                + '\n')
    print('snarf.py: %d excerpts' % len(seen))


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
# ------------------------------ END snarf.py ------------------------------- #
# --------------------------------------------------------------------------- #
