#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    srclinks.py                                                              #
#                                                                             #
#    Makes every file of the repository named in the manual a link to that    #
#    file on GitHub. A span `path`, `path:n` or `path:a-b` of the chapters    #
#    and of index.md (outside the blocks of code) that names a file or a      #
#    directory of the tree, by its whole path or by a unique final part of    #
#    it (a basename, a name without .cpp for an executable), becomes the      #
#    reference link [`path`][name], defined at the end of the chapter as      #
#    URL/blob/main/<path>#L<a>-L<b> (/tree/main/ for a directory); a line     #
#    made longer than 79 columns gives its last words to the next one.        #
#    With --check nothing is written: the exit status is 1 if a span that     #
#    looks like a path is not a link, or if a link names a path missing in    #
#    the tree or a line past the end of its file; the tables written by       #
#    tables.py in manual/generated are checked as well.                       #
#                                                                             #
#        python3 manual/build/srclinks.py [--check]                           #
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

BUILD = os.path.dirname(os.path.abspath(__file__))
MANUAL = os.path.dirname(BUILD)
ROOT = os.path.dirname(MANUAL)
URL = 'https://github.com/dmeoli/PBeFS'

EXT = ('h', 'hpp', 'cpp', 'cc', 'txt', 'md', 'py', 'sh', 'yml', 'yaml',
       'cmake', 'tex', 'csv', 'lua', 'js', 'css', 'json', 'in')
TOP = ('include/', 'src/', 'test/', 'tools/', 'examples/', 'data/',
       'docs/', 'manual/', 'cmake/', '.github/')
PATHLIKE = re.compile(
    r'^(?P<path>[\w.*/+-]+?)(?::(?P<a>\d+)(?:-(?P<b>\d+))?)?$')
SPAN = re.compile(r'(?<![`\\])`(?P<x>[^`\n]+)`(?!`)')
LINK = re.compile(re.escape(URL) +
                  r'/(?P<kind>blob|tree)/main/(?P<path>[^)#\s]*)'
                  r'(?:#L(?P<a>\d+)(?:-L(?P<b>\d+))?)?\)')
FENCE = re.compile(r'^\s*(```|~~~)')

FILES = subprocess.run(['git', 'ls-files'], cwd=ROOT, check=True,
                       capture_output=True, text=True).stdout.split()
FILESET = set(FILES)
DIRS = {os.path.dirname(f) for f in FILES}
DIRS |= {d.rsplit('/', k)[0] for d in list(DIRS) for k in range(1, 6)
         if d.count('/') >= k}
LINES = {}


def lines(path):
    if path not in LINES:
        with open(os.path.join(ROOT, path), encoding='utf-8',
                  errors='replace') as f:
            LINES[path] = sum(1 for _ in f)
    return LINES[path]


def pathlike(x):
    m = PATHLIKE.match(x)
    if not m:
        return None
    p = m.group('path')
    if p.startswith(TOP) or re.search(r'\.(%s)$' % '|'.join(EXT), p):
        return m
    return None


def resolve(p):
    """(kind, path) of the tree named by p, or None"""
    if '*' in p:
        d = os.path.dirname(p)
        return ('tree', d) if d in DIRS else None
    q = p.rstrip('/')
    if q in FILESET:
        return 'blob', q
    if q in DIRS:
        return 'tree', q
    if q + '.cpp' in FILESET:
        return 'blob', q + '.cpp'
    for name in (q, q + '.cpp'):
        c = [f for f in FILES if f.endswith('/' + name)]
        if c:
            c.sort(key=lambda f: (f.count('/'), f))
            if len(c) > 1 and c[0].count('/') == c[1].count('/'):
                return None
            return 'blob', c[0]
    c = sorted((d for d in DIRS if d.endswith('/' + q)),
               key=lambda d: d.count('/'))
    if len(c) == 1:
        return 'tree', c[0]
    return None


def target(m):
    r = resolve(m.group('path'))
    if r is None:
        return None, 'no file of the tree is %s' % m.group('path')
    kind, path = r
    a, b = m.group('a'), m.group('b')
    url = '%s/%s/main/%s' % (URL, kind, path)
    if a:
        if kind != 'blob' or int(b or a) > lines(path) or \
           (b and int(b) < int(a)):
            return None, 'line %s past the end of %s' % (b or a, path)
        url += '#L' + a + ('-L' + b if b else '')
    return url, None


MARK = '<!-- the links to the sources, written by srclinks.py -->'


def reflow(lines, changed):
    """gives the last words of a line made too long to the next line of
    its paragraph, or to a new line when it is the last one"""
    def scan(line, math):
        # the spaces outside code spans and math, and the state of the
        # math at the end of the line
        code = False
        spaces = []
        for i, ch in enumerate(line):
            if ch == '`':
                code = not code
            elif ch == '$' and not code:
                math = not math
            elif ch == ' ' and not code and not math and \
                    line[:i].strip():
                spaces.append(i)
        return spaces, math

    def split(line, math):
        # the last such space within 79 columns, else the first one
        spaces = scan(line, math)[0]
        inside = [i for i in spaces if i <= 79]
        return inside[-1] if inside else (spaces[0] if spaces else None)

    def plain(line):
        t = line.lstrip()
        return t != '' and not t.startswith(('|', '#', '<', '--8<--',
                                               '![', '$$', '[^', ':::',
                                               '!!!', '???'))

    item = re.compile(r'(\s*)([-*+]|\d+\.)(\s+)')
    out = []
    math = False
    i = 0
    while i < len(lines):
        line = lines[i]
        if i in changed and len(line) > 79 and plain(line):
            k = split(line, math)
            if k is not None:
                rest = line[k + 1:]
                head = line[:k]
                nxt = lines[i + 1] if i + 1 < len(lines) else ''
                if plain(nxt) and not item.match(nxt):
                    indent = re.match(r'\s*', nxt).group(0)
                    lines[i + 1] = indent + rest + ' ' + nxt.lstrip()
                else:
                    m = item.match(line)
                    indent = ' ' * len(m.group(0)) if m else \
                        re.match(r'\s*', line).group(0)
                    lines.insert(i + 1, indent + rest)
                    changed = {j + 1 if j > i else j for j in changed}
                changed.add(i + 1)
                line = head
        out.append(line)
        math = scan(line, math)[1] if line.strip() and \
            '$$' not in line else False
        i += 1
    return out


def process(text, name, errors, count):
    if MARK in text:
        text = text[:text.index(MARK)].rstrip('\n') + '\n'
    out = []
    defs = {}
    labels = {}
    changed = set()
    fenced = False
    for no, line in enumerate(text.split('\n'), 1):
        if FENCE.match(line):
            fenced = not fenced
        if fenced or FENCE.match(line):
            out.append(line)
            continue
        for m in LINK.finditer(line):
            p = m.group('path').rstrip('/')
            ok = p in FILESET if m.group('kind') == 'blob' else p in DIRS
            if not ok:
                errors.append('%s:%d: the link names %s, missing in the '
                              'tree' % (name, no, p))

        def sub(s):
            x = s.group('x')
            pm = pathlike(x)
            before = line[s.start() - 1:s.start()]
            after = line[s.end():s.end() + 2]
            if pm is None or (before == '[' and after == ']('):
                return s.group(0)
            url, err = target(pm)
            if url is None:
                errors.append('%s:%d: `%s`: %s' % (name, no, x, err))
                return s.group(0)
            if before == '[' and after == '][':
                # a link of an earlier run: its label is kept
                label = line[s.end() + 2:line.index(']', s.end() + 2)]
                defs[label.lower()] = (label, url)
                labels[label.lower()] = url
                count[0] += 1
                return s.group(0)
            # the shortest label that names this link only
            for label in (os.path.basename(x.rstrip('/')) or x, x,
                          url[len(URL) + 1:]):
                key = label.lower()
                if labels.get(key, url) == url:
                    break
            labels[key] = url
            defs[key] = (label, url)
            count[0] += 1
            return '[`%s`][%s]' % (x, label)
        new = SPAN.sub(sub, line)
        if new != line:
            changed.add(len(out))
        out.append(new)
    out = reflow(out, changed)
    text = '\n'.join(out).rstrip('\n') + '\n'
    if defs:
        text += '\n' + MARK + '\n\n' + ''.join(
            '[%s]: %s\n' % d for d in sorted(defs.values(),
                                             key=lambda d: d[0].lower()))
    return text


def main():
    check = '--check' in sys.argv[1:]
    errors = []
    files = sorted(glob.glob(os.path.join(MANUAL, 'chapters', '*.md')))
    files.append(os.path.join(MANUAL, 'index.md'))
    # the tables that tables.py writes are checked only: their links are
    # written there, inline
    tables = sorted(glob.glob(os.path.join(MANUAL, 'generated', 'tables',
                                           '*.md')))
    for f in files + tables:
        with open(f, encoding='utf-8') as h:
            text = h.read()
        count = [0]
        rel = os.path.relpath(f, ROOT)
        new = process(text, rel, errors, count)
        if new != text:
            if check or f in tables:
                errors.append('%s: some files of the tree are not links '
                              '(run manual/build/srclinks.py)' % rel)
            else:
                with open(f, 'w', encoding='utf-8') as h:
                    h.write(new)
        if not check and count[0]:
            print('%4d  %s' % (count[0], os.path.relpath(f, MANUAL)))
    for e in errors:
        print(e)
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())

# ----------------------------- END srclinks.py ----------------------------- #
