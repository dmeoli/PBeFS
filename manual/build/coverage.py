#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    coverage.py                                                              #
#                                                                             #
#    Checks that the manual accounts for the public interface of the          #
#    library: every header of include/PBeFS, and every class, struct,         #
#    concept and function declared at namespace scope in one (outside the     #
#    namespaces detail and of the implementations in structures/), must be    #
#    named in a chapter, as `name` in code or in a [Source: ...] line, or be  #
#    listed in coverage-omitted.txt with the reason why the manual leaves it  #
#    out. The names not accounted for are printed and the exit status is 1.   #
#                                                                             #
#        python3 manual/build/coverage.py                                     #
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
MANUAL = os.path.dirname(BUILD)
ROOT = os.path.dirname(MANUAL)
OMITTED = os.path.join(BUILD, 'coverage-omitted.txt')

DECL = re.compile(
    r'^(?:template\s*<[^;{]*?>\s*)?'
    r'(?:class|struct|concept)\s+([A-Za-z_]\w*)\b(?!\s*;)', re.M)
FUNC = re.compile(r'(?<![\w:.>~])([a-z_]\w*)\s*\(')
KEYWORDS = {'if', 'for', 'while', 'switch', 'return', 'operator', 'sizeof',
            'decltype', 'requires', 'noexcept', 'alignas', 'alignof',
            'static_assert', 'explicit', 'defined', 'template', 'typename',
            'const', 'constexpr', 'consteval', 'inline', 'static', 'void',
            'bool', 'int', 'unsigned', 'double', 'float', 'char', 'auto',
            'long', 'short', 'size_t', 'nodiscard', 'maybe_unused'}


def strip(text):
    """the text without its comments and without the bodies of the
    namespaces detail"""
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    text = re.sub(r'//[^\n]*', '', text)
    out, depth, skip = [], 0, None
    for line in text.split('\n'):
        if skip is None and re.match(r'\s*namespace\s+detail\b', line):
            skip = depth
        if skip is None:
            out.append(line)
        depth += line.count('{') - line.count('}')
        if skip is not None and depth <= skip and '}' in line:
            skip = None
    return '\n'.join(out)


def toplevel(text):
    """the text at the scope of a namespace, i.e., not inside a class, an
    enumeration or a function (the bodies are dropped)"""
    out, stack, pending = [], [], False
    for tok in re.split(r'(\bnamespace\b|[{};])', text):
        if tok == 'namespace':
            pending = True
        elif tok == '{':
            stack.append('ns' if pending else 'block')
            pending = False
            if stack[-1] == 'ns' or 'block' in stack[:-1]:
                continue
            out.append(tok)
            continue
        elif tok == '}':
            if stack:
                stack.pop()
            if 'block' in stack:
                continue
        elif tok == ';':
            pending = False
        if 'block' not in stack:
            out.append(tok)
    return ''.join(out)


def names():
    """{ name : header } of the public interface"""
    out = {}
    for path in sorted(glob.glob(os.path.join(ROOT, 'include', 'PBeFS', '**',
                                              '*.h'), recursive=True)):
        rel = os.path.relpath(path, ROOT)
        if '/structures/' in rel:
            continue
        out.setdefault(os.path.basename(path), rel)
        text = toplevel(strip(open(path, encoding='utf-8').read()))
        for m in DECL.finditer(text):
            out.setdefault(m.group(1), rel)
        text = re.sub(r'^\s*#.*$', '', text, flags=re.M)
        for m in FUNC.finditer(text):
            if m.group(1) not in KEYWORDS:
                out.setdefault(m.group(1), rel)
    return out


def main():
    text = ''
    for path in glob.glob(os.path.join(MANUAL, 'chapters', '*.md')):
        text += open(path, encoding='utf-8').read()
    omitted = {}
    if os.path.exists(OMITTED):
        for line in open(OMITTED, encoding='utf-8'):
            if line.strip() and not line.startswith('#'):
                name, why = line.rstrip('\n').split('\t', 1)
                omitted[name] = why
    missing = []
    found = names()
    for name, rel in sorted(found.items(), key=lambda x: (x[1], x[0])):
        if name in omitted:
            continue
        if not re.search(r'(?<![\w])%s(?![\w])' % re.escape(name), text):
            missing.append((rel, name))
    for rel, name in missing:
        print('coverage: %s (%s) is not in the manual' % (name, rel))
    print('coverage: %d names, %d omitted with a reason, %d unaccounted'
          % (len(found), len(omitted), len(missing)))
    sys.exit(1 if missing else 0)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------- #
# ----------------------------- END coverage.py ----------------------------- #
# --------------------------------------------------------------------------- #
