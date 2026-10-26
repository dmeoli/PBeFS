#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    apiref.py                                                                #
#                                                                             #
#    Writes the API reference of the site into docs/api, from the XML that    #
#    Doxygen extracts from the headers into build-docs/xml, as a tree that    #
#    mirrors the library: a page with its three layers (core, astar and       #
#    bnb), a page per layer, in it a page per group of headers (the           #
#    searches, the domains, the runtimes, ...), a page per header and a page  #
#    per class, each rendered by Breathe. The alternative sets of data        #
#    structures of A* (astar/structures/<set>) are left out of the XML,       #
#    since their classes have the names of those of the library, which        #
#    Doxygen would merge; Structures.h describes them. Doxygen writes its     #
#    warnings into build-docs/doxygen-warnings.log, which docs.yml wants      #
#    empty. pages.py calls it, before Sphinx.                                 #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

import collections
import glob
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INCLUDE = os.path.join(ROOT, 'include', 'PBeFS')
BUILD = os.path.join(ROOT, 'build-docs')
XML = os.path.join(BUILD, 'xml')
WARNINGS = os.path.join(BUILD, 'doxygen-warnings.log')
OUT = os.path.join(ROOT, 'docs', 'api')
NAMESPACE = 'PBeFS_di_unipi_it::'
BLOB = 'https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/'

DOXYFILE = '\n'.join([
    'INPUT = "%s"' % INCLUDE,
    'EXCLUDE = "%s"' % os.path.join(INCLUDE, 'astar', 'structures'),
    'FILE_PATTERNS = *.h',
    'RECURSIVE = YES',
    'EXCLUDE_SYMBOLS = detail::* PBeFS_di_unipi_it::detail::*',
    'EXTRACT_ALL = YES',
    'JAVADOC_AUTOBRIEF = YES',
    'BUILTIN_STL_SUPPORT = YES',
    'FULL_PATH_NAMES = YES',
    'STRIP_FROM_PATH = "%s"' % os.path.join(ROOT, 'include'),
    'GENERATE_HTML = NO',
    'GENERATE_LATEX = NO',
    'GENERATE_XML = YES',
    'XML_OUTPUT = "%s"' % XML,
    'QUIET = YES',
    'WARNINGS = YES',
    'WARN_LOGFILE = "%s"' % WARNINGS,
])

# the layers, in the order of the pages, with their line
LAYERS = [
    ('core', 'The runtimes of the workers, the messages between them and the '
             'helpers that do not depend on the kind of search.'),
    ('astar', 'The graph searches of A*, sequential and parallel, with their '
              'domains, work distributions and data structures.'),
    ('bnb', 'The tree searches of branch-and-bound, sequential and parallel, '
            'and their choice by scheme.'),
]

# the group of each header of a layer, in the order of the pages; a header
# not listed goes to the group of the helpers
GROUPS = {
    'core': [
        ('Runtimes', ['ThreadRuntime.h', 'OpenMPRuntime.h',
                      'FastFlowRuntime.h']),
        ('Messages', ['ThreadMailbox.h', 'MPIMailbox.h']),
    ],
    'astar': [
        ('Searches', ['AStar.h', 'HDAStar.h', 'SPAStar.h', 'SafePBNF.h',
                      'PASE.h', 'FarmAStar.h', 'HDAStarMPI.h']),
        ('Work distributions', ['Distribution.h']),
        ('Domains', ['Domain.h', 'SlidingTile.h', 'Grid.h',
                     'ExplicitGraph.h', 'SlowDomain.h', 'SequenceScore.h']),
        ('Data structures', ['Node.h', 'NodeTable.h', 'OpenList.h',
                             'Structures.h']),
        ('Results', ['Solution.h']),
    ],
    'bnb': [
        ('Searches', ['BranchAndBound.h', 'ParallelBranchAndBound.h',
                      'DeterministicBranchAndBound.h',
                      'StealingBranchAndBound.h', 'SelfSplitBranchAndBound.h',
                      'FarmBranchAndBound.h', 'MPIBranchAndBound.h']),
        ('Schemes', ['BranchAndBoundScheme.h']),
    ],
}
HELPERS = 'Helpers'
UMBRELLA = 'Headers of the layer'


# ----- reading the XML ----------------------------------------------------- #

def text(node):
    """the text of a node of the XML, on one line"""
    if node is None:
        return ''
    return ' '.join(''.join(node.itertext()).split())


def first_sentence(line):
    m = re.search(r'(?<!\be\.g)(?<!\bi\.e)(?<!\bcf)\.(\s|$)', line)
    return line[:m.start() + 1] if m else line


def description(cd):
    """the brief description of a compound, else the first sentence of the
    first paragraph of its detailed one"""
    brief = text(cd.find('briefdescription'))
    if brief:
        return brief
    detail = cd.find('detaileddescription')
    if detail is None:
        return ''
    for para in detail.findall('para'):
        para = ET.fromstring(ET.tostring(para))
        for sub in list(para):
            if sub.tag in ('simplesect', 'parameterlist', 'xrefsect'):
                para.remove(sub)
        line = text(para)
        if line:
            return first_sentence(line)
    return ''


def short(name):
    return name[len(NAMESPACE):] if name.startswith(NAMESPACE) else name


class Compound:
    def __init__(self, cd):
        self.refid = cd.get('id')
        self.kind = cd.get('kind')
        self.name = cd.findtext('compoundname')
        self.short = short(self.name)
        loc = cd.find('location')
        self.file = loc.get('file') if loc is not None else ''
        self.line = int(loc.get('line') or 0) if loc is not None else 0
        self.brief = description(cd)
        self.inner = [i.get('refid') for i in cd.findall('innerclass')]


def read_xml():
    """the compounds, and the members outside the classes, which lie in the
    compounds of the namespaces and of the files"""
    compounds, members = {}, {}
    for path in glob.glob(os.path.join(XML, '*.xml')):
        if os.path.basename(path) in ('index.xml', 'Doxyfile.xml'):
            continue
        for cd in ET.parse(path).getroot().findall('compounddef'):
            c = Compound(cd)
            compounds[c.refid] = c
            if c.kind not in ('namespace', 'file'):
                continue
            for m in cd.iter('memberdef'):
                loc = m.find('location')
                if loc is None or m.get('id') in members:
                    continue
                members[m.get('id')] = {
                    'kind': m.get('kind'), 'file': loc.get('file'),
                    'line': int(loc.get('line') or 0),
                    'name': m.findtext('name'),
                    'qualified': m.findtext('qualifiedname')
                    or m.findtext('name'),
                    'args': m.findtext('argsstring') or ''}
    return compounds, members


# ----- the headers --------------------------------------------------------- #

class Header:
    def __init__(self, f, path, classes):
        self.path = path                     # relative to include/PBeFS
        self.name = os.path.basename(path)
        self.stem = os.path.splitext(self.name)[0]
        self.refid = f.refid
        self.brief = f.brief
        self.classes = classes
        self.members = []
        # the path, so that astar.h and astar/AStar.h stay apart on a file
        # system that ignores the case
        self.page = 'header_' + re.sub(r'\W', '_', path)
        layer = path.split('/', 1)[0]
        self.layer = layer if '/' in path else self.stem
        self.group = UMBRELLA if '/' not in path else HELPERS
        for g, names in GROUPS.get(self.layer, []):
            if '/' in path and self.name in names:
                self.group = g
                self.order = names.index(self.name)
        if not hasattr(self, 'order'):
            self.order = 0 if self.group == UMBRELLA else 99


def headers_of(compounds, members):
    nested = {i for c in compounds.values() for i in c.inner
              if c.kind in ('class', 'struct', 'union')}
    by_file = collections.defaultdict(list)
    for c in compounds.values():
        if c.kind in ('class', 'struct', 'union') and c.refid not in nested:
            by_file[c.file].append(c)
    headers = []
    for f in compounds.values():
        if f.kind != 'file' or not f.file.endswith('.h'):
            continue
        path = os.path.relpath(f.file, INCLUDE) if os.path.isabs(f.file) \
            else re.sub(r'^PBeFS/', '', f.file)
        h = Header(f, path, sorted(by_file.get(f.file, []),
                                   key=lambda c: (c.short != f.name[:-2],
                                                  c.line)))
        h.members = sorted((m for m in members.values()
                            if m['file'] == f.file
                            and '<' not in m['qualified']),
                           key=lambda m: m['line'])
        headers.append(h)
    return headers


# ----- the pages ----------------------------------------------------------- #

def escape(line):
    return re.sub(r'([*_|`\\\[\]<>])', r'\\\1', line)


def capital(line):
    return line[:1].upper() + line[1:]


def title(line, char='='):
    return '%s\n%s\n\n' % (line, char * len(line))


def rubric(line):
    # a heading that stays out of the table of contents
    return '.. rubric:: %s\n\n' % line


def toctree(entries):
    return ('.. toctree::\n   :maxdepth: 1\n   :hidden:\n\n' +
            ''.join('   %s\n' % e for e in entries) + '\n')


def table(rows, header):
    if not rows:
        return ''
    out = '.. list-table::\n   :header-rows: 1\n   :widths: 30 70\n\n'
    for row in [header] + rows:
        out += '   * - %s\n     - %s\n' % (row[0], row[1] or ' ')
    return out + '\n'


def summary(h):
    """the line of a header in the tables: the brief description of its
    main class, else the first sentence of the comment of the file without
    its "Header file for" """
    main = [c for c in h.classes if c.short == h.stem and c.brief]
    if main:
        return capital(main[0].brief)
    return capital(re.sub(r'^Header file (for|that includes) (the )?', '',
                          first_sentence(h.brief)))


def class_page(c, h):
    kind = 'struct' if c.kind == 'struct' else 'class'
    out = title(escape(c.short))
    out += 'Declared in :doc:`%s <%s>`.\n\n' % (h.path, h.page)
    out += ('.. doxygen%s:: %s\n   :members:\n   :protected-members:\n'
            '   :undoc-members:\n' % (kind, c.name))
    return out


def member_directive(m):
    name = m['qualified']
    if m['kind'] == 'function':
        args = re.match(r'\(.*\)', m['args'])
        return '.. doxygenfunction:: %s%s\n' % (name,
                                                 args.group(0) if args else '')
    kind = {'typedef': 'typedef', 'enum': 'enum', 'variable': 'variable',
            'define': 'define'}.get(m['kind'])
    return '.. doxygen%s:: %s\n' % (kind, name) if kind else ''


def header_page(h, pages):
    out = title(escape(h.path))
    out += ('Header of the layer :doc:`%s <index>`, in the group :doc:`%s '
            '<%s>`; its source is on `GitHub <%s%s>`__.\n\n'
            % (h.layer, h.group, group_page(h.group), BLOB, h.path))
    if h.classes:
        out += rubric('Classes')
        out += table([(':doc:`%s <%s>`' % (escape(c.short), pages[c.refid]),
                       escape(capital(c.brief))) for c in h.classes],
                     ('Class', 'Description'))
        out += toctree([pages[c.refid] for c in h.classes])
    kinds = [('function', 'Functions'), ('typedef', 'Types'),
             ('enum', 'Enumerations'), ('variable', 'Variables'),
             ('define', 'Macros')]
    for kind, name in kinds:
        ms = [m for m in h.members if m['kind'] == kind]
        if kind == 'define':
            # the include guards
            ms = [m for m in ms if not m['name'].startswith('__')]
        if ms:
            out += rubric(name)
            out += ''.join(member_directive(m) + '\n' for m in ms)
    out += rubric('Description of the header')
    out += ('.. doxygenfile:: %s\n   :sections: briefdescription '
            'detaileddescription\n' % h.name)
    return out


def group_page(group):
    return 'group_' + re.sub(r'\W+', '_', group.lower())


def write(page, content):
    full = os.path.join(OUT, page + '.rst')
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, 'w') as f:
        f.write(content)


def doxygen():
    """the XML of the headers, and the number of the warnings of Doxygen"""
    if os.path.isdir(XML):
        shutil.rmtree(XML)
    os.makedirs(BUILD, exist_ok=True)
    try:
        subprocess.run(['doxygen', '-'], input=DOXYFILE, text=True,
                       cwd=ROOT, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    with open(WARNINGS) as f:
        return sum(1 for line in f if line.strip())


def main():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    warnings = doxygen()
    if warnings is None:
        write('index', title('API Reference') + 'Doxygen has not run.\n')
        return None
    compounds, members = read_xml()
    headers = headers_of(compounds, members)
    pages = {}
    for h in headers:
        for c in h.classes:
            pages[c.refid] = 'class_' + re.sub(r'\W+', '_',
                                               c.short).strip('_')

    index = title('API Reference')
    index += ('The library is in three layers, each with the header that '
              'includes it (``core.h``, ``astar.h`` and ``bnb.h``, all of '
              'them ``PBeFS.h``, save the parts that need OpenMP, FastFlow '
              'or MPI): the headers of each layer are gathered by what they '
              'define, and each class has a page of its own.\n\n')
    index += table([(':doc:`%s <%s/index>`' % (layer, layer), line)
                    for layer, line in LAYERS], ('Layer', 'Description'))
    index += toctree(['%s/index' % layer for layer, _ in LAYERS])
    write('index', index)

    for layer, line in LAYERS:
        hs = [h for h in headers if h.layer == layer]
        names = [UMBRELLA] + [g for g, _ in GROUPS[layer]] + [HELPERS]
        groups = collections.OrderedDict(
            (g, sorted((h for h in hs if h.group == g),
                       key=lambda h: (h.order, h.name)))
            for g in names)
        groups = collections.OrderedDict((g, v) for g, v in groups.items()
                                         if v)
        out = title(layer) + line + '\n\n'
        for group, members_of in groups.items():
            out += rubric(':doc:`%s <%s>`' % (group, group_page(group)))
            out += table([(':doc:`%s <%s>`' % (escape(h.path), h.page),
                           escape(summary(h))) for h in members_of],
                         ('Header', 'Description'))
            gp = title(group)
            gp += '*headers of the layer* :doc:`%s <index>`\n\n' % layer
            for h in members_of:
                gp += rubric(':doc:`%s <%s>`' % (escape(h.path), h.page))
                gp += escape(h.brief) + '\n\n'
                gp += table([(':doc:`%s <%s>`' % (escape(c.short),
                                                   pages[c.refid]),
                              escape(capital(c.brief)))
                             for c in h.classes], ('Class', 'Description'))
                write('%s/%s' % (layer, h.page), header_page(h, pages))
                for c in h.classes:
                    write('%s/%s' % (layer, pages[c.refid]), class_page(c, h))
            gp += toctree([h.page for h in members_of])
            write('%s/%s' % (layer, group_page(group)), gp)
        out += toctree([group_page(g) for g in groups])
        write('%s/index' % layer, out)
    return warnings


if __name__ == '__main__':
    n = main()
    print('Doxygen: %s' % ('not run' if n is None else '%d warnings' % n))

# --------------------------------------------------------------------------- #
#                             End of apiref.py                                #
# --------------------------------------------------------------------------- #
