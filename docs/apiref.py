#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    apiref.py                                                                #
#                                                                             #
#    Writes the API reference of the site into docs/api, from the XML that    #
#    docs/doxygen/Doxyfile has Doxygen extract from the headers into          #
#    build-docs/xml, as a tree that mirrors the library: a page with PBeFS.h  #
#    and the three layers (core, astar and bnb), a page per layer, in it a    #
#    page per group of headers (the searches, the domains, the runtimes,      #
#    ...), a page per header and a page per class, each rendered by Breathe   #
#    and linked to its page on the site of Doxygen. Doxygen writes its        #
#    warnings into build-docs/doxygen-warnings.log, which docs.yml wants      #
#    empty. pages.py calls it, before Sphinx.                                 #
#                                                                             #
#    The same tree goes to the Topics of the site of Doxygen: apiref.py       #
#    writes build-docs/modules.dox, a group per layer, in it a group per      #
#    group of headers and in that a group per header, which holds the header  #
#    and its classes and concepts; the file has only external documentation   #
#    blocks, so no header changes, and the second run of Doxygen in docs.yml  #
#    reads it.                                                                #
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
DOXYFILE = os.path.join(ROOT, 'docs', 'doxygen', 'Doxyfile')
BUILD = os.path.join(ROOT, 'build-docs')
XML = os.path.join(BUILD, 'xml')
WARNINGS = os.path.join(BUILD, 'doxygen-warnings.log')
TOPICS = os.path.join(BUILD, 'modules.dox')
PAGES = os.path.join(BUILD, 'pages.cfg')
OUT = os.path.join(ROOT, 'docs', 'api')
NAMESPACE = 'PBeFS_di_unipi_it::'
BLOB = 'https://github.com/dmeoli/PBeFS/blob/main/include/PBeFS/'

# where the pages of the site of Doxygen are, seen from a page of a layer,
# which lies in api/<layer>: the site of Doxygen is served in the
# directory doxygen of the site of Sphinx
DOXYGEN = '../../doxygen/'

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
                     'ExplicitGraph.h', 'SlowDomain.h', 'SequenceScore.h',
                     'Romania.h', 'NQueens.h', 'VacuumWorld.h',
                     'MissionariesCannibals.h', 'WaterJugs.h',
                     'TravelingSalesman.h']),
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
        ('Problems', ['Knapsack.h']),
    ],
}
LAYER_NAMES = [layer for layer, _ in LAYERS]
HELPERS = 'Helpers'

# the pages of a layer, from the Markdown files in docs, each with its title,
# its page on the site of Sphinx [see pages.py] and its line
DOCUMENTATION = 'Documentation'
DOCS = {
    'astar': [('theory/best-first-search.md', 'Best-first search and A*',
               '/generated/theory/best-first-search',
               'The search problem, A* with its invariants and guarantees, '
               'and the cost of its operations on the data structures of '
               'the library.'),
              ('theory/parallel-astar.md', 'Parallel A*',
               '/generated/theory/parallel-astar',
               'The overheads of a parallel best-first search, and what '
               'SPA*, HDA* with its distributions, SafePBNF, PA*SE and the '
               'farm guarantee and pay.')],
    'bnb': [('theory/branch-and-bound-search.md',
             'Branch-and-bound as best-first search',
             '/generated/theory/branch-and-bound-search',
             'The incumbent, the bounds, the pruning and the selection '
             'rules, and A* as a special case.'),
            ('theory/parallel-branch-and-bound.md',
             'Parallel branch-and-bound',
             '/generated/theory/parallel-branch-and-bound',
             'The shared, distributed, master-worker, deterministic, '
             'SelfSplit and stealing schemes: correctness, termination and '
             'determinism.'),
            ('branch-and-bound.md', 'Best-first search with an incumbent: '
             'design note', '/generated/user-guide/branch-and-bound-design',
             'The extension of the searches to branch-and-bound: the '
             'concept, the sequential and the parallel searches, the '
             'reoptimization.')],
}
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
        self.bases = [b.get('refid') for b in cd.findall('basecompoundref')
                      if b.get('refid')]
        self.derived = [b.get('refid')
                        for b in cd.findall('derivedcompoundref')
                        if b.get('refid')]


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

def html_page(refid):
    """the page of a compound on the site of Doxygen"""
    return refid + '.html'


class Header:
    def __init__(self, f, path, classes, concepts):
        self.path = path                     # relative to include/PBeFS
        self.name = os.path.basename(path)
        self.stem = os.path.splitext(self.name)[0]
        self.refid = f.refid
        self.brief = f.brief
        self.classes = classes
        self.concepts = concepts
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
    concepts = collections.defaultdict(list)
    for c in compounds.values():
        if c.kind == 'concept':
            concepts[c.file].append(c)
    headers = []
    for f in compounds.values():
        if f.kind != 'file' or not f.file.endswith('.h'):
            continue
        path = os.path.relpath(f.file, INCLUDE) if os.path.isabs(f.file) \
            else re.sub(r'^PBeFS/', '', f.file)
        h = Header(f, path, sorted(by_file.get(f.file, []),
                                   key=lambda c: (c.short != f.name[:-2],
                                                  c.line)),
                   sorted(concepts.get(f.file, []), key=lambda c: c.line))
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
    out += ('Declared in :doc:`%s <%s>`; the same %s on the `site of Doxygen '
            '<%s%s>`__.\n\n' % (h.path, h.page, kind, DOXYGEN,
                                  html_page(c.refid)))
    if c.bases or c.derived:
        # the diagram of the hierarchy that the site of Doxygen draws for
        # the class, which leads to its page there
        out += ('.. raw:: html\n\n   <p><a href="%s%s"><img alt="inheritance'
                ' diagram of %s" src="%s%s"></a></p>\n\n'
                % (DOXYGEN, html_page(c.refid), c.short, DOXYGEN,
                   c.refid + '.png'))
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


def header_page(h, pages, doxygen=DOXYGEN):
    out = title(escape(h.path))
    if h.layer in LAYER_NAMES:
        out += ('Header of the layer :doc:`%s <index>`, in the group '
                ':doc:`%s <%s>`; ' % (h.layer, h.group, group_page(h.group)))
    else:
        out += 'The header of the whole library; '
    out += ('the same header on the `site of Doxygen <%s%s>`__, its source '
            'on `GitHub <%s%s>`__.\n\n'
            % (doxygen, html_page(h.refid), BLOB, h.path))
    if h.classes:
        out += rubric('Classes')
        out += table([(':doc:`%s <%s>`' % (escape(c.short), pages[c.refid]),
                       escape(capital(c.brief))) for c in h.classes],
                     ('Class', 'Description'))
        out += toctree([pages[c.refid] for c in h.classes])
    if h.concepts:
        out += rubric('Concepts')
        out += ''.join('.. doxygenconcept:: %s\n\n' % c.name
                       for c in h.concepts)
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
    with open(DOXYFILE) as f:
        config = f.read() + '\nGENERATE_HTML = NO\n'
    try:
        subprocess.run(['doxygen', '-'], input=config, text=True,
                       cwd=ROOT, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    with open(WARNINGS) as f:
        return sum(1 for line in f if line.strip())


def tree(headers):
    """the headers of each layer by group, the groups in the order of the
    pages and the empty ones left out"""
    layers = collections.OrderedDict()
    for layer, _ in LAYERS:
        hs = [h for h in headers if h.layer == layer]
        names = [UMBRELLA] + [g for g, _ in GROUPS[layer]] + [HELPERS]
        groups = collections.OrderedDict(
            (g, sorted((h for h in hs if h.group == g),
                       key=lambda h: (h.order, h.name)))
            for g in names)
        layers[layer] = collections.OrderedDict(
            (g, v) for g, v in groups.items() if v)
    return layers


# ----- the Topics of the site of Doxygen ----------------------------------- #

def dox(line):
    """a line of text in a comment of Doxygen, its commands escaped"""
    return re.sub(r'([\\@&$#<>%])', r'\\\1', line)


def ident(*words):
    """the name of a group of Doxygen, which takes only letters, digits
    and underscores"""
    return '_'.join(re.sub(r'\W+', '_', w).strip('_') for w in words)


def block(*lines):
    """a documentation block of Doxygen with the given lines"""
    return '/** %s\n%s */\n' % (lines[0], ''.join(
        (' *  ' + line).rstrip() + '\n' for line in lines[1:]))


def header_topic(h, parent):
    """the group of a header, which holds the header, its classes and its
    concepts"""
    hdr = ident(h.page)
    lines = ['\\defgroup %s %s' % (hdr, h.path)]
    if parent:
        lines.append('\\ingroup %s' % parent)
    if summary(h):
        lines.append('\\brief %s' % dox(summary(h)))
    out = [block(*lines),
           block('\\file %s' % os.path.join('PBeFS', h.path),
                 '\\ingroup %s' % hdr)]
    for c in h.classes:
        out.append(block('\\%s %s' % ('struct' if c.kind == 'struct'
                                        else 'class', c.name),
                         '\\ingroup %s' % hdr))
    for c in h.concepts:
        out.append(block('\\concept %s' % c.name, '\\ingroup %s' % hdr))
    return out


def doc_groups():
    """the group of the site of Doxygen of each page of a layer, by the
    path of its Markdown file [see DOCS]"""
    return {'docs/' + d[0]: ident('layer', layer, DOCUMENTATION.lower(),
                                  os.path.splitext(d[0])[0])
            for layer, ds in DOCS.items() for d in ds}


def topics(library, layers):
    """the external documentation blocks that put the tree of the pages in
    the Topics of the site of Doxygen: the group of PBeFS.h, then a group
    per layer, a group per group of headers in it, a group per header in
    that"""
    out = ['/* written by docs/apiref.py from the XML of Doxygen */\n']
    filters = []
    for h in library:
        out += header_topic(h, None)
    for layer, line in LAYERS:
        lay = ident('layer', layer)
        out.append(block('\\defgroup %s %s' % (lay, layer),
                         '\\brief %s' % dox(line)))
        for group, hs in layers[layer].items():
            grp = ident('layer', layer, group.lower())
            out.append(block('\\defgroup %s %s' % (grp, group),
                             '\\ingroup %s' % lay))
            for h in hs:
                out += header_topic(h, grp)
        if DOCS.get(layer):
            doc = ident('layer', layer, DOCUMENTATION.lower())
            out.append(block('\\defgroup %s %s' % (doc, DOCUMENTATION),
                             '\\ingroup %s' % lay,
                             '\\brief The pages on the theory and on the '
                             'design of the layer.'))
            # a group per page, titled as the page, whose documentation is
            # the page itself, by an \addtogroup that docs/doxygen/markdown.py
            # puts in place of its title as the input filter; the group has
            # no brief, which Doxygen would print twice
            for d in DOCS[layer]:
                page = doc_groups()['docs/' + d[0]]
                out.append(block('\\defgroup %s %s' % (page, dox(d[1])),
                                 '\\ingroup %s' % doc))
                filters.append('"*/%s=python3 docs/doxygen/markdown.py '
                               '--topics %s"' % (d[0], page))
    with open(TOPICS, 'w') as f:
        f.write('\n'.join(out))
    # the filter of the other Markdown files comes last, the first pattern
    # that matches a file being the one used
    filters.append('"*.md=python3 docs/doxygen/markdown.py --topics"')
    # the mentions of the pages in the comments of the headers go to their
    # groups, which only this run has [see docs/doxygen/markdown.py]
    filters.append('"*.h=python3 docs/doxygen/markdown.py --topics"')
    with open(PAGES, 'w') as f:
        f.write('# written by docs/apiref.py\nFILTER_PATTERNS = %s\n'
                % ' \\\n                  '.join(filters))


def main():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    warnings = doxygen()
    if warnings is None:
        write('index', title('API reference') + 'Doxygen has not run.\n')
        return None
    compounds, members = read_xml()
    headers = headers_of(compounds, members)
    pages = {}
    for h in headers:
        for c in h.classes:
            pages[c.refid] = 'class_' + re.sub(r'\W+', '_',
                                               c.short).strip('_')
    library = [h for h in headers if h.layer not in LAYER_NAMES]
    layers = tree(headers)

    index = title('API reference')
    index += ('The library is in three layers, each with the header that '
              'includes it (``core.h``, ``astar.h`` and ``bnb.h``, all of '
              'them ``PBeFS.h``, save the parts that need OpenMP, FastFlow '
              'or MPI): the headers of each layer are gathered by what they '
              'define, and each class has a page of its own. The same '
              'reference, with the same tree under its Topics and with the '
              'lists of all the classes, of all the concepts and of all the '
              'files, is on the `site of Doxygen <%s>`__.\n\n'
              % DOXYGEN[3:])
    index += rubric('The library')
    index += table([(':doc:`%s <%s>`' % (escape(h.path), h.page),
                     escape(summary(h))) for h in library],
                   ('Header', 'Description'))
    index += rubric('The layers')
    index += table([(':doc:`%s <%s/index>`' % (layer, layer), line)
                    for layer, line in LAYERS], ('Layer', 'Description'))
    index += toctree([h.page for h in library] +
                     ['%s/index' % layer for layer, _ in LAYERS])
    write('index', index)
    for h in library:
        write(h.page, header_page(h, pages, DOXYGEN[3:]))
    topics(library, layers)

    for layer, line in LAYERS:
        groups = layers[layer]
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
        if DOCS.get(layer):
            out += rubric(DOCUMENTATION)
            out += table([(':doc:`%s <%s>`' % (escape(t), doc), line)
                          for _, t, doc, line in DOCS[layer]],
                         ('Page', 'Description'))
        out += toctree([group_page(g) for g in groups])
        write('%s/index' % layer, out)
    return warnings


if __name__ == '__main__':
    n = main()
    print('Doxygen: %s' % ('not run' if n is None else '%d warnings' % n))

# --------------------------------------------------------------------------- #
#                             End of apiref.py                                #
# --------------------------------------------------------------------------- #
