#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
#    tables.py                                                                #
#                                                                             #
#    Writes the tables of the benchmarks of the manual into                   #
#    manual/generated/tables, one Markdown file per table, which a chapter    #
#    includes by the line --8<-- "tables/<name>.md" (see include.lua), from   #
#    the CSV that                                                             #
#    tools/bench/campaign.sh (A.csv to F.csv) and tools/bench/pdb.sh          #
#    (search.csv) write in the directory of the results, whose machine.txt    #
#    describes, in one line, the machine of the runs. The sums and the        #
#    medians are those of tools/bench/tables.py, imported from there. A       #
#    table whose CSV is not there is written as a placeholder that says so,   #
#    so that the manual builds before the runs are done: its numbers are to   #
#    come from an idle machine, and none is made up in the meantime.          #
#                                                                             #
#        python3 manual/build/tables.py [<directory of the results>]          #
#                                                                             #
#    The directory is manual/results by default.                              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

import collections
import os
import statistics
import sys

MANUAL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(MANUAL)
OUT = os.path.join(MANUAL, 'generated', 'tables')
sys.path.insert(0, os.path.join(ROOT, 'tools', 'bench'))
import tables as bench  # noqa: E402  (tools/bench/tables.py)


def tex(s):
    """a string of the CSV as text of a Markdown table"""
    return s.replace('|', '\\|').replace('*', '\\*')


def part(head, rows, note, label=None):
    """the Markdown of one table, with its header, its rows and a note below
    it, after a label in italics if given; an empty table marked as a
    placeholder if there are no rows"""
    out = '*%s*\n\n' % label if label else ''
    out += '| ' + ' | '.join(head) + ' |\n'
    # the dashes of a column one more than the characters of its longest
    # word, which a table wider than the line takes as the width of the
    # column in the PDF, so that no word of a cell overflows it
    body = rows or []
    wide = [max(len(w) for c in [head[i]] + [r[i] for r in body]
                for w in c.split() + ['---']) for i in range(len(head))]
    out += '|' + '|'.join(':' + '-' * w if i == 0 else
                          '-' * w + ':'
                          for i, w in enumerate(wide)) + '|\n'
    for r in rows or [[''] * len(head)]:
        out += '| ' + ' | '.join(r) + ' |\n'
    out += '\n'
    if not rows:
        out += ('*Placeholder: the results of this table are to be measured '
                'on an idle machine and written here by '
                '[`manual/build/tables.py`](https://github.com/dmeoli/PBeFS/'
                'blob/main/manual/build/tables.py); no number is given until '
                'then.*\n')
    elif note:
        out += '*%s*\n' % note.strip()
    return out


def write(name, parts):
    """the file of the table name, made of the Markdown of its parts (the
    caption, with the number of the table, is written by the chapter)"""
    with open(os.path.join(OUT, name + '.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))


def machine(results):
    """the line of machine.txt, or nothing"""
    p = os.path.join(results, 'machine.txt')
    if not os.path.exists(p):
        return ''
    with open(p) as f:
        return 'Machine: %s.' % f.read().strip().rstrip('.')


def noise(gaps):
    """the A/A noise of a list of relative gaps"""
    if not gaps:
        return ''
    return (' A/A noise: median %.1f%%, largest %.1f%% over %d pairs.'
            % (100 * statistics.median(gaps), 100 * max(gaps), len(gaps)))


def scaling(rows, keep, seq, mach, label=None, only=None):
    """the part of a table with the seconds of the searches selected by
    keep on each number of cores, with the speedup over the sequential
    search seq on 1 core, on the instances of the list only (all if None)"""
    rows = [r for r in rows if keep(bench.name(r)) or bench.name(r) == seq]
    if only is not None:
        rows = [r for r in rows if r['instance'] in only]
    tot, cnt, gaps = bench.summarize(rows)
    inst = max(cnt.values()) if cnt else 0
    # the cores of tools/bench/campaign.sh (THREADS) when nothing was run
    cores = sorted({c for (_, c) in tot}) or [1, 2, 4, 8, 16, 32, 48, 64, 96]
    searches = sorted({s for (s, _) in tot}, key=lambda s: (s != seq, s))
    base = tot.get((seq, 1))
    out = []
    for s in searches:
        cells = [tex(s)]
        for c in cores:
            if cnt.get((s, c)) != inst:
                cells.append('')
            elif base:
                cells.append('%.2f (%.1f)' % (tot[(s, c)],
                                              base / tot[(s, c)]))
            else:
                cells.append('%.2f' % tot[(s, c)])
        out.append(cells)
    note = ('Seconds summed over %d instances, each the median of its '
            'repetitions, and in parentheses the speedup over %s on 1 core; '
            'an empty cell is a search that did not solve every instance. '
            '%s%s' % (inst, tex(seq), mach, noise(gaps)))
    return part(['search'] + [str(c) for c in cores], out, note, label)


def sequential(name, rows, mach):
    """the sequential searches on the instances of Korf"""
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r['status'] == 'solved':
            by[bench.name(r)][r['instance']].append(r)
    out, gaps = [], []
    for s, inst in sorted(by.items()):
        sec = sum(statistics.median(float(r['seconds']) for r in v)
                  for v in inst.values())
        exp = sum(int(v[0]['expanded']) for v in inst.values())
        for v in inst.values():
            t = [float(r['seconds']) for r in v]
            if len(t) >= 2 and max(t) > 0:
                gaps.append((max(t) - min(t)) / max(t))
        out.append([tex(s), str(len(inst)), '%.1f' % sec, '%.3g' % exp,
                    '%.2f' % (exp / sec / 1e6 if sec else 0)])
    write(name, [part(['search', 'solved', 'seconds', 'expanded',
                       'M expansions/s'], out,
                      'Each instance the median of its repetitions. %s%s'
                      % (mach, noise(gaps)))])


def heuristics(name, rows, mach):
    """A* with each heuristic of pdb.sh on the instances of Korf"""
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r.get('status') == 'solved':
            k = '%s %s' % (r['heur'], r['enc'])
            by[k][r['id']].append(r)
    out = []
    for k, inst in sorted(by.items()):
        sec = sum(statistics.median(float(r['seconds']) for r in v)
                  for v in inst.values())
        exp = sum(int(v[0]['expanded']) for v in inst.values())
        out.append([tex(k), str(len(inst)), '%.3g' % exp, '%.1f' % sec])
    write(name, [part(['heuristic', 'solved', 'expanded', 'seconds'], out,
                      'Each instance the median of its repetitions. '
                      + mach)])


def examples(name, rows, mach):
    """the problems of the examples"""
    by, order = collections.defaultdict(list), []
    for r in rows:
        k = (r['example'], r['arguments'])
        if k not in by:
            order.append(k)
        by[k].append(r)
    out = []
    for k in order:
        v = [r for r in by[k] if r['seconds']]
        if v:
            ms = statistics.median(float(r['seconds']) for r in v) * 1000
            out.append([tex(k[0]), tex(k[1] or '-'), v[0]['cost'] or '-',
                        v[0]['expanded'], '%.3g' % ms])
    write(name, [part(['example', 'arguments', 'cost', 'expanded',
                       'milliseconds'], out,
                      'Each problem the median of its repetitions. '
                      + mach)])


def main(results):
    os.makedirs(OUT, exist_ok=True)
    mach = machine(results)

    def csv(f):
        return bench.load(os.path.join(results, f))

    b, c, d = csv('B.csv'), csv('C.csv'), csv('D.csv')
    f = [r for r in csv('F.csv') if r['dist'] in bench.SELECTIONS]
    sequential('best-first-korf', csv('A.csv'), mach)
    heuristics('heuristics-pdb', csv(os.path.join('pdb', 'search.csv')),
               mach)
    # campaign B runs SPA* and MQA* on 6 of the 12 instances of the other
    # searches: the table is on those 6
    six = {r['instance'] for r in b if r['algo'] in ('spa', 'mqa')}
    write('shared-memory-tiles', [scaling(
        b, lambda s: s.split()[0] in ('SPA*', 'MQA*', 'farm'), 'A* buckets',
        mach, only=six or None)])
    # campaign D: one part for each work per expansion
    slows = sorted({r['slow'] for r in d}, key=int) or ['1000', '10000']
    write('shared-memory-slow', [scaling(
        [r for r in d if r['slow'] == w],
        lambda s: s.split()[0] in ('SPA*', 'MQA*', 'farm', 'PA*SE', 'HDA*'),
        'A* buckets', mach, label='work %s per expansion' % w)
        for w in slows])
    write('hash-distributed-tiles', [scaling(
        b, lambda s: s.startswith('HDA*'), 'A* buckets', mach)])
    write('safepbnf-grids', [scaling(
        [r for r in c if r['domain'] == 'grid8'],
        lambda s: s.startswith(('Safe PBNF', 'HDA*')), 'A* buckets', mach)])
    # table F: one part for each problem of branch-and-bound
    probs = sorted({r['problem'] for r in f}) or ['knapsack', 'tsp-bb']
    write('branch-and-bound-schemes', [scaling(
        [r for r in f if r['problem'] == q], lambda s: s.startswith('B&B'),
        'B&B sequential', mach, label=q) for q in probs])
    examples('domains-examples', csv('E.csv'), mach)
    # the tables whose runs tools/bench/campaign.sh does not make yet:
    # placeholders, with the columns that the chapters describe
    for name, head in PLANNED:
        write(name, [part(head, [], '')])


# the tables of the chapters whose campaign is still to be written, by name,
# with their columns (the cells of the scaling tables are the seconds and,
# in parentheses, the speedup)
CORES = ['1', '2', '4', '8', '16', '32', '48', '64', '96']
PLANNED = [
    ('best-first-greedy', ['search', 'solved', 'cost', 'cost / optimal',
                           'expanded', 'evaluated', 'seconds']),
    ('memory-bounded-classic', ['search', 'instances', 'solved', 'seconds',
                                'expanded', 'stored',
                                'expansions per second (millions)']),
    ('bounded-suboptimal-sequential', ['search', 'w', 'expanded',
                                       'cost / C*', 'cost / bound',
                                       'seconds', 'seconds to a path']),
    ('bounded-suboptimal-incremental', ['search', 'workers', 'policy',
                                        'expanded', 'under', 'messages',
                                        'restarts', 'seconds']),
    ('bidirectional-memory', ['instances', 'search', 'expanded forward',
                              'expanded backward', 'stored', 'seconds']),
    ('bidirectional-external', ['search', 'workers', 'expanded', 'stored',
                                'bytes written', 'bytes on disk at most',
                                'seconds']),
    ('heuristics-permutation', ['puzzle', 'search', 'heuristic', 'instances',
                                'entries', 'MB', 'seconds to build',
                                'solved', 'expanded', 'seconds']),
    ('shared-memory-edges', ['search'] + CORES),
    ('auto-distributions-tiles', ['distribution', 'estCO', 'estLB',
                                  'effesti', 'CO', 'LB', 'LBstored', 'SO',
                                  'seconds (speedup), 8 cores',
                                  'seconds (speedup), 32 cores']),
    ('auto-distributions-planning', ['distribution', 'estCO', 'estLB',
                                     'effesti', 'CO', 'LB', 'LBstored', 'SO',
                                     'seconds (speedup), 8 cores',
                                     'seconds (speedup), 32 cores']),
    ('more-parallel-schemes-tiles', ['search'] + CORES),
    ('more-parallel-schemes-greedy', ['search'] + CORES),
    ('branch-and-bound-domains', ['search'] + CORES),
    ('multi-objective-roads', ['search', 'workers', 'objectives', 'solved',
                               'front', 'expanded', 'most by a worker',
                               'checks', 'seconds']),
    ('msa-instances-hwd', ['instance', 'distribution', 'workers', 'seconds',
                           'speedup', 'SO', 'CO', 'LB']),
    ('planning-heuristics', ['task', 'C*', 'blind', 'hmax', 'lmcut', 'pdb',
                             'mas']),
    ('planning-parallel', ['task', 'search', 'seconds', 'speedup', 'SO',
                           'CO', 'LB']),
]


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1
         else os.path.join(MANUAL, 'results'))

# --------------------------------------------------------------------------- #
# ------------------------------ END tables.py ------------------------------ #
# --------------------------------------------------------------------------- #
