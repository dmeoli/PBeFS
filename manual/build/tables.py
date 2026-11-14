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


def write(name, caption, head, rows, note):
    """the table name, with its header, its rows and a note below it; an
    empty table marked as a placeholder if there are no rows (the caption,
    with the number of the table, is written by the chapter)"""
    del caption
    path = os.path.join(OUT, name + '.md')
    with open(path, 'w', encoding='utf-8') as f:
        f.write('| ' + ' | '.join(head) + ' |\n')
        f.write('|' + '|'.join(':--' if i == 0 else '--:'
                               for i in range(len(head))) + '|\n')
        for r in rows or [[''] * len(head)]:
            f.write('| ' + ' | '.join(r) + ' |\n')
        f.write('\n')
        if not rows:
            f.write('*Placeholder: the results of this table are to be '
                    'measured on an idle machine and written here by '
                    '`manual/build/tables.py`; no number is given until '
                    'then.*\n')
        elif note:
            f.write('*%s*\n' % note.strip())


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


def scaling(name, caption, rows, keep, seq, mach):
    """the seconds of the searches selected by keep on each number of
    cores, with the speedup over the sequential search seq on 1 core"""
    rows = [r for r in rows if keep(bench.name(r)) or bench.name(r) == seq]
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
    write(name, caption, ['search'] + [str(c) for c in cores], out, note)


def sequential(name, caption, rows, mach):
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
    write(name, caption, ['search', 'solved', 'seconds', 'expanded',
                          'M expansions/s'], out,
          'Each instance the median of its repetitions. %s%s'
          % (mach, noise(gaps)))


def heuristics(name, caption, rows, mach):
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
    write(name, caption, ['heuristic', 'solved', 'expanded', 'seconds'],
          out, 'Each instance the median of its repetitions. ' + mach)


def examples(name, caption, rows, mach):
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
    write(name, caption, ['example', 'arguments', 'cost', 'expanded',
                          'milliseconds'], out,
          'Each problem the median of its repetitions. ' + mach)


def main(results):
    os.makedirs(OUT, exist_ok=True)
    mach = machine(results)

    def csv(f):
        return bench.load(os.path.join(results, f))

    b, c, d = csv('B.csv'), csv('C.csv'), csv('D.csv')
    f = [r for r in csv('F.csv') if r['dist'] in bench.SELECTIONS]
    sequential('best-first-korf', 'Sequential A$^*$ on the 100 instances '
               'of Korf.', csv('A.csv'), mach)
    heuristics('heuristics-pdb', 'A$^*$ with the Manhattan distance and '
               'with the pattern databases on the instances of Korf.',
               csv(os.path.join('pdb', 'search.csv')), mach)
    scaling('shared-memory-tiles', 'The searches on a shared OPEN on the '
            '15-puzzle.', b,
            lambda s: s.split()[0] in ('SPA*', 'MQA*', 'farm', 'PA*SE'),
            'A* buckets', mach)
    scaling('shared-memory-slow', 'PA$^*$SE and HDA$^*$ with slow '
            'expansions on the 15-puzzle.', d,
            lambda s: s.startswith(('PA*SE', 'HDA*')), 'A* buckets', mach)
    scaling('hash-distributed-tiles', 'HDA$^*$ with each distribution on '
            'the 15-puzzle.', b, lambda s: s.startswith('HDA*'),
            'A* buckets', mach)
    scaling('safepbnf-grids', 'Safe PBNF and HDA$^*$ on random grids with '
            '8 moves.', [r for r in c if r['domain'] == 'grid8'],
            lambda s: s.startswith(('Safe PBNF', 'HDA*')), 'A* buckets',
            mach)
    scaling('branch-and-bound-schemes', 'The parallel schemes of '
            'branch-and-bound on the 0-1 knapsack.', f,
            lambda s: s.startswith('B&B'), 'B&B sequential', mach)
    examples('domains-examples', 'The problems of the examples, solved by '
             'sequential A$^*$.', csv('E.csv'), mach)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1
         else os.path.join(MANUAL, 'results'))

# --------------------------------------------------------------------------- #
# ------------------------------ END tables.py ------------------------------ #
# --------------------------------------------------------------------------- #
