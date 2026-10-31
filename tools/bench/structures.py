#!/usr/bin/env python3
###############################################################################
#                                                                             #
#   structures.py: the tables of A* on the sets of data structures            #
#                                                                             #
#   Usage:  structures.py <out of structures.sh> [--html]                     #
#                                                                             #
#   For each problem, size and search it gives (sequential A*, then each      #
#   parallel search by its number of threads), set by set (ADS, std and       #
#   Boost, one under the other), the total seconds over the instances (the    #
#   mean of the two repetitions of each), the nodes expanded per second, and  #
#   the time of the set divided by that of ADS, as the median and the range   #
#   over the instances (above 1 when ADS is faster). It also checks that      #
#   every set finds the same cost on every instance, and gives the A/A noise, #
#   i.e., the gap between the two repetitions of the same run, each in a      #
#   process of its own: on each row the median over its instances, and at the #
#   end the median and the 90th percentile over all the pairs. With --html    #
#   the table is HTML, a problem and a size spanning their rows, as in the    #
#   README; otherwise it is Markdown.                                         #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

import collections
import csv
import statistics
import sys

SETS = [('ads', 'ADS'), ('std', 'std'), ('boost', 'Boost')]
PARALLEL = {'hda': 'HDA*', 'spa': 'SPA*', 'pbnf': 'SafePBNF', 'pase': 'PA*SE',
            'farm': 'farm'}


def problem(r):
    """the problem of a row, and the key that orders the sizes"""
    if r['domain'] == 'tiles15':
        return '15-puzzle', 'about %s nodes' % r['size'], float(r['size'])
    moves = r['domain'][-1]
    side = int(r['size'])
    return ('grid, %s moves' % moves, '%d x %d' % (side, side), side)


def search(r):
    """the key that orders the search of a row, and its name"""
    if r['algo'] in PARALLEL:
        t = int(r['threads'])
        return ((1, r['algo'], t),
                '%s %d thread%s' % (PARALLEL[r['algo']], t,
                                    '' if t == 1 else 's'))
    if r['domain'] == 'tiles15':
        return (0, r['open'], 1), 'A* ' + r['open']
    return (0, '', 1), 'A*'


def fmt(x):
    """a number with three significant digits"""
    if x >= 100:
        return '%.0f' % x
    if x >= 10:
        return '%.1f' % x
    return '%.2f' % x


def main(out, html):
    with open(out + '/structures.csv') as f:
        rows = [r for r in csv.DictReader(f) if r['status'] != 'node_limit']
    secs = collections.defaultdict(list)    # run -> seconds of the reps
    nodes = {}                              # run -> expanded
    cost = collections.defaultdict(set)     # instance -> costs of the sets
    for r in rows:
        p, size, order = problem(r)
        inst = (p, size, r['domain'], r['instance'])
        run = (p, size, order, search(r), r['set'], r['instance'])
        secs[run].append(float(r['seconds']))
        nodes[run] = int(r['expanded'])
        cost[inst].add(r['cost'])
    wrong = [i for i, c in cost.items() if len(c) > 1]

    gaps = [abs(s[0] - s[1]) / statistics.mean(s) for s in secs.values()
            if len(s) == 2 and statistics.mean(s) > 0]
    groups = collections.defaultdict(dict)   # ( p , size , search ) -> set
    for (p, size, order, se, s, inst), t in secs.items():
        groups[(p, order, size, se[0], se[1])].setdefault(s, {})[inst] = (
            statistics.mean(t), nodes[(p, size, order, se, s, inst)],
            abs(t[0] - t[1]) / statistics.mean(t)
            if len(t) == 2 and statistics.mean(t) > 0 else None)

    lines = []
    for (p, order, size, key, se) in sorted(groups):
        g = groups[(p, order, size, key, se)]
        if 'ads' not in g:
            continue
        cells = []
        for s, label in SETS:
            if s not in g:
                continue
            tot = sum(t for t, n, a in g[s].values())
            nps = sum(n for t, n, a in g[s].values()) / tot if tot else 0
            ratios = [g[s][i][0] / g['ads'][i][0] for i in g[s]
                      if i in g['ads'] and g['ads'][i][0] > 0]
            aa = [a for t, n, a in g[s].values() if a is not None]
            noise = ('%.1f%%' % (100 * statistics.median(aa)) if aa
                     else '')
            cells.append((label, len(g[s]), tot, nps, ratios, noise))
        lines.append((p, size, se, cells))

    if html:
        print('<table>')
        print('<tr><th>problem</th><th>size</th><th>search</th><th>set</th>'
              '<th>instances</th><th>seconds</th><th>nodes/s</th>'
              '<th>/ ADS median</th><th>range</th><th>A/A</th></tr>')
        for p, size, se, cells in lines:
            first = True
            for label, n, tot, nps, ratios, noise in cells:
                head = ''
                if first:
                    k = len(cells)
                    head = ('<td rowspan="%d">%s</td><td rowspan="%d">%s</td>'
                            '<td rowspan="%d">%s</td>' % (k, p, k, size, k,
                                                          se))
                    first = False
                rng = ('%s to %s' % (fmt(min(ratios)), fmt(max(ratios)))
                       if label != 'ADS' else '')
                med = fmt(statistics.median(ratios)) if label != 'ADS' else ''
                print('<tr>%s<td>%s</td><td align="right">%d</td>'
                      '<td align="right">%s</td><td align="right">%.2e</td>'
                      '<td align="right">%s</td><td align="right">%s</td>'
                      '<td align="right">%s</td></tr>'
                      % (head, label, n, fmt(tot), nps, med, rng, noise))
        print('</table>')
    else:
        print('| problem | size | search | set | instances | seconds | '
              'nodes/s | / ADS median | range | A/A |')
        print('|---|---|---|---|--:|--:|--:|--:|--:|--:|')
        for p, size, se, cells in lines:
            for label, n, tot, nps, ratios, noise in cells:
                med = fmt(statistics.median(ratios)) if label != 'ADS' else ''
                rng = ('%s to %s' % (fmt(min(ratios)), fmt(max(ratios)))
                       if label != 'ADS' else '')
                print('| %s | %s | %s | %s | %d | %s | %.2e | %s | %s | %s |'
                      % (p, size, se, label, n, fmt(tot), nps, med, rng,
                         noise))
    print()
    if gaps:
        gaps.sort()
        print('A/A: %d pairs of runs, median gap %.1f%%, 90th percentile '
              '%.1f%%' % (len(gaps), 100 * statistics.median(gaps),
                          100 * gaps[int(0.9 * (len(gaps) - 1))]))
    print('costs: %s' % ('the same with every set' if not wrong else
                         'DIFFERENT on %s' % wrong))


if __name__ == '__main__':
    main(sys.argv[1], '--html' in sys.argv[2:])
