#!/usr/bin/env python3
###############################################################################
#                                                                             #
#   dimacs_sample.py: the sample road graph of data/dimacs                    #
#                                                                             #
#   Usage:  dimacs_sample.py <dir>                                            #
#                                                                             #
#   Writes <dir>/sample.gr and <dir>/sample.co, a synthetic graph in the      #
#   format of the 9th DIMACS Implementation Challenge: a 30 x 30 grid of      #
#   nodes 900 millionths of a degree apart, near 10.4 E 43.72 N, moved by up  #
#   to 250 at random, each edge of the grid kept with probability 0.9 as two  #
#   arcs whose length is the great-circle distance in meters times a factor   #
#   from 1 to 1.5, rounded up; the seed is fixed, and the files are those of  #
#   the data.                                                                 #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

import math
import random
import sys

random.seed(2006)
W = H = 30
LON0, LAT0, STEP = 10400000, 43720000, 900
co = {}
for y in range(H):
    for x in range(W):
        co[y * W + x + 1] = (LON0 + x * STEP + random.randint(-250, 250),
                             LAT0 + y * STEP + random.randint(-250, 250))


def dist(a, b):
    r = math.pi / 180e6
    p1, p2 = co[a][1] * r, co[b][1] * r
    dl = (co[b][0] - co[a][0]) * r
    h = math.sin((p2 - p1) / 2) ** 2 + \
        math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371008.8 * math.asin(min(1, math.sqrt(h)))


arcs = []
for y in range(H):
    for x in range(W):
        u = y * W + x + 1
        for dx, dy in ((1, 0), (0, 1)):
            if x + dx < W and y + dy < H and random.random() < 0.9:
                v = (y + dy) * W + x + dx + 1
                w = math.ceil(dist(u, v) * random.uniform(1.0, 1.5))
                arcs += [(u, v, w), (v, u, w)]
d = sys.argv[1]
with open(d + '/sample.gr', 'w') as f:
    f.write('c a synthetic road graph in the format of the 9th DIMACS '
            'Implementation Challenge\nc (see NOTICE): a 30 x 30 perturbed '
            'grid, lengths in meters\n')
    f.write('p sp %d %d\n' % (W * H, len(arcs)))
    for a in arcs:
        f.write('a %d %d %d\n' % a)
with open(d + '/sample.co', 'w') as f:
    f.write('c the coordinates of sample.gr, in millionths of a degree\n')
    f.write('p aux sp co %d\n' % (W * H))
    for v in range(1, W * H + 1):
        f.write('v %d %d %d\n' % (v, co[v][0], co[v][1]))
