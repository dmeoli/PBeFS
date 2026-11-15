#!/bin/bash
###############################################################################
#                                                                             #
#   fetch_dimacs.sh: the road graphs of the 9th DIMACS Challenge              #
#                                                                             #
#   Usage:  fetch_dimacs.sh <dir> [<graph>...] [--times]                      #
#                                                                             #
#   Downloads into <dir> the files .gr (distances, USA-road-d) and .co        #
#   (coordinates) of the graphs given (NY, BAY, COL, FLA, NW, NE, CAL, LKS,   #
#   E, W, CTR, USA; NY by default) from the site of the challenge, and        #
#   unpacks them; AStarBench reads them with --domain road (see the README).  #
#   With --times it downloads also the files .gr of the times (USA-road-t),   #
#   the second objective of MOABench.                                         #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

set -eu

dir=$1
shift
times=0
graphs=
for a in "$@"; do
    if [ "$a" = --times ]; then times=1; else graphs="$graphs $a"; fi
done
graphs=${graphs:-NY}
url=http://users.diag.uniroma1.it/challenge9/data
mkdir -p "$dir"
fetch() {   # fetch <kind d or t> <file>
    [ -f "$dir/$2" ] && return
    curl -fL -o "$dir/$2.gz" "$url/USA-road-$1/$2.gz"
    gunzip "$dir/$2.gz"
}
for g in $graphs; do
    fetch d "USA-road-d.$g.gr"
    fetch d "USA-road-d.$g.co"
    [ $times = 1 ] && fetch t "USA-road-t.$g.gr"
done
exit 0
