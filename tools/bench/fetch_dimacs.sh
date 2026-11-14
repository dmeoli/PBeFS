#!/bin/bash
###############################################################################
#                                                                             #
#   fetch_dimacs.sh: the road graphs of the 9th DIMACS Challenge              #
#                                                                             #
#   Usage:  fetch_dimacs.sh <dir> [<graph>...]                                #
#                                                                             #
#   Downloads into <dir> the files .gr (distances, USA-road-d) and .co        #
#   (coordinates) of the graphs given (NY, BAY, COL, FLA, NW, NE, CAL, LKS,   #
#   E, W, CTR, USA; NY by default) from the site of the challenge, and        #
#   unpacks them; AStarBench reads them with --domain road (see the README).  #
#   The times (USA-road-t) are the same with "-t" for "-d" in the names.      #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

set -eu

dir=$1
shift
graphs=${*:-NY}
url=http://users.diag.uniroma1.it/challenge9/data/USA-road-d
mkdir -p "$dir"
for g in $graphs; do
    for ext in gr co; do
        f=USA-road-d.$g.$ext
        [ -f "$dir/$f" ] && continue
        curl -fL -o "$dir/$f.gz" "$url/$f.gz"
        gunzip "$dir/$f.gz"
    done
done
