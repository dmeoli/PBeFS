#!/bin/bash
###############################################################################
#                                                                             #
#   campaign.sh: the runs behind the tables of the README                     #
#                                                                             #
#   Usage:  campaign.sh <build> <out> [<first core> <NUMA node>] [<tables>]   #
#                                                                             #
#   Runs tools/AStarBench (and tools/AStarBenchMPI) of the build directory    #
#   <build> on the instances of data/, and writes one CSV per table in        #
#   <out>: A (sequential A* on the 100 instances of Korf, and the solver of   #
#   Burns et al. if BURNS_SOLVER names it), B (the parallel searches on 12    #
#   of them), C (random grids), D (slow expansions) and E (the examples).     #
#   Each line of a CSV starts with the table, the cores given to the run      #
#   and the artificial work per expansion, then the line of AStarBench.       #
#   The runs of a table are pinned to the cores <first core> , ... of the     #
#   NUMA node <NUMA node> (default 96 and 1), with the memory on that node    #
#   and the workers one per core (PBEFS_PIN). Each run is a process of its    #
#   own, and each parallel one (and sequential A* with buckets) is run twice, #
#   in two processes, so that the gap between the two is the A/A noise of     #
#   the machine. The campaign holds a lock (flock) on each of its cores,      #
#   in LOCKS (default /tmp/pbefs-bench-locks), so that a second campaign on   #
#   the same cores waits for it. BURNS_SOLVER is to be built with the same    #
#   flags as AStarBench [see tools/CMakeLists.txt], so that the comparison    #
#   does not depend on the alignment of the code. Linux only: numactl places  #
#   the runs; launch.sh starts it detached, with its PID saved.               #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

set -u
BUILD=$(realpath "$1")
OUT=$(realpath -m "$2")
C0=${3:-96}
NODE=${4:-1}
TABLES=${5:-ABCDE}
SRC=$(realpath "$(dirname "$0")/../..")
B=$BUILD/tools/AStarBench
M=$BUILD/tools/AStarBenchMPI
K=$SRC/data/15-puzzle/korf100.txt
THREADS=${THREADS:-"1 2 4 8 16 32 48 64 96"}
LOCKS=${LOCKS:-/tmp/pbefs-bench-locks}
mkdir -p "$OUT" "$LOCKS"
export PBEFS_PIN=1

# a lock on each core of the campaign, held until it ends
max=1
for t in $THREADS; do [ "$t" -gt "$max" ] && max=$t; done
for c in $(seq "$C0" $(( C0 + max - 1 ))); do
  exec {fd}> "$LOCKS/core$c"
  flock -n "$fd" || { echo "core $c is busy, waiting" >&2; flock "$fd"; }
done

# runs AStarBench on t cores, the line prefixed by table, cores and slow
run() { # table cores slow options...
  local tab=$1 t=$2 slow=$3
  shift 3
  numactl --physcpubind=$C0-$(( C0 + t - 1 )) --membind=$NODE \
          "$B" --slow "$slow" "$@" 2>> "$OUT/errors.log" |
    sed "s/^/$tab,$t,$slow,/" >> "$OUT/$tab.csv"
  }

# the same with MPI, one process per core
run_mpi() { # table cores slow options...
  local tab=$1 t=$2 slow=$3
  shift 3
  mpirun -np "$t" --bind-to core --cpu-set $C0-$(( C0 + t - 1 )) \
         -x PBEFS_PIN numactl --membind=$NODE "$M" --slow "$slow" "$@" \
         2>> "$OUT/errors.log" | sed "s/^/$tab,$t,$slow,/" >> "$OUT/$tab.csv"
  }

# a run twice, in two processes, for the A/A noise
twice() { # command...
  "$@"
  "$@"
  }

# the parallel searches of a domain on t cores; the farm runs its master
# and t - 1 workers
parallel() { # table t slow options...
  local tab=$1 t=$2 slow=$3
  shift 3
  local w=$(( t > 1 ? t - 1 : 1 ))
  for d in Z AZ P; do
    twice run "$tab" "$t" "$slow" "$@" --threads "$t" --algo hda \
          --runtime thread --dist "$d"
  done
  for r in openmp fastflow; do
    for d in Z AZ; do
      twice run "$tab" "$t" "$slow" "$@" --threads "$t" --algo hda \
            --runtime "$r" --dist "$d"
    done
  done
  twice run "$tab" "$t" "$slow" "$@" --threads "$t" --algo pbnf
  twice run "$tab" "$t" "$slow" "$@" --threads "$w" --algo farm
  for d in Z AZ; do
    for k in 1 2; do
      run_mpi "$tab" "$t" "$slow" "$@" --dist "$d"
    done
  done
  }

header="table,cores,slow,domain,instance,algo,runtime,dist,threads,batch,"
header+="flush,open,status,cost,expanded,reexpanded,generated,duplicates,"
header+="sent,messages,stored,LB,CO,seconds"
for t in $(echo "$TABLES" | grep -o .); do
  [ -f "$OUT/$t.csv" ] || echo "$header" > "$OUT/$t.csv"
done

# A: sequential A* on the 100 instances of Korf, one process per instance;
# the buckets run again after the heap, the A/A of the table
if [[ $TABLES == *A* ]]; then
  for open in buckets heap buckets; do
    for i in $(seq 1 100); do
      run A 1 0 --file "$K" --ids "$i" --open "$open"
    done
  done
  if [ -n "${BURNS_SOLVER:-}" ]; then
    i=0
    grep -v '^#' "$K" | grep -v '^\s*$' | while read -r -a b; do
      i=$(( i + 1 ))
      { echo "4 4"; echo "starting positions for each tile:"
        for t in $(seq 0 15); do
          for c in $(seq 0 15); do [ "${b[$c]}" = "$t" ] && echo "$c"; done
        done
        echo "goal positions:"; seq 0 15; } > "$OUT/burns.in"
      numactl --physcpubind=$C0 --membind=$NODE "$BURNS_SOLVER" astar \
              < "$OUT/burns.in" 2>&1 | awk -v i="$i" -F'"' '
        /final sol cost/ { c = $4 } /total nodes expanded/ { e = $4 }
        /total nodes generated/ { g = $4 } /total wall time/ { s = $4 }
        END { printf "A,1,0,tiles15,%d,burns,,,1,,,,solved,%d,%d,0,%d,0,0,0," \
                     "0,1,0,%s\n", i , c , e , g , s }' >> "$OUT/A.csv"
    done
    rm -f "$OUT/burns.in"
  fi
fi

# B: the parallel searches on 12 instances of Korf, SPA* on 6 of them
if [[ $TABLES == *B* ]]; then
  I=1,3,10,14,15,22,26,27,32,33,37,52
  for t in $THREADS; do
    twice run B "$t" 0 --file "$K" --ids "$I" --threads "$t" --algo astar
    parallel B "$t" 0 --file "$K" --ids "$I"
    twice run B "$t" 0 --file "$K" --ids 1,3,10,14,15,26 --threads "$t" \
          --algo spa
  done
fi

# C: random grids of 5000 x 5000 cells, 35% blocked, with 4 and 8 moves, as
# in Burns et al. (JAIR 2010); the seeds are those whose corners connect
if [[ $TABLES == *C* ]]; then
  for moves in 4 8; do
    for seed in 1 2 5; do
      g=(--domain grid$moves --random 5000x5000:0.35:$seed)
      twice run C 1 0 "${g[@]}" --algo astar
      for t in $THREADS; do
        parallel C "$t" 0 "${g[@]}"
      done
    done
  done
fi

# D: slow expansions on an easy instance of Korf, with PA*SE
if [[ $TABLES == *D* ]]; then
  for slow in 1000 10000; do
    for t in $THREADS; do
      twice run D "$t" "$slow" --file "$K" --ids 79 --threads "$t" \
            --algo astar
      twice run D "$t" "$slow" --file "$K" --ids 79 --threads "$t" \
            --algo pase
      twice run D "$t" "$slow" --file "$K" --ids 79 --threads "$t" \
            --algo hda --dist AZ
      twice run D "$t" "$slow" --file "$K" --ids 79 --threads "$t" \
            --algo spa
      twice run D "$t" "$slow" --file "$K" --ids 79 \
            --threads $(( t > 1 ? t - 1 : 1 )) --algo farm
    done
  done
fi

# E: the examples, ten runs each
if [[ $TABLES == *E* ]]; then
  X=$BUILD/examples
  echo "example,arguments,run,cost,expanded,seconds" > "$OUT/E.csv"
  while IFS='|' read -r name args; do
    for k in $(seq 1 10); do
      numactl --physcpubind=$C0 --membind=$NODE "$X/$name" $args |
        awk -v n="$name" -v a="$args" -v k="$k" '
          /^cost:/ { c = $2 } /^moves:/ { m = $2 } /^expanded/ { e = $2 }
          / km\)$/ { km = $(NF - 1); sub( /\(/ , "" , km ) }
          /^seconds/ { s = $2 }
          END { if( c == "" ) c = km != "" ? km : m
                printf "%s,%s,%d,%s,%s,%s\n", n , a , k , c , e , s }' \
          >> "$OUT/E.csv"
    done
  done <<'LIST'
8Puzzle|
8Puzzle|-h h1
8Puzzle|worst
15Puzzle|
RouteFinding|
RouteFinding|Lugoj Neamt
NQueens|12
WaterJugs|
TravelingSalesman|16 1
FindPath|-r 1000x1000:0.3:2
FindPath|-8 -r 1000x1000:0.3:2
LIST
fi

###############################################################################
############################# End of campaign.sh ##############################
###############################################################################
