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
#   of them), C (random grids), D (slow expansions), E (the examples) and F   #
#   (every parallel search on every problem of the examples, A* and           #
#   Branch-and-Bound, the latter by tools/BranchAndBoundBench). Each line of  #
#   a CSV starts with the table (and, in F, the problem), the cores given to  #
#   the run and the artificial work per expansion, then the line of the       #
#   benchmark; F_LIST may name a file that replaces the list of the problems  #
#   of F [see below].                                                         #
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
TABLES=${5:-ABCDEF}
SRC=$(realpath "$(dirname "$0")/../..")
B=$BUILD/tools/AStarBench
M=$BUILD/tools/AStarBenchMPI
BB=$BUILD/tools/BranchAndBoundBench
BM=$BUILD/tools/BranchAndBoundBenchMPI
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

# runs AStarBench (or BIN, if set) on t cores, the line prefixed by table,
# cores and slow; a table "F,<problem>" writes to F.csv
run() { # table cores slow options...
  local tab=$1 t=$2 slow=$3
  shift 3
  numactl --physcpubind=$C0-$(( C0 + t - 1 )) --membind=$NODE \
          "${BIN:-$B}" --slow "$slow" "$@" 2>> "$OUT/errors.log" |
    sed "s/^/$tab,$t,$slow,/" >> "$OUT/${tab%%,*}.csv"
  }

# the same with MPI (AStarBenchMPI, or BIN_MPI if set), one process per
# core
run_mpi() { # table cores slow options...
  local tab=$1 t=$2 slow=$3
  shift 3
  mpirun -np "$t" --bind-to core --cpu-set $C0-$(( C0 + t - 1 )) \
         -x PBEFS_PIN numactl --membind=$NODE "${BIN_MPI:-$M}" \
         --slow "$slow" "$@" 2>> "$OUT/errors.log" |
    sed "s/^/$tab,$t,$slow,/" >> "$OUT/${tab%%,*}.csv"
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
  if [ "$t" = F ]; then
    [ -f "$OUT/F.csv" ] ||
      echo "${header/table,/table,problem,}" > "$OUT/F.csv"
  else
    [ -f "$OUT/$t.csv" ] || echo "$header" > "$OUT/$t.csv"
  fi
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

# F: every parallel search on every problem of the examples, each at the
# size where the sequential search takes about 1 to 10 s on one core, each
# run twice, the sequential search on 1 core being the reference of its
# speedup. A line of the list is
#   problem|bench|slow|searches|options
# with bench astar (AStarBench) or bnb (BranchAndBoundBench), slow the work
# per expansion, searches the distributions of HDA* that the Domain gives,
# the first being the one of the slow runs, and pbnf and pase if it has
# nblocks and a pairwise heuristic, and the options of the bench (the lines
# with the same problem are its instances, and a line ending with a
# backslash goes on in the next one). A line with slow 0 runs HDA*
# with every distribution on the threads, OpenMP and FastFlow, SPA* and
# PBNF on the threads and OpenMP, the farm and HDAStarMPI (every
# distribution but R); with slow > 0 it runs PA*SE on the threads and
# OpenMP, and HDA*, SPA* and the farm on the threads, as table D. HDA*[R]
# detects a duplicate only when its copies reach the same worker, and it
# runs with at most R_LIMIT nodes in all (default 200000000), the cell
# being empty if they are not enough. A bnb line runs the sequential B&B
# as the reference, every scheme on every runtime, the farm and the B&B
# over MPI. At the end, the cost of every run solved is checked against the
# reference of its instance, and the runs whose cost differs are written to
# errors.log and to F-mismatch.csv.
if [[ $TABLES == *F* ]]; then
  P8=$OUT/8puzzle.txt                # the 8-puzzle of figure 3.28, and the
  printf '%s\n' "7 2 4 5 0 6 8 3 1" "8 0 6 5 4 7 2 3 1" > "$P8"  # hardest
  RL=${R_LIMIT:-200000000}
  if [ -n "${F_LIST:-}" ]; then
    LIST=$(cat "$F_LIST")
  else
    LIST=$(cat <<'LIST'
tiles15|astar|0|Z AZ A P H R pbnf|--file "$K" --ids 1,10,14
tiles8|astar|0|Z AZ A P H R pbnf|--domain tiles8 --file "$P8"
romania|astar|0|Z P H R|--domain romania
grid4|astar|0|Z AZ A P H R pbnf|--domain grid4 --random 5000x5000:0.25:1:5
grid4|astar|0|Z AZ A P H R pbnf|--domain grid4 --random 5000x5000:0.25:2:5
grid4|astar|0|Z AZ A P H R pbnf|--domain grid4 --random 5000x5000:0.25:4:5
grid8|astar|0|Z AZ A P H R pbnf|--domain grid8 --random 5000x5000:0.25:1:5
grid8|astar|0|Z AZ A P H R pbnf|--domain grid8 --random 5000x5000:0.25:2:5
grid8|astar|0|Z AZ A P H R pbnf|--domain grid8 --random 5000x5000:0.25:4:5
vacuum|astar|0|Z AZ A P H R pbnf|--domain vacuum --size 22 --ids 5,8,11
queens|astar|0|Z H R|--domain queens --size 28
queens|astar|0|Z H R|--domain queens --size 29
queens|astar|0|Z H R|--domain queens --size 31
missionaries|astar|0|Z P H R|--domain missionaries --size 65535:65535:100
missionaries|astar|0|Z P H R|--domain missionaries --size 65535:65535:200
jugs|astar|0|Z P H R|--domain jugs --size 1:211:199:197:193
tsp|astar|0|Z P H R|--domain tsp --size 32 --ids 3,6,7
tiles15-slow|astar|10000|AZ Z A P H R pbnf pase|--file "$K" --ids 79
tiles8-slow|astar|10000|AZ Z A P H R pbnf pase|--domain tiles8 --file "$P8"
grid4-slow|astar|1000|AZ Z A P H R pbnf pase|--domain grid4 \
  --random 1000x1000:0.25:1:5
grid8-slow|astar|1000|AZ Z A P H R pbnf pase|--domain grid8 \
  --random 1000x1000:0.25:1:5
vacuum-slow|astar|1000|AZ Z A P H R pbnf pase|--domain vacuum --size 16 --ids 8
knapsack|bnb|0||--domain knapsack --size 2:90:1000 --ids 3,8,13
tsp-bb|bnb|0||--domain tsp --size 24 --ids 2,3,4
LIST
)
  fi
  while IFS='|' read -u 3 p bench slow feats opts; do
    { [ -z "$p" ] || [ "${p:0:1}" = "#" ]; } && continue
    eval "o=($opts)"
    tab="F,$p"
    if [ "$bench" = bnb ]; then
      BIN=$BB twice run "$tab" 1 0 "${o[@]}" --algo sequential
      for t in $THREADS; do
        for s in shared distributed deterministic selfsplit stealing; do
          for r in thread openmp fastflow; do
            BIN=$BB twice run "$tab" "$t" 0 "${o[@]}" --threads "$t" \
                --algo "$s" --runtime "$r"
          done
        done
        BIN=$BB twice run "$tab" "$t" 0 "${o[@]}" \
            --threads $(( t > 1 ? t - 1 : 1 )) --algo farm
        BIN_MPI=$BM twice run_mpi "$tab" "$t" 0 "${o[@]}"
      done
      continue
    fi
    dists=()
    pbnf=0
    pase=0
    for x in $feats; do
      case $x in
        pbnf) pbnf=1 ;;
        pase) pase=1 ;;
        *) dists+=("$x") ;;
      esac
    done
    twice run "$tab" 1 "$slow" "${o[@]}" --algo astar
    for t in $THREADS; do
      w=$(( t > 1 ? t - 1 : 1 ))
      if [ "$slow" -gt 0 ]; then
        if [ "$pase" = 1 ]; then
          for r in thread openmp; do
            twice run "$tab" "$t" "$slow" "${o[@]}" --threads "$t" \
                  --algo pase --runtime "$r"
          done
        fi
        twice run "$tab" "$t" "$slow" "${o[@]}" --threads "$t" --algo hda \
              --dist "${dists[0]}"
        twice run "$tab" "$t" "$slow" "${o[@]}" --threads "$t" --algo spa
        twice run "$tab" "$t" "$slow" "${o[@]}" --threads "$w" --algo farm
        continue
      fi
      for r in thread openmp fastflow; do
        for d in "${dists[@]}"; do
          lim=0
          [ "$d" = R ] && lim=$(( RL / t ))
          twice run "$tab" "$t" 0 "${o[@]}" --threads "$t" --algo hda \
                --runtime "$r" --dist "$d" --limit "$lim"
        done
      done
      for r in thread openmp; do
        twice run "$tab" "$t" 0 "${o[@]}" --threads "$t" --algo spa \
              --runtime "$r"
        if [ "$pbnf" = 1 ]; then
          twice run "$tab" "$t" 0 "${o[@]}" --threads "$t" --algo pbnf \
                --runtime "$r"
        fi
      done
      twice run "$tab" "$t" 0 "${o[@]}" --threads "$w" --algo farm
      for d in "${dists[@]}"; do
        [ "$d" = R ] || twice run_mpi "$tab" "$t" 0 "${o[@]}" --dist "$d"
      done
    done
  done 3<<< "$LIST"
  # the cost of each run solved against the reference of its instance
  awk -F, 'NR == FNR { if( ( $14 == "solved" ) &&
                           ( ( $7 == "astar" ) || ( $7 == "sequential" ) ) )
                        ref[ $2 SUBSEP $4 SUBSEP $6 ] = $15
                      next }
           FNR > 1 && $14 == "solved" {
             k = $2 SUBSEP $4 SUBSEP $6
             d = $15 - ref[ k ]
             if( ! ( k in ref ) || ( d * d > 1e-16 * ( 1 + ref[ k ] ^ 2 ) ) )
               print }' "$OUT/F.csv" "$OUT/F.csv" > "$OUT/F-mismatch.csv"
  if [ -s "$OUT/F-mismatch.csv" ]; then
    echo "F: $(wc -l < "$OUT/F-mismatch.csv") runs differ from the" \
         "sequential cost, in F-mismatch.csv" | tee -a "$OUT/errors.log" >&2
  fi
fi

###############################################################################
############################# End of campaign.sh ##############################
###############################################################################
