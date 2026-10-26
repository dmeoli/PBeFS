#!/bin/bash
###############################################################################
#                                                                             #
#   structures.sh: A* on the data structures of PBeFS, std, Boost and ADS     #
#                                                                             #
#   Usage:  structures.sh <out> <build own> <build std> <build boost>         #
#                         <build ads> [<first core> <NUMA node>]              #
#                                                                             #
#   Runs tools/AStarBench of the four builds, made with PBEFS_STRUCTURES own, #
#   std, boost and ads, on instances of growing size: the 15-puzzle in four   #
#   bands of the instances of Korf, by the nodes that A* expands (about 1e5,  #
#   1e6, 1e7 and above 7e7), and random grids of 2000 to 20000 cells a side,  #
#   35% blocked, with 4 and 8 moves (the first three seeds whose corners      #
#   connect, which A* of the build own finds). Each instance is solved by     #
#   sequential A* (with both OPEN lists on the 15-puzzle, with the one of its #
#   cost on the grids) and by HDA* on 8 threads, twice with every set, the    #
#   sets taking turns so that a drift of the machine falls on all of them;    #
#   structures.py makes the tables. Each run is a process of its own,         #
#   pinned to the cores <first core> , ... of the NUMA node <NUMA node>       #
#   (default 96 and 1) with its memory on that node, and the two runs of a    #
#   set are its A/A. The campaign holds a lock (flock) on each of its cores,  #
#   in LOCKS (default /tmp/pbefs-bench-locks), so that a second campaign on   #
#   the same cores waits for it. BANDS and SIDES restrict the sizes. Linux    #
#   only: numactl places the runs; launch.sh starts it detached, with its     #
#   PID saved.                                                                #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

set -u
OUT=$(realpath -m "$1")
declare -A BUILD=( [own]=$(realpath "$2") [std]=$(realpath "$3")
                   [boost]=$(realpath "$4") [ads]=$(realpath "$5") )
C0=${6:-96}
NODE=${7:-1}
SRC=$(realpath "$(dirname "$0")/../..")
K=$SRC/data/15-puzzle/korf100.txt
BANDS=${BANDS:-"1e5 1e6 1e7 1e8"}
SIDES=${SIDES:-"2000 5000 10000 20000"}
SETS="own std boost ads"
LOCKS=${LOCKS:-/tmp/pbefs-bench-locks}
mkdir -p "$OUT" "$LOCKS"
export PBEFS_PIN=1

# a lock on each of the 8 cores of the campaign, held until it ends
for c in $(seq "$C0" $(( C0 + 7 ))); do
  exec {fd}> "$LOCKS/core$c"
  flock -n "$fd" || { echo "core $c is busy, waiting" >&2; flock "$fd"; }
done

# the instances of Korf in each band, by the nodes that A* expands
declare -A KORF=( [1e5]=48,55,19,85,47,73,97,31,86,30
                  [1e6]=65,6,39,96,18,74,23,46,62,83
                  [1e7]=98,84,27,52,22,1,32,10,91,33
                  [1e8]=66,49,17,60,82,88 )

header="set,rep,size,domain,instance,algo,runtime,dist,threads,batch,flush,"
header+="open,status,cost,expanded,reexpanded,generated,duplicates,sent,"
header+="messages,stored,LB,CO,seconds"
[ -f "$OUT/structures.csv" ] || echo "$header" > "$OUT/structures.csv"

# runs AStarBench of the set s on t cores, the line prefixed by set, rep
# and size
run() { # set rep size cores options...
  local s=$1 rep=$2 size=$3 t=$4
  shift 4
  numactl --physcpubind=$C0-$(( C0 + t - 1 )) --membind=$NODE \
          "${BUILD[$s]}/tools/AStarBench" "$@" 2>> "$OUT/errors.log" |
    sed "s/^/$s,$rep,$size,/" >> "$OUT/structures.csv"
  }

# an instance with every search, twice with every set in turn; HDA* takes
# the OPEN list of the cost, and the 15-puzzle runs A* with both
all() { # size options...
  local size=$1
  shift
  for rep in 1 2; do
    for s in $SETS; do
      if [ "$size" != "${size#1e}" ]; then
        run "$s" "$rep" "$size" 1 "$@" --algo astar --open buckets
        run "$s" "$rep" "$size" 1 "$@" --algo astar --open heap
      else
        run "$s" "$rep" "$size" 1 "$@" --algo astar
      fi
      run "$s" "$rep" "$size" 8 "$@" --algo hda --threads 8 --dist Z
    done
  done
  }

for b in $BANDS; do
  for i in ${KORF[$b]//,/ }; do
    all "$b" --file "$K" --ids "$i"
  done
done

# the first three seeds whose grid of the given side connects its corners
# with 4 moves, and hence with 8 too
seeds() { # side
  local n=0 seed
  for seed in $(seq 1 50); do
    if "${BUILD[own]}/tools/AStarBench" --domain grid4 --algo astar \
         --random "$1x$1:0.35:$seed" | grep -q ',solved,'; then
      echo "$seed"
      n=$(( n + 1 ))
      [ $n -eq 3 ] && return
    fi
  done
  }

for side in $SIDES; do
  S=$(seeds "$side")
  for moves in 4 8; do
    for seed in $S; do
      all "$side" --domain grid$moves --random "${side}x$side:0.35:$seed"
    done
  done
done
