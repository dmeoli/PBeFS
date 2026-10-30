#!/bin/bash
###############################################################################
#                                                                             #
#   pdb.sh: the pattern databases of the 15-puzzle and of the 24-puzzle       #
#                                                                             #
#   Usage:  pdb.sh <out> <build> [<ids>]                                      #
#                                                                             #
#   Runs tools/PDBBench of the build <build>: first the mode table on each    #
#   partition and encoding (the size of the tables, the time to build them,   #
#   the mean heuristic on random states and the time of a lookup), the        #
#   tables kept in <out>/cache, then A* on the instances <ids> of Korf (all   #
#   by default) with the Manhattan distance and with each heuristic of the    #
#   tables, twice, the second run being the A/A of the first. The             #
#   expansions do not depend on the machine; the times are comparable only    #
#   on an idle one. Last, IDA* with the tables 6-6-6-6 on 50 instances of     #
#   the 24-puzzle (or those of IDS24), each at most LIMIT24 expansions        #
#   (default 1e10): the random walks of WALK24 moves from the goal (default   #
#   150), or, with WALK24 empty, the instances of data/24-puzzle. JOBS runs   #
#   go at once (default 4), each one a process of its own; table.csv,         #
#   search.csv and search24.csv in <out> hold the lines.                      #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

set -u
OUT=$(realpath -m "$1")
BIN=$(realpath "$2")/tools/PDBBench
IDS=${3:-1-100}
JOBS=${JOBS:-4}
SRC=$(realpath "$(dirname "$0")/../..")
K=$SRC/data/15-puzzle/korf100.txt
mkdir -p "$OUT/cache" "$OUT/runs"

# the tables: partition, encoding, ranks per entry, rule of the blank,
# reflection
TABLES="555:byte:1:ignored:1 555:nibble:1:ignored:1 555:mod3:1:ignored:0
        555:bit:1:ignored:0 555:min:4:ignored:1 555:min:16:ignored:1
        555:nibble:1:minimized:1 663:nibble:1:ignored:1
        78:nibble:1:ignored:1 78:bit:1:ignored:0 78:min:8:ignored:1
        78:nibble:1:minimized:1 6666:nibble:1:ignored:1
        6666:bit:1:ignored:0"
if "$BIN" --mode table --enc sparse --heur 555 --samples 1 > /dev/null \
   2>&1; then
  TABLES="$TABLES 555:sparse:1:ignored:1 78:sparse:1:ignored:1"
fi

echo "heur,enc,compress,blank,reflect,entries,bytes,bits,source,seconds,"\
"mean_h,mean_manhattan,ns_lookup,sink" > "$OUT/table.csv"
for t in $TABLES; do
  IFS=: read -r h e c b r <<< "$t"
  "$BIN" --mode table --heur "$h" --enc "$e" --compress "$c" --blank "$b" \
         --reflect "$r" --cache "$OUT/cache" >> "$OUT/table.csv" ||
    echo "table $t failed" >&2
done

# the searches: each heuristic twice, the tables loaded from the cache
SEARCHES="manhattan:nibble:1:ignored:1 $TABLES"
run() {  # <rep> <table>
  IFS=: read -r h e c b r <<< "$2"
  "$BIN" --file "$K" --ids "$IDS" --heur "$h" --enc "$e" --compress "$c" \
         --blank "$b" --reflect "$r" --cache "$OUT/cache" |
    sed "s/^/$1,/" > "$OUT/runs/$1-$2.csv"
}
export -f run
export BIN K IDS OUT
for rep in 1 2; do
  for t in $SEARCHES; do
    case $t in 6666:*|*:min:4:*|*:byte:*|*:mod3:*) ;; *) echo "$rep $t";;
    esac
  done
done | xargs -P "$JOBS" -L 1 bash -c 'run "$0" "$1"'
echo "rep,instance,heur,enc,compress,blank,reflect,algo,threads,status,"\
"cost,expanded,reexpanded,generated,seconds" > "$OUT/search.csv"
cat "$OUT"/runs/*.csv >> "$OUT/search.csv"

# the 24-puzzle: IDA* with the tables 6-6-6-6, each instance of IDS24 (all
# 50 by default) once, as its expansions do not depend on the machine, and
# at most LIMIT24 expansions (default 1e10) each; the instances are the
# random walks of WALK24 moves, or those of the file if WALK24 is empty
WALK24=${WALK24-150}
if [ -n "$WALK24" ]; then
  SET24="--walk $WALK24"
else
  SET24="--file $SRC/data/24-puzzle/instances50.txt"
fi
IDS24=${IDS24:-1-50}
LIMIT24=${LIMIT24:-10000000000}
run24() {  # <instance>
  "$BIN" --puzzle 24 $SET24 --ids "$1" --heur 6666 --enc nibble \
         --algo ida --limit "$LIMIT24" --cache "$OUT/cache" \
    > "$OUT/runs24/$1.csv"
}
export -f run24
export SET24 LIMIT24
mkdir -p "$OUT/runs24"
"$BIN" --mode table --heur 6666 --enc nibble --samples 1 \
       --cache "$OUT/cache" > /dev/null                # the tables, once
for i in $(seq 1 50); do echo "$i"; done |
  awk -v s="$IDS24" 'BEGIN { n = split(s, p, ","); for (k = 1; k <= n; ++k) {
    m = split(p[k], r, "-"); a = r[1]; b = m > 1 ? r[2] : r[1];
    for (i = a; i <= b; ++i) ok[i] = 1 } } ok[$1]' |
  xargs -P "$JOBS" -I{} bash -c 'run24 {}'
echo "instance,heur,enc,compress,blank,reflect,algo,threads,status,cost,"\
"expanded,reexpanded,generated,seconds" > "$OUT/search24.csv"
cat "$OUT"/runs24/*.csv >> "$OUT/search24.csv"
