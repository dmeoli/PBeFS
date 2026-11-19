#!/bin/bash
###############################################################################
#                                                                             #
#   rubik.sh: IDA* and parallel IDA* on the ten instances of Korf (1997)      #
#                                                                             #
#   Usage:  rubik.sh <build> <out> [<first core> <NUMA node>]                 #
#                                                                             #
#   Runs tools/AStarBench of the build <build> on the scrambles of            #
#   data/rubik/korf10.txt with the tables of Korf (the corners and two sets   #
#   of 6 edges, their maximum): sequential IDA* on the instances IDS          #
#   (default 1-10), then the parallel IDA* on each number of workers of       #
#   THREADS (default 1 2 4 8 16 32). Each run is a process of its own,        #
#   pinned to the cores <first core> , ... of the NUMA node <NUMA node>       #
#   (default 96 and 1), with the memory on that node, at most MEMKB           #
#   kilobytes of virtual memory (default 16000000) and TIMEOUT of time        #
#   (default 7d), never under nice. The tables are built once in CACHE        #
#   (default <out>/cache), or copied from CACHE_FROM, and pdb.csv in <out>    #
#   holds the time and the peak memory of the build with the size and the     #
#   md5 of each file; each run loads them, and its time excludes that load,   #
#   which is a column of its own. Before the campaign, the runs of AA         #
#   (default "ida:1 pida:32", algorithm:workers, empty for none) on the       #
#   instances AA_IDS (default 1) are made twice, the gap between the two      #
#   being the A/A noise of the machine. results.csv holds a line per run,     #
#   labelled aa1, aa2 or run; a run already there is not made again, so that  #
#   the script resumes where it stopped.                                      #
#                                                                             #
#   The script stops (exit 3) if ~/TIMED-RUN-IN-PROGRESS exists, or if the    #
#   machine is not at rest: before the campaign, PRECHECKS checks (default    #
#   3) PREGAP seconds apart (default 600), and before each run one check,     #
#   must find the load average of a minute below LOAD (default 6) and the     #
#   CPUs busy less than BUSY percent (default 2) of a minute by /proc/stat;   #
#   while a run goes, if the CPUs busy exceed its workers by more than        #
#   MARGIN (default 4) for 3 minutes in a row, the run is stopped, left out   #
#   of results.csv, and the script stops. launch.sh starts it detached.       #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

set -u
BUILD=$(realpath "$1")
OUT=$(realpath -m "$2")
C0=${3:-96}
NODE=${4:-1}
SRC=$(realpath "$(dirname "$0")/../..")
B=$BUILD/tools/AStarBench
K=$SRC/data/rubik/korf10.txt
IDS=${IDS:-1-10}
THREADS=${THREADS:-"1 2 4 8 16 32"}
AA=${AA-"ida:1 pida:32"}
AA_IDS=${AA_IDS:-1}
TIMEOUT=${TIMEOUT:-7d}
MEMKB=${MEMKB:-16000000}
CACHE=${CACHE:-$OUT/cache}
LOAD=${LOAD:-6}
BUSY=${BUSY:-2}
MARGIN=${MARGIN:-4}
PRECHECKS=${PRECHECKS:-3}
PREGAP=${PREGAP:-600}
NCPU=$(nproc --all)
LOG=$OUT/rubik.log
RES=$OUT/results.csv
mkdir -p "$OUT/runs" "$CACHE"
export PBEFS_PIN=1

say() { echo "$(date -Iseconds) $*" >> "$LOG"; }
stop() { say "stop: $*"; exit 3; }
[ "$(nice)" = 0 ] || stop "running under nice $(nice)"
[ -x "$B" ] || stop "no $B"

# the indices of a list a-b,c,...
expand() {
  echo "$1" | tr ',' '\n' | while IFS=- read -r a b; do seq "$a" "${b:-$a}"
  done
}

# the CPUs busy over <seconds>, in hundredths of a CPU, from /proc/stat
busy() {
  local a b
  a=$(awk '/^cpu / { print $2+$3+$4+$7+$8+$9, $5+$6 }' /proc/stat)
  sleep "$1"
  b=$(awk '/^cpu / { print $2+$3+$4+$7+$8+$9, $5+$6 }' /proc/stat)
  echo "$a $b" | awk -v n="$NCPU" '{ w = $3 - $1; i = $4 - $2;
    print int(100 * n * w / (w + i)) }'
}

# true if the machine is at rest: the marker absent, the load average of a
# minute below LOAD (waited for up to 10 minutes, as that of the run just
# ended decays) and the CPUs busy less than BUSY percent of a minute
at_rest() {
  [ -e ~/TIMED-RUN-IN-PROGRESS ] && stop "$HOME/TIMED-RUN-IN-PROGRESS exists"
  local k l
  for k in $(seq 1 20); do
    l=$(cut -d' ' -f1 /proc/loadavg)
    awk -v l="$l" -v m="$LOAD" 'BEGIN { exit !(l < m) }' && break
    sleep 30
  done
  awk -v l="$l" -v m="$LOAD" 'BEGIN { exit !(l < m) }' ||
    { say "load average $l"; return 1; }
  k=$(busy 60)
  [ "$k" -lt $(( BUSY * NCPU )) ] || { say "busy $k/100 CPUs"; return 1; }
  REST="$l,$k"
}

# stops the run <pid> of <workers> workers if the CPUs busy exceed them by
# more than MARGIN for 3 minutes in a row
guard() {
  local n=0 k
  while kill -0 "$1" 2>/dev/null; do
    k=$(busy 60)
    if [ "$k" -gt $(( ( $2 + MARGIN ) * 100 )) ]; then n=$(( n + 1 ))
    else n=0; fi
    if [ "$n" -ge 3 ]; then
      say "busy $k/100 CPUs during the run $3"
      touch "$OUT/busy"
      kill "$1"
      return
    fi
  done
}

# one run: <label> <algo> <workers> <instance>
run() {
  local tag=$1-$2-$3-$4 pid wpid rc st
  grep -q "^$1,$4,$2,$3," "$RES" && return
  at_rest || stop "machine not at rest before $tag"
  say "start $tag"
  ( ulimit -v "$MEMKB"
    exec timeout "$TIMEOUT" numactl --physcpubind="$C0-$(( C0 + $3 - 1 ))" \
         --membind="$NODE" "$B" --domain rubik --file "$K" --ids "$4" \
         --pdb korf --cache "$CACHE" --pdbinfo 1 --algo "$2" \
         --threads "$3" --header ) > "$OUT/runs/$tag.out" \
    2> "$OUT/runs/$tag.err" &
  pid=$!
  guard "$pid" "$3" "$tag" &
  wpid=$!
  wait "$pid"
  rc=$?
  kill "$wpid" 2>/dev/null
  wait "$wpid" 2>/dev/null
  [ -e "$OUT/busy" ] && { rm -f "$OUT/busy"; stop "busy during $tag"; }
  case $rc in 0) st=ok;; 124) st=timeout;; *) st=error;; esac
  say "end $tag rc $rc"
  awk -F, -v lab="$1" -v id="$4" -v a="$2" -v t="$3" -v st="$st" \
      -v rc="$rc" -v rest="$REST" -v err="$OUT/runs/$tag.err" '
    BEGIN { while ((getline l < err) > 0) { split(l, p, ",")
              if (p[1] == "pdb" && p[2] == "all") pdb = p[3] } }
    NR == 1 { for (i = 1; i <= NF; ++i) c[$i] = i; next }
    { print lab, id, a, t, st == "ok" ? $c["status"] : st, $c["cost"],
            $c["expanded"], $c["generated"], $c["seconds"], pdb, rc, rest
      done = 1 }
    END { if (!done) print lab, id, a, t, st, "", "", "", "", pdb, rc, rest
        }' OFS=, "$OUT/runs/$tag.out" >> "$RES"
}

[ -e ~/TIMED-RUN-IN-PROGRESS ] && stop "$HOME/TIMED-RUN-IN-PROGRESS exists"
say "start on $(hostname), cores $C0-, node $NODE, $NCPU CPUs"
for k in $(seq 1 "$PRECHECKS"); do
  at_rest || stop "machine not at rest (check $k)"
  say "at rest, load and busy $REST"
  [ "$k" -lt "$PRECHECKS" ] && sleep "$PREGAP"
done

# the tables: copied, or built once with their time and peak memory
[ -f "$RES" ] || echo "label,instance,algo,threads,status,length,expanded,"\
"generated,seconds,pdb_seconds,rc,load,busy" > "$RES"
if [ -z "$(ls "$CACHE"/*.pdb 2>/dev/null)" ]; then
  if [ -n "${CACHE_FROM:-}" ]; then
    cp "$CACHE_FROM"/*.pdb "$CACHE"/
    say "tables copied from $CACHE_FROM"
  else
    t=$(( NCPU < 32 ? NCPU : 32 ))
    numactl --physcpubind="$C0-$(( C0 + t - 1 ))" --membind="$NODE" \
      /usr/bin/time -f "%e,%M" -o "$OUT/pdb-time.txt" "$B" --domain rubik \
      --size 1 --pdb korf --cache "$CACHE" --pdbinfo 1 --threads "$t" \
      --algo ida > /dev/null 2> "$OUT/pdb-build.err" ||
      stop "the tables were not built"
    say "tables built on $t workers"
  fi
  { echo "file,bytes,md5"
    for f in "$CACHE"/*.pdb; do
      m=$(md5sum < "$f" | cut -d' ' -f1)
      echo "$(basename "$f"),$(stat -c %s "$f"),$m"
    done
    [ -f "$OUT/pdb-time.txt" ] &&
      echo "# seconds,maxrss_kb $(cat "$OUT/pdb-time.txt")"
  } > "$OUT/pdb.csv"
fi

# the A/A: the same runs twice
for rep in aa1 aa2; do
  for x in $AA; do
    for i in $(expand "$AA_IDS"); do
      run "$rep" "${x%%:*}" "${x##*:}" "$i"
    done
  done
done

# the campaign: sequential IDA*, then the parallel one on 1, 2, ... workers
for i in $(expand "$IDS"); do run run ida 1 "$i"; done
for t in $THREADS; do
  for i in $(expand "$IDS"); do run run pida "$t" "$i"; done
done
say "done"

###############################################################################
############################### End of rubik.sh ###############################
###############################################################################
