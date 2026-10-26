#!/bin/bash
###############################################################################
#                                                                             #
#   launch.sh: a campaign started detached, with its PID saved                #
#                                                                             #
#   Usage:  launch.sh <out> <campaign> [<arguments>...]                       #
#                                                                             #
#   Starts <campaign> (campaign.sh or structures.sh, with its arguments) by   #
#   nohup and setsid, so that it survives the end of the session that         #
#   started it and none of its signals, with its output in <out>/launch.log   #
#   and its PID, which is also the group of all its processes, in             #
#   <out>/launch.pid; it is stopped by kill -- -<PID>, and by no kill of      #
#   processes by their name, which would take those of other users too.       #
#   The timings are to be taken on a machine at rest, with the cores and the  #
#   NUMA node given to the campaign not used by anything else, and never on   #
#   a virtual machine whose memory is given and taken back while it runs      #
#   (e.g., wolfe, a VM of Hyper-V whose dynamic memory stays about 2.4 GB     #
#   for hours, so that any large instance swaps): its timings are not         #
#   comparable with each other.                                               #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

set -eu
OUT=$(realpath -m "$1")
shift
mkdir -p "$OUT"
if [ -f "$OUT/launch.pid" ] && kill -0 "$(cat "$OUT/launch.pid")" 2>/dev/null
then
  echo "launch.sh: a campaign runs in $OUT, PID $(cat "$OUT/launch.pid")" >&2
  exit 1
fi
{
  echo "# host     $(hostname)"
  echo "# date     $(date -Iseconds)"
  echo "# load     $(cut -d' ' -f1-3 /proc/loadavg)"
  echo "# command  $*"
} > "$OUT/launch.log"
nohup setsid "$@" >> "$OUT/launch.log" 2>&1 < /dev/null &
echo $! > "$OUT/launch.pid"
echo "PID $(cat "$OUT/launch.pid"), log $OUT/launch.log"

###############################################################################
############################## End of launch.sh ###############################
###############################################################################
