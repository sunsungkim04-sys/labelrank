"""Scheduler files for clusters: a SLURM array script and a plain command list (GNU parallel, xargs).

Each task runs one run folder's cmd.sh and writes the same DONE marker `labelrank run` writes, so
`labelrank report` reads the result either way. Sequential stopping needs `labelrank run`, which decides what to
launch next; an array launches everything.
"""
import os
import time

from .runner import output_folders

SBATCH = """#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --array=1-{n}%{maxpar}
#SBATCH --cpus-per-task=1
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --output="{logdir}/%A_%a.log"
{extra}
# written by labelrank {version} on {date}; one line of {listname} per array task
d=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {listfile})
[ -n "$d" ] || {{ echo "no run folder for task $SLURM_ARRAY_TASK_ID"; exit 1; }}
exec bash {runone} "$d"
"""

# runs one folder and leaves the DONE marker `labelrank run` leaves (same fields, same format)
RUN_ONE = """#!/bin/bash
d="$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -e "$d/DONE" ] && { echo "already done: $d"; exit 0; }
# the run folder lock 'labelrank run' uses too ('host pid ...'). A stale one (no holder alive here) is replaced
# only under a second lock, $lk.takeover, so two tasks that both found it stale cannot delete each other's new lock
lk="$d/RUN_LOCK"; tk="$lk.takeover"; me="$(hostname) $$"
mk() { ( set -o noclobber; echo "$me" > "$1" ) 2>/dev/null; }
stale() {
  [ -e "$lk" ] || return 0
  read -r h p1 p2 < "$lk" 2>/dev/null
  t=$(stat -c %Y "$lk" 2>/dev/null || stat -f %m "$lk" 2>/dev/null || date +%s)
  [ -z "$h" ] && [ $(( $(date +%s) - t )) -lt 10 ] && return 1
  [ -n "$h" ] && { [ "$h" != "$(hostname)" ] || kill -0 "$p1" 2>/dev/null || { [ -n "$p2" ] && kill -0 "$p2" 2>/dev/null; }; } && return 1
  return 0
}
got=
for i in 1 2 3 4 5 6 7 8 9 10; do
  mk "$lk" && { got=1; break; }
  stale || break
  if mk "$tk"; then
    stale && rm -f "$lk"
    mk "$lk" && got=1
    [ "$(cat "$tk" 2>/dev/null)" = "$me" ] && rm -f "$tk"
    [ -n "$got" ] && break
  fi
  sleep 1
done
[ -n "$got" ] || { echo "busy: $d (see $lk)"; exit 0; }
release() { [ "$(cat "$lk" 2>/dev/null)" = "$me" ] && rm -f "$lk"; }
# finished by another process between the check above and the lock
[ -e "$d/DONE" ] && { release; echo "already done: $d"; exit 0; }
stop() {   # TERM the run's process group, then KILL what is left, so no engine outlives this script
  kill -TERM -- -$pid 2>/dev/null
  for i in 1 2 3 4 5; do kill -0 -- -$pid 2>/dev/null || break; sleep 1; done
  kill -KILL -- -$pid 2>/dev/null
  wait $pid 2>/dev/null
}
# a run cancelled or stopped by the scheduler gets no DONE, so the next submission starts it again
trap 'stop; release; echo "interrupted: $d (no DONE written)"; exit 143' TERM INT
rm -f "$d/relax.json" "$d/relax.out" "$d/errors.log" "$d/alt/mlc" "$d/null/mlc" "$d/alt/rst" "$d/null/rst"
s=$(date +%s)
set -m                      # job control: the run gets its own process group, which the trap can stop whole
bash "$d/cmd.sh" & pid=$!
me="$me $pid"; echo "$me" > "$lk"
wait $pid; rc=$?
if [ -s "$d/relax.json" ] || [ -s "$d/alt/mlc" ]; then st=ok; else st=FAIL; fi
[ $rc -ne 0 ] && st=FAIL
printf "%s\\t%d\\t%ds\\t%s\\n" "$st" "$rc" "$(( $(date +%s) - s ))" "$(date '+%Y-%m-%d %H:%M:%S')" > "$d/DONE"
release
exit $rc
"""

MAX_ARRAY = 1000   # SLURM's default MaxArraySize is 1001; larger lists are split into several arrays


def _q(x):
    return "'" + x.replace("'", "'\\''") + "'"


def write(path, partition=None, account=None, time_limit="48:00:00", mem="2G", max_parallel=100, log=print):
    from . import __version__
    import json
    path = os.path.abspath(path)
    dirs, seq = [], False
    for o in output_folders(path):
        with open(os.path.join(o, "setup.json")) as fh:
            seq = seq or bool(json.load(fh).get("sequential_h"))
        rd = os.path.join(o, "runs")
        dirs += sorted(os.path.join(rd, x) for x in os.listdir(rd)
                       if not os.path.exists(os.path.join(rd, x, "DONE")))
    if not dirs:
        log("nothing to run: every run folder has a DONE marker (failed ones too: delete a failed run's DONE to "
            "retry it); no scheduler files written")
        return None
    # SLURM reads %% in --output as a pattern and the path sits in double quotes there; refuse what would break it
    bad = [c for c in '"%\n$`\\' if c in path]
    if bad:
        raise ValueError("the folder path %r contains %s, which an #SBATCH line cannot hold; move the folder to a "
                         "plainer path" % (path, " ".join(repr(c) for c in bad)))
    sd = os.path.join(path, "slurm")
    os.makedirs(os.path.join(sd, "logs"), exist_ok=True)
    runone = os.path.join(sd, "run_one.sh")
    with open(runone, "w") as fh:
        fh.write(RUN_ONE)
    with open(os.path.join(sd, "commands.txt"), "w") as fh:
        fh.write("".join("bash %s %s\n" % (_q(runone), _q(d)) for d in dirs))
    extra = "".join("#SBATCH --%s=%s\n" % (k, v) for k, v in (("partition", partition), ("account", account)) if v)
    chunks = [dirs[i:i + MAX_ARRAY] for i in range(0, len(dirs), MAX_ARRAY)]
    names = []
    for i, ch in enumerate(chunks, 1):
        suffix = "" if len(chunks) == 1 else "_%03d" % i
        listfile = os.path.join(sd, "runs%s.txt" % suffix)
        with open(listfile, "w") as fh:
            fh.write("".join(d + "\n" for d in ch))
        name = "array%s.sbatch" % suffix
        with open(os.path.join(sd, name), "w") as fh:
            fh.write(SBATCH.format(name="labelrank", n=len(ch), maxpar=max_parallel, time=time_limit, mem=mem,
                                   logdir=os.path.join(sd, "logs"), extra=extra, version=__version__,
                                   date=time.strftime("%Y-%m-%d"), listfile=_q(listfile),
                                   listname=os.path.basename(listfile), runone=_q(runone)))
        names.append(name)
    log("%d unfinished run folder(s); wrote %s in %s (sbatch each) and commands.txt. The files hold absolute "
        "paths: rerun 'labelrank slurm' if the folder moves. At most %d tasks run at once (--max-parallel)." % (
            len(dirs), ", ".join(names), sd, max_parallel))
    if seq:
        log("note: this setup uses sequential stopping, which an array cannot do; every alternative will run. "
            "Use 'labelrank run' to stop early.")
    return sd
