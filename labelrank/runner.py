"""Run every run folder that has no DONE marker, a fixed number at a time, one thread each.

Works on one output folder or on a batch folder (one output folder per gene). With a sequential setup
(--sequential h), a gene's alternatives start only after its observed runs finish, are launched in the pre-drawn
order, and no new one is launched once h of them have reached the observed LRT (Besag-Clifford stopping).
Runs already started finish; the report ignores anything past the stopping point.
"""
import json
import os
import signal
import socket
import subprocess
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait

SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
# outputs of an earlier attempt; removed before a run starts so a stale file can never pass for this run's result
STALE = ("relax.json", "relax.out", "errors.log", "alt/mlc", "null/mlc", "alt/rst", "null/rst", "alt/codeml.log",
         "null/codeml.log")

_live = {}                     # run folder -> Popen of the engine's process group
_live_lock = threading.RLock()   # re-entrant: a second Ctrl-C arrives in the handler while it holds the lock
_stopping = threading.Event()  # set on Ctrl-C / SIGTERM / SIGHUP: no new runs, no DONE for the interrupted ones


def _kill_group(proc):
    """SIGTERM the run's process group, then SIGKILL whatever is left of it: an engine that ignores SIGTERM must
    not outlive its shell and write into a folder that is already marked done or being rerun."""
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except OSError:
        pass
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    for _ in range(50):
        proc.poll()        # reap the shell, so a zombie leader does not keep the group looking alive
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except OSError:
            break          # the group is empty
        time.sleep(0.1)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def run_one(d, timeout):
    if os.path.exists(os.path.join(d, "DONE")):
        return "skip", d
    if _stopping.is_set():
        return "interrupted", d
    # a run folder lock shared with the scheduler's run_one.sh: a folder another process is running is left alone
    lk, me = os.path.join(d, "RUN_LOCK"), "%s %d" % (socket.gethostname(), os.getpid())
    if not take(lk, me):
        return "busy", d
    try:
        if os.path.exists(os.path.join(d, "DONE")):   # finished by another process between the check and the lock
            return "skip", d
        for f in STALE:
            try:
                os.remove(os.path.join(d, f))
            except OSError:
                pass
        env = dict(os.environ, **SINGLE_THREAD)
        t0 = time.time()
        # its own process group, so a timeout or an interrupt stops the engine itself, not only the shell around it
        proc = subprocess.Popen(["bash", os.path.join(d, "cmd.sh")], env=env, start_new_session=True)
        with _live_lock:
            _live[d] = proc
        if _stopping.is_set():   # Ctrl-C between the check above and here: the handler did not see this run
            _kill_group(proc)
        with open(lk, "w") as fh:   # the engine's group too: if this runner is killed outright, it still holds
            fh.write("%s %d\n" % (me, proc.pid))
        me = "%s %d" % (me, proc.pid)
        st = None
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_group(proc)
            rc, st = 124, "TIMEOUT"
        finally:
            with _live_lock:
                _live.pop(d, None)
        return _finish(d, rc, st, t0)
    finally:
        release(lk, me)


def _finish(d, rc, st, t0):
    if _stopping.is_set() and st is None and rc != 0:
        return "interrupted", d   # no DONE: the next 'labelrank run' starts it again from scratch
    ok = any(os.path.exists(os.path.join(d, f)) and os.path.getsize(os.path.join(d, f)) > 0
             for f in ("relax.json", "alt/mlc"))
    st = st or ("ok" if ok and rc == 0 else "FAIL")
    with open(os.path.join(d, "DONE"), "w") as fh:
        fh.write("%s\t%d\t%ds\t%s\n" % (st, rc, time.time() - t0, time.strftime("%Y-%m-%d %H:%M:%S")))
    return st, d


def _alive(content):
    """Whether a lock's holder ('host pid [pid ...]') may still be working: any listed process alive on this host,
    or any other host (which cannot be checked from here)."""
    f = content.split()
    if len(f) < 2:
        return False
    if f[0] != socket.gethostname():
        return True
    for x in f[1:]:
        try:
            os.kill(int(x), 0)
            return True
        except (OSError, ValueError):
            pass
    return False


def _create(p, content):
    try:
        fd = os.open(p, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as fh:
        fh.write(content + "\n")
    return True


def _stale(p):
    """True if lock file p is gone or has no live holder."""
    try:
        with open(p) as fh:
            held = fh.read()
        # created but not yet written by its holder (O_EXCL creates it empty): held, not stale
        return not (_alive(held) or (not held.strip() and time.time() - os.path.getmtime(p) < 10))
    except OSError:
        return True


def take(p, content):
    """Create lock file p atomically; a stale one (no holder alive) is replaced. Returns False if it is held.
    Replacing is done under a second lock, p.takeover, so two processes that both found p stale cannot each
    delete it and then the other's new lock: only the takeover holder deletes, after checking again."""
    tk = p + ".takeover"
    for _ in range(20):
        if _create(p, content):
            return True
        if not _stale(p):
            return False
        if not _create(tk, content):
            try:   # left by a process that died in the middle of a takeover, which takes milliseconds
                if time.time() - os.path.getmtime(tk) > 30:
                    os.remove(tk)
            except OSError:
                pass
            time.sleep(0.05)
            continue
        try:
            if _stale(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
            if _create(p, content):
                return True
        finally:
            release(tk, content)
    return False


def release(p, content):
    """Remove lock file p only if it is still ours."""
    try:
        with open(p) as fh:
            if fh.read().strip() == content:
                os.remove(p)
    except OSError:
        pass


def _lock(folders):
    """One 'labelrank run' per output folder: two would launch every run twice into the same files."""
    me = "%s %d" % (socket.gethostname(), os.getpid())
    taken = []
    for o in folders:
        p = os.path.join(o, "RUNNING")
        if not take(p, me):
            for q in taken:
                release(q, me)
            with open(p) as fh:
                other = fh.read().split()
            raise RuntimeError("%s is already being run by %s (process %s); if that is wrong, delete %s" % (
                o, other[0] if other else "?", other[1] if len(other) > 1 else "?", p))
        taken.append(p)
    return [(p, me) for p in taken]


def _unlock(taken):
    for p, me in taken:
        release(p, me)


def _log(msg):
    print(msg, flush=True)


def output_folders(path):
    """One output folder, or every gene folder listed in a batch folder's batch.json."""
    if os.path.exists(os.path.join(path, "batch.json")):
        with open(os.path.join(path, "batch.json")) as fh:
            genes = json.load(fh)["genes"]
        return [os.path.join(path, g) for g in genes if os.path.exists(os.path.join(path, g, "setup.json"))]
    return [path]


class _Gene(object):
    """Launch queue of one output folder."""

    def __init__(self, out):
        from .report import load_manifest
        self.out = out
        with open(os.path.join(out, "setup.json")) as fh:
            self.meta = json.load(fh)
        rows = load_manifest(out)
        self.h = self.meta.get("sequential_h")
        d = lambda r: os.path.join(out, "runs", r["run"])
        self.obs = [d(r) for r in rows if r["label"] == "observed"]
        alts = [r for r in rows if r["label"] != "observed"]
        if self.h:
            alts.sort(key=lambda r: (int(r["order"]), r["run"]))
        self.alts = [d(r) for r in alts]
        self.stopped = False
        self.blocked = False

    def _done(self, d):
        return os.path.exists(os.path.join(d, "DONE"))

    def next_runs(self, running):
        """Folders of this gene that may start now."""
        todo_obs = [d for d in self.obs if not self._done(d) and d not in running]
        if not self.h:
            return todo_obs + [d for d in self.alts if not self._done(d) and d not in running]
        if todo_obs or any(d in running for d in self.obs):
            return todo_obs
        if self.stopped or self.blocked:
            return []
        state = self.sequential_state()
        if state is None:
            self.blocked = True
            return []
        if state.get("complete") and state.get("stopped"):
            self.stopped = True
            return []
        return [d for d in self.alts if not self._done(d) and d not in running]

    def sequential_state(self):
        from .report import build, observed_statistic, sequential_reading
        meta, rows, labels, excluded = build(self.out)
        g = observed_statistic(meta, rows, labels, meta.get("alpha", 0.05), meta.get("direction"))
        if not g["evaluable"]:
            return None
        return sequential_reading(meta, rows, labels, g["lobs_raw"])


def run_all(path, jobs, timeout=None, log=_log):
    folders = output_folders(path)
    taken = _lock(sorted(set(folders + [path])))
    _stopping.clear()

    def stop(signum, frame):
        _stopping.set()
        with _live_lock:
            procs = list(_live.values())
        for p in procs:
            _kill_group(p)

    old = {}
    if threading.current_thread() is threading.main_thread():
        for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            old[s] = signal.signal(s, stop)
    try:
        return _run_all(path, folders, jobs, timeout, log)
    finally:
        for s, h in old.items():
            signal.signal(s, h)
        _unlock(taken)


def _run_all(path, folders, jobs, timeout, log):
    genes = [_Gene(o) for o in folders]
    total = sum(len(g.obs) + len(g.alts) for g in genes)
    left = sum(1 for g in genes for d in g.obs + g.alts if not os.path.exists(os.path.join(d, "DONE")))
    log("%d output folder(s), %d run folders, %d to run, %d at a time%s" % (
        len(genes), total, left, jobs, "; sequential stopping on" if any(g.h for g in genes) else ""))
    running = {}
    busy = set()   # held by another process (RUN_LOCK): not started by this one
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        while True:
            # fill free slots, observed runs first across genes so every gate is known early
            cand = []
            for g in genes:
                nxt = g.next_runs(running)
                cand += [(0 if d in g.obs else 1, i, d) for i, d in enumerate(nxt)]
            cand.sort()
            for _, _, d in cand:
                if len(running) >= jobs or _stopping.is_set():
                    break
                if d not in running and d not in busy:
                    running[d] = ex.submit(run_one, d, timeout)
            if not running:
                break
            done, _ = wait(list(running.values()), return_when=FIRST_COMPLETED)
            for f in done:
                st, d = f.result()
                running.pop(d, None)
                if st == "busy":
                    busy.add(d)
                    continue
                log("%s %s %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), st,
                                  os.path.relpath(d, path) if d.startswith(path) else d))
    if busy:
        log("%d run(s) were left alone because another process holds their RUN_LOCK, e.g. %s; if nothing else is "
            "running them, delete those RUN_LOCK files" % (len(busy), os.path.relpath(sorted(busy)[0], path)))
    if _stopping.is_set():
        log("interrupted: the runs in progress were stopped and left without a DONE marker; 'labelrank run' starts "
            "them again from scratch")
        return -1
    for g in genes:
        if g.blocked:
            log("%s: the observed label is not evaluable (a run failed or was excluded); alternatives not started"
                % g.out)
        elif g.stopped:
            log("%s: sequential stopping reached; remaining alternatives not run" % g.out)
    failed = []
    for g in genes:
        for d in g.obs + g.alts:
            p = os.path.join(d, "DONE")
            if os.path.exists(p):
                with open(p) as fh:
                    if not fh.read().startswith("ok"):
                        failed.append(d)
    if failed:
        log("%d run(s) failed, e.g. %s (see its relax.out / errors.log or alt|null/codeml.log). They are excluded "
            "from the report. To retry one, delete its DONE file and run 'labelrank run' again (its old outputs are "
            "cleared when it starts)."
            % (len(failed), os.path.relpath(failed[0], path)))
    if busy:   # not finished: no ALL_DONE, and a non-zero exit
        return len(failed) + len(busy)
    with open(os.path.join(path, "ALL_DONE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
    return len(failed)
