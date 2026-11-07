"""PAML codeml branch-site model A (alternative) against its omega2 = 1 null.

Each run folder holds alt/ and null/ with their own control files: codeml writes fixed-name scratch files
(rst, rub, 2NG.*) into its working folder, so two runs sharing a folder corrupt each other without an error.
codeml has no random seed; a 'seed' here picks the starting omega and kappa, so several seeds play the same
role as HyPhy seeds: different starts, and the spread of logL among them shows whether the optimum was reached.
"""
import os
import re
import subprocess

TESTED = ("4.10.10",)   # releases whose mlc output this parser has been checked against
STARTS = [(0.5, 2.0), (1.5, 1.0), (3.0, 4.0), (0.2, 0.5), (5.0, 2.0), (0.8, 8.0)]   # (omega, kappa)

CTL = """      seqfile = aln.phy
     treefile = tree.nwk
      outfile = mlc
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 2
        icode = {icode}
    fix_kappa = 0
        kappa = {kappa}
    fix_omega = {fix_omega}
        omega = {omega}
    cleandata = 0
  fix_blength = 0
       method = 0
"""


def codeml_version(codeml):
    """'PAML X' from codeml's banner. codeml prints it only after reading a complete control file, so it is run
    on a stub whose sequence file does not exist, in a scratch folder (it stops there)."""
    import tempfile
    stub = ("seqfile = none.phy\ntreefile = none.nwk\noutfile = mlc\nnoisy = 0\nseqtype = 1\nmodel = 0\n"
            "NSsites = 0\nCodonFreq = 0\nicode = 0\nrunmode = 0\n")
    try:
        with tempfile.TemporaryDirectory() as td:
            with open(os.path.join(td, "v.ctl"), "w") as fh:
                fh.write(stub)
            out = subprocess.run([codeml, "v.ctl"], cwd=td, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, timeout=30).stdout
        m = re.search(r"paml version ([0-9][^,\s]*)", out.decode("utf-8", "replace"))
        return "PAML %s" % m.group(1) if m else "codeml (version not reported)"
    except (OSError, subprocess.SubprocessError):
        return "codeml (version not reported)"


def write_phylip(records, path):
    with open(path, "w") as fh:
        fh.write(" %d %d\n" % (len(records), len(records[0][1])))
        for n, s in records:
            fh.write("%s  %s\n" % (n, s.replace("U", "T").replace(".", "-").replace("?", "N")))


def write_run(d, records, tree_text, seed, codeml, icode=0):
    if not 1 <= seed <= len(STARTS):
        raise ValueError("codeml seeds pick one of %d starting points: use 1-%d" % (len(STARTS), len(STARTS)))
    omega, kappa = STARTS[(seed - 1) % len(STARTS)]
    os.makedirs(d)
    for sub, fix, om in (("alt", 0, omega), ("null", 1, 1.0)):
        sd = os.path.join(d, sub)
        os.makedirs(sd)
        write_phylip(records, os.path.join(sd, "aln.phy"))
        with open(os.path.join(sd, "tree.nwk"), "w") as fh:
            fh.write(tree_text)
        with open(os.path.join(sd, "codeml.ctl"), "w") as fh:
            fh.write(CTL.format(kappa=kappa, fix_omega=fix, omega=om, icode=icode))
    with open(os.path.join(d, "cmd.sh"), "w") as fh:
        q = "'" + codeml.replace("'", "'\\''") + "'"
        # codeml 4.10.10 exits with status 1 after writing a complete mlc (seen with strace: exit_group(1), no
        # message), so success is judged by both mlc files holding their lnL line, not by the exit status
        fh.write("cd \"$(dirname \"$0\")\"\n(cd alt && %s codeml.ctl > codeml.log 2>&1)\n"
                 "(cd null && %s codeml.ctl > codeml.log 2>&1)\n"
                 "grep -q '^lnL(ntime:' alt/mlc && grep -q '^lnL(ntime:' null/mlc\n" % (q, q))


_LNL = re.compile(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s*(-?[0-9.]+)")


def parse_mlc(text):
    m = _LNL.search(text)
    if not m:
        return None
    r = {"logL": float(m.group(3)), "np": int(m.group(2))}
    fw = re.search(r"foreground w\s+([0-9.\s]+)", text)
    if fw:
        vals = [float(x) for x in fw.group(1).split()]
        r["foreground_w2"] = vals[-1] if vals else None
    return r


def chi2_1_sf(x):
    import math
    return math.erfc(math.sqrt(max(x, 0.0) / 2.0))


def read_run_dir(d):
    out = {"fits": [], "warnings": [], "sets": {}, "deleted": frozenset(), "LRT": None, "p": None, "K": None,
           "best_logL": None, "null_logL": None, "reported_LRT": None, "reported_p": None, "reported_K": None}
    a = _read(os.path.join(d, "alt", "mlc"))
    n = _read(os.path.join(d, "null", "mlc"))
    pa = parse_mlc(a) if a else None
    pn = parse_mlc(n) if n else None
    if not pa or not pn:
        out["warnings"].append("incomplete: mlc missing or without lnL")
        return out
    if pa["np"] != pn["np"] + 1:
        out["warnings"].append("parameter counts %d (alt) and %d (null) do not differ by one" % (pa["np"], pn["np"]))
    lrt = 2.0 * (pa["logL"] - pn["logL"])
    if lrt < -1e-3:
        out["warnings"].append("alternative logL below null by %.4f: rerun with other starts" % (-lrt / 2))
    vm = re.search(r"paml version (\S+?),", a)
    out["engine_version"] = "PAML %s" % vm.group(1) if vm else None
    out.update(LRT=max(lrt, 0.0), p=0.5 * chi2_1_sf(lrt) if lrt > 0 else 1.0, K=pa.get("foreground_w2"),
               best_logL=pa["logL"], null_logL=pn["logL"], fits=[dict(pa, kind="alternative"), dict(pn, kind="null")],
               best_kind="alternative")
    return out


def _read(p):
    if not os.path.exists(p):
        return None
    with open(p, "rb") as fh:
        return fh.read().decode("utf-8", "replace")
