"""HyPhy RELAX: run folders, and reading what a run printed.

Reading is by log-likelihood, not by the reported verdict. When RELAX finds a flat likelihood surface it refits K
from the other side of 1 ("Attempt to fit an alternative direction of K"), so one output can hold two alternative
fits. In at least one published output the reported P and direction came from the lower of the two. This module
reads both fits and the null, and computes LRT = 2 (max alternative logL - null logL) itself.
"""
import json
import math
import os
import re
import shutil
import subprocess

SEED_LINE = "SetParameter (RANDOM_SEED, %d, 0);\n"

# ---------- locating HyPhy ----------


TESTED = ("2.5.101",)   # releases (with RELAX.bf 4.7) whose output this parser has been checked against

def find_relax_bf(hyphy):
    exe = shutil.which(hyphy) or hyphy
    exe = os.path.realpath(exe)
    prefix = os.path.dirname(os.path.dirname(exe))
    cand = os.path.join(prefix, "share", "hyphy", "TemplateBatchFiles", "SelectionAnalyses", "RELAX.bf")
    if not os.path.exists(cand):
        raise FileNotFoundError("RELAX.bf not found next to %s; pass --relax-bf" % exe)
    return cand


def hyphy_version(hyphy):
    try:
        out = subprocess.run([hyphy, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=60).stdout.decode("utf-8", "replace")
    except (OSError, subprocess.TimeoutExpired) as e:
        return "unknown (%s)" % e
    return out.strip().splitlines()[0] if out.strip() else "unknown"


# ---------- writing a run ----------


def write_run(d, aln_src, tree_text, seed, test, refarg, relax_bf, hyphy, code="Universal"):
    os.makedirs(d)
    shutil.copyfile(aln_src, os.path.join(d, "aln.fasta"))   # byte copy: the alignment is never rewritten
    with open(os.path.join(d, "tree.nwk"), "w") as fh:
        fh.write(tree_text)
    with open(relax_bf) as fh:
        src = fh.read()
    if src.startswith("SetParameter"):
        raise ValueError("%s already starts with a SetParameter line" % relax_bf)
    # the seeded copy must sit where RELAX.bf's relative includes resolve
    with open(os.path.join(d, "RELAX_seeded.bf"), "w") as fh:
        fh.write(SEED_LINE % seed + src)
    mod = os.path.join(os.path.dirname(relax_bf), "modules")
    os.symlink(mod, os.path.join(d, "modules"))
    with open(os.path.join(d, "cmd.sh"), "w") as fh:
        fh.write("cd \"$(dirname \"$0\")\" && %s CPU=1 RELAX_seeded.bf --alignment aln.fasta --tree tree.nwk "
                 "--code %s --test %s --reference %s --output relax.json > relax.out 2>&1\n" % (
                     _sh(hyphy), _sh(code), _sh(test), _sh(refarg)))


def _sh(x):
    return "'" + x.replace("'", "'\\''") + "'"


# ---------- tolerant JSON ----------


def loads_tolerant(text):
    """json.loads that accepts the bare inf, -inf, nan, -nan HyPhy writes. Strings are left untouched."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if c == '"':
            j = i + 1
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == '"':
                    break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        m = re.match(r"(-?)(inf(?:inity)?|nan)\b", text[i:i + 10], re.I)
        if m and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] in "._")):
            out.append("NaN" if m.group(2).lower() == "nan" else m.group(1) + "Infinity")
            i += m.end()
            continue
        out.append(c)
        i += 1
    return json.loads("".join(out))


# ---------- reading relax.out ----------

_SEL = re.compile(r"Selected (\d+) branches as the _([A-Za-z]+)_ set: `([^`]*)`")
_NUI = re.compile(r"(\d+) branches are in the unclassified \(nuisance\) set: `([^`]*)`")
_DEL = re.compile(r"Deleted (\d+) zero-length internal branch(?:es)?: `([^`]*)`")
_LOGL = re.compile(r"Log\(L\)\s*=\s*(-?[0-9.]+(?:[eE][+-]?\d+)?|-?inf|nan)")
_K = re.compile(r"Relaxation/intensification parameter \(K\)\s*=\s*(-?[0-9.]+(?:[eE][+-]?\d+)?|inf|nan)")
_P = re.compile(r"Likelihood ratio test \*\*p\s*=\s*([0-9.]+(?:[eE][+-]?\d+)?)\*\*")


def _names(s):
    return frozenset(x.strip() for x in s.split(",") if x.strip())


def parse_out(text):
    r = {"sets": {}, "deleted": frozenset(), "fits": [], "reported_p": None, "reported_evidence": None,
         "reset_note": "Relaxation parameter reset to opposite mode" in text,
         # after a negative LRT RELAX refits up to 5 more times; those passes print no section headers, so the
         # fits below are from the abandoned first pass and only the JSON holds the final ones
         "refit": "Detected convergence issues (negative LRT)" in text}
    for m in _SEL.finditer(text):
        names = _names(m.group(3))
        if len(names) != int(m.group(1)):
            raise ValueError("branch-set line lists %d names but says %s" % (len(names), m.group(1)))
        r["sets"][m.group(2)] = names
    m = _NUI.search(text)
    if m:
        r["sets"]["nuisance"] = _names(m.group(2))
    elif "Test" in r["sets"] and "Reference" in r["sets"]:
        r["sets"]["nuisance"] = frozenset()   # RELAX prints this line only when the set is not empty
    dels = set()
    for m in _DEL.finditer(text):
        dels |= _names(m.group(2))
    r["deleted"] = frozenset(dels)
    # fits: walk the section headers in order and take the first Log(L) / K after each
    heads = [(m.start(), m.group(0)) for m in re.finditer(r"^#{2,3} .*$", text, re.M)]
    for i, (pos, h) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else len(text)
        body = text[pos:end]
        if "Fitting the alternative model" in h:
            kind = "alternative"
        elif "Attempt to fit an alternative direction of K" in h:
            kind = "alternative_other_direction"
        elif "Fitting the null (K := 1) model" in h:
            kind = "null"
        else:
            continue
        lm, km = _LOGL.search(body), _K.search(body)
        r["fits"].append({"kind": kind, "logL": _f(lm.group(1)) if lm else None,
                          "K": _f(km.group(1)) if km else None})
    m = _P.search(text)
    if m:
        r["reported_p"] = float(m.group(1))
    m = re.search(r">(No significant evidence|Evidence) for \*?(relaxation|intensification)", text)
    if m:
        r["reported_evidence"] = "none" if m.group(1).startswith("No") else m.group(2)
    return r


def _f(x):
    x = x.lower()
    return float("nan") if x == "nan" else float(x)


# ---------- reading relax.json ----------


def parse_json(text):
    j = loads_tolerant(text)
    if not isinstance(j, dict):
        raise ValueError("relax.json is not a JSON object")
    fits = j.get("fits") or {}
    tr = j.get("test results") or {}
    st = (j.get("analysis") or {}).get("settings") or {}
    alt = fits.get("RELAX alternative") or {}
    null = fits.get("RELAX null") or {}
    k = tr.get("relaxation or intensification parameter")
    if k is None:
        k = (alt.get("Rate Distributions") or {}).get("relaxation or intensification parameter")
    take2 = st.get("convergence-unstable-alernative") or st.get("convergence-unstable-alternative")
    t2 = None
    if take2:
        row = take2[0] if isinstance(take2[0], list) else take2
        t2 = {"K": float(row[0]), "logL": float(row[1])}
    # written when the other-direction fit won: [[K, logL] of the winner, [K, logL] of the first fit]
    unstable = [{"K": float(x[0]), "logL": float(x[1])} for x in (st.get("convergence-unstable") or [])
                if isinstance(x, list) and len(x) == 2]
    return {"negative_lrt_passes": st.get("convergence-negative-lrt"), "unstable": unstable,"alt_logL": _num(alt.get("Log Likelihood")), "null_logL": _num(null.get("Log Likelihood")),
            "reported_LRT": _num(tr.get("LRT")), "reported_p": _num(tr.get("p-value")), "reported_K": _num(k),
            "other_direction": t2, "version": (j.get("analysis") or {}).get("version")}


def _num(x):
    return None if x is None else float(x)


# ---------- one run, read by logL ----------


def chi2_1_sf(x):
    """Upper tail of chi-square with 1 df (RELAX's reference distribution)."""
    if x is None or x != x:
        return float("nan")
    return math.erfc(math.sqrt(max(x, 0.0) / 2.0))


def _finite(x):
    return x is not None and x == x and abs(x) != float("inf")


def side(k):
    """'relax' (K < 1), 'intens' (K > 1), or None when K is 1 within print rounding (relax.out gives 2 decimals)."""
    if k is None or k != k or abs(k - 1.0) < 0.006:
        return None
    return "relax" if k < 1 else "intens"


def read_result(out_text=None, json_text=None):
    """Combine relax.out and relax.json (either may be missing) into one reading.

    Returns dict: LRT, p, K (from the higher-logL alternative fit), the reported values, every fit seen, and
    warnings. JSON values are preferred (full precision); relax.out fills in when JSON is absent.
    """
    o = parse_out(out_text) if out_text else None
    js = parse_json(json_text) if json_text else None
    fits, null, notes = [], None, []
    if o and o["refit"]:
        if not js:
            return _unreadable(o, "RELAX refitted after a negative LRT and relax.json is missing; relax.out alone "
                                  "holds only the abandoned first pass")
        notes.append("RELAX refitted after a negative LRT (%s passes); read from relax.json" % js["negative_lrt_passes"])
        o_fits = []
    else:
        o_fits = o["fits"] if o else []
    if o_fits:
        # relax.out pairs each fit's logL with its K (2 decimals); a JSON value replaces it only when both match
        precise = []
        if js:
            od = js["other_direction"]
            precise = [(js["alt_logL"], js["reported_K"])] + ([(od["logL"], od["K"])] if od else []) + [
                (u["logL"], u["K"]) for u in js["unstable"]]
            precise = [(x, k) for x, k in precise if x is not None]
        used = set()
        for f in o_fits:
            if f["logL"] is None:
                continue
            g = dict(f, src="out")
            if f["kind"] == "null":
                if null is None:
                    null = g["logL"]
                    if js and js["null_logL"] is not None and abs(js["null_logL"] - null) <= 0.006:
                        null = js["null_logL"]
                continue
            for i, (x, kx) in enumerate(precise):
                if i in used or abs(x - f["logL"]) > 0.006:
                    continue
                if kx is not None and f["K"] is not None and abs(kx - f["K"]) > 0.006:
                    continue
                g["logL"], g["src"] = x, "json"
                if kx is not None:
                    g["K"] = kx
                used.add(i)
                break
            fits.append(g)
        if js and js["other_direction"] and not any(f["kind"] == "alternative_other_direction" for f in fits):
            fits.append(dict(js["other_direction"], kind="alternative_other_direction", src="json"))
    elif js and js["alt_logL"] is not None:
        fits.append({"kind": "alternative", "logL": js["alt_logL"], "K": js["reported_K"], "src": "json"})
        od = js["other_direction"]
        if od and not (abs(od["logL"] - js["alt_logL"]) < 1e-9 and od["K"] == js["reported_K"]):
            fits.append(dict(od, kind="alternative_other_direction", src="json"))
        for u in js["unstable"]:
            if not any(abs(u["logL"] - f["logL"]) < 1e-9 for f in fits):
                fits.append(dict(u, kind="alternative_first_fit", src="json"))
        null = js["null_logL"]
    warnings = notes
    res = {"fits": fits, "null_logL": null, "warnings": warnings,
           "reported_LRT": js["reported_LRT"] if js else None,
           "reported_p": js["reported_p"] if js else (o["reported_p"] if o else None),
           "reported_K": js["reported_K"] if js else None,
           "reported_evidence": o["reported_evidence"] if o else None,
           "sets": o["sets"] if o else {}, "deleted": o["deleted"] if o else frozenset()}
    if not fits or null is None:
        res.update(LRT=None, p=None, K=None, best_logL=None)
        warnings.append("incomplete: no alternative or null fit found")
        return res
    finite = [f for f in fits if _finite(f["logL"])]
    if not finite or not _finite(null):
        # a non-finite logL would turn into an LRT of nan or inf that compares silently in the rank
        return _unreadable(o, "non-finite log-likelihood (alternative %s, null %s)" % (
            [f["logL"] for f in fits], null))
    best = max(finite, key=lambda f: f["logL"])
    lrt = 2.0 * (best["logL"] - null)
    res.update(LRT=lrt, p=chi2_1_sf(lrt), K=best["K"], best_logL=best["logL"], best_kind=best["kind"])
    if lrt < -1e-3:
        warnings.append("alternative logL below null logL by %.4f: the optimizer did not converge" % (-lrt / 2))
    if len(finite) > 1:
        lo = min(finite, key=lambda f: f["logL"])
        if side(lo["K"]) and side(best["K"]) and side(lo["K"]) != side(best["K"]):
            warnings.append("two alternative fits on opposite sides of K=1 (logL %.3f with K %.3g; %.3f with K %.3g);"
                            " the higher one is used" % (best["logL"], best["K"], lo["logL"], lo["K"]))
    rep_lrt = res["reported_LRT"]
    if rep_lrt is not None and abs(rep_lrt - lrt) > 0.01:
        warnings.append("reported LRT %.4f differs from 2*(max alternative logL - null logL) = %.4f" % (rep_lrt, lrt))
    if side(res["reported_K"]) and side(best["K"]) and side(res["reported_K"]) != side(best["K"]):
        warnings.append("reported K %.3g is on the other side of 1 from the higher-logL fit (K %.3g)" % (
            res["reported_K"], best["K"]))
    ev = res["reported_evidence"]
    if ev in ("relaxation", "intensification") and side(best["K"]) and (
            (ev == "relaxation") != (side(best["K"]) == "relax")):
        warnings.append("reported direction '%s' disagrees with the higher-logL fit (K %.3g)" % (ev, best["K"]))
    if res["reported_p"] is not None and abs(res["reported_p"] - res["p"]) > max(1e-3, 0.05 * res["p"]):
        warnings.append("reported P %.4g differs from the P of the higher-logL fit %.4g" % (res["reported_p"], res["p"]))
    if best["K"] is not None and best["K"] >= 49.99:
        warnings.append("K at its upper bound (50)")
    return res


def _unreadable(o, why):
    return {"LRT": None, "p": None, "K": None, "fits": [], "sets": o["sets"] if o else {},
            "deleted": o["deleted"] if o else frozenset(), "warnings": [why], "null_logL": None, "best_logL": None,
            "reported_LRT": None, "reported_p": None, "reported_K": None, "reported_evidence": None, "unreadable": why}


def read_run_dir(d):
    out = _read(os.path.join(d, "relax.out"))
    js = _read(os.path.join(d, "relax.json"))
    try:
        return read_result(out, js)
    except Exception as e:   # a JSON or output that cannot be read is reported, never dropped silently
        return _unreadable(None, "unreadable output: %s: %s" % (type(e).__name__, e))


def _read(p):
    if not os.path.exists(p):
        return None
    with open(p, "rb") as fh:
        return fh.read().decode("utf-8", "replace")
