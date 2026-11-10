"""Where the observed label ranks among same-shaped labels on the same tree.

Rules (fixed here, not tuned to any result):
- Instrument check per run, by branch name: the Test, Reference and nuisance sets the program printed equal the
  sets this tool wrote. A run that fails the check, failed, or cannot be read is excluded and listed; a run without
  a DONE marker is 'not run' and never read. A run in which the program deleted a test branch as zero-length is
  excluded too (its test set changed); deleted reference branches are counted and reported.
- Every fit is read by log-likelihood (hyphy_relax.read_result): the higher of RELAX's two alternative fits is used.
- Label statistic, the same for every label: seeds are optimizer restarts, so LRT = 2 (best alternative logL over
  the label's seeds - best null logL over its seeds). Every label is set up with the same seeds by default.
- Gate: the observed label's LRT gives P <= alpha (with a direction, also K of its best fit on that side of 1).
- Co-primary readings, non-directional: rank = 1 + number of alternatives with LRT >= observed, out of N + 1, over
  (a) all alternatives and (b) the length-matched ones (test-set length within the setup's ratio of the observed
  one). A rank-based P = rank / (N + 1) is given for a reading only when its N + 1 >= 20 and every written
  alternative has finished.
- Secondary: directional (RELAX: alternatives at or above the observed LRT with K on the observed side), and the
  share of alternatives the chi-square reading calls significant.
- Length: the share of alternatives shorter than the observed test set is always printed; short lineages carry few
  substitutions and give small LRTs whatever their biology, so a high share is flagged.
- Optimizer check: for every label with more than one seed, the spread of the alternative and of the null logL
  across seeds; a spread above --logl-tolerance is flagged.
"""
import json
import os
import statistics

LRT_FLOOR = 0.01   # LRTs below this are optimizer noise around 0 and read as 0


def at_or_above(v, lobs_raw):
    """An alternative counts as at or above the observed label when its LRT is within LRT_FLOOR of it: noise of
    that size cannot order two labels, at 0 or anywhere else, and a tie counts against the observed label."""
    return max(v["lrt_raw"], 0.0) >= lobs_raw - LRT_FLOOR


def _side(k):
    if k is None or k != k or abs(k - 1.0) < 0.006:
        return None
    return "relax" if k < 1 else "intens"


def _jload(p):
    with open(p) as fh:
        return json.load(fh)


def _text(p):
    with open(p) as fh:
        return fh.read()


def load_manifest(out):
    with open(os.path.join(out, "runs.tsv")) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh]
    head = lines[0]
    return [dict(zip(head, r)) for r in lines[1:]]


def read_run(out, engine, run):
    d = os.path.join(out, "runs", run)
    if engine == "relax":
        from .hyphy_relax import read_run_dir
    else:
        from .codeml import read_run_dir
    res = read_run_dir(d)
    exp_p = os.path.join(d, "expected_sets.json")
    exp = _jload(exp_p) if os.path.exists(exp_p) else None
    problem = None
    done = _text(os.path.join(d, "DONE")).split("\t")[0] if os.path.exists(os.path.join(d, "DONE")) else None
    if done is None:
        # no DONE marker: not started or still running; partial output is never read
        problem = "not run"
    elif done != "ok":
        # the run's own exit status says it failed: its output, even if parseable, is not used
        problem = "the program failed (status %s; see the run's log files)" % done
    elif res.get("unreadable"):
        problem = res["unreadable"]
    elif res.get("LRT") is None:
        problem = "incomplete (status %s)" % done
    elif exp is None:
        problem = "expected_sets.json missing: the run cannot be checked"
    elif engine == "relax":
        got = res.get("sets") or {}
        if not got or any(k not in got for k in ("Test", "Reference", "nuisance")):
            problem = "branch sets not printed"
        else:
            shown = set().union(*[set(got[k]) for k in ("Test", "Reference", "nuisance")])
            unshown = set().union(*[set(exp[k]) for k in ("Test", "Reference", "nuisance")]) - shown
            if len(unshown) == 1 and unshown <= set(exp.get("root_children", [])):
                # the root branch HyPhy merged away; whatever set it was in, it is not compared
                merged = sorted(unshown)[0]
                if merged in exp["Test"]:
                    problem = "test branch %s was merged away at the root by HyPhy (test set changed)" % merged
                exp = {k: sorted(set(v) - unshown) if isinstance(v, list) else v for k, v in exp.items()}
                res["merged_root_branch"] = merged
            for k in ("Test", "Reference", "nuisance"):
                if problem:
                    break
                if set(got[k]) != set(exp[k]):
                    a, b = set(got[k]) - set(exp[k]), set(exp[k]) - set(got[k])
                    problem = "%s set differs from the one written (+%d, -%d; e.g. %s)" % (
                        k, len(a), len(b), sorted(a | b)[0])
                    break
        if problem is None and set(res.get("deleted") or ()) & set(exp["Test"]):
            problem = "test branch(es) deleted as zero-length: %s" % ", ".join(
                sorted(set(res["deleted"]) & set(exp["Test"])))
    res["problem"] = problem
    res["deleted_reference"] = sorted(set(res.get("deleted") or ()) & set(exp["Reference"])) if exp else []
    return res


def build(out):
    meta = _jload(os.path.join(out, "setup.json"))
    engine = meta["engine"]
    rows = load_manifest(out)
    per = {}
    excluded = []
    for r in rows:
        res = read_run(out, engine, r["run"])
        r["res"] = res
        if res["problem"]:
            if res["problem"] != "not run":
                excluded.append((r["run"], res["problem"]))
            continue
        per.setdefault(r["label"], []).append(r)
    labels = {}
    written = {}
    for r in rows:
        written[r["label"]] = written.get(r["label"], 0) + 1
    for lab, rs in per.items():
        bl = [x["res"]["best_logL"] for x in rs]
        nl = [x["res"]["null_logL"] for x in rs]
        sp_a = (max(bl) - min(bl)) if len(bl) > 1 else None
        sp_n = (max(nl) - min(nl)) if len(nl) > 1 else None
        # seeds are optimizer restarts of one analysis: each model's best fit over the seeds is its estimate,
        # and the LRT is taken between those (the same seeds are run for every label)
        best = max(rs, key=lambda x: x["res"]["best_logL"])
        lrt_raw = 2.0 * (best["res"]["best_logL"] - max(nl))
        # below LRT_FLOOR (negative included: an alternative fit below its own null did not converge) the LRT is
        # optimizer noise around 0: it reads as 0 everywhere, so noise ties instead of ordering labels, and a tie
        # counts against the observed label
        lrt = lrt_raw if lrt_raw >= LRT_FLOOR else 0.0
        labels[lab] = {"runs": rs, "lrt": lrt, "lrt_raw": lrt_raw, "k": best["res"]["K"],
                       "seeds_used": len(rs), "seeds_written": written[lab],
                       # both fits move with the optimizer: the spread is the larger of the two across seeds
                       "logl_spread": max(x for x in (sp_a, sp_n, -1.0) if x is not None) if len(rs) > 1 else None,
                       "logl_spread_alt": sp_a, "logl_spread_null": sp_n,
                       "length": float(rs[0]["length"]), "length_matched": rs[0]["length_matched"] == "yes",
                       "branches": rs[0]["branches"]}
    return meta, rows, labels, excluded


def pvalue_fn(engine, model=None):
    from .hyphy_relax import chi2_1_sf
    if engine == "relax" or model == "two-ratio":
        # RELAX, and codeml two-ratio (nested at the interior point omega_fg = omega_bg): chi-square(1)
        return lambda x: chi2_1_sf(x) if x > 0 else 1.0
    # branch-site model A: 50:50 mixture of a point mass at 0 and chi-square(1)
    return lambda x: 0.5 * chi2_1_sf(x) if x > 0 else 1.0


def observed_statistic(meta, rows, labels, alpha=0.05, direction=None):
    """Gate and observed LRT (also used by the sequential runner).

    The observed label is evaluable when every one of its seeds passed the instrument check. Its LRT is the
    label statistic (best alternative fit against best null fit over the seeds, as for every label). Gate: that
    LRT gives P <= alpha (and, with a direction, K of the best fit on that side of 1).
    """
    engine = meta["engine"]
    obs = [r for r in rows if r["label"] == "observed"]
    obs_ok = [r for r in obs if not r["res"]["problem"]]
    n_seed = len(meta["seeds"])
    g = {"obs": obs, "obs_ok": obs_ok, "n_seed": n_seed, "evaluable": len(obs_ok) >= n_seed and "observed" in labels}
    if not g["evaluable"]:
        return g
    ob = labels["observed"]
    pval = pvalue_fn(engine, meta.get("codeml_model"))
    side = _side(ob["k"]) if engine == "relax" else None
    gate = pval(ob["lrt"]) <= alpha and (not direction or side == direction)
    seed_pass = [r for r in obs_ok if r["res"]["p"] <= alpha]
    g.update(gate=gate, lobs=ob["lrt"], lobs_raw=ob["lrt_raw"], p_obs=pval(ob["lrt"]), obs_side=side,
             passing=seed_pass)
    return g


def sequential_reading(meta, rows, labels, lobs_raw):
    """Besag-Clifford sequential reading over alternatives in their pre-drawn order.

    Walk the order label by label. A label with any run not finished (no DONE) blocks the walk: the reading is
    incomplete. A label whose finished runs were all excluded is not a draw and is skipped. Otherwise the label is
    a draw, read by its label statistic (best fits over its usable seeds). Stop at the first L where h draws are at or above the observed
    LRT: P = h / L. If the order runs out first: P = (g + 1) / (n + 1). Labels past the stop are not read.
    """
    h = meta.get("sequential_h")
    runs_of = {}
    order = {}
    for r in rows:
        if r["label"] == "observed":
            continue
        runs_of.setdefault(r["label"], []).append(r)
        if r.get("order"):
            order[r["label"]] = int(r["order"])
    g, read = 0, []
    for lab in sorted(order, key=lambda x: (order[x], x)):
        if any(r["res"]["problem"] == "not run" for r in runs_of[lab]):
            return {"h": h, "complete": False, "read": len(read), "ge": g, "labels": read}
        if lab not in labels:
            continue
        read.append(lab)
        if at_or_above(labels[lab], lobs_raw):
            g += 1
            if g >= h:
                return {"h": h, "complete": True, "stopped": True, "read": len(read), "ge": g, "labels": read,
                        "P": float(h) / len(read)}
    return {"h": h, "complete": True, "stopped": False, "read": len(read), "ge": g, "labels": read,
            "P": (g + 1.0) / (len(read) + 1)}


def _ranks(xs):
    """Average ranks (ties share the mean of their positions)."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def spearman(xs, ys):
    """Spearman's rank correlation, or None when it is undefined (fewer than 3 pairs, or a constant)."""
    if len(xs) < 3:
        return None
    a, b = _ranks(xs), _ranks(ys)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    sab = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    saa = sum((x - ma) ** 2 for x in a)
    sbb = sum((y - mb) ** 2 for y in b)
    if saa == 0 or sbb == 0:
        return None
    return sab / (saa * sbb) ** 0.5


def length_rho(lengths, lrts):
    """(Spearman rho of length against LRT, number of LRTs at or above LRT_FLOOR). LRTs below the floor are
    optimizer noise around 0 (codeml leaves differences of 1e-6); ranked as they are they would turn noise into a
    correlation, so they tie at 0, and rho is None unless 3 alternatives clear the floor."""
    floored = [x if x >= LRT_FLOOR else 0.0 for x in lrts]
    above = sum(1 for x in floored if x > 0)
    return (spearman(lengths, floored) if above >= 3 else None), above


def summarize(out, alpha=None, direction=None, logl_tolerance=0.5):
    """Every number the report, the methods paragraph and the figure use.

    alpha and direction default to the values fixed at setup (0.05 and none when setup fixed nothing).
    """
    meta, rows, labels, excluded = build(out)
    alpha = alpha if alpha is not None else meta.get("alpha", 0.05)
    direction = direction if direction is not None else meta.get("direction")
    direction = None if direction == "none" or meta["engine"] != "relax" else direction
    engine = meta["engine"]
    pval = pvalue_fn(engine, meta.get("codeml_model"))
    S = {"meta": meta, "rows": rows, "labels": labels, "excluded": excluded, "alpha": alpha,
         "direction": direction, "logl_tolerance": logl_tolerance, "pval": pval}
    g = observed_statistic(meta, rows, labels, alpha, direction)
    S.update(g)
    S["engine_version_read"] = next((r["res"].get("engine_version") for r in rows
                                     if r["res"].get("engine_version")), None)
    if engine == "relax":
        from .hyphy_relax import TESTED
    else:
        from .codeml import TESTED
    ver = S["engine_version_read"] or meta.get("engine_version") or ""
    from .gencode import release
    S["untested_version"] = None if release(ver) in TESTED else (ver or "unknown", TESTED)
    if not g["evaluable"]:
        return S
    lobs = g["lobs"]
    alts = [(k, v) for k, v in labels.items() if k != "observed"]
    if meta.get("sequential_h"):
        sq = sequential_reading(meta, rows, labels, g["lobs_raw"])
        if meta["alternatives_written"] + 1 < 20:
            sq["P"] = None
        S["sequential"] = sq
        # every reading uses the alternatives walked in the drawn order, never runs that finished past the stop
        keep = set(sq["labels"])
        alts = [(k, v) for k, v in alts if k in keep]
    N = len(alts)
    ge = [(k, v) for k, v in alts if at_or_above(v, g["lobs_raw"])]
    S["not_run"] = [r["run"] for r in rows if r["res"]["problem"] == "not run"]
    # outside sequential mode, an unfinished alternative leaves the rank provisional: no P is given
    S["incomplete"] = (not meta.get("sequential_h")) and any(
        r["label"] != "observed" for r in rows if r["res"]["problem"] == "not run")
    # still running in either mode: N = 0 then means 'nothing read yet', not 'every run failed'
    S["pending"] = S["incomplete"] or bool(S.get("sequential") and not S["sequential"]["complete"])
    rank_p = ((1.0 + len(ge)) / (N + 1) if N + 1 >= 20 and not meta.get("sequential_h") and not S["incomplete"]
              else None)
    S.update(alts=alts, N=N, ge=ge, rank=1 + len(ge), rank_p=rank_p)
    if engine == "relax":
        S["directional"] = [k for k, v in ge if g["obs_side"] and _side(v["k"]) == g["obs_side"]]
    lm = [(k, v) for k, v in alts if v["length_matched"]]
    S["length_matched"] = lm
    S["length_matched_ge"] = [k for k, v in lm if at_or_above(v, g["lobs_raw"])]
    S["lm_rank"] = 1 + len(S["length_matched_ge"])
    S["lm_rank_p"] = ((1.0 + len(S["length_matched_ge"])) / (len(lm) + 1)
                      if len(lm) + 1 >= 20 and not meta.get("sequential_h") and not S["incomplete"] else None)
    olen = float("%.9g" % meta["observed_length"])   # runs.tsv carries lengths at this precision
    S["shorter_share"] = (sum(1 for k, v in alts if v["length"] < olen) / float(N)) if N else None
    S["length_rho"], S["above_floor"] = length_rho([v["length"] for k, v in alts], [v["lrt"] for k, v in alts])
    S["asymmetric_seeds"] = sorted(meta["seeds"]) != sorted(meta["alt_seeds"])
    S["partial_seeds"] = [k for k, v in labels.items() if v["seeds_used"] < v["seeds_written"]]
    S["significant"] = [(k, v) for k, v in alts if pval(v["lrt"]) <= alpha]
    S["below_null"] = sorted((k, v["lrt_raw"]) for k, v in labels.items() if v["lrt_raw"] < -LRT_FLOOR)
    S["spread"] = [(k, v["logl_spread"]) for k, v in alts
                   if v["logl_spread"] is not None and v["logl_spread"] > logl_tolerance]
    return S


def report(out, alpha=None, direction=None, logl_tolerance=0.5, S=None):
    S = S or summarize(out, alpha, direction, logl_tolerance)
    meta, rows, labels, excluded = S["meta"], S["rows"], S["labels"], S["excluded"]
    engine = meta["engine"]
    L = []
    p = L.append
    p("# labelrank report")
    p("engine: %s (%s); test label %r, %s; alternatives written %d of %s%s" % (
        engine + (" " + meta["codeml_model"] if meta.get("codeml_model") else ""),
        S.get("engine_version_read") or meta["engine_version"], meta["test"],
        "reference %r" % meta["reference"] if engine == "relax" else "background = every unmarked branch",
        meta["alternatives_written"], "%d found" % meta["alternatives_found"] if meta["alternatives_found"] is not None
        else "an unknown number (drawn at random)", " (sampled)" if meta["sampled"] else ""))
    if meta.get("genetic_code", 1) != 1:
        p("genetic code: NCBI table %d (%s)" % (meta["genetic_code"], meta.get("genetic_code_name")))
    if meta.get("exclude"):
        p("lineages kept out of the alternatives: %s; %d reference branch(es)" % (
            "; ".join(meta["exclude"]), meta.get("excluded_reference_branches", 0)))
    if S.get("untested_version"):
        p("WARNING: %s is not a release this tool was checked against (%s); check a few runs by hand" % (
            S["untested_version"][0], ", ".join(S["untested_version"][1])))
    for run, why in excluded:
        p("excluded: %s - %s" % (run, why))
    nr = [r["run"] for r in rows if r["res"]["problem"] == "not run"]
    if nr and (S.get("sequential") or {}).get("stopped"):
        p("not launched: %d run(s) past the sequential stopping point" % len(nr))
    elif nr:
        p("not run (no DONE marker): %d run(s)%s" % (len(nr), ", e.g. " + nr[0] if nr else ""))
    obs = S["obs"]
    p("")
    p("## observed label (%d seeds; %s)" % (len(obs), obs[0]["branches"] if obs else "?"))
    for r in obs:
        x = r["res"]
        if x["problem"]:
            continue
        p("  %s  LRT %.4f  P %.4g  K %.4g  best logL %.4f (%s)" % (
            r["run"], x["LRT"], x["p"], x["K"] if x["K"] is not None else float("nan"), x["best_logL"],
            x.get("best_kind", "")))
        for w in x["warnings"]:
            p("    warning: " + w)
        if x.get("deleted_reference"):
            p("    note: %d reference branch(es) deleted as zero-length (names as in the run's tree.nwk): %s" % (
                len(x["deleted_reference"]), ", ".join(x["deleted_reference"])))
    if not S["evaluable"]:
        nr = sum(1 for r in obs if r["res"]["problem"] == "not run")
        p("  gate not evaluable: %d of %d observed runs are usable (%d not run yet, %d excluded above)" % (
            len(S["obs_ok"]), S["n_seed"], nr, S["n_seed"] - len(S["obs_ok"]) - nr))
        return "\n".join(L) + "\n"
    logl_tolerance, alpha, direction, pval = S["logl_tolerance"], S["alpha"], S["direction"], S["pval"]
    ob = labels.get("observed")
    if ob["logl_spread"] is not None:
        p("  logL spread across seeds: alternative %.4f, null %.4f%s" % (
            ob["logl_spread_alt"] or 0.0, ob["logl_spread_null"] or 0.0, " - ABOVE TOLERANCE %.2f: the optimizer may "
                                                  "not have converged; add seeds" % logl_tolerance
                                                  if ob["logl_spread"] > logl_tolerance else ""))
    lobs, gate = S["lobs"], S["gate"]
    p("  observed LRT %.4f (best alternative vs best null fit over %d seeds), P %.4g%s" % (
        lobs, ob["seeds_used"], S["p_obs"], ", K %.4g" % ob["k"] if engine == "relax" and ob["k"] is not None else ""))
    p("  gate: P <= %g%s -> %s (seeds individually: %d of %d at P <= %g)" % (
        alpha, " with K on the %s side" % direction if direction else "",
        "PASSED" if gate else "NOT PASSED (the observed label is not significant here; the ranks below are still "
        "computed)", len(S["passing"]), S["n_seed"], alpha))
    if S["asymmetric_seeds"]:
        p("  warning: the observed label used seeds %s but alternatives used %s; the labels did not get the same "
          "optimizer restarts, so their LRTs are not read on equal terms (use 'labelrank add-seeds')" % (
              meta["seeds"], meta["alt_seeds"]))
    alts, N, ge = S["alts"], S["N"], S["ge"]
    nr_alt = sum(1 for r in rows if r["label"] != "observed" and r["res"]["problem"] == "not run")
    p("")
    p("## rank among same-shaped labels")
    if S.get("incomplete"):
        p("  INCOMPLETE: %d alternative run(s) have not finished; the ranks below are provisional and no P is given"
          % nr_alt)
    seqnote = "; in the drawn order up to the stopping point" if S.get("sequential") else ""
    if N == 0 and S.get("pending"):
        p("  NO RANK YET: no alternative has finished; run 'labelrank run' to continue")
    elif N == 0:
        p("  NO RANK: none of the %d alternatives written could be read (every run failed or was excluded; see "
          "below); the lines that follow have nothing to rank against" % meta["alternatives_written"])
    p("  (a) all alternatives: %d read (of %d written%s); LRT >= observed: %d; rank %d of %d%s" % (
        N, meta["alternatives_written"], seqnote, len(ge), 1 + len(ge), N + 1,
        "; rank-based P = %d/%d = %.4f" % (1 + len(ge), N + 1, S["rank_p"]) if S["rank_p"] is not None else ""))
    lm = S["length_matched"]
    p("  (b) length-matched (test-set length within x%g of %.5f): %d alternatives; LRT >= observed: %d; rank %d of "
      "%d%s" % (meta["length_ratio"], meta["observed_length"], len(lm), len(S["length_matched_ge"]), S["lm_rank"],
                len(lm) + 1, "; rank-based P = %d/%d = %.4f" % (S["lm_rank"], len(lm) + 1, S["lm_rank_p"])
                if S["lm_rank_p"] is not None else ""))
    if S.get("sequential"):
        p("  primary P: the sequential P below (a rank-based P over a stopped walk would be biased)")
    elif not S.get("incomplete") and (N + 1 < 20 or len(lm) + 1 < 20):
        p("  a reading with N + 1 < 20 gives its rank only; no P value is claimed for it")
    if S["shorter_share"] is not None:
        p("  length: %.0f%% of the alternatives are shorter than the observed test set%s" % (
            100 * S["shorter_share"], " - the observed label is longer than most; short lineages give small LRTs "
            "whatever their biology, so reading (a) favours it; weigh reading (b)" if S["shorter_share"] >= 0.75
            else ""))
    if S["length_rho"] is not None:
        p("  length vs LRT among the alternatives: Spearman rho %.2f (n = %d)%s" % (
            S["length_rho"], N, "; longer lineages give larger LRTs on this tree, so reading (b) matters"
            if S["length_rho"] >= 0.3 else ""))
    elif N >= 3 and S["above_floor"] < 3:
        p("  length vs LRT among the alternatives: not defined (%d of %d alternatives have LRT >= %g; the rest tie "
          "at 0)" % (S["above_floor"], N, LRT_FLOOR))
    elif N >= 3:
        p("  length vs LRT among the alternatives: not defined (their lengths, or their LRTs, are all equal)")
    p("  (ranks among this tree's lineages of the same shape; lineages differ in history and are not independent, "
      "so these are calibrations, not exact null distributions)")
    sq = S.get("sequential")
    if sq:
        if not sq["complete"]:
            p("  sequential reading (h = %d): incomplete; %d alternatives read in order, %d at or above" % (
                sq["h"], sq["read"], sq["ge"]))
        else:
            p("  sequential reading (Besag-Clifford, h = %d): %s after %d alternatives in the drawn order, %d at or "
              "above the observed LRT%s" % (
                  sq["h"], "stopped" if sq["stopped"] else "not stopped", sq["read"], sq["ge"],
                  "; sequential P = %.4f" % sq["P"] if sq["P"] is not None else "; rank only (N + 1 < 20)"))
    if engine == "relax":
        if S["obs_side"]:
            p("  secondary (directional): alternatives with LRT >= observed and K on the observed side (%s): %d of %d"
              % (S["obs_side"], len(S["directional"]), N))
        else:
            p("  secondary (directional): not applicable (the observed K is 1)")
    sig = S["significant"]
    p("  alternatives significant by the chi-square reading (P <= %g): %d of %d (%.1f%%)%s" % (
        alpha, len(sig), N, 100.0 * len(sig) / max(1, N),
        ", %d on the observed side" % sum(1 for k, v in sig if _side(v["k"]) == S["obs_side"])
        if engine == "relax" else ""))
    if S["partial_seeds"]:
        p("  warning: %d label(s) were read from fewer seeds than written (some seeds excluded): %s" % (
            len(S["partial_seeds"]), ", ".join(sorted(S["partial_seeds"]))))
    if S["below_null"]:
        p("  warning: %d label(s) whose best alternative fit is below the best null fit (not converged; read as LRT "
          "0, add seeds): %s" % (len(S["below_null"]), ", ".join("%s LRT %.3f" % s for s in S["below_null"])))
    if S["spread"]:
        p("  warning: %d alternative(s) with logL spread above %.2f across seeds: %s" % (
            len(S["spread"]), logl_tolerance, ", ".join("%s %.3f" % s for s in S["spread"])))
    read = set(k for k, _ in alts)
    aw = [(r["run"], w) for r in rows if r["label"] in read and not r["res"]["problem"] for w in r["res"]["warnings"]]
    for run, w in aw:
        p("  warning (%s): %s" % (run, w))
    p("")
    p("## per alternative (label, length, length-matched, LRT, P, %s, branches), by LRT" % (
        "K" if engine == "relax" else
        "foreground omega" if meta.get("codeml_model") == "two-ratio" else "foreground omega2"))
    for k, v in sorted(alts, key=lambda kv: -kv[1]["lrt"]):
        p("  %s  %.5f  %s  %.3f  %.3g  %s  %s" % (k, v["length"], "y" if v["length_matched"] else "n",
                                                  v["lrt"], pval(v["lrt"]),
                                                  "%.3f" % v["k"] if v["k"] is not None else "-", v["branches"]))
    return "\n".join(L) + "\n"


def table(out):
    """Per-run TSV of what was read, for independent recounting."""
    meta, rows, labels, excluded = build(out)
    head = "run\tlabel\tseed\tproblem\tLRT\tp\tK\tbest_logL\tnull_logL\treported_LRT\treported_p\treported_K\tn_fits\twarnings"
    lines = [head]
    for r in rows:
        x = r["res"]
        f = lambda v: "" if v is None else repr(v)
        lines.append("\t".join([r["run"], r["label"], r["seed"], x["problem"] or "", f(x.get("LRT")), f(x.get("p")),
                                f(x.get("K")), f(x.get("best_logL")), f(x.get("null_logL")), f(x.get("reported_LRT")),
                                f(x.get("reported_p")), f(x.get("reported_K")), str(len(x.get("fits") or [])),
                                " | ".join(x.get("warnings") or [])]))
    return "\n".join(lines) + "\n"
