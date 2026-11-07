"""Synthetic RELAX outputs in the layout HyPhy 2.5.101 prints, for testing the reader and the report.

Numbers here are invented. They exercise the instrument (parsing, set checks, logL reading, rank arithmetic);
no test reads a real result that a verdict depends on.
"""
import json
import os


def out_text(test, ref, nuis, alt=(-1000.0, 0.6), other=None, null=-1003.0, p=None, evidence=None,
             deleted=(), reset=False, refit=None):
    L = ["### Branch sets for RELAX analysis",
         "* Selected %d branches as the _Reference_ set: `%s`" % (len(ref), ", ".join(sorted(ref))),
         "* Selected %d branches as the _Test_ set: `%s`" % (len(test), ", ".join(sorted(test))),
         ] + (["* %d branches are in the unclassified (nuisance) set: `%s`" % (len(nuis), ", ".join(sorted(nuis)))]
              if nuis else []) + [   # RELAX prints the nuisance line only when the set is not empty
         "",
         "### Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model"]
    if deleted:
        L.append("### Deleted %d zero-length internal branches: `%s`" % (len(deleted), ", ".join(deleted)))
    L += ["* Log(L) = -1200.00, AIC-c = 2500.00 (20 estimated parameters)",
          "### Fitting the general descriptive (separate k per branch) model",
          "### * Log(L) = -990.00, AIC-c = 2100.00 (40 estimated parameters)",
          "### Fitting the alternative model to test K != 1",
          "* Log(L) = %.2f, AIC-c = 2050.00 (25 estimated parameters)" % alt[0],
          "* Relaxation/intensification parameter (K) =     %.2f" % alt[1]]
    if other is not None:
        L += ["### * Potential convergence issues due to flat likelihood surfaces; checking to see whether K > 1 "
              "or K < 1 is robustly inferred",
              "### * Attempt to fit an alternative direction of K",
              "* Log(L) = %.2f, AIC-c = 2050.00 (25 estimated parameters)" % other[0],
              "* Relaxation/intensification parameter (K) =     %.2f" % other[1]]
        if reset:
            L += ["### Potential for highly unreliable K inference due to multiple local maxima in the likelihood "
                  "function, treat results with caution ",
                  "> Relaxation parameter reset to opposite mode of evolution from that obtained in the initial "
                  "optimization."]
    L += ["### Fitting the null (K := 1) model",
          "* Log(L) = %.2f, AIC-c = 2051.00 (24 estimated parameters)" % null]
    if refit:
        # later passes after a negative LRT: same messages, no section headers (HyPhy prints a header once)
        L += ["* Detected convergence issues (negative LRT). Refitting the alterative/null model pair from a new "
              "starting point",
              "* Log(L) = %.2f, AIC-c = 2050.00 (25 estimated parameters)" % refit[0],
              "* Relaxation/intensification parameter (K) =     %.2f" % refit[1],
              "* Log(L) = %.2f, AIC-c = 2051.00 (24 estimated parameters)" % refit[2]]
    L += [
          "----",
          "## Test for relaxation (or intensification) of selection [RELAX]"]
    if p is not None:
        L.append("Likelihood ratio test **p =   %.4f**." % p)
        if evidence:
            L.append(">Evidence for *%s of selection* among **test** branches _relative_ to the **reference** "
                     "branches at P<=0.05" % evidence)
        else:
            L.append(">No significant evidence for relaxation (or intensification) of selection among **test** "
                     "branches _relative_ to the **reference** branches at P<=0.05")
    L += ["----", "### Fitting the partitioned descriptive model (completely separate rate distributions for branch "
          "sets)", "* Log(L) = -999.00, AIC-c = 2060.00 (29 estimated parameters)"]
    return "\n".join(L) + "\n"


def json_text(alt_logL, null_logL, lrt, p, k, other=None, raw_inf=True, unstable=None, negative_lrt=None):
    d = {"analysis": {"version": "4.7", "settings": {"convergence-flat-surface": -1.0}},
         "fits": {"RELAX alternative": {"Log Likelihood": alt_logL, "estimated parameters": 25},
                  "RELAX null": {"Log Likelihood": null_logL, "estimated parameters": 24}},
         "test results": {"LRT": lrt, "p-value": p, "relaxation or intensification parameter": k},
         "branch attributes": {"0": {"x": {"k (general descriptive)": "__INF__"}}}}
    if other is not None:
        d["analysis"]["settings"]["convergence-unstable-alernative"] = [[other[1], other[0]]]
    if unstable is not None:
        d["analysis"]["settings"]["convergence-unstable"] = [[k_, l_] for l_, k_ in unstable]
    if negative_lrt is not None:
        d["analysis"]["settings"]["convergence-negative-lrt"] = negative_lrt
    t = json.dumps(d, indent=1)
    return t.replace('"__INF__"', "inf" if raw_inf else "1e308")


def write_run(d, exp, **kw):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "expected_sets.json"), "w") as fh:
        json.dump(exp, fh)
    got_test = kw.pop("got_test", exp["Test"])
    got_ref = kw.pop("got_ref", exp["Reference"])
    got_nui = kw.pop("got_nuis", exp["nuisance"])
    js = kw.pop("json", None)
    with open(os.path.join(d, "relax.out"), "w") as fh:
        fh.write(out_text(got_test, got_ref, got_nui, **kw))
    if js is not None:
        with open(os.path.join(d, "relax.json"), "w") as fh:
            fh.write(js)
    with open(os.path.join(d, "DONE"), "w") as fh:
        fh.write("ok\t0\t1s\tnow\n")
