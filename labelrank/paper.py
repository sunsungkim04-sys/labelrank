"""Text and a figure for a manuscript: a methods paragraph, a results sentence, a caption, and an SVG strip plot.

The figure follows common journal practice: single-column width (85 mm), Arial 7.5 pt at print size, Okabe-Ito
colours, no grid, no explanatory text inside the plot (the caption carries it).
"""
from . import __version__, CITATION
from .report import LRT_FLOOR

GREY, VERMILLION, BLUE = "#7f7f7f", "#D55E00", "#0072B2"


def _shape_words(meta):
    out = []
    for kind, topo in meta["shapes"]:
        if kind == "tip":
            out.append("a single terminal branch")
        else:
            n = topo.count("t")
            out.append(("a %d-tip clade with all its branches" if kind == "clade" else
                        "the stem branch of a %d-tip clade") % n)
    return " plus ".join(out)


def _short_version(v):
    """'HYPHY 2.5.101(MP) for Linux on x86_64 ...' -> 'HyPhy 2.5.101(MP)'."""
    import re
    m = re.match(r"HYPHY\s+(\S+)", v or "", re.I)
    return "HyPhy " + m.group(1) if m else (v or "HyPhy")


def methods(S):
    m = S["meta"]
    relax = m["engine"] == "relax"
    two = m.get("codeml_model") == "two-ratio"
    eng = ("HyPhy RELAX" if relax else
           "codeml two-ratio branch model (model = 2 against the one-ratio model = 0, NSsites = 0; PAML)" if two else
           "codeml branch-site test (model A against its omega2 = 1 null; PAML)")
    if m.get("drawn_at_random"):
        which = "%d sets drawn at random (sampling seed %d) from the possible sets" % (
            m["alternatives_before_length_filter"], m["sample_seed"])
    elif m["sampled"]:
        which = "a uniform random sample (sampling seed %d) of %d of the %d possible sets" % (
            m["sample_seed"], m["alternatives_before_length_filter"], m["alternatives_found"])
    else:
        which = "each of the %d possible sets" % m["alternatives_found"]
    if m.get("length_matched_only"):
        which += (", keeping the %d whose total branch length was within a factor of %g of the test set's" % (
            m["alternatives_written"], m["length_ratio"]))
    ref = ("the branches labelled %s" % m["reference"] if m["reference"] != "Unlabeled branches"
           else "the unlabelled branches")
    moved = ("moved to a set that is neither test nor reference, so the command line was unchanged" if relax
             else "returned to the background")
    txt = [
        "To calibrate the trait-labelled result against comparable lineages of the same tree, we compared its %s "
        "result with every same-shaped label, using labelrank v%s (%s). The test set (%s) was moved in turn to %s of "
        "branches with the same shape lying entirely within the reference set (%s)%s; in each, the original test "
        "branches were %s." % (
            eng, __version__, CITATION, _shape_words(m), which, ref,
            ", excluding the two branches below the root" if m.get("bifurcating_root", True) else "", moved)]
    if m.get("exclude"):
        txt.append("Lineages named in advance were not used as alternatives (%s; each with every branch below it)."
                   % "; ".join(m["exclude"]))
    if m.get("genetic_code", 1) != 1:
        txt.append("Codons were read with NCBI genetic code table %d (%s)." % (
            m["genetic_code"], m.get("genetic_code_name")))
    seeds = ", ".join(map(str, m["seeds"]))
    same = sorted(m["seeds"]) == sorted(m["alt_seeds"])
    if relax:
        txt.append(
            "Every label was analysed with the same command (%s, RELAX with the random seed fixed in the batch file, "
            "CPU=1) under seeds %s%s. For each label, the alternative and the null model were each taken at their "
            "highest log-likelihood over the seeds (for the alternative, over both of RELAX's fits when it refitted K "
            "from the other side of 1), and LRT = 2(logL alternative - logL null), with P from chi-square(1)." % (
                _short_version(m["engine_version"]), seeds,
                "" if same else " for the observed label and %s for the alternatives" % ", ".join(
                    map(str, m["alt_seeds"]))))
    else:
        from .codeml import STARTS
        starts = "; ".join("(%g, %g)" % STARTS[s - 1] for s in sorted(set(m["seeds"]) | set(m["alt_seeds"])))
        txt.append(
            "Every label was analysed with the same codeml settings (%s; CodonFreq = 2, cleandata = 0, each fit in its "
            "own folder) from the starting (omega, kappa) values %s%s. For each label, each model was taken at its "
            "highest log-likelihood over the starts, LRT = 2(logL alternative - logL null), and P from %s." % (
                S.get("engine_version_read") or "PAML", starts,
                "" if same else " (starts %s for the observed label and %s for the alternatives)" % (
                    ", ".join(map(str, m["seeds"])), ", ".join(map(str, m["alt_seeds"]))),
                "chi-square(1) (the two models are nested at the interior point omega foreground = omega background)"
                if two else "the 50:50 mixture of 0 and chi-square(1) (the boundary null; the plain chi-square(1) is "
                "the conservative alternative)"))
    if relax:
        txt.append("Runs whose branch sets, as printed by HyPhy, differed by name from those written, and runs that "
                   "failed, were excluded.")
    else:
        txt.append("Runs that failed were excluded.")
    side = {"relax": " with K < 1", "intens": " with K > 1"}.get(S.get("direction") or "", "")
    txt.append(
        "The observed label's rank is 1 + the number of alternatives with LRT at or above it, out of N + 1, read over "
        "all alternatives and over the length-matched ones (total length within a factor of %g of the test set's), "
        "since short lineages carry few substitutions and give small LRTs regardless of selection; the share of "
        "alternatives shorter than the test set and Spearman's correlation between length and LRT among the "
        "alternatives (LRTs below %g counted as 0; given only when at least 3 alternatives reach it) are reported. "
        "The observed "
        "result counted as significant by the usual reading when its P <= %g%s." % (
            m["length_ratio"], LRT_FLOOR, S["alpha"], side))
    if m.get("sequential_h"):
        txt.append("Alternatives were run in a pre-drawn random order until %d of them reached the observed LRT, "
                   "giving the sequential P = h / L after L alternatives, or (g + 1) / (n + 1) if all n were run with "
                   "g at or above it (Besag and Clifford 1991); alternatives whose runs all failed were skipped."
                   % m["sequential_h"])
    else:
        txt.append("A rank-based P, rank / (N + 1), is reported only when N + 1 >= 20.")
    txt.append("Because the lineages of one tree differ in history and are not independent, the rank is a "
               "calibration against this tree's lineages, not an exact null distribution.")
    return " ".join(txt)


def sq_mode(S):
    return bool(S.get("sequential"))


def results(S):
    if not S.get("evaluable"):
        return "The observed label could not be evaluated: %d of %d of its runs were usable." % (
            len(S["obs_ok"]), S["n_seed"])
    N, lm = S["N"], S["length_matched"]
    s = "The observed label gave LRT %.2f (P %s by the usual reading). " % (
        S["lobs"], ("%.2g" % S["p_obs"]))
    if N == 0 and S.get("pending"):
        return s + "No alternative label had finished yet, so it was not ranked."
    if N == 0:
        return s + "No alternative label could be read (every alternative run failed or was excluded), so it was not ranked."
    if S.get("incomplete"):
        s += "Not every alternative had finished, so the ranks below are provisional. "
    s += "It ranked %d of %d among the same-shaped labels %s (%d alternative%s at or above it)" % (
        S["rank"], N + 1, "read before the sequential stop" if sq_mode(S) else "read",
        len(S["ge"]), "" if len(S["ge"]) == 1 else "s")
    sq = S.get("sequential")
    if sq:
        if not sq.get("complete"):
            s += "; the sequential reading is incomplete"
        elif sq.get("P") is not None:
            s += "; sequential P = %.3f after %d alternatives%s" % (
                sq["P"], sq["read"], "" if sq.get("stopped") else " (all run)")
    elif S["rank_p"] is not None:
        s += "; rank-based P = %.3f" % S["rank_p"]
    s += ", and %d of %d among the %d length-matched ones" % (S["lm_rank"], len(lm) + 1, len(lm))
    if S["lm_rank_p"] is not None:
        s += " (rank-based P = %.3f)" % S["lm_rank_p"]
    if (N + 1 < 20 or len(lm) + 1 < 20) and not sq:
        s += "; readings with N + 1 < 20 give a rank only"
    s += ". %.0f%% of the alternatives were shorter than the test set" % (100 * (S["shorter_share"] or 0))
    if S.get("length_rho") is not None:
        s += "; across all %d alternatives, length and LRT had Spearman's rho %.2f." % (N, S["length_rho"])
    elif N >= 3 and S.get("above_floor", 0) < 3:
        s += ("; the length-LRT correlation was not computed (%d of %d alternatives reached LRT %g)." % (
            S["above_floor"], N, LRT_FLOOR))
    else:
        s += "."
    s += (" %d of %d alternatives (%.0f%%) would be called significant by the chi-square reading at P <= %g." % (
        len(S["significant"]), N, 100.0 * len(S["significant"]) / max(1, N), S["alpha"]))
    if S.get("partial_seeds"):
        s += " %d label%s were read from fewer seeds than written." % (
            len(S["partial_seeds"]), "" if len(S["partial_seeds"]) == 1 else "s")
    if S["excluded"]:
        s += " %d run%s %s excluded (failed, or failed the instrument check)." % (
            len(S["excluded"]), "" if len(S["excluded"]) == 1 else "s", "was" if len(S["excluded"]) == 1 else "were")
    return s


def _plain_chi2(meta):
    return meta["engine"] == "relax" or meta.get("codeml_model") == "two-ratio"


def caption(S):
    crit = "3.84" if _plain_chi2(S["meta"]) else "2.71"
    return ("Likelihood-ratio statistic of the observed trait label (vermillion) and of every same-shaped "
            "alternative label read on the same tree (N = %d; filled: total branch length within a factor of %g of the "
            "observed test set's, open: outside it). Dashed line: LRT = %s, P = 0.05 by the %s reading." % (
                S.get("N", 0), S["meta"]["length_ratio"], crit,
                "chi-square(1)" if _plain_chi2(S["meta"]) else "50:50 mixture"))


def figure_svg(S):
    """Single-column strip plot of LRTs, observed marked. Returns SVG text."""
    W, H = 241.0, 92.0          # points; 85 mm wide
    left, right, base, top = 14.0, 8.0, 64.0, 12.0
    fs = 7.5
    altv = [(v["lrt"], v["length_matched"]) for _, v in S.get("alts", [])]
    alts = [x for x, _ in altv]
    lobs = S.get("lobs")
    crit = 3.841458820694124 if _plain_chi2(S["meta"]) else 2.705543454095404
    vals = alts + ([lobs] if lobs is not None else []) + [crit]
    xmax = max(vals) * 1.08 if vals else 1.0
    xmax = max(xmax, 1.0)
    X = lambda v: left + (W - left - right) * max(0.0, v) / xmax
    # stack dots in bins so equal values do not hide each other
    r = 2.2
    count = {}
    for v in alts:
        count[int(X(v) / (2 * r))] = count.get(int(X(v) / (2 * r)), 0) + 1
    # tall stacks (many near-zero alternatives) are compressed to fit the panel
    step = min(2 * r + 0.4, (base - top - fs - 8) / max(1, max(count.values()) if count else 1))
    bins = {}
    dots = []
    for v, matched in sorted(altv):
        b = int(X(v) / (2 * r))
        k = bins.get(b, 0)
        bins[b] = k + 1
        dots.append((X(v), base - 4 - k * step, matched))
    ticks = _ticks(xmax)
    o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%gpt" height="%gpt" viewBox="0 0 %g %g" '
         'font-family="Arial, Helvetica, sans-serif" font-size="%g">' % (W, H, W, H, fs),
         '<rect width="100%" height="100%" fill="white"/>',
         '<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="black" stroke-width="0.6"/>' % (left, base, W - right, base)]
    for t in ticks:
        o.append('<line x1="%.2f" y1="%g" x2="%.2f" y2="%g" stroke="black" stroke-width="0.6"/>' % (
            X(t), base, X(t), base + 2.5))
        o.append('<text x="%.2f" y="%g" text-anchor="middle">%s</text>' % (X(t), base + 2.5 + fs, _fmt(t)))
    o.append('<text x="%g" y="%g" text-anchor="middle">LRT</text>' % ((left + W - right) / 2, H - 3))
    o.append('<line x1="%.2f" y1="%g" x2="%.2f" y2="%g" stroke="%s" stroke-width="0.6" stroke-dasharray="2,1.5"/>' % (
        X(crit), top, X(crit), base, BLUE))
    o.append('<text x="%.2f" y="%g" text-anchor="middle" fill="%s" font-size="%g">P = 0.05</text>' % (
        X(crit), top - 1.5, BLUE, fs - 1))
    for x, y, matched in dots:
        if matched:
            o.append('<circle cx="%.2f" cy="%.2f" r="%g" fill="%s"/>' % (x, y, r, GREY))
        else:
            o.append('<circle cx="%.2f" cy="%.2f" r="%g" fill="white" stroke="%s" stroke-width="0.7"/>' % (
                x, y, r - 0.35, GREY))
    if lobs is not None:
        o.append('<line x1="%.2f" y1="%g" x2="%.2f" y2="%g" stroke="%s" stroke-width="1.2"/>' % (
            X(lobs), top + fs + 2, X(lobs), base, VERMILLION))
        anchor = "end" if X(lobs) > W * 0.7 else "start"
        o.append('<text x="%.2f" y="%g" text-anchor="%s" fill="%s">observed</text>' % (
            X(lobs) + (-1.5 if anchor == "end" else 1.5), top + fs, anchor, VERMILLION))
    o.append("</svg>")
    return "\n".join(o) + "\n"


def _ticks(xmax):
    import math
    mag = 10 ** math.floor(math.log10(max(xmax, 1e-9) / 6.0))
    step = next(m * mag for m in (1, 2, 5, 10) if xmax / (m * mag) <= 6)
    n = int(xmax / step) + 1
    return [i * step for i in range(n) if i * step <= xmax]


def _fmt(t):
    return "%g" % round(t, 10)
