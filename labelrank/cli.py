"""labelrank command line: check, setup, run, report, inspect, batch-setup, batch-report, slurm."""
import argparse
import os
import sys

from . import __version__


def _write_atomic(path, text):
    """Write to a temporary file beside path, then rename: a failed run never leaves a half or wrong file."""
    d = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(d):
        os.makedirs(d)
    tmp = "%s.tmp%d" % (path, os.getpid())
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def _names(s):
    names = [x.strip() for x in s.split(",")]
    if not all(names):
        raise argparse.ArgumentTypeError("empty name in %r" % s)
    return names


def _ints(s):
    return [int(x) for x in s.split(",") if x.strip()]


def main(argv=None):
    ap = argparse.ArgumentParser(prog="labelrank", description=(
        "Rank a trait-labelled branch-selection test (HyPhy RELAX, or codeml branch-site model A or the two-ratio "
        "branch model) among "
        "every same-shaped label on the same tree."))
    ap.add_argument("--version", action="version", version="labelrank " + __version__)
    sub = ap.add_subparsers(dest="cmd")

    def inputs(p, files=True):
        if files:
            p.add_argument("--alignment", required=True, help="codon alignment, FASTA")
            p.add_argument("--tree", required=True, help="Newick tree with the test branches labelled")
        p.add_argument("--test", required=files, help="test label: HyPhy label (e.g. 'Foreground') or '#1' for codeml")
        p.add_argument("--reference", help="RELAX reference label; default 'Unlabeled branches' when such branches exist")
        p.add_argument("--engine", choices=("relax", "codeml"), default="relax")
        p.add_argument("--code", type=int, default=1,
                       help="genetic code as an NCBI transl_table number (default 1, universal; 2 = vertebrate "
                            "mitochondrial); sets the stop-codon check, HyPhy's --code and codeml's icode")
        p.add_argument("--exclude", action="append", type=_names, default=None, metavar="NAMES",
                       help="keep a lineage out of the alternatives, e.g. one you expect to share the trait: a tip "
                            "or node name, or tips joined by '+' for their clade; comma-separated or repeated")

    c = sub.add_parser("check", help="check the inputs and list the alternatives; writes nothing")
    inputs(c)
    c.add_argument("--max-alternatives", type=int, default=None)
    c.add_argument("--sample-seed", type=int, default=1)
    c.add_argument("--length-ratio", type=float, default=2.0)

    def setup_options(s):
        s.add_argument("--out", required=True, help="new, empty output folder")
        s.add_argument("--seeds", type=_ints, default=[1, 2, 3],
                       help="seeds (optimizer restarts) for every label (default 1,2,3)")
        s.add_argument("--alt-seeds", type=_ints, default=None,
                       help="different seeds for the alternatives (default: the same as --seeds; fewer favour the "
                            "observed label, and the report says so)")
        s.add_argument("--max-alternatives", type=int, default=200,
                       help="uniform sample when more alternatives exist (default 200; 0 = all)")
        s.add_argument("--sample-seed", type=int, default=1)
        s.add_argument("--length-ratio", type=float, default=2.0,
                       help="length-matched = test-set length within this factor of the observed one (default 2)")
        s.add_argument("--length-matched-only", action="store_true", help="write only length-matched alternatives")
        s.add_argument("--sequential", type=int, metavar="H",
                       help="stop running alternatives once H of them reach the observed LRT (Besag-Clifford); "
                            "needs 'labelrank run'")
        s.add_argument("--alpha", type=float, default=0.05, help="gate level, fixed here (default 0.05)")
        s.add_argument("--direction", choices=("relax", "intens"),
                       help="RELAX: the gate also requires K on this side of 1 (e.g. the published direction)")
        s.add_argument("--hyphy", default="hyphy", help="hyphy executable")
        s.add_argument("--relax-bf", help="RELAX.bf to seed (default: the one shipped with --hyphy)")
        s.add_argument("--codeml", default="codeml", help="codeml executable")
        s.add_argument("--model", choices=("branch-site", "two-ratio"), default=None,
                       help="codeml model (default branch-site): branch-site = model A against omega2 = 1, P from "
                            "the 50:50 mixture; two-ratio = model 2 against one ratio (model 0), NSsites = 0, P from "
                            "chi-square(1)")

    s = sub.add_parser("setup", help="write one run folder per label and seed")
    inputs(s)
    setup_options(s)

    bs = sub.add_parser("batch-setup", help="setup for many genes from a manifest (gene, alignment, tree[, test, reference, exclude])")
    bs.add_argument("--manifest", required=True, help="TSV with a header line; paths relative to the manifest's folder")
    inputs(bs, files=False)
    setup_options(bs)

    br = sub.add_parser("batch-report", help="one row per gene, and the pooled reading")
    br.add_argument("out")
    br.add_argument("--alpha", type=float)
    br.add_argument("--direction", choices=("relax", "intens", "none"))
    br.add_argument("--logl-tolerance", type=float, default=0.5)
    br.add_argument("--table", help="write the per-gene TSV here (default <out>/genes.tsv)")

    ad = sub.add_parser("add-seeds", help="add run folders for more seeds to an existing setup (e.g. to make "
                                          "alternatives match the observed label's seeds)")
    ad.add_argument("out")
    ad.add_argument("--seeds", type=_ints, required=True, help="seeds to add")
    ad.add_argument("--which", choices=("alternatives", "observed", "all"), default="alternatives")

    sl = sub.add_parser("slurm", help="write a SLURM array script and a command list for unfinished runs")
    sl.add_argument("out")
    sl.add_argument("--partition")
    sl.add_argument("--account")
    sl.add_argument("--time", default="48:00:00")
    sl.add_argument("--mem", default="2G")
    sl.add_argument("--max-parallel", type=int, default=100)

    r = sub.add_parser("run", help="run every unfinished run folder")
    r.add_argument("out")
    r.add_argument("--jobs", type=int, default=1, help="runs at a time, one CPU each (default 1)")
    r.add_argument("--timeout-hours", type=float, default=48)

    p = sub.add_parser("report", help="instrument check, gate and rank")
    p.add_argument("out")
    p.add_argument("--alpha", type=float, help="default: the value fixed at setup")
    p.add_argument("--direction", choices=("relax", "intens", "none"), help="default: the value fixed at setup")
    p.add_argument("--logl-tolerance", type=float, default=0.5)
    p.add_argument("--save", help="also write the report text to this file")
    p.add_argument("--table", help="also write the per-run TSV here")
    p.add_argument("--methods", help="also write a methods paragraph, results sentence and figure caption (text)")
    p.add_argument("--figure", help="also write an SVG strip plot of the LRTs (single column, 85 mm)")

    i = sub.add_parser("inspect", help="read one RELAX output (relax.out and/or relax.json) by logL")
    i.add_argument("paths", nargs="+")

    a = ap.parse_args(argv)
    if getattr(a, "exclude", None):
        a.exclude = [x for group in a.exclude for x in group]
    if a.cmd is None:
        ap.print_help()
        return 1
    if a.cmd == "check":
        from .setup import prepare, check_options, SetupError as _SE
        try:
            check_options(a, setup=False)
        except _SE as e:
            print("ERROR: %s" % e)
            return 2
        from .alternatives import enumerate_alternatives, describe, total_length, label_counts
        from .alternatives import no_alternatives_message
        from . import newick
        try:
            root, records, ref, errors, warnings, info = prepare(a.alignment, a.tree, a.test, a.reference, a.engine,
                                                                 a.code)
        except (ValueError, OSError) as e:   # NewickError is a ValueError
            print("ERROR: cannot read the inputs: %s" % e)
            return 2
        for w in warnings:
            print("warning: " + w)
        for e in errors:
            print("ERROR: " + e)
        if errors:
            return 2
        from .alternatives import ShapeError
        from .alternatives import resolve_exclude
        try:
            excl, excl_desc = resolve_exclude(root, a.exclude, a.test, ref)
            newick.name_internal_nodes(root)
            alt = enumerate_alternatives(root, a.test, ref, a.max_alternatives or None, a.sample_seed, exclude=excl)
        except ShapeError as e:
            print("ERROR: " + str(e))
            return 2
        print("input ok: %d sequences, %d codons; %s" % (info["sequences"], info["codons"], "background = every unmarked branch"
              if a.engine == "codeml" else "reference %r" % ref))
        print("branch labels: %s" % label_counts(root))
        print("observed: %s (length %.5f)" % (describe(alt["test"]), total_length(alt["test"])))
        print("alternatives: %d%s" % (len(alt["alternatives"]), " (drawn at random from %d raw combinations)" % alt[
            "combinations"] if alt["drawn"] else " (sampled from %d)" % alt["n_found"] if alt["sampled"] else ""))
        if excl:
            print("excluded: %d reference branch(es) in %s" % (alt["excluded_reference"], "; ".join(excl_desc)))
        for k, m in enumerate(alt["alternatives"]):
            print("  alt%03d  %.5f  %s" % (k + 1, total_length(m), describe(m)))
        n = len(alt["alternatives"])
        if not n:
            print("ERROR: " + no_alternatives_message(root, len(alt["reference"]), ref, len(alt["test"])))
            return 2
        if a.max_alternatives is None and n > 200:
            print("note: setup writes a uniform sample of 200 of these by default (--max-alternatives); the lines "
                  "below are for all %d" % n)
        olen, r = total_length(alt["test"]), a.length_ratio
        nlm = sum(1 for m in alt["alternatives"]
                  if ((olen / r <= total_length(m) <= olen * r) if olen > 0 else total_length(m) <= 0))
        for what, k in (("all alternatives", n), ("length-matched (within x%g)" % r, nlm)):
            print("resolution, %s: %d, so the smallest rank-based P is 1/%d = %.3g%s" % (
                what, k, k + 1, 1.0 / (k + 1), "; N + 1 < 20: the report gives the rank only" if k + 1 < 20 else ""))
        return 0
    if a.cmd in ("setup", "batch-setup"):
        from .setup import setup, SetupError
        from .batch import batch_setup
        if a.max_alternatives == 0:
            a.max_alternatives = None
        if a.sequential is not None and a.sequential < 1:
            print("ERROR: --sequential needs H >= 1", file=sys.stderr)
            return 2
        try:
            (setup if a.cmd == "setup" else batch_setup)(a)
        except SetupError as e:
            print("ERROR: %s" % e, file=sys.stderr)
            return 2
        return 0
    if a.cmd in ("batch-report", "add-seeds", "slurm", "run", "report"):
        want = ("batch.json",) if a.cmd == "batch-report" else ("setup.json",) if a.cmd in ("add-seeds", "report") \
            else ("setup.json", "batch.json")
        if not any(os.path.isfile(os.path.join(a.out, w)) for w in want):
            print("ERROR: %s is not a labelrank %s folder (no %s); make one with 'labelrank %s'" % (
                a.out, "batch output" if want == ("batch.json",) else "output", " or ".join(want),
                "batch-setup" if want == ("batch.json",) else "setup"), file=sys.stderr)
            return 2
    if a.cmd in ("report", "batch-report"):
        if a.alpha is not None and not 0 < a.alpha < 1:   # also refuses nan
            print("ERROR: --alpha must be between 0 and 1, got %g" % a.alpha, file=sys.stderr)
            return 2
        if not a.logl_tolerance >= 0:
            print("ERROR: --logl-tolerance must be 0 or more, got %g" % a.logl_tolerance, file=sys.stderr)
            return 2
    if a.cmd == "batch-report":
        from .batch import batch_report
        txt, tsv = batch_report(a.out, a.alpha, a.direction, a.logl_tolerance)
        _write_atomic(a.table or os.path.join(a.out, "genes.tsv"), tsv)
        sys.stdout.write(txt)
        return 0
    if a.cmd == "add-seeds":
        from .setup import add_seeds, SetupError
        try:
            add_seeds(a.out, a.seeds, a.which)
        except SetupError as e:
            print("ERROR: %s" % e, file=sys.stderr)
            return 2
        return 0
    if a.cmd == "slurm":
        import re as _re
        if a.max_parallel < 1:
            print("ERROR: --max-parallel must be at least 1", file=sys.stderr)
            return 2
        for nm in ("partition", "account", "time", "mem"):
            v = getattr(a, nm)
            if v is not None and not _re.match(r"^[A-Za-z0-9:_.,+-]+$", v):
                print("ERROR: --%s %r: only letters, digits and : _ . , + - are allowed (it goes into an #SBATCH "
                      "line)" % (nm, v), file=sys.stderr)
                return 2
        from .slurm import write
        try:
            write(a.out, a.partition, a.account, a.time, a.mem, a.max_parallel)
        except ValueError as e:
            print("ERROR: %s" % e, file=sys.stderr)
            return 2
        return 0
    if a.cmd == "run":
        if a.jobs < 1:
            print("ERROR: --jobs must be at least 1", file=sys.stderr)
            return 2
        if not a.timeout_hours >= 0:   # also refuses nan
            print("ERROR: --timeout-hours must be positive (0 = no limit), got %g" % a.timeout_hours, file=sys.stderr)
            return 2
        from .runner import run_all
        try:
            nfail = run_all(os.path.abspath(a.out), a.jobs, a.timeout_hours * 3600 if a.timeout_hours else None)
        except RuntimeError as e:
            print("ERROR: %s" % e, file=sys.stderr)
            return 2
        return 130 if nfail < 0 else 1 if nfail else 0
    if a.cmd == "report":
        if os.path.exists(os.path.join(a.out, "batch.json")):
            print("ERROR: %s is a batch folder; use 'labelrank batch-report'" % a.out, file=sys.stderr)
            return 2
        from .report import report, table, summarize
        S = summarize(a.out, a.alpha, a.direction, a.logl_tolerance)
        if all(r["res"]["problem"] == "not run" for r in S["rows"]):
            print("ERROR: nothing has run yet in %s; run 'labelrank run %s' first" % (a.out, a.out), file=sys.stderr)
            return 2
        txt = report(a.out, S=S)
        # every output is computed in full before any file is written, then each is written atomically
        outs = []
        if a.save:
            outs.append((a.save, txt))
        if a.table:
            outs.append((a.table, table(a.out)))
        if a.methods or a.figure:
            from . import paper
            if a.methods:
                outs.append((a.methods, "Methods\n\n%s\n\nResults\n\n%s\n\nFigure caption\n\n%s\n" % (
                    paper.methods(S), paper.results(S), paper.caption(S))))
            if a.figure:
                outs.append((a.figure, paper.figure_svg(S)))
        for path, text in outs:
            _write_atomic(path, text)
        sys.stdout.write(txt)
        return 0
    if a.cmd == "inspect":
        from .hyphy_relax import read_result
        out = js = None
        for pth in a.paths:
            with open(pth, "rb") as fh:
                t = fh.read().decode("utf-8", "replace")
            if "CODONML" in t[:2000] or "lnL(ntime:" in t:
                print("ERROR: %s is a codeml output; inspect reads HyPhy RELAX outputs only" % pth, file=sys.stderr)
                return 2
            if t.lstrip().startswith("{"):
                js = t
            else:
                out = t
        if out and "RELAX" not in out and js is None:
            print("ERROR: %s does not look like a RELAX output" % a.paths[0], file=sys.stderr)
            return 2
        res = read_result(out, js)
        g = lambda x: "-" if x is None else "%.6g" % x
        for f in res["fits"]:
            print("fit %-28s logL %s  K %s  (from %s)" % (f["kind"], g(f["logL"]), g(f["K"]), f.get("src")))
        print("null logL %s" % g(res["null_logL"]))
        if res.get("LRT") is not None:
            print("by logL (higher alternative fit): LRT %.4f  P %.4g  K %s (%s)" % (
                res["LRT"], res["p"], g(res["K"]), res.get("best_kind")))
        print("as reported by RELAX: LRT %s  P %s  K %s  verdict %s" % (
            g(res["reported_LRT"]), g(res["reported_p"]), g(res["reported_K"]),
            {"none": "no significant evidence", None: "-"}.get(res["reported_evidence"], res["reported_evidence"])))
        for w in res["warnings"]:
            print("warning: " + w)
        return 0


if __name__ == "__main__":
    sys.exit(main())
