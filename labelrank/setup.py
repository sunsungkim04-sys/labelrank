"""Write one run folder per (label, seed): the observed label under several seeds, each alternative under one."""
import hashlib
import re
import json
import os

import random

from . import __version__, alignment, alternatives, gencode, newick

ENGINES = ("relax", "codeml")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class SetupError(Exception):
    pass


def prepare(aln_path, tree_path, test, reference=None, engine="relax", code=1):
    """Parse and check inputs. Returns (root, records, refarg, errors, warnings, info)."""
    with open(tree_path) as fh:
        text = fh.read()
    root = newick.parse(text)
    records = alignment.read_fasta(aln_path)
    _, _, stops = gencode.get(code, engine)
    errors, warnings, info = alignment.check(records, [t.name for t in newick.tips(root)], stops)
    br = newick.branches(root)
    if engine == "relax":
        if reference is None:
            if any(not n.label for n in br):
                reference = alternatives.UNLABELED
            else:
                errors.append("every branch carries a label: pass --reference with the reference label "
                              "(HyPhy 2.5.101 refuses 'Unlabeled branches' on such a tree)")
        if any(n.label.startswith("#") for n in br):
            errors.append("PAML marks (#n) in the tree: for codeml add --engine codeml; for RELAX label branches "
                          "as name{label}")
    else:
        if reference not in (None, alternatives.UNLABELED):
            errors.append("codeml has no reference label: the background is every unmarked branch; drop --reference")
        reference = alternatives.UNLABELED
        bad = sorted(set(n.label for n in br if n.label and n.label != test))
        if bad:
            errors.append("codeml takes one foreground mark (#1); other marks found: %s" % ", ".join(bad))
    if len(root.kids) == 2 and root.kids[0].label != root.kids[1].label:
        warnings.append("the two branches below a bifurcating root carry different labels (%r, %r); HyPhy merges "
                        "them into one branch" % (root.kids[0].label, root.kids[1].label))
    badint = [n.name for n in newick.nodes(root) if n.kids and n.name and not re.match(r"^[A-Za-z0-9_]+$", n.name)]
    if badint:
        errors.append("%d internal node name(s) contain characters HyPhy rewrites, e.g. %r" % (len(badint), badint[0]))
    if reference == test:
        errors.append("--reference is the same label as --test")
    if len(root.kids) == 2 and any(k.label == test for k in root.kids):
        errors.append("a test branch sits directly below a bifurcating root, which HyPhy merges with its sister "
                      "branch; reroot the tree (e.g. on an outgroup) so the test lineage is not root-adjacent")
    lens = [n.length for n in br]
    if all(not x for x in lens):
        warnings.append("the tree has no branch lengths (or all are 0): every label then has length 0, so the "
                        "length-matched reading and the length checks say nothing; give a tree with lengths")
    elif any(x is None for x in lens) and any(x is not None for x in lens):
        warnings.append("%d of %d branches have no length; they count as 0 in the length-matched reading" % (
            sum(1 for x in lens if x is None), len(lens)))
    badlen = [n for n in br if n.length is not None and (n.length < 0 or n.length != n.length
                                                         or n.length == float("inf"))]
    if badlen:
        errors.append("%d branch length(s) are negative or not finite, e.g. %r" % (len(badlen), badlen[0].name or "?"))
    zero = [n for n in br if n.length is not None and n.length == 0.0]
    if zero:
        warnings.append("%d branch(es) have length 0 in the input tree; HyPhy may delete zero-length internal "
                        "branches, from any set (reference included); deletions are reported per run" % len(zero))
    return root, records, reference, errors, warnings, info


def check_options(args, setup=True):
    """Refuse option values that would crash later or quietly make a reading meaningless. Shared with 'check'."""
    import math
    r = args.length_ratio
    if not (isinstance(r, float) and math.isfinite(r) and r >= 1):
        raise SetupError("--length-ratio is a factor of at least 1 (2 = within half to double), got %r" % r)
    if args.max_alternatives is not None and args.max_alternatives < 0:
        raise SetupError("--max-alternatives must be 0 (all) or more, got %d" % args.max_alternatives)
    if getattr(args, "exclude", None) is not None and any(not x.strip() for x in args.exclude):
        raise SetupError("--exclude got an empty name")
    if args.engine == "codeml" and args.test is not None and args.test != "#1":
        raise SetupError("for codeml the test label is the foreground mark '#1' (got %r)" % args.test)
    model = getattr(args, "model", None)
    if model is not None and args.engine != "codeml":
        raise SetupError("--model chooses the codeml model (branch-site or two-ratio); RELAX has none")
    if not setup:
        return
    if not 0 < args.alpha < 1:
        raise SetupError("--alpha must be between 0 and 1, got %g" % args.alpha)
    if getattr(args, "sequential", None) is not None and args.sequential < 1:
        raise SetupError("--sequential H must be at least 1, got %d" % args.sequential)
    for nm, ss in (("--seeds", args.seeds), ("--alt-seeds", args.alt_seeds)):
        if ss is None:
            continue
        if not ss or len(set(ss)) != len(ss):
            raise SetupError("%s must list distinct seeds, got %s" % (nm, ss))
        if args.engine == "codeml" and not all(1 <= x <= 6 for x in ss):
            raise SetupError("codeml seeds pick one of 6 starting points: use 1-6, got %s" % ss)
        if min(ss) < 0:
            raise SetupError("%s must be 0 or more, got %s" % (nm, ss))
    # cmd.sh runs inside each run folder, so a relative program or file path would point nowhere there
    for opt in ("hyphy", "codeml", "relax_bf"):
        v = getattr(args, opt, None)
        if v and (os.sep in v or os.path.isfile(v)):   # a bare name is looked up on PATH, even beside a folder of that name
            setattr(args, opt, os.path.abspath(v))


def setup(args, log=print):
    out = os.path.abspath(args.out)
    if os.path.exists(out) and not os.path.isdir(out):
        raise SetupError("%s exists and is a file; --out names a new folder" % out)
    if os.path.exists(out) and os.listdir(out):
        raise SetupError("%s exists and is not empty; refusing to write over earlier runs" % out)
    check_options(args)
    if args.engine == "codeml" and getattr(args, "direction", None):
        raise SetupError("--direction is a side of RELAX's K; it is not implemented for codeml (branch-site is "
                         "one-sided already; two-ratio ranks departures of the foreground omega either way)")
    model = (getattr(args, "model", None) or "branch-site") if args.engine == "codeml" else None
    if args.alt_seeds is None:
        args.alt_seeds = list(args.seeds)
    elif sorted(args.alt_seeds) != sorted(args.seeds):
        log("warning: alternatives get seeds %s, the observed label %s; unequal restarts favour the label with more"
            % (args.alt_seeds, args.seeds))
    code = getattr(args, "code", 1)
    try:
        code_name, icode, _ = gencode.get(code, args.engine)
    except ValueError as e:
        raise SetupError(str(e))
    try:
        root, records, refarg, errors, warnings, info = prepare(args.alignment, args.tree, args.test, args.reference,
                                                                args.engine, code)
    except (newick.NewickError, ValueError, OSError) as e:
        raise SetupError("cannot read the inputs: %s" % e)
    for w in warnings:
        log("warning: " + w)
    if errors:
        raise SetupError("input check failed:\n  - " + "\n  - ".join(errors))
    log("input: %d sequences, %d codons" % (info["sequences"], info["codons"]))
    try:
        excl, excl_desc = alternatives.resolve_exclude(root, getattr(args, "exclude", None), args.test, refarg)
    except alternatives.ShapeError as e:
        raise SetupError(str(e))
    named = newick.name_internal_nodes(root)
    try:
        alt = alternatives.enumerate_alternatives(root, args.test, refarg, args.max_alternatives, args.sample_seed,
                                                  exclude=excl)
    except alternatives.ShapeError as e:
        raise SetupError(str(e))
    if excl:
        log("excluded from the alternatives: %d reference branch(es) in %s" % (
            alt["excluded_reference"], "; ".join(excl_desc)))
    obs_len = alternatives.total_length(alt["test"])
    ratio = args.length_ratio

    def matched(m):
        L = alternatives.total_length(m)
        if obs_len <= 0:
            return L <= 0
        return obs_len / ratio <= L <= obs_len * ratio

    alts = alt["alternatives"]
    if args.length_matched_only:
        alts = [a for a in alts if matched(a)]
    taken = set(n.label for n in newick.branches(root))
    orig = "orig"
    while orig in taken:
        orig += "_"
    if args.engine == "codeml":
        orig = None   # codeml has two classes only: the published foreground returns to background
    shapes = ", ".join(alternatives.shape_words(k, t) for k, t in alt["shapes"])
    log("test set: %d branch(es) in %d group(s) [%s]; %s: %d branches; other: %d" % (
        len(alt["test"]), len(alt["components"]), shapes,
        "reference %s" % refarg if args.engine == "relax" else "background", len(alt["reference"]),
        len(alt["other"])))
    if alt["drawn"]:
        log("alternatives: %d drawn at random from %d raw combinations (too many to list; the number of valid "
            "ones is not known); %d written" % (len(alt["alternatives"]), alt["combinations"], len(alts)))
    else:
        log("alternatives: %d found; %d written%s%s" % (
            alt["n_found"], len(alts), " (uniform sample, seed %d)" % args.sample_seed if alt["sampled"] else "",
            " (length-matched only)" if args.length_matched_only else ""))
    h = getattr(args, "sequential", None)
    if h and h > len(alts):
        log("warning: --sequential %d exceeds the %d alternatives; every alternative will run" % (h, len(alts)))
    if len(alts) + 1 < 20:
        log("note: N + 1 = %d < 20; the report will give the rank only, no P value" % (len(alts) + 1))
    if not alts:   # a tree labelled all over except one branch gives this with the default reference
        raise SetupError(alternatives.no_alternatives_message(root, len(alt["reference"]), refarg, len(alt["test"]),
                                                              bool(args.length_matched_only and alt["alternatives"])))

    if args.engine == "relax":
        from . import hyphy_relax as eng
        import shutil as _sh
        if not (_sh.which(args.hyphy) or os.path.isfile(args.hyphy)):
            raise SetupError("hyphy executable %r not found on PATH; pass --hyphy /path/to/hyphy" % args.hyphy)
        try:
            relax_bf = args.relax_bf or eng.find_relax_bf(args.hyphy)
        except FileNotFoundError as e:
            raise SetupError(str(e))
        version = eng.hyphy_version(args.hyphy)
        if not version.upper().startswith("HYPHY"):
            log("warning: '%s --version' did not report a HyPhy version: %s" % (args.hyphy, version[:120]))
    else:
        from . import codeml as eng
        relax_bf, version = None, eng.codeml_version(args.codeml)
    if gencode.release(version) not in eng.TESTED:
        log("warning: %s is not a release this tool was checked against (%s); output formats change between "
            "releases, so check a few runs with 'labelrank inspect' or by hand before trusting the report" % (
                version, ", ".join(eng.TESTED)))

    created = not os.path.exists(out)
    os.makedirs(os.path.join(out, "runs"))
    try:
        return _write(args, out, root, records, refarg, alt, alts, orig, matched, ratio, obs_len, named, version,
                      relax_bf, code, code_name, icode, excl_desc, warnings, eng, log, model)
    except BaseException as e:
        # a half-written folder would block the retry ("exists and is not empty")
        import shutil
        shutil.rmtree(out if created else os.path.join(out, "runs"), ignore_errors=True)
        if isinstance(e, (ValueError, OSError)):
            raise SetupError("could not write the run folders (nothing was kept): %s" % e)
        raise


def _write(args, out, root, records, refarg, alt, alts, orig, matched, ratio, obs_len, named, version, relax_bf,
           code, code_name, icode, excl_desc, warnings, eng, log, model=None):
    test_set = set(alt["test"])
    rows = ["run\trole\tlabel\tseed\tn_branches\tlength\tlength_matched\tbranches\torder"]
    labels = [("observed", None)] + [("alt%03d" % (i + 1), a) for i, a in enumerate(alts)]
    # a random running order, drawn now: the sequential (early-stopping) reading walks alternatives in this order
    perm = list(range(1, len(alts) + 1))
    random.Random("order-%d" % args.sample_seed).shuffle(perm)
    order_of = {"alt%03d" % (i + 1): perm[i] for i in range(len(alts))}
    for lab_id, members in labels:
        relabel = {}
        if members is not None:
            for n in test_set:
                relabel[n] = orig if orig else ""
            for n in members:
                relabel[n] = args.test
            seeds = args.alt_seeds
        else:
            members = alt["test"]
            seeds = args.seeds
        lab_of = lambda n: relabel.get(n, n.label)
        brs = newick.branches(root)
        exp = {"Test": sorted(n.name for n in brs if lab_of(n) == args.test)}
        if refarg == alternatives.UNLABELED:
            exp["Reference"] = sorted(n.name for n in brs if not lab_of(n))
        else:
            exp["Reference"] = sorted(n.name for n in brs if lab_of(n) == refarg)
        inset = set(exp["Test"]) | set(exp["Reference"])
        exp["nuisance"] = sorted(n.name for n in brs if n.name not in inset)
        # HyPhy unroots a bifurcating root by merging its two child branches, so one of them is never printed
        exp["root_children"] = [k.name for k in root.kids] if len(root.kids) == 2 else []
        for seed in seeds:
            run = "%s_s%d" % (lab_id, seed)
            d = os.path.join(out, "runs", run)
            if args.engine == "relax":
                eng.write_run(d, args.alignment, newick.write_hyphy(root, relabel), seed, args.test, refarg,
                              relax_bf, args.hyphy, code_name)
            else:
                eng.write_run(d, records, newick.write_paml(root, relabel), seed, args.codeml, icode, model)
            with open(os.path.join(d, "expected_sets.json"), "w") as fh:
                json.dump(exp, fh, indent=1)
            rows.append("%s\t%s\t%s\t%d\t%d\t%.9g\t%s\t%s\t%s" % (
                run, "observed" if lab_id == "observed" else "alternative", lab_id, seed, len(members),
                alternatives.total_length(members), "yes" if matched(members) else "no",
                alternatives.describe(members), order_of.get(lab_id, "")))
    with open(os.path.join(out, "runs.tsv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    meta = {"tool": "labelrank", "engine": args.engine, "codeml_model": model,
            "alignment": os.path.abspath(args.alignment),
            "alignment_sha256": sha256(args.alignment), "tree": os.path.abspath(args.tree),
            "tree_sha256": sha256(args.tree), "test": args.test, "reference": refarg, "orig_label": orig,
            "seeds": args.seeds, "alt_seeds": args.alt_seeds, "length_ratio": ratio,
            "length_matched_only": args.length_matched_only, "max_alternatives": args.max_alternatives,
            "sample_seed": args.sample_seed, "alternatives_found": alt["n_found"], "alternatives_written": len(alts),
            "alternatives_before_length_filter": len(alt["alternatives"]), "drawn_at_random": alt["drawn"],
            "raw_combinations": alt["combinations"],
            "sampled": alt["sampled"], "shapes": alt["shapes"], "internal_nodes_named": named,
            "observed_length": obs_len, "engine_version": version,
            "relax_bf": relax_bf, "relax_bf_sha256": sha256(relax_bf) if relax_bf else None,
            "genetic_code": code, "genetic_code_name": code_name,
            "exclude": excl_desc, "excluded_reference_branches": alt["excluded_reference"],
            "input_warnings": warnings, "sequential_h": getattr(args, "sequential", None) or None,
            "alpha": getattr(args, "alpha", 0.05), "direction": getattr(args, "direction", None),
            "bifurcating_root": len(root.kids) == 2,
            "tool_version": __version__}
    with open(os.path.join(out, "setup.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    log("wrote %d run folders under %s/runs (engine %s%s: %s)" % (
        len(rows) - 1, out, args.engine, " " + model if model else "", version))
    return out


def add_seeds(out, seeds, which="alternatives", log=print):
    """Copy each existing label's first run folder to new seeds: same tree, alignment and command; only the seed
    line (RELAX) or the starting values (codeml) change. runs.tsv and setup.json are updated."""
    import shutil
    out = os.path.abspath(out)
    with open(os.path.join(out, "setup.json")) as fh:
        meta = json.load(fh)
    with open(os.path.join(out, "runs.tsv")) as fh:
        lines = fh.read().splitlines()
    head, rows = lines[0].split("\t"), [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:]]
    first = {}
    for r in rows:
        first.setdefault(r["label"], r)
    have = set((r["label"], int(r["seed"])) for r in rows)
    # every check before anything is copied, so a refusal leaves the folder as it was
    if len(set(seeds)) != len(seeds) or not seeds:
        raise SetupError("--seeds must list distinct seeds, got %s" % seeds)
    if min(seeds) < 0:
        raise SetupError("--seeds must be 0 or more, got %s" % seeds)
    if meta["engine"] == "codeml":
        from .codeml import STARTS
        if any(not 1 <= sd <= len(STARTS) for sd in seeds):
            raise SetupError("codeml seeds pick one of %d starting points: use 1-%d" % (len(STARTS), len(STARTS)))
    for r in rows:
        d = os.path.join(out, "runs", r["run"])
        started = os.path.exists(os.path.join(d, "relax.out")) or os.path.exists(os.path.join(d, "alt", "mlc"))
        if started and not os.path.exists(os.path.join(d, "DONE")):
            raise SetupError("run %s has output but no DONE marker: a run may be in progress; add seeds only when "
                             "'labelrank run' is not running on this folder" % r["run"])
    todo = [(lab, r, sd) for lab, r in first.items() for sd in seeds
            if not (which == "alternatives" and lab == "observed" or which == "observed" and lab != "observed")
            and (lab, sd) not in have]
    for lab, r, sd in todo:
        if os.path.exists(os.path.join(out, "runs", "%s_s%d" % (lab, sd))):
            raise SetupError("%s_s%d exists but is not in runs.tsv; remove it first" % (lab, sd))
    new = []
    for lab, r, sd in todo:
        src = os.path.join(out, "runs", r["run"])
        run = "%s_s%d" % (lab, sd)
        dst = os.path.join(out, "runs", run)
        if os.path.exists(dst):
            raise SetupError("%s already exists" % dst)
        shutil.copytree(src, dst, symlinks=True, ignore=shutil.ignore_patterns(
            "DONE", "relax.out", "relax.json", "errors.log", "*.log", "mlc", "rst*", "rub", "lnf", "2NG.*",
            "*.cache", "*.json.tmp*"))
        if meta["engine"] == "relax":
            bf = os.path.join(dst, "RELAX_seeded.bf")
            with open(bf) as fh:
                body = fh.read().split("\n", 1)
            if not body[0].startswith("SetParameter (RANDOM_SEED,"):
                raise SetupError("%s does not start with the seed line" % bf)
            with open(bf, "w") as fh:
                fh.write("SetParameter (RANDOM_SEED, %d, 0);\n" % sd + body[1])
        else:
            from .codeml import STARTS
            omega, kappa = STARTS[sd - 1]
            ctl = os.path.join(dst, "alt", "codeml.ctl")
            with open(ctl) as fh:
                t = fh.read()
            t = re.sub(r"(?m)^(\s*omega = ).*$", r"\g<1>%s" % omega, t)
            t = re.sub(r"(?m)^(\s*kappa = ).*$", r"\g<1>%s" % kappa, t)
            with open(ctl, "w") as fh:
                fh.write(t)
            nctl = os.path.join(dst, "null", "codeml.ctl")
            with open(nctl) as fh:
                t = fh.read()
            t = re.sub(r"(?m)^(\s*kappa = ).*$", r"\g<1>%s" % kappa, t)
            if re.search(r"(?m)^\s*fix_omega = 0\s*$", t):   # two-ratio: the one-ratio null starts from omega too
                t = re.sub(r"(?m)^(\s*omega = ).*$", r"\g<1>%s" % omega, t)
            with open(nctl, "w") as fh:
                fh.write(t)
        row = dict(r, run=run, seed=str(sd))
        new.append("\t".join(row.get(h, "") for h in head))
    if not new:
        log("nothing to add")
        return
    key = {"alternatives": ["alt_seeds"], "observed": ["seeds"], "all": ["seeds", "alt_seeds"]}[which]
    for k in key:
        meta[k] = sorted(set(meta[k]) | set(seeds))
    meta.setdefault("seed_additions", []).append({"seeds": seeds, "which": which})
    # both files are written to temporaries first, then swapped in one after the other
    with open(os.path.join(out, "runs.tsv.tmp"), "w") as fh:
        fh.write("\n".join(lines + new) + "\n")
    with open(os.path.join(out, "setup.json.tmp"), "w") as fh:
        json.dump(meta, fh, indent=1)
    os.replace(os.path.join(out, "setup.json.tmp"), os.path.join(out, "setup.json"))
    os.replace(os.path.join(out, "runs.tsv.tmp"), os.path.join(out, "runs.tsv"))
    log("added %d run folder(s); run 'labelrank run %s' to run them" % (len(new), out))
