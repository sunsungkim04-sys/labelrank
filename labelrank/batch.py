"""Many genes, one label: one output folder per gene, one table, and a pooled reading.

Manifest: a TSV with a header line and columns gene, alignment, tree, and optionally test, reference and exclude
(which override the command-line values for that row; exclude is comma-separated, as --exclude, since tip names
usually differ between genes). Paths are relative to the manifest's folder.

Pooled reading: if the trait label were exchangeable with its alternatives in every gene, its rank in a gene with
N alternatives would be uniform on 1..N+1. The table therefore gives, over every gene with a complete reading
(gate passed or not), the number ranked first against the number expected, sum of 1/(N+1). Genes on one tree share the
lineage being tested, so they are not independent replicates; this is a description, not a test.
"""
import copy
import json
import os

from .setup import SetupError, setup


RESERVED = {"batch.json", "genes.tsv", "batch_report.txt", "slurm", "all_done", "running"}


def read_manifest(path):
    base = os.path.dirname(os.path.abspath(path))
    with open(path) as fh:   # (file line number, fields); blank and # lines are skipped but still counted
        lines = [(i, l.rstrip("\r\n").split("\t")) for i, l in enumerate(fh, 1)
                 if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        raise SetupError("manifest %s has no gene lines (it needs a header gene, alignment, tree and one line per "
                         "gene)" % path)
    head = [h.strip().lower() for h in lines[0][1]]
    for need in ("gene", "alignment", "tree"):
        if need not in head:
            raise SetupError("manifest needs a header with columns gene, alignment, tree (missing %r)" % need)
    rows = []
    seen = set()
    for i, l in lines[1:]:
        r = dict(zip(head, [x.strip() for x in l]))
        if not r.get("alignment") or not r.get("tree"):
            raise SetupError("manifest line %d has no alignment or tree" % i)
        # unique even ignoring case: on macOS and Windows folders 'g1' and 'G1' are the same folder
        if not r.get("gene") or r["gene"].lower() in seen or "/" in r["gene"] or r["gene"].startswith("."):
            raise SetupError("manifest line %d: gene names must be unique (ignoring case), non-empty, not start with "
                             "'.' and contain no '/': %r" % (i, r.get("gene")))
        if r["gene"].lower() in RESERVED:
            raise SetupError("manifest line %d: %r is a name the batch folder uses itself" % (i, r["gene"]))
        seen.add(r["gene"].lower())
        for k in ("alignment", "tree"):
            if not os.path.isabs(r[k]):
                r[k] = os.path.join(base, r[k])
        rows.append(r)
    return rows


def batch_setup(args, log=print):
    out = os.path.abspath(args.out)
    if os.path.exists(out) and os.listdir(out):
        raise SetupError("%s exists and is not empty; refusing to write over earlier runs" % out)
    rows = read_manifest(args.manifest)
    os.makedirs(out, exist_ok=True)
    genes, failed = [], {}
    for r in rows:
        a = copy.copy(args)
        a.alignment, a.tree, a.out = r["alignment"], r["tree"], os.path.join(out, r["gene"])
        a.test = r.get("test") or args.test
        a.reference = r.get("reference") or args.reference
        if r.get("exclude"):
            a.exclude = [x.strip() for x in r["exclude"].split(",") if x.strip()]
        log("== %s" % r["gene"])
        try:
            setup(a, log=lambda m: log("   " + m))
            genes.append(r["gene"])
        except (SetupError, ValueError, OSError) as e:
            failed[r["gene"]] = str(e)
            log("   SKIPPED: %s" % e)
    with open(os.path.join(out, "batch.json"), "w") as fh:
        json.dump({"manifest": os.path.abspath(args.manifest), "genes": genes, "failed": failed}, fh, indent=1)
    log("batch: %d gene(s) set up, %d skipped (see batch.json)" % (len(genes), len(failed)))
    if not genes:
        raise SetupError("no gene could be set up (paths in the manifest are relative to the manifest's folder); "
                         "remove %s before trying again" % out)
    return out


def batch_report(path, alpha=None, direction=None, logl_tolerance=0.5):
    from .report import summarize
    with open(os.path.join(path, "batch.json")) as fh:
        b = json.load(fh)
    head = ("gene\tstatus\tgate\tobserved_LRT\tN\trank\trank_based_P\tsequential_P\tN_length_matched"
            "\trank_length_matched\trank_based_P_length_matched\tshare_alternatives_shorter"
            "\tchi2_significant_alternatives\texcluded_runs\twarnings")
    lines, L = [head], []
    ranked = []
    seq_skipped = []
    for g in b["genes"]:
        try:
            S = summarize(os.path.join(path, g), alpha, direction, logl_tolerance)
        except Exception as e:   # one unreadable gene must not hide the others
            lines.append("%s\tunreadable: %s" % (g, e) + "\t" * 13)
            continue
        nw = sum(len(r["res"]["warnings"]) for r in S["rows"] if not r["res"]["problem"])
        if not S["evaluable"]:
            lines.append("%s\tgate not evaluable" % g + "\t" * 11 + "\t%d\t%d" % (len(S["excluded"]), nw))
            continue
        if S["N"] == 0 and not S.get("pending"):
            lines.append("%s\tno alternative could be read" % g + "\t" * 11 + "\t%d\t%d" % (len(S["excluded"]), nw))
            continue
        sq = S.get("sequential") or {}
        lines.append("\t".join(str(x) for x in (
            g, "ok", "passed" if S["gate"] else "not passed", "%.4f" % S["lobs"], S["N"], S["rank"],
            "%.4f" % S["rank_p"] if S["rank_p"] is not None else "",
            "%.4f" % sq["P"] if sq.get("P") is not None else "", len(S["length_matched"]), S["lm_rank"],
            "%.4f" % S["lm_rank_p"] if S["lm_rank_p"] is not None else "",
            "%.3f" % S["shorter_share"] if S["shorter_share"] is not None else "",
            len(S["significant"]), len(S["excluded"]), nw)))
        sq = S.get("sequential")
        if S.get("incomplete") or (sq and not sq.get("complete")):
            lines[-1] = lines[-1].replace("\tok\t", "\tincomplete (alternatives not run)\t", 1)
        elif sq:
            # a stopped walk is not a full rank: its 1/(N+1) expectation does not hold, so it is left out
            seq_skipped.append(g)
        else:
            # every evaluable gene counts, whether or not it passed the gate: filtering on the gate selects high
            # observed LRTs, and the chance of ranking first would no longer be 1/(N+1)
            ranked.append((g, S["rank"], S["N"], S["gate"]))
    for g, why in sorted(b.get("failed", {}).items()):
        lines.append("%s\tskipped at setup: %s" % (g, why.replace("\n", " ")) + "\t" * 13)
    L.append("# labelrank batch report: %d gene(s) set up, %d skipped at setup" % (len(b["genes"]), len(b.get("failed", {}))))
    L.append("genes with a complete reading: %d; observed label passed the gate in %d of them" % (
        len(ranked), sum(1 for x in ranked if x[3])))
    if ranked:
        first = sum(1 for _, r, _, _ in ranked if r == 1)
        exp = sum(1.0 / (n + 1) for _, _, n, _ in ranked)
        L.append("ranked first among same-shaped labels in %d of the %d; expected if the label were exchangeable "
                 "with its alternatives: %.2f" % (first, len(ranked), exp))
        L.append("(over all genes with a complete reading, gate passed or not; genes share the tested lineage and "
                 "tree, so they are not independent replicates)")
    if seq_skipped:
        L.append("%d gene(s) set up with sequential stopping are not counted above (their reading is the sequential P; "
                 "a walk that could stop is not a full rank, even when it ran to the end)" % len(seq_skipped))
    return "\n".join(L) + "\n", "\n".join(lines) + "\n"
