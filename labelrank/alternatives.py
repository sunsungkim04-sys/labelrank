"""Same-shaped alternative labels on the same tree.

A test set is split into connected components (branches joined parent-to-child). Each component must be one of
  tip    - a single terminal branch
  clade  - an internal branch together with every branch below it
  stem   - a single internal branch (the branches below it are not in the set)
and its shape is the kind plus the rooted topology of the clade it sits on. An alternative gives the test label to
the same number of components with the same shapes, every branch of which is in the reference set, and that stay
separate components once labelled. The published test branches move to a set that is neither test nor reference,
so the command line (test label, reference label) is the same in every run.
"""
import itertools
import random

from . import newick

UNLABELED = "Unlabeled branches"


class ShapeError(ValueError):
    pass


def components(branch_set):
    """Connected components of a set of branches (nodes), each as (top, members)."""
    bs = set(branch_set)
    tops = [n for n in bs if n.parent not in bs]
    out = []
    for t in tops:
        mem, stack = [], [t]
        while stack:
            x = stack.pop()
            mem.append(x)
            stack.extend(k for k in x.kids if k in bs)
        out.append((t, mem))
    out.sort(key=lambda c: min(_order_key(x) for x in c[1]))
    return out


_ORDER = {}


def _order_key(n):
    return _ORDER.get(id(n), 0)


def index_order(root):
    _ORDER.clear()
    for i, n in enumerate(newick.nodes(root)):
        _ORDER[id(n)] = i


def shape(top, members):
    mem = set(members)
    if top.is_tip:
        return ("tip", "t")
    if mem == set(newick.subtree(top)):
        return ("clade", newick.topology(top))
    if mem == {top}:
        return ("stem", newick.topology(top))
    raise ShapeError("a part of the test set (the branch above %s) is neither a tip, a whole clade, nor a single "
                     "stem branch: %d of the %d branches of that clade carry the test label. Label the whole clade, "
                     "or only its stem" % ("+".join(sorted(t.name for t in newick.tips(top))[:4]) +
                                          ("+..." if len(newick.tips(top)) > 4 else ""), len(mem),
                                          len(newick.subtree(top))))


def label_counts(root):
    """'back 96, whale 35, ..., unlabelled 1': how the tree's branches are labelled, largest first."""
    counts = {}
    for n in newick.branches(root):
        k = n.label or "unlabelled"
        counts[k] = counts.get(k, 0) + 1
    return ", ".join("%s %d" % kv for kv in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def no_alternatives_message(root, n_reference, reference, n_test, length_matched_only=False):
    """Why there is nothing to rank against: a reference too small to hold the test set's shape (typically the
    default reference on a tree labelled all over), or a shape - usually a large clade - that occurs only once."""
    msg = "no same-shaped alternative lies within the reference (%d branch(es)), so there is nothing to rank against" % (
        n_reference)
    if length_matched_only:
        return msg + "; same-shaped sets exist but none is within the length ratio (drop --length-matched-only)"
    if reference == UNLABELED and n_reference < 2 * n_test:
        return msg + (". Branch labels in the tree: %s. If the branches you mean as reference carry a label, pass "
                      "--reference <label>" % label_counts(root))
    return msg + (": no other lineage of the reference has the test set's shape (a clade's exact topology rarely "
                  "recurs on one tree); this tool calibrates labels whose shape repeats, such as a tip or a small "
                  "clade")


def branch_sets(root, test, reference):
    """Test, reference and other branches of the tree as labelled."""
    br = newick.branches(root)
    t = [n for n in br if n.label == test]
    if not t:
        found = sorted(set(n.label for n in br if n.label))
        raise ShapeError("no branch carries the test label %r; labels in the tree: %s" % (
            test, ", ".join(found) if found else "none"))
    if reference == UNLABELED:
        r = [n for n in br if not n.label]
        if not r:
            raise ShapeError("every branch carries a label, and HyPhy 2.5.101 refuses --reference "
                             "'Unlabeled branches' then; name the reference label explicitly")
    else:
        r = [n for n in br if n.label == reference]
        if not r:
            found = sorted(set(n.label for n in br if n.label))
            raise ShapeError("no branch carries the reference label %r; labels in the tree: %s" % (
                reference, ", ".join(found) if found else "none"))
    ts, rs = set(t), set(r)
    other = [n for n in br if n not in ts and n not in rs]
    return t, r, other


def candidates(root, shp, ref_set):
    kind, topo = shp
    out = []
    # HyPhy merges the two branches below a bifurcating root into one; a set using either could lose a branch
    root_kids = set(root.kids) if len(root.kids) == 2 else set()
    for n in newick.branches(root):
        if n not in ref_set or n in root_kids:
            continue
        if kind == "tip":
            if n.is_tip:
                out.append((n,))
        elif n.kids and newick.topology(n) == topo:
            mem = tuple(newick.subtree(n)) if kind == "clade" else (n,)
            if all(x in ref_set for x in mem):
                out.append(mem)
    return out


def resolve_exclude(root, specs, test=None, reference=None):
    """Branches to keep out of every alternative, from --exclude items.

    An item is a tip or node name, or tip names joined by '+' for the clade of their most recent common ancestor;
    the lineage is excluded whole: that branch and every branch below it. Excluded branches stay in their set for
    the fits (they are not relabelled); they are only never given the test label. Call it before internal nodes
    are named, so only names in the user's tree resolve. With test and reference given, an item that covers no
    reference branch is refused (it would be reported as a lineage left out while changing nothing). Returns (set
    of nodes, list of readable descriptions).
    """
    if not specs:
        return set(), []
    ref = set(branch_sets(root, test, reference)[1]) if test is not None else None
    seen = []
    for spec in specs:
        if spec not in seen:
            seen.append(spec)
    specs = seen
    by_name = {}
    for n in newick.nodes(root):
        if n.name:
            by_name.setdefault(n.name, []).append(n)
    out, desc = set(), []
    for spec in specs:
        names = [x.strip() for x in spec.split("+")]
        if not all(names):
            raise ShapeError("--exclude %r has an empty name" % spec)
        hits = []
        for nm in names:
            if nm not in by_name:
                raise ShapeError("--exclude: no tip or node named %r in the tree" % nm)
            if len(by_name[nm]) > 1:
                raise ShapeError("--exclude: the name %r occurs more than once in the tree" % nm)
            hits.append(by_name[nm][0])
        top = hits[0]
        if len(hits) > 1:
            paths = []
            for h in hits:
                p, x = [], h
                while x is not None:
                    p.append(x)
                    x = x.parent
                paths.append(p)
            common = set(paths[0]).intersection(*[set(p) for p in paths[1:]])
            top = next(x for x in paths[0] if x in common)
        if top.parent is None:
            raise ShapeError("--exclude %r covers the whole tree (its common ancestor is the root)" % spec)
        sub = newick.subtree(top)
        n_ref = len(ref & set(sub)) if ref is not None else len(sub)
        if ref is not None and not n_ref:
            raise ShapeError("--exclude %r covers no reference branch, so it would leave nothing out" % spec)
        out.update(sub)
        desc.append("%s (%d %sbranch%s)" % (spec, n_ref, "reference " if ref is not None else "",
                                             "" if n_ref == 1 else "es"))
    return out, desc


def total_length(members):
    return sum(x.length or 0.0 for x in members)


def enumerate_alternatives(root, test, reference, max_alternatives=None, sample_seed=1, enumerate_limit=200000,
                           exclude=None):
    """Returns dict with shapes, the observed set, and a list of alternatives (each a tuple of member nodes).

    All alternatives are kept when there are at most max_alternatives; otherwise a uniform sample of that size,
    drawn with sample_seed, and 'sampled' says so. Branches in exclude (a set of nodes) are never part of one.
    """
    index_order(root)
    t, r, other = branch_sets(root, test, reference)
    comps = components(t)
    shapes = [shape(top, mem) for top, mem in comps]
    rset = set(r)
    n_excl = 0
    if exclude:
        n_excl = len(rset & set(exclude))
        rset -= set(exclude)
    pools = [candidates(root, s, rset) for s in shapes]
    want = sorted(shapes)
    total = 1
    for p in pools:
        total *= max(1, len(p))
    found = set()
    alts = []

    def accept(combo):
        mem = [x for part in combo for x in part]
        if len(set(mem)) != len(mem):
            return
        key = frozenset(mem)
        if key in found:
            return
        cs = components(mem)
        try:
            got = sorted(shape(top, m) for top, m in cs)
        except ShapeError:
            return
        if got != want:
            return
        found.add(key)
        alts.append(tuple(sorted(mem, key=_order_key)))

    sampled = False
    drawn = False
    if any(not p for p in pools):
        pass
    elif total <= enumerate_limit:
        for combo in itertools.product(*pools):
            accept(combo)
    else:
        # too many combinations to list: draw combinations at random; the number that exist is then unknown
        # (only the number of raw combinations, before overlap checks), and the report says so
        if max_alternatives is None:
            raise ShapeError("%d combinations are too many to list; give --max-alternatives (a sample size)" % total)
        rng = random.Random(sample_seed)
        tries = 0
        while len(alts) < max_alternatives and tries < 50 * max_alternatives:
            tries += 1
            accept(tuple(rng.choice(p) for p in pools))
        sampled = drawn = True
    alts.sort(key=lambda a: [_order_key(x) for x in a])
    n_all = len(alts)
    if max_alternatives is not None and len(alts) > max_alternatives:
        rng = random.Random(sample_seed)
        alts = sorted(rng.sample(alts, max_alternatives), key=lambda a: [_order_key(x) for x in a])
        sampled = True
    return {"test": t, "reference": r, "other": other, "components": comps, "shapes": shapes,
            "alternatives": alts, "n_found": None if drawn else n_all, "sampled": sampled, "drawn": drawn,
            "combinations": total, "excluded_reference": n_excl}


def shape_words(kind, topo):
    """('clade', '((t,t),t)') -> 'clade of 3 tips ((t,t),t)'."""
    if kind == "tip":
        return "single tip"
    return "%s of %d tips %s" % ("whole clade" if kind == "clade" else "stem of a clade", topo.count("t"), topo)


def describe(members):
    """Readable name of a branch set: tips by name, an internal branch as (tip1+tip2+...)."""
    out = []
    for x in members:
        if x.is_tip:
            out.append(x.name)
        else:
            out.append("(" + "+".join(sorted(t.name for t in newick.tips(x))) + ")")
    return ",".join(sorted(out))
