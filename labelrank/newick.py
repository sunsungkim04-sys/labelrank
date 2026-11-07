"""Newick trees with branch labels in HyPhy ({label}) or PAML (#1, $1) notation.

Each node stands for the branch above it (the root has no branch). Internal nodes without a name get one
(LRn1, LRn2, ...) when the tree is written for a run, so that the branch sets HyPhy prints can be checked by name
rather than by count: HyPhy's own Node<k> numbering is not something this tool can predict.
"""
import re
import sys


class Node(object):
    __slots__ = ("name", "label", "support", "length", "kids", "parent")

    def __init__(self):
        self.name = ""
        self.label = ""      # "" = unlabelled; HyPhy "{x}" -> "x"; PAML "#1" -> "#1"
        self.support = ""    # a numeric internal "name" (bootstrap etc.) is kept here, never used
        self.length = None   # float or None
        self.kids = []
        self.parent = None

    @property
    def is_tip(self):
        return not self.kids


_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _strip_comments(text):
    """Remove [comments] outside quoted names."""
    out, q, depth = [], False, 0
    for c in text:
        if c == "'" and depth == 0:
            q = not q
        if not q and c == "[":
            depth += 1
            continue
        if not q and c == "]" and depth:
            depth -= 1
            continue
        if depth == 0:
            out.append(c)
    return "".join(out)


class NewickError(ValueError):
    pass


def parse(text):
    """Parse one Newick tree. Accepts quoted names, [comments], {label}, #n / $n, support values and :lengths."""
    s = _strip_comments(text.lstrip("\ufeff")).strip()
    # the parser, topology() and the writers recurse once per tree level; deep (caterpillar) trees need room
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * s.count("(") + 1000))
    if not s.endswith(";"):
        raise NewickError("tree does not end with ';'")
    s = s[:-1]
    pos = [0]
    clade_marks = []   # (node, mark) for PAML '$n' clade labels, expanded after parsing

    def skip_ws():
        while pos[0] < len(s) and s[pos[0]].isspace():
            pos[0] += 1

    def read_name():
        skip_ws()
        if pos[0] < len(s) and s[pos[0]] == "'":
            j = pos[0] + 1
            while True:   # '' inside a quoted name is a literal quote
                j = s.index("'", j)
                if j + 1 < len(s) and s[j + 1] == "'":
                    j += 2
                    continue
                break
            nm = s[pos[0] + 1:j].replace("''", "'")
            pos[0] = j + 1
            return nm
        j = pos[0]
        while j < len(s) and s[j] not in "(),:;{}#$" and not s[j].isspace():
            j += 1
        nm = s[pos[0]:j]
        pos[0] = j
        return nm

    def read_marks(n):
        while True:
            skip_ws()
            if pos[0] >= len(s):
                return
            c = s[pos[0]]
            if c == "{":
                j = s.index("}", pos[0])
                if n.label:
                    raise NewickError("branch %r carries two labels" % n.name)
                n.label = s[pos[0] + 1:j].strip()
                pos[0] = j + 1
            elif c in "#$":
                m = re.match(r"[#$]\d+", s[pos[0]:])
                if not m:
                    raise NewickError("bad PAML mark near position %d" % pos[0])
                if c == "#":
                    if n.label:
                        raise NewickError("branch %r carries two labels" % n.name)
                    n.label = m.group(0)
                else:
                    clade_marks.append((n, "#" + m.group(0)[1:]))
                pos[0] += m.end()
            elif c == ":":
                if n.length is not None:
                    raise NewickError("two branch lengths for one branch near position %d" % pos[0])
                m = re.match(r":\s*([0-9.eE+-]+)", s[pos[0]:])
                if not m:
                    raise NewickError("bad branch length near position %d" % pos[0])
                try:
                    n.length = float(m.group(1))
                except ValueError:
                    raise NewickError("bad branch length %r near position %d" % (m.group(1), pos[0]))
                pos[0] += m.end()
            elif n.kids and not n.support and (c.isdigit() or c == "."):
                # a support value written after the label, as in "){clade}0.726:0.0005"
                m = re.match(r"[0-9.eE+-]+", s[pos[0]:])
                n.support = m.group(0)
                pos[0] += m.end()
            else:
                return

    def node():
        n = Node()
        skip_ws()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            while True:
                k = node()
                k.parent = n
                n.kids.append(k)
                skip_ws()
                c = s[pos[0]]
                pos[0] += 1
                if c == ")":
                    break
                if c != ",":
                    raise NewickError("unexpected %r at position %d" % (c, pos[0] - 1))
        nm = read_name()
        if n.kids and nm and _NUM.match(nm):
            n.support = nm
        else:
            n.name = nm
        read_marks(n)
        return n

    try:
        root = node()
    except IndexError:
        raise NewickError("unbalanced parentheses or truncated tree")
    except ValueError as e:   # an unclosed quote or brace
        if isinstance(e, NewickError):
            raise
        raise NewickError("unclosed quote or brace")
    skip_ws()
    if pos[0] != len(s):
        raise NewickError("trailing text after the tree at position %d" % pos[0])
    for n, mark in clade_marks:
        for x in subtree(n):
            if x.label and x.label != mark:
                raise NewickError("clade mark %s conflicts with label %s on %r" % (mark, x.label, x.name))
            x.label = mark
    root.label = ""   # the root has no branch
    names = [x.name for x in nodes(root) if x.name]
    dup = sorted(set(x for x in names if names.count(x) > 1))
    if dup:
        raise NewickError("duplicate node names: %s" % ", ".join(dup[:5]))
    for t in tips(root):
        if not t.name:
            raise NewickError("a tip has no name")
    return root


def nodes(root):
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(reversed(n.kids))
    return out


def branches(root):
    """Every node except the root, in pre-order: one per branch."""
    return [n for n in nodes(root) if n.parent is not None]


def subtree(n):
    return nodes(n)


def tips(n):
    return [x for x in nodes(n) if x.is_tip]


def name_internal_nodes(root, prefix="LRn"):
    """Give every unnamed internal node a unique name; returns how many were named."""
    taken = set(x.name for x in nodes(root))
    k, named = 0, 0
    for x in nodes(root):
        if x.kids and not x.name:
            while True:
                k += 1
                cand = "%s%d" % (prefix, k)
                if cand not in taken:
                    break
            x.name = cand
            taken.add(cand)
            named += 1
    return named


def topology(n):
    """Canonical string of the rooted, unlabelled shape below (and including) n."""
    if n.is_tip:
        return "t"
    return "(" + ",".join(sorted(topology(k) for k in n.kids)) + ")"


def _fmt_len(x):
    return repr(x) if isinstance(x, float) else str(x)


def write_hyphy(root, labels=None):
    """HyPhy Newick: name{label}:length. labels: optional dict node->label overriding node.label."""
    def lab(n):
        return labels.get(n, n.label) if labels is not None else n.label

    def w(n):
        inner = "(" + ",".join(w(k) for k in n.kids) + ")" if n.kids else ""
        lb = lab(n) if n.parent is not None else ""
        if lb.startswith("#"):
            raise NewickError("PAML label %s cannot be written for HyPhy" % lb)
        return inner + _q(n.name) + ("{%s}" % lb if lb else "") + (
            ":" + _fmt_len(n.length) if n.length is not None else "")
    return w(root) + ";\n"


def write_paml(root, labels=None):
    """PAML Newick without branch lengths: name #1. Unlabelled branches are background (#0)."""
    def lab(n):
        return labels.get(n, n.label) if labels is not None else n.label

    def w(n):
        inner = "(" + ",".join(w(k) for k in n.kids) + ")" if n.kids else ""
        lb = lab(n) if n.parent is not None else ""
        if lb and not lb.startswith("#"):
            raise NewickError("label %r is not a PAML mark (#n)" % lb)
        nm = n.name if n.is_tip else ""
        return inner + nm + (" " + lb if lb and lb != "#0" else "")
    return w(root) + ";\n"


def _q(name):
    if re.search(r"[\s(),:;{}\[\]']", name):
        return "'" + name.replace("'", "''") + "'"
    return name
