"""Codon alignment checks, run before anything is launched.

The checks are the ones that otherwise surface hours later as dead runs: an in-frame stop codon or a length that is
not a multiple of three makes every HyPhy run fail (in one real case, all 16 runs of a test), and names HyPhy would
rewrite break the branch-set check.
"""
import re

STOPS = {"TAA", "TAG", "TGA"}   # universal code; other codes pass their own set (gencode.py)
_OK_CHARS = set("ACGTUNRYKMSWBDHV-?.")
_HYPHY_SAFE = re.compile(r"^[A-Za-z0-9_]+$")


_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG",
          "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _expand(c):
    return [a + b + d for a in _IUPAC[c[0]] for b in _IUPAC[c[1]] for d in _IUPAC[c[2]]]


def read_fasta(path):
    seqs, order, name = {}, [], None
    with open(path, "rb") as fh:
        raw = fh.read().decode("utf-8", errors="replace").lstrip("\ufeff")
    for line in raw.splitlines():
        line = line.strip().replace("\r", "")
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            if name in seqs:
                raise ValueError("duplicate sequence name %r" % name)
            seqs[name] = []
            order.append(name)
        else:
            if name is None:
                raise ValueError("sequence data before the first '>' header")
            seqs[name].append(line.replace(" ", ""))
    return [(n, "".join(seqs[n]).upper()) for n in order]


def check(records, tip_names, stops=STOPS):
    """Returns (errors, warnings, info). Errors stop the tool before any run is written."""
    errors, warnings = [], []
    if not records:
        return ["the alignment has no sequences"], [], {}
    lengths = set(len(s) for _, s in records)
    if len(lengths) != 1:
        errors.append("sequences differ in length (%s): not an alignment" % ", ".join(str(x) for x in sorted(lengths)))
    L = max(lengths)
    if L % 3:
        errors.append("alignment length %d is not a multiple of 3" % L)
    bad_chars = {}
    stop_hits = []
    for n, s in records:
        odd = set(s) - _OK_CHARS
        if odd:
            bad_chars[n] = "".join(sorted(odd))
        s2 = s.replace("U", "T")
        ncod = len(s2) // 3
        for i in range(ncod):
            c = s2[3 * i:3 * i + 3]
            if c in stops:
                stop_hits.append((n, i + 1, c, i == ncod - 1 or set(s2[3 * i + 3:]) <= set("-?.")))
        if set(s) <= set("-?.N"):
            warnings.append("sequence %r has no determined bases" % n)
    amb = []
    for n, s in records:
        s2 = s.replace("U", "T")
        for i in range(len(s2) // 3):
            c = s2[3 * i:3 * i + 3]
            if c not in stops and set(c) <= set("ACGTRYKMSWBDHVN") and set(c) - set("ACGT") and all(
                    x in stops for x in _expand(c)):
                amb.append((n, i + 1, c))
    if amb:
        warnings.append("%d ambiguous codon(s) can only resolve to stop codons, e.g. %s codon %d %s; HyPhy may "
                        "reject or mask them" % (len(amb), amb[0][0], amb[0][1], amb[0][2]))
    gapped = sum(1 for n, s in records for i in range(len(s) // 3)
                 if 0 < sum(ch in "-." for ch in s[3 * i:3 * i + 3]) < 3)
    if gapped:
        warnings.append("%d codon(s) are partly gapped (e.g. 'A-G'); check that the alignment is in frame" % gapped)
    if bad_chars:
        errors.append("non-nucleotide characters in %d sequence(s), e.g. %s: %s" % (
            len(bad_chars), next(iter(bad_chars)), next(iter(bad_chars.values()))))
    if stop_hits:
        seqs_hit = sorted(set(h[0] for h in stop_hits))
        terminal = sum(1 for h in stop_hits if h[3])
        ex = "; ".join("%s codon %d %s%s" % (n, i, c, " (terminal)" if t else "") for n, i, c, t in stop_hits[:5])
        errors.append("%d stop codon(s) in %d of %d sequences (%d terminal); HyPhy rejects these and codeml "
                      "silently turns their columns into missing data. Remove terminal stops, check the reading "
                      "frame and the genetic code (--code). First: %s" % (
                          len(stop_hits), len(seqs_hit), len(records), terminal, ex))
    names = [n for n, _ in records]
    unsafe = [n for n in names if not _HYPHY_SAFE.match(n)]
    if unsafe:
        errors.append("%d name(s) contain characters HyPhy rewrites (only letters, digits and _ are kept), "
                      "e.g. %r; rename them in both alignment and tree" % (len(unsafe), unsafe[0]))
    a, t = set(names), set(tip_names)
    if a - t:
        errors.append("%d sequence(s) not in the tree, e.g. %r" % (len(a - t), sorted(a - t)[0]))
    if t - a:
        errors.append("%d tree tip(s) not in the alignment, e.g. %r" % (len(t - a), sorted(t - a)[0]))
    seen = {}
    for n, s in records:
        seen.setdefault(s, []).append(n)
    ndup = sum(len(v) - 1 for v in seen.values() if len(v) > 1)
    if ndup:
        warnings.append("%d sequence(s) identical to another (kept; they add no information)" % ndup)
    if L // 3 == 0:
        errors.append("the alignment has no codons")
    elif all(set(s) <= set("-?.N") for _, s in records):
        errors.append("no sequence has a determined base")
    info = {"sequences": len(records), "length_nt": L, "codons": L // 3}
    return errors, warnings, info
