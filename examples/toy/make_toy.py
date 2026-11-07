"""Regenerate the toy inputs (deterministic). A crude codon simulation with the same purifying selection on every
branch: 85% of changes keep the amino acid. No lineage is treated differently, so no label should stand out.

    python3 make_toy.py      # writes aln.fasta, tree_relax.nwk, tree_codeml.nwk next to this file
"""
import os
import random

BASES = "ACGT"
SENSE = [a + b + c for a in BASES for b in BASES for c in BASES if a + b + c not in ("TAA", "TAG", "TGA")]
TREE = "((((t1:0.08,t2:0.06):0.03,(t3:0.07,t4:0.09):0.02):0.04,((t5:0.10,t6:0.08):0.05,t7:0.12):0.03):0.02," \
       "((t8:0.05,t9:0.07):0.04,t10:0.11):0.02);"
NCOD = 300
AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")   # standard code, TCAG order
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate("TCAG") for j, b in enumerate("TCAG")
        for k, c in enumerate("TCAG")}
SYN = {c: [d for d in SENSE if CODE[d] == CODE[c] and d != c] for c in SENSE}


def evolve(seq, length, rng):
    seq = list(seq)
    for _ in range(sum(1 for _ in range(len(seq)) if rng.random() < length)):
        i = rng.randrange(len(seq))
        if SYN[seq[i]] and rng.random() < 0.85:
            seq[i] = rng.choice(SYN[seq[i]])
        else:
            seq[i] = rng.choice(SENSE)
    return seq


def main():
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from labelrank import newick
    rng = random.Random(20261026)
    root = newick.parse(TREE)
    seqs = {}

    def down(n, s):
        if n.is_tip:
            seqs[n.name] = "".join(s)
        for k in n.kids:
            down(k, evolve(s, k.length, rng))
    down(root, [rng.choice(SENSE) for _ in range(NCOD)])
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "aln.fasta"), "w") as fh:
        for n in sorted(seqs, key=lambda x: int(x[1:])):
            fh.write(">%s\n%s\n" % (n, seqs[n]))
    fg = TREE.replace("t5:", "t5{Foreground}:").replace("t6:0.08)", "t6{Foreground}:0.08){Foreground}")
    with open(os.path.join(here, "tree_relax.nwk"), "w") as fh:
        fh.write(fg + "\n")
    with open(os.path.join(here, "tree_codeml.nwk"), "w") as fh:
        fh.write(fg.replace("{Foreground}", " #1") + "\n")


if __name__ == "__main__":
    main()
