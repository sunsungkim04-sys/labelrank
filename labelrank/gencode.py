"""Genetic codes, by NCBI transl_table number: the stop codons the alignment check uses, HyPhy's --code name and
codeml's icode.

The stop sets are NCBI's. HyPhy's names are the ones in its chooseGeneticCode.def (2.5.101). codeml's icode
numbers follow PAML's documented list (icode 6 is NCBI 9, 7 is 10, 8 is 12, 9 is 13, 10 is 15); feeding codeml
4.10.10 every codon under every icode gave exactly these stop sets (which cannot tell apart the codes whose stops
are only TAA and TAG). Table 11 (bacterial, archaeal,
plastid) assigns codons as table 1 does and differs only in start codons, so it runs as the universal code.
codeml has no flatworm mitochondrial code (14).
"""

#  table: (HyPhy name, codeml icode or None, stop codons)
TABLES = {
    1: ("Universal", 0, ("TAA", "TAG", "TGA")),
    2: ("Vertebrate-mtDNA", 1, ("TAA", "TAG", "AGA", "AGG")),
    3: ("Yeast-mtDNA", 2, ("TAA", "TAG")),
    4: ("Mold-Protozoan-mtDNA", 3, ("TAA", "TAG")),
    5: ("Invertebrate-mtDNA", 4, ("TAA", "TAG")),
    6: ("Ciliate-Nuclear", 5, ("TGA",)),
    9: ("Echinoderm-mtDNA", 6, ("TAA", "TAG")),
    10: ("Euplotid-Nuclear", 7, ("TAA", "TAG")),
    11: ("Universal", 0, ("TAA", "TAG", "TGA")),
    12: ("Alt-Yeast-Nuclear", 8, ("TAA", "TAG", "TGA")),
    13: ("Ascidian-mtDNA", 9, ("TAA", "TAG")),
    14: ("Flatworm-mtDNA", None, ("TAG",)),
    15: ("Blepharisma-Nuclear", 10, ("TAA", "TGA")),
}


def release(version_string):
    """The release number from an engine's own banner ('HYPHY 2.5.101(MP) ...', 'PAML 4.10.10'), else None.
    Only the banner is read, never a path, so an install folder named after a release cannot pass for it."""
    import re
    m = re.match(r"\s*(?:HYPHY|PAML)\s+(\d+(?:\.\d+)+)", version_string or "", re.I)
    return m.group(1) if m else None


def get(table, engine):
    """(HyPhy name, icode, stop set) for an NCBI table; ValueError when the engine has no such code."""
    if table not in TABLES:
        raise ValueError("genetic code %r is not supported; use an NCBI transl_table number: %s" % (
            table, ", ".join(str(t) for t in sorted(TABLES))))
    name, icode, stops = TABLES[table]
    if engine == "codeml" and icode is None:
        raise ValueError("codeml has no genetic code for NCBI table %d (%s)" % (table, name))
    return name, icode, set(stops)
