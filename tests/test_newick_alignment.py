import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from labelrank import alignment, newick  # noqa: E402


class TestNewick(unittest.TestCase):
    def test_hyphy_labels_support_and_lengths(self):
        t = newick.parse("((A{fg}:0.1,B{fg}:0.2){fg}0.95:0.05,(C:0.3,D{bg}:0.4)88:0.01);")
        lab = {n.name or n.support: n.label for n in newick.branches(t)}
        self.assertEqual(lab["A"], "fg")
        self.assertEqual(lab["0.95"], "fg")   # support written after the label
        self.assertEqual(lab["C"], "")
        self.assertEqual(lab["88"], "")
        self.assertAlmostEqual(sum(n.length for n in newick.branches(t)), 1.06)

    def test_paml_marks_and_clade_mark(self):
        t = newick.parse("((A,B)$1,(C #1,D));")
        lab = {x.name: x.label for x in newick.nodes(t) if x.name}
        self.assertEqual(lab["A"], "#1")
        self.assertEqual(lab["B"], "#1")
        self.assertEqual(lab["C"], "#1")
        self.assertEqual(lab["D"], "")
        self.assertEqual(newick.write_paml(t), "((A #1,B #1) #1,(C #1,D));\n")

    def test_round_trip_keeps_labels_and_lengths(self):
        src = "(('sp one'{x}:1e-09,B:0.000000006)n1{x}:0.5,C{y}:0.25);"
        t = newick.parse(src)
        t2 = newick.parse(newick.write_hyphy(t))
        a = [(n.name, n.label, n.length) for n in newick.nodes(t)]
        b = [(n.name, n.label, n.length) for n in newick.nodes(t2)]
        self.assertEqual(a, b)

    def test_name_internal_nodes_avoids_existing_names(self):
        t = newick.parse("((A,B),(C,D)LRn1);")
        self.assertEqual(newick.name_internal_nodes(t), 2)
        names = [n.name for n in newick.nodes(t)]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn("LRn1", names)

    def test_topology_is_order_free(self):
        a = newick.parse("((A,(B,C)),D);")
        b = newick.parse("(D,((C,B),A));")
        self.assertEqual(newick.topology(a), newick.topology(b))
        self.assertNotEqual(newick.topology(a.kids[0]), newick.topology(newick.parse("((A,B),(C,D));").kids[0]))

    def test_errors(self):
        for bad in ("(A,B)", "(A,(B,C);", "(A,A);", "(A{x}{y},B);"):
            with self.assertRaises(newick.NewickError):
                newick.parse(bad)


def _fa(recs):
    f = tempfile.NamedTemporaryFile("wb", suffix=".fa", delete=False)
    f.write("".join(">%s\r\n%s\r\n" % r for r in recs).encode())
    f.close()
    return f.name


class TestAlignment(unittest.TestCase):
    def test_clean(self):
        recs = alignment.read_fasta(_fa([("A", "ATGAAA---"), ("B", "ATGAAG---")]))
        err, warn, info = alignment.check(recs, ["A", "B"])
        self.assertEqual(err, [])
        self.assertEqual(info["codons"], 3)

    def test_crlf_is_stripped(self):
        recs = alignment.read_fasta(_fa([("A", "ATGAAA")]))
        self.assertEqual(recs, [("A", "ATGAAA")])

    def test_internal_and_terminal_stops_are_errors(self):
        recs = alignment.read_fasta(_fa([("A", "ATGTGAAAA"), ("B", "ATGAAATAG")]))
        err, _, _ = alignment.check(recs, ["A", "B"])
        self.assertEqual(len(err), 1)
        self.assertIn("2 stop codon(s) in 2 of 2 sequences (1 terminal)", err[0])

    def test_stop_out_of_frame_is_not_a_stop(self):
        recs = alignment.read_fasta(_fa([("A", "ATAAAA")]))   # TAA at offset 1 is not in frame
        self.assertEqual(alignment.check(recs, ["A"])[0], [])

    def test_length_not_multiple_of_three(self):
        recs = alignment.read_fasta(_fa([("A", "ATGAA"), ("B", "ATGAA")]))
        self.assertTrue(any("not a multiple of 3" in e for e in alignment.check(recs, ["A", "B"])[0]))

    def test_names(self):
        recs = alignment.read_fasta(_fa([("A|1", "ATG"), ("B", "ATG")]))
        err = alignment.check(recs, ["A|1", "B", "C"])[0]
        self.assertTrue(any("HyPhy rewrites" in e for e in err))
        self.assertTrue(any("not in the alignment" in e for e in err))


if __name__ == "__main__":
    unittest.main()
