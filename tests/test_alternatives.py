import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from labelrank import alternatives as A  # noqa: E402
from labelrank import newick  # noqa: E402

# 8 tips; (G,H) is the labelled cherry; one 3-tip clade (A,(B,C)) with the same number of branches as a cherry
TREE = "(((A:1,(B:1,C:1):1):1,(D:1,E:1):1):1,((F:1,X:1):1,(G{fg}:1,H{fg}:1){fg}:1):1);"


def alts(tree, test="fg", ref=A.UNLABELED, **kw):
    t = newick.parse(tree)
    newick.name_internal_nodes(t)
    return t, A.enumerate_alternatives(t, test, ref, **kw)


class TestShapes(unittest.TestCase):
    def test_cherry_alternatives_are_cherries_only(self):
        _, r = alts(TREE)
        self.assertEqual(r["shapes"], [("clade", "(t,t)")])
        got = sorted(A.describe(a) for a in r["alternatives"])
        self.assertEqual(got, ["(B+C),B,C", "(D+E),D,E", "(F+X),F,X"])   # not the 3-tip clade

    def test_tip(self):
        _, r = alts("((A,B{fg}),(C,D));")
        self.assertEqual(sorted(A.describe(a) for a in r["alternatives"]), ["A", "C", "D"])

    def test_alternatives_stay_inside_the_reference_set(self):
        _, r = alts("((A{bg},B{bg}){bg},((C{o},D{bg}){bg},E{fg}){bg});", ref="bg")
        self.assertEqual(sorted(A.describe(a) for a in r["alternatives"]), ["A", "B", "D"])   # C is labelled 'o'

    def test_stem(self):
        _, r = alts("(((A,B){fg},(C,D)),((E,F),G));")
        self.assertEqual(r["shapes"], [("stem", "(t,t)")])
        self.assertEqual(sorted(A.describe(a) for a in r["alternatives"]), ["(C+D)", "(E+F)"])

    def test_two_tips_in_two_components(self):
        _, r = alts("((A{fg},B),(C{fg},(D,E)));")
        # pairs of reference tips: B, D, E -> 3 pairs
        self.assertEqual(len(r["alternatives"]), 3)
        self.assertTrue(all(len(a) == 2 for a in r["alternatives"]))

    def test_adjacent_pick_that_would_merge_is_refused(self):
        # two stems: the stem of (C,D) and the stem of ((C,D),E) touch, so together they are one component
        _, r = alts("((((A,B){fg},X){fg2}),(((C,D),E),(F,G)));")
        for a in r["alternatives"]:
            self.assertEqual(len(A.components(a)), 1)

    def test_partial_clade_is_refused(self):
        t = newick.parse("(((A{fg},B){fg},C),D);")
        with self.assertRaises(A.ShapeError):
            A.enumerate_alternatives(t, "fg", A.UNLABELED)

    def test_unlabeled_reference_on_fully_labelled_tree(self):
        t = newick.parse("((A{fg},B{bg}){bg},(C{bg},D{bg}){bg});")
        with self.assertRaises(A.ShapeError):
            A.enumerate_alternatives(t, "fg", A.UNLABELED)

    def test_sampling_is_reproducible_and_capped(self):
        tree = "(" + ",".join("T%d" % i for i in range(30)) + ",Q{fg});"
        _, r1 = alts(tree, max_alternatives=10, sample_seed=7)
        _, r2 = alts(tree, max_alternatives=10, sample_seed=7)
        self.assertEqual(len(r1["alternatives"]), 10)
        self.assertTrue(r1["sampled"])
        self.assertEqual(r1["n_found"], 30)
        self.assertEqual([A.describe(a) for a in r1["alternatives"]], [A.describe(a) for a in r2["alternatives"]])


if __name__ == "__main__":
    unittest.main()
