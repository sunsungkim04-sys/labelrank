import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import fake_relax as F  # noqa: E402
from labelrank import hyphy_relax as H  # noqa: E402
from labelrank import report as R  # noqa: E402

T, REF, NUI = ["a", "b", "n1"], ["c", "d", "e", "n2"], ["f"]


class TestTolerantJson(unittest.TestCase):
    def test_bare_inf_nan(self):
        j = H.loads_tolerant('{"a": inf, "b": -inf, "c": [nan, -nan, 1], "d": "inf nan", "e": 2}')
        self.assertEqual(j["a"], float("inf"))
        self.assertEqual(j["b"], float("-inf"))
        self.assertTrue(j["c"][0] != j["c"][0])
        self.assertEqual(j["d"], "inf nan")   # strings untouched
        self.assertEqual(j["e"], 2)

    def test_key_names_containing_inf_are_untouched(self):
        j = H.loads_tolerant('{"info": 1, "x_inf": inf, "nan_key": "ok"}')
        self.assertEqual(j["info"], 1)
        self.assertEqual(j["x_inf"], float("inf"))

    def test_relax_json_with_inf_is_read_not_dropped(self):
        js = F.json_text(-1000.0, -1003.0, 6.0, 0.0143, 0.6)
        r = H.read_result(F.out_text(T, REF, NUI), js)
        self.assertAlmostEqual(r["LRT"], 6.0)


class TestReadByLogL(unittest.TestCase):
    def test_reported_direction_from_lower_fit_is_caught(self):
        # the shape of the published case: K > 1 fit reported, the other direction fits better
        o = F.out_text(T, REF, NUI, alt=(-19895.85, 8.90), other=(-19895.66, 0.05), null=-19898.14, p=0.0323,
                       evidence="intensification")
        r = H.read_result(o, None)
        self.assertAlmostEqual(r["LRT"], 4.96, places=6)
        self.assertEqual(H.side(r["K"]), "relax")
        w = " ".join(r["warnings"])
        self.assertIn("opposite sides of K=1", w)
        self.assertIn("reported direction 'intensification' disagrees", w)
        self.assertIn("reported P 0.0323 differs", w)

    def test_other_direction_at_k1_is_not_a_conflict(self):
        o = F.out_text(T, REF, NUI, alt=(-1000.0, 0.65), other=(-1003.61, 1.00), null=-1003.61, p=0.007,
                       evidence="relaxation")
        r = H.read_result(o, None)
        self.assertEqual(r["warnings"], [])
        self.assertAlmostEqual(r["LRT"], 7.22, places=6)

    def test_json_precision_replaces_matching_out_values(self):
        o = F.out_text(T, REF, NUI, alt=(-1000.0, 0.65), other=(-1003.61, 1.00), null=-1003.61)
        js = F.json_text(-1000.00123, -1003.6061, 7.20974, 0.00725, 0.6459, other=(-1003.6088, 1.0))
        r = H.read_result(o, js)
        self.assertAlmostEqual(r["LRT"], 2 * (-1000.00123 + 1003.6061), places=9)
        self.assertAlmostEqual(r["K"], 0.6459)
        self.assertEqual(r["warnings"], [])

    def test_alternative_below_null(self):
        r = H.read_result(F.out_text(T, REF, NUI, alt=(-1005.0, 1.2), null=-1003.0), None)
        self.assertTrue(any("did not converge" in w for w in r["warnings"]))

    def test_sets_and_deleted(self):
        o = H.parse_out(F.out_text(T, REF, NUI, deleted=("n2", "x9")))
        self.assertEqual(o["sets"]["Test"], frozenset(T))
        self.assertEqual(o["deleted"], frozenset({"n2", "x9"}))

    def test_incomplete(self):
        r = H.read_result(None, None)
        self.assertIsNone(r["LRT"])

    def test_refit_after_negative_lrt_is_read_from_json(self):
        o = F.out_text(T, REF, NUI, alt=(-1002.0, 0.6), null=-1000.0, refit=(-996.0, 3.0, -1000.0))
        js = F.json_text(-996.0, -1000.0, 8.0, 0.0047, 3.0, negative_lrt=1)
        r = H.read_result(o, js)
        self.assertAlmostEqual(r["LRT"], 8.0)
        self.assertAlmostEqual(r["K"], 3.0)
        self.assertTrue(any("refitted after a negative LRT" in w for w in r["warnings"]))

    def test_refit_without_json_is_unreadable(self):
        o = F.out_text(T, REF, NUI, alt=(-1002.0, 0.6), null=-1000.0, refit=(-996.0, 3.0, -1000.0))
        self.assertIn("unreadable", H.read_result(o, None))

    def test_near_tie_pairs_logL_with_its_own_k(self):
        # both fits print -1000.00; the other direction (K 1.8) is higher by 0.004 in the JSON
        o = F.out_text(T, REF, NUI, alt=(-1000.0, 0.5), other=(-1000.0, 1.8), null=-1003.0, reset=True)
        js = F.json_text(-999.998, -1003.0, 6.004, 0.0143, 1.8, other=(-999.998, 1.8),
                         unstable=[(-999.998, 1.8), (-1000.002, 0.5)])
        r = H.read_result(o, js)
        self.assertAlmostEqual(r["K"], 1.8)
        self.assertAlmostEqual(r["best_logL"], -999.998)

    def test_chi2(self):
        self.assertAlmostEqual(H.chi2_1_sf(3.841458820694124), 0.05, places=9)


def make_out(tmp, obs_lrts, alt_lrts, n_seeds=3, alt_lengths=None, tweak=None, root_children=("n2", "n9"),
             engine="relax", alt_seeds=1):
    """A fake output folder: observed under n_seeds seeds, one seed per alternative. LRT = 2 * (alt - null)."""
    out = os.path.join(tmp, "o")
    os.makedirs(os.path.join(out, "runs"))
    rows = ["run\trole\tlabel\tseed\tn_branches\tlength\tlength_matched\tbranches"]
    exp = {"Test": T, "Reference": REF, "nuisance": NUI, "root_children": list(root_children)}
    specs = [("observed_s%d" % (i + 1), "observed", i + 1, l, 0.05) for i, l in enumerate(obs_lrts[:n_seeds])]
    for i, l in enumerate(alt_lrts):
        L = alt_lengths[i] if alt_lengths else 0.05
        for sd in range(1, alt_seeds + 1):
            specs.append(("alt%03d_s%d" % (i + 1, sd), "alt%03d" % (i + 1), sd, l, L))
    for run, lab, seed, lrt, L in specs:
        kw = dict(alt=(-1000.0 + lrt / 2.0, 0.5), null=-1000.0)
        if tweak and run in tweak:
            kw.update(tweak[run])
        F.write_run(os.path.join(out, "runs", run), exp, **kw)
        rows.append("%s\t%s\t%s\t%d\t3\t%g\t%s\tx" % (run, "observed" if lab == "observed" else "alternative",
                                                     lab, seed, L, "yes" if 0.025 <= L <= 0.1 else "no"))
    with open(os.path.join(out, "runs.tsv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    with open(os.path.join(out, "setup.json"), "w") as fh:
        json.dump({"engine": engine, "engine_version": "fake", "test": "fg", "reference": "bg",
               "alternatives_written": len(alt_lrts), "alternatives_found": len(alt_lrts), "sampled": False,
                   "seeds": list(range(1, n_seeds + 1)), "alt_seeds": list(range(1, alt_seeds + 1)),
                   "length_ratio": 2.0, "observed_length": 0.05}, fh)
    return out


def _fake_codeml(out, obs, alts, two_ratio=False):
    runs = ["observed_s%d" % (i + 1) for i in range(len(obs))] + ["alt%03d_s1" % (i + 1) for i in range(len(alts))]
    for run, lrt in zip(runs, list(obs) + list(alts)):
        for sub, ll, npar in (("alt", -1000.0 + lrt / 2, 16), ("null", -1000.0, 15)):
            d = os.path.join(out, "runs", run, sub)
            os.makedirs(d)
            with open(os.path.join(d, "mlc"), "w") as fh:
                # the lnL line as PAML 4.10.10 prints it once np reaches 100 (no space after 'np:')
                fh.write("lnL(ntime: 12  np:%d):  %.6f   +0.000000\n" % (npar, ll))
                if two_ratio:
                    fh.write("w (dN/dS) for branches:  0.12000 0.25000\n" if sub == "alt"
                             else "omega (dN/dS) =  0.14000\n")
            if two_ratio:
                with open(os.path.join(d, "codeml.ctl"), "w") as fh:
                    fh.write("        model = %d\n      NSsites = 0\n" % (2 if sub == "alt" else 0))
        with open(os.path.join(out, "runs", run, "expected_sets.json"), "w") as fh:
            json.dump({"Test": T, "Reference": REF, "nuisance": NUI}, fh)
        with open(os.path.join(out, "runs", run, "DONE"), "w") as fh:
            fh.write("ok\t0\t1s\tnow\n")


class TestReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_rank_only_below_20(self):
        out = make_out(self.tmp, [7.2, 7.0, 7.3], [1.0, 0.5, 8.0, 0.2], alt_seeds=3)
        txt = R.report(out)
        self.assertIn("observed LRT 7.3000 (best alternative vs best null fit over 3 seeds)", txt)
        self.assertIn("gate: P <= 0.05 -> PASSED (seeds individually: 3 of 3 at P <= 0.05)", txt)
        self.assertIn("(a) all alternatives: 4 read (of 4 written); LRT >= observed: 1; rank 2 of 5", txt)
        self.assertIn("a reading with N + 1 < 20 gives its rank only", txt)
        self.assertNotIn("rank-based P", txt)
        self.assertNotIn("warning: the observed label used seeds", txt)

    def test_rank_based_p_at_20_or_more(self):
        out = make_out(self.tmp, [9.0, 9.0, 9.0], [0.1 * i for i in range(19)] + [9.5])
        txt = R.report(out)
        self.assertIn("rank 2 of 21; rank-based P = 2/21 = 0.0952", txt)

    def test_ties_count_against_the_observed_label(self):
        out = make_out(self.tmp, [4.0, 4.0, 4.0], [4.0, 1.0])
        self.assertIn("LRT >= observed: 1; rank 2 of 3", R.report(out))

    def test_gate_uses_the_best_fits_over_seeds(self):
        txt = R.report(make_out(self.tmp, [5.0, 1.0, 1.0], [0.5]))
        self.assertIn("observed LRT 5.0000", txt)
        self.assertIn("-> PASSED (seeds individually: 1 of 3", txt)

    def test_gate_not_passed(self):
        txt = R.report(make_out(self.tmp, [1.0, 1.0, 1.0], [0.5]))
        self.assertIn("NOT PASSED", txt)

    def test_direction_requirement(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [0.5])   # K = 0.5 in every fake run
        self.assertIn("NOT PASSED", R.report(out, direction="intens"))
        self.assertIn("-> PASSED", R.report(out, direction="relax"))

    def test_asymmetric_seeds_are_flagged(self):
        txt = R.report(make_out(self.tmp, [7.0, 7.0, 7.0], [0.5], alt_seeds=1))
        self.assertIn("warning: the observed label used seeds [1, 2, 3] but alternatives used [1]", txt)

    def test_length_share_flag(self):
        txt = R.report(make_out(self.tmp, [7.0, 7.0, 7.0], [1.0, 1.0, 1.0, 1.0], alt_lengths=[0.01, 0.01, 0.01, 0.2]))
        self.assertIn("75% of the alternatives are shorter than the observed test set - the observed label is longer", txt)

    def test_unfinished_alternative_withholds_p(self):
        out = make_out(self.tmp, [9.0, 9.0, 9.0], [0.1 * i for i in range(20)])
        os.remove(os.path.join(out, "runs", "alt003_s1", "DONE"))
        txt = R.report(out)
        self.assertIn("INCOMPLETE: 1 alternative run(s) have not finished", txt)
        self.assertNotIn("rank-based P", txt)

    def test_nonfinite_null_is_excluded(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0, 2.0], tweak={"alt001_s1": {"null": float("nan")}})
        txt = R.report(out)
        self.assertIn("excluded: alt001_s1 - non-finite log-likelihood", txt)

    def test_wrong_branch_set_is_excluded(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0, 2.0],
                       tweak={"alt001_s1": {"got_test": ["a", "b"], "got_nuis": ["f", "n1"]}})
        txt = R.report(out)
        self.assertIn("excluded: alt001_s1 - Test set differs", txt)
        self.assertIn("(a) all alternatives: 1 read (of 2 written)", txt)

    def test_merged_root_branch_is_tolerated(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0],
                       tweak={"alt001_s1": {"got_ref": ["c", "d", "e"]}})   # n2 (a root child) never printed
        self.assertNotIn("excluded", R.report(out))

    def test_merged_root_branch_in_the_test_set_excludes(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0], root_children=("n1", "n9"),
                       tweak={"alt001_s1": {"got_test": ["a", "b"]}})
        self.assertIn("excluded: alt001_s1 - test branch n1 was merged away at the root", R.report(out))

    def test_empty_nuisance_set_prints_no_line_and_is_not_excluded(self):
        global NUI
        saved = NUI
        NUI = []
        try:
            out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0])
        finally:
            NUI = saved
        txt = R.report(out)
        self.assertNotIn("excluded", txt)
        self.assertIn("-> PASSED", txt)

    def test_codeml_secondary_uses_the_mixture_p(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [3.0], engine="codeml")
        # chi-square(1) P of LRT 3.0 is 0.083; the 50:50 mixture gives 0.042, significant at 0.05
        os.rename(os.path.join(out, "runs"), os.path.join(out, "runs_relax"))
        _fake_codeml(out, [7.0, 7.0, 7.0], [3.0])
        txt = R.report(out)
        self.assertIn("significant by the chi-square reading (P <= 0.05): 1 of 1", txt)

    def test_codeml_two_ratio_uses_plain_chi_square(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [3.0], engine="codeml")
        os.rename(os.path.join(out, "runs"), os.path.join(out, "runs_relax"))
        _fake_codeml(out, [7.0, 7.0, 7.0], [3.0], two_ratio=True)
        with open(os.path.join(out, "setup.json")) as fh:
            m = json.load(fh)
        m["codeml_model"] = "two-ratio"
        with open(os.path.join(out, "setup.json"), "w") as fh:
            json.dump(m, fh)
        txt = R.report(out)
        # chi-square(1) P of LRT 3.0 is 0.083: not significant (the mixture would have called it 0.042)
        self.assertIn("significant by the chi-square reading (P <= 0.05): 0 of 1", txt)
        self.assertIn("engine: codeml two-ratio", txt)
        self.assertIn("foreground omega,", txt)
        self.assertIn("  alt001  0.05000  y  3.000  0.0833  0.250  x", txt)
        from labelrank import codeml as C
        r = C.read_run_dir(os.path.join(out, "runs", "alt001_s1"))
        self.assertEqual((r["model"], r["K"]), ("two-ratio", 0.25))
        self.assertAlmostEqual(r["p"], 0.0832645166635504, places=9)

    def test_deleted_test_branch_excludes_deleted_reference_is_noted(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0],
                       tweak={"alt001_s1": {"deleted": ("a",)}, "observed_s1": {"deleted": ("c",)}})
        txt = R.report(out)
        self.assertIn("excluded: alt001_s1 - test branch(es) deleted as zero-length: a", txt)
        self.assertIn("1 reference branch(es) deleted as zero-length (names as in the run's tree.nwk): c", txt)

    def test_missing_gate_run_makes_gate_not_evaluable(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0])
        os.remove(os.path.join(out, "runs", "observed_s2", "relax.out"))
        self.assertIn("gate not evaluable: 2 of 3 observed runs are usable (0 not run yet, 1 excluded above)", R.report(out))

    def test_logl_spread_warning(self):
        out = make_out(self.tmp, [7.0, 7.0, 9.0], [1.0])
        self.assertIn("ABOVE TOLERANCE", R.report(out))
        self.assertNotIn("ABOVE TOLERANCE", R.report(out, logl_tolerance=2.0))

    def test_length_matched_coprimary(self):
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [8.0, 8.0, 1.0], alt_lengths=[0.2, 0.05, 0.05])
        self.assertIn("2 alternatives; LRT >= observed: 1; rank 2 of 3", R.report(out))

    def test_lower_logL_direction_is_not_used_for_rank(self):
        # alternative's reported fit is weak, its other-direction fit is strong: the strong one counts
        out = make_out(self.tmp, [7.0, 7.0, 7.0], [1.0],
                       tweak={"alt001_s1": {"alt": (-999.5, 3.0), "other": (-996.0, 0.2), "null": -1000.0}})
        txt = R.report(out)
        self.assertIn("LRT >= observed: 1", txt)
        self.assertIn("opposite sides of K=1", txt)


if __name__ == "__main__":
    unittest.main()
