"""Genetic codes, lineage exclusion, the length-LRT correlation, the version warning and the SLURM array script
run as SLURM would run one task."""
import io
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from labelrank import alignment, alternatives as A, cli, gencode, newick  # noqa: E402
from labelrank.report import spearman  # noqa: E402
from test_end_to_end import _caterpillar, _fasta, _tips  # noqa: E402

TREE = "(((A:1,(B:1,C:1):1):1,(D:1,E:1):1):1,((F:1,X:1):1,(G{fg}:1,H{fg}:1){fg}:1):1);"


class GeneticCode(unittest.TestCase):
    def test_vertebrate_mito_stops(self):
        name, icode, stops = gencode.get(2, "relax")
        self.assertEqual((name, icode), ("Vertebrate-mtDNA", 1))
        self.assertEqual(stops, {"TAA", "TAG", "AGA", "AGG"})

    def test_tga_is_sense_in_mito_and_aga_is_stop(self):
        recs = [("a", "ATGTGAAAA"), ("b", "ATGTGAAAA")]
        self.assertEqual(alignment.check(recs, ["a", "b"], gencode.get(2, "relax")[2])[0], [])
        self.assertTrue(alignment.check(recs, ["a", "b"])[0])          # universal: TGA is a stop
        recs = [("a", "ATGAGAAAA"), ("b", "ATGAAAAAA")]
        self.assertTrue(alignment.check(recs, ["a", "b"], gencode.get(2, "relax")[2])[0])

    def test_codeml_has_no_flatworm_code_and_unknown_tables_fail(self):
        self.assertRaises(ValueError, gencode.get, 14, "codeml")
        self.assertEqual(gencode.get(14, "relax")[0], "Flatworm-mtDNA")
        self.assertRaises(ValueError, gencode.get, 7, "relax")


class Exclude(unittest.TestCase):
    def setUp(self):
        self.t = newick.parse(TREE)
        newick.name_internal_nodes(self.t)

    def test_clade_by_tips(self):
        ex, desc = A.resolve_exclude(self.t, ["B+C", "F"])
        r = A.enumerate_alternatives(self.t, "fg", A.UNLABELED, exclude=ex)
        self.assertEqual(sorted(A.describe(a) for a in r["alternatives"]), ["(D+E),D,E"])
        self.assertEqual(r["excluded_reference"], 4)   # B, C, their stem, and F
        self.assertEqual(desc, ["B+C (3 branches)", "F (1 branch)"])
        ex, desc = A.resolve_exclude(self.t, ["B+C", "B+C"], "fg", A.UNLABELED)
        self.assertEqual(desc, ["B+C (3 reference branches)"])

    def test_a_wider_clade_excludes_everything_below(self):
        ex, _ = A.resolve_exclude(self.t, ["A+E"])
        r = A.enumerate_alternatives(self.t, "fg", A.UNLABELED, exclude=ex)
        self.assertEqual(sorted(A.describe(a) for a in r["alternatives"]), ["(F+X),F,X"])

    def test_bad_names(self):
        self.assertRaises(A.ShapeError, A.resolve_exclude, self.t, ["nope"])
        self.assertRaises(A.ShapeError, A.resolve_exclude, self.t, ["A+H"])   # the root
        self.assertRaises(A.ShapeError, A.resolve_exclude, self.t, ["B+"])
        self.assertRaises(A.ShapeError, A.resolve_exclude, self.t, ["G+H"], "fg", A.UNLABELED)   # test lineage
        t = newick.parse(TREE)   # before naming, as setup calls it: the tool's own LRn names do not resolve
        self.assertRaises(A.ShapeError, A.resolve_exclude, t, ["LRn3"])


class Versions(unittest.TestCase):
    def test_release_comes_from_the_banner_only(self):
        from labelrank.gencode import release
        self.assertEqual(release("HYPHY 2.5.101(MP) for Linux"), "2.5.101")
        self.assertEqual(release("PAML 4.10.10"), "4.10.10")
        self.assertIsNone(release("codeml (version not reported)"))
        self.assertIsNone(release("unknown (/opt/hyphy-2.5.101/bin/hyphy: not found)"))

    def test_codeml_banner_is_read_and_a_silent_codeml_is_not_trusted(self):
        from labelrank.codeml import codeml_version
        tmp = tempfile.mkdtemp(prefix="paml4.10.10-")
        try:
            fake = os.path.join(tmp, "codeml")
            with open(fake, "w") as fh:
                fh.write("#!/bin/sh\necho 'CODONML in paml version 4.10.10, 29 Jan 2026'\n")
            os.chmod(fake, 0o755)
            self.assertEqual(codeml_version(fake), "PAML 4.10.10")
            with open(fake, "w") as fh:
                fh.write("#!/bin/sh\necho nothing\n")
            self.assertEqual(codeml_version(fake), "codeml (version not reported)")
        finally:
            shutil.rmtree(tmp)


class Spearman(unittest.TestCase):
    def test_values(self):
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertAlmostEqual(spearman([1, 2, 2, 3], [1, 2, 3, 4]), 0.9486832980505138)
        self.assertIsNone(spearman([1, 1, 1], [1, 2, 3]))
        self.assertIsNone(spearman([1, 2], [1, 2]))

    def test_noise_level_lrts_tie_at_zero(self):
        from labelrank.report import length_rho
        # the codeml A_siren case: one real LRT, the rest 0 or optimizer noise of 1e-6 to 1e-4
        lengths = [0.12, 0.09, 0.037, 0.020, 0.039, 0.0025, 0.028, 0.0045, 0.011]
        lrts = [2.4e-4, 2e-6, 6e-6, 0.0, 0.758, 0.0, 8e-6, 0.0, 0.0]
        self.assertGreater(spearman(lengths, lrts), 0.7)       # what ranking the noise gives
        self.assertEqual(length_rho(lengths, lrts), (None, 1))
        rho, above = length_rho([1, 2, 3, 4], [0.005, 0.5, 1.0, 2.0])
        self.assertEqual(above, 3)
        self.assertAlmostEqual(rho, 1.0)


class _Fixture(unittest.TestCase):
    """A 25-tip caterpillar with tip T12 labelled, a stand-in hyphy and RELAX.bf."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(suffix=" with space")
        self.hy = os.path.join(self.tmp, "hyphy")
        with open(self.hy, "w") as fh:
            fh.write("#!/bin/sh\nexec '%s' '%s' \"$@\"\n" % (sys.executable, os.path.join(HERE, "fake_hyphy.py")))
        os.chmod(self.hy, os.stat(self.hy).st_mode | stat.S_IEXEC)
        self.bf = os.path.join(self.tmp, "bf", "RELAX.bf")
        os.makedirs(os.path.join(self.tmp, "bf", "modules"))
        with open(self.bf, "w") as fh:
            fh.write("// stand-in\n")
        self.names = _tips(25)
        self.tree = os.path.join(self.tmp, "t.nwk")
        with open(self.tree, "w") as fh:
            fh.write(_caterpillar(self.names, "T12"))
        self.aln = os.path.join(self.tmp, "a.fa")
        _fasta(self.names, self.aln)
        os.environ["FAKE_HIGH"] = "T12"

    def tearDown(self):
        shutil.rmtree(self.tmp)
        for k in ("FAKE_HIGH", "FAKE_VERSION", "FAKE_LRT"):
            os.environ.pop(k, None)

    def cli(self, *args):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cli.main(list(args))
        return rc, buf.getvalue()

    def setup(self, out, *extra):
        return self.cli("setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", out,
                        "--hyphy", self.hy, "--relax-bf", self.bf, "--seeds", "1", *extra)



class Pipeline(_Fixture):
    def test_noise_does_not_order_labels(self):
        # observed 0.0002 above alternatives at 0 (or at -0.02, a fit below its null): all are LRT 0, tied
        for default in ("0", "-0.02"):
            os.environ["FAKE_LRT"] = "T12=0.0002,*=%s" % default
            out = os.path.join(self.tmp, "n" + default)
            self.setup(out)
            self.cli("run", out, "--jobs", "4")
            txt = self.cli("report", out)[1]
            self.assertIn("LRT >= observed: 23; rank 24 of 24; rank-based P = 24/24 = 1.0000", txt)
            if default != "0":
                self.assertIn("23 label(s) whose best alternative fit is below the best null fit", txt)
        os.environ.pop("FAKE_LRT")

    def test_noise_just_above_the_floor_still_ties(self):
        os.environ["FAKE_LRT"] = "T12=0.012,*=0.009"
        out = os.path.join(self.tmp, "b")
        self.setup(out)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("LRT >= observed: 23; rank 24 of 24", self.cli("report", out)[1])
        os.environ.pop("FAKE_LRT")

    def test_equal_lengths_are_not_counted_shorter(self):
        # runs.tsv carries lengths to 9 significant digits; a tip of the same length must not count as shorter
        with open(self.tree, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace(":0.1", ":0.1234567890123"))
        out = os.path.join(self.tmp, "e")
        self.setup(out)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("length: 0% of the alternatives are shorter", self.cli("report", out)[1])

    def test_genes_tsv_rows_match_the_header_when_the_gate_is_not_evaluable(self):
        man = os.path.join(self.tmp, "genes.tsv")
        with open(man, "w") as fh:
            fh.write("gene\talignment\ttree\ng1\ta.fa\tt.nwk\ng2\ta.fa\tt.nwk\n")
        out = os.path.join(self.tmp, "b")
        self.cli("batch-setup", "--manifest", man, "--test", "FG", "--out", out, "--hyphy", self.hy,
                 "--relax-bf", self.bf, "--seeds", "1")
        self.cli("run", out, "--jobs", "4")
        with open(os.path.join(out, "g1", "runs", "observed_s1", "DONE"), "w") as fh:
            fh.write("FAIL stand-in failure\n")
        rc, txt = self.cli("batch-report", out)
        with open(os.path.join(out, "genes.tsv")) as fh:
            rows = fh.read().splitlines()
        self.assertTrue(any(r.startswith("g1\tgate not evaluable") for r in rows), rows)
        self.assertEqual(set(len(r.split("\t")) for r in rows), {15})

    def test_no_alternatives_is_an_error(self):
        # the published SCNN1 trees label every branch but one: the default reference is then that one branch
        import re
        txt = _caterpillar(self.names, "T12").replace("):0.05", "){back}:0.05").replace("){back}", ")", 1)
        txt = re.sub(r"(T\d+):0\.1", lambda m: m.group(0) if m.group(1) == "T12" else
                     m.group(1) + "{back}:0.1", txt)   # one internal branch stays unlabelled; the test is a tip
        with open(self.tree, "w") as fh:
            fh.write(txt)
        rc, out = self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "FG")
        self.assertEqual(rc, 2, out)
        self.assertIn("branch labels: back ", out)
        self.assertIn("unlabelled 1", out)
        self.assertIn("ERROR: no same-shaped alternative lies within the reference (1 branch(es))", out)
        o = os.path.join(self.tmp, "o")
        from contextlib import redirect_stderr
        err = io.StringIO()
        with redirect_stderr(err):
            rc, out = self.setup(o)
        self.assertEqual(rc, 2, out)
        self.assertIn("pass --reference <label>", err.getvalue())
        self.assertFalse(os.path.exists(os.path.join(o, "runs")))
        rc, out = self.setup(o, "--reference", "back")
        self.assertEqual(rc, 0, out)

    def test_shape_that_occurs_once(self):
        # a caterpillar has one cherry: label it as a clade and no other lineage has that shape
        n = self.names
        cherry = "%s:0.1,%s:0.1):0.05" % (n[0], n[1])
        txt = _caterpillar(n, None).replace(cherry, "%s{FG}:0.1,%s{FG}:0.1){FG}:0.05" % (n[0], n[1]), 1)
        with open(self.tree, "w") as fh:
            fh.write(txt)
        rc, out = self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "FG")
        self.assertEqual(rc, 2, out)
        self.assertIn("no other lineage of the reference has the test set's shape", out)
        self.assertNotIn("--reference <label>", out)

    def test_code_reaches_hyphy_and_setup_json(self):
        out = os.path.join(self.tmp, "o")
        rc, txt = self.setup(out, "--code", "2")
        self.assertEqual(rc, 0, txt)
        with open(os.path.join(out, "runs", "observed_s1", "cmd.sh")) as fh:
            self.assertIn("--code 'Vertebrate-mtDNA'", fh.read())
        with open(os.path.join(out, "setup.json")) as fh:
            self.assertEqual(json.load(fh)["genetic_code"], 2)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("genetic code: NCBI table 2 (Vertebrate-mtDNA)", self.cli("report", out)[1])

    def test_exclude_on_the_command_line(self):
        out = os.path.join(self.tmp, "o")
        rc, txt = self.setup(out, "--exclude", "T3,T4", "--exclude", "T5")
        self.assertEqual(rc, 0, txt)
        self.assertIn("excluded from the alternatives: 3 reference branch(es)", txt)
        with open(os.path.join(out, "setup.json")) as fh:
            m = json.load(fh)
        self.assertEqual(m["alternatives_written"], 20)
        self.assertEqual(m["exclude"], ["T3 (1 reference branch)", "T4 (1 reference branch)",
                                        "T5 (1 reference branch)"])
        rc, txt = self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--exclude",
                           "T3,T4,T5")
        self.assertIn("alternatives: 20", txt)
        self.assertIn("resolution, all alternatives: 20, so the smallest rank-based P is 1/21", txt)

    def test_report_gives_length_correlation(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("length vs LRT among the alternatives: not defined", self.cli("report", out)[1])
        with open(self.tree, "w") as fh:   # unequal tip lengths: now it is defined
            fh.write(_caterpillar(self.names, "T12").replace("T3:0.1", "T3:0.3").replace("T7:0.1", "T7:0.02"))
        out = os.path.join(self.tmp, "o2")
        self.setup(out)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("length vs LRT among the alternatives: Spearman rho", self.cli("report", out)[1])

    def test_codeml_icode_and_add_seeds_keep_the_code(self):
        tree = os.path.join(self.tmp, "m.nwk")
        with open(tree, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace("{FG}", " #1"))
        out = os.path.join(self.tmp, "cm")
        rc, txt = self.cli("setup", "--engine", "codeml", "--alignment", self.aln, "--tree", tree, "--test", "#1",
                           "--out", out, "--codeml", os.path.join(self.tmp, "no-codeml"), "--seeds", "1",
                           "--code", "2")
        self.assertEqual(rc, 0, txt)
        self.assertIn("not a release this tool was checked against", txt)
        self.cli("add-seeds", out, "--seeds", "2", "--which", "all")
        for run in ("observed_s1", "observed_s2", "alt001_s2"):
            for sub in ("alt", "null"):
                with open(os.path.join(out, "runs", run, sub, "codeml.ctl")) as fh:
                    self.assertIn("icode = 1\n", fh.read())
        out = os.path.join(self.tmp, "o")
        self.setup(out, "--code", "2")
        self.cli("add-seeds", out, "--seeds", "2")
        with open(os.path.join(out, "runs", "alt001_s2", "cmd.sh")) as fh:
            self.assertIn("--code 'Vertebrate-mtDNA'", fh.read())

    def test_manifest_exclude_column(self):
        man = os.path.join(self.tmp, "genes.tsv")
        with open(man, "w") as fh:
            fh.write("gene\talignment\ttree\texclude\ng1\ta.fa\tt.nwk\tT3,T4\ng2\ta.fa\tt.nwk\t\n")
        out = os.path.join(self.tmp, "b")
        rc, txt = self.cli("batch-setup", "--manifest", man, "--test", "FG", "--out", out, "--hyphy", self.hy,
                           "--relax-bf", self.bf, "--seeds", "1")
        self.assertEqual(rc, 0, txt)
        got = {}
        for g in ("g1", "g2"):
            with open(os.path.join(out, g, "setup.json")) as fh:
                got[g] = json.load(fh)["alternatives_written"]
        self.assertEqual(got, {"g1": 21, "g2": 23})

    def test_untested_version_is_flagged(self):
        os.environ["FAKE_VERSION"] = "2.5.50"
        out = os.path.join(self.tmp, "o")
        rc, txt = self.setup(out)
        self.assertIn("is not a release this tool was checked against (2.5.101)", txt)
        self.cli("run", out, "--jobs", "4")
        self.assertIn("WARNING: HYPHY 2.5.50(MP)", self.cli("report", out)[1])

    def test_slurm_array_task_runs_as_slurm_would(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        self.cli("slurm", out)
        sb = os.path.join(out, "slurm", "array.sbatch")
        env = dict(os.environ, SLURM_ARRAY_TASK_ID="2")
        self.assertEqual(subprocess.call(["bash", sb], env=env, stdout=subprocess.DEVNULL), 0)
        with open(os.path.join(out, "slurm", "runs.txt")) as fh:
            d = fh.read().splitlines()[1]
        with open(os.path.join(d, "DONE")) as fh:
            self.assertTrue(fh.read().startswith("ok\t0\t"))
        env["SLURM_ARRAY_TASK_ID"] = "999"
        self.assertNotEqual(subprocess.call(["bash", sb], env=env, stdout=subprocess.DEVNULL), 0)
        for i in range(1, 25):   # every task, then the report reads what the array left
            env["SLURM_ARRAY_TASK_ID"] = str(i)
            subprocess.call(["bash", sb], env=env, stdout=subprocess.DEVNULL)
        txt = self.cli("report", out)[1]
        self.assertIn("(a) all alternatives: 23 read (of 23 written); LRT >= observed: 0; rank 1 of 24", txt)


class Refusals(_Fixture):
    """Inputs that used to crash, leave a half-written folder, or give a meaningless reading."""

    def _setup_err(self, out, *extra):
        from contextlib import redirect_stderr
        err = io.StringIO()
        with redirect_stderr(err):
            rc, txt = self.setup(out, *extra)
        return rc, txt + err.getvalue()

    def test_bad_options_refused_before_anything_is_written(self):
        for extra, msg in ((("--length-ratio", "0"), "--length-ratio"), (("--length-ratio", "0.5"), "--length-ratio"),
                           (("--length-ratio", "nan"), "--length-ratio"), (("--max-alternatives", "-3"), "--max-alt"),
                           (("--alpha", "1.5"), "--alpha"), (("--seeds=-1,1",), "--seeds"),
                           (("--sequential", "0"), "--sequential")):
            out = os.path.join(self.tmp, "x")
            rc, txt = self._setup_err(out, *extra)
            self.assertEqual(rc, 2, (extra, txt))
            self.assertIn(msg, txt)
            self.assertFalse(os.path.exists(out), extra)
        rc, txt = self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--length-ratio", "0")
        self.assertEqual(rc, 2, txt)

    def test_codeml_mark_and_seeds(self):
        tree = os.path.join(self.tmp, "m.nwk")
        with open(tree, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace("{FG}", " #1"))
        for test, seeds in (("#0", "1"), ("#2", "1"), ("#1", "0,1")):
            out = os.path.join(self.tmp, "cm")
            rc, txt = self.cli("setup", "--engine", "codeml", "--alignment", self.aln, "--tree", tree, "--test", test,
                               "--out", out, "--codeml", os.path.join(self.tmp, "no-codeml"), "--seeds", seeds)
            self.assertEqual(rc, 2, (test, seeds))
            self.assertFalse(os.path.exists(out))

    def test_codeml_two_ratio_control_files_and_add_seeds(self):
        tree = os.path.join(self.tmp, "m.nwk")
        with open(tree, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace("{FG}", " #1"))
        out = os.path.join(self.tmp, "tr")
        rc, txt = self.cli("setup", "--engine", "codeml", "--model", "two-ratio", "--alignment", self.aln, "--tree",
                           tree, "--test", "#1", "--out", out, "--codeml", os.path.join(self.tmp, "no-codeml"),
                           "--seeds", "1")
        self.assertEqual(rc, 0, txt)
        self.assertIn("engine codeml two-ratio", txt)
        with open(os.path.join(out, "setup.json")) as fh:
            self.assertEqual(json.load(fh)["codeml_model"], "two-ratio")
        self.cli("add-seeds", out, "--seeds", "2", "--which", "all")

        def ctl(run, sub):
            with open(os.path.join(out, "runs", run, sub, "codeml.ctl")) as fh:
                return dict(l.split(" = ") for l in (x.strip() for x in fh) if " = " in l)
        for run, om, ka in (("observed_s1", "0.5", "2.0"), ("observed_s2", "1.5", "1.0"), ("alt001_s2", "1.5", "1.0")):
            a, n = ctl(run, "alt"), ctl(run, "null")
            self.assertEqual((a["model"], a["NSsites"], a["fix_omega"], a["omega"], a["kappa"]),
                             ("2", "0", "0", om, ka), run)
            # the one-ratio null has a free omega, so it starts from the seed's omega too
            self.assertEqual((n["model"], n["NSsites"], n["fix_omega"], n["omega"], n["kappa"]),
                             ("0", "0", "0", om, ka), run)
            self.assertEqual((a["CodonFreq"], a["cleandata"]), ("2", "0"))
        # branch-site stays as it was: null fixes omega2 at 1, and add-seeds leaves that 1 alone
        out = os.path.join(self.tmp, "bs")
        rc, txt = self.cli("setup", "--engine", "codeml", "--alignment", self.aln, "--tree", tree, "--test", "#1",
                           "--out", out, "--codeml", os.path.join(self.tmp, "no-codeml"), "--seeds", "1")
        self.assertEqual(rc, 0, txt)
        self.cli("add-seeds", out, "--seeds", "3", "--which", "all")
        n = ctl("observed_s3", "null")
        self.assertEqual((n["model"], n["NSsites"], n["fix_omega"], n["omega"]), ("2", "2", "1", "1.0"))
        with open(os.path.join(out, "setup.json")) as fh:
            self.assertEqual(json.load(fh)["codeml_model"], "branch-site")

    def test_model_is_refused_for_relax(self):
        out = os.path.join(self.tmp, "o")
        rc, txt = self.cli("setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", out,
                           "--hyphy", self.hy, "--model", "two-ratio")
        self.assertEqual(rc, 2, txt)   # the message goes to stderr: "--model chooses the codeml model ..."
        self.assertFalse(os.path.exists(out))

    def test_bare_program_name_beside_a_folder_of_that_name(self):
        bindir = os.path.join(self.tmp, "bin")
        os.makedirs(bindir)
        shutil.copy(self.hy, os.path.join(bindir, "hyphy"))
        os.makedirs(os.path.join(self.tmp, "w", "hyphy"))
        cwd, path = os.getcwd(), os.environ["PATH"]
        os.chdir(os.path.join(self.tmp, "w"))
        os.environ["PATH"] = bindir + os.pathsep + path
        try:
            rc, txt = self.cli("setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", "o",
                               "--relax-bf", self.bf, "--seeds", "1")
        finally:
            os.chdir(cwd)
            os.environ["PATH"] = path
        self.assertEqual(rc, 0, txt)
        with open(os.path.join(self.tmp, "w", "o", "runs", "observed_s1", "cmd.sh")) as fh:
            self.assertIn("&& 'hyphy' CPU=1", fh.read())

    def test_empty_exclude_name(self):
        from contextlib import redirect_stderr
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--exclude", "")

    def test_relative_program_paths_become_absolute(self):
        cwd = os.getcwd()
        os.chdir(self.tmp)
        try:
            out = os.path.join(self.tmp, "o")
            rc, txt = self.cli("setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", out,
                               "--hyphy", "./hyphy", "--relax-bf", os.path.join("bf", "RELAX.bf"), "--seeds", "1")
            self.assertEqual(rc, 0, txt)
        finally:
            os.chdir(cwd)
        with open(os.path.join(out, "runs", "observed_s1", "cmd.sh")) as fh:
            self.assertIn(self.hy, fh.read())
        self.assertTrue(os.path.isdir(os.path.join(out, "runs", "observed_s1", "modules")))

    def test_inputs_without_content(self):
        empty = os.path.join(self.tmp, "e.fa")
        with open(empty, "w") as fh:
            fh.write("".join(">%s\n---\n" % n for n in self.names))
        rc, txt = self.cli("check", "--alignment", empty, "--tree", self.tree, "--test", "FG")
        self.assertEqual(rc, 2)
        self.assertIn("no sequence has a determined base", txt)
        nolen = os.path.join(self.tmp, "n.nwk")
        with open(nolen, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace(":0.1", "").replace(":0.05", ""))
        rc, txt = self.cli("check", "--alignment", self.aln, "--tree", nolen, "--test", "FG")
        self.assertIn("the tree has no branch lengths", txt)
        two = os.path.join(self.tmp, "two.nwk")
        with open(two, "w") as fh:
            fh.write(_caterpillar(self.names, "T12").replace("T3:0.1", "T3:0.1:0.5"))
        rc, txt = self.cli("check", "--alignment", self.aln, "--tree", two, "--test", "FG")
        self.assertEqual(rc, 2)
        self.assertIn("two branch lengths", txt)

    def test_all_alternatives_failed_gives_no_rank(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        rd = os.path.join(out, "runs")
        for x in os.listdir(rd):
            if x.startswith("alt"):
                with open(os.path.join(rd, x, "cmd.sh"), "w") as fh:
                    fh.write("exit 1\n")
        self.cli("run", out, "--jobs", "4")
        txt = self.cli("report", out)[1]
        self.assertIn("NO RANK: none of the 23 alternatives written could be read", txt)
        from labelrank import paper, report
        self.assertIn("so it was not ranked", paper.results(report.summarize(out)))

    def test_folders_and_scheduler_values(self):
        from contextlib import redirect_stderr
        for cmd in (("report", os.path.join(self.tmp, "nope")), ("run", os.path.join(self.tmp, "nope")),
                    ("batch-report", self.tmp)):
            with redirect_stderr(io.StringIO()):
                self.assertEqual(self.cli(*cmd)[0], 2, cmd)
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        with redirect_stderr(io.StringIO()):
            self.assertEqual(self.cli("slurm", out, "--max-parallel", "0")[0], 2)
            self.assertEqual(self.cli("slurm", out, "--time", "1:00:00\n#SBATCH --x=y")[0], 2)
        self.cli("run", out, "--jobs", "4")
        fig = os.path.join(self.tmp, "new dir", "f.svg")
        self.assertEqual(self.cli("report", out, "--figure", fig)[0], 0)
        self.assertTrue(os.path.exists(fig))
        f = os.path.join(self.tmp, "afile")
        open(f, "w").close()
        self.assertEqual(self._setup_err(f)[0], 2)
        with redirect_stderr(io.StringIO()):
            for cmd in (("run", out, "--timeout-hours", "-1"), ("run", out, "--timeout-hours", "nan"),
                        ("report", out, "--alpha", "7"), ("report", out, "--alpha", "nan"),
                        ("add-seeds", out, "--seeds=-5")):
                self.assertEqual(self.cli(*cmd)[0], 2, cmd)
        self.assertFalse(os.path.exists(os.path.join(out, "runs", "observed_s-5")))

    def test_scheduler_paths_that_a_shell_or_slurm_would_misread(self):
        from contextlib import redirect_stderr
        out = os.path.join(self.tmp, "p$HOME q", "o")
        os.makedirs(os.path.dirname(out))
        self.setup(out)
        with redirect_stderr(io.StringIO()) as err:
            self.assertEqual(self.cli("slurm", out)[0], 2)
        self.assertIn("'$'", err.getvalue())
        out = os.path.join(self.tmp, "it's o")   # a quote is fine: every path is quoted for the shell
        self.setup(out)
        self.cli("slurm", out)
        env = dict(os.environ, SLURM_ARRAY_TASK_ID="1")
        sb = os.path.join(out, "slurm", "array.sbatch")
        self.assertEqual(subprocess.call(["bash", sb], env=env, stdout=subprocess.DEVNULL), 0)
        self.assertTrue(os.path.exists(os.path.join(out, "runs", "alt001_s1", "DONE")))

    def test_sequential_before_any_alternative_has_run(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out, "--sequential", "3")
        from labelrank import runner, paper, report
        self.assertEqual(runner.run_one(os.path.join(out, "runs", "observed_s1"), None)[0], "ok")
        txt = self.cli("report", out)[1]
        self.assertIn("NO RANK YET", txt)
        self.assertNotIn("could be read", txt)
        self.assertIn("had finished yet", paper.results(report.summarize(out)))

    def test_manifest_problems(self):
        from labelrank.batch import read_manifest
        from labelrank.setup import SetupError
        m = os.path.join(self.tmp, "m.tsv")
        for body, msg in (("", "no gene lines"), ("gene\talignment\ttree\n", "no gene lines"),
                          ("gene\talignment\ttree\n\n\ng1\ta\tt\nG1\ta\tt\n", "line 5"),
                          ("gene\talignment\ttree\nbatch.json\ta\tt\n", "the batch folder uses")):
            with open(m, "w") as fh:
                fh.write(body)
            with self.assertRaises(SetupError) as cm:
                read_manifest(m)
            self.assertIn(msg, str(cm.exception))


def _shq(x):
    return "'" + x.replace("'", "'\\''") + "'"


def _contend(p, bar, bar2, q):
    import socket
    from labelrank import runner
    bar.wait()
    q.put(runner.take(p, "%s %d" % (socket.gethostname(), os.getpid())))
    bar2.wait()   # holders stay alive until every contender has tried


class Runner(_Fixture):
    """Interrupts, timeouts and a second runner on the same folder, with a slow stand-in engine."""

    def _run(self, out, *extra, **env):
        e = dict(os.environ, PYTHONPATH=os.path.join(HERE, ".."), PYTHONDONTWRITEBYTECODE="1", **env)
        return subprocess.Popen([sys.executable, "-m", "labelrank", "run", out] + list(extra), env=e,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def _done(self, out):
        rd = os.path.join(out, "runs")
        res = {}
        for x in os.listdir(rd):
            if os.path.exists(os.path.join(rd, x, "DONE")):
                with open(os.path.join(rd, x, "DONE")) as fh:
                    res[x] = fh.read().split("\t")[0]
        return res

    def test_interrupt_leaves_no_done_and_resumes(self):
        import signal
        import time
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        p = self._run(out, "--jobs", "3", FAKE_SLEEP="3")
        time.sleep(1.5)
        p.send_signal(signal.SIGINT)
        txt = p.communicate(timeout=30)[0].decode()
        self.assertEqual(p.returncode, 130, txt)
        self.assertIn("interrupted", txt)
        self.assertEqual(self._done(out), {})                 # nothing marked FAIL for good
        self.assertFalse(os.path.exists(os.path.join(out, "RUNNING")))
        time.sleep(4)                                          # a surviving engine would write by now
        self.assertFalse([x for x in os.listdir(os.path.join(out, "runs"))
                          if os.path.exists(os.path.join(out, "runs", x, "relax.json"))])
        p = self._run(out, "--jobs", "4")
        p.communicate(timeout=120)
        self.assertEqual(p.returncode, 0)
        self.assertEqual(set(self._done(out).values()), {"ok"})

    def test_timeout_stops_the_engine(self):
        import time
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        p = self._run(out, "--jobs", "2", "--timeout-hours", "0.0003", FAKE_SLEEP="4")
        time.sleep(4)
        p.kill()   # only the first runs matter
        p.communicate()
        done = self._done(out)
        self.assertTrue(done and set(done.values()) == {"TIMEOUT"}, done)
        time.sleep(3)
        for x in done:
            self.assertFalse(os.path.exists(os.path.join(out, "runs", x, "relax.json")), x)

    def test_engine_that_ignores_sigterm_does_not_outlive_a_timeout(self):
        import time
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        d = os.path.join(out, "runs", "observed_s1")
        late = os.path.join(d, "late.txt")
        with open(os.path.join(d, "cmd.sh"), "w") as fh:   # the shell dies on TERM, its child ignores it
            fh.write("bash -c 'trap \"\" TERM; sleep 3; echo x > \"$1\"' _ %s &\nwait\n" % _shq(late))
        from labelrank import runner
        self.assertEqual(runner.run_one(d, 0.5)[0], "TIMEOUT")
        time.sleep(4)
        self.assertFalse(os.path.exists(late))

    def test_a_run_finished_by_another_process_while_waiting_for_its_lock_is_not_redone(self):
        from labelrank import runner
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        d = os.path.join(out, "runs", "observed_s1")
        take = runner.take

        def take_after_other_finished(p, content):
            with open(os.path.join(d, "DONE"), "w") as fh:
                fh.write("ok\t0\t1s\tother\n")
            return take(p, content)
        runner.take = take_after_other_finished
        try:
            self.assertEqual(runner.run_one(d, None)[0], "skip")
        finally:
            runner.take = take
        with open(os.path.join(d, "DONE")) as fh:
            self.assertIn("other", fh.read())
        self.assertFalse(os.path.exists(os.path.join(d, "RUN_LOCK")))

    def test_stale_lock_taken_over_by_exactly_one_of_many(self):
        import multiprocessing as mp
        import socket
        p = os.path.join(self.tmp, "RUN_LOCK")
        for _ in range(20):
            with open(p, "w") as fh:
                fh.write("%s 999999\n" % socket.gethostname())   # left by a dead process
            ctx = mp.get_context("fork")
            bar, bar2, q = ctx.Barrier(8), ctx.Barrier(8), ctx.Queue()
            ps = [ctx.Process(target=_contend, args=(p, bar, bar2, q)) for _ in range(8)]
            for x in ps:
                x.start()
            for x in ps:
                x.join()
            self.assertEqual(sum(q.get() for _ in range(8)), 1)
            os.remove(p)

    def test_second_runner_is_refused(self):
        import socket
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        with open(os.path.join(out, "RUNNING"), "w") as fh:
            fh.write("%s %d\n" % (socket.gethostname(), os.getpid()))
        rc, txt = self.cli("run", out)
        self.assertEqual(rc, 2)
        self.assertEqual(self._done(out), {})
        with open(os.path.join(out, "RUNNING"), "w") as fh:   # a stale lock (dead process) is taken over
            fh.write("%s 999999\n" % socket.gethostname())
        self.assertEqual(self.cli("run", out, "--jobs", "4")[0], 0)

    def test_two_runners_started_together_run_each_folder_once(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        count = os.path.join(self.tmp, "count.log")
        ps = [self._run(out, "--jobs", "4", FAKE_COUNT=count, FAKE_SLEEP="1") for _ in range(2)]
        for p in ps:
            p.communicate(timeout=120)
        self.assertEqual(sorted(p.returncode for p in ps), [0, 2])
        with open(count) as fh:
            self.assertEqual(len(fh.read().splitlines()), 24)
        self.assertFalse(os.path.exists(os.path.join(out, "RUNNING")))

    def test_a_scheduler_task_in_progress_is_left_alone(self):
        import time
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        self.cli("slurm", out)
        count = os.path.join(self.tmp, "count.log")
        d = os.path.join(out, "runs", "observed_s1")
        env = dict(os.environ, FAKE_SLEEP="3", FAKE_COUNT=count)
        task = subprocess.Popen(["bash", os.path.join(out, "slurm", "run_one.sh"), d], env=env,
                                stdout=subprocess.DEVNULL)
        time.sleep(1)
        p = self._run(out, "--jobs", "4", FAKE_COUNT=count)
        txt = p.communicate(timeout=120)[0].decode()
        self.assertIn("1 run(s) were left alone because another process holds their RUN_LOCK", txt)
        self.assertEqual(p.returncode, 1)                      # not finished: no ALL_DONE, non-zero exit
        self.assertFalse(os.path.exists(os.path.join(out, "ALL_DONE")))
        task.wait(timeout=30)
        with open(count) as fh:
            starts = fh.read().splitlines()
        self.assertEqual(len(starts), 24)                      # observed_s1 started once, by the task
        self.assertEqual(self._done(out)["observed_s1"], "ok")
        self.assertFalse(os.path.exists(os.path.join(d, "RUN_LOCK")))

    def test_stale_output_cannot_pass_as_a_result(self):
        out = os.path.join(self.tmp, "o")
        self.setup(out)
        d = os.path.join(out, "runs", "observed_s1")
        with open(os.path.join(d, "relax.json"), "w") as fh:
            fh.write("{}")
        with open(os.path.join(d, "cmd.sh"), "w") as fh:
            fh.write("exit 0\n")   # an engine that writes nothing and exits 0
        self.cli("run", out, "--jobs", "4")
        self.assertEqual(self._done(out)["observed_s1"], "FAIL")


if __name__ == "__main__":
    unittest.main()
