"""setup -> run -> report with a stand-in for HyPhy (tests/fake_hyphy.py): the pipeline, not any statistic."""
import io
import os
import shutil
import stat
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from labelrank import cli  # noqa: E402


def _rows(out):
    with open(os.path.join(out, "runs.tsv")) as fh:
        return [l.split("\t") for l in fh.read().splitlines()[1:]]


def shlex_last(cmd):
    import shlex
    return shlex.split(cmd)[-1]


def _tips(n):
    return ["T%d" % i for i in range(n)]


def _caterpillar(names, fg):
    t = names[0] + ":0.1"
    for x in names[1:]:
        t = "(%s,%s%s:0.1):0.05" % (t, x, "{FG}" if x == fg else "")
    return t + ";"


def _fasta(names, path, stop=False):
    with open(path, "w") as fh:
        for i, n in enumerate(names):
            fh.write(">%s\nATGAAA%s\n" % (n, "TAA" if stop and i == 0 else "CCC"))


class EndToEnd(unittest.TestCase):
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

    def tearDown(self):
        shutil.rmtree(self.tmp)
        os.environ.pop("FAKE_HIGH", None)

    def cli(self, *args):
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cli.main(list(args))
        return rc, buf.getvalue()

    def setup_args(self, out, *extra):
        return ["setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", out,
                "--hyphy", self.hy, "--relax-bf", self.bf] + list(extra)

    def test_full_pipeline_with_paper_outputs(self):
        os.environ["FAKE_HIGH"] = "T12"
        out = os.path.join(self.tmp, "o")
        self.assertEqual(self.cli(*self.setup_args(out))[0], 0)
        self.assertEqual(self.cli("run", out, "--jobs", "4")[0], 0)
        m, f = os.path.join(self.tmp, "m.txt"), os.path.join(self.tmp, "f.svg")
        rc, txt = self.cli("report", out, "--methods", m, "--figure", f, "--table", os.path.join(self.tmp, "t.tsv"))
        self.assertNotIn("excluded", txt)
        # 25 tips: 24 alternatives (the root's two branches are never used; one of them is a tip here)
        self.assertIn("-> PASSED", txt)
        self.assertIn("LRT >= observed: 0; rank 1 of 24; rank-based P = 1/24", txt)
        with open(m) as fh:
            self.assertIn("labelrank v", fh.read())
        with open(f) as fh:
            self.assertTrue(fh.read().startswith("<svg"))

    def test_sequential_stops_early(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names)   # every label is 'high': stop after h alternatives
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--sequential", "3"))
        self.cli("run", out, "--jobs", "1")
        done = [d for d in os.listdir(os.path.join(out, "runs")) if os.path.exists(os.path.join(out, "runs", d, "DONE"))]
        self.assertEqual(len(done), 3 + 3 * 3)   # 3 observed seeds, then 3 alternatives with 3 seeds each
        txt = self.cli("report", out)[1]
        self.assertIn("stopped after 3 alternatives in the drawn order, 3 at or above", txt)
        self.assertIn("sequential P = 1.0000", txt)

    def test_sequential_with_parallel_runs_ignores_runs_past_the_stop(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names)
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--sequential", "3"))
        self.cli("run", out, "--jobs", "4")
        txt = self.cli("report", out)[1]
        self.assertIn("stopped after 3 alternatives in the drawn order", txt)

    def _rank_line(self, txt):
        return [l for l in txt.splitlines() if "(a) all alternatives" in l][0]

    def test_sequential_rank_does_not_depend_on_jobs(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names[::2])
        lines = []
        for jobs in ("1", "6"):
            out = os.path.join(self.tmp, "o" + jobs)
            self.cli(*self.setup_args(out, "--sequential", "2"))
            self.cli("run", out, "--jobs", jobs)
            lines.append(self._rank_line(self.cli("report", out)[1]))
        self.assertEqual(lines[0], lines[1])

    def test_unfinished_run_in_the_order_blocks_the_sequential_reading(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names)
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--sequential", "3"))
        self.cli("run", out, "--jobs", "1")
        rows = _rows(out)
        first = [r[0] for r in rows if r[8] == "1"][0]
        os.remove(os.path.join(out, "runs", first, "DONE"))   # as if still running
        txt = self.cli("report", out)[1]
        self.assertIn("sequential reading (h = 3): incomplete", txt)
        self.assertIn("not run (no DONE marker): ", txt)

    def test_never_launched_runs_are_not_called_failures(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names)
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--sequential", "3"))
        self.cli("run", out, "--jobs", "1")
        m = os.path.join(self.tmp, "m.txt")
        txt = self.cli("report", out, "--methods", m)[1]
        self.assertNotIn("excluded", txt)
        with open(m) as fh:
            self.assertNotIn("failed the instrument check", fh.read())

    def test_alternative_with_one_failed_seed_is_still_a_draw(self):
        os.environ["FAKE_HIGH"] = ",".join(self.names)
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--sequential", "2", "--alt-seeds", "1,2"))
        self.cli("run", out, "--jobs", "1")
        rows = _rows(out)
        first = [r[0] for r in rows if r[8] == "1" and r[3] == "2"][0]
        with open(os.path.join(out, "runs", first, "DONE"), "w") as fh:
            fh.write("FAIL\t1\t1s\tnow\n")
        os.remove(os.path.join(out, "runs", first, "relax.json"))
        txt = self.cli("report", out)[1]
        self.assertIn("stopped after 2 alternatives in the drawn order, 2 at or above", txt)

    def test_commands_txt_runs_with_spaces_in_paths(self):
        os.environ["FAKE_HIGH"] = "T12"
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out))
        self.cli("slurm", out)
        import subprocess
        with open(os.path.join(out, "slurm", "commands.txt")) as fh:
            first = fh.readline()
        self.assertEqual(subprocess.call(first, shell=True), 0)
        d = shlex_last(first)
        with open(os.path.join(d, "DONE")) as fh:
            self.assertTrue(fh.read().startswith("ok\t0\t"))

    def test_slurm_with_nothing_left_writes_nothing(self):
        os.environ["FAKE_HIGH"] = "T12"
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out))
        self.cli("run", out, "--jobs", "4")
        txt = self.cli("slurm", out)[1]
        self.assertIn("nothing to run", txt)
        self.assertFalse(os.path.exists(os.path.join(out, "slurm")))

    def test_batch_skips_a_bad_gene_and_reports_the_rest(self):
        os.environ["FAKE_HIGH"] = "T12"
        bad = os.path.join(self.tmp, "bad.fa")
        _fasta(self.names, bad, stop=True)
        man = os.path.join(self.tmp, "genes.tsv")
        with open(man, "w") as fh:
            fh.write("gene\talignment\ttree\ng1\ta.fa\tt.nwk\ng2\tbad.fa\tt.nwk\n")
        out = os.path.join(self.tmp, "b")
        rc, txt = self.cli("batch-setup", "--manifest", man, "--test", "FG", "--out", out, "--hyphy", self.hy,
                           "--relax-bf", self.bf)
        self.assertIn("1 gene(s) set up, 1 skipped", txt)
        self.cli("run", out, "--jobs", "4")
        rc, txt = self.cli("batch-report", out)
        self.assertIn("genes with a complete reading: 1; observed label passed the gate in 1 of them", txt)
        self.assertIn("ranked first among same-shaped labels in 1 of the 1; expected if the label were exchangeable "
                      "with its alternatives: 0.04", txt)
        with open(os.path.join(out, "genes.tsv")) as fh:
            rows = fh.read()
        self.assertIn("g2\tskipped at setup", rows)
        self.assertEqual(set(len(l.split("\t")) for l in rows.splitlines()), {15})

    def test_slurm_lists_unfinished_runs(self):
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out))
        self.cli("slurm", out, "--partition", "short")
        sd = os.path.join(out, "slurm")
        with open(os.path.join(sd, "runs.txt")) as fh:
            self.assertEqual(len(fh.read().splitlines()), 3 + 23 * 3)
        with open(os.path.join(sd, "array.sbatch")) as fh:
            sb = fh.read()
        self.assertIn("#SBATCH --array=1-72%100", sb)
        self.assertIn("#SBATCH --partition=short", sb)

    def test_add_seeds_makes_alternatives_match(self):
        os.environ["FAKE_HIGH"] = "T12"
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--alt-seeds", "1"))
        self.cli("run", out, "--jobs", "4")
        self.assertIn("warning: the observed label used seeds", self.cli("report", out)[1])
        self.cli("add-seeds", out, "--seeds", "2,3")
        self.assertIn("INCOMPLETE", self.cli("report", out)[1])
        self.cli("run", out, "--jobs", "4")
        txt = self.cli("report", out)[1]
        self.assertNotIn("warning: the observed label used seeds", txt)
        self.assertNotIn("INCOMPLETE", txt)
        with open(os.path.join(out, "runs", "alt001_s3", "RELAX_seeded.bf")) as fh:
            self.assertEqual(fh.readline(), "SetParameter (RANDOM_SEED, 3, 0);\n")

    def test_user_mistakes_give_clean_errors(self):
        from contextlib import redirect_stderr
        out = os.path.join(self.tmp, "o")
        rc, txt = self.cli("check", "--alignment", self.aln, "--tree", self.tree, "--test", "fg")
        self.assertEqual(rc, 2)
        self.assertIn("labels in the tree: FG", txt)
        err = io.StringIO()
        with redirect_stderr(err):
            rc = self.cli("setup", "--alignment", self.aln, "--tree", self.tree, "--test", "FG", "--out", out,
                          "--hyphy", "no_such_hyphy_here")[0]
        self.assertEqual(rc, 2)
        self.assertIn("not found on PATH", err.getvalue())
        self.cli(*self.setup_args(out))
        err = io.StringIO()
        with redirect_stderr(err):
            rc = self.cli("report", out)[0]
        self.assertEqual(rc, 2)
        self.assertIn("nothing has run yet", err.getvalue())

    def test_add_seeds_refuses_without_leaving_folders(self):
        from contextlib import redirect_stderr
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out, "--alt-seeds", "1"))
        before = sorted(os.listdir(os.path.join(out, "runs")))
        with redirect_stderr(io.StringIO()):
            self.assertEqual(self.cli("add-seeds", out, "--seeds", "3,3")[0], 2)
        self.assertEqual(sorted(os.listdir(os.path.join(out, "runs"))), before)
        # a run with output but no DONE looks in progress: refused
        d = os.path.join(out, "runs", "alt001_s1")
        open(os.path.join(d, "relax.out"), "w").close()
        with redirect_stderr(io.StringIO()):
            self.assertEqual(self.cli("add-seeds", out, "--seeds", "2")[0], 2)
        self.assertEqual(sorted(os.listdir(os.path.join(out, "runs"))), before)

    def test_existing_output_is_refused(self):
        out = os.path.join(self.tmp, "o")
        self.cli(*self.setup_args(out))
        buf = io.StringIO()
        from contextlib import redirect_stderr
        with redirect_stderr(buf):
            rc = self.cli(*self.setup_args(out))[0]
        self.assertEqual(rc, 2)
        self.assertIn("refusing to write over", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
