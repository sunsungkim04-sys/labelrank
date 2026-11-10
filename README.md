# labelrank

Rank a trait-labelled branch-selection test among every same-shaped label on the same tree.

A branch-specific selection test (HyPhy RELAX, or PAML codeml branch-site model A or the two-ratio branch model)
run on a trait label is
usually read against its chi-square reference: P <= 0.05 and the trait lineage is said to be under relaxed,
intensified or positive selection. On real trees, other lineages that share nothing with the trait often give
"significant" results too. A chi-square P alone cannot separate a trait-linked signal from that background.

`labelrank` gives the missing comparison. It takes your alignment and your labelled tree, builds every other label
of the same shape on the same tree (or a uniform sample of them), runs each with the same command, and reports
where your label ranks among them.

## What "same shape" means

The test set is split into connected groups of branches. Each group must be one of:

| shape | test set | alternatives |
|---|---|---|
| tip | one terminal branch | every other terminal branch in the reference set |
| clade | an internal branch and every branch below it (e.g. a cherry: 3 branches) | every clade with the same rooted topology, all of whose branches are in the reference set |
| stem | one internal branch only | every internal branch above a clade of the same topology |

Several separate groups (for example three independent lineages) are matched group by group; the alternatives
are combinations of non-overlapping groups of the same shapes, sampled at random when there are too many.

In each alternative the alternative's branches take the test label. The original test branches move, for RELAX,
to a label that is neither test nor reference (they join the unclassified, nuisance branches, and the command line
is unchanged); for codeml, which has only foreground and background, they return to the background. The two
branches directly below a bifurcating root are never used (HyPhy merges them), and a test branch there is refused:
reroot the tree.

Length matters: a short lineage carries few substitutions and gives a small LRT whatever its biology. Every
alternative's total test-set length is recorded, and the report gives two co-primary ranks, one over all
alternatives and one over the length-matched ones (total length within a factor of 2 of yours,
`--length-ratio`). `--length-matched-only` writes only those. It also gives Spearman's correlation between length
and LRT among the alternatives, so you can see whether length drives the LRT on your tree at all; with fewer than
3 alternatives at LRT 0.01 or more it is reported as not defined.

## Install

Python 3.8 or later, no dependencies.

```bash
pip install .            # or run in place: python3 -m labelrank ...
```

HyPhy (tested with 2.5.101 MP) or PAML codeml (tested with 4.10.10) must be installed separately. Give their
paths with `--hyphy` / `--codeml` unless they are on your PATH. Other releases run, but setup and the report warn:
the output layout this tool parses can change between releases, so check a few runs by hand
(`labelrank inspect` for RELAX) before trusting the report.

## Use

```bash
# 1. check inputs, list the alternatives and the smallest P they allow (writes nothing)
labelrank check --alignment aln.fasta --tree labelled.nwk --test Foreground

# 2. write one run folder per (label, seed)
labelrank setup --alignment aln.fasta --tree labelled.nwk --test Foreground --out lr_out \
    --hyphy /path/to/hyphy

# 3. run them (one CPU each); or submit each runs/*/cmd.sh to your scheduler
labelrank run lr_out --jobs 8

# 4. read; optionally write manuscript text and a figure
labelrank report lr_out --table lr_out/runs_read.tsv --methods lr_out/methods.txt --figure lr_out/rank.svg
```

`--methods` writes a methods paragraph (with tool version and citation), a results sentence and a figure caption,
filled in from this run. `--figure` writes a single-column (85 mm) SVG strip plot: every alternative's LRT (filled
dots length-matched, open dots not), the observed label marked, and the P = 0.05 line (chi-square(1) for RELAX, the 50:50 mixture for codeml). `--save` writes the
report text. Every output file is written only after the whole report has been computed.

### Many genes at once

```bash
# genes.tsv: header 'gene<TAB>alignment<TAB>tree', optional 'test', 'reference' and 'exclude' columns per
# gene (exclude: comma-separated, as --exclude); relative paths are read relative to the folder holding genes.tsv
labelrank batch-setup --manifest genes.tsv --test Foreground --out scan --hyphy /path/to/hyphy
labelrank run scan --jobs 32                 # one pool across all genes
labelrank batch-report scan                  # writes scan/genes.tsv
```

A gene whose inputs fail the checks is skipped and listed; the others go ahead (if none can be set up, the
command fails). The batch report gives one row
per gene (both ranks, the share of shorter alternatives, and status) and, over every gene with a complete
reading, how many rank first against the number expected if the label were exchangeable with its alternatives
(the sum of 1/(N+1)). The count is not restricted to genes where the observed label was significant: that filter
would select large LRTs and the expectation would no longer hold. Genes on one tree share the tested lineage, so
this is a description, not a test.

### Stopping early

`--sequential H` at setup draws a random running order for the alternatives. `labelrank run` then runs the
observed label first and the alternatives in that order, and launches no more once H of them reach the observed
LRT. The report gives the Besag-Clifford sequential P = H / L (L alternatives read when the H-th was reached), or
(g + 1) / (N + 1) when the order runs out first. Every reading (rank, secondaries, figure) then uses only the
alternatives walked up to the stopping point, so the answer does not depend on `--jobs`; runs that finished past
it are ignored, and a run without a `DONE` marker (still running) blocks the reading rather than being read
half-way. Labels that clearly do not stand out stop after a few runs;
labels that do stand out still run every alternative. H = 10 is a common choice.

### Clusters

```bash
labelrank slurm lr_out --partition long --time 48:00:00 --mem 2G
sbatch lr_out/slurm/array.sbatch
```

This writes a SLURM array over the unfinished run folders (split into arrays of at most 1000 tasks) and a plain
`commands.txt` (one `bash run_one.sh <folder>` per line) for GNU parallel or another scheduler. Each task leaves
the same `DONE` marker as `labelrank run`; a task the scheduler cancels or stops at its time limit gets no
`DONE`, so the next submission runs it again (raise `--time` if one keeps hitting the limit). The files hold absolute paths (rerun `labelrank slurm` after moving
the folder), and at most 100 tasks run at once (`--max-parallel`). An array cannot stop early; use `labelrank run`
for `--sequential`. The tests run the generated array script task by task the way SLURM does (with
`SLURM_ARRAY_TASK_ID` set) and read the result; it has not yet been submitted on a real cluster.

### codeml

Mark the foreground with `#1` and pass `--engine codeml` to every command that reads the tree:

```bash
labelrank check --engine codeml --alignment aln.fasta --tree marked.nwk --test '#1'
labelrank setup --engine codeml --alignment aln.fasta --tree marked.nwk --test '#1' --out lr_cm --codeml /path/to/codeml
```

`--model` picks the codeml model at setup (default `branch-site`); later commands read it from the setup:

| `--model` | alternative (`alt/`) | null (`null/`) | P | report's K column |
|---|---|---|---|---|
| `branch-site` | model A: model = 2, NSsites = 2 | the same with fix_omega = 1, omega = 1 | 50:50 mixture of 0 and chi-square(1) | foreground omega2 |
| `two-ratio` | model = 2, NSsites = 0: the `#1` branches get their own omega | model = 0, NSsites = 0: one omega | chi-square(1) | foreground omega |

```bash
labelrank setup --engine codeml --model two-ratio --alignment aln.fasta --tree marked.nwk --test '#1' --out lr_tr --codeml /path/to/codeml
```

Each fit runs in its own working folder (codeml writes fixed-name scratch files), with CodonFreq = 2 (F3x4),
cleandata = 0 and an estimated kappa. codeml has no random seed: a seed number picks one of six starting
(omega, kappa) pairs, (0.5, 2), (1.5, 1), (3, 4), (0.2, 0.5), (5, 2), (0.8, 8). The branch-site null fixes
omega2 at 1, so only kappa takes the start there; the one-ratio null has a free omega and takes both.

The branch-site null lies on the boundary of model A, so P uses the 50:50 mixture of 0 and chi-square(1); the
plain chi-square(1) is the conservative choice the PAML documentation also mentions. The two-ratio and one-ratio
models are nested at the interior point omega foreground = omega background, so P is chi-square(1). The two-ratio
test is two-sided: a foreground omega below the background one gives as large an LRT as one above it, and the
rank counts both. Read the K column for the direction. codeml 4.10.10 exits with status 1 after a complete run,
so a run counts as finished when both `mlc` files hold their likelihood line. codeml does not stop at a stop
codon: it turns that column into missing data and carries on, so the stop-codon check before setup matters
here too.

### Options fixed at setup

- **Seeds.** Every label, observed and alternative, is run under the same seeds (default 1, 2, 3: `--seeds`).
  Seeds are optimizer restarts of one analysis, not replicates. `--alt-seeds` can give the alternatives fewer;
  the report then warns, because labels given different restarts are not compared on equal terms. `labelrank add-seeds`
  adds seeds to an existing setup (not while `labelrank run` is working on it).
- **Gate level and direction.** `--alpha` (default 0.05) and, for RELAX, `--direction relax|intens` (K must be on
  that side of 1) are recorded in `setup.json`; `report` uses them unless told otherwise.
- **Genetic code.** `--code` takes an NCBI translation table number (default 1; 2 = vertebrate mitochondrial, 5 =
  invertebrate mitochondrial, and 3, 4, 6, 9-15). It sets the stop codons the input check looks for, HyPhy's
  `--code` and codeml's `icode` (codeml has no table 14). Table 11 runs as table 1: they assign codons the same way.
- **Lineages left out.** `--exclude` keeps lineages out of the alternatives, for example other lineages you expect
  to share the trait, which would otherwise count against your label. Give a tip or node name, or tips joined by
  `+` for the clade of their common ancestor (that branch and everything below it):
  `--exclude t3,t8+t9` (tip t3, and the clade of t8 and t9). Excluded branches stay in their set for the fits;
  they are only never given the test label. Name them before you see any result.
- **How many alternatives.** At most 200, sampled uniformly with `--sample-seed` when more exist
  (`--max-alternatives 0` = all). When the combinations of a multi-group test set are too many to list, they are
  drawn at random and the report says the number of valid sets is unknown.

The reference set: `--reference <label>`, or, when omitted, the unlabelled branches. If every branch carries a
label, give the reference label explicitly: HyPhy 2.5.101 refuses `--reference "Unlabeled branches"` on such a tree.

## What the report says

From `validation/scnn1_A_siren/` (published data, 12 alternatives):

```
observed LRT 7.2261 (best alternative vs best null fit over 3 seeds), P 0.007185, K 0.6428
gate: P <= 0.05 -> PASSED (seeds individually: 3 of 3 at P <= 0.05)
(a) all alternatives: 12 read (of 12 written); LRT >= observed: 0; rank 1 of 13
(b) length-matched (test-set length within x2 of 0.05071): 2 alternatives; LRT >= observed: 0; rank 1 of 3
length: 92% of the alternatives are shorter than the observed test set - the observed label is longer than most ...
length vs LRT among the alternatives: Spearman rho 0.15 (n = 12)
```

- **Label statistic.** For every label, the alternative and the null model are each taken at their best
  log-likelihood over the seeds, and LRT = 2 (logL alternative - logL null). An LRT below 0.01, negative ones
  included, is optimizer noise around 0 and is read as 0 everywhere (gate, ranks, correlation), so noise cannot
  order labels; a label whose alternative fit lies more than 0.01 below its null is flagged as not converged.
- **Ranks.** 1 + the number of alternatives whose LRT is at least yours, out of N + 1, over (a) all alternatives
  and (b) the length-matched ones. Ties count against your label, and an LRT within 0.01 of yours is a tie. Read them together: (a) has more alternatives,
  (b) removes most of the length effect. The report prints the share of alternatives shorter than your test set
  and flags it at 75% or more, with the length-LRT rank correlation among the alternatives.
- **Rank-based P** = rank / (N + 1), given for a reading only when its N + 1 >= 20 and every written alternative
  has finished. Below 20, the smallest possible value is already above 0.05, and only the rank is reported.
- **These are ranks among this tree's lineages, not exact null distributions.** The comparison assumes that,
  apart from the trait, your lineage is exchangeable with lineages of its shape on the same tree. They differ in
  history, are not independent of one another, and every alternative lies in your run's reference set while your
  lineage never does. A low rank says the trait lineage stands out among lineages of its shape on this tree; a
  middling rank says the test, by itself, does not separate it from them. Neither says the biology is right or
  wrong. Runs that fail are excluded, and failures may not be random with respect to LRT; the report lists them.
- **Gate.** The usual reading of your label: its P <= alpha (with `--direction`, also K on that side of 1). The
  ranks are computed either way.
- **Secondary readings:** directional (RELAX: alternatives at or above yours with K on the same side of 1) and how
  many alternatives the chi-square reading would call significant.
- **Sequential P** (with `--sequential`): see above.

## Checks built in, and why

Each of these comes from something that went wrong in real use.

1. **Read by log-likelihood, not by the reported verdict.** When RELAX meets a flat likelihood surface it refits K
   from the other side of 1 ("Attempt to fit an alternative direction of K"), and the output can hold two
   alternative fits. In one published output, the reported P and direction came from the lower-likelihood fit;
   the higher one pointed the other way. `labelrank` reads both fits and the null, uses
   LRT = 2 (max alternative logL - null logL), and warns when the reported LRT, P, K or direction disagree with it.
   `labelrank inspect relax.out [relax.json]` does this for any single RELAX output, including published ones
(RELAX only).
   Fits are paired by logL and K together, so near-tied fits on opposite sides of 1 are not mixed up. When RELAX
   refits after a negative LRT, the later passes print no section headers; the final values are read from
   `relax.json`, and a run with such a refit but no JSON is excluded.
2. **Several seeds, and the spread of logL among them.** The same command can give different optima; the report
   prints the spread of the alternative and of the null logL per label and flags it above `--logl-tolerance`
   (0.5 by default).
3. **Input checks before anything is written.** In-frame stop codons under the chosen genetic code (internal or
   terminal), a length that is not a multiple of 3, non-nucleotide characters, names HyPhy would rewrite, and
   alignment/tree name mismatches stop the tool. (In one case all 16 runs of a test died inside HyPhy on stop
   codons; codeml would instead have run on, with those columns silently treated as missing.)
4. **HyPhy JSON with bare `inf` / `nan`** is parsed, not dropped. A run whose output cannot be read is listed as
   excluded with the reason, never silently skipped.
5. **Zero-length branches.** HyPhy deletes zero-length internal branches after its first fit; they can be in the
   reference set as well as in the test set. A run that lost a test branch is excluded (its test set changed);
   deleted reference branches are counted and reported.
6. **Explicit reference names** on fully labelled trees (see above).
7. **Branch sets checked by name (RELAX).** Unnamed internal nodes are named (`LRn1`, `LRn2`, ...) in the tree each run
   receives, and the Test, Reference and nuisance sets HyPhy prints are compared, name by name, with the ones
   written (RELAX omits the nuisance line when that set is empty; that is read as empty). HyPhy merges the two
   branches below a bifurcating root into one and never prints one of them; that one branch is allowed to be
   missing unless it was a test branch. Alternatives never use either root branch. Any other difference excludes
   the run.

8. **Unfinished and failed runs.** A run without a `DONE` marker is never read, even if it has partial output;
   while any alternative is unfinished the report gives provisional ranks and no P. A run whose program failed
   is excluded even if part of its output parses. A non-finite log-likelihood excludes the run.
9. **Engine releases.** The release is read from the program's own banner at setup (and, for codeml, again from
   each `mlc`), never from its path; one this tool has not been checked against, or one that reports no release,
   is flagged in setup and in the report.

For RELAX the alignment is copied byte for byte into each run folder (codeml gets a PHYLIP copy with U, '.' and
'?' written as T, '-' and N). RELAX is run through a copy of the installed
`RELAX.bf` with `SetParameter (RANDOM_SEED, <seed>, 0);` as its first line, with `CPU=1`.

## Output folder

```
lr_out/
  setup.json        inputs (with sha256), options, engine version, number of alternatives, sampling
  runs.tsv          one row per run: label, seed, test-set size and length, length-matched, branches
  runs/<label>_s<seed>/
      aln.fasta  tree.nwk  RELAX_seeded.bf  cmd.sh  expected_sets.json
      modules -> HyPhy's SelectionAnalyses/modules (a link RELAX_seeded.bf needs)
      relax.out  relax.json  DONE        (after running)
  ALL_DONE          written when 'labelrank run' has finished every run folder
```

With `--engine codeml` a run folder holds `alt/` and `null/` (each with its control file, tree and `mlc`) in
place of the RELAX files.

`labelrank report --table` writes what was read from every run (LRT, P, K, best and null logL, the reported values,
warnings), so the reading can be recounted independently of this tool.

## Example

`examples/toy/` holds a 10-taxon, 300-codon simulated alignment (the same purifying selection on every branch, so
no label should stand out) with a labelled cherry, in both notations (`tree_relax.nwk`, `tree_codeml.nwk`). It has
3 same-shaped alternatives, so it demonstrates the mechanics, not a useful rank. `make_toy.py` regenerates it.

```bash
labelrank setup --alignment examples/toy/aln.fasta --tree examples/toy/tree_relax.nwk --test Foreground \
    --out toy_relax --hyphy /path/to/hyphy
labelrank run toy_relax --jobs 4          # 12 runs; a few minutes
labelrank report toy_relax
```

With HyPhy 2.5.101 every run finishes. The observed label gives LRT 2.43 (P 0.12, not significant) and ranks
1 of 4; the report flags that all 3 alternatives are shorter than the labelled cherry, and warns that one
alternative's logL differs by 39 units between seeds (an optimizer that stopped early on one restart: the reason
every label gets several seeds).

The same toy with codeml (`--engine codeml --tree examples/toy/tree_codeml.nwk --test '#1'`) shows what the tool is
for. The simulator is crude and follows neither program's model, and codeml's branch-site test calls the labelled
cherry under positive selection (LRT 37.8, P < 1e-9), but it calls every one of the 3 alternative cherries
significant as well (P from below 1e-9 to 2e-4), although no lineage was simulated differently. A chi-square P alone
would have reported a trait-linked result.

`validation/scnn1_A_siren/` is a full run on published data: the SCNN1A alignment and siren-labelled tree of
PMC12227717 (Source Data 6, Zenodo 10.5281/zenodo.15255678). The published RELAX result reproduces in 3 of 3
seeds and ranks first of 13 same-shaped labels (first of 3 length-matched); the folder holds the report, the
per-run table, the methods text and the figure (the inputs themselves are in the original deposit).

`validation/scnn1_otter_batch/` runs two more published tests from the same deposit (the sea otter tip in SCNN1A
and SCNN1G) as one batch with `--sequential 3`, one seed and 19 sampled alternatives each. Every run equals, to
the last digit, an earlier run of the same label and seed made without this tool. SCNN1A (LRT 0.73, not
significant here) stopped after 5 alternatives, 3 of them at or above it (sequential P 0.60). SCNN1G (LRT 4.79,
P 0.029) ran all 19: rank 2 of 20 overall, 1 of 8 length-matched, sequential P 0.10. On that tree length and LRT
correlate strongly among the alternatives (Spearman rho 0.74), which is the case the length-matched reading is
for.

`validation/scnn1_A_siren_codeml/` runs the same siren label through codeml branch-site model A on the real
alignment (PAML 4.10.10, one start, 9 sampled alternatives of 27; each fit pair took 17 to 23 hours on one
core). It is an end-to-end check of the codeml path on real data, not a second test of the result: the
observed label gives LRT 0 (P 1), and so do 8 of the 9 alternatives, which differ from 0 only by optimizer
noise (up to 2e-4); one gives LRT 0.76 (P 0.19). Ties count against the label, so it ranks 10 of 10. The run
folder was set up with v0.3.0 and read with v0.5.0.

`validation/scnn1_B_siren/` runs the dugong tip (*Dugong dugong*) in SCNN1B (published RELAX P 0.138, not significant) with two
seeds and 19 of 46 sampled alternatives, set up with v0.4.0 and read with v0.5.0. The two seeds give LRT 2.46 and 3.61; the
best fits give 3.53 (P 0.060), against 2.22 in the published output. The published null fit is 0.20 log-likelihood
units higher than either seed's, so here two seeds were not enough, and the report flags the spread across
seeds (0.54 for the alternative model) as above tolerance. Either way the label is not significant, and it ranks
4 of 20 (3 of 6 length-matched); 3 of the 19 alternatives would be called significant by the chi-square reading.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests use invented numbers in HyPhy's output layout, and a stand-in for the hyphy executable
(`tests/fake_hyphy.py`) for the end-to-end runs; they check the instrument (parsing, set checks, logL reading,
rank arithmetic, batch, stopping, scheduler files), not any biological result.

## When a run fails

A run whose program exits with an error (for example HyPhy's occasional "Internal error in _SimplexMethod") is
marked `FAIL` in its `DONE` file and excluded, with the reason, even if part of its output can be read; one that
passes `--timeout-hours` (default 48) is stopped, engine included, and marked `TIMEOUT`. To retry it, delete that
folder's `DONE` and run `labelrank run` again: a run clears its old outputs when it starts, and the other runs are
kept. To add seeds, use `labelrank add-seeds` rather than editing `RELAX_seeded.bf`.

Stopping `labelrank run` (Ctrl-C, `kill`, a closed terminal) stops the runs in progress, engines included, and
leaves them without a `DONE` marker; running it again starts them from scratch and skips the finished ones. Only
one `labelrank run` works on a folder at a time (a `RUNNING` file holds the host and process; a stale one, whose
process is gone, is taken over), and each run folder being worked on holds a `RUN_LOCK`, which `labelrank run`
and the scheduler's `run_one.sh` both respect, so a folder is never run twice at once. If the runner itself is
killed outright (`kill -9`), its lock names the engine's process too, so the folder stays locked until that
engine ends. If any run was left alone because another process held its `RUN_LOCK`, `labelrank run` exits
non-zero and writes no `ALL_DONE`. A `.takeover` file next to a lock is the short-lived lock used while a stale
lock is replaced; one older than 30 seconds is cleared automatically.

## Limits

- One test set per run, with tip, whole-clade and single-stem groups. A partly labelled clade is refused.
- RELAX with two branch sets (test, reference); other labelled branches are nuisance. codeml: branch-site model A
  or the two-ratio branch model, one foreground mark (`#1`); no `--direction` for codeml.
- Genetic codes: NCBI tables 1-6 and 9-15 (codeml: not 14).
- Checked against HyPhy 2.5.101 (RELAX.bf 4.7) and codeml 4.10.10 only; other releases are flagged.
- The rank is a calibration against this tree's other lineages (see "What the report says"). Few alternatives
  (small trees, large clades) give coarse ranks; the report says when N + 1 < 20. A large clade's exact topology
  usually occurs once on a tree, and then there is nothing to rank against: on the SCNN1 trees of the validation
  deposit the phocid, otariid and whale clades all have 0 alternatives, and `check` and `setup` stop with that
  message. The tool is for labels whose shape repeats (tips, small clades). Several separate groups reuse
  the same few lineages in many combinations, so their effective number is smaller than the count.
- Alternatives exclude any lineage not wholly inside the reference set (for example other labelled clades), and
  any lineage named with `--exclude`. An alternative lineage may itself be under selection; that makes the
  comparison conservative for your label, not anti-conservative, but it can hide a real trait signal. Leaving out
  lineages you expect to share the trait, named in advance, is the remedy.

## Related work

Calibrating a lineage-specific signal against other lineages of the same tree has precedents: control lineages
in tests for adaptive convergence (Zou and Zhang 2015, Mol Biol Evol; Thomas and Hahn 2015, Mol Biol Evol) and phylogenetic permulations for trait-associated evolutionary rates in RERconverge (Kowalczyk et al.
2019, Bioinformatics; Saputra et al. 2021, Mol Biol Evol). Permulations randomise the trait over
species while keeping its phylogenetic structure, for a rate-correlation statistic across many genes;
`labelrank` instead moves one label to every same-shaped place on the tree and reruns the branch-specific test
itself (RELAX, Wertheim et al. 2015, Mol Biol Evol; branch-site model A, Zhang et al. 2005, Mol Biol Evol;
two-ratio branch model, Yang 1998, Mol Biol Evol). Branch-site tests can return false positives on lineages with no positive selection, for example
through multinucleotide mutations (Venkat et al. 2018, Nat Ecol Evol). Sequential stopping
follows Besag and Clifford (1991, Biometrika).

## Citation

See `CITATION.cff`, which GitHub shows as "Cite this repository"; the archived release and its DOI are on
Zenodo. The methods paragraph cites the tool as "(Kim 2026)"; replace it with the full reference, including the
version and DOI you used.

## License

MIT
