# Changelog

## 0.6.0 (2026-11-10)

- **codeml two-ratio branch model.** `setup --engine codeml --model two-ratio` fits the two-ratio branch model
  (model = 2, NSsites = 0; the `#1` branches get their own omega) against the one-ratio model (model = 0,
  NSsites = 0). P is chi-square(1), since the two models are nested at the interior point omega foreground =
  omega background. The report's K column holds the foreground omega. The rest follows the codeml engine of
  0.5.0: CodonFreq = 2, cleandata = 0, each fit in its own folder, a run finished when both `mlc` files hold their
  likelihood line, the release warning. The one-ratio null takes the seed's starting omega as well as its kappa.
- `--model branch-site` is the default and writes the same control files as 0.5.0. Setups written before 0.6.0
  have no model in `setup.json` and are read as branch-site, so their reports do not change.
- `--model` is refused with `--engine relax`.

## 0.5.0 (2026-11-07)

First public release (Zenodo 10.5281/zenodo.24035697): HyPhy RELAX and codeml branch-site model A, rank over all
and over length-matched alternatives, Besag-Clifford sequential stopping, batch use and SLURM arrays.
