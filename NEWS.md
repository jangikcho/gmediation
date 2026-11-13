# gmediation 2.0.0

Version 2.0.0 is a rewrite of the internals of `gmediate()` and `gmedsens()`
(same arguments, same statistical method) together with bug fixes, S3
classes, and a test suite. Results are unchanged from 1.0.0 except where a
bug fix below says otherwise. Development notes, the file layout and the
verification against 1.0.0 are in `DEVELOPMENT.md` (source repository only).

## Behaviour changes you may notice

* **Printing.** `gmediate()` and `gmedsens()` now return objects of class
  `"gmediate"` and `"gmedsens"`.
  * `res <- gmediate(...)` no longer prints. Type `res` or `print(res)` to
    see the settings and tables, as before.
  * `gmedsens()` is the exception: its report (settings, the paths of the CSV
    tables and of the plots PDF) is always printed when the analysis
    finishes, even if the result is assigned, because the report says where
    the results were written. Use the new argument `print_report = FALSE` to
    suppress it; `print(res)` shows it again.
* **Return value of `gmedsens()`.** It used to return `NULL`. It now returns
  a `"gmedsens"` list with `path`, `first.stage` and `second.stage` (the
  tables written to the CSV files), `curves` (every plotted curve, in long
  format), `settings` and `files` (paths written). `gmediate()` results gain
  a `settings` element; the other elements are unchanged.
* **`plot(res)` for `gmedsens` objects** redraws the sensitivity plots on the
  current device.
* **`ref = NULL`** is accepted as the entire sample (it was the default
  before 1.0.0), in both functions; a `ref` that is neither `"entire sample"`,
  `NULL` nor `c(variable, value, ...)` is reported with a clear message.
* **`gmedsens()` rejects Poisson and negative-binomial first-stage
  mediators** with a clear error. Version 1 failed for them with the obscure
  error "the condition has length > 1". Binomial and Gaussian first-stage
  mediators, and all families for second-stage mediators and the outcome,
  are supported as before.
* **`gmedsens(tolx = )`** is deprecated and ignored (the iterative search it
  controlled was replaced by a direct computation); `tol` and `maxcnt` are
  kept.
* **Several assumption pairs share bootstrap samples.** With
  `assumvars = "ALL"`, version 1 redrew the bootstrap samples and refitted
  every model for each pair. Version 2 draws the samples once and fits each
  model once per sample, so all pairs use common resamples and the run is
  faster. Results for the second and later pairs therefore differ from
  version 1 at Monte Carlo noise level (about 0.01 to 0.02 at `boots = 3`,
  smaller for realistic `boots`). A single pair gives the same numbers as
  version 1 unless a bug fix below applies.
* **Missing values** are removed only for the variables that are used (outcome
  and mediators with their covariates, exposure, cluster and reference
  variables). Version 1 removed every row with a missing value in *any* column
  of `data` (one unused all-`NA` column left 5 of 205 rows). Sample sizes, and
  hence estimates, change when unused columns contain `NA`.
* **The user's random seed is restored** after each call (`seed` still makes
  the analysis reproducible).

## Bug fixes (numbers change)

* `gmediate(scale = "RR")`: all path effects were exactly 0 (the code took
  `log(E / E)`); they are now `log(E1 / E0)`.
* `gmediate()`, negative-binomial models: the dispersion used for the draws
  was wrong; it is now the model's `theta`. Count-valued (Poisson,
  negative-binomial) first-stage mediators crashed in version 1 with an NA
  error and a negative-binomial second-stage mediator crashed with
  "subscript out of bounds"; both now run. The slow search over a fixed
  1..100 grid that also truncated counts above 100 was replaced by the exact
  inverse CDF (`qpois` / `qnbinom`).
* `gmedsens()`, binary second-stage mediator: the correction was inverted,
  so that at `rho = 0` the mediator was drawn as 1 with probability
  `1 - mu` instead of `mu` (mean 0.70 instead of 0.30 at `mu = 0.3`).
  Sensitivity curves for pairs involving a binary second-stage mediator
  change, and the first-stage curves for different pairs now agree at
  `rho = 0`, as they should.
* `gmedsens()`: mediator proportions (`pm1`, `pm2`) were taken from the last
  partner model instead of the partner of each pair; the standard deviation
  for a Poisson second-stage mediator was computed from stale data. Both are
  fixed. `scale = "OR"` and `"RR"` no longer stop with "the condition has
  length > 1".
* The convergence check looped over the last mediator only; it now checks
  every mediator model.
* `gmedsens()` output folders: an empty folder name now gives a warning, and
  file paths no longer treat `{` `}` in a folder name specially (`glue` was
  replaced by `file.path()`).
* Invalid `boots`, `CIPCT`, `scale`, `decomp`, `refmult`, `xcorrv`, `deltsi`
  and `g` are reported with clear messages (an invalid `g` gives a warning and
  the default 1 is used). A model formula that uses a variable missing from
  `data` is reported by name.
* Plotting: `ggplot2::aes_string()` (deprecated) is no longer used. Drawing
  the arrow in plot titles on a device without the glyph no longer warns;
  `plot(res, ascii_arrows = TRUE)` or `options(gmediation.ascii_arrows =
  TRUE)` writes `->` instead, and the base `pdf()` fallback does so
  automatically. A stray empty PDF is no longer left behind when an analysis
  fails.

## Internal changes

* `gmediate()` (about 1,450 lines) and `gmedsens()` (about 2,270 lines) are
  now short functions that call about 75 internal helpers in 15 files
  (`setup.R`, `simulate.R`, `summarise.R`, `engine-*.R`, `tables-*.R`,
  `plots-gmedsens.R`, `print-*.R`, `export.R`). Code shared by both functions
  (model specification, reference data, bootstrap, draws, summaries) exists
  once. The random-number stream is unchanged, so results match version 1.
* Dependencies: removed the imports `plyr`, `VGAM`, `dplyr`, `kableExtra`,
  `glue` and `stringr`; `knitr` moved to Suggests (vignettes only). Added
  `testthat (>= 3.0.0)` and `covr` to Suggests. Imports are now MASS,
  ggplot2, gridExtra, grid, writexl, grDevices, utils and stats.
* Source files are ASCII only (non-ASCII characters are written as `\u`
  escapes).
* New testthat 3e suite (389 expectations, about 97% line coverage):
  unit tests for every helper, print snapshots, and regression tests against
  frozen version 1 results (`tests/testthat/fixtures/golden.rds`).

## Known issues, deliberately not changed in 2.0.0

These behave as in version 1 and are flagged for a later decision.

* Cluster bootstrap samples the cluster ids of individual rows rather than
  the unique clusters.
* `refmult` is ignored for the entire-sample reference group in `gmediate()`
  but applied in `gmedsens()`.
* With `scale = "OR"` / `"RR"`, `gmedsens()` first- and second-stage tables
  are not exponentiated (path tables are).
* The binomial-outcome correction uses a power of 2 for outcome versus
  first-stage mediator and 1 elsewhere.
* The Poisson-outcome loop in `gmedsens()` starts at the first second-stage
  index.

# gmediation 1.0.0

## New features

* Added `gmedsens()`, a probabilistic sensitivity analysis for `gmediate()`
  path effects. `gmedsens()` re-evaluates non-identifiable path effects
  across a user-specified grid of assumed correlations (`xcorrv`) and
  reports both tabular summaries and sensitivity plots, with optional
  CSV export via `table_out` and optional export of all plots to a single
  PDF via `plot_out`. Supports cluster bootstrap resampling via `cluster`.
* `gmedsens()` evaluates every pair of assumption variables when
  `assumvars = "ALL"` and exports all pairs (identified by the new
  `Rho.Between` column) to the CSV files and the plots PDF. Output
  folders given in `plot_out` / `table_out` are created if missing.
* `gmediate()` gains a `table_out` argument to write its summary tables to
  an Excel workbook, and its printed output now uses base-R formatted tables.
* New example data set `dental`.
* New vignette `vignette("gmedsens")` documenting `gmedsens()`; the
  `gmediate` vignette now points to it.

## Breaking changes

* In `gmediate()`, the default value of `ref` changed from `NULL` to
  `"entire sample"`. Behavior is unchanged: as before, the default is to
  use the whole complete-case sample as the reference group. Code that
  explicitly passed `ref = NULL` should be updated to omit `ref` or pass
  `ref = "entire sample"`.
* Several `gmediate()` arguments were renamed relative to the previously
  archived CRAN release (0.1.0): `bootsims` -> `boots`, `sig.level` ->
  `CIPCT`, `sens.par` -> `xcorrv`. The default `seed` also changed.

* The package now depends on R (>= 4.1.0) (native pipe `|>`), and imports
  `knitr`, `stringr` and `writexl`.

## Other changes

* This is the first CRAN submission since the package was archived on
  2020-05-19.
