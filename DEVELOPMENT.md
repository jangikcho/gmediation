# Development notes for gmediation 2.0.0

(Not part of the built package; see `NEWS.md` for the user-facing changes.)

## Goals of the 2.0.0 refactor

Readable, maintainable code with shared helpers for `gmediate()` and
`gmedsens()`; a unit test for every helper; a structure that passes
`R CMD check` for CRAN; identical results to 1.0.0 except for documented bug
fixes. Changes were made in three stages, each verified before the next:

1. **Stage A** - restructure only; results identical to 1.0.0 (verified
   bit-for-bit on 30 fixed-seed configurations; the only difference was a stray
   empty PDF that 1.0.0 left behind on an error).
2. **Stage B** - bug fixes and input validation (missing values scoped to model
   variables, convergence check over all mediators, RR scale, negative-binomial
   dispersion, vectorised discrete draws, restoring the random seed, ...).
3. **Stage C** - `gmedsens()` shares bootstrap samples and model fits across
   assumption pairs; binary second-stage correction fixed.

## File layout

| File | Contents |
|---|---|
| `R/gmediate.R` | exported `gmediate()` (orchestrator, about 60 lines of code) |
| `R/gmedsens.R` | exported `gmedsens()` (orchestrator) |
| `R/setup.R` | `gm_count_mediators`, `gm_check_convergence`, `gm_family_name`, `gm_model_spec`, `gm_check_spec`, `gm_decomp_index`, `gm_normalise_ref`, `gm_complete_data`, `gm_check_args`, `gm_set_seed`, `gm_reference_data`, `gm_prediction_data`, `gm_boot_samples` |
| `R/simulate.R` | `gm_refit`, `gm_draw_m1`, `gm_draw_m2`, `gm_dispersion`, `gm_discrepant_draw` (Gaussian copula for non-identifiable paths) |
| `R/summarise.R` | `gm_boot_summary`, `gm_boot_pvalue`, `gm_ci_probs`, `gm_path_effect` (MeanDiff / OR / RR in one place), `gm_total_effect` |
| `R/engine-gmediate.R` | `gm_fit_boot_models`, `gm_boot_iteration`, `gm_run_bootstrap` |
| `R/engine-gmedsens.R` | `gs_context`, `gs_new_state`, `gs_new_arrays`, `gs_fit_boot_models`, `gs_draw_first_stage`, `gs_path_shown`, `gs_is_discrepant`, `gs_path_predictions`, `gs_correction`, `gs_outcome_mean`, `gs_expected_outcome`, `gs_scale_effect`, `gs_boot_iteration`, `gs_run_pairs` |
| `R/tables-gmediate.R`, `R/tables-gmedsens.R` | summary tables and sensitivity curves |
| `R/plots-gmedsens.R` | one ggplot builder, one paginator, PDF device with Cairo fallback, arrow-glyph handling |
| `R/print-*.R` | `print.gmediate`, `print.gmedsens`, `plot.gmedsens`, text-table utilities |
| `R/export.R` | xlsx / CSV writers, output-folder checks, assumption pairs |

Helpers are internal (`@noRd`, not exported) and are tested through
`gmediation:::`.

## Sensitivity grid (gmedsens)

`rho = -1 + i * deltsi` for `i = 1 .. round(2 / deltsi) - 1`. Assumption pairs
`(i1, i2)`: `i1 = 1` is the outcome, `i1 = k + 1` is second-stage mediator `k`;
`i2 <= n1` is a first-stage mediator, otherwise a second-stage mediator.

## Tests

* `tests/testthat/test-*.R`: unit tests for the helpers of every file, plus
  integration tests of `gmediate()` / `gmedsens()`.
* `test-print.R`: snapshot tests of printed output (skipped on CRAN).
* `test-golden.R`: reproduces frozen version 1 results
  (`fixtures/golden.rds`: gmediate decomp 1 and 2, Gaussian outcome, cluster
  bootstrap with two mediators per stage, gmedsens single pair) and the
  corrected results for the RR scale and for a second `gmedsens()` pair.
  Run with `NOT_CRAN=true` (skipped on CRAN).
* Coverage with `covr::package_coverage()`: about 97%.
* Golden files were produced in a C locale (`unset LANG LC_ALL`); the helper
  `compare_gmedsens` normalises the rho column name for other locales.

Run: `devtools::test()`; full: `Sys.setenv(NOT_CRAN = "true"); devtools::test()`.

## Verification against 1.0.0

Harness of 30 fixed-seed configurations (families gaussian, binomial, Poisson,
negative binomial; 1x1 and 2x2 mediators; cluster bootstrap; subset reference;
`decomp` 1 and 2; scales).

* After Stage A: `gmediate` 15/15 identical, `gmedsens` 14/15 identical.
* After Stages B and C: `gmediate` 12/15 identical; the 3 differences are the
  RR fix and two cases that crashed in 1.0.0. `gmedsens`: single-pair cases
  identical where no bug applied; others differ for the documented reasons
  (binary second-stage correction, per-partner proportions, Poisson
  second-stage standard deviation, shared resamples).

## Known issues not changed (see NEWS.md)

Cluster bootstrap by row; `refmult` handling differs between functions; OR/RR
not exponentiated in `gmedsens` mediator tables; `binom_power = 2` for
outcome versus first-stage mediator; start index of the Poisson outcome loop.

## Release checklist

1. `roxygen2::roxygenise()` (use a UTF-8 locale, e.g. `LANG=C.UTF-8`).
2. `devtools::check(args = "--as-cran")`, `devtools::check_win_devel()`.
3. Fill in the `[[...]]` placeholders in `cran-comments.md`; confirm the
   maintainer e-mail.
