# CRAN comments

## Resubmission after archival

`gmediation` (version 0.1.0) was previously on CRAN and was archived on
2020-05-19 because check problems were not corrected in time. This
submission (version 2.0.0) addresses [[NOTE: fill in specifically what the
original check problems were and how they were resolved, once you have
reproduced them with the current R CMD check output]].

## Summary of changes since the archived version

* Added a new exported function, `gmedsens()`, providing a probabilistic
  sensitivity analysis for path effects estimated by `gmediate()`.
* In `gmediate()`, the default value of `ref` changed from `NULL` to the
  string `"entire sample"`; behavior for the default case is unchanged
  (whole complete-case sample used as the reference group).
* `gmedsens()` gains `plot_out` and `table_out` options, and `gmediate()`
  gains `table_out`.
* Renamed several `gmediate()` arguments for consistency with `gmedsens()`:
  `bootsims` -> `boots`, `sig.level` -> `CIPCT`, `sens.par` -> `xcorrv`.

* Version 2.0.0 refactors `gmediate()` and `gmedsens()` into small internal
  functions with unit tests, adds S3 classes with `print()`/`plot()` methods
  (assigned results no longer print), drops six imports, and fixes several
  bugs (RR scale, negative-binomial dispersion, binary second-stage
  mediators in `gmedsens()`, missing-value handling); see NEWS.md for the full
  list and for the behaviour changes.

## Test environments

* Linux (Ubuntu 24.04), R 4.3.3: `R CMD check` (see results below).
* [[add: local macOS, win-builder, R-hub, etc. once run]]

## R CMD check results

* Linux, R 4.3.3: 0 errors | 0 warnings | 1 note. The only note was that the
  Suggested package `mediation` was not installed on the build machine.
* [[add results from win-builder / macOS / R-hub]]

## Downstream dependencies

* None (this package has no reverse dependencies on CRAN, as it was
  previously archived).
