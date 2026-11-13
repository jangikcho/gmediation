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
