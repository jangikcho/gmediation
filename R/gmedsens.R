#' Probabilistic Sensitivity Analysis for Generalized Causal Mediation
#'
#' Extends the mediation/path effect framework of \code{\link{gmediate}} with
#' a probabilistic sensitivity analysis. Path effects that involve
#' non-identifiable components (paths whose identification requires an
#' assumption about the correlation between residuals of two mediators, or
#' between a mediator and the outcome, that cannot be estimated from the
#' observed data) are re-evaluated across a range of assumed correlations, so
#' that the robustness of the substantive conclusions to this un-testable
#' assumption can be assessed.
#'
#' As in \code{\link{gmediate}}, the user supplies already-fit models for the
#' outcome and for each first- and second-stage mediator; these are re-fit
#' internally on the complete-case data, and bootstrap resampling is used to
#' construct confidence intervals. \code{gmedsens} additionally varies the
#' assumed correlation(s) between residuals for non-identifiable paths over
#' the grid of values given in \code{xcorrv}, and summarizes how the path
#' estimates and their confidence intervals change across that grid, both in
#' tabular form and in sensitivity plots (rendered using \pkg{ggplot2} and
#' \pkg{gridExtra}).
#'
#' As with \code{\link{gmediate}}, the default reference group, \code{ref =
#' "entire sample"}, uses the whole complete-case sample; a different
#' reference group may instead be specified as a vector giving the reference
#' group variable name and value(s).
#'
#' \code{gmedsens} returns an object of class \code{"gmedsens"} holding the
#' sensitivity tables and curves. The report (which files were written and
#' where, and the settings used) is printed when the analysis finishes, even
#' when the result is assigned to an object, unless \code{print_report =
#' FALSE}; use \code{print()} to show it again and \code{plot()} to draw the
#' sensitivity plots. When \code{plot_out[1] =
#' "true"} all sensitivity plots are also written to a single PDF in the
#' directory given in \code{plot_out[2]}, and, when \code{table_out[1] =
#' "true"}, the path-specific, first-stage, and second-stage sensitivity tables
#' are written as CSV files to the directory given in \code{table_out[2]}.
#'
#' @param data a data frame containing all variables used in the outcome and
#'   mediator models.
#' @param model.m1 a list of the model fit(s) for the first-stage mediator(s).
#' @param model.m2 a list of the model fit(s) for the second-stage
#'   mediator(s).
#' @param model.y an object with the model fit for the response (outcome)
#'   variable.
#' @param expos a vector with the name of the exposure variable and the
#'   exposure value, e.g. \code{c("exposure variable name", exposure value)}.
#' @param ref reference group specification. Either the default
#'   \code{"entire sample"}, meaning the whole complete-case sample is used
#'   as the reference group, or a vector with the name of the reference group
#'   variable and the reference group value(s), e.g. \code{c("reference group
#'   variable name", c(reference group value))}.
#' @param refmult numeric multiplier used to increase the number of Monte
#'   Carlo draws taken from the reference group, to reduce Monte Carlo error.
#' @param seed seed used to generate bootstrap samples.
#' @param boots number of bootstrap samples to generate.
#' @param cluster the cluster variable used for cluster bootstrap resampling,
#'   if applicable. Default is \code{NULL} (individual-level resampling).
#' @param decomp integer (1 or 2) specifying which of two available
#'   decompositions of the total exposure effect is used to define path
#'   effects.
#' @param CIPCT confidence level (as a percentage) used for the bootstrap
#'   confidence intervals.
#' @param xcorrv a vector of sensitivity parameter values (assumed
#'   correlations between residuals) evaluated for non-identifiable paths.
#' @param scale character string specifying the scale on which path effects
#'   are summarized (e.g. \code{"MeanDiff"}).
#' @param deltsi step of the grid of the sensitivity parameter \eqn{\rho}:
#'   \eqn{\rho} runs over \code{-1 + i * deltsi} (default 0.05, i.e. -0.95
#'   to 0.95). Must be in (0, 1) and give at least one grid point.
#' @param tol numeric convergence tolerance used in the root-finding /
#'   probability-matching steps of the sensitivity calculations.
#' @param tolx deprecated and ignored; kept so that existing calls still work.
#' @param assumvars character vector naming the mediator/outcome variables
#'   for which non-identifiable paths should be evaluated across
#'   \code{xcorrv}; \code{"ALL"} (the default) evaluates all applicable
#'   variables.
#' @param g numeric factor applied to scale the assumed between-mediator (or
#'   mediator-outcome) correlation grid.
#' @param n_rows number of rows of panels used when laying out sensitivity
#'   plots on a page.
#' @param n_cols number of columns of panels used when laying out sensitivity
#'   plots on a page.
#' @param ndig number of digits used when rounding values displayed in the
#'   sensitivity tables.
#' @param rows_per_page maximum number of table rows to display per printed
#'   page/table when results span many rows.
#' @param plot_out a length-2 character vector controlling whether the
#'   sensitivity plots are written to disk. The first element is
#'   \code{"true"} or \code{"false"} (the default), and the second is the
#'   directory in which to write the PDF \code{gmedsens_all_plots.pdf}
#'   (created if it does not exist) when the first element is \code{"true"}.
#'   The PDF contains the plots for every pair of assumption variables
#'   evaluated. If the Cairo graphics device is unavailable, a standard PDF
#'   device is used instead (in which case arrow glyphs in plot titles may
#'   render as dots).
#' @param table_out a length-2 character vector controlling whether the
#'   sensitivity summary tables are written to disk. The first element is
#'   \code{"true"} or \code{"false"} (the default), and the second is the
#'   directory (created if it does not exist) in which to write the CSV files
#'   \code{gmedsens_pathspec_out.csv}, \code{gmedsens_firststage_out.csv} and
#'   \code{gmedsens_secondstage_out.csv}. When several pairs of assumption
#'   variables are evaluated (\code{assumvars = "ALL"}), the files contain the
#'   results for every pair, identified by the \code{Rho.Between} column.
#' @param print_report print the report when the analysis finishes (default
#'   \code{TRUE}). It lists the output files and their location, so it is shown
#'   even when the result is assigned. Use \code{FALSE} to suppress it.
#' @param maxcnt maximum number of iterations allowed in the internal
#'   root-finding routines used for the sensitivity calculations.
#'
#' @return An object of class \code{"gmedsens"}, a list with elements
#'   \describe{
#'     \item{\code{path}}{data frame with the path-specific sensitivity
#'       results (columns \code{Path}, \code{Xcorr}, the sensitivity
#'       parameter \eqn{\rho}, \code{Estimate}, \code{Lower.CI},
#'       \code{Upper.CI} and \code{Rho.Between}); \code{NULL} if empty.}
#'     \item{\code{first.stage}}{data frame with the first-stage mediator
#'       results; \code{NULL} if empty.}
#'     \item{\code{second.stage}}{data frame with the second-stage mediator
#'       results; \code{NULL} if empty.}
#'     \item{\code{curves}}{long data frame with every sensitivity curve
#'       (estimate and confidence limits by \eqn{\rho}) used for the plots.}
#'     \item{\code{settings}}{list of the analysis settings.}
#'     \item{\code{files}}{paths of the files written (plots PDF, CSV
#'       tables), if any.}
#'   }
#'   Use \code{print()} for the report and \code{plot()} for the plots.
#'
#' @seealso \code{\link{gmediate}} for the underlying causal mediation and
#'   path analysis this function builds on.
#'
#' @examples
#' \donttest{
#' model.y2 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + visit +
#'                    seal + ohi, family = binomial, data = dental)
#' model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
#'                   family = binomial, data = dental)
#' model.m12 <- glm(visit ~ ses + grpvlb + grpbpd + race + SEX,
#'                   family = binomial, data = dental)
#' model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
#'                   family = binomial, data = dental)
#' model.m22 <- glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
#'                   family = gaussian, data = dental)
#'
#' gmedsens(data = dental, model.m1 = list(model.m11, model.m12),
#'          model.m2 = list(model.m21, model.m22), model.y = model.y2,
#'          expos = c("ses", 1), seed = 443095, boots = 100,
#'          cluster = "momid", decomp = 2, CIPCT = 95,
#'          xcorrv = c(0, 0.5, 0.9), scale = "MeanDiff", assumvars = c("dmftd",
#'          "ohi"), g = 1, ndig = 4, n_rows = 1, n_cols = 1,
#'          rows_per_page = 18, maxcnt = 100,
#'          plot_out = c("true", tempdir()), table_out = c("true", tempdir()))
#' }
#'
#' @importFrom stats predict quantile getCall pnorm pbinom ppois
#'   pnbinom rbinom rnbinom rnorm rpois runif dnorm qnorm complete.cases
#' @importFrom MASS glm.nb rnegbin
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs theme_bw theme
#'   scale_x_continuous element_text .data
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom utils write.csv
#' @export
gmedsens <- function(data, model.m1, model.m2, model.y, expos, ref = "entire sample",
                     refmult = 10, seed = 443095, boots = 499, cluster = NULL,
                     decomp = 2, CIPCT = 95, xcorrv = c(-0.9, 0, 0.5, 0.9),
                     scale = "MeanDiff", deltsi = 0.05, tol = 0.0001, tolx = 0.0001,
                     assumvars = "ALL", g = 1, n_rows = 1, n_cols = 1, ndig = 4,
                     rows_per_page = 16, plot_out = c("false", ""),
                     table_out = c("false", ""), maxcnt = 100,
                     print_report = TRUE) {
  envir <- parent.frame()

  # ---- set-up and validation -------------------------------------------
  n1 <- gm_count_mediators(model.m1)
  n2 <- gm_count_mediators(model.m2)
  gm_check_convergence(model.y, model.m1, model.m2, n1, n2)
  ref <- gm_normalise_ref(ref)
  gm_check_args(scale, boots, CIPCT, decomp, refmult, xcorrv)
  if (!(is.numeric(deltsi) && length(deltsi) == 1 && deltsi > 0 && deltsi <= 1 &&
        abs(2 / deltsi - round(2 / deltsi)) < 1e-8)) {
    stop("deltsi must be a number in (0, 1] for which 2 / deltsi is a whole number")
  }

  # g must be in (0, 1]; default to 1 if out of range
  if (g <= 0 | g > 1) {
    warning("Invalid g value. g must be in the range (0, 1]. Defaulting to g = 1.")
    g <- 1
  }

  restore_seed <- gm_set_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  spec <- gm_model_spec(model.y, model.m1, model.m2, n1, n2)
  gm_check_spec(spec, expos)
  if (!all(unlist(spec$family.m1) %in% c("gaussian", "binomial"))) {
    stop("gmedsens() supports only Gaussian and binomial first-stage mediators")
  }
  gdata <- gm_complete_data(data, spec, expos, cluster, ref)

  # ---- reference group, decomposition, assumption pairs -----------------
  data.cr <- gm_reference_data(gdata, ref, refmult, repeat_entire = TRUE)
  ctx <- gs_context(
    spec, expos, gm_prediction_data(data.cr, gdata, spec, expos, with_m2o = TRUE),
    data.cr, xcorrv, deltsi, g, tol, maxcnt, scale,
    gm_decomp_index(n1, n2, decomp), envir
  )
  pairs <- gs_assumption_pairs(assumvars, spec)
  plot_dir <- gs_output_dir(plot_out, "plot_out[2]", "to Save the Plots PDF")
  table_dir <- gs_output_dir(table_out, "table_out[2]", "to Save the Output Files")

  # ---- sensitivity analysis for each pair of assumptions ---------------
  st <- gs_new_state(ctx)
  probs <- gm_ci_probs(CIPCT)
  bootsamples <- gm_boot_samples(gdata, boots, cluster)
  accs <- gs_run_pairs(pairs, bootsamples, ctx, st, boots)
  results <- lapply(seq_len(nrow(pairs)), function(p) {
    gs_pair_outputs(accs[[p]], ctx, boots, probs, ndig, assumvars, pairs$i1[p], pairs$i2[p])
  })
  out <- gs_combine_pairs(results)

  # ---- files -------------------------------------------------------------
  files <- list()
  if (!is.null(table_dir)) {
    files$tables <- gs_write_csv(out, table_dir)
  }
  if (!is.null(plot_dir)) {
    plots_file <- file.path(plot_dir, "gmedsens_all_plots.pdf")
    written <- gs_write_pdf(out$curves, plots_file, CIPCT, n_rows, n_cols)
    if (is.null(written)) files$plots_failed <- plots_file else files$plots <- written
  }

  result <- structure(
    list(path = out$path, first.stage = out$first, second.stage = out$second,
         curves = out$curves,
         settings = list(expos = expos, y = spec$y, n = nrow(gdata), boots = boots,
                         ref = ref, refmult = refmult, scale = scale, cluster = cluster,
                         CIPCT = CIPCT, g = g, xcorrv = xcorrv, deltsi = deltsi,
                         assumvars = assumvars, n_rows = n_rows, n_cols = n_cols),
         files = files),
    class = "gmedsens"
  )
  if (isTRUE(print_report)) print(result)
  invisible(result)
}
