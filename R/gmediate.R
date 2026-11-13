#' Generalized Causal Mediation and Path Analysis
#'
#' Estimates causal mediation and path effects for models with up to two
#' stages of mediation. The function allows an arbitrary number of mediators
#' at each stage, computes individual path effects, first- and second-stage
#' mediator effects, and the total exposure effect, and derives bootstrap
#' confidence intervals for each.
#'
#' This function carries out inference for mediation (path) effects based on
#' user-specified models. The function allows for two stages of mediation (up
#' to three links from exposure/treatment, X, to final outcome, Y), and an
#' arbitrary number of mediators at each stage. The user specifies an allowed
#' generalized regression model for each outcome (Y and each mediator).
#' Allowed distributions (link functions) include: normal (identity link),
#' Bernoulli (logit link), Poisson (log link), and negative binomial (log
#' link); canonical links are used for each distribution. Unsaturated models
#' (i.e., leaving out earlier stage mediators) are allowed. Baseline
#' (pre-exposure) covariates may be included in each model. The user-specified
#' models are re-fit in the function using the complete cases (based on all
#' the model variables).
#'
#' The user may also specify a reference group, effectively defining the
#' empirical multivariate baseline covariate distribution used in the
#' mediation formula; this may be viewed as defining the sub-population to
#' which inference is performed. The default reference group, \code{ref =
#' "entire sample"}, is the whole (complete case) sample. The reference group
#' may also be "cloned" by a specified integer multiple (\code{refmult}) to
#' reduce Monte Carlo error (at the possible cost of longer computing time).
#'
#' Bootstrap resampling is used to obtain bootstrap percentile method
#' confidence intervals (at the specified confidence level, \code{CIPCT}).
#' The option \code{cluster} can be used to specify a unit (aside from the
#' individual, that is, each data record, which is the default) for which
#' cluster bootstrap resampling is performed.
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
#' @param xcorrv a vector of sensitivity parameter values (correlations, on
#'   the latent/linear predictor scale) used when a path is non-identifiable.
#' @param scale character string specifying the scale on which path effects
#'   are summarized (e.g. \code{"MeanDiff"}).
#' @param CIPCT confidence level (as a percentage) used for the bootstrap
#'   confidence intervals.
#' @param decomp integer (1 or 2) specifying which of two available
#'   decompositions of the total exposure effect is used to define path
#'   effects; the ordering of mediators in the specified models also
#'   generally affects the decomposition.
#' @param table_out a length-2 character vector controlling whether the
#'   summary tables are written to disk. The first element is \code{"true"}
#'   or \code{"false"} (the default), and the second is the directory in
#'   which to write the Excel workbook \code{gmediate_summary_tables.xlsx}
#'   (created if it does not exist) when the first element is \code{"true"}.
#'
#' @return An object of class \code{"gmediate"}: a list with the following
#'   components. Nothing is printed when the result is assigned; use
#'   \code{print()} (or type the object's name) to display the analysis
#'   settings and the summary tables:
#' \describe{
#'   \item{indiv.path}{estimates and confidence intervals for individual path
#'     effects.}
#'   \item{total.path}{estimate and confidence interval for the total
#'     exposure effect.}
#'   \item{first.path}{estimates and confidence intervals for first-stage
#'     mediator effects.}
#'   \item{second.path}{estimates and confidence intervals for second-stage
#'     mediator effects.}
#'   \item{prop.indiv.path}{estimated proportions of the total effect
#'     mediated through each individual path.}
#'   \item{prop.first.path}{estimated proportions of the total effect
#'     mediated through each first-stage mediator.}
#'   \item{prop.second.path}{estimated proportions of the total effect
#'     mediated through each second-stage mediator.}
#'   \item{boot.est}{path effect estimates from all bootstrap samples.}
#'   \item{model.fit}{model fits from the complete case data used to fit
#'     each mediator and outcome model.}
#'   \item{settings}{the analysis settings used by the print method.}
#' }
#'
#' @seealso \code{\link{gmedsens}} for a probabilistic sensitivity analysis
#'   built on the same modeling framework.
#'
#' @examples
#' \donttest{
#' model.y1 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + seal,
#'                  family = binomial, data = dental)
#' model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
#'                   family = binomial, data = dental)
#' model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush,
#'                   family = binomial, data = dental)
#'
#' single <- gmediate(data = dental, model.m1 = list(model.m11),
#'                     model.m2 = list(model.m21), model.y = model.y1,
#'                     expos = c("ses", 1), seed = 1234, boots = 100,
#'                     cluster = NULL, decomp = 1, CIPCT = 95,
#'                     xcorrv = c(0, 0.9), refmult = 1)
#' single
#' names(single)
#' }
#'
#' @importFrom stats predict quantile ecdf getCall pnorm pbinom ppois
#'   pnbinom rbinom rnbinom rnorm rpois runif qnorm qpois qnbinom dnorm
#'   complete.cases
#' @importFrom MASS glm.nb rnegbin
#' @export
gmediate <- function(data, model.m1, model.m2, model.y, expos, ref = "entire sample",
                     refmult = 100, seed = 443095, boots = 499, cluster = NULL,
                     xcorrv = c(0, 0.5, 0.9), scale = "MeanDiff", CIPCT = 95,
                     decomp = 2, table_out = c("false", "")) {
  envir <- parent.frame()

  # ---- set-up and validation -------------------------------------------
  n1 <- gm_count_mediators(model.m1)
  n2 <- gm_count_mediators(model.m2)
  gm_check_convergence(model.y, model.m1, model.m2, n1, n2)
  ref <- gm_normalise_ref(ref)
  gm_check_args(scale, boots, CIPCT, decomp, refmult, xcorrv)
  restore_seed <- gm_set_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  spec <- gm_model_spec(model.y, model.m1, model.m2, n1, n2)
  gm_check_spec(spec, expos)
  gdata <- gm_complete_data(data, spec, expos, cluster, ref)

  # ---- reference group, decomposition, bootstrap samples ---------------
  data.cr <- gm_reference_data(gdata, ref, refmult)
  ctx <- list(
    spec = spec, expos = expos, xcorrv = xcorrv, nrho = length(xcorrv),
    scale = scale, CIPCT = CIPCT, envir = envir, nref = nrow(data.cr),
    pred = gm_prediction_data(data.cr, gdata, spec, expos),
    indexDM = gm_decomp_index(n1, n2, decomp)
  )
  bootsamples <- gm_boot_samples(gdata, boots, cluster)

  # ---- bootstrap estimation --------------------------------------------
  res <- gm_run_bootstrap(bootsamples, ctx)

  # ---- summary tables --------------------------------------------------
  probs <- gm_ci_probs(CIPCT)
  ESTT <- round(res$EDT[boots + 1, ], 7)
  tables <- list(
    indiv.path = gm_path_table(res, ctx, boots, probs, ESTT),
    total.path = gm_total_table(res$EDT, boots, probs, CIPCT),
    first.path = gm_first_stage_table(res, ctx, boots, probs, ESTT),
    second.path = gm_second_stage_table(res, ctx, boots, probs, ESTT),
    prop.indiv.path = gm_path_table(res, ctx, boots, probs, ESTT, proportion = TRUE),
    prop.first.path = gm_first_stage_table(res, ctx, boots, probs, ESTT, proportion = TRUE),
    prop.second.path = gm_second_stage_table(res, ctx, boots, probs, ESTT, proportion = TRUE)
  )

  boot_est <- lapply(res$EDM, function(m) m[-(boots + 1), ])
  model_fit <- gm_model_fit_summary(res$last$fit, spec)
  gm_write_tables_xlsx(tables, table_out)

  structure(
    c(tables,
      list(boot.est = boot_est,
           model.fit = model_fit,
           settings = list(expos = expos, y = spec$y, n = nrow(gdata), boots = boots,
                           ref = ref, refmult = refmult, scale = scale,
                           cluster = cluster, CIPCT = CIPCT))),
    class = "gmediate"
  )
}
