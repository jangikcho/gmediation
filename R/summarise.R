# Internal helpers: summarising bootstrap distributions and effect scales.

#' Percentile CI and two-sided bootstrap p-value of a bootstrap vector
#'
#' @param x numeric vector of bootstrap estimates (original sample excluded).
#' @param lo,hi lower / upper quantile probabilities.
#' @param na.rm passed to quantile().
#' @param back function mapping the scale of `x` to the reported scale
#'   (`exp` for OR / RR effects that are bootstrapped on the log scale). It is
#'   applied to the CI limits only: the p-value is always computed on the scale
#'   of `x`, where "no effect" is 0.
#' @return list(lower, upper, p): CI limits rounded to 4 decimals and the
#'   p-value as a string with 4 decimals.
#' @noRd
gm_boot_summary <- function(x, lo, hi, na.rm = TRUE, back = identity) {
  list(
    lower = unname(round(back(quantile(x, lo, na.rm = na.rm)), 4)),
    upper = unname(round(back(quantile(x, hi, na.rm = na.rm)), 4)),
    p = gm_boot_pvalue(x)
  )
}

#' Function that maps log-scale OR / RR effects to the reported scale
#' @param proportion proportions of the total effect stay on the bootstrap
#'   (log) scale because they are ratios of additive effects.
#' @noRd
gm_scale_back <- function(scale, proportion = FALSE) {
  if (scale %in% c("OR", "RR") && !proportion) exp else identity
}

#' Name of the effect column: "Estimate", or "OR" / "RR" for ratio scales
#' @noRd
gm_estimate_label <- function(scale, proportion = FALSE) {
  if (proportion) return("Est Prop")
  if (scale %in% c("OR", "RR")) scale else "Estimate"
}

#' Two-sided bootstrap p-value for H0: effect = 0
#' @noRd
gm_boot_pvalue <- function(x) {
  format(round(2 * min(ecdf(x)(0), ecdf(-x)(0)), 4), nsmall = 4)
}

#' Confidence-level quantile probabilities
#' @return c(lower, upper) probabilities for a CIPCT percent interval.
#' @noRd
gm_ci_probs <- function(CIPCT) {
  alpha <- 1 - (CIPCT / 100)
  c(alpha / 2, 1 - alpha / 2)
}

#' Path effect from two expected potential outcomes
#'
#' @param e_hi,e_lo expected outcomes (vectors or matrices) before and after
#'   switching the path.
#' @param scale "MeanDiff" (or ""), "OR" or "RR".
#' @noRd
gm_path_effect <- function(e_hi, e_lo, scale) {
  if (scale == "MeanDiff" || scale == "") return(e_hi - e_lo)
  if (scale == "OR") return(log((e_hi / (1 - e_hi)) / (e_lo / (1 - e_lo))))
  if (scale == "RR") return(log(e_hi / e_lo))
  e_hi * 0
}

#' Total effect from the first and last expected potential outcomes
#' @noRd
gm_total_effect <- function(e_all_exposed, e_all_control, scale) {
  if (scale == "MeanDiff" || scale == "") return(e_all_exposed - e_all_control)
  if (scale == "OR") {
    return(log((e_all_exposed / (1 - e_all_exposed)) /
                 (e_all_control / (1 - e_all_control))))
  }
  if (scale == "RR") return(log(e_all_exposed / e_all_control))
  e_all_exposed * 0
}
