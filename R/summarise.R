# Internal helpers: summarising bootstrap distributions and effect scales.

#' Percentile CI and two-sided bootstrap p-value of a bootstrap vector
#'
#' @param x numeric vector of bootstrap estimates (original sample excluded).
#' @param lo,hi lower / upper quantile probabilities.
#' @param na.rm passed to quantile().
#' @return list(lower, upper, p): CI limits rounded to 4 decimals and the
#'   p-value as a string with 4 decimals.
#' @noRd
gm_boot_summary <- function(x, lo, hi, na.rm = TRUE) {
  list(
    lower = unname(round(quantile(x, lo, na.rm = na.rm), 4)),
    upper = unname(round(quantile(x, hi, na.rm = na.rm), 4)),
    p = gm_boot_pvalue(x)
  )
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
