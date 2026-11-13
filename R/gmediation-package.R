#' gmediation: Mediation Analysis for Multiple and Multi-Stage Mediators
#'
#' Conducts causal mediation and path analysis for multiple mediators
#' arranged in up to two stages.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{gmediate}}}{estimates individual path effects,
#'     first- and second-stage mediator effects, and the total exposure
#'     effect, with bootstrap confidence intervals.}
#'   \item{\code{\link{gmedsens}}}{extends \code{gmediate} with a
#'     probabilistic sensitivity analysis for non-identifiable paths.}
#' }
#'
#' Both functions default to using the entire (complete case) sample as the
#' reference group (\code{ref = "entire sample"}) unless a different
#' reference group is specified.
#'
#' See the package vignette, \code{vignette("gmediate", package =
#' "gmediation")}, for a worked example using the built-in \code{\link{dental}}
#' data.
#'
#' @keywords internal
"_PACKAGE"

# Column names used via non-standard evaluation inside gmedsens()
utils::globalVariables(c("Mediator", "Path", "Xcorr"))
