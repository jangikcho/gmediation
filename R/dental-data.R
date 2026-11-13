#' Dental Data: gmediate/gmedsens Example Data
#'
#' A dataset of 205 children used throughout the package examples and
#' vignette to illustrate causal mediation and path analysis with
#' \code{\link{gmediate}} and \code{\link{gmedsens}}. The exposure of
#' interest is low socioeconomic status (\code{ses}); the outcome is dental
#' decay status (\code{dmftd}); \code{brush} and \code{visit} are first-stage
#' mediators; \code{seal} and \code{ohi} are second-stage mediators; and
#' \code{grpvlb}, \code{grpbpd}, \code{race}, and \code{SEX} are baseline
#' covariates. \code{momid} identifies the mother of each child and can be
#' used as a cluster variable for cluster bootstrap resampling.
#'
#' @format A data frame with 205 rows and 13 variables:
#' \describe{
#'   \item{ID}{character. Permuted subject ID.}
#'   \item{dmft}{numeric. Count of decayed, missing, or filled teeth.}
#'   \item{SEX}{binary (0/1). 1 if male; 0 if female.}
#'   \item{grpbpd}{binary (0/1). 1 if BPD (bronchopulmonary dysplasia); 0
#'     otherwise.}
#'   \item{grpvlb}{binary (0/1). 1 if VLBW (very low birth weight); 0
#'     otherwise.}
#'   \item{brush}{binary (0/1). 1 if brushes at least twice per day; 0
#'     otherwise. First-stage mediator.}
#'   \item{visit}{binary (0/1). 1 if the child is taken to the dentist at
#'     most once per year; 0 otherwise. First-stage mediator.}
#'   \item{dmftd}{binary (0/1). 1 if the dmft count is greater than 0; 0 if
#'     it equals 0. Outcome variable.}
#'   \item{seal}{binary (0/1). 1 if has at least one dental sealant; 0
#'     otherwise. Second-stage mediator.}
#'   \item{ses}{binary (0/1). 1 if low socioeconomic status; 0 otherwise.
#'     Exposure variable.}
#'   \item{race}{binary (0/1). 1 if African American; 0 otherwise.}
#'   \item{ohi}{numeric. Oral Hygiene Index (between 0 and 3, higher is
#'     worse). Second-stage mediator.}
#'   \item{momid}{character. Permuted maternal ID; can be used as a cluster
#'     variable.}
#' }
#'
#' @source Simulated/de-identified example data distributed with the
#'   gmediation package.
"dental"
