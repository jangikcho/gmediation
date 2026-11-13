# Internal helpers: refitting models on a bootstrap sample and simulating
# mediator values from their fitted distributions.

#' Refit a user model on a (bootstrap) data set
#'
#' glm.nb models are refit with MASS::glm.nb(); every other model is refit by
#' evaluating the stored call in `envir` (the caller's frame) so that the
#' family and link are resolved exactly as in the user's original call.
#' @noRd
gm_refit <- function(call, family, data, envir) {
  call$data <- data
  if (family == "negbin") {
    glm.nb(call$formula, call$data)
  } else {
    eval(call, envir)
  }
}

#' Simulate first-stage mediator values (one draw per reference-group row)
#'
#' Negative binomial draws use MASS::rnegbin(), as in version 1.
#' @noRd
gm_draw_m1 <- function(family, n, mu, disp) {
  if (family == "gaussian") return(rnorm(n, mean = mu, sd = sqrt(disp)))
  if (family == "poisson")  return(rpois(n, lambda = mu))
  if (family == "negbin")   return(rnegbin(n, mu = mu, theta = disp))
  if (family == "binomial") return(rbinom(n, size = 1, prob = mu))
  NULL
}

#' Simulate second-stage mediator values
#'
#' Negative binomial draws use stats::rnbinom(), as in version 1.
#' @noRd
gm_draw_m2 <- function(family, n, mu, disp) {
  if (family == "gaussian") return(rnorm(n, mean = mu, sd = sqrt(disp)))
  if (family == "poisson")  return(rpois(n, lambda = mu))
  if (family == "binomial") return(rbinom(n, 1, prob = mu))
  if (family == "negbin")   return(rnbinom(n, size = disp, mu = mu))
  NULL
}

#' Draw a first-stage mediator under the opposite exposure, correlated with
#' its factual draw (Gaussian copula) -- gmediate() version
#'
#' @param m1j factual draws; `mu_fact`, `mu_cf` factual / counterfactual means.
#' @return numeric vector of counterfactual draws.
#' @noRd
gm_discrepant_draw <- function(family, m1j, mu_fact, mu_cf, disp, rho, n) {
  if (family == "gaussian") {
    ZM1 <- (m1j - mu_fact) / sqrt(disp)
    ErrM1 <- rnorm(n, mean = 0, sd = 1)
    ZM12 <- (rho * ZM1) + sqrt(1 - rho**2) * ErrM1
    return((sqrt(disp) * ZM12) + mu_cf)
  }
  if (family == "binomial") {
    probjl <- pbinom(m1j - 1, 1, mu_fact)
    probj <- pbinom(m1j, 1, mu_fact)
  } else if (family == "poisson") {
    probjl <- ppois(m1j - 1, mu_fact)
    probj <- ppois(m1j, mu_fact)
  } else if (family == "negbin") {
    probjl <- pnbinom(m1j - 1, size = disp, mu = mu_fact)
    probj <- pnbinom(m1j, size = disp, mu = mu_fact)
  }
  # randomised probability integral transform of the factual draw
  UM1 <- runif(n, min = 0, max = 1)
  probM1 <- probjl + UM1 * (probj - probjl)
  ZM1 <- qnorm(probM1)
  ErrM1 <- rnorm(n, mean = 0, sd = 1)
  ZM12 <- (rho * ZM1) + sqrt(1 - rho**2) * ErrM1
  probM12 <- pnorm(ZM12)

  # quantile of the counterfactual distribution
  eps <- .Machine$double.eps
  probM12 <- pmin(pmax(probM12, eps), 1 - eps)
  if (family == "binomial") return(as.numeric(probM12 > pbinom(0, 1, mu_cf)))
  if (family == "poisson") return(qpois(probM12, mu_cf))
  qnbinom(probM12, size = disp, mu = mu_cf)
}

#' Dispersion of a fitted model
#'
#' Residual variance for Gaussian models, the shape parameter `theta` for
#' negative binomial models and `NULL` for the other families.
#' @noRd
gm_dispersion <- function(fit, family) {
  if (family == "gaussian") return(summary(fit)$dispersion)
  if (family == "negbin") return(fit$theta)
  NULL
}
