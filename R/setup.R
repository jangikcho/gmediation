# Internal helpers shared by gmediate() and gmedsens(): counting mediators,
# checking models, extracting the model specification, building the reference
# group and bootstrap samples. None of these functions is exported.

#' Number of mediators supplied for one stage
#'
#' @param model_m a list of fitted models (one per mediator).
#' @return integer count.
#' @noRd
gm_count_mediators <- function(model_m) {
  if (identical(class(model_m)[1], "list")) {
    if (is.null(names(model_m))) {
      length(model_m)
    } else {
      sum(names(model_m) == "y")
    }
  } else {
    1L
  }
}

#' Stop when any supplied model reports a convergence warning
#' @noRd
gm_check_convergence <- function(model.y, model.m1, model.m2, n1, n2) {
  msg <- "please check your model convergence"
  has_warn <- function(m) "th.warn" %in% names(m)

  if (has_warn(model.y)) stop(msg)

  check_stage <- function(model_m, n) {
    if (identical(class(model_m)[1], "list")) {
      if (is.null(names(model_m))) {
        for (i in seq_len(n)) if (has_warn(model_m[[i]])) stop(msg)
      } else if (has_warn(model_m)) {
        stop(msg)
      }
    } else if (has_warn(model_m)) {
      stop(msg)
    }
  }
  check_stage(model.m1, n1)
  check_stage(model.m2, n2)
  invisible(TRUE)
}

#' Distribution family of a fitted glm / glm.nb, with "negbin" for glm.nb
#' @noRd
gm_family_name <- function(model) {
  fam <- model$family$family
  if (length(grep("Negative Binomial", substr(fam, 1, 17))) == 1) fam <- "negbin"
  fam
}

#' Extract everything the algorithms need from the user's fitted models
#'
#' @return a list with the variable names, covariate sets, distribution
#'   families and calls of the outcome and mediator models.
#' @noRd
gm_model_spec <- function(model.y, model.m1, model.m2, n1, n2) {
  yvar <- names(model.y$model)
  y <- yvar[1]
  y.allcov <- yvar[-1]

  m2 <- vector("list", n2)
  m2.eachcov <- vector("list", n2)
  for (i in seq_len(n2)) {
    m2[[i]] <- names(model.m2[[i]]$model)[1]
    m2.eachcov[[i]] <- names(model.m2[[i]]$model)[-1]
  }
  m1 <- vector("list", n1)
  m1.eachcov <- vector("list", n1)
  for (i in seq_len(n1)) {
    m1[[i]] <- names(model.m1[[i]]$model)[1]
    m1.eachcov[[i]] <- names(model.m1[[i]]$model)[-1]
  }

  m2.allcov <- unique(unlist(m2.eachcov))
  m1.allcov <- unique(unlist(m1.eachcov))

  # covariates of the Y model that are not mediators; same for the M2 models
  y.cov <- y.allcov[!(y.allcov %in% c(m2, m1))]
  m2.cov <- m2.allcov[!(m2.allcov %in% m1)]

  list(
    n1 = n1, n2 = n2, npath = (n1 + 1) * (n2 + 1),
    y = y, y.allcov = y.allcov, y.cov = y.cov,
    m1 = m1, m1.eachcov = m1.eachcov, m1.allcov = m1.allcov,
    m2 = m2, m2.eachcov = m2.eachcov, m2.allcov = m2.allcov, m2.cov = m2.cov,
    family.y = gm_family_name(model.y),
    family.m1 = lapply(model.m1[seq_len(n1)], gm_family_name),
    family.m2 = lapply(model.m2[seq_len(n2)], gm_family_name),
    call.y = getCall(model.y),
    call.m1 = lapply(model.m1[seq_len(n1)], getCall),
    call.m2 = lapply(model.m2[seq_len(n2)], getCall)
  )
}

#' Validate the mediation model structure (same checks and messages as v1)
#' @noRd
gm_check_spec <- function(spec, expos) {
  n1 <- spec$n1
  n2 <- spec$n2
  y <- spec$y
  `%ni%` <- Negate(`%in%`)

  if (n1 * n2 == 0) stop("Need at least 1 mediator for each stage")

  for (k in seq_len(n2)) {
    if (expos[1] == spec$m2[[k]]) stop("Exposure cannot be mediator")
    if (y == spec$m2[[k]]) stop("Final outcome cannot be mediator")
    if (!any(unlist(spec$m1) %in% spec$m2.eachcov[[k]])) {
      stop("At least one M1 must be in each M2 model")
    }
    if (spec$m2[[k]] %ni% spec$y.allcov) stop("Each M2 must be in Y model")
  }

  for (j in seq_len(n1)) {
    if (expos[1] %ni% spec$m1.eachcov[[j]]) stop("Exposure must be in each M1 model")
    if (expos[1] == spec$m1[[j]]) stop("Exposure cannot be mediator")
    if (y == spec$m1[[j]]) stop("Final outcome cannot be mediator")
  }

  for (j in seq_len(n1)) {
    in_some_m2 <- any(vapply(spec$m2.eachcov, function(cv) spec$m1[[j]] %in% cv, logical(1)))
    if (!(in_some_m2 || spec$m1[[j]] %in% spec$y.cov)) {
      stop("Each M1 must be a covariate in the Y or some M2 model")
    }
  }
  invisible(TRUE)
}

#' Order in which the path effects are peeled off (decomp = 1 or 2)
#'
#' @return a 2 x npath integer matrix; column c holds the (M1, M2) cell that is
#'   switched from exposed to control when path effect c is evaluated.
#' @noRd
gm_decomp_index <- function(n1, n2, decomp) {
  npath <- (n1 + 1) * (n2 + 1)
  idx <- matrix(NA, 2, npath)
  cnt <- 1
  for (j in 1:(n1 + 1)) {
    for (k in 1:(n2 + 1)) {
      idx[, cnt] <- c(j, k)
      cnt <- cnt + 1
    }
  }
  if (decomp == 2) idx <- idx[, npath:1]
  idx
}

#' Reference-group data set (optionally cloned `refmult` times)
#'
#' @param repeat_entire if FALSE (gmediate), `refmult` is ignored when
#'   `ref = "entire sample"`, as in version 1.
#' @noRd
gm_reference_data <- function(gdata, ref, refmult, repeat_entire = FALSE) {
  ref <- gm_normalise_ref(ref)
  if (identical(ref, "entire sample")) {
    if (!repeat_entire) return(gdata)
    data.cr0 <- gdata
  } else {
    data.cr0 <- c()
    for (i in 2:length(ref)) {
      data.cr0 <- rbind(data.cr0, gdata[gdata[, ref[1]] == ref[i], ])
    }
  }
  if (refmult == 1) return(data.cr0)
  data.cr <- c()
  for (i in 1:refmult) data.cr <- rbind(data.cr0, data.cr)
  data.cr
}

#' Prediction data sets for the reference group
#'
#' @param with_m2o also build the data for M2 under the all-exposed setting
#'   (needed by gmedsens()).
#' @return list(pred.dat.y, pred.dat.m2, pred.dat.m2o, pred.dat.m1.e,
#'   pred.dat.m1.c, exposed, control).
#' @noRd
gm_prediction_data <- function(data.cr, gdata, spec, expos, with_m2o = FALSE) {
  pred.dat.y <- data.cr[, spec$y.allcov]
  pred.dat.m2 <- vector("list", spec$n2)
  pred.dat.m2o <- vector("list", spec$n2)
  for (k in seq_len(spec$n2)) {
    pred.dat.m2[[k]] <- data.cr[, spec$m2.eachcov[[k]]]
    if (with_m2o) pred.dat.m2o[[k]] <- data.cr[, spec$m2.eachcov[[k]]]
  }

  lexp <- levels(as.factor(gdata[, expos[1]]))
  exposed <- as.numeric(expos[2])
  control <- as.numeric(lexp[lexp != expos[2]])

  pred.dat.m1.e <- vector("list", spec$n1)
  pred.dat.m1.c <- vector("list", spec$n1)
  for (j in seq_len(spec$n1)) {
    pred.dat.m1.e[[j]] <- pred.dat.m1.c[[j]] <- data.frame(data.cr[, spec$m1.eachcov[[j]]])
    colnames(pred.dat.m1.e[[j]]) <- c(spec$m1.eachcov[[j]])
    colnames(pred.dat.m1.c[[j]]) <- c(spec$m1.eachcov[[j]])
    pred.dat.m1.e[[j]][, expos[1]] <- exposed
    pred.dat.m1.c[[j]][, expos[1]] <- control
  }
  list(pred.dat.y = pred.dat.y, pred.dat.m2 = pred.dat.m2, pred.dat.m2o = pred.dat.m2o,
       pred.dat.m1.e = pred.dat.m1.e, pred.dat.m1.c = pred.dat.m1.c,
       exposed = exposed, control = control)
}

#' Bootstrap samples (the last element is the original sample)
#'
#' Consumes the random-number stream exactly as version 1 did.
#' @return list of `boots + 1` data frames.
#' @noRd
gm_boot_samples <- function(gdata, boots, cluster = NULL) {
  bootsamples <- vector("list", boots + 1)
  if (is.null(cluster)) {
    nclus <- nrow(gdata)
    gdata$numsamp <- seq_len(nclus)
    for (i in 1:boots) {
      samp <- sample(gdata$numsamp, size = nclus, replace = TRUE)
      bootsamples[[i]] <- gdata[samp, ]
    }
  } else {
    nclus <- length(unique(gdata[, cluster]))
    cl <- gdata[, cluster]
    for (i in 1:boots) {
      samp <- sample(cl, size = nclus, replace = TRUE)
      rows <- unlist(lapply(samp, function(s) which(cl == s)), use.names = FALSE)
      bootsamples[[i]] <- gdata[rows, ]
    }
  }
  bootsamples[[boots + 1]] <- gdata
  bootsamples
}

#' Complete cases on the variables the analysis uses
#'
#' Only the variables of the models, the exposure, the cluster variable and the
#' reference-group variable need to be observed; missing values in other
#' columns of `data` no longer remove rows.
#' @noRd
gm_complete_data <- function(data, spec, expos, cluster = NULL, ref = "entire sample") {
  ref <- gm_normalise_ref(ref)
  vars <- unique(c(spec$y, spec$y.allcov, unlist(spec$m1), unlist(spec$m1.eachcov),
                   unlist(spec$m2), unlist(spec$m2.eachcov), expos[1], cluster,
                   if (!identical(ref, "entire sample")) ref[1]))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("variables not found in `data`: ", paste(missing_vars, collapse = ", "))
  }
  data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
}

#' Check the scalar arguments shared by gmediate() and gmedsens()
#' @noRd
gm_check_args <- function(scale, boots, CIPCT, decomp, refmult, xcorrv) {
  if (!(is.character(scale) && length(scale) == 1 && scale %in% c("MeanDiff", "OR", "RR", ""))) {
    stop("scale must be one of \"MeanDiff\", \"OR\" or \"RR\"")
  }
  if (!(is.numeric(boots) && length(boots) == 1 && boots >= 1 && boots == round(boots))) {
    stop("boots must be a positive whole number")
  }
  if (!(is.numeric(CIPCT) && length(CIPCT) == 1 && CIPCT > 0 && CIPCT < 100)) {
    stop("CIPCT must be a number between 0 and 100")
  }
  if (!(length(decomp) == 1 && decomp %in% c(1, 2))) stop("decomp must be 1 or 2")
  if (!(is.numeric(refmult) && length(refmult) == 1 && refmult >= 1 && refmult == round(refmult))) {
    stop("refmult must be a positive whole number")
  }
  if (!(is.numeric(xcorrv) && length(xcorrv) >= 1 && all(abs(xcorrv) <= 1))) {
    stop("xcorrv must be numeric values between -1 and 1")
  }
  invisible(TRUE)
}

#' Seed the random number generator and restore the user's state afterwards
#'
#' Call as `on.exit(gm_set_seed(seed), add = TRUE)`-style: the returned
#' function restores `.Random.seed` as it was before the call.
#' @return a function that restores the previous random number state.
#' @noRd
gm_set_seed <- function(seed) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE)
  set.seed(seed)
  function() {
    if (had) {
      assign(".Random.seed", old, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }
}

#' Normalise the reference-group argument
#'
#' `NULL` (the version 1.0.0 default of gmediate) means the entire sample.
#' Otherwise `ref` must be c(variable, value, ...).
#' @noRd
gm_normalise_ref <- function(ref) {
  if (is.null(ref) || identical(ref, "entire sample")) return("entire sample")
  if (length(ref) < 2) {
    stop("`ref` must be \"entire sample\" (or NULL) or c(variable, value1, ...)",
         call. = FALSE)
  }
  ref
}
