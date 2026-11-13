# Internal engine of gmedsens(): for one pair of assumptions (i1, i2) it runs
# the bootstrap and, for every path, every correlation `xcorr` and every value
# of the sensitivity parameter rho, computes the expected potential outcome of
# Y. Path effects are differences (or log ratios) of consecutive expected
# potential outcomes.
#
# Assumption pair (i1, i2):
#   i1 = 1      : no unmeasured confounding between Y and ...
#   i1 = k + 1  : no unmeasured confounding between the k-th M2 and ...
#   i2 <= n1    : ... the i2-th first-stage mediator
#   i2 > n1     : ... the (i2 - n1)-th second-stage mediator (only with i1 = 1)
#
# Result arrays (all for one assumption pair; dimensions
# boots + 1 bootstrap samples (the last is the original sample) x nrcnt values
# of rho x nxcorr values of the cross-world correlation):
#   EYDM[[p]]   expected potential outcomes, p = 1 .. npath + 1
#   EDM_SC[[c]] path effects on the requested scale, c = 1 .. npath
#   EDM1_SC[[j]] total effect through first-stage mediator j
#   EDM2_SC[[k]] total effect through second-stage mediator k

# ---- context and state ---------------------------------------------------

gs_context <- function(spec, expos, pred, data.cr, xcorrv, deltsi, g, tol,
                       maxcnt, scale, indexDM, envir) {
  nrhosi <- 2 / deltsi # goes from -1 to 1
  list(
    spec = spec, expos = expos, pred = pred, nref = dim(data.cr)[1],
    exposed = pred$exposed, control = pred$control,
    xcorrv = xcorrv, nxcorr = length(xcorrv),
    deltsi = deltsi, nrcnt = nrhosi - 1,
    g = g, tol = tol, maxcnt = maxcnt, scale = scale,
    indexDM = indexDM, npath = spec$npath, envir = envir
  )
}

# Mutable state: the prediction data sets are modified in place along the paths
# (every path resets all columns it uses before reading them).
gs_new_state <- function(ctx) {
  st <- new.env(parent = emptyenv())
  st$pred.dat.y <- ctx$pred$pred.dat.y
  st$pred.dat.m2 <- ctx$pred$pred.dat.m2
  st$pred.dat.m2o <- ctx$pred$pred.dat.m2o
  st$xcorr <- c()
  st
}

gs_new_arrays <- function(ctx, boots) {
  n1 <- ctx$spec$n1
  n2 <- ctx$spec$n2
  npath <- ctx$npath
  nrcnt <- ctx$nrcnt
  nxcorr <- ctx$nxcorr
  dim3 <- c(boots + 1, nrcnt, nxcorr)
  acc <- new.env(parent = emptyenv())
  acc$EYDM <- lapply(seq_len(npath + 1), function(i) array(0, dim3))
  acc$EDM_SC <- lapply(seq_len(npath), function(i) array(0, dim3))
  acc$EDM1_SC <- lapply(seq_len(n1), function(i) matrix(0, boots + 1, nrcnt))
  acc$EDM2_SC <- lapply(seq_len(n2), function(i) array(0, dim3))
  acc$EDMM1 <- array(0, c(npath, nrcnt, boots + 1))
  acc$EDMM2 <- array(0, c(npath, nrcnt, boots + 1, nxcorr))
  acc$pathpr <- rep(1, npath)
  acc$nullv <- rep(0, npath)
  acc$xcrep <- rep(1, npath + 1)
  acc
}

gs_get <- function(lst, i) if (i <= length(lst)) lst[[i]] else NULL

# ---- model fitting -------------------------------------------------------

# Refit the Y, M1 and M2 models on one bootstrap sample and record the
# standard deviations / dispersions used by the sensitivity corrections. The
# fits do not depend on the assumption pair, so they are made once per sample.
gs_fit_boot_models <- function(bdata, ctx) {
  sp <- ctx$spec
  fit <- list(m1 = vector("list", sp$n1), m2 = vector("list", sp$n2),
              mu.m1.e = vector("list", sp$n1), mu.m1.c = vector("list", sp$n1),
              disp.m1 = list(), disp.m2 = list(), stdv.m1 = list(), stdv.m2 = list(),
              pm1 = list(), pm2 = list())

  for (k in seq_len(sp$n2)) {
    fam <- sp$family.m2[[k]]
    m <- gm_refit(sp$call.m2[[k]], fam, bdata[, c(sp$m2[[k]], sp$m2.eachcov[[k]])], ctx$envir)
    fit$m2[[k]] <- m
    if (fam == "binomial") {
      # crude marginal proportion
      pm2 <- mean(bdata[, c(sp$m2[[k]])])
      fit$pm2[[k]] <- pm2
      fit$stdv.m2[[k]] <- sqrt(pm2 * (1 - pm2))
    } else if (fam == "gaussian" || fam == "negbin") {
      fit$disp.m2[[k]] <- gm_dispersion(m, fam)
      fit$stdv.m2[[k]] <- if (fam == "gaussian") sqrt(fit$disp.m2[[k]]) else fit$disp.m2[[k]]
    } else if (fam == "poisson") {
      mu <- predict(m, type = "response", newdata = ctx$pred$pred.dat.m2[[k]])
      fit$stdv.m2[[k]] <- sqrt(mu)
    }
  }

  for (j in seq_len(sp$n1)) {
    fam <- sp$family.m1[[j]]
    m <- gm_refit(sp$call.m1[[j]], fam, bdata[, c(sp$m1[[j]], sp$m1.eachcov[[j]])], ctx$envir)
    fit$m1[[j]] <- m
    if (fam == "binomial") {
      pm1 <- mean(bdata[, c(sp$m1[[j]])])
      fit$pm1[[j]] <- pm1
      fit$stdv.m1[[j]] <- sqrt(pm1 * (1 - pm1))
    } else if (fam == "gaussian") {
      fit$disp.m1[[j]] <- gm_dispersion(m, fam)
      fit$stdv.m1[[j]] <- sqrt(fit$disp.m1[[j]])
    }
    fit$mu.m1.e[[j]] <- predict(m, type = "response", newdata = ctx$pred$pred.dat.m1.e[[j]])
    fit$mu.m1.c[[j]] <- predict(m, type = "response", newdata = ctx$pred$pred.dat.m1.c[[j]])
  }

  fit$y <- gm_refit(sp$call.y, sp$family.y, bdata[, c(sp$y, sp$y.allcov)], ctx$envir)
  if (sp$family.y == "negbin") {
    fit$disp.y <- gm_dispersion(fit$y, "negbin")
  } else if (sp$family.y == "gaussian") {
    fit$disp.y <- gm_dispersion(fit$y, "gaussian")
    fit$stdv.y <- sqrt(fit$disp.y)
  } else if (sp$family.y == "binomial") {
    py <- mean(bdata[, c(sp$y)])
    fit$stdv.y <- sqrt(py * (1 - py))
  }
  fit
}

# Draw first-stage mediator values for everyone in the reference group, once
# as if all were exposed and once as if all were control.
gs_draw_first_stage <- function(fit, ctx) {
  sp <- ctx$spec
  e <- vector("list", sp$n1)
  cc <- vector("list", sp$n1)
  for (j in seq_len(sp$n1)) {
    fam <- sp$family.m1[[j]]
    disp <- if (fam == "gaussian") fit$disp.m1[[j]] else NULL
    e[[j]] <- gm_draw_m1(fam, ctx$nref, fit$mu.m1.e[[j]], disp)
    cc[[j]] <- gm_draw_m1(fam, ctx$nref, fit$mu.m1.c[[j]], disp)
  }
  list(e = e, c = cc)
}

# ---- one path ------------------------------------------------------------

# Should a path be shown for assumption (i1, i2)? Paths whose two decompositions
# are both equal (or both unequal) are not affected by the assumption: the
# effects cancel exactly or approximately.
gs_path_shown <- function(DM, PrevDM, i1, i2, n1) {
  same <- function(a, b, a0, b0) {
    (a == b & a0 == b0) | (a != b & a0 != b0)
  }
  if (i1 == 1 && i2 < n1 + 1) {
    if (same(DM[1, 1], DM[i2 + 1, 1], PrevDM[1, 1], PrevDM[i2 + 1, 1])) return(0)
  }
  if (i1 == 1 && i2 > n1) {
    col <- i2 - n1 + 1
    if (same(DM[1, 1], DM[1, col], PrevDM[1, 1], PrevDM[1, col])) return(0)
  }
  if (i1 > 1) {
    if (same(DM[1, i1], DM[i2 + 1, i1], PrevDM[1, i1], PrevDM[i2 + 1, i1])) return(0)
  }
  1
}

# Is any first-stage mediator "discrepant" in DM (set differently for the
# second-stage mediator than for the exposure -> M1 link)?
gs_is_discrepant <- function(DM, ctx) {
  sp <- ctx$spec
  for (jj in 2:(sp$n1 + 1)) {
    for (kk in 2:(sp$n2 + 1)) {
      if (sp$m1[[jj - 1]] %in% sp$m2.eachcov[[kk - 1]] & DM[jj, kk] != DM[jj, 1]) {
        return(TRUE)
      }
    }
  }
  FALSE
}

# Mediator means / draws for the current decomposition matrix DM, and the
# updated prediction data for Y and the M2 models.
gs_path_predictions <- function(fit, draws, DM, jd, ctx, st, i1, i2) {
  sp <- ctx$spec
  expos <- ctx$expos
  nref <- ctx$nref
  n1 <- sp$n1
  n2 <- sp$n2
  mu.m1 <- pred.m1 <- pred.m1c <- mu.m12 <- mu.m1.d00 <- vector("list", n1)

  for (j in seq_len(n1)) {
    j1 <- j + 1
    mu.m1[[j]] <- DM[j1, 1] * fit$mu.m1.e[[j]] + (1 - DM[j1, 1]) * fit$mu.m1.c[[j]]
    pred.m1[[j]] <- DM[j1, 1] * draws$e[[j]] + (1 - DM[j1, 1]) * draws$c[[j]]
    if (i1 > 1) {
      pred.m1c[[j]] <- DM[i2 + 1, i1] * draws$e[[j]] + (1 - DM[i2 + 1, i1]) * draws$c[[j]]
    }
    mu.m12[[j]] <- (1 - DM[j1, 1]) * fit$mu.m1.e[[j]] + DM[j1, 1] * fit$mu.m1.c[[j]]
    mu.m1.d00[[j]] <- DM[1, 1] * fit$mu.m1.e[[j]] + (1 - DM[1, 1]) * fit$mu.m1.c[[j]]
    st$pred.dat.y[, sp$m1[[j]]] <- pred.m1[[j]]
  }

  mu.m2 <- mu.m2o <- pred.m2 <- vector("list", n2)
  for (k in seq_len(n2)) {
    k1 <- k + 1
    if (expos[1] %in% sp$m2.eachcov[[k]]) {
      st$pred.dat.m2[[k]][, expos[1]] <- DM[1, k1] * ctx$exposed + (1 - DM[1, k1]) * ctx$control
      # M2 means for the exposure value of the first decomposition (d00)
      st$pred.dat.m2o[[k]][, expos[1]] <- DM[1, 1] * ctx$exposed + (1 - DM[1, 1]) * ctx$control
    }

    for (j in seq_len(n1)) {
      if (!(sp$m1[[j]] %in% sp$m2.eachcov[[k]])) next
      j1 <- j + 1
      m1j <- pred.m1[[j]]
      if (DM[j1, k1] == DM[j1, 1]) {
        newval <- m1j # not discrepant
      } else {
        # discrepant case: non-identifiable path, uses the cross-world correlation
        fam <- sp$family.m1[[j]]
        disp <- if (fam %in% c("gaussian", "negbin")) fit$disp.m1[[j]] else NULL
        newval <- gm_discrepant_draw(fam, m1j, mu.m1[[j]], mu.m12[[j]], disp,
                                     st$xcorr, nref)
      }
      st$pred.dat.m2[[k]][, sp$m1[[j]]] <- newval
      st$pred.dat.m2o[[k]][, sp$m1[[j]]] <- newval
    }

    mu.m2[[k]] <- predict(fit$m2[[k]], type = "response", newdata = st$pred.dat.m2[[k]])
    mu.m2o[[k]] <- predict(fit$m2[[k]], type = "response", newdata = st$pred.dat.m2o[[k]])
    fam2 <- sp$family.m2[[k]]
    pred.m2[[k]] <- gm_draw_m2(fam2, nref, mu.m2[[k]],
                               if (fam2 %in% c("gaussian", "negbin")) fit$disp.m2[[k]] else NULL)
  }

  st$pred.dat.y[, expos[1]] <- DM[1, 1] * ctx$exposed + (1 - DM[1, 1]) * ctx$control
  for (j in seq_len(n1)) st$pred.dat.y[, sp$m1[[j]]] <- pred.m1[[j]]
  for (k in seq_len(n2)) st$pred.dat.y[, sp$m2[[k]]] <- pred.m2[[k]]

  list(mu.m1 = mu.m1, pred.m1 = pred.m1, pred.m1c = pred.m1c, mu.m1.d00 = mu.m1.d00,
       mu.m2 = mu.m2, mu.m2o = mu.m2o, pred.m2 = pred.m2)
}

# ---- sensitivity corrections ---------------------------------------------

# Correction term for a given value `rho` of the sensitivity parameter.
#   stdv        standard deviation of the partner variable
#   fam_partner family of the partner variable (mediator in the assumption)
#   fam_y       family of the outcome
#   p           marginal proportion of the partner (binomial case)
#   resid_term  (m * (g - 1) + mu_a - g * mu_b)
#   binom_power power of the density ratio when both are binomial
#               (2 for Y vs M1, as in version 1; 1 otherwise)
gs_correction <- function(rho, stdv, fam_partner, fam_y, p, resid_term,
                          binom_power = 1) {
  if (fam_partner != "binomial") {
    (rho / stdv) * (1 / sqrt(1 - rho**2)) * resid_term
  } else if (fam_y != "binomial") {
    (rho / stdv) * (0.5 + (1 + (dnorm(qnorm(p)) / stdv)) / 4) * resid_term
  } else {
    (rho / stdv) * ((dnorm(qnorm(p)) / stdv)**binom_power) * resid_term
  }
}

# Draw M2 values when the assumption involves that M2 (i1 > 1): mean is
# shifted by the correction, then values are generated from a common uniform.
gs_corrected_m2 <- function(fam, mu, stdv, disp, crctm2, rho, u, nref) {
  if (fam == "gaussian") {
    return(qnorm(u, mean = mu + stdv * crctm2, sd = stdv))
  }
  if (fam == "poisson" || fam == "negbin") {
    j1 <- 1
    PMj <- rep(0, nref)
    chkcum <- rep(0, nref)
    while (any(u > PMj)) {
      j <- j1 - 1
      j1 <- j1 + 1
      PM0j <- if (fam == "poisson") ppois(j, mu) else pnbinom(j, mu = mu, size = disp)
      PMj <- pnorm(qnorm(PM0j) - crctm2 / sqrt(1 - rho**2))
      chkcum <- chkcum + (u > PMj)
    }
    return(chkcum)
  }
  if (fam == "binomial") {
    # P(M2 = 0) shifted by the correction; version 1 had an extra `1 -` here, which
    # drew M2 = 1 with probability 1 - mu instead of mu at rho = 0.
    PMj <- pnorm(qnorm(1 - mu) - crctm2 / sqrt(1 - rho**2))
    return(1 * (u > PMj))
  }
  NULL
}

# Expected value of Y given the (corrected) prediction data.
gs_outcome_mean <- function(muY0, crct, rho, fam, stdv.y, disp.y, ctx, j_start) {
  nref <- ctx$nref
  if (fam == "gaussian") return(muY0 + stdv.y * crct)
  if (fam == "binomial") return(1 - pnorm(qnorm(1 - muY0) - crct))
  if (fam == "poisson" || fam == "negbin") {
    j1 <- 0
    j <- j_start # version 1 started from the last loop index; see NEWS
    Pj <- matrix(0, 1, nref)
    sumPj <- matrix(0, 1, nref)
    PjT <- matrix(1 - ctx$tol, 1, nref) # tolerance should be close to 1
    chk <- 0
    while (chk < nref & j < ctx$maxcnt) {
      j <- j1 - 1
      j1 <- j1 + 1
      P0j <- if (fam == "poisson") ppois(j, muY0) else pnbinom(j, mu = muY0, size = disp.y)
      sumPj <- sumPj + Pj
      Pj <- pnorm(qnorm(P0j) - crct / sqrt(1 - rho**2))
      chk <- sum(Pj > PjT)
    }
    maxY <- j
    return((maxY * Pj) - sumPj)
  }
  NULL
}

# Expected potential outcome of Y for every value of the sensitivity parameter.
gs_expected_outcome <- function(fit, pp, DM, ctx, st, i1, i2) {
  sp <- ctx$spec
  n1 <- sp$n1
  nref <- ctx$nref
  g <- ctx$g
  i2a <- i2 - n1
  fam.y <- sp$family.y
  out <- numeric(ctx$nrcnt)
  rho_vec <- -1 + seq_len(ctx$nrcnt) * ctx$deltsi

  pred.m2u <- if (i1 > 1) runif(nref) else NULL
  muY0 <- predict(fit$y, type = "response", newdata = st$pred.dat.y)

  if (i1 > 1) {
    i1m2 <- i1 - 1
    mu.d0k <- DM[1, i1] * fit$mu.m1.e[[i2]] + (1 - DM[1, i1]) * fit$mu.m1.c[[i2]]
    mu.djk <- DM[i2 + 1, i1] * fit$mu.m1.e[[i2]] + (1 - DM[i2 + 1, i1]) * fit$mu.m1.c[[i2]]
    m2c <- pp$pred.m1c[[i2]]
    resid_m2 <- m2c * (g - 1) + mu.d0k - (g * mu.djk)
  }

  for (ri in seq_len(ctx$nrcnt)) {
    rhosi <- rho_vec[ri]
    crct <- 0
    pred.dat.y <- st$pred.dat.y

    if (i1 == 1) {
      if (i2 < n1 + 1) { # assumption involves a first-stage mediator
        crct <- gs_correction(
          rhosi, fit$stdv.m1[[i2]], sp$family.m1[[i2]], fam.y, fit$pm1[[i2]],
          pp$pred.m1[[i2]] * (g - 1) + pp$mu.m1.d00[[i2]] - (g * pp$mu.m1[[i2]]),
          binom_power = 2
        )
      } else { # assumption involves a second-stage mediator
        crct <- gs_correction(
          rhosi, fit$stdv.m2[[i2a]], sp$family.m2[[i2a]], fam.y, fit$pm2[[i2a]],
          pp$pred.m2[[i2a]] * (g - 1) + pp$mu.m2o[[i2a]] - (g * pp$mu.m2[[i2a]])
        )
      }
      muY0_r <- muY0
    } else {
      # assumption involves M2 as outcome: no direct correction of mean Y,
      # but the M2 values are corrected
      crctm2 <- gs_correction(rhosi, fit$stdv.m1[[i2]], sp$family.m1[[i2]], fam.y,
                              fit$pm1[[i2]], resid_m2)
      fam2 <- sp$family.m2[[i1m2]]
      pm2new <- gs_corrected_m2(fam2, pp$mu.m2[[i1m2]], gs_get(fit$stdv.m2, i1m2),
                                gs_get(fit$disp.m2, i1m2), crctm2, rhosi, pred.m2u, nref)
      pred.dat.y[, sp$m2[[i1m2]]] <- pm2new
      muY0_r <- predict(fit$y, type = "response", newdata = pred.dat.y)
    }

    out[ri] <- mean(gs_outcome_mean(muY0_r, crct, rhosi, fam.y, fit$stdv.y, fit$disp.y,
                                    ctx, j_start = n1))
  }
  out
}

# ---- scales ----------------------------------------------------------------

# Path effect from consecutive expected potential outcomes on the chosen scale.
gs_scale_effect <- function(e_hi, e_lo, scale) {
  if (scale == "MeanDiff" || scale == "") return(e_hi - e_lo)
  if (scale == "OR") {
    if (!all(e_hi > 0 & e_lo > 0 & e_hi < 1 & e_lo < 1)) {
      stop("Expected value for Y not in (0,1) - consider a different scale")
    }
    return(log((e_hi / (1 - e_hi)) / (e_lo / (1 - e_lo))))
  }
  if (scale == "RR") {
    if (!all(e_hi > 0 & e_lo > 0)) {
      stop("Expected value for Y is negative  - use scale = 'MeanDiff'")
    }
    return(log(e_hi / e_lo))
  }
  rep(0, length(e_hi))
}

# ---- one bootstrap sample ----------------------------------------------------

gs_boot_iteration <- function(b, fit, ctx, acc, st, i1, i2) {
  sp <- ctx$spec
  n1 <- sp$n1
  n2 <- sp$n2
  npath <- ctx$npath
  nxcorr <- ctx$nxcorr

  draws <- gs_draw_first_stage(fit, ctx)
  st$xcorr <- c()

  DM <- matrix(1, n1 + 1, n2 + 1)
  PrevDM <- NULL
  nullv <- rep(0, npath)
  xcrep <- rep(1, npath + 1)
  m1d <- 0

  for (c1 in seq_len(npath + 1)) {
    c <- c1 - 1 # path index; c1 = 1 is the all-exposed potential outcome
    m1chk <- "Null"
    m2chk <- "Null"

    if (c1 > 1) {
      idx <- ctx$indexDM[, c1 - 1]
      PrevDM <- DM
      DM[idx[1], idx[2]] <- 0
      m1d <- idx[1] - 1
      m2d <- idx[2] - 1
      m1chk <- if (m1d > 0) sp$m1[[m1d]] else ctx$expos[1]
      m2chk <- if (m2d > 0) sp$m2.eachcov[[m2d]] else sp$y.allcov
      xcrep[c1] <- if (gs_is_discrepant(DM, ctx)) nxcorr else 1
      if (b == 1) acc$pathpr[c] <- gs_path_shown(DM, PrevDM, i1, i2, n1)
    }

    for (jd in seq_len(xcrep[c1])) {
      if (c1 == 1 || m1chk %in% m2chk) {
        if (m1d > 0) st$xcorr <- ctx$xcorrv[jd]
        pp <- gs_path_predictions(fit, draws, DM, jd, ctx, st, i1, i2)
        acc$EYDM[[c1]][b, , jd] <- gs_expected_outcome(fit, pp, DM, ctx, st, i1, i2)
      } else { # null path: same expected outcome as the previous one
        acc$EYDM[[c1]][b, , jd] <- acc$EYDM[[c1 - 1]][b, , jd]
        nullv[c] <- 1
        xcrep[c1] <- 1
      }
    }

    # extend arrays for paths not affected by the cross-world correlation
    if (xcrep[c1] < nxcorr) {
      for (lt in 2:nxcorr) {
        acc$EYDM[[c1]][b, , lt] <-
          if (c1 == 1 || nullv[c] == 0) acc$EYDM[[c1]][b, , 1] else acc$EYDM[[c1 - 1]][b, , lt]
      }
    }

    if (c1 > 1) {
      for (jd in seq_len(nxcorr)) {
        eff <- gs_scale_effect(acc$EYDM[[c1 - 1]][b, , jd], acc$EYDM[[c1]][b, , jd], ctx$scale)
        acc$EDM_SC[[c]][b, , jd] <- eff
        acc$EDMM2[c, , b, jd] <- eff
      }
      acc$EDMM1[c, , b] <- acc$EDM_SC[[c]][b, , 1]
    }
  }
  acc$nullv <- nullv
  acc$xcrep <- xcrep

  # total effect through each mediator: sum of the path effects that change it
  idx1 <- ctx$indexDM[1, ]
  idx2 <- ctx$indexDM[2, ]
  for (j in seq_len(n1)) {
    acc$EDM1_SC[[j]][b, ] <- vapply(seq_len(ctx$nrcnt), function(r) {
      sum(acc$EDMM1[, r, b] * (idx1 == j + 1), na.rm = TRUE)
    }, numeric(1))
  }
  for (k in seq_len(n2)) {
    for (jd in seq_len(nxcorr)) {
      acc$EDM2_SC[[k]][b, , jd] <- vapply(seq_len(ctx$nrcnt), function(r) {
        sum(acc$EDMM2[, r, b, jd] * (idx2 == k + 1), na.rm = TRUE)
      }, numeric(1))
    }
  }
  invisible(acc)
}

# Run the bootstrap for all assumption pairs. The bootstrap samples are drawn
# once and the models are refit once per sample, so every pair is evaluated on
# the same resamples. Returns one list of result arrays per pair.
gs_run_pairs <- function(pairs, bootsamples, ctx, st, boots) {
  accs <- lapply(seq_len(nrow(pairs)), function(p) gs_new_arrays(ctx, boots))
  for (b in seq_len(boots + 1)) {
    fit <- gs_fit_boot_models(bootsamples[[b]], ctx)
    for (p in seq_len(nrow(pairs))) {
      gs_boot_iteration(b, fit, ctx, accs[[p]], st, pairs$i1[p], pairs$i2[p])
    }
  }
  accs
}
