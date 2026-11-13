# Internal engine of gmediate(): one bootstrap iteration (model refits, draws
# of the potential mediators, expected potential outcomes) and the assembly
# of effects across iterations.

#' Refit all models on one bootstrap sample and compute reference-group means
#' @noRd
gm_fit_boot_models <- function(bdata, ctx) {
  spec <- ctx$spec
  fit <- list(m1 = list(), m2 = list(), disp.m1 = list(), disp.m2 = list(),
              mu.m1.e = list(), mu.m1.c = list())

  for (k in seq_len(spec$n2)) {
    fit$m2[[k]] <- gm_refit(spec$call.m2[[k]], spec$family.m2[[k]],
                            bdata[, c(spec$m2[[k]], spec$m2.eachcov[[k]])], ctx$envir)
    fit$disp.m2[[k]] <- gm_dispersion(fit$m2[[k]], spec$family.m2[[k]])
  }
  for (j in seq_len(spec$n1)) {
    fit$m1[[j]] <- gm_refit(spec$call.m1[[j]], spec$family.m1[[j]],
                            bdata[, c(spec$m1[[j]], spec$m1.eachcov[[j]])], ctx$envir)
    fit$disp.m1[[j]] <- gm_dispersion(fit$m1[[j]], spec$family.m1[[j]])
    fit$mu.m1.e[[j]] <- predict(fit$m1[[j]], type = "response", newdata = ctx$pred$pred.dat.m1.e[[j]])
    fit$mu.m1.c[[j]] <- predict(fit$m1[[j]], type = "response", newdata = ctx$pred$pred.dat.m1.c[[j]])
  }
  fit$y <- gm_refit(spec$call.y, spec$family.y,
                    bdata[, c(spec$y, spec$y.allcov)], ctx$envir)
  fit
}

#' One bootstrap iteration: expected potential outcomes for every path
#'
#' @param state environment holding `nonnullm1`, which (as in version 1)
#'   accumulates across iterations.
#' @return list(EY, nullv, IDENT, rhrepd, fit).
#' @noRd
gm_boot_iteration <- function(bdata, ctx, state) {
  spec <- ctx$spec
  n1 <- spec$n1; n2 <- spec$n2; npath <- spec$npath; nrho <- ctx$nrho
  m1 <- spec$m1; m2 <- spec$m2; expos <- ctx$expos; nref <- ctx$nref
  exposed <- ctx$pred$exposed; control <- ctx$pred$control
  indexDM <- ctx$indexDM

  fit <- gm_fit_boot_models(bdata, ctx)
  pred.dat.y <- ctx$pred$pred.dat.y
  pred.dat.m2 <- ctx$pred$pred.dat.m2

  # factual / counterfactual first-stage mediator draws (everyone exposed,
  # everyone control)
  pred.m1.e <- list(); pred.m1.c <- list()
  for (j in seq_len(n1)) {
    pred.m1.e[[j]] <- gm_draw_m1(spec$family.m1[[j]], nref, fit$mu.m1.e[[j]], fit$disp.m1[[j]])
    pred.m1.c[[j]] <- gm_draw_m1(spec$family.m1[[j]], nref, fit$mu.m1.c[[j]], fit$disp.m1[[j]])
  }

  DM <- matrix(1, n1 + 1, n2 + 1)
  nullv <- matrix(0, 1, npath)
  DISC <- matrix(0, 1, npath + 1)
  EY <- vector("list", npath + 1)
  IDENT <- list(); rhrepd <- list(); rhrep <- list()
  mu.m1 <- list(); pred.m1 <- list(); mu.m12 <- list(); pred.m12 <- list()
  mu.m2 <- list(); pred.m2 <- list()
  rho <- NULL
  rhrep[[1]] <- 1
  DISC[1, 1] <- 0

  for (c1 in 1:(npath + 1)) {
    c <- c1 - 1
    m1chk <- "Null"
    m2chk <- "Null"
    if (c1 == 1) m1d <- 0
    if (c1 > 1) {
      indexDMV <- matrix(indexDM[, c1 - 1], 2, 1)
      DM[indexDMV[1], indexDMV[2]] <- 0
      m1d <- indexDMV[1] - 1
      m2d <- indexDMV[2] - 1
      m1chk <- if (m1d > 0) m1[[m1d]] else expos[1]
      m2chk <- if (m2d > 0) spec$m2.eachcov[[m2d]] else spec$y.allcov

      # a path is "discrepant" when an M1 is switched on a link into an M2
      for (jj in 2:(n1 + 1)) {
        for (kk in 2:(n2 + 1)) {
          if (m1[[jj - 1]] %in% spec$m2.eachcov[[kk - 1]] & DM[jj, kk] != DM[jj, 1]) DISC[1, c1] <- 1
        }
      }
      rhrep[[c1]] <- if (DISC[1, c1] == 1) nrho else 1
    }

    EY[[c1]] <- numeric(nrho)
    for (jd in 1:rhrep[[c1]]) {
      if (c1 == 1 | m1chk %in% m2chk) {
        if (m1d > 0) {
          state$nonnullm1[[m1d]][1, m2d + 1] <- 1
          rho <- ctx$xcorrv[jd]
        }

        for (j in seq_len(n1)) {
          j1 <- j + 1
          mu.m1[[j]] <- DM[j1, 1] * fit$mu.m1.e[[j]] + (1 - DM[j1, 1]) * fit$mu.m1.c[[j]]
          pred.m1[[j]] <- DM[j1, 1] * pred.m1.e[[j]] + (1 - DM[j1, 1]) * pred.m1.c[[j]]
          mu.m12[[j]] <- (1 - DM[j1, 1]) * fit$mu.m1.e[[j]] + DM[j1, 1] * fit$mu.m1.c[[j]]
          pred.m12[[j]] <- (1 - DM[j1, 1]) * pred.m1.e[[j]] + DM[j1, 1] * pred.m1.c[[j]]
        }

        for (k in seq_len(n2)) {
          k1 <- k + 1
          if (expos[1] %in% spec$m2.eachcov[[k]]) {
            pred.dat.m2[[k]][, expos[1]] <- DM[1, k1] * exposed + (1 - DM[1, k1]) * control
          }
          for (j in seq_len(n1)) {
            if (m1[[j]] %in% spec$m2.eachcov[[k]]) {
              j1 <- j + 1
              m1j <- pred.m1[[j]]
              if (DM[j1, k1] == DM[j1, 1]) {
                pred.dat.m2[[k]][, m1[[j]]] <- m1j
              } else {
                pred.m12[[j]] <- gm_discrepant_draw(spec$family.m1[[j]], m1j, mu.m1[[j]],
                                                    mu.m12[[j]], fit$disp.m1[[j]], rho, nref)
                pred.dat.m2[[k]][, m1[[j]]] <- pred.m12[[j]]
              }
            }
          }
          mu.m2[[k]] <- predict(fit$m2[[k]], type = "response", newdata = pred.dat.m2[[k]])
          pred.m2[[k]] <- gm_draw_m2(spec$family.m2[[k]], nref, mu.m2[[k]], fit$disp.m2[[k]])
        }

        pred.dat.y[, expos[1]] <- DM[1, 1] * exposed + (1 - DM[1, 1]) * control
        for (j in seq_len(n1)) pred.dat.y[, m1[[j]]] <- pred.m1[[j]]
        for (k in seq_len(n2)) pred.dat.y[, m2[[k]]] <- pred.m2[[k]]
        muY <- predict(fit$y, type = "response", newdata = pred.dat.y)
        EY[[c1]][jd] <- mean(muY)
      } else {
        EY[[c1]][1] <- EY[[c1 - 1]][1]
        nullv[1, c] <- 1
        rhrep[[c1]] <- 1
      }
    }

    # paths unaffected by the correlation: repeat the single estimate
    if (rhrep[[c1]] < nrho) {
      for (lt in 2:nrho) {
        if (c1 == 1) {
          EY[[c1]][lt] <- EY[[c1]][1]
        } else if (nullv[1, c] == 0) {
          EY[[c1]][lt] <- EY[[c1]][1]
        } else {
          EY[[c1]][lt] <- EY[[c1 - 1]][lt]
        }
      }
    }

    if (c1 > 1) {
      IDENT[[c]] <- 0
      if (m1d == 0) {
        IDENT[[c]] <- 1
      } else if (m1d > 0) {
        if (nullv[1, c] == 1 | sum(state$nonnullm1[[m1d]]) <= 1) IDENT[[c]] <- 1
      }
      rhrepd[[c]] <- if (IDENT[[c]] == 1) 1 else nrho
    }
  }
  list(EY = EY, nullv = nullv, IDENT = IDENT, rhrepd = rhrepd, fit = fit)
}

#' Run all bootstrap iterations and assemble the effect arrays
#'
#' @return list with the path-effect matrices `EDM` (one `(boots+1) x nrho`
#'   matrix per path), total effects `EDT`, and the first- / second-stage
#'   mediator effect sums `EDB1`, `EDB2`, plus the last iteration's `nullv`,
#'   `IDENT`, `rhrepd` and model fits.
#' @noRd
gm_run_bootstrap <- function(bootsamples, ctx) {
  spec <- ctx$spec
  npath <- spec$npath; nrho <- ctx$nrho; B <- length(bootsamples); boots <- B - 1
  state <- new.env(parent = emptyenv())
  state$nonnullm1 <- lapply(seq_len(spec$n1), function(i) matrix(0, 1, spec$n2 + 1))

  EYDM <- lapply(seq_len(npath + 1), function(p) matrix(0, B, nrho))
  last <- NULL
  for (b in seq_len(B)) {
    last <- gm_boot_iteration(bootsamples[[b]], ctx, state)
    for (c1 in seq_len(npath + 1)) EYDM[[c1]][b, ] <- last$EY[[c1]]
  }

  scale <- ctx$scale
  EDM <- lapply(seq_len(npath), function(c) {
    gm_path_effect(EYDM[[c]], EYDM[[c + 1]], scale)
  })
  EDT <- matrix(0, B, 1)
  EDT[, 1] <- gm_total_effect(EYDM[[1]][, 1], EYDM[[npath + 1]][, 1], scale)

  # effects of each mediator = sums of the path effects that switch it
  EDMM1 <- matrix(0, npath, B)
  EDMM2 <- lapply(seq_len(nrho), function(jd) matrix(0, npath, B))
  for (c in seq_len(npath)) {
    EDMM1[c, ] <- EDM[[c]][, 1]
    for (jd in seq_len(nrho)) EDMM2[[jd]][c, ] <- EDM[[c]][, jd]
  }
  exp_scale <- scale == "OR" || scale == "RR"
  EDB1 <- lapply(seq_len(spec$n1), function(j) {
    distinguish <- sapply(1:ncol(EDMM1), function(x) EDMM1[, x] * (ctx$indexDM[1, ] == j + 1))
    s <- sapply(1:ncol(distinguish), function(x) sum(distinguish[, x], na.rm = TRUE))
    if (exp_scale) exp(s) else s
  })
  EDB2 <- lapply(seq_len(spec$n2), function(k) {
    m <- matrix(0, nrho, B)
    for (jd in seq_len(nrho)) {
      distinguish <- sapply(1:ncol(EDMM2[[jd]]), function(x) EDMM2[[jd]][, x] * (ctx$indexDM[2, ] == k + 1))
      m[jd, ] <- sapply(1:ncol(distinguish), function(x) sum(distinguish[, x], na.rm = TRUE))
      if (exp_scale) m[jd, ] <- exp(m[jd, ])
    }
    m
  })
  list(EDM = EDM, EDT = EDT, EDB1 = EDB1, EDB2 = EDB2, last = last)
}
