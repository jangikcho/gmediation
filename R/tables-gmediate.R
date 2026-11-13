# Internal helpers: turning bootstrap effect arrays into the seven summary
# tables returned (and printed) by gmediate().

#' Column names used in all summary tables
#' @noRd
gm_ci_names <- function(CIPCT) {
  c(paste(CIPCT, "% CI Lower", sep = ""), paste(CIPCT, "% CI Upper", sep = ""))
}

#' Summarise one path's bootstrap matrix over its correlation columns
#'
#' @param boot_mat `(boots+1) x ncol` matrix of effects (last row = original).
#' @param ncols number of columns (correlation values) that matter.
#' @param ratio_by optional vector to divide each bootstrap row by (used for
#'   proportions); the original-sample estimate is passed in `est_orig`.
#' @return list(est, lower, upper, p) vectors of length `ncols`.
#' @noRd
gm_summarise_columns <- function(boot_mat, ncols, boots, probs, ratio_by = NULL) {
  out <- list(est = numeric(ncols), lower = numeric(ncols), upper = numeric(ncols),
              p = character(ncols))
  for (jd in seq_len(ncols)) {
    x <- boot_mat[1:boots, jd]
    if (!is.null(ratio_by)) x <- x / ratio_by[1:boots]
    s <- gm_boot_summary(x, probs[1], probs[2], na.rm = TRUE)
    out$lower[jd] <- s$lower
    out$upper[jd] <- s$upper
    out$p[jd] <- s$p
  }
  out
}

#' Table of path-specific effects or of their proportions of the total effect
#'
#' @param proportion if TRUE build the proportion table (effect / total).
#' @return data.frame with columns Path, Xcorr, Estimate (or Est Prop), CI
#'   limits and p-value; null paths are marked "*".
#' @noRd
gm_path_table <- function(res, ctx, boots, probs, ESTT, proportion = FALSE) {
  spec <- ctx$spec
  npath <- spec$npath
  last <- res$last
  nullv <- last$nullv
  m1ex <- c("   ", spec$m1)
  m2ex <- c("   ", spec$m2)

  e1 <- e2 <- rho <- character(0)
  est <- lb <- ub <- pv <- list()
  for (c in seq_len(npath)) {
    ncols <- last$rhrepd[[c]]
    if (nullv[c] == 1) {
      cells <- list(est = rep("*", ncols), lower = rep("*", ncols),
                    upper = rep("*", ncols), p = rep("*", ncols))
    } else {
      cells <- gm_summarise_columns(res$EDM[[c]], ncols, boots, probs,
                                    ratio_by = if (proportion) res$EDT[, 1] else NULL)
      cells$est <- round(res$EDM[[c]][boots + 1, seq_len(ncols)], 7)
      if (proportion) cells$est <- round(cells$est / ESTT, 4)
    }
    for (jd in seq_len(ncols)) {
      e1[length(e1) + 1] <- paste(m1ex[ctx$indexDM[1, c]])
      e2[length(e2) + 1] <- paste(m2ex[ctx$indexDM[2, c]])
      rho[length(rho) + 1] <-
        if (nullv[c] == 0 & last$IDENT[[c]] == 1 || nullv[c] == 1) "-" else paste(ctx$xcorrv[jd])
    }
    est[[c]] <- cells$est; lb[[c]] <- cells$lower
    ub[[c]] <- cells$upper; pv[[c]] <- cells$p
  }

  Path <- mapply(function(a, b) {
    parts <- c(ctx$expos[1], a, b, spec$y)
    parts <- parts[!is.na(parts) & nzchar(trimws(parts))]
    parts <- substr(parts, 1, 7)
    paste(parts, collapse = " -> ")
  }, e1, e2, USE.NAMES = FALSE)

  tab <- data.frame(Path = Path, Xcorr = rho,
                    Estimate = unlist(est), `CI Lower` = unlist(lb),
                    `CI Upper` = unlist(ub), `p-value` = unlist(pv),
                    check.names = FALSE)
  colnames(tab) <- c("Path", "Xcorr", if (proportion) "Estimate" else "Estimate",
                     gm_ci_names(ctx$CIPCT), "p-value")
  tab
}

#' Total-effect table (one row named "Total")
#' @noRd
gm_total_table <- function(EDT, boots, probs, CIPCT) {
  ESTT <- round(EDT[boots + 1, ], 7)
  s <- gm_boot_summary(EDT[1:boots, ], probs[1], probs[2], na.rm = TRUE)
  tab <- data.frame(ESTT, s$lower, s$upper, s$p)
  colnames(tab) <- c("Estimate", gm_ci_names(CIPCT), "p-value")
  rownames(tab) <- "Total"
  tab
}

#' First-stage mediator table (effects or proportions)
#' @noRd
gm_first_stage_table <- function(res, ctx, boots, probs, ESTT, proportion = FALSE) {
  n1 <- ctx$spec$n1
  est <- lb <- ub <- pv <- numeric(0)
  for (j in seq_len(n1)) {
    e <- round(res$EDB1[[j]][boots + 1], 7)
    x <- res$EDB1[[j]][1:boots]
    if (proportion) {
      e <- round(e / ESTT, 4)
      x <- x / res$EDT[1:boots, 1]
    }
    s <- gm_boot_summary(x, probs[1], probs[2], na.rm = TRUE)
    est[j] <- e; lb[j] <- s$lower; ub[j] <- s$upper; pv[j] <- s$p
  }
  tab <- data.frame(unlist(ctx$spec$m1), est, lb, ub, pv)
  colnames(tab) <- c("Mediator", if (proportion) "Est Prop" else "Estimate",
                     gm_ci_names(ctx$CIPCT), "p-value")
  tab
}

#' Second-stage mediator table (effects or proportions)
#' @noRd
gm_second_stage_table <- function(res, ctx, boots, probs, ESTT, proportion = FALSE) {
  n2 <- ctx$spec$n2
  nrho <- ctx$nrho
  names_ <- est <- lb <- ub <- pv <- list()
  for (k in seq_len(n2)) {
    for (jd in seq_len(nrho)) {
      e <- round(res$EDB2[[k]][jd, boots + 1], 7)
      x <- res$EDB2[[k]][jd, 1:boots]
      if (proportion) {
        e <- round(e / ESTT, 4)
        x <- x / res$EDT[1:boots, 1]
      }
      s <- gm_boot_summary(x, probs[1], probs[2], na.rm = proportion)
      names_[[length(names_) + 1]] <- paste(ctx$spec$m2[[k]])
      est[[length(est) + 1]] <- e; lb[[length(lb) + 1]] <- s$lower
      ub[[length(ub) + 1]] <- s$upper; pv[[length(pv) + 1]] <- s$p
    }
  }
  rho <- sprintf("%.1f", rep(ctx$xcorrv, n2))
  tab <- data.frame(unlist(names_), rho, unlist(est), unlist(lb), unlist(ub), unlist(pv))
  colnames(tab) <- c("Mediator", "Xcorr", if (proportion) "Est Prop" else "Estimate",
                     gm_ci_names(ctx$CIPCT), "p-value")
  tab
}
