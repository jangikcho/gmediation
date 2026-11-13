# Internal helpers: turn the result arrays of one assumption pair into
# sensitivity curves (used for plots) and summary tables (written to CSV).

# Estimate and percentile confidence limits over the grid of rho values.
# `boot_mat` has one row per bootstrap sample (the last row is the original
# sample) and one column per value of rho.
gs_curve_data <- function(boot_mat, probs, boots, rho_plot, transform = identity) {
  boot_mat <- matrix(boot_mat, nrow = boots + 1)
  q <- apply(boot_mat, 2, quantile, probs)
  data.frame(
    rho = rho_plot,
    estimate = transform(boot_mat[boots + 1, ]),
    lower = transform(q[1, ]),
    upper = transform(q[2, ])
  )
}

# Table rows for one curve, rounded as in version 1.
gs_curve_rows <- function(curve, rho_tab, ndig, ...) {
  # in a non-UTF-8 locale the column name for rho cannot be shown natively
  suppressWarnings(data.frame(
    ...,
    "\u03c1" = rho_tab,
    Estimate = round(curve$estimate, ndig),
    `Lower CI` = round(curve$lower, ndig),
    `Upper CI` = round(curve$upper, ndig)
  ))
}

gs_curve_record <- function(kind, title, ylab, data) {
  list(kind = kind, title = title, ylab = ylab, data = data)
}

# Names of the pair of variables for which no unmeasured confounding is assumed.
gs_pair_names <- function(spec, i1, i2) {
  av1 <- if (i1 == 1) spec$y else spec$m2[[i1 - 1]]
  av2 <- if (i2 <= spec$n1) spec$m1[[i2]] else spec$m2[[i2 - spec$n1]]
  c(av1, av2)
}

# Sensitivity curves and tables of all paths and mediators for one pair.
gs_pair_outputs <- function(acc, ctx, boots, probs, ndig, assumvars, i1, i2) {
  sp <- ctx$spec
  nxcorr <- ctx$nxcorr
  rho_plot <- seq(from = -1 + ctx$deltsi, to = 1 - ctx$deltsi, by = ctx$deltsi)
  rho_tab <- round(-1 + seq_len(ctx$nrcnt) * ctx$deltsi, 2)
  scale_back <- if (ctx$scale %in% c("OR", "RR")) exp else identity
  all_assum <- "ALL" %in% toupper(assumvars)
  curves <- list()
  path_tabs <- list()
  first_tabs <- list()
  second_tabs <- list()

  # -- path-specific effects: non-null paths that the assumption can affect
  for (pathnum in seq_len(ctx$npath)) {
    if (acc$nullv[pathnum] != 0 || acc$pathpr[pathnum] != 1) next
    pathv <- ctx$indexDM[, pathnum]
    m1name <- if (pathv[1] > 1) trimws(sp$m1[[pathv[1] - 1]]) else ""
    m2name <- if (pathv[2] > 1) trimws(sp$m2[[pathv[2] - 1]]) else ""
    nxcp <- if (acc$xcrep[pathnum + 1] == 1 && acc$xcrep[pathnum] == 1) 1 else nxcorr
    mids <- c(ctx$expos[1], m1name[nzchar(m1name)], m2name[nzchar(m2name)])

    for (xci in seq_len(nxcp)) {
      last <- if (nzchar(m1name)) sprintf("%s, xcorr = %s", sp$y, ctx$xcorrv[[xci]]) else sp$y
      title <- paste("Path:", paste(c(mids, last), collapse = " \u2192 "))
      curve <- gs_curve_data(acc$EDM_SC[[pathnum]][, , xci], probs, boots, rho_plot, scale_back)
      curves[[length(curves) + 1]] <- gs_curve_record("path", title, "Estimated PSE", curve)

      rows <- gs_curve_rows(curve, rho_tab, ndig, Path = title)
      if (any(complete.cases(rows))) {
        xcorr <- if (nzchar(m1name)) as.numeric(as.character(ctx$xcorrv[[xci]])) else NA_real_
        rows$Path <- paste(c(mids, sp$y), collapse = " \u2192 ")
        path_tabs[[length(path_tabs) + 1]] <-
          data.frame(rows[1], Xcorr = xcorr, rows[-1], check.names = FALSE)
      }
    }
  }

  # -- total effect through each first-stage mediator
  for (j in seq_len(sp$n1)) {
    if (!all_assum && !(sp$m1[[j]] %in% assumvars)) next
    curve <- gs_curve_data(acc$EDM1_SC[[j]], probs, boots, rho_plot)
    curves[[length(curves) + 1]] <-
      gs_curve_record("first", paste0("Mediator = ", sp$m1[[j]]), "Estimated Mediator Effect", curve)
    first_tabs[[length(first_tabs) + 1]] <-
      gs_curve_rows(curve, rho_tab, ndig, Mediator = sp$m1[[j]])
  }

  # -- total effect through each second-stage mediator
  for (k in seq_len(sp$n2)) {
    if (!all_assum && !(sp$m2[[k]] %in% assumvars)) next
    for (xci in seq_len(nxcorr)) {
      curve <- gs_curve_data(acc$EDM2_SC[[k]][, , xci], probs, boots, rho_plot)
      title <- paste0("Mediator = ", sp$m2[[k]], ", xcorr = ", ctx$xcorrv[[xci]])
      curves[[length(curves) + 1]] <-
        gs_curve_record("second", title, "Estimated Mediator Effect", curve)
      rows <- gs_curve_rows(
        curve, rho_tab, ndig, Mediator = rep(sp$m2[[k]], nrow(curve)),
        Xcorr = rep(as.numeric(as.character(ctx$xcorrv[[xci]])), nrow(curve))
      )
      if (any(complete.cases(rows))) second_tabs[[length(second_tabs) + 1]] <- rows
    }
  }

  list(curves = curves,
       path = gs_bind_tables(path_tabs),
       first = gs_bind_tables(first_tabs),
       second = gs_bind_tables(second_tabs),
       assum = gs_pair_names(sp, i1, i2))
}

gs_bind_tables <- function(tabs) {
  if (length(tabs) == 0) return(NULL)
  out <- do.call(rbind, tabs)
  rownames(out) <- NULL
  out
}

# Combine the per-pair results into the tables and the long curve data frame
# stored in the "gmedsens" object.
gs_combine_pairs <- function(pairs) {
  tag <- function(tab, label) {
    if (is.data.frame(tab) && nrow(tab) > 0) {
      tab$Rho.Between <- label
      tab
    }
  }
  collect <- function(which) {
    gs_bind_tables(Filter(Negate(is.null), lapply(pairs, function(p) {
      tag(p[[which]], paste0(p$assum[1], " & ", p$assum[2]))
    })))
  }

  curve_rows <- list()
  id <- 0
  for (p in pairs) {
    for (cv in p$curves) {
      id <- id + 1
      curve_rows[[id]] <- data.frame(
        curve_id = id, assum1 = p$assum[1], assum2 = p$assum[2],
        kind = cv$kind, title = cv$title, ylab = cv$ylab, cv$data,
        stringsAsFactors = FALSE
      )
    }
  }
  list(path = collect("path"), first = collect("first"), second = collect("second"),
       curves = gs_bind_tables(curve_rows))
}
