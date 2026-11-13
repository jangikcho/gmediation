# Internal helpers: writing results to disk and summarising model fits.

#' Write the gmediate summary tables to one Excel workbook
#'
#' @param tables named list of data frames.
#' @param table_out length-2 character vector: "true"/"false", then folder.
#' @noRd
gm_write_tables_xlsx <- function(tables, table_out) {
  if (table_out[1] != "true") return(invisible(NULL))
  if (nchar(table_out[2]) == 0) {
    warning("Please specify a folder in table_out[2] to save the output files.",
            call. = FALSE)
    return(invisible(NULL))
  }
  sheets <- list(
    "Individual Path" = tables$indiv.path,
    "Total Effect" = tables$total.path,
    "First Stage" = tables$first.path,
    "Second Stage" = tables$second.path,
    "Prop Individual Path" = tables$prop.indiv.path,
    "Prop First Stage" = tables$prop.first.path,
    "Prop Second Stage" = tables$prop.second.path
  )
  if (!dir.exists(table_out[2])) dir.create(table_out[2], recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path = file.path(table_out[2], "gmediate_summary_tables.xlsx"))
  invisible(NULL)
}

# Formula of a fitted model; glm.nb() fits do not store it, so fall back to
# the formula of the original call.
gm_model_formula <- function(fit, call) {
  f <- fit$formula
  if (is.null(f)) call$formula else f
}

#' Formula and coefficient table of every refit model
#'
#' @param fit list returned by gm_fit_boot_models() for the original sample.
#' @return list with one element per model (M1s, then M2s, then Y).
#' @noRd
gm_model_fit_summary <- function(fit, spec) {
  n1 <- spec$n1
  n2 <- spec$n2
  mfit <- vector("list", n1 + n2 + 1)
  for (j in seq_len(n1)) {
    mfit[[j]][[1]] <- gm_model_formula(fit$m1[[j]], spec$call.m1[[j]])
    mfit[[j]][[2]] <- summary(fit$m1[[j]])$coefficient
  }
  for (k in seq_len(n2)) {
    mfit[[n1 + k]][[1]] <- gm_model_formula(fit$m2[[k]], spec$call.m2[[k]])
    mfit[[n1 + k]][[2]] <- summary(fit$m2[[k]])$coefficient
  }
  mfit[[n1 + n2 + 1]][[1]] <- gm_model_formula(fit$y, spec$call.y)
  mfit[[n1 + n2 + 1]][[2]] <- summary(fit$y)$coefficient
  mfit
}

# ---- gmedsens exports --------------------------------------------------------

# Output folder requested through plot_out / table_out (c("true", folder)):
# the folder (created if needed), or NULL when no output is requested or no
# folder was given.
gs_output_dir <- function(out, arg, what) {
  if (!identical(out[1], "true")) return(NULL)
  if (nchar(out[2]) == 0) {
    warning("Please specify a folder in ", arg, " ", what, call. = FALSE)
    return(NULL)
  }
  if (!dir.exists(out[2])) dir.create(out[2], recursive = TRUE, showWarnings = FALSE)
  out[2]
}

# Write the sensitivity tables as CSV files; returns the paths written.
gs_write_csv <- function(out, dir) {
  targets <- c(path = "gmedsens_pathspec_out.csv", first = "gmedsens_firststage_out.csv",
               second = "gmedsens_secondstage_out.csv")
  written <- character()
  for (nm in names(targets)) {
    tab <- out[[nm]]
    if (!is.null(tab) && nrow(tab) > 0) {
      f <- file.path(dir, targets[[nm]])
      write.csv(tab, f, row.names = FALSE)
      written <- c(written, f)
    }
  }
  written
}

# Pairs (i1, i2) of variables for which "no unmeasured confounding" is assumed.
gs_assumption_pairs <- function(assumvars, spec) {
  n1 <- spec$n1
  n2 <- spec$n2
  if (all(assumvars == "ALL")) {
    pairs <- do.call(rbind, lapply(seq_len(n2 + 1), function(i1) {
      i2end <- if (i1 == 1) n1 + n2 else n1
      data.frame(i1 = i1, i2 = seq_len(i2end))
    }))
    return(pairs)
  }
  if (length(assumvars) != 2 || any(assumvars == "ALL")) {
    stop("assumvars must be \"ALL\" or the names of two variables")
  }
  m1 <- unlist(spec$m1)
  m2 <- unlist(spec$m2)
  if (assumvars[1] == spec$y) {
    i1 <- 1
    i2 <- c(match(assumvars[2], m1), n1 + match(assumvars[2], m2))
  } else {
    i1 <- 1 + match(assumvars[1], m2)
    i2 <- match(assumvars[2], m1)
  }
  i2 <- i2[!is.na(i2)]
  if (is.na(i1) || length(i2) == 0) {
    stop("assumvars must name the outcome or a second-stage mediator, and a mediator")
  }
  data.frame(i1 = i1, i2 = i2[length(i2)])
}
