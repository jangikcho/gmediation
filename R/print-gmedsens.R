#' Print a gmedsens object
#'
#' Reproduces the report of [gmedsens()]: for every page of the plots PDF the
#' sensitivity plots it holds, followed by a summary of the analysis settings
#' and the location of the saved files. The sensitivity tables are in the
#' elements `path`, `first.stage` and `second.stage` of the object.
#'
#' @param x an object of class `"gmedsens"`, as returned by [gmedsens()].
#' @param ... unused.
#' @return `x`, invisibly.
#' @seealso [gmedsens()], [plot.gmedsens()]
#' @export
print.gmedsens <- function(x, ...) {
  s <- x$settings
  if (!is.null(x$files$plots)) {
    for (pg in gs_pages(x$curves, s$n_rows, s$n_cols)) {
      cat("\n----------------------------------------\n")
      cat("PDF Page ", pg$page, "\n", sep = "")
      cat(gs_page_heading(pg$kind), "\n",
          "Assumption: No Unmeasured Confounding between '", pg$assum1, "' and '",
          pg$assum2, "'", " \n", sep = "")
      for (tt in pg$titles) cat("  ", tt, "\n", sep = "")
      cat("----------------------------------------\n")
    }
  }

  cat("----------------------------------------\n")
  cat("Exposure ('exposed' value):", paste0(s$expos[1], "(", s$expos[2], ")", "\n"))
  cat("Final Outcome:", s$y, "\n")
  cat("Sample Size Used:", s$n, "\n")
  cat("Number of Bootstrap Samples / Number Used:", sprintf("%d / %d\n", s$boots, s$boots))
  if (!identical(s$ref, "entire sample")) {
    cat("Reference Group:", s$ref[1], "= (", paste(s$ref[-1], collapse = ", "), "), ",
        "Reference Group Multiplier (refmult):", s$refmult, "\n")
  } else {
    cat("Reference Group: Whole Sample, refmult:", s$refmult, "\n")
  }
  cat("Scale:", s$scale, "\n")
  if (!is.null(s$cluster)) cat("Cluster:", s$cluster, "\n")
  cat("Confidence Interval Percentage:", s$CIPCT, "\n")
  cat("g (between-correlation factor):", s$g, "\n")
  cat("----------------------------------------\n")

  if (length(x$files$tables) > 0) {
    cat("Tables saved to:\n")
    for (f in x$files$tables) cat("  ", f, "\n", sep = "")
  }
  if (!is.null(x$files$plots)) {
    cat("Plots saved to: ", x$files$plots, "\n", sep = "")
  } else if (!is.null(x$files$plots_failed)) {
    cat("WARNING: The plots PDF was not created:", x$files$plots_failed, "\n")
  }
  invisible(x)
}

#' Plot the sensitivity curves of a gmedsens object
#'
#' Draws the sensitivity plots (estimated path-specific and mediator effects
#' against the sensitivity parameter, with confidence bands), one page per
#' `n_rows * n_cols` plots, on the current graphics device.
#'
#' The arrow in the path titles needs a graphics device with that glyph.
#' Warnings about its substitution are silenced; set `ascii_arrows = TRUE`
#' (or `options(gmediation.ascii_arrows = TRUE)`) to write `->` instead.
#'
#' @param x an object of class `"gmedsens"`.
#' @param ascii_arrows write `->` instead of the arrow in plot titles.
#' @param ... unused.
#' @return `x`, invisibly.
#' @seealso [gmedsens()]
#' @export
plot.gmedsens <- function(x, ascii_arrows = getOption("gmediation.ascii_arrows", FALSE), ...) {
  s <- x$settings
  gs_draw_pages(x$curves, s$CIPCT, s$n_rows, s$n_cols, ascii = ascii_arrows)
  invisible(x)
}
