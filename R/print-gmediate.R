#' Print a gmediate result
#'
#' Prints the analysis settings followed by the seven summary tables.
#'
#' @param x an object of class \code{"gmediate"}, as returned by
#'   \code{\link{gmediate}}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{gmediate}}
#' @export
print.gmediate <- function(x, ...) {
  s <- x$settings
  cat("----------------------------------------\n")
  cat("Causal Mediation Sensitivity Analysis", "\n")
  cat("Exposure ('exposed' value):", paste0(s$expos[1], "(", s$expos[2], ")", "\n"))
  cat("Final Outcome:", s$y, "\n")
  cat("Sample Size Used:", s$n, "\n")
  cat("Number of Bootstrap Samples / Number Used:", s$boots, "/", s$boots, "\n")
  if (!identical(s$ref, "entire sample")) {
    cat("Reference Group:", s$ref[1], "= (", s$ref[2:length(s$ref)], "),",
        "Reference Group Multiplier (refmult) =", s$refmult, "\n")
  } else {
    cat("Reference Group: Whole Sample,", "Reference Group Multiplier (refmult)=", s$refmult, "\n")
  }
  cat("Scale:", s$scale, "\n")
  if (!is.null(s$cluster)) cat("Cluster:", s$cluster, "\n")
  cat("Confidence Interval Level:", s$CIPCT, "\n")
  cat("----------------------------------------\n")

  print_table_block("Path-Specific Effects", x$indiv.path, first_col_left = TRUE)
  print_table_block("Total Effect", x$total.path)
  print_table_block("First Stage Mediators", x$first.path)
  print_table_block("Second Stage Mediators", x$second.path)
  print_table_block("Proportion: Path-Specific Effects", x$prop.indiv.path, first_col_left = TRUE)
  print_table_block("Proportion: First Stage Mediators", x$prop.first.path)
  print_table_block("Proportion: Second Stage Mediators", x$prop.second.path)
  invisible(x)
}
