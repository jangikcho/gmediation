# Internal helpers for console output: centred plain-text tables.

#' Centre text in a field of a given width
#' @noRd
center_text <- function(x, width) {
  n <- nchar(x)
  left <- floor((width - n) / 2)
  right <- width - n - left
  paste0(strrep(" ", left), x, strrep(" ", right))
}

#' Left-align text in a field of a given width
#' @noRd
left_text <- function(x, width) {
  paste0(x, strrep(" ", width - nchar(x)))
}

#' Print a data frame as a centred plain-text table
#'
#' @param df data frame; numeric columns are shown with `digits` decimals.
#' @param digits number of decimals for numeric columns.
#' @param first_col_left left-align the first column.
#' @noRd
print_centered_table <- function(df, digits = 4, first_col_left = FALSE) {
  df_char <- as.data.frame(
    lapply(df, function(x) {
      if (is.numeric(x)) format(round(x, digits), nsmall = digits, trim = TRUE) else as.character(x)
    }),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  widths <- sapply(seq_along(df_char), function(j) {
    max(nchar(names(df_char)[j]), nchar(df_char[[j]]), na.rm = TRUE)
  })
  headers <- mapply(center_text, names(df_char), widths, USE.NAMES = FALSE)
  cat(paste(headers, collapse = "  "), "\n")
  for (i in seq_len(nrow(df_char))) {
    row_text <- sapply(seq_along(df_char), function(j) {
      value <- as.character(df_char[i, j])
      if (j == 1 && first_col_left) left_text(value, widths[j]) else center_text(value, widths[j])
    })
    cat(paste(row_text, collapse = "  "), "\n")
  }
}

#' Print a titled table block
#' @noRd
print_table_block <- function(title, df, first_col_left = FALSE) {
  cat("\n====================================\n")
  cat(title, "\n", sep = "")
  cat("====================================\n")
  print_centered_table(df, first_col_left = first_col_left)
}
