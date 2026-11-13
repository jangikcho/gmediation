# Internal helpers: sensitivity plots and their layout on pages.

# One sensitivity curve: estimate (line) and confidence band against rho.
gs_curve_plot <- function(data, title, ylab) {
  ggplot(data, aes(x = .data$rho, y = .data$estimate)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2) +
    scale_x_continuous(breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                       labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                       limits = c(-1.0, 1.0)) +
    labs(title = title, x = expression(rho), y = ylab) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10, face = "bold")
    )
}

# Title text for plots: optionally replace the arrow by "->" so that devices
# without the arrow glyph need no substitution.
gs_plot_text <- function(x, ascii = FALSE) {
  if (isTRUE(ascii)) gsub("\u2192", "->", x, fixed = TRUE) else x
}

# Run `expr`, silencing only the warnings R gives when a graphics device
# replaces a character it cannot draw (such as the arrow) by an ASCII one.
gs_mute_glyph_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("mbcsToSbcs|substituted for", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

gs_page_heading <- function(kind) {
  if (kind == "path") {
    "Estimated Path-Specific Effect (PSE) vs Sensitivity Parameter (rho)"
  } else {
    "Estimated Total Mediation Effect vs Sensitivity Parameter (rho)"
  }
}

# Split the curves into pages: for each assumption pair, the path-specific
# curves, then the first-stage and the second-stage mediator curves, with
# `n_rows * n_cols` curves per page. Page numbers continue across pairs.
gs_pages <- function(curves, n_rows, n_cols, start_page = 1) {
  pages <- list()
  if (is.null(curves) || nrow(curves) == 0) return(pages)
  per_page <- n_rows * n_cols
  cu <- curves[!duplicated(curves$curve_id), c("curve_id", "assum1", "assum2", "kind", "title")]
  pair_key <- paste(cu$assum1, cu$assum2, sep = "\r")
  page <- start_page
  for (pk in unique(pair_key)) {
    for (kind in c("path", "first", "second")) {
      sel <- cu[pair_key == pk & cu$kind == kind, , drop = FALSE]
      n <- nrow(sel)
      if (n == 0) next
      for (i in seq(1, n, by = per_page)) {
        rows <- sel[i:min(i + per_page - 1, n), , drop = FALSE]
        pages[[length(pages) + 1]] <- list(
          page = page, kind = kind, assum1 = rows$assum1[1], assum2 = rows$assum2[1],
          ids = rows$curve_id, titles = rows$title
        )
        page <- page + 1
      }
    }
  }
  pages
}

# Draw one page of curves on the current graphics device.
gs_draw_page <- function(pg, curves, CIPCT, n_rows, n_cols, ascii = FALSE) {
  plots <- lapply(pg$ids, function(id) {
    d <- curves[curves$curve_id == id, ]
    gs_curve_plot(d, gs_plot_text(d$title[1], ascii), d$ylab[1])
  })
  heading <- if (pg$kind == "path") {
    expression("Estimated Path-Specific Effect (PSE) vs Sensitivity Parameter (" * rho * ")")
  } else {
    expression("Estimated Total Mediation Effect vs Sensitivity Parameter (" * rho * ")")
  }
  bold <- grid::gpar(fontsize = 10, fontface = "bold")
  title_grob <- gridExtra::arrangeGrob(
    grid::textGrob(heading, x = 0.5, just = "center", gp = bold),
    grid::textGrob(
      paste0("Assumption: No Unmeasured Confounding between '", pg$assum1, "' and '",
             pg$assum2, "'"),
      x = 0.5, just = "center", gp = bold
    ),
    ncol = 1
  )
  subtitle_grob <- grid::textGrob(
    paste0("Solid line: estimates; Shaded gray area: ", CIPCT, "% confidence interval"),
    x = 0.5, just = "center", gp = grid::gpar(fontsize = 8)
  )
  gridExtra::grid.arrange(
    title_grob, subtitle_grob,
    do.call(gridExtra::arrangeGrob, c(plots, ncol = n_cols, nrow = n_rows)),
    ncol = 1, heights = c(0.12, 0.05, 0.83)
  )
}

gs_draw_pages <- function(curves, CIPCT, n_rows, n_cols,
                          ascii = getOption("gmediation.ascii_arrows", FALSE)) {
  gs_mute_glyph_warnings({
    for (pg in gs_pages(curves, n_rows, n_cols)) {
      gs_draw_page(pg, curves, CIPCT, n_rows, n_cols, ascii)
    }
  })
  invisible(NULL)
}

# Open a PDF device. Prefer cairo_pdf() (renders the arrows and rho correctly).
# On some Macs it only warns ("failed to load cairo DLL") and opens no device,
# so a warning, an error or "no new device" all count as failure and the
# base pdf() device is used instead.
gs_open_pdf <- function(path) {
  devs_before <- grDevices::dev.list()
  ok <- tryCatch({
    grDevices::cairo_pdf(path, width = 8, height = 6, onefile = TRUE)
    length(grDevices::dev.list()) > length(devs_before)
  }, error = function(e) {
    message("cairo_pdf() failed: ", conditionMessage(e))
    FALSE
  }, warning = function(w) {
    message("cairo_pdf() unavailable: ", conditionMessage(w))
    FALSE
  })
  if (!ok) {
    message("Using base pdf() instead of cairo_pdf(); arrows in plot titles will be written as '->'.")
    grDevices::pdf(path, width = 8, height = 6, onefile = TRUE)
  }
  structure(grDevices::dev.cur(), cairo = ok)
}

# Write all pages to one PDF file; returns the normalised path, or NULL (with a
# warning) when no file was produced.
gs_write_pdf <- function(curves, path, CIPCT, n_rows, n_cols) {
  previous <- grDevices::dev.cur()
  dev <- gs_open_pdf(path)
  # the base pdf() device has no arrow glyph: use "->" there
  ascii <- !isTRUE(attr(dev, "cairo")) || isTRUE(getOption("gmediation.ascii_arrows", FALSE))
  on.exit({
    if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
    if (previous > 1 && previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)
  gs_draw_pages(curves, CIPCT, n_rows, n_cols, ascii)
  grDevices::dev.off(dev)
  if (file.exists(path) && file.info(path)$size > 0) {
    normalizePath(path)
  } else {
    warning("The plots PDF was not created: ", path, call. = FALSE)
    NULL
  }
}
