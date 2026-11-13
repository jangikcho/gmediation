test_that("gs_curve_data summarises bootstrap columns and applies the transform", {
  set.seed(5)
  boot <- matrix(rnorm(6 * 4), 6, 4)
  d <- gmediation:::gs_curve_data(boot, c(0.025, 0.975), 5, c(-0.5, 0, 0.5, 0.9))
  expect_equal(names(d), c("rho", "estimate", "lower", "upper"))
  expect_equal(d$estimate, boot[6, ])
  expect_equal(d$lower, unname(apply(boot, 2, quantile, 0.025)))
  expect_true(all(d$lower <= d$upper))
  e <- gmediation:::gs_curve_data(boot, c(0.025, 0.975), 5, 1:4, transform = exp)
  expect_equal(e$estimate, exp(boot[6, ]))
  # a single value of rho
  one <- gmediation:::gs_curve_data(boot[, 1], c(0.025, 0.975), 5, 0)
  expect_equal(nrow(one), 1)
})

test_that("gs_curve_rows rounds and names the table columns", {
  curve <- data.frame(rho = c(-0.5, 0.5), estimate = c(0.123456, 0.2),
                      lower = c(0.1, 0.15), upper = c(0.2, 0.3))
  tab <- suppressWarnings(gmediation:::gs_curve_rows(curve, c(-0.5, 0.5), 2, Mediator = "m"))
  expect_equal(tab$Estimate, c(0.12, 0.2))
  expect_equal(tab$Mediator, c("m", "m"))
  expect_true(all(c("Lower.CI", "Upper.CI") %in% names(tab)))
  expect_equal(ncol(tab), 5)
})

test_that("gs_pages numbers pages across pairs and kinds", {
  mk <- function(id, a1, a2, kind, title) {
    data.frame(curve_id = id, assum1 = a1, assum2 = a2, kind = kind, title = title,
               ylab = "y", rho = c(-0.5, 0.5), estimate = 1, lower = 0, upper = 2)
  }
  curves <- rbind(
    mk(1, "y", "m1", "path", "p1"), mk(2, "y", "m1", "path", "p2"),
    mk(3, "y", "m1", "path", "p3"), mk(4, "y", "m1", "first", "f1"),
    mk(5, "y", "m2", "path", "p4"), mk(6, "y", "m2", "second", "s1")
  )
  pg <- gmediation:::gs_pages(curves, 1, 2)
  # pair 1: path (3 curves -> 2 pages), first (1); pair 2: path (1), second (1)
  expect_equal(vapply(pg, `[[`, numeric(1), "page"), 1:5)
  expect_equal(vapply(pg, `[[`, character(1), "kind"),
               c("path", "path", "first", "path", "second"))
  expect_equal(pg[[2]]$titles, "p3")
  expect_equal(gmediation:::gs_pages(curves, 1, 2, start_page = 10)[[1]]$page, 10)
  expect_length(gmediation:::gs_pages(NULL, 1, 1), 0)
  expect_equal(gmediation:::gs_page_heading("path"),
               "Estimated Path-Specific Effect (PSE) vs Sensitivity Parameter (rho)")
  expect_match(gmediation:::gs_page_heading("first"), "Total Mediation Effect")
})

test_that("gs_curve_plot builds a ggplot with the requested labels", {
  d <- data.frame(rho = c(-0.5, 0, 0.5), estimate = c(0.1, 0.2, 0.3),
                  lower = c(0, 0.1, 0.2), upper = c(0.2, 0.3, 0.4))
  p <- gmediation:::gs_curve_plot(d, "Title", "Estimated PSE")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Title")
  expect_equal(p$labels$y, "Estimated PSE")
})

test_that("gs_output_dir creates folders and warns when none is given", {
  expect_null(gmediation:::gs_output_dir(c("false", ""), "x", "y"))
  expect_warning(r <- gmediation:::gs_output_dir(c("true", ""), "plot_out[2]", "to save"), "plot_out")
  expect_null(r)
  d <- file.path(tempfile(), "a", "b")
  expect_equal(gmediation:::gs_output_dir(c("true", d), "x", "y"), d)
  expect_true(dir.exists(d))
})

test_that("gs_write_csv writes only the non-empty tables", {
  d <- tempfile(); dir.create(d)
  out <- list(path = data.frame(a = 1), first = NULL, second = data.frame(b = 2))
  files <- gmediation:::gs_write_csv(out, d)
  expect_equal(basename(files), c("gmedsens_pathspec_out.csv", "gmedsens_secondstage_out.csv"))
  expect_true(all(file.exists(files)))
  expect_equal(read.csv(files[1])$a, 1)
})

test_that("gs_plot_text replaces the arrow only when asked", {
  x <- "Path: ses → y"
  expect_identical(gmediation:::gs_plot_text(x, FALSE), x)
  expect_identical(gmediation:::gs_plot_text(x, TRUE), "Path: ses -> y")
})

test_that("gs_mute_glyph_warnings silences only glyph substitution warnings", {
  expect_no_warning(gmediation:::gs_mute_glyph_warnings(
    warning("in 'mbcsToSbcs': -> substituted for →")))
  expect_warning(gmediation:::gs_mute_glyph_warnings(warning("something else")),
                 "something else")
})
