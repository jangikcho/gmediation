test_that("gmedsens returns a classed object and prints nothing when print_report = FALSE", {
  expect_silent(res <- run_gmedsens(bin_models()))
  expect_s3_class(res, "gmedsens")
  expect_true(all(c("path", "first.stage", "second.stage", "curves", "settings", "files") %in%
                    names(res)))
  expect_equal(names(res$path)[c(1, 2)], c("Path", "Xcorr"))
  expect_true("Rho.Between" %in% names(res$path))
  expect_equal(unique(res$curves$kind), c("path", "first", "second"))
  expect_true(all(res$curves$lower <= res$curves$upper + 1e-12))
  expect_length(res$files, 0)
})

test_that("gmedsens evaluates every assumption pair for ALL", {
  res <- run_gmedsens(bin_models())
  pairs <- unique(res$curves[, c("assum1", "assum2")])
  # Y vs M1, Y vs M2, M2 vs M1
  expect_equal(nrow(pairs), 3)
  expect_setequal(unique(res$path$Rho.Between), c("dmftd & brush", "dmftd & seal", "seal & brush"))
  # one rho grid per curve
  n_rho <- 2 / 0.2 - 1
  expect_equal(nrow(res$curves) %% n_rho, 0)
})

test_that("a specific pair gives the same curves as the matching part of ALL", {
  single <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"))
  expect_equal(unique(single$path$Rho.Between), "dmftd & brush")
  expect_true(all(single$first.stage$Mediator == "brush"))
  expect_null(single$second.stage)
})

test_that("gmedsens is reproducible and keeps the user's random numbers", {
  m <- bin_models()
  set.seed(8)
  before <- .Random.seed
  a <- run_gmedsens(m, assumvars = c("dmftd", "brush"))
  expect_identical(.Random.seed, before)
  b <- run_gmedsens(m, assumvars = c("dmftd", "brush"))
  expect_identical(a$curves, b$curves)
})

test_that("sensitivity curves agree across assumption pairs at rho = 0", {
  skip_on_cran()
  m <- bin_models()
  res <- run_gmedsens(m, boots = 2, xcorrv = 0, refmult = 60, deltsi = 0.5, seed = 7)
  x <- res$curves[res$curves$kind == "first" & abs(res$curves$rho) < 1e-9, ]
  expect_equal(nrow(x), 3)
  expect_lt(diff(range(x$estimate)), 0.01)
})

test_that("gmedsens supports the OR and RR scales", {
  m <- bin_models()
  or <- run_gmedsens(m, scale = "OR", assumvars = c("dmftd", "brush"))
  rr <- run_gmedsens(m, scale = "RR", assumvars = c("dmftd", "brush"))
  expect_s3_class(or, "gmedsens")
  expect_s3_class(rr, "gmedsens")
  expect_true(all(or$curves$estimate[or$curves$kind == "path"] > 0))
})

test_that("gmedsens works with mixed mediators, clusters and a Gaussian outcome", {
  skip_on_cran()
  mix <- run_gmedsens(mix_models(), cluster = "momid", boots = 2, xcorrv = c(0, 0.9),
                      assumvars = c("dmftd", "ohi"))
  expect_s3_class(mix, "gmedsens")
  expect_equal(unique(mix$path$Rho.Between), "dmftd & ohi")
  g <- run_gmedsens(gaussian_y_models(), assumvars = c("ohi", "seal"))
  expect_s3_class(g, "gmedsens")
})

test_that("gmedsens writes CSV tables and the plots PDF", {
  d <- file.path(tempfile(), "out")
  res <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"),
                      plot_out = c("true", d), table_out = c("true", d))
  expect_setequal(basename(res$files$tables),
                  c("gmedsens_pathspec_out.csv", "gmedsens_firststage_out.csv"))
  expect_true(file.exists(res$files$plots))
  expect_gt(file.info(res$files$plots)$size, 0)
  tab <- read.csv(file.path(d, "gmedsens_pathspec_out.csv"))
  expect_equal(nrow(tab), nrow(res$path))
  expect_equal(tab$Rho.Between[1], "dmftd & brush")
})

test_that("plot and print methods work", {
  res <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"), n_cols = 2)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(suppressWarnings(plot(res)))
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl("Final Outcome: dmftd", out)))
  expect_true(any(grepl("Sample Size Used:", out)))
  expect_false(any(grepl("PDF Page", out)))
})

test_that("gmedsens checks its inputs", {
  m <- bin_models()
  expect_error(run_gmedsens(m, deltsi = 0.3), "deltsi")
  expect_error(run_gmedsens(m, scale = "x"), "scale")
  expect_error(run_gmedsens(m, assumvars = "dmftd"), "assumvars")
  expect_error(run_gmedsens(m, assumvars = c("dmftd", "nothing")), "assumvars")
  expect_warning(res <- run_gmedsens(m, g = 2, assumvars = c("dmftd", "brush")), "g value")
  expect_equal(res$settings$g, 1)
  expect_warning(run_gmedsens(m, assumvars = c("dmftd", "brush"), plot_out = c("true", "")),
                 "plot_out")
})

test_that("gmedsens rejects first-stage mediators it cannot handle", {
  m <- list(
    y = glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + dmft + seal, family = gaussian, data = dental),
    m1 = list(glm(dmft ~ ses + grpvlb + grpbpd + race + SEX, family = poisson, data = dental)),
    m2 = list(glm(seal ~ ses + grpvlb + grpbpd + race + SEX + dmft, family = binomial, data = dental))
  )
  expect_error(run_gmedsens(m), "first-stage")
})

test_that("missing values outside the model variables do not drop observations", {
  d <- dental
  d$unused <- NA_real_
  res <- run_gmedsens(bin_models(), data = d, assumvars = c("dmftd", "brush"))
  expect_equal(res$settings$n, nrow(na.omit(dental)))
})

test_that("gmedsens prints its report even when the result is assigned", {
  m <- bin_models()
  d <- file.path(tempdir(), "gs_report")
  out <- capture.output(res <- gmedsens(
    dental, m$m1, m$m2, m$y, expos = c("ses", 1), seed = 1, boots = 3, xcorrv = 0,
    refmult = 1, decomp = 1, deltsi = 0.5, table_out = c("true", d)))
  expect_s3_class(res, "gmedsens")
  expect_true(any(grepl("Exposure", out)))
  expect_true(any(grepl("Tables saved to", out)))
  expect_true(any(grepl(basename(d), out, fixed = TRUE)))
})
