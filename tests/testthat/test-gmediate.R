test_that("dental data loads with the expected shape", {
  expect_s3_class(dental, "data.frame")
  expect_equal(nrow(dental), 205)
  expect_true(all(c("ID", "dmftd", "ses", "brush", "seal", "visit", "ohi",
                    "grpvlb", "grpbpd", "race", "SEX", "momid") %in% names(dental)))
})

test_that("gmediate returns a classed object and prints nothing by itself", {
  expect_silent(res <- run_gmediate(bin_models()))
  expect_s3_class(res, "gmediate")
  expect_true(all(c("indiv.path", "total.path", "first.path", "second.path",
                    "prop.indiv.path", "prop.first.path", "prop.second.path",
                    "boot.est", "model.fit", "settings") %in% names(res)))
  expect_equal(res$settings$boots, 4)
  expect_equal(res$settings$n, nrow(na.omit(dental)))
  expect_equal(nrow(res$total.path), 1)
  expect_equal(rownames(res$total.path), "Total")
  expect_length(res$model.fit, 3)
  # printing is explicit
  expect_output(print(res), "Individual Path|Total")
})

test_that("gmediate is reproducible and leaves the user's random numbers alone", {
  m <- bin_models()
  set.seed(42)
  before <- .Random.seed
  a <- run_gmediate(m)
  expect_identical(.Random.seed, before)
  b <- run_gmediate(m)
  expect_identical(a$indiv.path, b$indiv.path)
  expect_identical(a$boot.est, b$boot.est)
  c2 <- run_gmediate(m, seed = 99)
  expect_false(identical(a$boot.est, c2$boot.est))
})

test_that("gmediate path effects add up to the total effect", {
  res <- run_gmediate(bin_models(), xcorrv = 0, boots = 3)
  est <- res$indiv.path$Estimate
  est <- suppressWarnings(as.numeric(est))
  expect_equal(sum(est, na.rm = TRUE), res$total.path$Estimate, tolerance = 1e-3)
})

test_that("gmediate supports decomposition 2, subsets of the reference group and clusters", {
  m <- bin_models()
  d1 <- run_gmediate(m, decomp = 1)
  d2 <- run_gmediate(m, decomp = 2)
  expect_s3_class(d2, "gmediate")
  expect_equal(d1$total.path$Estimate, d2$total.path$Estimate)
  r <- run_gmediate(m, ref = c("race", 0), refmult = 2)
  expect_s3_class(r, "gmediate")
  cl <- run_gmediate(mix_models(), cluster = "momid", xcorrv = c(0, 0.5), boots = 3)
  expect_s3_class(cl, "gmediate")
  expect_equal(cl$settings$cluster, "momid")
})

test_that("the risk ratio scale returns non-zero path effects (regression)", {
  res <- run_gmediate(bin_models(), scale = "RR", xcorrv = 0)
  est <- suppressWarnings(as.numeric(res$indiv.path$Estimate))
  expect_true(any(abs(est) > 1e-8, na.rm = TRUE))
  or <- run_gmediate(bin_models(), scale = "OR", xcorrv = 0)
  expect_s3_class(or, "gmediate")
})

test_that("a Gaussian outcome and negative binomial mediators are supported", {
  expect_s3_class(run_gmediate(gaussian_y_models()), "gmediate")
  nb <- list(
    y = glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + dmft,
            family = gaussian, data = dental),
    m1 = list(glm(brush ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental)),
    m2 = list(MASS::glm.nb(dmft ~ ses + grpvlb + grpbpd + race + SEX + brush, data = dental))
  )
  res <- run_gmediate(nb, boots = 3)
  expect_s3_class(res, "gmediate")
  expect_length(res$model.fit, 3)
})

test_that("missing values outside the model variables do not drop observations", {
  d <- dental
  d$unused <- NA_real_
  m <- bin_models()
  res <- run_gmediate(m, data = d)
  expect_equal(res$settings$n, nrow(na.omit(dental)))
})

test_that("gmediate checks its inputs", {
  m <- bin_models()
  expect_error(run_gmediate(m, scale = "bad"), "scale")
  expect_error(run_gmediate(m, boots = 0), "boots")
  expect_error(run_gmediate(m, expos = c("brush", 1)), "Exposure")
  bad <- m$y; bad$th.warn <- "x"
  expect_error(run_gmediate(list(y = bad, m1 = m$m1, m2 = m$m2)), "convergence")
})

test_that("gmediate writes the summary tables to an Excel workbook", {
  skip_if_not_installed("writexl")
  d <- file.path(tempfile(), "out")
  res <- run_gmediate(bin_models(), table_out = c("true", d))
  f <- file.path(d, "gmediate_summary_tables.xlsx")
  expect_true(file.exists(f))
  expect_s3_class(res, "gmediate")
  expect_warning(run_gmediate(bin_models(), table_out = c("true", "")), "table_out")
})

test_that("ref = NULL behaves like the entire sample", {
  m <- bin_models()
  a <- run_gmediate(m, boots = 3, xcorrv = c(0, 0.9))
  b <- gmediate(dental, m$m1, m$m2, m$y, expos = c("ses", 1), ref = NULL, seed = 1234,
                boots = 3, xcorrv = c(0, 0.9), refmult = 1, decomp = 1)
  expect_equal(unclass(a$indiv.path), unclass(b$indiv.path))
})
