test_that("gm_ci_probs gives the quantile probabilities of a CI", {
  expect_equal(gmediation:::gm_ci_probs(95), c(0.025, 0.975))
  expect_equal(gmediation:::gm_ci_probs(90), c(0.05, 0.95))
})

test_that("gm_boot_pvalue is a two-sided empirical p-value", {
  expect_equal(gmediation:::gm_boot_pvalue(1:100), "0.0000")
  x <- c(-2, -1, 1, 2, 3, 4, 5, 6, 7, 8)
  expect_equal(gmediation:::gm_boot_pvalue(x), format(0.4, nsmall = 4))
  expect_type(gmediation:::gm_boot_pvalue(rnorm(50)), "character")
})

test_that("gm_boot_summary returns rounded percentile limits and p-value", {
  x <- seq(-1, 3, length.out = 101)
  s <- gmediation:::gm_boot_summary(x, 0.025, 0.975)
  expect_equal(s$lower, round(unname(quantile(x, 0.025)), 4))
  expect_equal(s$upper, round(unname(quantile(x, 0.975)), 4))
  expect_equal(s$p, gmediation:::gm_boot_pvalue(x))
  s_na <- gmediation:::gm_boot_summary(c(x, NA), 0.025, 0.975, na.rm = TRUE)
  expect_equal(s_na$lower, s$lower)
})

test_that("gm_path_effect computes differences, odds ratios and risk ratios", {
  expect_equal(gmediation:::gm_path_effect(0.6, 0.3, "MeanDiff"), 0.3)
  expect_equal(gmediation:::gm_path_effect(0.6, 0.3, ""), 0.3)
  expect_equal(gmediation:::gm_path_effect(0.6, 0.3, "OR"), log((0.6 / 0.4) / (0.3 / 0.7)))
  # regression test: version 1 returned log(e_hi / e_hi) = 0 for RR
  expect_equal(gmediation:::gm_path_effect(0.6, 0.3, "RR"), log(2))
  expect_equal(gmediation:::gm_path_effect(c(0.6, 0.5), c(0.3, 0.5), "RR"), c(log(2), 0))
})

test_that("gm_total_effect uses the same scales", {
  expect_equal(gmediation:::gm_total_effect(0.6, 0.3, "MeanDiff"), 0.3)
  expect_equal(gmediation:::gm_total_effect(0.6, 0.3, "OR"), log((0.6 / 0.4) / (0.3 / 0.7)))
  expect_equal(gmediation:::gm_total_effect(0.6, 0.3, "RR"), log(2))
})

test_that("gm_boot_summary back-transforms the limits but tests on the original scale", {
  set.seed(2)
  x <- rnorm(200, mean = 0.5, sd = 0.2)
  a <- gmediation:::gm_boot_summary(x, 0.025, 0.975)
  b <- gmediation:::gm_boot_summary(x, 0.025, 0.975, back = exp)
  expect_equal(b$lower, round(exp(unname(quantile(x, 0.025))), 4))
  expect_equal(b$upper, round(exp(unname(quantile(x, 0.975))), 4))
  expect_identical(a$p, b$p)
  # a log-ratio centred at 0 is not "significant" (p is not 0)
  y <- rnorm(200, 0, 0.1)
  expect_gt(as.numeric(gmediation:::gm_boot_summary(y, 0.025, 0.975, back = exp)$p), 0.2)
})

test_that("gm_scale_back and gm_estimate_label depend on the scale", {
  expect_identical(gmediation:::gm_scale_back("RR"), exp)
  expect_identical(gmediation:::gm_scale_back("OR"), exp)
  expect_identical(gmediation:::gm_scale_back("MeanDiff"), identity)
  expect_identical(gmediation:::gm_scale_back("RR", proportion = TRUE), identity)
  expect_equal(gmediation:::gm_estimate_label("RR"), "RR")
  expect_equal(gmediation:::gm_estimate_label("OR"), "OR")
  expect_equal(gmediation:::gm_estimate_label("MeanDiff"), "Estimate")
  expect_equal(gmediation:::gm_estimate_label("RR", proportion = TRUE), "Est Prop")
})
