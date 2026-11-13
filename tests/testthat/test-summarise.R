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
