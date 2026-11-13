test_that("gs_correction is zero at rho = 0 and has the documented forms", {
  expect_equal(gmediation:::gs_correction(0, 1, "gaussian", "gaussian", NA, 5), 0)
  expect_equal(gmediation:::gs_correction(0, 1, "binomial", "gaussian", 0.3, 5), 0)

  # non-binomial partner
  expect_equal(gmediation:::gs_correction(0.5, 2, "gaussian", "binomial", NA, 3),
               (0.5 / 2) * (1 / sqrt(1 - 0.25)) * 3)
  # binomial partner, non-binomial outcome
  p <- 0.3
  expect_equal(gmediation:::gs_correction(0.5, 2, "binomial", "gaussian", p, 3),
               (0.5 / 2) * (0.5 + (1 + dnorm(qnorm(p)) / 2) / 4) * 3)
  # both binomial; the power is 2 for Y vs M1 and 1 otherwise
  expect_equal(gmediation:::gs_correction(0.5, 2, "binomial", "binomial", p, 3, binom_power = 2),
               (0.5 / 2) * (dnorm(qnorm(p)) / 2)^2 * 3)
  expect_equal(gmediation:::gs_correction(0.5, 2, "binomial", "binomial", p, 3),
               (0.5 / 2) * (dnorm(qnorm(p)) / 2) * 3)
  # odd in rho
  expect_equal(gmediation:::gs_correction(-0.5, 2, "gaussian", "gaussian", NA, 3),
               -gmediation:::gs_correction(0.5, 2, "gaussian", "gaussian", NA, 3))
})

test_that("gs_corrected_m2 reproduces the M2 distribution at rho = 0 (regression)", {
  set.seed(1)
  n <- 50000
  u <- runif(n)
  f <- gmediation:::gs_corrected_m2
  # version 1 drew M2 = 1 with probability 1 - mu in the binomial case
  expect_equal(mean(f("binomial", rep(0.3, n), NULL, NULL, 0, 0, u, n)), 0.3, tolerance = 0.02)
  expect_equal(mean(f("poisson", rep(2, n), NULL, NULL, 0, 0, u, n)), 2, tolerance = 0.02)
  expect_equal(mean(f("negbin", rep(2, n), NULL, 3, 0, 0, u, n)), 2, tolerance = 0.03)
  expect_equal(mean(f("gaussian", rep(2, n), 1.5, NULL, 0, 0, u, n)), 2, tolerance = 0.02)
  # a positive correction increases M2 in every family
  up <- rep(0.8, n)
  expect_gt(mean(f("binomial", rep(0.3, n), NULL, NULL, up, 0.5, u, n)), 0.3)
  expect_gt(mean(f("poisson", rep(2, n), NULL, NULL, up, 0.5, u, n)), 2)
  expect_gt(mean(f("gaussian", rep(2, n), 1.5, NULL, up, 0.5, u, n)), 2)
  expect_null(f("gamma", 1, 1, 1, 0, 0, 0.5, 1))
})

test_that("gs_outcome_mean equals the unperturbed mean at crct = 0", {
  ctx <- list(nref = 5, tol = 1e-4, maxcnt = 100)
  mu <- c(0.2, 0.4, 0.6, 0.8, 0.5)
  expect_equal(gmediation:::gs_outcome_mean(mu, 0, 0, "binomial", NA, NA, ctx, 1), mu)
  expect_equal(gmediation:::gs_outcome_mean(c(1, 2, 3), 0, 0, "gaussian", 2, NA, ctx, 1), c(1, 2, 3))
  ctx$nref <- 2000
  set.seed(4)
  mu_c <- rep(2, 2000)
  expect_equal(mean(gmediation:::gs_outcome_mean(mu_c, rep(0, 2000), 0, "poisson", NA, NA, ctx, 1)),
               2, tolerance = 0.01)
  expect_gt(mean(gmediation:::gs_outcome_mean(mu_c, rep(1, 2000), 0.5, "poisson", NA, NA, ctx, 1)), 2)
  expect_null(gmediation:::gs_outcome_mean(mu_c, 0, 0, "gamma", NA, NA, ctx, 1))
})

test_that("gs_scale_effect handles the scales and checks the range", {
  f <- gmediation:::gs_scale_effect
  expect_equal(f(c(0.6, 0.5), c(0.3, 0.2), "MeanDiff"), c(0.3, 0.3))
  expect_equal(f(0.6, 0.3, ""), 0.3)
  expect_equal(f(c(0.6, 0.5), c(0.3, 0.2), "OR"),
               log((c(0.6, 0.5) / c(0.4, 0.5)) / (c(0.3, 0.2) / c(0.7, 0.8))))
  expect_equal(f(c(0.6, 0.5), c(0.3, 0.25), "RR"), log(c(2, 2)))
  expect_error(f(c(0.6, 1.2), c(0.3, 0.2), "OR"), "not in \\(0,1\\)")
  expect_error(f(c(0.6, 0.5), c(0.3, -0.1), "RR"), "negative")
  expect_equal(f(c(1, 2), c(1, 1), "other"), c(0, 0))
})

test_that("gs_path_shown hides paths the assumption cannot affect", {
  eq <- matrix(0, 2, 2)
  expect_equal(gmediation:::gs_path_shown(eq, eq, 1, 1, 1), 0)
  DM <- matrix(1, 2, 2); DM[1, 1] <- 0
  expect_equal(gmediation:::gs_path_shown(DM, matrix(1, 2, 2), 1, 1, 1), 1)
  # Y vs M2 uses the second-stage column
  expect_equal(gmediation:::gs_path_shown(eq, eq, 1, 2, 1), 0)
  expect_equal(gmediation:::gs_path_shown(DM, matrix(1, 2, 2), 1, 2, 1), 1)
  # M2 vs M1
  expect_equal(gmediation:::gs_path_shown(eq, eq, 2, 1, 1), 0)
  DM2 <- matrix(1, 2, 2); DM2[1, 2] <- 0
  expect_equal(gmediation:::gs_path_shown(DM2, matrix(1, 2, 2), 2, 1, 1), 1)
})

test_that("gs_is_discrepant detects a mediator set differently across stages", {
  sp <- bin_spec()
  ctx <- list(spec = sp)
  DM <- matrix(1, 2, 2)
  expect_false(gmediation:::gs_is_discrepant(DM, ctx))
  DM[2, 2] <- 0 # M1 -> M2 link switched, exposure -> M1 link not
  expect_true(gmediation:::gs_is_discrepant(DM, ctx))
})

test_that("gs_assumption_pairs lists the pairs for ALL and checks specific ones", {
  sp <- gmediation:::gm_model_spec(mix_models()$y, mix_models()$m1, mix_models()$m2, 2, 2)
  all <- gmediation:::gs_assumption_pairs("ALL", sp)
  # (Y vs 4 mediators) + 2 x (M2 vs 2 M1)
  expect_equal(nrow(all), 4 + 2 * 2)
  expect_equal(all$i1, c(1, 1, 1, 1, 2, 2, 3, 3))
  expect_equal(all$i2, c(1, 2, 3, 4, 1, 2, 1, 2))

  p <- gmediation:::gs_assumption_pairs(c("dmftd", "ohi"), sp)
  expect_equal(c(p$i1, p$i2), c(1, 4))
  p <- gmediation:::gs_assumption_pairs(c("dmftd", "visit"), sp)
  expect_equal(c(p$i1, p$i2), c(1, 2))
  p <- gmediation:::gs_assumption_pairs(c("ohi", "brush"), sp)
  expect_equal(c(p$i1, p$i2), c(3, 1))
  expect_error(gmediation:::gs_assumption_pairs(c("dmftd", "nothing"), sp), "assumvars")
  expect_error(gmediation:::gs_assumption_pairs(c("nothing", "brush"), sp), "assumvars")
  expect_error(gmediation:::gs_assumption_pairs("dmftd", sp), "assumvars")
})

test_that("gs_get returns NULL outside the list", {
  expect_equal(gmediation:::gs_get(list(1, 2), 2), 2)
  expect_null(gmediation:::gs_get(list(1, 2), 3))
})
