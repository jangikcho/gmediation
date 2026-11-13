test_that("gm_refit refits glm and glm.nb models on new data", {
  m <- bin_models()
  cl <- stats::getCall(m$m1[[1]])
  d <- dental[1:150, c("brush", "ses", "grpvlb", "grpbpd", "race", "SEX")]
  fit <- gmediation:::gm_refit(cl, "binomial", d, environment())
  expect_s3_class(fit, "glm")
  expect_equal(nobs(fit), nrow(na.omit(d)))

  cl_nb <- quote(MASS::glm.nb(formula = dmft ~ ses + brush, data = dental))
  nb <- gmediation:::gm_refit(cl_nb, "negbin", dental[, c("dmft", "ses", "brush")], environment())
  expect_s3_class(nb, "negbin")
})

test_that("gm_dispersion returns variance, theta or NULL by family", {
  m <- gaussian_y_models()
  expect_equal(gmediation:::gm_dispersion(m$y, "gaussian"), summary(m$y)$dispersion)
  nb <- MASS::glm.nb(dmft ~ ses + brush, data = dental)
  expect_equal(gmediation:::gm_dispersion(nb, "negbin"), nb$theta)
  expect_null(gmediation:::gm_dispersion(m$m1[[1]], "binomial"))
})

test_that("gm_draw_m1 and gm_draw_m2 reproduce the fitted means", {
  set.seed(10)
  n <- 20000
  mu <- rep(0.4, n)
  expect_equal(mean(gmediation:::gm_draw_m1("binomial", n, mu, NULL)), 0.4, tolerance = 0.02)
  expect_equal(mean(gmediation:::gm_draw_m1("poisson", n, rep(3, n), NULL)), 3, tolerance = 0.02)
  expect_equal(mean(gmediation:::gm_draw_m1("gaussian", n, rep(2, n), 4)), 2, tolerance = 0.05)
  expect_equal(mean(gmediation:::gm_draw_m1("negbin", n, rep(3, n), 2)), 3, tolerance = 0.05)
  expect_equal(mean(gmediation:::gm_draw_m2("binomial", n, mu, NULL)), 0.4, tolerance = 0.02)
  expect_equal(mean(gmediation:::gm_draw_m2("poisson", n, rep(3, n), NULL)), 3, tolerance = 0.02)
  expect_equal(sd(gmediation:::gm_draw_m2("gaussian", n, rep(2, n), 4)), 2, tolerance = 0.05)
  expect_equal(mean(gmediation:::gm_draw_m2("negbin", n, rep(3, n), 2)), 3, tolerance = 0.05)
  expect_null(gmediation:::gm_draw_m1("gamma", 3, 1, 1))
  expect_null(gmediation:::gm_draw_m2("gamma", 3, 1, 1))
})

test_that("gm_discrepant_draw keeps the counterfactual marginal and respects rho", {
  set.seed(3)
  n <- 20000
  mu <- rep(0.3, n)
  m1j <- rbinom(n, 1, mu)
  cf <- rep(0.6, n)
  d0 <- gmediation:::gm_discrepant_draw("binomial", m1j, mu, cf, NULL, 0, n)
  expect_true(all(d0 %in% c(0, 1)))
  expect_equal(mean(d0), 0.6, tolerance = 0.03)
  # perfect correlation: the counterfactual is a monotone function of the factual draw
  d1 <- gmediation:::gm_discrepant_draw("binomial", m1j, mu, cf, NULL, 1, n)
  expect_true(mean(d1[m1j == 1]) >= mean(d1[m1j == 0]))
  expect_gt(mean(d1[m1j == 1]) - mean(d1[m1j == 0]), 0.5)

  mp <- rep(3, n); xp <- rpois(n, mp)
  dp <- gmediation:::gm_discrepant_draw("poisson", xp, mp, rep(5, n), NULL, 0.5, n)
  expect_equal(mean(dp), 5, tolerance = 0.05)
  expect_true(all(dp == round(dp)))
  dn <- gmediation:::gm_discrepant_draw("negbin", xp, mp, rep(5, n), 2, 0.5, n)
  expect_equal(mean(dn), 5, tolerance = 0.1)

  xg <- rnorm(n, 1, 2)
  dg <- gmediation:::gm_discrepant_draw("gaussian", xg, rep(1, n), rep(3, n), 4, 0.8, n)
  expect_equal(mean(dg), 3, tolerance = 0.05)
  expect_equal(cor(dg, xg), 0.8, tolerance = 0.05)
})
