test_that("gm_count_mediators handles a single model, a list and a named list", {
  m <- bin_models()
  expect_equal(gmediation:::gm_count_mediators(m$m1[[1]]), 1)
  expect_equal(gmediation:::gm_count_mediators(m$m1), 1)
  expect_equal(gmediation:::gm_count_mediators(list(m$m1[[1]], m$m1[[1]])), 2)
  expect_equal(gmediation:::gm_count_mediators(list(y = 1, y = 2, z = 3)), 2)
})

test_that("gm_check_convergence stops on a non-converged model in any position", {
  m <- bin_models()
  warned <- m$m1[[1]]
  warned$th.warn <- "iteration limit reached"
  ok <- gmediation:::gm_check_convergence(m$y, m$m1, m$m2, 1, 1)
  expect_true(ok)
  expect_error(gmediation:::gm_check_convergence(m$y, list(warned), m$m2, 1, 1),
               "convergence")
  y_warned <- m$y
  y_warned$th.warn <- "x"
  expect_error(gmediation:::gm_check_convergence(y_warned, m$m1, m$m2, 1, 1),
               "convergence")
  # the first of two mediators is checked as well (version 1 only checked the last)
  expect_error(
    gmediation:::gm_check_convergence(m$y, list(warned, m$m1[[1]]), m$m2, 2, 1),
    "convergence"
  )
})

test_that("gm_family_name recognises negative binomial fits", {
  m <- bin_models()
  expect_equal(gmediation:::gm_family_name(m$y), "binomial")
  nb <- MASS::glm.nb(dmft ~ ses + brush, data = dental)
  expect_equal(gmediation:::gm_family_name(nb), "negbin")
})

test_that("gm_model_spec extracts variables and drops mediators from covariate sets", {
  sp <- bin_spec()
  expect_equal(sp$y, "dmftd")
  expect_equal(sp$m1[[1]], "brush")
  expect_equal(sp$m2[[1]], "seal")
  expect_equal(sp$npath, 4)
  expect_setequal(sp$y.cov, c("ses", "grpvlb", "grpbpd", "race", "SEX"))
  expect_false("brush" %in% sp$m2.cov)
  expect_equal(sp$family.y, "binomial")
  expect_equal(unlist(sp$family.m1), "binomial")
})

test_that("gm_check_spec reports each kind of invalid specification", {
  sp <- bin_spec()
  expect_true(gmediation:::gm_check_spec(sp, c("ses", 1)))

  bad <- sp; bad$n1 <- 0
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "at least 1 mediator")

  expect_error(gmediation:::gm_check_spec(sp, c("brush", 1)), "Exposure")

  bad <- sp; bad$m2.eachcov[[1]] <- setdiff(bad$m2.eachcov[[1]], "brush")
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "At least one M1")

  bad <- sp; bad$y.allcov <- setdiff(bad$y.allcov, "seal"); bad$y.cov <- bad$y.allcov
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Each M2 must be in Y")

  bad <- sp; bad$m1.eachcov[[1]] <- setdiff(bad$m1.eachcov[[1]], "ses")
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Exposure must be in each M1")

  bad <- sp; bad$m2[[1]] <- "ses"
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Exposure cannot be mediator")
  bad <- sp; bad$m2[[1]] <- "dmftd"
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Final outcome cannot be mediator")

  bad <- sp; bad$m1[[1]] <- "ses"
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Exposure cannot be mediator")
  bad <- sp; bad$m1[[1]] <- "dmftd"; bad$m2.eachcov[[1]] <- c(bad$m2.eachcov[[1]], "dmftd")
  expect_error(gmediation:::gm_check_spec(bad, c("ses", 1)), "Final outcome cannot be mediator")

  bad <- sp; bad$y.cov <- bad$y.cov; bad$m1[[1]] <- "race"
  bad$m2.eachcov[[1]] <- setdiff(bad$m2.eachcov[[1]], "race"); bad$m2.eachcov[[1]] <- c(bad$m2.eachcov[[1]], "race")
  expect_true(gmediation:::gm_check_spec(bad, c("ses", 1)))
})

test_that("gm_decomp_index orders the paths for both decompositions", {
  d1 <- gmediation:::gm_decomp_index(1, 1, 1)
  expect_equal(d1, matrix(c(1, 1, 1, 2, 2, 1, 2, 2), 2))
  d2 <- gmediation:::gm_decomp_index(1, 1, 2)
  expect_equal(d2, d1[, 4:1])
  expect_equal(dim(gmediation:::gm_decomp_index(2, 3, 1)), c(2, 12))
})

test_that("gm_reference_data selects and repeats the reference group", {
  g <- na.omit(dental)
  expect_equal(nrow(gmediation:::gm_reference_data(g, "entire sample", 5)), nrow(g))
  expect_equal(nrow(gmediation:::gm_reference_data(g, "entire sample", 3, repeat_entire = TRUE)),
               3 * nrow(g))
  sub <- gmediation:::gm_reference_data(g, c("race", 0), 1)
  expect_true(all(sub$race == 0))
  sub2 <- gmediation:::gm_reference_data(g, c("race", 0), 2)
  expect_equal(nrow(sub2), 2 * nrow(sub))
})

test_that("gm_prediction_data sets exposure under both values", {
  g <- na.omit(dental)
  sp <- bin_spec()
  pd <- gmediation:::gm_prediction_data(g, g, sp, c("ses", 1), with_m2o = TRUE)
  expect_equal(pd$exposed, 1)
  expect_equal(pd$control, 0)
  expect_true(all(pd$pred.dat.m1.e[[1]]$ses == 1))
  expect_true(all(pd$pred.dat.m1.c[[1]]$ses == 0))
  expect_equal(names(pd$pred.dat.y), sp$y.allcov)
  expect_equal(nrow(pd$pred.dat.m2o[[1]]), nrow(g))
})

test_that("gm_boot_samples draws reproducible resamples with the original last", {
  g <- na.omit(dental)
  set.seed(1); a <- gmediation:::gm_boot_samples(g, 5)
  set.seed(1); b <- gmediation:::gm_boot_samples(g, 5)
  expect_length(a, 6)
  expect_identical(a, b)
  expect_equal(nrow(a[[1]]), nrow(g))
  expect_equal(a[[6]][, names(g)], g)
  expect_false(identical(a[[1]]$ID, a[[2]]$ID))
})

test_that("cluster resampling keeps whole clusters together", {
  g <- na.omit(dental)
  set.seed(2)
  bs <- gmediation:::gm_boot_samples(g, 3, cluster = "momid")
  size <- table(g$momid)
  for (s in bs[1:3]) {
    counts <- table(s$momid)
    # each cluster present is repeated a whole number of times
    expect_true(all(counts %% size[names(counts)] == 0))
  }
})

test_that("gm_complete_data only drops rows with missing model variables", {
  sp <- bin_spec()
  full <- nrow(gmediation:::gm_complete_data(dental, sp, c("ses", 1)))

  d <- dental
  d$unused <- NA_real_ # an unused all-NA column must not remove any row
  expect_equal(nrow(gmediation:::gm_complete_data(d, sp, c("ses", 1))), full)

  d2 <- dental
  d2$brush[which(complete.cases(dental[, c("brush", "seal")]))[1:3]] <- NA
  expect_equal(nrow(gmediation:::gm_complete_data(d2, sp, c("ses", 1))), full - 3)

  d3 <- dental[, setdiff(names(dental), "seal")]
  expect_error(gmediation:::gm_complete_data(d3, sp, c("ses", 1)), "not found")

  d4 <- dental
  d4$momid[1] <- NA
  expect_equal(nrow(gmediation:::gm_complete_data(d4, sp, c("ses", 1), cluster = "momid")),
               nrow(gmediation:::gm_complete_data(d4[-1, ], sp, c("ses", 1), cluster = "momid")))
})

test_that("gm_check_args validates scalar arguments", {
  ok <- function(...) gmediation:::gm_check_args(
    scale = "MeanDiff", boots = 10, CIPCT = 95, decomp = 1, refmult = 1, xcorrv = c(0, 0.5), ...)
  expect_true(ok())
  chk <- function(...) {
    a <- list(scale = "MeanDiff", boots = 10, CIPCT = 95, decomp = 1, refmult = 1, xcorrv = 0)
    a[names(list(...))] <- list(...)
    do.call(gmediation:::gm_check_args, a)
  }
  expect_error(chk(scale = "foo"), "scale")
  expect_error(chk(boots = 0), "boots")
  expect_error(chk(boots = 2.5), "boots")
  expect_error(chk(CIPCT = 100), "CIPCT")
  expect_error(chk(decomp = 3), "decomp")
  expect_error(chk(refmult = 0), "refmult")
  expect_error(chk(xcorrv = 1.5), "xcorrv")
})

test_that("gm_set_seed seeds the generator and restores the previous state", {
  set.seed(99)
  before <- .Random.seed
  restore <- gmediation:::gm_set_seed(1)
  x <- runif(1)
  restore()
  expect_identical(.Random.seed, before)
  set.seed(1)
  expect_equal(x, runif(1))

  # no seed existed before: it is removed again
  rm(".Random.seed", envir = globalenv())
  restore <- gmediation:::gm_set_seed(1)
  restore()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("gm_normalise_ref treats NULL as the entire sample", {
  expect_identical(gm_normalise_ref(NULL), "entire sample")
  expect_identical(gm_normalise_ref("entire sample"), "entire sample")
  expect_identical(gm_normalise_ref(c("SEX", 1)), c("SEX", 1))
  expect_error(gm_normalise_ref("SEX"), "ref")
})
