test_that("dental data loads with the expected shape", {
  expect_s3_class(dental, "data.frame")
  expect_equal(nrow(dental), 205)
  expect_true(all(c("ID", "dmftd", "ses", "brush", "seal", "visit", "ohi",
                     "grpvlb", "grpbpd", "race", "SEX", "momid") %in% names(dental)))
})

test_that("gmediate runs with the default (entire sample) reference group", {
  model.y1 <- glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + seal,
                   family = binomial, data = dental)
  model.m11 <- glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
                    family = binomial, data = dental)
  model.m21 <- glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush,
                    family = binomial, data = dental)

  # small number of bootstrap samples to keep the test fast
  res <- gmediate(data = dental, model.m1 = list(model.m11),
                   model.m2 = list(model.m21), model.y = model.y1,
                   expos = c("ses", 1), seed = 1234, boots = 5,
                   cluster = NULL, decomp = 1, CIPCT = 95,
                   xcorrv = c(0, 0.9), refmult = 1)

  expect_type(res, "list")
  expect_true(all(c("indiv.path", "total.path", "first.path", "second.path",
                     "prop.indiv.path", "prop.first.path", "prop.second.path",
                     "boot.est", "model.fit") %in% names(res)))
})
