# Models on the `dental` data used throughout the tests. Formulas are written
# out in full so that the stored calls can be re-evaluated on bootstrap samples.

bin_models <- function() {
  list(
    y = glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + seal,
            family = binomial, data = dental),
    m1 = list(glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
                  family = binomial, data = dental)),
    m2 = list(glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush,
                  family = binomial, data = dental))
  )
}

mix_models <- function() {
  list(
    y = glm(dmftd ~ ses + grpvlb + grpbpd + race + SEX + brush + visit + seal + ohi,
            family = binomial, data = dental),
    m1 = list(
      glm(brush ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental),
      glm(visit ~ ses + grpvlb + grpbpd + race + SEX, family = binomial, data = dental)
    ),
    m2 = list(
      glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
          family = binomial, data = dental),
      glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + visit,
          family = gaussian, data = dental)
    )
  )
}

gaussian_y_models <- function() {
  list(
    y = glm(ohi ~ ses + grpvlb + grpbpd + race + SEX + brush + seal,
            family = gaussian, data = dental),
    m1 = list(glm(brush ~ ses + grpvlb + grpbpd + race + SEX,
                  family = binomial, data = dental)),
    m2 = list(glm(seal ~ ses + grpvlb + grpbpd + race + SEX + brush,
                  family = binomial, data = dental))
  )
}

run_gmediate <- function(m, ...) {
  args <- list(data = dental, model.m1 = m$m1, model.m2 = m$m2, model.y = m$y,
               expos = c("ses", 1), seed = 1234, boots = 4, xcorrv = c(0, 0.9),
               refmult = 1, decomp = 1)
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(gmediate, args)
}

run_gmedsens <- function(m, ...) {
  args <- list(data = dental, model.m1 = m$m1, model.m2 = m$m2, model.y = m$y,
               expos = c("ses", 1), seed = 1234, boots = 3, xcorrv = c(0, 0.5),
               refmult = 1, decomp = 1, deltsi = 0.2, print_report = FALSE)
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(gmedsens, args)
}

# spec / context for unit tests of the internal helpers
bin_spec <- function() {
  m <- bin_models()
  gmediation:::gm_model_spec(m$y, m$m1, m$m2, 1, 1)
}
