# Regression tests against frozen results.
#  * g_bin11_d1, g_bin11_d2, g_mix22_clu, g_ygau, s_bin11_spec, s_mix22_spec:
#    output of gmediation 1.0.0 -- version 2 must reproduce it.
#  * g_bin11_RR, s_bin11_spec2: output of version 2 where version 1 was wrong
#    (risk-ratio path effects were all zero; the binary second-stage mediator
#    correction in gmedsens() had the wrong direction).

golden <- function() readRDS(test_path("fixtures", "golden.rds"))

gmediate_parts <- c("indiv.path", "total.path", "first.path", "second.path",
                    "prop.indiv.path", "prop.first.path", "prop.second.path", "boot.est")

compare_gmediate <- function(res, gold) {
  for (nm in gmediate_parts) {
    expect_equal(unclass(res[[nm]]), unclass(gold[[nm]]), tolerance = 1e-6, info = nm)
  }
}

compare_gmedsens <- function(res, gold) {
  tabs <- list(gmedsens_pathspec_out.csv = res$path,
               gmedsens_firststage_out.csv = res$first.stage,
               gmedsens_secondstage_out.csv = res$second.stage)
  tabs <- Filter(Negate(is.null), tabs)
  expect_setequal(names(tabs), names(gold))
  for (nm in names(gold)) {
    got <- tabs[[nm]]
    names(got)[names(got) %in% c("X.U.03C1.", "\u03c1")] <- "rho"
    expect_equal(dim(got), dim(gold[[nm]]), info = nm)
    for (col in names(gold[[nm]])) {
      exp <- gold[[nm]][[col]]
      if (is.logical(exp) && all(is.na(exp))) exp <- as.numeric(exp)
      expect_equal(got[[col]], exp, tolerance = 1e-6, info = paste(nm, col))
    }
  }
}

test_that("gmediate reproduces version 1 results", {
  skip_on_cran()
  gold <- golden()
  m <- bin_models()
  compare_gmediate(run_gmediate(m, boots = 8, xcorrv = c(0, 0.9), decomp = 1), gold$g_bin11_d1)
  compare_gmediate(run_gmediate(m, boots = 8, xcorrv = c(0, 0.9), decomp = 2), gold$g_bin11_d2)
  compare_gmediate(run_gmediate(gaussian_y_models(), boots = 8, xcorrv = c(0, 0.9)), gold$g_ygau)
  compare_gmediate(
    run_gmediate(mix_models(), boots = 6, xcorrv = c(0, 0.5, 0.9), cluster = "momid", decomp = 2),
    gold$g_mix22_clu
  )
})

test_that("gmedsens reproduces version 1 results for single assumption pairs", {
  skip_on_cran()
  gold <- golden()
  compare_gmedsens(
    run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"), deltsi = 0.05),
    gold$s_bin11_spec
  )
  compare_gmedsens(
    run_gmedsens(mix_models(), boots = 2, xcorrv = c(0, 0.9), assumvars = c("dmftd", "ohi"),
                 deltsi = 0.2),
    gold$s_mix22_spec
  )
})

test_that("version 2 results where version 1 was incorrect are stable", {
  skip_on_cran()
  gold <- golden()
  compare_gmediate(run_gmediate(bin_models(), boots = 8, xcorrv = c(0, 0.9), scale = "RR"),
                   gold$g_bin11_RR)
  compare_gmedsens(
    run_gmedsens(bin_models(), assumvars = c("seal", "brush"), deltsi = 0.05),
    gold$s_bin11_spec2
  )
})
