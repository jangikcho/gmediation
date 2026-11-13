test_that("print.gmediate reports the tables and settings", {
  res <- run_gmediate(bin_models(), xcorrv = 0.5, boots = 3)
  out <- capture.output(r <- print(res))
  expect_identical(r, res)
  expect_true(any(grepl("Path-Specific", out)))
  expect_true(any(grepl("Total", out)))
  expect_true(any(grepl("Bootstrap", out)))
})

test_that("print.gmedsens lists PDF pages and the saved plots when a PDF was written", {
  d <- file.path(tempfile(), "plots")
  res <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"), plot_out = c("true", d))
  out <- capture.output(print(res))
  expect_true(any(grepl("^PDF Page 1$", out)))
  expect_true(any(grepl("Assumption: No Unmeasured Confounding between 'dmftd' and 'brush'", out)))
  expect_true(any(grepl("^Plots saved to: ", out)))
  expect_true(any(grepl("^  Path: ses", out)))
  expect_true(any(grepl("Scale: MeanDiff", out)))
  expect_true(any(grepl("Reference Group: Whole Sample", out)))
})

test_that("print.gmedsens describes a subset reference group and clusters", {
  res <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"), ref = c("race", 0),
                      refmult = 2, boots = 2)
  out <- capture.output(print(res))
  expect_true(any(grepl("Reference Group: race = \\( 0 \\)", out)))
})

test_that("print helpers centre and align text", {
  expect_equal(gmediation:::center_text("ab", 6), "  ab  ")
  expect_equal(gmediation:::center_text("abc", 6), " abc  ")
  expect_equal(gmediation:::left_text("ab", 5), "ab   ")
  out <- capture.output(gmediation:::print_centered_table(
    data.frame(Path = c("a", "bb"), Est = c(1, 2.5)), digits = 2, first_col_left = TRUE))
  expect_length(out, 3)
  expect_match(out[2], "^a")
  expect_match(out[3], "2.50")
  blk <- capture.output(gmediation:::print_table_block("Title", data.frame(x = 1)))
  expect_true("Title" %in% blk)
})

test_that("print output is stable (snapshot)", {
  skip_on_cran()
  res <- run_gmedsens(bin_models(), assumvars = c("dmftd", "brush"))
  expect_snapshot(print(res))
  expect_snapshot(print(run_gmediate(bin_models(), xcorrv = 0.5, boots = 3)))
})
