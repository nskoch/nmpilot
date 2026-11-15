test_that("nm_extract_run_table builds the parameter table and the Word file", {
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  d <- file.path(tempfile("rt"), "run1.reg")
  dir.create(d, recursive = TRUE)
  writeLines(c(
    "TABLE NO.     1: First Order Conditional Estimation with Interaction",
    " ITERATION    THETA1       THETA2       THETA3       SIGMA(1,1)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OBJ",
    ext_row(-1000000000, c(4.2, 45, 0.15, 1, 0.09, 0.01, 0.04, 200)),
    ext_row(-1000000001, c(0.42, 9, 0.015, 1e10, 0.018, 0.005, 0.01, 0)),
    ext_row(-1000000006, c(0, 0, 0, 1, 0, 0, 0, 0))),
    file.path(d, "run1.ext"))
  res <- NULL
  invisible(capture.output(suppressMessages(res <- nm_extract_run_table(run_path = d))))
  expect_true(file.exists(res$output_file))
  p <- res$params
  expect_equal(p$Estimate[p$Type == "THETA"], c(4.2, 45, 0.15))
  # ETA SD derived from the variance, correlation from the covariance.
  expect_equal(p$Estimate[p$Type == "ETA_SD"], sqrt(c(0.09, 0.04)))
  expect_equal(p$Estimate[p$Type == "ETA_COR"], 0.01 / sqrt(0.09 * 0.04))
  expect_true(p$Fixed[p$Type == "SIGMA"])
})
