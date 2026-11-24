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
  invisible(capture.output(suppressMessages(res <- nm_extract_run_table(d))))
  expect_true(file.exists(res$output_file))
  p <- res$params
  expect_equal(p$Estimate[p$Type == "THETA"], c(4.2, 45, 0.15))
  # ETA SD derived from the variance, correlation from the covariance.
  expect_equal(p$Estimate[p$Type == "ETA_SD"], sqrt(c(0.09, 0.04)))
  expect_equal(p$Estimate[p$Type == "ETA_COR"], 0.01 / sqrt(0.09 * 0.04))
  expect_true(p$Fixed[p$Type == "SIGMA"])
})

test_that("nm_extract_run_table finds runs without .reg/.nm7 directories", {
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  ext_lines <- c(
    "TABLE NO.     1: First Order Conditional Estimation with Interaction",
    " ITERATION    THETA1       SIGMA(1,1)   OMEGA(1,1)   OBJ",
    ext_row(-1000000000, c(4.2, 1, 0.09, 200)),
    ext_row(-1000000001, c(0.42, 0.1, 0.018, 0)))
  # Flat layout: run1.ext next to run1.ctl
  root <- tempfile("rtflat"); dir.create(root)
  writeLines(ext_lines, file.path(root, "run1.ext"))
  res <- NULL
  invisible(capture.output(suppressMessages(
    res <- nm_extract_run_table(file.path(root, "run1")))))
  expect_equal(res$params$Estimate[res$params$Type == "THETA"], 4.2)
  expect_equal(normalizePath(res$run_dir), normalizePath(root))
  invisible(capture.output(suppressMessages(
    res2 <- nm_extract_run_table(file.path(root, "run1.ctl")))))
  expect_equal(res2$extension, "flat")
  # Plain folder named after the run
  d <- file.path(tempfile("rtdir"), "run2"); dir.create(d, recursive = TRUE)
  writeLines(ext_lines, file.path(d, "run2.ext"))
  invisible(capture.output(suppressMessages(
    res3 <- nm_extract_run_table(d))))
  expect_equal(res3$extension, "dir")
})
